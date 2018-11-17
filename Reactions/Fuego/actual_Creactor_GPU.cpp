#include <actual_Creactor_GPU.h> 

/**********************************/
/* Global Variables */
  N_Vector y         = NULL;
  SUNLinearSolver LS = NULL;
  int NEQ       = 0;
  int NCELLS    = 0;
  int iDense_Creact = 1;
  int iJac_Creact   = 0;
  int iE_Creact     = 1;
  int iverbose      = 2;
  void *cvode_mem   = NULL;
  double *rhoe_init = NULL;
  double *rhoh_init = NULL;
  double *rhoesrc_ext = NULL;
  double *rhohsrc_ext = NULL;
  double *rYsrc = NULL;
  double *dt_save;
  int *ncells_d;
  int *nspec;
  int *flagP;

/**********************************/
/* Definitions */
int extern_cInit(const int* cvode_meth,const int* cvode_itmeth, 
		const int* cvode_iJac, const int* cvode_iE,
		const int* cvode_iDense, const int* Ncells){

	int flag, neq_tot ;
	realtype reltol, time;
	N_Vector atol;
	realtype *ratol;
	double rwrk;
	int iwrk;
	int mm, ii, nfit;

	ckindx_(&iwrk,&rwrk,&mm,&NEQ,&ii,&nfit);
        if (iverbose > 1) {
	    printf("Nb of spec is %d \n", NEQ);
	}

	iDense_Creact  = *cvode_iDense;
	iJac_Creact    = *cvode_iJac;
	iE_Creact      = *cvode_iE;
	NCELLS         = *Ncells;
        neq_tot        = (NEQ + 1) * NCELLS;

        if (iverbose > 1) {
	    printf("Ncells in one solve ? %d\n",NCELLS);
	}

	/* Definition of main vector */
	y = N_VNew_Cuda(neq_tot);
	if(check_flag((void*)y, "N_VNew_Cuda", 0)) return(1);

	/* Call CVodeCreate to create the solver memory and specify the
	 * Backward Differentiation Formula and the use of a Newton iteration */
	if ((*cvode_meth == 2) && (*cvode_itmeth == 2))
	{
	    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
	    if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);
	} else {
	    amrex::Abort("\n--> Weird inputs to CVodeCreate. Viable options are CV_BDF (=2), CV_NEWTON (=2)\n");
	}

        time = 0.0e+0;
	/* Call CVodeInit to initialize the integrator memory and specify the
	 * user's right hand side function, the inital time, and 
	 * initial dependent variable vector y. */
	flag = CVodeInit(cvode_mem, cF_RHS, time, y);
	if (check_flag(&flag, "CVodeInit", 1)) return(1);
	
	/* Definition of tolerances: one for each species */
	reltol = 1.0e-05;
        atol  = N_VNew_Cuda(neq_tot);
	ratol = N_VGetHostArrayPointer_Cuda(atol);
        for (int i=0; i<neq_tot; i++) {
            ratol[i] = 1.0e-10;
        }
	N_VCopyToDevice_Cuda(atol);
	/* Call CVodeSVtolerances to specify the scalar relative tolerance
	 * and vector absolute tolerances */
	flag = CVodeSVtolerances(cvode_mem, reltol, atol);
	if (check_flag(&flag, "CVodeSVtolerances", 1)) return(1);

	if (iDense_Creact == 1) {
	    amrex::Abort("\n--> No Direct Dense Solver on GPU\n");
	} else if (iDense_Creact == 99) {
            printf("\n--> Using an Iterative Solver \n");

            /* Create the linear solver object */
	    LS = SUNSPGMR(y, PREC_NONE, 0);
	    if(check_flag((void *)LS, "SUNDenseLinearSolver", 0)) return(1);

	    /* Set CVSpils linear solver to LS */
	    flag = CVSpilsSetLinearSolver(cvode_mem, LS);
	    if(check_flag(&flag, "CVSpilsSetLinearSolver", 1)) return(1);
	} else {
	    amrex::Abort("\n--> Direct Sparse / Iterative Solvers not yet implemented ...\n");
	}

	if (iJac_Creact == 0) {
            printf("\n--> Without Analytical J\n");
	} else {
	    amrex::Abort("\n--> With Analytical J: not yet implemented");
	}

	flag = CVodeSetMaxNumSteps(cvode_mem, 10000);
	if(check_flag(&flag, "CVodeSetMaxNumSteps", 1)) return(1);

	/* Define vectors to be used later in creact */
	// GPU stuff: might want to rethink this and put everything in userdata
	cudaMallocManaged(&ncells_d, 1*sizeof(int));
	ncells_d[0] = NCELLS; 
	cudaMallocManaged(&nspec, 1*sizeof(int));
	*nspec = NEQ;
	cudaMallocManaged(&flagP, 1*sizeof(int));
        *flagP = iE_Creact; 
	cudaMallocManaged(&dt_save, 1*sizeof(double));
	if (iE_Creact == 1) { 
	    cudaMallocManaged(&rhoe_init, NCELLS*sizeof(double));
	    cudaMallocManaged(&rhoesrc_ext, NCELLS*sizeof(double));
	} else {
	    cudaMallocManaged(&rhoh_init, NCELLS*sizeof(double));
	    cudaMallocManaged(&rhohsrc_ext, NCELLS*sizeof(double));
	}
	cudaMallocManaged(&rYsrc, (NCELLS*NEQ)*sizeof(double));

	N_VDestroy(atol);          /* Free the atol vector */

	/* Ok we're done ...*/
        if (iverbose > 1) {
	    printf(" --> DONE WITH INITIALIZATION (GPU) %d \n", iE_Creact);
	}

	return(0);
}

/* Main routine for external looping */
int actual_cReact(realtype *rY_in, realtype *rY_src_in, 
		realtype *rX_in, realtype *rX_src_in,
		realtype *P_in, realtype *dt_react, realtype *time){

	realtype time_init, time_out ;
	int flag;

	cudaError_t cuda_status = cudaSuccess;

        time_init = *time;
	time_out  = *time + *dt_react;

	/* Get Device MemCpy of in arrays */
	/* Get Device pointer of solution vector */
	realtype *yvec_d      = N_VGetDeviceArrayPointer_Cuda(y);
	// rhoY,T
	cudaMemcpy(yvec_d, rY_in, sizeof(realtype) * ((NEQ+1)*NCELLS), cudaMemcpyHostToDevice);
	// rhoY_src_ext
	cudaMemcpy(rYsrc, rY_src_in, (NEQ*NCELLS)*sizeof(double), cudaMemcpyHostToDevice);
	// rhoE/rhoH
	if (iE_Creact == 1) { 
	    cudaMemcpy(rhoe_init, rX_in, sizeof(realtype) * NCELLS, cudaMemcpyHostToDevice);
	    cudaMemcpy(rhoesrc_ext, rX_src_in, sizeof(realtype) * NCELLS, cudaMemcpyHostToDevice);
	} else {
	    cudaMemcpy(rhoh_init, rX_in, sizeof(realtype) * NCELLS, cudaMemcpyHostToDevice);
	    cudaMemcpy(rhohsrc_ext, rX_src_in, sizeof(realtype) * NCELLS, cudaMemcpyHostToDevice);
	}

	/* Call CVODE: ReInit for convergence */
	CVodeReInit(cvode_mem, time_init, y);
	flag = CVode(cvode_mem, time_out, y, &time_init, CV_NORMAL);
	if (check_flag(&flag, "CVode", 1)) return(1);

	/* Pack data to return in main routine external */
	cudaMemcpy(rY_in, yvec_d, ((NEQ+1)*NCELLS)*sizeof(realtype), cudaMemcpyDeviceToHost);
	for  (int i = 0; i < NCELLS; i++) {
            rX_in[i] = rX_in[i] + (*dt_react) * rX_src_in[i];
	}

	/* tests HP PUT THIS ON DEVICE ! */
        for (int tid = 0; tid < NCELLS; tid ++) {
	    double rhov, energy, temp;
	    double MF[NEQ];
            int * iwrk;
            double *  rwrk;
	    int  lierr;
	    rhov = 0.0;
            int offset = tid * (NEQ + 1); 
            for (int k = 0; k < NEQ; k ++) {
		rhov =  rhov + rY_in[offset + k];
	    }
            for (int k = 0; k < NEQ; k ++) {
		MF[k] = rY_in[offset + k]/rhov;
	    }
	    energy = rX_in[tid]/rhov ;
	    if (iE_Creact == 1) { 
	        get_t_given_ey_(&energy, MF, iwrk, rwrk, &temp, &lierr);
	    } else {
	        get_t_given_hy_(&energy, MF, iwrk, rwrk, &temp, &lierr);
	    }
	    rY_in[offset + NEQ] =  temp;
	}

	/* If in debug mode: print stats */
        if (iverbose > 2) {
           int ierr;
	   double htmp;
           printf("\n......cvode done:\n");
           ierr = CVodeGetLastStep(cvode_mem, &htmp);
           printf(" -- last successful step size = %4.8e \n",htmp);
           ierr = CVodeGetCurrentStep(cvode_mem, &htmp);
           printf(" -- next step to try = %4.8e \n",htmp);
           ierr = CVodeGetCurrentTime(cvode_mem, &htmp);
           printf(" -- integrated time reached = %4.8e \n",htmp);
	   long int itmp, itmp2; 
           ierr = CVodeGetNumSteps(cvode_mem, &itmp);
           printf(" -- number of time steps (nst) = %-6ld \n",itmp);
           ierr = CVodeGetNumRhsEvals(cvode_mem, &itmp);
           ierr = CVSpilsGetNumRhsEvals(cvode_mem, &itmp2);
           itmp = itmp + itmp2;
           printf(" -- number of fRHS EVAL (nfe+nfels) = %-6ld \n", itmp);
           ierr = CVDlsGetNumJacEvals(cvode_mem, &itmp);
           printf(" -- number of Jac EVAL = %-6ld \n", itmp);
	   int itmp3; 
           ierr = CVodeGetLastOrder(cvode_mem, &itmp3);
           printf(" -- method order last used = %d \n", itmp3);
           ierr = CVodeGetCurrentOrder(cvode_mem, &itmp3);
           printf(" -- method order to be used = %d \n", itmp3);
           ierr = CVodeGetNumLinSolvSetups(cvode_mem, &itmp);
	   printf(" -- number of linear solver setups %-6ld \n", itmp);
           ierr = CVodeGetNumErrTestFails(cvode_mem,  &itmp);
	   printf(" -- number of err test fails (netf) %-6ld \n", itmp);
           ierr = CVodeGetNumNonlinSolvIters(cvode_mem, &itmp);
           printf(" -- number of Newton iterations (nni) %-6ld \n", itmp);
           ierr = CVodeGetNumNonlinSolvConvFails(cvode_mem, &itmp);
	   printf(" -- number of Newton failures %-6ld \n", itmp);
           ierr = CVodeGetNumGEvals(cvode_mem, &itmp);
	   printf(" -- nge ? %-6ld \n", itmp);
	}
	return(0);
}

/* RHS routine used in CVODE */
static int cF_RHS(realtype t, N_Vector y_in, N_Vector ydot_in, 
		void *user_data){

        //std::chrono::time_point<std::chrono::system_clock> start, end;		
	//start = std::chrono::system_clock::now();
	cudaError_t cuda_status = cudaSuccess;

	/* Get Device pointers for Kernel call */
	realtype *yvec_d      = N_VGetDeviceArrayPointer_Cuda(y_in);
	realtype *ydot_d      = N_VGetDeviceArrayPointer_Cuda(ydot_in);
	// UV/HP
	cudaMemcpy(dt_save, &t, sizeof(double), cudaMemcpyHostToDevice);

        if (iE_Creact == 1) {
	    /* GPU tests */
	    unsigned block = 32;
	    unsigned grid = NCELLS/32 + 1;
	    fKernelSpec<<<grid,block>>>(dt_save, ncells_d, nspec, 
			    yvec_d, ydot_d, 
			    rhoe_init, rhoesrc_ext, rYsrc, flagP);
	    cuda_status = cudaDeviceSynchronize();
	    //std::cout << "In fun_rhs, got cudaDeviceSynchronize error of: " << cudaGetErrorString(cuda_status) << std::endl;
	    assert(cuda_status == cudaSuccess);
	} else {
	    unsigned block = 32;
	    unsigned grid = NCELLS/32 + 1;
	    fKernelSpec<<<grid,block>>>(dt_save, ncells_d, nspec,
			    yvec_d, ydot_d, 
			    rhoh_init, rhohsrc_ext, rYsrc, flagP);
	    cuda_status = cudaDeviceSynchronize();
	    //std::cout << "In fun_rhs, got cudaDeviceSynchronize error of: " << cudaGetErrorString(cuda_status) << std::endl;
	    assert(cuda_status == cudaSuccess);
	}
	//end = std::chrono::system_clock::now();
	//std::chrono::duration<double> elapsed_seconds = end - start;
	//std::cout << " RHS duration " << elapsed_seconds.count() << std::endl;
	return(0);
}

/*
 * CUDA kernels
 */
__global__ void fKernelSpec(realtype *dt, int *ncells_d, int *nspec, 
		            realtype *yvec_d, realtype *ydot_d,  
		            double *rhoX_init, double *rhoXsrc_ext, double *rYs, int *flagDoP)
{
  int tid;

  tid = blockDim.x * blockIdx.x + threadIdx.x;

  if (tid < ncells_d[0]) {
      //How can I pass the 29 here I have no idea;
      realtype massfrac[56],activity[56];
      realtype Xi[56], cXi[56];
      realtype cdot[56], molecular_weight[56];
      realtype temp, energy;
      int lierr;

      int offset = tid * (*nspec + 1); 

      /* MW CGS */
      molecularWeight_d(molecular_weight);
      /* rho */ 
      realtype rho = 0.0;
      for (int i = 0; i < *nspec; i++){
          rho = rho + yvec_d[offset + i];
      }
      /* temp */
      temp = yvec_d[offset + *nspec];
      /* Yks, C CGS*/
      for (int i = 0; i < *nspec; i++){
          massfrac[i] = yvec_d[offset + i] / rho;
	  activity[i] = yvec_d[offset + i]/(molecular_weight[i]);
      }
      /* NRG CGS */
      energy = (rhoX_init[tid] + rhoXsrc_ext[tid]*(*dt)) /rho;

      /* Fuego calls on device */
      if (*flagDoP == 1){
          get_t_given_ey_d_(&energy, massfrac, &temp, &lierr);
          ckums_d(&temp, Xi);
          ckcvms_d(&temp, cXi);
      } else {
          get_t_given_hy_d_(&energy, massfrac, &temp, &lierr);
          ckhms_d(&temp, Xi);
          ckcpms_d(&temp, cXi);
      }
      ckwc_d(&temp, activity, cdot);
      int cX = 0.0;
      for (int i = 0; i < *nspec; i++){
          cX = cX + massfrac[i] * cXi[i];
      }

      /* Fill ydot vect */
      ydot_d[offset + *nspec] = rhoXsrc_ext[tid];
      for (int i = 0; i < *nspec; i++){
          ydot_d[offset + i] = cdot[i] * molecular_weight[i] + rYs[tid * (*nspec) + i];
          ydot_d[offset + *nspec] = ydot_d[offset + *nspec]  - ydot_d[offset + i] * Xi[i];
      }
      ydot_d[offset + *nspec] = ydot_d[offset + *nspec] /(rho * cX);
  }
}

 /* Free and destroy memory */
void extern_cFree(){

	printf("In cFree\n");
	SUNLinSolFree(LS);
	N_VDestroy(y);          /* Free the y vector */
	CVodeFree(&cvode_mem);
	if (iE_Creact == 1) { 
	    cudaFree(rhoe_init);
	    cudaFree(rhoesrc_ext);
	} else {
	    cudaFree(rhoh_init);
	    cudaFree(rhohsrc_ext);
	}
	cudaFree(rYsrc);
	cudaFree(ncells_d);
	cudaFree(nspec);
}

/* Get and print some final statistics */
static void PrintFinalStats(void *cvodeMem)
{
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
  int flag;

  flag = CVodeGetNumSteps(cvodeMem, &nst);
  check_flag(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetNumRhsEvals(cvodeMem, &nfe);
  check_flag(&flag, "CVodeGetNumRhsEvals", 1);
  flag = CVodeGetNumLinSolvSetups(cvodeMem, &nsetups);
  check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
  flag = CVodeGetNumErrTestFails(cvodeMem, &netf);
  check_flag(&flag, "CVodeGetNumErrTestFails", 1);
  flag = CVodeGetNumNonlinSolvIters(cvodeMem, &nni);
  check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvodeMem, &ncfn);
  check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);

  flag = CVDlsGetNumJacEvals(cvodeMem, &nje);
  check_flag(&flag, "CVDlsGetNumJacEvals", 1);
  flag = CVDlsGetNumRhsEvals(cvodeMem, &nfeLS);
  check_flag(&flag, "CVDlsGetNumRhsEvals", 1);

  flag = CVodeGetNumGEvals(cvodeMem, &nge);
  check_flag(&flag, "CVodeGetNumGEvals", 1);

  printf("\nFinal Statistics:\n");
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n",
	 nst, nfe, nsetups, nfeLS, nje);
  printf("nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n \n",
	 nni, ncfn, netf, nge);
}

/* Check function return value...
     opt == 0 means SUNDIALS function allocates memory so check if
              returned NULL pointer
     opt == 1 means SUNDIALS function returns a flag so check if
              flag >= 0
     opt == 2 means function allocates memory so check if returned
              NULL pointer */
static int check_flag(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
              funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

  return(0);
}

/*
 * CUDA device functions
 */
/*save inv molecular weights into array */
__device__ void imolecularWeight_d(double * iwt)
{
    iwt[0] = 0.992093; /*H */
    iwt[1] = 0.062502; /*O */
    iwt[2] = 0.058798; /*OH */
    iwt[3] = 0.030297; /*HO2 */
    iwt[4] = 0.496047; /*H2 */
    iwt[5] = 0.055508; /*H2O */
    iwt[6] = 0.029399; /*H2O2 */
    iwt[7] = 0.031251; /*O2 */
    iwt[8] = 0.076810; /*CH */
    iwt[9] = 0.071291; /*CH2 */
    iwt[10] = 0.071291; /*CH2* */
    iwt[11] = 0.066511; /*CH3 */
    iwt[12] = 0.062332; /*CH4 */
    iwt[13] = 0.034461; /*HCO */
    iwt[14] = 0.033304; /*CH2O */
    iwt[15] = 0.032222; /*CH3O */
    iwt[16] = 0.032222; /*CH2OH */
    iwt[17] = 0.031209; /*CH3OH */
    iwt[18] = 0.035701; /*CO */
    iwt[19] = 0.022722; /*CO2 */
    iwt[20] = 0.039952; /*C2H */
    iwt[21] = 0.038405; /*C2H2 */
    iwt[22] = 0.036974; /*C2H3 */
    iwt[23] = 0.035645; /*C2H4 */
    iwt[24] = 0.034409; /*C2H5 */
    iwt[25] = 0.033256; /*C2H6 */
    iwt[26] = 0.024373; /*HCCO */
    iwt[27] = 0.023788; /*CH2CO */
    iwt[28] = 0.023231; /*CH3CO */
    iwt[29] = 0.023231; /*CH2CHO */
    iwt[30] = 0.022700; /*CH3CHO */
    iwt[31] = 0.025603; /*C3H3 */
    iwt[32] = 0.024959; /*pC3H4 */
    iwt[33] = 0.024959; /*aC3H4 */
    iwt[34] = 0.024347; /*aC3H5 */
    iwt[35] = 0.024347; /*CH3CCH2 */
    iwt[36] = 0.023764; /*C3H6 */
    iwt[37] = 0.023208; /*nC3H7 */
    iwt[38] = 0.023208; /*iC3H7 */
    iwt[39] = 0.017837; /*C2H3CHO */
    iwt[40] = 0.019976; /*C4H2 */
    iwt[41] = 0.019582; /*iC4H3 */
    iwt[42] = 0.019203; /*C4H4 */
    iwt[43] = 0.018838; /*iC4H5 */
    iwt[44] = 0.018838; /*C4H5-2 */
    iwt[45] = 0.018487; /*C4H6 */
    iwt[46] = 0.018487; /*C4H612 */
    iwt[47] = 0.018487; /*C4H6-2 */
    iwt[48] = 0.018149; /*C4H7 */
    iwt[49] = 0.017823; /*C4H81 */
    iwt[50] = 0.017508; /*pC4H9 */
    iwt[51] = 0.005871; /*NC12H26 */
    iwt[52] = 0.011882; /*C6H12 */
    iwt[53] = 0.012026; /*C6H11 */
    iwt[54] = 0.014258; /*C5H10 */
    iwt[55] = 0.035697; /*N2 */

    return;
}


/* Returns R, Rc, Patm */
__device__ void ckrp_d( double * ru, double * ruc, double * pa)
{
     *ru  = 8.31451e+07; 
     *ruc = 1.98721558317399615845; 
     *pa  = 1.01325e+06; 
}


/*Compute P = rhoRT/W(x) */
__device__ void ckpx_d(double * rho, double * T, double * x, double * P)
{
    double XW = 0;/* To hold mean molecular wt */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*15.999400; /*O */
    XW += x[2]*17.007370; /*OH */
    XW += x[3]*33.006770; /*HO2 */
    XW += x[4]*2.015940; /*H2 */
    XW += x[5]*18.015340; /*H2O */
    XW += x[6]*34.014740; /*H2O2 */
    XW += x[7]*31.998800; /*O2 */
    XW += x[8]*13.019120; /*CH */
    XW += x[9]*14.027090; /*CH2 */
    XW += x[10]*14.027090; /*CH2* */
    XW += x[11]*15.035060; /*CH3 */
    XW += x[12]*16.043030; /*CH4 */
    XW += x[13]*29.018520; /*HCO */
    XW += x[14]*30.026490; /*CH2O */
    XW += x[15]*31.034460; /*CH3O */
    XW += x[16]*31.034460; /*CH2OH */
    XW += x[17]*32.042430; /*CH3OH */
    XW += x[18]*28.010550; /*CO */
    XW += x[19]*44.009950; /*CO2 */
    XW += x[20]*25.030270; /*C2H */
    XW += x[21]*26.038240; /*C2H2 */
    XW += x[22]*27.046210; /*C2H3 */
    XW += x[23]*28.054180; /*C2H4 */
    XW += x[24]*29.062150; /*C2H5 */
    XW += x[25]*30.070120; /*C2H6 */
    XW += x[26]*41.029670; /*HCCO */
    XW += x[27]*42.037640; /*CH2CO */
    XW += x[28]*43.045610; /*CH3CO */
    XW += x[29]*43.045610; /*CH2CHO */
    XW += x[30]*44.053580; /*CH3CHO */
    XW += x[31]*39.057360; /*C3H3 */
    XW += x[32]*40.065330; /*pC3H4 */
    XW += x[33]*40.065330; /*aC3H4 */
    XW += x[34]*41.073300; /*aC3H5 */
    XW += x[35]*41.073300; /*CH3CCH2 */
    XW += x[36]*42.081270; /*C3H6 */
    XW += x[37]*43.089240; /*nC3H7 */
    XW += x[38]*43.089240; /*iC3H7 */
    XW += x[39]*56.064730; /*C2H3CHO */
    XW += x[40]*50.060540; /*C4H2 */
    XW += x[41]*51.068510; /*iC4H3 */
    XW += x[42]*52.076480; /*C4H4 */
    XW += x[43]*53.084450; /*iC4H5 */
    XW += x[44]*53.084450; /*C4H5-2 */
    XW += x[45]*54.092420; /*C4H6 */
    XW += x[46]*54.092420; /*C4H612 */
    XW += x[47]*54.092420; /*C4H6-2 */
    XW += x[48]*55.100390; /*C4H7 */
    XW += x[49]*56.108360; /*C4H81 */
    XW += x[50]*57.116330; /*pC4H9 */
    XW += x[51]*170.341020; /*NC12H26 */
    XW += x[52]*84.162540; /*C6H12 */
    XW += x[53]*83.154570; /*C6H11 */
    XW += x[54]*70.135450; /*C5H10 */
    XW += x[55]*28.013400; /*N2 */
    *P = *rho * 8.31451e+07 * (*T) / XW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(y) */
__device__ void ckpy_d(double * rho, double * T, double * y_wk, double * P)
{
    double imw[56];/* inv molecular weight array */
    imolecularWeight_d(imw);
    double YOW = 0;/* for computing mean MW */
    YOW += y_wk[0]*imw[0]; /*H */
    YOW += y_wk[1]*imw[1]; /*O */
    YOW += y_wk[2]*imw[2]; /*OH */
    YOW += y_wk[3]*imw[3]; /*HO2 */
    YOW += y_wk[4]*imw[4]; /*H2 */
    YOW += y_wk[5]*imw[5]; /*H2O */
    YOW += y_wk[6]*imw[6]; /*H2O2 */
    YOW += y_wk[7]*imw[7]; /*O2 */
    YOW += y_wk[8]*imw[8]; /*CH */
    YOW += y_wk[9]*imw[9]; /*CH2 */
    YOW += y_wk[10]*imw[10]; /*CH2* */
    YOW += y_wk[11]*imw[11]; /*CH3 */
    YOW += y_wk[12]*imw[12]; /*CH4 */
    YOW += y_wk[13]*imw[13]; /*HCO */
    YOW += y_wk[14]*imw[14]; /*CH2O */
    YOW += y_wk[15]*imw[15]; /*CH3O */
    YOW += y_wk[16]*imw[16]; /*CH2OH */
    YOW += y_wk[17]*imw[17]; /*CH3OH */
    YOW += y_wk[18]*imw[18]; /*CO */
    YOW += y_wk[19]*imw[19]; /*CO2 */
    YOW += y_wk[20]*imw[20]; /*C2H */
    YOW += y_wk[21]*imw[21]; /*C2H2 */
    YOW += y_wk[22]*imw[22]; /*C2H3 */
    YOW += y_wk[23]*imw[23]; /*C2H4 */
    YOW += y_wk[24]*imw[24]; /*C2H5 */
    YOW += y_wk[25]*imw[25]; /*C2H6 */
    YOW += y_wk[26]*imw[26]; /*HCCO */
    YOW += y_wk[27]*imw[27]; /*CH2CO */
    YOW += y_wk[28]*imw[28]; /*CH3CO */
    YOW += y_wk[29]*imw[29]; /*CH2CHO */
    YOW += y_wk[30]*imw[30]; /*CH3CHO */
    YOW += y_wk[31]*imw[31]; /*C3H3 */
    YOW += y_wk[32]*imw[32]; /*pC3H4 */
    YOW += y_wk[33]*imw[33]; /*aC3H4 */
    YOW += y_wk[34]*imw[34]; /*aC3H5 */
    YOW += y_wk[35]*imw[35]; /*CH3CCH2 */
    YOW += y_wk[36]*imw[36]; /*C3H6 */
    YOW += y_wk[37]*imw[37]; /*nC3H7 */
    YOW += y_wk[38]*imw[38]; /*iC3H7 */
    YOW += y_wk[39]*imw[39]; /*C2H3CHO */
    YOW += y_wk[40]*imw[40]; /*C4H2 */
    YOW += y_wk[41]*imw[41]; /*iC4H3 */
    YOW += y_wk[42]*imw[42]; /*C4H4 */
    YOW += y_wk[43]*imw[43]; /*iC4H5 */
    YOW += y_wk[44]*imw[44]; /*C4H5-2 */
    YOW += y_wk[45]*imw[45]; /*C4H6 */
    YOW += y_wk[46]*imw[46]; /*C4H612 */
    YOW += y_wk[47]*imw[47]; /*C4H6-2 */
    YOW += y_wk[48]*imw[48]; /*C4H7 */
    YOW += y_wk[49]*imw[49]; /*C4H81 */
    YOW += y_wk[50]*imw[50]; /*pC4H9 */
    YOW += y_wk[51]*imw[51]; /*NC12H26 */
    YOW += y_wk[52]*imw[52]; /*C6H12 */
    YOW += y_wk[53]*imw[53]; /*C6H11 */
    YOW += y_wk[54]*imw[54]; /*C5H10 */
    YOW += y_wk[55]*imw[55]; /*N2 */
    *P = *rho * 8.31451e+07 * (*T) * YOW; /*P = rho*R*T/W */

    return;
}


/*Compute rho = P*W(y)/RT */
__device__ void ckrhoy_d(double * P, double * T, double * y_wk, double * rho)
{
    double YOW = 0;
    double imw[56];/* inv molecular weight array */
    imolecularWeight_d(imw);
    double tmp[56];

    for (int i = 0; i < 56; i++)
    {
        tmp[i] = y_wk[i]*imw[i];
    }
    for (int i = 0; i < 56; i++)
    {
        YOW += tmp[i];
    }

    *rho = *P / (8.31451e+07 * (*T) * YOW);/*rho = P*W/(R*T) */
    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
__device__ void ckytcr_d(double * rho, double * T, double * y_wk, double * c)
{
    double imw[56];/* inv molecular weight array */
    imolecularWeight_d(imw);
    for (int i = 0; i < 56; i++)
    {
        c[i] = (*rho)  * y_wk[i] * imw[i];
    }
}


/*Returns the specific heats at constant volume */
/*in mass units (Eq. 29) */
__device__ void ckcvms_d(double * T, double * cvms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cv_R_d(cvms, tc);
    /*multiply by R/molecularweight */
    cvms[0] *= 8.248767324424338e+07; /*H */
    cvms[1] *= 5.196763628636074e+06; /*O */
    cvms[2] *= 4.888768810227566e+06; /*OH */
    cvms[3] *= 2.519031701678171e+06; /*HO2 */
    cvms[4] *= 4.124383662212169e+07; /*H2 */
    cvms[5] *= 4.615239012974499e+06; /*H2O */
    cvms[6] *= 2.444384405113783e+06; /*H2O2 */
    cvms[7] *= 2.598381814318037e+06; /*O2 */
    cvms[8] *= 6.386384025955671e+06; /*CH */
    cvms[9] *= 5.927466067445207e+06; /*CH2 */
    cvms[10] *= 5.927466067445207e+06; /*CH2* */
    cvms[11] *= 5.530081023953346e+06; /*CH3 */
    cvms[12] *= 5.182630712527496e+06; /*CH4 */
    cvms[13] *= 2.865242610581105e+06; /*HCO */
    cvms[14] *= 2.769058254894261e+06; /*CH2O */
    cvms[15] *= 2.679121853578248e+06; /*CH3O */
    cvms[16] *= 2.679121853578248e+06; /*CH2OH */
    cvms[17] *= 2.594843774332970e+06; /*CH3OH */
    cvms[18] *= 2.968349425484326e+06; /*CO */
    cvms[19] *= 1.889234139098090e+06; /*CO2 */
    cvms[20] *= 3.321781986370902e+06; /*C2H */
    cvms[21] *= 3.193192012977835e+06; /*C2H2 */
    cvms[22] *= 3.074186734481467e+06; /*C2H3 */
    cvms[23] *= 2.963733033722604e+06; /*C2H4 */
    cvms[24] *= 2.860941121011349e+06; /*C2H5 */
    cvms[25] *= 2.765040511976673e+06; /*C2H6 */
    cvms[26] *= 2.026462801187531e+06; /*HCCO */
    cvms[27] *= 1.977872687429646e+06; /*CH2CO */
    cvms[28] *= 1.931558177477331e+06; /*CH3CO */
    cvms[29] *= 1.931558177477332e+06; /*CH2CHO */
    cvms[30] *= 1.887363070152301e+06; /*CH3CHO */
    cvms[31] *= 2.128794675318557e+06; /*C3H3 */
    cvms[32] *= 2.075238117344846e+06; /*pC3H4 */
    cvms[33] *= 2.075238117344846e+06; /*aC3H4 */
    cvms[34] *= 2.024310196648431e+06; /*aC3H5 */
    cvms[35] *= 2.024310196648431e+06; /*CH3CCH2 */
    cvms[36] *= 1.975822022481736e+06; /*C3H6 */
    cvms[37] *= 1.929602378691293e+06; /*nC3H7 */
    cvms[38] *= 1.929602378691293e+06; /*iC3H7 */
    cvms[39] *= 1.483019716673923e+06; /*C2H3CHO */
    cvms[40] *= 1.660890993185451e+06; /*C4H2 */
    cvms[41] *= 1.628108985361038e+06; /*iC4H3 */
    cvms[42] *= 1.596596006488918e+06; /*C4H4 */
    cvms[43] *= 1.566279767427184e+06; /*iC4H5 */
    cvms[44] *= 1.566279767427184e+06; /*C4H5-2 */
    cvms[45] *= 1.537093367240733e+06; /*C4H6 */
    cvms[46] *= 1.537093367240733e+06; /*C4H612 */
    cvms[47] *= 1.537093367240733e+06; /*C4H6-2 */
    cvms[48] *= 1.508974800359852e+06; /*C4H7 */
    cvms[49] *= 1.481866516861302e+06; /*C4H81 */
    cvms[50] *= 1.455715029309481e+06; /*pC4H9 */
    cvms[51] *= 4.881096755203180e+05; /*NC12H26 */
    cvms[52] *= 9.879110112408679e+05; /*C6H12 */
    cvms[53] *= 9.998861156999548e+05; /*C6H11 */
    cvms[54] *= 1.185493213489041e+06; /*C5H10 */
    cvms[55] *= 2.968047434442088e+06; /*N2 */
}


/*Returns the specific heats at constant pressure */
/*in mass units (Eq. 26) */
__device__ void ckcpms_d(double * T, double * cpms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R_d(cpms, tc);
    /*multiply by R/molecularweight */
    cpms[0] *= 8.248767324424338e+07; /*H */
    cpms[1] *= 5.196763628636074e+06; /*O */
    cpms[2] *= 4.888768810227566e+06; /*OH */
    cpms[3] *= 2.519031701678171e+06; /*HO2 */
    cpms[4] *= 4.124383662212169e+07; /*H2 */
    cpms[5] *= 4.615239012974499e+06; /*H2O */
    cpms[6] *= 2.444384405113783e+06; /*H2O2 */
    cpms[7] *= 2.598381814318037e+06; /*O2 */
    cpms[8] *= 6.386384025955671e+06; /*CH */
    cpms[9] *= 5.927466067445207e+06; /*CH2 */
    cpms[10] *= 5.927466067445207e+06; /*CH2* */
    cpms[11] *= 5.530081023953346e+06; /*CH3 */
    cpms[12] *= 5.182630712527496e+06; /*CH4 */
    cpms[13] *= 2.865242610581105e+06; /*HCO */
    cpms[14] *= 2.769058254894261e+06; /*CH2O */
    cpms[15] *= 2.679121853578248e+06; /*CH3O */
    cpms[16] *= 2.679121853578248e+06; /*CH2OH */
    cpms[17] *= 2.594843774332970e+06; /*CH3OH */
    cpms[18] *= 2.968349425484326e+06; /*CO */
    cpms[19] *= 1.889234139098090e+06; /*CO2 */
    cpms[20] *= 3.321781986370902e+06; /*C2H */
    cpms[21] *= 3.193192012977835e+06; /*C2H2 */
    cpms[22] *= 3.074186734481467e+06; /*C2H3 */
    cpms[23] *= 2.963733033722604e+06; /*C2H4 */
    cpms[24] *= 2.860941121011349e+06; /*C2H5 */
    cpms[25] *= 2.765040511976673e+06; /*C2H6 */
    cpms[26] *= 2.026462801187531e+06; /*HCCO */
    cpms[27] *= 1.977872687429646e+06; /*CH2CO */
    cpms[28] *= 1.931558177477331e+06; /*CH3CO */
    cpms[29] *= 1.931558177477332e+06; /*CH2CHO */
    cpms[30] *= 1.887363070152301e+06; /*CH3CHO */
    cpms[31] *= 2.128794675318557e+06; /*C3H3 */
    cpms[32] *= 2.075238117344846e+06; /*pC3H4 */
    cpms[33] *= 2.075238117344846e+06; /*aC3H4 */
    cpms[34] *= 2.024310196648431e+06; /*aC3H5 */
    cpms[35] *= 2.024310196648431e+06; /*CH3CCH2 */
    cpms[36] *= 1.975822022481736e+06; /*C3H6 */
    cpms[37] *= 1.929602378691293e+06; /*nC3H7 */
    cpms[38] *= 1.929602378691293e+06; /*iC3H7 */
    cpms[39] *= 1.483019716673923e+06; /*C2H3CHO */
    cpms[40] *= 1.660890993185451e+06; /*C4H2 */
    cpms[41] *= 1.628108985361038e+06; /*iC4H3 */
    cpms[42] *= 1.596596006488918e+06; /*C4H4 */
    cpms[43] *= 1.566279767427184e+06; /*iC4H5 */
    cpms[44] *= 1.566279767427184e+06; /*C4H5-2 */
    cpms[45] *= 1.537093367240733e+06; /*C4H6 */
    cpms[46] *= 1.537093367240733e+06; /*C4H612 */
    cpms[47] *= 1.537093367240733e+06; /*C4H6-2 */
    cpms[48] *= 1.508974800359852e+06; /*C4H7 */
    cpms[49] *= 1.481866516861302e+06; /*C4H81 */
    cpms[50] *= 1.455715029309481e+06; /*pC4H9 */
    cpms[51] *= 4.881096755203180e+05; /*NC12H26 */
    cpms[52] *= 9.879110112408679e+05; /*C6H12 */
    cpms[53] *= 9.998861156999548e+05; /*C6H11 */
    cpms[54] *= 1.185493213489041e+06; /*C5H10 */
    cpms[55] *= 2.968047434442088e+06; /*N2 */
}


/*Returns internal energy in mass units (Eq 30.) */
__device__ void ckums_d(double * T, double * ums)
{
    double imw[56];/* inv molecular weight array */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    imolecularWeight_d(imw);
    speciesInternalEnergy_d(ums, tc);
    for (int i = 0; i < 56; i++)
    {
        ums[i] *= RT*imw[i];
    }
}


/*Returns enthalpy in mass units (Eq 27.) */
__device__ void ckhms_d(double * T, double * hms)
{
    double imw[56];/* inv molecular weight array */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    imolecularWeight_d(imw);
    speciesEnthalpy_d(hms, tc);
    for (int i = 0; i < 56; i++)
    {
        hms[i] *= RT*imw[i];
    }
}


/*Returns the mean specific heat at CP (Eq. 34) */
__device__ void ckcpbs_d(double * T, double * y_wk, double * cpbs)
{
    double result = 0; 
    double imw[56];/* inv molecular weight array */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[56], tresult[56]; /* temporary storage */
    imolecularWeight_d(imw);
    cp_R_d(cpor, tc);
    for (int i = 0; i < 56; i++)
    {
        tresult[i] = cpor[i]*y_wk[i]*imw[i];

    }
    for (int i = 0; i < 56; i++)
    {
        result += tresult[i];
    }

    *cpbs = result * 8.31451e+07;
}


/*Returns the mean specific heat at CV (Eq. 36) */
__device__ void ckcvbs_d(double * T, double * y_wk, double * cvbs)
{
    double result = 0; 
    double imw[56];/* inv molecular weight array */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[56]; /* temporary storage */
    imolecularWeight_d(imw);
    cv_R_d(cvor, tc);
    /*multiply by y/molecularweight */
    result += cvor[0]*y_wk[0]*imw[0]; /*H */
    result += cvor[1]*y_wk[1]*imw[1]; /*O */
    result += cvor[2]*y_wk[2]*imw[2]; /*OH */
    result += cvor[3]*y_wk[3]*imw[3]; /*HO2 */
    result += cvor[4]*y_wk[4]*imw[4]; /*H2 */
    result += cvor[5]*y_wk[5]*imw[5]; /*H2O */
    result += cvor[6]*y_wk[6]*imw[6]; /*H2O2 */
    result += cvor[7]*y_wk[7]*imw[7]; /*O2 */
    result += cvor[8]*y_wk[8]*imw[8]; /*CH */
    result += cvor[9]*y_wk[9]*imw[9]; /*CH2 */
    result += cvor[10]*y_wk[10]*imw[10]; /*CH2* */
    result += cvor[11]*y_wk[11]*imw[11]; /*CH3 */
    result += cvor[12]*y_wk[12]*imw[12]; /*CH4 */
    result += cvor[13]*y_wk[13]*imw[13]; /*HCO */
    result += cvor[14]*y_wk[14]*imw[14]; /*CH2O */
    result += cvor[15]*y_wk[15]*imw[15]; /*CH3O */
    result += cvor[16]*y_wk[16]*imw[16]; /*CH2OH */
    result += cvor[17]*y_wk[17]*imw[17]; /*CH3OH */
    result += cvor[18]*y_wk[18]*imw[18]; /*CO */
    result += cvor[19]*y_wk[19]*imw[19]; /*CO2 */
    result += cvor[20]*y_wk[20]*imw[20]; /*C2H */
    result += cvor[21]*y_wk[21]*imw[21]; /*C2H2 */
    result += cvor[22]*y_wk[22]*imw[22]; /*C2H3 */
    result += cvor[23]*y_wk[23]*imw[23]; /*C2H4 */
    result += cvor[24]*y_wk[24]*imw[24]; /*C2H5 */
    result += cvor[25]*y_wk[25]*imw[25]; /*C2H6 */
    result += cvor[26]*y_wk[26]*imw[26]; /*HCCO */
    result += cvor[27]*y_wk[27]*imw[27]; /*CH2CO */
    result += cvor[28]*y_wk[28]*imw[28]; /*CH3CO */
    result += cvor[29]*y_wk[29]*imw[29]; /*CH2CHO */
    result += cvor[30]*y_wk[30]*imw[30]; /*CH3CHO */
    result += cvor[31]*y_wk[31]*imw[31]; /*C3H3 */
    result += cvor[32]*y_wk[32]*imw[32]; /*pC3H4 */
    result += cvor[33]*y_wk[33]*imw[33]; /*aC3H4 */
    result += cvor[34]*y_wk[34]*imw[34]; /*aC3H5 */
    result += cvor[35]*y_wk[35]*imw[35]; /*CH3CCH2 */
    result += cvor[36]*y_wk[36]*imw[36]; /*C3H6 */
    result += cvor[37]*y_wk[37]*imw[37]; /*nC3H7 */
    result += cvor[38]*y_wk[38]*imw[38]; /*iC3H7 */
    result += cvor[39]*y_wk[39]*imw[39]; /*C2H3CHO */
    result += cvor[40]*y_wk[40]*imw[40]; /*C4H2 */
    result += cvor[41]*y_wk[41]*imw[41]; /*iC4H3 */
    result += cvor[42]*y_wk[42]*imw[42]; /*C4H4 */
    result += cvor[43]*y_wk[43]*imw[43]; /*iC4H5 */
    result += cvor[44]*y_wk[44]*imw[44]; /*C4H5-2 */
    result += cvor[45]*y_wk[45]*imw[45]; /*C4H6 */
    result += cvor[46]*y_wk[46]*imw[46]; /*C4H612 */
    result += cvor[47]*y_wk[47]*imw[47]; /*C4H6-2 */
    result += cvor[48]*y_wk[48]*imw[48]; /*C4H7 */
    result += cvor[49]*y_wk[49]*imw[49]; /*C4H81 */
    result += cvor[50]*y_wk[50]*imw[50]; /*pC4H9 */
    result += cvor[51]*y_wk[51]*imw[51]; /*NC12H26 */
    result += cvor[52]*y_wk[52]*imw[52]; /*C6H12 */
    result += cvor[53]*y_wk[53]*imw[53]; /*C6H11 */
    result += cvor[54]*y_wk[54]*imw[54]; /*C5H10 */
    result += cvor[55]*y_wk[55]*imw[55]; /*N2 */

    *cvbs = result * 8.31451e+07;
}


/*Returns mean enthalpy of mixture in mass units */
__device__ void ckhbms_d(double * T, double * y_wk, double * hbms)
{
    double result = 0;
    double imw[56];/* inv molecular weight array */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[56], tmp[56]; /* temporary storage */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy_d(hml, tc);
    imolecularWeight_d(imw);
    int id;
    for (id = 0; id < 56; ++id) {
        tmp[id] = y_wk[id]*hml[id]*imw[id];
    }
    for (id = 0; id < 56; ++id) {
        result += tmp[id];
    }

    *hbms = result * RT;
}


/*get mean internal energy in mass units */
__device__ void ckubms_d(double * T, double * y_wk, double * ubms)
{
    double result = 0;
    double imw[56];/* inv molecular weight array */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double ums[56]; /* temporary energy array */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy_d(ums, tc);
    imolecularWeight_d(imw);
    /*perform dot product + scaling by wt */
    result += y_wk[0]*ums[0]*imw[0]; /*H */
    result += y_wk[1]*ums[1]*imw[1]; /*O */
    result += y_wk[2]*ums[2]*imw[2]; /*OH */
    result += y_wk[3]*ums[3]*imw[3]; /*HO2 */
    result += y_wk[4]*ums[4]*imw[4]; /*H2 */
    result += y_wk[5]*ums[5]*imw[5]; /*H2O */
    result += y_wk[6]*ums[6]*imw[6]; /*H2O2 */
    result += y_wk[7]*ums[7]*imw[7]; /*O2 */
    result += y_wk[8]*ums[8]*imw[8]; /*CH */
    result += y_wk[9]*ums[9]*imw[9]; /*CH2 */
    result += y_wk[10]*ums[10]*imw[10]; /*CH2* */
    result += y_wk[11]*ums[11]*imw[11]; /*CH3 */
    result += y_wk[12]*ums[12]*imw[12]; /*CH4 */
    result += y_wk[13]*ums[13]*imw[13]; /*HCO */
    result += y_wk[14]*ums[14]*imw[14]; /*CH2O */
    result += y_wk[15]*ums[15]*imw[15]; /*CH3O */
    result += y_wk[16]*ums[16]*imw[16]; /*CH2OH */
    result += y_wk[17]*ums[17]*imw[17]; /*CH3OH */
    result += y_wk[18]*ums[18]*imw[18]; /*CO */
    result += y_wk[19]*ums[19]*imw[19]; /*CO2 */
    result += y_wk[20]*ums[20]*imw[20]; /*C2H */
    result += y_wk[21]*ums[21]*imw[21]; /*C2H2 */
    result += y_wk[22]*ums[22]*imw[22]; /*C2H3 */
    result += y_wk[23]*ums[23]*imw[23]; /*C2H4 */
    result += y_wk[24]*ums[24]*imw[24]; /*C2H5 */
    result += y_wk[25]*ums[25]*imw[25]; /*C2H6 */
    result += y_wk[26]*ums[26]*imw[26]; /*HCCO */
    result += y_wk[27]*ums[27]*imw[27]; /*CH2CO */
    result += y_wk[28]*ums[28]*imw[28]; /*CH3CO */
    result += y_wk[29]*ums[29]*imw[29]; /*CH2CHO */
    result += y_wk[30]*ums[30]*imw[30]; /*CH3CHO */
    result += y_wk[31]*ums[31]*imw[31]; /*C3H3 */
    result += y_wk[32]*ums[32]*imw[32]; /*pC3H4 */
    result += y_wk[33]*ums[33]*imw[33]; /*aC3H4 */
    result += y_wk[34]*ums[34]*imw[34]; /*aC3H5 */
    result += y_wk[35]*ums[35]*imw[35]; /*CH3CCH2 */
    result += y_wk[36]*ums[36]*imw[36]; /*C3H6 */
    result += y_wk[37]*ums[37]*imw[37]; /*nC3H7 */
    result += y_wk[38]*ums[38]*imw[38]; /*iC3H7 */
    result += y_wk[39]*ums[39]*imw[39]; /*C2H3CHO */
    result += y_wk[40]*ums[40]*imw[40]; /*C4H2 */
    result += y_wk[41]*ums[41]*imw[41]; /*iC4H3 */
    result += y_wk[42]*ums[42]*imw[42]; /*C4H4 */
    result += y_wk[43]*ums[43]*imw[43]; /*iC4H5 */
    result += y_wk[44]*ums[44]*imw[44]; /*C4H5-2 */
    result += y_wk[45]*ums[45]*imw[45]; /*C4H6 */
    result += y_wk[46]*ums[46]*imw[46]; /*C4H612 */
    result += y_wk[47]*ums[47]*imw[47]; /*C4H6-2 */
    result += y_wk[48]*ums[48]*imw[48]; /*C4H7 */
    result += y_wk[49]*ums[49]*imw[49]; /*C4H81 */
    result += y_wk[50]*ums[50]*imw[50]; /*pC4H9 */
    result += y_wk[51]*ums[51]*imw[51]; /*NC12H26 */
    result += y_wk[52]*ums[52]*imw[52]; /*C6H12 */
    result += y_wk[53]*ums[53]*imw[53]; /*C6H11 */
    result += y_wk[54]*ums[54]*imw[54]; /*C5H10 */
    result += y_wk[55]*ums[55]*imw[55]; /*N2 */

    *ubms = result * RT;
}


/*compute the production rate for each species */
__device__ void ckwc_d(double * T, double * C, double * wdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 56; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    productionRate_d(wdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 56; ++id) {
        C[id] *= 1.0e-6;
        wdot[id] *= 1.0e-6;
    }
}

/*compute the production rate for each species */
__device__ void productionRate_d(double * wdot, double * sc, double T)
{
    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];

    double qdot, q_f[289], q_r[289];
    comp_qfqr_d(q_f, q_r, sc, tc, invT);

    for (int i = 0; i < 56; ++i) {
        wdot[i] = 0.0;
    }

    qdot = q_f[0]-q_r[0];
    wdot[0] -= qdot;
    wdot[3] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[1]-q_r[1];
    wdot[2] -= 2 * qdot;
    wdot[6] += qdot;

    qdot = q_f[2]-q_r[2];
    wdot[9] -= qdot;
    wdot[18] -= qdot;
    wdot[27] += qdot;

    qdot = q_f[3]-q_r[3];
    wdot[5] -= qdot;
    wdot[10] -= qdot;
    wdot[17] += qdot;

    qdot = q_f[4]-q_r[4];
    wdot[0] -= qdot;
    wdot[14] -= qdot;
    wdot[16] += qdot;

    qdot = q_f[5]-q_r[5];
    wdot[0] -= qdot;
    wdot[14] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[6]-q_r[6];
    wdot[0] -= qdot;
    wdot[11] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[7]-q_r[7];
    wdot[2] -= qdot;
    wdot[11] -= qdot;
    wdot[17] += qdot;

    qdot = q_f[8]-q_r[8];
    wdot[11] -= 2 * qdot;
    wdot[25] += qdot;

    qdot = q_f[9]-q_r[9];
    wdot[0] += qdot;
    wdot[21] += qdot;
    wdot[22] -= qdot;

    qdot = q_f[10]-q_r[10];
    wdot[0] -= qdot;
    wdot[27] -= qdot;
    wdot[29] += qdot;

    qdot = q_f[11]-q_r[11];
    wdot[0] -= qdot;
    wdot[22] -= qdot;
    wdot[23] += qdot;

    qdot = q_f[12]-q_r[12];
    wdot[11] -= qdot;
    wdot[22] -= qdot;
    wdot[36] += qdot;

    qdot = q_f[13]-q_r[13];
    wdot[11] -= qdot;
    wdot[18] -= qdot;
    wdot[28] += qdot;

    qdot = q_f[14]-q_r[14];
    wdot[11] -= qdot;
    wdot[13] -= qdot;
    wdot[30] += qdot;

    qdot = q_f[15]-q_r[15];
    wdot[4] += qdot;
    wdot[21] += qdot;
    wdot[23] -= qdot;

    qdot = q_f[16]-q_r[16];
    wdot[0] -= qdot;
    wdot[23] -= qdot;
    wdot[24] += qdot;

    qdot = q_f[17]-q_r[17];
    wdot[0] -= qdot;
    wdot[24] -= qdot;
    wdot[25] += qdot;

    qdot = q_f[18]-q_r[18];
    wdot[22] -= qdot;
    wdot[24] -= qdot;
    wdot[49] += qdot;

    qdot = q_f[19]-q_r[19];
    wdot[0] -= qdot;
    wdot[34] -= qdot;
    wdot[36] += qdot;

    qdot = q_f[20]-q_r[20];
    wdot[11] -= qdot;
    wdot[34] -= qdot;
    wdot[49] += qdot;

    qdot = q_f[21]-q_r[21];
    wdot[0] -= qdot;
    wdot[36] -= qdot;
    wdot[37] += qdot;

    qdot = q_f[22]-q_r[22];
    wdot[0] -= qdot;
    wdot[36] -= qdot;
    wdot[38] += qdot;

    qdot = q_f[23]-q_r[23];
    wdot[0] -= qdot;
    wdot[36] += qdot;
    wdot[37] += qdot;
    wdot[52] -= qdot;

    qdot = q_f[24]-q_r[24];
    wdot[0] -= qdot;
    wdot[24] += qdot;
    wdot[36] += qdot;
    wdot[54] -= qdot;

    qdot = q_f[25]-q_r[25];
    wdot[1] -= qdot;
    wdot[18] -= qdot;
    wdot[19] += qdot;

    qdot = q_f[26]-q_r[26];
    wdot[0] -= 2 * qdot;
    wdot[4] += qdot;

    qdot = q_f[27]-q_r[27];
    wdot[0] -= qdot;
    wdot[2] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[28]-q_r[28];
    wdot[0] -= qdot;
    wdot[1] -= qdot;
    wdot[2] += qdot;

    qdot = q_f[29]-q_r[29];
    wdot[1] -= 2 * qdot;
    wdot[7] += qdot;

    qdot = q_f[30]-q_r[30];
    wdot[0] += qdot;
    wdot[13] -= qdot;
    wdot[18] += qdot;

    qdot = q_f[31]-q_r[31];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[2] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[32]-q_r[32];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[4] -= qdot;

    qdot = q_f[33]-q_r[33];
    wdot[0] += qdot;
    wdot[2] -= qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[34]-q_r[34];
    wdot[1] += qdot;
    wdot[2] -= 2 * qdot;
    wdot[5] += qdot;

    qdot = q_f[35]-q_r[35];
    wdot[0] -= 2 * qdot;
    wdot[4] += qdot;
    wdot[5] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[36]-q_r[36];
    wdot[0] += qdot;
    wdot[3] += qdot;
    wdot[4] -= qdot;
    wdot[7] -= qdot;

    qdot = q_f[37]-q_r[37];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[3] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[38]-q_r[38];
    wdot[0] -= qdot;
    wdot[2] += 2 * qdot;
    wdot[3] -= qdot;

    qdot = q_f[39]-q_r[39];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[7] += qdot;

    qdot = q_f[40]-q_r[40];
    wdot[2] -= qdot;
    wdot[3] -= qdot;
    wdot[5] += qdot;
    wdot[7] += qdot;

    qdot = q_f[41]-q_r[41];
    wdot[2] -= qdot;
    wdot[3] -= qdot;
    wdot[5] += qdot;
    wdot[7] += qdot;

    qdot = q_f[42]-q_r[42];
    wdot[3] -= 2 * qdot;
    wdot[6] += qdot;
    wdot[7] += qdot;

    qdot = q_f[43]-q_r[43];
    wdot[3] -= 2 * qdot;
    wdot[6] += qdot;
    wdot[7] += qdot;

    qdot = q_f[44]-q_r[44];
    wdot[0] -= qdot;
    wdot[3] += qdot;
    wdot[4] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[45]-q_r[45];
    wdot[0] -= qdot;
    wdot[2] += qdot;
    wdot[5] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[46]-q_r[46];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[3] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[47]-q_r[47];
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[5] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[48]-q_r[48];
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[5] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[49]-q_r[49];
    wdot[0] += qdot;
    wdot[2] -= qdot;
    wdot[18] -= qdot;
    wdot[19] += qdot;

    qdot = q_f[50]-q_r[50];
    wdot[0] += qdot;
    wdot[2] -= qdot;
    wdot[18] -= qdot;
    wdot[19] += qdot;

    qdot = q_f[51]-q_r[51];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[18] -= qdot;
    wdot[19] += qdot;

    qdot = q_f[52]-q_r[52];
    wdot[0] -= qdot;
    wdot[4] += qdot;
    wdot[13] -= qdot;
    wdot[18] += qdot;

    qdot = q_f[53]-q_r[53];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[13] -= qdot;
    wdot[18] += qdot;

    qdot = q_f[54]-q_r[54];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[13] -= qdot;
    wdot[19] += qdot;

    qdot = q_f[55]-q_r[55];
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[13] -= qdot;
    wdot[18] += qdot;

    qdot = q_f[56]-q_r[56];
    wdot[0] += qdot;
    wdot[5] -= qdot;
    wdot[5] += qdot;
    wdot[13] -= qdot;
    wdot[18] += qdot;

    qdot = q_f[57]-q_r[57];
    wdot[3] += qdot;
    wdot[7] -= qdot;
    wdot[13] -= qdot;
    wdot[18] += qdot;

    qdot = q_f[58]-q_r[58];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[8] -= qdot;
    wdot[18] += qdot;

    qdot = q_f[59]-q_r[59];
    wdot[0] += qdot;
    wdot[2] -= qdot;
    wdot[8] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[60]-q_r[60];
    wdot[0] += qdot;
    wdot[4] -= qdot;
    wdot[8] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[61]-q_r[61];
    wdot[0] += qdot;
    wdot[5] -= qdot;
    wdot[8] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[62]-q_r[62];
    wdot[1] += qdot;
    wdot[7] -= qdot;
    wdot[8] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[63]-q_r[63];
    wdot[8] -= qdot;
    wdot[13] += qdot;
    wdot[18] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[64]-q_r[64];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[9] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[65]-q_r[65];
    wdot[0] += qdot;
    wdot[2] -= qdot;
    wdot[9] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[66]-q_r[66];
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[8] += qdot;
    wdot[9] -= qdot;

    qdot = q_f[67]-q_r[67];
    wdot[0] += qdot;
    wdot[4] -= qdot;
    wdot[9] -= qdot;
    wdot[11] += qdot;

    qdot = q_f[68]-q_r[68];
    wdot[2] += qdot;
    wdot[7] -= qdot;
    wdot[9] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[69]-q_r[69];
    wdot[0] += 2 * qdot;
    wdot[7] -= qdot;
    wdot[9] -= qdot;
    wdot[19] += qdot;

    qdot = q_f[70]-q_r[70];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[9] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[71]-q_r[71];
    wdot[9] += qdot;
    wdot[10] -= qdot;
    wdot[55] -= qdot;
    wdot[55] += qdot;

    qdot = q_f[72]-q_r[72];
    wdot[0] -= qdot;
    wdot[4] += qdot;
    wdot[8] += qdot;
    wdot[10] -= qdot;

    qdot = q_f[73]-q_r[73];
    wdot[0] += qdot;
    wdot[2] -= qdot;
    wdot[10] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[74]-q_r[74];
    wdot[0] += qdot;
    wdot[4] -= qdot;
    wdot[10] -= qdot;
    wdot[11] += qdot;

    qdot = q_f[75]-q_r[75];
    wdot[0] += qdot;
    wdot[2] += qdot;
    wdot[7] -= qdot;
    wdot[10] -= qdot;
    wdot[18] += qdot;

    qdot = q_f[76]-q_r[76];
    wdot[5] += qdot;
    wdot[7] -= qdot;
    wdot[10] -= qdot;
    wdot[18] += qdot;

    qdot = q_f[77]-q_r[77];
    wdot[5] -= qdot;
    wdot[5] += qdot;
    wdot[9] += qdot;
    wdot[10] -= qdot;

    qdot = q_f[78]-q_r[78];
    wdot[9] += qdot;
    wdot[10] -= qdot;
    wdot[18] -= qdot;
    wdot[18] += qdot;

    qdot = q_f[79]-q_r[79];
    wdot[9] += qdot;
    wdot[10] -= qdot;
    wdot[19] -= qdot;
    wdot[19] += qdot;

    qdot = q_f[80]-q_r[80];
    wdot[10] -= qdot;
    wdot[14] += qdot;
    wdot[18] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[81]-q_r[81];
    wdot[0] -= qdot;
    wdot[4] += qdot;
    wdot[13] += qdot;
    wdot[14] -= qdot;

    qdot = q_f[82]-q_r[82];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[13] += qdot;
    wdot[14] -= qdot;

    qdot = q_f[83]-q_r[83];
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[13] += qdot;
    wdot[14] -= qdot;

    qdot = q_f[84]-q_r[84];
    wdot[3] += qdot;
    wdot[7] -= qdot;
    wdot[13] += qdot;
    wdot[14] -= qdot;

    qdot = q_f[85]-q_r[85];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[13] += qdot;
    wdot[14] -= qdot;

    qdot = q_f[86]-q_r[86];
    wdot[0] += qdot;
    wdot[8] -= qdot;
    wdot[14] -= qdot;
    wdot[27] += qdot;

    qdot = q_f[87]-q_r[87];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[11] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[88]-q_r[88];
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[9] += qdot;
    wdot[11] -= qdot;

    qdot = q_f[89]-q_r[89];
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[10] += qdot;
    wdot[11] -= qdot;

    qdot = q_f[90]-q_r[90];
    wdot[1] += qdot;
    wdot[7] -= qdot;
    wdot[11] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[91]-q_r[91];
    wdot[2] += qdot;
    wdot[7] -= qdot;
    wdot[11] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[92]-q_r[92];
    wdot[3] -= qdot;
    wdot[7] += qdot;
    wdot[11] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[93]-q_r[93];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[11] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[94]-q_r[94];
    wdot[0] += qdot;
    wdot[8] -= qdot;
    wdot[11] -= qdot;
    wdot[22] += qdot;

    qdot = q_f[95]-q_r[95];
    wdot[11] -= qdot;
    wdot[12] += qdot;
    wdot[13] -= qdot;
    wdot[18] += qdot;

    qdot = q_f[96]-q_r[96];
    wdot[11] -= qdot;
    wdot[12] += qdot;
    wdot[13] += qdot;
    wdot[14] -= qdot;

    qdot = q_f[97]-q_r[97];
    wdot[0] += qdot;
    wdot[9] -= qdot;
    wdot[11] -= qdot;
    wdot[23] += qdot;

    qdot = q_f[98]-q_r[98];
    wdot[0] += qdot;
    wdot[11] -= 2 * qdot;
    wdot[24] += qdot;

    qdot = q_f[99]-q_r[99];
    wdot[11] -= qdot;
    wdot[18] += qdot;
    wdot[23] += qdot;
    wdot[26] -= qdot;

    qdot = q_f[100]-q_r[100];
    wdot[0] -= qdot;
    wdot[4] += qdot;
    wdot[14] += qdot;
    wdot[15] -= qdot;

    qdot = q_f[101]-q_r[101];
    wdot[0] -= qdot;
    wdot[2] += qdot;
    wdot[11] += qdot;
    wdot[15] -= qdot;

    qdot = q_f[102]-q_r[102];
    wdot[0] -= qdot;
    wdot[5] += qdot;
    wdot[10] += qdot;
    wdot[15] -= qdot;

    qdot = q_f[103]-q_r[103];
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[14] += qdot;
    wdot[15] -= qdot;

    qdot = q_f[104]-q_r[104];
    wdot[3] += qdot;
    wdot[7] -= qdot;
    wdot[14] += qdot;
    wdot[15] -= qdot;

    qdot = q_f[105]-q_r[105];
    wdot[0] -= qdot;
    wdot[4] += qdot;
    wdot[14] += qdot;
    wdot[16] -= qdot;

    qdot = q_f[106]-q_r[106];
    wdot[0] -= qdot;
    wdot[2] += qdot;
    wdot[11] += qdot;
    wdot[16] -= qdot;

    qdot = q_f[107]-q_r[107];
    wdot[0] -= qdot;
    wdot[5] += qdot;
    wdot[10] += qdot;
    wdot[16] -= qdot;

    qdot = q_f[108]-q_r[108];
    wdot[3] += qdot;
    wdot[7] -= qdot;
    wdot[14] += qdot;
    wdot[16] -= qdot;

    qdot = q_f[109]-q_r[109];
    wdot[0] -= qdot;
    wdot[4] += qdot;
    wdot[11] += qdot;
    wdot[12] -= qdot;

    qdot = q_f[110]-q_r[110];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[11] += qdot;
    wdot[12] -= qdot;

    qdot = q_f[111]-q_r[111];
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[11] += qdot;
    wdot[12] -= qdot;

    qdot = q_f[112]-q_r[112];
    wdot[0] += qdot;
    wdot[8] -= qdot;
    wdot[12] -= qdot;
    wdot[23] += qdot;

    qdot = q_f[113]-q_r[113];
    wdot[9] -= qdot;
    wdot[11] += 2 * qdot;
    wdot[12] -= qdot;

    qdot = q_f[114]-q_r[114];
    wdot[10] -= qdot;
    wdot[11] += 2 * qdot;
    wdot[12] -= qdot;

    qdot = q_f[115]-q_r[115];
    wdot[0] -= qdot;
    wdot[4] += qdot;
    wdot[16] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[116]-q_r[116];
    wdot[0] -= qdot;
    wdot[4] += qdot;
    wdot[15] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[117]-q_r[117];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[16] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[118]-q_r[118];
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[16] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[119]-q_r[119];
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[15] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[120]-q_r[120];
    wdot[1] -= qdot;
    wdot[8] += qdot;
    wdot[18] += qdot;
    wdot[20] -= qdot;

    qdot = q_f[121]-q_r[121];
    wdot[0] += qdot;
    wdot[2] -= qdot;
    wdot[20] -= qdot;
    wdot[26] += qdot;

    qdot = q_f[122]-q_r[122];
    wdot[7] -= qdot;
    wdot[13] += qdot;
    wdot[18] += qdot;
    wdot[20] -= qdot;

    qdot = q_f[123]-q_r[123];
    wdot[0] += qdot;
    wdot[4] -= qdot;
    wdot[20] -= qdot;
    wdot[21] += qdot;

    qdot = q_f[124]-q_r[124];
    wdot[0] -= qdot;
    wdot[10] += qdot;
    wdot[18] += qdot;
    wdot[26] -= qdot;

    qdot = q_f[125]-q_r[125];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[18] += 2 * qdot;
    wdot[26] -= qdot;

    qdot = q_f[126]-q_r[126];
    wdot[2] += qdot;
    wdot[7] -= qdot;
    wdot[18] += 2 * qdot;
    wdot[26] -= qdot;

    qdot = q_f[127]-q_r[127];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[20] += qdot;
    wdot[21] -= qdot;

    qdot = q_f[128]-q_r[128];
    wdot[1] -= qdot;
    wdot[9] += qdot;
    wdot[18] += qdot;
    wdot[21] -= qdot;

    qdot = q_f[129]-q_r[129];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[21] -= qdot;
    wdot[26] += qdot;

    qdot = q_f[130]-q_r[130];
    wdot[0] += qdot;
    wdot[2] -= qdot;
    wdot[21] -= qdot;
    wdot[27] += qdot;

    qdot = q_f[131]-q_r[131];
    wdot[0] += qdot;
    wdot[2] -= qdot;
    wdot[21] -= qdot;
    wdot[27] += qdot;

    qdot = q_f[132]-q_r[132];
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[20] += qdot;
    wdot[21] -= qdot;

    qdot = q_f[133]-q_r[133];
    wdot[13] -= qdot;
    wdot[18] += qdot;
    wdot[21] -= qdot;
    wdot[22] += qdot;

    qdot = q_f[134]-q_r[134];
    wdot[0] += qdot;
    wdot[9] -= qdot;
    wdot[21] -= qdot;
    wdot[31] += qdot;

    qdot = q_f[135]-q_r[135];
    wdot[0] += qdot;
    wdot[10] -= qdot;
    wdot[21] -= qdot;
    wdot[31] += qdot;

    qdot = q_f[136]-q_r[136];
    wdot[0] += qdot;
    wdot[20] -= qdot;
    wdot[21] -= qdot;
    wdot[40] += qdot;

    qdot = q_f[137]-q_r[137];
    wdot[0] += qdot;
    wdot[11] -= qdot;
    wdot[21] -= qdot;
    wdot[32] += qdot;

    qdot = q_f[138]-q_r[138];
    wdot[0] += qdot;
    wdot[11] -= qdot;
    wdot[21] -= qdot;
    wdot[33] += qdot;

    qdot = q_f[139]-q_r[139];
    wdot[0] -= qdot;
    wdot[4] += qdot;
    wdot[26] += qdot;
    wdot[27] -= qdot;

    qdot = q_f[140]-q_r[140];
    wdot[0] -= qdot;
    wdot[11] += qdot;
    wdot[18] += qdot;
    wdot[27] -= qdot;

    qdot = q_f[141]-q_r[141];
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[26] += qdot;
    wdot[27] -= qdot;

    qdot = q_f[142]-q_r[142];
    wdot[0] -= qdot;
    wdot[4] += qdot;
    wdot[21] += qdot;
    wdot[22] -= qdot;

    qdot = q_f[143]-q_r[143];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[22] -= qdot;
    wdot[27] += qdot;

    qdot = q_f[144]-q_r[144];
    wdot[1] -= qdot;
    wdot[11] += qdot;
    wdot[18] += qdot;
    wdot[22] -= qdot;

    qdot = q_f[145]-q_r[145];
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[21] += qdot;
    wdot[22] -= qdot;

    qdot = q_f[146]-q_r[146];
    wdot[3] += qdot;
    wdot[7] -= qdot;
    wdot[21] += qdot;
    wdot[22] -= qdot;

    qdot = q_f[147]-q_r[147];
    wdot[1] += qdot;
    wdot[7] -= qdot;
    wdot[22] -= qdot;
    wdot[29] += qdot;

    qdot = q_f[148]-q_r[148];
    wdot[7] -= qdot;
    wdot[13] += qdot;
    wdot[14] += qdot;
    wdot[22] -= qdot;

    qdot = q_f[149]-q_r[149];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[22] -= qdot;
    wdot[29] += qdot;

    qdot = q_f[150]-q_r[150];
    wdot[13] -= qdot;
    wdot[18] += qdot;
    wdot[22] -= qdot;
    wdot[23] += qdot;

    qdot = q_f[151]-q_r[151];
    wdot[13] -= qdot;
    wdot[22] -= qdot;
    wdot[39] += qdot;

    qdot = q_f[152]-q_r[152];
    wdot[0] += qdot;
    wdot[11] -= qdot;
    wdot[22] -= qdot;
    wdot[34] += qdot;

    qdot = q_f[153]-q_r[153];
    wdot[11] += qdot;
    wdot[18] += qdot;
    wdot[29] -= qdot;

    qdot = q_f[154]-q_r[154];
    wdot[0] -= qdot;
    wdot[0] += qdot;
    wdot[28] += qdot;
    wdot[29] -= qdot;

    qdot = q_f[155]-q_r[155];
    wdot[0] -= qdot;
    wdot[11] += qdot;
    wdot[13] += qdot;
    wdot[29] -= qdot;

    qdot = q_f[156]-q_r[156];
    wdot[0] -= qdot;
    wdot[4] += qdot;
    wdot[27] += qdot;
    wdot[29] -= qdot;

    qdot = q_f[157]-q_r[157];
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[27] += qdot;
    wdot[29] -= qdot;

    qdot = q_f[158]-q_r[158];
    wdot[3] += qdot;
    wdot[7] -= qdot;
    wdot[27] += qdot;
    wdot[29] -= qdot;

    qdot = q_f[159]-q_r[159];
    wdot[2] += qdot;
    wdot[7] -= qdot;
    wdot[14] += qdot;
    wdot[18] += qdot;
    wdot[29] -= qdot;

    qdot = q_f[160]-q_r[160];
    wdot[0] -= qdot;
    wdot[11] += qdot;
    wdot[13] += qdot;
    wdot[28] -= qdot;

    qdot = q_f[161]-q_r[161];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[11] += qdot;
    wdot[19] += qdot;
    wdot[28] -= qdot;

    qdot = q_f[162]-q_r[162];
    wdot[0] -= qdot;
    wdot[4] += qdot;
    wdot[28] += qdot;
    wdot[30] -= qdot;

    qdot = q_f[163]-q_r[163];
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[28] += qdot;
    wdot[30] -= qdot;

    qdot = q_f[164]-q_r[164];
    wdot[11] -= qdot;
    wdot[12] += qdot;
    wdot[28] += qdot;
    wdot[30] -= qdot;

    qdot = q_f[165]-q_r[165];
    wdot[3] += qdot;
    wdot[7] -= qdot;
    wdot[28] += qdot;
    wdot[30] -= qdot;

    qdot = q_f[166]-q_r[166];
    wdot[0] -= qdot;
    wdot[4] += qdot;
    wdot[22] += qdot;
    wdot[23] -= qdot;

    qdot = q_f[167]-q_r[167];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[22] += qdot;
    wdot[23] -= qdot;

    qdot = q_f[168]-q_r[168];
    wdot[1] -= qdot;
    wdot[11] += qdot;
    wdot[13] += qdot;
    wdot[23] -= qdot;

    qdot = q_f[169]-q_r[169];
    wdot[1] -= qdot;
    wdot[9] += qdot;
    wdot[14] += qdot;
    wdot[23] -= qdot;

    qdot = q_f[170]-q_r[170];
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[22] += qdot;
    wdot[23] -= qdot;

    qdot = q_f[171]-q_r[171];
    wdot[13] -= qdot;
    wdot[18] += qdot;
    wdot[23] -= qdot;
    wdot[24] += qdot;

    qdot = q_f[172]-q_r[172];
    wdot[0] += qdot;
    wdot[8] -= qdot;
    wdot[23] -= qdot;
    wdot[33] += qdot;

    qdot = q_f[173]-q_r[173];
    wdot[0] += qdot;
    wdot[8] -= qdot;
    wdot[23] -= qdot;
    wdot[32] += qdot;

    qdot = q_f[174]-q_r[174];
    wdot[0] += qdot;
    wdot[9] -= qdot;
    wdot[23] -= qdot;
    wdot[34] += qdot;

    qdot = q_f[175]-q_r[175];
    wdot[0] += qdot;
    wdot[10] -= qdot;
    wdot[23] -= qdot;
    wdot[34] += qdot;

    qdot = q_f[176]-q_r[176];
    wdot[11] -= qdot;
    wdot[12] += qdot;
    wdot[22] += qdot;
    wdot[23] -= qdot;

    qdot = q_f[177]-q_r[177];
    wdot[11] -= qdot;
    wdot[23] -= qdot;
    wdot[37] += qdot;

    qdot = q_f[178]-q_r[178];
    wdot[22] -= qdot;
    wdot[23] -= qdot;
    wdot[48] += qdot;

    qdot = q_f[179]-q_r[179];
    wdot[0] -= qdot;
    wdot[4] += qdot;
    wdot[23] += qdot;
    wdot[24] -= qdot;

    qdot = q_f[180]-q_r[180];
    wdot[1] -= qdot;
    wdot[11] += qdot;
    wdot[14] += qdot;
    wdot[24] -= qdot;

    qdot = q_f[181]-q_r[181];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[24] -= qdot;
    wdot[30] += qdot;

    qdot = q_f[182]-q_r[182];
    wdot[3] += qdot;
    wdot[7] -= qdot;
    wdot[23] += qdot;
    wdot[24] -= qdot;

    qdot = q_f[183]-q_r[183];
    wdot[3] -= qdot;
    wdot[7] += qdot;
    wdot[24] -= qdot;
    wdot[25] += qdot;

    qdot = q_f[184]-q_r[184];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[23] += qdot;
    wdot[24] -= qdot;

    qdot = q_f[185]-q_r[185];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[11] += qdot;
    wdot[14] += qdot;
    wdot[24] -= qdot;

    qdot = q_f[186]-q_r[186];
    wdot[11] += qdot;
    wdot[22] -= qdot;
    wdot[24] -= qdot;
    wdot[34] += qdot;

    qdot = q_f[187]-q_r[187];
    wdot[0] -= qdot;
    wdot[4] += qdot;
    wdot[24] += qdot;
    wdot[25] -= qdot;

    qdot = q_f[188]-q_r[188];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[24] += qdot;
    wdot[25] -= qdot;

    qdot = q_f[189]-q_r[189];
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[24] += qdot;
    wdot[25] -= qdot;

    qdot = q_f[190]-q_r[190];
    wdot[10] -= qdot;
    wdot[11] += qdot;
    wdot[24] += qdot;
    wdot[25] -= qdot;

    qdot = q_f[191]-q_r[191];
    wdot[11] -= qdot;
    wdot[12] += qdot;
    wdot[24] += qdot;
    wdot[25] -= qdot;

    qdot = q_f[192]-q_r[192];
    wdot[0] -= qdot;
    wdot[31] -= qdot;
    wdot[32] += qdot;

    qdot = q_f[193]-q_r[193];
    wdot[1] -= qdot;
    wdot[14] += qdot;
    wdot[20] += qdot;
    wdot[31] -= qdot;

    qdot = q_f[194]-q_r[194];
    wdot[7] -= qdot;
    wdot[13] += qdot;
    wdot[27] += qdot;
    wdot[31] -= qdot;

    qdot = q_f[195]-q_r[195];
    wdot[3] -= qdot;
    wdot[7] += qdot;
    wdot[31] -= qdot;
    wdot[32] += qdot;

    qdot = q_f[196]-q_r[196];
    wdot[0] -= qdot;
    wdot[33] -= qdot;
    wdot[35] += qdot;

    qdot = q_f[197]-q_r[197];
    wdot[0] -= qdot;
    wdot[33] -= qdot;
    wdot[34] += qdot;

    qdot = q_f[198]-q_r[198];
    wdot[1] -= qdot;
    wdot[18] += qdot;
    wdot[23] += qdot;
    wdot[33] -= qdot;

    qdot = q_f[199]-q_r[199];
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[31] += qdot;
    wdot[33] -= qdot;

    qdot = q_f[200]-q_r[200];
    wdot[32] -= qdot;
    wdot[33] += qdot;

    qdot = q_f[201]-q_r[201];
    wdot[0] -= qdot;
    wdot[0] += qdot;
    wdot[32] -= qdot;
    wdot[33] += qdot;

    qdot = q_f[202]-q_r[202];
    wdot[0] -= qdot;
    wdot[32] -= qdot;
    wdot[35] += qdot;

    qdot = q_f[203]-q_r[203];
    wdot[1] -= qdot;
    wdot[18] += qdot;
    wdot[23] += qdot;
    wdot[32] -= qdot;

    qdot = q_f[204]-q_r[204];
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[31] += qdot;
    wdot[32] -= qdot;

    qdot = q_f[205]-q_r[205];
    wdot[0] -= qdot;
    wdot[4] += qdot;
    wdot[33] += qdot;
    wdot[34] -= qdot;

    qdot = q_f[206]-q_r[206];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[34] -= qdot;
    wdot[39] += qdot;

    qdot = q_f[207]-q_r[207];
    wdot[0] += 2 * qdot;
    wdot[2] -= qdot;
    wdot[34] -= qdot;
    wdot[39] += qdot;

    qdot = q_f[208]-q_r[208];
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[33] += qdot;
    wdot[34] -= qdot;

    qdot = q_f[209]-q_r[209];
    wdot[3] -= qdot;
    wdot[7] += qdot;
    wdot[34] -= qdot;
    wdot[36] += qdot;

    qdot = q_f[210]-q_r[210];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[14] += qdot;
    wdot[22] += qdot;
    wdot[34] -= qdot;

    qdot = q_f[211]-q_r[211];
    wdot[13] -= qdot;
    wdot[18] += qdot;
    wdot[34] -= qdot;
    wdot[36] += qdot;

    qdot = q_f[212]-q_r[212];
    wdot[11] -= qdot;
    wdot[12] += qdot;
    wdot[33] += qdot;
    wdot[34] -= qdot;

    qdot = q_f[213]-q_r[213];
    wdot[7] -= qdot;
    wdot[14] += qdot;
    wdot[28] += qdot;
    wdot[35] -= qdot;

    qdot = q_f[214]-q_r[214];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[11] += qdot;
    wdot[27] += qdot;
    wdot[35] -= qdot;

    qdot = q_f[215]-q_r[215];
    wdot[0] -= qdot;
    wdot[11] += qdot;
    wdot[23] += qdot;
    wdot[36] -= qdot;

    qdot = q_f[216]-q_r[216];
    wdot[0] -= qdot;
    wdot[4] += qdot;
    wdot[34] += qdot;
    wdot[36] -= qdot;

    qdot = q_f[217]-q_r[217];
    wdot[0] -= qdot;
    wdot[4] += qdot;
    wdot[35] += qdot;
    wdot[36] -= qdot;

    qdot = q_f[218]-q_r[218];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[11] += qdot;
    wdot[27] += qdot;
    wdot[36] -= qdot;

    qdot = q_f[219]-q_r[219];
    wdot[0] += 2 * qdot;
    wdot[1] -= qdot;
    wdot[36] -= qdot;
    wdot[39] += qdot;

    qdot = q_f[220]-q_r[220];
    wdot[1] -= qdot;
    wdot[13] += qdot;
    wdot[24] += qdot;
    wdot[36] -= qdot;

    qdot = q_f[221]-q_r[221];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[34] += qdot;
    wdot[36] -= qdot;

    qdot = q_f[222]-q_r[222];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[35] += qdot;
    wdot[36] -= qdot;

    qdot = q_f[223]-q_r[223];
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[34] += qdot;
    wdot[36] -= qdot;

    qdot = q_f[224]-q_r[224];
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[35] += qdot;
    wdot[36] -= qdot;

    qdot = q_f[225]-q_r[225];
    wdot[11] -= qdot;
    wdot[12] += qdot;
    wdot[34] += qdot;
    wdot[36] -= qdot;

    qdot = q_f[226]-q_r[226];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[18] += qdot;
    wdot[22] += qdot;
    wdot[39] -= qdot;

    qdot = q_f[227]-q_r[227];
    wdot[1] -= qdot;
    wdot[14] += qdot;
    wdot[27] += qdot;
    wdot[39] -= qdot;

    qdot = q_f[228]-q_r[228];
    wdot[0] -= qdot;
    wdot[11] += qdot;
    wdot[24] += qdot;
    wdot[38] -= qdot;

    qdot = q_f[229]-q_r[229];
    wdot[1] -= qdot;
    wdot[11] += qdot;
    wdot[30] += qdot;
    wdot[38] -= qdot;

    qdot = q_f[230]-q_r[230];
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[36] += qdot;
    wdot[38] -= qdot;

    qdot = q_f[231]-q_r[231];
    wdot[3] += qdot;
    wdot[7] -= qdot;
    wdot[36] += qdot;
    wdot[38] -= qdot;

    qdot = q_f[232]-q_r[232];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[11] += qdot;
    wdot[30] += qdot;
    wdot[38] -= qdot;

    qdot = q_f[233]-q_r[233];
    wdot[11] -= qdot;
    wdot[12] += qdot;
    wdot[36] += qdot;
    wdot[38] -= qdot;

    qdot = q_f[234]-q_r[234];
    wdot[0] -= qdot;
    wdot[11] += qdot;
    wdot[24] += qdot;
    wdot[37] -= qdot;

    qdot = q_f[235]-q_r[235];
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[36] += qdot;
    wdot[37] -= qdot;

    qdot = q_f[236]-q_r[236];
    wdot[3] += qdot;
    wdot[7] -= qdot;
    wdot[36] += qdot;
    wdot[37] -= qdot;

    qdot = q_f[237]-q_r[237];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[14] += qdot;
    wdot[24] += qdot;
    wdot[37] -= qdot;

    qdot = q_f[238]-q_r[238];
    wdot[11] -= qdot;
    wdot[12] += qdot;
    wdot[36] += qdot;
    wdot[37] -= qdot;

    qdot = q_f[239]-q_r[239];
    wdot[0] -= qdot;
    wdot[40] -= qdot;
    wdot[41] += qdot;

    qdot = q_f[240]-q_r[240];
    wdot[0] -= qdot;
    wdot[4] += qdot;
    wdot[40] += qdot;
    wdot[41] -= qdot;

    qdot = q_f[241]-q_r[241];
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[41] += qdot;
    wdot[42] -= qdot;

    qdot = q_f[242]-q_r[242];
    wdot[43] += qdot;
    wdot[44] -= qdot;

    qdot = q_f[243]-q_r[243];
    wdot[0] -= qdot;
    wdot[22] += qdot;
    wdot[23] += qdot;
    wdot[45] -= qdot;

    qdot = q_f[244]-q_r[244];
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[43] += qdot;
    wdot[45] -= qdot;

    qdot = q_f[245]-q_r[245];
    wdot[0] += qdot;
    wdot[43] += qdot;
    wdot[46] -= qdot;

    qdot = q_f[246]-q_r[246];
    wdot[0] += qdot;
    wdot[44] += qdot;
    wdot[47] -= qdot;

    qdot = q_f[247]-q_r[247];
    wdot[0] += qdot;
    wdot[45] += qdot;
    wdot[48] -= qdot;

    qdot = q_f[248]-q_r[248];
    wdot[3] += qdot;
    wdot[7] -= qdot;
    wdot[45] += qdot;
    wdot[48] -= qdot;

    qdot = q_f[249]-q_r[249];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[14] += qdot;
    wdot[34] += qdot;
    wdot[48] -= qdot;

    qdot = q_f[250]-q_r[250];
    wdot[0] -= qdot;
    wdot[23] += qdot;
    wdot[24] += qdot;
    wdot[49] -= qdot;

    qdot = q_f[251]-q_r[251];
    wdot[0] -= qdot;
    wdot[11] += qdot;
    wdot[36] += qdot;
    wdot[49] -= qdot;

    qdot = q_f[252]-q_r[252];
    wdot[0] -= qdot;
    wdot[4] += qdot;
    wdot[48] += qdot;
    wdot[49] -= qdot;

    qdot = q_f[253]-q_r[253];
    wdot[1] -= qdot;
    wdot[13] += qdot;
    wdot[37] += qdot;
    wdot[49] -= qdot;

    qdot = q_f[254]-q_r[254];
    wdot[23] -= qdot;
    wdot[24] -= qdot;
    wdot[50] += qdot;

    qdot = q_f[255]-q_r[255];
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[49] += qdot;
    wdot[50] -= qdot;

    qdot = q_f[256]-q_r[256];
    wdot[3] += qdot;
    wdot[7] -= qdot;
    wdot[49] += qdot;
    wdot[50] -= qdot;

    qdot = q_f[257]-q_r[257];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[14] += qdot;
    wdot[37] += qdot;
    wdot[50] -= qdot;

    qdot = q_f[258]-q_r[258];
    wdot[11] -= qdot;
    wdot[12] += qdot;
    wdot[49] += qdot;
    wdot[50] -= qdot;

    qdot = q_f[259]-q_r[259];
    wdot[23] += 3 * qdot;
    wdot[37] += 2 * qdot;
    wdot[51] -= qdot;

    qdot = q_f[260]-q_r[260];
    wdot[23] += 2 * qdot;
    wdot[50] += 2 * qdot;
    wdot[51] -= qdot;

    qdot = q_f[261]-q_r[261];
    wdot[0] -= qdot;
    wdot[4] += qdot;
    wdot[23] += 4 * qdot;
    wdot[50] += qdot;
    wdot[51] -= qdot;

    qdot = q_f[262]-q_r[262];
    wdot[0] -= qdot;
    wdot[4] += qdot;
    wdot[23] += 2 * qdot;
    wdot[49] += qdot;
    wdot[50] += qdot;
    wdot[51] -= qdot;

    qdot = q_f[263]-q_r[263];
    wdot[0] -= qdot;
    wdot[4] += qdot;
    wdot[36] += qdot;
    wdot[37] += qdot;
    wdot[51] -= qdot;
    wdot[52] += qdot;

    qdot = q_f[264]-q_r[264];
    wdot[0] -= qdot;
    wdot[4] += qdot;
    wdot[23] += 2 * qdot;
    wdot[37] += qdot;
    wdot[51] -= qdot;
    wdot[54] += qdot;

    qdot = q_f[265]-q_r[265];
    wdot[0] -= qdot;
    wdot[4] += qdot;
    wdot[23] += qdot;
    wdot[50] += qdot;
    wdot[51] -= qdot;
    wdot[52] += qdot;

    qdot = q_f[266]-q_r[266];
    wdot[11] -= qdot;
    wdot[12] += qdot;
    wdot[23] += 4 * qdot;
    wdot[50] += qdot;
    wdot[51] -= qdot;

    qdot = q_f[267]-q_r[267];
    wdot[11] -= qdot;
    wdot[12] += qdot;
    wdot[23] += 2 * qdot;
    wdot[49] += qdot;
    wdot[50] += qdot;
    wdot[51] -= qdot;

    qdot = q_f[268]-q_r[268];
    wdot[11] -= qdot;
    wdot[12] += qdot;
    wdot[36] += qdot;
    wdot[37] += qdot;
    wdot[51] -= qdot;
    wdot[52] += qdot;

    qdot = q_f[269]-q_r[269];
    wdot[11] -= qdot;
    wdot[12] += qdot;
    wdot[23] += 2 * qdot;
    wdot[37] += qdot;
    wdot[51] -= qdot;
    wdot[54] += qdot;

    qdot = q_f[270]-q_r[270];
    wdot[11] -= qdot;
    wdot[12] += qdot;
    wdot[23] += qdot;
    wdot[50] += qdot;
    wdot[51] -= qdot;
    wdot[52] += qdot;

    qdot = q_f[271]-q_r[271];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[23] += 4 * qdot;
    wdot[50] += qdot;
    wdot[51] -= qdot;

    qdot = q_f[272]-q_r[272];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[23] += 2 * qdot;
    wdot[49] += qdot;
    wdot[50] += qdot;
    wdot[51] -= qdot;

    qdot = q_f[273]-q_r[273];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[36] += qdot;
    wdot[37] += qdot;
    wdot[51] -= qdot;
    wdot[52] += qdot;

    qdot = q_f[274]-q_r[274];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[23] += 2 * qdot;
    wdot[37] += qdot;
    wdot[51] -= qdot;
    wdot[54] += qdot;

    qdot = q_f[275]-q_r[275];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[23] += qdot;
    wdot[50] += qdot;
    wdot[51] -= qdot;
    wdot[52] += qdot;

    qdot = q_f[276]-q_r[276];
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[23] += 4 * qdot;
    wdot[50] += qdot;
    wdot[51] -= qdot;

    qdot = q_f[277]-q_r[277];
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[23] += 2 * qdot;
    wdot[49] += qdot;
    wdot[50] += qdot;
    wdot[51] -= qdot;

    qdot = q_f[278]-q_r[278];
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[36] += qdot;
    wdot[37] += qdot;
    wdot[51] -= qdot;
    wdot[52] += qdot;

    qdot = q_f[279]-q_r[279];
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[23] += 2 * qdot;
    wdot[37] += qdot;
    wdot[51] -= qdot;
    wdot[54] += qdot;

    qdot = q_f[280]-q_r[280];
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[23] += qdot;
    wdot[50] += qdot;
    wdot[51] -= qdot;
    wdot[52] += qdot;

    qdot = q_f[281]-q_r[281];
    wdot[0] -= qdot;
    wdot[23] += qdot;
    wdot[50] += qdot;
    wdot[52] -= qdot;

    qdot = q_f[282]-q_r[282];
    wdot[0] -= qdot;
    wdot[4] += qdot;
    wdot[52] -= qdot;
    wdot[53] += qdot;

    qdot = q_f[283]-q_r[283];
    wdot[0] -= qdot;
    wdot[23] += qdot;
    wdot[37] += qdot;
    wdot[54] -= qdot;

    qdot = q_f[284]-q_r[284];
    wdot[0] -= qdot;
    wdot[4] += qdot;
    wdot[23] += qdot;
    wdot[34] += qdot;
    wdot[54] -= qdot;

    qdot = q_f[285]-q_r[285];
    wdot[0] -= qdot;
    wdot[11] += qdot;
    wdot[23] += qdot;
    wdot[34] += qdot;
    wdot[53] -= qdot;

    qdot = q_f[286]-q_r[286];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[14] += qdot;
    wdot[23] += qdot;
    wdot[34] += qdot;
    wdot[53] -= qdot;

    qdot = q_f[287]-q_r[287];
    wdot[1] -= qdot;
    wdot[13] += qdot;
    wdot[23] += qdot;
    wdot[37] += qdot;
    wdot[52] -= qdot;

    qdot = q_f[288]-q_r[288];
    wdot[1] -= qdot;
    wdot[13] += qdot;
    wdot[50] += qdot;
    wdot[54] -= qdot;

    return;
}

__device__ void comp_k_f_d(double * tc, double invT, double * k_f, double * Corr, double * sc)
{
    double fwd_A[289], fwd_beta[289], fwd_Ea[289];
    double low_A[289], low_beta[289], low_Ea[289];
    double rev_A[289], rev_beta[289], rev_Ea[289];
    double troe_a[289],troe_Ts[289], troe_Tss[289], troe_Tsss[289];
    double sri_a[289], sri_b[289], sri_c[289], sri_d[289], sri_e[289];
    double activation_units[289], prefactor_units[289], phase_units[289];
    int is_PD[289], troe_len[289], sri_len[289], nTB[289], *TBid[289];
    double *TB[289];

    // (0):  H + O2 (+M) <=> HO2 (+M)
    fwd_A[0]     = 5116000000000;
    fwd_beta[0]  = 0.44;
    fwd_Ea[0]    = 0;
    low_A[0]     = 6.328e+19;
    low_beta[0]  = -1.3999999999999999;
    low_Ea[0]    = 0;
    troe_a[0]    = 0.5;
    troe_Tsss[0] = 1.0000000000000001e-30;
    troe_Ts[0]   = 1e+30;
    troe_len[0]  = 3;
    prefactor_units[0]  = 1.0000000000000002e-06;
    activation_units[0] = 0.50321666580471969;
    phase_units[0]      = 1e-12;
    is_PD[0] = 1;
    nTB[0] = 4;
    TB[0] = (double *) malloc(4 * sizeof(double));
    TBid[0] = (int *) malloc(4 * sizeof(int));
    TBid[0][0] = 5; TB[0][0] = 11.890000000000001; // H2O
    TBid[0][1] = 7; TB[0][1] = 0.84999999999999998; // O2
    TBid[0][2] = 18; TB[0][2] = 1.0900000000000001; // CO
    TBid[0][3] = 19; TB[0][3] = 2.1800000000000002; // CO2

    // (1):  2 OH (+M) <=> H2O2 (+M)
    fwd_A[1]     = 111000000000000;
    fwd_beta[1]  = -0.37;
    fwd_Ea[1]    = 0;
    low_A[1]     = 2.01e+17;
    low_beta[1]  = -0.58399999999999996;
    low_Ea[1]    = -2293;
    troe_a[1]    = 0.73460000000000003;
    troe_Tsss[1] = 94;
    troe_Ts[1]   = 1756;
    troe_Tss[1]  = 5182;
    troe_len[1]  = 4;
    prefactor_units[1]  = 1.0000000000000002e-06;
    activation_units[1] = 0.50321666580471969;
    phase_units[1]      = 1e-12;
    is_PD[1] = 1;
    nTB[1] = 4;
    TB[1] = (double *) malloc(4 * sizeof(double));
    TBid[1] = (int *) malloc(4 * sizeof(int));
    TBid[1][0] = 4; TB[1][0] = 2; // H2
    TBid[1][1] = 5; TB[1][1] = 6; // H2O
    TBid[1][2] = 18; TB[1][2] = 1.75; // CO
    TBid[1][3] = 19; TB[1][3] = 3.6000000000000001; // CO2

    // (2):  CH2 + CO (+M) <=> CH2CO (+M)
    fwd_A[2]     = 810000000000;
    fwd_beta[2]  = 0.5;
    fwd_Ea[2]    = 4510;
    low_A[2]     = 2.69e+33;
    low_beta[2]  = -5.1100000000000003;
    low_Ea[2]    = 7095;
    troe_a[2]    = 0.5907;
    troe_Tsss[2] = 275;
    troe_Ts[2]   = 1226;
    troe_Tss[2]  = 5185;
    troe_len[2]  = 4;
    prefactor_units[2]  = 1.0000000000000002e-06;
    activation_units[2] = 0.50321666580471969;
    phase_units[2]      = 1e-12;
    is_PD[2] = 1;
    nTB[2] = 6;
    TB[2] = (double *) malloc(6 * sizeof(double));
    TBid[2] = (int *) malloc(6 * sizeof(int));
    TBid[2][0] = 4; TB[2][0] = 2; // H2
    TBid[2][1] = 5; TB[2][1] = 6; // H2O
    TBid[2][2] = 12; TB[2][2] = 2; // CH4
    TBid[2][3] = 18; TB[2][3] = 1.5; // CO
    TBid[2][4] = 19; TB[2][4] = 2; // CO2
    TBid[2][5] = 25; TB[2][5] = 3; // C2H6

    // (3):  CH2* + H2O (+M) <=> CH3OH (+M)
    fwd_A[3]     = 20000000000000;
    fwd_beta[3]  = 0;
    fwd_Ea[3]    = 0;
    low_A[3]     = 2.7e+38;
    low_beta[3]  = -6.2999999999999998;
    low_Ea[3]    = 3100;
    troe_a[3]    = 0.1507;
    troe_Tsss[3] = 134;
    troe_Ts[3]   = 2383;
    troe_Tss[3]  = 7265;
    troe_len[3]  = 4;
    prefactor_units[3]  = 1.0000000000000002e-06;
    activation_units[3] = 0.50321666580471969;
    phase_units[3]      = 1e-12;
    is_PD[3] = 1;
    nTB[3] = 6;
    TB[3] = (double *) malloc(6 * sizeof(double));
    TBid[3] = (int *) malloc(6 * sizeof(int));
    TBid[3][0] = 4; TB[3][0] = 2; // H2
    TBid[3][1] = 5; TB[3][1] = 6; // H2O
    TBid[3][2] = 12; TB[3][2] = 2; // CH4
    TBid[3][3] = 18; TB[3][3] = 1.5; // CO
    TBid[3][4] = 19; TB[3][4] = 2; // CO2
    TBid[3][5] = 25; TB[3][5] = 3; // C2H6

    // (4):  CH2O + H (+M) <=> CH2OH (+M)
    fwd_A[4]     = 540000000000;
    fwd_beta[4]  = 0.45400000000000001;
    fwd_Ea[4]    = 3600;
    low_A[4]     = 1.27e+32;
    low_beta[4]  = -4.8200000000000003;
    low_Ea[4]    = 6530;
    troe_a[4]    = 0.71870000000000001;
    troe_Tsss[4] = 103;
    troe_Ts[4]   = 1291;
    troe_Tss[4]  = 4160;
    troe_len[4]  = 4;
    prefactor_units[4]  = 1.0000000000000002e-06;
    activation_units[4] = 0.50321666580471969;
    phase_units[4]      = 1e-12;
    is_PD[4] = 1;
    nTB[4] = 6;
    TB[4] = (double *) malloc(6 * sizeof(double));
    TBid[4] = (int *) malloc(6 * sizeof(int));
    TBid[4][0] = 4; TB[4][0] = 2; // H2
    TBid[4][1] = 5; TB[4][1] = 6; // H2O
    TBid[4][2] = 12; TB[4][2] = 2; // CH4
    TBid[4][3] = 18; TB[4][3] = 1.5; // CO
    TBid[4][4] = 19; TB[4][4] = 2; // CO2
    TBid[4][5] = 25; TB[4][5] = 3; // C2H6

    // (5):  CH2O + H (+M) <=> CH3O (+M)
    fwd_A[5]     = 540000000000;
    fwd_beta[5]  = 0.45400000000000001;
    fwd_Ea[5]    = 2600;
    low_A[5]     = 2.2e+30;
    low_beta[5]  = -4.7999999999999998;
    low_Ea[5]    = 5560;
    troe_a[5]    = 0.75800000000000001;
    troe_Tsss[5] = 94;
    troe_Ts[5]   = 1555;
    troe_Tss[5]  = 4200;
    troe_len[5]  = 4;
    prefactor_units[5]  = 1.0000000000000002e-06;
    activation_units[5] = 0.50321666580471969;
    phase_units[5]      = 1e-12;
    is_PD[5] = 1;
    nTB[5] = 6;
    TB[5] = (double *) malloc(6 * sizeof(double));
    TBid[5] = (int *) malloc(6 * sizeof(int));
    TBid[5][0] = 4; TB[5][0] = 2; // H2
    TBid[5][1] = 5; TB[5][1] = 6; // H2O
    TBid[5][2] = 12; TB[5][2] = 2; // CH4
    TBid[5][3] = 18; TB[5][3] = 1.5; // CO
    TBid[5][4] = 19; TB[5][4] = 2; // CO2
    TBid[5][5] = 25; TB[5][5] = 3; // C2H6

    // (6):  CH3 + H (+M) <=> CH4 (+M)
    fwd_A[6]     = 12700000000000000;
    fwd_beta[6]  = -0.63;
    fwd_Ea[6]    = 383;
    low_A[6]     = 2.4769999999999999e+33;
    low_beta[6]  = -4.7599999999999998;
    low_Ea[6]    = 2440;
    troe_a[6]    = 0.78300000000000003;
    troe_Tsss[6] = 74;
    troe_Ts[6]   = 2941;
    troe_Tss[6]  = 6964;
    troe_len[6]  = 4;
    prefactor_units[6]  = 1.0000000000000002e-06;
    activation_units[6] = 0.50321666580471969;
    phase_units[6]      = 1e-12;
    is_PD[6] = 1;
    nTB[6] = 6;
    TB[6] = (double *) malloc(6 * sizeof(double));
    TBid[6] = (int *) malloc(6 * sizeof(int));
    TBid[6][0] = 4; TB[6][0] = 2; // H2
    TBid[6][1] = 5; TB[6][1] = 6; // H2O
    TBid[6][2] = 12; TB[6][2] = 2; // CH4
    TBid[6][3] = 18; TB[6][3] = 1.5; // CO
    TBid[6][4] = 19; TB[6][4] = 2; // CO2
    TBid[6][5] = 25; TB[6][5] = 3; // C2H6

    // (7):  CH3 + OH (+M) <=> CH3OH (+M)
    fwd_A[7]     = 63000000000000;
    fwd_beta[7]  = 0;
    fwd_Ea[7]    = 0;
    low_A[7]     = 2.7e+38;
    low_beta[7]  = -6.2999999999999998;
    low_Ea[7]    = 3100;
    troe_a[7]    = 0.21049999999999999;
    troe_Tsss[7] = 83.5;
    troe_Ts[7]   = 5398;
    troe_Tss[7]  = 8370;
    troe_len[7]  = 4;
    prefactor_units[7]  = 1.0000000000000002e-06;
    activation_units[7] = 0.50321666580471969;
    phase_units[7]      = 1e-12;
    is_PD[7] = 1;
    nTB[7] = 6;
    TB[7] = (double *) malloc(6 * sizeof(double));
    TBid[7] = (int *) malloc(6 * sizeof(int));
    TBid[7][0] = 4; TB[7][0] = 2; // H2
    TBid[7][1] = 5; TB[7][1] = 6; // H2O
    TBid[7][2] = 12; TB[7][2] = 2; // CH4
    TBid[7][3] = 18; TB[7][3] = 1.5; // CO
    TBid[7][4] = 19; TB[7][4] = 2; // CO2
    TBid[7][5] = 25; TB[7][5] = 3; // C2H6

    // (8):  2 CH3 (+M) <=> C2H6 (+M)
    fwd_A[8]     = 21200000000000000;
    fwd_beta[8]  = -0.96999999999999997;
    fwd_Ea[8]    = 620;
    low_A[8]     = 1.7700000000000001e+50;
    low_beta[8]  = -9.6699999999999999;
    low_Ea[8]    = 6220;
    troe_a[8]    = 0.53249999999999997;
    troe_Tsss[8] = 151;
    troe_Ts[8]   = 1038;
    troe_Tss[8]  = 4970;
    troe_len[8]  = 4;
    prefactor_units[8]  = 1.0000000000000002e-06;
    activation_units[8] = 0.50321666580471969;
    phase_units[8]      = 1e-12;
    is_PD[8] = 1;
    nTB[8] = 6;
    TB[8] = (double *) malloc(6 * sizeof(double));
    TBid[8] = (int *) malloc(6 * sizeof(int));
    TBid[8][0] = 4; TB[8][0] = 2; // H2
    TBid[8][1] = 5; TB[8][1] = 6; // H2O
    TBid[8][2] = 12; TB[8][2] = 2; // CH4
    TBid[8][3] = 18; TB[8][3] = 1.5; // CO
    TBid[8][4] = 19; TB[8][4] = 2; // CO2
    TBid[8][5] = 25; TB[8][5] = 3; // C2H6

    // (9):  C2H3 (+M) <=> C2H2 + H (+M)
    fwd_A[9]     = 386000000;
    fwd_beta[9]  = 1.6200000000000001;
    fwd_Ea[9]    = 37048.199999999997;
    low_A[9]     = 2.5650000000000001e+27;
    low_beta[9]  = -3.3999999999999999;
    low_Ea[9]    = 35798.720000000001;
    troe_a[9]    = 1.9816;
    troe_Tsss[9] = 5383.6999999999998;
    troe_Ts[9]   = 4.2999999999999998;
    troe_Tss[9]  = -0.10000000000000001;
    troe_len[9]  = 4;
    prefactor_units[9]  = 1;
    activation_units[9] = 0.50321666580471969;
    phase_units[9]      = 1e-6;
    is_PD[9] = 1;
    nTB[9] = 8;
    TB[9] = (double *) malloc(8 * sizeof(double));
    TBid[9] = (int *) malloc(8 * sizeof(int));
    TBid[9][0] = 4; TB[9][0] = 2; // H2
    TBid[9][1] = 5; TB[9][1] = 6; // H2O
    TBid[9][2] = 12; TB[9][2] = 2; // CH4
    TBid[9][3] = 18; TB[9][3] = 1.5; // CO
    TBid[9][4] = 19; TB[9][4] = 2; // CO2
    TBid[9][5] = 21; TB[9][5] = 3; // C2H2
    TBid[9][6] = 23; TB[9][6] = 3; // C2H4
    TBid[9][7] = 25; TB[9][7] = 3; // C2H6

    // (10):  CH2CO + H (+M) <=> CH2CHO (+M)
    fwd_A[10]     = 330000000000000;
    fwd_beta[10]  = -0.059999999999999998;
    fwd_Ea[10]    = 8500;
    low_A[10]     = 3.8000000000000001e+41;
    low_beta[10]  = -7.6399999999999997;
    low_Ea[10]    = 11900;
    troe_a[10]    = 0.33700000000000002;
    troe_Tsss[10] = 1707;
    troe_Ts[10]   = 3200;
    troe_Tss[10]  = 4131;
    troe_len[10]  = 4;
    prefactor_units[10]  = 1.0000000000000002e-06;
    activation_units[10] = 0.50321666580471969;
    phase_units[10]      = 1e-12;
    is_PD[10] = 1;
    nTB[10] = 8;
    TB[10] = (double *) malloc(8 * sizeof(double));
    TBid[10] = (int *) malloc(8 * sizeof(int));
    TBid[10][0] = 4; TB[10][0] = 2; // H2
    TBid[10][1] = 5; TB[10][1] = 6; // H2O
    TBid[10][2] = 12; TB[10][2] = 2; // CH4
    TBid[10][3] = 18; TB[10][3] = 1.5; // CO
    TBid[10][4] = 19; TB[10][4] = 2; // CO2
    TBid[10][5] = 21; TB[10][5] = 3; // C2H2
    TBid[10][6] = 23; TB[10][6] = 3; // C2H4
    TBid[10][7] = 25; TB[10][7] = 3; // C2H6

    // (11):  C2H3 + H (+M) <=> C2H4 (+M)
    fwd_A[11]     = 6080000000000;
    fwd_beta[11]  = 0.27000000000000002;
    fwd_Ea[11]    = 280;
    low_A[11]     = 1.3999999999999999e+30;
    low_beta[11]  = -3.8599999999999999;
    low_Ea[11]    = 3320;
    troe_a[11]    = 0.78200000000000003;
    troe_Tsss[11] = 207.5;
    troe_Ts[11]   = 2663;
    troe_Tss[11]  = 6095;
    troe_len[11]  = 4;
    prefactor_units[11]  = 1.0000000000000002e-06;
    activation_units[11] = 0.50321666580471969;
    phase_units[11]      = 1e-12;
    is_PD[11] = 1;
    nTB[11] = 8;
    TB[11] = (double *) malloc(8 * sizeof(double));
    TBid[11] = (int *) malloc(8 * sizeof(int));
    TBid[11][0] = 4; TB[11][0] = 2; // H2
    TBid[11][1] = 5; TB[11][1] = 6; // H2O
    TBid[11][2] = 12; TB[11][2] = 2; // CH4
    TBid[11][3] = 18; TB[11][3] = 1.5; // CO
    TBid[11][4] = 19; TB[11][4] = 2; // CO2
    TBid[11][5] = 21; TB[11][5] = 3; // C2H2
    TBid[11][6] = 23; TB[11][6] = 3; // C2H4
    TBid[11][7] = 25; TB[11][7] = 3; // C2H6

    // (12):  C2H3 + CH3 (+M) <=> C3H6 (+M)
    fwd_A[12]     = 25000000000000;
    fwd_beta[12]  = 0;
    fwd_Ea[12]    = 0;
    low_A[12]     = 4.2699999999999999e+58;
    low_beta[12]  = -11.94;
    low_Ea[12]    = 9769.7999999999993;
    troe_a[12]    = 0.17499999999999999;
    troe_Tsss[12] = 1340.5999999999999;
    troe_Ts[12]   = 60000;
    troe_Tss[12]  = 10139.799999999999;
    troe_len[12]  = 4;
    prefactor_units[12]  = 1.0000000000000002e-06;
    activation_units[12] = 0.50321666580471969;
    phase_units[12]      = 1e-12;
    is_PD[12] = 1;
    nTB[12] = 7;
    TB[12] = (double *) malloc(7 * sizeof(double));
    TBid[12] = (int *) malloc(7 * sizeof(int));
    TBid[12][0] = 4; TB[12][0] = 2; // H2
    TBid[12][1] = 5; TB[12][1] = 6; // H2O
    TBid[12][2] = 12; TB[12][2] = 2; // CH4
    TBid[12][3] = 18; TB[12][3] = 1.5; // CO
    TBid[12][4] = 19; TB[12][4] = 2; // CO2
    TBid[12][5] = 23; TB[12][5] = 3; // C2H4
    TBid[12][6] = 25; TB[12][6] = 3; // C2H6

    // (13):  CH3 + CO (+M) <=> CH3CO (+M)
    fwd_A[13]     = 48500000;
    fwd_beta[13]  = 1.6499999999999999;
    fwd_Ea[13]    = 6150;
    low_A[13]     = 7.8000000000000002e+30;
    low_beta[13]  = -5.3949999999999996;
    low_Ea[13]    = 8600;
    troe_a[13]    = 0.25800000000000001;
    troe_Tsss[13] = 598;
    troe_Ts[13]   = 21002;
    troe_Tss[13]  = 1773;
    troe_len[13]  = 4;
    prefactor_units[13]  = 1.0000000000000002e-06;
    activation_units[13] = 0.50321666580471969;
    phase_units[13]      = 1e-12;
    is_PD[13] = 1;
    nTB[13] = 8;
    TB[13] = (double *) malloc(8 * sizeof(double));
    TBid[13] = (int *) malloc(8 * sizeof(int));
    TBid[13][0] = 4; TB[13][0] = 2; // H2
    TBid[13][1] = 5; TB[13][1] = 6; // H2O
    TBid[13][2] = 12; TB[13][2] = 2; // CH4
    TBid[13][3] = 18; TB[13][3] = 1.5; // CO
    TBid[13][4] = 19; TB[13][4] = 2; // CO2
    TBid[13][5] = 21; TB[13][5] = 3; // C2H2
    TBid[13][6] = 23; TB[13][6] = 3; // C2H4
    TBid[13][7] = 25; TB[13][7] = 3; // C2H6

    // (14):  CH3 + HCO (+M) <=> CH3CHO (+M)
    fwd_A[14]     = 18000000000000;
    fwd_beta[14]  = 0;
    fwd_Ea[14]    = 0;
    low_A[14]     = 2.1999999999999999e+48;
    low_beta[14]  = -9.5879999999999992;
    low_Ea[14]    = 5100;
    troe_a[14]    = 0.61729999999999996;
    troe_Tsss[14] = 13.1;
    troe_Ts[14]   = 2078;
    troe_Tss[14]  = 5093;
    troe_len[14]  = 4;
    prefactor_units[14]  = 1.0000000000000002e-06;
    activation_units[14] = 0.50321666580471969;
    phase_units[14]      = 1e-12;
    is_PD[14] = 1;
    nTB[14] = 8;
    TB[14] = (double *) malloc(8 * sizeof(double));
    TBid[14] = (int *) malloc(8 * sizeof(int));
    TBid[14][0] = 4; TB[14][0] = 2; // H2
    TBid[14][1] = 5; TB[14][1] = 6; // H2O
    TBid[14][2] = 12; TB[14][2] = 2; // CH4
    TBid[14][3] = 18; TB[14][3] = 1.5; // CO
    TBid[14][4] = 19; TB[14][4] = 2; // CO2
    TBid[14][5] = 21; TB[14][5] = 3; // C2H2
    TBid[14][6] = 23; TB[14][6] = 3; // C2H4
    TBid[14][7] = 25; TB[14][7] = 3; // C2H6

    // (15):  C2H4 (+M) <=> H2 + C2H2 (+M)
    fwd_A[15]     = 8000000000000;
    fwd_beta[15]  = 0.44;
    fwd_Ea[15]    = 88770;
    low_A[15]     = 7.0000000000000001e+50;
    low_beta[15]  = -9.3100000000000005;
    low_Ea[15]    = 99860;
    troe_a[15]    = 0.73450000000000004;
    troe_Tsss[15] = 180;
    troe_Ts[15]   = 1035;
    troe_Tss[15]  = 5417;
    troe_len[15]  = 4;
    prefactor_units[15]  = 1;
    activation_units[15] = 0.50321666580471969;
    phase_units[15]      = 1e-6;
    is_PD[15] = 1;
    nTB[15] = 6;
    TB[15] = (double *) malloc(6 * sizeof(double));
    TBid[15] = (int *) malloc(6 * sizeof(int));
    TBid[15][0] = 4; TB[15][0] = 2; // H2
    TBid[15][1] = 5; TB[15][1] = 6; // H2O
    TBid[15][2] = 12; TB[15][2] = 2; // CH4
    TBid[15][3] = 18; TB[15][3] = 1.5; // CO
    TBid[15][4] = 19; TB[15][4] = 2; // CO2
    TBid[15][5] = 25; TB[15][5] = 3; // C2H6

    // (16):  C2H4 + H (+M) <=> C2H5 (+M)
    fwd_A[16]     = 1367000000;
    fwd_beta[16]  = 1.4630000000000001;
    fwd_Ea[16]    = 1355;
    low_A[16]     = 2.0269999999999999e+39;
    low_beta[16]  = -6.6420000000000003;
    low_Ea[16]    = 5769;
    troe_a[16]    = -0.56899999999999995;
    troe_Tsss[16] = 299;
    troe_Ts[16]   = 9147;
    troe_Tss[16]  = -152.40000000000001;
    troe_len[16]  = 4;
    prefactor_units[16]  = 1.0000000000000002e-06;
    activation_units[16] = 0.50321666580471969;
    phase_units[16]      = 1e-12;
    is_PD[16] = 1;
    nTB[16] = 0;
    TB[16] = (double *) malloc(0 * sizeof(double));
    TBid[16] = (int *) malloc(0 * sizeof(int));

    // (17):  C2H5 + H (+M) <=> C2H6 (+M)
    fwd_A[17]     = 5.21e+17;
    fwd_beta[17]  = -0.98999999999999999;
    fwd_Ea[17]    = 1580;
    low_A[17]     = 1.9900000000000001e+41;
    low_beta[17]  = -7.0800000000000001;
    low_Ea[17]    = 6685;
    troe_a[17]    = 0.84219999999999995;
    troe_Tsss[17] = 125;
    troe_Ts[17]   = 2219;
    troe_Tss[17]  = 6882;
    troe_len[17]  = 4;
    prefactor_units[17]  = 1.0000000000000002e-06;
    activation_units[17] = 0.50321666580471969;
    phase_units[17]      = 1e-12;
    is_PD[17] = 1;
    nTB[17] = 6;
    TB[17] = (double *) malloc(6 * sizeof(double));
    TBid[17] = (int *) malloc(6 * sizeof(int));
    TBid[17][0] = 4; TB[17][0] = 2; // H2
    TBid[17][1] = 5; TB[17][1] = 6; // H2O
    TBid[17][2] = 12; TB[17][2] = 2; // CH4
    TBid[17][3] = 18; TB[17][3] = 1.5; // CO
    TBid[17][4] = 19; TB[17][4] = 2; // CO2
    TBid[17][5] = 25; TB[17][5] = 3; // C2H6

    // (18):  C2H5 + C2H3 (+M) <=> C4H81 (+M)
    fwd_A[18]     = 15000000000000;
    fwd_beta[18]  = 0;
    fwd_Ea[18]    = 0;
    low_A[18]     = 1.55e+56;
    low_beta[18]  = -11.789999999999999;
    low_Ea[18]    = 8984.5;
    troe_a[18]    = 0.19800000000000001;
    troe_Tsss[18] = 2277.9000000000001;
    troe_Ts[18]   = 60000;
    troe_Tss[18]  = 5723.1999999999998;
    troe_len[18]  = 4;
    prefactor_units[18]  = 1.0000000000000002e-06;
    activation_units[18] = 0.50321666580471969;
    phase_units[18]      = 1e-12;
    is_PD[18] = 1;
    nTB[18] = 6;
    TB[18] = (double *) malloc(6 * sizeof(double));
    TBid[18] = (int *) malloc(6 * sizeof(int));
    TBid[18][0] = 4; TB[18][0] = 2; // H2
    TBid[18][1] = 5; TB[18][1] = 6; // H2O
    TBid[18][2] = 12; TB[18][2] = 2; // CH4
    TBid[18][3] = 18; TB[18][3] = 1.5; // CO
    TBid[18][4] = 19; TB[18][4] = 2; // CO2
    TBid[18][5] = 25; TB[18][5] = 3; // C2H6

    // (19):  aC3H5 + H (+M) <=> C3H6 (+M)
    fwd_A[19]     = 200000000000000;
    fwd_beta[19]  = 0;
    fwd_Ea[19]    = 0;
    low_A[19]     = 1.33e+60;
    low_beta[19]  = -12;
    low_Ea[19]    = 5967.8000000000002;
    troe_a[19]    = 0.02;
    troe_Tsss[19] = 1096.5999999999999;
    troe_Ts[19]   = 1096.5999999999999;
    troe_Tss[19]  = 6859.5;
    troe_len[19]  = 4;
    prefactor_units[19]  = 1.0000000000000002e-06;
    activation_units[19] = 0.50321666580471969;
    phase_units[19]      = 1e-12;
    is_PD[19] = 1;
    nTB[19] = 6;
    TB[19] = (double *) malloc(6 * sizeof(double));
    TBid[19] = (int *) malloc(6 * sizeof(int));
    TBid[19][0] = 4; TB[19][0] = 2; // H2
    TBid[19][1] = 5; TB[19][1] = 6; // H2O
    TBid[19][2] = 12; TB[19][2] = 2; // CH4
    TBid[19][3] = 18; TB[19][3] = 1.5; // CO
    TBid[19][4] = 19; TB[19][4] = 2; // CO2
    TBid[19][5] = 25; TB[19][5] = 3; // C2H6

    // (20):  aC3H5 + CH3 (+M) <=> C4H81 (+M)
    fwd_A[20]     = 100000000000000;
    fwd_beta[20]  = -0.32000000000000001;
    fwd_Ea[20]    = -262.30000000000001;
    low_A[20]     = 3.9100000000000002e+60;
    low_beta[20]  = -12.81;
    low_Ea[20]    = 6250;
    troe_a[20]    = 0.104;
    troe_Tsss[20] = 1606;
    troe_Ts[20]   = 60000;
    troe_Tss[20]  = 6118.3999999999996;
    troe_len[20]  = 4;
    prefactor_units[20]  = 1.0000000000000002e-06;
    activation_units[20] = 0.50321666580471969;
    phase_units[20]      = 1e-12;
    is_PD[20] = 1;
    nTB[20] = 6;
    TB[20] = (double *) malloc(6 * sizeof(double));
    TBid[20] = (int *) malloc(6 * sizeof(int));
    TBid[20][0] = 4; TB[20][0] = 2; // H2
    TBid[20][1] = 5; TB[20][1] = 6; // H2O
    TBid[20][2] = 12; TB[20][2] = 2; // CH4
    TBid[20][3] = 18; TB[20][3] = 1.5; // CO
    TBid[20][4] = 19; TB[20][4] = 2; // CO2
    TBid[20][5] = 25; TB[20][5] = 3; // C2H6

    // (21):  C3H6 + H (+M) <=> nC3H7 (+M)
    fwd_A[21]     = 13300000000000;
    fwd_beta[21]  = 0;
    fwd_Ea[21]    = 3260.6999999999998;
    low_A[21]     = 6.2599999999999999e+38;
    low_beta[21]  = -6.6600000000000001;
    low_Ea[21]    = 7000;
    troe_a[21]    = 1;
    troe_Tsss[21] = 1000;
    troe_Ts[21]   = 1310;
    troe_Tss[21]  = 48097;
    troe_len[21]  = 4;
    prefactor_units[21]  = 1.0000000000000002e-06;
    activation_units[21] = 0.50321666580471969;
    phase_units[21]      = 1e-12;
    is_PD[21] = 1;
    nTB[21] = 6;
    TB[21] = (double *) malloc(6 * sizeof(double));
    TBid[21] = (int *) malloc(6 * sizeof(int));
    TBid[21][0] = 4; TB[21][0] = 2; // H2
    TBid[21][1] = 5; TB[21][1] = 6; // H2O
    TBid[21][2] = 12; TB[21][2] = 2; // CH4
    TBid[21][3] = 18; TB[21][3] = 1.5; // CO
    TBid[21][4] = 19; TB[21][4] = 2; // CO2
    TBid[21][5] = 25; TB[21][5] = 3; // C2H6

    // (22):  C3H6 + H (+M) <=> iC3H7 (+M)
    fwd_A[22]     = 13300000000000;
    fwd_beta[22]  = 0;
    fwd_Ea[22]    = 1559.8;
    low_A[22]     = 8.7e+42;
    low_beta[22]  = -7.5;
    low_Ea[22]    = 4721.8000000000002;
    troe_a[22]    = 1;
    troe_Tsss[22] = 1000;
    troe_Ts[22]   = 645.39999999999998;
    troe_Tss[22]  = 6844.3000000000002;
    troe_len[22]  = 4;
    prefactor_units[22]  = 1.0000000000000002e-06;
    activation_units[22] = 0.50321666580471969;
    phase_units[22]      = 1e-12;
    is_PD[22] = 1;
    nTB[22] = 6;
    TB[22] = (double *) malloc(6 * sizeof(double));
    TBid[22] = (int *) malloc(6 * sizeof(int));
    TBid[22][0] = 4; TB[22][0] = 2; // H2
    TBid[22][1] = 5; TB[22][1] = 6; // H2O
    TBid[22][2] = 12; TB[22][2] = 2; // CH4
    TBid[22][3] = 18; TB[22][3] = 1.5; // CO
    TBid[22][4] = 19; TB[22][4] = 2; // CO2
    TBid[22][5] = 25; TB[22][5] = 3; // C2H6

    // (23):  C6H12 + H (+M) <=> C3H6 + nC3H7 (+M)
    fwd_A[23]     = 13300000000000;
    fwd_beta[23]  = 0;
    fwd_Ea[23]    = 1559.8;
    low_A[23]     = 8.7e+42;
    low_beta[23]  = -7.5;
    low_Ea[23]    = 4721.8000000000002;
    troe_a[23]    = 1;
    troe_Tsss[23] = 1000;
    troe_Ts[23]   = 645.39999999999998;
    troe_Tss[23]  = 6844.3000000000002;
    troe_len[23]  = 4;
    prefactor_units[23]  = 1.0000000000000002e-06;
    activation_units[23] = 0.50321666580471969;
    phase_units[23]      = 1e-12;
    is_PD[23] = 1;
    nTB[23] = 6;
    TB[23] = (double *) malloc(6 * sizeof(double));
    TBid[23] = (int *) malloc(6 * sizeof(int));
    TBid[23][0] = 4; TB[23][0] = 2; // H2
    TBid[23][1] = 5; TB[23][1] = 6; // H2O
    TBid[23][2] = 12; TB[23][2] = 2; // CH4
    TBid[23][3] = 18; TB[23][3] = 1.5; // CO
    TBid[23][4] = 19; TB[23][4] = 2; // CO2
    TBid[23][5] = 25; TB[23][5] = 3; // C2H6

    // (24):  C5H10 + H (+M) <=> C3H6 + C2H5 (+M)
    fwd_A[24]     = 13300000000000;
    fwd_beta[24]  = 0;
    fwd_Ea[24]    = 1559.8;
    low_A[24]     = 8.7e+42;
    low_beta[24]  = -7.5;
    low_Ea[24]    = 4721.8000000000002;
    troe_a[24]    = 1;
    troe_Tsss[24] = 1000;
    troe_Ts[24]   = 645.39999999999998;
    troe_Tss[24]  = 6844.3000000000002;
    troe_len[24]  = 4;
    prefactor_units[24]  = 1.0000000000000002e-06;
    activation_units[24] = 0.50321666580471969;
    phase_units[24]      = 1e-12;
    is_PD[24] = 1;
    nTB[24] = 6;
    TB[24] = (double *) malloc(6 * sizeof(double));
    TBid[24] = (int *) malloc(6 * sizeof(int));
    TBid[24][0] = 4; TB[24][0] = 2; // H2
    TBid[24][1] = 5; TB[24][1] = 6; // H2O
    TBid[24][2] = 12; TB[24][2] = 2; // CH4
    TBid[24][3] = 18; TB[24][3] = 1.5; // CO
    TBid[24][4] = 19; TB[24][4] = 2; // CO2
    TBid[24][5] = 25; TB[24][5] = 3; // C2H6

    // (25):  CO + O (+M) <=> CO2 (+M)
    fwd_A[25]     = 13620000000;
    fwd_beta[25]  = 0;
    fwd_Ea[25]    = 2384;
    low_A[25]     = 1.1729999999999999e+24;
    low_beta[25]  = -2.79;
    low_Ea[25]    = 4191;
    prefactor_units[25]  = 1.0000000000000002e-06;
    activation_units[25] = 0.50321666580471969;
    phase_units[25]      = 1e-12;
    is_PD[25] = 1;
    nTB[25] = 4;
    TB[25] = (double *) malloc(4 * sizeof(double));
    TBid[25] = (int *) malloc(4 * sizeof(int));
    TBid[25][0] = 4; TB[25][0] = 2; // H2
    TBid[25][1] = 5; TB[25][1] = 12; // H2O
    TBid[25][2] = 18; TB[25][2] = 1.75; // CO
    TBid[25][3] = 19; TB[25][3] = 3.6000000000000001; // CO2

    // (26):  2 H + M <=> H2 + M
    fwd_A[26]     = 1.78e+18;
    fwd_beta[26]  = -1;
    fwd_Ea[26]    = 0;
    prefactor_units[26]  = 1.0000000000000002e-12;
    activation_units[26] = 0.50321666580471969;
    phase_units[26]      = 1e-12;
    is_PD[26] = 0;
    nTB[26] = 3;
    TB[26] = (double *) malloc(3 * sizeof(double));
    TBid[26] = (int *) malloc(3 * sizeof(int));
    TBid[26][0] = 4; TB[26][0] = 0; // H2
    TBid[26][1] = 5; TB[26][1] = 0; // H2O
    TBid[26][2] = 19; TB[26][2] = 0; // CO2

    // (27):  H + OH + M <=> H2O + M
    fwd_A[27]     = 4.4e+22;
    fwd_beta[27]  = -2;
    fwd_Ea[27]    = 0;
    prefactor_units[27]  = 1.0000000000000002e-12;
    activation_units[27] = 0.50321666580471969;
    phase_units[27]      = 1e-12;
    is_PD[27] = 0;
    nTB[27] = 4;
    TB[27] = (double *) malloc(4 * sizeof(double));
    TBid[27] = (int *) malloc(4 * sizeof(int));
    TBid[27][0] = 4; TB[27][0] = 2; // H2
    TBid[27][1] = 5; TB[27][1] = 6.2999999999999998; // H2O
    TBid[27][2] = 18; TB[27][2] = 1.75; // CO
    TBid[27][3] = 19; TB[27][3] = 3.6000000000000001; // CO2

    // (28):  O + H + M <=> OH + M
    fwd_A[28]     = 9.428e+18;
    fwd_beta[28]  = -1;
    fwd_Ea[28]    = 0;
    prefactor_units[28]  = 1.0000000000000002e-12;
    activation_units[28] = 0.50321666580471969;
    phase_units[28]      = 1e-12;
    is_PD[28] = 0;
    nTB[28] = 4;
    TB[28] = (double *) malloc(4 * sizeof(double));
    TBid[28] = (int *) malloc(4 * sizeof(int));
    TBid[28][0] = 4; TB[28][0] = 2; // H2
    TBid[28][1] = 5; TB[28][1] = 12; // H2O
    TBid[28][2] = 18; TB[28][2] = 1.75; // CO
    TBid[28][3] = 19; TB[28][3] = 3.6000000000000001; // CO2

    // (29):  2 O + M <=> O2 + M
    fwd_A[29]     = 1.2e+17;
    fwd_beta[29]  = -1;
    fwd_Ea[29]    = 0;
    prefactor_units[29]  = 1.0000000000000002e-12;
    activation_units[29] = 0.50321666580471969;
    phase_units[29]      = 1e-12;
    is_PD[29] = 0;
    nTB[29] = 4;
    TB[29] = (double *) malloc(4 * sizeof(double));
    TBid[29] = (int *) malloc(4 * sizeof(int));
    TBid[29][0] = 4; TB[29][0] = 2.3999999999999999; // H2
    TBid[29][1] = 5; TB[29][1] = 15.4; // H2O
    TBid[29][2] = 18; TB[29][2] = 1.75; // CO
    TBid[29][3] = 19; TB[29][3] = 3.6000000000000001; // CO2

    // (30):  HCO + M <=> CO + H + M
    fwd_A[30]     = 1.87e+17;
    fwd_beta[30]  = -1;
    fwd_Ea[30]    = 17000;
    prefactor_units[30]  = 1.0000000000000002e-06;
    activation_units[30] = 0.50321666580471969;
    phase_units[30]      = 1e-6;
    is_PD[30] = 0;
    nTB[30] = 4;
    TB[30] = (double *) malloc(4 * sizeof(double));
    TBid[30] = (int *) malloc(4 * sizeof(int));
    TBid[30][0] = 4; TB[30][0] = 2; // H2
    TBid[30][1] = 5; TB[30][1] = 0; // H2O
    TBid[30][2] = 18; TB[30][2] = 1.75; // CO
    TBid[30][3] = 19; TB[30][3] = 3.6000000000000001; // CO2

    // (31):  H + O2 <=> O + OH
    fwd_A[31]     = 26440000000000000;
    fwd_beta[31]  = -0.67100000000000004;
    fwd_Ea[31]    = 17041;
    prefactor_units[31]  = 1.0000000000000002e-06;
    activation_units[31] = 0.50321666580471969;
    phase_units[31]      = 1e-12;
    is_PD[31] = 0;
    nTB[31] = 0;

    // (32):  O + H2 <=> H + OH
    fwd_A[32]     = 45890;
    fwd_beta[32]  = 2.7000000000000002;
    fwd_Ea[32]    = 6260;
    prefactor_units[32]  = 1.0000000000000002e-06;
    activation_units[32] = 0.50321666580471969;
    phase_units[32]      = 1e-12;
    is_PD[32] = 0;
    nTB[32] = 0;

    // (33):  OH + H2 <=> H + H2O
    fwd_A[33]     = 173400000;
    fwd_beta[33]  = 1.51;
    fwd_Ea[33]    = 3430;
    prefactor_units[33]  = 1.0000000000000002e-06;
    activation_units[33] = 0.50321666580471969;
    phase_units[33]      = 1e-12;
    is_PD[33] = 0;
    nTB[33] = 0;

    // (34):  2 OH <=> O + H2O
    fwd_A[34]     = 39730;
    fwd_beta[34]  = 2.3999999999999999;
    fwd_Ea[34]    = -2110;
    prefactor_units[34]  = 1.0000000000000002e-06;
    activation_units[34] = 0.50321666580471969;
    phase_units[34]      = 1e-12;
    is_PD[34] = 0;
    nTB[34] = 0;

    // (35):  2 H + H2O <=> H2 + H2O
    fwd_A[35]     = 5.624e+19;
    fwd_beta[35]  = -1.25;
    fwd_Ea[35]    = 0;
    prefactor_units[35]  = 1.0000000000000002e-12;
    activation_units[35] = 0.50321666580471969;
    phase_units[35]      = 1e-18;
    is_PD[35] = 0;
    nTB[35] = 0;

    // (36):  H2 + O2 <=> HO2 + H
    fwd_A[36]     = 591600;
    fwd_beta[36]  = 2.4329999999999998;
    fwd_Ea[36]    = 53502;
    prefactor_units[36]  = 1.0000000000000002e-06;
    activation_units[36] = 0.50321666580471969;
    phase_units[36]      = 1e-12;
    is_PD[36] = 0;
    nTB[36] = 0;

    // (37):  HO2 + H <=> O + H2O
    fwd_A[37]     = 3970000000000;
    fwd_beta[37]  = 0;
    fwd_Ea[37]    = 671;
    prefactor_units[37]  = 1.0000000000000002e-06;
    activation_units[37] = 0.50321666580471969;
    phase_units[37]      = 1e-12;
    is_PD[37] = 0;
    nTB[37] = 0;

    // (38):  HO2 + H <=> 2 OH
    fwd_A[38]     = 74850000000000;
    fwd_beta[38]  = 0;
    fwd_Ea[38]    = 295;
    prefactor_units[38]  = 1.0000000000000002e-06;
    activation_units[38] = 0.50321666580471969;
    phase_units[38]      = 1e-12;
    is_PD[38] = 0;
    nTB[38] = 0;

    // (39):  HO2 + O <=> OH + O2
    fwd_A[39]     = 40000000000000;
    fwd_beta[39]  = 0;
    fwd_Ea[39]    = 0;
    prefactor_units[39]  = 1.0000000000000002e-06;
    activation_units[39] = 0.50321666580471969;
    phase_units[39]      = 1e-12;
    is_PD[39] = 0;
    nTB[39] = 0;

    // (40):  HO2 + OH <=> O2 + H2O
    fwd_A[40]     = 23750000000000;
    fwd_beta[40]  = 0;
    fwd_Ea[40]    = -500;
    prefactor_units[40]  = 1.0000000000000002e-06;
    activation_units[40] = 0.50321666580471969;
    phase_units[40]      = 1e-12;
    is_PD[40] = 0;
    nTB[40] = 0;

    // (41):  HO2 + OH <=> O2 + H2O
    fwd_A[41]     = 10000000000000000;
    fwd_beta[41]  = 0;
    fwd_Ea[41]    = 17330;
    prefactor_units[41]  = 1.0000000000000002e-06;
    activation_units[41] = 0.50321666580471969;
    phase_units[41]      = 1e-12;
    is_PD[41] = 0;
    nTB[41] = 0;

    // (42):  2 HO2 <=> O2 + H2O2
    fwd_A[42]     = 130000000000;
    fwd_beta[42]  = 0;
    fwd_Ea[42]    = -1630;
    prefactor_units[42]  = 1.0000000000000002e-06;
    activation_units[42] = 0.50321666580471969;
    phase_units[42]      = 1e-12;
    is_PD[42] = 0;
    nTB[42] = 0;

    // (43):  2 HO2 <=> O2 + H2O2
    fwd_A[43]     = 365800000000000;
    fwd_beta[43]  = 0;
    fwd_Ea[43]    = 12000;
    prefactor_units[43]  = 1.0000000000000002e-06;
    activation_units[43] = 0.50321666580471969;
    phase_units[43]      = 1e-12;
    is_PD[43] = 0;
    nTB[43] = 0;

    // (44):  H2O2 + H <=> HO2 + H2
    fwd_A[44]     = 6050000;
    fwd_beta[44]  = 2;
    fwd_Ea[44]    = 5200;
    prefactor_units[44]  = 1.0000000000000002e-06;
    activation_units[44] = 0.50321666580471969;
    phase_units[44]      = 1e-12;
    is_PD[44] = 0;
    nTB[44] = 0;

    // (45):  H2O2 + H <=> OH + H2O
    fwd_A[45]     = 24100000000000;
    fwd_beta[45]  = 0;
    fwd_Ea[45]    = 3970;
    prefactor_units[45]  = 1.0000000000000002e-06;
    activation_units[45] = 0.50321666580471969;
    phase_units[45]      = 1e-12;
    is_PD[45] = 0;
    nTB[45] = 0;

    // (46):  H2O2 + O <=> OH + HO2
    fwd_A[46]     = 9630000;
    fwd_beta[46]  = 2;
    fwd_Ea[46]    = 3970;
    prefactor_units[46]  = 1.0000000000000002e-06;
    activation_units[46] = 0.50321666580471969;
    phase_units[46]      = 1e-12;
    is_PD[46] = 0;
    nTB[46] = 0;

    // (47):  H2O2 + OH <=> HO2 + H2O
    fwd_A[47]     = 2000000000000;
    fwd_beta[47]  = 0;
    fwd_Ea[47]    = 427;
    prefactor_units[47]  = 1.0000000000000002e-06;
    activation_units[47] = 0.50321666580471969;
    phase_units[47]      = 1e-12;
    is_PD[47] = 0;
    nTB[47] = 0;

    // (48):  H2O2 + OH <=> HO2 + H2O
    fwd_A[48]     = 2.6700000000000001e+41;
    fwd_beta[48]  = -7;
    fwd_Ea[48]    = 37600;
    prefactor_units[48]  = 1.0000000000000002e-06;
    activation_units[48] = 0.50321666580471969;
    phase_units[48]      = 1e-12;
    is_PD[48] = 0;
    nTB[48] = 0;

    // (49):  CO + OH <=> CO2 + H
    fwd_A[49]     = 800000000000;
    fwd_beta[49]  = 0.14000000000000001;
    fwd_Ea[49]    = 7352;
    prefactor_units[49]  = 1.0000000000000002e-06;
    activation_units[49] = 0.50321666580471969;
    phase_units[49]      = 1e-12;
    is_PD[49] = 0;
    nTB[49] = 0;

    // (50):  CO + OH <=> CO2 + H
    fwd_A[50]     = 87840000000;
    fwd_beta[50]  = 0.029999999999999999;
    fwd_Ea[50]    = -16;
    prefactor_units[50]  = 1.0000000000000002e-06;
    activation_units[50] = 0.50321666580471969;
    phase_units[50]      = 1e-12;
    is_PD[50] = 0;
    nTB[50] = 0;

    // (51):  CO + HO2 <=> CO2 + OH
    fwd_A[51]     = 30100000000000;
    fwd_beta[51]  = 0;
    fwd_Ea[51]    = 23000;
    prefactor_units[51]  = 1.0000000000000002e-06;
    activation_units[51] = 0.50321666580471969;
    phase_units[51]      = 1e-12;
    is_PD[51] = 0;
    nTB[51] = 0;

    // (52):  HCO + H <=> CO + H2
    fwd_A[52]     = 120000000000000;
    fwd_beta[52]  = 0;
    fwd_Ea[52]    = 0;
    prefactor_units[52]  = 1.0000000000000002e-06;
    activation_units[52] = 0.50321666580471969;
    phase_units[52]      = 1e-12;
    is_PD[52] = 0;
    nTB[52] = 0;

    // (53):  HCO + O <=> CO + OH
    fwd_A[53]     = 30000000000000;
    fwd_beta[53]  = 0;
    fwd_Ea[53]    = 0;
    prefactor_units[53]  = 1.0000000000000002e-06;
    activation_units[53] = 0.50321666580471969;
    phase_units[53]      = 1e-12;
    is_PD[53] = 0;
    nTB[53] = 0;

    // (54):  HCO + O <=> CO2 + H
    fwd_A[54]     = 30000000000000;
    fwd_beta[54]  = 0;
    fwd_Ea[54]    = 0;
    prefactor_units[54]  = 1.0000000000000002e-06;
    activation_units[54] = 0.50321666580471969;
    phase_units[54]      = 1e-12;
    is_PD[54] = 0;
    nTB[54] = 0;

    // (55):  HCO + OH <=> CO + H2O
    fwd_A[55]     = 30200000000000;
    fwd_beta[55]  = 0;
    fwd_Ea[55]    = 0;
    prefactor_units[55]  = 1.0000000000000002e-06;
    activation_units[55] = 0.50321666580471969;
    phase_units[55]      = 1e-12;
    is_PD[55] = 0;
    nTB[55] = 0;

    // (56):  HCO + H2O <=> CO + H + H2O
    fwd_A[56]     = 2.244e+18;
    fwd_beta[56]  = -1;
    fwd_Ea[56]    = 17000;
    prefactor_units[56]  = 1.0000000000000002e-06;
    activation_units[56] = 0.50321666580471969;
    phase_units[56]      = 1e-12;
    is_PD[56] = 0;
    nTB[56] = 0;

    // (57):  HCO + O2 <=> CO + HO2
    fwd_A[57]     = 12040000000;
    fwd_beta[57]  = 0.80700000000000005;
    fwd_Ea[57]    = -727;
    prefactor_units[57]  = 1.0000000000000002e-06;
    activation_units[57] = 0.50321666580471969;
    phase_units[57]      = 1e-12;
    is_PD[57] = 0;
    nTB[57] = 0;

    // (58):  CH + O <=> CO + H
    fwd_A[58]     = 57000000000000;
    fwd_beta[58]  = 0;
    fwd_Ea[58]    = 0;
    prefactor_units[58]  = 1.0000000000000002e-06;
    activation_units[58] = 0.50321666580471969;
    phase_units[58]      = 1e-12;
    is_PD[58] = 0;
    nTB[58] = 0;

    // (59):  CH + OH <=> HCO + H
    fwd_A[59]     = 30000000000000;
    fwd_beta[59]  = 0;
    fwd_Ea[59]    = 0;
    prefactor_units[59]  = 1.0000000000000002e-06;
    activation_units[59] = 0.50321666580471969;
    phase_units[59]      = 1e-12;
    is_PD[59] = 0;
    nTB[59] = 0;

    // (60):  CH + H2 <=> CH2 + H
    fwd_A[60]     = 110700000;
    fwd_beta[60]  = 1.79;
    fwd_Ea[60]    = 1670;
    prefactor_units[60]  = 1.0000000000000002e-06;
    activation_units[60] = 0.50321666580471969;
    phase_units[60]      = 1e-12;
    is_PD[60] = 0;
    nTB[60] = 0;

    // (61):  CH + H2O <=> CH2O + H
    fwd_A[61]     = 5710000000000;
    fwd_beta[61]  = 0;
    fwd_Ea[61]    = -755;
    prefactor_units[61]  = 1.0000000000000002e-06;
    activation_units[61] = 0.50321666580471969;
    phase_units[61]      = 1e-12;
    is_PD[61] = 0;
    nTB[61] = 0;

    // (62):  CH + O2 <=> HCO + O
    fwd_A[62]     = 33000000000000;
    fwd_beta[62]  = 0;
    fwd_Ea[62]    = 0;
    prefactor_units[62]  = 1.0000000000000002e-06;
    activation_units[62] = 0.50321666580471969;
    phase_units[62]      = 1e-12;
    is_PD[62] = 0;
    nTB[62] = 0;

    // (63):  CH + CO2 <=> HCO + CO
    fwd_A[63]     = 3400000000000;
    fwd_beta[63]  = 0;
    fwd_Ea[63]    = 690;
    prefactor_units[63]  = 1.0000000000000002e-06;
    activation_units[63] = 0.50321666580471969;
    phase_units[63]      = 1e-12;
    is_PD[63] = 0;
    nTB[63] = 0;

    // (64):  CH2 + O <=> HCO + H
    fwd_A[64]     = 80000000000000;
    fwd_beta[64]  = 0;
    fwd_Ea[64]    = 0;
    prefactor_units[64]  = 1.0000000000000002e-06;
    activation_units[64] = 0.50321666580471969;
    phase_units[64]      = 1e-12;
    is_PD[64] = 0;
    nTB[64] = 0;

    // (65):  CH2 + OH <=> CH2O + H
    fwd_A[65]     = 20000000000000;
    fwd_beta[65]  = 0;
    fwd_Ea[65]    = 0;
    prefactor_units[65]  = 1.0000000000000002e-06;
    activation_units[65] = 0.50321666580471969;
    phase_units[65]      = 1e-12;
    is_PD[65] = 0;
    nTB[65] = 0;

    // (66):  CH2 + OH <=> CH + H2O
    fwd_A[66]     = 11300000;
    fwd_beta[66]  = 2;
    fwd_Ea[66]    = 3000;
    prefactor_units[66]  = 1.0000000000000002e-06;
    activation_units[66] = 0.50321666580471969;
    phase_units[66]      = 1e-12;
    is_PD[66] = 0;
    nTB[66] = 0;

    // (67):  CH2 + H2 <=> H + CH3
    fwd_A[67]     = 500000;
    fwd_beta[67]  = 2;
    fwd_Ea[67]    = 7230;
    prefactor_units[67]  = 1.0000000000000002e-06;
    activation_units[67] = 0.50321666580471969;
    phase_units[67]      = 1e-12;
    is_PD[67] = 0;
    nTB[67] = 0;

    // (68):  CH2 + O2 <=> HCO + OH
    fwd_A[68]     = 10600000000000;
    fwd_beta[68]  = 0;
    fwd_Ea[68]    = 1500;
    prefactor_units[68]  = 1.0000000000000002e-06;
    activation_units[68] = 0.50321666580471969;
    phase_units[68]      = 1e-12;
    is_PD[68] = 0;
    nTB[68] = 0;

    // (69):  CH2 + O2 <=> CO2 + 2 H
    fwd_A[69]     = 2640000000000;
    fwd_beta[69]  = 0;
    fwd_Ea[69]    = 1500;
    prefactor_units[69]  = 1.0000000000000002e-06;
    activation_units[69] = 0.50321666580471969;
    phase_units[69]      = 1e-12;
    is_PD[69] = 0;
    nTB[69] = 0;

    // (70):  CH2 + HO2 <=> CH2O + OH
    fwd_A[70]     = 20000000000000;
    fwd_beta[70]  = 0;
    fwd_Ea[70]    = 0;
    prefactor_units[70]  = 1.0000000000000002e-06;
    activation_units[70] = 0.50321666580471969;
    phase_units[70]      = 1e-12;
    is_PD[70] = 0;
    nTB[70] = 0;

    // (71):  CH2* + N2 <=> CH2 + N2
    fwd_A[71]     = 15000000000000;
    fwd_beta[71]  = 0;
    fwd_Ea[71]    = 600;
    prefactor_units[71]  = 1.0000000000000002e-06;
    activation_units[71] = 0.50321666580471969;
    phase_units[71]      = 1e-12;
    is_PD[71] = 0;
    nTB[71] = 0;

    // (72):  CH2* + H <=> CH + H2
    fwd_A[72]     = 30000000000000;
    fwd_beta[72]  = 0;
    fwd_Ea[72]    = 0;
    prefactor_units[72]  = 1.0000000000000002e-06;
    activation_units[72] = 0.50321666580471969;
    phase_units[72]      = 1e-12;
    is_PD[72] = 0;
    nTB[72] = 0;

    // (73):  CH2* + OH <=> CH2O + H
    fwd_A[73]     = 30000000000000;
    fwd_beta[73]  = 0;
    fwd_Ea[73]    = 0;
    prefactor_units[73]  = 1.0000000000000002e-06;
    activation_units[73] = 0.50321666580471969;
    phase_units[73]      = 1e-12;
    is_PD[73] = 0;
    nTB[73] = 0;

    // (74):  CH2* + H2 <=> CH3 + H
    fwd_A[74]     = 70000000000000;
    fwd_beta[74]  = 0;
    fwd_Ea[74]    = 0;
    prefactor_units[74]  = 1.0000000000000002e-06;
    activation_units[74] = 0.50321666580471969;
    phase_units[74]      = 1e-12;
    is_PD[74] = 0;
    nTB[74] = 0;

    // (75):  CH2* + O2 <=> H + OH + CO
    fwd_A[75]     = 28000000000000;
    fwd_beta[75]  = 0;
    fwd_Ea[75]    = 0;
    prefactor_units[75]  = 1.0000000000000002e-06;
    activation_units[75] = 0.50321666580471969;
    phase_units[75]      = 1e-12;
    is_PD[75] = 0;
    nTB[75] = 0;

    // (76):  CH2* + O2 <=> CO + H2O
    fwd_A[76]     = 12000000000000;
    fwd_beta[76]  = 0;
    fwd_Ea[76]    = 0;
    prefactor_units[76]  = 1.0000000000000002e-06;
    activation_units[76] = 0.50321666580471969;
    phase_units[76]      = 1e-12;
    is_PD[76] = 0;
    nTB[76] = 0;

    // (77):  CH2* + H2O <=> CH2 + H2O
    fwd_A[77]     = 30000000000000;
    fwd_beta[77]  = 0;
    fwd_Ea[77]    = 0;
    prefactor_units[77]  = 1.0000000000000002e-06;
    activation_units[77] = 0.50321666580471969;
    phase_units[77]      = 1e-12;
    is_PD[77] = 0;
    nTB[77] = 0;

    // (78):  CH2* + CO <=> CH2 + CO
    fwd_A[78]     = 9000000000000;
    fwd_beta[78]  = 0;
    fwd_Ea[78]    = 0;
    prefactor_units[78]  = 1.0000000000000002e-06;
    activation_units[78] = 0.50321666580471969;
    phase_units[78]      = 1e-12;
    is_PD[78] = 0;
    nTB[78] = 0;

    // (79):  CH2* + CO2 <=> CH2 + CO2
    fwd_A[79]     = 7000000000000;
    fwd_beta[79]  = 0;
    fwd_Ea[79]    = 0;
    prefactor_units[79]  = 1.0000000000000002e-06;
    activation_units[79] = 0.50321666580471969;
    phase_units[79]      = 1e-12;
    is_PD[79] = 0;
    nTB[79] = 0;

    // (80):  CH2* + CO2 <=> CH2O + CO
    fwd_A[80]     = 14000000000000;
    fwd_beta[80]  = 0;
    fwd_Ea[80]    = 0;
    prefactor_units[80]  = 1.0000000000000002e-06;
    activation_units[80] = 0.50321666580471969;
    phase_units[80]      = 1e-12;
    is_PD[80] = 0;
    nTB[80] = 0;

    // (81):  CH2O + H <=> HCO + H2
    fwd_A[81]     = 23000000000;
    fwd_beta[81]  = 1.05;
    fwd_Ea[81]    = 3275;
    prefactor_units[81]  = 1.0000000000000002e-06;
    activation_units[81] = 0.50321666580471969;
    phase_units[81]      = 1e-12;
    is_PD[81] = 0;
    nTB[81] = 0;

    // (82):  CH2O + O <=> HCO + OH
    fwd_A[82]     = 39000000000000;
    fwd_beta[82]  = 0;
    fwd_Ea[82]    = 3540;
    prefactor_units[82]  = 1.0000000000000002e-06;
    activation_units[82] = 0.50321666580471969;
    phase_units[82]      = 1e-12;
    is_PD[82] = 0;
    nTB[82] = 0;

    // (83):  CH2O + OH <=> HCO + H2O
    fwd_A[83]     = 3430000000;
    fwd_beta[83]  = 1.1799999999999999;
    fwd_Ea[83]    = -447;
    prefactor_units[83]  = 1.0000000000000002e-06;
    activation_units[83] = 0.50321666580471969;
    phase_units[83]      = 1e-12;
    is_PD[83] = 0;
    nTB[83] = 0;

    // (84):  CH2O + O2 <=> HCO + HO2
    fwd_A[84]     = 100000000000000;
    fwd_beta[84]  = 0;
    fwd_Ea[84]    = 40000;
    prefactor_units[84]  = 1.0000000000000002e-06;
    activation_units[84] = 0.50321666580471969;
    phase_units[84]      = 1e-12;
    is_PD[84] = 0;
    nTB[84] = 0;

    // (85):  CH2O + HO2 <=> HCO + H2O2
    fwd_A[85]     = 1000000000000;
    fwd_beta[85]  = 0;
    fwd_Ea[85]    = 8000;
    prefactor_units[85]  = 1.0000000000000002e-06;
    activation_units[85] = 0.50321666580471969;
    phase_units[85]      = 1e-12;
    is_PD[85] = 0;
    nTB[85] = 0;

    // (86):  CH2O + CH <=> CH2CO + H
    fwd_A[86]     = 94600000000000;
    fwd_beta[86]  = 0;
    fwd_Ea[86]    = -515;
    prefactor_units[86]  = 1.0000000000000002e-06;
    activation_units[86] = 0.50321666580471969;
    phase_units[86]      = 1e-12;
    is_PD[86] = 0;
    nTB[86] = 0;

    // (87):  CH3 + O <=> CH2O + H
    fwd_A[87]     = 84300000000000;
    fwd_beta[87]  = 0;
    fwd_Ea[87]    = 0;
    prefactor_units[87]  = 1.0000000000000002e-06;
    activation_units[87] = 0.50321666580471969;
    phase_units[87]      = 1e-12;
    is_PD[87] = 0;
    nTB[87] = 0;

    // (88):  CH3 + OH <=> CH2 + H2O
    fwd_A[88]     = 56000000;
    fwd_beta[88]  = 1.6000000000000001;
    fwd_Ea[88]    = 5420;
    prefactor_units[88]  = 1.0000000000000002e-06;
    activation_units[88] = 0.50321666580471969;
    phase_units[88]      = 1e-12;
    is_PD[88] = 0;
    nTB[88] = 0;

    // (89):  CH3 + OH <=> CH2* + H2O
    fwd_A[89]     = 25010000000000;
    fwd_beta[89]  = 0;
    fwd_Ea[89]    = 0;
    prefactor_units[89]  = 1.0000000000000002e-06;
    activation_units[89] = 0.50321666580471969;
    phase_units[89]      = 1e-12;
    is_PD[89] = 0;
    nTB[89] = 0;

    // (90):  CH3 + O2 <=> O + CH3O
    fwd_A[90]     = 30830000000000;
    fwd_beta[90]  = 0;
    fwd_Ea[90]    = 28800;
    prefactor_units[90]  = 1.0000000000000002e-06;
    activation_units[90] = 0.50321666580471969;
    phase_units[90]      = 1e-12;
    is_PD[90] = 0;
    nTB[90] = 0;

    // (91):  CH3 + O2 <=> OH + CH2O
    fwd_A[91]     = 36000000000;
    fwd_beta[91]  = 0;
    fwd_Ea[91]    = 8940;
    prefactor_units[91]  = 1.0000000000000002e-06;
    activation_units[91] = 0.50321666580471969;
    phase_units[91]      = 1e-12;
    is_PD[91] = 0;
    nTB[91] = 0;

    // (92):  CH3 + HO2 <=> CH4 + O2
    fwd_A[92]     = 1000000000000;
    fwd_beta[92]  = 0;
    fwd_Ea[92]    = 0;
    prefactor_units[92]  = 1.0000000000000002e-06;
    activation_units[92] = 0.50321666580471969;
    phase_units[92]      = 1e-12;
    is_PD[92] = 0;
    nTB[92] = 0;

    // (93):  CH3 + HO2 <=> CH3O + OH
    fwd_A[93]     = 13400000000000;
    fwd_beta[93]  = 0;
    fwd_Ea[93]    = 0;
    prefactor_units[93]  = 1.0000000000000002e-06;
    activation_units[93] = 0.50321666580471969;
    phase_units[93]      = 1e-12;
    is_PD[93] = 0;
    nTB[93] = 0;

    // (94):  CH3 + CH <=> C2H3 + H
    fwd_A[94]     = 30000000000000;
    fwd_beta[94]  = 0;
    fwd_Ea[94]    = 0;
    prefactor_units[94]  = 1.0000000000000002e-06;
    activation_units[94] = 0.50321666580471969;
    phase_units[94]      = 1e-12;
    is_PD[94] = 0;
    nTB[94] = 0;

    // (95):  CH3 + HCO <=> CH4 + CO
    fwd_A[95]     = 8480000000000;
    fwd_beta[95]  = 0;
    fwd_Ea[95]    = 0;
    prefactor_units[95]  = 1.0000000000000002e-06;
    activation_units[95] = 0.50321666580471969;
    phase_units[95]      = 1e-12;
    is_PD[95] = 0;
    nTB[95] = 0;

    // (96):  CH3 + CH2O <=> CH4 + HCO
    fwd_A[96]     = 3320;
    fwd_beta[96]  = 2.8100000000000001;
    fwd_Ea[96]    = 5860;
    prefactor_units[96]  = 1.0000000000000002e-06;
    activation_units[96] = 0.50321666580471969;
    phase_units[96]      = 1e-12;
    is_PD[96] = 0;
    nTB[96] = 0;

    // (97):  CH3 + CH2 <=> C2H4 + H
    fwd_A[97]     = 40000000000000;
    fwd_beta[97]  = 0;
    fwd_Ea[97]    = 0;
    prefactor_units[97]  = 1.0000000000000002e-06;
    activation_units[97] = 0.50321666580471969;
    phase_units[97]      = 1e-12;
    is_PD[97] = 0;
    nTB[97] = 0;

    // (98):  2 CH3 <=> H + C2H5
    fwd_A[98]     = 4990000000000;
    fwd_beta[98]  = 0.10000000000000001;
    fwd_Ea[98]    = 10600;
    prefactor_units[98]  = 1.0000000000000002e-06;
    activation_units[98] = 0.50321666580471969;
    phase_units[98]      = 1e-12;
    is_PD[98] = 0;
    nTB[98] = 0;

    // (99):  CH3 + HCCO <=> C2H4 + CO
    fwd_A[99]     = 50000000000000;
    fwd_beta[99]  = 0;
    fwd_Ea[99]    = 0;
    prefactor_units[99]  = 1.0000000000000002e-06;
    activation_units[99] = 0.50321666580471969;
    phase_units[99]      = 1e-12;
    is_PD[99] = 0;
    nTB[99] = 0;

    // (100):  CH3O + H <=> CH2O + H2
    fwd_A[100]     = 20000000000000;
    fwd_beta[100]  = 0;
    fwd_Ea[100]    = 0;
    prefactor_units[100]  = 1.0000000000000002e-06;
    activation_units[100] = 0.50321666580471969;
    phase_units[100]      = 1e-12;
    is_PD[100] = 0;
    nTB[100] = 0;

    // (101):  CH3O + H <=> CH3 + OH
    fwd_A[101]     = 32000000000000;
    fwd_beta[101]  = 0;
    fwd_Ea[101]    = 0;
    prefactor_units[101]  = 1.0000000000000002e-06;
    activation_units[101] = 0.50321666580471969;
    phase_units[101]      = 1e-12;
    is_PD[101] = 0;
    nTB[101] = 0;

    // (102):  CH3O + H <=> CH2* + H2O
    fwd_A[102]     = 16000000000000;
    fwd_beta[102]  = 0;
    fwd_Ea[102]    = 0;
    prefactor_units[102]  = 1.0000000000000002e-06;
    activation_units[102] = 0.50321666580471969;
    phase_units[102]      = 1e-12;
    is_PD[102] = 0;
    nTB[102] = 0;

    // (103):  CH3O + OH <=> CH2O + H2O
    fwd_A[103]     = 5000000000000;
    fwd_beta[103]  = 0;
    fwd_Ea[103]    = 0;
    prefactor_units[103]  = 1.0000000000000002e-06;
    activation_units[103] = 0.50321666580471969;
    phase_units[103]      = 1e-12;
    is_PD[103] = 0;
    nTB[103] = 0;

    // (104):  CH3O + O2 <=> CH2O + HO2
    fwd_A[104]     = 4.2799999999999999e-13;
    fwd_beta[104]  = 7.5999999999999996;
    fwd_Ea[104]    = -3530;
    prefactor_units[104]  = 1.0000000000000002e-06;
    activation_units[104] = 0.50321666580471969;
    phase_units[104]      = 1e-12;
    is_PD[104] = 0;
    nTB[104] = 0;

    // (105):  CH2OH + H <=> CH2O + H2
    fwd_A[105]     = 20000000000000;
    fwd_beta[105]  = 0;
    fwd_Ea[105]    = 0;
    prefactor_units[105]  = 1.0000000000000002e-06;
    activation_units[105] = 0.50321666580471969;
    phase_units[105]      = 1e-12;
    is_PD[105] = 0;
    nTB[105] = 0;

    // (106):  CH2OH + H <=> CH3 + OH
    fwd_A[106]     = 12000000000000;
    fwd_beta[106]  = 0;
    fwd_Ea[106]    = 0;
    prefactor_units[106]  = 1.0000000000000002e-06;
    activation_units[106] = 0.50321666580471969;
    phase_units[106]      = 1e-12;
    is_PD[106] = 0;
    nTB[106] = 0;

    // (107):  CH2OH + H <=> CH2* + H2O
    fwd_A[107]     = 6000000000000;
    fwd_beta[107]  = 0;
    fwd_Ea[107]    = 0;
    prefactor_units[107]  = 1.0000000000000002e-06;
    activation_units[107] = 0.50321666580471969;
    phase_units[107]      = 1e-12;
    is_PD[107] = 0;
    nTB[107] = 0;

    // (108):  CH2OH + O2 <=> CH2O + HO2
    fwd_A[108]     = 18000000000000;
    fwd_beta[108]  = 0;
    fwd_Ea[108]    = 900;
    prefactor_units[108]  = 1.0000000000000002e-06;
    activation_units[108] = 0.50321666580471969;
    phase_units[108]      = 1e-12;
    is_PD[108] = 0;
    nTB[108] = 0;

    // (109):  CH4 + H <=> CH3 + H2
    fwd_A[109]     = 660000000;
    fwd_beta[109]  = 1.6200000000000001;
    fwd_Ea[109]    = 10840;
    prefactor_units[109]  = 1.0000000000000002e-06;
    activation_units[109] = 0.50321666580471969;
    phase_units[109]      = 1e-12;
    is_PD[109] = 0;
    nTB[109] = 0;

    // (110):  CH4 + O <=> CH3 + OH
    fwd_A[110]     = 1020000000;
    fwd_beta[110]  = 1.5;
    fwd_Ea[110]    = 8600;
    prefactor_units[110]  = 1.0000000000000002e-06;
    activation_units[110] = 0.50321666580471969;
    phase_units[110]      = 1e-12;
    is_PD[110] = 0;
    nTB[110] = 0;

    // (111):  CH4 + OH <=> CH3 + H2O
    fwd_A[111]     = 100000000;
    fwd_beta[111]  = 1.6000000000000001;
    fwd_Ea[111]    = 3120;
    prefactor_units[111]  = 1.0000000000000002e-06;
    activation_units[111] = 0.50321666580471969;
    phase_units[111]      = 1e-12;
    is_PD[111] = 0;
    nTB[111] = 0;

    // (112):  CH4 + CH <=> C2H4 + H
    fwd_A[112]     = 60000000000000;
    fwd_beta[112]  = 0;
    fwd_Ea[112]    = 0;
    prefactor_units[112]  = 1.0000000000000002e-06;
    activation_units[112] = 0.50321666580471969;
    phase_units[112]      = 1e-12;
    is_PD[112] = 0;
    nTB[112] = 0;

    // (113):  CH4 + CH2 <=> 2 CH3
    fwd_A[113]     = 2460000;
    fwd_beta[113]  = 2;
    fwd_Ea[113]    = 8270;
    prefactor_units[113]  = 1.0000000000000002e-06;
    activation_units[113] = 0.50321666580471969;
    phase_units[113]      = 1e-12;
    is_PD[113] = 0;
    nTB[113] = 0;

    // (114):  CH4 + CH2* <=> 2 CH3
    fwd_A[114]     = 16000000000000;
    fwd_beta[114]  = 0;
    fwd_Ea[114]    = -570;
    prefactor_units[114]  = 1.0000000000000002e-06;
    activation_units[114] = 0.50321666580471969;
    phase_units[114]      = 1e-12;
    is_PD[114] = 0;
    nTB[114] = 0;

    // (115):  CH3OH + H <=> CH2OH + H2
    fwd_A[115]     = 17000000;
    fwd_beta[115]  = 2.1000000000000001;
    fwd_Ea[115]    = 4870;
    prefactor_units[115]  = 1.0000000000000002e-06;
    activation_units[115] = 0.50321666580471969;
    phase_units[115]      = 1e-12;
    is_PD[115] = 0;
    nTB[115] = 0;

    // (116):  CH3OH + H <=> CH3O + H2
    fwd_A[116]     = 4200000;
    fwd_beta[116]  = 2.1000000000000001;
    fwd_Ea[116]    = 4870;
    prefactor_units[116]  = 1.0000000000000002e-06;
    activation_units[116] = 0.50321666580471969;
    phase_units[116]      = 1e-12;
    is_PD[116] = 0;
    nTB[116] = 0;

    // (117):  CH3OH + O <=> CH2OH + OH
    fwd_A[117]     = 388000;
    fwd_beta[117]  = 2.5;
    fwd_Ea[117]    = 3100;
    prefactor_units[117]  = 1.0000000000000002e-06;
    activation_units[117] = 0.50321666580471969;
    phase_units[117]      = 1e-12;
    is_PD[117] = 0;
    nTB[117] = 0;

    // (118):  CH3OH + OH <=> CH2OH + H2O
    fwd_A[118]     = 1440000;
    fwd_beta[118]  = 2;
    fwd_Ea[118]    = -840;
    prefactor_units[118]  = 1.0000000000000002e-06;
    activation_units[118] = 0.50321666580471969;
    phase_units[118]      = 1e-12;
    is_PD[118] = 0;
    nTB[118] = 0;

    // (119):  CH3OH + OH <=> CH3O + H2O
    fwd_A[119]     = 6300000;
    fwd_beta[119]  = 2;
    fwd_Ea[119]    = 1500;
    prefactor_units[119]  = 1.0000000000000002e-06;
    activation_units[119] = 0.50321666580471969;
    phase_units[119]      = 1e-12;
    is_PD[119] = 0;
    nTB[119] = 0;

    // (120):  C2H + O <=> CH + CO
    fwd_A[120]     = 50000000000000;
    fwd_beta[120]  = 0;
    fwd_Ea[120]    = 0;
    prefactor_units[120]  = 1.0000000000000002e-06;
    activation_units[120] = 0.50321666580471969;
    phase_units[120]      = 1e-12;
    is_PD[120] = 0;
    nTB[120] = 0;

    // (121):  C2H + OH <=> H + HCCO
    fwd_A[121]     = 20000000000000;
    fwd_beta[121]  = 0;
    fwd_Ea[121]    = 0;
    prefactor_units[121]  = 1.0000000000000002e-06;
    activation_units[121] = 0.50321666580471969;
    phase_units[121]      = 1e-12;
    is_PD[121] = 0;
    nTB[121] = 0;

    // (122):  C2H + O2 <=> HCO + CO
    fwd_A[122]     = 50000000000000;
    fwd_beta[122]  = 0;
    fwd_Ea[122]    = 1500;
    prefactor_units[122]  = 1.0000000000000002e-06;
    activation_units[122] = 0.50321666580471969;
    phase_units[122]      = 1e-12;
    is_PD[122] = 0;
    nTB[122] = 0;

    // (123):  C2H + H2 <=> H + C2H2
    fwd_A[123]     = 490000;
    fwd_beta[123]  = 2.5;
    fwd_Ea[123]    = 560;
    prefactor_units[123]  = 1.0000000000000002e-06;
    activation_units[123] = 0.50321666580471969;
    phase_units[123]      = 1e-12;
    is_PD[123] = 0;
    nTB[123] = 0;

    // (124):  HCCO + H <=> CH2* + CO
    fwd_A[124]     = 100000000000000;
    fwd_beta[124]  = 0;
    fwd_Ea[124]    = 0;
    prefactor_units[124]  = 1.0000000000000002e-06;
    activation_units[124] = 0.50321666580471969;
    phase_units[124]      = 1e-12;
    is_PD[124] = 0;
    nTB[124] = 0;

    // (125):  HCCO + O <=> H + 2 CO
    fwd_A[125]     = 100000000000000;
    fwd_beta[125]  = 0;
    fwd_Ea[125]    = 0;
    prefactor_units[125]  = 1.0000000000000002e-06;
    activation_units[125] = 0.50321666580471969;
    phase_units[125]      = 1e-12;
    is_PD[125] = 0;
    nTB[125] = 0;

    // (126):  HCCO + O2 <=> OH + 2 CO
    fwd_A[126]     = 1600000000000;
    fwd_beta[126]  = 0;
    fwd_Ea[126]    = 854;
    prefactor_units[126]  = 1.0000000000000002e-06;
    activation_units[126] = 0.50321666580471969;
    phase_units[126]      = 1e-12;
    is_PD[126] = 0;
    nTB[126] = 0;

    // (127):  C2H2 + O <=> C2H + OH
    fwd_A[127]     = 4.6e+19;
    fwd_beta[127]  = -1.4099999999999999;
    fwd_Ea[127]    = 28950;
    prefactor_units[127]  = 1.0000000000000002e-06;
    activation_units[127] = 0.50321666580471969;
    phase_units[127]      = 1e-12;
    is_PD[127] = 0;
    nTB[127] = 0;

    // (128):  C2H2 + O <=> CH2 + CO
    fwd_A[128]     = 4080000;
    fwd_beta[128]  = 2;
    fwd_Ea[128]    = 1900;
    prefactor_units[128]  = 1.0000000000000002e-06;
    activation_units[128] = 0.50321666580471969;
    phase_units[128]      = 1e-12;
    is_PD[128] = 0;
    nTB[128] = 0;

    // (129):  C2H2 + O <=> HCCO + H
    fwd_A[129]     = 16320000;
    fwd_beta[129]  = 2;
    fwd_Ea[129]    = 1900;
    prefactor_units[129]  = 1.0000000000000002e-06;
    activation_units[129] = 0.50321666580471969;
    phase_units[129]      = 1e-12;
    is_PD[129] = 0;
    nTB[129] = 0;

    // (130):  C2H2 + OH <=> CH2CO + H
    fwd_A[130]     = 0.00021800000000000001;
    fwd_beta[130]  = 4.5;
    fwd_Ea[130]    = -1000;
    prefactor_units[130]  = 1.0000000000000002e-06;
    activation_units[130] = 0.50321666580471969;
    phase_units[130]      = 1e-12;
    is_PD[130] = 0;
    nTB[130] = 0;

    // (131):  C2H2 + OH <=> CH2CO + H
    fwd_A[131]     = 504000;
    fwd_beta[131]  = 2.2999999999999998;
    fwd_Ea[131]    = 13500;
    prefactor_units[131]  = 1.0000000000000002e-06;
    activation_units[131] = 0.50321666580471969;
    phase_units[131]      = 1e-12;
    is_PD[131] = 0;
    nTB[131] = 0;

    // (132):  C2H2 + OH <=> C2H + H2O
    fwd_A[132]     = 33700000;
    fwd_beta[132]  = 2;
    fwd_Ea[132]    = 14000;
    prefactor_units[132]  = 1.0000000000000002e-06;
    activation_units[132] = 0.50321666580471969;
    phase_units[132]      = 1e-12;
    is_PD[132] = 0;
    nTB[132] = 0;

    // (133):  C2H2 + HCO <=> C2H3 + CO
    fwd_A[133]     = 10000000;
    fwd_beta[133]  = 2;
    fwd_Ea[133]    = 6000;
    prefactor_units[133]  = 1.0000000000000002e-06;
    activation_units[133] = 0.50321666580471969;
    phase_units[133]      = 1e-12;
    is_PD[133] = 0;
    nTB[133] = 0;

    // (134):  C2H2 + CH2 <=> C3H3 + H
    fwd_A[134]     = 12000000000000;
    fwd_beta[134]  = 0;
    fwd_Ea[134]    = 6620;
    prefactor_units[134]  = 1.0000000000000002e-06;
    activation_units[134] = 0.50321666580471969;
    phase_units[134]      = 1e-12;
    is_PD[134] = 0;
    nTB[134] = 0;

    // (135):  C2H2 + CH2* <=> C3H3 + H
    fwd_A[135]     = 20000000000000;
    fwd_beta[135]  = 0;
    fwd_Ea[135]    = 0;
    prefactor_units[135]  = 1.0000000000000002e-06;
    activation_units[135] = 0.50321666580471969;
    phase_units[135]      = 1e-12;
    is_PD[135] = 0;
    nTB[135] = 0;

    // (136):  C2H2 + C2H <=> C4H2 + H
    fwd_A[136]     = 96000000000000;
    fwd_beta[136]  = 0;
    fwd_Ea[136]    = 0;
    prefactor_units[136]  = 1.0000000000000002e-06;
    activation_units[136] = 0.50321666580471969;
    phase_units[136]      = 1e-12;
    is_PD[136] = 0;
    nTB[136] = 0;

    // (137):  C2H2 + CH3 <=> pC3H4 + H
    fwd_A[137]     = 2560000000;
    fwd_beta[137]  = 1.1000000000000001;
    fwd_Ea[137]    = 13644;
    prefactor_units[137]  = 1.0000000000000002e-06;
    activation_units[137] = 0.50321666580471969;
    phase_units[137]      = 1e-12;
    is_PD[137] = 0;
    nTB[137] = 0;

    // (138):  C2H2 + CH3 <=> aC3H4 + H
    fwd_A[138]     = 5140000000;
    fwd_beta[138]  = 0.85999999999999999;
    fwd_Ea[138]    = 22153;
    prefactor_units[138]  = 1.0000000000000002e-06;
    activation_units[138] = 0.50321666580471969;
    phase_units[138]      = 1e-12;
    is_PD[138] = 0;
    nTB[138] = 0;

    // (139):  CH2CO + H <=> HCCO + H2
    fwd_A[139]     = 50000000000000;
    fwd_beta[139]  = 0;
    fwd_Ea[139]    = 8000;
    prefactor_units[139]  = 1.0000000000000002e-06;
    activation_units[139] = 0.50321666580471969;
    phase_units[139]      = 1e-12;
    is_PD[139] = 0;
    nTB[139] = 0;

    // (140):  CH2CO + H <=> CH3 + CO
    fwd_A[140]     = 1500000000;
    fwd_beta[140]  = 1.4299999999999999;
    fwd_Ea[140]    = 2690;
    prefactor_units[140]  = 1.0000000000000002e-06;
    activation_units[140] = 0.50321666580471969;
    phase_units[140]      = 1e-12;
    is_PD[140] = 0;
    nTB[140] = 0;

    // (141):  CH2CO + OH <=> HCCO + H2O
    fwd_A[141]     = 7500000000000;
    fwd_beta[141]  = 0;
    fwd_Ea[141]    = 2000;
    prefactor_units[141]  = 1.0000000000000002e-06;
    activation_units[141] = 0.50321666580471969;
    phase_units[141]      = 1e-12;
    is_PD[141] = 0;
    nTB[141] = 0;

    // (142):  C2H3 + H <=> C2H2 + H2
    fwd_A[142]     = 90000000000000;
    fwd_beta[142]  = 0;
    fwd_Ea[142]    = 0;
    prefactor_units[142]  = 1.0000000000000002e-06;
    activation_units[142] = 0.50321666580471969;
    phase_units[142]      = 1e-12;
    is_PD[142] = 0;
    nTB[142] = 0;

    // (143):  C2H3 + O <=> CH2CO + H
    fwd_A[143]     = 48000000000000;
    fwd_beta[143]  = 0;
    fwd_Ea[143]    = 0;
    prefactor_units[143]  = 1.0000000000000002e-06;
    activation_units[143] = 0.50321666580471969;
    phase_units[143]      = 1e-12;
    is_PD[143] = 0;
    nTB[143] = 0;

    // (144):  C2H3 + O <=> CH3 + CO
    fwd_A[144]     = 48000000000000;
    fwd_beta[144]  = 0;
    fwd_Ea[144]    = 0;
    prefactor_units[144]  = 1.0000000000000002e-06;
    activation_units[144] = 0.50321666580471969;
    phase_units[144]      = 1e-12;
    is_PD[144] = 0;
    nTB[144] = 0;

    // (145):  C2H3 + OH <=> C2H2 + H2O
    fwd_A[145]     = 30110000000000;
    fwd_beta[145]  = 0;
    fwd_Ea[145]    = 0;
    prefactor_units[145]  = 1.0000000000000002e-06;
    activation_units[145] = 0.50321666580471969;
    phase_units[145]      = 1e-12;
    is_PD[145] = 0;
    nTB[145] = 0;

    // (146):  C2H3 + O2 <=> C2H2 + HO2
    fwd_A[146]     = 1340000;
    fwd_beta[146]  = 1.6100000000000001;
    fwd_Ea[146]    = -383.39999999999998;
    prefactor_units[146]  = 1.0000000000000002e-06;
    activation_units[146] = 0.50321666580471969;
    phase_units[146]      = 1e-12;
    is_PD[146] = 0;
    nTB[146] = 0;

    // (147):  C2H3 + O2 <=> CH2CHO + O
    fwd_A[147]     = 300000000000;
    fwd_beta[147]  = 0.28999999999999998;
    fwd_Ea[147]    = 11;
    prefactor_units[147]  = 1.0000000000000002e-06;
    activation_units[147] = 0.50321666580471969;
    phase_units[147]      = 1e-12;
    is_PD[147] = 0;
    nTB[147] = 0;

    // (148):  C2H3 + O2 <=> HCO + CH2O
    fwd_A[148]     = 46000000000000000;
    fwd_beta[148]  = -1.3899999999999999;
    fwd_Ea[148]    = 1010;
    prefactor_units[148]  = 1.0000000000000002e-06;
    activation_units[148] = 0.50321666580471969;
    phase_units[148]      = 1e-12;
    is_PD[148] = 0;
    nTB[148] = 0;

    // (149):  C2H3 + HO2 <=> CH2CHO + OH
    fwd_A[149]     = 10000000000000;
    fwd_beta[149]  = 0;
    fwd_Ea[149]    = 0;
    prefactor_units[149]  = 1.0000000000000002e-06;
    activation_units[149] = 0.50321666580471969;
    phase_units[149]      = 1e-12;
    is_PD[149] = 0;
    nTB[149] = 0;

    // (150):  C2H3 + HCO <=> C2H4 + CO
    fwd_A[150]     = 90330000000000;
    fwd_beta[150]  = 0;
    fwd_Ea[150]    = 0;
    prefactor_units[150]  = 1.0000000000000002e-06;
    activation_units[150] = 0.50321666580471969;
    phase_units[150]      = 1e-12;
    is_PD[150] = 0;
    nTB[150] = 0;

    // (151):  C2H3 + HCO <=> C2H3CHO
    fwd_A[151]     = 18000000000000;
    fwd_beta[151]  = 0;
    fwd_Ea[151]    = 0;
    prefactor_units[151]  = 1.0000000000000002e-06;
    activation_units[151] = 0.50321666580471969;
    phase_units[151]      = 1e-12;
    is_PD[151] = 0;
    nTB[151] = 0;

    // (152):  C2H3 + CH3 <=> aC3H5 + H
    fwd_A[152]     = 1.5e+24;
    fwd_beta[152]  = -2.8300000000000001;
    fwd_Ea[152]    = 18618;
    prefactor_units[152]  = 1.0000000000000002e-06;
    activation_units[152] = 0.50321666580471969;
    phase_units[152]      = 1e-12;
    is_PD[152] = 0;
    nTB[152] = 0;

    // (153):  CH2CHO <=> CH3 + CO
    fwd_A[153]     = 7.7999999999999994e+41;
    fwd_beta[153]  = -9.1470000000000002;
    fwd_Ea[153]    = 46900;
    prefactor_units[153]  = 1;
    activation_units[153] = 0.50321666580471969;
    phase_units[153]      = 1e-6;
    is_PD[153] = 0;
    nTB[153] = 0;

    // (154):  CH2CHO + H <=> CH3CO + H
    fwd_A[154]     = 5000000000000;
    fwd_beta[154]  = 0;
    fwd_Ea[154]    = 0;
    prefactor_units[154]  = 1.0000000000000002e-06;
    activation_units[154] = 0.50321666580471969;
    phase_units[154]      = 1e-12;
    is_PD[154] = 0;
    nTB[154] = 0;

    // (155):  CH2CHO + H <=> CH3 + HCO
    fwd_A[155]     = 90000000000000;
    fwd_beta[155]  = 0;
    fwd_Ea[155]    = 0;
    prefactor_units[155]  = 1.0000000000000002e-06;
    activation_units[155] = 0.50321666580471969;
    phase_units[155]      = 1e-12;
    is_PD[155] = 0;
    nTB[155] = 0;

    // (156):  CH2CHO + H <=> CH2CO + H2
    fwd_A[156]     = 20000000000000;
    fwd_beta[156]  = 0;
    fwd_Ea[156]    = 4000;
    prefactor_units[156]  = 1.0000000000000002e-06;
    activation_units[156] = 0.50321666580471969;
    phase_units[156]      = 1e-12;
    is_PD[156] = 0;
    nTB[156] = 0;

    // (157):  CH2CHO + OH <=> CH2CO + H2O
    fwd_A[157]     = 10000000000000;
    fwd_beta[157]  = 0;
    fwd_Ea[157]    = 2000;
    prefactor_units[157]  = 1.0000000000000002e-06;
    activation_units[157] = 0.50321666580471969;
    phase_units[157]      = 1e-12;
    is_PD[157] = 0;
    nTB[157] = 0;

    // (158):  CH2CHO + O2 <=> CH2CO + HO2
    fwd_A[158]     = 140000000000;
    fwd_beta[158]  = 0;
    fwd_Ea[158]    = 0;
    prefactor_units[158]  = 1.0000000000000002e-06;
    activation_units[158] = 0.50321666580471969;
    phase_units[158]      = 1e-12;
    is_PD[158] = 0;
    nTB[158] = 0;

    // (159):  CH2CHO + O2 <=> CH2O + CO + OH
    fwd_A[159]     = 18000000000;
    fwd_beta[159]  = 0;
    fwd_Ea[159]    = 0;
    prefactor_units[159]  = 1.0000000000000002e-06;
    activation_units[159] = 0.50321666580471969;
    phase_units[159]      = 1e-12;
    is_PD[159] = 0;
    nTB[159] = 0;

    // (160):  CH3CO + H <=> CH3 + HCO
    fwd_A[160]     = 96000000000000;
    fwd_beta[160]  = 0;
    fwd_Ea[160]    = 0;
    prefactor_units[160]  = 1.0000000000000002e-06;
    activation_units[160] = 0.50321666580471969;
    phase_units[160]      = 1e-12;
    is_PD[160] = 0;
    nTB[160] = 0;

    // (161):  CH3CO + HO2 <=> CH3 + CO2 + OH
    fwd_A[161]     = 30000000000000;
    fwd_beta[161]  = 0;
    fwd_Ea[161]    = 0;
    prefactor_units[161]  = 1.0000000000000002e-06;
    activation_units[161] = 0.50321666580471969;
    phase_units[161]      = 1e-12;
    is_PD[161] = 0;
    nTB[161] = 0;

    // (162):  CH3CHO + H <=> CH3CO + H2
    fwd_A[162]     = 4100000000;
    fwd_beta[162]  = 1.1599999999999999;
    fwd_Ea[162]    = 2400;
    prefactor_units[162]  = 1.0000000000000002e-06;
    activation_units[162] = 0.50321666580471969;
    phase_units[162]      = 1e-12;
    is_PD[162] = 0;
    nTB[162] = 0;

    // (163):  CH3CHO + OH <=> CH3CO + H2O
    fwd_A[163]     = 23500000000;
    fwd_beta[163]  = 0.72999999999999998;
    fwd_Ea[163]    = -1110;
    prefactor_units[163]  = 1.0000000000000002e-06;
    activation_units[163] = 0.50321666580471969;
    phase_units[163]      = 1e-12;
    is_PD[163] = 0;
    nTB[163] = 0;

    // (164):  CH3CHO + CH3 <=> CH3CO + CH4
    fwd_A[164]     = 1.9999999999999999e-06;
    fwd_beta[164]  = 5.5999999999999996;
    fwd_Ea[164]    = 2460;
    prefactor_units[164]  = 1.0000000000000002e-06;
    activation_units[164] = 0.50321666580471969;
    phase_units[164]      = 1e-12;
    is_PD[164] = 0;
    nTB[164] = 0;

    // (165):  CH3CHO + O2 <=> CH3CO + HO2
    fwd_A[165]     = 30000000000000;
    fwd_beta[165]  = 0;
    fwd_Ea[165]    = 39100;
    prefactor_units[165]  = 1.0000000000000002e-06;
    activation_units[165] = 0.50321666580471969;
    phase_units[165]      = 1e-12;
    is_PD[165] = 0;
    nTB[165] = 0;

    // (166):  C2H4 + H <=> C2H3 + H2
    fwd_A[166]     = 50700000;
    fwd_beta[166]  = 1.8999999999999999;
    fwd_Ea[166]    = 12950;
    prefactor_units[166]  = 1.0000000000000002e-06;
    activation_units[166] = 0.50321666580471969;
    phase_units[166]      = 1e-12;
    is_PD[166] = 0;
    nTB[166] = 0;

    // (167):  C2H4 + O <=> C2H3 + OH
    fwd_A[167]     = 15100000;
    fwd_beta[167]  = 1.8999999999999999;
    fwd_Ea[167]    = 3740;
    prefactor_units[167]  = 1.0000000000000002e-06;
    activation_units[167] = 0.50321666580471969;
    phase_units[167]      = 1e-12;
    is_PD[167] = 0;
    nTB[167] = 0;

    // (168):  C2H4 + O <=> CH3 + HCO
    fwd_A[168]     = 19200000;
    fwd_beta[168]  = 1.8300000000000001;
    fwd_Ea[168]    = 220;
    prefactor_units[168]  = 1.0000000000000002e-06;
    activation_units[168] = 0.50321666580471969;
    phase_units[168]      = 1e-12;
    is_PD[168] = 0;
    nTB[168] = 0;

    // (169):  C2H4 + O <=> CH2 + CH2O
    fwd_A[169]     = 384000;
    fwd_beta[169]  = 1.8300000000000001;
    fwd_Ea[169]    = 220;
    prefactor_units[169]  = 1.0000000000000002e-06;
    activation_units[169] = 0.50321666580471969;
    phase_units[169]      = 1e-12;
    is_PD[169] = 0;
    nTB[169] = 0;

    // (170):  C2H4 + OH <=> C2H3 + H2O
    fwd_A[170]     = 3600000;
    fwd_beta[170]  = 2;
    fwd_Ea[170]    = 2500;
    prefactor_units[170]  = 1.0000000000000002e-06;
    activation_units[170] = 0.50321666580471969;
    phase_units[170]      = 1e-12;
    is_PD[170] = 0;
    nTB[170] = 0;

    // (171):  C2H4 + HCO <=> C2H5 + CO
    fwd_A[171]     = 10000000;
    fwd_beta[171]  = 2;
    fwd_Ea[171]    = 8000;
    prefactor_units[171]  = 1.0000000000000002e-06;
    activation_units[171] = 0.50321666580471969;
    phase_units[171]      = 1e-12;
    is_PD[171] = 0;
    nTB[171] = 0;

    // (172):  C2H4 + CH <=> aC3H4 + H
    fwd_A[172]     = 30000000000000;
    fwd_beta[172]  = 0;
    fwd_Ea[172]    = 0;
    prefactor_units[172]  = 1.0000000000000002e-06;
    activation_units[172] = 0.50321666580471969;
    phase_units[172]      = 1e-12;
    is_PD[172] = 0;
    nTB[172] = 0;

    // (173):  C2H4 + CH <=> pC3H4 + H
    fwd_A[173]     = 30000000000000;
    fwd_beta[173]  = 0;
    fwd_Ea[173]    = 0;
    prefactor_units[173]  = 1.0000000000000002e-06;
    activation_units[173] = 0.50321666580471969;
    phase_units[173]      = 1e-12;
    is_PD[173] = 0;
    nTB[173] = 0;

    // (174):  C2H4 + CH2 <=> aC3H5 + H
    fwd_A[174]     = 20000000000000;
    fwd_beta[174]  = 0;
    fwd_Ea[174]    = 6000;
    prefactor_units[174]  = 1.0000000000000002e-06;
    activation_units[174] = 0.50321666580471969;
    phase_units[174]      = 1e-12;
    is_PD[174] = 0;
    nTB[174] = 0;

    // (175):  C2H4 + CH2* <=> aC3H5 + H
    fwd_A[175]     = 50000000000000;
    fwd_beta[175]  = 0;
    fwd_Ea[175]    = 0;
    prefactor_units[175]  = 1.0000000000000002e-06;
    activation_units[175] = 0.50321666580471969;
    phase_units[175]      = 1e-12;
    is_PD[175] = 0;
    nTB[175] = 0;

    // (176):  C2H4 + CH3 <=> C2H3 + CH4
    fwd_A[176]     = 227000;
    fwd_beta[176]  = 2;
    fwd_Ea[176]    = 9200;
    prefactor_units[176]  = 1.0000000000000002e-06;
    activation_units[176] = 0.50321666580471969;
    phase_units[176]      = 1e-12;
    is_PD[176] = 0;
    nTB[176] = 0;

    // (177):  C2H4 + CH3 <=> nC3H7
    fwd_A[177]     = 330000000000;
    fwd_beta[177]  = 0;
    fwd_Ea[177]    = 7700;
    prefactor_units[177]  = 1.0000000000000002e-06;
    activation_units[177] = 0.50321666580471969;
    phase_units[177]      = 1e-12;
    is_PD[177] = 0;
    nTB[177] = 0;

    // (178):  C2H4 + C2H3 <=> C4H7
    fwd_A[178]     = 7.9300000000000007e+38;
    fwd_beta[178]  = -8.4700000000000006;
    fwd_Ea[178]    = 14220;
    prefactor_units[178]  = 1.0000000000000002e-06;
    activation_units[178] = 0.50321666580471969;
    phase_units[178]      = 1e-12;
    is_PD[178] = 0;
    nTB[178] = 0;

    // (179):  C2H5 + H <=> C2H4 + H2
    fwd_A[179]     = 2000000000000;
    fwd_beta[179]  = 0;
    fwd_Ea[179]    = 0;
    prefactor_units[179]  = 1.0000000000000002e-06;
    activation_units[179] = 0.50321666580471969;
    phase_units[179]      = 1e-12;
    is_PD[179] = 0;
    nTB[179] = 0;

    // (180):  C2H5 + O <=> CH3 + CH2O
    fwd_A[180]     = 16040000000000;
    fwd_beta[180]  = 0;
    fwd_Ea[180]    = 0;
    prefactor_units[180]  = 1.0000000000000002e-06;
    activation_units[180] = 0.50321666580471969;
    phase_units[180]      = 1e-12;
    is_PD[180] = 0;
    nTB[180] = 0;

    // (181):  C2H5 + O <=> CH3CHO + H
    fwd_A[181]     = 80200000000000;
    fwd_beta[181]  = 0;
    fwd_Ea[181]    = 0;
    prefactor_units[181]  = 1.0000000000000002e-06;
    activation_units[181] = 0.50321666580471969;
    phase_units[181]      = 1e-12;
    is_PD[181] = 0;
    nTB[181] = 0;

    // (182):  C2H5 + O2 <=> C2H4 + HO2
    fwd_A[182]     = 20000000000;
    fwd_beta[182]  = 0;
    fwd_Ea[182]    = 0;
    prefactor_units[182]  = 1.0000000000000002e-06;
    activation_units[182] = 0.50321666580471969;
    phase_units[182]      = 1e-12;
    is_PD[182] = 0;
    nTB[182] = 0;

    // (183):  C2H5 + HO2 <=> C2H6 + O2
    fwd_A[183]     = 300000000000;
    fwd_beta[183]  = 0;
    fwd_Ea[183]    = 0;
    prefactor_units[183]  = 1.0000000000000002e-06;
    activation_units[183] = 0.50321666580471969;
    phase_units[183]      = 1e-12;
    is_PD[183] = 0;
    nTB[183] = 0;

    // (184):  C2H5 + HO2 <=> C2H4 + H2O2
    fwd_A[184]     = 300000000000;
    fwd_beta[184]  = 0;
    fwd_Ea[184]    = 0;
    prefactor_units[184]  = 1.0000000000000002e-06;
    activation_units[184] = 0.50321666580471969;
    phase_units[184]      = 1e-12;
    is_PD[184] = 0;
    nTB[184] = 0;

    // (185):  C2H5 + HO2 <=> CH3 + CH2O + OH
    fwd_A[185]     = 24000000000000;
    fwd_beta[185]  = 0;
    fwd_Ea[185]    = 0;
    prefactor_units[185]  = 1.0000000000000002e-06;
    activation_units[185] = 0.50321666580471969;
    phase_units[185]      = 1e-12;
    is_PD[185] = 0;
    nTB[185] = 0;

    // (186):  C2H5 + C2H3 <=> aC3H5 + CH3
    fwd_A[186]     = 3.8999999999999999e+32;
    fwd_beta[186]  = -5.2199999999999998;
    fwd_Ea[186]    = 19747;
    prefactor_units[186]  = 1.0000000000000002e-06;
    activation_units[186] = 0.50321666580471969;
    phase_units[186]      = 1e-12;
    is_PD[186] = 0;
    nTB[186] = 0;

    // (187):  C2H6 + H <=> C2H5 + H2
    fwd_A[187]     = 115000000;
    fwd_beta[187]  = 1.8999999999999999;
    fwd_Ea[187]    = 7530;
    prefactor_units[187]  = 1.0000000000000002e-06;
    activation_units[187] = 0.50321666580471969;
    phase_units[187]      = 1e-12;
    is_PD[187] = 0;
    nTB[187] = 0;

    // (188):  C2H6 + O <=> C2H5 + OH
    fwd_A[188]     = 89800000;
    fwd_beta[188]  = 1.9199999999999999;
    fwd_Ea[188]    = 5690;
    prefactor_units[188]  = 1.0000000000000002e-06;
    activation_units[188] = 0.50321666580471969;
    phase_units[188]      = 1e-12;
    is_PD[188] = 0;
    nTB[188] = 0;

    // (189):  C2H6 + OH <=> C2H5 + H2O
    fwd_A[189]     = 3540000;
    fwd_beta[189]  = 2.1200000000000001;
    fwd_Ea[189]    = 870;
    prefactor_units[189]  = 1.0000000000000002e-06;
    activation_units[189] = 0.50321666580471969;
    phase_units[189]      = 1e-12;
    is_PD[189] = 0;
    nTB[189] = 0;

    // (190):  C2H6 + CH2* <=> C2H5 + CH3
    fwd_A[190]     = 40000000000000;
    fwd_beta[190]  = 0;
    fwd_Ea[190]    = -550;
    prefactor_units[190]  = 1.0000000000000002e-06;
    activation_units[190] = 0.50321666580471969;
    phase_units[190]      = 1e-12;
    is_PD[190] = 0;
    nTB[190] = 0;

    // (191):  C2H6 + CH3 <=> C2H5 + CH4
    fwd_A[191]     = 6140000;
    fwd_beta[191]  = 1.74;
    fwd_Ea[191]    = 10450;
    prefactor_units[191]  = 1.0000000000000002e-06;
    activation_units[191] = 0.50321666580471969;
    phase_units[191]      = 1e-12;
    is_PD[191] = 0;
    nTB[191] = 0;

    // (192):  C3H3 + H <=> pC3H4
    fwd_A[192]     = 15000000000000;
    fwd_beta[192]  = 0;
    fwd_Ea[192]    = 0;
    prefactor_units[192]  = 1.0000000000000002e-06;
    activation_units[192] = 0.50321666580471969;
    phase_units[192]      = 1e-12;
    is_PD[192] = 0;
    nTB[192] = 0;

    // (193):  C3H3 + O <=> CH2O + C2H
    fwd_A[193]     = 20000000000000;
    fwd_beta[193]  = 0;
    fwd_Ea[193]    = 0;
    prefactor_units[193]  = 1.0000000000000002e-06;
    activation_units[193] = 0.50321666580471969;
    phase_units[193]      = 1e-12;
    is_PD[193] = 0;
    nTB[193] = 0;

    // (194):  C3H3 + O2 <=> CH2CO + HCO
    fwd_A[194]     = 30000000000;
    fwd_beta[194]  = 0;
    fwd_Ea[194]    = 2868;
    prefactor_units[194]  = 1.0000000000000002e-06;
    activation_units[194] = 0.50321666580471969;
    phase_units[194]      = 1e-12;
    is_PD[194] = 0;
    nTB[194] = 0;

    // (195):  C3H3 + HO2 <=> pC3H4 + O2
    fwd_A[195]     = 2500000000000;
    fwd_beta[195]  = 0;
    fwd_Ea[195]    = 0;
    prefactor_units[195]  = 1.0000000000000002e-06;
    activation_units[195] = 0.50321666580471969;
    phase_units[195]      = 1e-12;
    is_PD[195] = 0;
    nTB[195] = 0;

    // (196):  aC3H4 + H <=> CH3CCH2
    fwd_A[196]     = 9.4599999999999998e+42;
    fwd_beta[196]  = -9.4299999999999997;
    fwd_Ea[196]    = 11190;
    prefactor_units[196]  = 1.0000000000000002e-06;
    activation_units[196] = 0.50321666580471969;
    phase_units[196]      = 1e-12;
    is_PD[196] = 0;
    nTB[196] = 0;

    // (197):  aC3H4 + H <=> aC3H5
    fwd_A[197]     = 1.5199999999999999e+59;
    fwd_beta[197]  = -13.539999999999999;
    fwd_Ea[197]    = 26949;
    prefactor_units[197]  = 1.0000000000000002e-06;
    activation_units[197] = 0.50321666580471969;
    phase_units[197]      = 1e-12;
    is_PD[197] = 0;
    nTB[197] = 0;

    // (198):  aC3H4 + O <=> C2H4 + CO
    fwd_A[198]     = 20000000;
    fwd_beta[198]  = 1.8;
    fwd_Ea[198]    = 1000;
    prefactor_units[198]  = 1.0000000000000002e-06;
    activation_units[198] = 0.50321666580471969;
    phase_units[198]      = 1e-12;
    is_PD[198] = 0;
    nTB[198] = 0;

    // (199):  aC3H4 + OH <=> C3H3 + H2O
    fwd_A[199]     = 5300000;
    fwd_beta[199]  = 2;
    fwd_Ea[199]    = 2000;
    prefactor_units[199]  = 1.0000000000000002e-06;
    activation_units[199] = 0.50321666580471969;
    phase_units[199]      = 1e-12;
    is_PD[199] = 0;
    nTB[199] = 0;

    // (200):  pC3H4 <=> aC3H4
    fwd_A[200]     = 5.1500000000000001e+60;
    fwd_beta[200]  = -13.93;
    fwd_Ea[200]    = 91117;
    prefactor_units[200]  = 1;
    activation_units[200] = 0.50321666580471969;
    phase_units[200]      = 1e-6;
    is_PD[200] = 0;
    nTB[200] = 0;

    // (201):  pC3H4 + H <=> aC3H4 + H
    fwd_A[201]     = 6.27e+17;
    fwd_beta[201]  = -0.91000000000000003;
    fwd_Ea[201]    = 10079;
    prefactor_units[201]  = 1.0000000000000002e-06;
    activation_units[201] = 0.50321666580471969;
    phase_units[201]      = 1e-12;
    is_PD[201] = 0;
    nTB[201] = 0;

    // (202):  pC3H4 + H <=> CH3CCH2
    fwd_A[202]     = 1.6599999999999999e+47;
    fwd_beta[202]  = -10.58;
    fwd_Ea[202]    = 13690;
    prefactor_units[202]  = 1.0000000000000002e-06;
    activation_units[202] = 0.50321666580471969;
    phase_units[202]      = 1e-12;
    is_PD[202] = 0;
    nTB[202] = 0;

    // (203):  pC3H4 + O <=> C2H4 + CO
    fwd_A[203]     = 10000000000000;
    fwd_beta[203]  = 0;
    fwd_Ea[203]    = 2250;
    prefactor_units[203]  = 1.0000000000000002e-06;
    activation_units[203] = 0.50321666580471969;
    phase_units[203]      = 1e-12;
    is_PD[203] = 0;
    nTB[203] = 0;

    // (204):  pC3H4 + OH <=> C3H3 + H2O
    fwd_A[204]     = 1000000;
    fwd_beta[204]  = 2;
    fwd_Ea[204]    = 100;
    prefactor_units[204]  = 1.0000000000000002e-06;
    activation_units[204] = 0.50321666580471969;
    phase_units[204]      = 1e-12;
    is_PD[204] = 0;
    nTB[204] = 0;

    // (205):  aC3H5 + H <=> aC3H4 + H2
    fwd_A[205]     = 18000000000000;
    fwd_beta[205]  = 0;
    fwd_Ea[205]    = 0;
    prefactor_units[205]  = 1.0000000000000002e-06;
    activation_units[205] = 0.50321666580471969;
    phase_units[205]      = 1e-12;
    is_PD[205] = 0;
    nTB[205] = 0;

    // (206):  aC3H5 + O <=> C2H3CHO + H
    fwd_A[206]     = 60000000000000;
    fwd_beta[206]  = 0;
    fwd_Ea[206]    = 0;
    prefactor_units[206]  = 1.0000000000000002e-06;
    activation_units[206] = 0.50321666580471969;
    phase_units[206]      = 1e-12;
    is_PD[206] = 0;
    nTB[206] = 0;

    // (207):  aC3H5 + OH <=> C2H3CHO + 2 H
    fwd_A[207]     = 4.2000000000000001e+32;
    fwd_beta[207]  = -5.1600000000000001;
    fwd_Ea[207]    = 30126;
    prefactor_units[207]  = 1.0000000000000002e-06;
    activation_units[207] = 0.50321666580471969;
    phase_units[207]      = 1e-12;
    is_PD[207] = 0;
    nTB[207] = 0;

    // (208):  aC3H5 + OH <=> aC3H4 + H2O
    fwd_A[208]     = 6000000000000;
    fwd_beta[208]  = 0;
    fwd_Ea[208]    = 0;
    prefactor_units[208]  = 1.0000000000000002e-06;
    activation_units[208] = 0.50321666580471969;
    phase_units[208]      = 1e-12;
    is_PD[208] = 0;
    nTB[208] = 0;

    // (209):  aC3H5 + HO2 <=> C3H6 + O2
    fwd_A[209]     = 2660000000000;
    fwd_beta[209]  = 0;
    fwd_Ea[209]    = 0;
    prefactor_units[209]  = 1.0000000000000002e-06;
    activation_units[209] = 0.50321666580471969;
    phase_units[209]      = 1e-12;
    is_PD[209] = 0;
    nTB[209] = 0;

    // (210):  aC3H5 + HO2 <=> OH + C2H3 + CH2O
    fwd_A[210]     = 6600000000000;
    fwd_beta[210]  = 0;
    fwd_Ea[210]    = 0;
    prefactor_units[210]  = 1.0000000000000002e-06;
    activation_units[210] = 0.50321666580471969;
    phase_units[210]      = 1e-12;
    is_PD[210] = 0;
    nTB[210] = 0;

    // (211):  aC3H5 + HCO <=> C3H6 + CO
    fwd_A[211]     = 60000000000000;
    fwd_beta[211]  = 0;
    fwd_Ea[211]    = 0;
    prefactor_units[211]  = 1.0000000000000002e-06;
    activation_units[211] = 0.50321666580471969;
    phase_units[211]      = 1e-12;
    is_PD[211] = 0;
    nTB[211] = 0;

    // (212):  aC3H5 + CH3 <=> aC3H4 + CH4
    fwd_A[212]     = 3000000000000;
    fwd_beta[212]  = -0.32000000000000001;
    fwd_Ea[212]    = -131;
    prefactor_units[212]  = 1.0000000000000002e-06;
    activation_units[212] = 0.50321666580471969;
    phase_units[212]      = 1e-12;
    is_PD[212] = 0;
    nTB[212] = 0;

    // (213):  CH3CCH2 + O2 <=> CH3CO + CH2O
    fwd_A[213]     = 100000000000;
    fwd_beta[213]  = 0;
    fwd_Ea[213]    = 0;
    prefactor_units[213]  = 1.0000000000000002e-06;
    activation_units[213] = 0.50321666580471969;
    phase_units[213]      = 1e-12;
    is_PD[213] = 0;
    nTB[213] = 0;

    // (214):  CH3CCH2 + HO2 <=> CH3 + CH2CO + OH
    fwd_A[214]     = 20000000000000;
    fwd_beta[214]  = 0;
    fwd_Ea[214]    = 0;
    prefactor_units[214]  = 1.0000000000000002e-06;
    activation_units[214] = 0.50321666580471969;
    phase_units[214]      = 1e-12;
    is_PD[214] = 0;
    nTB[214] = 0;

    // (215):  C3H6 + H <=> C2H4 + CH3
    fwd_A[215]     = 8e+21;
    fwd_beta[215]  = -2.3900000000000001;
    fwd_Ea[215]    = 11180;
    prefactor_units[215]  = 1.0000000000000002e-06;
    activation_units[215] = 0.50321666580471969;
    phase_units[215]      = 1e-12;
    is_PD[215] = 0;
    nTB[215] = 0;

    // (216):  C3H6 + H <=> aC3H5 + H2
    fwd_A[216]     = 173000;
    fwd_beta[216]  = 2.5;
    fwd_Ea[216]    = 2490;
    prefactor_units[216]  = 1.0000000000000002e-06;
    activation_units[216] = 0.50321666580471969;
    phase_units[216]      = 1e-12;
    is_PD[216] = 0;
    nTB[216] = 0;

    // (217):  C3H6 + H <=> CH3CCH2 + H2
    fwd_A[217]     = 400000;
    fwd_beta[217]  = 2.5;
    fwd_Ea[217]    = 9790;
    prefactor_units[217]  = 1.0000000000000002e-06;
    activation_units[217] = 0.50321666580471969;
    phase_units[217]      = 1e-12;
    is_PD[217] = 0;
    nTB[217] = 0;

    // (218):  C3H6 + O <=> CH2CO + CH3 + H
    fwd_A[218]     = 80000000;
    fwd_beta[218]  = 1.6499999999999999;
    fwd_Ea[218]    = 327;
    prefactor_units[218]  = 1.0000000000000002e-06;
    activation_units[218] = 0.50321666580471969;
    phase_units[218]      = 1e-12;
    is_PD[218] = 0;
    nTB[218] = 0;

    // (219):  C3H6 + O <=> C2H3CHO + 2 H
    fwd_A[219]     = 40000000;
    fwd_beta[219]  = 1.6499999999999999;
    fwd_Ea[219]    = 327;
    prefactor_units[219]  = 1.0000000000000002e-06;
    activation_units[219] = 0.50321666580471969;
    phase_units[219]      = 1e-12;
    is_PD[219] = 0;
    nTB[219] = 0;

    // (220):  C3H6 + O <=> C2H5 + HCO
    fwd_A[220]     = 35000000;
    fwd_beta[220]  = 1.6499999999999999;
    fwd_Ea[220]    = -972;
    prefactor_units[220]  = 1.0000000000000002e-06;
    activation_units[220] = 0.50321666580471969;
    phase_units[220]      = 1e-12;
    is_PD[220] = 0;
    nTB[220] = 0;

    // (221):  C3H6 + O <=> aC3H5 + OH
    fwd_A[221]     = 180000000000;
    fwd_beta[221]  = 0.69999999999999996;
    fwd_Ea[221]    = 5880;
    prefactor_units[221]  = 1.0000000000000002e-06;
    activation_units[221] = 0.50321666580471969;
    phase_units[221]      = 1e-12;
    is_PD[221] = 0;
    nTB[221] = 0;

    // (222):  C3H6 + O <=> CH3CCH2 + OH
    fwd_A[222]     = 60000000000;
    fwd_beta[222]  = 0.69999999999999996;
    fwd_Ea[222]    = 7630;
    prefactor_units[222]  = 1.0000000000000002e-06;
    activation_units[222] = 0.50321666580471969;
    phase_units[222]      = 1e-12;
    is_PD[222] = 0;
    nTB[222] = 0;

    // (223):  C3H6 + OH <=> aC3H5 + H2O
    fwd_A[223]     = 3100000;
    fwd_beta[223]  = 2;
    fwd_Ea[223]    = -298;
    prefactor_units[223]  = 1.0000000000000002e-06;
    activation_units[223] = 0.50321666580471969;
    phase_units[223]      = 1e-12;
    is_PD[223] = 0;
    nTB[223] = 0;

    // (224):  C3H6 + OH <=> CH3CCH2 + H2O
    fwd_A[224]     = 1100000;
    fwd_beta[224]  = 2;
    fwd_Ea[224]    = 1450;
    prefactor_units[224]  = 1.0000000000000002e-06;
    activation_units[224] = 0.50321666580471969;
    phase_units[224]      = 1e-12;
    is_PD[224] = 0;
    nTB[224] = 0;

    // (225):  C3H6 + CH3 <=> aC3H5 + CH4
    fwd_A[225]     = 2.2000000000000002;
    fwd_beta[225]  = 3.5;
    fwd_Ea[225]    = 5675;
    prefactor_units[225]  = 1.0000000000000002e-06;
    activation_units[225] = 0.50321666580471969;
    phase_units[225]      = 1e-12;
    is_PD[225] = 0;
    nTB[225] = 0;

    // (226):  C2H3CHO + O <=> C2H3 + OH + CO
    fwd_A[226]     = 30000000000000;
    fwd_beta[226]  = 0;
    fwd_Ea[226]    = 3540;
    prefactor_units[226]  = 1.0000000000000002e-06;
    activation_units[226] = 0.50321666580471969;
    phase_units[226]      = 1e-12;
    is_PD[226] = 0;
    nTB[226] = 0;

    // (227):  C2H3CHO + O <=> CH2O + CH2CO
    fwd_A[227]     = 19000000;
    fwd_beta[227]  = 1.8;
    fwd_Ea[227]    = 220;
    prefactor_units[227]  = 1.0000000000000002e-06;
    activation_units[227] = 0.50321666580471969;
    phase_units[227]      = 1e-12;
    is_PD[227] = 0;
    nTB[227] = 0;

    // (228):  iC3H7 + H <=> CH3 + C2H5
    fwd_A[228]     = 1.4000000000000001e+28;
    fwd_beta[228]  = -3.9399999999999999;
    fwd_Ea[228]    = 15916;
    prefactor_units[228]  = 1.0000000000000002e-06;
    activation_units[228] = 0.50321666580471969;
    phase_units[228]      = 1e-12;
    is_PD[228] = 0;
    nTB[228] = 0;

    // (229):  iC3H7 + O <=> CH3CHO + CH3
    fwd_A[229]     = 96000000000000;
    fwd_beta[229]  = 0;
    fwd_Ea[229]    = 0;
    prefactor_units[229]  = 1.0000000000000002e-06;
    activation_units[229] = 0.50321666580471969;
    phase_units[229]      = 1e-12;
    is_PD[229] = 0;
    nTB[229] = 0;

    // (230):  iC3H7 + OH <=> C3H6 + H2O
    fwd_A[230]     = 24000000000000;
    fwd_beta[230]  = 0;
    fwd_Ea[230]    = 0;
    prefactor_units[230]  = 1.0000000000000002e-06;
    activation_units[230] = 0.50321666580471969;
    phase_units[230]      = 1e-12;
    is_PD[230] = 0;
    nTB[230] = 0;

    // (231):  iC3H7 + O2 <=> C3H6 + HO2
    fwd_A[231]     = 130000000000;
    fwd_beta[231]  = 0;
    fwd_Ea[231]    = 0;
    prefactor_units[231]  = 1.0000000000000002e-06;
    activation_units[231] = 0.50321666580471969;
    phase_units[231]      = 1e-12;
    is_PD[231] = 0;
    nTB[231] = 0;

    // (232):  iC3H7 + HO2 <=> CH3CHO + CH3 + OH
    fwd_A[232]     = 24000000000000;
    fwd_beta[232]  = 0;
    fwd_Ea[232]    = 0;
    prefactor_units[232]  = 1.0000000000000002e-06;
    activation_units[232] = 0.50321666580471969;
    phase_units[232]      = 1e-12;
    is_PD[232] = 0;
    nTB[232] = 0;

    // (233):  iC3H7 + CH3 <=> CH4 + C3H6
    fwd_A[233]     = 220000000000000;
    fwd_beta[233]  = -0.68000000000000005;
    fwd_Ea[233]    = 0;
    prefactor_units[233]  = 1.0000000000000002e-06;
    activation_units[233] = 0.50321666580471969;
    phase_units[233]      = 1e-12;
    is_PD[233] = 0;
    nTB[233] = 0;

    // (234):  nC3H7 + H <=> C2H5 + CH3
    fwd_A[234]     = 3.7e+24;
    fwd_beta[234]  = -2.9199999999999999;
    fwd_Ea[234]    = 12505;
    prefactor_units[234]  = 1.0000000000000002e-06;
    activation_units[234] = 0.50321666580471969;
    phase_units[234]      = 1e-12;
    is_PD[234] = 0;
    nTB[234] = 0;

    // (235):  nC3H7 + OH <=> C3H6 + H2O
    fwd_A[235]     = 24000000000000;
    fwd_beta[235]  = 0;
    fwd_Ea[235]    = 0;
    prefactor_units[235]  = 1.0000000000000002e-06;
    activation_units[235] = 0.50321666580471969;
    phase_units[235]      = 1e-12;
    is_PD[235] = 0;
    nTB[235] = 0;

    // (236):  nC3H7 + O2 <=> C3H6 + HO2
    fwd_A[236]     = 90000000000;
    fwd_beta[236]  = 0;
    fwd_Ea[236]    = 0;
    prefactor_units[236]  = 1.0000000000000002e-06;
    activation_units[236] = 0.50321666580471969;
    phase_units[236]      = 1e-12;
    is_PD[236] = 0;
    nTB[236] = 0;

    // (237):  nC3H7 + HO2 <=> C2H5 + OH + CH2O
    fwd_A[237]     = 24000000000000;
    fwd_beta[237]  = 0;
    fwd_Ea[237]    = 0;
    prefactor_units[237]  = 1.0000000000000002e-06;
    activation_units[237] = 0.50321666580471969;
    phase_units[237]      = 1e-12;
    is_PD[237] = 0;
    nTB[237] = 0;

    // (238):  nC3H7 + CH3 <=> CH4 + C3H6
    fwd_A[238]     = 11000000000000;
    fwd_beta[238]  = 0;
    fwd_Ea[238]    = 0;
    prefactor_units[238]  = 1.0000000000000002e-06;
    activation_units[238] = 0.50321666580471969;
    phase_units[238]      = 1e-12;
    is_PD[238] = 0;
    nTB[238] = 0;

    // (239):  C4H2 + H <=> iC4H3
    fwd_A[239]     = 1.1e+30;
    fwd_beta[239]  = -4.9199999999999999;
    fwd_Ea[239]    = 10800;
    prefactor_units[239]  = 1.0000000000000002e-06;
    activation_units[239] = 0.50321666580471969;
    phase_units[239]      = 1e-12;
    is_PD[239] = 0;
    nTB[239] = 0;

    // (240):  iC4H3 + H <=> C4H2 + H2
    fwd_A[240]     = 60000000000000;
    fwd_beta[240]  = 0;
    fwd_Ea[240]    = 0;
    prefactor_units[240]  = 1.0000000000000002e-06;
    activation_units[240] = 0.50321666580471969;
    phase_units[240]      = 1e-12;
    is_PD[240] = 0;
    nTB[240] = 0;

    // (241):  C4H4 + OH <=> iC4H3 + H2O
    fwd_A[241]     = 15500000;
    fwd_beta[241]  = 2;
    fwd_Ea[241]    = 430;
    prefactor_units[241]  = 1.0000000000000002e-06;
    activation_units[241] = 0.50321666580471969;
    phase_units[241]      = 1e-12;
    is_PD[241] = 0;
    nTB[241] = 0;

    // (242):  C4H5-2 <=> iC4H5
    fwd_A[242]     = 1.5e+67;
    fwd_beta[242]  = -16.890000000000001;
    fwd_Ea[242]    = 59100;
    prefactor_units[242]  = 1;
    activation_units[242] = 0.50321666580471969;
    phase_units[242]      = 1e-6;
    is_PD[242] = 0;
    nTB[242] = 0;

    // (243):  C4H6 + H <=> C2H4 + C2H3
    fwd_A[243]     = 1.4599999999999999e+30;
    fwd_beta[243]  = -4.3399999999999999;
    fwd_Ea[243]    = 21647;
    prefactor_units[243]  = 1.0000000000000002e-06;
    activation_units[243] = 0.50321666580471969;
    phase_units[243]      = 1e-12;
    is_PD[243] = 0;
    nTB[243] = 0;

    // (244):  C4H6 + OH <=> iC4H5 + H2O
    fwd_A[244]     = 3100000;
    fwd_beta[244]  = 2;
    fwd_Ea[244]    = 430;
    prefactor_units[244]  = 1.0000000000000002e-06;
    activation_units[244] = 0.50321666580471969;
    phase_units[244]      = 1e-12;
    is_PD[244] = 0;
    nTB[244] = 0;

    // (245):  C4H612 <=> iC4H5 + H
    fwd_A[245]     = 4200000000000000;
    fwd_beta[245]  = 0;
    fwd_Ea[245]    = 92600;
    prefactor_units[245]  = 1;
    activation_units[245] = 0.50321666580471969;
    phase_units[245]      = 1e-6;
    is_PD[245] = 0;
    nTB[245] = 0;

    // (246):  C4H6-2 <=> H + C4H5-2
    fwd_A[246]     = 5000000000000000;
    fwd_beta[246]  = 0;
    fwd_Ea[246]    = 87300;
    prefactor_units[246]  = 1;
    activation_units[246] = 0.50321666580471969;
    phase_units[246]      = 1e-6;
    is_PD[246] = 0;
    nTB[246] = 0;

    // (247):  C4H7 <=> C4H6 + H
    fwd_A[247]     = 2.48e+53;
    fwd_beta[247]  = -12.300000000000001;
    fwd_Ea[247]    = 52000;
    prefactor_units[247]  = 1;
    activation_units[247] = 0.50321666580471969;
    phase_units[247]      = 1e-6;
    is_PD[247] = 0;
    nTB[247] = 0;

    // (248):  C4H7 + O2 <=> C4H6 + HO2
    fwd_A[248]     = 100000000000;
    fwd_beta[248]  = 0;
    fwd_Ea[248]    = 0;
    prefactor_units[248]  = 1.0000000000000002e-06;
    activation_units[248] = 0.50321666580471969;
    phase_units[248]      = 1e-12;
    is_PD[248] = 0;
    nTB[248] = 0;

    // (249):  C4H7 + HO2 <=> CH2O + OH + aC3H5
    fwd_A[249]     = 24000000000000;
    fwd_beta[249]  = 0;
    fwd_Ea[249]    = 0;
    prefactor_units[249]  = 1.0000000000000002e-06;
    activation_units[249] = 0.50321666580471969;
    phase_units[249]      = 1e-12;
    is_PD[249] = 0;
    nTB[249] = 0;

    // (250):  C4H81 + H <=> C2H4 + C2H5
    fwd_A[250]     = 1.6e+22;
    fwd_beta[250]  = -2.3900000000000001;
    fwd_Ea[250]    = 11180;
    prefactor_units[250]  = 1.0000000000000002e-06;
    activation_units[250] = 0.50321666580471969;
    phase_units[250]      = 1e-12;
    is_PD[250] = 0;
    nTB[250] = 0;

    // (251):  C4H81 + H <=> C3H6 + CH3
    fwd_A[251]     = 3.2e+22;
    fwd_beta[251]  = -2.3900000000000001;
    fwd_Ea[251]    = 11180;
    prefactor_units[251]  = 1.0000000000000002e-06;
    activation_units[251] = 0.50321666580471969;
    phase_units[251]      = 1e-12;
    is_PD[251] = 0;
    nTB[251] = 0;

    // (252):  C4H81 + H <=> C4H7 + H2
    fwd_A[252]     = 650000;
    fwd_beta[252]  = 2.54;
    fwd_Ea[252]    = 6756;
    prefactor_units[252]  = 1.0000000000000002e-06;
    activation_units[252] = 0.50321666580471969;
    phase_units[252]      = 1e-12;
    is_PD[252] = 0;
    nTB[252] = 0;

    // (253):  C4H81 + O <=> nC3H7 + HCO
    fwd_A[253]     = 330000000;
    fwd_beta[253]  = 1.45;
    fwd_Ea[253]    = -402;
    prefactor_units[253]  = 1.0000000000000002e-06;
    activation_units[253] = 0.50321666580471969;
    phase_units[253]      = 1e-12;
    is_PD[253] = 0;
    nTB[253] = 0;

    // (254):  C2H4 + C2H5 <=> pC4H9
    fwd_A[254]     = 150000000000;
    fwd_beta[254]  = 0;
    fwd_Ea[254]    = 7300;
    prefactor_units[254]  = 1.0000000000000002e-06;
    activation_units[254] = 0.50321666580471969;
    phase_units[254]      = 1e-12;
    is_PD[254] = 0;
    nTB[254] = 0;

    // (255):  pC4H9 + OH <=> C4H81 + H2O
    fwd_A[255]     = 24000000000000;
    fwd_beta[255]  = 0;
    fwd_Ea[255]    = 0;
    prefactor_units[255]  = 1.0000000000000002e-06;
    activation_units[255] = 0.50321666580471969;
    phase_units[255]      = 1e-12;
    is_PD[255] = 0;
    nTB[255] = 0;

    // (256):  pC4H9 + O2 <=> C4H81 + HO2
    fwd_A[256]     = 270000000000;
    fwd_beta[256]  = 0;
    fwd_Ea[256]    = 0;
    prefactor_units[256]  = 1.0000000000000002e-06;
    activation_units[256] = 0.50321666580471969;
    phase_units[256]      = 1e-12;
    is_PD[256] = 0;
    nTB[256] = 0;

    // (257):  pC4H9 + HO2 <=> nC3H7 + OH + CH2O
    fwd_A[257]     = 24000000000000;
    fwd_beta[257]  = 0;
    fwd_Ea[257]    = 0;
    prefactor_units[257]  = 1.0000000000000002e-06;
    activation_units[257] = 0.50321666580471969;
    phase_units[257]      = 1e-12;
    is_PD[257] = 0;
    nTB[257] = 0;

    // (258):  pC4H9 + CH3 <=> C4H81 + CH4
    fwd_A[258]     = 11000000000000;
    fwd_beta[258]  = 0;
    fwd_Ea[258]    = 0;
    prefactor_units[258]  = 1.0000000000000002e-06;
    activation_units[258] = 0.50321666580471969;
    phase_units[258]      = 1e-12;
    is_PD[258] = 0;
    nTB[258] = 0;

    // (259):  NC12H26 => 3 C2H4 + 2 nC3H7
    fwd_A[259]     = 5.6399999999999997e+26;
    fwd_beta[259]  = -2.6800000000000002;
    fwd_Ea[259]    = 88171;
    prefactor_units[259]  = 1;
    activation_units[259] = 0.50321666580471969;
    phase_units[259]      = 1e-6;
    is_PD[259] = 0;
    nTB[259] = 0;

    // (260):  NC12H26 => 2 C2H4 + 2 pC4H9
    fwd_A[260]     = 5.11e+25;
    fwd_beta[260]  = -2.5099999999999998;
    fwd_Ea[260]    = 88117;
    prefactor_units[260]  = 1;
    activation_units[260] = 0.50321666580471969;
    phase_units[260]      = 1e-6;
    is_PD[260] = 0;
    nTB[260] = 0;

    // (261):  NC12H26 + H => 4 C2H4 + pC4H9 + H2
    fwd_A[261]     = 1300000;
    fwd_beta[261]  = 2.54;
    fwd_Ea[261]    = 6756;
    prefactor_units[261]  = 1.0000000000000002e-06;
    activation_units[261] = 0.50321666580471969;
    phase_units[261]      = 1e-12;
    is_PD[261] = 0;
    nTB[261] = 0;

    // (262):  NC12H26 + H => C4H81 + 2 C2H4 + pC4H9 + H2
    fwd_A[262]     = 1300000;
    fwd_beta[262]  = 2.3999999999999999;
    fwd_Ea[262]    = 4471;
    prefactor_units[262]  = 1.0000000000000002e-06;
    activation_units[262] = 0.50321666580471969;
    phase_units[262]      = 1e-12;
    is_PD[262] = 0;
    nTB[262] = 0;

    // (263):  NC12H26 + H => C3H6 + C6H12 + nC3H7 + H2
    fwd_A[263]     = 1300000;
    fwd_beta[263]  = 2.3999999999999999;
    fwd_Ea[263]    = 4471;
    prefactor_units[263]  = 1.0000000000000002e-06;
    activation_units[263] = 0.50321666580471969;
    phase_units[263]      = 1e-12;
    is_PD[263] = 0;
    nTB[263] = 0;

    // (264):  NC12H26 + H => C5H10 + 2 C2H4 + nC3H7 + H2
    fwd_A[264]     = 1300000;
    fwd_beta[264]  = 2.3999999999999999;
    fwd_Ea[264]    = 4471;
    prefactor_units[264]  = 1.0000000000000002e-06;
    activation_units[264] = 0.50321666580471969;
    phase_units[264]      = 1e-12;
    is_PD[264] = 0;
    nTB[264] = 0;

    // (265):  NC12H26 + H => C6H12 + C2H4 + pC4H9 + H2
    fwd_A[265]     = 2600000;
    fwd_beta[265]  = 2.3999999999999999;
    fwd_Ea[265]    = 4471;
    prefactor_units[265]  = 1.0000000000000002e-06;
    activation_units[265] = 0.50321666580471969;
    phase_units[265]      = 1e-12;
    is_PD[265] = 0;
    nTB[265] = 0;

    // (266):  NC12H26 + CH3 => 4 C2H4 + pC4H9 + CH4
    fwd_A[266]     = 1.8100000000000001;
    fwd_beta[266]  = 3.6499999999999999;
    fwd_Ea[266]    = 7153;
    prefactor_units[266]  = 1.0000000000000002e-06;
    activation_units[266] = 0.50321666580471969;
    phase_units[266]      = 1e-12;
    is_PD[266] = 0;
    nTB[266] = 0;

    // (267):  NC12H26 + CH3 => C4H81 + 2 C2H4 + pC4H9 + CH4
    fwd_A[267]     = 3;
    fwd_beta[267]  = 3.46;
    fwd_Ea[267]    = 5480;
    prefactor_units[267]  = 1.0000000000000002e-06;
    activation_units[267] = 0.50321666580471969;
    phase_units[267]      = 1e-12;
    is_PD[267] = 0;
    nTB[267] = 0;

    // (268):  NC12H26 + CH3 => C3H6 + C6H12 + nC3H7 + CH4
    fwd_A[268]     = 3;
    fwd_beta[268]  = 3.46;
    fwd_Ea[268]    = 5480;
    prefactor_units[268]  = 1.0000000000000002e-06;
    activation_units[268] = 0.50321666580471969;
    phase_units[268]      = 1e-12;
    is_PD[268] = 0;
    nTB[268] = 0;

    // (269):  NC12H26 + CH3 => C5H10 + 2 C2H4 + nC3H7 + CH4
    fwd_A[269]     = 3;
    fwd_beta[269]  = 3.46;
    fwd_Ea[269]    = 5480;
    prefactor_units[269]  = 1.0000000000000002e-06;
    activation_units[269] = 0.50321666580471969;
    phase_units[269]      = 1e-12;
    is_PD[269] = 0;
    nTB[269] = 0;

    // (270):  NC12H26 + CH3 => C6H12 + C2H4 + pC4H9 + CH4
    fwd_A[270]     = 6;
    fwd_beta[270]  = 3.46;
    fwd_Ea[270]    = 5480;
    prefactor_units[270]  = 1.0000000000000002e-06;
    activation_units[270] = 0.50321666580471969;
    phase_units[270]      = 1e-12;
    is_PD[270] = 0;
    nTB[270] = 0;

    // (271):  NC12H26 + O => 4 C2H4 + pC4H9 + OH
    fwd_A[271]     = 190000;
    fwd_beta[271]  = 2.6800000000000002;
    fwd_Ea[271]    = 3716;
    prefactor_units[271]  = 1.0000000000000002e-06;
    activation_units[271] = 0.50321666580471969;
    phase_units[271]      = 1e-12;
    is_PD[271] = 0;
    nTB[271] = 0;

    // (272):  NC12H26 + O => C4H81 + 2 C2H4 + pC4H9 + OH
    fwd_A[272]     = 47600;
    fwd_beta[272]  = 2.71;
    fwd_Ea[272]    = 2106;
    prefactor_units[272]  = 1.0000000000000002e-06;
    activation_units[272] = 0.50321666580471969;
    phase_units[272]      = 1e-12;
    is_PD[272] = 0;
    nTB[272] = 0;

    // (273):  NC12H26 + O => C3H6 + C6H12 + nC3H7 + OH
    fwd_A[273]     = 47600;
    fwd_beta[273]  = 2.71;
    fwd_Ea[273]    = 2106;
    prefactor_units[273]  = 1.0000000000000002e-06;
    activation_units[273] = 0.50321666580471969;
    phase_units[273]      = 1e-12;
    is_PD[273] = 0;
    nTB[273] = 0;

    // (274):  NC12H26 + O => C5H10 + 2 C2H4 + nC3H7 + OH
    fwd_A[274]     = 47600;
    fwd_beta[274]  = 2.71;
    fwd_Ea[274]    = 2106;
    prefactor_units[274]  = 1.0000000000000002e-06;
    activation_units[274] = 0.50321666580471969;
    phase_units[274]      = 1e-12;
    is_PD[274] = 0;
    nTB[274] = 0;

    // (275):  NC12H26 + O => C6H12 + C2H4 + pC4H9 + OH
    fwd_A[275]     = 95200;
    fwd_beta[275]  = 2.71;
    fwd_Ea[275]    = 2106;
    prefactor_units[275]  = 1.0000000000000002e-06;
    activation_units[275] = 0.50321666580471969;
    phase_units[275]      = 1e-12;
    is_PD[275] = 0;
    nTB[275] = 0;

    // (276):  NC12H26 + OH => 4 C2H4 + pC4H9 + H2O
    fwd_A[276]     = 1400;
    fwd_beta[276]  = 2.6600000000000001;
    fwd_Ea[276]    = 527;
    prefactor_units[276]  = 1.0000000000000002e-06;
    activation_units[276] = 0.50321666580471969;
    phase_units[276]      = 1e-12;
    is_PD[276] = 0;
    nTB[276] = 0;

    // (277):  NC12H26 + OH => C4H81 + 2 C2H4 + pC4H9 + H2O
    fwd_A[277]     = 27000;
    fwd_beta[277]  = 2.3900000000000001;
    fwd_Ea[277]    = 393;
    prefactor_units[277]  = 1.0000000000000002e-06;
    activation_units[277] = 0.50321666580471969;
    phase_units[277]      = 1e-12;
    is_PD[277] = 0;
    nTB[277] = 0;

    // (278):  NC12H26 + OH => C3H6 + C6H12 + nC3H7 + H2O
    fwd_A[278]     = 27000;
    fwd_beta[278]  = 2.3900000000000001;
    fwd_Ea[278]    = 393;
    prefactor_units[278]  = 1.0000000000000002e-06;
    activation_units[278] = 0.50321666580471969;
    phase_units[278]      = 1e-12;
    is_PD[278] = 0;
    nTB[278] = 0;

    // (279):  NC12H26 + OH => C5H10 + 2 C2H4 + nC3H7 + H2O
    fwd_A[279]     = 27000;
    fwd_beta[279]  = 2.3900000000000001;
    fwd_Ea[279]    = 393;
    prefactor_units[279]  = 1.0000000000000002e-06;
    activation_units[279] = 0.50321666580471969;
    phase_units[279]      = 1e-12;
    is_PD[279] = 0;
    nTB[279] = 0;

    // (280):  NC12H26 + OH => C6H12 + C2H4 + pC4H9 + H2O
    fwd_A[280]     = 54000;
    fwd_beta[280]  = 2.3900000000000001;
    fwd_Ea[280]    = 393;
    prefactor_units[280]  = 1.0000000000000002e-06;
    activation_units[280] = 0.50321666580471969;
    phase_units[280]      = 1e-12;
    is_PD[280] = 0;
    nTB[280] = 0;

    // (281):  C6H12 + H <=> C2H4 + pC4H9
    fwd_A[281]     = 8e+21;
    fwd_beta[281]  = -2.3900000000000001;
    fwd_Ea[281]    = 11180;
    prefactor_units[281]  = 1.0000000000000002e-06;
    activation_units[281] = 0.50321666580471969;
    phase_units[281]      = 1e-12;
    is_PD[281] = 0;
    nTB[281] = 0;

    // (282):  C6H12 + H <=> C6H11 + H2
    fwd_A[282]     = 650000;
    fwd_beta[282]  = 2.54;
    fwd_Ea[282]    = 6756;
    prefactor_units[282]  = 1.0000000000000002e-06;
    activation_units[282] = 0.50321666580471969;
    phase_units[282]      = 1e-12;
    is_PD[282] = 0;
    nTB[282] = 0;

    // (283):  C5H10 + H <=> C2H4 + nC3H7
    fwd_A[283]     = 8e+21;
    fwd_beta[283]  = -2.3900000000000001;
    fwd_Ea[283]    = 11180;
    prefactor_units[283]  = 1.0000000000000002e-06;
    activation_units[283] = 0.50321666580471969;
    phase_units[283]      = 1e-12;
    is_PD[283] = 0;
    nTB[283] = 0;

    // (284):  C5H10 + H <=> C2H4 + aC3H5 + H2
    fwd_A[284]     = 650000;
    fwd_beta[284]  = 2.54;
    fwd_Ea[284]    = 6756;
    prefactor_units[284]  = 1.0000000000000002e-06;
    activation_units[284] = 0.50321666580471969;
    phase_units[284]      = 1e-12;
    is_PD[284] = 0;
    nTB[284] = 0;

    // (285):  C6H11 + H <=> CH3 + C2H4 + aC3H5
    fwd_A[285]     = 2e+21;
    fwd_beta[285]  = -2;
    fwd_Ea[285]    = 11000;
    prefactor_units[285]  = 1.0000000000000002e-06;
    activation_units[285] = 0.50321666580471969;
    phase_units[285]      = 1e-12;
    is_PD[285] = 0;
    nTB[285] = 0;

    // (286):  C6H11 + HO2 => CH2O + OH + aC3H5 + C2H4
    fwd_A[286]     = 24000000000000;
    fwd_beta[286]  = 0;
    fwd_Ea[286]    = 0;
    prefactor_units[286]  = 1.0000000000000002e-06;
    activation_units[286] = 0.50321666580471969;
    phase_units[286]      = 1e-12;
    is_PD[286] = 0;
    nTB[286] = 0;

    // (287):  C6H12 + O <=> C2H4 + nC3H7 + HCO
    fwd_A[287]     = 330000000;
    fwd_beta[287]  = 1.45;
    fwd_Ea[287]    = -402;
    prefactor_units[287]  = 1.0000000000000002e-06;
    activation_units[287] = 0.50321666580471969;
    phase_units[287]      = 1e-12;
    is_PD[287] = 0;
    nTB[287] = 0;

    // (288):  C5H10 + O <=> pC4H9 + HCO
    fwd_A[288]     = 330000000;
    fwd_beta[288]  = 1.45;
    fwd_Ea[288]    = -402;
    prefactor_units[288]  = 1.0000000000000002e-06;
    activation_units[288] = 0.50321666580471969;
    phase_units[288]      = 1e-12;
    is_PD[288] = 0;
    nTB[288] = 0;

    for (int i=0; i<289; ++i) {
        k_f[i] = prefactor_units[i] * fwd_A[i]
                    * exp(fwd_beta[i] * tc[0] - activation_units[i] * fwd_Ea[i] * invT);
    };
    /*compute the mixture concentration */
    double mixture = 0.0;
    for (int i = 0; i < 56; ++i) {
        mixture += sc[i];
    }

    /* troe */
    {
        double alpha[25];
        alpha[0] = mixture + (TB[0][0] - 1)*sc[5] + (TB[0][1] - 1)*sc[7] + (TB[0][2] - 1)*sc[18] + (TB[0][3] - 1)*sc[19];
        alpha[1] = mixture + (TB[1][0] - 1)*sc[4] + (TB[1][1] - 1)*sc[5] + (TB[1][2] - 1)*sc[18] + (TB[1][3] - 1)*sc[19];
        alpha[2] = mixture + (TB[2][0] - 1)*sc[4] + (TB[2][1] - 1)*sc[5] + (TB[2][2] - 1)*sc[12] + (TB[2][3] - 1)*sc[18] + (TB[2][4] - 1)*sc[19] + (TB[2][5] - 1)*sc[25];
        alpha[3] = mixture + (TB[3][0] - 1)*sc[4] + (TB[3][1] - 1)*sc[5] + (TB[3][2] - 1)*sc[12] + (TB[3][3] - 1)*sc[18] + (TB[3][4] - 1)*sc[19] + (TB[3][5] - 1)*sc[25];
        alpha[4] = mixture + (TB[4][0] - 1)*sc[4] + (TB[4][1] - 1)*sc[5] + (TB[4][2] - 1)*sc[12] + (TB[4][3] - 1)*sc[18] + (TB[4][4] - 1)*sc[19] + (TB[4][5] - 1)*sc[25];
        alpha[5] = mixture + (TB[5][0] - 1)*sc[4] + (TB[5][1] - 1)*sc[5] + (TB[5][2] - 1)*sc[12] + (TB[5][3] - 1)*sc[18] + (TB[5][4] - 1)*sc[19] + (TB[5][5] - 1)*sc[25];
        alpha[6] = mixture + (TB[6][0] - 1)*sc[4] + (TB[6][1] - 1)*sc[5] + (TB[6][2] - 1)*sc[12] + (TB[6][3] - 1)*sc[18] + (TB[6][4] - 1)*sc[19] + (TB[6][5] - 1)*sc[25];
        alpha[7] = mixture + (TB[7][0] - 1)*sc[4] + (TB[7][1] - 1)*sc[5] + (TB[7][2] - 1)*sc[12] + (TB[7][3] - 1)*sc[18] + (TB[7][4] - 1)*sc[19] + (TB[7][5] - 1)*sc[25];
        alpha[8] = mixture + (TB[8][0] - 1)*sc[4] + (TB[8][1] - 1)*sc[5] + (TB[8][2] - 1)*sc[12] + (TB[8][3] - 1)*sc[18] + (TB[8][4] - 1)*sc[19] + (TB[8][5] - 1)*sc[25];
        alpha[9] = mixture + (TB[9][0] - 1)*sc[4] + (TB[9][1] - 1)*sc[5] + (TB[9][2] - 1)*sc[12] + (TB[9][3] - 1)*sc[18] + (TB[9][4] - 1)*sc[19] + (TB[9][5] - 1)*sc[21] + (TB[9][6] - 1)*sc[23] + (TB[9][7] - 1)*sc[25];
        alpha[10] = mixture + (TB[10][0] - 1)*sc[4] + (TB[10][1] - 1)*sc[5] + (TB[10][2] - 1)*sc[12] + (TB[10][3] - 1)*sc[18] + (TB[10][4] - 1)*sc[19] + (TB[10][5] - 1)*sc[21] + (TB[10][6] - 1)*sc[23] + (TB[10][7] - 1)*sc[25];
        alpha[11] = mixture + (TB[11][0] - 1)*sc[4] + (TB[11][1] - 1)*sc[5] + (TB[11][2] - 1)*sc[12] + (TB[11][3] - 1)*sc[18] + (TB[11][4] - 1)*sc[19] + (TB[11][5] - 1)*sc[21] + (TB[11][6] - 1)*sc[23] + (TB[11][7] - 1)*sc[25];
        alpha[12] = mixture + (TB[12][0] - 1)*sc[4] + (TB[12][1] - 1)*sc[5] + (TB[12][2] - 1)*sc[12] + (TB[12][3] - 1)*sc[18] + (TB[12][4] - 1)*sc[19] + (TB[12][5] - 1)*sc[23] + (TB[12][6] - 1)*sc[25];
        alpha[13] = mixture + (TB[13][0] - 1)*sc[4] + (TB[13][1] - 1)*sc[5] + (TB[13][2] - 1)*sc[12] + (TB[13][3] - 1)*sc[18] + (TB[13][4] - 1)*sc[19] + (TB[13][5] - 1)*sc[21] + (TB[13][6] - 1)*sc[23] + (TB[13][7] - 1)*sc[25];
        alpha[14] = mixture + (TB[14][0] - 1)*sc[4] + (TB[14][1] - 1)*sc[5] + (TB[14][2] - 1)*sc[12] + (TB[14][3] - 1)*sc[18] + (TB[14][4] - 1)*sc[19] + (TB[14][5] - 1)*sc[21] + (TB[14][6] - 1)*sc[23] + (TB[14][7] - 1)*sc[25];
        alpha[15] = mixture + (TB[15][0] - 1)*sc[4] + (TB[15][1] - 1)*sc[5] + (TB[15][2] - 1)*sc[12] + (TB[15][3] - 1)*sc[18] + (TB[15][4] - 1)*sc[19] + (TB[15][5] - 1)*sc[25];
        alpha[16] = mixture;
        alpha[17] = mixture + (TB[17][0] - 1)*sc[4] + (TB[17][1] - 1)*sc[5] + (TB[17][2] - 1)*sc[12] + (TB[17][3] - 1)*sc[18] + (TB[17][4] - 1)*sc[19] + (TB[17][5] - 1)*sc[25];
        alpha[18] = mixture + (TB[18][0] - 1)*sc[4] + (TB[18][1] - 1)*sc[5] + (TB[18][2] - 1)*sc[12] + (TB[18][3] - 1)*sc[18] + (TB[18][4] - 1)*sc[19] + (TB[18][5] - 1)*sc[25];
        alpha[19] = mixture + (TB[19][0] - 1)*sc[4] + (TB[19][1] - 1)*sc[5] + (TB[19][2] - 1)*sc[12] + (TB[19][3] - 1)*sc[18] + (TB[19][4] - 1)*sc[19] + (TB[19][5] - 1)*sc[25];
        alpha[20] = mixture + (TB[20][0] - 1)*sc[4] + (TB[20][1] - 1)*sc[5] + (TB[20][2] - 1)*sc[12] + (TB[20][3] - 1)*sc[18] + (TB[20][4] - 1)*sc[19] + (TB[20][5] - 1)*sc[25];
        alpha[21] = mixture + (TB[21][0] - 1)*sc[4] + (TB[21][1] - 1)*sc[5] + (TB[21][2] - 1)*sc[12] + (TB[21][3] - 1)*sc[18] + (TB[21][4] - 1)*sc[19] + (TB[21][5] - 1)*sc[25];
        alpha[22] = mixture + (TB[22][0] - 1)*sc[4] + (TB[22][1] - 1)*sc[5] + (TB[22][2] - 1)*sc[12] + (TB[22][3] - 1)*sc[18] + (TB[22][4] - 1)*sc[19] + (TB[22][5] - 1)*sc[25];
        alpha[23] = mixture + (TB[23][0] - 1)*sc[4] + (TB[23][1] - 1)*sc[5] + (TB[23][2] - 1)*sc[12] + (TB[23][3] - 1)*sc[18] + (TB[23][4] - 1)*sc[19] + (TB[23][5] - 1)*sc[25];
        alpha[24] = mixture + (TB[24][0] - 1)*sc[4] + (TB[24][1] - 1)*sc[5] + (TB[24][2] - 1)*sc[12] + (TB[24][3] - 1)*sc[18] + (TB[24][4] - 1)*sc[19] + (TB[24][5] - 1)*sc[25];
        for (int i=0; i<25; i++)
        {
            double redP, F, logPred, logFcent, troe_c, troe_n, troe, F_troe;
            redP = alpha[i-0] / k_f[i] * phase_units[i] * low_A[i] * exp(low_beta[i] * tc[0] - activation_units[i] * low_Ea[i] *invT);
            F = redP / (1.0 + redP);
            logPred = log10(redP);
            logFcent = log10(
                (fabs(troe_Tsss[i]) > 1.e-100 ? (1.-troe_a[i])*exp(-tc[1]/troe_Tsss[i]) : 0.) 
                + (fabs(troe_Ts[i]) > 1.e-100 ? troe_a[i] * exp(-tc[1]/troe_Ts[i]) : 0.) 
                + (troe_len[i] == 4 ? exp(-troe_Tss[i] * invT) : 0.) );
            troe_c = -.4 - .67 * logFcent;
            troe_n = .75 - 1.27 * logFcent;
            troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
            F_troe = pow(10., logFcent / (1.0 + troe*troe));
            Corr[i] = F * F_troe;
        }
    }

    /* Lindemann */
    {
        double alpha;
        alpha = mixture + (TB[25][0] - 1)*sc[4] + (TB[25][1] - 1)*sc[5] + (TB[25][2] - 1)*sc[18] + (TB[25][3] - 1)*sc[19];
        double redP = alpha / k_f[25] * phase_units[25] * low_A[25] * exp(low_beta[25] * tc[0] - activation_units[25] * low_Ea[25] * invT);
        Corr[25] = redP / (1. + redP);
    }

    /* simple three-body correction */
    {
        double alpha;
        alpha = mixture + (TB[26][0] - 1)*sc[4] + (TB[26][1] - 1)*sc[5] + (TB[26][2] - 1)*sc[19];
        Corr[26] = alpha;
        alpha = mixture + (TB[27][0] - 1)*sc[4] + (TB[27][1] - 1)*sc[5] + (TB[27][2] - 1)*sc[18] + (TB[27][3] - 1)*sc[19];
        Corr[27] = alpha;
        alpha = mixture + (TB[28][0] - 1)*sc[4] + (TB[28][1] - 1)*sc[5] + (TB[28][2] - 1)*sc[18] + (TB[28][3] - 1)*sc[19];
        Corr[28] = alpha;
        alpha = mixture + (TB[29][0] - 1)*sc[4] + (TB[29][1] - 1)*sc[5] + (TB[29][2] - 1)*sc[18] + (TB[29][3] - 1)*sc[19];
        Corr[29] = alpha;
        alpha = mixture + (TB[30][0] - 1)*sc[4] + (TB[30][1] - 1)*sc[5] + (TB[30][2] - 1)*sc[18] + (TB[30][3] - 1)*sc[19];
        Corr[30] = alpha;
    }

    for (int i=0; i<289; i++) {
        if (nTB[i] != 0) {
            nTB[i] = 0;
            free(TB[i]);
            free(TBid[i]);
        }
    }
    return;
}

__device__ void comp_Kc_d(double * tc, double invT, double * Kc)
{
    /*compute the Gibbs free energy */
    double g_RT[56];
    gibbs_d(g_RT, tc);

    Kc[0] = g_RT[0] - g_RT[3] + g_RT[7];
    Kc[1] = 2*g_RT[2] - g_RT[6];
    Kc[2] = g_RT[9] + g_RT[18] - g_RT[27];
    Kc[3] = g_RT[5] + g_RT[10] - g_RT[17];
    Kc[4] = g_RT[0] + g_RT[14] - g_RT[16];
    Kc[5] = g_RT[0] + g_RT[14] - g_RT[15];
    Kc[6] = g_RT[0] + g_RT[11] - g_RT[12];
    Kc[7] = g_RT[2] + g_RT[11] - g_RT[17];
    Kc[8] = 2*g_RT[11] - g_RT[25];
    Kc[9] = -g_RT[0] - g_RT[21] + g_RT[22];
    Kc[10] = g_RT[0] + g_RT[27] - g_RT[29];
    Kc[11] = g_RT[0] + g_RT[22] - g_RT[23];
    Kc[12] = g_RT[11] + g_RT[22] - g_RT[36];
    Kc[13] = g_RT[11] + g_RT[18] - g_RT[28];
    Kc[14] = g_RT[11] + g_RT[13] - g_RT[30];
    Kc[15] = -g_RT[4] - g_RT[21] + g_RT[23];
    Kc[16] = g_RT[0] + g_RT[23] - g_RT[24];
    Kc[17] = g_RT[0] + g_RT[24] - g_RT[25];
    Kc[18] = g_RT[22] + g_RT[24] - g_RT[49];
    Kc[19] = g_RT[0] + g_RT[34] - g_RT[36];
    Kc[20] = g_RT[11] + g_RT[34] - g_RT[49];
    Kc[21] = g_RT[0] + g_RT[36] - g_RT[37];
    Kc[22] = g_RT[0] + g_RT[36] - g_RT[38];
    Kc[23] = g_RT[0] - g_RT[36] - g_RT[37] + g_RT[52];
    Kc[24] = g_RT[0] - g_RT[24] - g_RT[36] + g_RT[54];
    Kc[25] = g_RT[1] + g_RT[18] - g_RT[19];
    Kc[26] = 2*g_RT[0] - g_RT[4];
    Kc[27] = g_RT[0] + g_RT[2] - g_RT[5];
    Kc[28] = g_RT[0] + g_RT[1] - g_RT[2];
    Kc[29] = 2*g_RT[1] - g_RT[7];
    Kc[30] = -g_RT[0] + g_RT[13] - g_RT[18];
    Kc[31] = g_RT[0] - g_RT[1] - g_RT[2] + g_RT[7];
    Kc[32] = -g_RT[0] + g_RT[1] - g_RT[2] + g_RT[4];
    Kc[33] = -g_RT[0] + g_RT[2] + g_RT[4] - g_RT[5];
    Kc[34] = -g_RT[1] + 2*g_RT[2] - g_RT[5];
    Kc[35] = 2*g_RT[0] - g_RT[4] + g_RT[5] - g_RT[5];
    Kc[36] = -g_RT[0] - g_RT[3] + g_RT[4] + g_RT[7];
    Kc[37] = g_RT[0] - g_RT[1] + g_RT[3] - g_RT[5];
    Kc[38] = g_RT[0] - 2*g_RT[2] + g_RT[3];
    Kc[39] = g_RT[1] - g_RT[2] + g_RT[3] - g_RT[7];
    Kc[40] = g_RT[2] + g_RT[3] - g_RT[5] - g_RT[7];
    Kc[41] = g_RT[2] + g_RT[3] - g_RT[5] - g_RT[7];
    Kc[42] = 2*g_RT[3] - g_RT[6] - g_RT[7];
    Kc[43] = 2*g_RT[3] - g_RT[6] - g_RT[7];
    Kc[44] = g_RT[0] - g_RT[3] - g_RT[4] + g_RT[6];
    Kc[45] = g_RT[0] - g_RT[2] - g_RT[5] + g_RT[6];
    Kc[46] = g_RT[1] - g_RT[2] - g_RT[3] + g_RT[6];
    Kc[47] = g_RT[2] - g_RT[3] - g_RT[5] + g_RT[6];
    Kc[48] = g_RT[2] - g_RT[3] - g_RT[5] + g_RT[6];
    Kc[49] = -g_RT[0] + g_RT[2] + g_RT[18] - g_RT[19];
    Kc[50] = -g_RT[0] + g_RT[2] + g_RT[18] - g_RT[19];
    Kc[51] = -g_RT[2] + g_RT[3] + g_RT[18] - g_RT[19];
    Kc[52] = g_RT[0] - g_RT[4] + g_RT[13] - g_RT[18];
    Kc[53] = g_RT[1] - g_RT[2] + g_RT[13] - g_RT[18];
    Kc[54] = -g_RT[0] + g_RT[1] + g_RT[13] - g_RT[19];
    Kc[55] = g_RT[2] - g_RT[5] + g_RT[13] - g_RT[18];
    Kc[56] = -g_RT[0] + g_RT[5] - g_RT[5] + g_RT[13] - g_RT[18];
    Kc[57] = -g_RT[3] + g_RT[7] + g_RT[13] - g_RT[18];
    Kc[58] = -g_RT[0] + g_RT[1] + g_RT[8] - g_RT[18];
    Kc[59] = -g_RT[0] + g_RT[2] + g_RT[8] - g_RT[13];
    Kc[60] = -g_RT[0] + g_RT[4] + g_RT[8] - g_RT[9];
    Kc[61] = -g_RT[0] + g_RT[5] + g_RT[8] - g_RT[14];
    Kc[62] = -g_RT[1] + g_RT[7] + g_RT[8] - g_RT[13];
    Kc[63] = g_RT[8] - g_RT[13] - g_RT[18] + g_RT[19];
    Kc[64] = -g_RT[0] + g_RT[1] + g_RT[9] - g_RT[13];
    Kc[65] = -g_RT[0] + g_RT[2] + g_RT[9] - g_RT[14];
    Kc[66] = g_RT[2] - g_RT[5] - g_RT[8] + g_RT[9];
    Kc[67] = -g_RT[0] + g_RT[4] + g_RT[9] - g_RT[11];
    Kc[68] = -g_RT[2] + g_RT[7] + g_RT[9] - g_RT[13];
    Kc[69] = -2*g_RT[0] + g_RT[7] + g_RT[9] - g_RT[19];
    Kc[70] = -g_RT[2] + g_RT[3] + g_RT[9] - g_RT[14];
    Kc[71] = -g_RT[9] + g_RT[10] + g_RT[55] - g_RT[55];
    Kc[72] = g_RT[0] - g_RT[4] - g_RT[8] + g_RT[10];
    Kc[73] = -g_RT[0] + g_RT[2] + g_RT[10] - g_RT[14];
    Kc[74] = -g_RT[0] + g_RT[4] + g_RT[10] - g_RT[11];
    Kc[75] = -g_RT[0] - g_RT[2] + g_RT[7] + g_RT[10] - g_RT[18];
    Kc[76] = -g_RT[5] + g_RT[7] + g_RT[10] - g_RT[18];
    Kc[77] = g_RT[5] - g_RT[5] - g_RT[9] + g_RT[10];
    Kc[78] = -g_RT[9] + g_RT[10] + g_RT[18] - g_RT[18];
    Kc[79] = -g_RT[9] + g_RT[10] + g_RT[19] - g_RT[19];
    Kc[80] = g_RT[10] - g_RT[14] - g_RT[18] + g_RT[19];
    Kc[81] = g_RT[0] - g_RT[4] - g_RT[13] + g_RT[14];
    Kc[82] = g_RT[1] - g_RT[2] - g_RT[13] + g_RT[14];
    Kc[83] = g_RT[2] - g_RT[5] - g_RT[13] + g_RT[14];
    Kc[84] = -g_RT[3] + g_RT[7] - g_RT[13] + g_RT[14];
    Kc[85] = g_RT[3] - g_RT[6] - g_RT[13] + g_RT[14];
    Kc[86] = -g_RT[0] + g_RT[8] + g_RT[14] - g_RT[27];
    Kc[87] = -g_RT[0] + g_RT[1] + g_RT[11] - g_RT[14];
    Kc[88] = g_RT[2] - g_RT[5] - g_RT[9] + g_RT[11];
    Kc[89] = g_RT[2] - g_RT[5] - g_RT[10] + g_RT[11];
    Kc[90] = -g_RT[1] + g_RT[7] + g_RT[11] - g_RT[15];
    Kc[91] = -g_RT[2] + g_RT[7] + g_RT[11] - g_RT[14];
    Kc[92] = g_RT[3] - g_RT[7] + g_RT[11] - g_RT[12];
    Kc[93] = -g_RT[2] + g_RT[3] + g_RT[11] - g_RT[15];
    Kc[94] = -g_RT[0] + g_RT[8] + g_RT[11] - g_RT[22];
    Kc[95] = g_RT[11] - g_RT[12] + g_RT[13] - g_RT[18];
    Kc[96] = g_RT[11] - g_RT[12] - g_RT[13] + g_RT[14];
    Kc[97] = -g_RT[0] + g_RT[9] + g_RT[11] - g_RT[23];
    Kc[98] = -g_RT[0] + 2*g_RT[11] - g_RT[24];
    Kc[99] = g_RT[11] - g_RT[18] - g_RT[23] + g_RT[26];
    Kc[100] = g_RT[0] - g_RT[4] - g_RT[14] + g_RT[15];
    Kc[101] = g_RT[0] - g_RT[2] - g_RT[11] + g_RT[15];
    Kc[102] = g_RT[0] - g_RT[5] - g_RT[10] + g_RT[15];
    Kc[103] = g_RT[2] - g_RT[5] - g_RT[14] + g_RT[15];
    Kc[104] = -g_RT[3] + g_RT[7] - g_RT[14] + g_RT[15];
    Kc[105] = g_RT[0] - g_RT[4] - g_RT[14] + g_RT[16];
    Kc[106] = g_RT[0] - g_RT[2] - g_RT[11] + g_RT[16];
    Kc[107] = g_RT[0] - g_RT[5] - g_RT[10] + g_RT[16];
    Kc[108] = -g_RT[3] + g_RT[7] - g_RT[14] + g_RT[16];
    Kc[109] = g_RT[0] - g_RT[4] - g_RT[11] + g_RT[12];
    Kc[110] = g_RT[1] - g_RT[2] - g_RT[11] + g_RT[12];
    Kc[111] = g_RT[2] - g_RT[5] - g_RT[11] + g_RT[12];
    Kc[112] = -g_RT[0] + g_RT[8] + g_RT[12] - g_RT[23];
    Kc[113] = g_RT[9] - 2*g_RT[11] + g_RT[12];
    Kc[114] = g_RT[10] - 2*g_RT[11] + g_RT[12];
    Kc[115] = g_RT[0] - g_RT[4] - g_RT[16] + g_RT[17];
    Kc[116] = g_RT[0] - g_RT[4] - g_RT[15] + g_RT[17];
    Kc[117] = g_RT[1] - g_RT[2] - g_RT[16] + g_RT[17];
    Kc[118] = g_RT[2] - g_RT[5] - g_RT[16] + g_RT[17];
    Kc[119] = g_RT[2] - g_RT[5] - g_RT[15] + g_RT[17];
    Kc[120] = g_RT[1] - g_RT[8] - g_RT[18] + g_RT[20];
    Kc[121] = -g_RT[0] + g_RT[2] + g_RT[20] - g_RT[26];
    Kc[122] = g_RT[7] - g_RT[13] - g_RT[18] + g_RT[20];
    Kc[123] = -g_RT[0] + g_RT[4] + g_RT[20] - g_RT[21];
    Kc[124] = g_RT[0] - g_RT[10] - g_RT[18] + g_RT[26];
    Kc[125] = -g_RT[0] + g_RT[1] - 2*g_RT[18] + g_RT[26];
    Kc[126] = -g_RT[2] + g_RT[7] - 2*g_RT[18] + g_RT[26];
    Kc[127] = g_RT[1] - g_RT[2] - g_RT[20] + g_RT[21];
    Kc[128] = g_RT[1] - g_RT[9] - g_RT[18] + g_RT[21];
    Kc[129] = -g_RT[0] + g_RT[1] + g_RT[21] - g_RT[26];
    Kc[130] = -g_RT[0] + g_RT[2] + g_RT[21] - g_RT[27];
    Kc[131] = -g_RT[0] + g_RT[2] + g_RT[21] - g_RT[27];
    Kc[132] = g_RT[2] - g_RT[5] - g_RT[20] + g_RT[21];
    Kc[133] = g_RT[13] - g_RT[18] + g_RT[21] - g_RT[22];
    Kc[134] = -g_RT[0] + g_RT[9] + g_RT[21] - g_RT[31];
    Kc[135] = -g_RT[0] + g_RT[10] + g_RT[21] - g_RT[31];
    Kc[136] = -g_RT[0] + g_RT[20] + g_RT[21] - g_RT[40];
    Kc[137] = -g_RT[0] + g_RT[11] + g_RT[21] - g_RT[32];
    Kc[138] = -g_RT[0] + g_RT[11] + g_RT[21] - g_RT[33];
    Kc[139] = g_RT[0] - g_RT[4] - g_RT[26] + g_RT[27];
    Kc[140] = g_RT[0] - g_RT[11] - g_RT[18] + g_RT[27];
    Kc[141] = g_RT[2] - g_RT[5] - g_RT[26] + g_RT[27];
    Kc[142] = g_RT[0] - g_RT[4] - g_RT[21] + g_RT[22];
    Kc[143] = -g_RT[0] + g_RT[1] + g_RT[22] - g_RT[27];
    Kc[144] = g_RT[1] - g_RT[11] - g_RT[18] + g_RT[22];
    Kc[145] = g_RT[2] - g_RT[5] - g_RT[21] + g_RT[22];
    Kc[146] = -g_RT[3] + g_RT[7] - g_RT[21] + g_RT[22];
    Kc[147] = -g_RT[1] + g_RT[7] + g_RT[22] - g_RT[29];
    Kc[148] = g_RT[7] - g_RT[13] - g_RT[14] + g_RT[22];
    Kc[149] = -g_RT[2] + g_RT[3] + g_RT[22] - g_RT[29];
    Kc[150] = g_RT[13] - g_RT[18] + g_RT[22] - g_RT[23];
    Kc[151] = g_RT[13] + g_RT[22] - g_RT[39];
    Kc[152] = -g_RT[0] + g_RT[11] + g_RT[22] - g_RT[34];
    Kc[153] = -g_RT[11] - g_RT[18] + g_RT[29];
    Kc[154] = g_RT[0] - g_RT[0] - g_RT[28] + g_RT[29];
    Kc[155] = g_RT[0] - g_RT[11] - g_RT[13] + g_RT[29];
    Kc[156] = g_RT[0] - g_RT[4] - g_RT[27] + g_RT[29];
    Kc[157] = g_RT[2] - g_RT[5] - g_RT[27] + g_RT[29];
    Kc[158] = -g_RT[3] + g_RT[7] - g_RT[27] + g_RT[29];
    Kc[159] = -g_RT[2] + g_RT[7] - g_RT[14] - g_RT[18] + g_RT[29];
    Kc[160] = g_RT[0] - g_RT[11] - g_RT[13] + g_RT[28];
    Kc[161] = -g_RT[2] + g_RT[3] - g_RT[11] - g_RT[19] + g_RT[28];
    Kc[162] = g_RT[0] - g_RT[4] - g_RT[28] + g_RT[30];
    Kc[163] = g_RT[2] - g_RT[5] - g_RT[28] + g_RT[30];
    Kc[164] = g_RT[11] - g_RT[12] - g_RT[28] + g_RT[30];
    Kc[165] = -g_RT[3] + g_RT[7] - g_RT[28] + g_RT[30];
    Kc[166] = g_RT[0] - g_RT[4] - g_RT[22] + g_RT[23];
    Kc[167] = g_RT[1] - g_RT[2] - g_RT[22] + g_RT[23];
    Kc[168] = g_RT[1] - g_RT[11] - g_RT[13] + g_RT[23];
    Kc[169] = g_RT[1] - g_RT[9] - g_RT[14] + g_RT[23];
    Kc[170] = g_RT[2] - g_RT[5] - g_RT[22] + g_RT[23];
    Kc[171] = g_RT[13] - g_RT[18] + g_RT[23] - g_RT[24];
    Kc[172] = -g_RT[0] + g_RT[8] + g_RT[23] - g_RT[33];
    Kc[173] = -g_RT[0] + g_RT[8] + g_RT[23] - g_RT[32];
    Kc[174] = -g_RT[0] + g_RT[9] + g_RT[23] - g_RT[34];
    Kc[175] = -g_RT[0] + g_RT[10] + g_RT[23] - g_RT[34];
    Kc[176] = g_RT[11] - g_RT[12] - g_RT[22] + g_RT[23];
    Kc[177] = g_RT[11] + g_RT[23] - g_RT[37];
    Kc[178] = g_RT[22] + g_RT[23] - g_RT[48];
    Kc[179] = g_RT[0] - g_RT[4] - g_RT[23] + g_RT[24];
    Kc[180] = g_RT[1] - g_RT[11] - g_RT[14] + g_RT[24];
    Kc[181] = -g_RT[0] + g_RT[1] + g_RT[24] - g_RT[30];
    Kc[182] = -g_RT[3] + g_RT[7] - g_RT[23] + g_RT[24];
    Kc[183] = g_RT[3] - g_RT[7] + g_RT[24] - g_RT[25];
    Kc[184] = g_RT[3] - g_RT[6] - g_RT[23] + g_RT[24];
    Kc[185] = -g_RT[2] + g_RT[3] - g_RT[11] - g_RT[14] + g_RT[24];
    Kc[186] = -g_RT[11] + g_RT[22] + g_RT[24] - g_RT[34];
    Kc[187] = g_RT[0] - g_RT[4] - g_RT[24] + g_RT[25];
    Kc[188] = g_RT[1] - g_RT[2] - g_RT[24] + g_RT[25];
    Kc[189] = g_RT[2] - g_RT[5] - g_RT[24] + g_RT[25];
    Kc[190] = g_RT[10] - g_RT[11] - g_RT[24] + g_RT[25];
    Kc[191] = g_RT[11] - g_RT[12] - g_RT[24] + g_RT[25];
    Kc[192] = g_RT[0] + g_RT[31] - g_RT[32];
    Kc[193] = g_RT[1] - g_RT[14] - g_RT[20] + g_RT[31];
    Kc[194] = g_RT[7] - g_RT[13] - g_RT[27] + g_RT[31];
    Kc[195] = g_RT[3] - g_RT[7] + g_RT[31] - g_RT[32];
    Kc[196] = g_RT[0] + g_RT[33] - g_RT[35];
    Kc[197] = g_RT[0] + g_RT[33] - g_RT[34];
    Kc[198] = g_RT[1] - g_RT[18] - g_RT[23] + g_RT[33];
    Kc[199] = g_RT[2] - g_RT[5] - g_RT[31] + g_RT[33];
    Kc[200] = g_RT[32] - g_RT[33];
    Kc[201] = g_RT[0] - g_RT[0] + g_RT[32] - g_RT[33];
    Kc[202] = g_RT[0] + g_RT[32] - g_RT[35];
    Kc[203] = g_RT[1] - g_RT[18] - g_RT[23] + g_RT[32];
    Kc[204] = g_RT[2] - g_RT[5] - g_RT[31] + g_RT[32];
    Kc[205] = g_RT[0] - g_RT[4] - g_RT[33] + g_RT[34];
    Kc[206] = -g_RT[0] + g_RT[1] + g_RT[34] - g_RT[39];
    Kc[207] = -2*g_RT[0] + g_RT[2] + g_RT[34] - g_RT[39];
    Kc[208] = g_RT[2] - g_RT[5] - g_RT[33] + g_RT[34];
    Kc[209] = g_RT[3] - g_RT[7] + g_RT[34] - g_RT[36];
    Kc[210] = -g_RT[2] + g_RT[3] - g_RT[14] - g_RT[22] + g_RT[34];
    Kc[211] = g_RT[13] - g_RT[18] + g_RT[34] - g_RT[36];
    Kc[212] = g_RT[11] - g_RT[12] - g_RT[33] + g_RT[34];
    Kc[213] = g_RT[7] - g_RT[14] - g_RT[28] + g_RT[35];
    Kc[214] = -g_RT[2] + g_RT[3] - g_RT[11] - g_RT[27] + g_RT[35];
    Kc[215] = g_RT[0] - g_RT[11] - g_RT[23] + g_RT[36];
    Kc[216] = g_RT[0] - g_RT[4] - g_RT[34] + g_RT[36];
    Kc[217] = g_RT[0] - g_RT[4] - g_RT[35] + g_RT[36];
    Kc[218] = -g_RT[0] + g_RT[1] - g_RT[11] - g_RT[27] + g_RT[36];
    Kc[219] = -2*g_RT[0] + g_RT[1] + g_RT[36] - g_RT[39];
    Kc[220] = g_RT[1] - g_RT[13] - g_RT[24] + g_RT[36];
    Kc[221] = g_RT[1] - g_RT[2] - g_RT[34] + g_RT[36];
    Kc[222] = g_RT[1] - g_RT[2] - g_RT[35] + g_RT[36];
    Kc[223] = g_RT[2] - g_RT[5] - g_RT[34] + g_RT[36];
    Kc[224] = g_RT[2] - g_RT[5] - g_RT[35] + g_RT[36];
    Kc[225] = g_RT[11] - g_RT[12] - g_RT[34] + g_RT[36];
    Kc[226] = g_RT[1] - g_RT[2] - g_RT[18] - g_RT[22] + g_RT[39];
    Kc[227] = g_RT[1] - g_RT[14] - g_RT[27] + g_RT[39];
    Kc[228] = g_RT[0] - g_RT[11] - g_RT[24] + g_RT[38];
    Kc[229] = g_RT[1] - g_RT[11] - g_RT[30] + g_RT[38];
    Kc[230] = g_RT[2] - g_RT[5] - g_RT[36] + g_RT[38];
    Kc[231] = -g_RT[3] + g_RT[7] - g_RT[36] + g_RT[38];
    Kc[232] = -g_RT[2] + g_RT[3] - g_RT[11] - g_RT[30] + g_RT[38];
    Kc[233] = g_RT[11] - g_RT[12] - g_RT[36] + g_RT[38];
    Kc[234] = g_RT[0] - g_RT[11] - g_RT[24] + g_RT[37];
    Kc[235] = g_RT[2] - g_RT[5] - g_RT[36] + g_RT[37];
    Kc[236] = -g_RT[3] + g_RT[7] - g_RT[36] + g_RT[37];
    Kc[237] = -g_RT[2] + g_RT[3] - g_RT[14] - g_RT[24] + g_RT[37];
    Kc[238] = g_RT[11] - g_RT[12] - g_RT[36] + g_RT[37];
    Kc[239] = g_RT[0] + g_RT[40] - g_RT[41];
    Kc[240] = g_RT[0] - g_RT[4] - g_RT[40] + g_RT[41];
    Kc[241] = g_RT[2] - g_RT[5] - g_RT[41] + g_RT[42];
    Kc[242] = -g_RT[43] + g_RT[44];
    Kc[243] = g_RT[0] - g_RT[22] - g_RT[23] + g_RT[45];
    Kc[244] = g_RT[2] - g_RT[5] - g_RT[43] + g_RT[45];
    Kc[245] = -g_RT[0] - g_RT[43] + g_RT[46];
    Kc[246] = -g_RT[0] - g_RT[44] + g_RT[47];
    Kc[247] = -g_RT[0] - g_RT[45] + g_RT[48];
    Kc[248] = -g_RT[3] + g_RT[7] - g_RT[45] + g_RT[48];
    Kc[249] = -g_RT[2] + g_RT[3] - g_RT[14] - g_RT[34] + g_RT[48];
    Kc[250] = g_RT[0] - g_RT[23] - g_RT[24] + g_RT[49];
    Kc[251] = g_RT[0] - g_RT[11] - g_RT[36] + g_RT[49];
    Kc[252] = g_RT[0] - g_RT[4] - g_RT[48] + g_RT[49];
    Kc[253] = g_RT[1] - g_RT[13] - g_RT[37] + g_RT[49];
    Kc[254] = g_RT[23] + g_RT[24] - g_RT[50];
    Kc[255] = g_RT[2] - g_RT[5] - g_RT[49] + g_RT[50];
    Kc[256] = -g_RT[3] + g_RT[7] - g_RT[49] + g_RT[50];
    Kc[257] = -g_RT[2] + g_RT[3] - g_RT[14] - g_RT[37] + g_RT[50];
    Kc[258] = g_RT[11] - g_RT[12] - g_RT[49] + g_RT[50];
    Kc[259] = -3*g_RT[23] - 2*g_RT[37] + g_RT[51];
    Kc[260] = -2*g_RT[23] - 2*g_RT[50] + g_RT[51];
    Kc[261] = g_RT[0] - g_RT[4] - 4*g_RT[23] - g_RT[50] + g_RT[51];
    Kc[262] = g_RT[0] - g_RT[4] - 2*g_RT[23] - g_RT[49] - g_RT[50] + g_RT[51];
    Kc[263] = g_RT[0] - g_RT[4] - g_RT[36] - g_RT[37] + g_RT[51] - g_RT[52];
    Kc[264] = g_RT[0] - g_RT[4] - 2*g_RT[23] - g_RT[37] + g_RT[51] - g_RT[54];
    Kc[265] = g_RT[0] - g_RT[4] - g_RT[23] - g_RT[50] + g_RT[51] - g_RT[52];
    Kc[266] = g_RT[11] - g_RT[12] - 4*g_RT[23] - g_RT[50] + g_RT[51];
    Kc[267] = g_RT[11] - g_RT[12] - 2*g_RT[23] - g_RT[49] - g_RT[50] + g_RT[51];
    Kc[268] = g_RT[11] - g_RT[12] - g_RT[36] - g_RT[37] + g_RT[51] - g_RT[52];
    Kc[269] = g_RT[11] - g_RT[12] - 2*g_RT[23] - g_RT[37] + g_RT[51] - g_RT[54];
    Kc[270] = g_RT[11] - g_RT[12] - g_RT[23] - g_RT[50] + g_RT[51] - g_RT[52];
    Kc[271] = g_RT[1] - g_RT[2] - 4*g_RT[23] - g_RT[50] + g_RT[51];
    Kc[272] = g_RT[1] - g_RT[2] - 2*g_RT[23] - g_RT[49] - g_RT[50] + g_RT[51];
    Kc[273] = g_RT[1] - g_RT[2] - g_RT[36] - g_RT[37] + g_RT[51] - g_RT[52];
    Kc[274] = g_RT[1] - g_RT[2] - 2*g_RT[23] - g_RT[37] + g_RT[51] - g_RT[54];
    Kc[275] = g_RT[1] - g_RT[2] - g_RT[23] - g_RT[50] + g_RT[51] - g_RT[52];
    Kc[276] = g_RT[2] - g_RT[5] - 4*g_RT[23] - g_RT[50] + g_RT[51];
    Kc[277] = g_RT[2] - g_RT[5] - 2*g_RT[23] - g_RT[49] - g_RT[50] + g_RT[51];
    Kc[278] = g_RT[2] - g_RT[5] - g_RT[36] - g_RT[37] + g_RT[51] - g_RT[52];
    Kc[279] = g_RT[2] - g_RT[5] - 2*g_RT[23] - g_RT[37] + g_RT[51] - g_RT[54];
    Kc[280] = g_RT[2] - g_RT[5] - g_RT[23] - g_RT[50] + g_RT[51] - g_RT[52];
    Kc[281] = g_RT[0] - g_RT[23] - g_RT[50] + g_RT[52];
    Kc[282] = g_RT[0] - g_RT[4] + g_RT[52] - g_RT[53];
    Kc[283] = g_RT[0] - g_RT[23] - g_RT[37] + g_RT[54];
    Kc[284] = g_RT[0] - g_RT[4] - g_RT[23] - g_RT[34] + g_RT[54];
    Kc[285] = g_RT[0] - g_RT[11] - g_RT[23] - g_RT[34] + g_RT[53];
    Kc[286] = -g_RT[2] + g_RT[3] - g_RT[14] - g_RT[23] - g_RT[34] + g_RT[53];
    Kc[287] = g_RT[1] - g_RT[13] - g_RT[23] - g_RT[37] + g_RT[52];
    Kc[288] = g_RT[1] - g_RT[13] - g_RT[50] + g_RT[54];

    for (int i=0; i<289; ++i) {
        Kc[i] = exp(Kc[i]);
    };

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 * invT;
    double refCinv = 1 / refC;

    Kc[0] *= refCinv;
    Kc[1] *= refCinv;
    Kc[2] *= refCinv;
    Kc[3] *= refCinv;
    Kc[4] *= refCinv;
    Kc[5] *= refCinv;
    Kc[6] *= refCinv;
    Kc[7] *= refCinv;
    Kc[8] *= refCinv;
    Kc[9] *= refC;
    Kc[10] *= refCinv;
    Kc[11] *= refCinv;
    Kc[12] *= refCinv;
    Kc[13] *= refCinv;
    Kc[14] *= refCinv;
    Kc[15] *= refC;
    Kc[16] *= refCinv;
    Kc[17] *= refCinv;
    Kc[18] *= refCinv;
    Kc[19] *= refCinv;
    Kc[20] *= refCinv;
    Kc[21] *= refCinv;
    Kc[22] *= refCinv;
    Kc[25] *= refCinv;
    Kc[26] *= refCinv;
    Kc[27] *= refCinv;
    Kc[28] *= refCinv;
    Kc[29] *= refCinv;
    Kc[30] *= refC;
    Kc[35] *= refCinv;
    Kc[56] *= refC;
    Kc[69] *= refC;
    Kc[75] *= refC;
    Kc[125] *= refC;
    Kc[126] *= refC;
    Kc[151] *= refCinv;
    Kc[153] *= refC;
    Kc[159] *= refC;
    Kc[161] *= refC;
    Kc[177] *= refCinv;
    Kc[178] *= refCinv;
    Kc[185] *= refC;
    Kc[192] *= refCinv;
    Kc[196] *= refCinv;
    Kc[197] *= refCinv;
    Kc[202] *= refCinv;
    Kc[207] *= refC;
    Kc[210] *= refC;
    Kc[214] *= refC;
    Kc[218] *= refC;
    Kc[219] *= refC;
    Kc[226] *= refC;
    Kc[232] *= refC;
    Kc[237] *= refC;
    Kc[239] *= refCinv;
    Kc[245] *= refC;
    Kc[246] *= refC;
    Kc[247] *= refC;
    Kc[249] *= refC;
    Kc[254] *= refCinv;
    Kc[257] *= refC;
    Kc[259] *= refC*refC*refC*refC;
    Kc[260] *= refC*refC*refC;
    Kc[261] *= refC*refC*refC*refC;
    Kc[262] *= refC*refC*refC;
    Kc[263] *= refC*refC;
    Kc[264] *= refC*refC*refC;
    Kc[265] *= refC*refC;
    Kc[266] *= refC*refC*refC*refC;
    Kc[267] *= refC*refC*refC;
    Kc[268] *= refC*refC;
    Kc[269] *= refC*refC*refC;
    Kc[270] *= refC*refC;
    Kc[271] *= refC*refC*refC*refC;
    Kc[272] *= refC*refC*refC;
    Kc[273] *= refC*refC;
    Kc[274] *= refC*refC*refC;
    Kc[275] *= refC*refC;
    Kc[276] *= refC*refC*refC*refC;
    Kc[277] *= refC*refC*refC;
    Kc[278] *= refC*refC;
    Kc[279] *= refC*refC*refC;
    Kc[280] *= refC*refC;
    Kc[284] *= refC;
    Kc[285] *= refC;
    Kc[286] *= refC*refC;
    Kc[287] *= refC;

    return;
}

__device__ void comp_qfqr_d(double *  qf, double * qr, double * sc, double * tc, double invT)
{

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    qf[0] = sc[0]*sc[7];
    qr[0] = sc[3];

    /*reaction 2: 2 OH (+M) <=> H2O2 (+M) */
    qf[1] = sc[2]*sc[2];
    qr[1] = sc[6];

    /*reaction 3: CH2 + CO (+M) <=> CH2CO (+M) */
    qf[2] = sc[9]*sc[18];
    qr[2] = sc[27];

    /*reaction 4: CH2* + H2O (+M) <=> CH3OH (+M) */
    qf[3] = sc[5]*sc[10];
    qr[3] = sc[17];

    /*reaction 5: CH2O + H (+M) <=> CH2OH (+M) */
    qf[4] = sc[0]*sc[14];
    qr[4] = sc[16];

    /*reaction 6: CH2O + H (+M) <=> CH3O (+M) */
    qf[5] = sc[0]*sc[14];
    qr[5] = sc[15];

    /*reaction 7: CH3 + H (+M) <=> CH4 (+M) */
    qf[6] = sc[0]*sc[11];
    qr[6] = sc[12];

    /*reaction 8: CH3 + OH (+M) <=> CH3OH (+M) */
    qf[7] = sc[2]*sc[11];
    qr[7] = sc[17];

    /*reaction 9: 2 CH3 (+M) <=> C2H6 (+M) */
    qf[8] = sc[11]*sc[11];
    qr[8] = sc[25];

    /*reaction 10: C2H3 (+M) <=> C2H2 + H (+M) */
    qf[9] = sc[22];
    qr[9] = sc[0]*sc[21];

    /*reaction 11: CH2CO + H (+M) <=> CH2CHO (+M) */
    qf[10] = sc[0]*sc[27];
    qr[10] = sc[29];

    /*reaction 12: C2H3 + H (+M) <=> C2H4 (+M) */
    qf[11] = sc[0]*sc[22];
    qr[11] = sc[23];

    /*reaction 13: C2H3 + CH3 (+M) <=> C3H6 (+M) */
    qf[12] = sc[11]*sc[22];
    qr[12] = sc[36];

    /*reaction 14: CH3 + CO (+M) <=> CH3CO (+M) */
    qf[13] = sc[11]*sc[18];
    qr[13] = sc[28];

    /*reaction 15: CH3 + HCO (+M) <=> CH3CHO (+M) */
    qf[14] = sc[11]*sc[13];
    qr[14] = sc[30];

    /*reaction 16: C2H4 (+M) <=> H2 + C2H2 (+M) */
    qf[15] = sc[23];
    qr[15] = sc[4]*sc[21];

    /*reaction 17: C2H4 + H (+M) <=> C2H5 (+M) */
    qf[16] = sc[0]*sc[23];
    qr[16] = sc[24];

    /*reaction 18: C2H5 + H (+M) <=> C2H6 (+M) */
    qf[17] = sc[0]*sc[24];
    qr[17] = sc[25];

    /*reaction 19: C2H5 + C2H3 (+M) <=> C4H81 (+M) */
    qf[18] = sc[22]*sc[24];
    qr[18] = sc[49];

    /*reaction 20: aC3H5 + H (+M) <=> C3H6 (+M) */
    qf[19] = sc[0]*sc[34];
    qr[19] = sc[36];

    /*reaction 21: aC3H5 + CH3 (+M) <=> C4H81 (+M) */
    qf[20] = sc[11]*sc[34];
    qr[20] = sc[49];

    /*reaction 22: C3H6 + H (+M) <=> nC3H7 (+M) */
    qf[21] = sc[0]*sc[36];
    qr[21] = sc[37];

    /*reaction 23: C3H6 + H (+M) <=> iC3H7 (+M) */
    qf[22] = sc[0]*sc[36];
    qr[22] = sc[38];

    /*reaction 24: C6H12 + H (+M) <=> C3H6 + nC3H7 (+M) */
    qf[23] = sc[0]*sc[52];
    qr[23] = sc[36]*sc[37];

    /*reaction 25: C5H10 + H (+M) <=> C3H6 + C2H5 (+M) */
    qf[24] = sc[0]*sc[54];
    qr[24] = sc[24]*sc[36];

    /*reaction 26: CO + O (+M) <=> CO2 (+M) */
    qf[25] = sc[1]*sc[18];
    qr[25] = sc[19];

    /*reaction 27: 2 H + M <=> H2 + M */
    qf[26] = sc[0]*sc[0];
    qr[26] = sc[4];

    /*reaction 28: H + OH + M <=> H2O + M */
    qf[27] = sc[0]*sc[2];
    qr[27] = sc[5];

    /*reaction 29: O + H + M <=> OH + M */
    qf[28] = sc[0]*sc[1];
    qr[28] = sc[2];

    /*reaction 30: 2 O + M <=> O2 + M */
    qf[29] = sc[1]*sc[1];
    qr[29] = sc[7];

    /*reaction 31: HCO + M <=> CO + H + M */
    qf[30] = sc[13];
    qr[30] = sc[0]*sc[18];

    /*reaction 32: H + O2 <=> O + OH */
    qf[31] = sc[0]*sc[7];
    qr[31] = sc[1]*sc[2];

    /*reaction 33: O + H2 <=> H + OH */
    qf[32] = sc[1]*sc[4];
    qr[32] = sc[0]*sc[2];

    /*reaction 34: OH + H2 <=> H + H2O */
    qf[33] = sc[2]*sc[4];
    qr[33] = sc[0]*sc[5];

    /*reaction 35: 2 OH <=> O + H2O */
    qf[34] = sc[2]*sc[2];
    qr[34] = sc[1]*sc[5];

    /*reaction 36: 2 H + H2O <=> H2 + H2O */
    qf[35] = sc[0]*sc[0]*sc[5];
    qr[35] = sc[4]*sc[5];

    /*reaction 37: H2 + O2 <=> HO2 + H */
    qf[36] = sc[4]*sc[7];
    qr[36] = sc[0]*sc[3];

    /*reaction 38: HO2 + H <=> O + H2O */
    qf[37] = sc[0]*sc[3];
    qr[37] = sc[1]*sc[5];

    /*reaction 39: HO2 + H <=> 2 OH */
    qf[38] = sc[0]*sc[3];
    qr[38] = sc[2]*sc[2];

    /*reaction 40: HO2 + O <=> OH + O2 */
    qf[39] = sc[1]*sc[3];
    qr[39] = sc[2]*sc[7];

    /*reaction 41: HO2 + OH <=> O2 + H2O */
    qf[40] = sc[2]*sc[3];
    qr[40] = sc[5]*sc[7];

    /*reaction 42: HO2 + OH <=> O2 + H2O */
    qf[41] = sc[2]*sc[3];
    qr[41] = sc[5]*sc[7];

    /*reaction 43: 2 HO2 <=> O2 + H2O2 */
    qf[42] = sc[3]*sc[3];
    qr[42] = sc[6]*sc[7];

    /*reaction 44: 2 HO2 <=> O2 + H2O2 */
    qf[43] = sc[3]*sc[3];
    qr[43] = sc[6]*sc[7];

    /*reaction 45: H2O2 + H <=> HO2 + H2 */
    qf[44] = sc[0]*sc[6];
    qr[44] = sc[3]*sc[4];

    /*reaction 46: H2O2 + H <=> OH + H2O */
    qf[45] = sc[0]*sc[6];
    qr[45] = sc[2]*sc[5];

    /*reaction 47: H2O2 + O <=> OH + HO2 */
    qf[46] = sc[1]*sc[6];
    qr[46] = sc[2]*sc[3];

    /*reaction 48: H2O2 + OH <=> HO2 + H2O */
    qf[47] = sc[2]*sc[6];
    qr[47] = sc[3]*sc[5];

    /*reaction 49: H2O2 + OH <=> HO2 + H2O */
    qf[48] = sc[2]*sc[6];
    qr[48] = sc[3]*sc[5];

    /*reaction 50: CO + OH <=> CO2 + H */
    qf[49] = sc[2]*sc[18];
    qr[49] = sc[0]*sc[19];

    /*reaction 51: CO + OH <=> CO2 + H */
    qf[50] = sc[2]*sc[18];
    qr[50] = sc[0]*sc[19];

    /*reaction 52: CO + HO2 <=> CO2 + OH */
    qf[51] = sc[3]*sc[18];
    qr[51] = sc[2]*sc[19];

    /*reaction 53: HCO + H <=> CO + H2 */
    qf[52] = sc[0]*sc[13];
    qr[52] = sc[4]*sc[18];

    /*reaction 54: HCO + O <=> CO + OH */
    qf[53] = sc[1]*sc[13];
    qr[53] = sc[2]*sc[18];

    /*reaction 55: HCO + O <=> CO2 + H */
    qf[54] = sc[1]*sc[13];
    qr[54] = sc[0]*sc[19];

    /*reaction 56: HCO + OH <=> CO + H2O */
    qf[55] = sc[2]*sc[13];
    qr[55] = sc[5]*sc[18];

    /*reaction 57: HCO + H2O <=> CO + H + H2O */
    qf[56] = sc[5]*sc[13];
    qr[56] = sc[0]*sc[5]*sc[18];

    /*reaction 58: HCO + O2 <=> CO + HO2 */
    qf[57] = sc[7]*sc[13];
    qr[57] = sc[3]*sc[18];

    /*reaction 59: CH + O <=> CO + H */
    qf[58] = sc[1]*sc[8];
    qr[58] = sc[0]*sc[18];

    /*reaction 60: CH + OH <=> HCO + H */
    qf[59] = sc[2]*sc[8];
    qr[59] = sc[0]*sc[13];

    /*reaction 61: CH + H2 <=> CH2 + H */
    qf[60] = sc[4]*sc[8];
    qr[60] = sc[0]*sc[9];

    /*reaction 62: CH + H2O <=> CH2O + H */
    qf[61] = sc[5]*sc[8];
    qr[61] = sc[0]*sc[14];

    /*reaction 63: CH + O2 <=> HCO + O */
    qf[62] = sc[7]*sc[8];
    qr[62] = sc[1]*sc[13];

    /*reaction 64: CH + CO2 <=> HCO + CO */
    qf[63] = sc[8]*sc[19];
    qr[63] = sc[13]*sc[18];

    /*reaction 65: CH2 + O <=> HCO + H */
    qf[64] = sc[1]*sc[9];
    qr[64] = sc[0]*sc[13];

    /*reaction 66: CH2 + OH <=> CH2O + H */
    qf[65] = sc[2]*sc[9];
    qr[65] = sc[0]*sc[14];

    /*reaction 67: CH2 + OH <=> CH + H2O */
    qf[66] = sc[2]*sc[9];
    qr[66] = sc[5]*sc[8];

    /*reaction 68: CH2 + H2 <=> H + CH3 */
    qf[67] = sc[4]*sc[9];
    qr[67] = sc[0]*sc[11];

    /*reaction 69: CH2 + O2 <=> HCO + OH */
    qf[68] = sc[7]*sc[9];
    qr[68] = sc[2]*sc[13];

    /*reaction 70: CH2 + O2 <=> CO2 + 2 H */
    qf[69] = sc[7]*sc[9];
    qr[69] = sc[0]*sc[0]*sc[19];

    /*reaction 71: CH2 + HO2 <=> CH2O + OH */
    qf[70] = sc[3]*sc[9];
    qr[70] = sc[2]*sc[14];

    /*reaction 72: CH2* + N2 <=> CH2 + N2 */
    qf[71] = sc[10]*sc[55];
    qr[71] = sc[9]*sc[55];

    /*reaction 73: CH2* + H <=> CH + H2 */
    qf[72] = sc[0]*sc[10];
    qr[72] = sc[4]*sc[8];

    /*reaction 74: CH2* + OH <=> CH2O + H */
    qf[73] = sc[2]*sc[10];
    qr[73] = sc[0]*sc[14];

    /*reaction 75: CH2* + H2 <=> CH3 + H */
    qf[74] = sc[4]*sc[10];
    qr[74] = sc[0]*sc[11];

    /*reaction 76: CH2* + O2 <=> H + OH + CO */
    qf[75] = sc[7]*sc[10];
    qr[75] = sc[0]*sc[2]*sc[18];

    /*reaction 77: CH2* + O2 <=> CO + H2O */
    qf[76] = sc[7]*sc[10];
    qr[76] = sc[5]*sc[18];

    /*reaction 78: CH2* + H2O <=> CH2 + H2O */
    qf[77] = sc[5]*sc[10];
    qr[77] = sc[5]*sc[9];

    /*reaction 79: CH2* + CO <=> CH2 + CO */
    qf[78] = sc[10]*sc[18];
    qr[78] = sc[9]*sc[18];

    /*reaction 80: CH2* + CO2 <=> CH2 + CO2 */
    qf[79] = sc[10]*sc[19];
    qr[79] = sc[9]*sc[19];

    /*reaction 81: CH2* + CO2 <=> CH2O + CO */
    qf[80] = sc[10]*sc[19];
    qr[80] = sc[14]*sc[18];

    /*reaction 82: CH2O + H <=> HCO + H2 */
    qf[81] = sc[0]*sc[14];
    qr[81] = sc[4]*sc[13];

    /*reaction 83: CH2O + O <=> HCO + OH */
    qf[82] = sc[1]*sc[14];
    qr[82] = sc[2]*sc[13];

    /*reaction 84: CH2O + OH <=> HCO + H2O */
    qf[83] = sc[2]*sc[14];
    qr[83] = sc[5]*sc[13];

    /*reaction 85: CH2O + O2 <=> HCO + HO2 */
    qf[84] = sc[7]*sc[14];
    qr[84] = sc[3]*sc[13];

    /*reaction 86: CH2O + HO2 <=> HCO + H2O2 */
    qf[85] = sc[3]*sc[14];
    qr[85] = sc[6]*sc[13];

    /*reaction 87: CH2O + CH <=> CH2CO + H */
    qf[86] = sc[8]*sc[14];
    qr[86] = sc[0]*sc[27];

    /*reaction 88: CH3 + O <=> CH2O + H */
    qf[87] = sc[1]*sc[11];
    qr[87] = sc[0]*sc[14];

    /*reaction 89: CH3 + OH <=> CH2 + H2O */
    qf[88] = sc[2]*sc[11];
    qr[88] = sc[5]*sc[9];

    /*reaction 90: CH3 + OH <=> CH2* + H2O */
    qf[89] = sc[2]*sc[11];
    qr[89] = sc[5]*sc[10];

    /*reaction 91: CH3 + O2 <=> O + CH3O */
    qf[90] = sc[7]*sc[11];
    qr[90] = sc[1]*sc[15];

    /*reaction 92: CH3 + O2 <=> OH + CH2O */
    qf[91] = sc[7]*sc[11];
    qr[91] = sc[2]*sc[14];

    /*reaction 93: CH3 + HO2 <=> CH4 + O2 */
    qf[92] = sc[3]*sc[11];
    qr[92] = sc[7]*sc[12];

    /*reaction 94: CH3 + HO2 <=> CH3O + OH */
    qf[93] = sc[3]*sc[11];
    qr[93] = sc[2]*sc[15];

    /*reaction 95: CH3 + CH <=> C2H3 + H */
    qf[94] = sc[8]*sc[11];
    qr[94] = sc[0]*sc[22];

    /*reaction 96: CH3 + HCO <=> CH4 + CO */
    qf[95] = sc[11]*sc[13];
    qr[95] = sc[12]*sc[18];

    /*reaction 97: CH3 + CH2O <=> CH4 + HCO */
    qf[96] = sc[11]*sc[14];
    qr[96] = sc[12]*sc[13];

    /*reaction 98: CH3 + CH2 <=> C2H4 + H */
    qf[97] = sc[9]*sc[11];
    qr[97] = sc[0]*sc[23];

    /*reaction 99: 2 CH3 <=> H + C2H5 */
    qf[98] = sc[11]*sc[11];
    qr[98] = sc[0]*sc[24];

    /*reaction 100: CH3 + HCCO <=> C2H4 + CO */
    qf[99] = sc[11]*sc[26];
    qr[99] = sc[18]*sc[23];

    /*reaction 101: CH3O + H <=> CH2O + H2 */
    qf[100] = sc[0]*sc[15];
    qr[100] = sc[4]*sc[14];

    /*reaction 102: CH3O + H <=> CH3 + OH */
    qf[101] = sc[0]*sc[15];
    qr[101] = sc[2]*sc[11];

    /*reaction 103: CH3O + H <=> CH2* + H2O */
    qf[102] = sc[0]*sc[15];
    qr[102] = sc[5]*sc[10];

    /*reaction 104: CH3O + OH <=> CH2O + H2O */
    qf[103] = sc[2]*sc[15];
    qr[103] = sc[5]*sc[14];

    /*reaction 105: CH3O + O2 <=> CH2O + HO2 */
    qf[104] = sc[7]*sc[15];
    qr[104] = sc[3]*sc[14];

    /*reaction 106: CH2OH + H <=> CH2O + H2 */
    qf[105] = sc[0]*sc[16];
    qr[105] = sc[4]*sc[14];

    /*reaction 107: CH2OH + H <=> CH3 + OH */
    qf[106] = sc[0]*sc[16];
    qr[106] = sc[2]*sc[11];

    /*reaction 108: CH2OH + H <=> CH2* + H2O */
    qf[107] = sc[0]*sc[16];
    qr[107] = sc[5]*sc[10];

    /*reaction 109: CH2OH + O2 <=> CH2O + HO2 */
    qf[108] = sc[7]*sc[16];
    qr[108] = sc[3]*sc[14];

    /*reaction 110: CH4 + H <=> CH3 + H2 */
    qf[109] = sc[0]*sc[12];
    qr[109] = sc[4]*sc[11];

    /*reaction 111: CH4 + O <=> CH3 + OH */
    qf[110] = sc[1]*sc[12];
    qr[110] = sc[2]*sc[11];

    /*reaction 112: CH4 + OH <=> CH3 + H2O */
    qf[111] = sc[2]*sc[12];
    qr[111] = sc[5]*sc[11];

    /*reaction 113: CH4 + CH <=> C2H4 + H */
    qf[112] = sc[8]*sc[12];
    qr[112] = sc[0]*sc[23];

    /*reaction 114: CH4 + CH2 <=> 2 CH3 */
    qf[113] = sc[9]*sc[12];
    qr[113] = sc[11]*sc[11];

    /*reaction 115: CH4 + CH2* <=> 2 CH3 */
    qf[114] = sc[10]*sc[12];
    qr[114] = sc[11]*sc[11];

    /*reaction 116: CH3OH + H <=> CH2OH + H2 */
    qf[115] = sc[0]*sc[17];
    qr[115] = sc[4]*sc[16];

    /*reaction 117: CH3OH + H <=> CH3O + H2 */
    qf[116] = sc[0]*sc[17];
    qr[116] = sc[4]*sc[15];

    /*reaction 118: CH3OH + O <=> CH2OH + OH */
    qf[117] = sc[1]*sc[17];
    qr[117] = sc[2]*sc[16];

    /*reaction 119: CH3OH + OH <=> CH2OH + H2O */
    qf[118] = sc[2]*sc[17];
    qr[118] = sc[5]*sc[16];

    /*reaction 120: CH3OH + OH <=> CH3O + H2O */
    qf[119] = sc[2]*sc[17];
    qr[119] = sc[5]*sc[15];

    /*reaction 121: C2H + O <=> CH + CO */
    qf[120] = sc[1]*sc[20];
    qr[120] = sc[8]*sc[18];

    /*reaction 122: C2H + OH <=> H + HCCO */
    qf[121] = sc[2]*sc[20];
    qr[121] = sc[0]*sc[26];

    /*reaction 123: C2H + O2 <=> HCO + CO */
    qf[122] = sc[7]*sc[20];
    qr[122] = sc[13]*sc[18];

    /*reaction 124: C2H + H2 <=> H + C2H2 */
    qf[123] = sc[4]*sc[20];
    qr[123] = sc[0]*sc[21];

    /*reaction 125: HCCO + H <=> CH2* + CO */
    qf[124] = sc[0]*sc[26];
    qr[124] = sc[10]*sc[18];

    /*reaction 126: HCCO + O <=> H + 2 CO */
    qf[125] = sc[1]*sc[26];
    qr[125] = sc[0]*sc[18]*sc[18];

    /*reaction 127: HCCO + O2 <=> OH + 2 CO */
    qf[126] = sc[7]*sc[26];
    qr[126] = sc[2]*sc[18]*sc[18];

    /*reaction 128: C2H2 + O <=> C2H + OH */
    qf[127] = sc[1]*sc[21];
    qr[127] = sc[2]*sc[20];

    /*reaction 129: C2H2 + O <=> CH2 + CO */
    qf[128] = sc[1]*sc[21];
    qr[128] = sc[9]*sc[18];

    /*reaction 130: C2H2 + O <=> HCCO + H */
    qf[129] = sc[1]*sc[21];
    qr[129] = sc[0]*sc[26];

    /*reaction 131: C2H2 + OH <=> CH2CO + H */
    qf[130] = sc[2]*sc[21];
    qr[130] = sc[0]*sc[27];

    /*reaction 132: C2H2 + OH <=> CH2CO + H */
    qf[131] = sc[2]*sc[21];
    qr[131] = sc[0]*sc[27];

    /*reaction 133: C2H2 + OH <=> C2H + H2O */
    qf[132] = sc[2]*sc[21];
    qr[132] = sc[5]*sc[20];

    /*reaction 134: C2H2 + HCO <=> C2H3 + CO */
    qf[133] = sc[13]*sc[21];
    qr[133] = sc[18]*sc[22];

    /*reaction 135: C2H2 + CH2 <=> C3H3 + H */
    qf[134] = sc[9]*sc[21];
    qr[134] = sc[0]*sc[31];

    /*reaction 136: C2H2 + CH2* <=> C3H3 + H */
    qf[135] = sc[10]*sc[21];
    qr[135] = sc[0]*sc[31];

    /*reaction 137: C2H2 + C2H <=> C4H2 + H */
    qf[136] = sc[20]*sc[21];
    qr[136] = sc[0]*sc[40];

    /*reaction 138: C2H2 + CH3 <=> pC3H4 + H */
    qf[137] = sc[11]*sc[21];
    qr[137] = sc[0]*sc[32];

    /*reaction 139: C2H2 + CH3 <=> aC3H4 + H */
    qf[138] = sc[11]*sc[21];
    qr[138] = sc[0]*sc[33];

    /*reaction 140: CH2CO + H <=> HCCO + H2 */
    qf[139] = sc[0]*sc[27];
    qr[139] = sc[4]*sc[26];

    /*reaction 141: CH2CO + H <=> CH3 + CO */
    qf[140] = sc[0]*sc[27];
    qr[140] = sc[11]*sc[18];

    /*reaction 142: CH2CO + OH <=> HCCO + H2O */
    qf[141] = sc[2]*sc[27];
    qr[141] = sc[5]*sc[26];

    /*reaction 143: C2H3 + H <=> C2H2 + H2 */
    qf[142] = sc[0]*sc[22];
    qr[142] = sc[4]*sc[21];

    /*reaction 144: C2H3 + O <=> CH2CO + H */
    qf[143] = sc[1]*sc[22];
    qr[143] = sc[0]*sc[27];

    /*reaction 145: C2H3 + O <=> CH3 + CO */
    qf[144] = sc[1]*sc[22];
    qr[144] = sc[11]*sc[18];

    /*reaction 146: C2H3 + OH <=> C2H2 + H2O */
    qf[145] = sc[2]*sc[22];
    qr[145] = sc[5]*sc[21];

    /*reaction 147: C2H3 + O2 <=> C2H2 + HO2 */
    qf[146] = sc[7]*sc[22];
    qr[146] = sc[3]*sc[21];

    /*reaction 148: C2H3 + O2 <=> CH2CHO + O */
    qf[147] = sc[7]*sc[22];
    qr[147] = sc[1]*sc[29];

    /*reaction 149: C2H3 + O2 <=> HCO + CH2O */
    qf[148] = sc[7]*sc[22];
    qr[148] = sc[13]*sc[14];

    /*reaction 150: C2H3 + HO2 <=> CH2CHO + OH */
    qf[149] = sc[3]*sc[22];
    qr[149] = sc[2]*sc[29];

    /*reaction 151: C2H3 + HCO <=> C2H4 + CO */
    qf[150] = sc[13]*sc[22];
    qr[150] = sc[18]*sc[23];

    /*reaction 152: C2H3 + HCO <=> C2H3CHO */
    qf[151] = sc[13]*sc[22];
    qr[151] = sc[39];

    /*reaction 153: C2H3 + CH3 <=> aC3H5 + H */
    qf[152] = sc[11]*sc[22];
    qr[152] = sc[0]*sc[34];

    /*reaction 154: CH2CHO <=> CH3 + CO */
    qf[153] = sc[29];
    qr[153] = sc[11]*sc[18];

    /*reaction 155: CH2CHO + H <=> CH3CO + H */
    qf[154] = sc[0]*sc[29];
    qr[154] = sc[0]*sc[28];

    /*reaction 156: CH2CHO + H <=> CH3 + HCO */
    qf[155] = sc[0]*sc[29];
    qr[155] = sc[11]*sc[13];

    /*reaction 157: CH2CHO + H <=> CH2CO + H2 */
    qf[156] = sc[0]*sc[29];
    qr[156] = sc[4]*sc[27];

    /*reaction 158: CH2CHO + OH <=> CH2CO + H2O */
    qf[157] = sc[2]*sc[29];
    qr[157] = sc[5]*sc[27];

    /*reaction 159: CH2CHO + O2 <=> CH2CO + HO2 */
    qf[158] = sc[7]*sc[29];
    qr[158] = sc[3]*sc[27];

    /*reaction 160: CH2CHO + O2 <=> CH2O + CO + OH */
    qf[159] = sc[7]*sc[29];
    qr[159] = sc[2]*sc[14]*sc[18];

    /*reaction 161: CH3CO + H <=> CH3 + HCO */
    qf[160] = sc[0]*sc[28];
    qr[160] = sc[11]*sc[13];

    /*reaction 162: CH3CO + HO2 <=> CH3 + CO2 + OH */
    qf[161] = sc[3]*sc[28];
    qr[161] = sc[2]*sc[11]*sc[19];

    /*reaction 163: CH3CHO + H <=> CH3CO + H2 */
    qf[162] = sc[0]*sc[30];
    qr[162] = sc[4]*sc[28];

    /*reaction 164: CH3CHO + OH <=> CH3CO + H2O */
    qf[163] = sc[2]*sc[30];
    qr[163] = sc[5]*sc[28];

    /*reaction 165: CH3CHO + CH3 <=> CH3CO + CH4 */
    qf[164] = sc[11]*sc[30];
    qr[164] = sc[12]*sc[28];

    /*reaction 166: CH3CHO + O2 <=> CH3CO + HO2 */
    qf[165] = sc[7]*sc[30];
    qr[165] = sc[3]*sc[28];

    /*reaction 167: C2H4 + H <=> C2H3 + H2 */
    qf[166] = sc[0]*sc[23];
    qr[166] = sc[4]*sc[22];

    /*reaction 168: C2H4 + O <=> C2H3 + OH */
    qf[167] = sc[1]*sc[23];
    qr[167] = sc[2]*sc[22];

    /*reaction 169: C2H4 + O <=> CH3 + HCO */
    qf[168] = sc[1]*sc[23];
    qr[168] = sc[11]*sc[13];

    /*reaction 170: C2H4 + O <=> CH2 + CH2O */
    qf[169] = sc[1]*sc[23];
    qr[169] = sc[9]*sc[14];

    /*reaction 171: C2H4 + OH <=> C2H3 + H2O */
    qf[170] = sc[2]*sc[23];
    qr[170] = sc[5]*sc[22];

    /*reaction 172: C2H4 + HCO <=> C2H5 + CO */
    qf[171] = sc[13]*sc[23];
    qr[171] = sc[18]*sc[24];

    /*reaction 173: C2H4 + CH <=> aC3H4 + H */
    qf[172] = sc[8]*sc[23];
    qr[172] = sc[0]*sc[33];

    /*reaction 174: C2H4 + CH <=> pC3H4 + H */
    qf[173] = sc[8]*sc[23];
    qr[173] = sc[0]*sc[32];

    /*reaction 175: C2H4 + CH2 <=> aC3H5 + H */
    qf[174] = sc[9]*sc[23];
    qr[174] = sc[0]*sc[34];

    /*reaction 176: C2H4 + CH2* <=> aC3H5 + H */
    qf[175] = sc[10]*sc[23];
    qr[175] = sc[0]*sc[34];

    /*reaction 177: C2H4 + CH3 <=> C2H3 + CH4 */
    qf[176] = sc[11]*sc[23];
    qr[176] = sc[12]*sc[22];

    /*reaction 178: C2H4 + CH3 <=> nC3H7 */
    qf[177] = sc[11]*sc[23];
    qr[177] = sc[37];

    /*reaction 179: C2H4 + C2H3 <=> C4H7 */
    qf[178] = sc[22]*sc[23];
    qr[178] = sc[48];

    /*reaction 180: C2H5 + H <=> C2H4 + H2 */
    qf[179] = sc[0]*sc[24];
    qr[179] = sc[4]*sc[23];

    /*reaction 181: C2H5 + O <=> CH3 + CH2O */
    qf[180] = sc[1]*sc[24];
    qr[180] = sc[11]*sc[14];

    /*reaction 182: C2H5 + O <=> CH3CHO + H */
    qf[181] = sc[1]*sc[24];
    qr[181] = sc[0]*sc[30];

    /*reaction 183: C2H5 + O2 <=> C2H4 + HO2 */
    qf[182] = sc[7]*sc[24];
    qr[182] = sc[3]*sc[23];

    /*reaction 184: C2H5 + HO2 <=> C2H6 + O2 */
    qf[183] = sc[3]*sc[24];
    qr[183] = sc[7]*sc[25];

    /*reaction 185: C2H5 + HO2 <=> C2H4 + H2O2 */
    qf[184] = sc[3]*sc[24];
    qr[184] = sc[6]*sc[23];

    /*reaction 186: C2H5 + HO2 <=> CH3 + CH2O + OH */
    qf[185] = sc[3]*sc[24];
    qr[185] = sc[2]*sc[11]*sc[14];

    /*reaction 187: C2H5 + C2H3 <=> aC3H5 + CH3 */
    qf[186] = sc[22]*sc[24];
    qr[186] = sc[11]*sc[34];

    /*reaction 188: C2H6 + H <=> C2H5 + H2 */
    qf[187] = sc[0]*sc[25];
    qr[187] = sc[4]*sc[24];

    /*reaction 189: C2H6 + O <=> C2H5 + OH */
    qf[188] = sc[1]*sc[25];
    qr[188] = sc[2]*sc[24];

    /*reaction 190: C2H6 + OH <=> C2H5 + H2O */
    qf[189] = sc[2]*sc[25];
    qr[189] = sc[5]*sc[24];

    /*reaction 191: C2H6 + CH2* <=> C2H5 + CH3 */
    qf[190] = sc[10]*sc[25];
    qr[190] = sc[11]*sc[24];

    /*reaction 192: C2H6 + CH3 <=> C2H5 + CH4 */
    qf[191] = sc[11]*sc[25];
    qr[191] = sc[12]*sc[24];

    /*reaction 193: C3H3 + H <=> pC3H4 */
    qf[192] = sc[0]*sc[31];
    qr[192] = sc[32];

    /*reaction 194: C3H3 + O <=> CH2O + C2H */
    qf[193] = sc[1]*sc[31];
    qr[193] = sc[14]*sc[20];

    /*reaction 195: C3H3 + O2 <=> CH2CO + HCO */
    qf[194] = sc[7]*sc[31];
    qr[194] = sc[13]*sc[27];

    /*reaction 196: C3H3 + HO2 <=> pC3H4 + O2 */
    qf[195] = sc[3]*sc[31];
    qr[195] = sc[7]*sc[32];

    /*reaction 197: aC3H4 + H <=> CH3CCH2 */
    qf[196] = sc[0]*sc[33];
    qr[196] = sc[35];

    /*reaction 198: aC3H4 + H <=> aC3H5 */
    qf[197] = sc[0]*sc[33];
    qr[197] = sc[34];

    /*reaction 199: aC3H4 + O <=> C2H4 + CO */
    qf[198] = sc[1]*sc[33];
    qr[198] = sc[18]*sc[23];

    /*reaction 200: aC3H4 + OH <=> C3H3 + H2O */
    qf[199] = sc[2]*sc[33];
    qr[199] = sc[5]*sc[31];

    /*reaction 201: pC3H4 <=> aC3H4 */
    qf[200] = sc[32];
    qr[200] = sc[33];

    /*reaction 202: pC3H4 + H <=> aC3H4 + H */
    qf[201] = sc[0]*sc[32];
    qr[201] = sc[0]*sc[33];

    /*reaction 203: pC3H4 + H <=> CH3CCH2 */
    qf[202] = sc[0]*sc[32];
    qr[202] = sc[35];

    /*reaction 204: pC3H4 + O <=> C2H4 + CO */
    qf[203] = sc[1]*sc[32];
    qr[203] = sc[18]*sc[23];

    /*reaction 205: pC3H4 + OH <=> C3H3 + H2O */
    qf[204] = sc[2]*sc[32];
    qr[204] = sc[5]*sc[31];

    /*reaction 206: aC3H5 + H <=> aC3H4 + H2 */
    qf[205] = sc[0]*sc[34];
    qr[205] = sc[4]*sc[33];

    /*reaction 207: aC3H5 + O <=> C2H3CHO + H */
    qf[206] = sc[1]*sc[34];
    qr[206] = sc[0]*sc[39];

    /*reaction 208: aC3H5 + OH <=> C2H3CHO + 2 H */
    qf[207] = sc[2]*sc[34];
    qr[207] = sc[0]*sc[0]*sc[39];

    /*reaction 209: aC3H5 + OH <=> aC3H4 + H2O */
    qf[208] = sc[2]*sc[34];
    qr[208] = sc[5]*sc[33];

    /*reaction 210: aC3H5 + HO2 <=> C3H6 + O2 */
    qf[209] = sc[3]*sc[34];
    qr[209] = sc[7]*sc[36];

    /*reaction 211: aC3H5 + HO2 <=> OH + C2H3 + CH2O */
    qf[210] = sc[3]*sc[34];
    qr[210] = sc[2]*sc[14]*sc[22];

    /*reaction 212: aC3H5 + HCO <=> C3H6 + CO */
    qf[211] = sc[13]*sc[34];
    qr[211] = sc[18]*sc[36];

    /*reaction 213: aC3H5 + CH3 <=> aC3H4 + CH4 */
    qf[212] = sc[11]*sc[34];
    qr[212] = sc[12]*sc[33];

    /*reaction 214: CH3CCH2 + O2 <=> CH3CO + CH2O */
    qf[213] = sc[7]*sc[35];
    qr[213] = sc[14]*sc[28];

    /*reaction 215: CH3CCH2 + HO2 <=> CH3 + CH2CO + OH */
    qf[214] = sc[3]*sc[35];
    qr[214] = sc[2]*sc[11]*sc[27];

    /*reaction 216: C3H6 + H <=> C2H4 + CH3 */
    qf[215] = sc[0]*sc[36];
    qr[215] = sc[11]*sc[23];

    /*reaction 217: C3H6 + H <=> aC3H5 + H2 */
    qf[216] = sc[0]*sc[36];
    qr[216] = sc[4]*sc[34];

    /*reaction 218: C3H6 + H <=> CH3CCH2 + H2 */
    qf[217] = sc[0]*sc[36];
    qr[217] = sc[4]*sc[35];

    /*reaction 219: C3H6 + O <=> CH2CO + CH3 + H */
    qf[218] = sc[1]*sc[36];
    qr[218] = sc[0]*sc[11]*sc[27];

    /*reaction 220: C3H6 + O <=> C2H3CHO + 2 H */
    qf[219] = sc[1]*sc[36];
    qr[219] = sc[0]*sc[0]*sc[39];

    /*reaction 221: C3H6 + O <=> C2H5 + HCO */
    qf[220] = sc[1]*sc[36];
    qr[220] = sc[13]*sc[24];

    /*reaction 222: C3H6 + O <=> aC3H5 + OH */
    qf[221] = sc[1]*sc[36];
    qr[221] = sc[2]*sc[34];

    /*reaction 223: C3H6 + O <=> CH3CCH2 + OH */
    qf[222] = sc[1]*sc[36];
    qr[222] = sc[2]*sc[35];

    /*reaction 224: C3H6 + OH <=> aC3H5 + H2O */
    qf[223] = sc[2]*sc[36];
    qr[223] = sc[5]*sc[34];

    /*reaction 225: C3H6 + OH <=> CH3CCH2 + H2O */
    qf[224] = sc[2]*sc[36];
    qr[224] = sc[5]*sc[35];

    /*reaction 226: C3H6 + CH3 <=> aC3H5 + CH4 */
    qf[225] = sc[11]*sc[36];
    qr[225] = sc[12]*sc[34];

    /*reaction 227: C2H3CHO + O <=> C2H3 + OH + CO */
    qf[226] = sc[1]*sc[39];
    qr[226] = sc[2]*sc[18]*sc[22];

    /*reaction 228: C2H3CHO + O <=> CH2O + CH2CO */
    qf[227] = sc[1]*sc[39];
    qr[227] = sc[14]*sc[27];

    /*reaction 229: iC3H7 + H <=> CH3 + C2H5 */
    qf[228] = sc[0]*sc[38];
    qr[228] = sc[11]*sc[24];

    /*reaction 230: iC3H7 + O <=> CH3CHO + CH3 */
    qf[229] = sc[1]*sc[38];
    qr[229] = sc[11]*sc[30];

    /*reaction 231: iC3H7 + OH <=> C3H6 + H2O */
    qf[230] = sc[2]*sc[38];
    qr[230] = sc[5]*sc[36];

    /*reaction 232: iC3H7 + O2 <=> C3H6 + HO2 */
    qf[231] = sc[7]*sc[38];
    qr[231] = sc[3]*sc[36];

    /*reaction 233: iC3H7 + HO2 <=> CH3CHO + CH3 + OH */
    qf[232] = sc[3]*sc[38];
    qr[232] = sc[2]*sc[11]*sc[30];

    /*reaction 234: iC3H7 + CH3 <=> CH4 + C3H6 */
    qf[233] = sc[11]*sc[38];
    qr[233] = sc[12]*sc[36];

    /*reaction 235: nC3H7 + H <=> C2H5 + CH3 */
    qf[234] = sc[0]*sc[37];
    qr[234] = sc[11]*sc[24];

    /*reaction 236: nC3H7 + OH <=> C3H6 + H2O */
    qf[235] = sc[2]*sc[37];
    qr[235] = sc[5]*sc[36];

    /*reaction 237: nC3H7 + O2 <=> C3H6 + HO2 */
    qf[236] = sc[7]*sc[37];
    qr[236] = sc[3]*sc[36];

    /*reaction 238: nC3H7 + HO2 <=> C2H5 + OH + CH2O */
    qf[237] = sc[3]*sc[37];
    qr[237] = sc[2]*sc[14]*sc[24];

    /*reaction 239: nC3H7 + CH3 <=> CH4 + C3H6 */
    qf[238] = sc[11]*sc[37];
    qr[238] = sc[12]*sc[36];

    /*reaction 240: C4H2 + H <=> iC4H3 */
    qf[239] = sc[0]*sc[40];
    qr[239] = sc[41];

    /*reaction 241: iC4H3 + H <=> C4H2 + H2 */
    qf[240] = sc[0]*sc[41];
    qr[240] = sc[4]*sc[40];

    /*reaction 242: C4H4 + OH <=> iC4H3 + H2O */
    qf[241] = sc[2]*sc[42];
    qr[241] = sc[5]*sc[41];

    /*reaction 243: C4H5-2 <=> iC4H5 */
    qf[242] = sc[44];
    qr[242] = sc[43];

    /*reaction 244: C4H6 + H <=> C2H4 + C2H3 */
    qf[243] = sc[0]*sc[45];
    qr[243] = sc[22]*sc[23];

    /*reaction 245: C4H6 + OH <=> iC4H5 + H2O */
    qf[244] = sc[2]*sc[45];
    qr[244] = sc[5]*sc[43];

    /*reaction 246: C4H612 <=> iC4H5 + H */
    qf[245] = sc[46];
    qr[245] = sc[0]*sc[43];

    /*reaction 247: C4H6-2 <=> H + C4H5-2 */
    qf[246] = sc[47];
    qr[246] = sc[0]*sc[44];

    /*reaction 248: C4H7 <=> C4H6 + H */
    qf[247] = sc[48];
    qr[247] = sc[0]*sc[45];

    /*reaction 249: C4H7 + O2 <=> C4H6 + HO2 */
    qf[248] = sc[7]*sc[48];
    qr[248] = sc[3]*sc[45];

    /*reaction 250: C4H7 + HO2 <=> CH2O + OH + aC3H5 */
    qf[249] = sc[3]*sc[48];
    qr[249] = sc[2]*sc[14]*sc[34];

    /*reaction 251: C4H81 + H <=> C2H4 + C2H5 */
    qf[250] = sc[0]*sc[49];
    qr[250] = sc[23]*sc[24];

    /*reaction 252: C4H81 + H <=> C3H6 + CH3 */
    qf[251] = sc[0]*sc[49];
    qr[251] = sc[11]*sc[36];

    /*reaction 253: C4H81 + H <=> C4H7 + H2 */
    qf[252] = sc[0]*sc[49];
    qr[252] = sc[4]*sc[48];

    /*reaction 254: C4H81 + O <=> nC3H7 + HCO */
    qf[253] = sc[1]*sc[49];
    qr[253] = sc[13]*sc[37];

    /*reaction 255: C2H4 + C2H5 <=> pC4H9 */
    qf[254] = sc[23]*sc[24];
    qr[254] = sc[50];

    /*reaction 256: pC4H9 + OH <=> C4H81 + H2O */
    qf[255] = sc[2]*sc[50];
    qr[255] = sc[5]*sc[49];

    /*reaction 257: pC4H9 + O2 <=> C4H81 + HO2 */
    qf[256] = sc[7]*sc[50];
    qr[256] = sc[3]*sc[49];

    /*reaction 258: pC4H9 + HO2 <=> nC3H7 + OH + CH2O */
    qf[257] = sc[3]*sc[50];
    qr[257] = sc[2]*sc[14]*sc[37];

    /*reaction 259: pC4H9 + CH3 <=> C4H81 + CH4 */
    qf[258] = sc[11]*sc[50];
    qr[258] = sc[12]*sc[49];

    /*reaction 260: NC12H26 => 3 C2H4 + 2 nC3H7 */
    qf[259] = sc[51];
    qr[259] = 0.0;

    /*reaction 261: NC12H26 => 2 C2H4 + 2 pC4H9 */
    qf[260] = sc[51];
    qr[260] = 0.0;

    /*reaction 262: NC12H26 + H => 4 C2H4 + pC4H9 + H2 */
    qf[261] = sc[0]*sc[51];
    qr[261] = 0.0;

    /*reaction 263: NC12H26 + H => C4H81 + 2 C2H4 + pC4H9 + H2 */
    qf[262] = sc[0]*sc[51];
    qr[262] = 0.0;

    /*reaction 264: NC12H26 + H => C3H6 + C6H12 + nC3H7 + H2 */
    qf[263] = sc[0]*sc[51];
    qr[263] = 0.0;

    /*reaction 265: NC12H26 + H => C5H10 + 2 C2H4 + nC3H7 + H2 */
    qf[264] = sc[0]*sc[51];
    qr[264] = 0.0;

    /*reaction 266: NC12H26 + H => C6H12 + C2H4 + pC4H9 + H2 */
    qf[265] = sc[0]*sc[51];
    qr[265] = 0.0;

    /*reaction 267: NC12H26 + CH3 => 4 C2H4 + pC4H9 + CH4 */
    qf[266] = sc[11]*sc[51];
    qr[266] = 0.0;

    /*reaction 268: NC12H26 + CH3 => C4H81 + 2 C2H4 + pC4H9 + CH4 */
    qf[267] = sc[11]*sc[51];
    qr[267] = 0.0;

    /*reaction 269: NC12H26 + CH3 => C3H6 + C6H12 + nC3H7 + CH4 */
    qf[268] = sc[11]*sc[51];
    qr[268] = 0.0;

    /*reaction 270: NC12H26 + CH3 => C5H10 + 2 C2H4 + nC3H7 + CH4 */
    qf[269] = sc[11]*sc[51];
    qr[269] = 0.0;

    /*reaction 271: NC12H26 + CH3 => C6H12 + C2H4 + pC4H9 + CH4 */
    qf[270] = sc[11]*sc[51];
    qr[270] = 0.0;

    /*reaction 272: NC12H26 + O => 4 C2H4 + pC4H9 + OH */
    qf[271] = sc[1]*sc[51];
    qr[271] = 0.0;

    /*reaction 273: NC12H26 + O => C4H81 + 2 C2H4 + pC4H9 + OH */
    qf[272] = sc[1]*sc[51];
    qr[272] = 0.0;

    /*reaction 274: NC12H26 + O => C3H6 + C6H12 + nC3H7 + OH */
    qf[273] = sc[1]*sc[51];
    qr[273] = 0.0;

    /*reaction 275: NC12H26 + O => C5H10 + 2 C2H4 + nC3H7 + OH */
    qf[274] = sc[1]*sc[51];
    qr[274] = 0.0;

    /*reaction 276: NC12H26 + O => C6H12 + C2H4 + pC4H9 + OH */
    qf[275] = sc[1]*sc[51];
    qr[275] = 0.0;

    /*reaction 277: NC12H26 + OH => 4 C2H4 + pC4H9 + H2O */
    qf[276] = sc[2]*sc[51];
    qr[276] = 0.0;

    /*reaction 278: NC12H26 + OH => C4H81 + 2 C2H4 + pC4H9 + H2O */
    qf[277] = sc[2]*sc[51];
    qr[277] = 0.0;

    /*reaction 279: NC12H26 + OH => C3H6 + C6H12 + nC3H7 + H2O */
    qf[278] = sc[2]*sc[51];
    qr[278] = 0.0;

    /*reaction 280: NC12H26 + OH => C5H10 + 2 C2H4 + nC3H7 + H2O */
    qf[279] = sc[2]*sc[51];
    qr[279] = 0.0;

    /*reaction 281: NC12H26 + OH => C6H12 + C2H4 + pC4H9 + H2O */
    qf[280] = sc[2]*sc[51];
    qr[280] = 0.0;

    /*reaction 282: C6H12 + H <=> C2H4 + pC4H9 */
    qf[281] = sc[0]*sc[52];
    qr[281] = sc[23]*sc[50];

    /*reaction 283: C6H12 + H <=> C6H11 + H2 */
    qf[282] = sc[0]*sc[52];
    qr[282] = sc[4]*sc[53];

    /*reaction 284: C5H10 + H <=> C2H4 + nC3H7 */
    qf[283] = sc[0]*sc[54];
    qr[283] = sc[23]*sc[37];

    /*reaction 285: C5H10 + H <=> C2H4 + aC3H5 + H2 */
    qf[284] = sc[0]*sc[54];
    qr[284] = sc[4]*sc[23]*sc[34];

    /*reaction 286: C6H11 + H <=> CH3 + C2H4 + aC3H5 */
    qf[285] = sc[0]*sc[53];
    qr[285] = sc[11]*sc[23]*sc[34];

    /*reaction 287: C6H11 + HO2 => CH2O + OH + aC3H5 + C2H4 */
    qf[286] = sc[3]*sc[53];
    qr[286] = 0.0;

    /*reaction 288: C6H12 + O <=> C2H4 + nC3H7 + HCO */
    qf[287] = sc[1]*sc[52];
    qr[287] = sc[13]*sc[23]*sc[37];

    /*reaction 289: C5H10 + O <=> pC4H9 + HCO */
    qf[288] = sc[1]*sc[54];
    qr[288] = sc[13]*sc[50];

    double T = tc[1];
    double Corr[289];
    for (int i = 0; i < 289; ++i) {
        Corr[i] = 1.0;
    }
    double k_f_save[289];
    double Kc_save[289];
    comp_k_f_d(tc,invT,k_f_save,Corr,sc);
    comp_Kc_d(tc,invT,Kc_save);


    for (int i=0; i<289; i++)
    {
        qf[i] *= Corr[i] * k_f_save[i];
        qr[i] *= Corr[i] * k_f_save[i] / Kc_save[i];
    }

    return;
}


/*compute the g/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
__device__ void gibbs_d(double * species, double *  tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +2.547365990000000e+04 * invT
            +2.946682853000000e+00
            -2.500000000000000e+00 * tc[0]
            -3.526664095000000e-13 * tc[1]
            +3.326532733333333e-16 * tc[2]
            -1.917346933333333e-19 * tc[3]
            +4.638661660000000e-23 * tc[4];
        /*species 1: O */
        species[1] =
            +2.912225920000000e+04 * invT
            +1.116333640000000e+00
            -3.168267100000000e+00 * tc[0]
            +1.639659420000000e-03 * tc[1]
            -1.107177326666667e-06 * tc[2]
            +5.106721866666666e-10 * tc[3]
            -1.056329855000000e-13 * tc[4];
        /*species 2: OH */
        species[2] =
            +3.381538120000000e+03 * invT
            +4.815738570000000e+00
            -4.125305610000000e+00 * tc[0]
            +1.612724695000000e-03 * tc[1]
            -1.087941151666667e-06 * tc[2]
            +4.832113691666666e-10 * tc[3]
            -1.031186895000000e-13 * tc[4];
        /*species 3: HO2 */
        species[3] =
            +2.948080400000000e+02 * invT
            +5.851355599999999e-01
            -4.301798010000000e+00 * tc[0]
            +2.374560255000000e-03 * tc[1]
            -3.526381516666666e-06 * tc[2]
            +2.023032450000000e-09 * tc[3]
            -4.646125620000001e-13 * tc[4];
        /*species 4: H2 */
        species[4] =
            -9.179351730000000e+02 * invT
            +1.661320882000000e+00
            -2.344331120000000e+00 * tc[0]
            -3.990260375000000e-03 * tc[1]
            +3.246358500000000e-06 * tc[2]
            -1.679767450000000e-09 * tc[3]
            +3.688058805000000e-13 * tc[4];
        /*species 5: H2O */
        species[5] =
            -3.029372670000000e+04 * invT
            +5.047672768000000e+00
            -4.198640560000000e+00 * tc[0]
            +1.018217050000000e-03 * tc[1]
            -1.086733685000000e-06 * tc[2]
            +4.573308850000000e-10 * tc[3]
            -8.859890850000000e-14 * tc[4];
        /*species 6: H2O2 */
        species[6] =
            -1.770258210000000e+04 * invT
            +8.410619499999998e-01
            -4.276112690000000e+00 * tc[0]
            +2.714112085000000e-04 * tc[1]
            -2.788928350000000e-06 * tc[2]
            +1.798090108333333e-09 * tc[3]
            -4.312271815000000e-13 * tc[4];
        /*species 7: O2 */
        species[7] =
            -1.063943560000000e+03 * invT
            +1.247806300000001e-01
            -3.782456360000000e+00 * tc[0]
            +1.498367080000000e-03 * tc[1]
            -1.641217001666667e-06 * tc[2]
            +8.067745908333334e-10 * tc[3]
            -1.621864185000000e-13 * tc[4];
        /*species 8: CH */
        species[8] =
            +7.079729340000000e+04 * invT
            +1.405805570000000e+00
            -3.489816650000000e+00 * tc[0]
            -1.619177705000000e-04 * tc[1]
            +2.814984416666667e-07 * tc[2]
            -2.635144391666666e-10 * tc[3]
            +7.030453350000001e-14 * tc[4];
        /*species 9: CH2 */
        species[9] =
            +4.600404010000000e+04 * invT
            +2.200146820000000e+00
            -3.762678670000000e+00 * tc[0]
            -4.844360715000000e-04 * tc[1]
            -4.658164016666667e-07 * tc[2]
            +3.209092941666667e-10 * tc[3]
            -8.437085950000000e-14 * tc[4];
        /*species 10: CH2* */
        species[10] =
            +5.049681630000000e+04 * invT
            +4.967723077000000e+00
            -4.198604110000000e+00 * tc[0]
            +1.183307095000000e-03 * tc[1]
            -1.372160366666667e-06 * tc[2]
            +5.573466508333334e-10 * tc[3]
            -9.715736850000000e-14 * tc[4];
        /*species 11: CH3 */
        species[11] =
            +1.644499880000000e+04 * invT
            +2.069026070000000e+00
            -3.673590400000000e+00 * tc[0]
            -1.005475875000000e-03 * tc[1]
            -9.550364266666668e-07 * tc[2]
            +5.725978541666666e-10 * tc[3]
            -1.271928670000000e-13 * tc[4];
        /*species 12: CH4 */
        species[12] =
            -1.024664760000000e+04 * invT
            +9.791179889999999e+00
            -5.149876130000000e+00 * tc[0]
            +6.835489400000000e-03 * tc[1]
            -8.196676650000000e-06 * tc[2]
            +4.039525216666667e-09 * tc[3]
            -8.334697800000000e-13 * tc[4];
        /*species 13: HCO */
        species[13] =
            +3.839564960000000e+03 * invT
            +8.268134100000002e-01
            -4.221185840000000e+00 * tc[0]
            +1.621962660000000e-03 * tc[1]
            -2.296657433333333e-06 * tc[2]
            +1.109534108333333e-09 * tc[3]
            -2.168844325000000e-13 * tc[4];
        /*species 14: CH2O */
        species[14] =
            -1.430895670000000e+04 * invT
            +4.190910250000000e+00
            -4.793723150000000e+00 * tc[0]
            +4.954166845000000e-03 * tc[1]
            -6.220333466666666e-06 * tc[2]
            +3.160710508333333e-09 * tc[3]
            -6.588632600000000e-13 * tc[4];
        /*species 15: CH3O */
        species[15] =
            +1.295697600000000e+03 * invT
            -2.860603620000000e+00
            -3.711805020000000e+00 * tc[0]
            +1.402316530000000e-03 * tc[1]
            -6.275849516666667e-06 * tc[2]
            +3.942267408333333e-09 * tc[3]
            -9.329421000000001e-13 * tc[4];
        /*species 16: CH2OH */
        species[16] =
            -3.500728900000000e+03 * invT
            +1.169208670000000e+00
            -4.478343670000000e+00 * tc[0]
            +6.753515500000000e-04 * tc[1]
            -4.641416333333334e-06 * tc[2]
            +3.040575500000000e-09 * tc[3]
            -7.395372499999999e-13 * tc[4];
        /*species 17: CH3OH */
        species[17] =
            -2.564276560000000e+04 * invT
            +7.219494050000001e+00
            -5.715395820000000e+00 * tc[0]
            +7.615456450000000e-03 * tc[1]
            -1.087401925000000e-05 * tc[2]
            +5.923390741666667e-09 * tc[3]
            -1.306763490000000e-12 * tc[4];
        /*species 18: CO */
        species[18] =
            -1.434408600000000e+04 * invT
            +7.112418999999992e-02
            -3.579533470000000e+00 * tc[0]
            +3.051768400000000e-04 * tc[1]
            -1.694690550000000e-07 * tc[2]
            -7.558382366666667e-11 * tc[3]
            +4.522122495000000e-14 * tc[4];
        /*species 19: CO2 */
        species[19] =
            -4.837196970000000e+04 * invT
            -7.544278700000000e+00
            -2.356773520000000e+00 * tc[0]
            -4.492298385000000e-03 * tc[1]
            +1.187260448333333e-06 * tc[2]
            -2.049325183333333e-10 * tc[3]
            +7.184977399999999e-15 * tc[4];
        /*species 20: C2H */
        species[20] =
            +6.683939320000001e+04 * invT
            -3.333307050000000e+00
            -2.889657330000000e+00 * tc[0]
            -6.704980550000000e-03 * tc[1]
            +4.746158350000000e-06 * tc[2]
            -2.456592041666667e-09 * tc[3]
            +5.466575550000000e-13 * tc[4];
        /*species 21: C2H2 */
        species[21] =
            +2.642898070000000e+04 * invT
            -1.313102400600000e+01
            -8.086810940000000e-01 * tc[0]
            -1.168078145000000e-02 * tc[1]
            +5.919530250000000e-06 * tc[2]
            -2.334603641666667e-09 * tc[3]
            +4.250364870000000e-13 * tc[4];
        /*species 22: C2H3 */
        species[22] =
            +3.485984680000000e+04 * invT
            -5.298073800000000e+00
            -3.212466450000000e+00 * tc[0]
            -7.573958100000000e-04 * tc[1]
            -4.320156866666666e-06 * tc[2]
            +2.980482058333333e-09 * tc[3]
            -7.357543650000000e-13 * tc[4];
        /*species 23: C2H4 */
        species[23] =
            +5.089775930000000e+03 * invT
            -1.381294799999999e-01
            -3.959201480000000e+00 * tc[0]
            +3.785261235000000e-03 * tc[1]
            -9.516504866666667e-06 * tc[2]
            +5.763239608333333e-09 * tc[3]
            -1.349421865000000e-12 * tc[4];
        /*species 24: C2H5 */
        species[24] =
            +1.284162650000000e+04 * invT
            -4.007435600000004e-01
            -4.306465680000000e+00 * tc[0]
            +2.093294460000000e-03 * tc[1]
            -8.285713450000000e-06 * tc[2]
            +4.992721716666666e-09 * tc[3]
            -1.152545020000000e-12 * tc[4];
        /*species 25: C2H6 */
        species[25] =
            -1.152220550000000e+04 * invT
            +1.624601760000000e+00
            -4.291424920000000e+00 * tc[0]
            +2.750771350000000e-03 * tc[1]
            -9.990638133333334e-06 * tc[2]
            +5.903885708333334e-09 * tc[3]
            -1.343428855000000e-12 * tc[4];
        /*species 26: HCCO */
        species[26] =
            +2.005944900000000e+04 * invT
            -1.023869560000000e+01
            -2.251721400000000e+00 * tc[0]
            -8.827510500000000e-03 * tc[1]
            +3.954850166666666e-06 * tc[2]
            -1.439646583333334e-09 * tc[3]
            +2.533240550000000e-13 * tc[4];
        /*species 27: CH2CO */
        species[27] =
            -7.270000000000000e+03 * invT
            -1.007981170000000e+01
            -2.135836300000000e+00 * tc[0]
            -9.059436050000000e-03 * tc[1]
            +2.899124566666666e-06 * tc[2]
            -7.786646400000000e-10 * tc[3]
            +1.007288075000000e-13 * tc[4];
        /*species 28: CH3CO */
        species[28] =
            -2.657452900000000e+03 * invT
            -3.183402300000000e+00
            -4.163425700000000e+00 * tc[0]
            +1.163080500000000e-04 * tc[1]
            -5.711303333333333e-06 * tc[2]
            +3.675435583333333e-09 * tc[3]
            -8.637806000000000e-13 * tc[4];
        /*species 29: CH2CHO */
        species[29] =
            +6.200000000000000e+01 * invT
            -6.162391100000001e+00
            -3.409062400000000e+00 * tc[0]
            -5.369287000000000e-03 * tc[1]
            -3.152487500000000e-07 * tc[2]
            +5.965485916666667e-10 * tc[3]
            -1.433692550000000e-13 * tc[4];
        /*species 30: CH3CHO */
        species[30] =
            -2.157287800000000e+04 * invT
            +6.264436000000000e-01
            -4.729459500000000e+00 * tc[0]
            +1.596642900000000e-03 * tc[1]
            -7.922486833333334e-06 * tc[2]
            +4.788217583333333e-09 * tc[3]
            -1.096555600000000e-12 * tc[4];
        /*species 31: C3H3 */
        species[31] =
            +4.010577830000000e+04 * invT
            -1.385478313000000e+01
            -1.351109270000000e+00 * tc[0]
            -1.637056115000000e-02 * tc[1]
            +7.897118916666667e-06 * tc[2]
            -3.135915066666667e-09 * tc[3]
            +5.927046150000000e-13 * tc[4];
        /*species 32: pC3H4 */
        species[32] =
            +2.080237400000000e+04 * invT
            -7.196548200000001e+00
            -2.680386900000000e+00 * tc[0]
            -7.899825500000001e-03 * tc[1]
            -4.178432666666667e-07 * tc[2]
            +1.138135250000000e-09 * tc[3]
            -3.307714250000000e-13 * tc[4];
        /*species 33: aC3H4 */
        species[33] =
            +2.154156700000000e+04 * invT
            -7.613094500000000e+00
            -2.613044500000000e+00 * tc[0]
            -6.061287500000000e-03 * tc[1]
            -3.089980000000000e-06 * tc[2]
            +2.877095750000000e-09 * tc[3]
            -7.667539500000001e-13 * tc[4];
        /*species 34: aC3H5 */
        species[34] =
            +1.924562900000000e+04 * invT
            -1.581003050000000e+01
            -1.363183500000000e+00 * tc[0]
            -9.906910499999999e-03 * tc[1]
            -2.082843333333333e-06 * tc[2]
            +2.779629583333333e-09 * tc[3]
            -7.923285500000000e-13 * tc[4];
        /*species 35: CH3CCH2 */
        species[35] =
            +2.904049800000000e+04 * invT
            -1.483595710000000e+01
            -1.732920900000000e+00 * tc[0]
            -1.119731000000000e-02 * tc[1]
            +8.581768500000001e-07 * tc[2]
            +5.633038833333333e-10 * tc[3]
            -1.912660550000000e-13 * tc[4];
        /*species 36: C3H6 */
        species[36] =
            +1.074826000000000e+03 * invT
            -1.465203300000000e+01
            -1.493307000000000e+00 * tc[0]
            -1.046259000000000e-02 * tc[1]
            -7.477990000000000e-07 * tc[2]
            +1.390760000000000e-09 * tc[3]
            -3.579073000000000e-13 * tc[4];
        /*species 37: nC3H7 */
        species[37] =
            +1.031234600000000e+04 * invT
            -2.008691670000000e+01
            -1.049117300000000e+00 * tc[0]
            -1.300448650000000e-02 * tc[1]
            -3.923752666666667e-07 * tc[2]
            +1.632927666666667e-09 * tc[3]
            -4.686010350000000e-13 * tc[4];
        /*species 38: iC3H7 */
        species[38] =
            +9.422372400000000e+03 * invT
            -1.867139710000000e+01
            -1.444919900000000e+00 * tc[0]
            -1.049955600000000e-02 * tc[1]
            -1.283937033333333e-06 * tc[2]
            +1.539687750000000e-09 * tc[3]
            -3.564148100000000e-13 * tc[4];
        /*species 39: C2H3CHO */
        species[39] =
            -9.335734399999999e+03 * invT
            -1.822672720000000e+01
            -1.271349800000000e+00 * tc[0]
            -1.311552700000000e-02 * tc[1]
            +1.548538416666667e-06 * tc[2]
            +3.986439333333333e-10 * tc[3]
            -1.674027150000000e-13 * tc[4];
        /*species 40: C4H2 */
        species[40] =
            +5.418521100000000e+04 * invT
            -1.381219320000000e+01
            -1.054397800000000e+00 * tc[0]
            -2.081348000000000e-02 * tc[1]
            +1.097863066666667e-05 * tc[2]
            -4.438089583333333e-09 * tc[3]
            +8.341581000000001e-13 * tc[4];
        /*species 41: iC4H3 */
        species[41] =
            +5.883712100000000e+04 * invT
            -3.841576300000000e+00
            -3.722148200000000e+00 * tc[0]
            -1.297877150000000e-02 * tc[1]
            +4.392723833333334e-06 * tc[2]
            -1.292410000000000e-09 * tc[3]
            +1.902028250000000e-13 * tc[4];
        /*species 42: C4H4 */
        species[42] =
            +3.335949200000000e+04 * invT
            -2.006931052000000e+01
            -5.885704800000000e-01 * tc[0]
            -1.827334250000000e-02 * tc[1]
            +5.684494666666667e-06 * tc[2]
            -1.387718250000000e-09 * tc[3]
            +1.503231150000000e-13 * tc[4];
        /*species 43: iC4H5 */
        species[43] =
            +3.638337100000000e+04 * invT
            -2.357937595000000e+01
            -1.130810500000000e-01 * tc[0]
            -2.047530750000000e-02 * tc[1]
            +5.902263500000000e-06 * tc[2]
            -1.294247416666666e-09 * tc[3]
            +1.167756100000000e-13 * tc[4];
        /*species 44: C4H5-2 */
        species[44] =
            +3.550331600000000e+04 * invT
            -9.066423000000000e+00
            -2.969628000000000e+00 * tc[0]
            -1.222112250000000e-02 * tc[1]
            +1.520857066666667e-06 * tc[2]
            +3.538905916666667e-19 * tc[3]
            -8.152364000000000e-23 * tc[4];
        /*species 45: C4H6 */
        species[45] =
            +1.180227000000000e+04 * invT
            -2.297715135000000e+01
            -1.128446500000000e-01 * tc[0]
            -1.718451100000000e-02 * tc[1]
            +1.851232000000000e-06 * tc[2]
            +7.675555000000000e-10 * tc[3]
            -3.103258950000000e-13 * tc[4];
        /*species 46: C4H612 */
        species[46] =
            +1.811799000000000e+04 * invT
            -1.872719300000000e+01
            -1.023467000000000e+00 * tc[0]
            -1.747959500000000e-02 * tc[1]
            +3.668175000000000e-06 * tc[2]
            -5.785226666666666e-10 * tc[3]
            +3.939593500000000e-14 * tc[4];
        /*species 47: C4H6-2 */
        species[47] =
            +1.571090200000000e+04 * invT
            -1.139209220000000e+01
            -2.137333800000000e+00 * tc[0]
            -1.324311450000000e-02 * tc[1]
            +1.509478516666667e-06 * tc[2]
            +4.615533083333334e-20 * tc[3]
            -1.064094200000000e-23 * tc[4];
        /*species 48: C4H7 */
        species[48] =
            +2.265332800000000e+04 * invT
            -2.269338368000000e+01
            -7.444943200000000e-01 * tc[0]
            -1.983942850000000e-02 * tc[1]
            +3.816347666666667e-06 * tc[2]
            -1.779414416666667e-10 * tc[3]
            -1.154818750000000e-13 * tc[4];
        /*species 49: C4H81 */
        species[49] =
            -1.790400400000000e+03 * invT
            -1.988133100000000e+01
            -1.181138000000000e+00 * tc[0]
            -1.542669000000000e-02 * tc[1]
            -8.477541166666667e-07 * tc[2]
            +2.054574000000000e-09 * tc[3]
            -5.555096499999999e-13 * tc[4];
        /*species 50: pC4H9 */
        species[50] =
            +7.322104000000000e+03 * invT
            -2.096056380000000e+01
            -1.208704200000000e+00 * tc[0]
            -1.914874850000000e-02 * tc[1]
            +1.211008483333333e-06 * tc[2]
            +1.285712250000000e-09 * tc[3]
            -4.342971750000000e-13 * tc[4];
        /*species 55: N2 */
        species[55] =
            -1.020899900000000e+03 * invT
            -6.516950000000001e-01
            -3.298677000000000e+00 * tc[0]
            -7.041202000000000e-04 * tc[1]
            +6.605369999999999e-07 * tc[2]
            -4.701262500000001e-10 * tc[3]
            +1.222427000000000e-13 * tc[4];
    } else {
        /*species 0: H */
        species[0] =
            +2.547365990000000e+04 * invT
            +2.946682924000000e+00
            -2.500000010000000e+00 * tc[0]
            +1.154214865000000e-11 * tc[1]
            -2.692699133333334e-15 * tc[2]
            +3.945960291666667e-19 * tc[3]
            -2.490986785000000e-23 * tc[4];
        /*species 1: O */
        species[1] =
            +2.921757910000000e+04 * invT
            -2.214917859999999e+00
            -2.569420780000000e+00 * tc[0]
            +4.298705685000000e-05 * tc[1]
            -6.991409816666667e-09 * tc[2]
            +8.348149916666666e-13 * tc[3]
            -6.141684549999999e-17 * tc[4];
        /*species 2: OH */
        species[2] =
            +3.718857740000000e+03 * invT
            -2.836911870000000e+00
            -2.864728860000000e+00 * tc[0]
            -5.282522400000000e-04 * tc[1]
            +4.318045966666667e-08 * tc[2]
            -2.543488950000000e-12 * tc[3]
            +6.659793800000000e-17 * tc[4];
        /*species 3: HO2 */
        species[3] =
            +1.118567130000000e+02 * invT
            +2.321087500000001e-01
            -4.017210900000000e+00 * tc[0]
            -1.119910065000000e-03 * tc[1]
            +1.056096916666667e-07 * tc[2]
            -9.520530833333334e-12 * tc[3]
            +5.395426750000000e-16 * tc[4];
        /*species 4: H2 */
        species[4] =
            -9.501589220000000e+02 * invT
            +6.542302510000000e+00
            -3.337279200000000e+00 * tc[0]
            +2.470123655000000e-05 * tc[1]
            -8.324279633333333e-08 * tc[2]
            +1.496386616666667e-11 * tc[3]
            -1.001276880000000e-15 * tc[4];
        /*species 5: H2O */
        species[5] =
            -3.000429710000000e+04 * invT
            -1.932777610000000e+00
            -3.033992490000000e+00 * tc[0]
            -1.088459020000000e-03 * tc[1]
            +2.734541966666666e-08 * tc[2]
            +8.086832250000000e-12 * tc[3]
            -8.410049600000000e-16 * tc[4];
        /*species 6: H2O2 */
        species[6] =
            -1.786178770000000e+04 * invT
            +1.248846229999999e+00
            -4.165002850000000e+00 * tc[0]
            -2.454158470000000e-03 * tc[1]
            +3.168987083333333e-07 * tc[2]
            -3.093216550000000e-11 * tc[3]
            +1.439541525000000e-15 * tc[4];
        /*species 7: O2 */
        species[7] =
            -1.088457720000000e+03 * invT
            -2.170693450000000e+00
            -3.282537840000000e+00 * tc[0]
            -7.415437700000000e-04 * tc[1]
            +1.263277781666667e-07 * tc[2]
            -1.745587958333333e-11 * tc[3]
            +1.083588970000000e-15 * tc[4];
        /*species 8: CH */
        species[8] =
            +7.101243640000001e+04 * invT
            -2.606515260000000e+00
            -2.878464730000000e+00 * tc[0]
            -4.854568405000000e-04 * tc[1]
            -2.407427583333333e-08 * tc[2]
            +1.089065408333333e-11 * tc[3]
            -8.803969149999999e-16 * tc[4];
        /*species 9: CH2 */
        species[9] =
            +4.626360400000000e+04 * invT
            -3.297092110000000e+00
            -2.874101130000000e+00 * tc[0]
            -1.828196460000000e-03 * tc[1]
            +2.348243283333333e-07 * tc[2]
            -2.168162908333333e-11 * tc[3]
            +9.386378350000000e-16 * tc[4];
        /*species 10: CH2* */
        species[10] =
            +5.092599970000000e+04 * invT
            -6.334463270000000e+00
            -2.292038420000000e+00 * tc[0]
            -2.327943185000000e-03 * tc[1]
            +3.353199116666667e-07 * tc[2]
            -3.482550000000000e-11 * tc[3]
            +1.698581825000000e-15 * tc[4];
        /*species 11: CH3 */
        species[11] =
            +1.677558430000000e+04 * invT
            -6.194354070000000e+00
            -2.285717720000000e+00 * tc[0]
            -3.619950185000000e-03 * tc[1]
            +4.978572466666667e-07 * tc[2]
            -4.964038700000000e-11 * tc[3]
            +2.335771970000000e-15 * tc[4];
        /*species 12: CH4 */
        species[12] =
            -9.468344590000001e+03 * invT
            -1.836246650500000e+01
            -7.485149500000000e-02 * tc[0]
            -6.695473350000000e-03 * tc[1]
            +9.554763483333333e-07 * tc[2]
            -1.019104458333333e-10 * tc[3]
            +5.090761500000000e-15 * tc[4];
        /*species 13: HCO */
        species[13] =
            +4.011918150000000e+03 * invT
            -7.026170540000000e+00
            -2.772174380000000e+00 * tc[0]
            -2.478477630000000e-03 * tc[1]
            +4.140760216666667e-07 * tc[2]
            -4.909681483333334e-11 * tc[3]
            +2.667543555000000e-15 * tc[4];
        /*species 14: CH2O */
        species[14] =
            -1.399583230000000e+04 * invT
            -1.189563292000000e+01
            -1.760690080000000e+00 * tc[0]
            -4.600000410000000e-03 * tc[1]
            +7.370980216666666e-07 * tc[2]
            -8.386767666666666e-11 * tc[3]
            +4.419278200000001e-15 * tc[4];
        /*species 15: CH3O */
        species[15] =
            +3.781119400000000e+02 * invT
            +6.724592660000000e+00
            -4.757792380000000e+00 * tc[0]
            -3.720712370000000e-03 * tc[1]
            +4.495086266666666e-07 * tc[2]
            -3.650754200000000e-11 * tc[3]
            +1.317685490000000e-15 * tc[4];
        /*species 16: CH2OH */
        species[16] =
            -4.034096400000000e+03 * invT
            +6.940058629999999e+00
            -5.093143700000000e+00 * tc[0]
            -2.973806300000000e-03 * tc[1]
            +3.441624333333333e-07 * tc[2]
            -2.691734775000000e-11 * tc[3]
            +9.406295100000001e-16 * tc[4];
        /*species 17: CH3OH */
        species[17] =
            -2.537487470000000e+04 * invT
            -1.271265439000000e+01
            -1.789707910000000e+00 * tc[0]
            -7.046914600000000e-03 * tc[1]
            +1.060834725000000e-06 * tc[2]
            -1.151425708333333e-10 * tc[3]
            +5.853011000000000e-15 * tc[4];
        /*species 18: CO */
        species[18] =
            -1.415187240000000e+04 * invT
            -5.103502110000000e+00
            -2.715185610000000e+00 * tc[0]
            -1.031263715000000e-03 * tc[1]
            +1.664709618333334e-07 * tc[2]
            -1.917108400000000e-11 * tc[3]
            +1.018238580000000e-15 * tc[4];
        /*species 19: CO2 */
        species[19] =
            -4.875916600000000e+04 * invT
            +1.585822230000000e+00
            -3.857460290000000e+00 * tc[0]
            -2.207185130000000e-03 * tc[1]
            +3.691356733333334e-07 * tc[2]
            -4.362418233333334e-11 * tc[3]
            +2.360420820000000e-15 * tc[4];
        /*species 20: C2H */
        species[20] =
            +6.712106500000000e+04 * invT
            -3.468088230000000e+00
            -3.167806520000000e+00 * tc[0]
            -2.376109510000000e-03 * tc[1]
            +3.063117950000000e-07 * tc[2]
            -2.534918766666666e-11 * tc[3]
            +8.861638500000000e-16 * tc[4];
        /*species 21: C2H2 */
        species[21] =
            +2.593599920000000e+04 * invT
            +5.377850850000001e+00
            -4.147569640000000e+00 * tc[0]
            -2.980833320000000e-03 * tc[1]
            +3.954914200000000e-07 * tc[2]
            -3.895101425000000e-11 * tc[3]
            +1.806176065000000e-15 * tc[4];
        /*species 22: C2H3 */
        species[22] =
            +3.461287390000000e+04 * invT
            -4.770599780000000e+00
            -3.016724000000000e+00 * tc[0]
            -5.165114600000000e-03 * tc[1]
            +7.801372483333333e-07 * tc[2]
            -8.480274000000000e-11 * tc[3]
            +4.313035205000000e-15 * tc[4];
        /*species 23: C2H4 */
        species[23] =
            +4.939886140000000e+03 * invT
            -8.269258140000002e+00
            -2.036111160000000e+00 * tc[0]
            -7.322707550000000e-03 * tc[1]
            +1.118463191666667e-06 * tc[2]
            -1.226857691666667e-10 * tc[3]
            +6.285303050000000e-15 * tc[4];
        /*species 24: C2H5 */
        species[24] =
            +1.285752000000000e+04 * invT
            -1.150777788000000e+01
            -1.954656420000000e+00 * tc[0]
            -8.698636100000001e-03 * tc[1]
            +1.330344446666667e-06 * tc[2]
            -1.460147408333333e-10 * tc[3]
            +7.482078800000000e-15 * tc[4];
        /*species 25: C2H6 */
        species[25] =
            -1.142639320000000e+04 * invT
            -1.404372920000000e+01
            -1.071881500000000e+00 * tc[0]
            -1.084263385000000e-02 * tc[1]
            +1.670934450000000e-06 * tc[2]
            -1.845100008333333e-10 * tc[3]
            +9.500144500000000e-15 * tc[4];
        /*species 26: HCCO */
        species[26] =
            +1.932721500000000e+04 * invT
            +9.558465300000000e+00
            -5.628205800000000e+00 * tc[0]
            -2.042670050000000e-03 * tc[1]
            +2.655757833333333e-07 * tc[2]
            -2.385504333333333e-11 * tc[3]
            +9.703915999999999e-16 * tc[4];
        /*species 27: CH2CO */
        species[27] =
            -7.778500000000000e+03 * invT
            +3.879050115000000e+00
            -4.511297320000000e+00 * tc[0]
            -4.501798725000000e-03 * tc[1]
            +6.948993916666666e-07 * tc[2]
            -7.694549016666667e-11 * tc[3]
            +3.974191005000000e-15 * tc[4];
        /*species 28: CH3CO */
        species[28] =
            -3.787307500000000e+03 * invT
            +1.095844820000000e+01
            -5.944773100000000e+00 * tc[0]
            -3.933360250000000e-03 * tc[1]
            +4.810980333333333e-07 * tc[2]
            -3.939239583333334e-11 * tc[3]
            +1.429993050000000e-15 * tc[4];
        /*species 29: CH2CHO */
        species[29] =
            -9.695000000000000e+02 * invT
            +1.100775780000000e+01
            -5.975669900000000e+00 * tc[0]
            -4.065295700000000e-03 * tc[1]
            +4.572707500000000e-07 * tc[2]
            -3.391920083333333e-11 * tc[3]
            +1.088008550000000e-15 * tc[4];
        /*species 30: CH3CHO */
        species[30] =
            -2.259312200000000e+04 * invT
            +8.884902499999999e+00
            -5.404110800000000e+00 * tc[0]
            -5.861529500000000e-03 * tc[1]
            +7.043856166666666e-07 * tc[2]
            -5.697704250000000e-11 * tc[3]
            +2.049243150000000e-15 * tc[4];
        /*species 31: C3H3 */
        species[31] =
            +3.890874270000000e+04 * invT
            +1.972706240000000e+01
            -7.142218800000000e+00 * tc[0]
            -3.809510025000000e-03 * tc[1]
            +4.457665833333333e-07 * tc[2]
            -3.540956675000000e-11 * tc[3]
            +1.257377075000000e-15 * tc[4];
        /*species 32: pC3H4 */
        species[32] =
            +1.962094200000000e+04 * invT
            +1.462961850000000e+01
            -6.025240000000000e+00 * tc[0]
            -5.668271000000000e-03 * tc[1]
            +6.703898500000001e-07 * tc[2]
            -5.364671916666667e-11 * tc[3]
            +1.914981750000000e-15 * tc[4];
        /*species 33: aC3H4 */
        species[33] =
            +2.011749500000000e+04 * invT
            +1.731263820000000e+01
            -6.316872200000000e+00 * tc[0]
            -5.566864000000000e-03 * tc[1]
            +6.604896333333332e-07 * tc[2]
            -5.297019833333334e-11 * tc[3]
            +1.893777000000000e-15 * tc[4];
        /*species 34: aC3H5 */
        species[34] =
            +1.748244900000000e+04 * invT
            +1.774383770000000e+01
            -6.500787700000000e+00 * tc[0]
            -7.162365500000000e-03 * tc[1]
            +9.463605333333332e-07 * tc[2]
            -9.234000833333333e-11 * tc[3]
            +4.518194349999999e-15 * tc[4];
        /*species 35: CH3CCH2 */
        species[35] =
            +2.784302700000000e+04 * invT
            +8.778271200000001e+00
            -5.425552800000000e+00 * tc[0]
            -7.755536000000000e-03 * tc[1]
            +9.446391666666667e-07 * tc[2]
            -6.602032333333333e-11 * tc[3]
            +8.439017000000001e-16 * tc[4];
        /*species 36: C3H6 */
        species[36] =
            -9.235703000000000e+02 * invT
            +2.004560700000000e+01
            -6.732257000000000e+00 * tc[0]
            -7.454170000000000e-03 * tc[1]
            +8.249831666666666e-07 * tc[2]
            -6.010018333333334e-11 * tc[3]
            +1.883102000000000e-15 * tc[4];
        /*species 37: nC3H7 */
        species[37] =
            +7.976223600000000e+03 * invT
            +2.322504490000000e+01
            -7.709747900000000e+00 * tc[0]
            -8.015742500000001e-03 * tc[1]
            +8.786706333333332e-07 * tc[2]
            -6.324029333333334e-11 * tc[3]
            +1.943135950000000e-15 * tc[4];
        /*species 38: iC3H7 */
        species[38] =
            +7.322719300000000e+03 * invT
            +1.560229560000000e+01
            -6.519274100000000e+00 * tc[0]
            -8.610052000000000e-03 * tc[1]
            +9.560702833333334e-07 * tc[2]
            -7.010894333333333e-11 * tc[3]
            +2.228295650000000e-15 * tc[4];
        /*species 39: C2H3CHO */
        species[39] =
            -1.078405400000000e+04 * invT
            +1.066998720000000e+01
            -5.811186800000000e+00 * tc[0]
            -8.557128000000001e-03 * tc[1]
            +1.247236016666667e-06 * tc[2]
            -1.187687416666667e-10 * tc[3]
            +4.587342050000000e-15 * tc[4];
        /*species 40: C4H2 */
        species[40] =
            +5.258803900000000e+04 * invT
            +3.286909280000000e+01
            -9.157632800000000e+00 * tc[0]
            -2.771525900000000e-03 * tc[1]
            +2.265267333333333e-07 * tc[2]
            -1.565006250000000e-12 * tc[3]
            -1.159476800000000e-15 * tc[4];
        /*species 41: iC4H3 */
        species[41] =
            +5.795436300000000e+04 * invT
            +1.941033080000000e+01
            -7.653854800000000e+00 * tc[0]
            -5.602027500000000e-03 * tc[1]
            +7.733557000000000e-07 * tc[2]
            -7.232219916666667e-11 * tc[3]
            +2.871528100000000e-15 * tc[4];
        /*species 42: C4H4 */
        species[42] =
            +3.176601600000000e+04 * invT
            +1.988348110000000e+01
            -7.253960100000000e+00 * tc[0]
            -6.957047000000000e-03 * tc[1]
            +8.822035666666667e-07 * tc[2]
            -6.956704166666667e-11 * tc[3]
            +1.759894100000000e-15 * tc[4];
        /*species 43: iC4H5 */
        species[43] =
            +3.472509800000000e+04 * invT
            +1.761392390000000e+01
            -6.964602900000000e+00 * tc[0]
            -9.137166500000000e-03 * tc[1]
            +1.302228916666667e-06 * tc[2]
            -1.274346166666667e-10 * tc[3]
            +5.460246500000000e-15 * tc[4];
        /*species 44: C4H5-2 */
        species[44] =
            +3.325909500000000e+04 * invT
            +5.990766800000000e+01
            -1.453817100000000e+01 * tc[0]
            +4.283852800000000e-03 * tc[1]
            -3.926587333333333e-06 * tc[2]
            +1.139698250000000e-09 * tc[3]
            -1.221846350000000e-13 * tc[4];
        /*species 45: C4H6 */
        species[45] =
            +9.133851600000000e+03 * invT
            +3.219548440000000e+01
            -8.867313400000000e+00 * tc[0]
            -7.459335000000000e-03 * tc[1]
            +5.258119333333333e-07 * tc[2]
            +3.486777500000000e-11 * tc[3]
            -7.880629000000001e-15 * tc[4];
        /*species 46: C4H612 */
        species[46] =
            +1.267342000000000e+04 * invT
            +8.764219000000000e+01
            -1.781557000000000e+01 * tc[0]
            +2.128751000000000e-03 * tc[1]
            -1.751975000000000e-06 * tc[2]
            +3.728203333333333e-10 * tc[3]
            -2.924069000000000e-14 * tc[4];
        /*species 47: C4H6-2 */
        species[47] =
            +1.433506800000000e+04 * invT
            +3.001957530000000e+01
            -9.033813300000000e+00 * tc[0]
            -4.106225500000000e-03 * tc[1]
            -1.195899200000000e-06 * tc[2]
            +4.902861166666666e-10 * tc[3]
            -5.171957500000000e-14 * tc[4];
        /*species 48: C4H7 */
        species[48] =
            +2.095500800000000e+04 * invT
            +1.590279150000000e+01
            -7.013483500000000e+00 * tc[0]
            -1.131727900000000e-02 * tc[1]
            +1.542424500000000e-06 * tc[2]
            -1.400660583333333e-10 * tc[3]
            +5.204308500000000e-15 * tc[4];
        /*species 49: C4H81 */
        species[49] =
            -2.139723100000000e+03 * invT
            -1.348961690000000e+01
            -2.053584100000000e+00 * tc[0]
            -1.717525350000000e-02 * tc[1]
            +2.647199500000000e-06 * tc[2]
            -2.757471833333334e-10 * tc[3]
            +1.268052250000000e-14 * tc[4];
        /*species 50: pC4H9 */
        species[50] =
            +4.964405800000000e+03 * invT
            +2.657398650000000e+01
            -8.682239500000000e+00 * tc[0]
            -1.184553550000000e-02 * tc[1]
            +1.265814416666667e-06 * tc[2]
            -5.535594666666666e-11 * tc[3]
            -2.742256800000000e-15 * tc[4];
        /*species 55: N2 */
        species[55] =
            -9.227977000000000e+02 * invT
            -3.053888000000000e+00
            -2.926640000000000e+00 * tc[0]
            -7.439884000000000e-04 * tc[1]
            +9.474600000000001e-08 * tc[2]
            -8.414198333333333e-12 * tc[3]
            +3.376675500000000e-16 * tc[4];
    }

    /*species with midpoint at T=1392 kelvin */
    if (T < 1392) {
        /*species 52: C6H12 */
        species[52] =
            -7.343686170000000e+03 * invT
            -3.666482115000000e+01
            +1.352752050000000e+00 * tc[0]
            -3.493277130000000e-02 * tc[1]
            +7.656800366666667e-06 * tc[2]
            -1.308061191666667e-09 * tc[3]
            +1.106480875000000e-13 * tc[4];
        /*species 54: C5H10 */
        species[54] =
            -4.465466660000000e+03 * invT
            -3.333621381000000e+01
            +1.062234810000000e+00 * tc[0]
            -2.871091470000000e-02 * tc[1]
            +6.241448166666667e-06 * tc[2]
            -1.061374908333333e-09 * tc[3]
            +8.980489449999999e-14 * tc[4];
    } else {
        /*species 52: C6H12 */
        species[52] =
            -1.420628600000000e+04 * invT
            +8.621563800000001e+01
            -1.783375290000000e+01 * tc[0]
            -1.336888290000000e-02 * tc[1]
            +1.516727955000000e-06 * tc[2]
            -1.173498066666667e-10 * tc[3]
            +4.075621220000000e-15 * tc[4];
        /*species 54: C5H10 */
        species[54] =
            -1.008982050000000e+04 * invT
            +6.695354750000000e+01
            -1.458515390000000e+01 * tc[0]
            -1.120362355000000e-02 * tc[1]
            +1.272246708333333e-06 * tc[2]
            -9.849080500000001e-11 * tc[3]
            +3.421925695000000e-15 * tc[4];
    }

    /*species with midpoint at T=1389 kelvin */
    if (T < 1389) {
        /*species 53: C6H11 */
        species[53] =
            +9.663166950000001e+03 * invT
            -3.670371523999999e+01
            +1.555449440000000e+00 * tc[0]
            -3.384328010000000e-02 * tc[1]
            +7.450810583333334e-06 * tc[2]
            -1.268638583333333e-09 * tc[3]
            +1.071731885000000e-13 * tc[4];
    } else {
        /*species 53: C6H11 */
        species[53] =
            +2.680171740000000e+03 * invT
            +8.708080579999999e+01
            -1.773365500000000e+01 * tc[0]
            -1.244673875000000e-02 * tc[1]
            +1.433319083333333e-06 * tc[2]
            -1.120106900000000e-10 * tc[3]
            +3.917378330000000e-15 * tc[4];
    }

    /*species with midpoint at T=1391 kelvin */
    if (T < 1391) {
        /*species 51: NC12H26 */
        species[51] =
            -4.006542530000000e+04 * invT
            -5.272127854000000e+01
            +2.621815940000000e+00 * tc[0]
            -7.361885550000000e-02 * tc[1]
            +1.573283785000000e-05 * tc[2]
            -2.562010566666667e-09 * tc[3]
            +2.018011150000000e-13 * tc[4];
    } else {
        /*species 51: NC12H26 */
        species[51] =
            -5.488434650000000e+04 * invT
            +2.111804257000000e+02
            -3.850950370000000e+01 * tc[0]
            -2.817750240000000e-02 * tc[1]
            +3.191553333333333e-06 * tc[2]
            -2.466873850000000e-10 * tc[3]
            +8.562207500000000e-15 * tc[4];
    }
    return;
}


/*compute Cv/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
__device__ void cv_R_d(double * species, double *  tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +1.50000000e+00
            +7.05332819e-13 * tc[1]
            -1.99591964e-15 * tc[2]
            +2.30081632e-18 * tc[3]
            -9.27732332e-22 * tc[4];
        /*species 1: O */
        species[1] =
            +2.16826710e+00
            -3.27931884e-03 * tc[1]
            +6.64306396e-06 * tc[2]
            -6.12806624e-09 * tc[3]
            +2.11265971e-12 * tc[4];
        /*species 2: OH */
        species[2] =
            +3.12530561e+00
            -3.22544939e-03 * tc[1]
            +6.52764691e-06 * tc[2]
            -5.79853643e-09 * tc[3]
            +2.06237379e-12 * tc[4];
        /*species 3: HO2 */
        species[3] =
            +3.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
        /*species 4: H2 */
        species[4] =
            +1.34433112e+00
            +7.98052075e-03 * tc[1]
            -1.94781510e-05 * tc[2]
            +2.01572094e-08 * tc[3]
            -7.37611761e-12 * tc[4];
        /*species 5: H2O */
        species[5] =
            +3.19864056e+00
            -2.03643410e-03 * tc[1]
            +6.52040211e-06 * tc[2]
            -5.48797062e-09 * tc[3]
            +1.77197817e-12 * tc[4];
        /*species 6: H2O2 */
        species[6] =
            +3.27611269e+00
            -5.42822417e-04 * tc[1]
            +1.67335701e-05 * tc[2]
            -2.15770813e-08 * tc[3]
            +8.62454363e-12 * tc[4];
        /*species 7: O2 */
        species[7] =
            +2.78245636e+00
            -2.99673416e-03 * tc[1]
            +9.84730201e-06 * tc[2]
            -9.68129509e-09 * tc[3]
            +3.24372837e-12 * tc[4];
        /*species 8: CH */
        species[8] =
            +2.48981665e+00
            +3.23835541e-04 * tc[1]
            -1.68899065e-06 * tc[2]
            +3.16217327e-09 * tc[3]
            -1.40609067e-12 * tc[4];
        /*species 9: CH2 */
        species[9] =
            +2.76267867e+00
            +9.68872143e-04 * tc[1]
            +2.79489841e-06 * tc[2]
            -3.85091153e-09 * tc[3]
            +1.68741719e-12 * tc[4];
        /*species 10: CH2* */
        species[10] =
            +3.19860411e+00
            -2.36661419e-03 * tc[1]
            +8.23296220e-06 * tc[2]
            -6.68815981e-09 * tc[3]
            +1.94314737e-12 * tc[4];
        /*species 11: CH3 */
        species[11] =
            +2.67359040e+00
            +2.01095175e-03 * tc[1]
            +5.73021856e-06 * tc[2]
            -6.87117425e-09 * tc[3]
            +2.54385734e-12 * tc[4];
        /*species 12: CH4 */
        species[12] =
            +4.14987613e+00
            -1.36709788e-02 * tc[1]
            +4.91800599e-05 * tc[2]
            -4.84743026e-08 * tc[3]
            +1.66693956e-11 * tc[4];
        /*species 13: HCO */
        species[13] =
            +3.22118584e+00
            -3.24392532e-03 * tc[1]
            +1.37799446e-05 * tc[2]
            -1.33144093e-08 * tc[3]
            +4.33768865e-12 * tc[4];
        /*species 14: CH2O */
        species[14] =
            +3.79372315e+00
            -9.90833369e-03 * tc[1]
            +3.73220008e-05 * tc[2]
            -3.79285261e-08 * tc[3]
            +1.31772652e-11 * tc[4];
        /*species 15: CH3O */
        species[15] =
            +2.71180502e+00
            -2.80463306e-03 * tc[1]
            +3.76550971e-05 * tc[2]
            -4.73072089e-08 * tc[3]
            +1.86588420e-11 * tc[4];
        /*species 16: CH2OH */
        species[16] =
            +3.47834367e+00
            -1.35070310e-03 * tc[1]
            +2.78484980e-05 * tc[2]
            -3.64869060e-08 * tc[3]
            +1.47907450e-11 * tc[4];
        /*species 17: CH3OH */
        species[17] =
            +4.71539582e+00
            -1.52309129e-02 * tc[1]
            +6.52441155e-05 * tc[2]
            -7.10806889e-08 * tc[3]
            +2.61352698e-11 * tc[4];
        /*species 18: CO */
        species[18] =
            +2.57953347e+00
            -6.10353680e-04 * tc[1]
            +1.01681433e-06 * tc[2]
            +9.07005884e-10 * tc[3]
            -9.04424499e-13 * tc[4];
        /*species 19: CO2 */
        species[19] =
            +1.35677352e+00
            +8.98459677e-03 * tc[1]
            -7.12356269e-06 * tc[2]
            +2.45919022e-09 * tc[3]
            -1.43699548e-13 * tc[4];
        /*species 20: C2H */
        species[20] =
            +1.88965733e+00
            +1.34099611e-02 * tc[1]
            -2.84769501e-05 * tc[2]
            +2.94791045e-08 * tc[3]
            -1.09331511e-11 * tc[4];
        /*species 21: C2H2 */
        species[21] =
            -1.91318906e-01
            +2.33615629e-02 * tc[1]
            -3.55171815e-05 * tc[2]
            +2.80152437e-08 * tc[3]
            -8.50072974e-12 * tc[4];
        /*species 22: C2H3 */
        species[22] =
            +2.21246645e+00
            +1.51479162e-03 * tc[1]
            +2.59209412e-05 * tc[2]
            -3.57657847e-08 * tc[3]
            +1.47150873e-11 * tc[4];
        /*species 23: C2H4 */
        species[23] =
            +2.95920148e+00
            -7.57052247e-03 * tc[1]
            +5.70990292e-05 * tc[2]
            -6.91588753e-08 * tc[3]
            +2.69884373e-11 * tc[4];
        /*species 24: C2H5 */
        species[24] =
            +3.30646568e+00
            -4.18658892e-03 * tc[1]
            +4.97142807e-05 * tc[2]
            -5.99126606e-08 * tc[3]
            +2.30509004e-11 * tc[4];
        /*species 25: C2H6 */
        species[25] =
            +3.29142492e+00
            -5.50154270e-03 * tc[1]
            +5.99438288e-05 * tc[2]
            -7.08466285e-08 * tc[3]
            +2.68685771e-11 * tc[4];
        /*species 26: HCCO */
        species[26] =
            +1.25172140e+00
            +1.76550210e-02 * tc[1]
            -2.37291010e-05 * tc[2]
            +1.72757590e-08 * tc[3]
            -5.06648110e-12 * tc[4];
        /*species 27: CH2CO */
        species[27] =
            +1.13583630e+00
            +1.81188721e-02 * tc[1]
            -1.73947474e-05 * tc[2]
            +9.34397568e-09 * tc[3]
            -2.01457615e-12 * tc[4];
        /*species 28: CH3CO */
        species[28] =
            +3.16342570e+00
            -2.32616100e-04 * tc[1]
            +3.42678200e-05 * tc[2]
            -4.41052270e-08 * tc[3]
            +1.72756120e-11 * tc[4];
        /*species 29: CH2CHO */
        species[29] =
            +2.40906240e+00
            +1.07385740e-02 * tc[1]
            +1.89149250e-06 * tc[2]
            -7.15858310e-09 * tc[3]
            +2.86738510e-12 * tc[4];
        /*species 30: CH3CHO */
        species[30] =
            +3.72945950e+00
            -3.19328580e-03 * tc[1]
            +4.75349210e-05 * tc[2]
            -5.74586110e-08 * tc[3]
            +2.19311120e-11 * tc[4];
        /*species 31: C3H3 */
        species[31] =
            +3.51109270e-01
            +3.27411223e-02 * tc[1]
            -4.73827135e-05 * tc[2]
            +3.76309808e-08 * tc[3]
            -1.18540923e-11 * tc[4];
        /*species 32: pC3H4 */
        species[32] =
            +1.68038690e+00
            +1.57996510e-02 * tc[1]
            +2.50705960e-06 * tc[2]
            -1.36576230e-08 * tc[3]
            +6.61542850e-12 * tc[4];
        /*species 33: aC3H4 */
        species[33] =
            +1.61304450e+00
            +1.21225750e-02 * tc[1]
            +1.85398800e-05 * tc[2]
            -3.45251490e-08 * tc[3]
            +1.53350790e-11 * tc[4];
        /*species 34: aC3H5 */
        species[34] =
            +3.63183500e-01
            +1.98138210e-02 * tc[1]
            +1.24970600e-05 * tc[2]
            -3.33555550e-08 * tc[3]
            +1.58465710e-11 * tc[4];
        /*species 35: CH3CCH2 */
        species[35] =
            +7.32920900e-01
            +2.23946200e-02 * tc[1]
            -5.14906110e-06 * tc[2]
            -6.75964660e-09 * tc[3]
            +3.82532110e-12 * tc[4];
        /*species 36: C3H6 */
        species[36] =
            +4.93307000e-01
            +2.09251800e-02 * tc[1]
            +4.48679400e-06 * tc[2]
            -1.66891200e-08 * tc[3]
            +7.15814600e-12 * tc[4];
        /*species 37: nC3H7 */
        species[37] =
            +4.91173000e-02
            +2.60089730e-02 * tc[1]
            +2.35425160e-06 * tc[2]
            -1.95951320e-08 * tc[3]
            +9.37202070e-12 * tc[4];
        /*species 38: iC3H7 */
        species[38] =
            +4.44919900e-01
            +2.09991120e-02 * tc[1]
            +7.70362220e-06 * tc[2]
            -1.84762530e-08 * tc[3]
            +7.12829620e-12 * tc[4];
        /*species 39: C2H3CHO */
        species[39] =
            +2.71349800e-01
            +2.62310540e-02 * tc[1]
            -9.29123050e-06 * tc[2]
            -4.78372720e-09 * tc[3]
            +3.34805430e-12 * tc[4];
        /*species 40: C4H2 */
        species[40] =
            +5.43978000e-02
            +4.16269600e-02 * tc[1]
            -6.58717840e-05 * tc[2]
            +5.32570750e-08 * tc[3]
            -1.66831620e-11 * tc[4];
        /*species 41: iC4H3 */
        species[41] =
            +2.72214820e+00
            +2.59575430e-02 * tc[1]
            -2.63563430e-05 * tc[2]
            +1.55089200e-08 * tc[3]
            -3.80405650e-12 * tc[4];
        /*species 42: C4H4 */
        species[42] =
            -4.11429520e-01
            +3.65466850e-02 * tc[1]
            -3.41069680e-05 * tc[2]
            +1.66526190e-08 * tc[3]
            -3.00646230e-12 * tc[4];
        /*species 43: iC4H5 */
        species[43] =
            -8.86918950e-01
            +4.09506150e-02 * tc[1]
            -3.54135810e-05 * tc[2]
            +1.55309690e-08 * tc[3]
            -2.33551220e-12 * tc[4];
        /*species 44: C4H5-2 */
        species[44] =
            +1.96962800e+00
            +2.44422450e-02 * tc[1]
            -9.12514240e-06 * tc[2]
            -4.24668710e-18 * tc[3]
            +1.63047280e-21 * tc[4];
        /*species 45: C4H6 */
        species[45] =
            -8.87155350e-01
            +3.43690220e-02 * tc[1]
            -1.11073920e-05 * tc[2]
            -9.21066600e-09 * tc[3]
            +6.20651790e-12 * tc[4];
        /*species 46: C4H612 */
        species[46] =
            +2.34670000e-02
            +3.49591900e-02 * tc[1]
            -2.20090500e-05 * tc[2]
            +6.94227200e-09 * tc[3]
            -7.87918700e-13 * tc[4];
        /*species 47: C4H6-2 */
        species[47] =
            +1.13733380e+00
            +2.64862290e-02 * tc[1]
            -9.05687110e-06 * tc[2]
            -5.53863970e-19 * tc[3]
            +2.12818840e-22 * tc[4];
        /*species 48: C4H7 */
        species[48] =
            -2.55505680e-01
            +3.96788570e-02 * tc[1]
            -2.28980860e-05 * tc[2]
            +2.13529730e-09 * tc[3]
            +2.30963750e-12 * tc[4];
        /*species 49: C4H81 */
        species[49] =
            +1.81138000e-01
            +3.08533800e-02 * tc[1]
            +5.08652470e-06 * tc[2]
            -2.46548880e-08 * tc[3]
            +1.11101930e-11 * tc[4];
        /*species 50: pC4H9 */
        species[50] =
            +2.08704200e-01
            +3.82974970e-02 * tc[1]
            -7.26605090e-06 * tc[2]
            -1.54285470e-08 * tc[3]
            +8.68594350e-12 * tc[4];
        /*species 55: N2 */
        species[55] =
            +2.29867700e+00
            +1.40824040e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485400e-12 * tc[4];
    } else {
        /*species 0: H */
        species[0] =
            +1.50000001e+00
            -2.30842973e-11 * tc[1]
            +1.61561948e-14 * tc[2]
            -4.73515235e-18 * tc[3]
            +4.98197357e-22 * tc[4];
        /*species 1: O */
        species[1] =
            +1.56942078e+00
            -8.59741137e-05 * tc[1]
            +4.19484589e-08 * tc[2]
            -1.00177799e-11 * tc[3]
            +1.22833691e-15 * tc[4];
        /*species 2: OH */
        species[2] =
            +1.86472886e+00
            +1.05650448e-03 * tc[1]
            -2.59082758e-07 * tc[2]
            +3.05218674e-11 * tc[3]
            -1.33195876e-15 * tc[4];
        /*species 3: HO2 */
        species[3] =
            +3.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
        /*species 4: H2 */
        species[4] =
            +2.33727920e+00
            -4.94024731e-05 * tc[1]
            +4.99456778e-07 * tc[2]
            -1.79566394e-10 * tc[3]
            +2.00255376e-14 * tc[4];
        /*species 5: H2O */
        species[5] =
            +2.03399249e+00
            +2.17691804e-03 * tc[1]
            -1.64072518e-07 * tc[2]
            -9.70419870e-11 * tc[3]
            +1.68200992e-14 * tc[4];
        /*species 6: H2O2 */
        species[6] =
            +3.16500285e+00
            +4.90831694e-03 * tc[1]
            -1.90139225e-06 * tc[2]
            +3.71185986e-10 * tc[3]
            -2.87908305e-14 * tc[4];
        /*species 7: O2 */
        species[7] =
            +2.28253784e+00
            +1.48308754e-03 * tc[1]
            -7.57966669e-07 * tc[2]
            +2.09470555e-10 * tc[3]
            -2.16717794e-14 * tc[4];
        /*species 8: CH */
        species[8] =
            +1.87846473e+00
            +9.70913681e-04 * tc[1]
            +1.44445655e-07 * tc[2]
            -1.30687849e-10 * tc[3]
            +1.76079383e-14 * tc[4];
        /*species 9: CH2 */
        species[9] =
            +1.87410113e+00
            +3.65639292e-03 * tc[1]
            -1.40894597e-06 * tc[2]
            +2.60179549e-10 * tc[3]
            -1.87727567e-14 * tc[4];
        /*species 10: CH2* */
        species[10] =
            +1.29203842e+00
            +4.65588637e-03 * tc[1]
            -2.01191947e-06 * tc[2]
            +4.17906000e-10 * tc[3]
            -3.39716365e-14 * tc[4];
        /*species 11: CH3 */
        species[11] =
            +1.28571772e+00
            +7.23990037e-03 * tc[1]
            -2.98714348e-06 * tc[2]
            +5.95684644e-10 * tc[3]
            -4.67154394e-14 * tc[4];
        /*species 12: CH4 */
        species[12] =
            -9.25148505e-01
            +1.33909467e-02 * tc[1]
            -5.73285809e-06 * tc[2]
            +1.22292535e-09 * tc[3]
            -1.01815230e-13 * tc[4];
        /*species 13: HCO */
        species[13] =
            +1.77217438e+00
            +4.95695526e-03 * tc[1]
            -2.48445613e-06 * tc[2]
            +5.89161778e-10 * tc[3]
            -5.33508711e-14 * tc[4];
        /*species 14: CH2O */
        species[14] =
            +7.60690080e-01
            +9.20000082e-03 * tc[1]
            -4.42258813e-06 * tc[2]
            +1.00641212e-09 * tc[3]
            -8.83855640e-14 * tc[4];
        /*species 15: CH3O */
        species[15] =
            +3.75779238e+00
            +7.44142474e-03 * tc[1]
            -2.69705176e-06 * tc[2]
            +4.38090504e-10 * tc[3]
            -2.63537098e-14 * tc[4];
        /*species 16: CH2OH */
        species[16] =
            +4.09314370e+00
            +5.94761260e-03 * tc[1]
            -2.06497460e-06 * tc[2]
            +3.23008173e-10 * tc[3]
            -1.88125902e-14 * tc[4];
        /*species 17: CH3OH */
        species[17] =
            +7.89707910e-01
            +1.40938292e-02 * tc[1]
            -6.36500835e-06 * tc[2]
            +1.38171085e-09 * tc[3]
            -1.17060220e-13 * tc[4];
        /*species 18: CO */
        species[18] =
            +1.71518561e+00
            +2.06252743e-03 * tc[1]
            -9.98825771e-07 * tc[2]
            +2.30053008e-10 * tc[3]
            -2.03647716e-14 * tc[4];
        /*species 19: CO2 */
        species[19] =
            +2.85746029e+00
            +4.41437026e-03 * tc[1]
            -2.21481404e-06 * tc[2]
            +5.23490188e-10 * tc[3]
            -4.72084164e-14 * tc[4];
        /*species 20: C2H */
        species[20] =
            +2.16780652e+00
            +4.75221902e-03 * tc[1]
            -1.83787077e-06 * tc[2]
            +3.04190252e-10 * tc[3]
            -1.77232770e-14 * tc[4];
        /*species 21: C2H2 */
        species[21] =
            +3.14756964e+00
            +5.96166664e-03 * tc[1]
            -2.37294852e-06 * tc[2]
            +4.67412171e-10 * tc[3]
            -3.61235213e-14 * tc[4];
        /*species 22: C2H3 */
        species[22] =
            +2.01672400e+00
            +1.03302292e-02 * tc[1]
            -4.68082349e-06 * tc[2]
            +1.01763288e-09 * tc[3]
            -8.62607041e-14 * tc[4];
        /*species 23: C2H4 */
        species[23] =
            +1.03611116e+00
            +1.46454151e-02 * tc[1]
            -6.71077915e-06 * tc[2]
            +1.47222923e-09 * tc[3]
            -1.25706061e-13 * tc[4];
        /*species 24: C2H5 */
        species[24] =
            +9.54656420e-01
            +1.73972722e-02 * tc[1]
            -7.98206668e-06 * tc[2]
            +1.75217689e-09 * tc[3]
            -1.49641576e-13 * tc[4];
        /*species 25: C2H6 */
        species[25] =
            +7.18815000e-02
            +2.16852677e-02 * tc[1]
            -1.00256067e-05 * tc[2]
            +2.21412001e-09 * tc[3]
            -1.90002890e-13 * tc[4];
        /*species 26: HCCO */
        species[26] =
            +4.62820580e+00
            +4.08534010e-03 * tc[1]
            -1.59345470e-06 * tc[2]
            +2.86260520e-10 * tc[3]
            -1.94078320e-14 * tc[4];
        /*species 27: CH2CO */
        species[27] =
            +3.51129732e+00
            +9.00359745e-03 * tc[1]
            -4.16939635e-06 * tc[2]
            +9.23345882e-10 * tc[3]
            -7.94838201e-14 * tc[4];
        /*species 28: CH3CO */
        species[28] =
            +4.94477310e+00
            +7.86672050e-03 * tc[1]
            -2.88658820e-06 * tc[2]
            +4.72708750e-10 * tc[3]
            -2.85998610e-14 * tc[4];
        /*species 29: CH2CHO */
        species[29] =
            +4.97566990e+00
            +8.13059140e-03 * tc[1]
            -2.74362450e-06 * tc[2]
            +4.07030410e-10 * tc[3]
            -2.17601710e-14 * tc[4];
        /*species 30: CH3CHO */
        species[30] =
            +4.40411080e+00
            +1.17230590e-02 * tc[1]
            -4.22631370e-06 * tc[2]
            +6.83724510e-10 * tc[3]
            -4.09848630e-14 * tc[4];
        /*species 31: C3H3 */
        species[31] =
            +6.14221880e+00
            +7.61902005e-03 * tc[1]
            -2.67459950e-06 * tc[2]
            +4.24914801e-10 * tc[3]
            -2.51475415e-14 * tc[4];
        /*species 32: pC3H4 */
        species[32] =
            +5.02524000e+00
            +1.13365420e-02 * tc[1]
            -4.02233910e-06 * tc[2]
            +6.43760630e-10 * tc[3]
            -3.82996350e-14 * tc[4];
        /*species 33: aC3H4 */
        species[33] =
            +5.31687220e+00
            +1.11337280e-02 * tc[1]
            -3.96293780e-06 * tc[2]
            +6.35642380e-10 * tc[3]
            -3.78755400e-14 * tc[4];
        /*species 34: aC3H5 */
        species[34] =
            +5.50078770e+00
            +1.43247310e-02 * tc[1]
            -5.67816320e-06 * tc[2]
            +1.10808010e-09 * tc[3]
            -9.03638870e-14 * tc[4];
        /*species 35: CH3CCH2 */
        species[35] =
            +4.42555280e+00
            +1.55110720e-02 * tc[1]
            -5.66783500e-06 * tc[2]
            +7.92243880e-10 * tc[3]
            -1.68780340e-14 * tc[4];
        /*species 36: C3H6 */
        species[36] =
            +5.73225700e+00
            +1.49083400e-02 * tc[1]
            -4.94989900e-06 * tc[2]
            +7.21202200e-10 * tc[3]
            -3.76620400e-14 * tc[4];
        /*species 37: nC3H7 */
        species[37] =
            +6.70974790e+00
            +1.60314850e-02 * tc[1]
            -5.27202380e-06 * tc[2]
            +7.58883520e-10 * tc[3]
            -3.88627190e-14 * tc[4];
        /*species 38: iC3H7 */
        species[38] =
            +5.51927410e+00
            +1.72201040e-02 * tc[1]
            -5.73642170e-06 * tc[2]
            +8.41307320e-10 * tc[3]
            -4.45659130e-14 * tc[4];
        /*species 39: C2H3CHO */
        species[39] =
            +4.81118680e+00
            +1.71142560e-02 * tc[1]
            -7.48341610e-06 * tc[2]
            +1.42522490e-09 * tc[3]
            -9.17468410e-14 * tc[4];
        /*species 40: C4H2 */
        species[40] =
            +8.15763280e+00
            +5.54305180e-03 * tc[1]
            -1.35916040e-06 * tc[2]
            +1.87800750e-11 * tc[3]
            +2.31895360e-14 * tc[4];
        /*species 41: iC4H3 */
        species[41] =
            +6.65385480e+00
            +1.12040550e-02 * tc[1]
            -4.64013420e-06 * tc[2]
            +8.67866390e-10 * tc[3]
            -5.74305620e-14 * tc[4];
        /*species 42: C4H4 */
        species[42] =
            +6.25396010e+00
            +1.39140940e-02 * tc[1]
            -5.29322140e-06 * tc[2]
            +8.34804500e-10 * tc[3]
            -3.51978820e-14 * tc[4];
        /*species 43: iC4H5 */
        species[43] =
            +5.96460290e+00
            +1.82743330e-02 * tc[1]
            -7.81337350e-06 * tc[2]
            +1.52921540e-09 * tc[3]
            -1.09204930e-13 * tc[4];
        /*species 44: C4H5-2 */
        species[44] =
            +1.35381710e+01
            -8.56770560e-03 * tc[1]
            +2.35595240e-05 * tc[2]
            -1.36763790e-08 * tc[3]
            +2.44369270e-12 * tc[4];
        /*species 45: C4H6 */
        species[45] =
            +7.86731340e+00
            +1.49186700e-02 * tc[1]
            -3.15487160e-06 * tc[2]
            -4.18413300e-10 * tc[3]
            +1.57612580e-13 * tc[4];
        /*species 46: C4H612 */
        species[46] =
            +1.68155700e+01
            -4.25750200e-03 * tc[1]
            +1.05118500e-05 * tc[2]
            -4.47384400e-09 * tc[3]
            +5.84813800e-13 * tc[4];
        /*species 47: C4H6-2 */
        species[47] =
            +8.03381330e+00
            +8.21245100e-03 * tc[1]
            +7.17539520e-06 * tc[2]
            -5.88343340e-09 * tc[3]
            +1.03439150e-12 * tc[4];
        /*species 48: C4H7 */
        species[48] =
            +6.01348350e+00
            +2.26345580e-02 * tc[1]
            -9.25454700e-06 * tc[2]
            +1.68079270e-09 * tc[3]
            -1.04086170e-13 * tc[4];
        /*species 49: C4H81 */
        species[49] =
            +1.05358410e+00
            +3.43505070e-02 * tc[1]
            -1.58831970e-05 * tc[2]
            +3.30896620e-09 * tc[3]
            -2.53610450e-13 * tc[4];
        /*species 50: pC4H9 */
        species[50] =
            +7.68223950e+00
            +2.36910710e-02 * tc[1]
            -7.59488650e-06 * tc[2]
            +6.64271360e-10 * tc[3]
            +5.48451360e-14 * tc[4];
        /*species 55: N2 */
        species[55] =
            +1.92664000e+00
            +1.48797680e-03 * tc[1]
            -5.68476000e-07 * tc[2]
            +1.00970380e-10 * tc[3]
            -6.75335100e-15 * tc[4];
    }

    /*species with midpoint at T=1392 kelvin */
    if (T < 1392) {
        /*species 52: C6H12 */
        species[52] =
            -2.35275205e+00
            +6.98655426e-02 * tc[1]
            -4.59408022e-05 * tc[2]
            +1.56967343e-08 * tc[3]
            -2.21296175e-12 * tc[4];
        /*species 54: C5H10 */
        species[54] =
            -2.06223481e+00
            +5.74218294e-02 * tc[1]
            -3.74486890e-05 * tc[2]
            +1.27364989e-08 * tc[3]
            -1.79609789e-12 * tc[4];
    } else {
        /*species 52: C6H12 */
        species[52] =
            +1.68337529e+01
            +2.67377658e-02 * tc[1]
            -9.10036773e-06 * tc[2]
            +1.40819768e-09 * tc[3]
            -8.15124244e-14 * tc[4];
        /*species 54: C5H10 */
        species[54] =
            +1.35851539e+01
            +2.24072471e-02 * tc[1]
            -7.63348025e-06 * tc[2]
            +1.18188966e-09 * tc[3]
            -6.84385139e-14 * tc[4];
    }

    /*species with midpoint at T=1389 kelvin */
    if (T < 1389) {
        /*species 53: C6H11 */
        species[53] =
            -2.55544944e+00
            +6.76865602e-02 * tc[1]
            -4.47048635e-05 * tc[2]
            +1.52236630e-08 * tc[3]
            -2.14346377e-12 * tc[4];
    } else {
        /*species 53: C6H11 */
        species[53] =
            +1.67336550e+01
            +2.48934775e-02 * tc[1]
            -8.59991450e-06 * tc[2]
            +1.34412828e-09 * tc[3]
            -7.83475666e-14 * tc[4];
    }

    /*species with midpoint at T=1391 kelvin */
    if (T < 1391) {
        /*species 51: NC12H26 */
        species[51] =
            -3.62181594e+00
            +1.47237711e-01 * tc[1]
            -9.43970271e-05 * tc[2]
            +3.07441268e-08 * tc[3]
            -4.03602230e-12 * tc[4];
    } else {
        /*species 51: NC12H26 */
        species[51] =
            +3.75095037e+01
            +5.63550048e-02 * tc[1]
            -1.91493200e-05 * tc[2]
            +2.96024862e-09 * tc[3]
            -1.71244150e-13 * tc[4];
    }
    return;
}


/*compute Cp/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
__device__ void cp_R_d(double * species, double *  tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +2.50000000e+00
            +7.05332819e-13 * tc[1]
            -1.99591964e-15 * tc[2]
            +2.30081632e-18 * tc[3]
            -9.27732332e-22 * tc[4];
        /*species 1: O */
        species[1] =
            +3.16826710e+00
            -3.27931884e-03 * tc[1]
            +6.64306396e-06 * tc[2]
            -6.12806624e-09 * tc[3]
            +2.11265971e-12 * tc[4];
        /*species 2: OH */
        species[2] =
            +4.12530561e+00
            -3.22544939e-03 * tc[1]
            +6.52764691e-06 * tc[2]
            -5.79853643e-09 * tc[3]
            +2.06237379e-12 * tc[4];
        /*species 3: HO2 */
        species[3] =
            +4.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
        /*species 4: H2 */
        species[4] =
            +2.34433112e+00
            +7.98052075e-03 * tc[1]
            -1.94781510e-05 * tc[2]
            +2.01572094e-08 * tc[3]
            -7.37611761e-12 * tc[4];
        /*species 5: H2O */
        species[5] =
            +4.19864056e+00
            -2.03643410e-03 * tc[1]
            +6.52040211e-06 * tc[2]
            -5.48797062e-09 * tc[3]
            +1.77197817e-12 * tc[4];
        /*species 6: H2O2 */
        species[6] =
            +4.27611269e+00
            -5.42822417e-04 * tc[1]
            +1.67335701e-05 * tc[2]
            -2.15770813e-08 * tc[3]
            +8.62454363e-12 * tc[4];
        /*species 7: O2 */
        species[7] =
            +3.78245636e+00
            -2.99673416e-03 * tc[1]
            +9.84730201e-06 * tc[2]
            -9.68129509e-09 * tc[3]
            +3.24372837e-12 * tc[4];
        /*species 8: CH */
        species[8] =
            +3.48981665e+00
            +3.23835541e-04 * tc[1]
            -1.68899065e-06 * tc[2]
            +3.16217327e-09 * tc[3]
            -1.40609067e-12 * tc[4];
        /*species 9: CH2 */
        species[9] =
            +3.76267867e+00
            +9.68872143e-04 * tc[1]
            +2.79489841e-06 * tc[2]
            -3.85091153e-09 * tc[3]
            +1.68741719e-12 * tc[4];
        /*species 10: CH2* */
        species[10] =
            +4.19860411e+00
            -2.36661419e-03 * tc[1]
            +8.23296220e-06 * tc[2]
            -6.68815981e-09 * tc[3]
            +1.94314737e-12 * tc[4];
        /*species 11: CH3 */
        species[11] =
            +3.67359040e+00
            +2.01095175e-03 * tc[1]
            +5.73021856e-06 * tc[2]
            -6.87117425e-09 * tc[3]
            +2.54385734e-12 * tc[4];
        /*species 12: CH4 */
        species[12] =
            +5.14987613e+00
            -1.36709788e-02 * tc[1]
            +4.91800599e-05 * tc[2]
            -4.84743026e-08 * tc[3]
            +1.66693956e-11 * tc[4];
        /*species 13: HCO */
        species[13] =
            +4.22118584e+00
            -3.24392532e-03 * tc[1]
            +1.37799446e-05 * tc[2]
            -1.33144093e-08 * tc[3]
            +4.33768865e-12 * tc[4];
        /*species 14: CH2O */
        species[14] =
            +4.79372315e+00
            -9.90833369e-03 * tc[1]
            +3.73220008e-05 * tc[2]
            -3.79285261e-08 * tc[3]
            +1.31772652e-11 * tc[4];
        /*species 15: CH3O */
        species[15] =
            +3.71180502e+00
            -2.80463306e-03 * tc[1]
            +3.76550971e-05 * tc[2]
            -4.73072089e-08 * tc[3]
            +1.86588420e-11 * tc[4];
        /*species 16: CH2OH */
        species[16] =
            +4.47834367e+00
            -1.35070310e-03 * tc[1]
            +2.78484980e-05 * tc[2]
            -3.64869060e-08 * tc[3]
            +1.47907450e-11 * tc[4];
        /*species 17: CH3OH */
        species[17] =
            +5.71539582e+00
            -1.52309129e-02 * tc[1]
            +6.52441155e-05 * tc[2]
            -7.10806889e-08 * tc[3]
            +2.61352698e-11 * tc[4];
        /*species 18: CO */
        species[18] =
            +3.57953347e+00
            -6.10353680e-04 * tc[1]
            +1.01681433e-06 * tc[2]
            +9.07005884e-10 * tc[3]
            -9.04424499e-13 * tc[4];
        /*species 19: CO2 */
        species[19] =
            +2.35677352e+00
            +8.98459677e-03 * tc[1]
            -7.12356269e-06 * tc[2]
            +2.45919022e-09 * tc[3]
            -1.43699548e-13 * tc[4];
        /*species 20: C2H */
        species[20] =
            +2.88965733e+00
            +1.34099611e-02 * tc[1]
            -2.84769501e-05 * tc[2]
            +2.94791045e-08 * tc[3]
            -1.09331511e-11 * tc[4];
        /*species 21: C2H2 */
        species[21] =
            +8.08681094e-01
            +2.33615629e-02 * tc[1]
            -3.55171815e-05 * tc[2]
            +2.80152437e-08 * tc[3]
            -8.50072974e-12 * tc[4];
        /*species 22: C2H3 */
        species[22] =
            +3.21246645e+00
            +1.51479162e-03 * tc[1]
            +2.59209412e-05 * tc[2]
            -3.57657847e-08 * tc[3]
            +1.47150873e-11 * tc[4];
        /*species 23: C2H4 */
        species[23] =
            +3.95920148e+00
            -7.57052247e-03 * tc[1]
            +5.70990292e-05 * tc[2]
            -6.91588753e-08 * tc[3]
            +2.69884373e-11 * tc[4];
        /*species 24: C2H5 */
        species[24] =
            +4.30646568e+00
            -4.18658892e-03 * tc[1]
            +4.97142807e-05 * tc[2]
            -5.99126606e-08 * tc[3]
            +2.30509004e-11 * tc[4];
        /*species 25: C2H6 */
        species[25] =
            +4.29142492e+00
            -5.50154270e-03 * tc[1]
            +5.99438288e-05 * tc[2]
            -7.08466285e-08 * tc[3]
            +2.68685771e-11 * tc[4];
        /*species 26: HCCO */
        species[26] =
            +2.25172140e+00
            +1.76550210e-02 * tc[1]
            -2.37291010e-05 * tc[2]
            +1.72757590e-08 * tc[3]
            -5.06648110e-12 * tc[4];
        /*species 27: CH2CO */
        species[27] =
            +2.13583630e+00
            +1.81188721e-02 * tc[1]
            -1.73947474e-05 * tc[2]
            +9.34397568e-09 * tc[3]
            -2.01457615e-12 * tc[4];
        /*species 28: CH3CO */
        species[28] =
            +4.16342570e+00
            -2.32616100e-04 * tc[1]
            +3.42678200e-05 * tc[2]
            -4.41052270e-08 * tc[3]
            +1.72756120e-11 * tc[4];
        /*species 29: CH2CHO */
        species[29] =
            +3.40906240e+00
            +1.07385740e-02 * tc[1]
            +1.89149250e-06 * tc[2]
            -7.15858310e-09 * tc[3]
            +2.86738510e-12 * tc[4];
        /*species 30: CH3CHO */
        species[30] =
            +4.72945950e+00
            -3.19328580e-03 * tc[1]
            +4.75349210e-05 * tc[2]
            -5.74586110e-08 * tc[3]
            +2.19311120e-11 * tc[4];
        /*species 31: C3H3 */
        species[31] =
            +1.35110927e+00
            +3.27411223e-02 * tc[1]
            -4.73827135e-05 * tc[2]
            +3.76309808e-08 * tc[3]
            -1.18540923e-11 * tc[4];
        /*species 32: pC3H4 */
        species[32] =
            +2.68038690e+00
            +1.57996510e-02 * tc[1]
            +2.50705960e-06 * tc[2]
            -1.36576230e-08 * tc[3]
            +6.61542850e-12 * tc[4];
        /*species 33: aC3H4 */
        species[33] =
            +2.61304450e+00
            +1.21225750e-02 * tc[1]
            +1.85398800e-05 * tc[2]
            -3.45251490e-08 * tc[3]
            +1.53350790e-11 * tc[4];
        /*species 34: aC3H5 */
        species[34] =
            +1.36318350e+00
            +1.98138210e-02 * tc[1]
            +1.24970600e-05 * tc[2]
            -3.33555550e-08 * tc[3]
            +1.58465710e-11 * tc[4];
        /*species 35: CH3CCH2 */
        species[35] =
            +1.73292090e+00
            +2.23946200e-02 * tc[1]
            -5.14906110e-06 * tc[2]
            -6.75964660e-09 * tc[3]
            +3.82532110e-12 * tc[4];
        /*species 36: C3H6 */
        species[36] =
            +1.49330700e+00
            +2.09251800e-02 * tc[1]
            +4.48679400e-06 * tc[2]
            -1.66891200e-08 * tc[3]
            +7.15814600e-12 * tc[4];
        /*species 37: nC3H7 */
        species[37] =
            +1.04911730e+00
            +2.60089730e-02 * tc[1]
            +2.35425160e-06 * tc[2]
            -1.95951320e-08 * tc[3]
            +9.37202070e-12 * tc[4];
        /*species 38: iC3H7 */
        species[38] =
            +1.44491990e+00
            +2.09991120e-02 * tc[1]
            +7.70362220e-06 * tc[2]
            -1.84762530e-08 * tc[3]
            +7.12829620e-12 * tc[4];
        /*species 39: C2H3CHO */
        species[39] =
            +1.27134980e+00
            +2.62310540e-02 * tc[1]
            -9.29123050e-06 * tc[2]
            -4.78372720e-09 * tc[3]
            +3.34805430e-12 * tc[4];
        /*species 40: C4H2 */
        species[40] =
            +1.05439780e+00
            +4.16269600e-02 * tc[1]
            -6.58717840e-05 * tc[2]
            +5.32570750e-08 * tc[3]
            -1.66831620e-11 * tc[4];
        /*species 41: iC4H3 */
        species[41] =
            +3.72214820e+00
            +2.59575430e-02 * tc[1]
            -2.63563430e-05 * tc[2]
            +1.55089200e-08 * tc[3]
            -3.80405650e-12 * tc[4];
        /*species 42: C4H4 */
        species[42] =
            +5.88570480e-01
            +3.65466850e-02 * tc[1]
            -3.41069680e-05 * tc[2]
            +1.66526190e-08 * tc[3]
            -3.00646230e-12 * tc[4];
        /*species 43: iC4H5 */
        species[43] =
            +1.13081050e-01
            +4.09506150e-02 * tc[1]
            -3.54135810e-05 * tc[2]
            +1.55309690e-08 * tc[3]
            -2.33551220e-12 * tc[4];
        /*species 44: C4H5-2 */
        species[44] =
            +2.96962800e+00
            +2.44422450e-02 * tc[1]
            -9.12514240e-06 * tc[2]
            -4.24668710e-18 * tc[3]
            +1.63047280e-21 * tc[4];
        /*species 45: C4H6 */
        species[45] =
            +1.12844650e-01
            +3.43690220e-02 * tc[1]
            -1.11073920e-05 * tc[2]
            -9.21066600e-09 * tc[3]
            +6.20651790e-12 * tc[4];
        /*species 46: C4H612 */
        species[46] =
            +1.02346700e+00
            +3.49591900e-02 * tc[1]
            -2.20090500e-05 * tc[2]
            +6.94227200e-09 * tc[3]
            -7.87918700e-13 * tc[4];
        /*species 47: C4H6-2 */
        species[47] =
            +2.13733380e+00
            +2.64862290e-02 * tc[1]
            -9.05687110e-06 * tc[2]
            -5.53863970e-19 * tc[3]
            +2.12818840e-22 * tc[4];
        /*species 48: C4H7 */
        species[48] =
            +7.44494320e-01
            +3.96788570e-02 * tc[1]
            -2.28980860e-05 * tc[2]
            +2.13529730e-09 * tc[3]
            +2.30963750e-12 * tc[4];
        /*species 49: C4H81 */
        species[49] =
            +1.18113800e+00
            +3.08533800e-02 * tc[1]
            +5.08652470e-06 * tc[2]
            -2.46548880e-08 * tc[3]
            +1.11101930e-11 * tc[4];
        /*species 50: pC4H9 */
        species[50] =
            +1.20870420e+00
            +3.82974970e-02 * tc[1]
            -7.26605090e-06 * tc[2]
            -1.54285470e-08 * tc[3]
            +8.68594350e-12 * tc[4];
        /*species 55: N2 */
        species[55] =
            +3.29867700e+00
            +1.40824040e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485400e-12 * tc[4];
    } else {
        /*species 0: H */
        species[0] =
            +2.50000001e+00
            -2.30842973e-11 * tc[1]
            +1.61561948e-14 * tc[2]
            -4.73515235e-18 * tc[3]
            +4.98197357e-22 * tc[4];
        /*species 1: O */
        species[1] =
            +2.56942078e+00
            -8.59741137e-05 * tc[1]
            +4.19484589e-08 * tc[2]
            -1.00177799e-11 * tc[3]
            +1.22833691e-15 * tc[4];
        /*species 2: OH */
        species[2] =
            +2.86472886e+00
            +1.05650448e-03 * tc[1]
            -2.59082758e-07 * tc[2]
            +3.05218674e-11 * tc[3]
            -1.33195876e-15 * tc[4];
        /*species 3: HO2 */
        species[3] =
            +4.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
        /*species 4: H2 */
        species[4] =
            +3.33727920e+00
            -4.94024731e-05 * tc[1]
            +4.99456778e-07 * tc[2]
            -1.79566394e-10 * tc[3]
            +2.00255376e-14 * tc[4];
        /*species 5: H2O */
        species[5] =
            +3.03399249e+00
            +2.17691804e-03 * tc[1]
            -1.64072518e-07 * tc[2]
            -9.70419870e-11 * tc[3]
            +1.68200992e-14 * tc[4];
        /*species 6: H2O2 */
        species[6] =
            +4.16500285e+00
            +4.90831694e-03 * tc[1]
            -1.90139225e-06 * tc[2]
            +3.71185986e-10 * tc[3]
            -2.87908305e-14 * tc[4];
        /*species 7: O2 */
        species[7] =
            +3.28253784e+00
            +1.48308754e-03 * tc[1]
            -7.57966669e-07 * tc[2]
            +2.09470555e-10 * tc[3]
            -2.16717794e-14 * tc[4];
        /*species 8: CH */
        species[8] =
            +2.87846473e+00
            +9.70913681e-04 * tc[1]
            +1.44445655e-07 * tc[2]
            -1.30687849e-10 * tc[3]
            +1.76079383e-14 * tc[4];
        /*species 9: CH2 */
        species[9] =
            +2.87410113e+00
            +3.65639292e-03 * tc[1]
            -1.40894597e-06 * tc[2]
            +2.60179549e-10 * tc[3]
            -1.87727567e-14 * tc[4];
        /*species 10: CH2* */
        species[10] =
            +2.29203842e+00
            +4.65588637e-03 * tc[1]
            -2.01191947e-06 * tc[2]
            +4.17906000e-10 * tc[3]
            -3.39716365e-14 * tc[4];
        /*species 11: CH3 */
        species[11] =
            +2.28571772e+00
            +7.23990037e-03 * tc[1]
            -2.98714348e-06 * tc[2]
            +5.95684644e-10 * tc[3]
            -4.67154394e-14 * tc[4];
        /*species 12: CH4 */
        species[12] =
            +7.48514950e-02
            +1.33909467e-02 * tc[1]
            -5.73285809e-06 * tc[2]
            +1.22292535e-09 * tc[3]
            -1.01815230e-13 * tc[4];
        /*species 13: HCO */
        species[13] =
            +2.77217438e+00
            +4.95695526e-03 * tc[1]
            -2.48445613e-06 * tc[2]
            +5.89161778e-10 * tc[3]
            -5.33508711e-14 * tc[4];
        /*species 14: CH2O */
        species[14] =
            +1.76069008e+00
            +9.20000082e-03 * tc[1]
            -4.42258813e-06 * tc[2]
            +1.00641212e-09 * tc[3]
            -8.83855640e-14 * tc[4];
        /*species 15: CH3O */
        species[15] =
            +4.75779238e+00
            +7.44142474e-03 * tc[1]
            -2.69705176e-06 * tc[2]
            +4.38090504e-10 * tc[3]
            -2.63537098e-14 * tc[4];
        /*species 16: CH2OH */
        species[16] =
            +5.09314370e+00
            +5.94761260e-03 * tc[1]
            -2.06497460e-06 * tc[2]
            +3.23008173e-10 * tc[3]
            -1.88125902e-14 * tc[4];
        /*species 17: CH3OH */
        species[17] =
            +1.78970791e+00
            +1.40938292e-02 * tc[1]
            -6.36500835e-06 * tc[2]
            +1.38171085e-09 * tc[3]
            -1.17060220e-13 * tc[4];
        /*species 18: CO */
        species[18] =
            +2.71518561e+00
            +2.06252743e-03 * tc[1]
            -9.98825771e-07 * tc[2]
            +2.30053008e-10 * tc[3]
            -2.03647716e-14 * tc[4];
        /*species 19: CO2 */
        species[19] =
            +3.85746029e+00
            +4.41437026e-03 * tc[1]
            -2.21481404e-06 * tc[2]
            +5.23490188e-10 * tc[3]
            -4.72084164e-14 * tc[4];
        /*species 20: C2H */
        species[20] =
            +3.16780652e+00
            +4.75221902e-03 * tc[1]
            -1.83787077e-06 * tc[2]
            +3.04190252e-10 * tc[3]
            -1.77232770e-14 * tc[4];
        /*species 21: C2H2 */
        species[21] =
            +4.14756964e+00
            +5.96166664e-03 * tc[1]
            -2.37294852e-06 * tc[2]
            +4.67412171e-10 * tc[3]
            -3.61235213e-14 * tc[4];
        /*species 22: C2H3 */
        species[22] =
            +3.01672400e+00
            +1.03302292e-02 * tc[1]
            -4.68082349e-06 * tc[2]
            +1.01763288e-09 * tc[3]
            -8.62607041e-14 * tc[4];
        /*species 23: C2H4 */
        species[23] =
            +2.03611116e+00
            +1.46454151e-02 * tc[1]
            -6.71077915e-06 * tc[2]
            +1.47222923e-09 * tc[3]
            -1.25706061e-13 * tc[4];
        /*species 24: C2H5 */
        species[24] =
            +1.95465642e+00
            +1.73972722e-02 * tc[1]
            -7.98206668e-06 * tc[2]
            +1.75217689e-09 * tc[3]
            -1.49641576e-13 * tc[4];
        /*species 25: C2H6 */
        species[25] =
            +1.07188150e+00
            +2.16852677e-02 * tc[1]
            -1.00256067e-05 * tc[2]
            +2.21412001e-09 * tc[3]
            -1.90002890e-13 * tc[4];
        /*species 26: HCCO */
        species[26] =
            +5.62820580e+00
            +4.08534010e-03 * tc[1]
            -1.59345470e-06 * tc[2]
            +2.86260520e-10 * tc[3]
            -1.94078320e-14 * tc[4];
        /*species 27: CH2CO */
        species[27] =
            +4.51129732e+00
            +9.00359745e-03 * tc[1]
            -4.16939635e-06 * tc[2]
            +9.23345882e-10 * tc[3]
            -7.94838201e-14 * tc[4];
        /*species 28: CH3CO */
        species[28] =
            +5.94477310e+00
            +7.86672050e-03 * tc[1]
            -2.88658820e-06 * tc[2]
            +4.72708750e-10 * tc[3]
            -2.85998610e-14 * tc[4];
        /*species 29: CH2CHO */
        species[29] =
            +5.97566990e+00
            +8.13059140e-03 * tc[1]
            -2.74362450e-06 * tc[2]
            +4.07030410e-10 * tc[3]
            -2.17601710e-14 * tc[4];
        /*species 30: CH3CHO */
        species[30] =
            +5.40411080e+00
            +1.17230590e-02 * tc[1]
            -4.22631370e-06 * tc[2]
            +6.83724510e-10 * tc[3]
            -4.09848630e-14 * tc[4];
        /*species 31: C3H3 */
        species[31] =
            +7.14221880e+00
            +7.61902005e-03 * tc[1]
            -2.67459950e-06 * tc[2]
            +4.24914801e-10 * tc[3]
            -2.51475415e-14 * tc[4];
        /*species 32: pC3H4 */
        species[32] =
            +6.02524000e+00
            +1.13365420e-02 * tc[1]
            -4.02233910e-06 * tc[2]
            +6.43760630e-10 * tc[3]
            -3.82996350e-14 * tc[4];
        /*species 33: aC3H4 */
        species[33] =
            +6.31687220e+00
            +1.11337280e-02 * tc[1]
            -3.96293780e-06 * tc[2]
            +6.35642380e-10 * tc[3]
            -3.78755400e-14 * tc[4];
        /*species 34: aC3H5 */
        species[34] =
            +6.50078770e+00
            +1.43247310e-02 * tc[1]
            -5.67816320e-06 * tc[2]
            +1.10808010e-09 * tc[3]
            -9.03638870e-14 * tc[4];
        /*species 35: CH3CCH2 */
        species[35] =
            +5.42555280e+00
            +1.55110720e-02 * tc[1]
            -5.66783500e-06 * tc[2]
            +7.92243880e-10 * tc[3]
            -1.68780340e-14 * tc[4];
        /*species 36: C3H6 */
        species[36] =
            +6.73225700e+00
            +1.49083400e-02 * tc[1]
            -4.94989900e-06 * tc[2]
            +7.21202200e-10 * tc[3]
            -3.76620400e-14 * tc[4];
        /*species 37: nC3H7 */
        species[37] =
            +7.70974790e+00
            +1.60314850e-02 * tc[1]
            -5.27202380e-06 * tc[2]
            +7.58883520e-10 * tc[3]
            -3.88627190e-14 * tc[4];
        /*species 38: iC3H7 */
        species[38] =
            +6.51927410e+00
            +1.72201040e-02 * tc[1]
            -5.73642170e-06 * tc[2]
            +8.41307320e-10 * tc[3]
            -4.45659130e-14 * tc[4];
        /*species 39: C2H3CHO */
        species[39] =
            +5.81118680e+00
            +1.71142560e-02 * tc[1]
            -7.48341610e-06 * tc[2]
            +1.42522490e-09 * tc[3]
            -9.17468410e-14 * tc[4];
        /*species 40: C4H2 */
        species[40] =
            +9.15763280e+00
            +5.54305180e-03 * tc[1]
            -1.35916040e-06 * tc[2]
            +1.87800750e-11 * tc[3]
            +2.31895360e-14 * tc[4];
        /*species 41: iC4H3 */
        species[41] =
            +7.65385480e+00
            +1.12040550e-02 * tc[1]
            -4.64013420e-06 * tc[2]
            +8.67866390e-10 * tc[3]
            -5.74305620e-14 * tc[4];
        /*species 42: C4H4 */
        species[42] =
            +7.25396010e+00
            +1.39140940e-02 * tc[1]
            -5.29322140e-06 * tc[2]
            +8.34804500e-10 * tc[3]
            -3.51978820e-14 * tc[4];
        /*species 43: iC4H5 */
        species[43] =
            +6.96460290e+00
            +1.82743330e-02 * tc[1]
            -7.81337350e-06 * tc[2]
            +1.52921540e-09 * tc[3]
            -1.09204930e-13 * tc[4];
        /*species 44: C4H5-2 */
        species[44] =
            +1.45381710e+01
            -8.56770560e-03 * tc[1]
            +2.35595240e-05 * tc[2]
            -1.36763790e-08 * tc[3]
            +2.44369270e-12 * tc[4];
        /*species 45: C4H6 */
        species[45] =
            +8.86731340e+00
            +1.49186700e-02 * tc[1]
            -3.15487160e-06 * tc[2]
            -4.18413300e-10 * tc[3]
            +1.57612580e-13 * tc[4];
        /*species 46: C4H612 */
        species[46] =
            +1.78155700e+01
            -4.25750200e-03 * tc[1]
            +1.05118500e-05 * tc[2]
            -4.47384400e-09 * tc[3]
            +5.84813800e-13 * tc[4];
        /*species 47: C4H6-2 */
        species[47] =
            +9.03381330e+00
            +8.21245100e-03 * tc[1]
            +7.17539520e-06 * tc[2]
            -5.88343340e-09 * tc[3]
            +1.03439150e-12 * tc[4];
        /*species 48: C4H7 */
        species[48] =
            +7.01348350e+00
            +2.26345580e-02 * tc[1]
            -9.25454700e-06 * tc[2]
            +1.68079270e-09 * tc[3]
            -1.04086170e-13 * tc[4];
        /*species 49: C4H81 */
        species[49] =
            +2.05358410e+00
            +3.43505070e-02 * tc[1]
            -1.58831970e-05 * tc[2]
            +3.30896620e-09 * tc[3]
            -2.53610450e-13 * tc[4];
        /*species 50: pC4H9 */
        species[50] =
            +8.68223950e+00
            +2.36910710e-02 * tc[1]
            -7.59488650e-06 * tc[2]
            +6.64271360e-10 * tc[3]
            +5.48451360e-14 * tc[4];
        /*species 55: N2 */
        species[55] =
            +2.92664000e+00
            +1.48797680e-03 * tc[1]
            -5.68476000e-07 * tc[2]
            +1.00970380e-10 * tc[3]
            -6.75335100e-15 * tc[4];
    }

    /*species with midpoint at T=1392 kelvin */
    if (T < 1392) {
        /*species 52: C6H12 */
        species[52] =
            -1.35275205e+00
            +6.98655426e-02 * tc[1]
            -4.59408022e-05 * tc[2]
            +1.56967343e-08 * tc[3]
            -2.21296175e-12 * tc[4];
        /*species 54: C5H10 */
        species[54] =
            -1.06223481e+00
            +5.74218294e-02 * tc[1]
            -3.74486890e-05 * tc[2]
            +1.27364989e-08 * tc[3]
            -1.79609789e-12 * tc[4];
    } else {
        /*species 52: C6H12 */
        species[52] =
            +1.78337529e+01
            +2.67377658e-02 * tc[1]
            -9.10036773e-06 * tc[2]
            +1.40819768e-09 * tc[3]
            -8.15124244e-14 * tc[4];
        /*species 54: C5H10 */
        species[54] =
            +1.45851539e+01
            +2.24072471e-02 * tc[1]
            -7.63348025e-06 * tc[2]
            +1.18188966e-09 * tc[3]
            -6.84385139e-14 * tc[4];
    }

    /*species with midpoint at T=1389 kelvin */
    if (T < 1389) {
        /*species 53: C6H11 */
        species[53] =
            -1.55544944e+00
            +6.76865602e-02 * tc[1]
            -4.47048635e-05 * tc[2]
            +1.52236630e-08 * tc[3]
            -2.14346377e-12 * tc[4];
    } else {
        /*species 53: C6H11 */
        species[53] =
            +1.77336550e+01
            +2.48934775e-02 * tc[1]
            -8.59991450e-06 * tc[2]
            +1.34412828e-09 * tc[3]
            -7.83475666e-14 * tc[4];
    }

    /*species with midpoint at T=1391 kelvin */
    if (T < 1391) {
        /*species 51: NC12H26 */
        species[51] =
            -2.62181594e+00
            +1.47237711e-01 * tc[1]
            -9.43970271e-05 * tc[2]
            +3.07441268e-08 * tc[3]
            -4.03602230e-12 * tc[4];
    } else {
        /*species 51: NC12H26 */
        species[51] =
            +3.85095037e+01
            +5.63550048e-02 * tc[1]
            -1.91493200e-05 * tc[2]
            +2.96024862e-09 * tc[3]
            -1.71244150e-13 * tc[4];
    }
    return;
}


/*compute the e/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
__device__ void speciesInternalEnergy_d(double * species, double *  tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +1.50000000e+00
            +3.52666409e-13 * tc[1]
            -6.65306547e-16 * tc[2]
            +5.75204080e-19 * tc[3]
            -1.85546466e-22 * tc[4]
            +2.54736599e+04 * invT;
        /*species 1: O */
        species[1] =
            +2.16826710e+00
            -1.63965942e-03 * tc[1]
            +2.21435465e-06 * tc[2]
            -1.53201656e-09 * tc[3]
            +4.22531942e-13 * tc[4]
            +2.91222592e+04 * invT;
        /*species 2: OH */
        species[2] =
            +3.12530561e+00
            -1.61272470e-03 * tc[1]
            +2.17588230e-06 * tc[2]
            -1.44963411e-09 * tc[3]
            +4.12474758e-13 * tc[4]
            +3.38153812e+03 * invT;
        /*species 3: HO2 */
        species[3] =
            +3.30179801e+00
            -2.37456025e-03 * tc[1]
            +7.05276303e-06 * tc[2]
            -6.06909735e-09 * tc[3]
            +1.85845025e-12 * tc[4]
            +2.94808040e+02 * invT;
        /*species 4: H2 */
        species[4] =
            +1.34433112e+00
            +3.99026037e-03 * tc[1]
            -6.49271700e-06 * tc[2]
            +5.03930235e-09 * tc[3]
            -1.47522352e-12 * tc[4]
            -9.17935173e+02 * invT;
        /*species 5: H2O */
        species[5] =
            +3.19864056e+00
            -1.01821705e-03 * tc[1]
            +2.17346737e-06 * tc[2]
            -1.37199266e-09 * tc[3]
            +3.54395634e-13 * tc[4]
            -3.02937267e+04 * invT;
        /*species 6: H2O2 */
        species[6] =
            +3.27611269e+00
            -2.71411208e-04 * tc[1]
            +5.57785670e-06 * tc[2]
            -5.39427032e-09 * tc[3]
            +1.72490873e-12 * tc[4]
            -1.77025821e+04 * invT;
        /*species 7: O2 */
        species[7] =
            +2.78245636e+00
            -1.49836708e-03 * tc[1]
            +3.28243400e-06 * tc[2]
            -2.42032377e-09 * tc[3]
            +6.48745674e-13 * tc[4]
            -1.06394356e+03 * invT;
        /*species 8: CH */
        species[8] =
            +2.48981665e+00
            +1.61917771e-04 * tc[1]
            -5.62996883e-07 * tc[2]
            +7.90543317e-10 * tc[3]
            -2.81218134e-13 * tc[4]
            +7.07972934e+04 * invT;
        /*species 9: CH2 */
        species[9] =
            +2.76267867e+00
            +4.84436072e-04 * tc[1]
            +9.31632803e-07 * tc[2]
            -9.62727883e-10 * tc[3]
            +3.37483438e-13 * tc[4]
            +4.60040401e+04 * invT;
        /*species 10: CH2* */
        species[10] =
            +3.19860411e+00
            -1.18330710e-03 * tc[1]
            +2.74432073e-06 * tc[2]
            -1.67203995e-09 * tc[3]
            +3.88629474e-13 * tc[4]
            +5.04968163e+04 * invT;
        /*species 11: CH3 */
        species[11] =
            +2.67359040e+00
            +1.00547588e-03 * tc[1]
            +1.91007285e-06 * tc[2]
            -1.71779356e-09 * tc[3]
            +5.08771468e-13 * tc[4]
            +1.64449988e+04 * invT;
        /*species 12: CH4 */
        species[12] =
            +4.14987613e+00
            -6.83548940e-03 * tc[1]
            +1.63933533e-05 * tc[2]
            -1.21185757e-08 * tc[3]
            +3.33387912e-12 * tc[4]
            -1.02466476e+04 * invT;
        /*species 13: HCO */
        species[13] =
            +3.22118584e+00
            -1.62196266e-03 * tc[1]
            +4.59331487e-06 * tc[2]
            -3.32860233e-09 * tc[3]
            +8.67537730e-13 * tc[4]
            +3.83956496e+03 * invT;
        /*species 14: CH2O */
        species[14] =
            +3.79372315e+00
            -4.95416684e-03 * tc[1]
            +1.24406669e-05 * tc[2]
            -9.48213152e-09 * tc[3]
            +2.63545304e-12 * tc[4]
            -1.43089567e+04 * invT;
        /*species 15: CH3O */
        species[15] =
            +2.71180502e+00
            -1.40231653e-03 * tc[1]
            +1.25516990e-05 * tc[2]
            -1.18268022e-08 * tc[3]
            +3.73176840e-12 * tc[4]
            +1.29569760e+03 * invT;
        /*species 16: CH2OH */
        species[16] =
            +3.47834367e+00
            -6.75351550e-04 * tc[1]
            +9.28283267e-06 * tc[2]
            -9.12172650e-09 * tc[3]
            +2.95814900e-12 * tc[4]
            -3.50072890e+03 * invT;
        /*species 17: CH3OH */
        species[17] =
            +4.71539582e+00
            -7.61545645e-03 * tc[1]
            +2.17480385e-05 * tc[2]
            -1.77701722e-08 * tc[3]
            +5.22705396e-12 * tc[4]
            -2.56427656e+04 * invT;
        /*species 18: CO */
        species[18] =
            +2.57953347e+00
            -3.05176840e-04 * tc[1]
            +3.38938110e-07 * tc[2]
            +2.26751471e-10 * tc[3]
            -1.80884900e-13 * tc[4]
            -1.43440860e+04 * invT;
        /*species 19: CO2 */
        species[19] =
            +1.35677352e+00
            +4.49229839e-03 * tc[1]
            -2.37452090e-06 * tc[2]
            +6.14797555e-10 * tc[3]
            -2.87399096e-14 * tc[4]
            -4.83719697e+04 * invT;
        /*species 20: C2H */
        species[20] =
            +1.88965733e+00
            +6.70498055e-03 * tc[1]
            -9.49231670e-06 * tc[2]
            +7.36977613e-09 * tc[3]
            -2.18663022e-12 * tc[4]
            +6.68393932e+04 * invT;
        /*species 21: C2H2 */
        species[21] =
            -1.91318906e-01
            +1.16807815e-02 * tc[1]
            -1.18390605e-05 * tc[2]
            +7.00381092e-09 * tc[3]
            -1.70014595e-12 * tc[4]
            +2.64289807e+04 * invT;
        /*species 22: C2H3 */
        species[22] =
            +2.21246645e+00
            +7.57395810e-04 * tc[1]
            +8.64031373e-06 * tc[2]
            -8.94144617e-09 * tc[3]
            +2.94301746e-12 * tc[4]
            +3.48598468e+04 * invT;
        /*species 23: C2H4 */
        species[23] =
            +2.95920148e+00
            -3.78526124e-03 * tc[1]
            +1.90330097e-05 * tc[2]
            -1.72897188e-08 * tc[3]
            +5.39768746e-12 * tc[4]
            +5.08977593e+03 * invT;
        /*species 24: C2H5 */
        species[24] =
            +3.30646568e+00
            -2.09329446e-03 * tc[1]
            +1.65714269e-05 * tc[2]
            -1.49781651e-08 * tc[3]
            +4.61018008e-12 * tc[4]
            +1.28416265e+04 * invT;
        /*species 25: C2H6 */
        species[25] =
            +3.29142492e+00
            -2.75077135e-03 * tc[1]
            +1.99812763e-05 * tc[2]
            -1.77116571e-08 * tc[3]
            +5.37371542e-12 * tc[4]
            -1.15222055e+04 * invT;
        /*species 26: HCCO */
        species[26] =
            +1.25172140e+00
            +8.82751050e-03 * tc[1]
            -7.90970033e-06 * tc[2]
            +4.31893975e-09 * tc[3]
            -1.01329622e-12 * tc[4]
            +2.00594490e+04 * invT;
        /*species 27: CH2CO */
        species[27] =
            +1.13583630e+00
            +9.05943605e-03 * tc[1]
            -5.79824913e-06 * tc[2]
            +2.33599392e-09 * tc[3]
            -4.02915230e-13 * tc[4]
            -7.27000000e+03 * invT;
        /*species 28: CH3CO */
        species[28] =
            +3.16342570e+00
            -1.16308050e-04 * tc[1]
            +1.14226067e-05 * tc[2]
            -1.10263067e-08 * tc[3]
            +3.45512240e-12 * tc[4]
            -2.65745290e+03 * invT;
        /*species 29: CH2CHO */
        species[29] =
            +2.40906240e+00
            +5.36928700e-03 * tc[1]
            +6.30497500e-07 * tc[2]
            -1.78964578e-09 * tc[3]
            +5.73477020e-13 * tc[4]
            +6.20000000e+01 * invT;
        /*species 30: CH3CHO */
        species[30] =
            +3.72945950e+00
            -1.59664290e-03 * tc[1]
            +1.58449737e-05 * tc[2]
            -1.43646527e-08 * tc[3]
            +4.38622240e-12 * tc[4]
            -2.15728780e+04 * invT;
        /*species 31: C3H3 */
        species[31] =
            +3.51109270e-01
            +1.63705612e-02 * tc[1]
            -1.57942378e-05 * tc[2]
            +9.40774520e-09 * tc[3]
            -2.37081846e-12 * tc[4]
            +4.01057783e+04 * invT;
        /*species 32: pC3H4 */
        species[32] =
            +1.68038690e+00
            +7.89982550e-03 * tc[1]
            +8.35686533e-07 * tc[2]
            -3.41440575e-09 * tc[3]
            +1.32308570e-12 * tc[4]
            +2.08023740e+04 * invT;
        /*species 33: aC3H4 */
        species[33] =
            +1.61304450e+00
            +6.06128750e-03 * tc[1]
            +6.17996000e-06 * tc[2]
            -8.63128725e-09 * tc[3]
            +3.06701580e-12 * tc[4]
            +2.15415670e+04 * invT;
        /*species 34: aC3H5 */
        species[34] =
            +3.63183500e-01
            +9.90691050e-03 * tc[1]
            +4.16568667e-06 * tc[2]
            -8.33888875e-09 * tc[3]
            +3.16931420e-12 * tc[4]
            +1.92456290e+04 * invT;
        /*species 35: CH3CCH2 */
        species[35] =
            +7.32920900e-01
            +1.11973100e-02 * tc[1]
            -1.71635370e-06 * tc[2]
            -1.68991165e-09 * tc[3]
            +7.65064220e-13 * tc[4]
            +2.90404980e+04 * invT;
        /*species 36: C3H6 */
        species[36] =
            +4.93307000e-01
            +1.04625900e-02 * tc[1]
            +1.49559800e-06 * tc[2]
            -4.17228000e-09 * tc[3]
            +1.43162920e-12 * tc[4]
            +1.07482600e+03 * invT;
        /*species 37: nC3H7 */
        species[37] =
            +4.91173000e-02
            +1.30044865e-02 * tc[1]
            +7.84750533e-07 * tc[2]
            -4.89878300e-09 * tc[3]
            +1.87440414e-12 * tc[4]
            +1.03123460e+04 * invT;
        /*species 38: iC3H7 */
        species[38] =
            +4.44919900e-01
            +1.04995560e-02 * tc[1]
            +2.56787407e-06 * tc[2]
            -4.61906325e-09 * tc[3]
            +1.42565924e-12 * tc[4]
            +9.42237240e+03 * invT;
        /*species 39: C2H3CHO */
        species[39] =
            +2.71349800e-01
            +1.31155270e-02 * tc[1]
            -3.09707683e-06 * tc[2]
            -1.19593180e-09 * tc[3]
            +6.69610860e-13 * tc[4]
            -9.33573440e+03 * invT;
        /*species 40: C4H2 */
        species[40] =
            +5.43978000e-02
            +2.08134800e-02 * tc[1]
            -2.19572613e-05 * tc[2]
            +1.33142687e-08 * tc[3]
            -3.33663240e-12 * tc[4]
            +5.41852110e+04 * invT;
        /*species 41: iC4H3 */
        species[41] =
            +2.72214820e+00
            +1.29787715e-02 * tc[1]
            -8.78544767e-06 * tc[2]
            +3.87723000e-09 * tc[3]
            -7.60811300e-13 * tc[4]
            +5.88371210e+04 * invT;
        /*species 42: C4H4 */
        species[42] =
            -4.11429520e-01
            +1.82733425e-02 * tc[1]
            -1.13689893e-05 * tc[2]
            +4.16315475e-09 * tc[3]
            -6.01292460e-13 * tc[4]
            +3.33594920e+04 * invT;
        /*species 43: iC4H5 */
        species[43] =
            -8.86918950e-01
            +2.04753075e-02 * tc[1]
            -1.18045270e-05 * tc[2]
            +3.88274225e-09 * tc[3]
            -4.67102440e-13 * tc[4]
            +3.63833710e+04 * invT;
        /*species 44: C4H5-2 */
        species[44] =
            +1.96962800e+00
            +1.22211225e-02 * tc[1]
            -3.04171413e-06 * tc[2]
            -1.06167177e-18 * tc[3]
            +3.26094560e-22 * tc[4]
            +3.55033160e+04 * invT;
        /*species 45: C4H6 */
        species[45] =
            -8.87155350e-01
            +1.71845110e-02 * tc[1]
            -3.70246400e-06 * tc[2]
            -2.30266650e-09 * tc[3]
            +1.24130358e-12 * tc[4]
            +1.18022700e+04 * invT;
        /*species 46: C4H612 */
        species[46] =
            +2.34670000e-02
            +1.74795950e-02 * tc[1]
            -7.33635000e-06 * tc[2]
            +1.73556800e-09 * tc[3]
            -1.57583740e-13 * tc[4]
            +1.81179900e+04 * invT;
        /*species 47: C4H6-2 */
        species[47] =
            +1.13733380e+00
            +1.32431145e-02 * tc[1]
            -3.01895703e-06 * tc[2]
            -1.38465992e-19 * tc[3]
            +4.25637680e-23 * tc[4]
            +1.57109020e+04 * invT;
        /*species 48: C4H7 */
        species[48] =
            -2.55505680e-01
            +1.98394285e-02 * tc[1]
            -7.63269533e-06 * tc[2]
            +5.33824325e-10 * tc[3]
            +4.61927500e-13 * tc[4]
            +2.26533280e+04 * invT;
        /*species 49: C4H81 */
        species[49] =
            +1.81138000e-01
            +1.54266900e-02 * tc[1]
            +1.69550823e-06 * tc[2]
            -6.16372200e-09 * tc[3]
            +2.22203860e-12 * tc[4]
            -1.79040040e+03 * invT;
        /*species 50: pC4H9 */
        species[50] =
            +2.08704200e-01
            +1.91487485e-02 * tc[1]
            -2.42201697e-06 * tc[2]
            -3.85713675e-09 * tc[3]
            +1.73718870e-12 * tc[4]
            +7.32210400e+03 * invT;
        /*species 55: N2 */
        species[55] =
            +2.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 * invT;
    } else {
        /*species 0: H */
        species[0] =
            +1.50000001e+00
            -1.15421486e-11 * tc[1]
            +5.38539827e-15 * tc[2]
            -1.18378809e-18 * tc[3]
            +9.96394714e-23 * tc[4]
            +2.54736599e+04 * invT;
        /*species 1: O */
        species[1] =
            +1.56942078e+00
            -4.29870569e-05 * tc[1]
            +1.39828196e-08 * tc[2]
            -2.50444497e-12 * tc[3]
            +2.45667382e-16 * tc[4]
            +2.92175791e+04 * invT;
        /*species 2: OH */
        species[2] =
            +1.86472886e+00
            +5.28252240e-04 * tc[1]
            -8.63609193e-08 * tc[2]
            +7.63046685e-12 * tc[3]
            -2.66391752e-16 * tc[4]
            +3.71885774e+03 * invT;
        /*species 3: HO2 */
        species[3] =
            +3.01721090e+00
            +1.11991006e-03 * tc[1]
            -2.11219383e-07 * tc[2]
            +2.85615925e-11 * tc[3]
            -2.15817070e-15 * tc[4]
            +1.11856713e+02 * invT;
        /*species 4: H2 */
        species[4] =
            +2.33727920e+00
            -2.47012365e-05 * tc[1]
            +1.66485593e-07 * tc[2]
            -4.48915985e-11 * tc[3]
            +4.00510752e-15 * tc[4]
            -9.50158922e+02 * invT;
        /*species 5: H2O */
        species[5] =
            +2.03399249e+00
            +1.08845902e-03 * tc[1]
            -5.46908393e-08 * tc[2]
            -2.42604967e-11 * tc[3]
            +3.36401984e-15 * tc[4]
            -3.00042971e+04 * invT;
        /*species 6: H2O2 */
        species[6] =
            +3.16500285e+00
            +2.45415847e-03 * tc[1]
            -6.33797417e-07 * tc[2]
            +9.27964965e-11 * tc[3]
            -5.75816610e-15 * tc[4]
            -1.78617877e+04 * invT;
        /*species 7: O2 */
        species[7] =
            +2.28253784e+00
            +7.41543770e-04 * tc[1]
            -2.52655556e-07 * tc[2]
            +5.23676387e-11 * tc[3]
            -4.33435588e-15 * tc[4]
            -1.08845772e+03 * invT;
        /*species 8: CH */
        species[8] =
            +1.87846473e+00
            +4.85456840e-04 * tc[1]
            +4.81485517e-08 * tc[2]
            -3.26719623e-11 * tc[3]
            +3.52158766e-15 * tc[4]
            +7.10124364e+04 * invT;
        /*species 9: CH2 */
        species[9] =
            +1.87410113e+00
            +1.82819646e-03 * tc[1]
            -4.69648657e-07 * tc[2]
            +6.50448872e-11 * tc[3]
            -3.75455134e-15 * tc[4]
            +4.62636040e+04 * invT;
        /*species 10: CH2* */
        species[10] =
            +1.29203842e+00
            +2.32794318e-03 * tc[1]
            -6.70639823e-07 * tc[2]
            +1.04476500e-10 * tc[3]
            -6.79432730e-15 * tc[4]
            +5.09259997e+04 * invT;
        /*species 11: CH3 */
        species[11] =
            +1.28571772e+00
            +3.61995018e-03 * tc[1]
            -9.95714493e-07 * tc[2]
            +1.48921161e-10 * tc[3]
            -9.34308788e-15 * tc[4]
            +1.67755843e+04 * invT;
        /*species 12: CH4 */
        species[12] =
            -9.25148505e-01
            +6.69547335e-03 * tc[1]
            -1.91095270e-06 * tc[2]
            +3.05731338e-10 * tc[3]
            -2.03630460e-14 * tc[4]
            -9.46834459e+03 * invT;
        /*species 13: HCO */
        species[13] =
            +1.77217438e+00
            +2.47847763e-03 * tc[1]
            -8.28152043e-07 * tc[2]
            +1.47290445e-10 * tc[3]
            -1.06701742e-14 * tc[4]
            +4.01191815e+03 * invT;
        /*species 14: CH2O */
        species[14] =
            +7.60690080e-01
            +4.60000041e-03 * tc[1]
            -1.47419604e-06 * tc[2]
            +2.51603030e-10 * tc[3]
            -1.76771128e-14 * tc[4]
            -1.39958323e+04 * invT;
        /*species 15: CH3O */
        species[15] =
            +3.75779238e+00
            +3.72071237e-03 * tc[1]
            -8.99017253e-07 * tc[2]
            +1.09522626e-10 * tc[3]
            -5.27074196e-15 * tc[4]
            +3.78111940e+02 * invT;
        /*species 16: CH2OH */
        species[16] =
            +4.09314370e+00
            +2.97380630e-03 * tc[1]
            -6.88324867e-07 * tc[2]
            +8.07520433e-11 * tc[3]
            -3.76251804e-15 * tc[4]
            -4.03409640e+03 * invT;
        /*species 17: CH3OH */
        species[17] =
            +7.89707910e-01
            +7.04691460e-03 * tc[1]
            -2.12166945e-06 * tc[2]
            +3.45427713e-10 * tc[3]
            -2.34120440e-14 * tc[4]
            -2.53748747e+04 * invT;
        /*species 18: CO */
        species[18] =
            +1.71518561e+00
            +1.03126372e-03 * tc[1]
            -3.32941924e-07 * tc[2]
            +5.75132520e-11 * tc[3]
            -4.07295432e-15 * tc[4]
            -1.41518724e+04 * invT;
        /*species 19: CO2 */
        species[19] =
            +2.85746029e+00
            +2.20718513e-03 * tc[1]
            -7.38271347e-07 * tc[2]
            +1.30872547e-10 * tc[3]
            -9.44168328e-15 * tc[4]
            -4.87591660e+04 * invT;
        /*species 20: C2H */
        species[20] =
            +2.16780652e+00
            +2.37610951e-03 * tc[1]
            -6.12623590e-07 * tc[2]
            +7.60475630e-11 * tc[3]
            -3.54465540e-15 * tc[4]
            +6.71210650e+04 * invT;
        /*species 21: C2H2 */
        species[21] =
            +3.14756964e+00
            +2.98083332e-03 * tc[1]
            -7.90982840e-07 * tc[2]
            +1.16853043e-10 * tc[3]
            -7.22470426e-15 * tc[4]
            +2.59359992e+04 * invT;
        /*species 22: C2H3 */
        species[22] =
            +2.01672400e+00
            +5.16511460e-03 * tc[1]
            -1.56027450e-06 * tc[2]
            +2.54408220e-10 * tc[3]
            -1.72521408e-14 * tc[4]
            +3.46128739e+04 * invT;
        /*species 23: C2H4 */
        species[23] =
            +1.03611116e+00
            +7.32270755e-03 * tc[1]
            -2.23692638e-06 * tc[2]
            +3.68057308e-10 * tc[3]
            -2.51412122e-14 * tc[4]
            +4.93988614e+03 * invT;
        /*species 24: C2H5 */
        species[24] =
            +9.54656420e-01
            +8.69863610e-03 * tc[1]
            -2.66068889e-06 * tc[2]
            +4.38044223e-10 * tc[3]
            -2.99283152e-14 * tc[4]
            +1.28575200e+04 * invT;
        /*species 25: C2H6 */
        species[25] =
            +7.18815000e-02
            +1.08426339e-02 * tc[1]
            -3.34186890e-06 * tc[2]
            +5.53530003e-10 * tc[3]
            -3.80005780e-14 * tc[4]
            -1.14263932e+04 * invT;
        /*species 26: HCCO */
        species[26] =
            +4.62820580e+00
            +2.04267005e-03 * tc[1]
            -5.31151567e-07 * tc[2]
            +7.15651300e-11 * tc[3]
            -3.88156640e-15 * tc[4]
            +1.93272150e+04 * invT;
        /*species 27: CH2CO */
        species[27] =
            +3.51129732e+00
            +4.50179872e-03 * tc[1]
            -1.38979878e-06 * tc[2]
            +2.30836470e-10 * tc[3]
            -1.58967640e-14 * tc[4]
            -7.77850000e+03 * invT;
        /*species 28: CH3CO */
        species[28] =
            +4.94477310e+00
            +3.93336025e-03 * tc[1]
            -9.62196067e-07 * tc[2]
            +1.18177188e-10 * tc[3]
            -5.71997220e-15 * tc[4]
            -3.78730750e+03 * invT;
        /*species 29: CH2CHO */
        species[29] =
            +4.97566990e+00
            +4.06529570e-03 * tc[1]
            -9.14541500e-07 * tc[2]
            +1.01757603e-10 * tc[3]
            -4.35203420e-15 * tc[4]
            -9.69500000e+02 * invT;
        /*species 30: CH3CHO */
        species[30] =
            +4.40411080e+00
            +5.86152950e-03 * tc[1]
            -1.40877123e-06 * tc[2]
            +1.70931128e-10 * tc[3]
            -8.19697260e-15 * tc[4]
            -2.25931220e+04 * invT;
        /*species 31: C3H3 */
        species[31] =
            +6.14221880e+00
            +3.80951002e-03 * tc[1]
            -8.91533167e-07 * tc[2]
            +1.06228700e-10 * tc[3]
            -5.02950830e-15 * tc[4]
            +3.89087427e+04 * invT;
        /*species 32: pC3H4 */
        species[32] =
            +5.02524000e+00
            +5.66827100e-03 * tc[1]
            -1.34077970e-06 * tc[2]
            +1.60940157e-10 * tc[3]
            -7.65992700e-15 * tc[4]
            +1.96209420e+04 * invT;
        /*species 33: aC3H4 */
        species[33] =
            +5.31687220e+00
            +5.56686400e-03 * tc[1]
            -1.32097927e-06 * tc[2]
            +1.58910595e-10 * tc[3]
            -7.57510800e-15 * tc[4]
            +2.01174950e+04 * invT;
        /*species 34: aC3H5 */
        species[34] =
            +5.50078770e+00
            +7.16236550e-03 * tc[1]
            -1.89272107e-06 * tc[2]
            +2.77020025e-10 * tc[3]
            -1.80727774e-14 * tc[4]
            +1.74824490e+04 * invT;
        /*species 35: CH3CCH2 */
        species[35] =
            +4.42555280e+00
            +7.75553600e-03 * tc[1]
            -1.88927833e-06 * tc[2]
            +1.98060970e-10 * tc[3]
            -3.37560680e-15 * tc[4]
            +2.78430270e+04 * invT;
        /*species 36: C3H6 */
        species[36] =
            +5.73225700e+00
            +7.45417000e-03 * tc[1]
            -1.64996633e-06 * tc[2]
            +1.80300550e-10 * tc[3]
            -7.53240800e-15 * tc[4]
            -9.23570300e+02 * invT;
        /*species 37: nC3H7 */
        species[37] =
            +6.70974790e+00
            +8.01574250e-03 * tc[1]
            -1.75734127e-06 * tc[2]
            +1.89720880e-10 * tc[3]
            -7.77254380e-15 * tc[4]
            +7.97622360e+03 * invT;
        /*species 38: iC3H7 */
        species[38] =
            +5.51927410e+00
            +8.61005200e-03 * tc[1]
            -1.91214057e-06 * tc[2]
            +2.10326830e-10 * tc[3]
            -8.91318260e-15 * tc[4]
            +7.32271930e+03 * invT;
        /*species 39: C2H3CHO */
        species[39] =
            +4.81118680e+00
            +8.55712800e-03 * tc[1]
            -2.49447203e-06 * tc[2]
            +3.56306225e-10 * tc[3]
            -1.83493682e-14 * tc[4]
            -1.07840540e+04 * invT;
        /*species 40: C4H2 */
        species[40] =
            +8.15763280e+00
            +2.77152590e-03 * tc[1]
            -4.53053467e-07 * tc[2]
            +4.69501875e-12 * tc[3]
            +4.63790720e-15 * tc[4]
            +5.25880390e+04 * invT;
        /*species 41: iC4H3 */
        species[41] =
            +6.65385480e+00
            +5.60202750e-03 * tc[1]
            -1.54671140e-06 * tc[2]
            +2.16966597e-10 * tc[3]
            -1.14861124e-14 * tc[4]
            +5.79543630e+04 * invT;
        /*species 42: C4H4 */
        species[42] =
            +6.25396010e+00
            +6.95704700e-03 * tc[1]
            -1.76440713e-06 * tc[2]
            +2.08701125e-10 * tc[3]
            -7.03957640e-15 * tc[4]
            +3.17660160e+04 * invT;
        /*species 43: iC4H5 */
        species[43] =
            +5.96460290e+00
            +9.13716650e-03 * tc[1]
            -2.60445783e-06 * tc[2]
            +3.82303850e-10 * tc[3]
            -2.18409860e-14 * tc[4]
            +3.47250980e+04 * invT;
        /*species 44: C4H5-2 */
        species[44] =
            +1.35381710e+01
            -4.28385280e-03 * tc[1]
            +7.85317467e-06 * tc[2]
            -3.41909475e-09 * tc[3]
            +4.88738540e-13 * tc[4]
            +3.32590950e+04 * invT;
        /*species 45: C4H6 */
        species[45] =
            +7.86731340e+00
            +7.45933500e-03 * tc[1]
            -1.05162387e-06 * tc[2]
            -1.04603325e-10 * tc[3]
            +3.15225160e-14 * tc[4]
            +9.13385160e+03 * invT;
        /*species 46: C4H612 */
        species[46] =
            +1.68155700e+01
            -2.12875100e-03 * tc[1]
            +3.50395000e-06 * tc[2]
            -1.11846100e-09 * tc[3]
            +1.16962760e-13 * tc[4]
            +1.26734200e+04 * invT;
        /*species 47: C4H6-2 */
        species[47] =
            +8.03381330e+00
            +4.10622550e-03 * tc[1]
            +2.39179840e-06 * tc[2]
            -1.47085835e-09 * tc[3]
            +2.06878300e-13 * tc[4]
            +1.43350680e+04 * invT;
        /*species 48: C4H7 */
        species[48] =
            +6.01348350e+00
            +1.13172790e-02 * tc[1]
            -3.08484900e-06 * tc[2]
            +4.20198175e-10 * tc[3]
            -2.08172340e-14 * tc[4]
            +2.09550080e+04 * invT;
        /*species 49: C4H81 */
        species[49] =
            +1.05358410e+00
            +1.71752535e-02 * tc[1]
            -5.29439900e-06 * tc[2]
            +8.27241550e-10 * tc[3]
            -5.07220900e-14 * tc[4]
            -2.13972310e+03 * invT;
        /*species 50: pC4H9 */
        species[50] =
            +7.68223950e+00
            +1.18455355e-02 * tc[1]
            -2.53162883e-06 * tc[2]
            +1.66067840e-10 * tc[3]
            +1.09690272e-14 * tc[4]
            +4.96440580e+03 * invT;
        /*species 55: N2 */
        species[55] =
            +1.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
    }

    /*species with midpoint at T=1392 kelvin */
    if (T < 1392) {
        /*species 52: C6H12 */
        species[52] =
            -2.35275205e+00
            +3.49327713e-02 * tc[1]
            -1.53136007e-05 * tc[2]
            +3.92418358e-09 * tc[3]
            -4.42592350e-13 * tc[4]
            -7.34368617e+03 * invT;
        /*species 54: C5H10 */
        species[54] =
            -2.06223481e+00
            +2.87109147e-02 * tc[1]
            -1.24828963e-05 * tc[2]
            +3.18412472e-09 * tc[3]
            -3.59219578e-13 * tc[4]
            -4.46546666e+03 * invT;
    } else {
        /*species 52: C6H12 */
        species[52] =
            +1.68337529e+01
            +1.33688829e-02 * tc[1]
            -3.03345591e-06 * tc[2]
            +3.52049420e-10 * tc[3]
            -1.63024849e-14 * tc[4]
            -1.42062860e+04 * invT;
        /*species 54: C5H10 */
        species[54] =
            +1.35851539e+01
            +1.12036235e-02 * tc[1]
            -2.54449342e-06 * tc[2]
            +2.95472415e-10 * tc[3]
            -1.36877028e-14 * tc[4]
            -1.00898205e+04 * invT;
    }

    /*species with midpoint at T=1389 kelvin */
    if (T < 1389) {
        /*species 53: C6H11 */
        species[53] =
            -2.55544944e+00
            +3.38432801e-02 * tc[1]
            -1.49016212e-05 * tc[2]
            +3.80591575e-09 * tc[3]
            -4.28692754e-13 * tc[4]
            +9.66316695e+03 * invT;
    } else {
        /*species 53: C6H11 */
        species[53] =
            +1.67336550e+01
            +1.24467388e-02 * tc[1]
            -2.86663817e-06 * tc[2]
            +3.36032070e-10 * tc[3]
            -1.56695133e-14 * tc[4]
            +2.68017174e+03 * invT;
    }

    /*species with midpoint at T=1391 kelvin */
    if (T < 1391) {
        /*species 51: NC12H26 */
        species[51] =
            -3.62181594e+00
            +7.36188555e-02 * tc[1]
            -3.14656757e-05 * tc[2]
            +7.68603170e-09 * tc[3]
            -8.07204460e-13 * tc[4]
            -4.00654253e+04 * invT;
    } else {
        /*species 51: NC12H26 */
        species[51] =
            +3.75095037e+01
            +2.81775024e-02 * tc[1]
            -6.38310667e-06 * tc[2]
            +7.40062155e-10 * tc[3]
            -3.42488300e-14 * tc[4]
            -5.48843465e+04 * invT;
    }
    return;
}


/*compute the h/(RT) at the given temperature (Eq 20) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
__device__ void speciesEnthalpy_d(double * species, double *  tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +2.50000000e+00
            +3.52666409e-13 * tc[1]
            -6.65306547e-16 * tc[2]
            +5.75204080e-19 * tc[3]
            -1.85546466e-22 * tc[4]
            +2.54736599e+04 * invT;
        /*species 1: O */
        species[1] =
            +3.16826710e+00
            -1.63965942e-03 * tc[1]
            +2.21435465e-06 * tc[2]
            -1.53201656e-09 * tc[3]
            +4.22531942e-13 * tc[4]
            +2.91222592e+04 * invT;
        /*species 2: OH */
        species[2] =
            +4.12530561e+00
            -1.61272470e-03 * tc[1]
            +2.17588230e-06 * tc[2]
            -1.44963411e-09 * tc[3]
            +4.12474758e-13 * tc[4]
            +3.38153812e+03 * invT;
        /*species 3: HO2 */
        species[3] =
            +4.30179801e+00
            -2.37456025e-03 * tc[1]
            +7.05276303e-06 * tc[2]
            -6.06909735e-09 * tc[3]
            +1.85845025e-12 * tc[4]
            +2.94808040e+02 * invT;
        /*species 4: H2 */
        species[4] =
            +2.34433112e+00
            +3.99026037e-03 * tc[1]
            -6.49271700e-06 * tc[2]
            +5.03930235e-09 * tc[3]
            -1.47522352e-12 * tc[4]
            -9.17935173e+02 * invT;
        /*species 5: H2O */
        species[5] =
            +4.19864056e+00
            -1.01821705e-03 * tc[1]
            +2.17346737e-06 * tc[2]
            -1.37199266e-09 * tc[3]
            +3.54395634e-13 * tc[4]
            -3.02937267e+04 * invT;
        /*species 6: H2O2 */
        species[6] =
            +4.27611269e+00
            -2.71411208e-04 * tc[1]
            +5.57785670e-06 * tc[2]
            -5.39427032e-09 * tc[3]
            +1.72490873e-12 * tc[4]
            -1.77025821e+04 * invT;
        /*species 7: O2 */
        species[7] =
            +3.78245636e+00
            -1.49836708e-03 * tc[1]
            +3.28243400e-06 * tc[2]
            -2.42032377e-09 * tc[3]
            +6.48745674e-13 * tc[4]
            -1.06394356e+03 * invT;
        /*species 8: CH */
        species[8] =
            +3.48981665e+00
            +1.61917771e-04 * tc[1]
            -5.62996883e-07 * tc[2]
            +7.90543317e-10 * tc[3]
            -2.81218134e-13 * tc[4]
            +7.07972934e+04 * invT;
        /*species 9: CH2 */
        species[9] =
            +3.76267867e+00
            +4.84436072e-04 * tc[1]
            +9.31632803e-07 * tc[2]
            -9.62727883e-10 * tc[3]
            +3.37483438e-13 * tc[4]
            +4.60040401e+04 * invT;
        /*species 10: CH2* */
        species[10] =
            +4.19860411e+00
            -1.18330710e-03 * tc[1]
            +2.74432073e-06 * tc[2]
            -1.67203995e-09 * tc[3]
            +3.88629474e-13 * tc[4]
            +5.04968163e+04 * invT;
        /*species 11: CH3 */
        species[11] =
            +3.67359040e+00
            +1.00547588e-03 * tc[1]
            +1.91007285e-06 * tc[2]
            -1.71779356e-09 * tc[3]
            +5.08771468e-13 * tc[4]
            +1.64449988e+04 * invT;
        /*species 12: CH4 */
        species[12] =
            +5.14987613e+00
            -6.83548940e-03 * tc[1]
            +1.63933533e-05 * tc[2]
            -1.21185757e-08 * tc[3]
            +3.33387912e-12 * tc[4]
            -1.02466476e+04 * invT;
        /*species 13: HCO */
        species[13] =
            +4.22118584e+00
            -1.62196266e-03 * tc[1]
            +4.59331487e-06 * tc[2]
            -3.32860233e-09 * tc[3]
            +8.67537730e-13 * tc[4]
            +3.83956496e+03 * invT;
        /*species 14: CH2O */
        species[14] =
            +4.79372315e+00
            -4.95416684e-03 * tc[1]
            +1.24406669e-05 * tc[2]
            -9.48213152e-09 * tc[3]
            +2.63545304e-12 * tc[4]
            -1.43089567e+04 * invT;
        /*species 15: CH3O */
        species[15] =
            +3.71180502e+00
            -1.40231653e-03 * tc[1]
            +1.25516990e-05 * tc[2]
            -1.18268022e-08 * tc[3]
            +3.73176840e-12 * tc[4]
            +1.29569760e+03 * invT;
        /*species 16: CH2OH */
        species[16] =
            +4.47834367e+00
            -6.75351550e-04 * tc[1]
            +9.28283267e-06 * tc[2]
            -9.12172650e-09 * tc[3]
            +2.95814900e-12 * tc[4]
            -3.50072890e+03 * invT;
        /*species 17: CH3OH */
        species[17] =
            +5.71539582e+00
            -7.61545645e-03 * tc[1]
            +2.17480385e-05 * tc[2]
            -1.77701722e-08 * tc[3]
            +5.22705396e-12 * tc[4]
            -2.56427656e+04 * invT;
        /*species 18: CO */
        species[18] =
            +3.57953347e+00
            -3.05176840e-04 * tc[1]
            +3.38938110e-07 * tc[2]
            +2.26751471e-10 * tc[3]
            -1.80884900e-13 * tc[4]
            -1.43440860e+04 * invT;
        /*species 19: CO2 */
        species[19] =
            +2.35677352e+00
            +4.49229839e-03 * tc[1]
            -2.37452090e-06 * tc[2]
            +6.14797555e-10 * tc[3]
            -2.87399096e-14 * tc[4]
            -4.83719697e+04 * invT;
        /*species 20: C2H */
        species[20] =
            +2.88965733e+00
            +6.70498055e-03 * tc[1]
            -9.49231670e-06 * tc[2]
            +7.36977613e-09 * tc[3]
            -2.18663022e-12 * tc[4]
            +6.68393932e+04 * invT;
        /*species 21: C2H2 */
        species[21] =
            +8.08681094e-01
            +1.16807815e-02 * tc[1]
            -1.18390605e-05 * tc[2]
            +7.00381092e-09 * tc[3]
            -1.70014595e-12 * tc[4]
            +2.64289807e+04 * invT;
        /*species 22: C2H3 */
        species[22] =
            +3.21246645e+00
            +7.57395810e-04 * tc[1]
            +8.64031373e-06 * tc[2]
            -8.94144617e-09 * tc[3]
            +2.94301746e-12 * tc[4]
            +3.48598468e+04 * invT;
        /*species 23: C2H4 */
        species[23] =
            +3.95920148e+00
            -3.78526124e-03 * tc[1]
            +1.90330097e-05 * tc[2]
            -1.72897188e-08 * tc[3]
            +5.39768746e-12 * tc[4]
            +5.08977593e+03 * invT;
        /*species 24: C2H5 */
        species[24] =
            +4.30646568e+00
            -2.09329446e-03 * tc[1]
            +1.65714269e-05 * tc[2]
            -1.49781651e-08 * tc[3]
            +4.61018008e-12 * tc[4]
            +1.28416265e+04 * invT;
        /*species 25: C2H6 */
        species[25] =
            +4.29142492e+00
            -2.75077135e-03 * tc[1]
            +1.99812763e-05 * tc[2]
            -1.77116571e-08 * tc[3]
            +5.37371542e-12 * tc[4]
            -1.15222055e+04 * invT;
        /*species 26: HCCO */
        species[26] =
            +2.25172140e+00
            +8.82751050e-03 * tc[1]
            -7.90970033e-06 * tc[2]
            +4.31893975e-09 * tc[3]
            -1.01329622e-12 * tc[4]
            +2.00594490e+04 * invT;
        /*species 27: CH2CO */
        species[27] =
            +2.13583630e+00
            +9.05943605e-03 * tc[1]
            -5.79824913e-06 * tc[2]
            +2.33599392e-09 * tc[3]
            -4.02915230e-13 * tc[4]
            -7.27000000e+03 * invT;
        /*species 28: CH3CO */
        species[28] =
            +4.16342570e+00
            -1.16308050e-04 * tc[1]
            +1.14226067e-05 * tc[2]
            -1.10263067e-08 * tc[3]
            +3.45512240e-12 * tc[4]
            -2.65745290e+03 * invT;
        /*species 29: CH2CHO */
        species[29] =
            +3.40906240e+00
            +5.36928700e-03 * tc[1]
            +6.30497500e-07 * tc[2]
            -1.78964578e-09 * tc[3]
            +5.73477020e-13 * tc[4]
            +6.20000000e+01 * invT;
        /*species 30: CH3CHO */
        species[30] =
            +4.72945950e+00
            -1.59664290e-03 * tc[1]
            +1.58449737e-05 * tc[2]
            -1.43646527e-08 * tc[3]
            +4.38622240e-12 * tc[4]
            -2.15728780e+04 * invT;
        /*species 31: C3H3 */
        species[31] =
            +1.35110927e+00
            +1.63705612e-02 * tc[1]
            -1.57942378e-05 * tc[2]
            +9.40774520e-09 * tc[3]
            -2.37081846e-12 * tc[4]
            +4.01057783e+04 * invT;
        /*species 32: pC3H4 */
        species[32] =
            +2.68038690e+00
            +7.89982550e-03 * tc[1]
            +8.35686533e-07 * tc[2]
            -3.41440575e-09 * tc[3]
            +1.32308570e-12 * tc[4]
            +2.08023740e+04 * invT;
        /*species 33: aC3H4 */
        species[33] =
            +2.61304450e+00
            +6.06128750e-03 * tc[1]
            +6.17996000e-06 * tc[2]
            -8.63128725e-09 * tc[3]
            +3.06701580e-12 * tc[4]
            +2.15415670e+04 * invT;
        /*species 34: aC3H5 */
        species[34] =
            +1.36318350e+00
            +9.90691050e-03 * tc[1]
            +4.16568667e-06 * tc[2]
            -8.33888875e-09 * tc[3]
            +3.16931420e-12 * tc[4]
            +1.92456290e+04 * invT;
        /*species 35: CH3CCH2 */
        species[35] =
            +1.73292090e+00
            +1.11973100e-02 * tc[1]
            -1.71635370e-06 * tc[2]
            -1.68991165e-09 * tc[3]
            +7.65064220e-13 * tc[4]
            +2.90404980e+04 * invT;
        /*species 36: C3H6 */
        species[36] =
            +1.49330700e+00
            +1.04625900e-02 * tc[1]
            +1.49559800e-06 * tc[2]
            -4.17228000e-09 * tc[3]
            +1.43162920e-12 * tc[4]
            +1.07482600e+03 * invT;
        /*species 37: nC3H7 */
        species[37] =
            +1.04911730e+00
            +1.30044865e-02 * tc[1]
            +7.84750533e-07 * tc[2]
            -4.89878300e-09 * tc[3]
            +1.87440414e-12 * tc[4]
            +1.03123460e+04 * invT;
        /*species 38: iC3H7 */
        species[38] =
            +1.44491990e+00
            +1.04995560e-02 * tc[1]
            +2.56787407e-06 * tc[2]
            -4.61906325e-09 * tc[3]
            +1.42565924e-12 * tc[4]
            +9.42237240e+03 * invT;
        /*species 39: C2H3CHO */
        species[39] =
            +1.27134980e+00
            +1.31155270e-02 * tc[1]
            -3.09707683e-06 * tc[2]
            -1.19593180e-09 * tc[3]
            +6.69610860e-13 * tc[4]
            -9.33573440e+03 * invT;
        /*species 40: C4H2 */
        species[40] =
            +1.05439780e+00
            +2.08134800e-02 * tc[1]
            -2.19572613e-05 * tc[2]
            +1.33142687e-08 * tc[3]
            -3.33663240e-12 * tc[4]
            +5.41852110e+04 * invT;
        /*species 41: iC4H3 */
        species[41] =
            +3.72214820e+00
            +1.29787715e-02 * tc[1]
            -8.78544767e-06 * tc[2]
            +3.87723000e-09 * tc[3]
            -7.60811300e-13 * tc[4]
            +5.88371210e+04 * invT;
        /*species 42: C4H4 */
        species[42] =
            +5.88570480e-01
            +1.82733425e-02 * tc[1]
            -1.13689893e-05 * tc[2]
            +4.16315475e-09 * tc[3]
            -6.01292460e-13 * tc[4]
            +3.33594920e+04 * invT;
        /*species 43: iC4H5 */
        species[43] =
            +1.13081050e-01
            +2.04753075e-02 * tc[1]
            -1.18045270e-05 * tc[2]
            +3.88274225e-09 * tc[3]
            -4.67102440e-13 * tc[4]
            +3.63833710e+04 * invT;
        /*species 44: C4H5-2 */
        species[44] =
            +2.96962800e+00
            +1.22211225e-02 * tc[1]
            -3.04171413e-06 * tc[2]
            -1.06167177e-18 * tc[3]
            +3.26094560e-22 * tc[4]
            +3.55033160e+04 * invT;
        /*species 45: C4H6 */
        species[45] =
            +1.12844650e-01
            +1.71845110e-02 * tc[1]
            -3.70246400e-06 * tc[2]
            -2.30266650e-09 * tc[3]
            +1.24130358e-12 * tc[4]
            +1.18022700e+04 * invT;
        /*species 46: C4H612 */
        species[46] =
            +1.02346700e+00
            +1.74795950e-02 * tc[1]
            -7.33635000e-06 * tc[2]
            +1.73556800e-09 * tc[3]
            -1.57583740e-13 * tc[4]
            +1.81179900e+04 * invT;
        /*species 47: C4H6-2 */
        species[47] =
            +2.13733380e+00
            +1.32431145e-02 * tc[1]
            -3.01895703e-06 * tc[2]
            -1.38465992e-19 * tc[3]
            +4.25637680e-23 * tc[4]
            +1.57109020e+04 * invT;
        /*species 48: C4H7 */
        species[48] =
            +7.44494320e-01
            +1.98394285e-02 * tc[1]
            -7.63269533e-06 * tc[2]
            +5.33824325e-10 * tc[3]
            +4.61927500e-13 * tc[4]
            +2.26533280e+04 * invT;
        /*species 49: C4H81 */
        species[49] =
            +1.18113800e+00
            +1.54266900e-02 * tc[1]
            +1.69550823e-06 * tc[2]
            -6.16372200e-09 * tc[3]
            +2.22203860e-12 * tc[4]
            -1.79040040e+03 * invT;
        /*species 50: pC4H9 */
        species[50] =
            +1.20870420e+00
            +1.91487485e-02 * tc[1]
            -2.42201697e-06 * tc[2]
            -3.85713675e-09 * tc[3]
            +1.73718870e-12 * tc[4]
            +7.32210400e+03 * invT;
        /*species 55: N2 */
        species[55] =
            +3.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 * invT;
    } else {
        /*species 0: H */
        species[0] =
            +2.50000001e+00
            -1.15421486e-11 * tc[1]
            +5.38539827e-15 * tc[2]
            -1.18378809e-18 * tc[3]
            +9.96394714e-23 * tc[4]
            +2.54736599e+04 * invT;
        /*species 1: O */
        species[1] =
            +2.56942078e+00
            -4.29870569e-05 * tc[1]
            +1.39828196e-08 * tc[2]
            -2.50444497e-12 * tc[3]
            +2.45667382e-16 * tc[4]
            +2.92175791e+04 * invT;
        /*species 2: OH */
        species[2] =
            +2.86472886e+00
            +5.28252240e-04 * tc[1]
            -8.63609193e-08 * tc[2]
            +7.63046685e-12 * tc[3]
            -2.66391752e-16 * tc[4]
            +3.71885774e+03 * invT;
        /*species 3: HO2 */
        species[3] =
            +4.01721090e+00
            +1.11991006e-03 * tc[1]
            -2.11219383e-07 * tc[2]
            +2.85615925e-11 * tc[3]
            -2.15817070e-15 * tc[4]
            +1.11856713e+02 * invT;
        /*species 4: H2 */
        species[4] =
            +3.33727920e+00
            -2.47012365e-05 * tc[1]
            +1.66485593e-07 * tc[2]
            -4.48915985e-11 * tc[3]
            +4.00510752e-15 * tc[4]
            -9.50158922e+02 * invT;
        /*species 5: H2O */
        species[5] =
            +3.03399249e+00
            +1.08845902e-03 * tc[1]
            -5.46908393e-08 * tc[2]
            -2.42604967e-11 * tc[3]
            +3.36401984e-15 * tc[4]
            -3.00042971e+04 * invT;
        /*species 6: H2O2 */
        species[6] =
            +4.16500285e+00
            +2.45415847e-03 * tc[1]
            -6.33797417e-07 * tc[2]
            +9.27964965e-11 * tc[3]
            -5.75816610e-15 * tc[4]
            -1.78617877e+04 * invT;
        /*species 7: O2 */
        species[7] =
            +3.28253784e+00
            +7.41543770e-04 * tc[1]
            -2.52655556e-07 * tc[2]
            +5.23676387e-11 * tc[3]
            -4.33435588e-15 * tc[4]
            -1.08845772e+03 * invT;
        /*species 8: CH */
        species[8] =
            +2.87846473e+00
            +4.85456840e-04 * tc[1]
            +4.81485517e-08 * tc[2]
            -3.26719623e-11 * tc[3]
            +3.52158766e-15 * tc[4]
            +7.10124364e+04 * invT;
        /*species 9: CH2 */
        species[9] =
            +2.87410113e+00
            +1.82819646e-03 * tc[1]
            -4.69648657e-07 * tc[2]
            +6.50448872e-11 * tc[3]
            -3.75455134e-15 * tc[4]
            +4.62636040e+04 * invT;
        /*species 10: CH2* */
        species[10] =
            +2.29203842e+00
            +2.32794318e-03 * tc[1]
            -6.70639823e-07 * tc[2]
            +1.04476500e-10 * tc[3]
            -6.79432730e-15 * tc[4]
            +5.09259997e+04 * invT;
        /*species 11: CH3 */
        species[11] =
            +2.28571772e+00
            +3.61995018e-03 * tc[1]
            -9.95714493e-07 * tc[2]
            +1.48921161e-10 * tc[3]
            -9.34308788e-15 * tc[4]
            +1.67755843e+04 * invT;
        /*species 12: CH4 */
        species[12] =
            +7.48514950e-02
            +6.69547335e-03 * tc[1]
            -1.91095270e-06 * tc[2]
            +3.05731338e-10 * tc[3]
            -2.03630460e-14 * tc[4]
            -9.46834459e+03 * invT;
        /*species 13: HCO */
        species[13] =
            +2.77217438e+00
            +2.47847763e-03 * tc[1]
            -8.28152043e-07 * tc[2]
            +1.47290445e-10 * tc[3]
            -1.06701742e-14 * tc[4]
            +4.01191815e+03 * invT;
        /*species 14: CH2O */
        species[14] =
            +1.76069008e+00
            +4.60000041e-03 * tc[1]
            -1.47419604e-06 * tc[2]
            +2.51603030e-10 * tc[3]
            -1.76771128e-14 * tc[4]
            -1.39958323e+04 * invT;
        /*species 15: CH3O */
        species[15] =
            +4.75779238e+00
            +3.72071237e-03 * tc[1]
            -8.99017253e-07 * tc[2]
            +1.09522626e-10 * tc[3]
            -5.27074196e-15 * tc[4]
            +3.78111940e+02 * invT;
        /*species 16: CH2OH */
        species[16] =
            +5.09314370e+00
            +2.97380630e-03 * tc[1]
            -6.88324867e-07 * tc[2]
            +8.07520433e-11 * tc[3]
            -3.76251804e-15 * tc[4]
            -4.03409640e+03 * invT;
        /*species 17: CH3OH */
        species[17] =
            +1.78970791e+00
            +7.04691460e-03 * tc[1]
            -2.12166945e-06 * tc[2]
            +3.45427713e-10 * tc[3]
            -2.34120440e-14 * tc[4]
            -2.53748747e+04 * invT;
        /*species 18: CO */
        species[18] =
            +2.71518561e+00
            +1.03126372e-03 * tc[1]
            -3.32941924e-07 * tc[2]
            +5.75132520e-11 * tc[3]
            -4.07295432e-15 * tc[4]
            -1.41518724e+04 * invT;
        /*species 19: CO2 */
        species[19] =
            +3.85746029e+00
            +2.20718513e-03 * tc[1]
            -7.38271347e-07 * tc[2]
            +1.30872547e-10 * tc[3]
            -9.44168328e-15 * tc[4]
            -4.87591660e+04 * invT;
        /*species 20: C2H */
        species[20] =
            +3.16780652e+00
            +2.37610951e-03 * tc[1]
            -6.12623590e-07 * tc[2]
            +7.60475630e-11 * tc[3]
            -3.54465540e-15 * tc[4]
            +6.71210650e+04 * invT;
        /*species 21: C2H2 */
        species[21] =
            +4.14756964e+00
            +2.98083332e-03 * tc[1]
            -7.90982840e-07 * tc[2]
            +1.16853043e-10 * tc[3]
            -7.22470426e-15 * tc[4]
            +2.59359992e+04 * invT;
        /*species 22: C2H3 */
        species[22] =
            +3.01672400e+00
            +5.16511460e-03 * tc[1]
            -1.56027450e-06 * tc[2]
            +2.54408220e-10 * tc[3]
            -1.72521408e-14 * tc[4]
            +3.46128739e+04 * invT;
        /*species 23: C2H4 */
        species[23] =
            +2.03611116e+00
            +7.32270755e-03 * tc[1]
            -2.23692638e-06 * tc[2]
            +3.68057308e-10 * tc[3]
            -2.51412122e-14 * tc[4]
            +4.93988614e+03 * invT;
        /*species 24: C2H5 */
        species[24] =
            +1.95465642e+00
            +8.69863610e-03 * tc[1]
            -2.66068889e-06 * tc[2]
            +4.38044223e-10 * tc[3]
            -2.99283152e-14 * tc[4]
            +1.28575200e+04 * invT;
        /*species 25: C2H6 */
        species[25] =
            +1.07188150e+00
            +1.08426339e-02 * tc[1]
            -3.34186890e-06 * tc[2]
            +5.53530003e-10 * tc[3]
            -3.80005780e-14 * tc[4]
            -1.14263932e+04 * invT;
        /*species 26: HCCO */
        species[26] =
            +5.62820580e+00
            +2.04267005e-03 * tc[1]
            -5.31151567e-07 * tc[2]
            +7.15651300e-11 * tc[3]
            -3.88156640e-15 * tc[4]
            +1.93272150e+04 * invT;
        /*species 27: CH2CO */
        species[27] =
            +4.51129732e+00
            +4.50179872e-03 * tc[1]
            -1.38979878e-06 * tc[2]
            +2.30836470e-10 * tc[3]
            -1.58967640e-14 * tc[4]
            -7.77850000e+03 * invT;
        /*species 28: CH3CO */
        species[28] =
            +5.94477310e+00
            +3.93336025e-03 * tc[1]
            -9.62196067e-07 * tc[2]
            +1.18177188e-10 * tc[3]
            -5.71997220e-15 * tc[4]
            -3.78730750e+03 * invT;
        /*species 29: CH2CHO */
        species[29] =
            +5.97566990e+00
            +4.06529570e-03 * tc[1]
            -9.14541500e-07 * tc[2]
            +1.01757603e-10 * tc[3]
            -4.35203420e-15 * tc[4]
            -9.69500000e+02 * invT;
        /*species 30: CH3CHO */
        species[30] =
            +5.40411080e+00
            +5.86152950e-03 * tc[1]
            -1.40877123e-06 * tc[2]
            +1.70931128e-10 * tc[3]
            -8.19697260e-15 * tc[4]
            -2.25931220e+04 * invT;
        /*species 31: C3H3 */
        species[31] =
            +7.14221880e+00
            +3.80951002e-03 * tc[1]
            -8.91533167e-07 * tc[2]
            +1.06228700e-10 * tc[3]
            -5.02950830e-15 * tc[4]
            +3.89087427e+04 * invT;
        /*species 32: pC3H4 */
        species[32] =
            +6.02524000e+00
            +5.66827100e-03 * tc[1]
            -1.34077970e-06 * tc[2]
            +1.60940157e-10 * tc[3]
            -7.65992700e-15 * tc[4]
            +1.96209420e+04 * invT;
        /*species 33: aC3H4 */
        species[33] =
            +6.31687220e+00
            +5.56686400e-03 * tc[1]
            -1.32097927e-06 * tc[2]
            +1.58910595e-10 * tc[3]
            -7.57510800e-15 * tc[4]
            +2.01174950e+04 * invT;
        /*species 34: aC3H5 */
        species[34] =
            +6.50078770e+00
            +7.16236550e-03 * tc[1]
            -1.89272107e-06 * tc[2]
            +2.77020025e-10 * tc[3]
            -1.80727774e-14 * tc[4]
            +1.74824490e+04 * invT;
        /*species 35: CH3CCH2 */
        species[35] =
            +5.42555280e+00
            +7.75553600e-03 * tc[1]
            -1.88927833e-06 * tc[2]
            +1.98060970e-10 * tc[3]
            -3.37560680e-15 * tc[4]
            +2.78430270e+04 * invT;
        /*species 36: C3H6 */
        species[36] =
            +6.73225700e+00
            +7.45417000e-03 * tc[1]
            -1.64996633e-06 * tc[2]
            +1.80300550e-10 * tc[3]
            -7.53240800e-15 * tc[4]
            -9.23570300e+02 * invT;
        /*species 37: nC3H7 */
        species[37] =
            +7.70974790e+00
            +8.01574250e-03 * tc[1]
            -1.75734127e-06 * tc[2]
            +1.89720880e-10 * tc[3]
            -7.77254380e-15 * tc[4]
            +7.97622360e+03 * invT;
        /*species 38: iC3H7 */
        species[38] =
            +6.51927410e+00
            +8.61005200e-03 * tc[1]
            -1.91214057e-06 * tc[2]
            +2.10326830e-10 * tc[3]
            -8.91318260e-15 * tc[4]
            +7.32271930e+03 * invT;
        /*species 39: C2H3CHO */
        species[39] =
            +5.81118680e+00
            +8.55712800e-03 * tc[1]
            -2.49447203e-06 * tc[2]
            +3.56306225e-10 * tc[3]
            -1.83493682e-14 * tc[4]
            -1.07840540e+04 * invT;
        /*species 40: C4H2 */
        species[40] =
            +9.15763280e+00
            +2.77152590e-03 * tc[1]
            -4.53053467e-07 * tc[2]
            +4.69501875e-12 * tc[3]
            +4.63790720e-15 * tc[4]
            +5.25880390e+04 * invT;
        /*species 41: iC4H3 */
        species[41] =
            +7.65385480e+00
            +5.60202750e-03 * tc[1]
            -1.54671140e-06 * tc[2]
            +2.16966597e-10 * tc[3]
            -1.14861124e-14 * tc[4]
            +5.79543630e+04 * invT;
        /*species 42: C4H4 */
        species[42] =
            +7.25396010e+00
            +6.95704700e-03 * tc[1]
            -1.76440713e-06 * tc[2]
            +2.08701125e-10 * tc[3]
            -7.03957640e-15 * tc[4]
            +3.17660160e+04 * invT;
        /*species 43: iC4H5 */
        species[43] =
            +6.96460290e+00
            +9.13716650e-03 * tc[1]
            -2.60445783e-06 * tc[2]
            +3.82303850e-10 * tc[3]
            -2.18409860e-14 * tc[4]
            +3.47250980e+04 * invT;
        /*species 44: C4H5-2 */
        species[44] =
            +1.45381710e+01
            -4.28385280e-03 * tc[1]
            +7.85317467e-06 * tc[2]
            -3.41909475e-09 * tc[3]
            +4.88738540e-13 * tc[4]
            +3.32590950e+04 * invT;
        /*species 45: C4H6 */
        species[45] =
            +8.86731340e+00
            +7.45933500e-03 * tc[1]
            -1.05162387e-06 * tc[2]
            -1.04603325e-10 * tc[3]
            +3.15225160e-14 * tc[4]
            +9.13385160e+03 * invT;
        /*species 46: C4H612 */
        species[46] =
            +1.78155700e+01
            -2.12875100e-03 * tc[1]
            +3.50395000e-06 * tc[2]
            -1.11846100e-09 * tc[3]
            +1.16962760e-13 * tc[4]
            +1.26734200e+04 * invT;
        /*species 47: C4H6-2 */
        species[47] =
            +9.03381330e+00
            +4.10622550e-03 * tc[1]
            +2.39179840e-06 * tc[2]
            -1.47085835e-09 * tc[3]
            +2.06878300e-13 * tc[4]
            +1.43350680e+04 * invT;
        /*species 48: C4H7 */
        species[48] =
            +7.01348350e+00
            +1.13172790e-02 * tc[1]
            -3.08484900e-06 * tc[2]
            +4.20198175e-10 * tc[3]
            -2.08172340e-14 * tc[4]
            +2.09550080e+04 * invT;
        /*species 49: C4H81 */
        species[49] =
            +2.05358410e+00
            +1.71752535e-02 * tc[1]
            -5.29439900e-06 * tc[2]
            +8.27241550e-10 * tc[3]
            -5.07220900e-14 * tc[4]
            -2.13972310e+03 * invT;
        /*species 50: pC4H9 */
        species[50] =
            +8.68223950e+00
            +1.18455355e-02 * tc[1]
            -2.53162883e-06 * tc[2]
            +1.66067840e-10 * tc[3]
            +1.09690272e-14 * tc[4]
            +4.96440580e+03 * invT;
        /*species 55: N2 */
        species[55] =
            +2.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
    }

    /*species with midpoint at T=1392 kelvin */
    if (T < 1392) {
        /*species 52: C6H12 */
        species[52] =
            -1.35275205e+00
            +3.49327713e-02 * tc[1]
            -1.53136007e-05 * tc[2]
            +3.92418358e-09 * tc[3]
            -4.42592350e-13 * tc[4]
            -7.34368617e+03 * invT;
        /*species 54: C5H10 */
        species[54] =
            -1.06223481e+00
            +2.87109147e-02 * tc[1]
            -1.24828963e-05 * tc[2]
            +3.18412472e-09 * tc[3]
            -3.59219578e-13 * tc[4]
            -4.46546666e+03 * invT;
    } else {
        /*species 52: C6H12 */
        species[52] =
            +1.78337529e+01
            +1.33688829e-02 * tc[1]
            -3.03345591e-06 * tc[2]
            +3.52049420e-10 * tc[3]
            -1.63024849e-14 * tc[4]
            -1.42062860e+04 * invT;
        /*species 54: C5H10 */
        species[54] =
            +1.45851539e+01
            +1.12036235e-02 * tc[1]
            -2.54449342e-06 * tc[2]
            +2.95472415e-10 * tc[3]
            -1.36877028e-14 * tc[4]
            -1.00898205e+04 * invT;
    }

    /*species with midpoint at T=1389 kelvin */
    if (T < 1389) {
        /*species 53: C6H11 */
        species[53] =
            -1.55544944e+00
            +3.38432801e-02 * tc[1]
            -1.49016212e-05 * tc[2]
            +3.80591575e-09 * tc[3]
            -4.28692754e-13 * tc[4]
            +9.66316695e+03 * invT;
    } else {
        /*species 53: C6H11 */
        species[53] =
            +1.77336550e+01
            +1.24467388e-02 * tc[1]
            -2.86663817e-06 * tc[2]
            +3.36032070e-10 * tc[3]
            -1.56695133e-14 * tc[4]
            +2.68017174e+03 * invT;
    }

    /*species with midpoint at T=1391 kelvin */
    if (T < 1391) {
        /*species 51: NC12H26 */
        species[51] =
            -2.62181594e+00
            +7.36188555e-02 * tc[1]
            -3.14656757e-05 * tc[2]
            +7.68603170e-09 * tc[3]
            -8.07204460e-13 * tc[4]
            -4.00654253e+04 * invT;
    } else {
        /*species 51: NC12H26 */
        species[51] =
            +3.85095037e+01
            +2.81775024e-02 * tc[1]
            -6.38310667e-06 * tc[2]
            +7.40062155e-10 * tc[3]
            -3.42488300e-14 * tc[4]
            -5.48843465e+04 * invT;
    }
    return;
}


/*save molecular weights into array */
__device__ void molecularWeight_d(double * wt)
{
    wt[0] = 1.007970; /*H */
    wt[1] = 15.999400; /*O */
    wt[2] = 17.007370; /*OH */
    wt[3] = 33.006770; /*HO2 */
    wt[4] = 2.015940; /*H2 */
    wt[5] = 18.015340; /*H2O */
    wt[6] = 34.014740; /*H2O2 */
    wt[7] = 31.998800; /*O2 */
    wt[8] = 13.019120; /*CH */
    wt[9] = 14.027090; /*CH2 */
    wt[10] = 14.027090; /*CH2* */
    wt[11] = 15.035060; /*CH3 */
    wt[12] = 16.043030; /*CH4 */
    wt[13] = 29.018520; /*HCO */
    wt[14] = 30.026490; /*CH2O */
    wt[15] = 31.034460; /*CH3O */
    wt[16] = 31.034460; /*CH2OH */
    wt[17] = 32.042430; /*CH3OH */
    wt[18] = 28.010550; /*CO */
    wt[19] = 44.009950; /*CO2 */
    wt[20] = 25.030270; /*C2H */
    wt[21] = 26.038240; /*C2H2 */
    wt[22] = 27.046210; /*C2H3 */
    wt[23] = 28.054180; /*C2H4 */
    wt[24] = 29.062150; /*C2H5 */
    wt[25] = 30.070120; /*C2H6 */
    wt[26] = 41.029670; /*HCCO */
    wt[27] = 42.037640; /*CH2CO */
    wt[28] = 43.045610; /*CH3CO */
    wt[29] = 43.045610; /*CH2CHO */
    wt[30] = 44.053580; /*CH3CHO */
    wt[31] = 39.057360; /*C3H3 */
    wt[32] = 40.065330; /*pC3H4 */
    wt[33] = 40.065330; /*aC3H4 */
    wt[34] = 41.073300; /*aC3H5 */
    wt[35] = 41.073300; /*CH3CCH2 */
    wt[36] = 42.081270; /*C3H6 */
    wt[37] = 43.089240; /*nC3H7 */
    wt[38] = 43.089240; /*iC3H7 */
    wt[39] = 56.064730; /*C2H3CHO */
    wt[40] = 50.060540; /*C4H2 */
    wt[41] = 51.068510; /*iC4H3 */
    wt[42] = 52.076480; /*C4H4 */
    wt[43] = 53.084450; /*iC4H5 */
    wt[44] = 53.084450; /*C4H5-2 */
    wt[45] = 54.092420; /*C4H6 */
    wt[46] = 54.092420; /*C4H612 */
    wt[47] = 54.092420; /*C4H6-2 */
    wt[48] = 55.100390; /*C4H7 */
    wt[49] = 56.108360; /*C4H81 */
    wt[50] = 57.116330; /*pC4H9 */
    wt[51] = 170.341020; /*NC12H26 */
    wt[52] = 84.162540; /*C6H12 */
    wt[53] = 83.154570; /*C6H11 */
    wt[54] = 70.135450; /*C5H10 */
    wt[55] = 28.013400; /*N2 */

    return;
}
/* get temperature given internal energy in mass units and mass fracs */
__device__ void get_t_given_ey_d_(double * e, double * y_wk, double * t, int * ierr)
{
    const int maxiter = 200;
    const double tol  = 1.e-6;
    double ein  = *e;
    double tmin = 90;/*max lower bound for thermo def */
    double tmax = 4000;/*min upper bound for thermo def */
    double e1,emin,emax,cv,t1,dt;
    int i;/* loop counter */
    ckubms_d(&tmin, y_wk, &emin);
    ckubms_d(&tmax, y_wk, &emax);
    if (ein < emin) {
        /*Linear Extrapolation below tmin */
        ckcvbs_d(&tmin, y_wk, &cv);
        *t = tmin - (emin-ein)/cv;
        *ierr = 1;
        return;
    }
    if (ein > emax) {
        /*Linear Extrapolation above tmax */
        ckcvbs_d(&tmax, y_wk,&cv);
        *t = tmax - (emax-ein)/cv;
        *ierr = 1;
        return;
    }
    t1 = *t;
    if (t1 < tmin || t1 > tmax) {
        t1 = tmin + (tmax-tmin)/(emax-emin)*(ein-emin);
    }
    for (i = 0; i < maxiter; ++i) {
        ckubms_d(&t1,y_wk,&e1);
        ckcvbs_d(&t1,y_wk,&cv);
        dt = (ein - e1) / cv;
        if (dt > 100.) { dt = 100.; }
        else if (dt < -100.) { dt = -100.; }
        else if (fabs(dt) < tol) break;
        else if (t1+dt == t1) break;
        t1 += dt;
    }
    *t = t1;
    *ierr = 0;
    return;
}
/* get temperature given enthalpy in mass units and mass fracs */
__device__ void get_t_given_hy_d_(double * h, double * y_wk, double * t, int * ierr)
{
    const int maxiter = 200;
    const double tol  = 1.e-6;
    double hin  = *h;
    double tmin = 90;/*max lower bound for thermo def */
    double tmax = 4000;/*min upper bound for thermo def */
    double h1,hmin,hmax,cp,t1,dt;
    int i;/* loop counter */
    ckhbms_d(&tmin, y_wk,  &hmin);
    ckhbms_d(&tmax, y_wk,  &hmax);
    if (hin < hmin) {
        /*Linear Extrapolation below tmin */
        ckcpbs_d(&tmin, y_wk, &cp);
        *t = tmin - (hmin-hin)/cp;
        *ierr = 1;
        return;
    }
    if (hin > hmax) {
        /*Linear Extrapolation above tmax */
        ckcpbs_d(&tmax, y_wk, &cp);
        *t = tmax - (hmax-hin)/cp;
        *ierr = 1;
        return;
    }
    t1 = *t;
    if (t1 < tmin || t1 > tmax) {
        t1 = tmin + (tmax-tmin)/(hmax-hmin)*(hin-hmin);
    }
    for (i = 0; i < maxiter; ++i) {
        ckhbms_d(&t1,y_wk,&h1);
        ckcpbs_d(&t1,y_wk,&cp);
        dt = (hin - h1) / cp;
        if (dt > 100.) { dt = 100.; }
        else if (dt < -100.) { dt = -100.; }
        else if (fabs(dt) < tol) break;
        else if (t1+dt == t1) break;
        t1 += dt;
    }
    *t = t1;
    *ierr = 0;
    return;
}

/* End of file  */




#if 0




\\
\\
\\  This is the mechanism file
\\
\\
!****************************************************************************
! Title: "A reduced high-temperature oxidation model for n-dodecane"
!
! Developed by Xiaoqing You and Hai Wang, 12/2007.
!
! Please contact Hai Wang at haiw@usc.edu for questions and comments. 
!
!****************************************************************************
ELEMENTS
O     H     C     N   
END
SPECIES
   H            O            OH           HO2          H2        
   H2O          H2O2         O2           CH           CH2          
   CH2*         CH3          CH4          HCO          CH2O         
   CH3O         CH2OH        CH3OH        CO           CO2          
   C2H          C2H2         C2H3         C2H4         C2H5         
   C2H6         HCCO         CH2CO        CH3CO        CH2CHO       
   CH3CHO       C3H3         pC3H4        aC3H4        aC3H5        
   CH3CCH2      C3H6         nC3H7        iC3H7        C2H3CHO      
   C4H2         iC4H3        C4H4         iC4H5        C4H5-2       
   C4H6         C4H612       C4H6-2       C4H7         C4H81     
   pC4H9        NC12H26      C6H12        C6H11        C5H10        
   N2
END
TRANS ALL
H                  0   145.000     2.050     0.000     0.000     0.000          
O                  0    80.000     2.750     0.000     0.000     0.000          
OH                 1    80.000     2.750     0.000     0.000     0.000          
HO2                2   107.400     3.458     0.000     0.000     1.000          
H2                 1    38.000     2.920     0.000     0.790   280.000          
H2O                2   572.400     2.605     1.844     0.000     4.000          
H2O2               2   107.400     3.458     0.000     0.000     3.800          
O2                 1   107.400     3.458     0.000     1.600     3.800          
CH                 1    80.000     2.750     0.000     0.000     0.000          
CH2                1   144.000     3.800     0.000     0.000     0.000          
CH2*               1   144.000     3.800     0.000     0.000     0.000          
CH3                1   144.000     3.800     0.000     0.000     0.000          
CH4                2   141.400     3.746     0.000     2.600    13.000          
HCO                2   498.000     3.590     0.000     0.000     0.000          
CH2O               2   498.000     3.590     0.000     0.000     2.000          
CH3O               2   417.000     3.690     1.700     0.000     2.000          
CH2OH              2   417.000     3.690     1.700     0.000     2.000          
CH3OH              2   481.800     3.626     0.000     0.000     1.000          
CO                 1    98.100     3.650     0.000     1.950     1.800          
CO2                1   244.000     3.763     0.000     2.650     2.100          
C2H                1   209.000     4.100     0.000     0.000     2.500          
C2H2               1   209.000     4.100     0.000     0.000     2.500          
C2H3               2   209.000     4.100     0.000     0.000     1.000          
C2H4               2   280.800     3.971     0.000     0.000     1.500          
C2H5               2   252.300     4.302     0.000     0.000     1.500          
C2H6               2   252.300     4.302     0.000     0.000     1.500          
HCCO               2   150.000     2.500     0.000     0.000     1.000          
CH2CO              2   436.000     3.970     0.000     0.000     2.000          
CH3CO              2   436.000     3.970     0.000     0.000     2.000          
CH2CHO             2   436.000     3.970     0.000     0.000     2.000          
CH3CHO             2   436.000     3.970     0.000     0.000     2.000          
C3H3               2   252.000     4.760     0.000     0.000     1.000          
pC3H4              1   252.000     4.760     0.000     0.000     1.000          
aC3H4              1   252.000     4.760     0.000     0.000     1.000          
aC3H5              2   266.800     4.982     0.000     0.000     1.000          
CH3CCH2            2   266.800     4.982     0.000     0.000     1.000          
C3H6               2   266.800     4.982     0.000     0.000     1.000          
nC3H7              2   266.800     4.982     0.000     0.000     1.000          
iC3H7              2   266.800     4.982     0.000     0.000     1.000          
C2H3CHO            2   357.000     5.176     0.000     0.000     1.000          
C4H2               1   357.000     5.180     0.000     0.000     1.000          
iC4H3              2   357.000     5.180     0.000     0.000     1.000          
C4H4               2   357.000     5.180     0.000     0.000     1.000          
iC4H5              2   357.000     5.176     0.000     0.000     1.000          
C4H5-2             2   357.000     5.180     0.000     0.000     1.000          
C4H6               2   357.000     5.176     0.000     0.000     1.000          
C4H612             2   357.000     5.180     0.000     0.000     1.000          
C4H6-2             2   357.000     5.180     0.000     0.000     1.000          
C4H7               2   357.000     5.176     0.000     0.000     1.000          
C4H81              2   357.000     5.176     0.000     0.000     1.000          
pC4H9              2   357.000     5.176     0.000     0.000     1.000          
NC12H26            2   789.980     7.047     0.000     0.000     1.000          
C6H12              2   504.629     5.628     0.000     0.000     1.000          
C6H11              2   504.629     5.628     0.000     0.000     1.000          
C5H10              2   448.508     5.342     0.000     0.000     1.000          
N2                 1    97.530     3.621     0.000     1.760     4.000          
END
REACTIONS
H+O2<=>O+OH                             2.6440E+16   -0.671  17041.0
O+H2<=>H+OH                             4.5890E+04    2.700   6260.0
OH+H2<=>H+H2O                           1.7340E+08    1.510   3430.0
2OH<=>O+H2O                             3.9730E+04    2.400  -2110.0
2H+M<=>H2+M                             1.7800E+18   -1.000      0.0
H2/ 0.00/ H2O/ 0.00/ 
CO2/ 0.00/ 
2H+H2O<=>H2+H2O                         5.6240E+19   -1.250      0.0
H+OH+M<=>H2O+M                          4.4000E+22   -2.000      0.0
H2/ 2.00/ H2O/ 6.30/ 
CO/ 1.75/ CO2/ 3.60/ 
O+H+M<=>OH+M                            9.4280E+18   -1.000      0.0
H2/ 2.00/ H2O/12.00/ 
CO/ 1.75/ CO2/ 3.60/ 
2O+M<=>O2+M                             1.2000E+17   -1.000      0.0
H2/ 2.40/ H2O/15.40/ 
CO/ 1.75/ CO2/ 3.60/ 
H+O2(+M)<=>HO2(+M)                      5.1160E+12    0.440      0.0
     LOW  / 6.3280E+19   -1.400      0.00/
     TROE/  0.50000E+00  0.10000E-29  0.10000E+31/
H2O/11.89/ O2/ 0.85/ 
CO/ 1.09/ CO2/ 2.18/ 
H2+O2<=>HO2+H                           5.9160E+05    2.433  53502.0
2OH(+M)<=>H2O2(+M)                      1.1100E+14   -0.370      0.0
     LOW  / 2.0100E+17   -0.584  -2293.00/
     TROE/  0.7346    94.0   1756.0   5182.0 /
H2/ 2.00/ H2O/ 6.00/ 
CO/ 1.75/ CO2/ 3.60/ 
HO2+H<=>O+H2O                           3.9700E+12    0.000    671.0
HO2+H<=>2OH                             7.4850E+13    0.000    295.0
HO2+O<=>OH+O2                           4.0000E+13    0.000      0.0
HO2+OH<=>O2+H2O                         2.3750E+13    0.000   -500.0
 DUPLICATE
HO2+OH<=>O2+H2O                         1.0000E+16    0.000  17330.0
 DUPLICATE
2HO2<=>O2+H2O2                          1.3000E+11    0.000  -1630.0
 DUPLICATE
2HO2<=>O2+H2O2                          3.6580E+14    0.000  12000.0
 DUPLICATE
H2O2+H<=>HO2+H2                         6.0500E+06    2.000   5200.0
H2O2+H<=>OH+H2O                         2.4100E+13    0.000   3970.0
H2O2+O<=>OH+HO2                         9.6300E+06    2.000   3970.0
H2O2+OH<=>HO2+H2O                       2.0000E+12    0.000    427.0
 DUPLICATE
H2O2+OH<=>HO2+H2O                       2.6700E+41   -7.000  37600.0
 DUPLICATE
CO+O(+M)<=>CO2(+M)                      1.3620E+10    0.000   2384.0
     LOW  / 1.1730E+24   -2.790   4191.00/
H2/ 2.00/ H2O/12.00/ 
CO/ 1.75/ CO2/ 3.60/ 
CO+OH<=>CO2+H                           8.0000E+11    0.140   7352.0
 DUPLICATE
CO+OH<=>CO2+H                           8.7840E+10    0.030    -16.0
 DUPLICATE
CO+HO2<=>CO2+OH                         3.0100E+13    0.000  23000.0
HCO+H<=>CO+H2                           1.2000E+14    0.000      0.0
HCO+O<=>CO+OH                           3.0000E+13    0.000      0.0
HCO+O<=>CO2+H                           3.0000E+13    0.000      0.0
HCO+OH<=>CO+H2O                         3.0200E+13    0.000      0.0
HCO+M<=>CO+H+M                          1.8700E+17   -1.000  17000.0
H2/ 2.00/ H2O/ 0.00/ CO/ 1.75/ 
CO2/ 3.60/ 
HCO+H2O<=>CO+H+H2O                      2.2440E+18   -1.000  17000.0
HCO+O2<=>CO+HO2                         1.2040E+10    0.807   -727.0
CH+O<=>CO+H                             5.7000E+13    0.000      0.0
CH+OH<=>HCO+H                           3.0000E+13    0.000      0.0
CH+H2<=>CH2+H                           1.1070E+08    1.790   1670.0
CH+H2O<=>CH2O+H                         5.7100E+12    0.000   -755.0
CH+O2<=>HCO+O                           3.3000E+13    0.000      0.0
CH+CO2<=>HCO+CO                         3.4000E+12    0.000    690.0
CH2+O<=>HCO+H                           8.0000E+13    0.000      0.0
CH2+OH<=>CH2O+H                         2.0000E+13    0.000      0.0
CH2+OH<=>CH+H2O                         1.1300E+07    2.000   3000.0
CH2+H2<=>H+CH3                          5.0000E+05    2.000   7230.0
CH2+O2<=>HCO+OH                         1.0600E+13    0.000   1500.0
CH2+O2<=>CO2+2H                         2.6400E+12    0.000   1500.0
CH2+HO2<=>CH2O+OH                       2.0000E+13    0.000      0.0
CH2+CO(+M)<=>CH2CO(+M)                  8.1000E+11    0.500   4510.0
     LOW  / 2.6900E+33   -5.110   7095.00/
     TROE/  0.5907   275.0   1226.0   5185.0 /
H2/ 2.00/ H2O/ 6.00/ 
CH4/ 2.00/ CO/ 1.50/ CO2/ 2.00/ C2H6/ 3.00/ 
CH2*+N2<=>CH2+N2                        1.5000E+13    0.000    600.0
CH2*+H<=>CH+H2                          3.0000E+13    0.000      0.0
CH2*+OH<=>CH2O+H                        3.0000E+13    0.000      0.0
CH2*+H2<=>CH3+H                         7.0000E+13    0.000      0.0
CH2*+O2<=>H+OH+CO                       2.8000E+13    0.000      0.0
CH2*+O2<=>CO+H2O                        1.2000E+13    0.000      0.0
CH2*+H2O(+M)<=>CH3OH(+M)                2.0000E+13    0.000      0.0
     LOW  / 2.7000E+38   -6.300   3100.00/
     TROE/  0.1507   134.0   2383.0   7265.0 /
H2/ 2.00/ H2O/ 6.00/ CH4/ 2.00/ 
CO/ 1.50/ CO2/ 2.00/ C2H6/ 3.00/ 
CH2*+H2O<=>CH2+H2O                      3.0000E+13    0.000      0.0
CH2*+CO<=>CH2+CO                        9.0000E+12    0.000      0.0
CH2*+CO2<=>CH2+CO2                      7.0000E+12    0.000      0.0
CH2*+CO2<=>CH2O+CO                      1.4000E+13    0.000      0.0
CH2O+H(+M)<=>CH2OH(+M)                  5.4000E+11    0.454   3600.0
     LOW  / 1.2700E+32   -4.820   6530.00/
     TROE/  0.7187   103.0   1291.0   4160.0 /
H2/ 2.00/ H2O/ 6.00/ CH4/ 2.00/ 
CO/ 1.50/ CO2/ 2.00/ C2H6/ 3.00/ 
CH2O+H(+M)<=>CH3O(+M)                   5.4000E+11    0.454   2600.0
     LOW  / 2.2000E+30   -4.800   5560.00/
     TROE/  0.7580    94.0   1555.0   4200.0 /
H2/ 2.00/ H2O/ 6.00/ CH4/ 2.00/ 
CO/ 1.50/ CO2/ 2.00/ C2H6/ 3.00/ 
CH2O+H<=>HCO+H2                         2.3000E+10    1.050   3275.0
CH2O+O<=>HCO+OH                         3.9000E+13    0.000   3540.0
CH2O+OH<=>HCO+H2O                       3.4300E+09    1.180   -447.0
CH2O+O2<=>HCO+HO2                       1.0000E+14    0.000  40000.0
CH2O+HO2<=>HCO+H2O2                     1.0000E+12    0.000   8000.0
CH2O+CH<=>CH2CO+H                       9.4600E+13    0.000   -515.0
CH3+H(+M)<=>CH4(+M)                     1.2700E+16   -0.630    383.0
     LOW  / 2.4770E+33   -4.760   2440.00/
     TROE/  0.7830    74.0   2941.0   6964.0 /
H2/ 2.00/ H2O/ 6.00/ 
CH4/ 2.00/ CO/ 1.50/ CO2/ 2.00/ C2H6/ 3.00/ 
CH3+O<=>CH2O+H                          8.4300E+13    0.000      0.0
CH3+OH(+M)<=>CH3OH(+M)                  6.3000E+13    0.000      0.0
     LOW  / 2.7000E+38   -6.300   3100.00/
     TROE/  0.2105    83.5   5398.0   8370.0 /
H2/ 2.00/ H2O/ 6.00/ CH4/ 2.00/ 
CO/ 1.50/ CO2/ 2.00/ C2H6/ 3.00/ 
CH3+OH<=>CH2+H2O                        5.6000E+07    1.600   5420.0
CH3+OH<=>CH2*+H2O                       2.5010E+13    0.000      0.0
CH3+O2<=>O+CH3O                         3.0830E+13    0.000  28800.0
CH3+O2<=>OH+CH2O                        3.6000E+10    0.000   8940.0
CH3+HO2<=>CH4+O2                        1.0000E+12    0.000      0.0
CH3+HO2<=>CH3O+OH                       1.3400E+13    0.000      0.0
CH3+CH<=>C2H3+H                         3.0000E+13    0.000      0.0
CH3+HCO<=>CH4+CO                        8.4800E+12    0.000      0.0
CH3+CH2O<=>CH4+HCO                      3.3200E+03    2.810   5860.0
CH3+CH2<=>C2H4+H                        4.0000E+13    0.000      0.0
2CH3(+M)<=>C2H6(+M)                     2.1200E+16   -0.970    620.0
     LOW  / 1.7700E+50   -9.670   6220.00/
     TROE/  0.5325   151.0   1038.0   4970.0 /
H2/ 2.00/ H2O/ 6.00/ 
CH4/ 2.00/ CO/ 1.50/ CO2/ 2.00/ C2H6/ 3.00/ 
2CH3<=>H+C2H5                           4.9900E+12    0.100  10600.0
CH3+HCCO<=>C2H4+CO                      5.0000E+13    0.000      0.0
CH3O+H<=>CH2O+H2                        2.0000E+13    0.000      0.0
CH3O+H<=>CH3+OH                         3.2000E+13    0.000      0.0
CH3O+H<=>CH2*+H2O                       1.6000E+13    0.000      0.0
CH3O+OH<=>CH2O+H2O                      5.0000E+12    0.000      0.0
CH3O+O2<=>CH2O+HO2                      4.2800E-13    7.600  -3530.0
CH2OH+H<=>CH2O+H2                       2.0000E+13    0.000      0.0
CH2OH+H<=>CH3+OH                        1.2000E+13    0.000      0.0
CH2OH+H<=>CH2*+H2O                      6.0000E+12    0.000      0.0
CH2OH+O2<=>CH2O+HO2                     1.8000E+13    0.000    900.0
CH4+H<=>CH3+H2                          6.6000E+08    1.620  10840.0
CH4+O<=>CH3+OH                          1.0200E+09    1.500   8600.0
CH4+OH<=>CH3+H2O                        1.0000E+08    1.600   3120.0
CH4+CH<=>C2H4+H                         6.0000E+13    0.000      0.0
CH4+CH2<=>2CH3                          2.4600E+06    2.000   8270.0
CH4+CH2*<=>2CH3                         1.6000E+13    0.000   -570.0
CH3OH+H<=>CH2OH+H2                      1.7000E+07    2.100   4870.0
CH3OH+H<=>CH3O+H2                       4.2000E+06    2.100   4870.0
CH3OH+O<=>CH2OH+OH                      3.8800E+05    2.500   3100.0
CH3OH+OH<=>CH2OH+H2O                    1.4400E+06    2.000   -840.0
CH3OH+OH<=>CH3O+H2O                     6.3000E+06    2.000   1500.0
C2H+O<=>CH+CO                           5.0000E+13    0.000      0.0
C2H+OH<=>H+HCCO                         2.0000E+13    0.000      0.0
C2H+O2<=>HCO+CO                         5.0000E+13    0.000   1500.0
C2H+H2<=>H+C2H2                         4.9000E+05    2.500    560.0
HCCO+H<=>CH2*+CO                        1.0000E+14    0.000      0.0
HCCO+O<=>H+2CO                          1.0000E+14    0.000      0.0
HCCO+O2<=>OH+2CO                        1.6000E+12    0.000    854.0
C2H3(+M)<=>C2H2+H(+M)                   3.8600E+08    1.620  37048.2
     LOW  / 2.5650E+27   -3.400  35798.72/
     TROE/  1.9816  5383.7      4.3     -0.1 /
H2/ 2.00/ H2O/ 6.00/ 
CH4/ 2.00/ CO/ 1.50/ CO2/ 2.00/ C2H2/ 3.00/ C2H4/ 3.00/ C2H6/ 3.00/ 
C2H2+O<=>C2H+OH                         4.6000E+19   -1.410  28950.0
C2H2+O<=>CH2+CO                         4.0800E+06    2.000   1900.0
C2H2+O<=>HCCO+H                         1.6320E+07    2.000   1900.0
C2H2+OH<=>CH2CO+H                       2.1800E-04    4.500  -1000.0
DUPLICATE
C2H2+OH<=>CH2CO+H                       5.0400E+05    2.300  13500.0
DUPLICATE
C2H2+OH<=>C2H+H2O                       3.3700E+07    2.000  14000.0
C2H2+HCO<=>C2H3+CO                      1.0000E+07    2.000   6000.0
C2H2+CH2<=>C3H3+H                       1.2000E+13    0.000   6620.0
C2H2+CH2*<=>C3H3+H                      2.0000E+13    0.000      0.0
C2H2+C2H<=>C4H2+H                       9.6000E+13    0.000      0.0
C2H2+CH3<=>pC3H4+H                      2.5600E+09    1.100  13644.0
C2H2+CH3<=>aC3H4+H                      5.1400E+09    0.860  22153.0
CH2CO+H(+M)<=>CH2CHO(+M)                3.3000E+14   -0.060   8500.0
     LOW  / 3.8000E+41   -7.640  11900.00/
     TROE/  0.3370  1707.0   3200.0   4131.0 /
H2/ 2.00/ H2O/ 6.00/ 
CH4/ 2.00/ CO/ 1.50/ CO2/ 2.00/ C2H2/ 3.00/ C2H4/ 3.00/ C2H6/ 3.00/ 
CH2CO+H<=>HCCO+H2                       5.0000E+13    0.000   8000.0
CH2CO+H<=>CH3+CO                        1.5000E+09    1.430   2690.0
CH2CO+OH<=>HCCO+H2O                     7.5000E+12    0.000   2000.0
C2H3+H(+M)<=>C2H4(+M)                   6.0800E+12    0.270    280.0
     LOW  / 1.4000E+30   -3.860   3320.00/
     TROE/  0.7820   207.5   2663.0   6095.0 /
H2/ 2.00/ H2O/ 6.00/ 
CH4/ 2.00/ CO/ 1.50/ CO2/ 2.00/ C2H2/ 3.00/ C2H4/ 3.00/ C2H6/ 3.00/ 
C2H3+H<=>C2H2+H2                        9.0000E+13    0.000      0.0
C2H3+O<=>CH2CO+H                        4.8000E+13    0.000      0.0
C2H3+O<=>CH3+CO                         4.8000E+13    0.000      0.0
C2H3+OH<=>C2H2+H2O                      3.0110E+13    0.000      0.0
C2H3+O2<=>C2H2+HO2                      1.3400E+06    1.610   -383.4
C2H3+O2<=>CH2CHO+O                      3.0000E+11    0.290     11.0
C2H3+O2<=>HCO+CH2O                      4.6000E+16   -1.390   1010.0
C2H3+HO2<=>CH2CHO+OH                    1.0000E+13    0.000      0.0
C2H3+HCO<=>C2H4+CO                      9.0330E+13    0.000      0.0
C2H3+HCO<=>C2H3CHO                      1.8000E+13    0.000      0.0
C2H3+CH3(+M)<=>C3H6(+M)                 2.5000E+13    0.000      0.0
     LOW  / 4.2700E+58  -11.940   9769.80/
     TROE/  0.1750  1340.6  60000.0  10139.8 /
H2/ 2.00/ H2O/ 6.00/ 
CH4/ 2.00/ CO/ 1.50/ CO2/ 2.00/ C2H4/ 3.00/ C2H6/ 3.00/ 
C2H3+CH3<=>aC3H5+H                      1.5000E+24   -2.830  18618.0
CH2CHO<=>CH3+CO                         7.8000E+41   -9.147  46900.0
CH2CHO+H<=>CH3CO+H                      5.0000E+12    0.000      0.0
CH2CHO+H<=>CH3+HCO                      9.0000E+13    0.000      0.0
CH2CHO+H<=>CH2CO+H2                     2.0000E+13    0.000   4000.0
CH2CHO+OH<=>CH2CO+H2O                   1.0000E+13    0.000   2000.0
CH2CHO+O2<=>CH2CO+HO2                   1.4000E+11    0.000      0.0
CH2CHO+O2<=>CH2O+CO+OH                  1.8000E+10    0.000      0.0
CH3+CO(+M)<=>CH3CO(+M)                  4.8500E+07    1.650   6150.0
     LOW  / 7.8000E+30   -5.395   8600.00/
     TROE/  0.2580   598.0  21002.0   1773.0 /
H2/ 2.00/ H2O/ 6.00/ 
CH4/ 2.00/ CO/ 1.50/ CO2/ 2.00/ C2H2/ 3.00/ C2H4/ 3.00/ C2H6/ 3.00/ 
CH3CO+H<=>CH3+HCO                       9.6000E+13    0.000      0.0
CH3CO+HO2<=>CH3+CO2+OH                  3.0000E+13    0.000      0.0
CH3+HCO(+M)<=>CH3CHO(+M)                1.8000E+13    0.000      0.0
     LOW  / 2.2000E+48   -9.588   5100.00/
     TROE/  0.6173    13.1   2078.0   5093.0 /
H2/ 2.00/ H2O/ 6.00/ CH4/ 2.00/ 
CO/ 1.50/ CO2/ 2.00/ C2H2/ 3.00/ C2H4/ 3.00/ C2H6/ 3.00/ 
CH3CHO+H<=>CH3CO+H2                     4.1000E+09    1.160   2400.0
CH3CHO+OH<=>CH3CO+H2O                   2.3500E+10    0.730  -1110.0
CH3CHO+CH3<=>CH3CO+CH4                  2.0000E-06    5.600   2460.0
CH3CHO+O2<=>CH3CO+HO2                   3.0000E+13    0.000  39100.0
C2H4(+M)<=>H2+C2H2(+M)                  8.0000E+12    0.440  88770.0
     LOW  / 7.0000E+50   -9.310  99860.00/
     TROE/  0.7345   180.0   1035.0   5417.0 /
H2/ 2.00/ H2O/ 6.00/ 
CH4/ 2.00/ CO/ 1.50/ CO2/ 2.00/ C2H6/ 3.00/ 
C2H4+H(+M)<=>C2H5(+M)                   1.3670E+09    1.463   1355.0
     LOW  / 2.0270E+39   -6.642   5769.00/
     TROE/ -0.5690   299.0   9147.0   -152.4 /
C2H4+H<=>C2H3+H2                        5.0700E+07    1.900  12950.0
C2H4+O<=>C2H3+OH                        1.5100E+07    1.900   3740.0
C2H4+O<=>CH3+HCO                        1.9200E+07    1.830    220.0
C2H4+O<=>CH2+CH2O                       3.8400E+05    1.830    220.0
C2H4+OH<=>C2H3+H2O                      3.6000E+06    2.000   2500.0
C2H4+HCO<=>C2H5+CO                      1.0000E+07    2.000   8000.0
C2H4+CH<=>aC3H4+H                       3.0000E+13    0.000      0.0
C2H4+CH<=>pC3H4+H                       3.0000E+13    0.000      0.0
C2H4+CH2<=>aC3H5+H                      2.0000E+13    0.000   6000.0
C2H4+CH2*<=>aC3H5+H                     5.0000E+13    0.000      0.0
C2H4+CH3<=>C2H3+CH4                     2.2700E+05    2.000   9200.0
C2H4+CH3<=>nC3H7                        3.3000E+11    0.000   7700.0
C2H4+C2H3<=>C4H7                        7.9300E+38   -8.470  14220.0
C2H5+H(+M)<=>C2H6(+M)                   5.2100E+17   -0.990   1580.0
     LOW  / 1.9900E+41   -7.080   6685.00/
     TROE/  0.8422   125.0   2219.0   6882.0 /
H2/ 2.00/ H2O/ 6.00/ 
CH4/ 2.00/ CO/ 1.50/ CO2/ 2.00/ C2H6/ 3.00/ 
C2H5+H<=>C2H4+H2                        2.0000E+12    0.000      0.0
C2H5+O<=>CH3+CH2O                       1.6040E+13    0.000      0.0
C2H5+O<=>CH3CHO+H                       8.0200E+13    0.000      0.0
C2H5+O2<=>C2H4+HO2                      2.0000E+10    0.000      0.0
C2H5+HO2<=>C2H6+O2                      3.0000E+11    0.000      0.0
C2H5+HO2<=>C2H4+H2O2                    3.0000E+11    0.000      0.0
C2H5+HO2<=>CH3+CH2O+OH                  2.4000E+13    0.000      0.0
C2H5+C2H3(+M)<=>C4H81(+M)               1.5000E+13    0.000      0.0
     LOW  / 1.5500E+56  -11.790   8984.50/
     TROE/  0.1980  2277.9  60000.0   5723.2 /
H2/ 2.00/ H2O/ 6.00/ 
CH4/ 2.00/ CO/ 1.50/ CO2/ 2.00/ C2H6/ 3.00/ 
C2H5+C2H3<=>aC3H5+CH3                   3.9000E+32   -5.220  19747.0
C2H6+H<=>C2H5+H2                        1.1500E+08    1.900   7530.0
C2H6+O<=>C2H5+OH                        8.9800E+07    1.920   5690.0
C2H6+OH<=>C2H5+H2O                      3.5400E+06    2.120    870.0
C2H6+CH2*<=>C2H5+CH3                    4.0000E+13    0.000   -550.0
C2H6+CH3<=>C2H5+CH4                     6.1400E+06    1.740  10450.0
C3H3+H<=>pC3H4                          1.5000E+13    0.000      0.0
C3H3+O<=>CH2O+C2H                       2.0000E+13    0.000      0.0
C3H3+O2<=>CH2CO+HCO                     3.0000E+10    0.000   2868.0
C3H3+HO2<=>pC3H4+O2                     2.5000E+12    0.000      0.0
aC3H4+H<=>CH3CCH2                       9.4600E+42   -9.430  11190.0
aC3H4+H<=>aC3H5                         1.5200E+59  -13.540  26949.0
aC3H4+O<=>C2H4+CO                       2.0000E+07    1.800   1000.0
aC3H4+OH<=>C3H3+H2O                     5.3000E+06    2.000   2000.0
pC3H4<=>aC3H4                           5.1500E+60  -13.930  91117.0
pC3H4+H<=>aC3H4+H                       6.2700E+17   -0.910  10079.0
pC3H4+H<=>CH3CCH2                       1.6600E+47  -10.580  13690.0
pC3H4+O<=>C2H4+CO                       1.0000E+13    0.000   2250.0
pC3H4+OH<=>C3H3+H2O                     1.0000E+06    2.000    100.0
aC3H5+H(+M)<=>C3H6(+M)                  2.0000E+14    0.000      0.0
     LOW  / 1.3300E+60  -12.000   5967.80/
     TROE/  0.0200  1096.6   1096.6   6859.5 /
H2/ 2.00/ H2O/ 6.00/ 
CH4/ 2.00/ CO/ 1.50/ CO2/ 2.00/ C2H6/ 3.00/ 
aC3H5+H<=>aC3H4+H2                      1.8000E+13    0.000      0.0
aC3H5+O<=>C2H3CHO+H                     6.0000E+13    0.000      0.0
aC3H5+OH<=>C2H3CHO+2H                   4.2000E+32   -5.160  30126.0
aC3H5+OH<=>aC3H4+H2O                    6.0000E+12    0.000      0.0
aC3H5+HO2<=>C3H6+O2                     2.6600E+12    0.000      0.0
aC3H5+HO2<=>OH+C2H3+CH2O                6.6000E+12    0.000      0.0
aC3H5+HCO<=>C3H6+CO                     6.0000E+13    0.000      0.0
aC3H5+CH3(+M)<=>C4H81(+M)               1.0000E+14   -0.320   -262.3
     LOW  / 3.9100E+60  -12.810   6250.00/
     TROE/  0.1040  1606.0  60000.0   6118.4 /
H2/ 2.00/ H2O/ 6.00/ 
CH4/ 2.00/ CO/ 1.50/ CO2/ 2.00/ C2H6/ 3.00/ 
aC3H5+CH3<=>aC3H4+CH4                   3.0000E+12   -0.320   -131.0
CH3CCH2+O2<=>CH3CO+CH2O                 1.0000E+11    0.000      0.0
CH3CCH2+HO2<=>CH3+CH2CO+OH              2.0000E+13    0.000      0.0
C3H6+H(+M)<=>nC3H7(+M)                  1.3300E+13    0.000   3260.7
     LOW  / 6.2600E+38   -6.660   7000.00/
     TROE/  1.0000  1000.0   1310.0  48097.0 /
H2/ 2.00/ H2O/ 6.00/ 
CH4/ 2.00/ CO/ 1.50/ CO2/ 2.00/ C2H6/ 3.00/ 
C3H6+H(+M)<=>iC3H7(+M)                  1.3300E+13    0.000   1559.8
     LOW  / 8.7000E+42   -7.500   4721.80/
     TROE/  1.0000  1000.0    645.4   6844.3 /
H2/ 2.00/ H2O/ 6.00/ 
CH4/ 2.00/ CO/ 1.50/ CO2/ 2.00/ C2H6/ 3.00/ 
C3H6+H<=>C2H4+CH3                       8.0000E+21   -2.390  11180.0
C3H6+H<=>aC3H5+H2                       1.7300E+05    2.500   2490.0
C3H6+H<=>CH3CCH2+H2                     4.0000E+05    2.500   9790.0
C3H6+O<=>CH2CO+CH3+H                    8.0000E+07    1.650    327.0
C3H6+O<=>C2H3CHO+2H                     4.0000E+07    1.650    327.0
C3H6+O<=>C2H5+HCO                       3.5000E+07    1.650   -972.0
C3H6+O<=>aC3H5+OH                       1.8000E+11    0.700   5880.0
C3H6+O<=>CH3CCH2+OH                     6.0000E+10    0.700   7630.0
C3H6+OH<=>aC3H5+H2O                     3.1000E+06    2.000   -298.0
C3H6+OH<=>CH3CCH2+H2O                   1.1000E+06    2.000   1450.0
C3H6+CH3<=>aC3H5+CH4                    2.2000E+00    3.500   5675.0
C2H3CHO+O<=>C2H3+OH+CO                  3.0000E+13    0.000   3540.0
C2H3CHO+O<=>CH2O+CH2CO                  1.9000E+07    1.800    220.0
iC3H7+H<=>CH3+C2H5                      1.4000E+28   -3.940  15916.0
iC3H7+O<=>CH3CHO+CH3                    9.6000E+13    0.000      0.0
iC3H7+OH<=>C3H6+H2O                     2.4000E+13    0.000      0.0
iC3H7+O2<=>C3H6+HO2                     1.3000E+11    0.000      0.0
iC3H7+HO2<=>CH3CHO+CH3+OH               2.4000E+13    0.000      0.0
iC3H7+CH3<=>CH4+C3H6                    2.2000E+14   -0.680      0.0
nC3H7+H<=>C2H5+CH3                      3.7000E+24   -2.920  12505.0
nC3H7+OH<=>C3H6+H2O                     2.4000E+13    0.000      0.0
nC3H7+O2<=>C3H6+HO2                     9.0000E+10    0.000      0.0
nC3H7+HO2<=>C2H5+OH+CH2O                2.4000E+13    0.000      0.0
nC3H7+CH3<=>CH4+C3H6                    1.1000E+13    0.000      0.0
C4H2+H<=>iC4H3                          1.1000E+30   -4.920  10800.0
iC4H3+H<=>C4H2+H2                       6.0000E+13    0.000      0.0
C4H4+OH<=>iC4H3+H2O                     1.5500E+07    2.000    430.0
C4H5-2<=>iC4H5                          1.5000E+67  -16.890  59100.0
C4H6+H<=>C2H4+C2H3                      1.4600E+30   -4.340  21647.0
C4H6+OH<=>iC4H5+H2O                     3.1000E+06    2.000    430.0
C4H612<=>iC4H5+H                        4.2000E+15    0.000  92600.0
C4H6-2<=>H+C4H5-2                       5.0000E+15    0.000  87300.0
C4H7<=>C4H6+H                           2.4800E+53  -12.300  52000.0
C4H7+O2<=>C4H6+HO2                      1.0000E+11    0.000      0.0
C4H7+HO2<=>CH2O+OH+aC3H5                2.4000E+13    0.000      0.0
C4H81+H<=>C2H4+C2H5                     1.6000E+22   -2.390  11180.0
C4H81+H<=>C3H6+CH3                      3.2000E+22   -2.390  11180.0
C4H81+H<=>C4H7+H2                       6.5000E+05    2.540   6756.0
C4H81+O<=>nC3H7+HCO                     3.3000E+08    1.450   -402.0
C2H4+C2H5<=>pC4H9                       1.5000E+11    0.000   7300.0
pC4H9+OH<=>C4H81+H2O                    2.4000E+13    0.000      0.0
pC4H9+O2<=>C4H81+HO2                    2.7000E+11    0.000      0.0
pC4H9+HO2<=>nC3H7+OH+CH2O               2.4000E+13    0.000      0.0
pC4H9+CH3<=>C4H81+CH4                   1.1000E+13    0.000      0.0
NC12H26 => 3C2H4+2nC3H7                 5.6400E+26   -2.680  88171.0      
NC12H26 => 2C2H4+2pC4H9                 5.1100E+25   -2.510  88117.0     
NC12H26+H => 4C2H4+pC4H9+H2             1.3000E+06    2.540   6756.0  
NC12H26+H => C4H81+2C2H4+pC4H9+H2       1.3000E+06    2.400   4471.0  
NC12H26+H => C3H6+C6H12+nC3H7+H2        1.3000E+06    2.400   4471.0  
NC12H26+H => C5H10+2C2H4+nC3H7+H2       1.3000E+06    2.400   4471.0  
NC12H26+H => C6H12+C2H4+pC4H9+H2        2.6000E+06    2.400   4471.0  
NC12H26+CH3=>4C2H4+pC4H9+CH4            1.8100E+00    3.650   7153.0
NC12H26+CH3=>C4H81+2C2H4+pC4H9+CH4      3.0000E+00    3.460   5480.0
NC12H26+CH3=>C3H6+C6H12+nC3H7+CH4       3.0000E+00    3.460   5480.0
NC12H26+CH3=>C5H10+2C2H4+nC3H7+CH4      3.0000E+00    3.460   5480.0
NC12H26+CH3=>C6H12+C2H4+pC4H9+CH4       6.0000E+00    3.460   5480.0
NC12H26+O=>4C2H4+pC4H9+OH               1.9000E+05    2.680   3716.0
NC12H26+O=>C4H81+2C2H4+pC4H9+OH         4.7600E+04    2.710   2106.0
NC12H26+O=>C3H6+C6H12+nC3H7+OH          4.7600E+04    2.710   2106.0
NC12H26+O=>C5H10+2C2H4+nC3H7+OH         4.7600E+04    2.710   2106.0
NC12H26+O=>C6H12+C2H4+pC4H9+OH          9.5200E+04    2.710   2106.0
NC12H26+OH=>4C2H4+pC4H9+H2O             1.4000E+03    2.660    527.0
NC12H26+OH=>C4H81+2C2H4+pC4H9+H2O       2.7000E+04    2.390    393.0
NC12H26+OH=>C3H6+C6H12+nC3H7+H2O        2.7000E+04    2.390    393.0
NC12H26+OH=>C5H10+2C2H4+nC3H7+H2O       2.7000E+04    2.390    393.0
NC12H26+OH=>C6H12+C2H4+pC4H9+H2O        5.4000E+04    2.390    393.0
C6H12+H(+M) = C3H6+nC3H7(+M)            1.3300E+13    0.000   1559.8 
                       LOW  / 8.70E+42  -7.50   4721.8      /
                       TROE / 1.000  1000.0   645.4  6844.3 /
    H2/2/ H2O/6/ CH4/2/ CO/1.5/ CO2/2/ C2H6/3/
C6H12+H = C2H4+pC4H9                    8.0000E+21   -2.390  11180.0 
C6H12+H = C6H11+H2                      6.5000E+05    2.540   6756.0 
C5H10+H(+M) = C3H6+C2H5(+M)             1.3300E+13    0.000   1559.8 
                       LOW  / 8.70E+42  -7.50   4721.8      /
                       TROE / 1.000  1000.0   645.4  6844.3 /
    H2/2/ H2O/6/ CH4/2/ CO/1.5/ CO2/2/ C2H6/3/ 
C5H10+H = C2H4+nC3H7                    8.0000E+21   -2.390  11180.0
C5H10+H = C2H4+aC3H5+H2                 6.5000E+05    2.540   6756.0
C6H11+H = CH3+C2H4+aC3H5                2.0000E+21   -2.000  11000.0 
C6H11+HO2 => CH2O+OH+aC3H5+C2H4         2.4000E+13    0.000      0.0 
C6H12+O = C2H4+nC3H7+HCO                3.3000E+08    1.450   -402.0 
C5H10+O = pC4H9+HCO                     3.3000E+08    1.450   -402.0 
END



\\
\\
\\  This is the therm file
\\
\\
THERMO
   298.000  1000.000  5000.000
N2                121286N   2               G  0300.00   5000.00  1000.00      1
 0.02926640E+02 0.14879768E-02-0.05684760E-05 0.10097038E-09-0.06753351E-13    2
-0.09227977E+04 0.05980528E+02 0.03298677E+02 0.14082404E-02-0.03963222E-04    3
 0.05641515E-07-0.02444854E-10-0.10208999E+04 0.03950372E+02                   4
AR                120186AR  1               G  0300.00   5000.00  1000.00      1
 0.02500000E+02 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
-0.07453750E+04 0.04366000E+02 0.02500000E+02 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00-0.07453750E+04 0.04366000E+02                   4
HE                L10/90HE  1    0    0    0G   200.000  6000.000 1000.        1
 2.50000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
-7.45375000E+02 9.28723974E-01 2.50000000E+00 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00-7.45375000E+02 9.28723974E-01 0.00000000E+00    4
NE                L10/92NE  1    0    0    0G   200.000  6000.000   1000.00    1
 0.25000000E+01 0.             0.             0.             0.                2
-0.74537500E+03 0.33553227E+01 0.25000000E+01 0.             0.                3
 0.             0.            -0.74537498E+03 0.33553227E+01 0.00000000E+00    4
O                 L 1/90O   1   00   00   00G   200.000  3500.000  1000.000    1
 2.56942078E+00-8.59741137E-05 4.19484589E-08-1.00177799E-11 1.22833691E-15    2
 2.92175791E+04 4.78433864E+00 3.16826710E+00-3.27931884E-03 6.64306396E-06    3
-6.12806624E-09 2.11265971E-12 2.91222592E+04 2.05193346E+00 6.72540300E+03    4
O2                TPIS89O   2   00   00   00G   200.000  3500.000  1000.000    1
 3.28253784E+00 1.48308754E-03-7.57966669E-07 2.09470555E-10-2.16717794E-14    2
-1.08845772E+03 5.45323129E+00 3.78245636E+00-2.99673416E-03 9.84730201E-06    3
-9.68129509E-09 3.24372837E-12-1.06394356E+03 3.65767573E+00 8.68010400E+03    4
H                 L 7/88H   1   00   00   00G   200.000  3500.000   1000.00    1
 2.50000001E+00-2.30842973E-11 1.61561948E-14-4.73515235E-18 4.98197357E-22    2
 2.54736599E+04-4.46682914E-01 2.50000000E+00 7.05332819E-13-1.99591964E-15    3
 2.30081632E-18-9.27732332E-22 2.54736599E+04-4.46682853E-01 6.19742800E+03    4
H2                TPIS78H   2   00   00   00G   200.000  3500.000   1000.00    1
 3.33727920E+00-4.94024731E-05 4.99456778E-07-1.79566394E-10 2.00255376E-14    2
-9.50158922E+02-3.20502331E+00 2.34433112E+00 7.98052075E-03-1.94781510E-05    3
 2.01572094E-08-7.37611761E-12-9.17935173E+02 6.83010238E-01 8.46810200E+03    4
OH                S 9/01O   1H   1    0    0G   200.000  6000.000   1000.00    1
 2.86472886E+00 1.05650448E-03-2.59082758E-07 3.05218674E-11-1.33195876E-15    2
 3.71885774E+03 5.70164073E+00 4.12530561E+00-3.22544939E-03 6.52764691E-06    3
-5.79853643E-09 2.06237379E-12 3.38153812E+03-6.90432960E-01 4.51532273E+03    4
H2O               L 8/89H   2O   1   00   00G   200.000  3500.000  1000.000    1
 3.03399249E+00 2.17691804E-03-1.64072518E-07-9.70419870E-11 1.68200992E-14    2
-3.00042971E+04 4.96677010E+00 4.19864056E+00-2.03643410E-03 6.52040211E-06    3
-5.48797062E-09 1.77197817E-12-3.02937267E+04-8.49032208E-01 9.90409200E+03    4
HO2               L 5/89H   1O   2   00   00G   200.000  3500.000  1000.000    1
 4.01721090E+00 2.23982013E-03-6.33658150E-07 1.14246370E-10-1.07908535E-14    2
 1.11856713E+02 3.78510215E+00 4.30179801E+00-4.74912051E-03 2.11582891E-05    3
-2.42763894E-08 9.29225124E-12 2.94808040E+02 3.71666245E+00 1.00021620E+04    4
H2O2              L 7/88H   2O   2   00   00G   200.000  3500.000  1000.000    1
 4.16500285E+00 4.90831694E-03-1.90139225E-06 3.71185986E-10-2.87908305E-14    2
-1.78617877E+04 2.91615662E+00 4.27611269E+00-5.42822417E-04 1.67335701E-05    3
-2.15770813E-08 8.62454363E-12-1.77025821E+04 3.43505074E+00 1.11588350E+04    4
C                 L11/88C   1   00   00   00G   200.000  3500.000  1000.000    1
 2.49266888E+00 4.79889284E-05-7.24335020E-08 3.74291029E-11-4.87277893E-15    2
 8.54512953E+04 4.80150373E+00 2.55423955E+00-3.21537724E-04 7.33792245E-07    3
-7.32234889E-10 2.66521446E-13 8.54438832E+04 4.53130848E+00 6.53589500E+03    4
CH                TPIS79C   1H   1   00   00G   200.000  3500.000  1000.000    1
 2.87846473E+00 9.70913681E-04 1.44445655E-07-1.30687849E-10 1.76079383E-14    2
 7.10124364E+04 5.48497999E+00 3.48981665E+00 3.23835541E-04-1.68899065E-06    3
 3.16217327E-09-1.40609067E-12 7.07972934E+04 2.08401108E+00 8.62500000E+03    4
CH2               L S/93C   1H   2   00   00G   200.000  3500.000  1000.000    1
 2.87410113E+00 3.65639292E-03-1.40894597E-06 2.60179549E-10-1.87727567E-14    2
 4.62636040E+04 6.17119324E+00 3.76267867E+00 9.68872143E-04 2.79489841E-06    3
-3.85091153E-09 1.68741719E-12 4.60040401E+04 1.56253185E+00 1.00274170E+04    4
CH2*              L S/93C   1H   2   00   00G   200.000  3500.000  1000.000    1
 2.29203842E+00 4.65588637E-03-2.01191947E-06 4.17906000E-10-3.39716365E-14    2
 5.09259997E+04 8.62650169E+00 4.19860411E+00-2.36661419E-03 8.23296220E-06    3
-6.68815981E-09 1.94314737E-12 5.04968163E+04-7.69118967E-01 9.93967200E+03    4
CH3               L11/89C   1H   3   00   00G   200.000  3500.000  1000.000    1
 2.28571772E+00 7.23990037E-03-2.98714348E-06 5.95684644E-10-4.67154394E-14    2
 1.67755843E+04 8.48007179E+00 3.67359040E+00 2.01095175E-03 5.73021856E-06    3
-6.87117425E-09 2.54385734E-12 1.64449988E+04 1.60456433E+00 1.03663400E+04    4
CH4               L 8/88C   1H   4   00   00G   200.000  3500.000  1000.000    1
 7.48514950E-02 1.33909467E-02-5.73285809E-06 1.22292535E-09-1.01815230E-13    2
-9.46834459E+03 1.84373180E+01 5.14987613E+00-1.36709788E-02 4.91800599E-05    3
-4.84743026E-08 1.66693956E-11-1.02466476E+04-4.64130376E+00 1.00161980E+04    4
CO                TPIS79C   1O   1   00   00G   200.000  3500.000  1000.000    1
 2.71518561E+00 2.06252743E-03-9.98825771E-07 2.30053008E-10-2.03647716E-14    2
-1.41518724E+04 7.81868772E+00 3.57953347E+00-6.10353680E-04 1.01681433E-06    3
 9.07005884E-10-9.04424499E-13-1.43440860E+04 3.50840928E+00 8.67100000E+03    4
CO2               L 7/88C   1O   2   00   00G   200.000  3500.000  1000.000    1
 3.85746029E+00 4.41437026E-03-2.21481404E-06 5.23490188E-10-4.72084164E-14    2
-4.87591660E+04 2.27163806E+00 2.35677352E+00 8.98459677E-03-7.12356269E-06    3
 2.45919022E-09-1.43699548E-13-4.83719697E+04 9.90105222E+00 9.36546900E+03    4
HCO               L12/89H   1C   1O   1   00G   200.000  3500.000  1000.000    1
 2.77217438E+00 4.95695526E-03-2.48445613E-06 5.89161778E-10-5.33508711E-14    2
 4.01191815E+03 9.79834492E+00 4.22118584E+00-3.24392532E-03 1.37799446E-05    3
-1.33144093E-08 4.33768865E-12 3.83956496E+03 3.39437243E+00 9.98945000E+03    4
CH2O              L 8/88H   2C   1O   1   00G   200.000  3500.000  1000.000    1
 1.76069008E+00 9.20000082E-03-4.42258813E-06 1.00641212E-09-8.83855640E-14    2
-1.39958323E+04 1.36563230E+01 4.79372315E+00-9.90833369E-03 3.73220008E-05    3
-3.79285261E-08 1.31772652E-11-1.43089567E+04 6.02812900E-01 1.00197170E+04    4
CH2OH             IU2/03C   1H   3O   1   00G   200.000  6000.00               1
 5.09314370E+00 5.94761260E-03-2.06497460E-06 3.23008173E-10-1.88125902E-14    2
-4.03409640E+03-1.84691493E+00 4.47834367E+00-1.35070310E-03 2.78484980E-05    3
-3.64869060E-08 1.47907450E-11-3.50072890E+03 3.30913500E+00-2.04462770E+03    4
CH3O              IU1/03C   1H   3O   1     G   200.000  6000.00               1
 4.75779238E+00 7.44142474E-03-2.69705176E-06 4.38090504E-10-2.63537098E-14    2
 3.78111940E+02-1.96680028E+00 3.71180502E+00-2.80463306E-03 3.76550971E-05    3
-4.73072089E-08 1.86588420E-11 1.29569760E+03 6.57240864E+00 2.52571660E+03    4
CH3OH             L 8/88C   1H   4O   1   00G   200.000  3500.000  1000.000    1
 1.78970791E+00 1.40938292E-02-6.36500835E-06 1.38171085E-09-1.17060220E-13    2
-2.53748747E+04 1.45023623E+01 5.71539582E+00-1.52309129E-02 6.52441155E-05    3
-7.10806889E-08 2.61352698E-11-2.56427656E+04-1.50409823E+00 1.14352770E+04    4
C2H               L 1/91C   2H   1   00   00G   200.000  3500.000  1000.000    1
 3.16780652E+00 4.75221902E-03-1.83787077E-06 3.04190252E-10-1.77232770E-14    2
 6.71210650E+04 6.63589475E+00 2.88965733E+00 1.34099611E-02-2.84769501E-05    3
 2.94791045E-08-1.09331511E-11 6.68393932E+04 6.22296438E+00 1.04544720E+04    4
C2H2              L 1/91C   2H   2   00   00G   200.000  3500.000  1000.000    1
 4.14756964E+00 5.96166664E-03-2.37294852E-06 4.67412171E-10-3.61235213E-14    2
 2.59359992E+04-1.23028121E+00 8.08681094E-01 2.33615629E-02-3.55171815E-05    3
 2.80152437E-08-8.50072974E-12 2.64289807E+04 1.39397051E+01 1.00058390E+04    4
H2CC              L12/89H   2C   2    0    0G   200.000  6000.000  1000.000    1
 0.42780340E+01 0.47562804E-02-0.16301009E-05 0.25462806E-09-0.14886379E-13    2
 0.48316688E+05 0.64023701E+00 0.32815483E+01 0.69764791E-02-0.23855244E-05    3
-0.12104432E-08 0.98189545E-12 0.48621794E+05 0.59203910E+01 0.49887266E+05    4
C2H3              L 2/92C   2H   3   00   00G   200.000  3500.000  1000.000    1
 3.01672400E+00 1.03302292E-02-4.68082349E-06 1.01763288E-09-8.62607041E-14    2
 3.46128739E+04 7.78732378E+00 3.21246645E+00 1.51479162E-03 2.59209412E-05    3
-3.57657847E-08 1.47150873E-11 3.48598468E+04 8.51054025E+00 1.05750490E+04    4
C2H4              L 1/91C   2H   4   00   00G   200.000  3500.000  1000.000    1
 2.03611116E+00 1.46454151E-02-6.71077915E-06 1.47222923E-09-1.25706061E-13    2
 4.93988614E+03 1.03053693E+01 3.95920148E+00-7.57052247E-03 5.70990292E-05    3
-6.91588753E-08 2.69884373E-11 5.08977593E+03 4.09733096E+00 1.05186890E+04    4
C2H5              L12/92C   2H   5   00   00G   200.000  3500.000  1000.000    1
 1.95465642E+00 1.73972722E-02-7.98206668E-06 1.75217689E-09-1.49641576E-13    2
 1.28575200E+04 1.34624343E+01 4.30646568E+00-4.18658892E-03 4.97142807E-05    3
-5.99126606E-08 2.30509004E-11 1.28416265E+04 4.70720924E+00 1.21852440E+04    4
C2H6              L 8/88C   2H   6   00   00G   200.000  3500.000  1000.000    1
 1.07188150E+00 2.16852677E-02-1.00256067E-05 2.21412001E-09-1.90002890E-13    2
-1.14263932E+04 1.51156107E+01 4.29142492E+00-5.50154270E-03 5.99438288E-05    3
-7.08466285E-08 2.68685771E-11-1.15222055E+04 2.66682316E+00 1.18915940E+04    4
CH2CO             D05/90C   2H   2O   1   00G   200.000  3500.000  1000.000    1
 4.51129732E+00 9.00359745E-03-4.16939635E-06 9.23345882E-10-7.94838201E-14    2
-7.77850000E+03 6.32247205E-01 2.13583630E+00 1.81188721E-02-1.73947474E-05    3
 9.34397568E-09-2.01457615E-12-7.27000000E+03 1.22156480E+01 1.17977430E+04    4
CH2CHO            D05/83O   1H   3C   2    0G   300.000  5000.000              1
 0.59756699E+01 0.81305914E-02-0.27436245E-05 0.40703041E-09-0.21760171E-13    2
-0.96950000E+03-0.50320879E+01 0.34090624E+01 0.10738574E-01 0.18914925E-05    3
-0.71585831E-08 0.28673851E-11 0.62000000E+02 0.95714535E+01 0.30474436E+04    4
CH2OCH            A12/04C   2H   3O   1    0G   298.150  3000.000   500.0      1
 0.44994054E+01 0.11552625E-01-0.48144129E-05 0.89234919E-09-0.56870585E-13    2
 0.17473963E+05 0.33925515E+00-0.38396084E+00 0.23879038E-01-0.12467587E-04    3
-0.17686411E-08 0.28142438E-11 0.18836203E+05 0.25741745E+02                   4
CH2OCH2           T 6/92C   2H   4O   1    0G   298.150  3000.0    1000.0      1
 0.54887641E+01 0.12046190E-01-0.43336931E-05 0.70028311E-09-0.41949088E-13    2
-0.91804251E+04-0.70799605E+01 0.37590532E+01-0.94412180E-02 0.80309721E-04    3
-0.10080788E-06 0.40039921E-10-0.75608143E+04 0.78497475E+01-0.63304657E+04    4
CH3CO             T 9/92C   2H   3O   1    0G   200.000  6000.0    1000.0      1
 0.59447731E+01 0.78667205E-02-0.28865882E-05 0.47270875E-09-0.28599861E-13    2
-0.37873075E+04-0.50136751E+01 0.41634257E+01-0.23261610E-03 0.34267820E-04    3
-0.44105227E-07 0.17275612E-10-0.26574529E+04 0.73468280E+01-0.12027167E+04    4
CH3CHO            L 8/88C   2H   4O   1    0G   200.000  6000.0    1000.0      1
 0.54041108E+01 0.11723059E-01-0.42263137E-05 0.68372451E-09-0.40984863E-13    2
-0.22593122E+05-0.34807917E+01 0.47294595E+01-0.31932858E-02 0.47534921E-04    3
-0.57458611E-07 0.21931112E-10-0.21572878E+05 0.41030159E+01-0.19987949E+05    4
HCCO              SRIC91H   1C   2O   1     G  0300.00   4000.00  1000.00      1
 0.56282058E+01 0.40853401E-02-0.15934547E-05 0.28626052E-09-0.19407832E-13    2
 0.19327215E+05-0.39302595E+01 0.22517214E+01 0.17655021E-01-0.23729101E-04    3
 0.17275759E-07-0.50664811E-11 0.20059449E+05 0.12490417E+02                   4
C2O               RUS 79C   2O   1    0    0G   200.000  6000.000              1
 0.51512722E+01 0.23726722E-02-0.76135971E-06 0.11706415E-09-0.70257804E-14    2
 0.33241888E+05-0.22183135E+01 0.28648610E+01 0.11990216E-01-0.18362448E-04    3
 0.15769739E-07-0.53897452E-11 0.33749932E+05 0.88867772E+01 0.35003406E+05    4
C3H2              T12/00C   3H   2    0    0G   200.000  6000.000              1
 0.73481207E+01 0.44476404E-02-0.12610332E-05 0.78131814E-10 0.13216298E-13    2
 0.62551656E+05-0.91040211E+01 0.45094776E+01 0.17438605E-01-0.24516321E-04    3
 0.18993967E-07-0.57996520E-11 0.63080191E+05 0.42892461E+01                   4
C3H2-2            S 4/01C   3H   2    0    0G   200.000  3000.000              1
 7.47247827E+00 4.57765160E-03-1.56482125E-06 2.43991965E-10-1.42462924E-14    2
 8.83321441E+04-1.27113314E+01 3.74356467E+00 2.51955211E-02-4.62608277E-05    3
 4.34360520E-08-1.53992558E-11 8.89297787E+04 4.22612394E+00 9.08356403E+04    4
cC3H2             121686C   3H   2          G  0300.00   5000.00  1000.00      1
 0.06530853E+02 0.05870316E-01-0.01720777E-04 0.02127498E-08-0.08291910E-13    2
 0.05115214E+06-0.01122728E+03 0.02691077E+02 0.01480366E+00-0.03250551E-04    3
-0.08644363E-07 0.05284878E-10 0.05219072E+06 0.08757391E+02                   4
C3H3              T 5/97C   3H   3    0    0G   200.000  6000.000              1
 7.14221880E+00 7.61902005E-03-2.67459950E-06 4.24914801E-10-2.51475415E-14    2
 3.89087427E+04-1.25848436E+01 1.35110927E+00 3.27411223E-02-4.73827135E-05    3
 3.76309808E-08-1.18540923E-11 4.01057783E+04 1.52058924E+01 4.16139977E+04    4
aC3H4             L 8/89C   3H   4    0    0G   200.000  6000.000              1
 0.63168722E+01 0.11133728E-01-0.39629378E-05 0.63564238E-09-0.37875540E-13    2
 0.20117495E+05-0.10995766E+02 0.26130445E+01 0.12122575E-01 0.18539880E-04    3
-0.34525149E-07 0.15335079E-10 0.21541567E+05 0.10226139E+02 0.22962267E+05    4
pC3H4             T 2/90H   4C   3    0    0G   200.000  6000.000              1
 0.60252400E+01 0.11336542E-01-0.40223391E-05 0.64376063E-09-0.38299635E-13    2
 0.19620942E+05-0.86043785E+01 0.26803869E+01 0.15799651E-01 0.25070596E-05    3
-0.13657623E-07 0.66154285E-11 0.20802374E+05 0.98769351E+01 0.22302059E+05    4
cC3H4             T12/81C   3H   4    0    0G   300.000  5000.000              1
 0.66999931E+01 0.10357372E-01-0.34551167E-05 0.50652949E-09-0.26682276E-13    2
 0.30199051E+05-0.13378770E+02-0.24621047E-01 0.23197215E-01-0.18474357E-05    3
-0.15927593E-07 0.86846155E-11 0.32334137E+05 0.22729762E+02 0.3332728 E+05    4
C3H8              P11/94C   3H   8    0    0G   300.000  3000.000              1
 0.75244152E+01 0.18898282E-01-0.62921041E-05 0.92161457E-09-0.48684478E-13    2
-0.16564394E+05-0.17838375E+02 0.92851093E+00 0.26460566E-01 0.60332446E-05    3
-0.21914953E-07 0.94961544E-11-0.14057907E+05 0.19225538E+02                   4
nC3H7             P11/94C   3H   7    0    0G   300.000  3000.000              1
 0.77097479E+01 0.16031485E-01-0.52720238E-05 0.75888352E-09-0.38862719E-13    2
 0.79762236E+04-0.15515297E+02 0.10491173E+01 0.26008973E-01 0.23542516E-05    3
-0.19595132E-07 0.93720207E-11 0.10312346E+05 0.21136034E+02                   4
iC3H7             P11/94C   3H   7    0    0G   300.000  3000.000              1
 0.65192741E+01 0.17220104E-01-0.57364217E-05 0.84130732E-09-0.44565913E-13    2
 0.73227193E+04-0.90830215E+01 0.14449199E+01 0.20999112E-01 0.77036222E-05    3
-0.18476253E-07 0.71282962E-11 0.94223724E+04 0.20116317E+02                   4
C3H6              120186C   3H   6          G  0300.00   5000.00  1000.00      1
 0.06732257E+02 0.01490834E+00-0.04949899E-04 0.07212022E-08-0.03766204E-12    2
-0.09235703E+04-0.01331335E+03 0.01493307E+02 0.02092518E+00 0.04486794E-04    3
-0.01668912E-06 0.07158146E-10 0.01074826E+05 0.01614534E+03                   4
CH2CHCO           T05/99C   3H   3O   1    0G   200.000  6000.0    1000.0      1
 6.95842227E+00 1.07193211E-02-3.85218494E-06 6.22009064E-10-3.72401640E-14    2
 5.64826498E+03-1.14745786E+01 3.21169467E+00 1.18422105E-02 1.67462582E-05    3
-3.06947176E-08 1.33048816E-11 7.12815750E+03 1.00881663E+01 8.70564832E+03    4
CH2CHCO           USC/07C   3H   3O   1    0G   300.000  5000.000              1
 0.73338666E+01 0.11401899E-01-0.45696443E-05 0.79430967E-09-0.44163078E-13    2
 0.83094941E+04-0.11019943E+02 0.23135836E+01 0.28253129E-01-0.25737754E-04    3
 0.12222654E-07-0.21353429E-11 0.95213496E+04 0.14105129E+02                   4
CH3CH2CHO         USC/07C   3H   6O   1    0G   300.000  5000.000              1
 0.62637410E+01 0.19976260E-01-0.76195147E-05 0.11687118E-08-0.41959993E-13    2
-0.25885953E+05-0.57786498E+01 0.27255676E+01 0.23236005E-01 0.29740656E-05    3
-0.16613415E-07 0.74250103E-11-0.24556711E+05 0.14166277E+02                   4
CH3COCH3          T 5/92C   3H   6O   1    0G   200.000  6000.000  1000.0      1
 0.72975991E+01 0.17566207E-01-0.63170456E-05 0.10203086E-08-0.61094016E-13    2
-0.29817680E+05-0.12756981E+02 0.55557943E+01-0.28365428E-02 0.70568945E-04    3
-0.87810488E-07 0.34028266E-10-0.28113337E+05 0.23226600E+01-0.26116945E+05    4
C2H3CHO           USC/07C   3H   4O   1    0G   300.000  5000.000              1
 0.58111868E+01 0.17114256E-01-0.74834161E-05 0.14252249E-08-0.91746841E-13    2
-0.10784054E+05-0.48588004E+01 0.12713498E+01 0.26231054E-01-0.92912305E-05    3
-0.47837272E-08 0.33480543E-11-0.93357344E+04 0.19498077E+02                   4
aC3H5             PD5/98C   3H   5    0    0G   300.000  3000.000              1
 0.65007877E+01 0.14324731E-01-0.56781632E-05 0.11080801E-08-0.90363887E-13    2
 0.17482449E+05-0.11243050E+02 0.13631835E+01 0.19813821E-01 0.12497060E-04    3
-0.33355555E-07 0.15846571E-10 0.19245629E+05 0.17173214E+02                   4
CH3CCH2           PD5/98C   3H   5    0    0G   300.000  3000.000              1
 0.54255528E+01 0.15511072E-01-0.56678350E-05 0.79224388E-09-0.16878034E-13    2
 0.27843027E+05-0.33527184E+01 0.17329209E+01 0.22394620E-01-0.51490611E-05    3
-0.67596466E-08 0.38253211E-11 0.29040498E+05 0.16568878E+02                   4
CH3CHCH           PD5/98C   3H   5    0    0G   300.000  3000.000              1
 0.53725281E+01 0.15780509E-01-0.59922850E-05 0.93089664E-09-0.36550966E-13    2
 0.29614760E+05-0.34186478E+01 0.91372931E+00 0.26432343E-01-0.11758950E-04    3
-0.23035678E-08 0.27715488E-11 0.30916867E+05 0.19989269E+02                   4
C4H               P 1/93C   4H   1    0    0G   300.000  3000.000              1
 0.77697593E+01 0.49829976E-02-0.17628546E-05 0.28144284E-09-0.16689869E-13    2
 0.94345900E+05-0.14165274E+02 0.13186295E+01 0.38582956E-01-0.71385623E-04    3
 0.65356359E-07-0.22617666E-10 0.95456106E+05 0.15567583E+02                   4
C4H2              D11/99C   4H   2    0    0G   300.000  3000.000              1
 0.91576328E+01 0.55430518E-02-0.13591604E-05 0.18780075E-10 0.23189536E-13    2
 0.52588039E+05-0.23711460E+02 0.10543978E+01 0.41626960E-01-0.65871784E-04    3
 0.53257075E-07-0.16683162E-10 0.54185211E+05 0.14866591E+02                   4
nC4H3             USC/07C   4H   3O   0    0G   300.000  5000.000              1
 0.78045716E+01 0.10712364E-01-0.41939124E-05 0.70446277E-09-0.36271326E-13    2
 0.62987805E+05-0.14129741E+02 0.81667686E+00 0.38716201E-01-0.48045651E-04    3
 0.32066808E-07-0.85628215E-11 0.64455754E+05 0.19740503E+02                   4
iC4H3             USC/07C   4H   3O   0    0G   300.000  5000.000              1
 0.76538548E+01 0.11204055E-01-0.46401342E-05 0.86786639E-09-0.57430562E-13    2
 0.57954363E+05-0.11756476E+02 0.37221482E+01 0.25957543E-01-0.26356343E-04    3
 0.15508920E-07-0.38040565E-11 0.58837121E+05 0.75637245E+01                   4
H2C4O             USC/07C   4H   2O   1    0G   300.000  5000.000              1
 0.84292183E+01 0.10502701E-01-0.42066836E-05 0.71184902E-09-0.35796602E-13    2
 0.22907807E+05-0.16511997E+02 0.31811900E+01 0.29840752E-01-0.32832409E-04    3
 0.20631813E-07-0.54200598E-11 0.24125576E+05 0.94210100E+01                   4
C4H4              USC/07C   4H   4O   0    0G   300.000  5000.000              1
 0.72539601E+01 0.13914094E-01-0.52932214E-05 0.83480450E-09-0.35197882E-13    2
 0.31766016E+05-0.12629521E+02 0.58857048E+00 0.36546685E-01-0.34106968E-04    3
 0.16652619E-07-0.30064623E-11 0.33359492E+05 0.20657881E+02                   4
nC4H5             USC/07C   4H   5O   0    0G   300.000  5000.000              1
 0.74087291E+01 0.17752748E-01-0.75601506E-05 0.14203795E-08-0.91100182E-13    2
 0.40438762E+05-0.13150027E+02 0.22611290E+00 0.36742371E-01-0.22120474E-04    3
 0.14390138E-08 0.26435809E-11 0.42428410E+05 0.24066401E+02                   4
iC4H5             USC/07C   4H   5O   0    0G   300.000  5000.000              1
 0.69646029E+01 0.18274333E-01-0.78133735E-05 0.15292154E-08-0.10920493E-12    2
 0.34725098E+05-0.10649321E+02 0.11308105E+00 0.40950615E-01-0.35413581E-04    3
 0.15530969E-07-0.23355122E-11 0.36383371E+05 0.23692457E+02                   4
C4H5-2            H6W/94C   4H   5    0    0G   300.000  3000.000              1
 1.45381710E+01-8.56770560E-03 2.35595240E-05-1.36763790E-08 2.44369270E-12    2
 3.32590950E+04-4.53694970E+01 2.96962800E+00 2.44422450E-02-9.12514240E-06    3
-4.24668710E-18 1.63047280E-21 3.55033160E+04 1.20360510E+01 3.73930550E+04    4
c-C4H5            PUPM3 C   4H   5    0    0G   300.000  3000.000              1
 0.67467155E+01 0.17283000E-01-0.65168579E-05 0.98917574E-09-0.34604908E-13    2
 0.32808359E+05-0.12912880E+02-0.26397593E+01 0.41549157E-01-0.21920954E-04    3
-0.46559014E-08 0.61348890E-11 0.35373828E+05 0.35701797E+02                   4
C4H7              USC/07C   4H   7O   0    0G   300.000  5000.000              1
 0.70134835E+01 0.22634558E-01-0.92545470E-05 0.16807927E-08-0.10408617E-12    2
 0.20955008E+05-0.88893080E+01 0.74449432E+00 0.39678857E-01-0.22898086E-04    3
 0.21352973E-08 0.23096375E-11 0.22653328E+05 0.23437878E+02                   4
C4H6              H6W/94C   4H   6    0    0G   300.000  3000.000              1
 0.88673134E+01 0.14918670E-01-0.31548716E-05-0.41841330E-09 0.15761258E-12    2
 0.91338516E+04-0.23328171E+02 0.11284465E+00 0.34369022E-01-0.11107392E-04    3
-0.92106660E-08 0.62065179E-11 0.11802270E+05 0.23089996E+02                   4
C4H612            A 8/83C   4H   6    0    0G   300.     3000.     1000.0      1
  0.1781557E+02 -0.4257502E-02  0.1051185E-04 -0.4473844E-08  0.5848138E-12    2
  0.1267342E+05 -0.6982662E+02  0.1023467E+01  0.3495919E-01 -0.2200905E-04    3
  0.6942272E-08 -0.7879187E-12  0.1811799E+05  0.1975066E+02  0.1950807E+05    4
C4H6-2            A 8/83C   4H   6    0    0G   300.     3000.     1000.0      1
  9.0338133E+00 8.2124510E-03   7.1753952E-06 -5.8834334E-09  1.0343915E-12    2
  1.4335068E+04 -2.0985762E+01  2.1373338E+00  2.6486229E-02 -9.0568711E-06    3
 -5.5386397E-19 2.1281884E-22   1.5710902E+04  1.3529426E+01  1.7488676E+04    4
C4H10             P11/94C   4H  10    0    0G   300.000  3000.000              1
 0.10526774E+02 0.23590738E-01-0.78522480E-05 0.11448408E-08-0.59827703E-13    2
-0.20479223E+05-0.32198579E+02 0.15685419E+01 0.34652278E-01 0.68168129E-05    3
-0.27995097E-07 0.12307742E-10-0.17129977E+05 0.17908045E+02                   4
iC4H10            P11/94C   4H  10    0    0G   300.000  3000.000              1
 0.10846169E+02 0.23338389E-01-0.77833962E-05 0.11393807E-08-0.59918289E-13    2
-0.21669854E+05-0.35870573E+02 0.54109489E+00 0.37860301E-01 0.55459804E-05    3
-0.30500110E-07 0.14033357E-10-0.17977644E+05 0.21150935E+02                   4
pC4H9             USC/07C   4H   9O   0    0G   300.000  5000.000              1
 0.86822395E+01 0.23691071E-01-0.75948865E-05 0.66427136E-09 0.54845136E-13    2
 0.49644058E+04-0.17891747E+02 0.12087042E+01 0.38297497E-01-0.72660509E-05    3
-0.15428547E-07 0.86859435E-11 0.73221040E+04 0.22169268E+02                   4
sC4H9             P11/94C   4H   9    0    0G   300.000  3000.000              1
 0.94263839E+01 0.21918998E-01-0.72868375E-05 0.10630334E-08-0.55649464E-13    2
 0.31965874E+04-0.22406051E+02 0.69428423E+00 0.33113346E-01 0.62942577E-05    3
-0.27025274E-07 0.11989315E-10 0.64175654E+04 0.26279789E+02                   4
tC4H9             P11/94C   4H   9    0    0G   300.000  3000.000              1
 0.76607261E+01 0.23879414E-01-0.80890353E-05 0.12057521E-08-0.65009814E-13    2
 0.16207623E+04-0.14800281E+02 0.96167553E+00 0.25735856E-01 0.15609033E-04    3
-0.26656519E-07 0.89418010E-11 0.46564412E+04 0.24805366E+02                   4
iC4H9             USC/07C   4H   9O   0    0G   300.000  5000.000              1
 0.84981728E+01 0.24689538E-01-0.86487589E-05 0.10779325E-08-0.64340570E-15    2
 0.44288174E+04-0.18441397E+02 0.97527862E+00 0.41613799E-01-0.14467331E-04    3
-0.93852393E-08 0.68797377E-11 0.66688267E+04 0.21277582E+02                   4
C4H81             T 6/83C   4H   8    0    0G   300.000  5000.000              1
 0.20535841E+01 0.34350507E-01-0.15883197E-04 0.33089662E-08-0.25361045E-12    2
-0.21397231E+04 0.15543201E+02 0.11811380E+01 0.30853380E-01 0.50865247E-05    3
-0.24654888E-07 0.11110193E-10-0.17904004E+04 0.21062469E+02                   4
C4H82             T 6/83C   4H   8    0    0G   300.000  5000.00               1
 0.82797676E+00 0.35864539E-01-0.16634498E-04 0.34732759E-08-0.26657398E-12    2
-0.30521033E+04 0.21342545E+02 0.12594252E+01 0.27808424E-01 0.87013932E-05    3
-0.24402205E-07 0.98977710E-11-0.29647742E+04 0.20501129E+02                   4
iC4H8             T 6/83H   8C   4    0    0G   300.000  5000.0                1
 0.44609470E+01 0.29611487E-01-0.13077129E-04 0.26571934E-08-0.20134713E-12    2
-0.50066758E+04 0.10671549E+01 0.26471405E+01 0.25902957E-01 0.81985354E-05    3
-0.22193259E-07 0.88958580E-11-0.40373069E+04 0.12676388E+02                   4
iC4H7             USC/07C   4H   7O   0    0G   300.000  5000.000              1
 0.71485939E+01 0.22189671E-01-0.84400172E-05 0.13133353E-08-0.51617927E-13    2
 0.12712294E+05-0.12131183E+02-0.10375890E+01 0.45566667E-01-0.30476231E-04    3
 0.71102568E-08 0.99685722E-12 0.14896458E+05 0.29863663E+02                   4
C2H3CHOCH2        A 8/83C   4H   6O   1    0G   300.     3000.     1000.0      1
-4.72093360E+00 3.91413780E-02-6.52872650E-06-7.68209500E-09 2.51473310E-12    2
 1.75352252E+03 5.17190420E+01 7.97985440E-01 3.44034320E-02-1.24598510E-05    3
-5.18062790E-18 1.99359540E-21-6.48927540E+02 2.18896980E+01 1.00654250E+03    4
CH3CHCHCHO        T 5/92C   4H   6O   1    0G   298.150  3000.0    1000.0      1
 1.98794540E+01-2.09130550E-02 4.45360508E-05-2.60374870E-08 4.86836120E-12    2
-1.95278768E+04-6.87200320E+01-1.55577660E+00 4.09640630E-02-1.69868810E-05    3
-6.00928140E-18 2.31368530E-21-1.41394920E+04 3.74707580E+01-1.29340710E+04    4
CH2CHCOCH3        T 3/97C   4H   6O   1    0G   200.000  3000.0    1000.0      1
 1.98794540E+01-2.09130550E-02 4.45360580E-05-2.60374870E-08 4.86836120E-12    2
-1.90786168E+04-6.97265750E+01-1.55577660E+00 4.09640630E-02-1.69868810E-05    3
-6.00928140E-18 2.31368530E-21-1.49447258E+04 3.64642160E+01-1.66079520E+04    4
C4H4O             T03/97C   4H   4O   1    0G   200.000  6000.0    1000.0      1
 9.38935003E+00 1.40291241E-02-5.07755110E-06 8.24137332E-10-4.95319963E-14    2
-8.68241814E+03-2.79162920E+01 8.47469463E-01 1.31773796E-02 5.99735901E-05    3
-9.71562904E-08 4.22733796E-11-5.36785445E+03 2.14945172E+01-4.17166616E+03    4
CH3CHCHCO         USC/07C   4H   5O   1    0G   300.000  5000.000              1
 0.77608204E+01 0.20031804E-01-0.80631016E-05 0.13361392E-08-0.62308408E-13    2
 0.45708291E+04-0.11095638E+02 0.53053460E+01 0.15749373E-01 0.21623913E-04    3
-0.36607769E-07 0.14932489E-10 0.57588633E+04 0.42043533E+01                   4
CH2CHCHCHO        USC/07C   4H   5O   1    0G   300.000  5000.000              1
 0.83010607E+01 0.19945331E-01-0.82903771E-05 0.15100753E-08-0.91581155E-13    2
 0.15788387E+03-0.16910566E+02 0.12108673E+01 0.35205878E-01-0.10939090E-04    3
-0.11720642E-07 0.76174908E-11 0.22665703E+04 0.20613544E+02                   4
CH2CHCH2CHO       T 5/92C   4H   6O   1    0G   298.150  3000.0    1000.0      1
 1.98794540E+01-2.09130550E-02 4.45360508E-05-2.60374870E-08 4.86836120E-12    2
-1.58539966E+04-6.71095639E+01-1.55577660E+00 4.09640630E-02-1.69868810E-05    3
-6.00928140E-18 2.31368530E-21-1.04656118E+04 3.90812260E+01-1.29340710E+04    4
C4H6O25           T 3/97C   4H   6O   1    0G   200.000  5000.000  1000.0      1
 8.60658242E+00 2.08310051E-02-8.42229481E-06 1.56717640E-09-1.09391202E-13    2
-1.76177415E+04-2.32464750E+01 2.67053463E+00 4.92586420E-03 8.86967406E-05    3
-1.26219194E-07 5.23991321E-11-1.46572472E+04 1.45722395E+01-1.30831522E+04    4
C4H6O23           T 3/97C   4H   6O   1    0G   200.000  5000.000  1000.0      1
 8.60658242E+00 2.08310051E-02-8.42229481E-06 1.56717640E-09-1.09391202E-13    2
-1.32392815E+04-2.32464750E+01 2.67053463E+00 4.92586420E-03 8.86967406E-05    3
-1.26219194E-07 5.23991321E-11-1.02787872E+04 1.45722395E+01-1.30831522E+04    4
sC4H9             T07/95C   4H   9    0    0G   200.000  6000.000  1000.0      1
 0.88057265E+01 0.23630381E-01-0.84564737E-05 0.13612584E-08-0.81313232E-13    2
 0.37941169E+04-0.19996770E+02 0.46457042E+01 0.79313214E-02 0.70027013E-04    3
-0.95973349E-07 0.38628890E-10 0.62341181E+04 0.79642808E+01 0.84190169E+04    4
C5H2               20587C   5H   2          G  0300.00   5000.00  1000.00      1
 0.01132917E+03 0.07424057E-01-0.02628189E-04 0.04082541E-08-0.02301333E-12    2
 0.07878706E+06-0.03617117E+03 0.03062322E+02 0.02709998E+00-0.01009170E-03    3
-0.01272745E-06 0.09167219E-10 0.08114969E+06 0.07071078E+02                   4
C5H3               20387C   5H   3          G  0300.00   5000.00  1000.00      1
 0.01078762E+03 0.09539619E-01-0.03206745E-04 0.04733323E-08-0.02512135E-12    2
 0.06392904E+06-0.03005444E+03 0.04328720E+02 0.02352480E+00-0.05856723E-04    3
-0.01215449E-06 0.07726478E-10 0.06588531E+06 0.04173259E+02                   4
C5H6              T 1/90C   5H 6      0    0G   200.000  6000.000              1
 0.99757848E+01 0.18905543E-01-0.68411461E-05 0.11099340E-08-0.66680236E-13    2
 0.11081693E+05-0.32209454E+02 0.86108957E+00 0.14804031E-01 0.72108895E-04    3
-0.11338055E-06 0.48689972E-10 0.14801755E+05 0.21353453E+02 0.16152485E+05    4
C5H5              T12/89C   5H 5      0    0G   300.00   2000.000 1000.00      1
 0.74743938E+01 0.16012733E-01-0.64823148E-08-0.35819703E-08 0.92365071E-12    2
 2.80860000E+04-0.16133000E+02 0.98349822E+00 0.33651476E-01-0.11054181E-06    3
-0.36743394E-07 0.23141184E-10 2.96260000E+04 0.16585519E+02                   4
cC5H8             T03/97C   5H   8O   0    0G   200.000  6000.000  1000.0      1
 0.77244792E+01 0.28322316E-01-0.11545236E-04 0.21540815E-08-0.15054178E-12    2
-0.78261573E+03-0.19769698E+02 0.26898140E+01 0.20954550E-02 0.11303687E-03    3
-0.15408070E-06 0.62763658E-10 0.23139663E+04 0.15294056E+02 0.39328836E+04    4
lC5H9             T03/97C   5H   9O   0    0G   200.000  6000.000  1000.0      1
 0.20313000E+02 0.10869880E-01-0.19063805E-05 0.00000000E+00 0.00000000E+00    2
 0.94061603E+04-0.82533815E+02 0.11430827E+01 0.44350789E-01-0.17825470E-04    3
 0.00000000E+00 0.00000000E+00 0.16967656E+05 0.24181940E+02 0.19122233E+05    4
cC5H9             T03/97C   5H   9O   0    0G   200.000  6000.000  1000.0      1
 0.11406802E+02 0.22563988E-01-0.70235595E-05 0.11321968E-08-0.73438204E-13    2
 0.75268769E+04-0.39636280E+02 0.29427128E+00 0.13823374E-01 0.90847653E-04    3
-0.13008694E-06 0.53051811E-10 0.12565712E+05 0.27389773E+02 0.13838458E+05    4
C5H4O             T 8/99C   5H   4O   1    0G   200.000  6000.000              1
 1.00806824E+01 1.61143465E-02-5.83314509E-06 9.46759320E-10-5.68972206E-14    2
 1.94364771E+03-2.94521623E+01 2.64576497E-01 3.34873827E-02 1.67738470E-06    3
-2.96207455E-08 1.54431476E-11 5.11159287E+03 2.35409513E+01 6.64245999E+03    4
C5H4OH            T 8/99C   5H   5O   1    0G   200.000  6000.000              1
 1.33741248E+01 1.51996469E-02-5.45685046E-06 8.80944866E-10-5.27493258E-14    2
 2.20358027E+03-4.59569069E+01-1.28398054E+00 4.90298511E-02-1.35844414E-05    3
-2.92983743E-08 1.90820619E-11 6.37364803E+03 3.08073591E+01 8.00114499E+03    4
C5H5O(2,4)        D 9/97C   5H   5O   1   00G   300.000  3000.000              1
 0.85405312E+01 0.22989510E-01-0.95437563E-05 0.17061612E-08-0.97459360E-13    2
 0.22263699E+05-0.20818825E+02-0.30777600E+01 0.52581679E-01-0.28856513E-04    3
-0.33885479E-08 0.63361399E-11 0.25510455E+05 0.39591522E+02 0.26570048E+05    4
C5H5O(1,3)        DU0997C   5H   5O   1   00G   300.000  3000.000   1000.00    1
 0.92431440E+01 0.22201257E-01-0.93105946E-05 0.17155222E-08-0.10613969E-12    2
 0.15908394E+04-0.24087738E+02-0.29566984E+01 0.55851892E-01-0.37241636E-04    3
 0.41624357E-08 0.39272010E-11 0.48573193E+04 0.38676682E+02                   4
C5H5OH            HWZD99C   5H   6O   1    0G   300.000  3000.000              1
 0.34893970E+01 0.38052600E-01-0.21654527E-04 0.59238574E-08-0.62763461E-12    2
-0.82131025E+04 0.71248055E+01-0.50430169E+01 0.71253479E-01-0.70918177E-04    3
 0.38680220E-07-0.87888264E-11-0.64167788E+04 0.48617100E+02                   4
bi-C5H4O          DU0997C   5H   4O   1   00G   300.000  3000.000   1000.00    1
 0.10514051E+02 0.16667502E-01-0.61001861E-05 0.81804008E-09-0.88743752E-14    2
 0.27501334E+05-0.30678673E+02-0.35879545E+01 0.59943721E-01-0.52969943E-04    3
 0.19971461E-07-0.14667430E-11 0.31091709E+05 0.40873169E+02                   4
lC5H6             HWZD99C   5H   6    0    0G   300.     3000.      1000.      1
 0.86914568E+01 0.21268958E-01-0.79818547E-05 0.11795985E-08-0.35253359E-13    2
 0.25763865E+05-0.19189083E+02 0.58391756E+00 0.42602919E-01-0.24962495E-04    3
 0.25815752E-08 0.23169964E-11 0.28043699E+05 0.22916780E+02                   4
lC5H7             HWZD99C   5H   7    0    0G   300.000  3000.000              1
 0.22246480E+01 0.39601296E-01-0.22345617E-04 0.60649676E-08-0.63840047E-12    2
 0.22303428E+05 0.14009951E+02-0.40974307E+01 0.61832044E-01-0.48770780E-04    3
 0.16696418E-07-0.75334899E-12 0.23683646E+05 0.45148109E+02                   4
C6H2              D11/99C   6H   2    0    0G   300.000  3000.000              1
 0.12893918E+02 0.79145068E-02-0.24027240E-05 0.24340149E-09 0.31383246E-14    2
 0.79832406E+05-0.40771996E+02 0.45099974E+00 0.67475192E-01-0.11809925E-03    3
 0.10367632E-06-0.34851039E-10 0.82173062E+05 0.17704124E+02                   4
C6H               P 1/93C   6H   1    0    0G   300.000  3000.000              1
 0.12370055E+02 0.52177699E-02-0.16885009E-05 0.25807149E-09-0.15472851E-13    2
 0.12158739E+06-0.34952797E+02-0.25630299E+00 0.63793827E-01-0.11440118E-03    3
 0.10136744E-06-0.34361855E-10 0.12408855E+06 0.24930750E+02                   4
l-C6H4            H6W/94C   6H   4    0    0G   300.000  3000.000              1
 0.12715182E+02 0.13839662E-01-0.43765440E-05 0.31541636E-09 0.46619026E-13    2
 0.57031148E+05-0.39464600E+02 0.29590225E+00 0.58053318E-01-0.67766756E-04    3
 0.43376762E-07-0.11418864E-10 0.60001371E+05 0.22318970E+02                   4
l-C6H6            H6W/94C   6H   6    0    0G   300.000  3000.000              1
 0.17584442E+02 0.64486600E-02 0.48933980E-05-0.34696221E-08 0.56150749E-12    2
 0.34111988E+05-0.66017838E+02-0.10170622E+01 0.61794821E-01-0.59461061E-04    3
 0.31873491E-07-0.71717693E-11 0.39202707E+05 0.29460373E+02                   4
c-C6H7            H6W/94C   6H   7    0    0G   300.000  3000.000              1
 0.19996841E+02 0.11189543E-02 0.11649756E-04-0.62779471E-08 0.94939508E-12    2
 0.16730059E+05-0.83746933E+02-0.30328493E+01 0.50804518E-01-0.69150292E-05    3
-0.29715974E-07 0.16296353E-10 0.23895383E+05 0.38909180E+02                   4
n-C6H7            H6W/94C   6H   7    0    0G   300.000  3000.000              1
 0.22577469E+02-0.30737517E-02 0.14225234E-04-0.69880848E-08 0.10232874E-11    2
 0.41228980E+05-0.91568619E+02 0.13248032E+00 0.57103366E-01-0.43712644E-04    3
 0.15538603E-07-0.12976356E-11 0.47730512E+05 0.25339081E+02                   4
C6H8              H6W/94C   6H   8    0    0G   300.000  3000.000              1
 0.28481979E+02-0.15702948E-01 0.26771697E-04-0.11780109E-07 0.16573427E-11    2
 0.93346445E+04-0.12500226E+03 0.15850439E+01 0.40215142E-01 0.78439543E-05    3
-0.38761325E-07 0.18545207E-10 0.17949613E+05 0.19112625E+02                   4
cC6H8             T03/97C   6H   8O   0    0G   200.000  6000.000  1000.0      1
 0.11779870E+02 0.25519980E-01-0.92666947E-05 0.15068122E-08-0.90658701E-13    2
 0.65486686E+04-0.41618805E+02 0.17265319E+01 0.14887612E-01 0.94809230E-04    3
-0.14083394E-06 0.58859873E-10 0.11021297E+05 0.19130886E+02 0.12784878E+05    4
lC6H9             T 2/92C   6H   9O    0   0G   200.000  3000.000 1000.0       1
 0.23165919E+02 0.10813608E-01-0.17638168E-05 0.00000000E+00 0.00000000E+00    2
 0.11162402E+05-0.98600332E+02 0.31671271E+00 0.52069818E-01-0.21965057E-04    3
 0.00000000E+00 0.00000000E+00 0.19926824E+05 0.27879902E+02 0.22141533E+05    4
cC6H9             T 2/92C   6H   9O    0   0G   200.000  3000.000 1000.0       1
 0.26295828E+02 0.86828857E-02-0.15770376E-05 0.00000000E+00 0.00000000E+00    2
 0.20863563E+04-0.12573825E+03-0.35714300E+01 0.61696043E-01-0.26928803E-04    3
 0.00000000E+00 0.00000000E+00 0.13657039E+05 0.39986250E+02 0.15096500E+05    4
cC6H10            T03/97C   6H  10O   0    0G   200.000  6000.000  1000.0      1
 0.11773904E+02 0.30947360E-01-0.11234330E-04 0.18262494E-08-0.10985119E-12    2
-0.72028376E+04-0.42658688E+02 0.23662378E+01 0.10681712E-01 0.11822112E-03    3
-0.16567854E-06 0.67612802E-10-0.24824973E+04 0.16769357E+02-0.55324968E+03    4
C6H3              H6W/94C   6H   3    0    0G   300.000  3000.000              1
 0.58188343E+01 0.27933408E-01-0.17825427E-04 0.53702536E-08-0.61707627E-12    2
 0.85188250E+05-0.92147827E+00 0.11790619E+01 0.55547360E-01-0.73076168E-04    3
 0.52076736E-07-0.15046964E-10 0.85647312E+05 0.19179199E+02                   4
i-C6H5            H6W/94C   6H   5    0    0G   300.000  3000.000              1
 0.22501663E+02-0.81009977E-02 0.15955695E-04-0.72310371E-08 0.10310424E-11    2
 0.58473410E+05-0.91224777E+02-0.78585434E+00 0.60221825E-01-0.62890264E-04    3
 0.36310730E-07-0.87000259E-11 0.64942270E+05 0.28658905E+02                   4
i-C6H7            H6W/94C   6H   7    0    0G   300.000  3000.000              1
 0.20481506E+02 0.79439697E-03 0.11450761E-04-0.60991177E-08 0.91756724E-12    2
 0.37728426E+05-0.81812073E+02-0.17099094E+01 0.62486034E-01-0.54290707E-04    3
 0.26959682E-07-0.58999090E-11 0.44086621E+05 0.33344772E+02                   4
o-C6H4            D11/99C   6H   4    0    0G   300.000  3000.000              1
 0.88432961E+01 0.20301474E-01-0.88674269E-05 0.17264292E-08-0.11786047E-12    2
 0.49317113E+05-0.24014301E+02-0.38454189E+01 0.58391564E-01-0.48644750E-04    3
 0.16770320E-07-0.78580680E-12 0.52592500E+05 0.40587132E+02                   4
m-C6H4            D11/99C   6H   4    0    0G   300.000  3000.000              1
 0.95307283E+01 0.19178549E-01-0.80941481E-05 0.14811132E-08-0.88632260E-13    2
 0.56865535E+05-0.27623203E+02-0.39450364E+01 0.59887171E-01-0.50811577E-04    3
 0.17603140E-07-0.72608743E-12 0.60323117E+05 0.40899506E+02                   4
p-C6H4            D11/99C   6H   4    0    0G   300.000  3000.000              1
 0.98300371E+01 0.18499156E-01-0.75165058E-05 0.12727610E-08-0.61767120E-13    2
 0.64446117E+05-0.29418266E+02-0.39744728E+01 0.58399867E-01-0.44950713E-04    3
 0.10307744E-07 0.22412619E-11 0.68058992E+05 0.41168865E+02                   4
l-C6H4Z           D11/99C   6H   4    0    0G   300.000  3000.000              1
 0.11186811E+02 0.17122138E-01-0.73898623E-05 0.14678845E-08-0.10733922E-12    2
 0.60743207E+05-0.29537384E+02 0.20895090E+01 0.53276263E-01-0.63299172E-04    3
 0.40811642E-07-0.10598600E-10 0.62662203E+05 0.14613283E+02                   4
nC6H5             D11/99C   6H   5    0    0G   300.000  3000.000              1
 0.11263281E+02 0.19379666E-01-0.76874276E-05 0.12866819E-08-0.63244650E-13    2
 0.68052773E+05-0.30487534E+02-0.27013230E+00 0.59389681E-01-0.60963321E-04    3
 0.33169378E-07-0.71466453E-11 0.70785828E+05 0.26953651E+02                   4
C6H5              D11/99C   6H   5    0    0G   300.000  3000.000              1
 0.85973110E+01 0.22241630E-01-0.87199978E-05 0.13788785E-08-0.53146056E-13    2
 0.36261047E+05-0.22954643E+02-0.36931453E+01 0.52178968E-01-0.25558427E-04    3
-0.70661121E-08 0.75833975E-11 0.39779590E+05 0.41332535E+02                   4
C6H6              D11/99C   6H   6    0    0G   300.000  3000.000              1
 0.91381245E+01 0.23854433E-01-0.88127726E-05 0.12099021E-08-0.18221503E-13    2
 0.52043462E+04-0.29115665E+02-0.48437734E+01 0.58427613E-01-0.29485855E-04    3
-0.69390440E-08 0.82125253E-11 0.91817773E+04 0.43889832E+02                   4
C5H5CH3           P 1/93C   6H   8    0    0G   300.000  2500.000              1
 0.14628364E+02 0.19849248E-01-0.50529134E-05 0.10556275E-10 0.11381723E-12    2
 0.55674092E+04-0.56114021E+02-0.45763016E-01 0.29978730E-01 0.61898092E-04    3
-0.11171783E-06 0.49435803E-10 0.10927480E+05 0.26558569E+02                   4
C5H4CH2           P 1/93C   6H   6    0    0G   300.000  2500.000              1
 0.75731055E+04-0.18843678E+02 0.17058320E-01-0.65980571E-05 0.93053393E-09    2
-0.22894220E+07-0.40003195E+05 0.78428810E+02-0.43919629E+00 0.13370259E-02    3
-0.14196110E-05 0.56357985E-09 0.22226365E+05-0.41005380E+03                   4
C6H5C6H5          HW /94C  12H  10    0    0G   300.000  3000.000              1
 0.50761871E+02-0.34501564E-01 0.50293413E-04-0.21559579E-07 0.30097192E-11    2
 0.21538867E+04-0.24670712E+03-0.10283234E+02 0.12428707E+00-0.95990268E-04    3
 0.32294793E-07-0.23045229E-11 0.20165258E+05 0.72707947E+02                   4
C6H5C2H           H6W/94C   8H   6    0    0G   300.000  3000.000              1
 0.24090759E+02 0.78232400E-03 0.11453964E-04-0.61620504E-08 0.93346685E-12    2
 0.27429445E+05-0.10499631E+03-0.52645016E+01 0.84511042E-01-0.76597848E-04    3
 0.33216978E-07-0.47673063E-11 0.35566242E+05 0.46378815E+02                   4
C6H5CH3           L 6/87C   7H   8    0    0G   200.000  6000.000              1
 0.12940034E+02 0.26691287E-01-0.96838505E-05 0.15738629E-08-0.94663601E-13    2
-0.69764908E+03-0.46728785E+02 0.16152663E+01 0.21099438E-01 0.85366018E-04    3
-0.13261066E-06 0.55956604E-10 0.40756300E+04 0.20282210E+02 0.60135835E+04    4
C6H5CH2           T08/90C   7H   7    0    0G   200.000  6000.000              1
 0.14043980E+02 0.23493873E-01-0.85375367E-05 0.13890841E-08-0.83614420E-13    2
 0.18564203E+05-0.51665589E+02 0.48111540E+00 0.38512832E-01 0.32861492E-04    3
-0.76972721E-07 0.35423068E-10 0.23307027E+05 0.23548820E+02 0.25317186E+05    4
C6H5C2H3          T12/94C   8H   8    0    0G   298.150  5000.000              1
 0.16139277E+02 0.24210847E-01-0.72678359E-05 0.11392276E-08-0.72984881E-13    2
 0.10249251E+05-0.61169437E+02-0.10717708E+02 0.12666725E+00-0.17762493E-03    3
 0.14344049E-06-0.47616577E-10 0.16597133E+05 0.71526331E+02 0.17723291E+05    4
C6H5CH2OH         L 7/87C   7H   8O   1    0G   200.000  6000.000 1000.00      1
 0.15281154E+02 0.27208501E-01-0.98584660E-05 0.16012183E-08-0.96278057E-13    2
-0.19700471E+05-0.59418673E+02 0.20642021E+01 0.22775140E-01 0.95972053E-04    3
-0.15085110E-06 0.64175832E-10-0.14285021E+05 0.18148312E+02-0.12077200E+05    4
C6H5CHO           L 3/86C   7H   6O   1    0G   298.150  5000.000 1000.00      1
 0.13650737E+02 0.25680419E-01-0.10466729E-04 0.19413430E-08-0.13483792E-12    2
-0.11019744E+05-0.47965796E+02-0.31627334E+01 0.66369245E-01-0.34816353E-04    3
-0.62999377E-08 0.85807101E-11-0.61169349E+04 0.40231735E+02-0.44259974E+04    4
C6H5CO    EST/BUR P 1/93C   7H   5O   1    0G   300.000  2500.000              1
 0.13374409E+02 0.23999289E-01-0.10465724E-04 0.21669131E-08-0.18007045E-12    2
 0.69147837E+04-0.44659218E+02-0.20251155E+01 0.61512541E-01-0.31603653E-04    3
-0.69724599E-08 0.79835149E-11 0.11255803E+05 0.35778175E+02                   4
C6H4O2            PUML96C   6H   4O   2    0G   300.000  5000.000 1000.000     1
 0.11730840E+02 0.23614995E-01-0.10234576E-04 0.19532174E-08-0.12746022E-12    2
-0.21085770E+05-0.36300453E+02-0.95193005E+00 0.57842445E-01-0.38214439E-04    3
 0.46312656E-08 0.36296651E-11-0.17611047E+05 0.29239513E+02                   4
C6H5O             T05/02C   6H   5O   1    0G   200.000  6000.000 1000.000     1
 1.37221720E+01 1.74688771E-02-6.35504520E-06 1.03492308E-09-6.23410504E-14    2
 2.87274751E+02-4.88181680E+01-4.66204455E-01 4.13443975E-02 1.32412991E-05    3
-5.72872769E-08 2.89763707E-11 4.77858391E+03 2.76990274E+01 6.49467016E+03    4
C6H5OH            L 4/84C   6H   6O   1    0G   300.000  5000.000              1
 0.14912073E+02 0.18378135E-01-0.61983128E-05 0.91983221E-09-0.49209565E-13    2
-0.18375199E+05-0.55924103E+02-0.16956539E+01 0.52271299E-01-0.72024050E-05    3
-0.35859603E-07 0.20449073E-10-0.13284121E+05 0.32542160E+02-0.11594207E+05    4
C6H5C2H5          A 6/83C   8H  10    0    0G   300.     3000.    1000.00      1
  0.3878978E+01  0.5810059E-01 -0.3196380E-04  0.8448993E-08 -0.8694825E-12    2
 -0.5024922E+03  0.3837099E+01 -0.7266845E+01  0.1003089E+00 -0.9651715E-04    3
  0.5565908E-07 -0.1453370E-10  0.1987290E+04  0.5857746E+02  0.3529492E+04    4
HOC6H4CH3 AVG CRESOL6/87C   7H   8O   1    0G   200.000  6000.000 1000.00      1
 0.15932987E+02 0.27011160E-01-0.99448722E-05 0.16296689E-08-0.98513298E-13    2
-0.23592065E+05-0.59732841E+02 0.42258267E+00 0.45551636E-01 0.32012513E-04    3
-0.81121959E-07 0.37665658E-10-0.18202621E+05 0.26032903E+02-0.15911701E+05    4
OC6H4CH3  EST/BUR P 1/93C   7H   7O   1    0G   300.000  2500.000              1
 0.22609371E+02 0.75646150E-02 0.65960894E-05-0.47150865E-08 0.80409063E-12    2
-0.82025244E+04-0.97292511E+02-0.28855777E+00 0.48003536E-01 0.18032993E-04    3
-0.61741488E-07 0.28852587E-10-0.68945581E+03 0.26720068E+02                   4
bi-C6H5CH2        A 6/83C  14H  14    0    0G   300.     3000.    1000.00      1
  0.7292035E+01  0.9250200E-01 -0.5168641E-04  0.1362709E-07 -0.1381148E-11    2
  0.1031673E+05 -0.1132738E+02 -0.1388958E+02  0.1720984E+00 -0.1700660E-03    3
  0.9601888E-07 -0.2373253E-10  0.1503234E+05  0.9270736E+02  0.1721641E+05    4
C10H8             H6W/94C  10H   8    0    0G   300.000  3000.000              1
 0.36468643E+02-0.15419513E-01 0.30160038E-04-0.13700120E-07 0.19582730E-11    2
 0.35091445E+04-0.17329489E+03-0.94505043E+01 0.11137849E+00-0.10345667E-03    3
 0.52800392E-07-0.11804439E-10 0.16695594E+05 0.65187668E+02                   4
C6H4CH3           P 1/93C   7H   7    0    0G   300.000  2500.000              1
 0.11615498E+02 0.27431838E-01-0.10899345E-04 0.18641830E-08-0.10191607E-12    2
 0.31209334E+05-0.38994637E+02-0.31415942E+01 0.56723077E-01-0.86885111E-05    3
-0.34249616E-07 0.19266902E-10 0.35738547E+05 0.39742840E+02                   4
NC5H12     1/ 2/ 7 THERMC   5H  12    0    0G   300.000  5000.000 1390.000    41
 1.57257603E+01 2.61086045E-02-8.90970996E-06 1.38102248E-09-8.00296536E-14    2
-2.60519543E+04-6.03365457E+01-7.36766122E-01 6.07200973E-02-3.57592761E-05    3
 1.04907042E-08-1.21487315E-12-1.98934934E+04 2.95358287E+01                   4
PXC5H11    1/ 2/ 7 THERMC   5H  11    0    0G   300.000  5000.000 1390.000    41
 1.52977446E+01 2.39735310E-02-8.18392948E-06 1.26883076E-09-7.35409055E-14    2
-9.80712307E+02-5.44829293E+01 5.24384081E-02 5.60796958E-02-3.31545803E-05    3
 9.77533781E-09-1.14009660E-12 4.71611460E+03 2.87238666E+01                   4
SXC5H11    1/ 2/ 7 THERMC   5H  11    0    0G   300.000  5000.000 1377.000    41
 1.50998007E+01 2.37225333E-02-8.01388900E-06 1.23431039E-09-7.12300125E-14    2
-2.33420039E+03-5.29613979E+01 4.98943592E-01 5.09850184E-02-2.40687488E-05    3
 3.59465211E-09 3.01383099E-13 3.40702366E+03 2.78600953E+01                   4
S2XC5H11   1/ 2/ 7 THERMC   5H  11    0    0G   300.000  5000.000 1377.000    41
 1.50998007E+01 2.37225333E-02-8.01388900E-06 1.23431039E-09-7.12300125E-14    2
-2.33420039E+03-5.29613979E+01 4.98943592E-01 5.09850184E-02-2.40687488E-05    3
 3.59465211E-09 3.01383099E-13 3.40702366E+03 2.78600953E+01                   4
C5H10      1/ 2/ 7 THERMC   5H  10    0    0G   300.000  5000.000 1392.000    31
 1.45851539E+01 2.24072471E-02-7.63348025E-06 1.18188966E-09-6.84385139E-14    2
-1.00898205E+04-5.23683936E+01-1.06223481E+00 5.74218294E-02-3.74486890E-05    3
 1.27364989E-08-1.79609789E-12-4.46546666E+03 3.22739790E+01                   4
C5H9       1/ 2/ 7 THERMC   5H   9    0    0G   300.000  5000.000 1386.000    31
 1.45546753E+01 2.04331611E-02-7.07391431E-06 1.10720105E-09-6.46021493E-14    2
 6.76813422E+03-5.37243740E+01-1.22720644E+00 5.49419289E-02-3.56284020E-05    3
 1.18245127E-08-1.61716004E-12 1.25391664E+04 3.19626373E+01                   4
NC6H14     1/ 2/ 7 THERMC   6H  14    0    0G   300.000  5000.000 1390.000    51
 1.89634117E+01 3.04480204E-02-1.03794829E-05 1.60775457E-09-9.31269728E-14    2
-3.01628739E+04-7.62839355E+01-9.69606184E-01 7.29085608E-02-4.38853919E-05    3
 1.32312807E-08-1.58437423E-12-2.27803862E+04 3.23069798E+01                   4
PXC6H13    1/ 2/ 7 THERMC   6H  13    0    0G   300.000  5000.000 1390.000    51
 1.85385470E+01 2.83107962E-02-9.65307246E-06 1.49547585E-09-8.66336064E-14    2
-5.09299041E+03-7.04490943E+01-2.04871465E-01 6.83801272E-02-4.14447912E-05    3
 1.26155802E-08-1.53120058E-12 1.83280393E+03 3.16075093E+01                   4
SXC6H13    1/ 2/ 7 THERMC   6H  13    0    0G   300.000  5000.000 1380.000    51
 1.83687363E+01 2.80268110E-02-9.47032396E-06 1.45888527E-09-8.42002461E-14    2
-6.46093974E+03-6.90934018E+01 2.29560149E-01 6.33327323E-02-3.24135431E-05    3
 6.46387687E-09-9.61420427E-14 5.25639156E+02 3.08006138E+01                   4
S2XC6H13   1/ 2/ 7 THERMC   6H  13    0    0G   300.000  5000.000 1380.000    51
 1.83687363E+01 2.80268110E-02-9.47032396E-06 1.45888527E-09-8.42002461E-14    2
-6.46093974E+03-6.90934018E+01 2.29560149E-01 6.33327323E-02-3.24135431E-05    3
 6.46387687E-09-9.61420427E-14 5.25639156E+02 3.08006138E+01                   4
C6H12      1/22/ 7 THERMC   6H  12    0    0G   300.000  5000.000 1392.000    41
 1.78337529E+01 2.67377658E-02-9.10036773E-06 1.40819768E-09-8.15124244E-14    2
-1.42062860E+04-6.83818851E+01-1.35275205E+00 6.98655426E-02-4.59408022E-05    3
 1.56967343E-08-2.21296175E-12-7.34368617E+03 3.53120691E+01                   4
C6H11      1/22/ 7 THERMC   6H  11    0    0G   300.000  5000.000 1389.000    41
 1.77336550E+01 2.48934775E-02-8.59991450E-06 1.34412828E-09-7.83475666E-14    2
 2.68017174E+03-6.93471508E+01-1.55544944E+00 6.76865602E-02-4.47048635E-05    3
 1.52236630E-08-2.14346377E-12 9.66316695E+03 3.51482658E+01                   4
NC7H16     1/ 2/ 7 THERMC   7H  16    0    0G   300.000  5000.000 1391.000    61
 2.22148969E+01 3.47675750E-02-1.18407129E-05 1.83298478E-09-1.06130266E-13    2
-3.42760081E+04-9.23040196E+01-1.26836187E+00 8.54355820E-02-5.25346786E-05    3
 1.62945721E-08-2.02394925E-12-2.56586565E+04 3.53732912E+01                   4
PXC7H15    1/ 2/ 7 THERMC   7H  15    0    0G   300.000  5000.000 1390.000    61
 2.17940709E+01 3.26280243E-02-1.11138244E-05 1.72067148E-09-9.96366999E-14    2
-9.20938221E+03-8.64954311E+01-4.99570406E-01 8.08826467E-02-5.00532754E-05    3
 1.56549308E-08-1.96616227E-12-1.04590223E+03 3.46564011E+01                   4
SXC7H15    1/ 2/ 7 THERMC   7H  15    0    0G   300.000  5000.000 1382.000    61
 2.16368842E+01 3.23324804E-02-1.09273807E-05 1.68357060E-09-9.71774091E-14    2
-1.05873616E+04-8.52209653E+01-3.79155767E-02 7.56726570E-02-4.07473634E-05    3
 9.32678943E-09-4.92360745E-13-2.35605303E+03 3.37321506E+01                   4
S2XC7H15   1/ 2/ 7 THERMC   7H  15    0    0G   300.000  5000.000 1382.000    61
 2.16368842E+01 3.23324804E-02-1.09273807E-05 1.68357060E-09-9.71774091E-14    2
-1.05873616E+04-8.52209653E+01-3.79155767E-02 7.56726570E-02-4.07473634E-05    3
 9.32678943E-09-4.92360745E-13-2.35605303E+03 3.37321506E+01                   4
S3XC7H15   1/ 2/ 7 THERMC   7H  15    0    0G   300.000  5000.000 1382.000    61
 2.16368842E+01 3.23324804E-02-1.09273807E-05 1.68357060E-09-9.71774091E-14    2
-1.05873616E+04-8.52209653E+01-3.79155767E-02 7.56726570E-02-4.07473634E-05    3
 9.32678943E-09-4.92360745E-13-2.35605303E+03 3.37321506E+01                   4
C7H14      1/ 2/ 7 THERMC   7H  14    0    0G   300.000  5000.000 1392.000    51
 2.10898039E+01 3.10607878E-02-1.05644793E-05 1.63405780E-09-9.45598219E-14    2
-1.83260065E+04-8.44391108E+01-1.67720549E+00 8.24611601E-02-5.46504108E-05    3
 1.87862303E-08-2.65737983E-12-1.02168601E+04 3.85068032E+01                   4
C7H13      1/ 2/ 7 THERMC   7H  13    0    0G   300.000  5000.000 1389.000    51
 2.09278134E+01 2.92841022E-02-1.00899640E-05 1.57425554E-09-9.16505991E-14    2
-1.41296217E+03-8.50447021E+01-1.66945935E+00 7.93202117E-02-5.20231516E-05    3
 1.74804211E-08-2.40912996E-12 6.75927466E+03 3.73759517E+01                   4
NC8H18     1/ 2/ 7 THERMC   8H  18    0    0G   300.000  5000.000 1391.000    71
 2.54710194E+01 3.90887037E-02-1.33038777E-05 2.05867527E-09-1.19167174E-13    2
-3.83962755E+04-1.08361094E+02-1.54218406E+00 9.78112063E-02-6.09318358E-05    3
 1.92005591E-08-2.42996250E-12-2.85395641E+04 3.83327978E+01                   4
PXC8H17    1/ 2/ 7 THERMC   8H  17    0    0G   300.000  5000.000 1390.000    71
 2.50510356E+01 3.69480162E-02-1.25765264E-05 1.94628409E-09-1.12668898E-13    2
-1.33300535E+04-1.02557384E+02-7.72759438E-01 9.32549705E-02-5.84447245E-05    3
 1.85570214E-08-2.37127483E-12-3.92689511E+03 3.76130631E+01                   4
SXC8H17    1/ 2/ 7 THERMC   8H  17    0    0G   300.000  5000.000 1383.000    71
 2.49043504E+01 3.66393792E-02-1.23849888E-05 1.90835394E-09-1.10160725E-13    2
-1.47135299E+04-1.01344740E+02-3.04233324E-01 8.80077253E-02-4.90742611E-05    3
 1.21857563E-08-8.87773198E-13-5.23792835E+03 3.66582632E+01                   4
S2XC8H17   1/ 2/ 7 THERMC   8H  17    0    0G   300.000  5000.000 1383.000    71
 2.49043504E+01 3.66393792E-02-1.23849888E-05 1.90835394E-09-1.10160725E-13    2
-1.47135299E+04-1.01344740E+02-3.04233324E-01 8.80077253E-02-4.90742611E-05    3
 1.21857563E-08-8.87773198E-13-5.23792835E+03 3.66582632E+01                   4
S3XC8H17   1/ 2/ 7 THERMC   8H  17    0    0G   300.000  5000.000 1383.000    71
 2.49043504E+01 3.66393792E-02-1.23849888E-05 1.90835394E-09-1.10160725E-13    2
-1.47135299E+04-1.01344740E+02-3.04233324E-01 8.80077253E-02-4.90742611E-05    3
 1.21857563E-08-8.87773198E-13-5.23792835E+03 3.66582632E+01                   4
C8H16      1/ 2/ 7 THERMC   8H  16    0    0G   300.000  5000.000 1392.000    61
 2.43540125E+01 3.53666462E-02-1.20208388E-05 1.85855053E-09-1.07522262E-13    2
-2.24485674E+04-1.00537716E+02-1.89226915E+00 9.46066357E-02-6.27385521E-05    3
 2.15158309E-08-3.02718683E-12-1.31074559E+04 4.11878981E+01                   4
C8H15      1/ 2/ 7 THERMC   8H  15    0    0G   300.000  5000.000 1390.000    61
 2.41485380E+01 3.36466429E-02-1.15692934E-05 1.80261649E-09-1.04847473E-13    2
-5.51631151E+03-1.00894875E+02-1.90098561E+00 9.15067740E-02-6.01588113E-05    3
 2.02337556E-08-2.78235289E-12 3.87213558E+03 4.01405031E+01                   4
NC9H20     1/ 2/ 7 THERMC   9H  20    0    0G   300.000  5000.000 1391.000    81
 2.87289600E+01 4.34074576E-02-1.47660985E-05 2.28420987E-09-1.32194796E-13    2
-4.25174479E+04-1.24428751E+02-1.81390458E+00 1.10176644E-01-6.93124463E-05    3
 2.20957601E-08-2.83355715E-12-3.14207716E+04 4.12827220E+01                   4
PXC9H19    1/ 2/ 7 THERMC   9H  19    0    0G   300.000  5000.000 1390.000    81
 2.83097514E+01 4.12657344E-02-1.40383289E-05 2.17174871E-09-1.25692307E-13    2
-1.74516030E+04-1.16837897E+02-1.04387292E+00 1.05617283E-01-6.68199971E-05    3
 2.14486166E-08-2.77404275E-12-6.80818512E+03 4.23518992E+01                   4
SXC9H19    1/ 2/ 7 THERMC   9H  19    0    0G   300.000  5000.000 1386.000    81
 2.80393256E+01 4.11440297E-02-1.39260043E-05 2.14745952E-09-1.24022172E-13    2
-1.87727855E+04-1.16696832E+02-5.76046059E-01 1.00511821E-01-5.77755119E-05    3
 1.53276702E-08-1.35056750E-12-8.12289806E+03 3.95795606E+01                   4
S2XC9H19   1/ 2/ 7 THERMC   9H  19    0    0G   300.000  5000.000 1386.000    81
 2.80393256E+01 4.11440297E-02-1.39260043E-05 2.14745952E-09-1.24022172E-13    2
-1.87727855E+04-1.16696832E+02-5.76046059E-01 1.00511821E-01-5.77755119E-05    3
 1.53276702E-08-1.35056750E-12-8.12289806E+03 3.95795606E+01                   4
S3XC9H19   1/ 2/ 7 THERMC   9H  19    0    0G   300.000  5000.000 1386.000    81
 2.80393256E+01 4.11440297E-02-1.39260043E-05 2.14745952E-09-1.24022172E-13    2
-1.87727855E+04-1.16696832E+02-5.76046059E-01 1.00511821E-01-5.77755119E-05    3
 1.53276702E-08-1.35056750E-12-8.12289806E+03 3.95795606E+01                   4
S4XC9H19   1/ 2/ 7 THERMC   9H  19    0    0G   300.000  5000.000 1386.000    81
 2.80393256E+01 4.11440297E-02-1.39260043E-05 2.14745952E-09-1.24022172E-13    2
-1.87727855E+04-1.16696832E+02-5.76046059E-01 1.00511821E-01-5.77755119E-05    3
 1.53276702E-08-1.35056750E-12-8.12289806E+03 3.95795606E+01                   4
C9H18      1/ 2/ 7 THERMC   9H  18    0    0G   300.000  5000.000 1392.000    71
 2.76142176E+01 3.96825287E-02-1.34819446E-05 2.08390452E-09-1.20539294E-13    2
-2.65709061E+04-1.16618623E+02-2.16108263E+00 1.06958297E-01-7.10973244E-05    3
 2.43971077E-08-3.42771547E-12-1.59890847E+04 4.41245128E+01                   4
C9H17      1/ 3/ 7 THERMC   9H  17    0    0G   300.000  5000.000 1390.000    71
 2.73846125E+01 3.79931250E-02-1.30425564E-05 2.02998652E-09-1.17985203E-13    2
-9.62653063E+03-1.17685762E+02-2.20061475E+00 1.03997426E-01-6.87316813E-05    3
 2.32482420E-08-3.21153346E-12 9.95158605E+02 4.23691203E+01                   4
NC10H22    1/ 2/ 7 THERMC  10H  22    0    0G   300.000  5000.000 1391.000    91
 3.19882239E+01 4.77244922E-02-1.62276391E-05 2.50963259E-09-1.45215772E-13    2
-4.66392840E+04-1.40504121E+02-2.08416969E+00 1.22535012E-01-7.76815739E-05    3
 2.49834877E-08-3.23548038E-12-3.43021863E+04 4.42260140E+01                   4
PXC10H21   1/ 2/ 7 THERMC  10H  21    0    0G   300.000  5000.000 1390.000    91
 3.15697160E+01 4.55818403E-02-1.54994965E-05 2.39710933E-09-1.38709559E-13    2
-2.15737832E+04-1.34708986E+02-1.31358348E+00 1.17972813E-01-7.51843079E-05    3
 2.43331106E-08-3.17522852E-12-9.68967550E+03 4.35010452E+01                   4
SXC10H21   1/ 2/ 7 THERMC  10H  21    0    0G   300.000  5000.000 1385.000    91
 3.14447580E+01 4.52778532E-02-1.53145696E-05 2.36072411E-09-1.36311835E-13    2
-2.29702700E+04-1.33634423E+02-9.30536886E-01 1.13137924E-01-6.64034118E-05    3
 1.83220872E-08-1.77128003E-12-1.09890165E+04 4.29335080E+01                   4
S2XC10H21  1/ 2/ 7 THERMC  10H  21    0    0G   300.000  5000.000 1385.000    91
 3.14447580E+01 4.52778532E-02-1.53145696E-05 2.36072411E-09-1.36311835E-13    2
-2.29702700E+04-1.33634423E+02-9.30536886E-01 1.13137924E-01-6.64034118E-05    3
 1.83220872E-08-1.77128003E-12-1.09890165E+04 4.29335080E+01                   4
S3XC10H21  1/ 2/ 7 THERMC  10H  21    0    0G   300.000  5000.000 1385.000    91
 3.14447580E+01 4.52778532E-02-1.53145696E-05 2.36072411E-09-1.36311835E-13    2
-2.29702700E+04-1.33634423E+02-9.30536886E-01 1.13137924E-01-6.64034118E-05    3
 1.83220872E-08-1.77128003E-12-1.09890165E+04 4.29335080E+01                   4
S4XC10H21  1/ 2/ 7 THERMC  10H  21    0    0G   300.000  5000.000 1385.000    91
 3.14447580E+01 4.52778532E-02-1.53145696E-05 2.36072411E-09-1.36311835E-13    2
-2.29702700E+04-1.33634423E+02-9.30536886E-01 1.13137924E-01-6.64034118E-05    3
 1.83220872E-08-1.77128003E-12-1.09890165E+04 4.29335080E+01                   4
C10H20     1/22/ 7 THERMC  10H  20    0    0G   300.000  5000.000 1392.000    81
 3.08753903E+01 4.39971526E-02-1.49425530E-05 2.30917678E-09-1.33551477E-13    2
-3.06937307E+04-1.32705172E+02-2.42901688E+00 1.19305598E-01-7.94489025E-05    3
 2.72736596E-08-3.82718373E-12-1.88708365E+04 4.70571383E+01                   4
C10H19     1/22/ 7 THERMC  10H  19    0    0G   300.000  5000.000 1390.000    81
 3.06442992E+01 4.23089356E-02-1.45032261E-05 2.25520722E-09-1.30991162E-13    2
-1.37453412E+04-1.32907157E+02-2.49340025E+00 1.16496175E-01-7.73332421E-05    3
 2.62845301E-08-3.64630937E-12-1.88400568E+03 4.62584816E+01                   4
NC11H24    1/22/ 7 THERMC  11H  24    0    0G   300.000  5000.000 1391.000    11
 3.52484813E+01 5.20402416E-02-1.76886732E-05 2.73497226E-09-1.58231832E-13    2
-5.07616214E+04-1.56585288E+02-2.35338447E+00 1.34888270E-01-8.60424000E-05    3
 2.78658195E-08-3.63619953E-12-3.71837502E+04 4.71645217E+01                   4
PXC11H23   1/22/ 7 THERMC  11H  23    0    0G   300.000  5000.000 1391.000    11
 3.48306023E+01 4.98967610E-02-1.69601993E-05 2.62239406E-09-1.51722334E-13    2
-2.56964317E+04-1.50793815E+02-1.58230042E+00 1.30323530E-01-8.35408417E-05    3
 2.72125730E-08-3.57529580E-12-1.25713075E+04 4.64373071E+01                   4
SXC11H23   1/22/ 7 THERMC  11H  23    0    0G   300.000  5000.000 1385.000    11
 3.47027943E+01 4.95633551E-02-1.67588574E-05 2.58285265E-09-1.49118584E-13    2
-2.70891758E+04-1.49690851E+02-1.10250355E+00 1.25021794E-01-7.40802024E-05    3
 2.07818900E-08-2.07855904E-12-1.38839647E+04 4.54310281E+01                   4
S2XC11H23  1/22/ 7 THERMC  11H  23    0    0G   300.000  5000.000 1385.000    11
 3.47027943E+01 4.95633551E-02-1.67588574E-05 2.58285265E-09-1.49118584E-13    2
-2.70891758E+04-1.49690851E+02-1.10250355E+00 1.25021794E-01-7.40802024E-05    3
 2.07818900E-08-2.07855904E-12-1.38839647E+04 4.54310281E+01                   4
S3XC11H23  1/22/ 7 THERMC  11H  23    0    0G   300.000  5000.000 1385.000    11
 3.47027943E+01 4.95633551E-02-1.67588574E-05 2.58285265E-09-1.49118584E-13    2
-2.70891758E+04-1.49690851E+02-1.10250355E+00 1.25021794E-01-7.40802024E-05    3
 2.07818900E-08-2.07855904E-12-1.38839647E+04 4.54310281E+01                   4
S4XC11H23  1/22/ 7 THERMC  11H  23    0    0G   300.000  5000.000 1385.000    11
 3.47027943E+01 4.95633551E-02-1.67588574E-05 2.58285265E-09-1.49118584E-13    2
-2.70891758E+04-1.49690851E+02-1.10250355E+00 1.25021794E-01-7.40802024E-05    3
 2.07818900E-08-2.07855904E-12-1.38839647E+04 4.54310281E+01                   4
S5XC11H23  1/22/ 7 THERMC  11H  23    0    0G   300.000  5000.000 1385.000    11
 3.47027943E+01 4.95633551E-02-1.67588574E-05 2.58285265E-09-1.49118584E-13    2
-2.70891758E+04-1.49690851E+02-1.10250355E+00 1.25021794E-01-7.40802024E-05    3
 2.07818900E-08-2.07855904E-12-1.38839647E+04 4.54310281E+01                   4
C11H22     1/ 2/ 7 THERMC  11H  22    0    0G   300.000  5000.000 1392.000    91
 3.41376800E+01 4.83102359E-02-1.64025319E-05 2.53434308E-09-1.46557255E-13    2
-3.48170932E+04-1.48798197E+02-2.69653994E+00 1.31650370E-01-8.77957172E-05    3
 3.01467965E-08-4.22584486E-12-2.17526343E+04 4.99879917E+01                   4
C11H21     1/ 2/ 7 THERMC  11H  21    0    0G   300.000  5000.000 1390.000    91
 3.38968067E+01 4.66347047E-02-1.59682128E-05 2.48119562E-09-1.44045675E-13    2
-1.78637749E+04-1.48943146E+02-2.77277326E+00 1.28897937E-01-8.57717335E-05    3
 2.92169310E-08-4.05817299E-12-4.76410432E+03 4.92434369E+01                   4
NC12H26    1/ 2/ 7 THERMC  12H  26    0    0G   300.000  5000.000 1391.000    11
 3.85095037E+01 5.63550048E-02-1.91493200E-05 2.96024862E-09-1.71244150E-13    2
-5.48843465E+04-1.72670922E+02-2.62181594E+00 1.47237711E-01-9.43970271E-05    3
 3.07441268E-08-4.03602230E-12-4.00654253E+04 5.00994626E+01                   4
PXC12H25   1/ 2/ 7 THERMC  12H  25    0    0G   300.000  5000.000 1390.000    11
 3.80921885E+01 5.42107848E-02-1.84205517E-05 2.84762173E-09-1.64731748E-13    2
-2.98194375E+04-1.66882734E+02-1.85028741E+00 1.42670708E-01-9.18916555E-05    3
 3.00883392E-08-3.97454300E-12-1.54530435E+04 4.93702421E+01                   4
SXC12H25   1/ 2/ 7 THERMC  12H  25    0    0G   300.000  5000.000 1385.000    11
 3.79688268E+01 5.38719464E-02-1.82171263E-05 2.80774503E-09-1.62108420E-13    2
-3.12144988E+04-1.65805933E+02-1.36787089E+00 1.37355348E-01-8.24076158E-05    3
 2.36421562E-08-2.47435932E-12-1.67660539E+04 4.83521895E+01                   4
S2XC12H25  1/ 2/ 7 THERMC  12H  25    0    0G   300.000  5000.000 1385.000    11
 3.79688268E+01 5.38719464E-02-1.82171263E-05 2.80774503E-09-1.62108420E-13    2
-3.12144988E+04-1.65805933E+02-1.36787089E+00 1.37355348E-01-8.24076158E-05    3
 2.36421562E-08-2.47435932E-12-1.67660539E+04 4.83521895E+01                   4
S3XC12H25  1/ 2/ 7 THERMC  12H  25    0    0G   300.000  5000.000 1385.000    11
 3.79688268E+01 5.38719464E-02-1.82171263E-05 2.80774503E-09-1.62108420E-13    2
-3.12144988E+04-1.65805933E+02-1.36787089E+00 1.37355348E-01-8.24076158E-05    3
 2.36421562E-08-2.47435932E-12-1.67660539E+04 4.83521895E+01                   4
S4XC12H25  1/ 2/ 7 THERMC  12H  25    0    0G   300.000  5000.000 1385.000    11
 3.79688268E+01 5.38719464E-02-1.82171263E-05 2.80774503E-09-1.62108420E-13    2
-3.12144988E+04-1.65805933E+02-1.36787089E+00 1.37355348E-01-8.24076158E-05    3
 2.36421562E-08-2.47435932E-12-1.67660539E+04 4.83521895E+01                   4
S5XC12H25  1/ 2/ 7 THERMC  12H  25    0    0G   300.000  5000.000 1385.000    11
 3.79688268E+01 5.38719464E-02-1.82171263E-05 2.80774503E-09-1.62108420E-13    2
-3.12144988E+04-1.65805933E+02-1.36787089E+00 1.37355348E-01-8.24076158E-05    3
 2.36421562E-08-2.47435932E-12-1.67660539E+04 4.83521895E+01                   4
C12H24     1/22/ 7 THERMC  12H  24    0    0G   300.000  5000.000 1391.000    11
 3.74002111E+01 5.26230753E-02-1.78624319E-05 2.75949863E-09-1.59562499E-13    2
-3.89405962E+04-1.64892663E+02-2.96342681E+00 1.43992360E-01-9.61384015E-05    3
 3.30174473E-08-4.62398190E-12-2.46345299E+04 5.29158870E+01                   4
C12H23     1/22/ 7 THERMC  12H  23    0    0G   300.000  5000.000 1391.000    11
 3.71516071E+01 5.09574750E-02-1.74320046E-05 2.70698624E-09-1.57088382E-13    2
-2.19833461E+04-1.64992468E+02-3.04897817E+00 1.41284540E-01-9.41858477E-05    3
 3.21335288E-08-4.46649750E-12-7.64465960E+03 5.22139103E+01                   4
C12H26O2   1/24/ 7 THERMC  12H  26O   2    0G   300.000  5000.000 1390.000    31
 4.35257442E+01 5.78016040E-02-1.99632594E-05 3.11946519E-09-1.81796787E-13    2
-6.68122451E+04-1.92119235E+02-5.33975911E-01 1.56648083E-01-1.05028971E-04    3
 3.67984049E-08-5.38662910E-12-5.09380146E+04 4.62325087E+01                   4
PC12H25O2  1/24/ 7 THERMC  12H  25O   2    0G   300.000  5000.000 1389.000    21
 4.10720559E+01 5.76223879E-02-1.98682057E-05 3.10098353E-09-1.80567383E-13    2
-4.84961160E+04-1.77176648E+02 4.13301762E-01 1.47893311E-01-9.69755376E-05    3
 3.35724323E-08-4.91454646E-12-3.36844450E+04 4.32258761E+01                   4
P12OOHX2   1/24/ 7 THERMC  12H  25O   2    0G   300.000  5000.000 1387.000    31
 4.29994696E+01 5.50833554E-02-1.88955870E-05 2.93978363E-09-1.70823564E-13    2
-4.30658136E+04-1.85856878E+02 1.60721186E-01 1.49325879E-01-9.67266789E-05    3
 3.17154092E-08-4.20449417E-12-2.75569298E+04 4.63656051E+01                   4
C12H26O4   1/24/ 7 THERMC  12H  26O   4    0G   300.000  5000.000 1390.000    51
 5.03368769E+01 5.67494383E-02-1.97115026E-05 3.09203723E-09-1.80689139E-13    2
-8.14466175E+04-2.23085533E+02-2.94753190E-01 1.75805891E-01-1.28643781E-04    3
 4.93995330E-08-7.88291511E-12-6.37896715E+04 4.88174790E+01                   4
SOO12OOH   1/24/ 7 THERMC  12H  25O   4    0G   300.000  5000.000 1389.000    41
 4.78237378E+01 5.66012764E-02-1.96226767E-05 3.07404604E-09-1.79468346E-13    2
-6.30929400E+04-2.07768890E+02 6.91692611E-01 1.67000994E-01-1.20729787E-04    3
 4.63412670E-08-7.45842745E-12-4.65446291E+04 4.56074372E+01                   4
PX212OOH   1/24/ 7 THERMC  12H  25O   4    0G   300.000  5000.000 1389.000    51
 4.89358966E+01 5.52991412E-02-1.91904385E-05 3.00850701E-09-1.75735606E-13    2
-5.73849295E+04-2.11902966E+02 1.58567624E+00 1.63721110E-01-1.14803987E-04    3
 4.18812947E-08-6.34154032E-12-4.05947115E+04 4.33893943E+01                   4
OC12OOH    1/24/ 7 THERMC  12H  24O   3    0G   300.000  5000.000 1392.000    21
 4.49664435E+01 5.45029400E-02-1.88364401E-05 2.94475639E-09-1.71672676E-13    2
-7.80784935E+04-1.96879294E+02-9.04780358E-01 1.62600770E-01-1.18036891E-04    3
 4.52507973E-08-7.22772167E-12-6.21001449E+04 4.93825152E+01                   4
PX212POO   1/24/ 7 THERMC  12H  23O   2    0G   300.000  5000.000 1394.000    01
 4.08114284E+01 5.32062949E-02-1.82595844E-05 2.84109969E-09-1.65082504E-13    2
-5.82425323E+04-1.76923634E+02-1.18205451E+00 1.51086901E-01-1.06393480E-04    3
 3.94388594E-08-6.07746274E-12-4.35503389E+04 4.88285955E+01                   4
ENDOFDATA


#endif
