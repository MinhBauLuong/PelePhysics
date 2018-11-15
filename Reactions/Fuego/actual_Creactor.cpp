#include <actual_Creactor.h> 

/**********************************/
/* Global Variables */
  N_Vector y         = NULL;
  SUNLinearSolver LS = NULL;
  SUNMatrix A = NULL;
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
	y = N_VNew_Serial(neq_tot);
	if(check_flag((void*)y, "N_VNew_Serial", 0)) return(1);

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
        atol  = N_VNew_Serial(neq_tot);
	ratol = N_VGetArrayPointer(atol);
        for (int i=0; i<neq_tot; i++) {
            ratol[i] = 1.0e-10;
        }
	/* Call CVodeSVtolerances to specify the scalar relative tolerance
	 * and vector absolute tolerances */
	flag = CVodeSVtolerances(cvode_mem, reltol, atol);
	if (check_flag(&flag, "CVodeSVtolerances", 1)) return(1);

	if (iDense_Creact == 1) {
            printf("\n--> Using a Direct Dense Solver \n");

            /* Create dense SUNMatrix for use in linear solves */
	    A = SUNDenseMatrix(neq_tot, neq_tot);
	    if(check_flag((void *)A, "SUNDenseMatrix", 0)) return(1);

	    /* Create dense SUNLinearSolver object for use by CVode */
	    LS = SUNDenseLinearSolver(y, A);
	    if(check_flag((void *)LS, "SUNDenseLinearSolver", 0)) return(1);

	    /* Call CVDlsSetLinearSolver to attach the matrix and linear solver to CVode */
	    flag = CVDlsSetLinearSolver(cvode_mem, LS, A);
	    if(check_flag(&flag, "CVDlsSetLinearSolver", 1)) return(1);

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
	    if (iDense_Creact == 1) {
                printf("\n--> With Analytical J\n");
                flag = CVDlsSetJacFn(cvode_mem, cJac);
	        if(check_flag(&flag, "CVDlsSetJacFn", 1)) return(1);
	    } else {
	        amrex::Abort("\n--> No Analytical J for Iterative Solvers");
	    }
	}

	flag = CVodeSetMaxNumSteps(cvode_mem, 10000);
	if(check_flag(&flag, "CVodeSetMaxNumSteps", 1)) return(1);

	/* Define vectors to be used later in creact */
	// GPU stuff: might want to rethink this and put everything in userdata
	if (iE_Creact == 1) { 
	    rhoe_init = (double *) malloc(NCELLS*sizeof(double));
	    rhoesrc_ext = (double *) malloc( NCELLS*sizeof(double));
	} else {
	    rhoh_init = (double *) malloc(NCELLS*sizeof(double));
	    rhohsrc_ext = (double *) malloc( NCELLS*sizeof(double));
	}
	rYsrc = (double *)  malloc((NCELLS*NEQ)*sizeof(double));

	N_VDestroy(atol);          /* Free the atol vector */

	/* Ok we're done ...*/
        if (iverbose > 1) {
	    printf(" --> DONE WITH INITIALIZATION, %d \n", iE_Creact);
	}

	return(0);
}

/* Main routine for external looping */
int actual_cReact(realtype *rY_in, realtype *rY_src_in, 
		realtype *rX_in, realtype *rX_src_in,
		realtype *P_in, realtype *dt_react, realtype *time){

	realtype time_init, time_out ;
	int flag;

        time_init = *time;
	time_out  = *time + *dt_react;

	/* Get Device MemCpy of in arrays */
	/* Get Device pointer of solution vector */
	realtype *yvec_d      = N_VGetArrayPointer(y);
	// rhoY,T
	std::memcpy(yvec_d, rY_in, sizeof(realtype) * ((NEQ+1)*NCELLS));
	// rhoY_src_ext
	std::memcpy(rYsrc, rY_src_in, (NEQ*NCELLS)*sizeof(double));
	// rhoE/rhoH
	if (iE_Creact == 1) { 
	    std::memcpy(rhoe_init, rX_in, sizeof(realtype) * NCELLS);
	    std::memcpy(rhoesrc_ext, rX_src_in, sizeof(realtype) * NCELLS);
	} else {
	    std::memcpy(rhoh_init, rX_in, sizeof(realtype) * NCELLS);
	    std::memcpy(rhohsrc_ext, rX_src_in, sizeof(realtype) * NCELLS);
	}

	/* Call CVODE: ReInit for convergence */
	CVodeReInit(cvode_mem, time_init, y);
	flag = CVode(cvode_mem, time_out, y, &time_init, CV_NORMAL);
	if (check_flag(&flag, "CVode", 1)) return(1);

	/* Pack data to return in main routine external */
	std::memcpy(rY_in, yvec_d, ((NEQ+1)*NCELLS)*sizeof(realtype));
	for  (int i = 0; i < NCELLS; i++) {
            rX_in[i] = rX_in[i] + (*dt_react) * rX_src_in[i];
	}

	/* tests HP */
        for (int tid = 0; tid < NCELLS; tid ++) {
	    double rhov, energy, temp;
	    double MF[NEQ];
            //realtype activity[NEQ], cdot[NEQ], molecular_weight[NEQ];
            int * iwrk;
            double *  rwrk;
	    int  lierr;
	    rhov = 0.0;
            int offset = tid * (NEQ + 1); 
            for (int k = 0; k < NEQ; k ++) {
		rhov =  rhov + rY_in[offset + k];
	    }
            //ckwt_(iwrk, rwrk, molecular_weight);
            for (int k = 0; k < NEQ; k ++) {
		MF[k] = rY_in[offset + k]/rhov;
	        //activity[k] = rY_in[offset + k]/(molecular_weight[k]);
	    }
	    energy = rX_in[tid]/rhov ;
	    if (iE_Creact == 1) { 
	        get_t_given_ey_(&energy, MF, iwrk, rwrk, &temp, &lierr);
	    } else {
	        get_t_given_hy_(&energy, MF, iwrk, rwrk, &temp, &lierr);
	    }
	    rY_in[offset + NEQ] =  temp;
	    // DEBUG CHEKS
            //ckwc_(&temp, activity, iwrk, rwrk, cdot);
	    //*P_in = cdot[2];
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

	realtype *y_d      = N_VGetArrayPointer(y_in);
	realtype *ydot_d   = N_VGetArrayPointer(ydot_in);

        if (iE_Creact == 1) {
	    fKernelSpec(&t, y_d, ydot_d, 
			    rhoe_init, rhoesrc_ext, rYsrc);
	} else {
	    fKernelSpec(&t, y_d, ydot_d, 
			    rhoh_init, rhohsrc_ext, rYsrc);
	}
	//end = std::chrono::system_clock::now();
	//std::chrono::duration<double> elapsed_seconds = end - start;
	//std::cout << " RHS duration " << elapsed_seconds.count() << std::endl;
	return(0);
}

/*
 * CUDA kernels
 */
void fKernelSpec(realtype *dt, realtype *yvec_d, realtype *ydot_d,  
		            double *rhoX_init, double *rhoXsrc_ext, double *rYs)
{
  int tid;
  for (tid = 0; tid < NCELLS; tid ++) {
      realtype massfrac[NEQ],activity[NEQ];
      realtype Xi[NEQ], cXi[NEQ];
      realtype cdot[NEQ], molecular_weight[NEQ];
      realtype temp, energy;
      int lierr;
      int * iwrk;
      double * rwrk;

      int offset = tid * (NEQ + 1); 

      /* MW CGS */
      ckwt_(iwrk, rwrk, molecular_weight);
      /* rho MKS */ 
      realtype rho = 0.0;
      for (int i = 0; i < NEQ; i++){
          rho = rho + yvec_d[offset + i];
      }
      /* temp */
      temp = yvec_d[offset + NEQ];
      /* Yks, C CGS*/
      for (int i = 0; i < NEQ; i++){
          massfrac[i] = yvec_d[offset + i] / rho;
	  activity[i] = yvec_d[offset + i]/(molecular_weight[i]);
      }
      /* NRG CGS */
      energy = (rhoX_init[tid] + rhoXsrc_ext[tid]*(*dt)) /rho;

      /* Fuego calls on device */
      if (iE_Creact == 1){
          get_t_given_ey_(&energy, massfrac, iwrk, rwrk, &temp, &lierr);
          ckums_(&temp, iwrk, rwrk, Xi);
          ckcvms_(&temp, iwrk, rwrk, cXi);
      } else {
          get_t_given_hy_(&energy, massfrac, iwrk, rwrk, &temp, &lierr);
          ckhms_(&temp, iwrk, rwrk, Xi);
          ckcpms_(&temp, iwrk, rwrk, cXi);
      }
      ckwc_(&temp, activity, iwrk, rwrk, cdot);
      int cX = 0.0;
      for (int i = 0; i < NEQ; i++){
          cX = cX + massfrac[i] * cXi[i];
      }

      /* Fill ydot vect */
      ydot_d[offset + NEQ] = rhoXsrc_ext[tid];
      for (int i = 0; i < NEQ; i++){
          ydot_d[offset + i] = cdot[i] * molecular_weight[i] + rYs[tid * (NEQ) + i];
          ydot_d[offset + NEQ] = ydot_d[offset + NEQ]  - ydot_d[offset + i] * Xi[i];
      }
      ydot_d[offset + NEQ] = ydot_d[offset + NEQ] /(rho * cX);
  }
}

/* Analytical Jacobian evaluation */
static int cJac(realtype tn, N_Vector u, N_Vector fu, SUNMatrix J,
		void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){

  realtype *ydata  = N_VGetArrayPointer(u);

  int tid;
  for (tid = 0; tid < NCELLS; tid ++) {
      realtype *J_col_k;
      realtype rho, temp, energy;
      realtype massfrac[NEQ], activity[NEQ], molecular_weight[NEQ];
      realtype Jmat_tmp[(NEQ+1)*(NEQ+1)];
      double * rwrk;
      int * iwrk;

      int offset = tid * (NEQ + 1); 
  
      /* MW CGS */
      ckwt_(iwrk, rwrk, molecular_weight);
      /* rho MKS */ 
      rho = 0.0;
      for (int i = 0; i < NEQ; i++){
          rho = rho + ydata[offset + i];
      }
      /* temp */
      temp = ydata[offset + NEQ];
      for (int i = 0; i < NEQ; i++){
          massfrac[i] = ydata[offset + i]/rho;
          activity[i] = ydata[offset + i]/(molecular_weight[i]);
      }
      /* NRG CGS */
      if (iE_Creact == 1) {
            energy = (rhoe_init[tid] + rhoesrc_ext[tid] * tn) /rho;
	    int consP = 0 ;
            dwdot_(Jmat_tmp, activity, &temp, &consP);
      } else {
            energy = (rhoh_init[tid] + rhohsrc_ext[tid] * tn) /rho;
	    int consP = 1 ;
            dwdot_(Jmat_tmp, activity, &temp, &consP);
      }
      /* fill the sunMat */
      for (int k = 0; k < NEQ; k++){
	    J_col_k = SM_COLUMN_D(J,offset + k);
	    for (int i = 0; i < NEQ; i++){
	        J_col_k[offset + i] = Jmat_tmp[k*(NEQ+1)+i] * molecular_weight[i] / molecular_weight[k]; 
	    }
	    J_col_k[offset + NEQ] = Jmat_tmp[k*(NEQ+1)+NEQ] *  rho / molecular_weight[k]; 
      }
      J_col_k = SM_COLUMN_D(J,offset + NEQ);
      for (int i = 0; i < NEQ; i++){
       	    J_col_k[offset + i] = Jmat_tmp[NEQ*(NEQ+1)+i] * molecular_weight[i] / rho; 
      }
  }

  return(0);

}


 /* Free and destroy memory */
void extern_cFree(){

	printf("In cFree\n");
	SUNLinSolFree(LS);
	N_VDestroy(y);          /* Free the y vector */
	CVodeFree(&cvode_mem);
	if (iDense_Creact == 1) {
	    SUNMatDestroy(A);
	}
	//Free(rYsrc);
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

