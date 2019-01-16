#include <actual_Creactor.h> 

/**********************************/

/* Global Variables */
  N_Vector y  = NULL;
  SUNLinearSolver LS = NULL;
  SUNMatrix A = NULL;
  long int nsetups_old = 0;
  long int nst_old =  0;
  long int nfe_old =  0;
  long int nni_old =  0;
  long int ncfn_old = 0;
  long int nfeLS_old = 0;
  long int nje_old = 0;
  long int npe_old = 0;
  long int nps_old = 0;
  int NEQ    = 0;
  int NCELLS    = 0;
  int iDense_Creact = 1;
  int iJac_Creact   = 0;
  int iE_Creact     = 1;
  int iverbose      = 3;
  void *cvode_mem   = NULL;
  double *rhoe_init = NULL;
  double *rhoh_init = NULL;
  double *rhoesrc_ext = NULL;
  double *rhohsrc_ext = NULL;
  double *rYsrc = NULL;
  double temp_old = 0;
  bool InitPartial = false;
  bool FirstTimePrecond = true;
  UserData data = NULL;
  //std::vector< double > arr_time_PSolve;
  std::chrono::duration<double> elapsed_seconds;
  std::chrono::duration<double> elapsed_seconds_RHS;
  std::chrono::duration<double> elapsed_seconds_Pcond;
  //booleantype jok_hack;

/**********************************/
/* Definitions */
/* Initialization routine, called once at the begining of the problem */
int extern_cInit(const int* cvode_meth,const int* cvode_itmeth, 
		const int* cvode_iJac, const int* cvode_iE,
		const int* cvode_iDense, const int* Ncells){

	int flag;
	realtype reltol, time;
	N_Vector atol;
	realtype *ratol;
	double rwrk;
	int iwrk;
	int mm, ii, nfit;
	int neq_tot;

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
	if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);

	/* Call CVodeCreate to create the solver memory and specify the
	 * Backward Differentiation Formula and the use of a Newton iteration */
	if ((*cvode_meth == 2) && (*cvode_itmeth == 2))
	{
	    cvode_mem = CVodeCreate(CV_BDF);
	    if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);
	} else {
	    amrex::Abort("\n--> Weird inputs to CVodeCreate. Viable options are CV_BDF (=2), CV_NEWTON (=2)\n");
	}

        /* Does not work for more than 1 cell right now */
	data = AllocUserData();
	if(check_flag((void *)data, "AllocUserData", 2)) return(1);

	/* Set the pointer to user-defined data */
	flag = CVodeSetUserData(cvode_mem, data);
	if(check_flag(&flag, "CVodeSetUserData", 1)) return(1);   

        time = 0.0e+0;
	/* Call CVodeInit to initialize the integrator memory and specify the
	 * user's right hand side function, the inital time, and 
	 * initial dependent variable vector y. */
	flag = CVodeInit(cvode_mem, cF_RHS, time, y);
	if (check_flag(&flag, "CVodeInit", 1)) return(1);
	
	/* Definition of tolerances: one for each species */
	reltol = 1.0e-10;
        atol  = N_VNew_Serial(neq_tot);
	ratol = N_VGetArrayPointer(atol);
        for (int i=0; i<neq_tot; i++) {
            ratol[i] = 1.0e-10;
        }
	/* Call CVodeSVtolerances to specify the scalar relative tolerance
	 * and vector absolute tolerances */
	flag = CVodeSVtolerances(cvode_mem, reltol, atol);
	if (check_flag(&flag, "CVodeSVtolerances", 1)) return(1);

	//flag = CVodeSetInitStep(cvode_mem, 1.0e-08);
	//if (check_flag(&flag, "CVodeSetInitStep", 1)) return(1);

	flag = CVodeSetNonlinConvCoef(cvode_mem, 1.0e-03);
	if (check_flag(&flag, "CVodeSetNonlinConvCoef", 1)) return(1);

	flag = CVodeSetMaxNonlinIters(cvode_mem, 15);
	if (check_flag(&flag, "CVodeSetMaxNonlinIters", 1)) return(1);

	flag = CVodeSetMaxErrTestFails(cvode_mem, 100);
	if (check_flag(&flag, "CVodeSetMaxErrTestFails", 1)) return(1);

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
	    if (iJac_Creact == 0) {
	        LS = SUNSPGMR(y, PREC_NONE, 0);
	        if(check_flag((void *)LS, "SUNDenseLinearSolver", 0)) return(1);
	    } else {
	        LS = SUNSPGMR(y, PREC_LEFT, 0);
	        if(check_flag((void *)LS, "SUNDenseLinearSolver", 0)) return(1);
	    }

	    /* Set CVSpils linear solver to LS */
	    flag = CVSpilsSetLinearSolver(cvode_mem, LS);
	    if(check_flag(&flag, "CVSpilsSetLinearSolver", 1)) return(1);
	} else {
	    int nJdata;
	    int HP = 0;
	    sparsity_info_(&nJdata, &HP);
            printf("--> SPARSE solver -- non zero entries %d represents %f %% sparsity pattern.", nJdata, nJdata/float((NEQ+1) * (NEQ+1)) *100.0);
	    amrex::Abort("\n--> Direct Sparse / Iterative Solvers not yet implemented ...\n");
	}

	if (iJac_Creact == 0) {
            printf("\n--> Without Analytical J\n");
	} else {
            printf("\n--> With Analytical J\n");
	    if (iDense_Creact == 99) {
                if (iverbose > 1) {
                    printf("\n    (99)\n");
		}
	        /* Set the JAcobian-times-vector function */
	        flag = CVSpilsSetJacTimes(cvode_mem, NULL, NULL);
	        if(check_flag(&flag, "CVSpilsSetJacTimes", 1)) return(1);
	        /* Set the preconditioner solve and setup functions */
#ifdef USE_KLU 
	        flag = CVSpilsSetPreconditioner(cvode_mem, Precond_sparse, PSolve_sparse);
	        if(check_flag(&flag, "CVSpilsSetPreconditioner", 1)) return(1);
#else
	        flag = CVSpilsSetPreconditioner(cvode_mem, Precond, PSolve);
	        if(check_flag(&flag, "CVSpilsSetPreconditioner", 1)) return(1);
#endif
	    } else {
                if (iverbose > 1) {
                    printf("\n    (1)\n");
		}
	        /* Set the user-supplied Jacobian routine Jac */
                flag = CVDlsSetJacFn(cvode_mem, cJac);
		if(check_flag(&flag, "CVDlsSetJacFn", 1)) return(1);
	    }
	}

        /* Set the max number of time steps */ 
	flag = CVodeSetMaxNumSteps(cvode_mem, 10000);
	if(check_flag(&flag, "CVodeSetMaxNumSteps", 1)) return(1);

        /* Set the max order */ 
        flag = CVodeSetMaxOrd(cvode_mem, 2);
	if(check_flag(&flag, "CVodeSetMaxOrd", 1)) return(1);

	/* Define vectors to be used later in creact */
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
	    printf(" --> DONE WITH INITIALIZATION (CPU) %d \n", iE_Creact);
	}

	return(0);
}

/* Main CVODE call routine */
int actual_cReact(realtype *rY_in, realtype *rY_src_in, 
		realtype *rX_in, realtype *rX_src_in,
		realtype *P_in, 
                realtype *dt_react, realtype *time, int *Init){

        std::chrono::time_point<std::chrono::system_clock> start, end;		
        std::chrono::duration<double> total_elapsed;
	realtype time_init, time_out, dummy_time, temperature_save ;
	int flag;

        //FirstTimePrecond = true;

        time_init = *time;
	time_out  = *time + (*dt_react);
	start = std::chrono::system_clock::now();
        elapsed_seconds = start - start;
        elapsed_seconds_RHS = start - start;
        elapsed_seconds_Pcond = start - start;

	printf("BEG : time curr is %14.6e and dt_react is %14.6e and final time should be %14.6e \n", time_init, *dt_react, time_out);

	/* Get Device MemCpy of in arrays */
	/* Get Device pointer of solution vector */
	realtype *yvec_d      = N_VGetArrayPointer(y);
	// rhoY,T
	std::memcpy(yvec_d, rY_in, sizeof(realtype) * ((NEQ+1)*NCELLS));
	temperature_save = rY_in[NEQ];
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
            printf("\n -------------------------------------\n");
	if (*Init == 1) {
            printf("ReInit always \n");
	    CVodeReInit(cvode_mem, time_init, y);
	    InitPartial = false;
	} else {
	    temp_old = abs(rY_in[NEQ] - temp_old);
	    // Sloppy but I can't think of anything better now
            if (temp_old > 50.0) {
	        printf("ReInit delta_T = %f \n", temp_old);
	        CVodeReInit(cvode_mem, time_init, y);
		InitPartial = false;
	    } else {
	        printf("ReInit Partial delta_T = %f \n", temp_old);
	        CVodeReInitPartial(cvode_mem, time_init, y);
		InitPartial = true;
	    }
	}
	flag = CVode(cvode_mem, time_out, y, &dummy_time, CV_NORMAL);
	//flag = CVode(cvode_mem, time_out, y, &dummy_time, CV_ONE_STEP);
	if (check_flag(&flag, "CVode", 1)) return(1);

	//CHECK THIS
	//CVodeGetCurrentTime(cvode_mem, time);
	//*time = time_out;
	*dt_react = dummy_time - time_init;
	printf("END : time curr is %14.6e and actual dt_react is %14.6e \n", dummy_time, *dt_react);

	/* Pack data to return in main routine external */
	std::memcpy(rY_in, yvec_d, ((NEQ+1)*NCELLS)*sizeof(realtype));
	for  (int i = 0; i < NCELLS; i++) {
            rX_in[i] = rX_in[i] + (*dt_react) * rX_src_in[i];
	}

        end = std::chrono::system_clock::now();
        total_elapsed = end - start;

	if (*Init != 1) {
	    temp_old = rY_in[NEQ];
	}

	/* VERBOSE */
        if (iverbose > 5) {
            for (int tid = 0; tid < NCELLS; tid ++) {
	        double rhov, energy, temp, energy2;
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
	            ckhbms_(&temp, MF, iwrk, rwrk, &energy2);
	            ckubms_(&temp, MF, iwrk, rwrk, &energy);
	    	ckpy_(&rhov, &temp, MF, iwrk, rwrk, P_in);
	            printf("e,h,p,rho ? %4.16e %4.16e %4.16e %4.16e \n",energy, energy2, *P_in, rhov);
	        } else {
	            get_t_given_hy_(&energy, MF, iwrk, rwrk, &temp, &lierr);
	            ckhbms_(&temp, MF, iwrk, rwrk, &energy);
	            ckubms_(&temp, MF, iwrk, rwrk, &energy2);
	    	ckpy_(&rhov, &temp, MF, iwrk, rwrk, P_in);
	            printf("e,h,p,rho ? %4.16e %4.16e %4.16e %4.16e\n",energy2, energy, *P_in, rhov);
	        }
	        //rY_in[offset + NEQ] =  temp;
	        // DEBUG CHEKS
                //ckwc_(&temp, activity, iwrk, rwrk, cdot);
	        // *P_in = cdot[2];
	    }

	    PrintFinalStats(cvode_mem, rY_in[NEQ], InitPartial);

	} else if (iverbose > 2) {
	    printf("\nAdditional verbose info --\n");
	    PrintFinalStats(cvode_mem, temperature_save, InitPartial);
	    std::cout << "Temp, chemistry solve, RHSeval, PSolve, Precond = " << rY_in[NEQ] << " " << total_elapsed.count() << " "<< elapsed_seconds_RHS.count() << " " << elapsed_seconds.count() << " " << elapsed_seconds_Pcond.count() << std::endl; 
	    std::cout << "Temp, RHSeval represnts, PSolve represents, Precond represents = " << rY_in[NEQ] << " " << elapsed_seconds_RHS.count() / total_elapsed.count() * 100.0 <<  " " << elapsed_seconds.count()/total_elapsed.count() * 100.0 << " " << elapsed_seconds_Pcond.count() /total_elapsed.count() * 100.0 << std::endl;
            printf(" -------------------------------------\n");
	}

        long int nfe;
	flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
	return nfe;
}

/* RHS routine used in CVODE */
static int cF_RHS(realtype t, N_Vector y_in, N_Vector ydot_in, 
		void *user_data){

        std::chrono::time_point<std::chrono::system_clock> start, end;		
	start = std::chrono::system_clock::now();

	realtype *y_d      = N_VGetArrayPointer(y_in);
	realtype *ydot_d   = N_VGetArrayPointer(ydot_in);

        if (iE_Creact == 1) {
	    fKernelSpec(&t, y_d, ydot_d, 
			    rhoe_init, rhoesrc_ext, rYsrc);
	} else {
	    fKernelSpec(&t, y_d, ydot_d, 
			    rhoh_init, rhohsrc_ext, rYsrc);
	}
	end = std::chrono::system_clock::now();
        elapsed_seconds_RHS = elapsed_seconds_RHS + end - start;
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
      realtype  temp; 
      realtype activity[NEQ], molecular_weight[NEQ];
      realtype Jmat_tmp[(NEQ+1)*(NEQ+1)];
      double rwrk;
      int iwrk;

      int offset = tid * (NEQ + 1); 

      /* MW CGS */
      ckwt_(&iwrk, &rwrk, molecular_weight);
      /* temp */
      temp = ydata[offset + NEQ];
      for (int i = 0; i < NEQ; i++){
          activity[i] = ydata[offset + i]/(molecular_weight[i]);
      }
      /* NRG CGS */
      int consP;
      if (iE_Creact == 1) {
	  consP = 0;
          dwdot_(Jmat_tmp, activity, &temp, &consP);
      } else {
          consP = 1;
          dwdot_(Jmat_tmp, activity, &temp, &consP);
      }
      /* fill the sunMat */
      for (int k = 0; k < NEQ; k++){
	  J_col_k = SM_COLUMN_D(J,offset + k);
	  for (int i = 0; i < NEQ; i++){
	        J_col_k[offset + i] = Jmat_tmp[k*(NEQ+1)+i] * molecular_weight[i] / molecular_weight[k]; 
          }
	  J_col_k[offset + NEQ] = Jmat_tmp[k*(NEQ+1)+NEQ] / molecular_weight[k]; 
      }
      J_col_k = SM_COLUMN_D(J,offset + NEQ);
      for (int i = 0; i < NEQ; i++){
          J_col_k[offset + i] = Jmat_tmp[NEQ*(NEQ+1)+i] * molecular_weight[i]; 
      }
  }

  return(0);

}

/* Jacobian-times-vector routine. */
static int jtv(N_Vector v, N_Vector Jv, realtype t, N_Vector u, N_Vector fu, 
			  void *user_data, N_Vector tmp)
{
  realtype *udata, *vdata, *Jvdata;
  vdata  = N_VGetArrayPointer(v); 
  udata  = N_VGetArrayPointer(u);
  Jvdata = N_VGetArrayPointer(Jv);

  int tid;
  for (tid = 0; tid < NCELLS; tid ++) {
	realtype temp;
	realtype activity[NEQ], molecular_weight[NEQ];
	realtype J[(NEQ+1)*(NEQ+1)];
	double rwrk;
	int iwrk;

        int offset = tid * (NEQ + 1); 

        /* MW CGS */
	ckwt_(&iwrk, &rwrk, molecular_weight);
	for (int i = 0; i < NEQ; i++){
            activity[i] = udata[offset + i]/(molecular_weight[i]);
	}
        /* temp */
	temp = udata[offset + NEQ];
        /* NRG CGS */
        if (iE_Creact == 1) {
            int consP = 0;
            dwdot_(J, activity, &temp, &consP);
        } else {
	    int consP = 1 ;
            dwdot_(J, activity, &temp, &consP);
        }

	/* PRINT JAC INFO: debug mode
	for (int i = 0; i < NEQ+1; i++){
	    for (int j = 0; j < NEQ+1; j++){
                printf(" %4.8e ", (J[j*(NEQ+1)+i]));
	    }
            printf(" \n");
	}
	amrex::Abort("\n--> ABORT\n");
	*/

	/* reinit */
	for (int i = 0; i < NEQ+1; i++){
	    Jvdata[offset + i] = 0.0;
	}
	/* COmpute */
	for (int i = 0; i < NEQ; i++){
            for (int j = 0; j < NEQ; j++){
                Jvdata[offset + i] = Jvdata[offset + i] + J[j*(NEQ+1)+i] * vdata[offset + j] * molecular_weight[i] / molecular_weight[j];
	    }
            Jvdata[offset + i] = Jvdata[offset + i] + J[NEQ*(NEQ+1)+i] * vdata[offset + NEQ] * molecular_weight[i];
        }
	for (int j = 0; j < NEQ; j++){
	    Jvdata[offset + NEQ] = Jvdata[offset + NEQ] + J[j*(NEQ+1)+NEQ] * vdata[offset + j] / molecular_weight[j];
	}
	Jvdata[offset + NEQ] = Jvdata[offset + NEQ]  + J[(NEQ+1)*(NEQ+1)-1] * vdata[offset + NEQ] ;  
  }

  return(0);
}

/* Preconditioner SPARSE setup routine. Generate and preprocess P. */
#ifdef USE_KLU 
static int Precond_sparse(realtype tn, N_Vector u, N_Vector fu, booleantype jok, 
		booleantype *jcurPtr, realtype gamma, void *user_data)
{
  std::chrono::time_point<std::chrono::system_clock> start, end;		
  std::chrono::time_point<std::chrono::system_clock> start_Jcomp, end_Jcomp;		
  std::chrono::time_point<std::chrono::system_clock> start_JNcomp, end_JNcomp;		
  std::chrono::time_point<std::chrono::system_clock> start_LUfac, end_LUfac;		
  
  //std::chrono::duration<double> elapsed_seconds_Jcomp;
  //std::chrono::duration<double> elapsed_seconds_JNcomp;
  //std::chrono::duration<double> elapsed_seconds_LUfac;
  //std::chrono::duration<double> elapsed_seconds_Pcond_prov;

  UserData data_wk;
  realtype **(*Jbd)[1];
  SUNMatrix PS;
  realtype *udata;
  realtype temp;
  realtype Jmat[(NEQ+1)*(NEQ+1)];
  realtype activity[NEQ], molecular_weight[NEQ];
  double rwrk;
  int iwrk, ok;

  start = std::chrono::system_clock::now();

  /* Make local copies of pointers in user_data, and of pointer to u's data */
  data_wk = (UserData) user_data;   
  PS = (data_wk->PS[0][0]);
  Jbd = (data_wk->Jbd);
  udata = N_VGetArrayPointer(u);

  //start_Jcomp = std::chrono::system_clock::now();
  if (jok) {
      /* jok = SUNTRUE: Copy Jbd to P */
      //printf("\n Reuse J ... ");
      *jcurPtr = SUNFALSE;
  } else {
	  //printf("\n Recompute J from scratch ... ");
          ckwt_(&iwrk, &rwrk, molecular_weight);
          for (int i = 0; i < NEQ; i++){
              activity[i] = udata[i]/(molecular_weight[i]);
          }
          temp = udata[NEQ];

          // C in mol/cm3
          int consP;
          if (iE_Creact == 1) { 
              consP = 0;
          } else {
              consP = 1;
          }
          dwdot_precond_(Jmat, activity, &temp, &consP);

          /* Compute Jacobian.  Load into P. */
          denseScale(0.0, Jbd[0][0], NEQ+1, NEQ+1);

          for (int i = 0; i < NEQ; i++) {
              for (int k = 0; k < NEQ; k++) {
                  (Jbd[0][0])[k][i] = Jmat[k*(NEQ+1) + i] * molecular_weight[i] / molecular_weight[k];
              }
              (Jbd[0][0])[i][NEQ] = Jmat[i*(NEQ+1) + NEQ] / molecular_weight[i];
          }
          for (int i = 0; i < NEQ; i++) {
              (Jbd[0][0])[NEQ][i] = Jmat[NEQ*(NEQ+1) + i] * molecular_weight[i];
          }
	  (Jbd[0][0])[NEQ][NEQ] = Jmat[(NEQ+1)*(NEQ+1)-1];

          *jcurPtr = SUNTRUE;
  }
  //end_Jcomp = std::chrono::system_clock::now();
  //elapsed_seconds_Jcomp = end_Jcomp - start_Jcomp;

  //start_JNcomp = std::chrono::system_clock::now();
  int nbVals;
  for (int i = 1; i < NEQ+2; i++) {
	  nbVals = data_wk->colPtrs[i]-data_wk->colPtrs[i-1];
	  //printf("Rows %d : we have %d nonzero values \n", i-1, nbVals);
	  for (int j = 0; j < nbVals; j++) {
		  int idx = data_wk->rowVals[ data_wk->colPtrs[i-1] + j ];
                  /* Scale by -gamma */
                  /* Add identity matrix */
		  if (idx == (i-1)) {
	              data_wk->Jdata[ data_wk->colPtrs[i-1] + j ] = 1.0 - gamma * (Jbd[0][0])[ i-1 ][ idx ]; 
		      //printf(" -- idx %d %4.16e %4.16e \n", idx, Jmat[(i-1)*(NEQ+1) + idx], (Jbd[0][0])[ i-1 ] [ idx ]);
		  } else {
	              data_wk->Jdata[ data_wk->colPtrs[i-1] + j ] = - gamma * (Jbd[0][0])[ i-1 ][ idx ]; 
		  }
	  }
  }
  //end_JNcomp = std::chrono::system_clock::now();
  //elapsed_seconds_JNcomp = end_JNcomp - start_JNcomp;

  //start_LUfac = std::chrono::system_clock::now();
  if (!FirstTimePrecond) {
  //if (jok) {
      //printf("and reuse pivots ...");
      ok = klu_refactor(data_wk->colPtrs, data_wk->rowVals, data_wk->Jdata, data_wk->Symbolic, data_wk->Numeric, &(data_wk->Common));
  } else {
      //printf("and compute pivots ...");
      data_wk->Numeric = klu_factor  (data_wk->colPtrs, data_wk->rowVals, data_wk->Jdata, data_wk->Symbolic, &(data_wk->Common)) ; 
      FirstTimePrecond = false;
  }
  //end_LUfac = std::chrono::system_clock::now();
  //elapsed_seconds_LUfac = end_LUfac - start_LUfac;

  end = std::chrono::system_clock::now();
  elapsed_seconds_Pcond = elapsed_seconds_Pcond + end - start;

  //elapsed_seconds_Pcond_prov = end - start;
  //std::cout << " stats (Jcomp,JNcomp,pivots,TOTAL) = " << elapsed_seconds_Jcomp.count() <<  " " << elapsed_seconds_JNcomp.count() << " " << elapsed_seconds_LUfac.count() << " " << elapsed_seconds_Pcond_prov.count() << std::endl;

  return(0);
}
#endif

/* Preconditioner setup routine. Generate and preprocess P. */
static int Precond(realtype tn, N_Vector u, N_Vector fu, booleantype jok, 
		booleantype *jcurPtr, realtype gamma, void *user_data)
{
  std::chrono::time_point<std::chrono::system_clock> start, end;		
  start = std::chrono::system_clock::now();

  UserData data_wk;
  realtype **(*P)[1], **(*Jbd)[1];
  sunindextype *(*pivot)[1], ierr;
  realtype *udata, **a, **j;
  realtype temp;
  realtype Jmat[(NEQ+1)*(NEQ+1)];
  realtype activity[NEQ], molecular_weight[NEQ];
  double rwrk;
  int iwrk;

  /* Make local copies of pointers in user_data, and of pointer to u's data */
  data_wk = (UserData) user_data;   
  P = (data_wk->P);
  Jbd = (data_wk->Jbd);
  pivot = (data_wk->pivot);
  udata = N_VGetArrayPointer(u);

  if (jok) {
      /* jok = SUNTRUE: Copy Jbd to P */
      denseCopy(Jbd[0][0], P[0][0], NEQ+1, NEQ+1);
      *jcurPtr = SUNFALSE;
  } else {
      /* jok = SUNFALSE: Generate Jbd from scratch and copy to P */
      /* Make local copies of problem variables, for efficiency. */
      ckwt_(&iwrk, &rwrk, molecular_weight);
      for (int i = 0; i < NEQ; i++){
          activity[i] = udata[i]/(molecular_weight[i]);
      }
      temp = udata[NEQ];

      // C in mol/cm3
      int consP;
      if (iE_Creact == 1) { 
          consP = 0;
      } else {
          consP = 1;
      }
      dwdot_precond_(Jmat, activity, &temp, &consP);
      //dwdot_(Jmat, activity, &temp, &consP);
      /* Compute Jacobian.  Load into P. */
      denseScale(0.0, Jbd[0][0], NEQ+1, NEQ+1);
      for (int i = 0; i < NEQ; i++) {
          for (int k = 0; k < NEQ; k++) {
              (Jbd[0][0])[k][i] = Jmat[k*(NEQ+1) + i] * molecular_weight[i] / molecular_weight[k];
              //(Jbd[0][0])[i][k] = Jmat[k*(NEQ+1) + i] * molecular_weight[i] / molecular_weight[k];
          }
          (Jbd[0][0])[i][NEQ] = Jmat[i*(NEQ+1) + NEQ] / molecular_weight[i];
          //(Jbd[0][0])[NEQ][i] = Jmat[i*(NEQ+1) + NEQ] / molecular_weight[i];
      }
      for (int i = 0; i < NEQ; i++) {
          (Jbd[0][0])[NEQ][i] = Jmat[NEQ*(NEQ+1) + i] * molecular_weight[i];
          //(Jbd[0][0])[i][NEQ] = Jmat[NEQ*(NEQ+1) + i] * molecular_weight[i];
      }
      (Jbd[0][0])[NEQ][NEQ] = Jmat[(NEQ+1)*(NEQ+1)-1];

      denseCopy(Jbd[0][0], P[0][0], NEQ+1, NEQ+1);

      *jcurPtr = SUNTRUE;
  }
  
  /* Scale by -gamma */
  denseScale(-gamma, P[0][0], NEQ+1, NEQ+1);
  //denseScale(0.0, P[0][0], NEQ+1, NEQ+1);
  
  /* Add identity matrix and do LU decompositions on blocks in place. */
  denseAddIdentity(P[0][0], NEQ+1);
  ierr = denseGETRF(P[0][0], NEQ+1, NEQ+1, pivot[0][0]);
  if (ierr != 0) return(1);
  
  end = std::chrono::system_clock::now();
  elapsed_seconds_Pcond = elapsed_seconds_Pcond + end - start;

  return(0);
}

#ifdef USE_KLU 
static int PSolve_sparse(realtype tn, N_Vector u, N_Vector fu, N_Vector r, N_Vector z,
                  realtype gamma, realtype delta, int lr, void *user_data)
{
  std::chrono::time_point<std::chrono::system_clock> start, end;		
  start = std::chrono::system_clock::now();

  realtype *zdata;
  UserData data_wk;

  data_wk = (UserData) user_data;

  N_VScale(1.0, r, z);

  zdata = N_VGetArrayPointer(z);
  
  /* Solve the block-diagonal system Pz = r using LU factors stored
     in P and pivot data in pivot, and return the solution in z. */

  klu_solve (data_wk->Symbolic, data_wk->Numeric, NEQ+1, 1, zdata, &(data_wk->Common)) ; 

  end = std::chrono::system_clock::now();
  elapsed_seconds = elapsed_seconds + end - start;

  return(0);
}
#endif

static int PSolve(realtype tn, N_Vector u, N_Vector fu, N_Vector r, N_Vector z,
                  realtype gamma, realtype delta, int lr, void *user_data)
{
  std::chrono::time_point<std::chrono::system_clock> start, end;		
  start = std::chrono::system_clock::now();

  realtype **(*P)[1];
  sunindextype *(*pivot)[1];
  realtype *zdata, *v;
  UserData data_wk;

  /* Extract the P and pivot arrays from user_data. */

  data_wk = (UserData) user_data;
  P = data_wk->P;
  pivot = data_wk->pivot;
  zdata = N_VGetArrayPointer(z);
  
  N_VScale(1.0, r, z);
  
  /* Solve the block-diagonal system Pz = r using LU factors stored
     in P and pivot data in pivot, and return the solution in z. */
  v = zdata;
  denseGETRS(P[0][0], NEQ+1, pivot[0][0], v);

  end = std::chrono::system_clock::now();
  elapsed_seconds = elapsed_seconds + end - start;
  //std::cout << " RHS duration " << elapsed_seconds.count() << std::endl;

  return(0);
}
 
void extern_cFree(){

	// Free y and abstol vectors 
	CVodeFree(&cvode_mem);
	SUNLinSolFree(LS);
	if (iDense_Creact == 1) {
	    SUNMatDestroy(A);
	}
	FreeUserData(data);
}

/* 
 * Get and print some final statistics
 */
static void PrintFinalStats(void *cvodeMem, realtype Temp, bool InitPartial)
{
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
  long int nli, npe, nps, ncfl, netfails;
  int flag;
  realtype hlast, hinused, hcur;

  flag = CVodeGetNumSteps(cvodeMem, &nst);
  check_flag(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetNumRhsEvals(cvodeMem, &nfe);
  check_flag(&flag, "CVodeGetNumRhsEvals", 1);
  flag = CVodeGetNumNonlinSolvIters(cvodeMem, &nni);
  check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumLinSolvSetups(cvodeMem, &nsetups);
  check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvodeMem, &ncfn);
  check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);
  flag = CVodeGetNumErrTestFails(cvodeMem, &netfails);
  check_flag(&flag, "CVodeGetNumErrTestFails", 1);
  flag = CVodeGetLastStep(cvodeMem, &hlast);
  check_flag(&flag, "CVodeGetLastStep", 1);
  flag = CVodeGetActualInitStep(cvodeMem, &hinused);
  check_flag(&flag, "CVodeGetActualInitStep", 1);
  flag = CVodeGetCurrentTime(cvodeMem, &hcur);
  check_flag(&flag, "CVodeGetCurrentTime", 1);

  if (iDense_Creact == 1){
      flag = CVDlsGetNumRhsEvals(cvodeMem, &nfeLS);
      check_flag(&flag, "CVDlsGetNumRhsEvals", 1);
      flag = CVDlsGetNumJacEvals(cvodeMem, &nje);
      check_flag(&flag, "CVDlsGetNumJacEvals", 1);
  } else if (iDense_Creact == 99){
      flag = CVSpilsGetNumRhsEvals(cvode_mem, &nfeLS);
      check_flag(&flag, "CVSpilsGetNumRhsEvals", 1);
      flag = CVSpilsGetNumJtimesEvals(cvodeMem, &nje);
      //flag = CVSpilsGetNumJTSetupEvals(cvodeMem, &nje);
      check_flag(&flag, "CVSpilsGetNumJTSetupEvals", 1);
      flag = CVSpilsGetNumPrecEvals(cvodeMem, &npe);
      check_flag(&flag, "CVSpilsGetNumPrecEvals", 1);
      flag = CVSpilsGetNumPrecSolves(cvodeMem, &nps);
      check_flag(&flag, "CVSpilsGetNumPrecSolves", 1);
      flag = CVSpilsGetNumLinIters(cvodeMem, &nli);
      check_flag(&flag, "CVSpilsGetNumLinIters", 1);
      flag = CVSpilsGetNumConvFails(cvodeMem, &ncfl); 
      check_flag(&flag, "CVSpilsGetNumConvFails", 1);
  }

  printf("-- Final Statistics --\n");
  printf("NonLinear (Newton) related --\n");
  if (InitPartial) {
         printf("    DT(dt-dttrue), RHS, Iterations, ConvFails, LinSolvSetups = %f %-6ld(%14.6e) %-6ld %-6ld %-6ld %-6ld \n",
	 Temp, nst-nst_old, hlast-hinused, nfe-nfe_old, nni-nni_old, ncfn-ncfn_old, nsetups-nsetups_old);
  }else{
         printf("    DT(dt, dtcur), RHS, Iterations, ErrTestFails, LinSolvSetups = %f %-6ld(%14.6e %14.6e) %-6ld %-6ld %-6ld %-6ld \n",
	 Temp, nst, hlast, hcur, nfe, nni, netfails, nsetups);
  }
  /* RESET */
  nsetups_old = nsetups;
  nst_old = nst;
  nfe_old = nfe;
  nni_old = nni;
  ncfn_old = ncfn;

  if (iDense_Creact == 1){
      printf("Linear (Dense Direct Solve) related --\n");
      if (InitPartial) {
          printf("    FD RHS, NumJacEvals                           = %f %-6ld %-6ld \n", Temp, nfeLS-nfeLS_old, nje-nje_old);
      } else {
          printf("    FD RHS, NumJacEvals                           = %f %-6ld %-6ld \n", Temp, nfeLS, nje);
      }
      /* RESET */
      nfeLS_old = nfeLS;
      nje_old = nje;
  } else if (iDense_Creact == 99){
	  // LinSolvSetups actually reflects the number of time the LinSolver has been called. 
	  // NonLinIterations can be taken without the need for LinItes
          //printf("Temp, NumLinIters/LinSolvSetups, NumConvfails                      = %f %f %-6ld \n", Temp, float(nli)/float(nsetups), ncfl);
      printf("Linear (Krylov GMRES Solve) related --\n");
      if (InitPartial) {
          printf("    RHSeval, jtvEval, NumPrecEvals, NumPrecSolves = %f %-6ld %-6ld %-6ld %-6ld \n", 
            	      Temp, nfeLS-nfeLS_old, nje-nje_old, npe-npe_old, nps-nps_old);
      } else {
          printf("    RHSeval, jtvEval, NumPrecEvals, NumPrecSolves = %f %-6ld %-6ld %-6ld %-6ld \n", 
            	      Temp, nfeLS, nje, npe, nps);
          printf("    Iterations, ConvFails = %f %-6ld %-6ld \n", 
            	      Temp, nli, ncfl );
          //printf("Temp, NumLinIters/LinSolvSetups, NumConvfails                      = %f %f %-6ld \n", Temp, float(nli)/float(nsetups), ncfl);
      }
      /* RESET */
      nfeLS_old = nfeLS;
      nje_old = nje;
      npe_old = npe;
      nps_old = nps;
  }
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

static UserData AllocUserData(void)
{
  int nJdata;
  int NNZ;
  UserData data_wk;

  data_wk = (UserData) malloc(sizeof *data_wk);

  (data_wk->P)[0][0] = newDenseMat(NEQ+1, NEQ+1);
  (data_wk->Jbd)[0][0] = newDenseMat(NEQ+1, NEQ+1);
  (data_wk->pivot)[0][0] = newIndexArray(NEQ+1);

#ifdef USE_KLU 
  if (iDense_Creact == 99) {
      int HP;
      if (iE_Creact == 1) {
	  HP = 0;
      } else {
          HP = 1;
      }
      sparsity_info_(&nJdata,&HP);
      NNZ = nJdata;  
      (data_wk->PS)[0][0] = SUNSparseMatrix(NEQ+1, NEQ+1, NNZ, CSC_MAT);
      printf("--> SPARSE Preconditioner -- non zero entries %d represents %f %% sparsity pattern.", nJdata, nJdata/float((NEQ+1) * (NEQ+1)) *100.0);
      //data_wk->colPtrs = (int *) malloc((NEQ+1)*sizeof(int));
      data_wk->colPtrs = (int*) SUNSparseMatrix_IndexPointers((data_wk->PS)[0][0]); 
      //data_wk->rowVals = (int *) malloc((NNZ)*sizeof(int));
      data_wk->rowVals = (int*) SUNSparseMatrix_IndexValues((data_wk->PS)[0][0]);
      //data_wk->Jdata = (double *) malloc((NNZ)*sizeof(double));
      data_wk->Jdata = SUNSparseMatrix_Data((data_wk->PS)[0][0]);
      sparsity_preproc_(data_wk->rowVals,data_wk->colPtrs,&HP);
      klu_defaults (&(data_wk->Common));
      //data_wk->Common.btf = 0;
      data_wk->Common.maxwork = 15;
      //data_wk->Common.ordering = 1;
      data_wk->Symbolic = klu_analyze (NEQ+1, data_wk->colPtrs, data_wk->rowVals, &(data_wk->Common)) ; 
  }
#endif

  return(data_wk);
}

/* Free data memory */
static void FreeUserData(UserData data_wk)
{
  destroyMat((data_wk->P)[0][0]);
  destroyMat((data_wk->Jbd)[0][0]);
  destroyArray((data_wk->pivot)[0][0]);
  free(data_wk);
} 









