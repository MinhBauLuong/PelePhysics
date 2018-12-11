#include <actual_Creactor.h> 

/**********************************/

/* Global Variables */
  N_Vector y  = NULL;
  SUNLinearSolver LS = NULL;
  SUNMatrix A = NULL;
  int NEQ    = 0;
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
  UserData data = NULL;

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

        /* Does not work for more than 1 cell right now */
	data = AllocUserData();
	if(check_flag((void *)data, "AllocUserData", 2)) return(1);

	/* Call CVodeCreate to create the solver memory and specify the
	 * Backward Differentiation Formula and the use of a Newton iteration */
	if ((*cvode_meth == 2) && (*cvode_itmeth == 2))
	{
	    cvode_mem = CVodeCreate(CV_BDF);
	    if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);
	} else {
	    amrex::Abort("\n--> Weird inputs to CVodeCreate. Viable options are CV_BDF (=2), CV_NEWTON (=2)\n");
	}

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
	        flag = CVSpilsSetPreconditioner(cvode_mem, Precond, PSolve);
	        if(check_flag(&flag, "CVSpilsSetPreconditioner", 1)) return(1);
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
                realtype *dt_react, realtype *time){

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

	PrintFinalStats(cvode_mem);

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


/* Preconditioner setup routine. Generate and preprocess P. */
static int Precond(realtype tn, N_Vector u, N_Vector fu, booleantype jok, 
		booleantype *jcurPtr, realtype gamma, void *user_data)
{
  UserData data_wk;
  realtype **(*P)[1], **(*Jbd)[1];
  sunindextype *(*pivot)[1], ierr;
  realtype *udata, **a, **j;
  realtype temp;
  realtype massfrac[NEQ];
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
    realtype rho = 0.0;
    for (int i = 0; i < NEQ; i++){
        rho = rho + udata[i];
    }
    for (int i = 0; i < NEQ; i++){
        massfrac[i] = udata[i]/rho;
    }
    temp = udata[NEQ];

    ckytcr_(&rho, &temp, massfrac, &iwrk, &rwrk, activity);
    // C in mol/cm3
    int consP;
    if (iE_Creact == 1) { 
        consP = 0;
    } else {
        consP = 1;
    }
    dwdot_precond_(Jmat, activity, &temp, &consP);
    ckwt_(&iwrk, &rwrk, molecular_weight);
    /* Compute Jacobian.  Load into P. */
    //printf("GOT HERE %f \n", Jmat[0*(NEQ+1) + 0]);
    denseScale(0.0, Jbd[0][0], NEQ+1, NEQ+1);
    //j = Jbd[0][0];
    //printf("J ? %f \n", j[0][0]);
    //a = P[0][0];
    for (int i = 0; i < NEQ; i++) {
        for (int k = 0; k < NEQ; k++) {
            (Jbd[0][0])[k][i] = Jmat[k*(NEQ+1) + i] * molecular_weight[i] / molecular_weight[k];
            //j[i][k] = Jmat[k*(NEQ+1) + i] * molecular_weight[i] / molecular_weight[k];
        }
        (Jbd[0][0])[i][NEQ] = Jmat[i*(NEQ+1) + NEQ] / molecular_weight[i];
        //j[NEQ][i] = Jmat[i*(NEQ+1) + NEQ] / molecular_weight[i];
    }
    for (int i = 0; i < NEQ; i++) {
        (Jbd[0][0])[NEQ][i] = Jmat[NEQ*(NEQ+1) + i] * molecular_weight[i];
        //j[i][NEQ] = Jmat[NEQ*(NEQ+1) + i] * molecular_weight[i];
    }
      
    denseCopy(Jbd[0][0], P[0][0], NEQ+1, NEQ+1);
    *jcurPtr = SUNTRUE;
  }
  
  /* Scale by -gamma */
  denseScale(-gamma, P[0][0], NEQ+1, NEQ+1);
  
  /* Add identity matrix and do LU decompositions on blocks in place. */
  denseAddIdentity(P[0][0], NEQ+1);
  ierr = denseGETRF(P[0][0], NEQ+1, NEQ+1, pivot[0][0]);
  if (ierr != 0) return(1);

  return(0);
}


static int PSolve(realtype tn, N_Vector u, N_Vector fu, N_Vector r, N_Vector z,
                  realtype gamma, realtype delta, int lr, void *user_data)
{
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

  return(0);
}
 
void extern_cFree(){

	// Print some final statistics 
	PrintFinalStats(cvode_mem);

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

static void PrintFinalStats(void *cvodeMem)
{
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
  long int nli, npe, nps, ncfl;
  int flag;

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

  if (iDense_Creact == 1){
      flag = CVDlsGetNumRhsEvals(cvodeMem, &nfeLS);
      check_flag(&flag, "CVDlsGetNumRhsEvals", 1);
      flag = CVDlsGetNumJacEvals(cvodeMem, &nje);
      check_flag(&flag, "CVDlsGetNumJacEvals", 1);
  } else if (iDense_Creact == 99){
      flag = CVSpilsGetNumRhsEvals(cvode_mem, &nfeLS);
      check_flag(&flag, "CVSpilsGetNumRhsEvals", 1);
      //flag = CVSpilsGetNumJtimesEvals(cvodeMem, &nje);
      flag = CVSpilsGetNumJTSetupEvals(cvodeMem, &nje);
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

  printf("\nFinal Statistics:\n");
  printf("nb dt = %-6ld nb RHS eval = %-6ld nb NonlinSolvIters = %-6ld nb NonlinSolvConvFails = %-6ld \n",
	 nst, nfe, nni, ncfn);
  printf("nb LinSolvSetups = %-6ld \n",
	 nsetups);
  if (iDense_Creact == 1){
      printf("NumRhsEvals (for FD) = %-6ld NumJacEvals  = %-6ld \n", nfeLS, nje);
  } else if (iDense_Creact == 99){
      printf("NumRhsEvals (for FD jtv) = %-6ld NumJacEvals (not sure)  = %-6ld NumPrecEvals = %-6ld NumPrecSolves = %-6ld \n", 
		      nfeLS, nje, npe, nps);
      printf("NumLinIters = %-6ld NumConvfails = %-6ld \n", nli, ncfl);
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
  UserData data_wk;

  data_wk = (UserData) malloc(sizeof *data_wk);

  (data_wk->P)[0][0] = newDenseMat(NEQ+1, NEQ+1);
  (data_wk->Jbd)[0][0] = newDenseMat(NEQ+1, NEQ+1);
  (data_wk->pivot)[0][0] = newIndexArray(NEQ+1);

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









