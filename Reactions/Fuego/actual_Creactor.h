#include <math.h>
#include <iostream>
#include <cassert>
#include <chrono>
#include <ctime>
#include <cstring>

#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sunlinsol/sunlinsol_spgmr.h> /* access to SPGMR SUNLinearSolver     */
#include <cvode/cvode_direct.h>        /* access to CVDls interface            */
#include <cvode/cvode_spils.h>         /* access to CVSpils interface */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include <sundials/sundials_math.h>

#include <AMReX_Print.H>

/**********************************/
/* Functions Called by the Solver */
static int cF_RHS(realtype t, N_Vector y_in, N_Vector ydot, void *user_data);

static int cJac(realtype t, N_Vector y_in, N_Vector fy, SUNMatrix J,
		void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static int jtv(N_Vector v, N_Vector Jv, realtype t, N_Vector u,
		N_Vector fu, void *user_data, N_Vector tmp);

static int PSolve(realtype tn, N_Vector u, N_Vector fu, N_Vector r, N_Vector z,
                  realtype gamma, realtype delta, int lr, void *user_data);

static int Precond(realtype tn, N_Vector u, N_Vector fu, booleantype jok, 
		booleantype *jcurPtr, realtype gamma, void *user_data);

/**********************************/
/* Functions Called by the main program */
int extern_cInit(const int* cvode_meth,const int* cvode_itmeth, 
		const int* cvode_iJac, const int* cvode_iE,
		const int* cvode_iDense, const int* Ncells);

int actual_cReact(realtype *rY_in, realtype *rY_src_in, 
		realtype *rX_in, realtype *rX_src_in, 
		realtype *P_in, realtype *dt_react, realtype *time);

void extern_cFree();

/**********************************/
/* Additional useful functions */
static int check_flag(void *flagvalue, const char *funcname, int opt);

static void PrintFinalStats(void *cvode_mem);

/* Stuff that comes from Fuego on Host */
extern "C" {
    void get_t_given_ey_(double * e, double * y, int * iwrk, double * rwrk, double * t, int *ierr); 
    void get_t_given_hy_(double * h, double * y, int * iwrk, double * rwrk, double * t, int *ierr); 
    void ckpy_(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * P);
    void ckrhoy_(double * P, double * T, double * y, int * iwrk, double * rwrk, double * rho);
    void ckytcr_(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * c);
    void ckcvms_(double * T, int * iwrk, double *  rwrk, double * cvms);
    void ckcpms_(double * T, int * iwrk, double *  rwrk, double * cpms);
    void ckums_(double * T, int * iwrk, double * rwrk, double * ums);
    void ckhms_(double * T, int * iwrk, double * rwrk, double * hms);
    void ckwc_(double * T, double * C, int * iwrk, double * rwrk, double * wdot);
    void ckwt_(int * iwrk, double * rwrk, double * wt);
    void dwdot_(double * J, double * sc, double * Tp, int * consP);
    void ckindx_(int * iwrk, double * rwrk, int * mm, int * kk, int * ii, int * nfit); 
    void ajacobian_diag_(double * J, double * sc, double T, int consP);
}


/**********************************/
/* Main Kernel fct called in solver RHS */
void fKernelSpec(realtype *dt, realtype *yvec_d, realtype *ydot_d,  
		            double *rhoX_init, double *rhoXsrc_ext, double *rYs);
