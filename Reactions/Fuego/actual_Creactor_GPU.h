#include <math.h>
#include <iostream>
#include <cassert>
#include <chrono>
#include <ctime>
#include <assert.h>

#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sunlinsol/sunlinsol_spgmr.h> /* access to SPGMR SUNLinearSolver     */
#include <cvode/cvode_direct.h>        /* access to CVDls interface            */
#include <cvode/cvode_spils.h>         /* access to CVSpils interface */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include <sundials/sundials_math.h>

#include <nvector/nvector_cuda.h>

#include <cusolver/cvode_cusolver_spqr.h>

#include <AMReX_Print.H>

/**********************************/
typedef struct CVodeUserData {
    int ncells_d[1]; 
    int neqs_per_cell[1];
    int flagP;
    double fwd_A[289], fwd_beta[289], fwd_Ea[289];
    double low_A[289], low_beta[289], low_Ea[289];
    double rev_A[289], rev_beta[289], rev_Ea[289];
    double troe_a[289],troe_Ts[289], troe_Tss[289], troe_Tsss[289];
    double sri_a[289], sri_b[289], sri_c[289], sri_d[289], sri_e[289];
    double activation_units[289], prefactor_units[289], phase_units[289]; 
    int is_PD[289], troe_len[289], sri_len[289], nTB[289], *TBid[289];
    double *TB[289];
}* UserData;

/* Functions Called by the Solver */
static int cF_RHS(realtype t, N_Vector y_in, N_Vector ydot, void *user_data);

int fun_csr_jac(realtype t, N_Vector y_in, N_Vector fy_in,  
                CV_cuSolver_csr_sys csr_sys, void* user_data);


/**********************************/
/* Functions Called by the main program */
int extern_cInit(const int* cvode_meth,const int* cvode_itmeth, 
		const int* cvode_iJac, const int* cvode_iE,
		const int* cvode_iDense, const int* Ncells);

int actual_cReact(realtype *rY_in, realtype *rY_src_in, 
		realtype *rX_in, realtype *rX_src_in, 
		realtype *P_in, 
		realtype *dt_react, realtype *time, int *Init);

void extern_cFree();

/**********************************/
/* Additional useful functions */
static void initialize_chemistry_device(UserData user_data);

//static void finalize_chemistry_device(void *user_data);

static int check_flag(void *flagvalue, const char *funcname, int opt);

static void PrintFinalStats(void *cvode_mem);

/* Stuff that comes from Fuego on Host */
//extern "C" {
    //void imolecularWeight_(double * iwt);
    //void ckpy_(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * P);
    void ckindx_(int * iwrk, double * rwrk, int * mm, int * kk, int * ii, int * nfit); 

    //void get_t_given_ey_(double * e, double * y, int * iwrk, double * rwrk, double * t, int *ierr); 
    //void get_t_given_hy_(double * h, double * y, int * iwrk, double * rwrk, double * t, int *ierr); 
    //void ckrhoy_(double * P, double * T, double * y, int * iwrk, double * rwrk, double * rho);
    //void ckytcr_(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * c);
    //void ckcvms_(double * T, int * iwrk, double *  rwrk, double * cvms);
    //void ckcpms_(double * T, int * iwrk, double *  rwrk, double * cpms);
    //void ckums_(double * T, int * iwrk, double * rwrk, double * ums);
    //void ckhms_(double * T, int * iwrk, double * rwrk, double * hms);
    //void ckwc_(double * T, double * C, int * iwrk, double * rwrk, double * wdot);
    //void ckwt_(int * iwrk, double * rwrk, double * wt);
    //void dwdot_(double * J, double * sc, double * Tp, int * consP);
    //void ajacobian_diag_(double * J, double * sc, double T, int consP);
//}


/**********************************/
/* Device crap               */

/* Main Kernel fct called in solver RHS */
__global__ void fKernelSpec(realtype *dt, void *user_data, 
		            realtype *yvec_d, realtype *ydot_d,  
		            double *rhoX_init, double *rhoXsrc_ext, double *rYs);


__global__ void fKernelJacCSR(realtype t, void *user_data,
                                          realtype *yvec_d, realtype *ydot_d,
                                          realtype* csr_jac, 
                                          const int size, const int nnz, 
                                          const int nbatched);


/* FROM FUEGO */
/*save inv molecular weights into array */
__device__ void imolecularWeight_d(double * iwt);
/* Returns R, Rc, Patm */
__device__ void ckrp_d( double * ru, double * ruc, double * pa);
/*Compute P = rhoRT/W(x) */
__device__ void ckpx_d(double * rho, double * T, double * x, double * P);
/*Compute P = rhoRT/W(y) */
__device__ void ckpy_d(double * rho, double * T, double * y_wk, double * P);
/*Compute rho = P*W(y)/RT */
__device__ void ckrhoy_d(double * P, double * T, double * y_wk, double * rho);
/*convert y[species] (mass fracs) to c[species] (molar conc) */
__device__ void ckytcr_d(double * rho, double * T, double * y_wk, double * c);
/*Returns the specific heats at constant volume */
/*in mass units (Eq. 29) */
__device__ void ckcvms_d(double * T, double * cvms);
/*Returns the specific heats at constant pressure */
/*in mass units (Eq. 26) */
__device__ void ckcpms_d(double * T, double * cpms);
/*Returns internal energy in mass units (Eq 30.) */
__device__ void ckums_d(double * T, double * ums);
/*Returns enthalpy in mass units (Eq 27.) */
__device__ void ckhms_d(double * T, double * hms);
/*Returns the mean specific heat at CP (Eq. 34) */
__device__ void ckcpbs_d(double * T, double * y_wk, double * cpbs);
/*Returns the mean specific heat at CV (Eq. 36) */
__device__ void ckcvbs_d(double * T, double * y_wk, double * cvbs);
/*Returns mean enthalpy of mixture in mass units */
__device__ void ckhbms_d(double * T, double * y_wk, double * hbms);
/*get mean internal energy in mass units */
__device__ void ckubms_d(double * T, double * y_wk, double * ubms);
/*compute the production rate for each species */
__device__ void ckwc_d(double * T, double * C, double * wdot, void *user_data);
/*compute the production rate for each species */
__device__ void productionRate_d(double * wdot, double * sc, double T, void *user_data);
/*compute the reaction Jacobian */
__device__ void dwdot_d(double * J, double * sc, double * Tp, int * consP, void *user_data);
__device__ void ajacobian_d(double * J, double * sc, double T, int consP, void *user_data);
__device__ void comp_k_f_d(double * tc, double invT, double * k_f, double * Corr, double * sc, void *user_data);
__device__ void comp_Kc_d(double * tc, double invT, double * Kc);
__device__ void comp_qfqr_d(double *  qf, double * qr, double * sc, double * tc, double invT, void *user_data);
/*compute the g/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
__device__ void gibbs_d(double * species, double *  tc);
/*compute Cv/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
__device__ void cv_R_d(double * species, double *  tc);
/*compute Cp/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
__device__ void cp_R_d(double * species, double *  tc);
/*compute d(Cp/R)/dT and d(Cv/R)/dT at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
__device__ void dcvpRdT_d(double * species, double * tc);
/*compute the e/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
__device__ void speciesInternalEnergy_d(double * species, double *  tc);
/*compute the h/(RT) at the given temperature (Eq 20) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
__device__ void speciesEnthalpy_d(double * species, double *  tc);
/*save molecular weights into array */
__device__ void molecularWeight_d(double * wt);
/* get temperature given internal energy in mass units and mass fracs */
__device__ void get_t_given_ey_d_(double * e, double * y_wk, double * t, int * ierr);
/* get temperature given enthalpy in mass units and mass fracs */
__device__ void get_t_given_hy_d_(double * h, double * y_wk, double * t, int * ierr);
