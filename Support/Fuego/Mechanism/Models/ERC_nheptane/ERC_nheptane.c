
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#if defined(BL_FORT_USE_UPPERCASE)
#define CKINDX CKINDX
#define CKINIT CKINIT
#define CKFINALIZE CKFINALIZE
#define CKXNUM CKXNUM
#define CKSYME CKSYME
#define CKSYMS CKSYMS
#define CKRP CKRP
#define CKPX CKPX
#define CKPY CKPY
#define CKPC CKPC
#define CKRHOX CKRHOX
#define CKRHOY CKRHOY
#define CKRHOC CKRHOC
#define CKWT CKWT
#define CKAWT CKAWT
#define CKMMWY CKMMWY
#define CKMMWX CKMMWX
#define CKMMWC CKMMWC
#define CKYTX CKYTX
#define CKYTCP CKYTCP
#define CKYTCR CKYTCR
#define CKXTY CKXTY
#define CKXTCP CKXTCP
#define CKXTCR CKXTCR
#define CKCTX CKCTX
#define CKCTY CKCTY
#define CKCPOR CKCPOR
#define CKHORT CKHORT
#define CKSOR CKSOR
#define CKCVML CKCVML
#define CKCPML CKCPML
#define CKUML CKUML
#define CKHML CKHML
#define CKGML CKGML
#define CKAML CKAML
#define CKSML CKSML
#define CKCVMS CKCVMS
#define CKCPMS CKCPMS
#define CKUMS CKUMS
#define CKHMS CKHMS
#define CKGMS CKGMS
#define CKAMS CKAMS
#define CKSMS CKSMS
#define CKCPBL CKCPBL
#define CKCPBS CKCPBS
#define CKCVBL CKCVBL
#define CKCVBS CKCVBS
#define CKHBML CKHBML
#define CKHBMS CKHBMS
#define CKUBML CKUBML
#define CKUBMS CKUBMS
#define CKSBML CKSBML
#define CKSBMS CKSBMS
#define CKGBML CKGBML
#define CKGBMS CKGBMS
#define CKABML CKABML
#define CKABMS CKABMS
#define CKWC CKWC
#define VCKWC VCKWC
#define CKWYP CKWYP
#define CKWXP CKWXP
#define CKWYR CKWYR
#define CKWXR CKWXR
#define CKQC CKQC
#define CKKFKR CKKFKR
#define CKQYP CKQYP
#define CKQXP CKQXP
#define CKQYR CKQYR
#define CKQXR CKQXR
#define CKNU CKNU
#define CKNCF CKNCF
#define CKABE CKABE
#define CKEQC CKEQC
#define CKEQYP CKEQYP
#define CKEQXP CKEQXP
#define CKEQYR CKEQYR
#define CKEQXR CKEQXR
#define DWDOT DWDOT
#define SPARSITY_INFO SPARSITY_INFO
#define SPARSITY_PREPROC SPARSITY_PREPROC
#define VGET_T_GIVEN_EY  VGET_T_GIVEN_EY
#define VCKUBMS VCKUBMS
#define VCKYTCR VCKYTCR
#define VCKUMS VCKUMS
#define VCKCVMS VCKCVMS
#define VCKHMS VCKHMS
#define VCKPY VCKPY
#define VCKWYR VCKWYR
#define VCKYTX VCKYTX
#define GET_T_GIVEN_EY GET_T_GIVEN_EY
#define GET_T_GIVEN_HY GET_T_GIVEN_HY
#define GET_REACTION_MAP GET_REACTION_MAP
#define GET_CRITPARAMS GET_CRITPARAMS
#elif defined(BL_FORT_USE_LOWERCASE)
#define CKINDX ckindx
#define CKINIT ckinit
#define CKFINALIZE ckfinalize
#define CKXNUM ckxnum
#define CKSYME cksyme
#define CKSYMS cksyms
#define CKRP ckrp
#define CKPX ckpx
#define CKPY ckpy
#define CKPC ckpc
#define CKRHOX ckrhox
#define CKRHOY ckrhoy
#define CKRHOC ckrhoc
#define CKWT ckwt
#define CKAWT ckawt
#define CKMMWY ckmmwy
#define CKMMWX ckmmwx
#define CKMMWC ckmmwc
#define CKYTX ckytx
#define CKYTCP ckytcp
#define CKYTCR ckytcr
#define CKXTY ckxty
#define CKXTCP ckxtcp
#define CKXTCR ckxtcr
#define CKCTX ckctx
#define CKCTY ckcty
#define CKCPOR ckcpor
#define CKHORT ckhort
#define CKSOR cksor
#define CKCVML ckcvml
#define CKCPML ckcpml
#define CKUML ckuml
#define CKHML ckhml
#define CKGML ckgml
#define CKAML ckaml
#define CKSML cksml
#define CKCVMS ckcvms
#define CKCPMS ckcpms
#define CKUMS ckums
#define CKHMS ckhms
#define CKGMS ckgms
#define CKAMS ckams
#define CKSMS cksms
#define CKCPBL ckcpbl
#define CKCPBS ckcpbs
#define CKCVBL ckcvbl
#define CKCVBS ckcvbs
#define CKHBML ckhbml
#define CKHBMS ckhbms
#define CKUBML ckubml
#define CKUBMS ckubms
#define CKSBML cksbml
#define CKSBMS cksbms
#define CKGBML ckgbml
#define CKGBMS ckgbms
#define CKABML ckabml
#define CKABMS ckabms
#define CKWC ckwc
#define VCKWC vckwc
#define CKWYP ckwyp
#define CKWXP ckwxp
#define CKWYR ckwyr
#define CKWXR ckwxr
#define CKQC ckqc
#define CKKFKR ckkfkr
#define CKQYP ckqyp
#define CKQXP ckqxp
#define CKQYR ckqyr
#define CKQXR ckqxr
#define CKNU cknu
#define CKNCF ckncf
#define CKABE ckabe
#define CKEQC ckeqc
#define CKEQYP ckeqyp
#define CKEQXP ckeqxp
#define CKEQYR ckeqyr
#define CKEQXR ckeqxr
#define DWDOT dwdot
#define SPARSITY_INFO sparsity_info
#define SPARSITY_PREPROC sparsity_preproc
#define VGET_T_GIVEN_EY  vget_t_given_ey
#define VCKUBMS vckubms
#define VCKYTCR vckytcr
#define VCKUMS vckums
#define VCKCVMS vckcvms
#define VCKHMS vckhms
#define VCKPY vckpy
#define VCKWYR vckwyr
#define VCKYTX vckytx
#define GET_T_GIVEN_EY get_t_given_ey
#define GET_T_GIVEN_HY get_t_given_hy
#define GET_REACTION_MAP get_reaction_map
#define GET_CRITPARAMS get_critparams
#elif defined(BL_FORT_USE_UNDERSCORE)
#define CKINDX ckindx_
#define CKINIT ckinit_
#define CKFINALIZE ckfinalize_
#define CKXNUM ckxnum_
#define CKSYME cksyme_
#define CKSYMS cksyms_
#define CKRP ckrp_
#define CKPX ckpx_
#define CKPY ckpy_
#define CKPC ckpc_
#define CKRHOX ckrhox_
#define CKRHOY ckrhoy_
#define CKRHOC ckrhoc_
#define CKWT ckwt_
#define CKAWT ckawt_
#define CKMMWY ckmmwy_
#define CKMMWX ckmmwx_
#define CKMMWC ckmmwc_
#define CKYTX ckytx_
#define CKYTCP ckytcp_
#define CKYTCR ckytcr_
#define CKXTY ckxty_
#define CKXTCP ckxtcp_
#define CKXTCR ckxtcr_
#define CKCTX ckctx_
#define CKCTY ckcty_
#define CKCPOR ckcpor_
#define CKHORT ckhort_
#define CKSOR cksor_
#define CKCVML ckcvml_
#define CKCPML ckcpml_
#define CKUML ckuml_
#define CKHML ckhml_
#define CKGML ckgml_
#define CKAML ckaml_
#define CKSML cksml_
#define CKCVMS ckcvms_
#define CKCPMS ckcpms_
#define CKUMS ckums_
#define CKHMS ckhms_
#define CKGMS ckgms_
#define CKAMS ckams_
#define CKSMS cksms_
#define CKCPBL ckcpbl_
#define CKCPBS ckcpbs_
#define CKCVBL ckcvbl_
#define CKCVBS ckcvbs_
#define CKHBML ckhbml_
#define CKHBMS ckhbms_
#define CKUBML ckubml_
#define CKUBMS ckubms_
#define CKSBML cksbml_
#define CKSBMS cksbms_
#define CKGBML ckgbml_
#define CKGBMS ckgbms_
#define CKABML ckabml_
#define CKABMS ckabms_
#define CKWC ckwc_
#define VCKWC vckwc_
#define CKWYP ckwyp_
#define CKWXP ckwxp_
#define CKWYR ckwyr_
#define CKWXR ckwxr_
#define CKQC ckqc_
#define CKKFKR ckkfkr_
#define CKQYP ckqyp_
#define CKQXP ckqxp_
#define CKQYR ckqyr_
#define CKQXR ckqxr_
#define CKNU cknu_
#define CKNCF ckncf_
#define CKABE ckabe_
#define CKEQC ckeqc_
#define CKEQYP ckeqyp_
#define CKEQXP ckeqxp_
#define CKEQYR ckeqyr_
#define CKEQXR ckeqxr_
#define DWDOT dwdot_
#define SPARSITY_INFO sparsity_info_
#define SPARSITY_PREPROC sparsity_preproc_
#define VGET_T_GIVEN_EY  vget_t_given_ey_
#define VCKUBMS vckubms_
#define VCKYTCR vckytcr_
#define VCKUMS vckums_
#define VCKCVMS vckcvms_
#define VCKHMS vckhms_
#define VCKPY vckpy_
#define VCKWYR vckwyr_
#define VCKYTX vckytx_
#define GET_T_GIVEN_EY get_t_given_ey_
#define GET_T_GIVEN_HY get_t_given_hy_
#define GET_REACTION_MAP get_reaction_map_
#define GET_CRITPARAMS get_critparams_
#endif

/*function declarations */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetEPS EGTRANSETEPS
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetEPS egtranseteps
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetEPS egtranseteps_
#endif
void egtransetEPS(double *  EPS);
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetSIG EGTRANSETSIG
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetSIG egtransetsig
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetSIG egtransetsig_
#endif
void egtransetSIG(double* SIG);
void atomicWeight(double * restrict awt);
void molecularWeight(double * restrict wt);
void gibbs(double * restrict species, double * restrict tc);
void helmholtz(double * restrict species, double * restrict tc);
void speciesInternalEnergy(double * restrict species, double * restrict tc);
void speciesEnthalpy(double * restrict species, double * restrict tc);
void speciesEntropy(double * restrict species, double * restrict tc);
void cp_R(double * restrict species, double * restrict tc);
void cv_R(double * restrict species, double * restrict tc);
void equilibriumConstants(double * restrict kc, double * restrict g_RT, double T);
void productionRate(double * restrict wdot, double * restrict sc, double T);
void comp_k_f(double * restrict tc, double invT, double * restrict k_f);
void comp_Kc(double * restrict tc, double invT, double * restrict Kc);
void comp_qfqr(double * restrict q_f, double * restrict q_r, double * restrict sc, double * restrict tc, double invT);
void progressRate(double * restrict qdot, double * restrict speciesConc, double T);
void progressRateFR(double * restrict q_f, double * restrict q_r, double * restrict speciesConc, double T);
void CKINIT();
void CKFINALIZE();
void CKINDX(int * iwrk, double * restrict rwrk, int * mm, int * kk, int * ii, int * nfit );
void CKXNUM(char * line, int * nexp, int * lout, int * nval, double * restrict rval, int * kerr, int lenline);
void CKSNUM(char * line, int * nexp, int * lout, char * kray, int * nn, int * knum, int * nval, double * restrict rval, int * kerr, int lenline, int lenkray);
void CKSYME(int * kname, int * lenkname);
void CKSYMS(int * kname, int * lenkname);
void CKRP(int * ickwrk, double * restrict rckwrk, double * restrict ru, double * restrict ruc, double * restrict pa);
void CKPX(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict P);
void CKPY(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict P);
void CKPC(double * restrict rho, double * restrict T, double * restrict c, int * iwrk, double * restrict rwrk, double * restrict P);
void CKRHOX(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict rho);
void CKRHOY(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict rho);
void CKRHOC(double * restrict P, double * restrict T, double * restrict c, int * iwrk, double * restrict rwrk, double * restrict rho);
void CKWT(int * iwrk, double * restrict rwrk, double * restrict wt);
void CKAWT(int * iwrk, double * restrict rwrk, double * restrict awt);
void CKMMWY(double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wtm);
void CKMMWX(double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wtm);
void CKMMWC(double * restrict c, int * iwrk, double * restrict rwrk, double * restrict wtm);
void CKYTX(double * restrict y, int * iwrk, double * restrict rwrk, double * restrict x);
void CKYTCP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict c);
void CKYTCR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict c);
void CKXTY(double * restrict x, int * iwrk, double * restrict rwrk, double * restrict y);
void CKXTCP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict c);
void CKXTCR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict c);
void CKCTX(double * restrict c, int * iwrk, double * restrict rwrk, double * restrict x);
void CKCTY(double * restrict c, int * iwrk, double * restrict rwrk, double * restrict y);
void CKCPOR(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cpor);
void CKHORT(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict hort);
void CKSOR(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict sor);
void CKCVML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cvml);
void CKCPML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cvml);
void CKUML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict uml);
void CKHML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict uml);
void CKGML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict gml);
void CKAML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict aml);
void CKSML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict sml);
void CKCVMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cvms);
void CKCPMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cvms);
void CKUMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict ums);
void CKHMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict ums);
void CKGMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict gms);
void CKAMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict ams);
void CKSMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict sms);
void CKCPBL(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict cpbl);
void CKCPBS(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict cpbs);
void CKCVBL(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict cpbl);
void CKCVBS(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict cpbs);
void CKHBML(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict hbml);
void CKHBMS(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict hbms);
void CKUBML(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict ubml);
void CKUBMS(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict ubms);
void CKSBML(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict sbml);
void CKSBMS(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict sbms);
void CKGBML(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict gbml);
void CKGBMS(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict gbms);
void CKABML(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict abml);
void CKABMS(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict abms);
void CKWC(double * restrict T, double * restrict C, int * iwrk, double * restrict rwrk, double * restrict wdot);
void CKWYP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wdot);
void CKWXP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wdot);
void CKWYR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wdot);
void CKWXR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wdot);
void CKQC(double * restrict T, double * restrict C, int * iwrk, double * restrict rwrk, double * restrict qdot);
void CKKFKR(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict q_f, double * restrict q_r);
void CKQYP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict qdot);
void CKQXP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict qdot);
void CKQYR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict qdot);
void CKQXR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict qdot);
void CKNU(int * kdim, int * iwrk, double * restrict rwrk, int * nuki);
void CKNCF(int * mdim, int * iwrk, double * restrict rwrk, int * ncf);
void CKABE(int * iwrk, double * restrict rwrk, double * restrict a, double * restrict b, double * restrict e );
void CKEQC(double * restrict T, double * restrict C , int * iwrk, double * restrict rwrk, double * restrict eqcon );
void CKEQYP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict eqcon);
void CKEQXP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict eqcon);
void CKEQYR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict eqcon);
void CKEQXR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict eqcon);
void DWDOT(double * restrict J, double * restrict sc, double * restrict T, int * consP);
void aJacobian(double * restrict J, double * restrict sc, double T, int consP);
void SPARSITY_INFO(int * nJdata);
void SPARSITY_PREPROC(int * restrict rowVals, int * restrict colPtrs);
void dcvpRdT(double * restrict species, double * restrict tc);
void GET_T_GIVEN_EY(double * restrict e, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict t, int *ierr);
void GET_T_GIVEN_HY(double * restrict h, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict t, int *ierr);
void GET_REACTION_MAP(int * restrict rmap);
/*vector version */
// ANNE
void VCKWC(int * restrict np, double * restrict T, double * restrict C, int * iwrk, double * restrict rwrk, double * restrict wdot);
void vcv_R(int npt, double * restrict species, double * restrict T_in);
//void CKCVMS(                 double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cvms);
void VCKCVMS(int * restrict np, double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cvms);
//void CKUMS(                 double * restrict T, int * iwrk, double * restrict rwrk, double * restrict ums);
void VCKUMS(int * restrict np, double * restrict T, int * iwrk, double * restrict rwrk, double * restrict ums);
//void CKYTCR(d                 double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict c);
void VCKYTCR(int * restrict np, double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict c);
//void CKUBMS(                  double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict ubms);
void VCKUBMS(int * restrict np, double * restrict T, double * restrict y, double * restrict ubms);
void vspeciesInternalEnergy(int npt, double * restrict species, double * restrict tc);
void VGET_T_GIVEN_EY(int * restrict np, double * restrict e, double * restrict y, double * restrict T, int * ierr);
void vproductionRate(int npt, double * restrict wdot, double * restrict c, double * restrict T);
void VCKHMS(int * restrict np, double * restrict T, int * iwrk, double * restrict rwrk, double * restrict ums);
//void CKPY(                  double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict P);
void VCKPY(int * restrict np, double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict P);
void VCKWYR(int * restrict np, double * restrict rho, double * restrict T,
            double * restrict y, int * restrict iwrk, double * restrict rwrk,
            double * restrict wdot);
void VCKYTX(int * restrict np, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict x);
void vcomp_k_f(int npt, double * restrict k_f_s, double * restrict tc, double * restrict invT);
void vcomp_gibbs(int npt, double * restrict g_RT, double * restrict tc);
void vcomp_Kc(int npt, double * restrict Kc_s, double * restrict g_RT, double * restrict invT);
void GET_CRITPARAMS(double * restrict Tci, double * restrict ai, double * restrict bi, double * restrict acentric_i);
void vcomp_wdot_1_50(int npt, double * restrict wdot, double * restrict mixture, double * restrict sc,
                double * restrict k_f_s, double * restrict Kc_s,
                double * restrict tc, double * restrict invT, double * restrict T);
void vcomp_wdot_51_52(int npt, double * restrict wdot, double * restrict mixture, double * restrict sc,
                double * restrict k_f_s, double * restrict Kc_s,
                double * restrict tc, double * restrict invT, double * restrict T);

/* Inverse molecular weights */
static const double imw[29] = {
    1.0 / 100.205570,  /*NC7H16 */
    1.0 / 31.998800,  /*O2 */
    1.0 / 44.009950,  /*CO2 */
    1.0 / 18.015340,  /*H2O */
    1.0 / 28.010550,  /*CO */
    1.0 / 2.015940,  /*H2 */
    1.0 / 17.007370,  /*OH */
    1.0 / 34.014740,  /*H2O2 */
    1.0 / 33.006770,  /*HO2 */
    1.0 / 1.007970,  /*H */
    1.0 / 15.999400,  /*O */
    1.0 / 31.034460,  /*CH3O */
    1.0 / 30.026490,  /*CH2O */
    1.0 / 29.018520,  /*HCO */
    1.0 / 14.027090,  /*CH2 */
    1.0 / 15.035060,  /*CH3 */
    1.0 / 16.043030,  /*CH4 */
    1.0 / 27.046210,  /*C2H3 */
    1.0 / 28.054180,  /*C2H4 */
    1.0 / 29.062150,  /*C2H5 */
    1.0 / 40.065330,  /*C3H4 */
    1.0 / 41.073300,  /*C3H5 */
    1.0 / 42.081270,  /*C3H6 */
    1.0 / 43.089240,  /*C3H7 */
    1.0 / 99.197600,  /*C7H15-2 */
    1.0 / 131.196400,  /*C7H15O2 */
    1.0 / 146.187830,  /*C7KET12 */
    1.0 / 99.153970,  /*C5H11CO */
    1.0 / 28.013400};  /*N2 */



static double fwd_A[52], fwd_beta[52], fwd_Ea[52];
static double low_A[52], low_beta[52], low_Ea[52];
static double rev_A[52], rev_beta[52], rev_Ea[52];
static double troe_a[52],troe_Ts[52], troe_Tss[52], troe_Tsss[52];
static double sri_a[52], sri_b[52], sri_c[52], sri_d[52], sri_e[52];
static double activation_units[52], prefactor_units[52], phase_units[52];
static int is_PD[52], troe_len[52], sri_len[52], nTB[52], *TBid[52];
static double *TB[52];

static double fwd_A_DEF[52], fwd_beta_DEF[52], fwd_Ea_DEF[52];
static double low_A_DEF[52], low_beta_DEF[52], low_Ea_DEF[52];
static double rev_A_DEF[52], rev_beta_DEF[52], rev_Ea_DEF[52];
static double troe_a_DEF[52],troe_Ts_DEF[52], troe_Tss_DEF[52], troe_Tsss_DEF[52];
static double sri_a_DEF[52], sri_b_DEF[52], sri_c_DEF[52], sri_d_DEF[52], sri_e_DEF[52];
static double activation_units_DEF[52], prefactor_units_DEF[52], phase_units_DEF[52];
static int is_PD_DEF[52], troe_len_DEF[52], sri_len_DEF[52], nTB_DEF[52], *TBid_DEF[52];
static double *TB_DEF[52];
static int rxn_map[52] = {3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,0,1,32,33,34,35,36,2,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51};

void GET_REACTION_MAP(int *rmap)
{
    for (int i=0; i<52; ++i) {
        rmap[i] = rxn_map[i];
    }
}


#include <ReactionData.H>
double* GetParamPtr(int                reaction_id,
                    REACTION_PARAMETER param_id,
                    int                species_id,
                    int                get_default)
{
  double* ret = 0;
  if (reaction_id<0 || reaction_id>=52) {
    printf("Bad reaction id = %d",reaction_id);
    abort();
  };
  int mrid = rxn_map[reaction_id];

  if (param_id == THIRD_BODY) {
    if (species_id<0 || species_id>=29) {
      printf("GetParamPtr: Bad species id = %d",species_id);
      abort();
    }
    if (get_default) {
      for (int i=0; i<nTB_DEF[mrid]; ++i) {
        if (species_id == TBid_DEF[mrid][i]) {
          ret = &(TB_DEF[mrid][i]);
        }
      }
    }
    else {
      for (int i=0; i<nTB[mrid]; ++i) {
        if (species_id == TBid[mrid][i]) {
          ret = &(TB[mrid][i]);
        }
      }
    }
    if (ret == 0) {
      printf("GetParamPtr: No TB for reaction id = %d",reaction_id);
      abort();
    }
  }
  else {
    if (     param_id == FWD_A)     {ret = (get_default ? &(fwd_A_DEF[mrid]) : &(fwd_A[mrid]));}
      else if (param_id == FWD_BETA)  {ret = (get_default ? &(fwd_beta_DEF[mrid]) : &(fwd_beta[mrid]));}
      else if (param_id == FWD_EA)    {ret = (get_default ? &(fwd_Ea_DEF[mrid]) : &(fwd_Ea[mrid]));}
      else if (param_id == LOW_A)     {ret = (get_default ? &(low_A_DEF[mrid]) : &(low_A[mrid]));}
      else if (param_id == LOW_BETA)  {ret = (get_default ? &(low_beta_DEF[mrid]) : &(low_beta[mrid]));}
      else if (param_id == LOW_EA)    {ret = (get_default ? &(low_Ea_DEF[mrid]) : &(low_Ea[mrid]));}
      else if (param_id == REV_A)     {ret = (get_default ? &(rev_A_DEF[mrid]) : &(rev_A[mrid]));}
      else if (param_id == REV_BETA)  {ret = (get_default ? &(rev_beta_DEF[mrid]) : &(rev_beta[mrid]));}
      else if (param_id == REV_EA)    {ret = (get_default ? &(rev_Ea_DEF[mrid]) : &(rev_Ea[mrid]));}
      else if (param_id == TROE_A)    {ret = (get_default ? &(troe_a_DEF[mrid]) : &(troe_a[mrid]));}
      else if (param_id == TROE_TS)   {ret = (get_default ? &(troe_Ts_DEF[mrid]) : &(troe_Ts[mrid]));}
      else if (param_id == TROE_TSS)  {ret = (get_default ? &(troe_Tss_DEF[mrid]) : &(troe_Tss[mrid]));}
      else if (param_id == TROE_TSSS) {ret = (get_default ? &(troe_Tsss_DEF[mrid]) : &(troe_Tsss[mrid]));}
      else if (param_id == SRI_A)     {ret = (get_default ? &(sri_a_DEF[mrid]) : &(sri_a[mrid]));}
      else if (param_id == SRI_B)     {ret = (get_default ? &(sri_b_DEF[mrid]) : &(sri_b[mrid]));}
      else if (param_id == SRI_C)     {ret = (get_default ? &(sri_c_DEF[mrid]) : &(sri_c[mrid]));}
      else if (param_id == SRI_D)     {ret = (get_default ? &(sri_d_DEF[mrid]) : &(sri_d[mrid]));}
      else if (param_id == SRI_E)     {ret = (get_default ? &(sri_e_DEF[mrid]) : &(sri_e[mrid]));}
    else {
      printf("GetParamPtr: Unknown parameter id");
      abort();
    }
  }
  return ret;
}

void ResetAllParametersToDefault()
{
    for (int i=0; i<52; i++) {
        if (nTB[i] != 0) {
            nTB[i] = 0;
            free(TB[i]);
            free(TBid[i]);
        }

        fwd_A[i]    = fwd_A_DEF[i];
        fwd_beta[i] = fwd_beta_DEF[i];
        fwd_Ea[i]   = fwd_Ea_DEF[i];

        low_A[i]    = low_A_DEF[i];
        low_beta[i] = low_beta_DEF[i];
        low_Ea[i]   = low_Ea_DEF[i];

        rev_A[i]    = rev_A_DEF[i];
        rev_beta[i] = rev_beta_DEF[i];
        rev_Ea[i]   = rev_Ea_DEF[i];

        troe_a[i]    = troe_a_DEF[i];
        troe_Ts[i]   = troe_Ts_DEF[i];
        troe_Tss[i]  = troe_Tss_DEF[i];
        troe_Tsss[i] = troe_Tsss_DEF[i];

        sri_a[i] = sri_a_DEF[i];
        sri_b[i] = sri_b_DEF[i];
        sri_c[i] = sri_c_DEF[i];
        sri_d[i] = sri_d_DEF[i];
        sri_e[i] = sri_e_DEF[i];

        is_PD[i]    = is_PD_DEF[i];
        troe_len[i] = troe_len_DEF[i];
        sri_len[i]  = sri_len_DEF[i];

        activation_units[i] = activation_units_DEF[i];
        prefactor_units[i]  = prefactor_units_DEF[i];
        phase_units[i]      = phase_units_DEF[i];

        nTB[i]  = nTB_DEF[i];
        if (nTB[i] != 0) {
           TB[i] = (double *) malloc(sizeof(double) * nTB[i]);
           TBid[i] = (int *) malloc(sizeof(int) * nTB[i]);
           for (int j=0; j<nTB[i]; j++) {
             TB[i][j] = TB_DEF[i][j];
             TBid[i][j] = TBid_DEF[i][j];
           }
        }
    }
}

void SetAllDefaults()
{
    for (int i=0; i<52; i++) {
        if (nTB_DEF[i] != 0) {
            nTB_DEF[i] = 0;
            free(TB_DEF[i]);
            free(TBid_DEF[i]);
        }

        fwd_A_DEF[i]    = fwd_A[i];
        fwd_beta_DEF[i] = fwd_beta[i];
        fwd_Ea_DEF[i]   = fwd_Ea[i];

        low_A_DEF[i]    = low_A[i];
        low_beta_DEF[i] = low_beta[i];
        low_Ea_DEF[i]   = low_Ea[i];

        rev_A_DEF[i]    = rev_A[i];
        rev_beta_DEF[i] = rev_beta[i];
        rev_Ea_DEF[i]   = rev_Ea[i];

        troe_a_DEF[i]    = troe_a[i];
        troe_Ts_DEF[i]   = troe_Ts[i];
        troe_Tss_DEF[i]  = troe_Tss[i];
        troe_Tsss_DEF[i] = troe_Tsss[i];

        sri_a_DEF[i] = sri_a[i];
        sri_b_DEF[i] = sri_b[i];
        sri_c_DEF[i] = sri_c[i];
        sri_d_DEF[i] = sri_d[i];
        sri_e_DEF[i] = sri_e[i];

        is_PD_DEF[i]    = is_PD[i];
        troe_len_DEF[i] = troe_len[i];
        sri_len_DEF[i]  = sri_len[i];

        activation_units_DEF[i] = activation_units[i];
        prefactor_units_DEF[i]  = prefactor_units[i];
        phase_units_DEF[i]      = phase_units[i];

        nTB_DEF[i]  = nTB[i];
        if (nTB_DEF[i] != 0) {
           TB_DEF[i] = (double *) malloc(sizeof(double) * nTB_DEF[i]);
           TBid_DEF[i] = (int *) malloc(sizeof(int) * nTB_DEF[i]);
           for (int j=0; j<nTB_DEF[i]; j++) {
             TB_DEF[i][j] = TB[i][j];
             TBid_DEF[i][j] = TBid[i][j];
           }
        }
    }
}

/* Finalizes parameter database */
void CKFINALIZE()
{
  for (int i=0; i<52; ++i) {
    free(TB[i]); TB[i] = 0; 
    free(TBid[i]); TBid[i] = 0;
    nTB[i] = 0;

    free(TB_DEF[i]); TB_DEF[i] = 0; 
    free(TBid_DEF[i]); TBid_DEF[i] = 0;
    nTB_DEF[i] = 0;
  }
}

/* Initializes parameter database */
void CKINIT()
{
    // (0):  NC7H16 + H <=> H2 + C7H15-2
    fwd_A[3]     = 43800000;
    fwd_beta[3]  = 2;
    fwd_Ea[3]    = 4759.4840000000004;
    prefactor_units[3]  = 1.0000000000000002e-06;
    activation_units[3] = 0.50321666580471969;
    phase_units[3]      = 1e-12;
    is_PD[3] = 0;
    nTB[3] = 0;

    // (1):  NC7H16 + OH <=> H2O + C7H15-2
    fwd_A[4]     = 9700000000;
    fwd_beta[4]  = 1.3;
    fwd_Ea[4]    = 1689.817;
    prefactor_units[4]  = 1.0000000000000002e-06;
    activation_units[4] = 0.50321666580471969;
    phase_units[4]      = 1e-12;
    is_PD[4] = 0;
    nTB[4] = 0;

    // (2):  NC7H16 + HO2 <=> H2O2 + C7H15-2
    fwd_A[5]     = 16500000000000;
    fwd_beta[5]  = 0;
    fwd_Ea[5]    = 16948.16;
    prefactor_units[5]  = 1.0000000000000002e-06;
    activation_units[5] = 0.50321666580471969;
    phase_units[5]      = 1e-12;
    is_PD[5] = 0;
    nTB[5] = 0;

    // (3):  NC7H16 + O2 <=> HO2 + C7H15-2
    fwd_A[6]     = 2000000000000000;
    fwd_beta[6]  = 0;
    fwd_Ea[6]    = 47374.860000000001;
    prefactor_units[6]  = 1.0000000000000002e-06;
    activation_units[6] = 0.50321666580471969;
    phase_units[6]      = 1e-12;
    is_PD[6] = 0;
    nTB[6] = 0;

    // (4):  O2 + C7H15-2 <=> C7H15O2
    fwd_A[7]     = 1560000000000;
    fwd_beta[7]  = 0;
    fwd_Ea[7]    = 0;
    prefactor_units[7]  = 1.0000000000000002e-06;
    activation_units[7] = 0.50321666580471969;
    phase_units[7]      = 1e-12;
    is_PD[7] = 0;
    nTB[7] = 0;

    // (5):  O2 + C7H15O2 <=> OH + C7KET12
    fwd_A[8]     = 450000000000000;
    fwd_beta[8]  = 0;
    fwd_Ea[8]    = 18230.73;
    prefactor_units[8]  = 1.0000000000000002e-06;
    activation_units[8] = 0.50321666580471969;
    phase_units[8]      = 1e-12;
    is_PD[8] = 0;
    nTB[8] = 0;

    // (6):  C7KET12 <=> OH + CH2O + C5H11CO
    fwd_A[9]     = 953000000000000;
    fwd_beta[9]  = 0;
    fwd_Ea[9]    = 41095.540000000001;
    prefactor_units[9]  = 1;
    activation_units[9] = 0.50321666580471969;
    phase_units[9]      = 1e-6;
    is_PD[9] = 0;
    nTB[9] = 0;

    // (7):  C5H11CO <=> CO + C2H4 + C3H7
    fwd_A[10]     = 9840000000000000;
    fwd_beta[10]  = 0;
    fwd_Ea[10]    = 40195.639999999999;
    prefactor_units[10]  = 1;
    activation_units[10] = 0.50321666580471969;
    phase_units[10]      = 1e-6;
    is_PD[10] = 0;
    nTB[10] = 0;

    // (8):  C7H15-2 <=> C2H4 + C2H5 + C3H6
    fwd_A[11]     = 704500000000000;
    fwd_beta[11]  = 0;
    fwd_Ea[11]    = 34596.25;
    prefactor_units[11]  = 1;
    activation_units[11] = 0.50321666580471969;
    phase_units[11]      = 1e-6;
    is_PD[11] = 0;
    nTB[11] = 0;

    // (9):  C3H7 <=> CH3 + C2H4
    fwd_A[12]     = 96000000000000;
    fwd_beta[12]  = 0;
    fwd_Ea[12]    = 30946.639999999999;
    prefactor_units[12]  = 1;
    activation_units[12] = 0.50321666580471969;
    phase_units[12]      = 1e-6;
    is_PD[12] = 0;
    nTB[12] = 0;

    // (10):  C3H7 <=> H + C3H6
    fwd_A[13]     = 125000000000000;
    fwd_beta[13]  = 0;
    fwd_Ea[13]    = 36896;
    prefactor_units[13]  = 1;
    activation_units[13] = 0.50321666580471969;
    phase_units[13]      = 1e-6;
    is_PD[13] = 0;
    nTB[13] = 0;

    // (11):  CH3 + C3H6 <=> CH4 + C3H5
    fwd_A[14]     = 9000000000000;
    fwd_beta[14]  = 0;
    fwd_Ea[14]    = 8479.0799999999999;
    prefactor_units[14]  = 1.0000000000000002e-06;
    activation_units[14] = 0.50321666580471969;
    phase_units[14]      = 1e-12;
    is_PD[14] = 0;
    nTB[14] = 0;

    // (12):  O2 + C3H5 <=> HO2 + C3H4
    fwd_A[15]     = 600000000000;
    fwd_beta[15]  = 0;
    fwd_Ea[15]    = 9998.9150000000009;
    prefactor_units[15]  = 1.0000000000000002e-06;
    activation_units[15] = 0.50321666580471969;
    phase_units[15]      = 1e-12;
    is_PD[15] = 0;
    nTB[15] = 0;

    // (13):  OH + C3H4 <=> CH2O + C2H3
    fwd_A[16]     = 1000000000000;
    fwd_beta[16]  = 0;
    fwd_Ea[16]    = 0;
    prefactor_units[16]  = 1.0000000000000002e-06;
    activation_units[16] = 0.50321666580471969;
    phase_units[16]      = 1e-12;
    is_PD[16] = 0;
    nTB[16] = 0;

    // (14):  OH + C3H4 <=> HCO + C2H4
    fwd_A[17]     = 1000000000000;
    fwd_beta[17]  = 0;
    fwd_Ea[17]    = 0;
    prefactor_units[17]  = 1.0000000000000002e-06;
    activation_units[17] = 0.50321666580471969;
    phase_units[17]      = 1e-12;
    is_PD[17] = 0;
    nTB[17] = 0;

    // (15):  HO2 + CH3 <=> OH + CH3O
    fwd_A[18]     = 50000000000000;
    fwd_beta[18]  = 0;
    fwd_Ea[18]    = 0;
    prefactor_units[18]  = 1.0000000000000002e-06;
    activation_units[18] = 0.50321666580471969;
    phase_units[18]      = 1e-12;
    is_PD[18] = 0;
    nTB[18] = 0;

    // (16):  OH + CH3 <=> H2O + CH2
    fwd_A[19]     = 7500000;
    fwd_beta[19]  = 2;
    fwd_Ea[19]    = 4999.4579999999996;
    prefactor_units[19]  = 1.0000000000000002e-06;
    activation_units[19] = 0.50321666580471969;
    phase_units[19]      = 1e-12;
    is_PD[19] = 0;
    nTB[19] = 0;

    // (17):  OH + CH2 <=> H + CH2O
    fwd_A[20]     = 25000000000000;
    fwd_beta[20]  = 0;
    fwd_Ea[20]    = 0;
    prefactor_units[20]  = 1.0000000000000002e-06;
    activation_units[20] = 0.50321666580471969;
    phase_units[20]      = 1e-12;
    is_PD[20] = 0;
    nTB[20] = 0;

    // (18):  O2 + CH2 <=> OH + HCO
    fwd_A[21]     = 43000000000;
    fwd_beta[21]  = 0;
    fwd_Ea[21]    = -499.94580000000002;
    prefactor_units[21]  = 1.0000000000000002e-06;
    activation_units[21] = 0.50321666580471969;
    phase_units[21]      = 1e-12;
    is_PD[21] = 0;
    nTB[21] = 0;

    // (19):  O2 + CH2 <=> CO2 + H2
    fwd_A[22]     = 690000000000;
    fwd_beta[22]  = 0;
    fwd_Ea[22]    = 499.94580000000002;
    prefactor_units[22]  = 1.0000000000000002e-06;
    activation_units[22] = 0.50321666580471969;
    phase_units[22]      = 1e-12;
    is_PD[22] = 0;
    nTB[22] = 0;

    // (20):  O2 + CH2 <=> H2O + CO
    fwd_A[23]     = 20000000000;
    fwd_beta[23]  = 0;
    fwd_Ea[23]    = -999.89149999999995;
    prefactor_units[23]  = 1.0000000000000002e-06;
    activation_units[23] = 0.50321666580471969;
    phase_units[23]      = 1e-12;
    is_PD[23] = 0;
    nTB[23] = 0;

    // (21):  O2 + CH2 <=> O + CH2O
    fwd_A[24]     = 50000000000000;
    fwd_beta[24]  = 0;
    fwd_Ea[24]    = 8999.0239999999994;
    prefactor_units[24]  = 1.0000000000000002e-06;
    activation_units[24] = 0.50321666580471969;
    phase_units[24]      = 1e-12;
    is_PD[24] = 0;
    nTB[24] = 0;

    // (22):  O2 + CH2 <=> CO2 + 2 H
    fwd_A[25]     = 1600000000000;
    fwd_beta[25]  = 0;
    fwd_Ea[25]    = 999.89149999999995;
    prefactor_units[25]  = 1.0000000000000002e-06;
    activation_units[25] = 0.50321666580471969;
    phase_units[25]      = 1e-12;
    is_PD[25] = 0;
    nTB[25] = 0;

    // (23):  O2 + CH2 <=> CO + OH + H
    fwd_A[26]     = 86000000000;
    fwd_beta[26]  = 0;
    fwd_Ea[26]    = -499.94580000000002;
    prefactor_units[26]  = 1.0000000000000002e-06;
    activation_units[26] = 0.50321666580471969;
    phase_units[26]      = 1e-12;
    is_PD[26] = 0;
    nTB[26] = 0;

    // (24):  CO + CH3O <=> CO2 + CH3
    fwd_A[27]     = 157000000000000;
    fwd_beta[27]  = 0;
    fwd_Ea[27]    = 11798.719999999999;
    prefactor_units[27]  = 1.0000000000000002e-06;
    activation_units[27] = 0.50321666580471969;
    phase_units[27]      = 1e-12;
    is_PD[27] = 0;
    nTB[27] = 0;

    // (25):  CO + OH <=> CO2 + H
    fwd_A[28]     = 89870000;
    fwd_beta[28]  = 1.3799999999999999;
    fwd_Ea[28]    = 5232.3090000000002;
    prefactor_units[28]  = 1.0000000000000002e-06;
    activation_units[28] = 0.50321666580471969;
    phase_units[28]      = 1e-12;
    is_PD[28] = 0;
    nTB[28] = 0;

    // (26):  OH + O <=> O2 + H
    fwd_A[29]     = 400000000000000;
    fwd_beta[29]  = -0.5;
    fwd_Ea[29]    = 0;
    prefactor_units[29]  = 1.0000000000000002e-06;
    activation_units[29] = 0.50321666580471969;
    phase_units[29]      = 1e-12;
    is_PD[29] = 0;
    nTB[29] = 0;

    // (27):  HO2 + H <=> 2 OH
    fwd_A[30]     = 170000000000000;
    fwd_beta[30]  = 0;
    fwd_Ea[30]    = 874.90509999999995;
    prefactor_units[30]  = 1.0000000000000002e-06;
    activation_units[30] = 0.50321666580471969;
    phase_units[30]      = 1e-12;
    is_PD[30] = 0;
    nTB[30] = 0;

    // (28):  2 OH <=> H2O + O
    fwd_A[31]     = 600000000;
    fwd_beta[31]  = 1.3;
    fwd_Ea[31]    = 0;
    prefactor_units[31]  = 1.0000000000000002e-06;
    activation_units[31] = 0.50321666580471969;
    phase_units[31]      = 1e-12;
    is_PD[31] = 0;
    nTB[31] = 0;

    // (29):  O2 + H + M <=> HO2 + M
    fwd_A[0]     = 3.6e+17;
    fwd_beta[0]  = -0.71999999999999997;
    fwd_Ea[0]    = 0;
    prefactor_units[0]  = 1.0000000000000002e-12;
    activation_units[0] = 0.50321666580471969;
    phase_units[0]      = 1e-12;
    is_PD[0] = 0;
    nTB[0] = 4;
    TB[0] = (double *) malloc(4 * sizeof(double));
    TBid[0] = (int *) malloc(4 * sizeof(int));
    TBid[0][0] = 2; TB[0][0] = 5; // CO2
    TBid[0][1] = 3; TB[0][1] = 21; // H2O
    TBid[0][2] = 4; TB[0][2] = 2; // CO
    TBid[0][3] = 5; TB[0][3] = 3.2999999999999998; // H2

    // (30):  H2O2 + M <=> 2 OH + M
    fwd_A[1]     = 10000000000000000;
    fwd_beta[1]  = 0;
    fwd_Ea[1]    = 45495.059999999998;
    prefactor_units[1]  = 1.0000000000000002e-06;
    activation_units[1] = 0.50321666580471969;
    phase_units[1]      = 1e-6;
    is_PD[1] = 0;
    nTB[1] = 4;
    TB[1] = (double *) malloc(4 * sizeof(double));
    TBid[1] = (int *) malloc(4 * sizeof(int));
    TBid[1][0] = 2; TB[1][0] = 5; // CO2
    TBid[1][1] = 3; TB[1][1] = 21; // H2O
    TBid[1][2] = 4; TB[1][2] = 2; // CO
    TBid[1][3] = 5; TB[1][3] = 3.2999999999999998; // H2

    // (31):  H2 + OH <=> H2O + H
    fwd_A[32]     = 1170000000;
    fwd_beta[32]  = 1.3;
    fwd_Ea[32]    = 3625.607;
    prefactor_units[32]  = 1.0000000000000002e-06;
    activation_units[32] = 0.50321666580471969;
    phase_units[32]      = 1e-12;
    is_PD[32] = 0;
    nTB[32] = 0;

    // (32):  2 HO2 <=> O2 + H2O2
    fwd_A[33]     = 3000000000000;
    fwd_beta[33]  = 0;
    fwd_Ea[33]    = 0;
    prefactor_units[33]  = 1.0000000000000002e-06;
    activation_units[33] = 0.50321666580471969;
    phase_units[33]      = 1e-12;
    is_PD[33] = 0;
    nTB[33] = 0;

    // (33):  OH + CH2O <=> H2O + HCO
    fwd_A[34]     = 55630000000;
    fwd_beta[34]  = 1.0900000000000001;
    fwd_Ea[34]    = -76.508700000000005;
    prefactor_units[34]  = 1.0000000000000002e-06;
    activation_units[34] = 0.50321666580471969;
    phase_units[34]      = 1e-12;
    is_PD[34] = 0;
    nTB[34] = 0;

    // (34):  HO2 + CH2O <=> H2O2 + HCO
    fwd_A[35]     = 3000000000000;
    fwd_beta[35]  = 0;
    fwd_Ea[35]    = 7999.1319999999996;
    prefactor_units[35]  = 1.0000000000000002e-06;
    activation_units[35] = 0.50321666580471969;
    phase_units[35]      = 1e-12;
    is_PD[35] = 0;
    nTB[35] = 0;

    // (35):  O2 + HCO <=> CO + HO2
    fwd_A[36]     = 33000000000000;
    fwd_beta[36]  = -0.40000000000000002;
    fwd_Ea[36]    = 0;
    prefactor_units[36]  = 1.0000000000000002e-06;
    activation_units[36] = 0.50321666580471969;
    phase_units[36]      = 1e-12;
    is_PD[36] = 0;
    nTB[36] = 0;

    // (36):  HCO + M <=> CO + H + M
    fwd_A[2]     = 1.591e+18;
    fwd_beta[2]  = 0.94999999999999996;
    fwd_Ea[2]    = 56706.18;
    prefactor_units[2]  = 1.0000000000000002e-06;
    activation_units[2] = 0.50321666580471969;
    phase_units[2]      = 1e-6;
    is_PD[2] = 0;
    nTB[2] = 0;
    //TB[2] = (double *) malloc(0 * sizeof(double));
    //TBid[2] = (int *) malloc(0 * sizeof(int));

    // (37):  CH3O + CH3 <=> CH2O + CH4
    fwd_A[37]     = 430000000000000;
    fwd_beta[37]  = 0;
    fwd_Ea[37]    = 0;
    prefactor_units[37]  = 1.0000000000000002e-06;
    activation_units[37] = 0.50321666580471969;
    phase_units[37]      = 1e-12;
    is_PD[37] = 0;
    nTB[37] = 0;

    // (38):  OH + C2H4 <=> CH2O + CH3
    fwd_A[38]     = 60000000000000;
    fwd_beta[38]  = 0;
    fwd_Ea[38]    = 959.89589999999998;
    prefactor_units[38]  = 1.0000000000000002e-06;
    activation_units[38] = 0.50321666580471969;
    phase_units[38]      = 1e-12;
    is_PD[38] = 0;
    nTB[38] = 0;

    // (39):  OH + C2H4 <=> H2O + C2H3
    fwd_A[39]     = 80200000000000;
    fwd_beta[39]  = 0;
    fwd_Ea[39]    = 5954.3540000000003;
    prefactor_units[39]  = 1.0000000000000002e-06;
    activation_units[39] = 0.50321666580471969;
    phase_units[39]      = 1e-12;
    is_PD[39] = 0;
    nTB[39] = 0;

    // (40):  O2 + C2H3 <=> CH2O + HCO
    fwd_A[40]     = 4000000000000;
    fwd_beta[40]  = 0;
    fwd_Ea[40]    = -249.97290000000001;
    prefactor_units[40]  = 1.0000000000000002e-06;
    activation_units[40] = 0.50321666580471969;
    phase_units[40]      = 1e-12;
    is_PD[40] = 0;
    nTB[40] = 0;

    // (41):  HCO + C2H3 <=> CO + C2H4
    fwd_A[41]     = 60340000000000;
    fwd_beta[41]  = 0;
    fwd_Ea[41]    = 0;
    prefactor_units[41]  = 1.0000000000000002e-06;
    activation_units[41] = 0.50321666580471969;
    phase_units[41]      = 1e-12;
    is_PD[41] = 0;
    nTB[41] = 0;

    // (42):  O2 + C2H5 <=> HO2 + C2H4
    fwd_A[42]     = 20000000000;
    fwd_beta[42]  = 0;
    fwd_Ea[42]    = -2199.761;
    prefactor_units[42]  = 1.0000000000000002e-06;
    activation_units[42] = 0.50321666580471969;
    phase_units[42]      = 1e-12;
    is_PD[42] = 0;
    nTB[42] = 0;

    // (43):  O2 + CH4 <=> HO2 + CH3
    fwd_A[43]     = 79000000000000;
    fwd_beta[43]  = 0;
    fwd_Ea[43]    = 55993.93;
    prefactor_units[43]  = 1.0000000000000002e-06;
    activation_units[43] = 0.50321666580471969;
    phase_units[43]      = 1e-12;
    is_PD[43] = 0;
    nTB[43] = 0;

    // (44):  OH + HO2 <=> O2 + H2O
    fwd_A[44]     = 7500000000000;
    fwd_beta[44]  = 0;
    fwd_Ea[44]    = 0;
    prefactor_units[44]  = 1.0000000000000002e-06;
    activation_units[44] = 0.50321666580471969;
    phase_units[44]      = 1e-12;
    is_PD[44] = 0;
    nTB[44] = 0;

    // (45):  O2 + CH3 <=> OH + CH2O
    fwd_A[45]     = 380000000000;
    fwd_beta[45]  = 0;
    fwd_Ea[45]    = 8999.0239999999994;
    prefactor_units[45]  = 1.0000000000000002e-06;
    activation_units[45] = 0.50321666580471969;
    phase_units[45]      = 1e-12;
    is_PD[45] = 0;
    nTB[45] = 0;

    // (46):  H + CH4 <=> H2 + CH3
    fwd_A[46]     = 660000000;
    fwd_beta[46]  = 1.6000000000000001;
    fwd_Ea[46]    = 10838.82;
    prefactor_units[46]  = 1.0000000000000002e-06;
    activation_units[46] = 0.50321666580471969;
    phase_units[46]      = 1e-12;
    is_PD[46] = 0;
    nTB[46] = 0;

    // (47):  OH + CH4 <=> H2O + CH3
    fwd_A[47]     = 1600000;
    fwd_beta[47]  = 2.1000000000000001;
    fwd_Ea[47]    = 2459.7330000000002;
    prefactor_units[47]  = 1.0000000000000002e-06;
    activation_units[47] = 0.50321666580471969;
    phase_units[47]      = 1e-12;
    is_PD[47] = 0;
    nTB[47] = 0;

    // (48):  O + CH4 <=> OH + CH3
    fwd_A[48]     = 1020000000;
    fwd_beta[48]  = 1.5;
    fwd_Ea[48]    = 8603.0669999999991;
    prefactor_units[48]  = 1.0000000000000002e-06;
    activation_units[48] = 0.50321666580471969;
    phase_units[48]      = 1e-12;
    is_PD[48] = 0;
    nTB[48] = 0;

    // (49):  HO2 + CH4 <=> H2O2 + CH3
    fwd_A[49]     = 900000000000;
    fwd_beta[49]  = 0;
    fwd_Ea[49]    = 18697.970000000001;
    prefactor_units[49]  = 1.0000000000000002e-06;
    activation_units[49] = 0.50321666580471969;
    phase_units[49]      = 1e-12;
    is_PD[49] = 0;
    nTB[49] = 0;

    // (50):  CH2 + CH4 <=> 2 CH3
    fwd_A[50]     = 4000000000000;
    fwd_beta[50]  = 0;
    fwd_Ea[50]    = -569.93820000000005;
    prefactor_units[50]  = 1.0000000000000002e-06;
    activation_units[50] = 0.50321666580471969;
    phase_units[50]      = 1e-12;
    is_PD[50] = 0;
    nTB[50] = 0;

    // (51):  C3H6 <=> CH3 + C2H3
    fwd_A[51]     = 3150000000000000;
    fwd_beta[51]  = 0;
    fwd_Ea[51]    = 85490.720000000001;
    prefactor_units[51]  = 1;
    activation_units[51] = 0.50321666580471969;
    phase_units[51]      = 1e-6;
    is_PD[51] = 0;
    nTB[51] = 0;

    SetAllDefaults();
}



/*A few mechanism parameters */
void CKINDX(int * iwrk, double * restrict rwrk, int * mm, int * kk, int * ii, int * nfit)
{
    *mm = 4;
    *kk = 29;
    *ii = 52;
    *nfit = -1; /*Why do you need this anyway ?  */
}



/* ckxnum... for parsing strings  */
void CKXNUM(char * line, int * nexp, int * lout, int * nval, double * restrict rval, int * kerr, int lenline )
{
    int n,i; /*Loop Counters */
    char cstr[1000];
    char *saveptr;
    char *p; /*String Tokens */
    /* Strip Comments  */
    for (i=0; i<lenline; ++i) {
        if (line[i]=='!') {
            break;
        }
        cstr[i] = line[i];
    }
    cstr[i] = '\0';

    p = strtok_r(cstr," ", &saveptr);
    if (!p) {
        *nval = 0;
        *kerr = 1;
        return;
    }
    for (n=0; n<*nexp; ++n) {
        rval[n] = atof(p);
        p = strtok_r(NULL, " ", &saveptr);
        if (!p) break;
    }
    *nval = n+1;
    if (*nval < *nexp) *kerr = 1;
    return;
}


/* cksnum... for parsing strings  */
void CKSNUM(char * line, int * nexp, int * lout, char * kray, int * nn, int * knum, int * nval, double * restrict rval, int * kerr, int lenline, int lenkray)
{
    /*Not done yet ... */
}


/* Returns the char strings of element names */
void CKSYME(int * kname, int * plenkname )
{
    int i; /*Loop Counter */
    int lenkname = *plenkname;
    /*clear kname */
    for (i=0; i<lenkname*4; i++) {
        kname[i] = ' ';
    }

    /* H  */
    kname[ 0*lenkname + 0 ] = 'H';
    kname[ 0*lenkname + 1 ] = ' ';

    /* C  */
    kname[ 1*lenkname + 0 ] = 'C';
    kname[ 1*lenkname + 1 ] = ' ';

    /* O  */
    kname[ 2*lenkname + 0 ] = 'O';
    kname[ 2*lenkname + 1 ] = ' ';

    /* N  */
    kname[ 3*lenkname + 0 ] = 'N';
    kname[ 3*lenkname + 1 ] = ' ';

}


/* Returns the char strings of species names */
void CKSYMS(int * kname, int * plenkname )
{
    int i; /*Loop Counter */
    int lenkname = *plenkname;
    /*clear kname */
    for (i=0; i<lenkname*29; i++) {
        kname[i] = ' ';
    }

    /* NC7H16  */
    kname[ 0*lenkname + 0 ] = 'N';
    kname[ 0*lenkname + 1 ] = 'C';
    kname[ 0*lenkname + 2 ] = '7';
    kname[ 0*lenkname + 3 ] = 'H';
    kname[ 0*lenkname + 4 ] = '1';
    kname[ 0*lenkname + 5 ] = '6';
    kname[ 0*lenkname + 6 ] = ' ';

    /* O2  */
    kname[ 1*lenkname + 0 ] = 'O';
    kname[ 1*lenkname + 1 ] = '2';
    kname[ 1*lenkname + 2 ] = ' ';

    /* CO2  */
    kname[ 2*lenkname + 0 ] = 'C';
    kname[ 2*lenkname + 1 ] = 'O';
    kname[ 2*lenkname + 2 ] = '2';
    kname[ 2*lenkname + 3 ] = ' ';

    /* H2O  */
    kname[ 3*lenkname + 0 ] = 'H';
    kname[ 3*lenkname + 1 ] = '2';
    kname[ 3*lenkname + 2 ] = 'O';
    kname[ 3*lenkname + 3 ] = ' ';

    /* CO  */
    kname[ 4*lenkname + 0 ] = 'C';
    kname[ 4*lenkname + 1 ] = 'O';
    kname[ 4*lenkname + 2 ] = ' ';

    /* H2  */
    kname[ 5*lenkname + 0 ] = 'H';
    kname[ 5*lenkname + 1 ] = '2';
    kname[ 5*lenkname + 2 ] = ' ';

    /* OH  */
    kname[ 6*lenkname + 0 ] = 'O';
    kname[ 6*lenkname + 1 ] = 'H';
    kname[ 6*lenkname + 2 ] = ' ';

    /* H2O2  */
    kname[ 7*lenkname + 0 ] = 'H';
    kname[ 7*lenkname + 1 ] = '2';
    kname[ 7*lenkname + 2 ] = 'O';
    kname[ 7*lenkname + 3 ] = '2';
    kname[ 7*lenkname + 4 ] = ' ';

    /* HO2  */
    kname[ 8*lenkname + 0 ] = 'H';
    kname[ 8*lenkname + 1 ] = 'O';
    kname[ 8*lenkname + 2 ] = '2';
    kname[ 8*lenkname + 3 ] = ' ';

    /* H  */
    kname[ 9*lenkname + 0 ] = 'H';
    kname[ 9*lenkname + 1 ] = ' ';

    /* O  */
    kname[ 10*lenkname + 0 ] = 'O';
    kname[ 10*lenkname + 1 ] = ' ';

    /* CH3O  */
    kname[ 11*lenkname + 0 ] = 'C';
    kname[ 11*lenkname + 1 ] = 'H';
    kname[ 11*lenkname + 2 ] = '3';
    kname[ 11*lenkname + 3 ] = 'O';
    kname[ 11*lenkname + 4 ] = ' ';

    /* CH2O  */
    kname[ 12*lenkname + 0 ] = 'C';
    kname[ 12*lenkname + 1 ] = 'H';
    kname[ 12*lenkname + 2 ] = '2';
    kname[ 12*lenkname + 3 ] = 'O';
    kname[ 12*lenkname + 4 ] = ' ';

    /* HCO  */
    kname[ 13*lenkname + 0 ] = 'H';
    kname[ 13*lenkname + 1 ] = 'C';
    kname[ 13*lenkname + 2 ] = 'O';
    kname[ 13*lenkname + 3 ] = ' ';

    /* CH2  */
    kname[ 14*lenkname + 0 ] = 'C';
    kname[ 14*lenkname + 1 ] = 'H';
    kname[ 14*lenkname + 2 ] = '2';
    kname[ 14*lenkname + 3 ] = ' ';

    /* CH3  */
    kname[ 15*lenkname + 0 ] = 'C';
    kname[ 15*lenkname + 1 ] = 'H';
    kname[ 15*lenkname + 2 ] = '3';
    kname[ 15*lenkname + 3 ] = ' ';

    /* CH4  */
    kname[ 16*lenkname + 0 ] = 'C';
    kname[ 16*lenkname + 1 ] = 'H';
    kname[ 16*lenkname + 2 ] = '4';
    kname[ 16*lenkname + 3 ] = ' ';

    /* C2H3  */
    kname[ 17*lenkname + 0 ] = 'C';
    kname[ 17*lenkname + 1 ] = '2';
    kname[ 17*lenkname + 2 ] = 'H';
    kname[ 17*lenkname + 3 ] = '3';
    kname[ 17*lenkname + 4 ] = ' ';

    /* C2H4  */
    kname[ 18*lenkname + 0 ] = 'C';
    kname[ 18*lenkname + 1 ] = '2';
    kname[ 18*lenkname + 2 ] = 'H';
    kname[ 18*lenkname + 3 ] = '4';
    kname[ 18*lenkname + 4 ] = ' ';

    /* C2H5  */
    kname[ 19*lenkname + 0 ] = 'C';
    kname[ 19*lenkname + 1 ] = '2';
    kname[ 19*lenkname + 2 ] = 'H';
    kname[ 19*lenkname + 3 ] = '5';
    kname[ 19*lenkname + 4 ] = ' ';

    /* C3H4  */
    kname[ 20*lenkname + 0 ] = 'C';
    kname[ 20*lenkname + 1 ] = '3';
    kname[ 20*lenkname + 2 ] = 'H';
    kname[ 20*lenkname + 3 ] = '4';
    kname[ 20*lenkname + 4 ] = ' ';

    /* C3H5  */
    kname[ 21*lenkname + 0 ] = 'C';
    kname[ 21*lenkname + 1 ] = '3';
    kname[ 21*lenkname + 2 ] = 'H';
    kname[ 21*lenkname + 3 ] = '5';
    kname[ 21*lenkname + 4 ] = ' ';

    /* C3H6  */
    kname[ 22*lenkname + 0 ] = 'C';
    kname[ 22*lenkname + 1 ] = '3';
    kname[ 22*lenkname + 2 ] = 'H';
    kname[ 22*lenkname + 3 ] = '6';
    kname[ 22*lenkname + 4 ] = ' ';

    /* C3H7  */
    kname[ 23*lenkname + 0 ] = 'C';
    kname[ 23*lenkname + 1 ] = '3';
    kname[ 23*lenkname + 2 ] = 'H';
    kname[ 23*lenkname + 3 ] = '7';
    kname[ 23*lenkname + 4 ] = ' ';

    /* C7H15-2  */
    kname[ 24*lenkname + 0 ] = 'C';
    kname[ 24*lenkname + 1 ] = '7';
    kname[ 24*lenkname + 2 ] = 'H';
    kname[ 24*lenkname + 3 ] = '1';
    kname[ 24*lenkname + 4 ] = '5';
    kname[ 24*lenkname + 5 ] = '-';
    kname[ 24*lenkname + 6 ] = '2';
    kname[ 24*lenkname + 7 ] = ' ';

    /* C7H15O2  */
    kname[ 25*lenkname + 0 ] = 'C';
    kname[ 25*lenkname + 1 ] = '7';
    kname[ 25*lenkname + 2 ] = 'H';
    kname[ 25*lenkname + 3 ] = '1';
    kname[ 25*lenkname + 4 ] = '5';
    kname[ 25*lenkname + 5 ] = 'O';
    kname[ 25*lenkname + 6 ] = '2';
    kname[ 25*lenkname + 7 ] = ' ';

    /* C7KET12  */
    kname[ 26*lenkname + 0 ] = 'C';
    kname[ 26*lenkname + 1 ] = '7';
    kname[ 26*lenkname + 2 ] = 'K';
    kname[ 26*lenkname + 3 ] = 'E';
    kname[ 26*lenkname + 4 ] = 'T';
    kname[ 26*lenkname + 5 ] = '1';
    kname[ 26*lenkname + 6 ] = '2';
    kname[ 26*lenkname + 7 ] = ' ';

    /* C5H11CO  */
    kname[ 27*lenkname + 0 ] = 'C';
    kname[ 27*lenkname + 1 ] = '5';
    kname[ 27*lenkname + 2 ] = 'H';
    kname[ 27*lenkname + 3 ] = '1';
    kname[ 27*lenkname + 4 ] = '1';
    kname[ 27*lenkname + 5 ] = 'C';
    kname[ 27*lenkname + 6 ] = 'O';
    kname[ 27*lenkname + 7 ] = ' ';

    /* N2  */
    kname[ 28*lenkname + 0 ] = 'N';
    kname[ 28*lenkname + 1 ] = '2';
    kname[ 28*lenkname + 2 ] = ' ';

}


/* Returns R, Rc, Patm */
void CKRP(int * ickwrk, double * restrict rckwrk, double * restrict ru, double * restrict ruc, double * restrict pa)
{
     *ru  = 8.31451e+07; 
     *ruc = 1.98721558317399615845; 
     *pa  = 1.01325e+06; 
}


/*Compute P = rhoRT/W(x) */
void CKPX(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict P)
{
    double XW = 0;/* To hold mean molecular wt */
    XW += x[0]*100.205570; /*NC7H16 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*44.009950; /*CO2 */
    XW += x[3]*18.015340; /*H2O */
    XW += x[4]*28.010550; /*CO */
    XW += x[5]*2.015940; /*H2 */
    XW += x[6]*17.007370; /*OH */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*33.006770; /*HO2 */
    XW += x[9]*1.007970; /*H */
    XW += x[10]*15.999400; /*O */
    XW += x[11]*31.034460; /*CH3O */
    XW += x[12]*30.026490; /*CH2O */
    XW += x[13]*29.018520; /*HCO */
    XW += x[14]*14.027090; /*CH2 */
    XW += x[15]*15.035060; /*CH3 */
    XW += x[16]*16.043030; /*CH4 */
    XW += x[17]*27.046210; /*C2H3 */
    XW += x[18]*28.054180; /*C2H4 */
    XW += x[19]*29.062150; /*C2H5 */
    XW += x[20]*40.065330; /*C3H4 */
    XW += x[21]*41.073300; /*C3H5 */
    XW += x[22]*42.081270; /*C3H6 */
    XW += x[23]*43.089240; /*C3H7 */
    XW += x[24]*99.197600; /*C7H15-2 */
    XW += x[25]*131.196400; /*C7H15O2 */
    XW += x[26]*146.187830; /*C7KET12 */
    XW += x[27]*99.153970; /*C5H11CO */
    XW += x[28]*28.013400; /*N2 */
    *P = *rho * 8.31451e+07 * (*T) / XW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(y) */
void CKPY(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict P)
{
    double YOW = 0;/* for computing mean MW */
    YOW += y[0]*imw[0]; /*NC7H16 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*CO2 */
    YOW += y[3]*imw[3]; /*H2O */
    YOW += y[4]*imw[4]; /*CO */
    YOW += y[5]*imw[5]; /*H2 */
    YOW += y[6]*imw[6]; /*OH */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*HO2 */
    YOW += y[9]*imw[9]; /*H */
    YOW += y[10]*imw[10]; /*O */
    YOW += y[11]*imw[11]; /*CH3O */
    YOW += y[12]*imw[12]; /*CH2O */
    YOW += y[13]*imw[13]; /*HCO */
    YOW += y[14]*imw[14]; /*CH2 */
    YOW += y[15]*imw[15]; /*CH3 */
    YOW += y[16]*imw[16]; /*CH4 */
    YOW += y[17]*imw[17]; /*C2H3 */
    YOW += y[18]*imw[18]; /*C2H4 */
    YOW += y[19]*imw[19]; /*C2H5 */
    YOW += y[20]*imw[20]; /*C3H4 */
    YOW += y[21]*imw[21]; /*C3H5 */
    YOW += y[22]*imw[22]; /*C3H6 */
    YOW += y[23]*imw[23]; /*C3H7 */
    YOW += y[24]*imw[24]; /*C7H15-2 */
    YOW += y[25]*imw[25]; /*C7H15O2 */
    YOW += y[26]*imw[26]; /*C7KET12 */
    YOW += y[27]*imw[27]; /*C5H11CO */
    YOW += y[28]*imw[28]; /*N2 */
    *P = *rho * 8.31451e+07 * (*T) * YOW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(y) */
void VCKPY(int * restrict np, double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict P)
{
    double YOW[*np];
    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

    for (int n=0; n<29; n++) {
        for (int i=0; i<(*np); i++) {
            YOW[i] += y[n*(*np)+i] * imw[n];
        }
    }

    for (int i=0; i<(*np); i++) {
        P[i] = rho[i] * 8.31451e+07 * T[i] * YOW[i]; /*P = rho*R*T/W */
    }

    return;
}


/*Compute P = rhoRT/W(c) */
void CKPC(double * restrict rho, double * restrict T, double * restrict c, int * iwrk, double * restrict rwrk, double * restrict P)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    double W = 0;
    double sumC = 0;
    W += c[0]*100.205570; /*NC7H16 */
    W += c[1]*31.998800; /*O2 */
    W += c[2]*44.009950; /*CO2 */
    W += c[3]*18.015340; /*H2O */
    W += c[4]*28.010550; /*CO */
    W += c[5]*2.015940; /*H2 */
    W += c[6]*17.007370; /*OH */
    W += c[7]*34.014740; /*H2O2 */
    W += c[8]*33.006770; /*HO2 */
    W += c[9]*1.007970; /*H */
    W += c[10]*15.999400; /*O */
    W += c[11]*31.034460; /*CH3O */
    W += c[12]*30.026490; /*CH2O */
    W += c[13]*29.018520; /*HCO */
    W += c[14]*14.027090; /*CH2 */
    W += c[15]*15.035060; /*CH3 */
    W += c[16]*16.043030; /*CH4 */
    W += c[17]*27.046210; /*C2H3 */
    W += c[18]*28.054180; /*C2H4 */
    W += c[19]*29.062150; /*C2H5 */
    W += c[20]*40.065330; /*C3H4 */
    W += c[21]*41.073300; /*C3H5 */
    W += c[22]*42.081270; /*C3H6 */
    W += c[23]*43.089240; /*C3H7 */
    W += c[24]*99.197600; /*C7H15-2 */
    W += c[25]*131.196400; /*C7H15O2 */
    W += c[26]*146.187830; /*C7KET12 */
    W += c[27]*99.153970; /*C5H11CO */
    W += c[28]*28.013400; /*N2 */

    for (id = 0; id < 29; ++id) {
        sumC += c[id];
    }
    *P = *rho * 8.31451e+07 * (*T) * sumC / W; /*P = rho*R*T/W */

    return;
}


/*Compute rho = PW(x)/RT */
void CKRHOX(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict rho)
{
    double XW = 0;/* To hold mean molecular wt */
    XW += x[0]*100.205570; /*NC7H16 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*44.009950; /*CO2 */
    XW += x[3]*18.015340; /*H2O */
    XW += x[4]*28.010550; /*CO */
    XW += x[5]*2.015940; /*H2 */
    XW += x[6]*17.007370; /*OH */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*33.006770; /*HO2 */
    XW += x[9]*1.007970; /*H */
    XW += x[10]*15.999400; /*O */
    XW += x[11]*31.034460; /*CH3O */
    XW += x[12]*30.026490; /*CH2O */
    XW += x[13]*29.018520; /*HCO */
    XW += x[14]*14.027090; /*CH2 */
    XW += x[15]*15.035060; /*CH3 */
    XW += x[16]*16.043030; /*CH4 */
    XW += x[17]*27.046210; /*C2H3 */
    XW += x[18]*28.054180; /*C2H4 */
    XW += x[19]*29.062150; /*C2H5 */
    XW += x[20]*40.065330; /*C3H4 */
    XW += x[21]*41.073300; /*C3H5 */
    XW += x[22]*42.081270; /*C3H6 */
    XW += x[23]*43.089240; /*C3H7 */
    XW += x[24]*99.197600; /*C7H15-2 */
    XW += x[25]*131.196400; /*C7H15O2 */
    XW += x[26]*146.187830; /*C7KET12 */
    XW += x[27]*99.153970; /*C5H11CO */
    XW += x[28]*28.013400; /*N2 */
    *rho = *P * XW / (8.31451e+07 * (*T)); /*rho = P*W/(R*T) */

    return;
}


/*Compute rho = P*W(y)/RT */
void CKRHOY(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict rho)
{
    double YOW = 0;
    double tmp[29];

    for (int i = 0; i < 29; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 29; i++)
    {
        YOW += tmp[i];
    }

    *rho = *P / (8.31451e+07 * (*T) * YOW);/*rho = P*W/(R*T) */
    return;
}


/*Compute rho = P*W(c)/(R*T) */
void CKRHOC(double * restrict P, double * restrict T, double * restrict c, int * iwrk, double * restrict rwrk, double * restrict rho)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    double W = 0;
    double sumC = 0;
    W += c[0]*100.205570; /*NC7H16 */
    W += c[1]*31.998800; /*O2 */
    W += c[2]*44.009950; /*CO2 */
    W += c[3]*18.015340; /*H2O */
    W += c[4]*28.010550; /*CO */
    W += c[5]*2.015940; /*H2 */
    W += c[6]*17.007370; /*OH */
    W += c[7]*34.014740; /*H2O2 */
    W += c[8]*33.006770; /*HO2 */
    W += c[9]*1.007970; /*H */
    W += c[10]*15.999400; /*O */
    W += c[11]*31.034460; /*CH3O */
    W += c[12]*30.026490; /*CH2O */
    W += c[13]*29.018520; /*HCO */
    W += c[14]*14.027090; /*CH2 */
    W += c[15]*15.035060; /*CH3 */
    W += c[16]*16.043030; /*CH4 */
    W += c[17]*27.046210; /*C2H3 */
    W += c[18]*28.054180; /*C2H4 */
    W += c[19]*29.062150; /*C2H5 */
    W += c[20]*40.065330; /*C3H4 */
    W += c[21]*41.073300; /*C3H5 */
    W += c[22]*42.081270; /*C3H6 */
    W += c[23]*43.089240; /*C3H7 */
    W += c[24]*99.197600; /*C7H15-2 */
    W += c[25]*131.196400; /*C7H15O2 */
    W += c[26]*146.187830; /*C7KET12 */
    W += c[27]*99.153970; /*C5H11CO */
    W += c[28]*28.013400; /*N2 */

    for (id = 0; id < 29; ++id) {
        sumC += c[id];
    }
    *rho = *P * W / (sumC * (*T) * 8.31451e+07); /*rho = PW/(R*T) */

    return;
}


/*get molecular weight for all species */
void CKWT(int * iwrk, double * restrict rwrk, double * restrict wt)
{
    molecularWeight(wt);
}


/*get atomic weight for all elements */
void CKAWT(int * iwrk, double * restrict rwrk, double * restrict awt)
{
    atomicWeight(awt);
}


/*given y[species]: mass fractions */
/*returns mean molecular weight (gm/mole) */
void CKMMWY(double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wtm)
{
    double YOW = 0;
    double tmp[29];

    for (int i = 0; i < 29; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 29; i++)
    {
        YOW += tmp[i];
    }

    *wtm = 1.0 / YOW;
    return;
}


/*given x[species]: mole fractions */
/*returns mean molecular weight (gm/mole) */
void CKMMWX(double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wtm)
{
    double XW = 0;/* see Eq 4 in CK Manual */
    XW += x[0]*100.205570; /*NC7H16 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*44.009950; /*CO2 */
    XW += x[3]*18.015340; /*H2O */
    XW += x[4]*28.010550; /*CO */
    XW += x[5]*2.015940; /*H2 */
    XW += x[6]*17.007370; /*OH */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*33.006770; /*HO2 */
    XW += x[9]*1.007970; /*H */
    XW += x[10]*15.999400; /*O */
    XW += x[11]*31.034460; /*CH3O */
    XW += x[12]*30.026490; /*CH2O */
    XW += x[13]*29.018520; /*HCO */
    XW += x[14]*14.027090; /*CH2 */
    XW += x[15]*15.035060; /*CH3 */
    XW += x[16]*16.043030; /*CH4 */
    XW += x[17]*27.046210; /*C2H3 */
    XW += x[18]*28.054180; /*C2H4 */
    XW += x[19]*29.062150; /*C2H5 */
    XW += x[20]*40.065330; /*C3H4 */
    XW += x[21]*41.073300; /*C3H5 */
    XW += x[22]*42.081270; /*C3H6 */
    XW += x[23]*43.089240; /*C3H7 */
    XW += x[24]*99.197600; /*C7H15-2 */
    XW += x[25]*131.196400; /*C7H15O2 */
    XW += x[26]*146.187830; /*C7KET12 */
    XW += x[27]*99.153970; /*C5H11CO */
    XW += x[28]*28.013400; /*N2 */
    *wtm = XW;

    return;
}


/*given c[species]: molar concentration */
/*returns mean molecular weight (gm/mole) */
void CKMMWC(double * restrict c, int * iwrk, double * restrict rwrk, double * restrict wtm)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    double W = 0;
    double sumC = 0;
    W += c[0]*100.205570; /*NC7H16 */
    W += c[1]*31.998800; /*O2 */
    W += c[2]*44.009950; /*CO2 */
    W += c[3]*18.015340; /*H2O */
    W += c[4]*28.010550; /*CO */
    W += c[5]*2.015940; /*H2 */
    W += c[6]*17.007370; /*OH */
    W += c[7]*34.014740; /*H2O2 */
    W += c[8]*33.006770; /*HO2 */
    W += c[9]*1.007970; /*H */
    W += c[10]*15.999400; /*O */
    W += c[11]*31.034460; /*CH3O */
    W += c[12]*30.026490; /*CH2O */
    W += c[13]*29.018520; /*HCO */
    W += c[14]*14.027090; /*CH2 */
    W += c[15]*15.035060; /*CH3 */
    W += c[16]*16.043030; /*CH4 */
    W += c[17]*27.046210; /*C2H3 */
    W += c[18]*28.054180; /*C2H4 */
    W += c[19]*29.062150; /*C2H5 */
    W += c[20]*40.065330; /*C3H4 */
    W += c[21]*41.073300; /*C3H5 */
    W += c[22]*42.081270; /*C3H6 */
    W += c[23]*43.089240; /*C3H7 */
    W += c[24]*99.197600; /*C7H15-2 */
    W += c[25]*131.196400; /*C7H15O2 */
    W += c[26]*146.187830; /*C7KET12 */
    W += c[27]*99.153970; /*C5H11CO */
    W += c[28]*28.013400; /*N2 */

    for (id = 0; id < 29; ++id) {
        sumC += c[id];
    }
    /* CK provides no guard against divison by zero */
    *wtm = W/sumC;

    return;
}


/*convert y[species] (mass fracs) to x[species] (mole fracs) */
void CKYTX(double * restrict y, int * iwrk, double * restrict rwrk, double * restrict x)
{
    double YOW = 0;
    double tmp[29];

    for (int i = 0; i < 29; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 29; i++)
    {
        YOW += tmp[i];
    }

    double YOWINV = 1.0/YOW;

    for (int i = 0; i < 29; i++)
    {
        x[i] = y[i]*imw[i]*YOWINV;
    }
    return;
}


/*convert y[npoints*species] (mass fracs) to x[npoints*species] (mole fracs) */
void VCKYTX(int * restrict np, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict x)
{
    double YOW[*np];
    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

    for (int n=0; n<29; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] = y[n*(*np)+i] * imw[n];
            YOW[i] += x[n*(*np)+i];
        }
    }

    for (int i=0; i<(*np); i++) {
        YOW[i] = 1.0/YOW[i];
    }

    for (int n=0; n<29; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] *=  YOW[i];
        }
    }
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
void CKYTCP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict c)
{
    double YOW = 0;
    double PWORT;

    /*Compute inverse of mean molecular wt first */
    for (int i = 0; i < 29; i++)
    {
        c[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 29; i++)
    {
        YOW += c[i];
    }

    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31451e+07 * (*T)); 
    /*Now compute conversion */

    for (int i = 0; i < 29; i++)
    {
        c[i] = PWORT * y[i] * imw[i];
    }
    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
void CKYTCR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict c)
{
    for (int i = 0; i < 29; i++)
    {
        c[i] = (*rho)  * y[i] * imw[i];
    }
}

/*convert y[species] (mass fracs) to c[species] (molar conc) */
void VCKYTCR(int * restrict np, double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict c)
{
    for (int n=0; n<29; n++) {
        for (int i=0; i<(*np); i++) {
            c[i*(29)+n] = y[i*(29)+n] * rho[i] * imw[n];
        }
    }
}


/*convert x[species] (mole fracs) to y[species] (mass fracs) */
void CKXTY(double * restrict x, int * iwrk, double * restrict rwrk, double * restrict y)
{
    double XW = 0; /*See Eq 4, 9 in CK Manual */
    /*Compute mean molecular wt first */
    XW += x[0]*100.205570; /*NC7H16 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*44.009950; /*CO2 */
    XW += x[3]*18.015340; /*H2O */
    XW += x[4]*28.010550; /*CO */
    XW += x[5]*2.015940; /*H2 */
    XW += x[6]*17.007370; /*OH */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*33.006770; /*HO2 */
    XW += x[9]*1.007970; /*H */
    XW += x[10]*15.999400; /*O */
    XW += x[11]*31.034460; /*CH3O */
    XW += x[12]*30.026490; /*CH2O */
    XW += x[13]*29.018520; /*HCO */
    XW += x[14]*14.027090; /*CH2 */
    XW += x[15]*15.035060; /*CH3 */
    XW += x[16]*16.043030; /*CH4 */
    XW += x[17]*27.046210; /*C2H3 */
    XW += x[18]*28.054180; /*C2H4 */
    XW += x[19]*29.062150; /*C2H5 */
    XW += x[20]*40.065330; /*C3H4 */
    XW += x[21]*41.073300; /*C3H5 */
    XW += x[22]*42.081270; /*C3H6 */
    XW += x[23]*43.089240; /*C3H7 */
    XW += x[24]*99.197600; /*C7H15-2 */
    XW += x[25]*131.196400; /*C7H15O2 */
    XW += x[26]*146.187830; /*C7KET12 */
    XW += x[27]*99.153970; /*C5H11CO */
    XW += x[28]*28.013400; /*N2 */
    /*Now compute conversion */
    double XWinv = 1.0/XW;
    y[0] = x[0]*100.205570*XWinv; 
    y[1] = x[1]*31.998800*XWinv; 
    y[2] = x[2]*44.009950*XWinv; 
    y[3] = x[3]*18.015340*XWinv; 
    y[4] = x[4]*28.010550*XWinv; 
    y[5] = x[5]*2.015940*XWinv; 
    y[6] = x[6]*17.007370*XWinv; 
    y[7] = x[7]*34.014740*XWinv; 
    y[8] = x[8]*33.006770*XWinv; 
    y[9] = x[9]*1.007970*XWinv; 
    y[10] = x[10]*15.999400*XWinv; 
    y[11] = x[11]*31.034460*XWinv; 
    y[12] = x[12]*30.026490*XWinv; 
    y[13] = x[13]*29.018520*XWinv; 
    y[14] = x[14]*14.027090*XWinv; 
    y[15] = x[15]*15.035060*XWinv; 
    y[16] = x[16]*16.043030*XWinv; 
    y[17] = x[17]*27.046210*XWinv; 
    y[18] = x[18]*28.054180*XWinv; 
    y[19] = x[19]*29.062150*XWinv; 
    y[20] = x[20]*40.065330*XWinv; 
    y[21] = x[21]*41.073300*XWinv; 
    y[22] = x[22]*42.081270*XWinv; 
    y[23] = x[23]*43.089240*XWinv; 
    y[24] = x[24]*99.197600*XWinv; 
    y[25] = x[25]*131.196400*XWinv; 
    y[26] = x[26]*146.187830*XWinv; 
    y[27] = x[27]*99.153970*XWinv; 
    y[28] = x[28]*28.013400*XWinv; 

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void CKXTCP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict c)
{
    int id; /*loop counter */
    double PORT = (*P)/(8.31451e+07 * (*T)); /*P/RT */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 29; ++id) {
        c[id] = x[id]*PORT;
    }

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void CKXTCR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict c)
{
    int id; /*loop counter */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*100.205570; /*NC7H16 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*44.009950; /*CO2 */
    XW += x[3]*18.015340; /*H2O */
    XW += x[4]*28.010550; /*CO */
    XW += x[5]*2.015940; /*H2 */
    XW += x[6]*17.007370; /*OH */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*33.006770; /*HO2 */
    XW += x[9]*1.007970; /*H */
    XW += x[10]*15.999400; /*O */
    XW += x[11]*31.034460; /*CH3O */
    XW += x[12]*30.026490; /*CH2O */
    XW += x[13]*29.018520; /*HCO */
    XW += x[14]*14.027090; /*CH2 */
    XW += x[15]*15.035060; /*CH3 */
    XW += x[16]*16.043030; /*CH4 */
    XW += x[17]*27.046210; /*C2H3 */
    XW += x[18]*28.054180; /*C2H4 */
    XW += x[19]*29.062150; /*C2H5 */
    XW += x[20]*40.065330; /*C3H4 */
    XW += x[21]*41.073300; /*C3H5 */
    XW += x[22]*42.081270; /*C3H6 */
    XW += x[23]*43.089240; /*C3H7 */
    XW += x[24]*99.197600; /*C7H15-2 */
    XW += x[25]*131.196400; /*C7H15O2 */
    XW += x[26]*146.187830; /*C7KET12 */
    XW += x[27]*99.153970; /*C5H11CO */
    XW += x[28]*28.013400; /*N2 */
    ROW = (*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 29; ++id) {
        c[id] = x[id]*ROW;
    }

    return;
}


/*convert c[species] (molar conc) to x[species] (mole fracs) */
void CKCTX(double * restrict c, int * iwrk, double * restrict rwrk, double * restrict x)
{
    int id; /*loop counter */
    double sumC = 0; 

    /*compute sum of c  */
    for (id = 0; id < 29; ++id) {
        sumC += c[id];
    }

    /* See Eq 13  */
    double sumCinv = 1.0/sumC;
    for (id = 0; id < 29; ++id) {
        x[id] = c[id]*sumCinv;
    }

    return;
}


/*convert c[species] (molar conc) to y[species] (mass fracs) */
void CKCTY(double * restrict c, int * iwrk, double * restrict rwrk, double * restrict y)
{
    double CW = 0; /*See Eq 12 in CK Manual */
    /*compute denominator in eq 12 first */
    CW += c[0]*100.205570; /*NC7H16 */
    CW += c[1]*31.998800; /*O2 */
    CW += c[2]*44.009950; /*CO2 */
    CW += c[3]*18.015340; /*H2O */
    CW += c[4]*28.010550; /*CO */
    CW += c[5]*2.015940; /*H2 */
    CW += c[6]*17.007370; /*OH */
    CW += c[7]*34.014740; /*H2O2 */
    CW += c[8]*33.006770; /*HO2 */
    CW += c[9]*1.007970; /*H */
    CW += c[10]*15.999400; /*O */
    CW += c[11]*31.034460; /*CH3O */
    CW += c[12]*30.026490; /*CH2O */
    CW += c[13]*29.018520; /*HCO */
    CW += c[14]*14.027090; /*CH2 */
    CW += c[15]*15.035060; /*CH3 */
    CW += c[16]*16.043030; /*CH4 */
    CW += c[17]*27.046210; /*C2H3 */
    CW += c[18]*28.054180; /*C2H4 */
    CW += c[19]*29.062150; /*C2H5 */
    CW += c[20]*40.065330; /*C3H4 */
    CW += c[21]*41.073300; /*C3H5 */
    CW += c[22]*42.081270; /*C3H6 */
    CW += c[23]*43.089240; /*C3H7 */
    CW += c[24]*99.197600; /*C7H15-2 */
    CW += c[25]*131.196400; /*C7H15O2 */
    CW += c[26]*146.187830; /*C7KET12 */
    CW += c[27]*99.153970; /*C5H11CO */
    CW += c[28]*28.013400; /*N2 */
    /*Now compute conversion */
    double CWinv = 1.0/CW;
    y[0] = c[0]*100.205570*CWinv; 
    y[1] = c[1]*31.998800*CWinv; 
    y[2] = c[2]*44.009950*CWinv; 
    y[3] = c[3]*18.015340*CWinv; 
    y[4] = c[4]*28.010550*CWinv; 
    y[5] = c[5]*2.015940*CWinv; 
    y[6] = c[6]*17.007370*CWinv; 
    y[7] = c[7]*34.014740*CWinv; 
    y[8] = c[8]*33.006770*CWinv; 
    y[9] = c[9]*1.007970*CWinv; 
    y[10] = c[10]*15.999400*CWinv; 
    y[11] = c[11]*31.034460*CWinv; 
    y[12] = c[12]*30.026490*CWinv; 
    y[13] = c[13]*29.018520*CWinv; 
    y[14] = c[14]*14.027090*CWinv; 
    y[15] = c[15]*15.035060*CWinv; 
    y[16] = c[16]*16.043030*CWinv; 
    y[17] = c[17]*27.046210*CWinv; 
    y[18] = c[18]*28.054180*CWinv; 
    y[19] = c[19]*29.062150*CWinv; 
    y[20] = c[20]*40.065330*CWinv; 
    y[21] = c[21]*41.073300*CWinv; 
    y[22] = c[22]*42.081270*CWinv; 
    y[23] = c[23]*43.089240*CWinv; 
    y[24] = c[24]*99.197600*CWinv; 
    y[25] = c[25]*131.196400*CWinv; 
    y[26] = c[26]*146.187830*CWinv; 
    y[27] = c[27]*99.153970*CWinv; 
    y[28] = c[28]*28.013400*CWinv; 

    return;
}


/*get Cp/R as a function of T  */
/*for all species (Eq 19) */
void CKCPOR(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cpor)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpor, tc);
}


/*get H/RT as a function of T  */
/*for all species (Eq 20) */
void CKHORT(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict hort)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEnthalpy(hort, tc);
}


/*get S/R as a function of T  */
/*for all species (Eq 21) */
void CKSOR(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict sor)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sor, tc);
}


/*get specific heat at constant volume as a function  */
/*of T for all species (molar units) */
void CKCVML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cvml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cv_R(cvml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 29; ++id) {
        cvml[id] *= 8.31451e+07;
    }
}


/*get specific heat at constant pressure as a  */
/*function of T for all species (molar units) */
void CKCPML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cpml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 29; ++id) {
        cpml[id] *= 8.31451e+07;
    }
}


/*get internal energy as a function  */
/*of T for all species (molar units) */
void CKUML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict uml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 29; ++id) {
        uml[id] *= RT;
    }
}


/*get enthalpy as a function  */
/*of T for all species (molar units) */
void CKHML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict hml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 29; ++id) {
        hml[id] *= RT;
    }
}


/*get standard-state Gibbs energy as a function  */
/*of T for all species (molar units) */
void CKGML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict gml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    gibbs(gml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 29; ++id) {
        gml[id] *= RT;
    }
}


/*get standard-state Helmholtz free energy as a  */
/*function of T for all species (molar units) */
void CKAML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict aml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    helmholtz(aml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 29; ++id) {
        aml[id] *= RT;
    }
}


/*Returns the standard-state entropies in molar units */
void CKSML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict sml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 29; ++id) {
        sml[id] *= 8.31451e+07;
    }
}


/*Returns the specific heats at constant volume */
/*in mass units (Eq. 29) */
void CKCVMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cvms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cv_R(cvms, tc);
    /*multiply by R/molecularweight */
    cvms[0] *= 8.297452926019980e+05; /*NC7H16 */
    cvms[1] *= 2.598381814318037e+06; /*O2 */
    cvms[2] *= 1.889234139098090e+06; /*CO2 */
    cvms[3] *= 4.615239012974499e+06; /*H2O */
    cvms[4] *= 2.968349425484326e+06; /*CO */
    cvms[5] *= 4.124383662212169e+07; /*H2 */
    cvms[6] *= 4.888768810227566e+06; /*OH */
    cvms[7] *= 2.444384405113783e+06; /*H2O2 */
    cvms[8] *= 2.519031701678171e+06; /*HO2 */
    cvms[9] *= 8.248767324424338e+07; /*H */
    cvms[10] *= 5.196763628636074e+06; /*O */
    cvms[11] *= 2.679121853578248e+06; /*CH3O */
    cvms[12] *= 2.769058254894261e+06; /*CH2O */
    cvms[13] *= 2.865242610581105e+06; /*HCO */
    cvms[14] *= 5.927466067445207e+06; /*CH2 */
    cvms[15] *= 5.530081023953346e+06; /*CH3 */
    cvms[16] *= 5.182630712527496e+06; /*CH4 */
    cvms[17] *= 3.074186734481467e+06; /*C2H3 */
    cvms[18] *= 2.963733033722604e+06; /*C2H4 */
    cvms[19] *= 2.860941121011349e+06; /*C2H5 */
    cvms[20] *= 2.075238117344846e+06; /*C3H4 */
    cvms[21] *= 2.024310196648431e+06; /*C3H5 */
    cvms[22] *= 1.975822022481736e+06; /*C3H6 */
    cvms[23] *= 1.929602378691293e+06; /*C3H7 */
    cvms[24] *= 8.381765284643982e+05; /*C7H15-2 */
    cvms[25] *= 6.337452856938147e+05; /*C7H15O2 */
    cvms[26] *= 5.687552787396871e+05; /*C7KET12 */
    cvms[27] *= 8.385453451838590e+05; /*C5H11CO */
    cvms[28] *= 2.968047434442088e+06; /*N2 */
}

/*Returns the specific heats at constant volume */
/*in mass units (Eq. 29) */
void VCKCVMS(int * restrict np, double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cvms)
{
    double cvms_tmp[29*(*np)]; /* temporary energy array */
    //double tT = *T; /*temporary temperature */
    //double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    vcv_R(*np, cvms_tmp, T);
    /*multiply by R/molecularweight */
    for (int i=0; i<(*np); i++) {
        cvms[i*(29)+0]  = cvms_tmp[0*(*np)+ i] * 8.297452926019980e+05; /*NC7H16 */
        cvms[i*(29)+1]  = cvms_tmp[1*(*np)+ i] * 2.598381814318037e+06; /*O2 */
        cvms[i*(29)+2]  = cvms_tmp[2*(*np)+ i] * 1.889234139098090e+06; /*CO2 */
        cvms[i*(29)+3]  = cvms_tmp[3*(*np)+ i] * 4.615239012974499e+06; /*H2O */
        cvms[i*(29)+4]  = cvms_tmp[4*(*np)+ i] * 2.968349425484326e+06; /*CO */
        cvms[i*(29)+5]  = cvms_tmp[5*(*np)+ i] * 4.124383662212169e+07; /*H2 */
        cvms[i*(29)+6]  = cvms_tmp[6*(*np)+ i] * 4.888768810227566e+06; /*OH */
        cvms[i*(29)+7]  = cvms_tmp[7*(*np)+ i] * 2.444384405113783e+06; /*H2O2 */
        cvms[i*(29)+8]  = cvms_tmp[8*(*np)+ i] * 2.519031701678171e+06; /*HO2 */
        cvms[i*(29)+9]  = cvms_tmp[9*(*np)+ i] * 8.248767324424338e+07; /*H */
        cvms[i*(29)+10] = cvms_tmp[10*(*np)+ i] * 5.196763628636074e+06; /*O */
        cvms[i*(29)+11] = cvms_tmp[11*(*np)+ i] * 2.679121853578248e+06; /*CH3O */
        cvms[i*(29)+12] = cvms_tmp[12*(*np)+ i] * 2.769058254894261e+06; /*CH2O */
        cvms[i*(29)+13] = cvms_tmp[13*(*np)+ i] * 2.865242610581105e+06; /*HCO */
        cvms[i*(29)+14] = cvms_tmp[14*(*np)+ i] * 5.927466067445207e+06; /*CH2 */
        cvms[i*(29)+15] = cvms_tmp[15*(*np)+ i] * 5.530081023953346e+06; /*CH3 */
        cvms[i*(29)+16] = cvms_tmp[16*(*np)+ i] * 5.182630712527496e+06; /*CH4 */
        cvms[i*(29)+17] = cvms_tmp[17*(*np)+ i] * 3.074186734481467e+06; /*C2H3 */
        cvms[i*(29)+18] = cvms_tmp[18*(*np)+ i] * 2.963733033722604e+06; /*C2H4 */
        cvms[i*(29)+19] = cvms_tmp[19*(*np)+ i] * 2.860941121011349e+06; /*C2H5 */
        cvms[i*(29)+20] = cvms_tmp[20*(*np)+ i] * 2.075238117344846e+06; /*C3H4 */
        cvms[i*(29)+21] = cvms_tmp[21*(*np)+ i] * 2.024310196648431e+06; /*C3H5 */
        cvms[i*(29)+22] = cvms_tmp[22*(*np)+ i] * 1.975822022481736e+06; /*C3H6 */
        cvms[i*(29)+23] = cvms_tmp[23*(*np)+ i] * 1.929602378691293e+06; /*C3H7 */
        cvms[i*(29)+24] = cvms_tmp[24*(*np)+ i] * 8.381765284643982e+05; /*C7H15-2 */
        cvms[i*(29)+25] = cvms_tmp[25*(*np)+ i] * 6.337452856938147e+05; /*C7H15O2 */
        cvms[i*(29)+26] = cvms_tmp[26*(*np)+ i] * 5.687552787396871e+05; /*C7KET12 */
        cvms[i*(29)+27] = cvms_tmp[27*(*np)+ i] * 8.385453451838590e+05; /*C5H11CO */
        cvms[i*(29)+28] = cvms_tmp[28*(*np)+ i] * 2.968047434442088e+06; /*N2 */
    }
}


/*Returns the specific heats at constant pressure */
/*in mass units (Eq. 26) */
void CKCPMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cpms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpms, tc);
    /*multiply by R/molecularweight */
    cpms[0] *= 8.297452926019980e+05; /*NC7H16 */
    cpms[1] *= 2.598381814318037e+06; /*O2 */
    cpms[2] *= 1.889234139098090e+06; /*CO2 */
    cpms[3] *= 4.615239012974499e+06; /*H2O */
    cpms[4] *= 2.968349425484326e+06; /*CO */
    cpms[5] *= 4.124383662212169e+07; /*H2 */
    cpms[6] *= 4.888768810227566e+06; /*OH */
    cpms[7] *= 2.444384405113783e+06; /*H2O2 */
    cpms[8] *= 2.519031701678171e+06; /*HO2 */
    cpms[9] *= 8.248767324424338e+07; /*H */
    cpms[10] *= 5.196763628636074e+06; /*O */
    cpms[11] *= 2.679121853578248e+06; /*CH3O */
    cpms[12] *= 2.769058254894261e+06; /*CH2O */
    cpms[13] *= 2.865242610581105e+06; /*HCO */
    cpms[14] *= 5.927466067445207e+06; /*CH2 */
    cpms[15] *= 5.530081023953346e+06; /*CH3 */
    cpms[16] *= 5.182630712527496e+06; /*CH4 */
    cpms[17] *= 3.074186734481467e+06; /*C2H3 */
    cpms[18] *= 2.963733033722604e+06; /*C2H4 */
    cpms[19] *= 2.860941121011349e+06; /*C2H5 */
    cpms[20] *= 2.075238117344846e+06; /*C3H4 */
    cpms[21] *= 2.024310196648431e+06; /*C3H5 */
    cpms[22] *= 1.975822022481736e+06; /*C3H6 */
    cpms[23] *= 1.929602378691293e+06; /*C3H7 */
    cpms[24] *= 8.381765284643982e+05; /*C7H15-2 */
    cpms[25] *= 6.337452856938147e+05; /*C7H15O2 */
    cpms[26] *= 5.687552787396871e+05; /*C7KET12 */
    cpms[27] *= 8.385453451838590e+05; /*C5H11CO */
    cpms[28] *= 2.968047434442088e+06; /*N2 */
}


/*Returns internal energy in mass units (Eq 30.) */
void CKUMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict ums)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    for (int i = 0; i < 29; i++)
    {
        ums[i] *= RT*imw[i];
    }
}

/*Returns internal energy in mass units (Eq 30.) */
void VCKUMS(int * restrict np, double * restrict T, int * iwrk, double * restrict rwrk, double * restrict ums)
{
    double ums_tmp[29*(*np)]; /* temporary energy array */
    vspeciesInternalEnergy(*np, ums_tmp, T);
    for (int n=0; n<29; n++) {
        for (int i=0; i<(*np); i++) {
            ums[i*(29)+n] = ums_tmp[n*(*np)+i] * 8.31451e+07 * T[i] * imw[n];
        }
    }
}


/*Returns enthalpy in mass units (Eq 27.) */
void CKHMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict hms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hms, tc);
    for (int i = 0; i < 29; i++)
    {
        hms[i] *= RT*imw[i];
    }
}


/*Returns enthalpy in mass units (Eq 27.) */
void VCKHMS(int * restrict np, double * restrict T, int * iwrk, double * restrict rwrk, double * restrict hms)
{
    double tc[5], h[29];

    for (int i=0; i<(*np); i++) {
        tc[0] = 0.0;
        tc[1] = T[i];
        tc[2] = T[i]*T[i];
        tc[3] = T[i]*T[i]*T[i];
        tc[4] = T[i]*T[i]*T[i]*T[i];

        speciesEnthalpy(h, tc);

        hms[0*(*np)+i] = h[0];
        hms[1*(*np)+i] = h[1];
        hms[2*(*np)+i] = h[2];
        hms[3*(*np)+i] = h[3];
        hms[4*(*np)+i] = h[4];
        hms[5*(*np)+i] = h[5];
        hms[6*(*np)+i] = h[6];
        hms[7*(*np)+i] = h[7];
        hms[8*(*np)+i] = h[8];
        hms[9*(*np)+i] = h[9];
        hms[10*(*np)+i] = h[10];
        hms[11*(*np)+i] = h[11];
        hms[12*(*np)+i] = h[12];
        hms[13*(*np)+i] = h[13];
        hms[14*(*np)+i] = h[14];
        hms[15*(*np)+i] = h[15];
        hms[16*(*np)+i] = h[16];
        hms[17*(*np)+i] = h[17];
        hms[18*(*np)+i] = h[18];
        hms[19*(*np)+i] = h[19];
        hms[20*(*np)+i] = h[20];
        hms[21*(*np)+i] = h[21];
        hms[22*(*np)+i] = h[22];
        hms[23*(*np)+i] = h[23];
        hms[24*(*np)+i] = h[24];
        hms[25*(*np)+i] = h[25];
        hms[26*(*np)+i] = h[26];
        hms[27*(*np)+i] = h[27];
        hms[28*(*np)+i] = h[28];
    }

    for (int n=0; n<29; n++) {
        for (int i=0; i<(*np); i++) {
            hms[n*(*np)+i] *= 8.31451e+07 * T[i] * imw[n];
        }
    }
}


/*Returns gibbs in mass units (Eq 31.) */
void CKGMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict gms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    gibbs(gms, tc);
    for (int i = 0; i < 29; i++)
    {
        gms[i] *= RT*imw[i];
    }
}


/*Returns helmholtz in mass units (Eq 32.) */
void CKAMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict ams)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    helmholtz(ams, tc);
    for (int i = 0; i < 29; i++)
    {
        ams[i] *= RT*imw[i];
    }
}


/*Returns the entropies in mass units (Eq 28.) */
void CKSMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict sms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sms, tc);
    /*multiply by R/molecularweight */
    sms[0] *= 8.297452926019980e+05; /*NC7H16 */
    sms[1] *= 2.598381814318037e+06; /*O2 */
    sms[2] *= 1.889234139098090e+06; /*CO2 */
    sms[3] *= 4.615239012974499e+06; /*H2O */
    sms[4] *= 2.968349425484326e+06; /*CO */
    sms[5] *= 4.124383662212169e+07; /*H2 */
    sms[6] *= 4.888768810227566e+06; /*OH */
    sms[7] *= 2.444384405113783e+06; /*H2O2 */
    sms[8] *= 2.519031701678171e+06; /*HO2 */
    sms[9] *= 8.248767324424338e+07; /*H */
    sms[10] *= 5.196763628636074e+06; /*O */
    sms[11] *= 2.679121853578248e+06; /*CH3O */
    sms[12] *= 2.769058254894261e+06; /*CH2O */
    sms[13] *= 2.865242610581105e+06; /*HCO */
    sms[14] *= 5.927466067445207e+06; /*CH2 */
    sms[15] *= 5.530081023953346e+06; /*CH3 */
    sms[16] *= 5.182630712527496e+06; /*CH4 */
    sms[17] *= 3.074186734481467e+06; /*C2H3 */
    sms[18] *= 2.963733033722604e+06; /*C2H4 */
    sms[19] *= 2.860941121011349e+06; /*C2H5 */
    sms[20] *= 2.075238117344846e+06; /*C3H4 */
    sms[21] *= 2.024310196648431e+06; /*C3H5 */
    sms[22] *= 1.975822022481736e+06; /*C3H6 */
    sms[23] *= 1.929602378691293e+06; /*C3H7 */
    sms[24] *= 8.381765284643982e+05; /*C7H15-2 */
    sms[25] *= 6.337452856938147e+05; /*C7H15O2 */
    sms[26] *= 5.687552787396871e+05; /*C7KET12 */
    sms[27] *= 8.385453451838590e+05; /*C5H11CO */
    sms[28] *= 2.968047434442088e+06; /*N2 */
}


/*Returns the mean specific heat at CP (Eq. 33) */
void CKCPBL(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict cpbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[29]; /* temporary storage */
    cp_R(cpor, tc);

    /*perform dot product */
    for (id = 0; id < 29; ++id) {
        result += x[id]*cpor[id];
    }

    *cpbl = result * 8.31451e+07;
}


/*Returns the mean specific heat at CP (Eq. 34) */
void CKCPBS(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict cpbs)
{
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[29], tresult[29]; /* temporary storage */
    cp_R(cpor, tc);
    for (int i = 0; i < 29; i++)
    {
        tresult[i] = cpor[i]*y[i]*imw[i];

    }
    for (int i = 0; i < 29; i++)
    {
        result += tresult[i];
    }

    *cpbs = result * 8.31451e+07;
}


/*Returns the mean specific heat at CV (Eq. 35) */
void CKCVBL(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict cvbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[29]; /* temporary storage */
    cv_R(cvor, tc);

    /*perform dot product */
    for (id = 0; id < 29; ++id) {
        result += x[id]*cvor[id];
    }

    *cvbl = result * 8.31451e+07;
}


/*Returns the mean specific heat at CV (Eq. 36) */
void CKCVBS(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict cvbs)
{
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[29]; /* temporary storage */
    cv_R(cvor, tc);
    /*multiply by y/molecularweight */
    result += cvor[0]*y[0]*imw[0]; /*NC7H16 */
    result += cvor[1]*y[1]*imw[1]; /*O2 */
    result += cvor[2]*y[2]*imw[2]; /*CO2 */
    result += cvor[3]*y[3]*imw[3]; /*H2O */
    result += cvor[4]*y[4]*imw[4]; /*CO */
    result += cvor[5]*y[5]*imw[5]; /*H2 */
    result += cvor[6]*y[6]*imw[6]; /*OH */
    result += cvor[7]*y[7]*imw[7]; /*H2O2 */
    result += cvor[8]*y[8]*imw[8]; /*HO2 */
    result += cvor[9]*y[9]*imw[9]; /*H */
    result += cvor[10]*y[10]*imw[10]; /*O */
    result += cvor[11]*y[11]*imw[11]; /*CH3O */
    result += cvor[12]*y[12]*imw[12]; /*CH2O */
    result += cvor[13]*y[13]*imw[13]; /*HCO */
    result += cvor[14]*y[14]*imw[14]; /*CH2 */
    result += cvor[15]*y[15]*imw[15]; /*CH3 */
    result += cvor[16]*y[16]*imw[16]; /*CH4 */
    result += cvor[17]*y[17]*imw[17]; /*C2H3 */
    result += cvor[18]*y[18]*imw[18]; /*C2H4 */
    result += cvor[19]*y[19]*imw[19]; /*C2H5 */
    result += cvor[20]*y[20]*imw[20]; /*C3H4 */
    result += cvor[21]*y[21]*imw[21]; /*C3H5 */
    result += cvor[22]*y[22]*imw[22]; /*C3H6 */
    result += cvor[23]*y[23]*imw[23]; /*C3H7 */
    result += cvor[24]*y[24]*imw[24]; /*C7H15-2 */
    result += cvor[25]*y[25]*imw[25]; /*C7H15O2 */
    result += cvor[26]*y[26]*imw[26]; /*C7KET12 */
    result += cvor[27]*y[27]*imw[27]; /*C5H11CO */
    result += cvor[28]*y[28]*imw[28]; /*N2 */

    *cvbs = result * 8.31451e+07;
}


/*Returns the mean enthalpy of the mixture in molar units */
void CKHBML(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict hbml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[29]; /* temporary storage */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*perform dot product */
    for (id = 0; id < 29; ++id) {
        result += x[id]*hml[id];
    }

    *hbml = result * RT;
}


/*Returns mean enthalpy of mixture in mass units */
void CKHBMS(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict hbms)
{
    double result = 0;
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[29], tmp[29]; /* temporary storage */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);
    int id;
    for (id = 0; id < 29; ++id) {
        tmp[id] = y[id]*hml[id]*imw[id];
    }
    for (id = 0; id < 29; ++id) {
        result += tmp[id];
    }

    *hbms = result * RT;
}


/*get mean internal energy in molar units */
void CKUBML(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict ubml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double uml[29]; /* temporary energy array */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*perform dot product */
    for (id = 0; id < 29; ++id) {
        result += x[id]*uml[id];
    }

    *ubml = result * RT;
}


/*get mean internal energy in mass units */
void CKUBMS(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict ubms)
{
    double result = 0;
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double ums[29]; /* temporary energy array */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    /* DEBUG
    for (int n=0; n<29; n++) {
        printf(" IN CKUBMS nspec=%d ums=%.12f \n", n,  ums[n]);
    }
    DEBUG */

    /*perform dot product + scaling by wt */
    result += y[0]*ums[0]*imw[0]; /*NC7H16 */
    result += y[1]*ums[1]*imw[1]; /*O2 */
    result += y[2]*ums[2]*imw[2]; /*CO2 */
    result += y[3]*ums[3]*imw[3]; /*H2O */
    result += y[4]*ums[4]*imw[4]; /*CO */
    result += y[5]*ums[5]*imw[5]; /*H2 */
    result += y[6]*ums[6]*imw[6]; /*OH */
    result += y[7]*ums[7]*imw[7]; /*H2O2 */
    result += y[8]*ums[8]*imw[8]; /*HO2 */
    result += y[9]*ums[9]*imw[9]; /*H */
    result += y[10]*ums[10]*imw[10]; /*O */
    result += y[11]*ums[11]*imw[11]; /*CH3O */
    result += y[12]*ums[12]*imw[12]; /*CH2O */
    result += y[13]*ums[13]*imw[13]; /*HCO */
    result += y[14]*ums[14]*imw[14]; /*CH2 */
    result += y[15]*ums[15]*imw[15]; /*CH3 */
    result += y[16]*ums[16]*imw[16]; /*CH4 */
    result += y[17]*ums[17]*imw[17]; /*C2H3 */
    result += y[18]*ums[18]*imw[18]; /*C2H4 */
    result += y[19]*ums[19]*imw[19]; /*C2H5 */
    result += y[20]*ums[20]*imw[20]; /*C3H4 */
    result += y[21]*ums[21]*imw[21]; /*C3H5 */
    result += y[22]*ums[22]*imw[22]; /*C3H6 */
    result += y[23]*ums[23]*imw[23]; /*C3H7 */
    result += y[24]*ums[24]*imw[24]; /*C7H15-2 */
    result += y[25]*ums[25]*imw[25]; /*C7H15O2 */
    result += y[26]*ums[26]*imw[26]; /*C7KET12 */
    result += y[27]*ums[27]*imw[27]; /*C5H11CO */
    result += y[28]*ums[28]*imw[28]; /*N2 */

    //printf(" IN CKUBMS result=%.12f \n",result);

    *ubms = result * RT;
}

/*get mean internal energy in mass units */
void VCKUBMS(int * restrict np, double * restrict T, double * restrict y, double * restrict ubms)
{
    double ums[29*(*np)]; /* temporary energy array */
    double YOW[*np];

    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

    vspeciesInternalEnergy(*np, ums, T);
    /* DEBUG 
    for (int i=0; i<(*np); i++) {
            for (int n=0; n<29; n++) {
                printf(" IN VCKUBMS nspec=%d ipt=%d ums=%.12f \n", n, i, ums[n*(*np)+i]);
	    }
    }
    DEBUG */

    for (int n=0; n<29; n++) {
        for (int i=0; i<(*np); i++) {
            YOW[i] += y[i*(29)+n] * ums[n*(*np)+i] * imw[n];
        }
    }

    /* DEBUG  
    for (int n=0; n<(*np); n++) {
        printf(" IN VCKUBMS ipt=%d result=%.12f \n", n, YOW[n]);
    }
     DEBUG */

    for (int i=0; i<(*np); i++) {
        ubms[i] = YOW[i] * 8.31451e+07 * T[i] ; 
    }

    return;
}


/*get mixture entropy in molar units */
void CKSBML(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict sbml)
{
    int id; /*loop counter */
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double sor[29]; /* temporary storage */
    speciesEntropy(sor, tc);

    /*Compute Eq 42 */
    for (id = 0; id < 29; ++id) {
        result += x[id]*(sor[id]-log((x[id]+1e-100))-logPratio);
    }

    *sbml = result * 8.31451e+07;
}


/*get mixture entropy in mass units */
void CKSBMS(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict sbms)
{
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double sor[29]; /* temporary storage */
    double x[29]; /* need a ytx conversion */
    double YOW = 0; /*See Eq 4, 6 in CK Manual */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*NC7H16 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*CO2 */
    YOW += y[3]*imw[3]; /*H2O */
    YOW += y[4]*imw[4]; /*CO */
    YOW += y[5]*imw[5]; /*H2 */
    YOW += y[6]*imw[6]; /*OH */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*HO2 */
    YOW += y[9]*imw[9]; /*H */
    YOW += y[10]*imw[10]; /*O */
    YOW += y[11]*imw[11]; /*CH3O */
    YOW += y[12]*imw[12]; /*CH2O */
    YOW += y[13]*imw[13]; /*HCO */
    YOW += y[14]*imw[14]; /*CH2 */
    YOW += y[15]*imw[15]; /*CH3 */
    YOW += y[16]*imw[16]; /*CH4 */
    YOW += y[17]*imw[17]; /*C2H3 */
    YOW += y[18]*imw[18]; /*C2H4 */
    YOW += y[19]*imw[19]; /*C2H5 */
    YOW += y[20]*imw[20]; /*C3H4 */
    YOW += y[21]*imw[21]; /*C3H5 */
    YOW += y[22]*imw[22]; /*C3H6 */
    YOW += y[23]*imw[23]; /*C3H7 */
    YOW += y[24]*imw[24]; /*C7H15-2 */
    YOW += y[25]*imw[25]; /*C7H15O2 */
    YOW += y[26]*imw[26]; /*C7KET12 */
    YOW += y[27]*imw[27]; /*C5H11CO */
    YOW += y[28]*imw[28]; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(100.205570*YOW); 
    x[1] = y[1]/(31.998800*YOW); 
    x[2] = y[2]/(44.009950*YOW); 
    x[3] = y[3]/(18.015340*YOW); 
    x[4] = y[4]/(28.010550*YOW); 
    x[5] = y[5]/(2.015940*YOW); 
    x[6] = y[6]/(17.007370*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(33.006770*YOW); 
    x[9] = y[9]/(1.007970*YOW); 
    x[10] = y[10]/(15.999400*YOW); 
    x[11] = y[11]/(31.034460*YOW); 
    x[12] = y[12]/(30.026490*YOW); 
    x[13] = y[13]/(29.018520*YOW); 
    x[14] = y[14]/(14.027090*YOW); 
    x[15] = y[15]/(15.035060*YOW); 
    x[16] = y[16]/(16.043030*YOW); 
    x[17] = y[17]/(27.046210*YOW); 
    x[18] = y[18]/(28.054180*YOW); 
    x[19] = y[19]/(29.062150*YOW); 
    x[20] = y[20]/(40.065330*YOW); 
    x[21] = y[21]/(41.073300*YOW); 
    x[22] = y[22]/(42.081270*YOW); 
    x[23] = y[23]/(43.089240*YOW); 
    x[24] = y[24]/(99.197600*YOW); 
    x[25] = y[25]/(131.196400*YOW); 
    x[26] = y[26]/(146.187830*YOW); 
    x[27] = y[27]/(99.153970*YOW); 
    x[28] = y[28]/(28.013400*YOW); 
    speciesEntropy(sor, tc);
    /*Perform computation in Eq 42 and 43 */
    result += x[0]*(sor[0]-log((x[0]+1e-100))-logPratio);
    result += x[1]*(sor[1]-log((x[1]+1e-100))-logPratio);
    result += x[2]*(sor[2]-log((x[2]+1e-100))-logPratio);
    result += x[3]*(sor[3]-log((x[3]+1e-100))-logPratio);
    result += x[4]*(sor[4]-log((x[4]+1e-100))-logPratio);
    result += x[5]*(sor[5]-log((x[5]+1e-100))-logPratio);
    result += x[6]*(sor[6]-log((x[6]+1e-100))-logPratio);
    result += x[7]*(sor[7]-log((x[7]+1e-100))-logPratio);
    result += x[8]*(sor[8]-log((x[8]+1e-100))-logPratio);
    result += x[9]*(sor[9]-log((x[9]+1e-100))-logPratio);
    result += x[10]*(sor[10]-log((x[10]+1e-100))-logPratio);
    result += x[11]*(sor[11]-log((x[11]+1e-100))-logPratio);
    result += x[12]*(sor[12]-log((x[12]+1e-100))-logPratio);
    result += x[13]*(sor[13]-log((x[13]+1e-100))-logPratio);
    result += x[14]*(sor[14]-log((x[14]+1e-100))-logPratio);
    result += x[15]*(sor[15]-log((x[15]+1e-100))-logPratio);
    result += x[16]*(sor[16]-log((x[16]+1e-100))-logPratio);
    result += x[17]*(sor[17]-log((x[17]+1e-100))-logPratio);
    result += x[18]*(sor[18]-log((x[18]+1e-100))-logPratio);
    result += x[19]*(sor[19]-log((x[19]+1e-100))-logPratio);
    result += x[20]*(sor[20]-log((x[20]+1e-100))-logPratio);
    result += x[21]*(sor[21]-log((x[21]+1e-100))-logPratio);
    result += x[22]*(sor[22]-log((x[22]+1e-100))-logPratio);
    result += x[23]*(sor[23]-log((x[23]+1e-100))-logPratio);
    result += x[24]*(sor[24]-log((x[24]+1e-100))-logPratio);
    result += x[25]*(sor[25]-log((x[25]+1e-100))-logPratio);
    result += x[26]*(sor[26]-log((x[26]+1e-100))-logPratio);
    result += x[27]*(sor[27]-log((x[27]+1e-100))-logPratio);
    result += x[28]*(sor[28]-log((x[28]+1e-100))-logPratio);
    /*Scale by R/W */
    *sbms = result * 8.31451e+07 * YOW;
}


/*Returns mean gibbs free energy in molar units */
void CKGBML(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict gbml)
{
    int id; /*loop counter */
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    double gort[29]; /* temporary storage */
    /*Compute g/RT */
    gibbs(gort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 29; ++id) {
        result += x[id]*(gort[id]+log((x[id]+1e-100))+logPratio);
    }

    *gbml = result * RT;
}


/*Returns mixture gibbs free energy in mass units */
void CKGBMS(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict gbms)
{
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    double gort[29]; /* temporary storage */
    double x[29]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*NC7H16 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*CO2 */
    YOW += y[3]*imw[3]; /*H2O */
    YOW += y[4]*imw[4]; /*CO */
    YOW += y[5]*imw[5]; /*H2 */
    YOW += y[6]*imw[6]; /*OH */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*HO2 */
    YOW += y[9]*imw[9]; /*H */
    YOW += y[10]*imw[10]; /*O */
    YOW += y[11]*imw[11]; /*CH3O */
    YOW += y[12]*imw[12]; /*CH2O */
    YOW += y[13]*imw[13]; /*HCO */
    YOW += y[14]*imw[14]; /*CH2 */
    YOW += y[15]*imw[15]; /*CH3 */
    YOW += y[16]*imw[16]; /*CH4 */
    YOW += y[17]*imw[17]; /*C2H3 */
    YOW += y[18]*imw[18]; /*C2H4 */
    YOW += y[19]*imw[19]; /*C2H5 */
    YOW += y[20]*imw[20]; /*C3H4 */
    YOW += y[21]*imw[21]; /*C3H5 */
    YOW += y[22]*imw[22]; /*C3H6 */
    YOW += y[23]*imw[23]; /*C3H7 */
    YOW += y[24]*imw[24]; /*C7H15-2 */
    YOW += y[25]*imw[25]; /*C7H15O2 */
    YOW += y[26]*imw[26]; /*C7KET12 */
    YOW += y[27]*imw[27]; /*C5H11CO */
    YOW += y[28]*imw[28]; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(100.205570*YOW); 
    x[1] = y[1]/(31.998800*YOW); 
    x[2] = y[2]/(44.009950*YOW); 
    x[3] = y[3]/(18.015340*YOW); 
    x[4] = y[4]/(28.010550*YOW); 
    x[5] = y[5]/(2.015940*YOW); 
    x[6] = y[6]/(17.007370*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(33.006770*YOW); 
    x[9] = y[9]/(1.007970*YOW); 
    x[10] = y[10]/(15.999400*YOW); 
    x[11] = y[11]/(31.034460*YOW); 
    x[12] = y[12]/(30.026490*YOW); 
    x[13] = y[13]/(29.018520*YOW); 
    x[14] = y[14]/(14.027090*YOW); 
    x[15] = y[15]/(15.035060*YOW); 
    x[16] = y[16]/(16.043030*YOW); 
    x[17] = y[17]/(27.046210*YOW); 
    x[18] = y[18]/(28.054180*YOW); 
    x[19] = y[19]/(29.062150*YOW); 
    x[20] = y[20]/(40.065330*YOW); 
    x[21] = y[21]/(41.073300*YOW); 
    x[22] = y[22]/(42.081270*YOW); 
    x[23] = y[23]/(43.089240*YOW); 
    x[24] = y[24]/(99.197600*YOW); 
    x[25] = y[25]/(131.196400*YOW); 
    x[26] = y[26]/(146.187830*YOW); 
    x[27] = y[27]/(99.153970*YOW); 
    x[28] = y[28]/(28.013400*YOW); 
    gibbs(gort, tc);
    /*Perform computation in Eq 44 */
    result += x[0]*(gort[0]+log((x[0]+1e-100))+logPratio);
    result += x[1]*(gort[1]+log((x[1]+1e-100))+logPratio);
    result += x[2]*(gort[2]+log((x[2]+1e-100))+logPratio);
    result += x[3]*(gort[3]+log((x[3]+1e-100))+logPratio);
    result += x[4]*(gort[4]+log((x[4]+1e-100))+logPratio);
    result += x[5]*(gort[5]+log((x[5]+1e-100))+logPratio);
    result += x[6]*(gort[6]+log((x[6]+1e-100))+logPratio);
    result += x[7]*(gort[7]+log((x[7]+1e-100))+logPratio);
    result += x[8]*(gort[8]+log((x[8]+1e-100))+logPratio);
    result += x[9]*(gort[9]+log((x[9]+1e-100))+logPratio);
    result += x[10]*(gort[10]+log((x[10]+1e-100))+logPratio);
    result += x[11]*(gort[11]+log((x[11]+1e-100))+logPratio);
    result += x[12]*(gort[12]+log((x[12]+1e-100))+logPratio);
    result += x[13]*(gort[13]+log((x[13]+1e-100))+logPratio);
    result += x[14]*(gort[14]+log((x[14]+1e-100))+logPratio);
    result += x[15]*(gort[15]+log((x[15]+1e-100))+logPratio);
    result += x[16]*(gort[16]+log((x[16]+1e-100))+logPratio);
    result += x[17]*(gort[17]+log((x[17]+1e-100))+logPratio);
    result += x[18]*(gort[18]+log((x[18]+1e-100))+logPratio);
    result += x[19]*(gort[19]+log((x[19]+1e-100))+logPratio);
    result += x[20]*(gort[20]+log((x[20]+1e-100))+logPratio);
    result += x[21]*(gort[21]+log((x[21]+1e-100))+logPratio);
    result += x[22]*(gort[22]+log((x[22]+1e-100))+logPratio);
    result += x[23]*(gort[23]+log((x[23]+1e-100))+logPratio);
    result += x[24]*(gort[24]+log((x[24]+1e-100))+logPratio);
    result += x[25]*(gort[25]+log((x[25]+1e-100))+logPratio);
    result += x[26]*(gort[26]+log((x[26]+1e-100))+logPratio);
    result += x[27]*(gort[27]+log((x[27]+1e-100))+logPratio);
    result += x[28]*(gort[28]+log((x[28]+1e-100))+logPratio);
    /*Scale by RT/W */
    *gbms = result * RT * YOW;
}


/*Returns mean helmholtz free energy in molar units */
void CKABML(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict abml)
{
    int id; /*loop counter */
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    double aort[29]; /* temporary storage */
    /*Compute g/RT */
    helmholtz(aort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 29; ++id) {
        result += x[id]*(aort[id]+log((x[id]+1e-100))+logPratio);
    }

    *abml = result * RT;
}


/*Returns mixture helmholtz free energy in mass units */
void CKABMS(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict abms)
{
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    double aort[29]; /* temporary storage */
    double x[29]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*NC7H16 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*CO2 */
    YOW += y[3]*imw[3]; /*H2O */
    YOW += y[4]*imw[4]; /*CO */
    YOW += y[5]*imw[5]; /*H2 */
    YOW += y[6]*imw[6]; /*OH */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*HO2 */
    YOW += y[9]*imw[9]; /*H */
    YOW += y[10]*imw[10]; /*O */
    YOW += y[11]*imw[11]; /*CH3O */
    YOW += y[12]*imw[12]; /*CH2O */
    YOW += y[13]*imw[13]; /*HCO */
    YOW += y[14]*imw[14]; /*CH2 */
    YOW += y[15]*imw[15]; /*CH3 */
    YOW += y[16]*imw[16]; /*CH4 */
    YOW += y[17]*imw[17]; /*C2H3 */
    YOW += y[18]*imw[18]; /*C2H4 */
    YOW += y[19]*imw[19]; /*C2H5 */
    YOW += y[20]*imw[20]; /*C3H4 */
    YOW += y[21]*imw[21]; /*C3H5 */
    YOW += y[22]*imw[22]; /*C3H6 */
    YOW += y[23]*imw[23]; /*C3H7 */
    YOW += y[24]*imw[24]; /*C7H15-2 */
    YOW += y[25]*imw[25]; /*C7H15O2 */
    YOW += y[26]*imw[26]; /*C7KET12 */
    YOW += y[27]*imw[27]; /*C5H11CO */
    YOW += y[28]*imw[28]; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(100.205570*YOW); 
    x[1] = y[1]/(31.998800*YOW); 
    x[2] = y[2]/(44.009950*YOW); 
    x[3] = y[3]/(18.015340*YOW); 
    x[4] = y[4]/(28.010550*YOW); 
    x[5] = y[5]/(2.015940*YOW); 
    x[6] = y[6]/(17.007370*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(33.006770*YOW); 
    x[9] = y[9]/(1.007970*YOW); 
    x[10] = y[10]/(15.999400*YOW); 
    x[11] = y[11]/(31.034460*YOW); 
    x[12] = y[12]/(30.026490*YOW); 
    x[13] = y[13]/(29.018520*YOW); 
    x[14] = y[14]/(14.027090*YOW); 
    x[15] = y[15]/(15.035060*YOW); 
    x[16] = y[16]/(16.043030*YOW); 
    x[17] = y[17]/(27.046210*YOW); 
    x[18] = y[18]/(28.054180*YOW); 
    x[19] = y[19]/(29.062150*YOW); 
    x[20] = y[20]/(40.065330*YOW); 
    x[21] = y[21]/(41.073300*YOW); 
    x[22] = y[22]/(42.081270*YOW); 
    x[23] = y[23]/(43.089240*YOW); 
    x[24] = y[24]/(99.197600*YOW); 
    x[25] = y[25]/(131.196400*YOW); 
    x[26] = y[26]/(146.187830*YOW); 
    x[27] = y[27]/(99.153970*YOW); 
    x[28] = y[28]/(28.013400*YOW); 
    helmholtz(aort, tc);
    /*Perform computation in Eq 44 */
    result += x[0]*(aort[0]+log((x[0]+1e-100))+logPratio);
    result += x[1]*(aort[1]+log((x[1]+1e-100))+logPratio);
    result += x[2]*(aort[2]+log((x[2]+1e-100))+logPratio);
    result += x[3]*(aort[3]+log((x[3]+1e-100))+logPratio);
    result += x[4]*(aort[4]+log((x[4]+1e-100))+logPratio);
    result += x[5]*(aort[5]+log((x[5]+1e-100))+logPratio);
    result += x[6]*(aort[6]+log((x[6]+1e-100))+logPratio);
    result += x[7]*(aort[7]+log((x[7]+1e-100))+logPratio);
    result += x[8]*(aort[8]+log((x[8]+1e-100))+logPratio);
    result += x[9]*(aort[9]+log((x[9]+1e-100))+logPratio);
    result += x[10]*(aort[10]+log((x[10]+1e-100))+logPratio);
    result += x[11]*(aort[11]+log((x[11]+1e-100))+logPratio);
    result += x[12]*(aort[12]+log((x[12]+1e-100))+logPratio);
    result += x[13]*(aort[13]+log((x[13]+1e-100))+logPratio);
    result += x[14]*(aort[14]+log((x[14]+1e-100))+logPratio);
    result += x[15]*(aort[15]+log((x[15]+1e-100))+logPratio);
    result += x[16]*(aort[16]+log((x[16]+1e-100))+logPratio);
    result += x[17]*(aort[17]+log((x[17]+1e-100))+logPratio);
    result += x[18]*(aort[18]+log((x[18]+1e-100))+logPratio);
    result += x[19]*(aort[19]+log((x[19]+1e-100))+logPratio);
    result += x[20]*(aort[20]+log((x[20]+1e-100))+logPratio);
    result += x[21]*(aort[21]+log((x[21]+1e-100))+logPratio);
    result += x[22]*(aort[22]+log((x[22]+1e-100))+logPratio);
    result += x[23]*(aort[23]+log((x[23]+1e-100))+logPratio);
    result += x[24]*(aort[24]+log((x[24]+1e-100))+logPratio);
    result += x[25]*(aort[25]+log((x[25]+1e-100))+logPratio);
    result += x[26]*(aort[26]+log((x[26]+1e-100))+logPratio);
    result += x[27]*(aort[27]+log((x[27]+1e-100))+logPratio);
    result += x[28]*(aort[28]+log((x[28]+1e-100))+logPratio);
    /*Scale by RT/W */
    *abms = result * RT * YOW;
}


/*compute the production rate for each species */
void CKWC(double * restrict T, double * restrict C, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 29; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    productionRate(wdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 29; ++id) {
        C[id] *= 1.0e-6;
        wdot[id] *= 1.0e-6;
    }
}

/*compute the production rate for each species */
void VCKWC(int * restrict np, double * restrict T, double * restrict C, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    double C_tmp[29*(*np)]; /* temporary energy array */
    double wdot_tmp[29*(*np)]; /* temporary energy array */

    /*convert to SI */
    for (int n=0; n<29; n++) {
        for (int i=0; i<(*np); i++) {
            C_tmp[n*(*np)+i] = C[i*(29)+n] * 1.0e6;
        }
    }

    /*convert to chemkin units */
    vproductionRate(*np, wdot_tmp, C_tmp, T);

    /*convert to chemkin units */
    for (int n=0; n<29; n++) {
        for (int i=0; i<(*np); i++) {
        //C[i*(29)+n] = C_tmp[n*(*np)+i]*1.0e-6;
        wdot[i*(29)+n] = wdot_tmp[n*(*np)+i]*1.0e-6;
        }
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mass fractions */
void CKWYP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */
    double c[29]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*NC7H16 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*CO2 */
    YOW += y[3]*imw[3]; /*H2O */
    YOW += y[4]*imw[4]; /*CO */
    YOW += y[5]*imw[5]; /*H2 */
    YOW += y[6]*imw[6]; /*OH */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*HO2 */
    YOW += y[9]*imw[9]; /*H */
    YOW += y[10]*imw[10]; /*O */
    YOW += y[11]*imw[11]; /*CH3O */
    YOW += y[12]*imw[12]; /*CH2O */
    YOW += y[13]*imw[13]; /*HCO */
    YOW += y[14]*imw[14]; /*CH2 */
    YOW += y[15]*imw[15]; /*CH3 */
    YOW += y[16]*imw[16]; /*CH4 */
    YOW += y[17]*imw[17]; /*C2H3 */
    YOW += y[18]*imw[18]; /*C2H4 */
    YOW += y[19]*imw[19]; /*C2H5 */
    YOW += y[20]*imw[20]; /*C3H4 */
    YOW += y[21]*imw[21]; /*C3H5 */
    YOW += y[22]*imw[22]; /*C3H6 */
    YOW += y[23]*imw[23]; /*C3H7 */
    YOW += y[24]*imw[24]; /*C7H15-2 */
    YOW += y[25]*imw[25]; /*C7H15O2 */
    YOW += y[26]*imw[26]; /*C7KET12 */
    YOW += y[27]*imw[27]; /*C5H11CO */
    YOW += y[28]*imw[28]; /*N2 */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31451e+07 * (*T)); 
    /*multiply by 1e6 so c goes to SI */
    PWORT *= 1e6; 
    /*Now compute conversion (and go to SI) */
    c[0] = PWORT * y[0]*imw[0]; 
    c[1] = PWORT * y[1]*imw[1]; 
    c[2] = PWORT * y[2]*imw[2]; 
    c[3] = PWORT * y[3]*imw[3]; 
    c[4] = PWORT * y[4]*imw[4]; 
    c[5] = PWORT * y[5]*imw[5]; 
    c[6] = PWORT * y[6]*imw[6]; 
    c[7] = PWORT * y[7]*imw[7]; 
    c[8] = PWORT * y[8]*imw[8]; 
    c[9] = PWORT * y[9]*imw[9]; 
    c[10] = PWORT * y[10]*imw[10]; 
    c[11] = PWORT * y[11]*imw[11]; 
    c[12] = PWORT * y[12]*imw[12]; 
    c[13] = PWORT * y[13]*imw[13]; 
    c[14] = PWORT * y[14]*imw[14]; 
    c[15] = PWORT * y[15]*imw[15]; 
    c[16] = PWORT * y[16]*imw[16]; 
    c[17] = PWORT * y[17]*imw[17]; 
    c[18] = PWORT * y[18]*imw[18]; 
    c[19] = PWORT * y[19]*imw[19]; 
    c[20] = PWORT * y[20]*imw[20]; 
    c[21] = PWORT * y[21]*imw[21]; 
    c[22] = PWORT * y[22]*imw[22]; 
    c[23] = PWORT * y[23]*imw[23]; 
    c[24] = PWORT * y[24]*imw[24]; 
    c[25] = PWORT * y[25]*imw[25]; 
    c[26] = PWORT * y[26]*imw[26]; 
    c[27] = PWORT * y[27]*imw[27]; 
    c[28] = PWORT * y[28]*imw[28]; 

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 29; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mole fractions */
void CKWXP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */
    double c[29]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 29; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 29; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void CKWYR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */
    double c[29]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    c[0] = 1e6 * (*rho) * y[0]*imw[0]; 
    c[1] = 1e6 * (*rho) * y[1]*imw[1]; 
    c[2] = 1e6 * (*rho) * y[2]*imw[2]; 
    c[3] = 1e6 * (*rho) * y[3]*imw[3]; 
    c[4] = 1e6 * (*rho) * y[4]*imw[4]; 
    c[5] = 1e6 * (*rho) * y[5]*imw[5]; 
    c[6] = 1e6 * (*rho) * y[6]*imw[6]; 
    c[7] = 1e6 * (*rho) * y[7]*imw[7]; 
    c[8] = 1e6 * (*rho) * y[8]*imw[8]; 
    c[9] = 1e6 * (*rho) * y[9]*imw[9]; 
    c[10] = 1e6 * (*rho) * y[10]*imw[10]; 
    c[11] = 1e6 * (*rho) * y[11]*imw[11]; 
    c[12] = 1e6 * (*rho) * y[12]*imw[12]; 
    c[13] = 1e6 * (*rho) * y[13]*imw[13]; 
    c[14] = 1e6 * (*rho) * y[14]*imw[14]; 
    c[15] = 1e6 * (*rho) * y[15]*imw[15]; 
    c[16] = 1e6 * (*rho) * y[16]*imw[16]; 
    c[17] = 1e6 * (*rho) * y[17]*imw[17]; 
    c[18] = 1e6 * (*rho) * y[18]*imw[18]; 
    c[19] = 1e6 * (*rho) * y[19]*imw[19]; 
    c[20] = 1e6 * (*rho) * y[20]*imw[20]; 
    c[21] = 1e6 * (*rho) * y[21]*imw[21]; 
    c[22] = 1e6 * (*rho) * y[22]*imw[22]; 
    c[23] = 1e6 * (*rho) * y[23]*imw[23]; 
    c[24] = 1e6 * (*rho) * y[24]*imw[24]; 
    c[25] = 1e6 * (*rho) * y[25]*imw[25]; 
    c[26] = 1e6 * (*rho) * y[26]*imw[26]; 
    c[27] = 1e6 * (*rho) * y[27]*imw[27]; 
    c[28] = 1e6 * (*rho) * y[28]*imw[28]; 

    /*call productionRate */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 29; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void VCKWYR(int * restrict np, double * restrict rho, double * restrict T,
	    double * restrict y, int * restrict iwrk, double * restrict rwrk,
	    double * restrict wdot)
{
    double c[29*(*np)]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    for (int n=0; n<29; n++) {
        for (int i=0; i<(*np); i++) {
            c[n*(*np)+i] = 1.0e6 * rho[i] * y[n*(*np)+i] * imw[n];
        }
    }

    /*call productionRate */
    vproductionRate(*np, wdot, c, T);

    /*convert to chemkin units */
    for (int i=0; i<29*(*np); i++) {
        wdot[i] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mole fractions */
void CKWXR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */
    double c[29]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*100.205570; /*NC7H16 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*44.009950; /*CO2 */
    XW += x[3]*18.015340; /*H2O */
    XW += x[4]*28.010550; /*CO */
    XW += x[5]*2.015940; /*H2 */
    XW += x[6]*17.007370; /*OH */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*33.006770; /*HO2 */
    XW += x[9]*1.007970; /*H */
    XW += x[10]*15.999400; /*O */
    XW += x[11]*31.034460; /*CH3O */
    XW += x[12]*30.026490; /*CH2O */
    XW += x[13]*29.018520; /*HCO */
    XW += x[14]*14.027090; /*CH2 */
    XW += x[15]*15.035060; /*CH3 */
    XW += x[16]*16.043030; /*CH4 */
    XW += x[17]*27.046210; /*C2H3 */
    XW += x[18]*28.054180; /*C2H4 */
    XW += x[19]*29.062150; /*C2H5 */
    XW += x[20]*40.065330; /*C3H4 */
    XW += x[21]*41.073300; /*C3H5 */
    XW += x[22]*42.081270; /*C3H6 */
    XW += x[23]*43.089240; /*C3H7 */
    XW += x[24]*99.197600; /*C7H15-2 */
    XW += x[25]*131.196400; /*C7H15O2 */
    XW += x[26]*146.187830; /*C7KET12 */
    XW += x[27]*99.153970; /*C5H11CO */
    XW += x[28]*28.013400; /*N2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 29; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 29; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the rate of progress for each reaction */
void CKQC(double * restrict T, double * restrict C, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 29; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    progressRate(qdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 29; ++id) {
        C[id] *= 1.0e-6;
    }

    for (id = 0; id < 52; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKKFKR(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict q_f, double * restrict q_r)
{
    int id; /*loop counter */
    double c[29]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 29; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRateFR(q_f, q_r, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 52; ++id) {
        q_f[id] *= 1.0e-6;
        q_r[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mass fractions */
void CKQYP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */
    double c[29]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*NC7H16 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*CO2 */
    YOW += y[3]*imw[3]; /*H2O */
    YOW += y[4]*imw[4]; /*CO */
    YOW += y[5]*imw[5]; /*H2 */
    YOW += y[6]*imw[6]; /*OH */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*HO2 */
    YOW += y[9]*imw[9]; /*H */
    YOW += y[10]*imw[10]; /*O */
    YOW += y[11]*imw[11]; /*CH3O */
    YOW += y[12]*imw[12]; /*CH2O */
    YOW += y[13]*imw[13]; /*HCO */
    YOW += y[14]*imw[14]; /*CH2 */
    YOW += y[15]*imw[15]; /*CH3 */
    YOW += y[16]*imw[16]; /*CH4 */
    YOW += y[17]*imw[17]; /*C2H3 */
    YOW += y[18]*imw[18]; /*C2H4 */
    YOW += y[19]*imw[19]; /*C2H5 */
    YOW += y[20]*imw[20]; /*C3H4 */
    YOW += y[21]*imw[21]; /*C3H5 */
    YOW += y[22]*imw[22]; /*C3H6 */
    YOW += y[23]*imw[23]; /*C3H7 */
    YOW += y[24]*imw[24]; /*C7H15-2 */
    YOW += y[25]*imw[25]; /*C7H15O2 */
    YOW += y[26]*imw[26]; /*C7KET12 */
    YOW += y[27]*imw[27]; /*C5H11CO */
    YOW += y[28]*imw[28]; /*N2 */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31451e+07 * (*T)); 
    /*multiply by 1e6 so c goes to SI */
    PWORT *= 1e6; 
    /*Now compute conversion (and go to SI) */
    c[0] = PWORT * y[0]*imw[0]; 
    c[1] = PWORT * y[1]*imw[1]; 
    c[2] = PWORT * y[2]*imw[2]; 
    c[3] = PWORT * y[3]*imw[3]; 
    c[4] = PWORT * y[4]*imw[4]; 
    c[5] = PWORT * y[5]*imw[5]; 
    c[6] = PWORT * y[6]*imw[6]; 
    c[7] = PWORT * y[7]*imw[7]; 
    c[8] = PWORT * y[8]*imw[8]; 
    c[9] = PWORT * y[9]*imw[9]; 
    c[10] = PWORT * y[10]*imw[10]; 
    c[11] = PWORT * y[11]*imw[11]; 
    c[12] = PWORT * y[12]*imw[12]; 
    c[13] = PWORT * y[13]*imw[13]; 
    c[14] = PWORT * y[14]*imw[14]; 
    c[15] = PWORT * y[15]*imw[15]; 
    c[16] = PWORT * y[16]*imw[16]; 
    c[17] = PWORT * y[17]*imw[17]; 
    c[18] = PWORT * y[18]*imw[18]; 
    c[19] = PWORT * y[19]*imw[19]; 
    c[20] = PWORT * y[20]*imw[20]; 
    c[21] = PWORT * y[21]*imw[21]; 
    c[22] = PWORT * y[22]*imw[22]; 
    c[23] = PWORT * y[23]*imw[23]; 
    c[24] = PWORT * y[24]*imw[24]; 
    c[25] = PWORT * y[25]*imw[25]; 
    c[26] = PWORT * y[26]*imw[26]; 
    c[27] = PWORT * y[27]*imw[27]; 
    c[28] = PWORT * y[28]*imw[28]; 

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 52; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKQXP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */
    double c[29]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 29; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 52; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mass fractions */
void CKQYR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */
    double c[29]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    c[0] = 1e6 * (*rho) * y[0]*imw[0]; 
    c[1] = 1e6 * (*rho) * y[1]*imw[1]; 
    c[2] = 1e6 * (*rho) * y[2]*imw[2]; 
    c[3] = 1e6 * (*rho) * y[3]*imw[3]; 
    c[4] = 1e6 * (*rho) * y[4]*imw[4]; 
    c[5] = 1e6 * (*rho) * y[5]*imw[5]; 
    c[6] = 1e6 * (*rho) * y[6]*imw[6]; 
    c[7] = 1e6 * (*rho) * y[7]*imw[7]; 
    c[8] = 1e6 * (*rho) * y[8]*imw[8]; 
    c[9] = 1e6 * (*rho) * y[9]*imw[9]; 
    c[10] = 1e6 * (*rho) * y[10]*imw[10]; 
    c[11] = 1e6 * (*rho) * y[11]*imw[11]; 
    c[12] = 1e6 * (*rho) * y[12]*imw[12]; 
    c[13] = 1e6 * (*rho) * y[13]*imw[13]; 
    c[14] = 1e6 * (*rho) * y[14]*imw[14]; 
    c[15] = 1e6 * (*rho) * y[15]*imw[15]; 
    c[16] = 1e6 * (*rho) * y[16]*imw[16]; 
    c[17] = 1e6 * (*rho) * y[17]*imw[17]; 
    c[18] = 1e6 * (*rho) * y[18]*imw[18]; 
    c[19] = 1e6 * (*rho) * y[19]*imw[19]; 
    c[20] = 1e6 * (*rho) * y[20]*imw[20]; 
    c[21] = 1e6 * (*rho) * y[21]*imw[21]; 
    c[22] = 1e6 * (*rho) * y[22]*imw[22]; 
    c[23] = 1e6 * (*rho) * y[23]*imw[23]; 
    c[24] = 1e6 * (*rho) * y[24]*imw[24]; 
    c[25] = 1e6 * (*rho) * y[25]*imw[25]; 
    c[26] = 1e6 * (*rho) * y[26]*imw[26]; 
    c[27] = 1e6 * (*rho) * y[27]*imw[27]; 
    c[28] = 1e6 * (*rho) * y[28]*imw[28]; 

    /*call progressRate */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 52; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mole fractions */
void CKQXR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */
    double c[29]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*100.205570; /*NC7H16 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*44.009950; /*CO2 */
    XW += x[3]*18.015340; /*H2O */
    XW += x[4]*28.010550; /*CO */
    XW += x[5]*2.015940; /*H2 */
    XW += x[6]*17.007370; /*OH */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*33.006770; /*HO2 */
    XW += x[9]*1.007970; /*H */
    XW += x[10]*15.999400; /*O */
    XW += x[11]*31.034460; /*CH3O */
    XW += x[12]*30.026490; /*CH2O */
    XW += x[13]*29.018520; /*HCO */
    XW += x[14]*14.027090; /*CH2 */
    XW += x[15]*15.035060; /*CH3 */
    XW += x[16]*16.043030; /*CH4 */
    XW += x[17]*27.046210; /*C2H3 */
    XW += x[18]*28.054180; /*C2H4 */
    XW += x[19]*29.062150; /*C2H5 */
    XW += x[20]*40.065330; /*C3H4 */
    XW += x[21]*41.073300; /*C3H5 */
    XW += x[22]*42.081270; /*C3H6 */
    XW += x[23]*43.089240; /*C3H7 */
    XW += x[24]*99.197600; /*C7H15-2 */
    XW += x[25]*131.196400; /*C7H15O2 */
    XW += x[26]*146.187830; /*C7KET12 */
    XW += x[27]*99.153970; /*C5H11CO */
    XW += x[28]*28.013400; /*N2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 29; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 52; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the stoichiometric coefficients */
/*of the reaction mechanism. (Eq 50) */
void CKNU(int * kdim, int * iwrk, double * restrict rwrk, int * nuki)
{
    int id; /*loop counter */
    int kd = (*kdim); 
    /*Zero nuki */
    for (id = 0; id < 29 * kd; ++ id) {
         nuki[id] = 0; 
    }

    /*reaction 1: O2 + H + M <=> HO2 + M */
    nuki[ 1 * kd + 0 ] += -1 ;
    nuki[ 9 * kd + 0 ] += -1 ;
    nuki[ 8 * kd + 0 ] += +1 ;

    /*reaction 2: H2O2 + M <=> 2 OH + M */
    nuki[ 7 * kd + 1 ] += -1 ;
    nuki[ 6 * kd + 1 ] += +2 ;

    /*reaction 3: HCO + M <=> CO + H + M */
    nuki[ 13 * kd + 2 ] += -1 ;
    nuki[ 4 * kd + 2 ] += +1 ;
    nuki[ 9 * kd + 2 ] += +1 ;

    /*reaction 4: NC7H16 + H <=> H2 + C7H15-2 */
    nuki[ 0 * kd + 3 ] += -1 ;
    nuki[ 9 * kd + 3 ] += -1 ;
    nuki[ 5 * kd + 3 ] += +1 ;
    nuki[ 24 * kd + 3 ] += +1 ;

    /*reaction 5: NC7H16 + OH <=> H2O + C7H15-2 */
    nuki[ 0 * kd + 4 ] += -1 ;
    nuki[ 6 * kd + 4 ] += -1 ;
    nuki[ 3 * kd + 4 ] += +1 ;
    nuki[ 24 * kd + 4 ] += +1 ;

    /*reaction 6: NC7H16 + HO2 <=> H2O2 + C7H15-2 */
    nuki[ 0 * kd + 5 ] += -1 ;
    nuki[ 8 * kd + 5 ] += -1 ;
    nuki[ 7 * kd + 5 ] += +1 ;
    nuki[ 24 * kd + 5 ] += +1 ;

    /*reaction 7: NC7H16 + O2 <=> HO2 + C7H15-2 */
    nuki[ 0 * kd + 6 ] += -1 ;
    nuki[ 1 * kd + 6 ] += -1 ;
    nuki[ 8 * kd + 6 ] += +1 ;
    nuki[ 24 * kd + 6 ] += +1 ;

    /*reaction 8: O2 + C7H15-2 <=> C7H15O2 */
    nuki[ 1 * kd + 7 ] += -1 ;
    nuki[ 24 * kd + 7 ] += -1 ;
    nuki[ 25 * kd + 7 ] += +1 ;

    /*reaction 9: O2 + C7H15O2 <=> OH + C7KET12 */
    nuki[ 1 * kd + 8 ] += -1 ;
    nuki[ 25 * kd + 8 ] += -1 ;
    nuki[ 6 * kd + 8 ] += +1 ;
    nuki[ 26 * kd + 8 ] += +1 ;

    /*reaction 10: C7KET12 <=> OH + CH2O + C5H11CO */
    nuki[ 26 * kd + 9 ] += -1 ;
    nuki[ 6 * kd + 9 ] += +1 ;
    nuki[ 12 * kd + 9 ] += +1 ;
    nuki[ 27 * kd + 9 ] += +1 ;

    /*reaction 11: C5H11CO <=> CO + C2H4 + C3H7 */
    nuki[ 27 * kd + 10 ] += -1 ;
    nuki[ 4 * kd + 10 ] += +1 ;
    nuki[ 18 * kd + 10 ] += +1 ;
    nuki[ 23 * kd + 10 ] += +1 ;

    /*reaction 12: C7H15-2 <=> C2H4 + C2H5 + C3H6 */
    nuki[ 24 * kd + 11 ] += -1 ;
    nuki[ 18 * kd + 11 ] += +1 ;
    nuki[ 19 * kd + 11 ] += +1 ;
    nuki[ 22 * kd + 11 ] += +1 ;

    /*reaction 13: C3H7 <=> CH3 + C2H4 */
    nuki[ 23 * kd + 12 ] += -1 ;
    nuki[ 15 * kd + 12 ] += +1 ;
    nuki[ 18 * kd + 12 ] += +1 ;

    /*reaction 14: C3H7 <=> H + C3H6 */
    nuki[ 23 * kd + 13 ] += -1 ;
    nuki[ 9 * kd + 13 ] += +1 ;
    nuki[ 22 * kd + 13 ] += +1 ;

    /*reaction 15: CH3 + C3H6 <=> CH4 + C3H5 */
    nuki[ 15 * kd + 14 ] += -1 ;
    nuki[ 22 * kd + 14 ] += -1 ;
    nuki[ 16 * kd + 14 ] += +1 ;
    nuki[ 21 * kd + 14 ] += +1 ;

    /*reaction 16: O2 + C3H5 <=> HO2 + C3H4 */
    nuki[ 1 * kd + 15 ] += -1 ;
    nuki[ 21 * kd + 15 ] += -1 ;
    nuki[ 8 * kd + 15 ] += +1 ;
    nuki[ 20 * kd + 15 ] += +1 ;

    /*reaction 17: OH + C3H4 <=> CH2O + C2H3 */
    nuki[ 6 * kd + 16 ] += -1 ;
    nuki[ 20 * kd + 16 ] += -1 ;
    nuki[ 12 * kd + 16 ] += +1 ;
    nuki[ 17 * kd + 16 ] += +1 ;

    /*reaction 18: OH + C3H4 <=> HCO + C2H4 */
    nuki[ 6 * kd + 17 ] += -1 ;
    nuki[ 20 * kd + 17 ] += -1 ;
    nuki[ 13 * kd + 17 ] += +1 ;
    nuki[ 18 * kd + 17 ] += +1 ;

    /*reaction 19: HO2 + CH3 <=> OH + CH3O */
    nuki[ 8 * kd + 18 ] += -1 ;
    nuki[ 15 * kd + 18 ] += -1 ;
    nuki[ 6 * kd + 18 ] += +1 ;
    nuki[ 11 * kd + 18 ] += +1 ;

    /*reaction 20: OH + CH3 <=> H2O + CH2 */
    nuki[ 6 * kd + 19 ] += -1 ;
    nuki[ 15 * kd + 19 ] += -1 ;
    nuki[ 3 * kd + 19 ] += +1 ;
    nuki[ 14 * kd + 19 ] += +1 ;

    /*reaction 21: OH + CH2 <=> H + CH2O */
    nuki[ 6 * kd + 20 ] += -1 ;
    nuki[ 14 * kd + 20 ] += -1 ;
    nuki[ 9 * kd + 20 ] += +1 ;
    nuki[ 12 * kd + 20 ] += +1 ;

    /*reaction 22: O2 + CH2 <=> OH + HCO */
    nuki[ 1 * kd + 21 ] += -1 ;
    nuki[ 14 * kd + 21 ] += -1 ;
    nuki[ 6 * kd + 21 ] += +1 ;
    nuki[ 13 * kd + 21 ] += +1 ;

    /*reaction 23: O2 + CH2 <=> CO2 + H2 */
    nuki[ 1 * kd + 22 ] += -1 ;
    nuki[ 14 * kd + 22 ] += -1 ;
    nuki[ 2 * kd + 22 ] += +1 ;
    nuki[ 5 * kd + 22 ] += +1 ;

    /*reaction 24: O2 + CH2 <=> H2O + CO */
    nuki[ 1 * kd + 23 ] += -1 ;
    nuki[ 14 * kd + 23 ] += -1 ;
    nuki[ 3 * kd + 23 ] += +1 ;
    nuki[ 4 * kd + 23 ] += +1 ;

    /*reaction 25: O2 + CH2 <=> O + CH2O */
    nuki[ 1 * kd + 24 ] += -1 ;
    nuki[ 14 * kd + 24 ] += -1 ;
    nuki[ 10 * kd + 24 ] += +1 ;
    nuki[ 12 * kd + 24 ] += +1 ;

    /*reaction 26: O2 + CH2 <=> CO2 + 2 H */
    nuki[ 1 * kd + 25 ] += -1 ;
    nuki[ 14 * kd + 25 ] += -1 ;
    nuki[ 2 * kd + 25 ] += +1 ;
    nuki[ 9 * kd + 25 ] += +2 ;

    /*reaction 27: O2 + CH2 <=> CO + OH + H */
    nuki[ 1 * kd + 26 ] += -1 ;
    nuki[ 14 * kd + 26 ] += -1 ;
    nuki[ 4 * kd + 26 ] += +1 ;
    nuki[ 6 * kd + 26 ] += +1 ;
    nuki[ 9 * kd + 26 ] += +1 ;

    /*reaction 28: CO + CH3O <=> CO2 + CH3 */
    nuki[ 4 * kd + 27 ] += -1 ;
    nuki[ 11 * kd + 27 ] += -1 ;
    nuki[ 2 * kd + 27 ] += +1 ;
    nuki[ 15 * kd + 27 ] += +1 ;

    /*reaction 29: CO + OH <=> CO2 + H */
    nuki[ 4 * kd + 28 ] += -1 ;
    nuki[ 6 * kd + 28 ] += -1 ;
    nuki[ 2 * kd + 28 ] += +1 ;
    nuki[ 9 * kd + 28 ] += +1 ;

    /*reaction 30: OH + O <=> O2 + H */
    nuki[ 6 * kd + 29 ] += -1 ;
    nuki[ 10 * kd + 29 ] += -1 ;
    nuki[ 1 * kd + 29 ] += +1 ;
    nuki[ 9 * kd + 29 ] += +1 ;

    /*reaction 31: HO2 + H <=> 2 OH */
    nuki[ 8 * kd + 30 ] += -1 ;
    nuki[ 9 * kd + 30 ] += -1 ;
    nuki[ 6 * kd + 30 ] += +2 ;

    /*reaction 32: 2 OH <=> H2O + O */
    nuki[ 6 * kd + 31 ] += -2 ;
    nuki[ 3 * kd + 31 ] += +1 ;
    nuki[ 10 * kd + 31 ] += +1 ;

    /*reaction 33: H2 + OH <=> H2O + H */
    nuki[ 5 * kd + 32 ] += -1 ;
    nuki[ 6 * kd + 32 ] += -1 ;
    nuki[ 3 * kd + 32 ] += +1 ;
    nuki[ 9 * kd + 32 ] += +1 ;

    /*reaction 34: 2 HO2 <=> O2 + H2O2 */
    nuki[ 8 * kd + 33 ] += -2 ;
    nuki[ 1 * kd + 33 ] += +1 ;
    nuki[ 7 * kd + 33 ] += +1 ;

    /*reaction 35: OH + CH2O <=> H2O + HCO */
    nuki[ 6 * kd + 34 ] += -1 ;
    nuki[ 12 * kd + 34 ] += -1 ;
    nuki[ 3 * kd + 34 ] += +1 ;
    nuki[ 13 * kd + 34 ] += +1 ;

    /*reaction 36: HO2 + CH2O <=> H2O2 + HCO */
    nuki[ 8 * kd + 35 ] += -1 ;
    nuki[ 12 * kd + 35 ] += -1 ;
    nuki[ 7 * kd + 35 ] += +1 ;
    nuki[ 13 * kd + 35 ] += +1 ;

    /*reaction 37: O2 + HCO <=> CO + HO2 */
    nuki[ 1 * kd + 36 ] += -1 ;
    nuki[ 13 * kd + 36 ] += -1 ;
    nuki[ 4 * kd + 36 ] += +1 ;
    nuki[ 8 * kd + 36 ] += +1 ;

    /*reaction 38: CH3O + CH3 <=> CH2O + CH4 */
    nuki[ 11 * kd + 37 ] += -1 ;
    nuki[ 15 * kd + 37 ] += -1 ;
    nuki[ 12 * kd + 37 ] += +1 ;
    nuki[ 16 * kd + 37 ] += +1 ;

    /*reaction 39: OH + C2H4 <=> CH2O + CH3 */
    nuki[ 6 * kd + 38 ] += -1 ;
    nuki[ 18 * kd + 38 ] += -1 ;
    nuki[ 12 * kd + 38 ] += +1 ;
    nuki[ 15 * kd + 38 ] += +1 ;

    /*reaction 40: OH + C2H4 <=> H2O + C2H3 */
    nuki[ 6 * kd + 39 ] += -1 ;
    nuki[ 18 * kd + 39 ] += -1 ;
    nuki[ 3 * kd + 39 ] += +1 ;
    nuki[ 17 * kd + 39 ] += +1 ;

    /*reaction 41: O2 + C2H3 <=> CH2O + HCO */
    nuki[ 1 * kd + 40 ] += -1 ;
    nuki[ 17 * kd + 40 ] += -1 ;
    nuki[ 12 * kd + 40 ] += +1 ;
    nuki[ 13 * kd + 40 ] += +1 ;

    /*reaction 42: HCO + C2H3 <=> CO + C2H4 */
    nuki[ 13 * kd + 41 ] += -1 ;
    nuki[ 17 * kd + 41 ] += -1 ;
    nuki[ 4 * kd + 41 ] += +1 ;
    nuki[ 18 * kd + 41 ] += +1 ;

    /*reaction 43: O2 + C2H5 <=> HO2 + C2H4 */
    nuki[ 1 * kd + 42 ] += -1 ;
    nuki[ 19 * kd + 42 ] += -1 ;
    nuki[ 8 * kd + 42 ] += +1 ;
    nuki[ 18 * kd + 42 ] += +1 ;

    /*reaction 44: O2 + CH4 <=> HO2 + CH3 */
    nuki[ 1 * kd + 43 ] += -1 ;
    nuki[ 16 * kd + 43 ] += -1 ;
    nuki[ 8 * kd + 43 ] += +1 ;
    nuki[ 15 * kd + 43 ] += +1 ;

    /*reaction 45: OH + HO2 <=> O2 + H2O */
    nuki[ 6 * kd + 44 ] += -1 ;
    nuki[ 8 * kd + 44 ] += -1 ;
    nuki[ 1 * kd + 44 ] += +1 ;
    nuki[ 3 * kd + 44 ] += +1 ;

    /*reaction 46: O2 + CH3 <=> OH + CH2O */
    nuki[ 1 * kd + 45 ] += -1 ;
    nuki[ 15 * kd + 45 ] += -1 ;
    nuki[ 6 * kd + 45 ] += +1 ;
    nuki[ 12 * kd + 45 ] += +1 ;

    /*reaction 47: H + CH4 <=> H2 + CH3 */
    nuki[ 9 * kd + 46 ] += -1 ;
    nuki[ 16 * kd + 46 ] += -1 ;
    nuki[ 5 * kd + 46 ] += +1 ;
    nuki[ 15 * kd + 46 ] += +1 ;

    /*reaction 48: OH + CH4 <=> H2O + CH3 */
    nuki[ 6 * kd + 47 ] += -1 ;
    nuki[ 16 * kd + 47 ] += -1 ;
    nuki[ 3 * kd + 47 ] += +1 ;
    nuki[ 15 * kd + 47 ] += +1 ;

    /*reaction 49: O + CH4 <=> OH + CH3 */
    nuki[ 10 * kd + 48 ] += -1 ;
    nuki[ 16 * kd + 48 ] += -1 ;
    nuki[ 6 * kd + 48 ] += +1 ;
    nuki[ 15 * kd + 48 ] += +1 ;

    /*reaction 50: HO2 + CH4 <=> H2O2 + CH3 */
    nuki[ 8 * kd + 49 ] += -1 ;
    nuki[ 16 * kd + 49 ] += -1 ;
    nuki[ 7 * kd + 49 ] += +1 ;
    nuki[ 15 * kd + 49 ] += +1 ;

    /*reaction 51: CH2 + CH4 <=> 2 CH3 */
    nuki[ 14 * kd + 50 ] += -1 ;
    nuki[ 16 * kd + 50 ] += -1 ;
    nuki[ 15 * kd + 50 ] += +2 ;

    /*reaction 52: C3H6 <=> CH3 + C2H3 */
    nuki[ 22 * kd + 51 ] += -1 ;
    nuki[ 15 * kd + 51 ] += +1 ;
    nuki[ 17 * kd + 51 ] += +1 ;
}


/*Returns the elemental composition  */
/*of the speciesi (mdim is num of elements) */
void CKNCF(int * mdim, int * iwrk, double * restrict rwrk, int * ncf)
{
    int id; /*loop counter */
    int kd = (*mdim); 
    /*Zero ncf */
    for (id = 0; id < kd * 29; ++ id) {
         ncf[id] = 0; 
    }

    /*NC7H16 */
    ncf[ 0 * kd + 0 ] = 16; /*H */
    ncf[ 0 * kd + 1 ] = 7; /*C */

    /*O2 */
    ncf[ 1 * kd + 2 ] = 2; /*O */

    /*CO2 */
    ncf[ 2 * kd + 1 ] = 1; /*C */
    ncf[ 2 * kd + 2 ] = 2; /*O */

    /*H2O */
    ncf[ 3 * kd + 0 ] = 2; /*H */
    ncf[ 3 * kd + 2 ] = 1; /*O */

    /*CO */
    ncf[ 4 * kd + 1 ] = 1; /*C */
    ncf[ 4 * kd + 2 ] = 1; /*O */

    /*H2 */
    ncf[ 5 * kd + 0 ] = 2; /*H */

    /*OH */
    ncf[ 6 * kd + 0 ] = 1; /*H */
    ncf[ 6 * kd + 2 ] = 1; /*O */

    /*H2O2 */
    ncf[ 7 * kd + 0 ] = 2; /*H */
    ncf[ 7 * kd + 2 ] = 2; /*O */

    /*HO2 */
    ncf[ 8 * kd + 0 ] = 1; /*H */
    ncf[ 8 * kd + 2 ] = 2; /*O */

    /*H */
    ncf[ 9 * kd + 0 ] = 1; /*H */

    /*O */
    ncf[ 10 * kd + 2 ] = 1; /*O */

    /*CH3O */
    ncf[ 11 * kd + 0 ] = 3; /*H */
    ncf[ 11 * kd + 1 ] = 1; /*C */
    ncf[ 11 * kd + 2 ] = 1; /*O */

    /*CH2O */
    ncf[ 12 * kd + 0 ] = 2; /*H */
    ncf[ 12 * kd + 1 ] = 1; /*C */
    ncf[ 12 * kd + 2 ] = 1; /*O */

    /*HCO */
    ncf[ 13 * kd + 0 ] = 1; /*H */
    ncf[ 13 * kd + 1 ] = 1; /*C */
    ncf[ 13 * kd + 2 ] = 1; /*O */

    /*CH2 */
    ncf[ 14 * kd + 0 ] = 2; /*H */
    ncf[ 14 * kd + 1 ] = 1; /*C */

    /*CH3 */
    ncf[ 15 * kd + 0 ] = 3; /*H */
    ncf[ 15 * kd + 1 ] = 1; /*C */

    /*CH4 */
    ncf[ 16 * kd + 0 ] = 4; /*H */
    ncf[ 16 * kd + 1 ] = 1; /*C */

    /*C2H3 */
    ncf[ 17 * kd + 0 ] = 3; /*H */
    ncf[ 17 * kd + 1 ] = 2; /*C */

    /*C2H4 */
    ncf[ 18 * kd + 0 ] = 4; /*H */
    ncf[ 18 * kd + 1 ] = 2; /*C */

    /*C2H5 */
    ncf[ 19 * kd + 0 ] = 5; /*H */
    ncf[ 19 * kd + 1 ] = 2; /*C */

    /*C3H4 */
    ncf[ 20 * kd + 0 ] = 4; /*H */
    ncf[ 20 * kd + 1 ] = 3; /*C */

    /*C3H5 */
    ncf[ 21 * kd + 0 ] = 5; /*H */
    ncf[ 21 * kd + 1 ] = 3; /*C */

    /*C3H6 */
    ncf[ 22 * kd + 0 ] = 6; /*H */
    ncf[ 22 * kd + 1 ] = 3; /*C */

    /*C3H7 */
    ncf[ 23 * kd + 0 ] = 7; /*H */
    ncf[ 23 * kd + 1 ] = 3; /*C */

    /*C7H15-2 */
    ncf[ 24 * kd + 0 ] = 15; /*H */
    ncf[ 24 * kd + 1 ] = 7; /*C */

    /*C7H15O2 */
    ncf[ 25 * kd + 0 ] = 15; /*H */
    ncf[ 25 * kd + 1 ] = 7; /*C */
    ncf[ 25 * kd + 2 ] = 2; /*O */

    /*C7KET12 */
    ncf[ 26 * kd + 0 ] = 14; /*H */
    ncf[ 26 * kd + 1 ] = 7; /*C */
    ncf[ 26 * kd + 2 ] = 3; /*O */

    /*C5H11CO */
    ncf[ 27 * kd + 0 ] = 11; /*H */
    ncf[ 27 * kd + 1 ] = 6; /*C */
    ncf[ 27 * kd + 2 ] = 1; /*O */

    /*N2 */
    ncf[ 28 * kd + 3 ] = 2; /*N */

}


/*Returns the arrehenius coefficients  */
/*for all reactions */
void CKABE(int * iwrk, double * restrict rwrk, double * restrict a, double * restrict b, double * restrict e)
{
    for (int i=0; i<52; ++i) {
        a[i] = fwd_A[i];
        b[i] = fwd_beta[i];
        e[i] = fwd_Ea[i];
    }

    return;
}


/*Returns the equil constants for each reaction */
void CKEQC(double * restrict T, double * restrict C, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[29]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: O2 + H + M <=> HO2 + M */
    eqcon[0] *= 1e+06; 

    /*reaction 2: H2O2 + M <=> 2 OH + M */
    eqcon[1] *= 1e-06; 

    /*reaction 3: HCO + M <=> CO + H + M */
    eqcon[2] *= 1e-06; 

    /*reaction 4: NC7H16 + H <=> H2 + C7H15-2 */
    /*eqcon[3] *= 1;  */

    /*reaction 5: NC7H16 + OH <=> H2O + C7H15-2 */
    /*eqcon[4] *= 1;  */

    /*reaction 6: NC7H16 + HO2 <=> H2O2 + C7H15-2 */
    /*eqcon[5] *= 1;  */

    /*reaction 7: NC7H16 + O2 <=> HO2 + C7H15-2 */
    /*eqcon[6] *= 1;  */

    /*reaction 8: O2 + C7H15-2 <=> C7H15O2 */
    eqcon[7] *= 1e+06; 

    /*reaction 9: O2 + C7H15O2 <=> OH + C7KET12 */
    /*eqcon[8] *= 1;  */

    /*reaction 10: C7KET12 <=> OH + CH2O + C5H11CO */
    eqcon[9] *= 1e-12; 

    /*reaction 11: C5H11CO <=> CO + C2H4 + C3H7 */
    eqcon[10] *= 1e-12; 

    /*reaction 12: C7H15-2 <=> C2H4 + C2H5 + C3H6 */
    eqcon[11] *= 1e-12; 

    /*reaction 13: C3H7 <=> CH3 + C2H4 */
    eqcon[12] *= 1e-06; 

    /*reaction 14: C3H7 <=> H + C3H6 */
    eqcon[13] *= 1e-06; 

    /*reaction 15: CH3 + C3H6 <=> CH4 + C3H5 */
    /*eqcon[14] *= 1;  */

    /*reaction 16: O2 + C3H5 <=> HO2 + C3H4 */
    /*eqcon[15] *= 1;  */

    /*reaction 17: OH + C3H4 <=> CH2O + C2H3 */
    /*eqcon[16] *= 1;  */

    /*reaction 18: OH + C3H4 <=> HCO + C2H4 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: HO2 + CH3 <=> OH + CH3O */
    /*eqcon[18] *= 1;  */

    /*reaction 20: OH + CH3 <=> H2O + CH2 */
    /*eqcon[19] *= 1;  */

    /*reaction 21: OH + CH2 <=> H + CH2O */
    /*eqcon[20] *= 1;  */

    /*reaction 22: O2 + CH2 <=> OH + HCO */
    /*eqcon[21] *= 1;  */

    /*reaction 23: O2 + CH2 <=> CO2 + H2 */
    /*eqcon[22] *= 1;  */

    /*reaction 24: O2 + CH2 <=> H2O + CO */
    /*eqcon[23] *= 1;  */

    /*reaction 25: O2 + CH2 <=> O + CH2O */
    /*eqcon[24] *= 1;  */

    /*reaction 26: O2 + CH2 <=> CO2 + 2 H */
    eqcon[25] *= 1e-06; 

    /*reaction 27: O2 + CH2 <=> CO + OH + H */
    eqcon[26] *= 1e-06; 

    /*reaction 28: CO + CH3O <=> CO2 + CH3 */
    /*eqcon[27] *= 1;  */

    /*reaction 29: CO + OH <=> CO2 + H */
    /*eqcon[28] *= 1;  */

    /*reaction 30: OH + O <=> O2 + H */
    /*eqcon[29] *= 1;  */

    /*reaction 31: HO2 + H <=> 2 OH */
    /*eqcon[30] *= 1;  */

    /*reaction 32: 2 OH <=> H2O + O */
    /*eqcon[31] *= 1;  */

    /*reaction 33: H2 + OH <=> H2O + H */
    /*eqcon[32] *= 1;  */

    /*reaction 34: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[33] *= 1;  */

    /*reaction 35: OH + CH2O <=> H2O + HCO */
    /*eqcon[34] *= 1;  */

    /*reaction 36: HO2 + CH2O <=> H2O2 + HCO */
    /*eqcon[35] *= 1;  */

    /*reaction 37: O2 + HCO <=> CO + HO2 */
    /*eqcon[36] *= 1;  */

    /*reaction 38: CH3O + CH3 <=> CH2O + CH4 */
    /*eqcon[37] *= 1;  */

    /*reaction 39: OH + C2H4 <=> CH2O + CH3 */
    /*eqcon[38] *= 1;  */

    /*reaction 40: OH + C2H4 <=> H2O + C2H3 */
    /*eqcon[39] *= 1;  */

    /*reaction 41: O2 + C2H3 <=> CH2O + HCO */
    /*eqcon[40] *= 1;  */

    /*reaction 42: HCO + C2H3 <=> CO + C2H4 */
    /*eqcon[41] *= 1;  */

    /*reaction 43: O2 + C2H5 <=> HO2 + C2H4 */
    /*eqcon[42] *= 1;  */

    /*reaction 44: O2 + CH4 <=> HO2 + CH3 */
    /*eqcon[43] *= 1;  */

    /*reaction 45: OH + HO2 <=> O2 + H2O */
    /*eqcon[44] *= 1;  */

    /*reaction 46: O2 + CH3 <=> OH + CH2O */
    /*eqcon[45] *= 1;  */

    /*reaction 47: H + CH4 <=> H2 + CH3 */
    /*eqcon[46] *= 1;  */

    /*reaction 48: OH + CH4 <=> H2O + CH3 */
    /*eqcon[47] *= 1;  */

    /*reaction 49: O + CH4 <=> OH + CH3 */
    /*eqcon[48] *= 1;  */

    /*reaction 50: HO2 + CH4 <=> H2O2 + CH3 */
    /*eqcon[49] *= 1;  */

    /*reaction 51: CH2 + CH4 <=> 2 CH3 */
    /*eqcon[50] *= 1;  */

    /*reaction 52: C3H6 <=> CH3 + C2H3 */
    eqcon[51] *= 1e-06; 
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mass fractions */
void CKEQYP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[29]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: O2 + H + M <=> HO2 + M */
    eqcon[0] *= 1e+06; 

    /*reaction 2: H2O2 + M <=> 2 OH + M */
    eqcon[1] *= 1e-06; 

    /*reaction 3: HCO + M <=> CO + H + M */
    eqcon[2] *= 1e-06; 

    /*reaction 4: NC7H16 + H <=> H2 + C7H15-2 */
    /*eqcon[3] *= 1;  */

    /*reaction 5: NC7H16 + OH <=> H2O + C7H15-2 */
    /*eqcon[4] *= 1;  */

    /*reaction 6: NC7H16 + HO2 <=> H2O2 + C7H15-2 */
    /*eqcon[5] *= 1;  */

    /*reaction 7: NC7H16 + O2 <=> HO2 + C7H15-2 */
    /*eqcon[6] *= 1;  */

    /*reaction 8: O2 + C7H15-2 <=> C7H15O2 */
    eqcon[7] *= 1e+06; 

    /*reaction 9: O2 + C7H15O2 <=> OH + C7KET12 */
    /*eqcon[8] *= 1;  */

    /*reaction 10: C7KET12 <=> OH + CH2O + C5H11CO */
    eqcon[9] *= 1e-12; 

    /*reaction 11: C5H11CO <=> CO + C2H4 + C3H7 */
    eqcon[10] *= 1e-12; 

    /*reaction 12: C7H15-2 <=> C2H4 + C2H5 + C3H6 */
    eqcon[11] *= 1e-12; 

    /*reaction 13: C3H7 <=> CH3 + C2H4 */
    eqcon[12] *= 1e-06; 

    /*reaction 14: C3H7 <=> H + C3H6 */
    eqcon[13] *= 1e-06; 

    /*reaction 15: CH3 + C3H6 <=> CH4 + C3H5 */
    /*eqcon[14] *= 1;  */

    /*reaction 16: O2 + C3H5 <=> HO2 + C3H4 */
    /*eqcon[15] *= 1;  */

    /*reaction 17: OH + C3H4 <=> CH2O + C2H3 */
    /*eqcon[16] *= 1;  */

    /*reaction 18: OH + C3H4 <=> HCO + C2H4 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: HO2 + CH3 <=> OH + CH3O */
    /*eqcon[18] *= 1;  */

    /*reaction 20: OH + CH3 <=> H2O + CH2 */
    /*eqcon[19] *= 1;  */

    /*reaction 21: OH + CH2 <=> H + CH2O */
    /*eqcon[20] *= 1;  */

    /*reaction 22: O2 + CH2 <=> OH + HCO */
    /*eqcon[21] *= 1;  */

    /*reaction 23: O2 + CH2 <=> CO2 + H2 */
    /*eqcon[22] *= 1;  */

    /*reaction 24: O2 + CH2 <=> H2O + CO */
    /*eqcon[23] *= 1;  */

    /*reaction 25: O2 + CH2 <=> O + CH2O */
    /*eqcon[24] *= 1;  */

    /*reaction 26: O2 + CH2 <=> CO2 + 2 H */
    eqcon[25] *= 1e-06; 

    /*reaction 27: O2 + CH2 <=> CO + OH + H */
    eqcon[26] *= 1e-06; 

    /*reaction 28: CO + CH3O <=> CO2 + CH3 */
    /*eqcon[27] *= 1;  */

    /*reaction 29: CO + OH <=> CO2 + H */
    /*eqcon[28] *= 1;  */

    /*reaction 30: OH + O <=> O2 + H */
    /*eqcon[29] *= 1;  */

    /*reaction 31: HO2 + H <=> 2 OH */
    /*eqcon[30] *= 1;  */

    /*reaction 32: 2 OH <=> H2O + O */
    /*eqcon[31] *= 1;  */

    /*reaction 33: H2 + OH <=> H2O + H */
    /*eqcon[32] *= 1;  */

    /*reaction 34: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[33] *= 1;  */

    /*reaction 35: OH + CH2O <=> H2O + HCO */
    /*eqcon[34] *= 1;  */

    /*reaction 36: HO2 + CH2O <=> H2O2 + HCO */
    /*eqcon[35] *= 1;  */

    /*reaction 37: O2 + HCO <=> CO + HO2 */
    /*eqcon[36] *= 1;  */

    /*reaction 38: CH3O + CH3 <=> CH2O + CH4 */
    /*eqcon[37] *= 1;  */

    /*reaction 39: OH + C2H4 <=> CH2O + CH3 */
    /*eqcon[38] *= 1;  */

    /*reaction 40: OH + C2H4 <=> H2O + C2H3 */
    /*eqcon[39] *= 1;  */

    /*reaction 41: O2 + C2H3 <=> CH2O + HCO */
    /*eqcon[40] *= 1;  */

    /*reaction 42: HCO + C2H3 <=> CO + C2H4 */
    /*eqcon[41] *= 1;  */

    /*reaction 43: O2 + C2H5 <=> HO2 + C2H4 */
    /*eqcon[42] *= 1;  */

    /*reaction 44: O2 + CH4 <=> HO2 + CH3 */
    /*eqcon[43] *= 1;  */

    /*reaction 45: OH + HO2 <=> O2 + H2O */
    /*eqcon[44] *= 1;  */

    /*reaction 46: O2 + CH3 <=> OH + CH2O */
    /*eqcon[45] *= 1;  */

    /*reaction 47: H + CH4 <=> H2 + CH3 */
    /*eqcon[46] *= 1;  */

    /*reaction 48: OH + CH4 <=> H2O + CH3 */
    /*eqcon[47] *= 1;  */

    /*reaction 49: O + CH4 <=> OH + CH3 */
    /*eqcon[48] *= 1;  */

    /*reaction 50: HO2 + CH4 <=> H2O2 + CH3 */
    /*eqcon[49] *= 1;  */

    /*reaction 51: CH2 + CH4 <=> 2 CH3 */
    /*eqcon[50] *= 1;  */

    /*reaction 52: C3H6 <=> CH3 + C2H3 */
    eqcon[51] *= 1e-06; 
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mole fractions */
void CKEQXP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[29]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: O2 + H + M <=> HO2 + M */
    eqcon[0] *= 1e+06; 

    /*reaction 2: H2O2 + M <=> 2 OH + M */
    eqcon[1] *= 1e-06; 

    /*reaction 3: HCO + M <=> CO + H + M */
    eqcon[2] *= 1e-06; 

    /*reaction 4: NC7H16 + H <=> H2 + C7H15-2 */
    /*eqcon[3] *= 1;  */

    /*reaction 5: NC7H16 + OH <=> H2O + C7H15-2 */
    /*eqcon[4] *= 1;  */

    /*reaction 6: NC7H16 + HO2 <=> H2O2 + C7H15-2 */
    /*eqcon[5] *= 1;  */

    /*reaction 7: NC7H16 + O2 <=> HO2 + C7H15-2 */
    /*eqcon[6] *= 1;  */

    /*reaction 8: O2 + C7H15-2 <=> C7H15O2 */
    eqcon[7] *= 1e+06; 

    /*reaction 9: O2 + C7H15O2 <=> OH + C7KET12 */
    /*eqcon[8] *= 1;  */

    /*reaction 10: C7KET12 <=> OH + CH2O + C5H11CO */
    eqcon[9] *= 1e-12; 

    /*reaction 11: C5H11CO <=> CO + C2H4 + C3H7 */
    eqcon[10] *= 1e-12; 

    /*reaction 12: C7H15-2 <=> C2H4 + C2H5 + C3H6 */
    eqcon[11] *= 1e-12; 

    /*reaction 13: C3H7 <=> CH3 + C2H4 */
    eqcon[12] *= 1e-06; 

    /*reaction 14: C3H7 <=> H + C3H6 */
    eqcon[13] *= 1e-06; 

    /*reaction 15: CH3 + C3H6 <=> CH4 + C3H5 */
    /*eqcon[14] *= 1;  */

    /*reaction 16: O2 + C3H5 <=> HO2 + C3H4 */
    /*eqcon[15] *= 1;  */

    /*reaction 17: OH + C3H4 <=> CH2O + C2H3 */
    /*eqcon[16] *= 1;  */

    /*reaction 18: OH + C3H4 <=> HCO + C2H4 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: HO2 + CH3 <=> OH + CH3O */
    /*eqcon[18] *= 1;  */

    /*reaction 20: OH + CH3 <=> H2O + CH2 */
    /*eqcon[19] *= 1;  */

    /*reaction 21: OH + CH2 <=> H + CH2O */
    /*eqcon[20] *= 1;  */

    /*reaction 22: O2 + CH2 <=> OH + HCO */
    /*eqcon[21] *= 1;  */

    /*reaction 23: O2 + CH2 <=> CO2 + H2 */
    /*eqcon[22] *= 1;  */

    /*reaction 24: O2 + CH2 <=> H2O + CO */
    /*eqcon[23] *= 1;  */

    /*reaction 25: O2 + CH2 <=> O + CH2O */
    /*eqcon[24] *= 1;  */

    /*reaction 26: O2 + CH2 <=> CO2 + 2 H */
    eqcon[25] *= 1e-06; 

    /*reaction 27: O2 + CH2 <=> CO + OH + H */
    eqcon[26] *= 1e-06; 

    /*reaction 28: CO + CH3O <=> CO2 + CH3 */
    /*eqcon[27] *= 1;  */

    /*reaction 29: CO + OH <=> CO2 + H */
    /*eqcon[28] *= 1;  */

    /*reaction 30: OH + O <=> O2 + H */
    /*eqcon[29] *= 1;  */

    /*reaction 31: HO2 + H <=> 2 OH */
    /*eqcon[30] *= 1;  */

    /*reaction 32: 2 OH <=> H2O + O */
    /*eqcon[31] *= 1;  */

    /*reaction 33: H2 + OH <=> H2O + H */
    /*eqcon[32] *= 1;  */

    /*reaction 34: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[33] *= 1;  */

    /*reaction 35: OH + CH2O <=> H2O + HCO */
    /*eqcon[34] *= 1;  */

    /*reaction 36: HO2 + CH2O <=> H2O2 + HCO */
    /*eqcon[35] *= 1;  */

    /*reaction 37: O2 + HCO <=> CO + HO2 */
    /*eqcon[36] *= 1;  */

    /*reaction 38: CH3O + CH3 <=> CH2O + CH4 */
    /*eqcon[37] *= 1;  */

    /*reaction 39: OH + C2H4 <=> CH2O + CH3 */
    /*eqcon[38] *= 1;  */

    /*reaction 40: OH + C2H4 <=> H2O + C2H3 */
    /*eqcon[39] *= 1;  */

    /*reaction 41: O2 + C2H3 <=> CH2O + HCO */
    /*eqcon[40] *= 1;  */

    /*reaction 42: HCO + C2H3 <=> CO + C2H4 */
    /*eqcon[41] *= 1;  */

    /*reaction 43: O2 + C2H5 <=> HO2 + C2H4 */
    /*eqcon[42] *= 1;  */

    /*reaction 44: O2 + CH4 <=> HO2 + CH3 */
    /*eqcon[43] *= 1;  */

    /*reaction 45: OH + HO2 <=> O2 + H2O */
    /*eqcon[44] *= 1;  */

    /*reaction 46: O2 + CH3 <=> OH + CH2O */
    /*eqcon[45] *= 1;  */

    /*reaction 47: H + CH4 <=> H2 + CH3 */
    /*eqcon[46] *= 1;  */

    /*reaction 48: OH + CH4 <=> H2O + CH3 */
    /*eqcon[47] *= 1;  */

    /*reaction 49: O + CH4 <=> OH + CH3 */
    /*eqcon[48] *= 1;  */

    /*reaction 50: HO2 + CH4 <=> H2O2 + CH3 */
    /*eqcon[49] *= 1;  */

    /*reaction 51: CH2 + CH4 <=> 2 CH3 */
    /*eqcon[50] *= 1;  */

    /*reaction 52: C3H6 <=> CH3 + C2H3 */
    eqcon[51] *= 1e-06; 
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mass fractions */
void CKEQYR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[29]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: O2 + H + M <=> HO2 + M */
    eqcon[0] *= 1e+06; 

    /*reaction 2: H2O2 + M <=> 2 OH + M */
    eqcon[1] *= 1e-06; 

    /*reaction 3: HCO + M <=> CO + H + M */
    eqcon[2] *= 1e-06; 

    /*reaction 4: NC7H16 + H <=> H2 + C7H15-2 */
    /*eqcon[3] *= 1;  */

    /*reaction 5: NC7H16 + OH <=> H2O + C7H15-2 */
    /*eqcon[4] *= 1;  */

    /*reaction 6: NC7H16 + HO2 <=> H2O2 + C7H15-2 */
    /*eqcon[5] *= 1;  */

    /*reaction 7: NC7H16 + O2 <=> HO2 + C7H15-2 */
    /*eqcon[6] *= 1;  */

    /*reaction 8: O2 + C7H15-2 <=> C7H15O2 */
    eqcon[7] *= 1e+06; 

    /*reaction 9: O2 + C7H15O2 <=> OH + C7KET12 */
    /*eqcon[8] *= 1;  */

    /*reaction 10: C7KET12 <=> OH + CH2O + C5H11CO */
    eqcon[9] *= 1e-12; 

    /*reaction 11: C5H11CO <=> CO + C2H4 + C3H7 */
    eqcon[10] *= 1e-12; 

    /*reaction 12: C7H15-2 <=> C2H4 + C2H5 + C3H6 */
    eqcon[11] *= 1e-12; 

    /*reaction 13: C3H7 <=> CH3 + C2H4 */
    eqcon[12] *= 1e-06; 

    /*reaction 14: C3H7 <=> H + C3H6 */
    eqcon[13] *= 1e-06; 

    /*reaction 15: CH3 + C3H6 <=> CH4 + C3H5 */
    /*eqcon[14] *= 1;  */

    /*reaction 16: O2 + C3H5 <=> HO2 + C3H4 */
    /*eqcon[15] *= 1;  */

    /*reaction 17: OH + C3H4 <=> CH2O + C2H3 */
    /*eqcon[16] *= 1;  */

    /*reaction 18: OH + C3H4 <=> HCO + C2H4 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: HO2 + CH3 <=> OH + CH3O */
    /*eqcon[18] *= 1;  */

    /*reaction 20: OH + CH3 <=> H2O + CH2 */
    /*eqcon[19] *= 1;  */

    /*reaction 21: OH + CH2 <=> H + CH2O */
    /*eqcon[20] *= 1;  */

    /*reaction 22: O2 + CH2 <=> OH + HCO */
    /*eqcon[21] *= 1;  */

    /*reaction 23: O2 + CH2 <=> CO2 + H2 */
    /*eqcon[22] *= 1;  */

    /*reaction 24: O2 + CH2 <=> H2O + CO */
    /*eqcon[23] *= 1;  */

    /*reaction 25: O2 + CH2 <=> O + CH2O */
    /*eqcon[24] *= 1;  */

    /*reaction 26: O2 + CH2 <=> CO2 + 2 H */
    eqcon[25] *= 1e-06; 

    /*reaction 27: O2 + CH2 <=> CO + OH + H */
    eqcon[26] *= 1e-06; 

    /*reaction 28: CO + CH3O <=> CO2 + CH3 */
    /*eqcon[27] *= 1;  */

    /*reaction 29: CO + OH <=> CO2 + H */
    /*eqcon[28] *= 1;  */

    /*reaction 30: OH + O <=> O2 + H */
    /*eqcon[29] *= 1;  */

    /*reaction 31: HO2 + H <=> 2 OH */
    /*eqcon[30] *= 1;  */

    /*reaction 32: 2 OH <=> H2O + O */
    /*eqcon[31] *= 1;  */

    /*reaction 33: H2 + OH <=> H2O + H */
    /*eqcon[32] *= 1;  */

    /*reaction 34: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[33] *= 1;  */

    /*reaction 35: OH + CH2O <=> H2O + HCO */
    /*eqcon[34] *= 1;  */

    /*reaction 36: HO2 + CH2O <=> H2O2 + HCO */
    /*eqcon[35] *= 1;  */

    /*reaction 37: O2 + HCO <=> CO + HO2 */
    /*eqcon[36] *= 1;  */

    /*reaction 38: CH3O + CH3 <=> CH2O + CH4 */
    /*eqcon[37] *= 1;  */

    /*reaction 39: OH + C2H4 <=> CH2O + CH3 */
    /*eqcon[38] *= 1;  */

    /*reaction 40: OH + C2H4 <=> H2O + C2H3 */
    /*eqcon[39] *= 1;  */

    /*reaction 41: O2 + C2H3 <=> CH2O + HCO */
    /*eqcon[40] *= 1;  */

    /*reaction 42: HCO + C2H3 <=> CO + C2H4 */
    /*eqcon[41] *= 1;  */

    /*reaction 43: O2 + C2H5 <=> HO2 + C2H4 */
    /*eqcon[42] *= 1;  */

    /*reaction 44: O2 + CH4 <=> HO2 + CH3 */
    /*eqcon[43] *= 1;  */

    /*reaction 45: OH + HO2 <=> O2 + H2O */
    /*eqcon[44] *= 1;  */

    /*reaction 46: O2 + CH3 <=> OH + CH2O */
    /*eqcon[45] *= 1;  */

    /*reaction 47: H + CH4 <=> H2 + CH3 */
    /*eqcon[46] *= 1;  */

    /*reaction 48: OH + CH4 <=> H2O + CH3 */
    /*eqcon[47] *= 1;  */

    /*reaction 49: O + CH4 <=> OH + CH3 */
    /*eqcon[48] *= 1;  */

    /*reaction 50: HO2 + CH4 <=> H2O2 + CH3 */
    /*eqcon[49] *= 1;  */

    /*reaction 51: CH2 + CH4 <=> 2 CH3 */
    /*eqcon[50] *= 1;  */

    /*reaction 52: C3H6 <=> CH3 + C2H3 */
    eqcon[51] *= 1e-06; 
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mole fractions */
void CKEQXR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[29]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: O2 + H + M <=> HO2 + M */
    eqcon[0] *= 1e+06; 

    /*reaction 2: H2O2 + M <=> 2 OH + M */
    eqcon[1] *= 1e-06; 

    /*reaction 3: HCO + M <=> CO + H + M */
    eqcon[2] *= 1e-06; 

    /*reaction 4: NC7H16 + H <=> H2 + C7H15-2 */
    /*eqcon[3] *= 1;  */

    /*reaction 5: NC7H16 + OH <=> H2O + C7H15-2 */
    /*eqcon[4] *= 1;  */

    /*reaction 6: NC7H16 + HO2 <=> H2O2 + C7H15-2 */
    /*eqcon[5] *= 1;  */

    /*reaction 7: NC7H16 + O2 <=> HO2 + C7H15-2 */
    /*eqcon[6] *= 1;  */

    /*reaction 8: O2 + C7H15-2 <=> C7H15O2 */
    eqcon[7] *= 1e+06; 

    /*reaction 9: O2 + C7H15O2 <=> OH + C7KET12 */
    /*eqcon[8] *= 1;  */

    /*reaction 10: C7KET12 <=> OH + CH2O + C5H11CO */
    eqcon[9] *= 1e-12; 

    /*reaction 11: C5H11CO <=> CO + C2H4 + C3H7 */
    eqcon[10] *= 1e-12; 

    /*reaction 12: C7H15-2 <=> C2H4 + C2H5 + C3H6 */
    eqcon[11] *= 1e-12; 

    /*reaction 13: C3H7 <=> CH3 + C2H4 */
    eqcon[12] *= 1e-06; 

    /*reaction 14: C3H7 <=> H + C3H6 */
    eqcon[13] *= 1e-06; 

    /*reaction 15: CH3 + C3H6 <=> CH4 + C3H5 */
    /*eqcon[14] *= 1;  */

    /*reaction 16: O2 + C3H5 <=> HO2 + C3H4 */
    /*eqcon[15] *= 1;  */

    /*reaction 17: OH + C3H4 <=> CH2O + C2H3 */
    /*eqcon[16] *= 1;  */

    /*reaction 18: OH + C3H4 <=> HCO + C2H4 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: HO2 + CH3 <=> OH + CH3O */
    /*eqcon[18] *= 1;  */

    /*reaction 20: OH + CH3 <=> H2O + CH2 */
    /*eqcon[19] *= 1;  */

    /*reaction 21: OH + CH2 <=> H + CH2O */
    /*eqcon[20] *= 1;  */

    /*reaction 22: O2 + CH2 <=> OH + HCO */
    /*eqcon[21] *= 1;  */

    /*reaction 23: O2 + CH2 <=> CO2 + H2 */
    /*eqcon[22] *= 1;  */

    /*reaction 24: O2 + CH2 <=> H2O + CO */
    /*eqcon[23] *= 1;  */

    /*reaction 25: O2 + CH2 <=> O + CH2O */
    /*eqcon[24] *= 1;  */

    /*reaction 26: O2 + CH2 <=> CO2 + 2 H */
    eqcon[25] *= 1e-06; 

    /*reaction 27: O2 + CH2 <=> CO + OH + H */
    eqcon[26] *= 1e-06; 

    /*reaction 28: CO + CH3O <=> CO2 + CH3 */
    /*eqcon[27] *= 1;  */

    /*reaction 29: CO + OH <=> CO2 + H */
    /*eqcon[28] *= 1;  */

    /*reaction 30: OH + O <=> O2 + H */
    /*eqcon[29] *= 1;  */

    /*reaction 31: HO2 + H <=> 2 OH */
    /*eqcon[30] *= 1;  */

    /*reaction 32: 2 OH <=> H2O + O */
    /*eqcon[31] *= 1;  */

    /*reaction 33: H2 + OH <=> H2O + H */
    /*eqcon[32] *= 1;  */

    /*reaction 34: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[33] *= 1;  */

    /*reaction 35: OH + CH2O <=> H2O + HCO */
    /*eqcon[34] *= 1;  */

    /*reaction 36: HO2 + CH2O <=> H2O2 + HCO */
    /*eqcon[35] *= 1;  */

    /*reaction 37: O2 + HCO <=> CO + HO2 */
    /*eqcon[36] *= 1;  */

    /*reaction 38: CH3O + CH3 <=> CH2O + CH4 */
    /*eqcon[37] *= 1;  */

    /*reaction 39: OH + C2H4 <=> CH2O + CH3 */
    /*eqcon[38] *= 1;  */

    /*reaction 40: OH + C2H4 <=> H2O + C2H3 */
    /*eqcon[39] *= 1;  */

    /*reaction 41: O2 + C2H3 <=> CH2O + HCO */
    /*eqcon[40] *= 1;  */

    /*reaction 42: HCO + C2H3 <=> CO + C2H4 */
    /*eqcon[41] *= 1;  */

    /*reaction 43: O2 + C2H5 <=> HO2 + C2H4 */
    /*eqcon[42] *= 1;  */

    /*reaction 44: O2 + CH4 <=> HO2 + CH3 */
    /*eqcon[43] *= 1;  */

    /*reaction 45: OH + HO2 <=> O2 + H2O */
    /*eqcon[44] *= 1;  */

    /*reaction 46: O2 + CH3 <=> OH + CH2O */
    /*eqcon[45] *= 1;  */

    /*reaction 47: H + CH4 <=> H2 + CH3 */
    /*eqcon[46] *= 1;  */

    /*reaction 48: OH + CH4 <=> H2O + CH3 */
    /*eqcon[47] *= 1;  */

    /*reaction 49: O + CH4 <=> OH + CH3 */
    /*eqcon[48] *= 1;  */

    /*reaction 50: HO2 + CH4 <=> H2O2 + CH3 */
    /*eqcon[49] *= 1;  */

    /*reaction 51: CH2 + CH4 <=> 2 CH3 */
    /*eqcon[50] *= 1;  */

    /*reaction 52: C3H6 <=> CH3 + C2H3 */
    eqcon[51] *= 1e-06; 
}

static double T_save = -1;
#ifdef _OPENMP
#pragma omp threadprivate(T_save)
#endif

static double k_f_save[52];
#ifdef _OPENMP
#pragma omp threadprivate(k_f_save)
#endif

static double Kc_save[52];
#ifdef _OPENMP
#pragma omp threadprivate(Kc_save)
#endif


/*compute the production rate for each species */
void productionRate(double * restrict wdot, double * restrict sc, double T)
{
    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];

    if (T != T_save)
    {
        T_save = T;
        comp_k_f(tc,invT,k_f_save);
        comp_Kc(tc,invT,Kc_save);
    }

    double qdot, q_f[52], q_r[52];
    comp_qfqr(q_f, q_r, sc, tc, invT);

    for (int i = 0; i < 29; ++i) {
        wdot[i] = 0.0;
    }

    qdot = q_f[0]-q_r[0];
    wdot[1] -= qdot;
    wdot[8] += qdot;
    wdot[9] -= qdot;

    qdot = q_f[1]-q_r[1];
    wdot[6] += 2 * qdot;
    wdot[7] -= qdot;

    qdot = q_f[2]-q_r[2];
    wdot[4] += qdot;
    wdot[9] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[3]-q_r[3];
    wdot[0] -= qdot;
    wdot[5] += qdot;
    wdot[9] -= qdot;
    wdot[24] += qdot;

    qdot = q_f[4]-q_r[4];
    wdot[0] -= qdot;
    wdot[3] += qdot;
    wdot[6] -= qdot;
    wdot[24] += qdot;

    qdot = q_f[5]-q_r[5];
    wdot[0] -= qdot;
    wdot[7] += qdot;
    wdot[8] -= qdot;
    wdot[24] += qdot;

    qdot = q_f[6]-q_r[6];
    wdot[0] -= qdot;
    wdot[1] -= qdot;
    wdot[8] += qdot;
    wdot[24] += qdot;

    qdot = q_f[7]-q_r[7];
    wdot[1] -= qdot;
    wdot[24] -= qdot;
    wdot[25] += qdot;

    qdot = q_f[8]-q_r[8];
    wdot[1] -= qdot;
    wdot[6] += qdot;
    wdot[25] -= qdot;
    wdot[26] += qdot;

    qdot = q_f[9]-q_r[9];
    wdot[6] += qdot;
    wdot[12] += qdot;
    wdot[26] -= qdot;
    wdot[27] += qdot;

    qdot = q_f[10]-q_r[10];
    wdot[4] += qdot;
    wdot[18] += qdot;
    wdot[23] += qdot;
    wdot[27] -= qdot;

    qdot = q_f[11]-q_r[11];
    wdot[18] += qdot;
    wdot[19] += qdot;
    wdot[22] += qdot;
    wdot[24] -= qdot;

    qdot = q_f[12]-q_r[12];
    wdot[15] += qdot;
    wdot[18] += qdot;
    wdot[23] -= qdot;

    qdot = q_f[13]-q_r[13];
    wdot[9] += qdot;
    wdot[22] += qdot;
    wdot[23] -= qdot;

    qdot = q_f[14]-q_r[14];
    wdot[15] -= qdot;
    wdot[16] += qdot;
    wdot[21] += qdot;
    wdot[22] -= qdot;

    qdot = q_f[15]-q_r[15];
    wdot[1] -= qdot;
    wdot[8] += qdot;
    wdot[20] += qdot;
    wdot[21] -= qdot;

    qdot = q_f[16]-q_r[16];
    wdot[6] -= qdot;
    wdot[12] += qdot;
    wdot[17] += qdot;
    wdot[20] -= qdot;

    qdot = q_f[17]-q_r[17];
    wdot[6] -= qdot;
    wdot[13] += qdot;
    wdot[18] += qdot;
    wdot[20] -= qdot;

    qdot = q_f[18]-q_r[18];
    wdot[6] += qdot;
    wdot[8] -= qdot;
    wdot[11] += qdot;
    wdot[15] -= qdot;

    qdot = q_f[19]-q_r[19];
    wdot[3] += qdot;
    wdot[6] -= qdot;
    wdot[14] += qdot;
    wdot[15] -= qdot;

    qdot = q_f[20]-q_r[20];
    wdot[6] -= qdot;
    wdot[9] += qdot;
    wdot[12] += qdot;
    wdot[14] -= qdot;

    qdot = q_f[21]-q_r[21];
    wdot[1] -= qdot;
    wdot[6] += qdot;
    wdot[13] += qdot;
    wdot[14] -= qdot;

    qdot = q_f[22]-q_r[22];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[5] += qdot;
    wdot[14] -= qdot;

    qdot = q_f[23]-q_r[23];
    wdot[1] -= qdot;
    wdot[3] += qdot;
    wdot[4] += qdot;
    wdot[14] -= qdot;

    qdot = q_f[24]-q_r[24];
    wdot[1] -= qdot;
    wdot[10] += qdot;
    wdot[12] += qdot;
    wdot[14] -= qdot;

    qdot = q_f[25]-q_r[25];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[9] += 2 * qdot;
    wdot[14] -= qdot;

    qdot = q_f[26]-q_r[26];
    wdot[1] -= qdot;
    wdot[4] += qdot;
    wdot[6] += qdot;
    wdot[9] += qdot;
    wdot[14] -= qdot;

    qdot = q_f[27]-q_r[27];
    wdot[2] += qdot;
    wdot[4] -= qdot;
    wdot[11] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[28]-q_r[28];
    wdot[2] += qdot;
    wdot[4] -= qdot;
    wdot[6] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[29]-q_r[29];
    wdot[1] += qdot;
    wdot[6] -= qdot;
    wdot[9] += qdot;
    wdot[10] -= qdot;

    qdot = q_f[30]-q_r[30];
    wdot[6] += 2 * qdot;
    wdot[8] -= qdot;
    wdot[9] -= qdot;

    qdot = q_f[31]-q_r[31];
    wdot[3] += qdot;
    wdot[6] -= 2 * qdot;
    wdot[10] += qdot;

    qdot = q_f[32]-q_r[32];
    wdot[3] += qdot;
    wdot[5] -= qdot;
    wdot[6] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[33]-q_r[33];
    wdot[1] += qdot;
    wdot[7] += qdot;
    wdot[8] -= 2 * qdot;

    qdot = q_f[34]-q_r[34];
    wdot[3] += qdot;
    wdot[6] -= qdot;
    wdot[12] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[35]-q_r[35];
    wdot[7] += qdot;
    wdot[8] -= qdot;
    wdot[12] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[36]-q_r[36];
    wdot[1] -= qdot;
    wdot[4] += qdot;
    wdot[8] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[37]-q_r[37];
    wdot[11] -= qdot;
    wdot[12] += qdot;
    wdot[15] -= qdot;
    wdot[16] += qdot;

    qdot = q_f[38]-q_r[38];
    wdot[6] -= qdot;
    wdot[12] += qdot;
    wdot[15] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[39]-q_r[39];
    wdot[3] += qdot;
    wdot[6] -= qdot;
    wdot[17] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[40]-q_r[40];
    wdot[1] -= qdot;
    wdot[12] += qdot;
    wdot[13] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[41]-q_r[41];
    wdot[4] += qdot;
    wdot[13] -= qdot;
    wdot[17] -= qdot;
    wdot[18] += qdot;

    qdot = q_f[42]-q_r[42];
    wdot[1] -= qdot;
    wdot[8] += qdot;
    wdot[18] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[43]-q_r[43];
    wdot[1] -= qdot;
    wdot[8] += qdot;
    wdot[15] += qdot;
    wdot[16] -= qdot;

    qdot = q_f[44]-q_r[44];
    wdot[1] += qdot;
    wdot[3] += qdot;
    wdot[6] -= qdot;
    wdot[8] -= qdot;

    qdot = q_f[45]-q_r[45];
    wdot[1] -= qdot;
    wdot[6] += qdot;
    wdot[12] += qdot;
    wdot[15] -= qdot;

    qdot = q_f[46]-q_r[46];
    wdot[5] += qdot;
    wdot[9] -= qdot;
    wdot[15] += qdot;
    wdot[16] -= qdot;

    qdot = q_f[47]-q_r[47];
    wdot[3] += qdot;
    wdot[6] -= qdot;
    wdot[15] += qdot;
    wdot[16] -= qdot;

    qdot = q_f[48]-q_r[48];
    wdot[6] += qdot;
    wdot[10] -= qdot;
    wdot[15] += qdot;
    wdot[16] -= qdot;

    qdot = q_f[49]-q_r[49];
    wdot[7] += qdot;
    wdot[8] -= qdot;
    wdot[15] += qdot;
    wdot[16] -= qdot;

    qdot = q_f[50]-q_r[50];
    wdot[14] -= qdot;
    wdot[15] += 2 * qdot;
    wdot[16] -= qdot;

    qdot = q_f[51]-q_r[51];
    wdot[15] += qdot;
    wdot[17] += qdot;
    wdot[22] -= qdot;

    return;
}

void comp_k_f(double * restrict tc, double invT, double * restrict k_f)
{
#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<52; ++i) {
        k_f[i] = prefactor_units[i] * fwd_A[i]
                    * exp(fwd_beta[i] * tc[0] - activation_units[i] * fwd_Ea[i] * invT);
    };
    return;
}

void comp_Kc(double * restrict tc, double invT, double * restrict Kc)
{
    /*compute the Gibbs free energy */
    double g_RT[29];
    gibbs(g_RT, tc);

    Kc[0] = g_RT[1] - g_RT[8] + g_RT[9];
    Kc[1] = -2*g_RT[6] + g_RT[7];
    Kc[2] = -g_RT[4] - g_RT[9] + g_RT[13];
    Kc[3] = g_RT[0] - g_RT[5] + g_RT[9] - g_RT[24];
    Kc[4] = g_RT[0] - g_RT[3] + g_RT[6] - g_RT[24];
    Kc[5] = g_RT[0] - g_RT[7] + g_RT[8] - g_RT[24];
    Kc[6] = g_RT[0] + g_RT[1] - g_RT[8] - g_RT[24];
    Kc[7] = g_RT[1] + g_RT[24] - g_RT[25];
    Kc[8] = g_RT[1] - g_RT[6] + g_RT[25] - g_RT[26];
    Kc[9] = -g_RT[6] - g_RT[12] + g_RT[26] - g_RT[27];
    Kc[10] = -g_RT[4] - g_RT[18] - g_RT[23] + g_RT[27];
    Kc[11] = -g_RT[18] - g_RT[19] - g_RT[22] + g_RT[24];
    Kc[12] = -g_RT[15] - g_RT[18] + g_RT[23];
    Kc[13] = -g_RT[9] - g_RT[22] + g_RT[23];
    Kc[14] = g_RT[15] - g_RT[16] - g_RT[21] + g_RT[22];
    Kc[15] = g_RT[1] - g_RT[8] - g_RT[20] + g_RT[21];
    Kc[16] = g_RT[6] - g_RT[12] - g_RT[17] + g_RT[20];
    Kc[17] = g_RT[6] - g_RT[13] - g_RT[18] + g_RT[20];
    Kc[18] = -g_RT[6] + g_RT[8] - g_RT[11] + g_RT[15];
    Kc[19] = -g_RT[3] + g_RT[6] - g_RT[14] + g_RT[15];
    Kc[20] = g_RT[6] - g_RT[9] - g_RT[12] + g_RT[14];
    Kc[21] = g_RT[1] - g_RT[6] - g_RT[13] + g_RT[14];
    Kc[22] = g_RT[1] - g_RT[2] - g_RT[5] + g_RT[14];
    Kc[23] = g_RT[1] - g_RT[3] - g_RT[4] + g_RT[14];
    Kc[24] = g_RT[1] - g_RT[10] - g_RT[12] + g_RT[14];
    Kc[25] = g_RT[1] - g_RT[2] - 2*g_RT[9] + g_RT[14];
    Kc[26] = g_RT[1] - g_RT[4] - g_RT[6] - g_RT[9] + g_RT[14];
    Kc[27] = -g_RT[2] + g_RT[4] + g_RT[11] - g_RT[15];
    Kc[28] = -g_RT[2] + g_RT[4] + g_RT[6] - g_RT[9];
    Kc[29] = -g_RT[1] + g_RT[6] - g_RT[9] + g_RT[10];
    Kc[30] = -2*g_RT[6] + g_RT[8] + g_RT[9];
    Kc[31] = -g_RT[3] + 2*g_RT[6] - g_RT[10];
    Kc[32] = -g_RT[3] + g_RT[5] + g_RT[6] - g_RT[9];
    Kc[33] = -g_RT[1] - g_RT[7] + 2*g_RT[8];
    Kc[34] = -g_RT[3] + g_RT[6] + g_RT[12] - g_RT[13];
    Kc[35] = -g_RT[7] + g_RT[8] + g_RT[12] - g_RT[13];
    Kc[36] = g_RT[1] - g_RT[4] - g_RT[8] + g_RT[13];
    Kc[37] = g_RT[11] - g_RT[12] + g_RT[15] - g_RT[16];
    Kc[38] = g_RT[6] - g_RT[12] - g_RT[15] + g_RT[18];
    Kc[39] = -g_RT[3] + g_RT[6] - g_RT[17] + g_RT[18];
    Kc[40] = g_RT[1] - g_RT[12] - g_RT[13] + g_RT[17];
    Kc[41] = -g_RT[4] + g_RT[13] + g_RT[17] - g_RT[18];
    Kc[42] = g_RT[1] - g_RT[8] - g_RT[18] + g_RT[19];
    Kc[43] = g_RT[1] - g_RT[8] - g_RT[15] + g_RT[16];
    Kc[44] = -g_RT[1] - g_RT[3] + g_RT[6] + g_RT[8];
    Kc[45] = g_RT[1] - g_RT[6] - g_RT[12] + g_RT[15];
    Kc[46] = -g_RT[5] + g_RT[9] - g_RT[15] + g_RT[16];
    Kc[47] = -g_RT[3] + g_RT[6] - g_RT[15] + g_RT[16];
    Kc[48] = -g_RT[6] + g_RT[10] - g_RT[15] + g_RT[16];
    Kc[49] = -g_RT[7] + g_RT[8] - g_RT[15] + g_RT[16];
    Kc[50] = g_RT[14] - 2*g_RT[15] + g_RT[16];
    Kc[51] = -g_RT[15] - g_RT[17] + g_RT[22];

#ifdef __INTEL_COMPILER
     #pragma simd
#endif
    for (int i=0; i<52; ++i) {
        Kc[i] = exp(Kc[i]);
    };

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 * invT;
    double refCinv = 1 / refC;

    Kc[0] *= refCinv;
    Kc[1] *= refC;
    Kc[2] *= refC;
    Kc[7] *= refCinv;
    Kc[9] *= refC*refC;
    Kc[10] *= refC*refC;
    Kc[11] *= refC*refC;
    Kc[12] *= refC;
    Kc[13] *= refC;
    Kc[25] *= refC;
    Kc[26] *= refC;
    Kc[51] *= refC;

    return;
}

void comp_qfqr(double * restrict qf, double * restrict qr, double * restrict sc, double * restrict tc, double invT)
{

    /*reaction 1: O2 + H + M <=> HO2 + M */
    qf[0] = sc[1]*sc[9];
    qr[0] = sc[8];

    /*reaction 2: H2O2 + M <=> 2 OH + M */
    qf[1] = sc[7];
    qr[1] = sc[6]*sc[6];

    /*reaction 3: HCO + M <=> CO + H + M */
    qf[2] = sc[13];
    qr[2] = sc[4]*sc[9];

    /*reaction 4: NC7H16 + H <=> H2 + C7H15-2 */
    qf[3] = sc[0]*sc[9];
    qr[3] = sc[5]*sc[24];

    /*reaction 5: NC7H16 + OH <=> H2O + C7H15-2 */
    qf[4] = sc[0]*sc[6];
    qr[4] = sc[3]*sc[24];

    /*reaction 6: NC7H16 + HO2 <=> H2O2 + C7H15-2 */
    qf[5] = sc[0]*sc[8];
    qr[5] = sc[7]*sc[24];

    /*reaction 7: NC7H16 + O2 <=> HO2 + C7H15-2 */
    qf[6] = sc[0]*sc[1];
    qr[6] = sc[8]*sc[24];

    /*reaction 8: O2 + C7H15-2 <=> C7H15O2 */
    qf[7] = sc[1]*sc[24];
    qr[7] = sc[25];

    /*reaction 9: O2 + C7H15O2 <=> OH + C7KET12 */
    qf[8] = sc[1]*sc[25];
    qr[8] = sc[6]*sc[26];

    /*reaction 10: C7KET12 <=> OH + CH2O + C5H11CO */
    qf[9] = sc[26];
    qr[9] = sc[6]*sc[12]*sc[27];

    /*reaction 11: C5H11CO <=> CO + C2H4 + C3H7 */
    qf[10] = sc[27];
    qr[10] = sc[4]*sc[18]*sc[23];

    /*reaction 12: C7H15-2 <=> C2H4 + C2H5 + C3H6 */
    qf[11] = sc[24];
    qr[11] = sc[18]*sc[19]*sc[22];

    /*reaction 13: C3H7 <=> CH3 + C2H4 */
    qf[12] = sc[23];
    qr[12] = sc[15]*sc[18];

    /*reaction 14: C3H7 <=> H + C3H6 */
    qf[13] = sc[23];
    qr[13] = sc[9]*sc[22];

    /*reaction 15: CH3 + C3H6 <=> CH4 + C3H5 */
    qf[14] = sc[15]*sc[22];
    qr[14] = sc[16]*sc[21];

    /*reaction 16: O2 + C3H5 <=> HO2 + C3H4 */
    qf[15] = sc[1]*sc[21];
    qr[15] = sc[8]*sc[20];

    /*reaction 17: OH + C3H4 <=> CH2O + C2H3 */
    qf[16] = sc[6]*sc[20];
    qr[16] = sc[12]*sc[17];

    /*reaction 18: OH + C3H4 <=> HCO + C2H4 */
    qf[17] = sc[6]*sc[20];
    qr[17] = sc[13]*sc[18];

    /*reaction 19: HO2 + CH3 <=> OH + CH3O */
    qf[18] = sc[8]*sc[15];
    qr[18] = sc[6]*sc[11];

    /*reaction 20: OH + CH3 <=> H2O + CH2 */
    qf[19] = sc[6]*sc[15];
    qr[19] = sc[3]*sc[14];

    /*reaction 21: OH + CH2 <=> H + CH2O */
    qf[20] = sc[6]*sc[14];
    qr[20] = sc[9]*sc[12];

    /*reaction 22: O2 + CH2 <=> OH + HCO */
    qf[21] = sc[1]*sc[14];
    qr[21] = sc[6]*sc[13];

    /*reaction 23: O2 + CH2 <=> CO2 + H2 */
    qf[22] = sc[1]*sc[14];
    qr[22] = sc[2]*sc[5];

    /*reaction 24: O2 + CH2 <=> H2O + CO */
    qf[23] = sc[1]*sc[14];
    qr[23] = sc[3]*sc[4];

    /*reaction 25: O2 + CH2 <=> O + CH2O */
    qf[24] = sc[1]*sc[14];
    qr[24] = sc[10]*sc[12];

    /*reaction 26: O2 + CH2 <=> CO2 + 2 H */
    qf[25] = sc[1]*sc[14];
    qr[25] = sc[2]*sc[9]*sc[9];

    /*reaction 27: O2 + CH2 <=> CO + OH + H */
    qf[26] = sc[1]*sc[14];
    qr[26] = sc[4]*sc[6]*sc[9];

    /*reaction 28: CO + CH3O <=> CO2 + CH3 */
    qf[27] = sc[4]*sc[11];
    qr[27] = sc[2]*sc[15];

    /*reaction 29: CO + OH <=> CO2 + H */
    qf[28] = sc[4]*sc[6];
    qr[28] = sc[2]*sc[9];

    /*reaction 30: OH + O <=> O2 + H */
    qf[29] = sc[6]*sc[10];
    qr[29] = sc[1]*sc[9];

    /*reaction 31: HO2 + H <=> 2 OH */
    qf[30] = sc[8]*sc[9];
    qr[30] = sc[6]*sc[6];

    /*reaction 32: 2 OH <=> H2O + O */
    qf[31] = sc[6]*sc[6];
    qr[31] = sc[3]*sc[10];

    /*reaction 33: H2 + OH <=> H2O + H */
    qf[32] = sc[5]*sc[6];
    qr[32] = sc[3]*sc[9];

    /*reaction 34: 2 HO2 <=> O2 + H2O2 */
    qf[33] = sc[8]*sc[8];
    qr[33] = sc[1]*sc[7];

    /*reaction 35: OH + CH2O <=> H2O + HCO */
    qf[34] = sc[6]*sc[12];
    qr[34] = sc[3]*sc[13];

    /*reaction 36: HO2 + CH2O <=> H2O2 + HCO */
    qf[35] = sc[8]*sc[12];
    qr[35] = sc[7]*sc[13];

    /*reaction 37: O2 + HCO <=> CO + HO2 */
    qf[36] = sc[1]*sc[13];
    qr[36] = sc[4]*sc[8];

    /*reaction 38: CH3O + CH3 <=> CH2O + CH4 */
    qf[37] = sc[11]*sc[15];
    qr[37] = sc[12]*sc[16];

    /*reaction 39: OH + C2H4 <=> CH2O + CH3 */
    qf[38] = sc[6]*sc[18];
    qr[38] = sc[12]*sc[15];

    /*reaction 40: OH + C2H4 <=> H2O + C2H3 */
    qf[39] = sc[6]*sc[18];
    qr[39] = sc[3]*sc[17];

    /*reaction 41: O2 + C2H3 <=> CH2O + HCO */
    qf[40] = sc[1]*sc[17];
    qr[40] = sc[12]*sc[13];

    /*reaction 42: HCO + C2H3 <=> CO + C2H4 */
    qf[41] = sc[13]*sc[17];
    qr[41] = sc[4]*sc[18];

    /*reaction 43: O2 + C2H5 <=> HO2 + C2H4 */
    qf[42] = sc[1]*sc[19];
    qr[42] = sc[8]*sc[18];

    /*reaction 44: O2 + CH4 <=> HO2 + CH3 */
    qf[43] = sc[1]*sc[16];
    qr[43] = sc[8]*sc[15];

    /*reaction 45: OH + HO2 <=> O2 + H2O */
    qf[44] = sc[6]*sc[8];
    qr[44] = sc[1]*sc[3];

    /*reaction 46: O2 + CH3 <=> OH + CH2O */
    qf[45] = sc[1]*sc[15];
    qr[45] = sc[6]*sc[12];

    /*reaction 47: H + CH4 <=> H2 + CH3 */
    qf[46] = sc[9]*sc[16];
    qr[46] = sc[5]*sc[15];

    /*reaction 48: OH + CH4 <=> H2O + CH3 */
    qf[47] = sc[6]*sc[16];
    qr[47] = sc[3]*sc[15];

    /*reaction 49: O + CH4 <=> OH + CH3 */
    qf[48] = sc[10]*sc[16];
    qr[48] = sc[6]*sc[15];

    /*reaction 50: HO2 + CH4 <=> H2O2 + CH3 */
    qf[49] = sc[8]*sc[16];
    qr[49] = sc[7]*sc[15];

    /*reaction 51: CH2 + CH4 <=> 2 CH3 */
    qf[50] = sc[14]*sc[16];
    qr[50] = sc[15]*sc[15];

    /*reaction 52: C3H6 <=> CH3 + C2H3 */
    qf[51] = sc[22];
    qr[51] = sc[15]*sc[17];

    double T = tc[1];

    /*compute the mixture concentration */
    double mixture = 0.0;
    for (int i = 0; i < 29; ++i) {
        mixture += sc[i];
    }

    double Corr[52];
    for (int i = 0; i < 52; ++i) {
        Corr[i] = 1.0;
    }

    /* simple three-body correction */
    {
        double alpha;
        alpha = mixture + (TB[0][0] - 1)*sc[2] + (TB[0][1] - 1)*sc[3] + (TB[0][2] - 1)*sc[4] + (TB[0][3] - 1)*sc[5];
        Corr[0] = alpha;
        alpha = mixture + (TB[1][0] - 1)*sc[2] + (TB[1][1] - 1)*sc[3] + (TB[1][2] - 1)*sc[4] + (TB[1][3] - 1)*sc[5];
        Corr[1] = alpha;
        alpha = mixture;
        Corr[2] = alpha;
    }

    for (int i=0; i<52; i++)
    {
        qf[i] *= Corr[i] * k_f_save[i];
        qr[i] *= Corr[i] * k_f_save[i] / Kc_save[i];
    }

    return;
}


/*compute the production rate for each species */
void vproductionRate(int npt, double * restrict wdot, double * restrict sc, double * restrict T)
{
    double k_f_s[52*npt], Kc_s[52*npt], mixture[npt], g_RT[29*npt];
    double tc[5*npt], invT[npt];

#ifdef __INTEL_COMPILER
     #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        tc[0*npt+i] = log(T[i]);
        tc[1*npt+i] = T[i];
        tc[2*npt+i] = T[i]*T[i];
        tc[3*npt+i] = T[i]*T[i]*T[i];
        tc[4*npt+i] = T[i]*T[i]*T[i]*T[i];
        invT[i] = 1.0 / T[i];
    }

    for (int i=0; i<npt; i++) {
        mixture[i] = 0.0;
    }

    for (int n=0; n<29; n++) {
        for (int i=0; i<npt; i++) {
            mixture[i] += sc[n*npt+i];
            wdot[n*npt+i] = 0.0;
        }
    }

    vcomp_k_f(npt, k_f_s, tc, invT);

    vcomp_gibbs(npt, g_RT, tc);

    vcomp_Kc(npt, Kc_s, g_RT, invT);

    vcomp_wdot_1_50(npt, wdot, mixture, sc, k_f_s, Kc_s, tc, invT, T);
    vcomp_wdot_51_52(npt, wdot, mixture, sc, k_f_s, Kc_s, tc, invT, T);
}

void vcomp_k_f(int npt, double * restrict k_f_s, double * restrict tc, double * restrict invT)
{
#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        k_f_s[0*npt+i] = prefactor_units[0] * fwd_A[0] * exp(fwd_beta[0] * tc[i] - activation_units[0] * fwd_Ea[0] * invT[i]);
        k_f_s[1*npt+i] = prefactor_units[1] * fwd_A[1] * exp(fwd_beta[1] * tc[i] - activation_units[1] * fwd_Ea[1] * invT[i]);
        k_f_s[2*npt+i] = prefactor_units[2] * fwd_A[2] * exp(fwd_beta[2] * tc[i] - activation_units[2] * fwd_Ea[2] * invT[i]);
        k_f_s[3*npt+i] = prefactor_units[3] * fwd_A[3] * exp(fwd_beta[3] * tc[i] - activation_units[3] * fwd_Ea[3] * invT[i]);
        k_f_s[4*npt+i] = prefactor_units[4] * fwd_A[4] * exp(fwd_beta[4] * tc[i] - activation_units[4] * fwd_Ea[4] * invT[i]);
        k_f_s[5*npt+i] = prefactor_units[5] * fwd_A[5] * exp(fwd_beta[5] * tc[i] - activation_units[5] * fwd_Ea[5] * invT[i]);
        k_f_s[6*npt+i] = prefactor_units[6] * fwd_A[6] * exp(fwd_beta[6] * tc[i] - activation_units[6] * fwd_Ea[6] * invT[i]);
        k_f_s[7*npt+i] = prefactor_units[7] * fwd_A[7] * exp(fwd_beta[7] * tc[i] - activation_units[7] * fwd_Ea[7] * invT[i]);
        k_f_s[8*npt+i] = prefactor_units[8] * fwd_A[8] * exp(fwd_beta[8] * tc[i] - activation_units[8] * fwd_Ea[8] * invT[i]);
        k_f_s[9*npt+i] = prefactor_units[9] * fwd_A[9] * exp(fwd_beta[9] * tc[i] - activation_units[9] * fwd_Ea[9] * invT[i]);
        k_f_s[10*npt+i] = prefactor_units[10] * fwd_A[10] * exp(fwd_beta[10] * tc[i] - activation_units[10] * fwd_Ea[10] * invT[i]);
        k_f_s[11*npt+i] = prefactor_units[11] * fwd_A[11] * exp(fwd_beta[11] * tc[i] - activation_units[11] * fwd_Ea[11] * invT[i]);
        k_f_s[12*npt+i] = prefactor_units[12] * fwd_A[12] * exp(fwd_beta[12] * tc[i] - activation_units[12] * fwd_Ea[12] * invT[i]);
        k_f_s[13*npt+i] = prefactor_units[13] * fwd_A[13] * exp(fwd_beta[13] * tc[i] - activation_units[13] * fwd_Ea[13] * invT[i]);
        k_f_s[14*npt+i] = prefactor_units[14] * fwd_A[14] * exp(fwd_beta[14] * tc[i] - activation_units[14] * fwd_Ea[14] * invT[i]);
        k_f_s[15*npt+i] = prefactor_units[15] * fwd_A[15] * exp(fwd_beta[15] * tc[i] - activation_units[15] * fwd_Ea[15] * invT[i]);
        k_f_s[16*npt+i] = prefactor_units[16] * fwd_A[16] * exp(fwd_beta[16] * tc[i] - activation_units[16] * fwd_Ea[16] * invT[i]);
        k_f_s[17*npt+i] = prefactor_units[17] * fwd_A[17] * exp(fwd_beta[17] * tc[i] - activation_units[17] * fwd_Ea[17] * invT[i]);
        k_f_s[18*npt+i] = prefactor_units[18] * fwd_A[18] * exp(fwd_beta[18] * tc[i] - activation_units[18] * fwd_Ea[18] * invT[i]);
        k_f_s[19*npt+i] = prefactor_units[19] * fwd_A[19] * exp(fwd_beta[19] * tc[i] - activation_units[19] * fwd_Ea[19] * invT[i]);
        k_f_s[20*npt+i] = prefactor_units[20] * fwd_A[20] * exp(fwd_beta[20] * tc[i] - activation_units[20] * fwd_Ea[20] * invT[i]);
        k_f_s[21*npt+i] = prefactor_units[21] * fwd_A[21] * exp(fwd_beta[21] * tc[i] - activation_units[21] * fwd_Ea[21] * invT[i]);
        k_f_s[22*npt+i] = prefactor_units[22] * fwd_A[22] * exp(fwd_beta[22] * tc[i] - activation_units[22] * fwd_Ea[22] * invT[i]);
        k_f_s[23*npt+i] = prefactor_units[23] * fwd_A[23] * exp(fwd_beta[23] * tc[i] - activation_units[23] * fwd_Ea[23] * invT[i]);
        k_f_s[24*npt+i] = prefactor_units[24] * fwd_A[24] * exp(fwd_beta[24] * tc[i] - activation_units[24] * fwd_Ea[24] * invT[i]);
        k_f_s[25*npt+i] = prefactor_units[25] * fwd_A[25] * exp(fwd_beta[25] * tc[i] - activation_units[25] * fwd_Ea[25] * invT[i]);
        k_f_s[26*npt+i] = prefactor_units[26] * fwd_A[26] * exp(fwd_beta[26] * tc[i] - activation_units[26] * fwd_Ea[26] * invT[i]);
        k_f_s[27*npt+i] = prefactor_units[27] * fwd_A[27] * exp(fwd_beta[27] * tc[i] - activation_units[27] * fwd_Ea[27] * invT[i]);
        k_f_s[28*npt+i] = prefactor_units[28] * fwd_A[28] * exp(fwd_beta[28] * tc[i] - activation_units[28] * fwd_Ea[28] * invT[i]);
        k_f_s[29*npt+i] = prefactor_units[29] * fwd_A[29] * exp(fwd_beta[29] * tc[i] - activation_units[29] * fwd_Ea[29] * invT[i]);
        k_f_s[30*npt+i] = prefactor_units[30] * fwd_A[30] * exp(fwd_beta[30] * tc[i] - activation_units[30] * fwd_Ea[30] * invT[i]);
        k_f_s[31*npt+i] = prefactor_units[31] * fwd_A[31] * exp(fwd_beta[31] * tc[i] - activation_units[31] * fwd_Ea[31] * invT[i]);
        k_f_s[32*npt+i] = prefactor_units[32] * fwd_A[32] * exp(fwd_beta[32] * tc[i] - activation_units[32] * fwd_Ea[32] * invT[i]);
        k_f_s[33*npt+i] = prefactor_units[33] * fwd_A[33] * exp(fwd_beta[33] * tc[i] - activation_units[33] * fwd_Ea[33] * invT[i]);
        k_f_s[34*npt+i] = prefactor_units[34] * fwd_A[34] * exp(fwd_beta[34] * tc[i] - activation_units[34] * fwd_Ea[34] * invT[i]);
        k_f_s[35*npt+i] = prefactor_units[35] * fwd_A[35] * exp(fwd_beta[35] * tc[i] - activation_units[35] * fwd_Ea[35] * invT[i]);
        k_f_s[36*npt+i] = prefactor_units[36] * fwd_A[36] * exp(fwd_beta[36] * tc[i] - activation_units[36] * fwd_Ea[36] * invT[i]);
        k_f_s[37*npt+i] = prefactor_units[37] * fwd_A[37] * exp(fwd_beta[37] * tc[i] - activation_units[37] * fwd_Ea[37] * invT[i]);
        k_f_s[38*npt+i] = prefactor_units[38] * fwd_A[38] * exp(fwd_beta[38] * tc[i] - activation_units[38] * fwd_Ea[38] * invT[i]);
        k_f_s[39*npt+i] = prefactor_units[39] * fwd_A[39] * exp(fwd_beta[39] * tc[i] - activation_units[39] * fwd_Ea[39] * invT[i]);
        k_f_s[40*npt+i] = prefactor_units[40] * fwd_A[40] * exp(fwd_beta[40] * tc[i] - activation_units[40] * fwd_Ea[40] * invT[i]);
        k_f_s[41*npt+i] = prefactor_units[41] * fwd_A[41] * exp(fwd_beta[41] * tc[i] - activation_units[41] * fwd_Ea[41] * invT[i]);
        k_f_s[42*npt+i] = prefactor_units[42] * fwd_A[42] * exp(fwd_beta[42] * tc[i] - activation_units[42] * fwd_Ea[42] * invT[i]);
        k_f_s[43*npt+i] = prefactor_units[43] * fwd_A[43] * exp(fwd_beta[43] * tc[i] - activation_units[43] * fwd_Ea[43] * invT[i]);
        k_f_s[44*npt+i] = prefactor_units[44] * fwd_A[44] * exp(fwd_beta[44] * tc[i] - activation_units[44] * fwd_Ea[44] * invT[i]);
        k_f_s[45*npt+i] = prefactor_units[45] * fwd_A[45] * exp(fwd_beta[45] * tc[i] - activation_units[45] * fwd_Ea[45] * invT[i]);
        k_f_s[46*npt+i] = prefactor_units[46] * fwd_A[46] * exp(fwd_beta[46] * tc[i] - activation_units[46] * fwd_Ea[46] * invT[i]);
        k_f_s[47*npt+i] = prefactor_units[47] * fwd_A[47] * exp(fwd_beta[47] * tc[i] - activation_units[47] * fwd_Ea[47] * invT[i]);
        k_f_s[48*npt+i] = prefactor_units[48] * fwd_A[48] * exp(fwd_beta[48] * tc[i] - activation_units[48] * fwd_Ea[48] * invT[i]);
        k_f_s[49*npt+i] = prefactor_units[49] * fwd_A[49] * exp(fwd_beta[49] * tc[i] - activation_units[49] * fwd_Ea[49] * invT[i]);
        k_f_s[50*npt+i] = prefactor_units[50] * fwd_A[50] * exp(fwd_beta[50] * tc[i] - activation_units[50] * fwd_Ea[50] * invT[i]);
        k_f_s[51*npt+i] = prefactor_units[51] * fwd_A[51] * exp(fwd_beta[51] * tc[i] - activation_units[51] * fwd_Ea[51] * invT[i]);
    }
}

void vcomp_gibbs(int npt, double * restrict g_RT, double * restrict tc)
{
    /*compute the Gibbs free energy */
    for (int i=0; i<npt; i++) {
        double tg[5], g[29];
        tg[0] = tc[0*npt+i];
        tg[1] = tc[1*npt+i];
        tg[2] = tc[2*npt+i];
        tg[3] = tc[3*npt+i];
        tg[4] = tc[4*npt+i];

        gibbs(g, tg);

        g_RT[0*npt+i] = g[0];
        g_RT[1*npt+i] = g[1];
        g_RT[2*npt+i] = g[2];
        g_RT[3*npt+i] = g[3];
        g_RT[4*npt+i] = g[4];
        g_RT[5*npt+i] = g[5];
        g_RT[6*npt+i] = g[6];
        g_RT[7*npt+i] = g[7];
        g_RT[8*npt+i] = g[8];
        g_RT[9*npt+i] = g[9];
        g_RT[10*npt+i] = g[10];
        g_RT[11*npt+i] = g[11];
        g_RT[12*npt+i] = g[12];
        g_RT[13*npt+i] = g[13];
        g_RT[14*npt+i] = g[14];
        g_RT[15*npt+i] = g[15];
        g_RT[16*npt+i] = g[16];
        g_RT[17*npt+i] = g[17];
        g_RT[18*npt+i] = g[18];
        g_RT[19*npt+i] = g[19];
        g_RT[20*npt+i] = g[20];
        g_RT[21*npt+i] = g[21];
        g_RT[22*npt+i] = g[22];
        g_RT[23*npt+i] = g[23];
        g_RT[24*npt+i] = g[24];
        g_RT[25*npt+i] = g[25];
        g_RT[26*npt+i] = g[26];
        g_RT[27*npt+i] = g[27];
        g_RT[28*npt+i] = g[28];
    }
}

void vcomp_Kc(int npt, double * restrict Kc_s, double * restrict g_RT, double * restrict invT)
{
#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
        double refC = (101325. / 8.31451) * invT[i];
        double refCinv = 1.0 / refC;

        Kc_s[0*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[9*npt+i]) - (g_RT[8*npt+i]));
        Kc_s[1*npt+i] = refC * exp((g_RT[7*npt+i]) - (2 * g_RT[6*npt+i]));
        Kc_s[2*npt+i] = refC * exp((g_RT[13*npt+i]) - (g_RT[4*npt+i] + g_RT[9*npt+i]));
        Kc_s[3*npt+i] = exp((g_RT[0*npt+i] + g_RT[9*npt+i]) - (g_RT[5*npt+i] + g_RT[24*npt+i]));
        Kc_s[4*npt+i] = exp((g_RT[0*npt+i] + g_RT[6*npt+i]) - (g_RT[3*npt+i] + g_RT[24*npt+i]));
        Kc_s[5*npt+i] = exp((g_RT[0*npt+i] + g_RT[8*npt+i]) - (g_RT[7*npt+i] + g_RT[24*npt+i]));
        Kc_s[6*npt+i] = exp((g_RT[0*npt+i] + g_RT[1*npt+i]) - (g_RT[8*npt+i] + g_RT[24*npt+i]));
        Kc_s[7*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[24*npt+i]) - (g_RT[25*npt+i]));
        Kc_s[8*npt+i] = exp((g_RT[1*npt+i] + g_RT[25*npt+i]) - (g_RT[6*npt+i] + g_RT[26*npt+i]));
        Kc_s[9*npt+i] = refC*refC * exp((g_RT[26*npt+i]) - (g_RT[6*npt+i] + g_RT[12*npt+i] + g_RT[27*npt+i]));
        Kc_s[10*npt+i] = refC*refC * exp((g_RT[27*npt+i]) - (g_RT[4*npt+i] + g_RT[18*npt+i] + g_RT[23*npt+i]));
        Kc_s[11*npt+i] = refC*refC * exp((g_RT[24*npt+i]) - (g_RT[18*npt+i] + g_RT[19*npt+i] + g_RT[22*npt+i]));
        Kc_s[12*npt+i] = refC * exp((g_RT[23*npt+i]) - (g_RT[15*npt+i] + g_RT[18*npt+i]));
        Kc_s[13*npt+i] = refC * exp((g_RT[23*npt+i]) - (g_RT[9*npt+i] + g_RT[22*npt+i]));
        Kc_s[14*npt+i] = exp((g_RT[15*npt+i] + g_RT[22*npt+i]) - (g_RT[16*npt+i] + g_RT[21*npt+i]));
        Kc_s[15*npt+i] = exp((g_RT[1*npt+i] + g_RT[21*npt+i]) - (g_RT[8*npt+i] + g_RT[20*npt+i]));
        Kc_s[16*npt+i] = exp((g_RT[6*npt+i] + g_RT[20*npt+i]) - (g_RT[12*npt+i] + g_RT[17*npt+i]));
        Kc_s[17*npt+i] = exp((g_RT[6*npt+i] + g_RT[20*npt+i]) - (g_RT[13*npt+i] + g_RT[18*npt+i]));
        Kc_s[18*npt+i] = exp((g_RT[8*npt+i] + g_RT[15*npt+i]) - (g_RT[6*npt+i] + g_RT[11*npt+i]));
        Kc_s[19*npt+i] = exp((g_RT[6*npt+i] + g_RT[15*npt+i]) - (g_RT[3*npt+i] + g_RT[14*npt+i]));
        Kc_s[20*npt+i] = exp((g_RT[6*npt+i] + g_RT[14*npt+i]) - (g_RT[9*npt+i] + g_RT[12*npt+i]));
        Kc_s[21*npt+i] = exp((g_RT[1*npt+i] + g_RT[14*npt+i]) - (g_RT[6*npt+i] + g_RT[13*npt+i]));
        Kc_s[22*npt+i] = exp((g_RT[1*npt+i] + g_RT[14*npt+i]) - (g_RT[2*npt+i] + g_RT[5*npt+i]));
        Kc_s[23*npt+i] = exp((g_RT[1*npt+i] + g_RT[14*npt+i]) - (g_RT[3*npt+i] + g_RT[4*npt+i]));
        Kc_s[24*npt+i] = exp((g_RT[1*npt+i] + g_RT[14*npt+i]) - (g_RT[10*npt+i] + g_RT[12*npt+i]));
        Kc_s[25*npt+i] = refC * exp((g_RT[1*npt+i] + g_RT[14*npt+i]) - (g_RT[2*npt+i] + 2 * g_RT[9*npt+i]));
        Kc_s[26*npt+i] = refC * exp((g_RT[1*npt+i] + g_RT[14*npt+i]) - (g_RT[4*npt+i] + g_RT[6*npt+i] + g_RT[9*npt+i]));
        Kc_s[27*npt+i] = exp((g_RT[4*npt+i] + g_RT[11*npt+i]) - (g_RT[2*npt+i] + g_RT[15*npt+i]));
        Kc_s[28*npt+i] = exp((g_RT[4*npt+i] + g_RT[6*npt+i]) - (g_RT[2*npt+i] + g_RT[9*npt+i]));
        Kc_s[29*npt+i] = exp((g_RT[6*npt+i] + g_RT[10*npt+i]) - (g_RT[1*npt+i] + g_RT[9*npt+i]));
        Kc_s[30*npt+i] = exp((g_RT[8*npt+i] + g_RT[9*npt+i]) - (2 * g_RT[6*npt+i]));
        Kc_s[31*npt+i] = exp((2 * g_RT[6*npt+i]) - (g_RT[3*npt+i] + g_RT[10*npt+i]));
        Kc_s[32*npt+i] = exp((g_RT[5*npt+i] + g_RT[6*npt+i]) - (g_RT[3*npt+i] + g_RT[9*npt+i]));
        Kc_s[33*npt+i] = exp((2 * g_RT[8*npt+i]) - (g_RT[1*npt+i] + g_RT[7*npt+i]));
        Kc_s[34*npt+i] = exp((g_RT[6*npt+i] + g_RT[12*npt+i]) - (g_RT[3*npt+i] + g_RT[13*npt+i]));
        Kc_s[35*npt+i] = exp((g_RT[8*npt+i] + g_RT[12*npt+i]) - (g_RT[7*npt+i] + g_RT[13*npt+i]));
        Kc_s[36*npt+i] = exp((g_RT[1*npt+i] + g_RT[13*npt+i]) - (g_RT[4*npt+i] + g_RT[8*npt+i]));
        Kc_s[37*npt+i] = exp((g_RT[11*npt+i] + g_RT[15*npt+i]) - (g_RT[12*npt+i] + g_RT[16*npt+i]));
        Kc_s[38*npt+i] = exp((g_RT[6*npt+i] + g_RT[18*npt+i]) - (g_RT[12*npt+i] + g_RT[15*npt+i]));
        Kc_s[39*npt+i] = exp((g_RT[6*npt+i] + g_RT[18*npt+i]) - (g_RT[3*npt+i] + g_RT[17*npt+i]));
        Kc_s[40*npt+i] = exp((g_RT[1*npt+i] + g_RT[17*npt+i]) - (g_RT[12*npt+i] + g_RT[13*npt+i]));
        Kc_s[41*npt+i] = exp((g_RT[13*npt+i] + g_RT[17*npt+i]) - (g_RT[4*npt+i] + g_RT[18*npt+i]));
        Kc_s[42*npt+i] = exp((g_RT[1*npt+i] + g_RT[19*npt+i]) - (g_RT[8*npt+i] + g_RT[18*npt+i]));
        Kc_s[43*npt+i] = exp((g_RT[1*npt+i] + g_RT[16*npt+i]) - (g_RT[8*npt+i] + g_RT[15*npt+i]));
        Kc_s[44*npt+i] = exp((g_RT[6*npt+i] + g_RT[8*npt+i]) - (g_RT[1*npt+i] + g_RT[3*npt+i]));
        Kc_s[45*npt+i] = exp((g_RT[1*npt+i] + g_RT[15*npt+i]) - (g_RT[6*npt+i] + g_RT[12*npt+i]));
        Kc_s[46*npt+i] = exp((g_RT[9*npt+i] + g_RT[16*npt+i]) - (g_RT[5*npt+i] + g_RT[15*npt+i]));
        Kc_s[47*npt+i] = exp((g_RT[6*npt+i] + g_RT[16*npt+i]) - (g_RT[3*npt+i] + g_RT[15*npt+i]));
        Kc_s[48*npt+i] = exp((g_RT[10*npt+i] + g_RT[16*npt+i]) - (g_RT[6*npt+i] + g_RT[15*npt+i]));
        Kc_s[49*npt+i] = exp((g_RT[8*npt+i] + g_RT[16*npt+i]) - (g_RT[7*npt+i] + g_RT[15*npt+i]));
        Kc_s[50*npt+i] = exp((g_RT[14*npt+i] + g_RT[16*npt+i]) - (2 * g_RT[15*npt+i]));
        Kc_s[51*npt+i] = refC * exp((g_RT[22*npt+i]) - (g_RT[15*npt+i] + g_RT[17*npt+i]));
    }
}

void vcomp_wdot_1_50(int npt, double * restrict wdot, double * restrict mixture, double * restrict sc,
		double * restrict k_f_s, double * restrict Kc_s,
		double * restrict tc, double * restrict invT, double * restrict T)
{
#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        double qdot, q_f, q_r, phi_f, phi_r, k_f, k_r, Kc;
        double alpha;

        /*reaction 1: O2 + H + M <=> HO2 + M */
        phi_f = sc[1*npt+i]*sc[9*npt+i];
        alpha = mixture[i] + (TB[0][0] - 1)*sc[2*npt+i] + (TB[0][1] - 1)*sc[3*npt+i] + (TB[0][2] - 1)*sc[4*npt+i] + (TB[0][3] - 1)*sc[5*npt+i];
        k_f = alpha * k_f_s[0*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i];
        Kc = Kc_s[0*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[9*npt+i] -= qdot;

        /*reaction 2: H2O2 + M <=> 2 OH + M */
        phi_f = sc[7*npt+i];
        alpha = mixture[i] + (TB[1][0] - 1)*sc[2*npt+i] + (TB[1][1] - 1)*sc[3*npt+i] + (TB[1][2] - 1)*sc[4*npt+i] + (TB[1][3] - 1)*sc[5*npt+i];
        k_f = alpha * k_f_s[1*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[6*npt+i];
        Kc = Kc_s[1*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] += 2 * qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 3: HCO + M <=> CO + H + M */
        phi_f = sc[13*npt+i];
        alpha = mixture[i];
        k_f = alpha * k_f_s[2*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[9*npt+i];
        Kc = Kc_s[2*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] += qdot;
        wdot[9*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;

        /*reaction 4: NC7H16 + H <=> H2 + C7H15-2 */
        phi_f = sc[0*npt+i]*sc[9*npt+i];
        k_f = k_f_s[3*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[24*npt+i];
        Kc = Kc_s[3*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[9*npt+i] -= qdot;
        wdot[24*npt+i] += qdot;

        /*reaction 5: NC7H16 + OH <=> H2O + C7H15-2 */
        phi_f = sc[0*npt+i]*sc[6*npt+i];
        k_f = k_f_s[4*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[24*npt+i];
        Kc = Kc_s[4*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[24*npt+i] += qdot;

        /*reaction 6: NC7H16 + HO2 <=> H2O2 + C7H15-2 */
        phi_f = sc[0*npt+i]*sc[8*npt+i];
        k_f = k_f_s[5*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[24*npt+i];
        Kc = Kc_s[5*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[24*npt+i] += qdot;

        /*reaction 7: NC7H16 + O2 <=> HO2 + C7H15-2 */
        phi_f = sc[0*npt+i]*sc[1*npt+i];
        k_f = k_f_s[6*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[24*npt+i];
        Kc = Kc_s[6*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[24*npt+i] += qdot;

        /*reaction 8: O2 + C7H15-2 <=> C7H15O2 */
        phi_f = sc[1*npt+i]*sc[24*npt+i];
        k_f = k_f_s[7*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[25*npt+i];
        Kc = Kc_s[7*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[24*npt+i] -= qdot;
        wdot[25*npt+i] += qdot;

        /*reaction 9: O2 + C7H15O2 <=> OH + C7KET12 */
        phi_f = sc[1*npt+i]*sc[25*npt+i];
        k_f = k_f_s[8*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[26*npt+i];
        Kc = Kc_s[8*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[25*npt+i] -= qdot;
        wdot[26*npt+i] += qdot;

        /*reaction 10: C7KET12 <=> OH + CH2O + C5H11CO */
        phi_f = sc[26*npt+i];
        k_f = k_f_s[9*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[12*npt+i]*sc[27*npt+i];
        Kc = Kc_s[9*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] += qdot;
        wdot[12*npt+i] += qdot;
        wdot[26*npt+i] -= qdot;
        wdot[27*npt+i] += qdot;

        /*reaction 11: C5H11CO <=> CO + C2H4 + C3H7 */
        phi_f = sc[27*npt+i];
        k_f = k_f_s[10*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[18*npt+i]*sc[23*npt+i];
        Kc = Kc_s[10*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] += qdot;
        wdot[18*npt+i] += qdot;
        wdot[23*npt+i] += qdot;
        wdot[27*npt+i] -= qdot;

        /*reaction 12: C7H15-2 <=> C2H4 + C2H5 + C3H6 */
        phi_f = sc[24*npt+i];
        k_f = k_f_s[11*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[18*npt+i]*sc[19*npt+i]*sc[22*npt+i];
        Kc = Kc_s[11*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[18*npt+i] += qdot;
        wdot[19*npt+i] += qdot;
        wdot[22*npt+i] += qdot;
        wdot[24*npt+i] -= qdot;

        /*reaction 13: C3H7 <=> CH3 + C2H4 */
        phi_f = sc[23*npt+i];
        k_f = k_f_s[12*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[15*npt+i]*sc[18*npt+i];
        Kc = Kc_s[12*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[15*npt+i] += qdot;
        wdot[18*npt+i] += qdot;
        wdot[23*npt+i] -= qdot;

        /*reaction 14: C3H7 <=> H + C3H6 */
        phi_f = sc[23*npt+i];
        k_f = k_f_s[13*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[9*npt+i]*sc[22*npt+i];
        Kc = Kc_s[13*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[9*npt+i] += qdot;
        wdot[22*npt+i] += qdot;
        wdot[23*npt+i] -= qdot;

        /*reaction 15: CH3 + C3H6 <=> CH4 + C3H5 */
        phi_f = sc[15*npt+i]*sc[22*npt+i];
        k_f = k_f_s[14*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[16*npt+i]*sc[21*npt+i];
        Kc = Kc_s[14*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[15*npt+i] -= qdot;
        wdot[16*npt+i] += qdot;
        wdot[21*npt+i] += qdot;
        wdot[22*npt+i] -= qdot;

        /*reaction 16: O2 + C3H5 <=> HO2 + C3H4 */
        phi_f = sc[1*npt+i]*sc[21*npt+i];
        k_f = k_f_s[15*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[20*npt+i];
        Kc = Kc_s[15*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[20*npt+i] += qdot;
        wdot[21*npt+i] -= qdot;

        /*reaction 17: OH + C3H4 <=> CH2O + C2H3 */
        phi_f = sc[6*npt+i]*sc[20*npt+i];
        k_f = k_f_s[16*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[12*npt+i]*sc[17*npt+i];
        Kc = Kc_s[16*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;
        wdot[17*npt+i] += qdot;
        wdot[20*npt+i] -= qdot;

        /*reaction 18: OH + C3H4 <=> HCO + C2H4 */
        phi_f = sc[6*npt+i]*sc[20*npt+i];
        k_f = k_f_s[17*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[13*npt+i]*sc[18*npt+i];
        Kc = Kc_s[17*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;
        wdot[18*npt+i] += qdot;
        wdot[20*npt+i] -= qdot;

        /*reaction 19: HO2 + CH3 <=> OH + CH3O */
        phi_f = sc[8*npt+i]*sc[15*npt+i];
        k_f = k_f_s[18*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[11*npt+i];
        Kc = Kc_s[18*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;
        wdot[15*npt+i] -= qdot;

        /*reaction 20: OH + CH3 <=> H2O + CH2 */
        phi_f = sc[6*npt+i]*sc[15*npt+i];
        k_f = k_f_s[19*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[14*npt+i];
        Kc = Kc_s[19*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;
        wdot[15*npt+i] -= qdot;

        /*reaction 21: OH + CH2 <=> H + CH2O */
        phi_f = sc[6*npt+i]*sc[14*npt+i];
        k_f = k_f_s[20*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[9*npt+i]*sc[12*npt+i];
        Kc = Kc_s[20*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;
        wdot[12*npt+i] += qdot;
        wdot[14*npt+i] -= qdot;

        /*reaction 22: O2 + CH2 <=> OH + HCO */
        phi_f = sc[1*npt+i]*sc[14*npt+i];
        k_f = k_f_s[21*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[13*npt+i];
        Kc = Kc_s[21*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[13*npt+i] += qdot;
        wdot[14*npt+i] -= qdot;

        /*reaction 23: O2 + CH2 <=> CO2 + H2 */
        phi_f = sc[1*npt+i]*sc[14*npt+i];
        k_f = k_f_s[22*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[5*npt+i];
        Kc = Kc_s[22*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[5*npt+i] += qdot;
        wdot[14*npt+i] -= qdot;

        /*reaction 24: O2 + CH2 <=> H2O + CO */
        phi_f = sc[1*npt+i]*sc[14*npt+i];
        k_f = k_f_s[23*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[4*npt+i];
        Kc = Kc_s[23*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[4*npt+i] += qdot;
        wdot[14*npt+i] -= qdot;

        /*reaction 25: O2 + CH2 <=> O + CH2O */
        phi_f = sc[1*npt+i]*sc[14*npt+i];
        k_f = k_f_s[24*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[10*npt+i]*sc[12*npt+i];
        Kc = Kc_s[24*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;
        wdot[12*npt+i] += qdot;
        wdot[14*npt+i] -= qdot;

        /*reaction 26: O2 + CH2 <=> CO2 + 2 H */
        phi_f = sc[1*npt+i]*sc[14*npt+i];
        k_f = k_f_s[25*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[9*npt+i]*sc[9*npt+i];
        Kc = Kc_s[25*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[9*npt+i] += 2 * qdot;
        wdot[14*npt+i] -= qdot;

        /*reaction 27: O2 + CH2 <=> CO + OH + H */
        phi_f = sc[1*npt+i]*sc[14*npt+i];
        k_f = k_f_s[26*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[6*npt+i]*sc[9*npt+i];
        Kc = Kc_s[26*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[6*npt+i] += qdot;
        wdot[9*npt+i] += qdot;
        wdot[14*npt+i] -= qdot;

        /*reaction 28: CO + CH3O <=> CO2 + CH3 */
        phi_f = sc[4*npt+i]*sc[11*npt+i];
        k_f = k_f_s[27*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[15*npt+i];
        Kc = Kc_s[27*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[11*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;

        /*reaction 29: CO + OH <=> CO2 + H */
        phi_f = sc[4*npt+i]*sc[6*npt+i];
        k_f = k_f_s[28*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[9*npt+i];
        Kc = Kc_s[28*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[6*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;

        /*reaction 30: OH + O <=> O2 + H */
        phi_f = sc[6*npt+i]*sc[10*npt+i];
        k_f = k_f_s[29*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[9*npt+i];
        Kc = Kc_s[29*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;

        /*reaction 31: HO2 + H <=> 2 OH */
        phi_f = sc[8*npt+i]*sc[9*npt+i];
        k_f = k_f_s[30*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[6*npt+i];
        Kc = Kc_s[30*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] += 2 * qdot;
        wdot[8*npt+i] -= qdot;
        wdot[9*npt+i] -= qdot;

        /*reaction 32: 2 OH <=> H2O + O */
        phi_f = sc[6*npt+i]*sc[6*npt+i];
        k_f = k_f_s[31*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[10*npt+i];
        Kc = Kc_s[31*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[6*npt+i] -= 2 * qdot;
        wdot[10*npt+i] += qdot;

        /*reaction 33: H2 + OH <=> H2O + H */
        phi_f = sc[5*npt+i]*sc[6*npt+i];
        k_f = k_f_s[32*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[9*npt+i];
        Kc = Kc_s[32*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;

        /*reaction 34: 2 HO2 <=> O2 + H2O2 */
        phi_f = sc[8*npt+i]*sc[8*npt+i];
        k_f = k_f_s[33*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[7*npt+i];
        Kc = Kc_s[33*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[7*npt+i] += qdot;
        wdot[8*npt+i] -= 2 * qdot;

        /*reaction 35: OH + CH2O <=> H2O + HCO */
        phi_f = sc[6*npt+i]*sc[12*npt+i];
        k_f = k_f_s[34*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[13*npt+i];
        Kc = Kc_s[34*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[12*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;

        /*reaction 36: HO2 + CH2O <=> H2O2 + HCO */
        phi_f = sc[8*npt+i]*sc[12*npt+i];
        k_f = k_f_s[35*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[13*npt+i];
        Kc = Kc_s[35*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[12*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;

        /*reaction 37: O2 + HCO <=> CO + HO2 */
        phi_f = sc[1*npt+i]*sc[13*npt+i];
        k_f = k_f_s[36*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[8*npt+i];
        Kc = Kc_s[36*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[8*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;

        /*reaction 38: CH3O + CH3 <=> CH2O + CH4 */
        phi_f = sc[11*npt+i]*sc[15*npt+i];
        k_f = k_f_s[37*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[12*npt+i]*sc[16*npt+i];
        Kc = Kc_s[37*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[11*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;
        wdot[15*npt+i] -= qdot;
        wdot[16*npt+i] += qdot;

        /*reaction 39: OH + C2H4 <=> CH2O + CH3 */
        phi_f = sc[6*npt+i]*sc[18*npt+i];
        k_f = k_f_s[38*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[12*npt+i]*sc[15*npt+i];
        Kc = Kc_s[38*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;
        wdot[15*npt+i] += qdot;
        wdot[18*npt+i] -= qdot;

        /*reaction 40: OH + C2H4 <=> H2O + C2H3 */
        phi_f = sc[6*npt+i]*sc[18*npt+i];
        k_f = k_f_s[39*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[17*npt+i];
        Kc = Kc_s[39*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[17*npt+i] += qdot;
        wdot[18*npt+i] -= qdot;

        /*reaction 41: O2 + C2H3 <=> CH2O + HCO */
        phi_f = sc[1*npt+i]*sc[17*npt+i];
        k_f = k_f_s[40*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[12*npt+i]*sc[13*npt+i];
        Kc = Kc_s[40*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;
        wdot[13*npt+i] += qdot;
        wdot[17*npt+i] -= qdot;

        /*reaction 42: HCO + C2H3 <=> CO + C2H4 */
        phi_f = sc[13*npt+i]*sc[17*npt+i];
        k_f = k_f_s[41*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[18*npt+i];
        Kc = Kc_s[41*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;
        wdot[17*npt+i] -= qdot;
        wdot[18*npt+i] += qdot;

        /*reaction 43: O2 + C2H5 <=> HO2 + C2H4 */
        phi_f = sc[1*npt+i]*sc[19*npt+i];
        k_f = k_f_s[42*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[18*npt+i];
        Kc = Kc_s[42*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[18*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 44: O2 + CH4 <=> HO2 + CH3 */
        phi_f = sc[1*npt+i]*sc[16*npt+i];
        k_f = k_f_s[43*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[15*npt+i];
        Kc = Kc_s[43*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[15*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;

        /*reaction 45: OH + HO2 <=> O2 + H2O */
        phi_f = sc[6*npt+i]*sc[8*npt+i];
        k_f = k_f_s[44*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[3*npt+i];
        Kc = Kc_s[44*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[3*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[8*npt+i] -= qdot;

        /*reaction 46: O2 + CH3 <=> OH + CH2O */
        phi_f = sc[1*npt+i]*sc[15*npt+i];
        k_f = k_f_s[45*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[12*npt+i];
        Kc = Kc_s[45*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[12*npt+i] += qdot;
        wdot[15*npt+i] -= qdot;

        /*reaction 47: H + CH4 <=> H2 + CH3 */
        phi_f = sc[9*npt+i]*sc[16*npt+i];
        k_f = k_f_s[46*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[15*npt+i];
        Kc = Kc_s[46*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] += qdot;
        wdot[9*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;

        /*reaction 48: OH + CH4 <=> H2O + CH3 */
        phi_f = sc[6*npt+i]*sc[16*npt+i];
        k_f = k_f_s[47*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[15*npt+i];
        Kc = Kc_s[47*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;

        /*reaction 49: O + CH4 <=> OH + CH3 */
        phi_f = sc[10*npt+i]*sc[16*npt+i];
        k_f = k_f_s[48*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[15*npt+i];
        Kc = Kc_s[48*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;

        /*reaction 50: HO2 + CH4 <=> H2O2 + CH3 */
        phi_f = sc[8*npt+i]*sc[16*npt+i];
        k_f = k_f_s[49*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[15*npt+i];
        Kc = Kc_s[49*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;
    }
}

void vcomp_wdot_51_52(int npt, double * restrict wdot, double * restrict mixture, double * restrict sc,
		double * restrict k_f_s, double * restrict Kc_s,
		double * restrict tc, double * restrict invT, double * restrict T)
{
#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        double qdot, q_f, q_r, phi_f, phi_r, k_f, k_r, Kc;

        /*reaction 51: CH2 + CH4 <=> 2 CH3 */
        phi_f = sc[14*npt+i]*sc[16*npt+i];
        k_f = k_f_s[50*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[15*npt+i]*sc[15*npt+i];
        Kc = Kc_s[50*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[14*npt+i] -= qdot;
        wdot[15*npt+i] += 2 * qdot;
        wdot[16*npt+i] -= qdot;

        /*reaction 52: C3H6 <=> CH3 + C2H3 */
        phi_f = sc[22*npt+i];
        k_f = k_f_s[51*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[15*npt+i]*sc[17*npt+i];
        Kc = Kc_s[51*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[15*npt+i] += qdot;
        wdot[17*npt+i] += qdot;
        wdot[22*npt+i] -= qdot;
    }
}


/*compute the reaction Jacobian */
void DWDOT(double * restrict J, double * restrict sc, double * restrict Tp, int * consP)
{
    double c[29];

    for (int k=0; k<29; k++) {
        c[k] = 1.e6 * sc[k];
    }

    aJacobian(J, c, *Tp, *consP);

    /* dwdot[k]/dT */
    for (int k=0; k<29; k++) {
        J[870+k] *= 1.e-6;
    }

    /* dTdot/d[X] */
    for (int k=0; k<29; k++) {
        J[k*30+29] *= 1.e6;
    }

    return;
}

/*compute the sparsity pattern Jacobian */
void SPARSITY_INFO( int * nJdata)
{
    double c[29];
    double J[900];

    for (int k=0; k<29; k++) {
        c[k] = 1.0/29.0;
    }

    aJacobian(J, c, 1500.0, 0);

    int nJdata_tmp = 0;
    for (int k=0; k<900; k++) {
        if(J[k] != 0.0){
	    nJdata_tmp = nJdata_tmp + 1;
	}
    }

    nJdata[0] = nJdata_tmp;
}

/*compute the sparsity pattern Jacobian */
void SPARSITY_PREPROC(int * restrict rowVals, int * restrict colPtrs)
{
    double c[29];
    double J[900];

    for (int k=0; k<29; k++) {
        c[k] = 1.0/29.0;
    }

    aJacobian(J, c, 1500.0, 0);

    colPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int k=0; k<30; k++) {
        for (int l=0; l<30; l++) {
            // Debug version
            if(J[30*k + l] != 0.0){
	        rowVals[nJdata_tmp] = l+1; 
	        nJdata_tmp = nJdata_tmp + 1;
	    }
	}
        colPtrs[k+1] = nJdata_tmp;
    }
}

/*compute the reaction Jacobian */
void aJacobian(double * restrict J, double * restrict sc, double T, int consP)
{
    for (int i=0; i<900; i++) {
        J[i] = 0.0;
    }

    double wdot[29];
    for (int k=0; k<29; k++) {
        wdot[k] = 0.0;
    }

    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];
    double invT2 = invT * invT;

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 / T;
    double refCinv = 1.0 / refC;

    /*compute the mixture concentration */
    double mixture = 0.0;
    for (int k = 0; k < 29; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    double g_RT[29];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    double h_RT[29];
    speciesEnthalpy(h_RT, tc);

    double phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    double dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    double dqdci, dcdc_fac, dqdc[29];
    double Pr, fPr, F, k_0, logPr;
    double logFcent, troe_c, troe_n, troePr_den, troePr, troe;
    double Fcent1, Fcent2, Fcent3, Fcent;
    double dlogFdc, dlogFdn, dlogFdcn_fac;
    double dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
    const double ln10 = log(10.0);
    const double log10e = 1.0/log(10.0);
    /*reaction 1: O2 + H + M <=> HO2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[0][0] - 1)*sc[2] + (TB[0][1] - 1)*sc[3] + (TB[0][2] - 1)*sc[4] + (TB[0][3] - 1)*sc[5];
    /* forward */
    phi_f = sc[1]*sc[9];
    k_f = prefactor_units[0] * fwd_A[0]
                * exp(fwd_beta[0] * tc[0] - activation_units[0] * fwd_Ea[0] * invT);
    dlnkfdT = fwd_beta[0] * invT + activation_units[0] * fwd_Ea[0] * invT2;
    /* reverse */
    phi_r = sc[8];
    Kc = refCinv * exp(g_RT[1] - g_RT[8] + g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[9]) + (h_RT[8]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    wdot[9] -= q; /* H */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[O2] */
        dqdci =  + k_f*sc[9];
        J[31] -= dqdci;               /* dwdot[O2]/d[O2] */
        J[38] += dqdci;               /* dwdot[HO2]/d[O2] */
        J[39] -= dqdci;               /* dwdot[H]/d[O2] */
        /* d()/d[CO2] */
        dqdci = (TB[0][0] - 1)*q_nocor;
        J[61] -= dqdci;               /* dwdot[O2]/d[CO2] */
        J[68] += dqdci;               /* dwdot[HO2]/d[CO2] */
        J[69] -= dqdci;               /* dwdot[H]/d[CO2] */
        /* d()/d[H2O] */
        dqdci = (TB[0][1] - 1)*q_nocor;
        J[91] -= dqdci;               /* dwdot[O2]/d[H2O] */
        J[98] += dqdci;               /* dwdot[HO2]/d[H2O] */
        J[99] -= dqdci;               /* dwdot[H]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[0][2] - 1)*q_nocor;
        J[121] -= dqdci;              /* dwdot[O2]/d[CO] */
        J[128] += dqdci;              /* dwdot[HO2]/d[CO] */
        J[129] -= dqdci;              /* dwdot[H]/d[CO] */
        /* d()/d[H2] */
        dqdci = (TB[0][3] - 1)*q_nocor;
        J[151] -= dqdci;              /* dwdot[O2]/d[H2] */
        J[158] += dqdci;              /* dwdot[HO2]/d[H2] */
        J[159] -= dqdci;              /* dwdot[H]/d[H2] */
        /* d()/d[HO2] */
        dqdci =  - k_r;
        J[241] -= dqdci;              /* dwdot[O2]/d[HO2] */
        J[248] += dqdci;              /* dwdot[HO2]/d[HO2] */
        J[249] -= dqdci;              /* dwdot[H]/d[HO2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[1];
        J[271] -= dqdci;              /* dwdot[O2]/d[H] */
        J[278] += dqdci;              /* dwdot[HO2]/d[H] */
        J[279] -= dqdci;              /* dwdot[H]/d[H] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[9];
        dqdc[2] = TB[0][0];
        dqdc[3] = TB[0][1];
        dqdc[4] = TB[0][2];
        dqdc[5] = TB[0][3];
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac - k_r;
        dqdc[9] = dcdc_fac + k_f*sc[1];
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = dcdc_fac;
        dqdc[23] = dcdc_fac;
        dqdc[24] = dcdc_fac;
        dqdc[25] = dcdc_fac;
        dqdc[26] = dcdc_fac;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        for (int k=0; k<29; k++) {
            J[30*k+1] -= dqdc[k];
            J[30*k+8] += dqdc[k];
            J[30*k+9] -= dqdc[k];
        }
    }
    J[871] -= dqdT; /* dwdot[O2]/dT */
    J[878] += dqdT; /* dwdot[HO2]/dT */
    J[879] -= dqdT; /* dwdot[H]/dT */

    /*reaction 2: H2O2 + M <=> 2 OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[1][0] - 1)*sc[2] + (TB[1][1] - 1)*sc[3] + (TB[1][2] - 1)*sc[4] + (TB[1][3] - 1)*sc[5];
    /* forward */
    phi_f = sc[7];
    k_f = prefactor_units[1] * fwd_A[1]
                * exp(fwd_beta[1] * tc[0] - activation_units[1] * fwd_Ea[1] * invT);
    dlnkfdT = fwd_beta[1] * invT + activation_units[1] * fwd_Ea[1] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[6];
    Kc = refC * exp(-2*g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7]) + (2*h_RT[6]) - 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[6] += 2 * q; /* OH */
    wdot[7] -= q; /* H2O2 */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[CO2] */
        dqdci = (TB[1][0] - 1)*q_nocor;
        J[66] += 2 * dqdci;           /* dwdot[OH]/d[CO2] */
        J[67] -= dqdci;               /* dwdot[H2O2]/d[CO2] */
        /* d()/d[H2O] */
        dqdci = (TB[1][1] - 1)*q_nocor;
        J[96] += 2 * dqdci;           /* dwdot[OH]/d[H2O] */
        J[97] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[1][2] - 1)*q_nocor;
        J[126] += 2 * dqdci;          /* dwdot[OH]/d[CO] */
        J[127] -= dqdci;              /* dwdot[H2O2]/d[CO] */
        /* d()/d[H2] */
        dqdci = (TB[1][3] - 1)*q_nocor;
        J[156] += 2 * dqdci;          /* dwdot[OH]/d[H2] */
        J[157] -= dqdci;              /* dwdot[H2O2]/d[H2] */
        /* d()/d[OH] */
        dqdci =  - k_r*2*sc[6];
        J[186] += 2 * dqdci;          /* dwdot[OH]/d[OH] */
        J[187] -= dqdci;              /* dwdot[H2O2]/d[OH] */
        /* d()/d[H2O2] */
        dqdci =  + k_f;
        J[216] += 2 * dqdci;          /* dwdot[OH]/d[H2O2] */
        J[217] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = TB[1][0];
        dqdc[3] = TB[1][1];
        dqdc[4] = TB[1][2];
        dqdc[5] = TB[1][3];
        dqdc[6] = dcdc_fac - k_r*2*sc[6];
        dqdc[7] = dcdc_fac + k_f;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = dcdc_fac;
        dqdc[23] = dcdc_fac;
        dqdc[24] = dcdc_fac;
        dqdc[25] = dcdc_fac;
        dqdc[26] = dcdc_fac;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        for (int k=0; k<29; k++) {
            J[30*k+6] += 2 * dqdc[k];
            J[30*k+7] -= dqdc[k];
        }
    }
    J[876] += 2 * dqdT; /* dwdot[OH]/dT */
    J[877] -= dqdT; /* dwdot[H2O2]/dT */

    /*reaction 3: HCO + M <=> CO + H + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture;
    /* forward */
    phi_f = sc[13];
    k_f = prefactor_units[2] * fwd_A[2]
                * exp(fwd_beta[2] * tc[0] - activation_units[2] * fwd_Ea[2] * invT);
    dlnkfdT = fwd_beta[2] * invT + activation_units[2] * fwd_Ea[2] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[9];
    Kc = refC * exp(-g_RT[4] - g_RT[9] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[13]) + (h_RT[4] + h_RT[9]) - 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* CO */
    wdot[9] += q; /* H */
    wdot[13] -= q; /* HCO */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[CO] */
        dqdci =  - k_r*sc[9];
        J[124] += dqdci;              /* dwdot[CO]/d[CO] */
        J[129] += dqdci;              /* dwdot[H]/d[CO] */
        J[133] -= dqdci;              /* dwdot[HCO]/d[CO] */
        /* d()/d[H] */
        dqdci =  - k_r*sc[4];
        J[274] += dqdci;              /* dwdot[CO]/d[H] */
        J[279] += dqdci;              /* dwdot[H]/d[H] */
        J[283] -= dqdci;              /* dwdot[HCO]/d[H] */
        /* d()/d[HCO] */
        dqdci =  + k_f;
        J[394] += dqdci;              /* dwdot[CO]/d[HCO] */
        J[399] += dqdci;              /* dwdot[H]/d[HCO] */
        J[403] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac - k_r*sc[9];
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac - k_r*sc[4];
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = dcdc_fac + k_f;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = dcdc_fac;
        dqdc[23] = dcdc_fac;
        dqdc[24] = dcdc_fac;
        dqdc[25] = dcdc_fac;
        dqdc[26] = dcdc_fac;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        for (int k=0; k<29; k++) {
            J[30*k+4] += dqdc[k];
            J[30*k+9] += dqdc[k];
            J[30*k+13] -= dqdc[k];
        }
    }
    J[874] += dqdT; /* dwdot[CO]/dT */
    J[879] += dqdT; /* dwdot[H]/dT */
    J[883] -= dqdT; /* dwdot[HCO]/dT */

    /*reaction 4: NC7H16 + H <=> H2 + C7H15-2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[9];
    k_f = prefactor_units[3] * fwd_A[3]
                * exp(fwd_beta[3] * tc[0] - activation_units[3] * fwd_Ea[3] * invT);
    dlnkfdT = fwd_beta[3] * invT + activation_units[3] * fwd_Ea[3] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[24];
    Kc = exp(g_RT[0] - g_RT[5] + g_RT[9] - g_RT[24]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[9]) + (h_RT[5] + h_RT[24]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* NC7H16 */
    wdot[5] += q; /* H2 */
    wdot[9] -= q; /* H */
    wdot[24] += q; /* C7H15-2 */
    /* d()/d[NC7H16] */
    dqdci =  + k_f*sc[9];
    J[0] -= dqdci;                /* dwdot[NC7H16]/d[NC7H16] */
    J[5] += dqdci;                /* dwdot[H2]/d[NC7H16] */
    J[9] -= dqdci;                /* dwdot[H]/d[NC7H16] */
    J[24] += dqdci;               /* dwdot[C7H15-2]/d[NC7H16] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[24];
    J[150] -= dqdci;              /* dwdot[NC7H16]/d[H2] */
    J[155] += dqdci;              /* dwdot[H2]/d[H2] */
    J[159] -= dqdci;              /* dwdot[H]/d[H2] */
    J[174] += dqdci;              /* dwdot[C7H15-2]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[0];
    J[270] -= dqdci;              /* dwdot[NC7H16]/d[H] */
    J[275] += dqdci;              /* dwdot[H2]/d[H] */
    J[279] -= dqdci;              /* dwdot[H]/d[H] */
    J[294] += dqdci;              /* dwdot[C7H15-2]/d[H] */
    /* d()/d[C7H15-2] */
    dqdci =  - k_r*sc[5];
    J[720] -= dqdci;              /* dwdot[NC7H16]/d[C7H15-2] */
    J[725] += dqdci;              /* dwdot[H2]/d[C7H15-2] */
    J[729] -= dqdci;              /* dwdot[H]/d[C7H15-2] */
    J[744] += dqdci;              /* dwdot[C7H15-2]/d[C7H15-2] */
    /* d()/dT */
    J[870] -= dqdT;               /* dwdot[NC7H16]/dT */
    J[875] += dqdT;               /* dwdot[H2]/dT */
    J[879] -= dqdT;               /* dwdot[H]/dT */
    J[894] += dqdT;               /* dwdot[C7H15-2]/dT */

    /*reaction 5: NC7H16 + OH <=> H2O + C7H15-2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[6];
    k_f = prefactor_units[4] * fwd_A[4]
                * exp(fwd_beta[4] * tc[0] - activation_units[4] * fwd_Ea[4] * invT);
    dlnkfdT = fwd_beta[4] * invT + activation_units[4] * fwd_Ea[4] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[24];
    Kc = exp(g_RT[0] - g_RT[3] + g_RT[6] - g_RT[24]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[6]) + (h_RT[3] + h_RT[24]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* NC7H16 */
    wdot[3] += q; /* H2O */
    wdot[6] -= q; /* OH */
    wdot[24] += q; /* C7H15-2 */
    /* d()/d[NC7H16] */
    dqdci =  + k_f*sc[6];
    J[0] -= dqdci;                /* dwdot[NC7H16]/d[NC7H16] */
    J[3] += dqdci;                /* dwdot[H2O]/d[NC7H16] */
    J[6] -= dqdci;                /* dwdot[OH]/d[NC7H16] */
    J[24] += dqdci;               /* dwdot[C7H15-2]/d[NC7H16] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[24];
    J[90] -= dqdci;               /* dwdot[NC7H16]/d[H2O] */
    J[93] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[96] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[114] += dqdci;              /* dwdot[C7H15-2]/d[H2O] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[0];
    J[180] -= dqdci;              /* dwdot[NC7H16]/d[OH] */
    J[183] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[186] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[204] += dqdci;              /* dwdot[C7H15-2]/d[OH] */
    /* d()/d[C7H15-2] */
    dqdci =  - k_r*sc[3];
    J[720] -= dqdci;              /* dwdot[NC7H16]/d[C7H15-2] */
    J[723] += dqdci;              /* dwdot[H2O]/d[C7H15-2] */
    J[726] -= dqdci;              /* dwdot[OH]/d[C7H15-2] */
    J[744] += dqdci;              /* dwdot[C7H15-2]/d[C7H15-2] */
    /* d()/dT */
    J[870] -= dqdT;               /* dwdot[NC7H16]/dT */
    J[873] += dqdT;               /* dwdot[H2O]/dT */
    J[876] -= dqdT;               /* dwdot[OH]/dT */
    J[894] += dqdT;               /* dwdot[C7H15-2]/dT */

    /*reaction 6: NC7H16 + HO2 <=> H2O2 + C7H15-2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[8];
    k_f = prefactor_units[5] * fwd_A[5]
                * exp(fwd_beta[5] * tc[0] - activation_units[5] * fwd_Ea[5] * invT);
    dlnkfdT = fwd_beta[5] * invT + activation_units[5] * fwd_Ea[5] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[24];
    Kc = exp(g_RT[0] - g_RT[7] + g_RT[8] - g_RT[24]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[8]) + (h_RT[7] + h_RT[24]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* NC7H16 */
    wdot[7] += q; /* H2O2 */
    wdot[8] -= q; /* HO2 */
    wdot[24] += q; /* C7H15-2 */
    /* d()/d[NC7H16] */
    dqdci =  + k_f*sc[8];
    J[0] -= dqdci;                /* dwdot[NC7H16]/d[NC7H16] */
    J[7] += dqdci;                /* dwdot[H2O2]/d[NC7H16] */
    J[8] -= dqdci;                /* dwdot[HO2]/d[NC7H16] */
    J[24] += dqdci;               /* dwdot[C7H15-2]/d[NC7H16] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[24];
    J[210] -= dqdci;              /* dwdot[NC7H16]/d[H2O2] */
    J[217] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
    J[218] -= dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[234] += dqdci;              /* dwdot[C7H15-2]/d[H2O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[0];
    J[240] -= dqdci;              /* dwdot[NC7H16]/d[HO2] */
    J[247] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    J[248] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[264] += dqdci;              /* dwdot[C7H15-2]/d[HO2] */
    /* d()/d[C7H15-2] */
    dqdci =  - k_r*sc[7];
    J[720] -= dqdci;              /* dwdot[NC7H16]/d[C7H15-2] */
    J[727] += dqdci;              /* dwdot[H2O2]/d[C7H15-2] */
    J[728] -= dqdci;              /* dwdot[HO2]/d[C7H15-2] */
    J[744] += dqdci;              /* dwdot[C7H15-2]/d[C7H15-2] */
    /* d()/dT */
    J[870] -= dqdT;               /* dwdot[NC7H16]/dT */
    J[877] += dqdT;               /* dwdot[H2O2]/dT */
    J[878] -= dqdT;               /* dwdot[HO2]/dT */
    J[894] += dqdT;               /* dwdot[C7H15-2]/dT */

    /*reaction 7: NC7H16 + O2 <=> HO2 + C7H15-2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[1];
    k_f = prefactor_units[6] * fwd_A[6]
                * exp(fwd_beta[6] * tc[0] - activation_units[6] * fwd_Ea[6] * invT);
    dlnkfdT = fwd_beta[6] * invT + activation_units[6] * fwd_Ea[6] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[24];
    Kc = exp(g_RT[0] + g_RT[1] - g_RT[8] - g_RT[24]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[1]) + (h_RT[8] + h_RT[24]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* NC7H16 */
    wdot[1] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    wdot[24] += q; /* C7H15-2 */
    /* d()/d[NC7H16] */
    dqdci =  + k_f*sc[1];
    J[0] -= dqdci;                /* dwdot[NC7H16]/d[NC7H16] */
    J[1] -= dqdci;                /* dwdot[O2]/d[NC7H16] */
    J[8] += dqdci;                /* dwdot[HO2]/d[NC7H16] */
    J[24] += dqdci;               /* dwdot[C7H15-2]/d[NC7H16] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[0];
    J[30] -= dqdci;               /* dwdot[NC7H16]/d[O2] */
    J[31] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[38] += dqdci;               /* dwdot[HO2]/d[O2] */
    J[54] += dqdci;               /* dwdot[C7H15-2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[24];
    J[240] -= dqdci;              /* dwdot[NC7H16]/d[HO2] */
    J[241] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[248] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[264] += dqdci;              /* dwdot[C7H15-2]/d[HO2] */
    /* d()/d[C7H15-2] */
    dqdci =  - k_r*sc[8];
    J[720] -= dqdci;              /* dwdot[NC7H16]/d[C7H15-2] */
    J[721] -= dqdci;              /* dwdot[O2]/d[C7H15-2] */
    J[728] += dqdci;              /* dwdot[HO2]/d[C7H15-2] */
    J[744] += dqdci;              /* dwdot[C7H15-2]/d[C7H15-2] */
    /* d()/dT */
    J[870] -= dqdT;               /* dwdot[NC7H16]/dT */
    J[871] -= dqdT;               /* dwdot[O2]/dT */
    J[878] += dqdT;               /* dwdot[HO2]/dT */
    J[894] += dqdT;               /* dwdot[C7H15-2]/dT */

    /*reaction 8: O2 + C7H15-2 <=> C7H15O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[24];
    k_f = prefactor_units[7] * fwd_A[7]
                * exp(fwd_beta[7] * tc[0] - activation_units[7] * fwd_Ea[7] * invT);
    dlnkfdT = fwd_beta[7] * invT + activation_units[7] * fwd_Ea[7] * invT2;
    /* reverse */
    phi_r = sc[25];
    Kc = refCinv * exp(g_RT[1] + g_RT[24] - g_RT[25]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[24]) + (h_RT[25]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* O2 */
    wdot[24] -= q; /* C7H15-2 */
    wdot[25] += q; /* C7H15O2 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[24];
    J[31] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[54] -= dqdci;               /* dwdot[C7H15-2]/d[O2] */
    J[55] += dqdci;               /* dwdot[C7H15O2]/d[O2] */
    /* d()/d[C7H15-2] */
    dqdci =  + k_f*sc[1];
    J[721] -= dqdci;              /* dwdot[O2]/d[C7H15-2] */
    J[744] -= dqdci;              /* dwdot[C7H15-2]/d[C7H15-2] */
    J[745] += dqdci;              /* dwdot[C7H15O2]/d[C7H15-2] */
    /* d()/d[C7H15O2] */
    dqdci =  - k_r;
    J[751] -= dqdci;              /* dwdot[O2]/d[C7H15O2] */
    J[774] -= dqdci;              /* dwdot[C7H15-2]/d[C7H15O2] */
    J[775] += dqdci;              /* dwdot[C7H15O2]/d[C7H15O2] */
    /* d()/dT */
    J[871] -= dqdT;               /* dwdot[O2]/dT */
    J[894] -= dqdT;               /* dwdot[C7H15-2]/dT */
    J[895] += dqdT;               /* dwdot[C7H15O2]/dT */

    /*reaction 9: O2 + C7H15O2 <=> OH + C7KET12 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[25];
    k_f = prefactor_units[8] * fwd_A[8]
                * exp(fwd_beta[8] * tc[0] - activation_units[8] * fwd_Ea[8] * invT);
    dlnkfdT = fwd_beta[8] * invT + activation_units[8] * fwd_Ea[8] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[26];
    Kc = exp(g_RT[1] - g_RT[6] + g_RT[25] - g_RT[26]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[25]) + (h_RT[6] + h_RT[26]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* O2 */
    wdot[6] += q; /* OH */
    wdot[25] -= q; /* C7H15O2 */
    wdot[26] += q; /* C7KET12 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[25];
    J[31] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[36] += dqdci;               /* dwdot[OH]/d[O2] */
    J[55] -= dqdci;               /* dwdot[C7H15O2]/d[O2] */
    J[56] += dqdci;               /* dwdot[C7KET12]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[26];
    J[181] -= dqdci;              /* dwdot[O2]/d[OH] */
    J[186] += dqdci;              /* dwdot[OH]/d[OH] */
    J[205] -= dqdci;              /* dwdot[C7H15O2]/d[OH] */
    J[206] += dqdci;              /* dwdot[C7KET12]/d[OH] */
    /* d()/d[C7H15O2] */
    dqdci =  + k_f*sc[1];
    J[751] -= dqdci;              /* dwdot[O2]/d[C7H15O2] */
    J[756] += dqdci;              /* dwdot[OH]/d[C7H15O2] */
    J[775] -= dqdci;              /* dwdot[C7H15O2]/d[C7H15O2] */
    J[776] += dqdci;              /* dwdot[C7KET12]/d[C7H15O2] */
    /* d()/d[C7KET12] */
    dqdci =  - k_r*sc[6];
    J[781] -= dqdci;              /* dwdot[O2]/d[C7KET12] */
    J[786] += dqdci;              /* dwdot[OH]/d[C7KET12] */
    J[805] -= dqdci;              /* dwdot[C7H15O2]/d[C7KET12] */
    J[806] += dqdci;              /* dwdot[C7KET12]/d[C7KET12] */
    /* d()/dT */
    J[871] -= dqdT;               /* dwdot[O2]/dT */
    J[876] += dqdT;               /* dwdot[OH]/dT */
    J[895] -= dqdT;               /* dwdot[C7H15O2]/dT */
    J[896] += dqdT;               /* dwdot[C7KET12]/dT */

    /*reaction 10: C7KET12 <=> OH + CH2O + C5H11CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[26];
    k_f = prefactor_units[9] * fwd_A[9]
                * exp(fwd_beta[9] * tc[0] - activation_units[9] * fwd_Ea[9] * invT);
    dlnkfdT = fwd_beta[9] * invT + activation_units[9] * fwd_Ea[9] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[12]*sc[27];
    Kc = refC*refC * exp(-g_RT[6] - g_RT[12] + g_RT[26] - g_RT[27]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[26]) + (h_RT[6] + h_RT[12] + h_RT[27]) - 2);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[6] += q; /* OH */
    wdot[12] += q; /* CH2O */
    wdot[26] -= q; /* C7KET12 */
    wdot[27] += q; /* C5H11CO */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[12]*sc[27];
    J[186] += dqdci;              /* dwdot[OH]/d[OH] */
    J[192] += dqdci;              /* dwdot[CH2O]/d[OH] */
    J[206] -= dqdci;              /* dwdot[C7KET12]/d[OH] */
    J[207] += dqdci;              /* dwdot[C5H11CO]/d[OH] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[6]*sc[27];
    J[366] += dqdci;              /* dwdot[OH]/d[CH2O] */
    J[372] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[386] -= dqdci;              /* dwdot[C7KET12]/d[CH2O] */
    J[387] += dqdci;              /* dwdot[C5H11CO]/d[CH2O] */
    /* d()/d[C7KET12] */
    dqdci =  + k_f;
    J[786] += dqdci;              /* dwdot[OH]/d[C7KET12] */
    J[792] += dqdci;              /* dwdot[CH2O]/d[C7KET12] */
    J[806] -= dqdci;              /* dwdot[C7KET12]/d[C7KET12] */
    J[807] += dqdci;              /* dwdot[C5H11CO]/d[C7KET12] */
    /* d()/d[C5H11CO] */
    dqdci =  - k_r*sc[6]*sc[12];
    J[816] += dqdci;              /* dwdot[OH]/d[C5H11CO] */
    J[822] += dqdci;              /* dwdot[CH2O]/d[C5H11CO] */
    J[836] -= dqdci;              /* dwdot[C7KET12]/d[C5H11CO] */
    J[837] += dqdci;              /* dwdot[C5H11CO]/d[C5H11CO] */
    /* d()/dT */
    J[876] += dqdT;               /* dwdot[OH]/dT */
    J[882] += dqdT;               /* dwdot[CH2O]/dT */
    J[896] -= dqdT;               /* dwdot[C7KET12]/dT */
    J[897] += dqdT;               /* dwdot[C5H11CO]/dT */

    /*reaction 11: C5H11CO <=> CO + C2H4 + C3H7 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[27];
    k_f = prefactor_units[10] * fwd_A[10]
                * exp(fwd_beta[10] * tc[0] - activation_units[10] * fwd_Ea[10] * invT);
    dlnkfdT = fwd_beta[10] * invT + activation_units[10] * fwd_Ea[10] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[18]*sc[23];
    Kc = refC*refC * exp(-g_RT[4] - g_RT[18] - g_RT[23] + g_RT[27]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[27]) + (h_RT[4] + h_RT[18] + h_RT[23]) - 2);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* CO */
    wdot[18] += q; /* C2H4 */
    wdot[23] += q; /* C3H7 */
    wdot[27] -= q; /* C5H11CO */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[18]*sc[23];
    J[124] += dqdci;              /* dwdot[CO]/d[CO] */
    J[138] += dqdci;              /* dwdot[C2H4]/d[CO] */
    J[143] += dqdci;              /* dwdot[C3H7]/d[CO] */
    J[147] -= dqdci;              /* dwdot[C5H11CO]/d[CO] */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[4]*sc[23];
    J[544] += dqdci;              /* dwdot[CO]/d[C2H4] */
    J[558] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[563] += dqdci;              /* dwdot[C3H7]/d[C2H4] */
    J[567] -= dqdci;              /* dwdot[C5H11CO]/d[C2H4] */
    /* d()/d[C3H7] */
    dqdci =  - k_r*sc[4]*sc[18];
    J[694] += dqdci;              /* dwdot[CO]/d[C3H7] */
    J[708] += dqdci;              /* dwdot[C2H4]/d[C3H7] */
    J[713] += dqdci;              /* dwdot[C3H7]/d[C3H7] */
    J[717] -= dqdci;              /* dwdot[C5H11CO]/d[C3H7] */
    /* d()/d[C5H11CO] */
    dqdci =  + k_f;
    J[814] += dqdci;              /* dwdot[CO]/d[C5H11CO] */
    J[828] += dqdci;              /* dwdot[C2H4]/d[C5H11CO] */
    J[833] += dqdci;              /* dwdot[C3H7]/d[C5H11CO] */
    J[837] -= dqdci;              /* dwdot[C5H11CO]/d[C5H11CO] */
    /* d()/dT */
    J[874] += dqdT;               /* dwdot[CO]/dT */
    J[888] += dqdT;               /* dwdot[C2H4]/dT */
    J[893] += dqdT;               /* dwdot[C3H7]/dT */
    J[897] -= dqdT;               /* dwdot[C5H11CO]/dT */

    /*reaction 12: C7H15-2 <=> C2H4 + C2H5 + C3H6 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[24];
    k_f = prefactor_units[11] * fwd_A[11]
                * exp(fwd_beta[11] * tc[0] - activation_units[11] * fwd_Ea[11] * invT);
    dlnkfdT = fwd_beta[11] * invT + activation_units[11] * fwd_Ea[11] * invT2;
    /* reverse */
    phi_r = sc[18]*sc[19]*sc[22];
    Kc = refC*refC * exp(-g_RT[18] - g_RT[19] - g_RT[22] + g_RT[24]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[24]) + (h_RT[18] + h_RT[19] + h_RT[22]) - 2);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[18] += q; /* C2H4 */
    wdot[19] += q; /* C2H5 */
    wdot[22] += q; /* C3H6 */
    wdot[24] -= q; /* C7H15-2 */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[19]*sc[22];
    J[558] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[559] += dqdci;              /* dwdot[C2H5]/d[C2H4] */
    J[562] += dqdci;              /* dwdot[C3H6]/d[C2H4] */
    J[564] -= dqdci;              /* dwdot[C7H15-2]/d[C2H4] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[18]*sc[22];
    J[588] += dqdci;              /* dwdot[C2H4]/d[C2H5] */
    J[589] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[592] += dqdci;              /* dwdot[C3H6]/d[C2H5] */
    J[594] -= dqdci;              /* dwdot[C7H15-2]/d[C2H5] */
    /* d()/d[C3H6] */
    dqdci =  - k_r*sc[18]*sc[19];
    J[678] += dqdci;              /* dwdot[C2H4]/d[C3H6] */
    J[679] += dqdci;              /* dwdot[C2H5]/d[C3H6] */
    J[682] += dqdci;              /* dwdot[C3H6]/d[C3H6] */
    J[684] -= dqdci;              /* dwdot[C7H15-2]/d[C3H6] */
    /* d()/d[C7H15-2] */
    dqdci =  + k_f;
    J[738] += dqdci;              /* dwdot[C2H4]/d[C7H15-2] */
    J[739] += dqdci;              /* dwdot[C2H5]/d[C7H15-2] */
    J[742] += dqdci;              /* dwdot[C3H6]/d[C7H15-2] */
    J[744] -= dqdci;              /* dwdot[C7H15-2]/d[C7H15-2] */
    /* d()/dT */
    J[888] += dqdT;               /* dwdot[C2H4]/dT */
    J[889] += dqdT;               /* dwdot[C2H5]/dT */
    J[892] += dqdT;               /* dwdot[C3H6]/dT */
    J[894] -= dqdT;               /* dwdot[C7H15-2]/dT */

    /*reaction 13: C3H7 <=> CH3 + C2H4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[23];
    k_f = prefactor_units[12] * fwd_A[12]
                * exp(fwd_beta[12] * tc[0] - activation_units[12] * fwd_Ea[12] * invT);
    dlnkfdT = fwd_beta[12] * invT + activation_units[12] * fwd_Ea[12] * invT2;
    /* reverse */
    phi_r = sc[15]*sc[18];
    Kc = refC * exp(-g_RT[15] - g_RT[18] + g_RT[23]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[23]) + (h_RT[15] + h_RT[18]) - 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[15] += q; /* CH3 */
    wdot[18] += q; /* C2H4 */
    wdot[23] -= q; /* C3H7 */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[18];
    J[465] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[468] += dqdci;              /* dwdot[C2H4]/d[CH3] */
    J[473] -= dqdci;              /* dwdot[C3H7]/d[CH3] */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[15];
    J[555] += dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[558] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[563] -= dqdci;              /* dwdot[C3H7]/d[C2H4] */
    /* d()/d[C3H7] */
    dqdci =  + k_f;
    J[705] += dqdci;              /* dwdot[CH3]/d[C3H7] */
    J[708] += dqdci;              /* dwdot[C2H4]/d[C3H7] */
    J[713] -= dqdci;              /* dwdot[C3H7]/d[C3H7] */
    /* d()/dT */
    J[885] += dqdT;               /* dwdot[CH3]/dT */
    J[888] += dqdT;               /* dwdot[C2H4]/dT */
    J[893] -= dqdT;               /* dwdot[C3H7]/dT */

    /*reaction 14: C3H7 <=> H + C3H6 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[23];
    k_f = prefactor_units[13] * fwd_A[13]
                * exp(fwd_beta[13] * tc[0] - activation_units[13] * fwd_Ea[13] * invT);
    dlnkfdT = fwd_beta[13] * invT + activation_units[13] * fwd_Ea[13] * invT2;
    /* reverse */
    phi_r = sc[9]*sc[22];
    Kc = refC * exp(-g_RT[9] - g_RT[22] + g_RT[23]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[23]) + (h_RT[9] + h_RT[22]) - 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[9] += q; /* H */
    wdot[22] += q; /* C3H6 */
    wdot[23] -= q; /* C3H7 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[22];
    J[279] += dqdci;              /* dwdot[H]/d[H] */
    J[292] += dqdci;              /* dwdot[C3H6]/d[H] */
    J[293] -= dqdci;              /* dwdot[C3H7]/d[H] */
    /* d()/d[C3H6] */
    dqdci =  - k_r*sc[9];
    J[669] += dqdci;              /* dwdot[H]/d[C3H6] */
    J[682] += dqdci;              /* dwdot[C3H6]/d[C3H6] */
    J[683] -= dqdci;              /* dwdot[C3H7]/d[C3H6] */
    /* d()/d[C3H7] */
    dqdci =  + k_f;
    J[699] += dqdci;              /* dwdot[H]/d[C3H7] */
    J[712] += dqdci;              /* dwdot[C3H6]/d[C3H7] */
    J[713] -= dqdci;              /* dwdot[C3H7]/d[C3H7] */
    /* d()/dT */
    J[879] += dqdT;               /* dwdot[H]/dT */
    J[892] += dqdT;               /* dwdot[C3H6]/dT */
    J[893] -= dqdT;               /* dwdot[C3H7]/dT */

    /*reaction 15: CH3 + C3H6 <=> CH4 + C3H5 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[15]*sc[22];
    k_f = prefactor_units[14] * fwd_A[14]
                * exp(fwd_beta[14] * tc[0] - activation_units[14] * fwd_Ea[14] * invT);
    dlnkfdT = fwd_beta[14] * invT + activation_units[14] * fwd_Ea[14] * invT2;
    /* reverse */
    phi_r = sc[16]*sc[21];
    Kc = exp(g_RT[15] - g_RT[16] - g_RT[21] + g_RT[22]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[15] + h_RT[22]) + (h_RT[16] + h_RT[21]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[15] -= q; /* CH3 */
    wdot[16] += q; /* CH4 */
    wdot[21] += q; /* C3H5 */
    wdot[22] -= q; /* C3H6 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[22];
    J[465] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[466] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[471] += dqdci;              /* dwdot[C3H5]/d[CH3] */
    J[472] -= dqdci;              /* dwdot[C3H6]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[21];
    J[495] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[496] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[501] += dqdci;              /* dwdot[C3H5]/d[CH4] */
    J[502] -= dqdci;              /* dwdot[C3H6]/d[CH4] */
    /* d()/d[C3H5] */
    dqdci =  - k_r*sc[16];
    J[645] -= dqdci;              /* dwdot[CH3]/d[C3H5] */
    J[646] += dqdci;              /* dwdot[CH4]/d[C3H5] */
    J[651] += dqdci;              /* dwdot[C3H5]/d[C3H5] */
    J[652] -= dqdci;              /* dwdot[C3H6]/d[C3H5] */
    /* d()/d[C3H6] */
    dqdci =  + k_f*sc[15];
    J[675] -= dqdci;              /* dwdot[CH3]/d[C3H6] */
    J[676] += dqdci;              /* dwdot[CH4]/d[C3H6] */
    J[681] += dqdci;              /* dwdot[C3H5]/d[C3H6] */
    J[682] -= dqdci;              /* dwdot[C3H6]/d[C3H6] */
    /* d()/dT */
    J[885] -= dqdT;               /* dwdot[CH3]/dT */
    J[886] += dqdT;               /* dwdot[CH4]/dT */
    J[891] += dqdT;               /* dwdot[C3H5]/dT */
    J[892] -= dqdT;               /* dwdot[C3H6]/dT */

    /*reaction 16: O2 + C3H5 <=> HO2 + C3H4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[21];
    k_f = prefactor_units[15] * fwd_A[15]
                * exp(fwd_beta[15] * tc[0] - activation_units[15] * fwd_Ea[15] * invT);
    dlnkfdT = fwd_beta[15] * invT + activation_units[15] * fwd_Ea[15] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[20];
    Kc = exp(g_RT[1] - g_RT[8] - g_RT[20] + g_RT[21]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[21]) + (h_RT[8] + h_RT[20]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    wdot[20] += q; /* C3H4 */
    wdot[21] -= q; /* C3H5 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[21];
    J[31] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[38] += dqdci;               /* dwdot[HO2]/d[O2] */
    J[50] += dqdci;               /* dwdot[C3H4]/d[O2] */
    J[51] -= dqdci;               /* dwdot[C3H5]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[20];
    J[241] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[248] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[260] += dqdci;              /* dwdot[C3H4]/d[HO2] */
    J[261] -= dqdci;              /* dwdot[C3H5]/d[HO2] */
    /* d()/d[C3H4] */
    dqdci =  - k_r*sc[8];
    J[601] -= dqdci;              /* dwdot[O2]/d[C3H4] */
    J[608] += dqdci;              /* dwdot[HO2]/d[C3H4] */
    J[620] += dqdci;              /* dwdot[C3H4]/d[C3H4] */
    J[621] -= dqdci;              /* dwdot[C3H5]/d[C3H4] */
    /* d()/d[C3H5] */
    dqdci =  + k_f*sc[1];
    J[631] -= dqdci;              /* dwdot[O2]/d[C3H5] */
    J[638] += dqdci;              /* dwdot[HO2]/d[C3H5] */
    J[650] += dqdci;              /* dwdot[C3H4]/d[C3H5] */
    J[651] -= dqdci;              /* dwdot[C3H5]/d[C3H5] */
    /* d()/dT */
    J[871] -= dqdT;               /* dwdot[O2]/dT */
    J[878] += dqdT;               /* dwdot[HO2]/dT */
    J[890] += dqdT;               /* dwdot[C3H4]/dT */
    J[891] -= dqdT;               /* dwdot[C3H5]/dT */

    /*reaction 17: OH + C3H4 <=> CH2O + C2H3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[20];
    k_f = prefactor_units[16] * fwd_A[16]
                * exp(fwd_beta[16] * tc[0] - activation_units[16] * fwd_Ea[16] * invT);
    dlnkfdT = fwd_beta[16] * invT + activation_units[16] * fwd_Ea[16] * invT2;
    /* reverse */
    phi_r = sc[12]*sc[17];
    Kc = exp(g_RT[6] - g_RT[12] - g_RT[17] + g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[20]) + (h_RT[12] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[6] -= q; /* OH */
    wdot[12] += q; /* CH2O */
    wdot[17] += q; /* C2H3 */
    wdot[20] -= q; /* C3H4 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[20];
    J[186] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[192] += dqdci;              /* dwdot[CH2O]/d[OH] */
    J[197] += dqdci;              /* dwdot[C2H3]/d[OH] */
    J[200] -= dqdci;              /* dwdot[C3H4]/d[OH] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[17];
    J[366] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[372] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[377] += dqdci;              /* dwdot[C2H3]/d[CH2O] */
    J[380] -= dqdci;              /* dwdot[C3H4]/d[CH2O] */
    /* d()/d[C2H3] */
    dqdci =  - k_r*sc[12];
    J[516] -= dqdci;              /* dwdot[OH]/d[C2H3] */
    J[522] += dqdci;              /* dwdot[CH2O]/d[C2H3] */
    J[527] += dqdci;              /* dwdot[C2H3]/d[C2H3] */
    J[530] -= dqdci;              /* dwdot[C3H4]/d[C2H3] */
    /* d()/d[C3H4] */
    dqdci =  + k_f*sc[6];
    J[606] -= dqdci;              /* dwdot[OH]/d[C3H4] */
    J[612] += dqdci;              /* dwdot[CH2O]/d[C3H4] */
    J[617] += dqdci;              /* dwdot[C2H3]/d[C3H4] */
    J[620] -= dqdci;              /* dwdot[C3H4]/d[C3H4] */
    /* d()/dT */
    J[876] -= dqdT;               /* dwdot[OH]/dT */
    J[882] += dqdT;               /* dwdot[CH2O]/dT */
    J[887] += dqdT;               /* dwdot[C2H3]/dT */
    J[890] -= dqdT;               /* dwdot[C3H4]/dT */

    /*reaction 18: OH + C3H4 <=> HCO + C2H4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[20];
    k_f = prefactor_units[17] * fwd_A[17]
                * exp(fwd_beta[17] * tc[0] - activation_units[17] * fwd_Ea[17] * invT);
    dlnkfdT = fwd_beta[17] * invT + activation_units[17] * fwd_Ea[17] * invT2;
    /* reverse */
    phi_r = sc[13]*sc[18];
    Kc = exp(g_RT[6] - g_RT[13] - g_RT[18] + g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[20]) + (h_RT[13] + h_RT[18]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[6] -= q; /* OH */
    wdot[13] += q; /* HCO */
    wdot[18] += q; /* C2H4 */
    wdot[20] -= q; /* C3H4 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[20];
    J[186] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[193] += dqdci;              /* dwdot[HCO]/d[OH] */
    J[198] += dqdci;              /* dwdot[C2H4]/d[OH] */
    J[200] -= dqdci;              /* dwdot[C3H4]/d[OH] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[18];
    J[396] -= dqdci;              /* dwdot[OH]/d[HCO] */
    J[403] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[408] += dqdci;              /* dwdot[C2H4]/d[HCO] */
    J[410] -= dqdci;              /* dwdot[C3H4]/d[HCO] */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[13];
    J[546] -= dqdci;              /* dwdot[OH]/d[C2H4] */
    J[553] += dqdci;              /* dwdot[HCO]/d[C2H4] */
    J[558] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[560] -= dqdci;              /* dwdot[C3H4]/d[C2H4] */
    /* d()/d[C3H4] */
    dqdci =  + k_f*sc[6];
    J[606] -= dqdci;              /* dwdot[OH]/d[C3H4] */
    J[613] += dqdci;              /* dwdot[HCO]/d[C3H4] */
    J[618] += dqdci;              /* dwdot[C2H4]/d[C3H4] */
    J[620] -= dqdci;              /* dwdot[C3H4]/d[C3H4] */
    /* d()/dT */
    J[876] -= dqdT;               /* dwdot[OH]/dT */
    J[883] += dqdT;               /* dwdot[HCO]/dT */
    J[888] += dqdT;               /* dwdot[C2H4]/dT */
    J[890] -= dqdT;               /* dwdot[C3H4]/dT */

    /*reaction 19: HO2 + CH3 <=> OH + CH3O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[15];
    k_f = prefactor_units[18] * fwd_A[18]
                * exp(fwd_beta[18] * tc[0] - activation_units[18] * fwd_Ea[18] * invT);
    dlnkfdT = fwd_beta[18] * invT + activation_units[18] * fwd_Ea[18] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[11];
    Kc = exp(-g_RT[6] + g_RT[8] - g_RT[11] + g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[15]) + (h_RT[6] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[6] += q; /* OH */
    wdot[8] -= q; /* HO2 */
    wdot[11] += q; /* CH3O */
    wdot[15] -= q; /* CH3 */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[11];
    J[186] += dqdci;              /* dwdot[OH]/d[OH] */
    J[188] -= dqdci;              /* dwdot[HO2]/d[OH] */
    J[191] += dqdci;              /* dwdot[CH3O]/d[OH] */
    J[195] -= dqdci;              /* dwdot[CH3]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[15];
    J[246] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[248] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[251] += dqdci;              /* dwdot[CH3O]/d[HO2] */
    J[255] -= dqdci;              /* dwdot[CH3]/d[HO2] */
    /* d()/d[CH3O] */
    dqdci =  - k_r*sc[6];
    J[336] += dqdci;              /* dwdot[OH]/d[CH3O] */
    J[338] -= dqdci;              /* dwdot[HO2]/d[CH3O] */
    J[341] += dqdci;              /* dwdot[CH3O]/d[CH3O] */
    J[345] -= dqdci;              /* dwdot[CH3]/d[CH3O] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[8];
    J[456] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[458] -= dqdci;              /* dwdot[HO2]/d[CH3] */
    J[461] += dqdci;              /* dwdot[CH3O]/d[CH3] */
    J[465] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[876] += dqdT;               /* dwdot[OH]/dT */
    J[878] -= dqdT;               /* dwdot[HO2]/dT */
    J[881] += dqdT;               /* dwdot[CH3O]/dT */
    J[885] -= dqdT;               /* dwdot[CH3]/dT */

    /*reaction 20: OH + CH3 <=> H2O + CH2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[15];
    k_f = prefactor_units[19] * fwd_A[19]
                * exp(fwd_beta[19] * tc[0] - activation_units[19] * fwd_Ea[19] * invT);
    dlnkfdT = fwd_beta[19] * invT + activation_units[19] * fwd_Ea[19] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[14];
    Kc = exp(-g_RT[3] + g_RT[6] - g_RT[14] + g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[15]) + (h_RT[3] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* H2O */
    wdot[6] -= q; /* OH */
    wdot[14] += q; /* CH2 */
    wdot[15] -= q; /* CH3 */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[14];
    J[93] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[96] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[104] += dqdci;              /* dwdot[CH2]/d[H2O] */
    J[105] -= dqdci;              /* dwdot[CH3]/d[H2O] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[15];
    J[183] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[186] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[194] += dqdci;              /* dwdot[CH2]/d[OH] */
    J[195] -= dqdci;              /* dwdot[CH3]/d[OH] */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[3];
    J[423] += dqdci;              /* dwdot[H2O]/d[CH2] */
    J[426] -= dqdci;              /* dwdot[OH]/d[CH2] */
    J[434] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[435] -= dqdci;              /* dwdot[CH3]/d[CH2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[6];
    J[453] += dqdci;              /* dwdot[H2O]/d[CH3] */
    J[456] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[464] += dqdci;              /* dwdot[CH2]/d[CH3] */
    J[465] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[873] += dqdT;               /* dwdot[H2O]/dT */
    J[876] -= dqdT;               /* dwdot[OH]/dT */
    J[884] += dqdT;               /* dwdot[CH2]/dT */
    J[885] -= dqdT;               /* dwdot[CH3]/dT */

    /*reaction 21: OH + CH2 <=> H + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[14];
    k_f = prefactor_units[20] * fwd_A[20]
                * exp(fwd_beta[20] * tc[0] - activation_units[20] * fwd_Ea[20] * invT);
    dlnkfdT = fwd_beta[20] * invT + activation_units[20] * fwd_Ea[20] * invT2;
    /* reverse */
    phi_r = sc[9]*sc[12];
    Kc = exp(g_RT[6] - g_RT[9] - g_RT[12] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[14]) + (h_RT[9] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[6] -= q; /* OH */
    wdot[9] += q; /* H */
    wdot[12] += q; /* CH2O */
    wdot[14] -= q; /* CH2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[14];
    J[186] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[189] += dqdci;              /* dwdot[H]/d[OH] */
    J[192] += dqdci;              /* dwdot[CH2O]/d[OH] */
    J[194] -= dqdci;              /* dwdot[CH2]/d[OH] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[12];
    J[276] -= dqdci;              /* dwdot[OH]/d[H] */
    J[279] += dqdci;              /* dwdot[H]/d[H] */
    J[282] += dqdci;              /* dwdot[CH2O]/d[H] */
    J[284] -= dqdci;              /* dwdot[CH2]/d[H] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[9];
    J[366] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[369] += dqdci;              /* dwdot[H]/d[CH2O] */
    J[372] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[374] -= dqdci;              /* dwdot[CH2]/d[CH2O] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[6];
    J[426] -= dqdci;              /* dwdot[OH]/d[CH2] */
    J[429] += dqdci;              /* dwdot[H]/d[CH2] */
    J[432] += dqdci;              /* dwdot[CH2O]/d[CH2] */
    J[434] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    /* d()/dT */
    J[876] -= dqdT;               /* dwdot[OH]/dT */
    J[879] += dqdT;               /* dwdot[H]/dT */
    J[882] += dqdT;               /* dwdot[CH2O]/dT */
    J[884] -= dqdT;               /* dwdot[CH2]/dT */

    /*reaction 22: O2 + CH2 <=> OH + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[14];
    k_f = prefactor_units[21] * fwd_A[21]
                * exp(fwd_beta[21] * tc[0] - activation_units[21] * fwd_Ea[21] * invT);
    dlnkfdT = fwd_beta[21] * invT + activation_units[21] * fwd_Ea[21] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[13];
    Kc = exp(g_RT[1] - g_RT[6] - g_RT[13] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[14]) + (h_RT[6] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* O2 */
    wdot[6] += q; /* OH */
    wdot[13] += q; /* HCO */
    wdot[14] -= q; /* CH2 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[14];
    J[31] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[36] += dqdci;               /* dwdot[OH]/d[O2] */
    J[43] += dqdci;               /* dwdot[HCO]/d[O2] */
    J[44] -= dqdci;               /* dwdot[CH2]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[13];
    J[181] -= dqdci;              /* dwdot[O2]/d[OH] */
    J[186] += dqdci;              /* dwdot[OH]/d[OH] */
    J[193] += dqdci;              /* dwdot[HCO]/d[OH] */
    J[194] -= dqdci;              /* dwdot[CH2]/d[OH] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[6];
    J[391] -= dqdci;              /* dwdot[O2]/d[HCO] */
    J[396] += dqdci;              /* dwdot[OH]/d[HCO] */
    J[403] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[404] -= dqdci;              /* dwdot[CH2]/d[HCO] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[1];
    J[421] -= dqdci;              /* dwdot[O2]/d[CH2] */
    J[426] += dqdci;              /* dwdot[OH]/d[CH2] */
    J[433] += dqdci;              /* dwdot[HCO]/d[CH2] */
    J[434] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    /* d()/dT */
    J[871] -= dqdT;               /* dwdot[O2]/dT */
    J[876] += dqdT;               /* dwdot[OH]/dT */
    J[883] += dqdT;               /* dwdot[HCO]/dT */
    J[884] -= dqdT;               /* dwdot[CH2]/dT */

    /*reaction 23: O2 + CH2 <=> CO2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[14];
    k_f = prefactor_units[22] * fwd_A[22]
                * exp(fwd_beta[22] * tc[0] - activation_units[22] * fwd_Ea[22] * invT);
    dlnkfdT = fwd_beta[22] * invT + activation_units[22] * fwd_Ea[22] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[5];
    Kc = exp(g_RT[1] - g_RT[2] - g_RT[5] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[14]) + (h_RT[2] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* O2 */
    wdot[2] += q; /* CO2 */
    wdot[5] += q; /* H2 */
    wdot[14] -= q; /* CH2 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[14];
    J[31] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[32] += dqdci;               /* dwdot[CO2]/d[O2] */
    J[35] += dqdci;               /* dwdot[H2]/d[O2] */
    J[44] -= dqdci;               /* dwdot[CH2]/d[O2] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[5];
    J[61] -= dqdci;               /* dwdot[O2]/d[CO2] */
    J[62] += dqdci;               /* dwdot[CO2]/d[CO2] */
    J[65] += dqdci;               /* dwdot[H2]/d[CO2] */
    J[74] -= dqdci;               /* dwdot[CH2]/d[CO2] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[2];
    J[151] -= dqdci;              /* dwdot[O2]/d[H2] */
    J[152] += dqdci;              /* dwdot[CO2]/d[H2] */
    J[155] += dqdci;              /* dwdot[H2]/d[H2] */
    J[164] -= dqdci;              /* dwdot[CH2]/d[H2] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[1];
    J[421] -= dqdci;              /* dwdot[O2]/d[CH2] */
    J[422] += dqdci;              /* dwdot[CO2]/d[CH2] */
    J[425] += dqdci;              /* dwdot[H2]/d[CH2] */
    J[434] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    /* d()/dT */
    J[871] -= dqdT;               /* dwdot[O2]/dT */
    J[872] += dqdT;               /* dwdot[CO2]/dT */
    J[875] += dqdT;               /* dwdot[H2]/dT */
    J[884] -= dqdT;               /* dwdot[CH2]/dT */

    /*reaction 24: O2 + CH2 <=> H2O + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[14];
    k_f = prefactor_units[23] * fwd_A[23]
                * exp(fwd_beta[23] * tc[0] - activation_units[23] * fwd_Ea[23] * invT);
    dlnkfdT = fwd_beta[23] * invT + activation_units[23] * fwd_Ea[23] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[4];
    Kc = exp(g_RT[1] - g_RT[3] - g_RT[4] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[14]) + (h_RT[3] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* O2 */
    wdot[3] += q; /* H2O */
    wdot[4] += q; /* CO */
    wdot[14] -= q; /* CH2 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[14];
    J[31] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[33] += dqdci;               /* dwdot[H2O]/d[O2] */
    J[34] += dqdci;               /* dwdot[CO]/d[O2] */
    J[44] -= dqdci;               /* dwdot[CH2]/d[O2] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[4];
    J[91] -= dqdci;               /* dwdot[O2]/d[H2O] */
    J[93] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[94] += dqdci;               /* dwdot[CO]/d[H2O] */
    J[104] -= dqdci;              /* dwdot[CH2]/d[H2O] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[3];
    J[121] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[123] += dqdci;              /* dwdot[H2O]/d[CO] */
    J[124] += dqdci;              /* dwdot[CO]/d[CO] */
    J[134] -= dqdci;              /* dwdot[CH2]/d[CO] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[1];
    J[421] -= dqdci;              /* dwdot[O2]/d[CH2] */
    J[423] += dqdci;              /* dwdot[H2O]/d[CH2] */
    J[424] += dqdci;              /* dwdot[CO]/d[CH2] */
    J[434] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    /* d()/dT */
    J[871] -= dqdT;               /* dwdot[O2]/dT */
    J[873] += dqdT;               /* dwdot[H2O]/dT */
    J[874] += dqdT;               /* dwdot[CO]/dT */
    J[884] -= dqdT;               /* dwdot[CH2]/dT */

    /*reaction 25: O2 + CH2 <=> O + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[14];
    k_f = prefactor_units[24] * fwd_A[24]
                * exp(fwd_beta[24] * tc[0] - activation_units[24] * fwd_Ea[24] * invT);
    dlnkfdT = fwd_beta[24] * invT + activation_units[24] * fwd_Ea[24] * invT2;
    /* reverse */
    phi_r = sc[10]*sc[12];
    Kc = exp(g_RT[1] - g_RT[10] - g_RT[12] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[14]) + (h_RT[10] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* O2 */
    wdot[10] += q; /* O */
    wdot[12] += q; /* CH2O */
    wdot[14] -= q; /* CH2 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[14];
    J[31] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[40] += dqdci;               /* dwdot[O]/d[O2] */
    J[42] += dqdci;               /* dwdot[CH2O]/d[O2] */
    J[44] -= dqdci;               /* dwdot[CH2]/d[O2] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[12];
    J[301] -= dqdci;              /* dwdot[O2]/d[O] */
    J[310] += dqdci;              /* dwdot[O]/d[O] */
    J[312] += dqdci;              /* dwdot[CH2O]/d[O] */
    J[314] -= dqdci;              /* dwdot[CH2]/d[O] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[10];
    J[361] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    J[370] += dqdci;              /* dwdot[O]/d[CH2O] */
    J[372] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[374] -= dqdci;              /* dwdot[CH2]/d[CH2O] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[1];
    J[421] -= dqdci;              /* dwdot[O2]/d[CH2] */
    J[430] += dqdci;              /* dwdot[O]/d[CH2] */
    J[432] += dqdci;              /* dwdot[CH2O]/d[CH2] */
    J[434] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    /* d()/dT */
    J[871] -= dqdT;               /* dwdot[O2]/dT */
    J[880] += dqdT;               /* dwdot[O]/dT */
    J[882] += dqdT;               /* dwdot[CH2O]/dT */
    J[884] -= dqdT;               /* dwdot[CH2]/dT */

    /*reaction 26: O2 + CH2 <=> CO2 + 2 H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[14];
    k_f = prefactor_units[25] * fwd_A[25]
                * exp(fwd_beta[25] * tc[0] - activation_units[25] * fwd_Ea[25] * invT);
    dlnkfdT = fwd_beta[25] * invT + activation_units[25] * fwd_Ea[25] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[9]*sc[9];
    Kc = refC * exp(g_RT[1] - g_RT[2] - 2*g_RT[9] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[14]) + (h_RT[2] + 2*h_RT[9]) - 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* O2 */
    wdot[2] += q; /* CO2 */
    wdot[9] += 2 * q; /* H */
    wdot[14] -= q; /* CH2 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[14];
    J[31] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[32] += dqdci;               /* dwdot[CO2]/d[O2] */
    J[39] += 2 * dqdci;           /* dwdot[H]/d[O2] */
    J[44] -= dqdci;               /* dwdot[CH2]/d[O2] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[9]*sc[9];
    J[61] -= dqdci;               /* dwdot[O2]/d[CO2] */
    J[62] += dqdci;               /* dwdot[CO2]/d[CO2] */
    J[69] += 2 * dqdci;           /* dwdot[H]/d[CO2] */
    J[74] -= dqdci;               /* dwdot[CH2]/d[CO2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[2]*2*sc[9];
    J[271] -= dqdci;              /* dwdot[O2]/d[H] */
    J[272] += dqdci;              /* dwdot[CO2]/d[H] */
    J[279] += 2 * dqdci;          /* dwdot[H]/d[H] */
    J[284] -= dqdci;              /* dwdot[CH2]/d[H] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[1];
    J[421] -= dqdci;              /* dwdot[O2]/d[CH2] */
    J[422] += dqdci;              /* dwdot[CO2]/d[CH2] */
    J[429] += 2 * dqdci;          /* dwdot[H]/d[CH2] */
    J[434] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    /* d()/dT */
    J[871] -= dqdT;               /* dwdot[O2]/dT */
    J[872] += dqdT;               /* dwdot[CO2]/dT */
    J[879] += 2 * dqdT;           /* dwdot[H]/dT */
    J[884] -= dqdT;               /* dwdot[CH2]/dT */

    /*reaction 27: O2 + CH2 <=> CO + OH + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[14];
    k_f = prefactor_units[26] * fwd_A[26]
                * exp(fwd_beta[26] * tc[0] - activation_units[26] * fwd_Ea[26] * invT);
    dlnkfdT = fwd_beta[26] * invT + activation_units[26] * fwd_Ea[26] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[6]*sc[9];
    Kc = refC * exp(g_RT[1] - g_RT[4] - g_RT[6] - g_RT[9] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[14]) + (h_RT[4] + h_RT[6] + h_RT[9]) - 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* O2 */
    wdot[4] += q; /* CO */
    wdot[6] += q; /* OH */
    wdot[9] += q; /* H */
    wdot[14] -= q; /* CH2 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[14];
    J[31] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[34] += dqdci;               /* dwdot[CO]/d[O2] */
    J[36] += dqdci;               /* dwdot[OH]/d[O2] */
    J[39] += dqdci;               /* dwdot[H]/d[O2] */
    J[44] -= dqdci;               /* dwdot[CH2]/d[O2] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[6]*sc[9];
    J[121] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[124] += dqdci;              /* dwdot[CO]/d[CO] */
    J[126] += dqdci;              /* dwdot[OH]/d[CO] */
    J[129] += dqdci;              /* dwdot[H]/d[CO] */
    J[134] -= dqdci;              /* dwdot[CH2]/d[CO] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[4]*sc[9];
    J[181] -= dqdci;              /* dwdot[O2]/d[OH] */
    J[184] += dqdci;              /* dwdot[CO]/d[OH] */
    J[186] += dqdci;              /* dwdot[OH]/d[OH] */
    J[189] += dqdci;              /* dwdot[H]/d[OH] */
    J[194] -= dqdci;              /* dwdot[CH2]/d[OH] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[4]*sc[6];
    J[271] -= dqdci;              /* dwdot[O2]/d[H] */
    J[274] += dqdci;              /* dwdot[CO]/d[H] */
    J[276] += dqdci;              /* dwdot[OH]/d[H] */
    J[279] += dqdci;              /* dwdot[H]/d[H] */
    J[284] -= dqdci;              /* dwdot[CH2]/d[H] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[1];
    J[421] -= dqdci;              /* dwdot[O2]/d[CH2] */
    J[424] += dqdci;              /* dwdot[CO]/d[CH2] */
    J[426] += dqdci;              /* dwdot[OH]/d[CH2] */
    J[429] += dqdci;              /* dwdot[H]/d[CH2] */
    J[434] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    /* d()/dT */
    J[871] -= dqdT;               /* dwdot[O2]/dT */
    J[874] += dqdT;               /* dwdot[CO]/dT */
    J[876] += dqdT;               /* dwdot[OH]/dT */
    J[879] += dqdT;               /* dwdot[H]/dT */
    J[884] -= dqdT;               /* dwdot[CH2]/dT */

    /*reaction 28: CO + CH3O <=> CO2 + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[11];
    k_f = prefactor_units[27] * fwd_A[27]
                * exp(fwd_beta[27] * tc[0] - activation_units[27] * fwd_Ea[27] * invT);
    dlnkfdT = fwd_beta[27] * invT + activation_units[27] * fwd_Ea[27] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[15];
    Kc = exp(-g_RT[2] + g_RT[4] + g_RT[11] - g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[11]) + (h_RT[2] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* CO2 */
    wdot[4] -= q; /* CO */
    wdot[11] -= q; /* CH3O */
    wdot[15] += q; /* CH3 */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[15];
    J[62] += dqdci;               /* dwdot[CO2]/d[CO2] */
    J[64] -= dqdci;               /* dwdot[CO]/d[CO2] */
    J[71] -= dqdci;               /* dwdot[CH3O]/d[CO2] */
    J[75] += dqdci;               /* dwdot[CH3]/d[CO2] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[11];
    J[122] += dqdci;              /* dwdot[CO2]/d[CO] */
    J[124] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[131] -= dqdci;              /* dwdot[CH3O]/d[CO] */
    J[135] += dqdci;              /* dwdot[CH3]/d[CO] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[4];
    J[332] += dqdci;              /* dwdot[CO2]/d[CH3O] */
    J[334] -= dqdci;              /* dwdot[CO]/d[CH3O] */
    J[341] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    J[345] += dqdci;              /* dwdot[CH3]/d[CH3O] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[2];
    J[452] += dqdci;              /* dwdot[CO2]/d[CH3] */
    J[454] -= dqdci;              /* dwdot[CO]/d[CH3] */
    J[461] -= dqdci;              /* dwdot[CH3O]/d[CH3] */
    J[465] += dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[872] += dqdT;               /* dwdot[CO2]/dT */
    J[874] -= dqdT;               /* dwdot[CO]/dT */
    J[881] -= dqdT;               /* dwdot[CH3O]/dT */
    J[885] += dqdT;               /* dwdot[CH3]/dT */

    /*reaction 29: CO + OH <=> CO2 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[6];
    k_f = prefactor_units[28] * fwd_A[28]
                * exp(fwd_beta[28] * tc[0] - activation_units[28] * fwd_Ea[28] * invT);
    dlnkfdT = fwd_beta[28] * invT + activation_units[28] * fwd_Ea[28] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[9];
    Kc = exp(-g_RT[2] + g_RT[4] + g_RT[6] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[6]) + (h_RT[2] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* CO2 */
    wdot[4] -= q; /* CO */
    wdot[6] -= q; /* OH */
    wdot[9] += q; /* H */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[9];
    J[62] += dqdci;               /* dwdot[CO2]/d[CO2] */
    J[64] -= dqdci;               /* dwdot[CO]/d[CO2] */
    J[66] -= dqdci;               /* dwdot[OH]/d[CO2] */
    J[69] += dqdci;               /* dwdot[H]/d[CO2] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[6];
    J[122] += dqdci;              /* dwdot[CO2]/d[CO] */
    J[124] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[126] -= dqdci;              /* dwdot[OH]/d[CO] */
    J[129] += dqdci;              /* dwdot[H]/d[CO] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[4];
    J[182] += dqdci;              /* dwdot[CO2]/d[OH] */
    J[184] -= dqdci;              /* dwdot[CO]/d[OH] */
    J[186] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[189] += dqdci;              /* dwdot[H]/d[OH] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[2];
    J[272] += dqdci;              /* dwdot[CO2]/d[H] */
    J[274] -= dqdci;              /* dwdot[CO]/d[H] */
    J[276] -= dqdci;              /* dwdot[OH]/d[H] */
    J[279] += dqdci;              /* dwdot[H]/d[H] */
    /* d()/dT */
    J[872] += dqdT;               /* dwdot[CO2]/dT */
    J[874] -= dqdT;               /* dwdot[CO]/dT */
    J[876] -= dqdT;               /* dwdot[OH]/dT */
    J[879] += dqdT;               /* dwdot[H]/dT */

    /*reaction 30: OH + O <=> O2 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[10];
    k_f = prefactor_units[29] * fwd_A[29]
                * exp(fwd_beta[29] * tc[0] - activation_units[29] * fwd_Ea[29] * invT);
    dlnkfdT = fwd_beta[29] * invT + activation_units[29] * fwd_Ea[29] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[9];
    Kc = exp(-g_RT[1] + g_RT[6] - g_RT[9] + g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[10]) + (h_RT[1] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* O2 */
    wdot[6] -= q; /* OH */
    wdot[9] += q; /* H */
    wdot[10] -= q; /* O */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[9];
    J[31] += dqdci;               /* dwdot[O2]/d[O2] */
    J[36] -= dqdci;               /* dwdot[OH]/d[O2] */
    J[39] += dqdci;               /* dwdot[H]/d[O2] */
    J[40] -= dqdci;               /* dwdot[O]/d[O2] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[10];
    J[181] += dqdci;              /* dwdot[O2]/d[OH] */
    J[186] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[189] += dqdci;              /* dwdot[H]/d[OH] */
    J[190] -= dqdci;              /* dwdot[O]/d[OH] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[1];
    J[271] += dqdci;              /* dwdot[O2]/d[H] */
    J[276] -= dqdci;              /* dwdot[OH]/d[H] */
    J[279] += dqdci;              /* dwdot[H]/d[H] */
    J[280] -= dqdci;              /* dwdot[O]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[6];
    J[301] += dqdci;              /* dwdot[O2]/d[O] */
    J[306] -= dqdci;              /* dwdot[OH]/d[O] */
    J[309] += dqdci;              /* dwdot[H]/d[O] */
    J[310] -= dqdci;              /* dwdot[O]/d[O] */
    /* d()/dT */
    J[871] += dqdT;               /* dwdot[O2]/dT */
    J[876] -= dqdT;               /* dwdot[OH]/dT */
    J[879] += dqdT;               /* dwdot[H]/dT */
    J[880] -= dqdT;               /* dwdot[O]/dT */

    /*reaction 31: HO2 + H <=> 2 OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[9];
    k_f = prefactor_units[30] * fwd_A[30]
                * exp(fwd_beta[30] * tc[0] - activation_units[30] * fwd_Ea[30] * invT);
    dlnkfdT = fwd_beta[30] * invT + activation_units[30] * fwd_Ea[30] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[6];
    Kc = exp(-2*g_RT[6] + g_RT[8] + g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[9]) + (2*h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[6] += 2 * q; /* OH */
    wdot[8] -= q; /* HO2 */
    wdot[9] -= q; /* H */
    /* d()/d[OH] */
    dqdci =  - k_r*2*sc[6];
    J[186] += 2 * dqdci;          /* dwdot[OH]/d[OH] */
    J[188] -= dqdci;              /* dwdot[HO2]/d[OH] */
    J[189] -= dqdci;              /* dwdot[H]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[9];
    J[246] += 2 * dqdci;          /* dwdot[OH]/d[HO2] */
    J[248] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[249] -= dqdci;              /* dwdot[H]/d[HO2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[8];
    J[276] += 2 * dqdci;          /* dwdot[OH]/d[H] */
    J[278] -= dqdci;              /* dwdot[HO2]/d[H] */
    J[279] -= dqdci;              /* dwdot[H]/d[H] */
    /* d()/dT */
    J[876] += 2 * dqdT;           /* dwdot[OH]/dT */
    J[878] -= dqdT;               /* dwdot[HO2]/dT */
    J[879] -= dqdT;               /* dwdot[H]/dT */

    /*reaction 32: 2 OH <=> H2O + O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[6];
    k_f = prefactor_units[31] * fwd_A[31]
                * exp(fwd_beta[31] * tc[0] - activation_units[31] * fwd_Ea[31] * invT);
    dlnkfdT = fwd_beta[31] * invT + activation_units[31] * fwd_Ea[31] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[10];
    Kc = exp(-g_RT[3] + 2*g_RT[6] - g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[6]) + (h_RT[3] + h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* H2O */
    wdot[6] -= 2 * q; /* OH */
    wdot[10] += q; /* O */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[10];
    J[93] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[96] += -2 * dqdci;          /* dwdot[OH]/d[H2O] */
    J[100] += dqdci;              /* dwdot[O]/d[H2O] */
    /* d()/d[OH] */
    dqdci =  + k_f*2*sc[6];
    J[183] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[186] += -2 * dqdci;         /* dwdot[OH]/d[OH] */
    J[190] += dqdci;              /* dwdot[O]/d[OH] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[3];
    J[303] += dqdci;              /* dwdot[H2O]/d[O] */
    J[306] += -2 * dqdci;         /* dwdot[OH]/d[O] */
    J[310] += dqdci;              /* dwdot[O]/d[O] */
    /* d()/dT */
    J[873] += dqdT;               /* dwdot[H2O]/dT */
    J[876] += -2 * dqdT;          /* dwdot[OH]/dT */
    J[880] += dqdT;               /* dwdot[O]/dT */

    /*reaction 33: H2 + OH <=> H2O + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[6];
    k_f = prefactor_units[32] * fwd_A[32]
                * exp(fwd_beta[32] * tc[0] - activation_units[32] * fwd_Ea[32] * invT);
    dlnkfdT = fwd_beta[32] * invT + activation_units[32] * fwd_Ea[32] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[9];
    Kc = exp(-g_RT[3] + g_RT[5] + g_RT[6] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[6]) + (h_RT[3] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* H2O */
    wdot[5] -= q; /* H2 */
    wdot[6] -= q; /* OH */
    wdot[9] += q; /* H */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[9];
    J[93] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[95] -= dqdci;               /* dwdot[H2]/d[H2O] */
    J[96] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[99] += dqdci;               /* dwdot[H]/d[H2O] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[6];
    J[153] += dqdci;              /* dwdot[H2O]/d[H2] */
    J[155] -= dqdci;              /* dwdot[H2]/d[H2] */
    J[156] -= dqdci;              /* dwdot[OH]/d[H2] */
    J[159] += dqdci;              /* dwdot[H]/d[H2] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[5];
    J[183] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[185] -= dqdci;              /* dwdot[H2]/d[OH] */
    J[186] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[189] += dqdci;              /* dwdot[H]/d[OH] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[3];
    J[273] += dqdci;              /* dwdot[H2O]/d[H] */
    J[275] -= dqdci;              /* dwdot[H2]/d[H] */
    J[276] -= dqdci;              /* dwdot[OH]/d[H] */
    J[279] += dqdci;              /* dwdot[H]/d[H] */
    /* d()/dT */
    J[873] += dqdT;               /* dwdot[H2O]/dT */
    J[875] -= dqdT;               /* dwdot[H2]/dT */
    J[876] -= dqdT;               /* dwdot[OH]/dT */
    J[879] += dqdT;               /* dwdot[H]/dT */

    /*reaction 34: 2 HO2 <=> O2 + H2O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[8];
    k_f = prefactor_units[33] * fwd_A[33]
                * exp(fwd_beta[33] * tc[0] - activation_units[33] * fwd_Ea[33] * invT);
    dlnkfdT = fwd_beta[33] * invT + activation_units[33] * fwd_Ea[33] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[7];
    Kc = exp(-g_RT[1] - g_RT[7] + 2*g_RT[8]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[8]) + (h_RT[1] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* O2 */
    wdot[7] += q; /* H2O2 */
    wdot[8] -= 2 * q; /* HO2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[7];
    J[31] += dqdci;               /* dwdot[O2]/d[O2] */
    J[37] += dqdci;               /* dwdot[H2O2]/d[O2] */
    J[38] += -2 * dqdci;          /* dwdot[HO2]/d[O2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[1];
    J[211] += dqdci;              /* dwdot[O2]/d[H2O2] */
    J[217] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
    J[218] += -2 * dqdci;         /* dwdot[HO2]/d[H2O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2*sc[8];
    J[241] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[247] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    J[248] += -2 * dqdci;         /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[871] += dqdT;               /* dwdot[O2]/dT */
    J[877] += dqdT;               /* dwdot[H2O2]/dT */
    J[878] += -2 * dqdT;          /* dwdot[HO2]/dT */

    /*reaction 35: OH + CH2O <=> H2O + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[12];
    k_f = prefactor_units[34] * fwd_A[34]
                * exp(fwd_beta[34] * tc[0] - activation_units[34] * fwd_Ea[34] * invT);
    dlnkfdT = fwd_beta[34] * invT + activation_units[34] * fwd_Ea[34] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[13];
    Kc = exp(-g_RT[3] + g_RT[6] + g_RT[12] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[12]) + (h_RT[3] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* H2O */
    wdot[6] -= q; /* OH */
    wdot[12] -= q; /* CH2O */
    wdot[13] += q; /* HCO */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[13];
    J[93] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[96] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[102] -= dqdci;              /* dwdot[CH2O]/d[H2O] */
    J[103] += dqdci;              /* dwdot[HCO]/d[H2O] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[12];
    J[183] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[186] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[192] -= dqdci;              /* dwdot[CH2O]/d[OH] */
    J[193] += dqdci;              /* dwdot[HCO]/d[OH] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[6];
    J[363] += dqdci;              /* dwdot[H2O]/d[CH2O] */
    J[366] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[372] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[373] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[3];
    J[393] += dqdci;              /* dwdot[H2O]/d[HCO] */
    J[396] -= dqdci;              /* dwdot[OH]/d[HCO] */
    J[402] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    J[403] += dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[873] += dqdT;               /* dwdot[H2O]/dT */
    J[876] -= dqdT;               /* dwdot[OH]/dT */
    J[882] -= dqdT;               /* dwdot[CH2O]/dT */
    J[883] += dqdT;               /* dwdot[HCO]/dT */

    /*reaction 36: HO2 + CH2O <=> H2O2 + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[12];
    k_f = prefactor_units[35] * fwd_A[35]
                * exp(fwd_beta[35] * tc[0] - activation_units[35] * fwd_Ea[35] * invT);
    dlnkfdT = fwd_beta[35] * invT + activation_units[35] * fwd_Ea[35] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[13];
    Kc = exp(-g_RT[7] + g_RT[8] + g_RT[12] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[12]) + (h_RT[7] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] += q; /* H2O2 */
    wdot[8] -= q; /* HO2 */
    wdot[12] -= q; /* CH2O */
    wdot[13] += q; /* HCO */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[13];
    J[217] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
    J[218] -= dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[222] -= dqdci;              /* dwdot[CH2O]/d[H2O2] */
    J[223] += dqdci;              /* dwdot[HCO]/d[H2O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[12];
    J[247] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    J[248] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[252] -= dqdci;              /* dwdot[CH2O]/d[HO2] */
    J[253] += dqdci;              /* dwdot[HCO]/d[HO2] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[8];
    J[367] += dqdci;              /* dwdot[H2O2]/d[CH2O] */
    J[368] -= dqdci;              /* dwdot[HO2]/d[CH2O] */
    J[372] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[373] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[7];
    J[397] += dqdci;              /* dwdot[H2O2]/d[HCO] */
    J[398] -= dqdci;              /* dwdot[HO2]/d[HCO] */
    J[402] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    J[403] += dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[877] += dqdT;               /* dwdot[H2O2]/dT */
    J[878] -= dqdT;               /* dwdot[HO2]/dT */
    J[882] -= dqdT;               /* dwdot[CH2O]/dT */
    J[883] += dqdT;               /* dwdot[HCO]/dT */

    /*reaction 37: O2 + HCO <=> CO + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[13];
    k_f = prefactor_units[36] * fwd_A[36]
                * exp(fwd_beta[36] * tc[0] - activation_units[36] * fwd_Ea[36] * invT);
    dlnkfdT = fwd_beta[36] * invT + activation_units[36] * fwd_Ea[36] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[8];
    Kc = exp(g_RT[1] - g_RT[4] - g_RT[8] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[13]) + (h_RT[4] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* O2 */
    wdot[4] += q; /* CO */
    wdot[8] += q; /* HO2 */
    wdot[13] -= q; /* HCO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[13];
    J[31] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[34] += dqdci;               /* dwdot[CO]/d[O2] */
    J[38] += dqdci;               /* dwdot[HO2]/d[O2] */
    J[43] -= dqdci;               /* dwdot[HCO]/d[O2] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[8];
    J[121] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[124] += dqdci;              /* dwdot[CO]/d[CO] */
    J[128] += dqdci;              /* dwdot[HO2]/d[CO] */
    J[133] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[4];
    J[241] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[244] += dqdci;              /* dwdot[CO]/d[HO2] */
    J[248] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[253] -= dqdci;              /* dwdot[HCO]/d[HO2] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[1];
    J[391] -= dqdci;              /* dwdot[O2]/d[HCO] */
    J[394] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[398] += dqdci;              /* dwdot[HO2]/d[HCO] */
    J[403] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[871] -= dqdT;               /* dwdot[O2]/dT */
    J[874] += dqdT;               /* dwdot[CO]/dT */
    J[878] += dqdT;               /* dwdot[HO2]/dT */
    J[883] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 38: CH3O + CH3 <=> CH2O + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[15];
    k_f = prefactor_units[37] * fwd_A[37]
                * exp(fwd_beta[37] * tc[0] - activation_units[37] * fwd_Ea[37] * invT);
    dlnkfdT = fwd_beta[37] * invT + activation_units[37] * fwd_Ea[37] * invT2;
    /* reverse */
    phi_r = sc[12]*sc[16];
    Kc = exp(g_RT[11] - g_RT[12] + g_RT[15] - g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[11] + h_RT[15]) + (h_RT[12] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[11] -= q; /* CH3O */
    wdot[12] += q; /* CH2O */
    wdot[15] -= q; /* CH3 */
    wdot[16] += q; /* CH4 */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[15];
    J[341] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    J[342] += dqdci;              /* dwdot[CH2O]/d[CH3O] */
    J[345] -= dqdci;              /* dwdot[CH3]/d[CH3O] */
    J[346] += dqdci;              /* dwdot[CH4]/d[CH3O] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[16];
    J[371] -= dqdci;              /* dwdot[CH3O]/d[CH2O] */
    J[372] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[375] -= dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[376] += dqdci;              /* dwdot[CH4]/d[CH2O] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[11];
    J[461] -= dqdci;              /* dwdot[CH3O]/d[CH3] */
    J[462] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    J[465] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[466] += dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[12];
    J[491] -= dqdci;              /* dwdot[CH3O]/d[CH4] */
    J[492] += dqdci;              /* dwdot[CH2O]/d[CH4] */
    J[495] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[496] += dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[881] -= dqdT;               /* dwdot[CH3O]/dT */
    J[882] += dqdT;               /* dwdot[CH2O]/dT */
    J[885] -= dqdT;               /* dwdot[CH3]/dT */
    J[886] += dqdT;               /* dwdot[CH4]/dT */

    /*reaction 39: OH + C2H4 <=> CH2O + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[18];
    k_f = prefactor_units[38] * fwd_A[38]
                * exp(fwd_beta[38] * tc[0] - activation_units[38] * fwd_Ea[38] * invT);
    dlnkfdT = fwd_beta[38] * invT + activation_units[38] * fwd_Ea[38] * invT2;
    /* reverse */
    phi_r = sc[12]*sc[15];
    Kc = exp(g_RT[6] - g_RT[12] - g_RT[15] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[18]) + (h_RT[12] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[6] -= q; /* OH */
    wdot[12] += q; /* CH2O */
    wdot[15] += q; /* CH3 */
    wdot[18] -= q; /* C2H4 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[18];
    J[186] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[192] += dqdci;              /* dwdot[CH2O]/d[OH] */
    J[195] += dqdci;              /* dwdot[CH3]/d[OH] */
    J[198] -= dqdci;              /* dwdot[C2H4]/d[OH] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[15];
    J[366] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[372] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[375] += dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[378] -= dqdci;              /* dwdot[C2H4]/d[CH2O] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[12];
    J[456] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[462] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    J[465] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[468] -= dqdci;              /* dwdot[C2H4]/d[CH3] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[6];
    J[546] -= dqdci;              /* dwdot[OH]/d[C2H4] */
    J[552] += dqdci;              /* dwdot[CH2O]/d[C2H4] */
    J[555] += dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[558] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    /* d()/dT */
    J[876] -= dqdT;               /* dwdot[OH]/dT */
    J[882] += dqdT;               /* dwdot[CH2O]/dT */
    J[885] += dqdT;               /* dwdot[CH3]/dT */
    J[888] -= dqdT;               /* dwdot[C2H4]/dT */

    /*reaction 40: OH + C2H4 <=> H2O + C2H3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[18];
    k_f = prefactor_units[39] * fwd_A[39]
                * exp(fwd_beta[39] * tc[0] - activation_units[39] * fwd_Ea[39] * invT);
    dlnkfdT = fwd_beta[39] * invT + activation_units[39] * fwd_Ea[39] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[17];
    Kc = exp(-g_RT[3] + g_RT[6] - g_RT[17] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[18]) + (h_RT[3] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* H2O */
    wdot[6] -= q; /* OH */
    wdot[17] += q; /* C2H3 */
    wdot[18] -= q; /* C2H4 */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[17];
    J[93] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[96] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[107] += dqdci;              /* dwdot[C2H3]/d[H2O] */
    J[108] -= dqdci;              /* dwdot[C2H4]/d[H2O] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[18];
    J[183] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[186] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[197] += dqdci;              /* dwdot[C2H3]/d[OH] */
    J[198] -= dqdci;              /* dwdot[C2H4]/d[OH] */
    /* d()/d[C2H3] */
    dqdci =  - k_r*sc[3];
    J[513] += dqdci;              /* dwdot[H2O]/d[C2H3] */
    J[516] -= dqdci;              /* dwdot[OH]/d[C2H3] */
    J[527] += dqdci;              /* dwdot[C2H3]/d[C2H3] */
    J[528] -= dqdci;              /* dwdot[C2H4]/d[C2H3] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[6];
    J[543] += dqdci;              /* dwdot[H2O]/d[C2H4] */
    J[546] -= dqdci;              /* dwdot[OH]/d[C2H4] */
    J[557] += dqdci;              /* dwdot[C2H3]/d[C2H4] */
    J[558] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    /* d()/dT */
    J[873] += dqdT;               /* dwdot[H2O]/dT */
    J[876] -= dqdT;               /* dwdot[OH]/dT */
    J[887] += dqdT;               /* dwdot[C2H3]/dT */
    J[888] -= dqdT;               /* dwdot[C2H4]/dT */

    /*reaction 41: O2 + C2H3 <=> CH2O + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[17];
    k_f = prefactor_units[40] * fwd_A[40]
                * exp(fwd_beta[40] * tc[0] - activation_units[40] * fwd_Ea[40] * invT);
    dlnkfdT = fwd_beta[40] * invT + activation_units[40] * fwd_Ea[40] * invT2;
    /* reverse */
    phi_r = sc[12]*sc[13];
    Kc = exp(g_RT[1] - g_RT[12] - g_RT[13] + g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[17]) + (h_RT[12] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* O2 */
    wdot[12] += q; /* CH2O */
    wdot[13] += q; /* HCO */
    wdot[17] -= q; /* C2H3 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[17];
    J[31] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[42] += dqdci;               /* dwdot[CH2O]/d[O2] */
    J[43] += dqdci;               /* dwdot[HCO]/d[O2] */
    J[47] -= dqdci;               /* dwdot[C2H3]/d[O2] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[13];
    J[361] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    J[372] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[373] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[377] -= dqdci;              /* dwdot[C2H3]/d[CH2O] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[12];
    J[391] -= dqdci;              /* dwdot[O2]/d[HCO] */
    J[402] += dqdci;              /* dwdot[CH2O]/d[HCO] */
    J[403] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[407] -= dqdci;              /* dwdot[C2H3]/d[HCO] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[1];
    J[511] -= dqdci;              /* dwdot[O2]/d[C2H3] */
    J[522] += dqdci;              /* dwdot[CH2O]/d[C2H3] */
    J[523] += dqdci;              /* dwdot[HCO]/d[C2H3] */
    J[527] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    /* d()/dT */
    J[871] -= dqdT;               /* dwdot[O2]/dT */
    J[882] += dqdT;               /* dwdot[CH2O]/dT */
    J[883] += dqdT;               /* dwdot[HCO]/dT */
    J[887] -= dqdT;               /* dwdot[C2H3]/dT */

    /*reaction 42: HCO + C2H3 <=> CO + C2H4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[13]*sc[17];
    k_f = prefactor_units[41] * fwd_A[41]
                * exp(fwd_beta[41] * tc[0] - activation_units[41] * fwd_Ea[41] * invT);
    dlnkfdT = fwd_beta[41] * invT + activation_units[41] * fwd_Ea[41] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[18];
    Kc = exp(-g_RT[4] + g_RT[13] + g_RT[17] - g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[13] + h_RT[17]) + (h_RT[4] + h_RT[18]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* CO */
    wdot[13] -= q; /* HCO */
    wdot[17] -= q; /* C2H3 */
    wdot[18] += q; /* C2H4 */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[18];
    J[124] += dqdci;              /* dwdot[CO]/d[CO] */
    J[133] -= dqdci;              /* dwdot[HCO]/d[CO] */
    J[137] -= dqdci;              /* dwdot[C2H3]/d[CO] */
    J[138] += dqdci;              /* dwdot[C2H4]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[17];
    J[394] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[403] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    J[407] -= dqdci;              /* dwdot[C2H3]/d[HCO] */
    J[408] += dqdci;              /* dwdot[C2H4]/d[HCO] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[13];
    J[514] += dqdci;              /* dwdot[CO]/d[C2H3] */
    J[523] -= dqdci;              /* dwdot[HCO]/d[C2H3] */
    J[527] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    J[528] += dqdci;              /* dwdot[C2H4]/d[C2H3] */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[4];
    J[544] += dqdci;              /* dwdot[CO]/d[C2H4] */
    J[553] -= dqdci;              /* dwdot[HCO]/d[C2H4] */
    J[557] -= dqdci;              /* dwdot[C2H3]/d[C2H4] */
    J[558] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    /* d()/dT */
    J[874] += dqdT;               /* dwdot[CO]/dT */
    J[883] -= dqdT;               /* dwdot[HCO]/dT */
    J[887] -= dqdT;               /* dwdot[C2H3]/dT */
    J[888] += dqdT;               /* dwdot[C2H4]/dT */

    /*reaction 43: O2 + C2H5 <=> HO2 + C2H4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[19];
    k_f = prefactor_units[42] * fwd_A[42]
                * exp(fwd_beta[42] * tc[0] - activation_units[42] * fwd_Ea[42] * invT);
    dlnkfdT = fwd_beta[42] * invT + activation_units[42] * fwd_Ea[42] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[18];
    Kc = exp(g_RT[1] - g_RT[8] - g_RT[18] + g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[19]) + (h_RT[8] + h_RT[18]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    wdot[18] += q; /* C2H4 */
    wdot[19] -= q; /* C2H5 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[19];
    J[31] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[38] += dqdci;               /* dwdot[HO2]/d[O2] */
    J[48] += dqdci;               /* dwdot[C2H4]/d[O2] */
    J[49] -= dqdci;               /* dwdot[C2H5]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[18];
    J[241] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[248] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[258] += dqdci;              /* dwdot[C2H4]/d[HO2] */
    J[259] -= dqdci;              /* dwdot[C2H5]/d[HO2] */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[8];
    J[541] -= dqdci;              /* dwdot[O2]/d[C2H4] */
    J[548] += dqdci;              /* dwdot[HO2]/d[C2H4] */
    J[558] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[559] -= dqdci;              /* dwdot[C2H5]/d[C2H4] */
    /* d()/d[C2H5] */
    dqdci =  + k_f*sc[1];
    J[571] -= dqdci;              /* dwdot[O2]/d[C2H5] */
    J[578] += dqdci;              /* dwdot[HO2]/d[C2H5] */
    J[588] += dqdci;              /* dwdot[C2H4]/d[C2H5] */
    J[589] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
    /* d()/dT */
    J[871] -= dqdT;               /* dwdot[O2]/dT */
    J[878] += dqdT;               /* dwdot[HO2]/dT */
    J[888] += dqdT;               /* dwdot[C2H4]/dT */
    J[889] -= dqdT;               /* dwdot[C2H5]/dT */

    /*reaction 44: O2 + CH4 <=> HO2 + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[16];
    k_f = prefactor_units[43] * fwd_A[43]
                * exp(fwd_beta[43] * tc[0] - activation_units[43] * fwd_Ea[43] * invT);
    dlnkfdT = fwd_beta[43] * invT + activation_units[43] * fwd_Ea[43] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[15];
    Kc = exp(g_RT[1] - g_RT[8] - g_RT[15] + g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[16]) + (h_RT[8] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    wdot[15] += q; /* CH3 */
    wdot[16] -= q; /* CH4 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[16];
    J[31] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[38] += dqdci;               /* dwdot[HO2]/d[O2] */
    J[45] += dqdci;               /* dwdot[CH3]/d[O2] */
    J[46] -= dqdci;               /* dwdot[CH4]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[15];
    J[241] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[248] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[255] += dqdci;              /* dwdot[CH3]/d[HO2] */
    J[256] -= dqdci;              /* dwdot[CH4]/d[HO2] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[8];
    J[451] -= dqdci;              /* dwdot[O2]/d[CH3] */
    J[458] += dqdci;              /* dwdot[HO2]/d[CH3] */
    J[465] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[466] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[1];
    J[481] -= dqdci;              /* dwdot[O2]/d[CH4] */
    J[488] += dqdci;              /* dwdot[HO2]/d[CH4] */
    J[495] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[496] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[871] -= dqdT;               /* dwdot[O2]/dT */
    J[878] += dqdT;               /* dwdot[HO2]/dT */
    J[885] += dqdT;               /* dwdot[CH3]/dT */
    J[886] -= dqdT;               /* dwdot[CH4]/dT */

    /*reaction 45: OH + HO2 <=> O2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[8];
    k_f = prefactor_units[44] * fwd_A[44]
                * exp(fwd_beta[44] * tc[0] - activation_units[44] * fwd_Ea[44] * invT);
    dlnkfdT = fwd_beta[44] * invT + activation_units[44] * fwd_Ea[44] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[3];
    Kc = exp(-g_RT[1] - g_RT[3] + g_RT[6] + g_RT[8]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[8]) + (h_RT[1] + h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* O2 */
    wdot[3] += q; /* H2O */
    wdot[6] -= q; /* OH */
    wdot[8] -= q; /* HO2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[3];
    J[31] += dqdci;               /* dwdot[O2]/d[O2] */
    J[33] += dqdci;               /* dwdot[H2O]/d[O2] */
    J[36] -= dqdci;               /* dwdot[OH]/d[O2] */
    J[38] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[1];
    J[91] += dqdci;               /* dwdot[O2]/d[H2O] */
    J[93] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[96] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[98] -= dqdci;               /* dwdot[HO2]/d[H2O] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[8];
    J[181] += dqdci;              /* dwdot[O2]/d[OH] */
    J[183] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[186] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[188] -= dqdci;              /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[6];
    J[241] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[243] += dqdci;              /* dwdot[H2O]/d[HO2] */
    J[246] -= dqdci;              /* dwdot[OH]/d[HO2] */
    J[248] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[871] += dqdT;               /* dwdot[O2]/dT */
    J[873] += dqdT;               /* dwdot[H2O]/dT */
    J[876] -= dqdT;               /* dwdot[OH]/dT */
    J[878] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 46: O2 + CH3 <=> OH + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[15];
    k_f = prefactor_units[45] * fwd_A[45]
                * exp(fwd_beta[45] * tc[0] - activation_units[45] * fwd_Ea[45] * invT);
    dlnkfdT = fwd_beta[45] * invT + activation_units[45] * fwd_Ea[45] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[12];
    Kc = exp(g_RT[1] - g_RT[6] - g_RT[12] + g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[15]) + (h_RT[6] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* O2 */
    wdot[6] += q; /* OH */
    wdot[12] += q; /* CH2O */
    wdot[15] -= q; /* CH3 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[15];
    J[31] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[36] += dqdci;               /* dwdot[OH]/d[O2] */
    J[42] += dqdci;               /* dwdot[CH2O]/d[O2] */
    J[45] -= dqdci;               /* dwdot[CH3]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[12];
    J[181] -= dqdci;              /* dwdot[O2]/d[OH] */
    J[186] += dqdci;              /* dwdot[OH]/d[OH] */
    J[192] += dqdci;              /* dwdot[CH2O]/d[OH] */
    J[195] -= dqdci;              /* dwdot[CH3]/d[OH] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[6];
    J[361] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    J[366] += dqdci;              /* dwdot[OH]/d[CH2O] */
    J[372] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[375] -= dqdci;              /* dwdot[CH3]/d[CH2O] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[1];
    J[451] -= dqdci;              /* dwdot[O2]/d[CH3] */
    J[456] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[462] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    J[465] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[871] -= dqdT;               /* dwdot[O2]/dT */
    J[876] += dqdT;               /* dwdot[OH]/dT */
    J[882] += dqdT;               /* dwdot[CH2O]/dT */
    J[885] -= dqdT;               /* dwdot[CH3]/dT */

    /*reaction 47: H + CH4 <=> H2 + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[16];
    k_f = prefactor_units[46] * fwd_A[46]
                * exp(fwd_beta[46] * tc[0] - activation_units[46] * fwd_Ea[46] * invT);
    dlnkfdT = fwd_beta[46] * invT + activation_units[46] * fwd_Ea[46] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[15];
    Kc = exp(-g_RT[5] + g_RT[9] - g_RT[15] + g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[9] + h_RT[16]) + (h_RT[5] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] += q; /* H2 */
    wdot[9] -= q; /* H */
    wdot[15] += q; /* CH3 */
    wdot[16] -= q; /* CH4 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[15];
    J[155] += dqdci;              /* dwdot[H2]/d[H2] */
    J[159] -= dqdci;              /* dwdot[H]/d[H2] */
    J[165] += dqdci;              /* dwdot[CH3]/d[H2] */
    J[166] -= dqdci;              /* dwdot[CH4]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[16];
    J[275] += dqdci;              /* dwdot[H2]/d[H] */
    J[279] -= dqdci;              /* dwdot[H]/d[H] */
    J[285] += dqdci;              /* dwdot[CH3]/d[H] */
    J[286] -= dqdci;              /* dwdot[CH4]/d[H] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[5];
    J[455] += dqdci;              /* dwdot[H2]/d[CH3] */
    J[459] -= dqdci;              /* dwdot[H]/d[CH3] */
    J[465] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[466] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[9];
    J[485] += dqdci;              /* dwdot[H2]/d[CH4] */
    J[489] -= dqdci;              /* dwdot[H]/d[CH4] */
    J[495] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[496] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[875] += dqdT;               /* dwdot[H2]/dT */
    J[879] -= dqdT;               /* dwdot[H]/dT */
    J[885] += dqdT;               /* dwdot[CH3]/dT */
    J[886] -= dqdT;               /* dwdot[CH4]/dT */

    /*reaction 48: OH + CH4 <=> H2O + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[16];
    k_f = prefactor_units[47] * fwd_A[47]
                * exp(fwd_beta[47] * tc[0] - activation_units[47] * fwd_Ea[47] * invT);
    dlnkfdT = fwd_beta[47] * invT + activation_units[47] * fwd_Ea[47] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[15];
    Kc = exp(-g_RT[3] + g_RT[6] - g_RT[15] + g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[16]) + (h_RT[3] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* H2O */
    wdot[6] -= q; /* OH */
    wdot[15] += q; /* CH3 */
    wdot[16] -= q; /* CH4 */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[15];
    J[93] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[96] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[105] += dqdci;              /* dwdot[CH3]/d[H2O] */
    J[106] -= dqdci;              /* dwdot[CH4]/d[H2O] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[16];
    J[183] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[186] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[195] += dqdci;              /* dwdot[CH3]/d[OH] */
    J[196] -= dqdci;              /* dwdot[CH4]/d[OH] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[3];
    J[453] += dqdci;              /* dwdot[H2O]/d[CH3] */
    J[456] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[465] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[466] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[6];
    J[483] += dqdci;              /* dwdot[H2O]/d[CH4] */
    J[486] -= dqdci;              /* dwdot[OH]/d[CH4] */
    J[495] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[496] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[873] += dqdT;               /* dwdot[H2O]/dT */
    J[876] -= dqdT;               /* dwdot[OH]/dT */
    J[885] += dqdT;               /* dwdot[CH3]/dT */
    J[886] -= dqdT;               /* dwdot[CH4]/dT */

    /*reaction 49: O + CH4 <=> OH + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[10]*sc[16];
    k_f = prefactor_units[48] * fwd_A[48]
                * exp(fwd_beta[48] * tc[0] - activation_units[48] * fwd_Ea[48] * invT);
    dlnkfdT = fwd_beta[48] * invT + activation_units[48] * fwd_Ea[48] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[15];
    Kc = exp(-g_RT[6] + g_RT[10] - g_RT[15] + g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[10] + h_RT[16]) + (h_RT[6] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[6] += q; /* OH */
    wdot[10] -= q; /* O */
    wdot[15] += q; /* CH3 */
    wdot[16] -= q; /* CH4 */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[15];
    J[186] += dqdci;              /* dwdot[OH]/d[OH] */
    J[190] -= dqdci;              /* dwdot[O]/d[OH] */
    J[195] += dqdci;              /* dwdot[CH3]/d[OH] */
    J[196] -= dqdci;              /* dwdot[CH4]/d[OH] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[16];
    J[306] += dqdci;              /* dwdot[OH]/d[O] */
    J[310] -= dqdci;              /* dwdot[O]/d[O] */
    J[315] += dqdci;              /* dwdot[CH3]/d[O] */
    J[316] -= dqdci;              /* dwdot[CH4]/d[O] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[6];
    J[456] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[460] -= dqdci;              /* dwdot[O]/d[CH3] */
    J[465] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[466] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[10];
    J[486] += dqdci;              /* dwdot[OH]/d[CH4] */
    J[490] -= dqdci;              /* dwdot[O]/d[CH4] */
    J[495] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[496] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[876] += dqdT;               /* dwdot[OH]/dT */
    J[880] -= dqdT;               /* dwdot[O]/dT */
    J[885] += dqdT;               /* dwdot[CH3]/dT */
    J[886] -= dqdT;               /* dwdot[CH4]/dT */

    /*reaction 50: HO2 + CH4 <=> H2O2 + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[16];
    k_f = prefactor_units[49] * fwd_A[49]
                * exp(fwd_beta[49] * tc[0] - activation_units[49] * fwd_Ea[49] * invT);
    dlnkfdT = fwd_beta[49] * invT + activation_units[49] * fwd_Ea[49] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[15];
    Kc = exp(-g_RT[7] + g_RT[8] - g_RT[15] + g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[16]) + (h_RT[7] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] += q; /* H2O2 */
    wdot[8] -= q; /* HO2 */
    wdot[15] += q; /* CH3 */
    wdot[16] -= q; /* CH4 */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[15];
    J[217] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
    J[218] -= dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[225] += dqdci;              /* dwdot[CH3]/d[H2O2] */
    J[226] -= dqdci;              /* dwdot[CH4]/d[H2O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[16];
    J[247] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    J[248] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[255] += dqdci;              /* dwdot[CH3]/d[HO2] */
    J[256] -= dqdci;              /* dwdot[CH4]/d[HO2] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[7];
    J[457] += dqdci;              /* dwdot[H2O2]/d[CH3] */
    J[458] -= dqdci;              /* dwdot[HO2]/d[CH3] */
    J[465] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[466] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[8];
    J[487] += dqdci;              /* dwdot[H2O2]/d[CH4] */
    J[488] -= dqdci;              /* dwdot[HO2]/d[CH4] */
    J[495] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[496] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[877] += dqdT;               /* dwdot[H2O2]/dT */
    J[878] -= dqdT;               /* dwdot[HO2]/dT */
    J[885] += dqdT;               /* dwdot[CH3]/dT */
    J[886] -= dqdT;               /* dwdot[CH4]/dT */

    /*reaction 51: CH2 + CH4 <=> 2 CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[14]*sc[16];
    k_f = prefactor_units[50] * fwd_A[50]
                * exp(fwd_beta[50] * tc[0] - activation_units[50] * fwd_Ea[50] * invT);
    dlnkfdT = fwd_beta[50] * invT + activation_units[50] * fwd_Ea[50] * invT2;
    /* reverse */
    phi_r = sc[15]*sc[15];
    Kc = exp(g_RT[14] - 2*g_RT[15] + g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[14] + h_RT[16]) + (2*h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[14] -= q; /* CH2 */
    wdot[15] += 2 * q; /* CH3 */
    wdot[16] -= q; /* CH4 */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[16];
    J[434] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[435] += 2 * dqdci;          /* dwdot[CH3]/d[CH2] */
    J[436] -= dqdci;              /* dwdot[CH4]/d[CH2] */
    /* d()/d[CH3] */
    dqdci =  - k_r*2*sc[15];
    J[464] -= dqdci;              /* dwdot[CH2]/d[CH3] */
    J[465] += 2 * dqdci;          /* dwdot[CH3]/d[CH3] */
    J[466] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[14];
    J[494] -= dqdci;              /* dwdot[CH2]/d[CH4] */
    J[495] += 2 * dqdci;          /* dwdot[CH3]/d[CH4] */
    J[496] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[884] -= dqdT;               /* dwdot[CH2]/dT */
    J[885] += 2 * dqdT;           /* dwdot[CH3]/dT */
    J[886] -= dqdT;               /* dwdot[CH4]/dT */

    /*reaction 52: C3H6 <=> CH3 + C2H3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[22];
    k_f = prefactor_units[51] * fwd_A[51]
                * exp(fwd_beta[51] * tc[0] - activation_units[51] * fwd_Ea[51] * invT);
    dlnkfdT = fwd_beta[51] * invT + activation_units[51] * fwd_Ea[51] * invT2;
    /* reverse */
    phi_r = sc[15]*sc[17];
    Kc = refC * exp(-g_RT[15] - g_RT[17] + g_RT[22]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[22]) + (h_RT[15] + h_RT[17]) - 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[15] += q; /* CH3 */
    wdot[17] += q; /* C2H3 */
    wdot[22] -= q; /* C3H6 */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[17];
    J[465] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[467] += dqdci;              /* dwdot[C2H3]/d[CH3] */
    J[472] -= dqdci;              /* dwdot[C3H6]/d[CH3] */
    /* d()/d[C2H3] */
    dqdci =  - k_r*sc[15];
    J[525] += dqdci;              /* dwdot[CH3]/d[C2H3] */
    J[527] += dqdci;              /* dwdot[C2H3]/d[C2H3] */
    J[532] -= dqdci;              /* dwdot[C3H6]/d[C2H3] */
    /* d()/d[C3H6] */
    dqdci =  + k_f;
    J[675] += dqdci;              /* dwdot[CH3]/d[C3H6] */
    J[677] += dqdci;              /* dwdot[C2H3]/d[C3H6] */
    J[682] -= dqdci;              /* dwdot[C3H6]/d[C3H6] */
    /* d()/dT */
    J[885] += dqdT;               /* dwdot[CH3]/dT */
    J[887] += dqdT;               /* dwdot[C2H3]/dT */
    J[892] -= dqdT;               /* dwdot[C3H6]/dT */

    double c_R[29], dcRdT[29], e_RT[29];
    double * eh_RT;
    if (consP) {
        cp_R(c_R, tc);
        dcvpRdT(dcRdT, tc);
        eh_RT = &h_RT[0];
    }
    else {
        cv_R(c_R, tc);
        dcvpRdT(dcRdT, tc);
        speciesInternalEnergy(e_RT, tc);
        eh_RT = &e_RT[0];
    }

    double cmix = 0.0, ehmix = 0.0, dcmixdT=0.0, dehmixdT=0.0;
    for (int k = 0; k < 29; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[870+k];
    }

    double cmixinv = 1.0/cmix;
    double tmp1 = ehmix*cmixinv;
    double tmp3 = cmixinv*T;
    double tmp2 = tmp1*tmp3;
    double dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 29; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 29; ++m) {
            dehmixdc += eh_RT[m]*J[k*30+m];
        }
        J[k*30+29] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[899] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;
}


/*compute d(Cp/R)/dT and d(Cv/R)/dT at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void dcvpRdT(double * restrict species, double * restrict tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1382 kelvin */
    if (T < 1382) {
        /*species 24: C7H15-2 */
        species[24] =
            +7.56726570e-02
            -8.14947268e-05 * tc[1]
            +2.79803683e-08 * tc[2]
            -1.96944298e-12 * tc[3];
    } else {
        /*species 24: C7H15-2 */
        species[24] =
            +3.23324804e-02
            -2.18547614e-05 * tc[1]
            +5.05071180e-09 * tc[2]
            -3.88709636e-13 * tc[3];
    }

    /*species with midpoint at T=1383 kelvin */
    if (T < 1383) {
        /*species 27: C5H11CO */
        species[27] =
            +6.17863563e-02
            -7.48269380e-05 * tc[1]
            +3.39851385e-08 * tc[2]
            -5.47670792e-12 * tc[3];
    } else {
        /*species 27: C5H11CO */
        species[27] =
            +2.50466029e-02
            -1.70972269e-05 * tc[1]
            +3.97673832e-09 * tc[2]
            -3.07401318e-13 * tc[3];
    }

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 1: O2 */
        species[1] =
            +1.12748600e-03
            -1.15123000e-06 * tc[1]
            +3.94163100e-09 * tc[2]
            -3.50742160e-12 * tc[3];
        /*species 2: CO2 */
        species[2] =
            +9.92207200e-03
            -2.08182200e-05 * tc[1]
            +2.06000610e-08 * tc[2]
            -8.46912000e-12 * tc[3];
        /*species 3: H2O */
        species[3] =
            +3.47498200e-03
            -1.27093920e-05 * tc[1]
            +2.09057430e-08 * tc[2]
            -1.00263520e-11 * tc[3];
        /*species 4: CO */
        species[4] =
            +1.51194100e-03
            -7.76351000e-06 * tc[1]
            +1.67458320e-08 * tc[2]
            -9.89980400e-12 * tc[3];
        /*species 5: H2 */
        species[5] =
            +8.24944200e-04
            -1.62860300e-06 * tc[1]
            -2.84263020e-10 * tc[2]
            +1.65394880e-12 * tc[3];
        /*species 6: OH */
        species[6] =
            +1.85091000e-04
            -3.35233000e-06 * tc[1]
            +7.16160900e-09 * tc[2]
            -3.37257680e-12 * tc[3];
        /*species 7: H2O2 */
        species[7] =
            +6.56922600e-03
            -2.97002600e-07 * tc[1]
            -1.38774180e-08 * tc[2]
            +9.88606000e-12 * tc[3];
        /*species 8: HO2 */
        species[8] =
            +4.99669700e-03
            -7.58199400e-06 * tc[1]
            +7.06257600e-09 * tc[2]
            -3.23560960e-12 * tc[3];
        /*species 9: H */
        species[9] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 10: O */
        species[10] =
            -1.63816600e-03
            +4.84206400e-06 * tc[1]
            -4.80852900e-09 * tc[2]
            +1.55627840e-12 * tc[3];
        /*species 11: CH3O */
        species[11] =
            +7.21659500e-03
            +1.06769440e-05 * tc[1]
            -2.21329080e-08 * tc[2]
            +8.30244400e-12 * tc[3];
        /*species 12: CH2O */
        species[12] =
            +1.26314400e-02
            -3.77633600e-05 * tc[1]
            +6.15009300e-08 * tc[2]
            -3.36529480e-11 * tc[3];
        /*species 13: HCO */
        species[13] =
            +6.19914700e-03
            -1.92461680e-05 * tc[1]
            +3.26947500e-08 * tc[2]
            -1.82995400e-11 * tc[3];
        /*species 14: CH2 */
        species[14] =
            +1.15981900e-03
            +4.97917000e-07 * tc[1]
            +2.64025080e-09 * tc[2]
            -2.93297400e-12 * tc[3];
        /*species 15: CH3 */
        species[15] =
            +1.11241000e-02
            -3.36044000e-05 * tc[1]
            +4.86548700e-08 * tc[2]
            -2.34598120e-11 * tc[3];
        /*species 16: CH4 */
        species[16] =
            +1.74766800e-02
            -5.56681800e-05 * tc[1]
            +9.14912400e-08 * tc[2]
            -4.89572400e-11 * tc[3];
        /*species 17: C2H3 */
        species[17] =
            +7.37147600e-03
            +4.21974600e-06 * tc[1]
            -3.96492600e-09 * tc[2]
            -4.73913600e-12 * tc[3];
        /*species 18: C2H4 */
        species[18] =
            +2.79616300e-02
            -6.77735400e-05 * tc[1]
            +8.35545600e-08 * tc[2]
            -3.89515160e-11 * tc[3];
        /*species 19: C2H5 */
        species[19] =
            +8.71913300e-03
            +8.83967800e-06 * tc[1]
            +2.80161090e-09 * tc[2]
            -1.57110920e-11 * tc[3];
        /*species 20: C3H4 */
        species[20] =
            +1.21223371e-02
            +3.70810800e-05 * tc[1]
            -1.03577542e-07 * tc[2]
            +6.13413556e-11 * tc[3];
        /*species 21: C3H5 */
        species[21] =
            +9.48414335e-03
            +4.84686736e-05 * tc[1]
            -1.09681203e-07 * tc[2]
            +5.94369424e-11 * tc[3];
        /*species 23: C3H7 */
        species[23] =
            +2.59919800e-02
            +4.76010800e-06 * tc[1]
            -5.88287070e-08 * tc[2]
            +3.74929880e-11 * tc[3];
        /*species 28: N2 */
        species[28] =
            +1.40824000e-03
            -7.92644400e-06 * tc[1]
            +1.69245450e-08 * tc[2]
            -9.77942000e-12 * tc[3];
    } else {
        /*species 1: O2 */
        species[1] =
            +6.13519700e-04
            -2.51768400e-07 * tc[1]
            +5.32584300e-11 * tc[2]
            -4.54574000e-15 * tc[3];
        /*species 2: CO2 */
        species[2] =
            +3.14016900e-03
            -2.55682200e-06 * tc[1]
            +7.18199100e-10 * tc[2]
            -6.67613200e-14 * tc[3];
        /*species 3: H2O */
        species[3] =
            +3.05629300e-03
            -1.74605200e-06 * tc[1]
            +3.60298800e-10 * tc[2]
            -2.55664720e-14 * tc[3];
        /*species 4: CO */
        species[4] =
            +1.44268900e-03
            -1.12616560e-06 * tc[1]
            +3.05574300e-10 * tc[2]
            -2.76438080e-14 * tc[3];
        /*species 5: H2 */
        species[5] =
            +7.00064400e-04
            -1.12676580e-07 * tc[1]
            -2.76947340e-11 * tc[2]
            +6.33100800e-15 * tc[3];
        /*species 6: OH */
        species[6] =
            +1.01397400e-03
            -4.55375400e-07 * tc[1]
            +6.52405200e-11 * tc[2]
            -2.05052200e-15 * tc[3];
        /*species 7: H2O2 */
        species[7] =
            +4.33613600e-03
            -2.94937800e-06 * tc[1]
            +7.04671200e-10 * tc[2]
            -5.72661600e-14 * tc[3];
        /*species 8: HO2 */
        species[8] =
            +2.13129600e-03
            -1.06162900e-06 * tc[1]
            +1.83368070e-10 * tc[2]
            -1.13646600e-14 * tc[3];
        /*species 9: H */
        species[9] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 10: O */
        species[10] =
            -2.75506200e-05
            -6.20560600e-09 * tc[1]
            +1.36532010e-11 * tc[2]
            -1.74722080e-15 * tc[3];
        /*species 11: CH3O */
        species[11] =
            +7.87149700e-03
            -5.31276800e-06 * tc[1]
            +1.18332930e-09 * tc[2]
            -8.45046400e-14 * tc[3];
        /*species 12: CH2O */
        species[12] =
            +6.68132100e-03
            -5.25791000e-06 * tc[1]
            +1.42114590e-09 * tc[2]
            -1.28500680e-13 * tc[3];
        /*species 13: HCO */
        species[13] =
            +3.34557300e-03
            -2.67001200e-06 * tc[1]
            +7.41171900e-10 * tc[2]
            -6.85540400e-14 * tc[3];
        /*species 14: CH2 */
        species[14] =
            +1.93305700e-03
            -3.37403200e-07 * tc[1]
            -3.02969700e-10 * tc[2]
            +7.23302400e-14 * tc[3];
        /*species 15: CH3 */
        species[15] =
            +6.13797400e-03
            -4.46069000e-06 * tc[1]
            +1.13554830e-09 * tc[2]
            -9.80863600e-14 * tc[3];
        /*species 16: CH4 */
        species[16] =
            +1.02372400e-02
            -7.75025800e-06 * tc[1]
            +2.03567550e-09 * tc[2]
            -1.80136920e-13 * tc[3];
        /*species 17: C2H3 */
        species[17] =
            +4.01774600e-03
            -7.93348000e-07 * tc[1]
            -4.32380100e-10 * tc[2]
            +9.51457600e-14 * tc[3];
        /*species 18: C2H4 */
        species[18] =
            +1.14851800e-02
            -8.83677000e-06 * tc[1]
            +2.35338030e-09 * tc[2]
            -2.10673920e-13 * tc[3];
        /*species 19: C2H5 */
        species[19] =
            +6.48407700e-03
            -1.28561300e-06 * tc[1]
            -7.04363700e-10 * tc[2]
            +1.55235080e-13 * tc[3];
        /*species 20: C3H4 */
        species[20] =
            +1.11336262e-02
            -7.92578036e-06 * tc[1]
            +1.90690132e-09 * tc[2]
            -1.51499954e-13 * tc[3];
        /*species 21: C3H5 */
        species[21] =
            +1.33152246e-02
            -9.56666200e-06 * tc[1]
            +2.31584944e-09 * tc[2]
            -1.84772323e-13 * tc[3];
        /*species 23: C3H7 */
        species[23] =
            +1.60442030e-02
            -1.05666440e-05 * tc[1]
            +2.28895770e-09 * tc[2]
            -1.57569136e-13 * tc[3];
        /*species 28: N2 */
        species[28] =
            +1.48797700e-03
            -1.13695220e-06 * tc[1]
            +3.02911200e-10 * tc[2]
            -2.70134040e-14 * tc[3];
    }

    /*species with midpoint at T=1388 kelvin */
    if (T < 1388) {
        /*species 22: C3H6 */
        species[22] =
            +2.89107662e-02
            -3.09773616e-05 * tc[1]
            +1.16644263e-08 * tc[2]
            -1.35156141e-12 * tc[3];
    } else {
        /*species 22: C3H6 */
        species[22] =
            +1.37023634e-02
            -9.32499466e-06 * tc[1]
            +2.16376321e-09 * tc[2]
            -1.66948050e-13 * tc[3];
    }

    /*species with midpoint at T=1390 kelvin */
    if (T < 1390) {
        /*species 25: C7H15O2 */
        species[25] =
            +8.34651906e-02
            -1.02779464e-04 * tc[1]
            +4.92652986e-08 * tc[2]
            -8.78020864e-12 * tc[3];
    } else {
        /*species 25: C7H15O2 */
        species[25] =
            +3.50716920e-02
            -2.40880612e-05 * tc[1]
            +5.62394466e-09 * tc[2]
            -4.35791164e-13 * tc[3];
    }

    /*species with midpoint at T=1391 kelvin */
    if (T < 1391) {
        /*species 0: NC7H16 */
        species[0] =
            +8.54355820e-02
            -1.05069357e-04 * tc[1]
            +4.88837163e-08 * tc[2]
            -8.09579700e-12 * tc[3];
    } else {
        /*species 0: NC7H16 */
        species[0] =
            +3.47675750e-02
            -2.36814258e-05 * tc[1]
            +5.49895434e-09 * tc[2]
            -4.24521064e-13 * tc[3];
    }

    /*species with midpoint at T=1396 kelvin */
    if (T < 1396) {
        /*species 26: C7KET12 */
        species[26] =
            +1.01207869e-01
            -1.53171199e-04 * tc[1]
            +9.02215818e-08 * tc[2]
            -1.93161117e-11 * tc[3];
    } else {
        /*species 26: C7KET12 */
        species[26] =
            +3.06622294e-02
            -2.11127180e-05 * tc[1]
            +4.93882029e-09 * tc[2]
            -3.83268670e-13 * tc[3];
    }
    return;
}


/*compute the progress rate for each reaction */
void progressRate(double * restrict qdot, double * restrict sc, double T)
{
    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];

    if (T != T_save)
    {
        T_save = T;
        comp_k_f(tc,invT,k_f_save);
        comp_Kc(tc,invT,Kc_save);
    }

    double q_f[52], q_r[52];
    comp_qfqr(q_f, q_r, sc, tc, invT);

    for (int i = 0; i < 52; ++i) {
        qdot[i] = q_f[i] - q_r[i];
    }

    return;
}


/*compute the progress rate for each reaction */
void progressRateFR(double * restrict q_f, double * restrict q_r, double * restrict sc, double T)
{
    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];

    if (T != T_save)
    {
        T_save = T;
        comp_k_f(tc,invT,k_f_save);
        comp_Kc(tc,invT,Kc_save);
    }

    comp_qfqr(q_f, q_r, sc, tc, invT);

    return;
}


/*compute the equilibrium constants for each reaction */
void equilibriumConstants(double * restrict kc, double * restrict g_RT, double T)
{
    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 / T;

    /*reaction 1: O2 + H + M <=> HO2 + M */
    kc[0] = 1.0 / (refC) * exp((g_RT[1] + g_RT[9]) - (g_RT[8]));

    /*reaction 2: H2O2 + M <=> 2 OH + M */
    kc[1] = refC * exp((g_RT[7]) - (2 * g_RT[6]));

    /*reaction 3: HCO + M <=> CO + H + M */
    kc[2] = refC * exp((g_RT[13]) - (g_RT[4] + g_RT[9]));

    /*reaction 4: NC7H16 + H <=> H2 + C7H15-2 */
    kc[3] = exp((g_RT[0] + g_RT[9]) - (g_RT[5] + g_RT[24]));

    /*reaction 5: NC7H16 + OH <=> H2O + C7H15-2 */
    kc[4] = exp((g_RT[0] + g_RT[6]) - (g_RT[3] + g_RT[24]));

    /*reaction 6: NC7H16 + HO2 <=> H2O2 + C7H15-2 */
    kc[5] = exp((g_RT[0] + g_RT[8]) - (g_RT[7] + g_RT[24]));

    /*reaction 7: NC7H16 + O2 <=> HO2 + C7H15-2 */
    kc[6] = exp((g_RT[0] + g_RT[1]) - (g_RT[8] + g_RT[24]));

    /*reaction 8: O2 + C7H15-2 <=> C7H15O2 */
    kc[7] = 1.0 / (refC) * exp((g_RT[1] + g_RT[24]) - (g_RT[25]));

    /*reaction 9: O2 + C7H15O2 <=> OH + C7KET12 */
    kc[8] = exp((g_RT[1] + g_RT[25]) - (g_RT[6] + g_RT[26]));

    /*reaction 10: C7KET12 <=> OH + CH2O + C5H11CO */
    kc[9] = refC*refC * exp((g_RT[26]) - (g_RT[6] + g_RT[12] + g_RT[27]));

    /*reaction 11: C5H11CO <=> CO + C2H4 + C3H7 */
    kc[10] = refC*refC * exp((g_RT[27]) - (g_RT[4] + g_RT[18] + g_RT[23]));

    /*reaction 12: C7H15-2 <=> C2H4 + C2H5 + C3H6 */
    kc[11] = refC*refC * exp((g_RT[24]) - (g_RT[18] + g_RT[19] + g_RT[22]));

    /*reaction 13: C3H7 <=> CH3 + C2H4 */
    kc[12] = refC * exp((g_RT[23]) - (g_RT[15] + g_RT[18]));

    /*reaction 14: C3H7 <=> H + C3H6 */
    kc[13] = refC * exp((g_RT[23]) - (g_RT[9] + g_RT[22]));

    /*reaction 15: CH3 + C3H6 <=> CH4 + C3H5 */
    kc[14] = exp((g_RT[15] + g_RT[22]) - (g_RT[16] + g_RT[21]));

    /*reaction 16: O2 + C3H5 <=> HO2 + C3H4 */
    kc[15] = exp((g_RT[1] + g_RT[21]) - (g_RT[8] + g_RT[20]));

    /*reaction 17: OH + C3H4 <=> CH2O + C2H3 */
    kc[16] = exp((g_RT[6] + g_RT[20]) - (g_RT[12] + g_RT[17]));

    /*reaction 18: OH + C3H4 <=> HCO + C2H4 */
    kc[17] = exp((g_RT[6] + g_RT[20]) - (g_RT[13] + g_RT[18]));

    /*reaction 19: HO2 + CH3 <=> OH + CH3O */
    kc[18] = exp((g_RT[8] + g_RT[15]) - (g_RT[6] + g_RT[11]));

    /*reaction 20: OH + CH3 <=> H2O + CH2 */
    kc[19] = exp((g_RT[6] + g_RT[15]) - (g_RT[3] + g_RT[14]));

    /*reaction 21: OH + CH2 <=> H + CH2O */
    kc[20] = exp((g_RT[6] + g_RT[14]) - (g_RT[9] + g_RT[12]));

    /*reaction 22: O2 + CH2 <=> OH + HCO */
    kc[21] = exp((g_RT[1] + g_RT[14]) - (g_RT[6] + g_RT[13]));

    /*reaction 23: O2 + CH2 <=> CO2 + H2 */
    kc[22] = exp((g_RT[1] + g_RT[14]) - (g_RT[2] + g_RT[5]));

    /*reaction 24: O2 + CH2 <=> H2O + CO */
    kc[23] = exp((g_RT[1] + g_RT[14]) - (g_RT[3] + g_RT[4]));

    /*reaction 25: O2 + CH2 <=> O + CH2O */
    kc[24] = exp((g_RT[1] + g_RT[14]) - (g_RT[10] + g_RT[12]));

    /*reaction 26: O2 + CH2 <=> CO2 + 2 H */
    kc[25] = refC * exp((g_RT[1] + g_RT[14]) - (g_RT[2] + 2 * g_RT[9]));

    /*reaction 27: O2 + CH2 <=> CO + OH + H */
    kc[26] = refC * exp((g_RT[1] + g_RT[14]) - (g_RT[4] + g_RT[6] + g_RT[9]));

    /*reaction 28: CO + CH3O <=> CO2 + CH3 */
    kc[27] = exp((g_RT[4] + g_RT[11]) - (g_RT[2] + g_RT[15]));

    /*reaction 29: CO + OH <=> CO2 + H */
    kc[28] = exp((g_RT[4] + g_RT[6]) - (g_RT[2] + g_RT[9]));

    /*reaction 30: OH + O <=> O2 + H */
    kc[29] = exp((g_RT[6] + g_RT[10]) - (g_RT[1] + g_RT[9]));

    /*reaction 31: HO2 + H <=> 2 OH */
    kc[30] = exp((g_RT[8] + g_RT[9]) - (2 * g_RT[6]));

    /*reaction 32: 2 OH <=> H2O + O */
    kc[31] = exp((2 * g_RT[6]) - (g_RT[3] + g_RT[10]));

    /*reaction 33: H2 + OH <=> H2O + H */
    kc[32] = exp((g_RT[5] + g_RT[6]) - (g_RT[3] + g_RT[9]));

    /*reaction 34: 2 HO2 <=> O2 + H2O2 */
    kc[33] = exp((2 * g_RT[8]) - (g_RT[1] + g_RT[7]));

    /*reaction 35: OH + CH2O <=> H2O + HCO */
    kc[34] = exp((g_RT[6] + g_RT[12]) - (g_RT[3] + g_RT[13]));

    /*reaction 36: HO2 + CH2O <=> H2O2 + HCO */
    kc[35] = exp((g_RT[8] + g_RT[12]) - (g_RT[7] + g_RT[13]));

    /*reaction 37: O2 + HCO <=> CO + HO2 */
    kc[36] = exp((g_RT[1] + g_RT[13]) - (g_RT[4] + g_RT[8]));

    /*reaction 38: CH3O + CH3 <=> CH2O + CH4 */
    kc[37] = exp((g_RT[11] + g_RT[15]) - (g_RT[12] + g_RT[16]));

    /*reaction 39: OH + C2H4 <=> CH2O + CH3 */
    kc[38] = exp((g_RT[6] + g_RT[18]) - (g_RT[12] + g_RT[15]));

    /*reaction 40: OH + C2H4 <=> H2O + C2H3 */
    kc[39] = exp((g_RT[6] + g_RT[18]) - (g_RT[3] + g_RT[17]));

    /*reaction 41: O2 + C2H3 <=> CH2O + HCO */
    kc[40] = exp((g_RT[1] + g_RT[17]) - (g_RT[12] + g_RT[13]));

    /*reaction 42: HCO + C2H3 <=> CO + C2H4 */
    kc[41] = exp((g_RT[13] + g_RT[17]) - (g_RT[4] + g_RT[18]));

    /*reaction 43: O2 + C2H5 <=> HO2 + C2H4 */
    kc[42] = exp((g_RT[1] + g_RT[19]) - (g_RT[8] + g_RT[18]));

    /*reaction 44: O2 + CH4 <=> HO2 + CH3 */
    kc[43] = exp((g_RT[1] + g_RT[16]) - (g_RT[8] + g_RT[15]));

    /*reaction 45: OH + HO2 <=> O2 + H2O */
    kc[44] = exp((g_RT[6] + g_RT[8]) - (g_RT[1] + g_RT[3]));

    /*reaction 46: O2 + CH3 <=> OH + CH2O */
    kc[45] = exp((g_RT[1] + g_RT[15]) - (g_RT[6] + g_RT[12]));

    /*reaction 47: H + CH4 <=> H2 + CH3 */
    kc[46] = exp((g_RT[9] + g_RT[16]) - (g_RT[5] + g_RT[15]));

    /*reaction 48: OH + CH4 <=> H2O + CH3 */
    kc[47] = exp((g_RT[6] + g_RT[16]) - (g_RT[3] + g_RT[15]));

    /*reaction 49: O + CH4 <=> OH + CH3 */
    kc[48] = exp((g_RT[10] + g_RT[16]) - (g_RT[6] + g_RT[15]));

    /*reaction 50: HO2 + CH4 <=> H2O2 + CH3 */
    kc[49] = exp((g_RT[8] + g_RT[16]) - (g_RT[7] + g_RT[15]));

    /*reaction 51: CH2 + CH4 <=> 2 CH3 */
    kc[50] = exp((g_RT[14] + g_RT[16]) - (2 * g_RT[15]));

    /*reaction 52: C3H6 <=> CH3 + C2H3 */
    kc[51] = refC * exp((g_RT[22]) - (g_RT[15] + g_RT[17]));

    return;
}


/*compute the g/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void gibbs(double * restrict species, double * restrict tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=1382 kelvin */
    if (T < 1382) {
        /*species 24: C7H15-2 */
        species[24] =
            -2.356053030000000e+03 * invT
            -3.377006617670000e+01
            +3.791557670000000e-02 * tc[0]
            -3.783632850000000e-02 * tc[1]
            +6.791227233333333e-06 * tc[2]
            -7.772324525000000e-10 * tc[3]
            +2.461803725000000e-14 * tc[4];
    } else {
        /*species 24: C7H15-2 */
        species[24] =
            -1.058736160000000e+04 * invT
            +1.068578495000000e+02
            -2.163688420000000e+01 * tc[0]
            -1.616624020000000e-02 * tc[1]
            +1.821230116666667e-06 * tc[2]
            -1.402975500000000e-10 * tc[3]
            +4.858870455000000e-15 * tc[4];
    }

    /*species with midpoint at T=1383 kelvin */
    if (T < 1383) {
        /*species 27: C5H11CO */
        species[27] =
            -1.434511720000000e+04 * invT
            -2.016801381000000e+01
            -2.144790690000000e+00 * tc[0]
            -3.089317815000000e-02 * tc[1]
            +6.235578166666666e-06 * tc[2]
            -9.440316250000000e-10 * tc[3]
            +6.845884899999999e-14 * tc[4];
    } else {
        /*species 27: C5H11CO */
        species[27] =
            -2.079239370000000e+04 * invT
            +9.167793899999999e+01
            -1.947838120000000e+01 * tc[0]
            -1.252330145000000e-02 * tc[1]
            +1.424768910000000e-06 * tc[2]
            -1.104649533333333e-10 * tc[3]
            +3.842516480000000e-15 * tc[4];
    }

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 1: O2 */
        species[1] =
            -1.005249000000000e+03 * invT
            -2.821802000000000e+00
            -3.212936000000000e+00 * tc[0]
            -5.637430000000000e-04 * tc[1]
            +9.593583333333333e-08 * tc[2]
            -1.094897500000000e-10 * tc[3]
            +4.384277000000000e-14 * tc[4];
        /*species 2: CO2 */
        species[2] =
            -4.837314000000000e+04 * invT
            -7.912765000000000e+00
            -2.275725000000000e+00 * tc[0]
            -4.961036000000000e-03 * tc[1]
            +1.734851666666667e-06 * tc[2]
            -5.722239166666667e-10 * tc[3]
            +1.058640000000000e-13 * tc[4];
        /*species 3: H2O */
        species[3] =
            -3.020811000000000e+04 * invT
            +7.966090000000001e-01
            -3.386842000000000e+00 * tc[0]
            -1.737491000000000e-03 * tc[1]
            +1.059116000000000e-06 * tc[2]
            -5.807150833333333e-10 * tc[3]
            +1.253294000000000e-13 * tc[4];
        /*species 4: CO */
        species[4] =
            -1.431054000000000e+04 * invT
            -1.586445000000000e+00
            -3.262452000000000e+00 * tc[0]
            -7.559705000000000e-04 * tc[1]
            +6.469591666666667e-07 * tc[2]
            -4.651620000000000e-10 * tc[3]
            +1.237475500000000e-13 * tc[4];
        /*species 5: H2 */
        species[5] =
            -1.012521000000000e+03 * invT
            +6.592218000000000e+00
            -3.298124000000000e+00 * tc[0]
            -4.124721000000000e-04 * tc[1]
            +1.357169166666667e-07 * tc[2]
            +7.896194999999999e-12 * tc[3]
            -2.067436000000000e-14 * tc[4];
        /*species 6: OH */
        species[6] =
            +3.606782000000000e+03 * invT
            +2.278406000000000e+00
            -3.637266000000000e+00 * tc[0]
            -9.254550000000000e-05 * tc[1]
            +2.793608333333333e-07 * tc[2]
            -1.989335833333333e-10 * tc[3]
            +4.215721000000000e-14 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.766315000000000e+04 * invT
            -3.396609000000000e+00
            -3.388754000000000e+00 * tc[0]
            -3.284613000000000e-03 * tc[1]
            +2.475021666666666e-08 * tc[2]
            +3.854838333333333e-10 * tc[3]
            -1.235757500000000e-13 * tc[4];
        /*species 8: HO2 */
        species[8] =
            +1.762274000000000e+02 * invT
            -6.242761000000000e+00
            -2.979963000000000e+00 * tc[0]
            -2.498348500000000e-03 * tc[1]
            +6.318328333333333e-07 * tc[2]
            -1.961826666666667e-10 * tc[3]
            +4.044512000000000e-14 * tc[4];
        /*species 9: H */
        species[9] =
            +2.547163000000000e+04 * invT
            +2.960117600000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 10: O */
        species[10] =
            +2.914764000000000e+04 * invT
            -1.756599999999997e-02
            -2.946429000000000e+00 * tc[0]
            +8.190830000000000e-04 * tc[1]
            -4.035053333333333e-07 * tc[2]
            +1.335702500000000e-10 * tc[3]
            -1.945348000000000e-14 * tc[4];
        /*species 11: CH3O */
        species[11] =
            +9.786011000000000e+02 * invT
            -1.104597600000000e+01
            -2.106204000000000e+00 * tc[0]
            -3.608297500000000e-03 * tc[1]
            -8.897453333333333e-07 * tc[2]
            +6.148030000000000e-10 * tc[3]
            -1.037805500000000e-13 * tc[4];
        /*species 12: CH2O */
        species[12] =
            -1.486540000000000e+04 * invT
            -1.213208900000000e+01
            -1.652731000000000e+00 * tc[0]
            -6.315720000000000e-03 * tc[1]
            +3.146946666666667e-06 * tc[2]
            -1.708359166666667e-09 * tc[3]
            +4.206618500000000e-13 * tc[4];
        /*species 13: HCO */
        species[13] =
            +4.159922000000000e+03 * invT
            -6.085284000000000e+00
            -2.898330000000000e+00 * tc[0]
            -3.099573500000000e-03 * tc[1]
            +1.603847333333333e-06 * tc[2]
            -9.081875000000000e-10 * tc[3]
            +2.287442500000000e-13 * tc[4];
        /*species 14: CH2 */
        species[14] =
            +4.536791000000000e+04 * invT
            +2.049659000000000e+00
            -3.762237000000000e+00 * tc[0]
            -5.799095000000000e-04 * tc[1]
            -4.149308333333333e-08 * tc[2]
            -7.334030000000001e-11 * tc[3]
            +3.666217500000000e-14 * tc[4];
        /*species 15: CH3 */
        species[15] =
            +1.642378000000000e+04 * invT
            -4.359351000000000e+00
            -2.430443000000000e+00 * tc[0]
            -5.562050000000000e-03 * tc[1]
            +2.800366666666666e-06 * tc[2]
            -1.351524166666667e-09 * tc[3]
            +2.932476500000000e-13 * tc[4];
        /*species 16: CH4 */
        species[16] =
            -9.825228999999999e+03 * invT
            -1.294344850000000e+01
            -7.787415000000000e-01 * tc[0]
            -8.738340000000001e-03 * tc[1]
            +4.639015000000000e-06 * tc[2]
            -2.541423333333333e-09 * tc[3]
            +6.119655000000000e-13 * tc[4];
        /*species 17: C2H3 */
        species[17] =
            +3.335225000000000e+04 * invT
            -9.096924000000001e+00
            -2.459276000000000e+00 * tc[0]
            -3.685738000000000e-03 * tc[1]
            -3.516455000000000e-07 * tc[2]
            +1.101368333333333e-10 * tc[3]
            +5.923920000000000e-14 * tc[4];
        /*species 18: C2H4 */
        species[18] =
            +5.573046000000000e+03 * invT
            -2.507297800000000e+01
            +8.614880000000000e-01 * tc[0]
            -1.398081500000000e-02 * tc[1]
            +5.647795000000000e-06 * tc[2]
            -2.320960000000000e-09 * tc[3]
            +4.868939500000000e-13 * tc[4];
        /*species 19: C2H5 */
        species[19] =
            +1.287040000000000e+04 * invT
            -9.447498000000000e+00
            -2.690702000000000e+00 * tc[0]
            -4.359566500000000e-03 * tc[1]
            -7.366398333333332e-07 * tc[2]
            -7.782252500000001e-11 * tc[3]
            +1.963886500000000e-13 * tc[4];
        /*species 20: C3H4 */
        species[20] =
            +2.154156420000000e+04 * invT
            -7.637257030000001e+00
            -2.613074870000000e+00 * tc[0]
            -6.061168550000000e-03 * tc[1]
            -3.090090000000000e-06 * tc[2]
            +2.877153958333333e-09 * tc[3]
            -7.667669450000000e-13 * tc[4];
        /*species 21: C3H5 */
        species[21] =
            +1.862612180000000e+04 * invT
            -4.040278060000000e+00
            -3.787946930000000e+00 * tc[0]
            -4.742071675000000e-03 * tc[1]
            -4.039056133333334e-06 * tc[2]
            +3.046700083333334e-09 * tc[3]
            -7.429617800000000e-13 * tc[4];
        /*species 23: C3H7 */
        species[23] =
            +1.063186300000000e+04 * invT
            -2.007100720000000e+01
            -1.051551800000000e+00 * tc[0]
            -1.299599000000000e-02 * tc[1]
            -3.966756666666666e-07 * tc[2]
            +1.634130750000000e-09 * tc[3]
            -4.686623500000000e-13 * tc[4];
        /*species 28: N2 */
        species[28] =
            -1.020900000000000e+03 * invT
            -6.516950000000001e-01
            -3.298677000000000e+00 * tc[0]
            -7.041200000000000e-04 * tc[1]
            +6.605369999999999e-07 * tc[2]
            -4.701262500000001e-10 * tc[3]
            +1.222427500000000e-13 * tc[4];
    } else {
        /*species 1: O2 */
        species[1] =
            -1.233930000000000e+03 * invT
            +5.084119999999999e-01
            -3.697578000000000e+00 * tc[0]
            -3.067598500000000e-04 * tc[1]
            +2.098070000000000e-08 * tc[2]
            -1.479400833333333e-12 * tc[3]
            +5.682175000000001e-17 * tc[4];
        /*species 2: CO2 */
        species[2] =
            -4.896696000000000e+04 * invT
            +5.409018900000000e+00
            -4.453623000000000e+00 * tc[0]
            -1.570084500000000e-03 * tc[1]
            +2.130685000000000e-07 * tc[2]
            -1.994997500000000e-11 * tc[3]
            +8.345165000000000e-16 * tc[4];
        /*species 3: H2O */
        species[3] =
            -2.989921000000000e+04 * invT
            -4.190671000000000e+00
            -2.672146000000000e+00 * tc[0]
            -1.528146500000000e-03 * tc[1]
            +1.455043333333333e-07 * tc[2]
            -1.000830000000000e-11 * tc[3]
            +3.195809000000000e-16 * tc[4];
        /*species 4: CO */
        species[4] =
            -1.426835000000000e+04 * invT
            -3.083140000000000e+00
            -3.025078000000000e+00 * tc[0]
            -7.213445000000000e-04 * tc[1]
            +9.384713333333334e-08 * tc[2]
            -8.488174999999999e-12 * tc[3]
            +3.455476000000000e-16 * tc[4];
        /*species 5: H2 */
        species[5] =
            -8.350340000000000e+02 * invT
            +4.346533000000000e+00
            -2.991423000000000e+00 * tc[0]
            -3.500322000000000e-04 * tc[1]
            +9.389715000000000e-09 * tc[2]
            +7.692981666666667e-13 * tc[3]
            -7.913760000000000e-17 * tc[4];
        /*species 6: OH */
        species[6] =
            +3.886888000000000e+03 * invT
            -2.712982000000000e+00
            -2.882730000000000e+00 * tc[0]
            -5.069870000000000e-04 * tc[1]
            +3.794795000000000e-08 * tc[2]
            -1.812236666666667e-12 * tc[3]
            +2.563152500000000e-17 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.800696000000000e+04 * invT
            +4.072030000000000e+00
            -4.573167000000000e+00 * tc[0]
            -2.168068000000000e-03 * tc[1]
            +2.457815000000000e-07 * tc[2]
            -1.957420000000000e-11 * tc[3]
            +7.158270000000000e-16 * tc[4];
        /*species 8: HO2 */
        species[8] =
            -1.579727000000000e+02 * invT
            +5.961620000000001e-01
            -4.072191000000000e+00 * tc[0]
            -1.065648000000000e-03 * tc[1]
            +8.846908333333335e-08 * tc[2]
            -5.093557500000001e-12 * tc[3]
            +1.420582500000000e-16 * tc[4];
        /*species 9: H */
        species[9] =
            +2.547163000000000e+04 * invT
            +2.960117600000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 10: O */
        species[10] =
            +2.923080000000000e+04 * invT
            -2.378248000000000e+00
            -2.542060000000000e+00 * tc[0]
            +1.377531000000000e-05 * tc[1]
            +5.171338333333333e-10 * tc[2]
            -3.792555833333334e-13 * tc[3]
            +2.184026000000000e-17 * tc[4];
        /*species 11: CH3O */
        species[11] =
            +1.278325000000000e+02 * invT
            +8.412250000000001e-01
            -3.770800000000000e+00 * tc[0]
            -3.935748500000000e-03 * tc[1]
            +4.427306666666667e-07 * tc[2]
            -3.287025833333333e-11 * tc[3]
            +1.056308000000000e-15 * tc[4];
        /*species 12: CH2O */
        species[12] =
            -1.532037000000000e+04 * invT
            -3.916966000000000e+00
            -2.995606000000000e+00 * tc[0]
            -3.340660500000000e-03 * tc[1]
            +4.381591666666666e-07 * tc[2]
            -3.947627500000000e-11 * tc[3]
            +1.606258500000000e-15 * tc[4];
        /*species 13: HCO */
        species[13] =
            +3.916324000000000e+03 * invT
            -1.995028000000000e+00
            -3.557271000000000e+00 * tc[0]
            -1.672786500000000e-03 * tc[1]
            +2.225010000000000e-07 * tc[2]
            -2.058810833333333e-11 * tc[3]
            +8.569255000000000e-16 * tc[4];
        /*species 14: CH2 */
        species[14] =
            +4.534134000000000e+04 * invT
            +1.479847000000000e+00
            -3.636408000000000e+00 * tc[0]
            -9.665285000000000e-04 * tc[1]
            +2.811693333333333e-08 * tc[2]
            +8.415825000000000e-12 * tc[3]
            -9.041279999999999e-16 * tc[4];
        /*species 15: CH3 */
        species[15] =
            +1.643781000000000e+04 * invT
            -2.608645000000000e+00
            -2.844052000000000e+00 * tc[0]
            -3.068987000000000e-03 * tc[1]
            +3.717241666666666e-07 * tc[2]
            -3.154300833333333e-11 * tc[3]
            +1.226079500000000e-15 * tc[4];
        /*species 16: CH4 */
        species[16] =
            -1.008079000000000e+04 * invT
            -7.939916000000000e+00
            -1.683479000000000e+00 * tc[0]
            -5.118620000000000e-03 * tc[1]
            +6.458548333333333e-07 * tc[2]
            -5.654654166666667e-11 * tc[3]
            +2.251711500000000e-15 * tc[4];
        /*species 17: C2H3 */
        species[17] =
            +3.185435000000000e+04 * invT
            +1.446378100000000e+01
            -5.933468000000000e+00 * tc[0]
            -2.008873000000000e-03 * tc[1]
            +6.611233333333333e-08 * tc[2]
            +1.201055833333333e-11 * tc[3]
            -1.189322000000000e-15 * tc[4];
        /*species 18: C2H4 */
        species[18] =
            +4.428289000000000e+03 * invT
            +1.298030000000000e+00
            -3.528419000000000e+00 * tc[0]
            -5.742590000000000e-03 * tc[1]
            +7.363975000000000e-07 * tc[2]
            -6.537167500000001e-11 * tc[3]
            +2.633424000000000e-15 * tc[4];
        /*species 19: C2H5 */
        species[19] =
            +1.067455000000000e+04 * invT
            +2.197137000000000e+01
            -7.190480000000000e+00 * tc[0]
            -3.242038500000000e-03 * tc[1]
            +1.071344166666667e-07 * tc[2]
            +1.956565833333333e-11 * tc[3]
            -1.940438500000000e-15 * tc[4];
        /*species 20: C3H4 */
        species[20] =
            +2.011746170000000e+04 * invT
            +1.728883489000000e+01
            -6.316948690000000e+00 * tc[0]
            -5.566813100000000e-03 * tc[1]
            +6.604816966666667e-07 * tc[2]
            -5.296948125000000e-11 * tc[3]
            +1.893749425000000e-15 * tc[4];
        /*species 21: C3H5 */
        species[21] =
            +1.727147070000000e+04 * invT
            +1.582247973000000e+01
            -6.547611320000000e+00 * tc[0]
            -6.657612300000000e-03 * tc[1]
            +7.972218333333334e-07 * tc[2]
            -6.432915116666666e-11 * tc[3]
            +2.309654040000000e-15 * tc[4];
        /*species 23: C3H7 */
        species[23] =
            +8.298433600000000e+03 * invT
            +2.318287870000000e+01
            -7.702698700000000e+00 * tc[0]
            -8.022101500000000e-03 * tc[1]
            +8.805536666666666e-07 * tc[2]
            -6.358215833333333e-11 * tc[3]
            +1.969614200000000e-15 * tc[4];
        /*species 28: N2 */
        species[28] =
            -9.227977000000000e+02 * invT
            -3.053888000000000e+00
            -2.926640000000000e+00 * tc[0]
            -7.439885000000000e-04 * tc[1]
            +9.474601666666666e-08 * tc[2]
            -8.414199999999999e-12 * tc[3]
            +3.376675500000000e-16 * tc[4];
    }

    /*species with midpoint at T=1388 kelvin */
    if (T < 1388) {
        /*species 22: C3H6 */
        species[22] =
            +1.066881640000000e+03 * invT
            -2.150575815600000e+01
            -3.946154440000000e-01 * tc[0]
            -1.445538310000000e-02 * tc[1]
            +2.581446800000000e-06 * tc[2]
            -3.240118408333333e-10 * tc[3]
            +1.689451760000000e-14 * tc[4];
    } else {
        /*species 22: C3H6 */
        species[22] =
            -1.878212710000000e+03 * invT
            +2.803202638000000e+01
            -8.015959580000001e+00 * tc[0]
            -6.851181700000000e-03 * tc[1]
            +7.770828883333333e-07 * tc[2]
            -6.010453350000000e-11 * tc[3]
            +2.086850630000000e-15 * tc[4];
    }

    /*species with midpoint at T=1390 kelvin */
    if (T < 1390) {
        /*species 25: C7H15O2 */
        species[25] =
            -1.992379610000000e+04 * invT
            -2.293174086000000e+01
            -2.374993340000000e+00 * tc[0]
            -4.173259530000000e-02 * tc[1]
            +8.564955333333334e-06 * tc[2]
            -1.368480516666667e-09 * tc[3]
            +1.097526080000000e-13 * tc[4];
    } else {
        /*species 25: C7H15O2 */
        species[25] =
            -2.829760500000000e+04 * invT
            +1.222947231000000e+02
            -2.490236890000000e+01 * tc[0]
            -1.753584600000000e-02 * tc[1]
            +2.007338433333334e-06 * tc[2]
            -1.562206850000000e-10 * tc[3]
            +5.447389550000000e-15 * tc[4];
    }

    /*species with midpoint at T=1391 kelvin */
    if (T < 1391) {
        /*species 0: NC7H16 */
        species[0] =
            -2.565865650000000e+04 * invT
            -3.664165307000000e+01
            +1.268361870000000e+00 * tc[0]
            -4.271779100000000e-02 * tc[1]
            +8.755779766666667e-06 * tc[2]
            -1.357881008333333e-09 * tc[3]
            +1.011974625000000e-13 * tc[4];
    } else {
        /*species 0: NC7H16 */
        species[0] =
            -3.427600810000000e+04 * invT
            +1.145189165000000e+02
            -2.221489690000000e+01 * tc[0]
            -1.738378750000000e-02 * tc[1]
            +1.973452150000000e-06 * tc[2]
            -1.527487316666667e-10 * tc[3]
            +5.306513300000000e-15 * tc[4];
    }

    /*species with midpoint at T=1396 kelvin */
    if (T < 1396) {
        /*species 26: C7KET12 */
        species[26] =
            -4.680544190000000e+04 * invT
            -3.275071120300000e+01
            -5.824336970000000e-01 * tc[0]
            -5.060393450000000e-02 * tc[1]
            +1.276426660000000e-05 * tc[2]
            -2.506155050000000e-09 * tc[3]
            +2.414513960000000e-13 * tc[4];
    } else {
        /*species 26: C7KET12 */
        species[26] =
            -5.668568280000000e+04 * invT
            +1.521797806000000e+02
            -2.974729060000000e+01 * tc[0]
            -1.533111470000000e-02 * tc[1]
            +1.759393166666667e-06 * tc[2]
            -1.371894525000000e-10 * tc[3]
            +4.790858375000000e-15 * tc[4];
    }
    return;
}


/*compute the a/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void helmholtz(double * restrict species, double * restrict tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=1382 kelvin */
    if (T < 1382) {
        /*species 24: C7H15-2 */
        species[24] =
            -2.35605303e+03 * invT
            -3.47700662e+01
            +3.79155767e-02 * tc[0]
            -3.78363285e-02 * tc[1]
            +6.79122723e-06 * tc[2]
            -7.77232453e-10 * tc[3]
            +2.46180373e-14 * tc[4];
    } else {
        /*species 24: C7H15-2 */
        species[24] =
            -1.05873616e+04 * invT
            +1.05857850e+02
            -2.16368842e+01 * tc[0]
            -1.61662402e-02 * tc[1]
            +1.82123012e-06 * tc[2]
            -1.40297550e-10 * tc[3]
            +4.85887045e-15 * tc[4];
    }

    /*species with midpoint at T=1383 kelvin */
    if (T < 1383) {
        /*species 27: C5H11CO */
        species[27] =
            -1.43451172e+04 * invT
            -2.11680138e+01
            -2.14479069e+00 * tc[0]
            -3.08931781e-02 * tc[1]
            +6.23557817e-06 * tc[2]
            -9.44031625e-10 * tc[3]
            +6.84588490e-14 * tc[4];
    } else {
        /*species 27: C5H11CO */
        species[27] =
            -2.07923937e+04 * invT
            +9.06779390e+01
            -1.94783812e+01 * tc[0]
            -1.25233014e-02 * tc[1]
            +1.42476891e-06 * tc[2]
            -1.10464953e-10 * tc[3]
            +3.84251648e-15 * tc[4];
    }

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 1: O2 */
        species[1] =
            -1.00524900e+03 * invT
            -3.82180200e+00
            -3.21293600e+00 * tc[0]
            -5.63743000e-04 * tc[1]
            +9.59358333e-08 * tc[2]
            -1.09489750e-10 * tc[3]
            +4.38427700e-14 * tc[4];
        /*species 2: CO2 */
        species[2] =
            -4.83731400e+04 * invT
            -8.91276500e+00
            -2.27572500e+00 * tc[0]
            -4.96103600e-03 * tc[1]
            +1.73485167e-06 * tc[2]
            -5.72223917e-10 * tc[3]
            +1.05864000e-13 * tc[4];
        /*species 3: H2O */
        species[3] =
            -3.02081100e+04 * invT
            -2.03391000e-01
            -3.38684200e+00 * tc[0]
            -1.73749100e-03 * tc[1]
            +1.05911600e-06 * tc[2]
            -5.80715083e-10 * tc[3]
            +1.25329400e-13 * tc[4];
        /*species 4: CO */
        species[4] =
            -1.43105400e+04 * invT
            -2.58644500e+00
            -3.26245200e+00 * tc[0]
            -7.55970500e-04 * tc[1]
            +6.46959167e-07 * tc[2]
            -4.65162000e-10 * tc[3]
            +1.23747550e-13 * tc[4];
        /*species 5: H2 */
        species[5] =
            -1.01252100e+03 * invT
            +5.59221800e+00
            -3.29812400e+00 * tc[0]
            -4.12472100e-04 * tc[1]
            +1.35716917e-07 * tc[2]
            +7.89619500e-12 * tc[3]
            -2.06743600e-14 * tc[4];
        /*species 6: OH */
        species[6] =
            +3.60678200e+03 * invT
            +1.27840600e+00
            -3.63726600e+00 * tc[0]
            -9.25455000e-05 * tc[1]
            +2.79360833e-07 * tc[2]
            -1.98933583e-10 * tc[3]
            +4.21572100e-14 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.76631500e+04 * invT
            -4.39660900e+00
            -3.38875400e+00 * tc[0]
            -3.28461300e-03 * tc[1]
            +2.47502167e-08 * tc[2]
            +3.85483833e-10 * tc[3]
            -1.23575750e-13 * tc[4];
        /*species 8: HO2 */
        species[8] =
            +1.76227400e+02 * invT
            -7.24276100e+00
            -2.97996300e+00 * tc[0]
            -2.49834850e-03 * tc[1]
            +6.31832833e-07 * tc[2]
            -1.96182667e-10 * tc[3]
            +4.04451200e-14 * tc[4];
        /*species 9: H */
        species[9] =
            +2.54716300e+04 * invT
            +1.96011760e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 10: O */
        species[10] =
            +2.91476400e+04 * invT
            -1.01756600e+00
            -2.94642900e+00 * tc[0]
            +8.19083000e-04 * tc[1]
            -4.03505333e-07 * tc[2]
            +1.33570250e-10 * tc[3]
            -1.94534800e-14 * tc[4];
        /*species 11: CH3O */
        species[11] =
            +9.78601100e+02 * invT
            -1.20459760e+01
            -2.10620400e+00 * tc[0]
            -3.60829750e-03 * tc[1]
            -8.89745333e-07 * tc[2]
            +6.14803000e-10 * tc[3]
            -1.03780550e-13 * tc[4];
        /*species 12: CH2O */
        species[12] =
            -1.48654000e+04 * invT
            -1.31320890e+01
            -1.65273100e+00 * tc[0]
            -6.31572000e-03 * tc[1]
            +3.14694667e-06 * tc[2]
            -1.70835917e-09 * tc[3]
            +4.20661850e-13 * tc[4];
        /*species 13: HCO */
        species[13] =
            +4.15992200e+03 * invT
            -7.08528400e+00
            -2.89833000e+00 * tc[0]
            -3.09957350e-03 * tc[1]
            +1.60384733e-06 * tc[2]
            -9.08187500e-10 * tc[3]
            +2.28744250e-13 * tc[4];
        /*species 14: CH2 */
        species[14] =
            +4.53679100e+04 * invT
            +1.04965900e+00
            -3.76223700e+00 * tc[0]
            -5.79909500e-04 * tc[1]
            -4.14930833e-08 * tc[2]
            -7.33403000e-11 * tc[3]
            +3.66621750e-14 * tc[4];
        /*species 15: CH3 */
        species[15] =
            +1.64237800e+04 * invT
            -5.35935100e+00
            -2.43044300e+00 * tc[0]
            -5.56205000e-03 * tc[1]
            +2.80036667e-06 * tc[2]
            -1.35152417e-09 * tc[3]
            +2.93247650e-13 * tc[4];
        /*species 16: CH4 */
        species[16] =
            -9.82522900e+03 * invT
            -1.39434485e+01
            -7.78741500e-01 * tc[0]
            -8.73834000e-03 * tc[1]
            +4.63901500e-06 * tc[2]
            -2.54142333e-09 * tc[3]
            +6.11965500e-13 * tc[4];
        /*species 17: C2H3 */
        species[17] =
            +3.33522500e+04 * invT
            -1.00969240e+01
            -2.45927600e+00 * tc[0]
            -3.68573800e-03 * tc[1]
            -3.51645500e-07 * tc[2]
            +1.10136833e-10 * tc[3]
            +5.92392000e-14 * tc[4];
        /*species 18: C2H4 */
        species[18] =
            +5.57304600e+03 * invT
            -2.60729780e+01
            +8.61488000e-01 * tc[0]
            -1.39808150e-02 * tc[1]
            +5.64779500e-06 * tc[2]
            -2.32096000e-09 * tc[3]
            +4.86893950e-13 * tc[4];
        /*species 19: C2H5 */
        species[19] =
            +1.28704000e+04 * invT
            -1.04474980e+01
            -2.69070200e+00 * tc[0]
            -4.35956650e-03 * tc[1]
            -7.36639833e-07 * tc[2]
            -7.78225250e-11 * tc[3]
            +1.96388650e-13 * tc[4];
        /*species 20: C3H4 */
        species[20] =
            +2.15415642e+04 * invT
            -8.63725703e+00
            -2.61307487e+00 * tc[0]
            -6.06116855e-03 * tc[1]
            -3.09009000e-06 * tc[2]
            +2.87715396e-09 * tc[3]
            -7.66766945e-13 * tc[4];
        /*species 21: C3H5 */
        species[21] =
            +1.86261218e+04 * invT
            -5.04027806e+00
            -3.78794693e+00 * tc[0]
            -4.74207167e-03 * tc[1]
            -4.03905613e-06 * tc[2]
            +3.04670008e-09 * tc[3]
            -7.42961780e-13 * tc[4];
        /*species 23: C3H7 */
        species[23] =
            +1.06318630e+04 * invT
            -2.10710072e+01
            -1.05155180e+00 * tc[0]
            -1.29959900e-02 * tc[1]
            -3.96675667e-07 * tc[2]
            +1.63413075e-09 * tc[3]
            -4.68662350e-13 * tc[4];
        /*species 28: N2 */
        species[28] =
            -1.02090000e+03 * invT
            -1.65169500e+00
            -3.29867700e+00 * tc[0]
            -7.04120000e-04 * tc[1]
            +6.60537000e-07 * tc[2]
            -4.70126250e-10 * tc[3]
            +1.22242750e-13 * tc[4];
    } else {
        /*species 1: O2 */
        species[1] =
            -1.23393000e+03 * invT
            -4.91588000e-01
            -3.69757800e+00 * tc[0]
            -3.06759850e-04 * tc[1]
            +2.09807000e-08 * tc[2]
            -1.47940083e-12 * tc[3]
            +5.68217500e-17 * tc[4];
        /*species 2: CO2 */
        species[2] =
            -4.89669600e+04 * invT
            +4.40901890e+00
            -4.45362300e+00 * tc[0]
            -1.57008450e-03 * tc[1]
            +2.13068500e-07 * tc[2]
            -1.99499750e-11 * tc[3]
            +8.34516500e-16 * tc[4];
        /*species 3: H2O */
        species[3] =
            -2.98992100e+04 * invT
            -5.19067100e+00
            -2.67214600e+00 * tc[0]
            -1.52814650e-03 * tc[1]
            +1.45504333e-07 * tc[2]
            -1.00083000e-11 * tc[3]
            +3.19580900e-16 * tc[4];
        /*species 4: CO */
        species[4] =
            -1.42683500e+04 * invT
            -4.08314000e+00
            -3.02507800e+00 * tc[0]
            -7.21344500e-04 * tc[1]
            +9.38471333e-08 * tc[2]
            -8.48817500e-12 * tc[3]
            +3.45547600e-16 * tc[4];
        /*species 5: H2 */
        species[5] =
            -8.35034000e+02 * invT
            +3.34653300e+00
            -2.99142300e+00 * tc[0]
            -3.50032200e-04 * tc[1]
            +9.38971500e-09 * tc[2]
            +7.69298167e-13 * tc[3]
            -7.91376000e-17 * tc[4];
        /*species 6: OH */
        species[6] =
            +3.88688800e+03 * invT
            -3.71298200e+00
            -2.88273000e+00 * tc[0]
            -5.06987000e-04 * tc[1]
            +3.79479500e-08 * tc[2]
            -1.81223667e-12 * tc[3]
            +2.56315250e-17 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.80069600e+04 * invT
            +3.07203000e+00
            -4.57316700e+00 * tc[0]
            -2.16806800e-03 * tc[1]
            +2.45781500e-07 * tc[2]
            -1.95742000e-11 * tc[3]
            +7.15827000e-16 * tc[4];
        /*species 8: HO2 */
        species[8] =
            -1.57972700e+02 * invT
            -4.03838000e-01
            -4.07219100e+00 * tc[0]
            -1.06564800e-03 * tc[1]
            +8.84690833e-08 * tc[2]
            -5.09355750e-12 * tc[3]
            +1.42058250e-16 * tc[4];
        /*species 9: H */
        species[9] =
            +2.54716300e+04 * invT
            +1.96011760e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 10: O */
        species[10] =
            +2.92308000e+04 * invT
            -3.37824800e+00
            -2.54206000e+00 * tc[0]
            +1.37753100e-05 * tc[1]
            +5.17133833e-10 * tc[2]
            -3.79255583e-13 * tc[3]
            +2.18402600e-17 * tc[4];
        /*species 11: CH3O */
        species[11] =
            +1.27832500e+02 * invT
            -1.58775000e-01
            -3.77080000e+00 * tc[0]
            -3.93574850e-03 * tc[1]
            +4.42730667e-07 * tc[2]
            -3.28702583e-11 * tc[3]
            +1.05630800e-15 * tc[4];
        /*species 12: CH2O */
        species[12] =
            -1.53203700e+04 * invT
            -4.91696600e+00
            -2.99560600e+00 * tc[0]
            -3.34066050e-03 * tc[1]
            +4.38159167e-07 * tc[2]
            -3.94762750e-11 * tc[3]
            +1.60625850e-15 * tc[4];
        /*species 13: HCO */
        species[13] =
            +3.91632400e+03 * invT
            -2.99502800e+00
            -3.55727100e+00 * tc[0]
            -1.67278650e-03 * tc[1]
            +2.22501000e-07 * tc[2]
            -2.05881083e-11 * tc[3]
            +8.56925500e-16 * tc[4];
        /*species 14: CH2 */
        species[14] =
            +4.53413400e+04 * invT
            +4.79847000e-01
            -3.63640800e+00 * tc[0]
            -9.66528500e-04 * tc[1]
            +2.81169333e-08 * tc[2]
            +8.41582500e-12 * tc[3]
            -9.04128000e-16 * tc[4];
        /*species 15: CH3 */
        species[15] =
            +1.64378100e+04 * invT
            -3.60864500e+00
            -2.84405200e+00 * tc[0]
            -3.06898700e-03 * tc[1]
            +3.71724167e-07 * tc[2]
            -3.15430083e-11 * tc[3]
            +1.22607950e-15 * tc[4];
        /*species 16: CH4 */
        species[16] =
            -1.00807900e+04 * invT
            -8.93991600e+00
            -1.68347900e+00 * tc[0]
            -5.11862000e-03 * tc[1]
            +6.45854833e-07 * tc[2]
            -5.65465417e-11 * tc[3]
            +2.25171150e-15 * tc[4];
        /*species 17: C2H3 */
        species[17] =
            +3.18543500e+04 * invT
            +1.34637810e+01
            -5.93346800e+00 * tc[0]
            -2.00887300e-03 * tc[1]
            +6.61123333e-08 * tc[2]
            +1.20105583e-11 * tc[3]
            -1.18932200e-15 * tc[4];
        /*species 18: C2H4 */
        species[18] =
            +4.42828900e+03 * invT
            +2.98030000e-01
            -3.52841900e+00 * tc[0]
            -5.74259000e-03 * tc[1]
            +7.36397500e-07 * tc[2]
            -6.53716750e-11 * tc[3]
            +2.63342400e-15 * tc[4];
        /*species 19: C2H5 */
        species[19] =
            +1.06745500e+04 * invT
            +2.09713700e+01
            -7.19048000e+00 * tc[0]
            -3.24203850e-03 * tc[1]
            +1.07134417e-07 * tc[2]
            +1.95656583e-11 * tc[3]
            -1.94043850e-15 * tc[4];
        /*species 20: C3H4 */
        species[20] =
            +2.01174617e+04 * invT
            +1.62888349e+01
            -6.31694869e+00 * tc[0]
            -5.56681310e-03 * tc[1]
            +6.60481697e-07 * tc[2]
            -5.29694812e-11 * tc[3]
            +1.89374943e-15 * tc[4];
        /*species 21: C3H5 */
        species[21] =
            +1.72714707e+04 * invT
            +1.48224797e+01
            -6.54761132e+00 * tc[0]
            -6.65761230e-03 * tc[1]
            +7.97221833e-07 * tc[2]
            -6.43291512e-11 * tc[3]
            +2.30965404e-15 * tc[4];
        /*species 23: C3H7 */
        species[23] =
            +8.29843360e+03 * invT
            +2.21828787e+01
            -7.70269870e+00 * tc[0]
            -8.02210150e-03 * tc[1]
            +8.80553667e-07 * tc[2]
            -6.35821583e-11 * tc[3]
            +1.96961420e-15 * tc[4];
        /*species 28: N2 */
        species[28] =
            -9.22797700e+02 * invT
            -4.05388800e+00
            -2.92664000e+00 * tc[0]
            -7.43988500e-04 * tc[1]
            +9.47460167e-08 * tc[2]
            -8.41420000e-12 * tc[3]
            +3.37667550e-16 * tc[4];
    }

    /*species with midpoint at T=1388 kelvin */
    if (T < 1388) {
        /*species 22: C3H6 */
        species[22] =
            +1.06688164e+03 * invT
            -2.25057582e+01
            -3.94615444e-01 * tc[0]
            -1.44553831e-02 * tc[1]
            +2.58144680e-06 * tc[2]
            -3.24011841e-10 * tc[3]
            +1.68945176e-14 * tc[4];
    } else {
        /*species 22: C3H6 */
        species[22] =
            -1.87821271e+03 * invT
            +2.70320264e+01
            -8.01595958e+00 * tc[0]
            -6.85118170e-03 * tc[1]
            +7.77082888e-07 * tc[2]
            -6.01045335e-11 * tc[3]
            +2.08685063e-15 * tc[4];
    }

    /*species with midpoint at T=1390 kelvin */
    if (T < 1390) {
        /*species 25: C7H15O2 */
        species[25] =
            -1.99237961e+04 * invT
            -2.39317409e+01
            -2.37499334e+00 * tc[0]
            -4.17325953e-02 * tc[1]
            +8.56495533e-06 * tc[2]
            -1.36848052e-09 * tc[3]
            +1.09752608e-13 * tc[4];
    } else {
        /*species 25: C7H15O2 */
        species[25] =
            -2.82976050e+04 * invT
            +1.21294723e+02
            -2.49023689e+01 * tc[0]
            -1.75358460e-02 * tc[1]
            +2.00733843e-06 * tc[2]
            -1.56220685e-10 * tc[3]
            +5.44738955e-15 * tc[4];
    }

    /*species with midpoint at T=1391 kelvin */
    if (T < 1391) {
        /*species 0: NC7H16 */
        species[0] =
            -2.56586565e+04 * invT
            -3.76416531e+01
            +1.26836187e+00 * tc[0]
            -4.27177910e-02 * tc[1]
            +8.75577977e-06 * tc[2]
            -1.35788101e-09 * tc[3]
            +1.01197462e-13 * tc[4];
    } else {
        /*species 0: NC7H16 */
        species[0] =
            -3.42760081e+04 * invT
            +1.13518917e+02
            -2.22148969e+01 * tc[0]
            -1.73837875e-02 * tc[1]
            +1.97345215e-06 * tc[2]
            -1.52748732e-10 * tc[3]
            +5.30651330e-15 * tc[4];
    }

    /*species with midpoint at T=1396 kelvin */
    if (T < 1396) {
        /*species 26: C7KET12 */
        species[26] =
            -4.68054419e+04 * invT
            -3.37507112e+01
            -5.82433697e-01 * tc[0]
            -5.06039345e-02 * tc[1]
            +1.27642666e-05 * tc[2]
            -2.50615505e-09 * tc[3]
            +2.41451396e-13 * tc[4];
    } else {
        /*species 26: C7KET12 */
        species[26] =
            -5.66856828e+04 * invT
            +1.51179781e+02
            -2.97472906e+01 * tc[0]
            -1.53311147e-02 * tc[1]
            +1.75939317e-06 * tc[2]
            -1.37189452e-10 * tc[3]
            +4.79085838e-15 * tc[4];
    }
    return;
}

/*compute Cv/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void vcv_R(int npt, double * restrict species, double * restrict T_in)
{
    for (int i=0; i<npt; i++) {
        /*temperature */
	double tc[5];
        tc[0] = 0; //tc_in[0*npt+i];
        tc[1] = T_in[i]; //tc_in[1*npt+i];
        tc[2] = T_in[i]*T_in[i]; //tc_in[2*npt+i];
        tc[3] = T_in[i]*T_in[i]*T_in[i]; //tc_in[3*npt+i];
        tc[4] = T_in[i]*T_in[i]*T_in[i]*T_in[i]; //tc_in[4*npt+i];
	double T = T_in[i];

        /*species with midpoint at T=1382 kelvin */
        if (T < 1382) {
            /*species 24: C7H15-2 */
            species[24*npt+i] =
                -1.03791558e+00
                +7.56726570e-02 * tc[1]
                -4.07473634e-05 * tc[2]
                +9.32678943e-09 * tc[3]
                -4.92360745e-13 * tc[4];
        } else {
            /*species 24: C7H15-2 */
            species[24*npt+i] =
                +2.06368842e+01
                +3.23324804e-02 * tc[1]
                -1.09273807e-05 * tc[2]
                +1.68357060e-09 * tc[3]
                -9.71774091e-14 * tc[4];
        }

        /*species with midpoint at T=1383 kelvin */
        if (T < 1383) {
            /*species 27: C5H11CO */
            species[27*npt+i] =
                +1.14479069e+00
                +6.17863563e-02 * tc[1]
                -3.74134690e-05 * tc[2]
                +1.13283795e-08 * tc[3]
                -1.36917698e-12 * tc[4];
        } else {
            /*species 27: C5H11CO */
            species[27*npt+i] =
                +1.84783812e+01
                +2.50466029e-02 * tc[1]
                -8.54861346e-06 * tc[2]
                +1.32557944e-09 * tc[3]
                -7.68503296e-14 * tc[4];
        }

        /*species with midpoint at T=1000 kelvin */
        if (T < 1000) {
            /*species 1: O2 */
            species[1*npt+i] =
                +2.21293600e+00
                +1.12748600e-03 * tc[1]
                -5.75615000e-07 * tc[2]
                +1.31387700e-09 * tc[3]
                -8.76855400e-13 * tc[4];
            /*species 2: CO2 */
            species[2*npt+i] =
                +1.27572500e+00
                +9.92207200e-03 * tc[1]
                -1.04091100e-05 * tc[2]
                +6.86668700e-09 * tc[3]
                -2.11728000e-12 * tc[4];
            /*species 3: H2O */
            species[3*npt+i] =
                +2.38684200e+00
                +3.47498200e-03 * tc[1]
                -6.35469600e-06 * tc[2]
                +6.96858100e-09 * tc[3]
                -2.50658800e-12 * tc[4];
            /*species 4: CO */
            species[4*npt+i] =
                +2.26245200e+00
                +1.51194100e-03 * tc[1]
                -3.88175500e-06 * tc[2]
                +5.58194400e-09 * tc[3]
                -2.47495100e-12 * tc[4];
            /*species 5: H2 */
            species[5*npt+i] =
                +2.29812400e+00
                +8.24944200e-04 * tc[1]
                -8.14301500e-07 * tc[2]
                -9.47543400e-11 * tc[3]
                +4.13487200e-13 * tc[4];
            /*species 6: OH */
            species[6*npt+i] =
                +2.63726600e+00
                +1.85091000e-04 * tc[1]
                -1.67616500e-06 * tc[2]
                +2.38720300e-09 * tc[3]
                -8.43144200e-13 * tc[4];
            /*species 7: H2O2 */
            species[7*npt+i] =
                +2.38875400e+00
                +6.56922600e-03 * tc[1]
                -1.48501300e-07 * tc[2]
                -4.62580600e-09 * tc[3]
                +2.47151500e-12 * tc[4];
            /*species 8: HO2 */
            species[8*npt+i] =
                +1.97996300e+00
                +4.99669700e-03 * tc[1]
                -3.79099700e-06 * tc[2]
                +2.35419200e-09 * tc[3]
                -8.08902400e-13 * tc[4];
            /*species 9: H */
            species[9*npt+i] =
                +1.50000000e+00
                +0.00000000e+00 * tc[1]
                +0.00000000e+00 * tc[2]
                +0.00000000e+00 * tc[3]
                +0.00000000e+00 * tc[4];
            /*species 10: O */
            species[10*npt+i] =
                +1.94642900e+00
                -1.63816600e-03 * tc[1]
                +2.42103200e-06 * tc[2]
                -1.60284300e-09 * tc[3]
                +3.89069600e-13 * tc[4];
            /*species 11: CH3O */
            species[11*npt+i] =
                +1.10620400e+00
                +7.21659500e-03 * tc[1]
                +5.33847200e-06 * tc[2]
                -7.37763600e-09 * tc[3]
                +2.07561100e-12 * tc[4];
            /*species 12: CH2O */
            species[12*npt+i] =
                +6.52731000e-01
                +1.26314400e-02 * tc[1]
                -1.88816800e-05 * tc[2]
                +2.05003100e-08 * tc[3]
                -8.41323700e-12 * tc[4];
            /*species 13: HCO */
            species[13*npt+i] =
                +1.89833000e+00
                +6.19914700e-03 * tc[1]
                -9.62308400e-06 * tc[2]
                +1.08982500e-08 * tc[3]
                -4.57488500e-12 * tc[4];
            /*species 14: CH2 */
            species[14*npt+i] =
                +2.76223700e+00
                +1.15981900e-03 * tc[1]
                +2.48958500e-07 * tc[2]
                +8.80083600e-10 * tc[3]
                -7.33243500e-13 * tc[4];
            /*species 15: CH3 */
            species[15*npt+i] =
                +1.43044300e+00
                +1.11241000e-02 * tc[1]
                -1.68022000e-05 * tc[2]
                +1.62182900e-08 * tc[3]
                -5.86495300e-12 * tc[4];
            /*species 16: CH4 */
            species[16*npt+i] =
                -2.21258500e-01
                +1.74766800e-02 * tc[1]
                -2.78340900e-05 * tc[2]
                +3.04970800e-08 * tc[3]
                -1.22393100e-11 * tc[4];
            /*species 17: C2H3 */
            species[17*npt+i] =
                +1.45927600e+00
                +7.37147600e-03 * tc[1]
                +2.10987300e-06 * tc[2]
                -1.32164200e-09 * tc[3]
                -1.18478400e-12 * tc[4];
            /*species 18: C2H4 */
            species[18*npt+i] =
                -1.86148800e+00
                +2.79616300e-02 * tc[1]
                -3.38867700e-05 * tc[2]
                +2.78515200e-08 * tc[3]
                -9.73787900e-12 * tc[4];
            /*species 19: C2H5 */
            species[19*npt+i] =
                +1.69070200e+00
                +8.71913300e-03 * tc[1]
                +4.41983900e-06 * tc[2]
                +9.33870300e-10 * tc[3]
                -3.92777300e-12 * tc[4];
            /*species 20: C3H4 */
            species[20*npt+i] =
                +1.61307487e+00
                +1.21223371e-02 * tc[1]
                +1.85405400e-05 * tc[2]
                -3.45258475e-08 * tc[3]
                +1.53353389e-11 * tc[4];
            /*species 21: C3H5 */
            species[21*npt+i] =
                +2.78794693e+00
                +9.48414335e-03 * tc[1]
                +2.42343368e-05 * tc[2]
                -3.65604010e-08 * tc[3]
                +1.48592356e-11 * tc[4];
            /*species 23: C3H7 */
            species[23*npt+i] =
                +5.15518000e-02
                +2.59919800e-02 * tc[1]
                +2.38005400e-06 * tc[2]
                -1.96095690e-08 * tc[3]
                +9.37324700e-12 * tc[4];
            /*species 28: N2 */
            species[28*npt+i] =
                +2.29867700e+00
                +1.40824000e-03 * tc[1]
                -3.96322200e-06 * tc[2]
                +5.64151500e-09 * tc[3]
                -2.44485500e-12 * tc[4];
        } else {
            /*species 1: O2 */
            species[1*npt+i] =
                +2.69757800e+00
                +6.13519700e-04 * tc[1]
                -1.25884200e-07 * tc[2]
                +1.77528100e-11 * tc[3]
                -1.13643500e-15 * tc[4];
            /*species 2: CO2 */
            species[2*npt+i] =
                +3.45362300e+00
                +3.14016900e-03 * tc[1]
                -1.27841100e-06 * tc[2]
                +2.39399700e-10 * tc[3]
                -1.66903300e-14 * tc[4];
            /*species 3: H2O */
            species[3*npt+i] =
                +1.67214600e+00
                +3.05629300e-03 * tc[1]
                -8.73026000e-07 * tc[2]
                +1.20099600e-10 * tc[3]
                -6.39161800e-15 * tc[4];
            /*species 4: CO */
            species[4*npt+i] =
                +2.02507800e+00
                +1.44268900e-03 * tc[1]
                -5.63082800e-07 * tc[2]
                +1.01858100e-10 * tc[3]
                -6.91095200e-15 * tc[4];
            /*species 5: H2 */
            species[5*npt+i] =
                +1.99142300e+00
                +7.00064400e-04 * tc[1]
                -5.63382900e-08 * tc[2]
                -9.23157800e-12 * tc[3]
                +1.58275200e-15 * tc[4];
            /*species 6: OH */
            species[6*npt+i] =
                +1.88273000e+00
                +1.01397400e-03 * tc[1]
                -2.27687700e-07 * tc[2]
                +2.17468400e-11 * tc[3]
                -5.12630500e-16 * tc[4];
            /*species 7: H2O2 */
            species[7*npt+i] =
                +3.57316700e+00
                +4.33613600e-03 * tc[1]
                -1.47468900e-06 * tc[2]
                +2.34890400e-10 * tc[3]
                -1.43165400e-14 * tc[4];
            /*species 8: HO2 */
            species[8*npt+i] =
                +3.07219100e+00
                +2.13129600e-03 * tc[1]
                -5.30814500e-07 * tc[2]
                +6.11226900e-11 * tc[3]
                -2.84116500e-15 * tc[4];
            /*species 9: H */
            species[9*npt+i] =
                +1.50000000e+00
                +0.00000000e+00 * tc[1]
                +0.00000000e+00 * tc[2]
                +0.00000000e+00 * tc[3]
                +0.00000000e+00 * tc[4];
            /*species 10: O */
            species[10*npt+i] =
                +1.54206000e+00
                -2.75506200e-05 * tc[1]
                -3.10280300e-09 * tc[2]
                +4.55106700e-12 * tc[3]
                -4.36805200e-16 * tc[4];
            /*species 11: CH3O */
            species[11*npt+i] =
                +2.77080000e+00
                +7.87149700e-03 * tc[1]
                -2.65638400e-06 * tc[2]
                +3.94443100e-10 * tc[3]
                -2.11261600e-14 * tc[4];
            /*species 12: CH2O */
            species[12*npt+i] =
                +1.99560600e+00
                +6.68132100e-03 * tc[1]
                -2.62895500e-06 * tc[2]
                +4.73715300e-10 * tc[3]
                -3.21251700e-14 * tc[4];
            /*species 13: HCO */
            species[13*npt+i] =
                +2.55727100e+00
                +3.34557300e-03 * tc[1]
                -1.33500600e-06 * tc[2]
                +2.47057300e-10 * tc[3]
                -1.71385100e-14 * tc[4];
            /*species 14: CH2 */
            species[14*npt+i] =
                +2.63640800e+00
                +1.93305700e-03 * tc[1]
                -1.68701600e-07 * tc[2]
                -1.00989900e-10 * tc[3]
                +1.80825600e-14 * tc[4];
            /*species 15: CH3 */
            species[15*npt+i] =
                +1.84405200e+00
                +6.13797400e-03 * tc[1]
                -2.23034500e-06 * tc[2]
                +3.78516100e-10 * tc[3]
                -2.45215900e-14 * tc[4];
            /*species 16: CH4 */
            species[16*npt+i] =
                +6.83479000e-01
                +1.02372400e-02 * tc[1]
                -3.87512900e-06 * tc[2]
                +6.78558500e-10 * tc[3]
                -4.50342300e-14 * tc[4];
            /*species 17: C2H3 */
            species[17*npt+i] =
                +4.93346800e+00
                +4.01774600e-03 * tc[1]
                -3.96674000e-07 * tc[2]
                -1.44126700e-10 * tc[3]
                +2.37864400e-14 * tc[4];
            /*species 18: C2H4 */
            species[18*npt+i] =
                +2.52841900e+00
                +1.14851800e-02 * tc[1]
                -4.41838500e-06 * tc[2]
                +7.84460100e-10 * tc[3]
                -5.26684800e-14 * tc[4];
            /*species 19: C2H5 */
            species[19*npt+i] =
                +6.19048000e+00
                +6.48407700e-03 * tc[1]
                -6.42806500e-07 * tc[2]
                -2.34787900e-10 * tc[3]
                +3.88087700e-14 * tc[4];
            /*species 20: C3H4 */
            species[20*npt+i] =
                +5.31694869e+00
                +1.11336262e-02 * tc[1]
                -3.96289018e-06 * tc[2]
                +6.35633775e-10 * tc[3]
                -3.78749885e-14 * tc[4];
            /*species 21: C3H5 */
            species[21*npt+i] =
                +5.54761132e+00
                +1.33152246e-02 * tc[1]
                -4.78333100e-06 * tc[2]
                +7.71949814e-10 * tc[3]
                -4.61930808e-14 * tc[4];
            /*species 23: C3H7 */
            species[23*npt+i] =
                +6.70269870e+00
                +1.60442030e-02 * tc[1]
                -5.28332200e-06 * tc[2]
                +7.62985900e-10 * tc[3]
                -3.93922840e-14 * tc[4];
            /*species 28: N2 */
            species[28*npt+i] =
                +1.92664000e+00
                +1.48797700e-03 * tc[1]
                -5.68476100e-07 * tc[2]
                +1.00970400e-10 * tc[3]
                -6.75335100e-15 * tc[4];
        }

        /*species with midpoint at T=1388 kelvin */
        if (T < 1388) {
            /*species 22: C3H6 */
            species[22*npt+i] =
                -6.05384556e-01
                +2.89107662e-02 * tc[1]
                -1.54886808e-05 * tc[2]
                +3.88814209e-09 * tc[3]
                -3.37890352e-13 * tc[4];
        } else {
            /*species 22: C3H6 */
            species[22*npt+i] =
                +7.01595958e+00
                +1.37023634e-02 * tc[1]
                -4.66249733e-06 * tc[2]
                +7.21254402e-10 * tc[3]
                -4.17370126e-14 * tc[4];
        }

        /*species with midpoint at T=1390 kelvin */
        if (T < 1390) {
            /*species 25: C7H15O2 */
            species[25*npt+i] =
                +1.37499334e+00
                +8.34651906e-02 * tc[1]
                -5.13897320e-05 * tc[2]
                +1.64217662e-08 * tc[3]
                -2.19505216e-12 * tc[4];
        } else {
            /*species 25: C7H15O2 */
            species[25*npt+i] =
                +2.39023689e+01
                +3.50716920e-02 * tc[1]
                -1.20440306e-05 * tc[2]
                +1.87464822e-09 * tc[3]
                -1.08947791e-13 * tc[4];
        }

        /*species with midpoint at T=1391 kelvin */
        if (T < 1391) {
            /*species 0: NC7H16 */
            species[0*npt+i] =
                -2.26836187e+00
                +8.54355820e-02 * tc[1]
                -5.25346786e-05 * tc[2]
                +1.62945721e-08 * tc[3]
                -2.02394925e-12 * tc[4];
        } else {
            /*species 0: NC7H16 */
            species[0*npt+i] =
                +2.12148969e+01
                +3.47675750e-02 * tc[1]
                -1.18407129e-05 * tc[2]
                +1.83298478e-09 * tc[3]
                -1.06130266e-13 * tc[4];
        }

        /*species with midpoint at T=1396 kelvin */
        if (T < 1396) {
            /*species 26: C7KET12 */
            species[26*npt+i] =
                -4.17566303e-01
                +1.01207869e-01 * tc[1]
                -7.65855996e-05 * tc[2]
                +3.00738606e-08 * tc[3]
                -4.82902792e-12 * tc[4];
        } else {
            /*species 26: C7KET12 */
            species[26*npt+i] =
                +2.87472906e+01
                +3.06622294e-02 * tc[1]
                -1.05563590e-05 * tc[2]
                +1.64627343e-09 * tc[3]
                -9.58171675e-14 * tc[4];
        }
    }
    return;
}

/*compute Cv/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void cv_R(double * restrict species, double * restrict tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1382 kelvin */
    if (T < 1382) {
        /*species 24: C7H15-2 */
        species[24] =
            -1.03791558e+00
            +7.56726570e-02 * tc[1]
            -4.07473634e-05 * tc[2]
            +9.32678943e-09 * tc[3]
            -4.92360745e-13 * tc[4];
    } else {
        /*species 24: C7H15-2 */
        species[24] =
            +2.06368842e+01
            +3.23324804e-02 * tc[1]
            -1.09273807e-05 * tc[2]
            +1.68357060e-09 * tc[3]
            -9.71774091e-14 * tc[4];
    }

    /*species with midpoint at T=1383 kelvin */
    if (T < 1383) {
        /*species 27: C5H11CO */
        species[27] =
            +1.14479069e+00
            +6.17863563e-02 * tc[1]
            -3.74134690e-05 * tc[2]
            +1.13283795e-08 * tc[3]
            -1.36917698e-12 * tc[4];
    } else {
        /*species 27: C5H11CO */
        species[27] =
            +1.84783812e+01
            +2.50466029e-02 * tc[1]
            -8.54861346e-06 * tc[2]
            +1.32557944e-09 * tc[3]
            -7.68503296e-14 * tc[4];
    }

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 1: O2 */
        species[1] =
            +2.21293600e+00
            +1.12748600e-03 * tc[1]
            -5.75615000e-07 * tc[2]
            +1.31387700e-09 * tc[3]
            -8.76855400e-13 * tc[4];
        /*species 2: CO2 */
        species[2] =
            +1.27572500e+00
            +9.92207200e-03 * tc[1]
            -1.04091100e-05 * tc[2]
            +6.86668700e-09 * tc[3]
            -2.11728000e-12 * tc[4];
        /*species 3: H2O */
        species[3] =
            +2.38684200e+00
            +3.47498200e-03 * tc[1]
            -6.35469600e-06 * tc[2]
            +6.96858100e-09 * tc[3]
            -2.50658800e-12 * tc[4];
        /*species 4: CO */
        species[4] =
            +2.26245200e+00
            +1.51194100e-03 * tc[1]
            -3.88175500e-06 * tc[2]
            +5.58194400e-09 * tc[3]
            -2.47495100e-12 * tc[4];
        /*species 5: H2 */
        species[5] =
            +2.29812400e+00
            +8.24944200e-04 * tc[1]
            -8.14301500e-07 * tc[2]
            -9.47543400e-11 * tc[3]
            +4.13487200e-13 * tc[4];
        /*species 6: OH */
        species[6] =
            +2.63726600e+00
            +1.85091000e-04 * tc[1]
            -1.67616500e-06 * tc[2]
            +2.38720300e-09 * tc[3]
            -8.43144200e-13 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +2.38875400e+00
            +6.56922600e-03 * tc[1]
            -1.48501300e-07 * tc[2]
            -4.62580600e-09 * tc[3]
            +2.47151500e-12 * tc[4];
        /*species 8: HO2 */
        species[8] =
            +1.97996300e+00
            +4.99669700e-03 * tc[1]
            -3.79099700e-06 * tc[2]
            +2.35419200e-09 * tc[3]
            -8.08902400e-13 * tc[4];
        /*species 9: H */
        species[9] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 10: O */
        species[10] =
            +1.94642900e+00
            -1.63816600e-03 * tc[1]
            +2.42103200e-06 * tc[2]
            -1.60284300e-09 * tc[3]
            +3.89069600e-13 * tc[4];
        /*species 11: CH3O */
        species[11] =
            +1.10620400e+00
            +7.21659500e-03 * tc[1]
            +5.33847200e-06 * tc[2]
            -7.37763600e-09 * tc[3]
            +2.07561100e-12 * tc[4];
        /*species 12: CH2O */
        species[12] =
            +6.52731000e-01
            +1.26314400e-02 * tc[1]
            -1.88816800e-05 * tc[2]
            +2.05003100e-08 * tc[3]
            -8.41323700e-12 * tc[4];
        /*species 13: HCO */
        species[13] =
            +1.89833000e+00
            +6.19914700e-03 * tc[1]
            -9.62308400e-06 * tc[2]
            +1.08982500e-08 * tc[3]
            -4.57488500e-12 * tc[4];
        /*species 14: CH2 */
        species[14] =
            +2.76223700e+00
            +1.15981900e-03 * tc[1]
            +2.48958500e-07 * tc[2]
            +8.80083600e-10 * tc[3]
            -7.33243500e-13 * tc[4];
        /*species 15: CH3 */
        species[15] =
            +1.43044300e+00
            +1.11241000e-02 * tc[1]
            -1.68022000e-05 * tc[2]
            +1.62182900e-08 * tc[3]
            -5.86495300e-12 * tc[4];
        /*species 16: CH4 */
        species[16] =
            -2.21258500e-01
            +1.74766800e-02 * tc[1]
            -2.78340900e-05 * tc[2]
            +3.04970800e-08 * tc[3]
            -1.22393100e-11 * tc[4];
        /*species 17: C2H3 */
        species[17] =
            +1.45927600e+00
            +7.37147600e-03 * tc[1]
            +2.10987300e-06 * tc[2]
            -1.32164200e-09 * tc[3]
            -1.18478400e-12 * tc[4];
        /*species 18: C2H4 */
        species[18] =
            -1.86148800e+00
            +2.79616300e-02 * tc[1]
            -3.38867700e-05 * tc[2]
            +2.78515200e-08 * tc[3]
            -9.73787900e-12 * tc[4];
        /*species 19: C2H5 */
        species[19] =
            +1.69070200e+00
            +8.71913300e-03 * tc[1]
            +4.41983900e-06 * tc[2]
            +9.33870300e-10 * tc[3]
            -3.92777300e-12 * tc[4];
        /*species 20: C3H4 */
        species[20] =
            +1.61307487e+00
            +1.21223371e-02 * tc[1]
            +1.85405400e-05 * tc[2]
            -3.45258475e-08 * tc[3]
            +1.53353389e-11 * tc[4];
        /*species 21: C3H5 */
        species[21] =
            +2.78794693e+00
            +9.48414335e-03 * tc[1]
            +2.42343368e-05 * tc[2]
            -3.65604010e-08 * tc[3]
            +1.48592356e-11 * tc[4];
        /*species 23: C3H7 */
        species[23] =
            +5.15518000e-02
            +2.59919800e-02 * tc[1]
            +2.38005400e-06 * tc[2]
            -1.96095690e-08 * tc[3]
            +9.37324700e-12 * tc[4];
        /*species 28: N2 */
        species[28] =
            +2.29867700e+00
            +1.40824000e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485500e-12 * tc[4];
    } else {
        /*species 1: O2 */
        species[1] =
            +2.69757800e+00
            +6.13519700e-04 * tc[1]
            -1.25884200e-07 * tc[2]
            +1.77528100e-11 * tc[3]
            -1.13643500e-15 * tc[4];
        /*species 2: CO2 */
        species[2] =
            +3.45362300e+00
            +3.14016900e-03 * tc[1]
            -1.27841100e-06 * tc[2]
            +2.39399700e-10 * tc[3]
            -1.66903300e-14 * tc[4];
        /*species 3: H2O */
        species[3] =
            +1.67214600e+00
            +3.05629300e-03 * tc[1]
            -8.73026000e-07 * tc[2]
            +1.20099600e-10 * tc[3]
            -6.39161800e-15 * tc[4];
        /*species 4: CO */
        species[4] =
            +2.02507800e+00
            +1.44268900e-03 * tc[1]
            -5.63082800e-07 * tc[2]
            +1.01858100e-10 * tc[3]
            -6.91095200e-15 * tc[4];
        /*species 5: H2 */
        species[5] =
            +1.99142300e+00
            +7.00064400e-04 * tc[1]
            -5.63382900e-08 * tc[2]
            -9.23157800e-12 * tc[3]
            +1.58275200e-15 * tc[4];
        /*species 6: OH */
        species[6] =
            +1.88273000e+00
            +1.01397400e-03 * tc[1]
            -2.27687700e-07 * tc[2]
            +2.17468400e-11 * tc[3]
            -5.12630500e-16 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +3.57316700e+00
            +4.33613600e-03 * tc[1]
            -1.47468900e-06 * tc[2]
            +2.34890400e-10 * tc[3]
            -1.43165400e-14 * tc[4];
        /*species 8: HO2 */
        species[8] =
            +3.07219100e+00
            +2.13129600e-03 * tc[1]
            -5.30814500e-07 * tc[2]
            +6.11226900e-11 * tc[3]
            -2.84116500e-15 * tc[4];
        /*species 9: H */
        species[9] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 10: O */
        species[10] =
            +1.54206000e+00
            -2.75506200e-05 * tc[1]
            -3.10280300e-09 * tc[2]
            +4.55106700e-12 * tc[3]
            -4.36805200e-16 * tc[4];
        /*species 11: CH3O */
        species[11] =
            +2.77080000e+00
            +7.87149700e-03 * tc[1]
            -2.65638400e-06 * tc[2]
            +3.94443100e-10 * tc[3]
            -2.11261600e-14 * tc[4];
        /*species 12: CH2O */
        species[12] =
            +1.99560600e+00
            +6.68132100e-03 * tc[1]
            -2.62895500e-06 * tc[2]
            +4.73715300e-10 * tc[3]
            -3.21251700e-14 * tc[4];
        /*species 13: HCO */
        species[13] =
            +2.55727100e+00
            +3.34557300e-03 * tc[1]
            -1.33500600e-06 * tc[2]
            +2.47057300e-10 * tc[3]
            -1.71385100e-14 * tc[4];
        /*species 14: CH2 */
        species[14] =
            +2.63640800e+00
            +1.93305700e-03 * tc[1]
            -1.68701600e-07 * tc[2]
            -1.00989900e-10 * tc[3]
            +1.80825600e-14 * tc[4];
        /*species 15: CH3 */
        species[15] =
            +1.84405200e+00
            +6.13797400e-03 * tc[1]
            -2.23034500e-06 * tc[2]
            +3.78516100e-10 * tc[3]
            -2.45215900e-14 * tc[4];
        /*species 16: CH4 */
        species[16] =
            +6.83479000e-01
            +1.02372400e-02 * tc[1]
            -3.87512900e-06 * tc[2]
            +6.78558500e-10 * tc[3]
            -4.50342300e-14 * tc[4];
        /*species 17: C2H3 */
        species[17] =
            +4.93346800e+00
            +4.01774600e-03 * tc[1]
            -3.96674000e-07 * tc[2]
            -1.44126700e-10 * tc[3]
            +2.37864400e-14 * tc[4];
        /*species 18: C2H4 */
        species[18] =
            +2.52841900e+00
            +1.14851800e-02 * tc[1]
            -4.41838500e-06 * tc[2]
            +7.84460100e-10 * tc[3]
            -5.26684800e-14 * tc[4];
        /*species 19: C2H5 */
        species[19] =
            +6.19048000e+00
            +6.48407700e-03 * tc[1]
            -6.42806500e-07 * tc[2]
            -2.34787900e-10 * tc[3]
            +3.88087700e-14 * tc[4];
        /*species 20: C3H4 */
        species[20] =
            +5.31694869e+00
            +1.11336262e-02 * tc[1]
            -3.96289018e-06 * tc[2]
            +6.35633775e-10 * tc[3]
            -3.78749885e-14 * tc[4];
        /*species 21: C3H5 */
        species[21] =
            +5.54761132e+00
            +1.33152246e-02 * tc[1]
            -4.78333100e-06 * tc[2]
            +7.71949814e-10 * tc[3]
            -4.61930808e-14 * tc[4];
        /*species 23: C3H7 */
        species[23] =
            +6.70269870e+00
            +1.60442030e-02 * tc[1]
            -5.28332200e-06 * tc[2]
            +7.62985900e-10 * tc[3]
            -3.93922840e-14 * tc[4];
        /*species 28: N2 */
        species[28] =
            +1.92664000e+00
            +1.48797700e-03 * tc[1]
            -5.68476100e-07 * tc[2]
            +1.00970400e-10 * tc[3]
            -6.75335100e-15 * tc[4];
    }

    /*species with midpoint at T=1388 kelvin */
    if (T < 1388) {
        /*species 22: C3H6 */
        species[22] =
            -6.05384556e-01
            +2.89107662e-02 * tc[1]
            -1.54886808e-05 * tc[2]
            +3.88814209e-09 * tc[3]
            -3.37890352e-13 * tc[4];
    } else {
        /*species 22: C3H6 */
        species[22] =
            +7.01595958e+00
            +1.37023634e-02 * tc[1]
            -4.66249733e-06 * tc[2]
            +7.21254402e-10 * tc[3]
            -4.17370126e-14 * tc[4];
    }

    /*species with midpoint at T=1390 kelvin */
    if (T < 1390) {
        /*species 25: C7H15O2 */
        species[25] =
            +1.37499334e+00
            +8.34651906e-02 * tc[1]
            -5.13897320e-05 * tc[2]
            +1.64217662e-08 * tc[3]
            -2.19505216e-12 * tc[4];
    } else {
        /*species 25: C7H15O2 */
        species[25] =
            +2.39023689e+01
            +3.50716920e-02 * tc[1]
            -1.20440306e-05 * tc[2]
            +1.87464822e-09 * tc[3]
            -1.08947791e-13 * tc[4];
    }

    /*species with midpoint at T=1391 kelvin */
    if (T < 1391) {
        /*species 0: NC7H16 */
        species[0] =
            -2.26836187e+00
            +8.54355820e-02 * tc[1]
            -5.25346786e-05 * tc[2]
            +1.62945721e-08 * tc[3]
            -2.02394925e-12 * tc[4];
    } else {
        /*species 0: NC7H16 */
        species[0] =
            +2.12148969e+01
            +3.47675750e-02 * tc[1]
            -1.18407129e-05 * tc[2]
            +1.83298478e-09 * tc[3]
            -1.06130266e-13 * tc[4];
    }

    /*species with midpoint at T=1396 kelvin */
    if (T < 1396) {
        /*species 26: C7KET12 */
        species[26] =
            -4.17566303e-01
            +1.01207869e-01 * tc[1]
            -7.65855996e-05 * tc[2]
            +3.00738606e-08 * tc[3]
            -4.82902792e-12 * tc[4];
    } else {
        /*species 26: C7KET12 */
        species[26] =
            +2.87472906e+01
            +3.06622294e-02 * tc[1]
            -1.05563590e-05 * tc[2]
            +1.64627343e-09 * tc[3]
            -9.58171675e-14 * tc[4];
    }
    return;
}


/*compute Cp/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void cp_R(double * restrict species, double * restrict tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1382 kelvin */
    if (T < 1382) {
        /*species 24: C7H15-2 */
        species[24] =
            -3.79155767e-02
            +7.56726570e-02 * tc[1]
            -4.07473634e-05 * tc[2]
            +9.32678943e-09 * tc[3]
            -4.92360745e-13 * tc[4];
    } else {
        /*species 24: C7H15-2 */
        species[24] =
            +2.16368842e+01
            +3.23324804e-02 * tc[1]
            -1.09273807e-05 * tc[2]
            +1.68357060e-09 * tc[3]
            -9.71774091e-14 * tc[4];
    }

    /*species with midpoint at T=1383 kelvin */
    if (T < 1383) {
        /*species 27: C5H11CO */
        species[27] =
            +2.14479069e+00
            +6.17863563e-02 * tc[1]
            -3.74134690e-05 * tc[2]
            +1.13283795e-08 * tc[3]
            -1.36917698e-12 * tc[4];
    } else {
        /*species 27: C5H11CO */
        species[27] =
            +1.94783812e+01
            +2.50466029e-02 * tc[1]
            -8.54861346e-06 * tc[2]
            +1.32557944e-09 * tc[3]
            -7.68503296e-14 * tc[4];
    }

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 1: O2 */
        species[1] =
            +3.21293600e+00
            +1.12748600e-03 * tc[1]
            -5.75615000e-07 * tc[2]
            +1.31387700e-09 * tc[3]
            -8.76855400e-13 * tc[4];
        /*species 2: CO2 */
        species[2] =
            +2.27572500e+00
            +9.92207200e-03 * tc[1]
            -1.04091100e-05 * tc[2]
            +6.86668700e-09 * tc[3]
            -2.11728000e-12 * tc[4];
        /*species 3: H2O */
        species[3] =
            +3.38684200e+00
            +3.47498200e-03 * tc[1]
            -6.35469600e-06 * tc[2]
            +6.96858100e-09 * tc[3]
            -2.50658800e-12 * tc[4];
        /*species 4: CO */
        species[4] =
            +3.26245200e+00
            +1.51194100e-03 * tc[1]
            -3.88175500e-06 * tc[2]
            +5.58194400e-09 * tc[3]
            -2.47495100e-12 * tc[4];
        /*species 5: H2 */
        species[5] =
            +3.29812400e+00
            +8.24944200e-04 * tc[1]
            -8.14301500e-07 * tc[2]
            -9.47543400e-11 * tc[3]
            +4.13487200e-13 * tc[4];
        /*species 6: OH */
        species[6] =
            +3.63726600e+00
            +1.85091000e-04 * tc[1]
            -1.67616500e-06 * tc[2]
            +2.38720300e-09 * tc[3]
            -8.43144200e-13 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +3.38875400e+00
            +6.56922600e-03 * tc[1]
            -1.48501300e-07 * tc[2]
            -4.62580600e-09 * tc[3]
            +2.47151500e-12 * tc[4];
        /*species 8: HO2 */
        species[8] =
            +2.97996300e+00
            +4.99669700e-03 * tc[1]
            -3.79099700e-06 * tc[2]
            +2.35419200e-09 * tc[3]
            -8.08902400e-13 * tc[4];
        /*species 9: H */
        species[9] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 10: O */
        species[10] =
            +2.94642900e+00
            -1.63816600e-03 * tc[1]
            +2.42103200e-06 * tc[2]
            -1.60284300e-09 * tc[3]
            +3.89069600e-13 * tc[4];
        /*species 11: CH3O */
        species[11] =
            +2.10620400e+00
            +7.21659500e-03 * tc[1]
            +5.33847200e-06 * tc[2]
            -7.37763600e-09 * tc[3]
            +2.07561100e-12 * tc[4];
        /*species 12: CH2O */
        species[12] =
            +1.65273100e+00
            +1.26314400e-02 * tc[1]
            -1.88816800e-05 * tc[2]
            +2.05003100e-08 * tc[3]
            -8.41323700e-12 * tc[4];
        /*species 13: HCO */
        species[13] =
            +2.89833000e+00
            +6.19914700e-03 * tc[1]
            -9.62308400e-06 * tc[2]
            +1.08982500e-08 * tc[3]
            -4.57488500e-12 * tc[4];
        /*species 14: CH2 */
        species[14] =
            +3.76223700e+00
            +1.15981900e-03 * tc[1]
            +2.48958500e-07 * tc[2]
            +8.80083600e-10 * tc[3]
            -7.33243500e-13 * tc[4];
        /*species 15: CH3 */
        species[15] =
            +2.43044300e+00
            +1.11241000e-02 * tc[1]
            -1.68022000e-05 * tc[2]
            +1.62182900e-08 * tc[3]
            -5.86495300e-12 * tc[4];
        /*species 16: CH4 */
        species[16] =
            +7.78741500e-01
            +1.74766800e-02 * tc[1]
            -2.78340900e-05 * tc[2]
            +3.04970800e-08 * tc[3]
            -1.22393100e-11 * tc[4];
        /*species 17: C2H3 */
        species[17] =
            +2.45927600e+00
            +7.37147600e-03 * tc[1]
            +2.10987300e-06 * tc[2]
            -1.32164200e-09 * tc[3]
            -1.18478400e-12 * tc[4];
        /*species 18: C2H4 */
        species[18] =
            -8.61488000e-01
            +2.79616300e-02 * tc[1]
            -3.38867700e-05 * tc[2]
            +2.78515200e-08 * tc[3]
            -9.73787900e-12 * tc[4];
        /*species 19: C2H5 */
        species[19] =
            +2.69070200e+00
            +8.71913300e-03 * tc[1]
            +4.41983900e-06 * tc[2]
            +9.33870300e-10 * tc[3]
            -3.92777300e-12 * tc[4];
        /*species 20: C3H4 */
        species[20] =
            +2.61307487e+00
            +1.21223371e-02 * tc[1]
            +1.85405400e-05 * tc[2]
            -3.45258475e-08 * tc[3]
            +1.53353389e-11 * tc[4];
        /*species 21: C3H5 */
        species[21] =
            +3.78794693e+00
            +9.48414335e-03 * tc[1]
            +2.42343368e-05 * tc[2]
            -3.65604010e-08 * tc[3]
            +1.48592356e-11 * tc[4];
        /*species 23: C3H7 */
        species[23] =
            +1.05155180e+00
            +2.59919800e-02 * tc[1]
            +2.38005400e-06 * tc[2]
            -1.96095690e-08 * tc[3]
            +9.37324700e-12 * tc[4];
        /*species 28: N2 */
        species[28] =
            +3.29867700e+00
            +1.40824000e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485500e-12 * tc[4];
    } else {
        /*species 1: O2 */
        species[1] =
            +3.69757800e+00
            +6.13519700e-04 * tc[1]
            -1.25884200e-07 * tc[2]
            +1.77528100e-11 * tc[3]
            -1.13643500e-15 * tc[4];
        /*species 2: CO2 */
        species[2] =
            +4.45362300e+00
            +3.14016900e-03 * tc[1]
            -1.27841100e-06 * tc[2]
            +2.39399700e-10 * tc[3]
            -1.66903300e-14 * tc[4];
        /*species 3: H2O */
        species[3] =
            +2.67214600e+00
            +3.05629300e-03 * tc[1]
            -8.73026000e-07 * tc[2]
            +1.20099600e-10 * tc[3]
            -6.39161800e-15 * tc[4];
        /*species 4: CO */
        species[4] =
            +3.02507800e+00
            +1.44268900e-03 * tc[1]
            -5.63082800e-07 * tc[2]
            +1.01858100e-10 * tc[3]
            -6.91095200e-15 * tc[4];
        /*species 5: H2 */
        species[5] =
            +2.99142300e+00
            +7.00064400e-04 * tc[1]
            -5.63382900e-08 * tc[2]
            -9.23157800e-12 * tc[3]
            +1.58275200e-15 * tc[4];
        /*species 6: OH */
        species[6] =
            +2.88273000e+00
            +1.01397400e-03 * tc[1]
            -2.27687700e-07 * tc[2]
            +2.17468400e-11 * tc[3]
            -5.12630500e-16 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +4.57316700e+00
            +4.33613600e-03 * tc[1]
            -1.47468900e-06 * tc[2]
            +2.34890400e-10 * tc[3]
            -1.43165400e-14 * tc[4];
        /*species 8: HO2 */
        species[8] =
            +4.07219100e+00
            +2.13129600e-03 * tc[1]
            -5.30814500e-07 * tc[2]
            +6.11226900e-11 * tc[3]
            -2.84116500e-15 * tc[4];
        /*species 9: H */
        species[9] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 10: O */
        species[10] =
            +2.54206000e+00
            -2.75506200e-05 * tc[1]
            -3.10280300e-09 * tc[2]
            +4.55106700e-12 * tc[3]
            -4.36805200e-16 * tc[4];
        /*species 11: CH3O */
        species[11] =
            +3.77080000e+00
            +7.87149700e-03 * tc[1]
            -2.65638400e-06 * tc[2]
            +3.94443100e-10 * tc[3]
            -2.11261600e-14 * tc[4];
        /*species 12: CH2O */
        species[12] =
            +2.99560600e+00
            +6.68132100e-03 * tc[1]
            -2.62895500e-06 * tc[2]
            +4.73715300e-10 * tc[3]
            -3.21251700e-14 * tc[4];
        /*species 13: HCO */
        species[13] =
            +3.55727100e+00
            +3.34557300e-03 * tc[1]
            -1.33500600e-06 * tc[2]
            +2.47057300e-10 * tc[3]
            -1.71385100e-14 * tc[4];
        /*species 14: CH2 */
        species[14] =
            +3.63640800e+00
            +1.93305700e-03 * tc[1]
            -1.68701600e-07 * tc[2]
            -1.00989900e-10 * tc[3]
            +1.80825600e-14 * tc[4];
        /*species 15: CH3 */
        species[15] =
            +2.84405200e+00
            +6.13797400e-03 * tc[1]
            -2.23034500e-06 * tc[2]
            +3.78516100e-10 * tc[3]
            -2.45215900e-14 * tc[4];
        /*species 16: CH4 */
        species[16] =
            +1.68347900e+00
            +1.02372400e-02 * tc[1]
            -3.87512900e-06 * tc[2]
            +6.78558500e-10 * tc[3]
            -4.50342300e-14 * tc[4];
        /*species 17: C2H3 */
        species[17] =
            +5.93346800e+00
            +4.01774600e-03 * tc[1]
            -3.96674000e-07 * tc[2]
            -1.44126700e-10 * tc[3]
            +2.37864400e-14 * tc[4];
        /*species 18: C2H4 */
        species[18] =
            +3.52841900e+00
            +1.14851800e-02 * tc[1]
            -4.41838500e-06 * tc[2]
            +7.84460100e-10 * tc[3]
            -5.26684800e-14 * tc[4];
        /*species 19: C2H5 */
        species[19] =
            +7.19048000e+00
            +6.48407700e-03 * tc[1]
            -6.42806500e-07 * tc[2]
            -2.34787900e-10 * tc[3]
            +3.88087700e-14 * tc[4];
        /*species 20: C3H4 */
        species[20] =
            +6.31694869e+00
            +1.11336262e-02 * tc[1]
            -3.96289018e-06 * tc[2]
            +6.35633775e-10 * tc[3]
            -3.78749885e-14 * tc[4];
        /*species 21: C3H5 */
        species[21] =
            +6.54761132e+00
            +1.33152246e-02 * tc[1]
            -4.78333100e-06 * tc[2]
            +7.71949814e-10 * tc[3]
            -4.61930808e-14 * tc[4];
        /*species 23: C3H7 */
        species[23] =
            +7.70269870e+00
            +1.60442030e-02 * tc[1]
            -5.28332200e-06 * tc[2]
            +7.62985900e-10 * tc[3]
            -3.93922840e-14 * tc[4];
        /*species 28: N2 */
        species[28] =
            +2.92664000e+00
            +1.48797700e-03 * tc[1]
            -5.68476100e-07 * tc[2]
            +1.00970400e-10 * tc[3]
            -6.75335100e-15 * tc[4];
    }

    /*species with midpoint at T=1388 kelvin */
    if (T < 1388) {
        /*species 22: C3H6 */
        species[22] =
            +3.94615444e-01
            +2.89107662e-02 * tc[1]
            -1.54886808e-05 * tc[2]
            +3.88814209e-09 * tc[3]
            -3.37890352e-13 * tc[4];
    } else {
        /*species 22: C3H6 */
        species[22] =
            +8.01595958e+00
            +1.37023634e-02 * tc[1]
            -4.66249733e-06 * tc[2]
            +7.21254402e-10 * tc[3]
            -4.17370126e-14 * tc[4];
    }

    /*species with midpoint at T=1390 kelvin */
    if (T < 1390) {
        /*species 25: C7H15O2 */
        species[25] =
            +2.37499334e+00
            +8.34651906e-02 * tc[1]
            -5.13897320e-05 * tc[2]
            +1.64217662e-08 * tc[3]
            -2.19505216e-12 * tc[4];
    } else {
        /*species 25: C7H15O2 */
        species[25] =
            +2.49023689e+01
            +3.50716920e-02 * tc[1]
            -1.20440306e-05 * tc[2]
            +1.87464822e-09 * tc[3]
            -1.08947791e-13 * tc[4];
    }

    /*species with midpoint at T=1391 kelvin */
    if (T < 1391) {
        /*species 0: NC7H16 */
        species[0] =
            -1.26836187e+00
            +8.54355820e-02 * tc[1]
            -5.25346786e-05 * tc[2]
            +1.62945721e-08 * tc[3]
            -2.02394925e-12 * tc[4];
    } else {
        /*species 0: NC7H16 */
        species[0] =
            +2.22148969e+01
            +3.47675750e-02 * tc[1]
            -1.18407129e-05 * tc[2]
            +1.83298478e-09 * tc[3]
            -1.06130266e-13 * tc[4];
    }

    /*species with midpoint at T=1396 kelvin */
    if (T < 1396) {
        /*species 26: C7KET12 */
        species[26] =
            +5.82433697e-01
            +1.01207869e-01 * tc[1]
            -7.65855996e-05 * tc[2]
            +3.00738606e-08 * tc[3]
            -4.82902792e-12 * tc[4];
    } else {
        /*species 26: C7KET12 */
        species[26] =
            +2.97472906e+01
            +3.06622294e-02 * tc[1]
            -1.05563590e-05 * tc[2]
            +1.64627343e-09 * tc[3]
            -9.58171675e-14 * tc[4];
    }
    return;
}

/*compute the e/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void vspeciesInternalEnergy(int npt, double * restrict species, double * restrict T_in)
{
    for (int i=0; i<npt; i++) {
        /*temperature */
	double tc[5];
        tc[0] = 0; //tc_in[0*npt+i];
        tc[1] = T_in[i]; //tc_in[1*npt+i];
        tc[2] = T_in[i]*T_in[i]; //tc_in[2*npt+i];
        tc[3] = T_in[i]*T_in[i]*T_in[i]; //tc_in[3*npt+i];
        tc[4] = T_in[i]*T_in[i]*T_in[i]*T_in[i]; //tc_in[4*npt+i];
	double T = T_in[i];
        double invT = 1 / T;

        //printf(" IN vspeciesInternalEnergy npt=%d i=%d T=%.12f ", npt, i, T);

        /*species with midpoint at T=1382 kelvin */
        if (T < 1382) {
            /*species 24: C7H15-2 */
            species[24*npt+i] =
                -1.03791558e+00
                +3.78363285e-02 * tc[1]
                -1.35824545e-05 * tc[2]
                +2.33169736e-09 * tc[3]
                -9.84721490e-14 * tc[4]
                -2.35605303e+03 * invT;
        } else {
            /*species 24: C7H15-2 */
            species[24*npt+i] =
                +2.06368842e+01
                +1.61662402e-02 * tc[1]
                -3.64246023e-06 * tc[2]
                +4.20892650e-10 * tc[3]
                -1.94354818e-14 * tc[4]
                -1.05873616e+04 * invT;
        }

        /*species with midpoint at T=1383 kelvin */
        if (T < 1383) {
            /*species 27: C5H11CO */
            species[27*npt+i] =
                +1.14479069e+00
                +3.08931781e-02 * tc[1]
                -1.24711563e-05 * tc[2]
                +2.83209488e-09 * tc[3]
                -2.73835396e-13 * tc[4]
                -1.43451172e+04 * invT;
        } else {
            /*species 27: C5H11CO */
            species[27*npt+i] =
                +1.84783812e+01
                +1.25233014e-02 * tc[1]
                -2.84953782e-06 * tc[2]
                +3.31394860e-10 * tc[3]
                -1.53700659e-14 * tc[4]
                -2.07923937e+04 * invT;
        }

        /*species with midpoint at T=1000 kelvin */
        if (T < 1000) {
            /*species 1: O2 */
            species[1*npt+i] =
                +2.21293600e+00
                +5.63743000e-04 * tc[1]
                -1.91871667e-07 * tc[2]
                +3.28469250e-10 * tc[3]
                -1.75371080e-13 * tc[4]
                -1.00524900e+03 * invT;
            /*species 2: CO2 */
            species[2*npt+i] =
                +1.27572500e+00
                +4.96103600e-03 * tc[1]
                -3.46970333e-06 * tc[2]
                +1.71667175e-09 * tc[3]
                -4.23456000e-13 * tc[4]
                -4.83731400e+04 * invT;
            /*species 3: H2O */
            species[3*npt+i] =
                +2.38684200e+00
                +1.73749100e-03 * tc[1]
                -2.11823200e-06 * tc[2]
                +1.74214525e-09 * tc[3]
                -5.01317600e-13 * tc[4]
                -3.02081100e+04 * invT;
            /*species 4: CO */
            species[4*npt+i] =
                +2.26245200e+00
                +7.55970500e-04 * tc[1]
                -1.29391833e-06 * tc[2]
                +1.39548600e-09 * tc[3]
                -4.94990200e-13 * tc[4]
                -1.43105400e+04 * invT;
            /*species 5: H2 */
            species[5*npt+i] =
                +2.29812400e+00
                +4.12472100e-04 * tc[1]
                -2.71433833e-07 * tc[2]
                -2.36885850e-11 * tc[3]
                +8.26974400e-14 * tc[4]
                -1.01252100e+03 * invT;
            /*species 6: OH */
            species[6*npt+i] =
                +2.63726600e+00
                +9.25455000e-05 * tc[1]
                -5.58721667e-07 * tc[2]
                +5.96800750e-10 * tc[3]
                -1.68628840e-13 * tc[4]
                +3.60678200e+03 * invT;
            /*species 7: H2O2 */
            species[7*npt+i] =
                +2.38875400e+00
                +3.28461300e-03 * tc[1]
                -4.95004333e-08 * tc[2]
                -1.15645150e-09 * tc[3]
                +4.94303000e-13 * tc[4]
                -1.76631500e+04 * invT;
            /*species 8: HO2 */
            species[8*npt+i] =
                +1.97996300e+00
                +2.49834850e-03 * tc[1]
                -1.26366567e-06 * tc[2]
                +5.88548000e-10 * tc[3]
                -1.61780480e-13 * tc[4]
                +1.76227400e+02 * invT;
            /*species 9: H */
            species[9*npt+i] =
                +1.50000000e+00
                +0.00000000e+00 * tc[1]
                +0.00000000e+00 * tc[2]
                +0.00000000e+00 * tc[3]
                +0.00000000e+00 * tc[4]
                +2.54716300e+04 * invT;
            /*species 10: O */
            species[10*npt+i] =
                +1.94642900e+00
                -8.19083000e-04 * tc[1]
                +8.07010667e-07 * tc[2]
                -4.00710750e-10 * tc[3]
                +7.78139200e-14 * tc[4]
                +2.91476400e+04 * invT;
            /*species 11: CH3O */
            species[11*npt+i] =
                +1.10620400e+00
                +3.60829750e-03 * tc[1]
                +1.77949067e-06 * tc[2]
                -1.84440900e-09 * tc[3]
                +4.15122200e-13 * tc[4]
                +9.78601100e+02 * invT;
            /*species 12: CH2O */
            species[12*npt+i] =
                +6.52731000e-01
                +6.31572000e-03 * tc[1]
                -6.29389333e-06 * tc[2]
                +5.12507750e-09 * tc[3]
                -1.68264740e-12 * tc[4]
                -1.48654000e+04 * invT;
            /*species 13: HCO */
            species[13*npt+i] =
                +1.89833000e+00
                +3.09957350e-03 * tc[1]
                -3.20769467e-06 * tc[2]
                +2.72456250e-09 * tc[3]
                -9.14977000e-13 * tc[4]
                +4.15992200e+03 * invT;
            /*species 14: CH2 */
            species[14*npt+i] =
                +2.76223700e+00
                +5.79909500e-04 * tc[1]
                +8.29861667e-08 * tc[2]
                +2.20020900e-10 * tc[3]
                -1.46648700e-13 * tc[4]
                +4.53679100e+04 * invT;
            /*species 15: CH3 */
            species[15*npt+i] =
                +1.43044300e+00
                +5.56205000e-03 * tc[1]
                -5.60073333e-06 * tc[2]
                +4.05457250e-09 * tc[3]
                -1.17299060e-12 * tc[4]
                +1.64237800e+04 * invT;
            /*species 16: CH4 */
            species[16*npt+i] =
                -2.21258500e-01
                +8.73834000e-03 * tc[1]
                -9.27803000e-06 * tc[2]
                +7.62427000e-09 * tc[3]
                -2.44786200e-12 * tc[4]
                -9.82522900e+03 * invT;
            /*species 17: C2H3 */
            species[17*npt+i] =
                +1.45927600e+00
                +3.68573800e-03 * tc[1]
                +7.03291000e-07 * tc[2]
                -3.30410500e-10 * tc[3]
                -2.36956800e-13 * tc[4]
                +3.33522500e+04 * invT;
            /*species 18: C2H4 */
            species[18*npt+i] =
                -1.86148800e+00
                +1.39808150e-02 * tc[1]
                -1.12955900e-05 * tc[2]
                +6.96288000e-09 * tc[3]
                -1.94757580e-12 * tc[4]
                +5.57304600e+03 * invT;
            /*species 19: C2H5 */
            species[19*npt+i] =
                +1.69070200e+00
                +4.35956650e-03 * tc[1]
                +1.47327967e-06 * tc[2]
                +2.33467575e-10 * tc[3]
                -7.85554600e-13 * tc[4]
                +1.28704000e+04 * invT;
            /*species 20: C3H4 */
            species[20*npt+i] =
                +1.61307487e+00
                +6.06116855e-03 * tc[1]
                +6.18018000e-06 * tc[2]
                -8.63146187e-09 * tc[3]
                +3.06706778e-12 * tc[4]
                +2.15415642e+04 * invT;
            /*species 21: C3H5 */
            species[21*npt+i] =
                +2.78794693e+00
                +4.74207167e-03 * tc[1]
                +8.07811227e-06 * tc[2]
                -9.14010025e-09 * tc[3]
                +2.97184712e-12 * tc[4]
                +1.86261218e+04 * invT;
            /*species 23: C3H7 */
            species[23*npt+i] =
                +5.15518000e-02
                +1.29959900e-02 * tc[1]
                +7.93351333e-07 * tc[2]
                -4.90239225e-09 * tc[3]
                +1.87464940e-12 * tc[4]
                +1.06318630e+04 * invT;
            /*species 28: N2 */
            species[28*npt+i] =
                +2.29867700e+00
                +7.04120000e-04 * tc[1]
                -1.32107400e-06 * tc[2]
                +1.41037875e-09 * tc[3]
                -4.88971000e-13 * tc[4]
                -1.02090000e+03 * invT;
        } else {
            /*species 1: O2 */
            species[1*npt+i] =
                +2.69757800e+00
                +3.06759850e-04 * tc[1]
                -4.19614000e-08 * tc[2]
                +4.43820250e-12 * tc[3]
                -2.27287000e-16 * tc[4]
                -1.23393000e+03 * invT;
            /*species 2: CO2 */
            species[2*npt+i] =
                +3.45362300e+00
                +1.57008450e-03 * tc[1]
                -4.26137000e-07 * tc[2]
                +5.98499250e-11 * tc[3]
                -3.33806600e-15 * tc[4]
                -4.89669600e+04 * invT;
            /*species 3: H2O */
            species[3*npt+i] =
                +1.67214600e+00
                +1.52814650e-03 * tc[1]
                -2.91008667e-07 * tc[2]
                +3.00249000e-11 * tc[3]
                -1.27832360e-15 * tc[4]
                -2.98992100e+04 * invT;
            /*species 4: CO */
            species[4*npt+i] =
                +2.02507800e+00
                +7.21344500e-04 * tc[1]
                -1.87694267e-07 * tc[2]
                +2.54645250e-11 * tc[3]
                -1.38219040e-15 * tc[4]
                -1.42683500e+04 * invT;
            /*species 5: H2 */
            species[5*npt+i] =
                +1.99142300e+00
                +3.50032200e-04 * tc[1]
                -1.87794300e-08 * tc[2]
                -2.30789450e-12 * tc[3]
                +3.16550400e-16 * tc[4]
                -8.35034000e+02 * invT;
            /*species 6: OH */
            species[6*npt+i] =
                +1.88273000e+00
                +5.06987000e-04 * tc[1]
                -7.58959000e-08 * tc[2]
                +5.43671000e-12 * tc[3]
                -1.02526100e-16 * tc[4]
                +3.88688800e+03 * invT;
            /*species 7: H2O2 */
            species[7*npt+i] =
                +3.57316700e+00
                +2.16806800e-03 * tc[1]
                -4.91563000e-07 * tc[2]
                +5.87226000e-11 * tc[3]
                -2.86330800e-15 * tc[4]
                -1.80069600e+04 * invT;
            /*species 8: HO2 */
            species[8*npt+i] =
                +3.07219100e+00
                +1.06564800e-03 * tc[1]
                -1.76938167e-07 * tc[2]
                +1.52806725e-11 * tc[3]
                -5.68233000e-16 * tc[4]
                -1.57972700e+02 * invT;
            /*species 9: H */
            species[9*npt+i] =
                +1.50000000e+00
                +0.00000000e+00 * tc[1]
                +0.00000000e+00 * tc[2]
                +0.00000000e+00 * tc[3]
                +0.00000000e+00 * tc[4]
                +2.54716300e+04 * invT;
            /*species 10: O */
            species[10*npt+i] =
                +1.54206000e+00
                -1.37753100e-05 * tc[1]
                -1.03426767e-09 * tc[2]
                +1.13776675e-12 * tc[3]
                -8.73610400e-17 * tc[4]
                +2.92308000e+04 * invT;
            /*species 11: CH3O */
            species[11*npt+i] =
                +2.77080000e+00
                +3.93574850e-03 * tc[1]
                -8.85461333e-07 * tc[2]
                +9.86107750e-11 * tc[3]
                -4.22523200e-15 * tc[4]
                +1.27832500e+02 * invT;
            /*species 12: CH2O */
            species[12*npt+i] =
                +1.99560600e+00
                +3.34066050e-03 * tc[1]
                -8.76318333e-07 * tc[2]
                +1.18428825e-10 * tc[3]
                -6.42503400e-15 * tc[4]
                -1.53203700e+04 * invT;
            /*species 13: HCO */
            species[13*npt+i] =
                +2.55727100e+00
                +1.67278650e-03 * tc[1]
                -4.45002000e-07 * tc[2]
                +6.17643250e-11 * tc[3]
                -3.42770200e-15 * tc[4]
                +3.91632400e+03 * invT;
            /*species 14: CH2 */
            species[14*npt+i] =
                +2.63640800e+00
                +9.66528500e-04 * tc[1]
                -5.62338667e-08 * tc[2]
                -2.52474750e-11 * tc[3]
                +3.61651200e-15 * tc[4]
                +4.53413400e+04 * invT;
            /*species 15: CH3 */
            species[15*npt+i] =
                +1.84405200e+00
                +3.06898700e-03 * tc[1]
                -7.43448333e-07 * tc[2]
                +9.46290250e-11 * tc[3]
                -4.90431800e-15 * tc[4]
                +1.64378100e+04 * invT;
            /*species 16: CH4 */
            species[16*npt+i] =
                +6.83479000e-01
                +5.11862000e-03 * tc[1]
                -1.29170967e-06 * tc[2]
                +1.69639625e-10 * tc[3]
                -9.00684600e-15 * tc[4]
                -1.00807900e+04 * invT;
            /*species 17: C2H3 */
            species[17*npt+i] =
                +4.93346800e+00
                +2.00887300e-03 * tc[1]
                -1.32224667e-07 * tc[2]
                -3.60316750e-11 * tc[3]
                +4.75728800e-15 * tc[4]
                +3.18543500e+04 * invT;
            /*species 18: C2H4 */
            species[18*npt+i] =
                +2.52841900e+00
                +5.74259000e-03 * tc[1]
                -1.47279500e-06 * tc[2]
                +1.96115025e-10 * tc[3]
                -1.05336960e-14 * tc[4]
                +4.42828900e+03 * invT;
            /*species 19: C2H5 */
            species[19*npt+i] =
                +6.19048000e+00
                +3.24203850e-03 * tc[1]
                -2.14268833e-07 * tc[2]
                -5.86969750e-11 * tc[3]
                +7.76175400e-15 * tc[4]
                +1.06745500e+04 * invT;
            /*species 20: C3H4 */
            species[20*npt+i] =
                +5.31694869e+00
                +5.56681310e-03 * tc[1]
                -1.32096339e-06 * tc[2]
                +1.58908444e-10 * tc[3]
                -7.57499770e-15 * tc[4]
                +2.01174617e+04 * invT;
            /*species 21: C3H5 */
            species[21*npt+i] =
                +5.54761132e+00
                +6.65761230e-03 * tc[1]
                -1.59444367e-06 * tc[2]
                +1.92987453e-10 * tc[3]
                -9.23861616e-15 * tc[4]
                +1.72714707e+04 * invT;
            /*species 23: C3H7 */
            species[23*npt+i] =
                +6.70269870e+00
                +8.02210150e-03 * tc[1]
                -1.76110733e-06 * tc[2]
                +1.90746475e-10 * tc[3]
                -7.87845680e-15 * tc[4]
                +8.29843360e+03 * invT;
            /*species 28: N2 */
            species[28*npt+i] =
                +1.92664000e+00
                +7.43988500e-04 * tc[1]
                -1.89492033e-07 * tc[2]
                +2.52426000e-11 * tc[3]
                -1.35067020e-15 * tc[4]
                -9.22797700e+02 * invT;
        }

        /*species with midpoint at T=1388 kelvin */
        if (T < 1388) {
            /*species 22: C3H6 */
            species[22*npt+i] =
                -6.05384556e-01
                +1.44553831e-02 * tc[1]
                -5.16289360e-06 * tc[2]
                +9.72035522e-10 * tc[3]
                -6.75780704e-14 * tc[4]
                +1.06688164e+03 * invT;
        } else {
            /*species 22: C3H6 */
            species[22*npt+i] =
                +7.01595958e+00
                +6.85118170e-03 * tc[1]
                -1.55416578e-06 * tc[2]
                +1.80313601e-10 * tc[3]
                -8.34740252e-15 * tc[4]
                -1.87821271e+03 * invT;
        }

        /*species with midpoint at T=1390 kelvin */
        if (T < 1390) {
            /*species 25: C7H15O2 */
            species[25*npt+i] =
                +1.37499334e+00
                +4.17325953e-02 * tc[1]
                -1.71299107e-05 * tc[2]
                +4.10544155e-09 * tc[3]
                -4.39010432e-13 * tc[4]
                -1.99237961e+04 * invT;
        } else {
            /*species 25: C7H15O2 */
            species[25*npt+i] =
                +2.39023689e+01
                +1.75358460e-02 * tc[1]
                -4.01467687e-06 * tc[2]
                +4.68662055e-10 * tc[3]
                -2.17895582e-14 * tc[4]
                -2.82976050e+04 * invT;
        }

        /*species with midpoint at T=1391 kelvin */
        if (T < 1391) {
            /*species 0: NC7H16 */
            species[0*npt+i] =
                -2.26836187e+00
                +4.27177910e-02 * tc[1]
                -1.75115595e-05 * tc[2]
                +4.07364302e-09 * tc[3]
                -4.04789850e-13 * tc[4]
                -2.56586565e+04 * invT;
        } else {
            /*species 0: NC7H16 */
            species[0*npt+i] =
                +2.12148969e+01
                +1.73837875e-02 * tc[1]
                -3.94690430e-06 * tc[2]
                +4.58246195e-10 * tc[3]
                -2.12260532e-14 * tc[4]
                -3.42760081e+04 * invT;
        }

        /*species with midpoint at T=1396 kelvin */
        if (T < 1396) {
            /*species 26: C7KET12 */
            species[26*npt+i] =
                -4.17566303e-01
                +5.06039345e-02 * tc[1]
                -2.55285332e-05 * tc[2]
                +7.51846515e-09 * tc[3]
                -9.65805584e-13 * tc[4]
                -4.68054419e+04 * invT;
        } else {
            /*species 26: C7KET12 */
            species[26*npt+i] =
                +2.87472906e+01
                +1.53311147e-02 * tc[1]
                -3.51878633e-06 * tc[2]
                +4.11568357e-10 * tc[3]
                -1.91634335e-14 * tc[4]
                -5.66856828e+04 * invT;
        }
        //printf(" spec[%d]=%.12f \n", 26*npt+i, species[26*npt+i]);
    }
    return;
}

/*compute the e/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void speciesInternalEnergy(double * restrict species, double * restrict tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    //printf(" IN speciesInternalEnergy T=%.12f ", T);

    /*species with midpoint at T=1382 kelvin */
    if (T < 1382) {
        /*species 24: C7H15-2 */
        species[24] =
            -1.03791558e+00
            +3.78363285e-02 * tc[1]
            -1.35824545e-05 * tc[2]
            +2.33169736e-09 * tc[3]
            -9.84721490e-14 * tc[4]
            -2.35605303e+03 * invT;
    } else {
        /*species 24: C7H15-2 */
        species[24] =
            +2.06368842e+01
            +1.61662402e-02 * tc[1]
            -3.64246023e-06 * tc[2]
            +4.20892650e-10 * tc[3]
            -1.94354818e-14 * tc[4]
            -1.05873616e+04 * invT;
    }

    /*species with midpoint at T=1383 kelvin */
    if (T < 1383) {
        /*species 27: C5H11CO */
        species[27] =
            +1.14479069e+00
            +3.08931781e-02 * tc[1]
            -1.24711563e-05 * tc[2]
            +2.83209488e-09 * tc[3]
            -2.73835396e-13 * tc[4]
            -1.43451172e+04 * invT;
    } else {
        /*species 27: C5H11CO */
        species[27] =
            +1.84783812e+01
            +1.25233014e-02 * tc[1]
            -2.84953782e-06 * tc[2]
            +3.31394860e-10 * tc[3]
            -1.53700659e-14 * tc[4]
            -2.07923937e+04 * invT;
    }

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 1: O2 */
        species[1] =
            +2.21293600e+00
            +5.63743000e-04 * tc[1]
            -1.91871667e-07 * tc[2]
            +3.28469250e-10 * tc[3]
            -1.75371080e-13 * tc[4]
            -1.00524900e+03 * invT;
        /*species 2: CO2 */
        species[2] =
            +1.27572500e+00
            +4.96103600e-03 * tc[1]
            -3.46970333e-06 * tc[2]
            +1.71667175e-09 * tc[3]
            -4.23456000e-13 * tc[4]
            -4.83731400e+04 * invT;
        /*species 3: H2O */
        species[3] =
            +2.38684200e+00
            +1.73749100e-03 * tc[1]
            -2.11823200e-06 * tc[2]
            +1.74214525e-09 * tc[3]
            -5.01317600e-13 * tc[4]
            -3.02081100e+04 * invT;
        /*species 4: CO */
        species[4] =
            +2.26245200e+00
            +7.55970500e-04 * tc[1]
            -1.29391833e-06 * tc[2]
            +1.39548600e-09 * tc[3]
            -4.94990200e-13 * tc[4]
            -1.43105400e+04 * invT;
        /*species 5: H2 */
        species[5] =
            +2.29812400e+00
            +4.12472100e-04 * tc[1]
            -2.71433833e-07 * tc[2]
            -2.36885850e-11 * tc[3]
            +8.26974400e-14 * tc[4]
            -1.01252100e+03 * invT;
        /*species 6: OH */
        species[6] =
            +2.63726600e+00
            +9.25455000e-05 * tc[1]
            -5.58721667e-07 * tc[2]
            +5.96800750e-10 * tc[3]
            -1.68628840e-13 * tc[4]
            +3.60678200e+03 * invT;
        /*species 7: H2O2 */
        species[7] =
            +2.38875400e+00
            +3.28461300e-03 * tc[1]
            -4.95004333e-08 * tc[2]
            -1.15645150e-09 * tc[3]
            +4.94303000e-13 * tc[4]
            -1.76631500e+04 * invT;
        /*species 8: HO2 */
        species[8] =
            +1.97996300e+00
            +2.49834850e-03 * tc[1]
            -1.26366567e-06 * tc[2]
            +5.88548000e-10 * tc[3]
            -1.61780480e-13 * tc[4]
            +1.76227400e+02 * invT;
        /*species 9: H */
        species[9] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716300e+04 * invT;
        /*species 10: O */
        species[10] =
            +1.94642900e+00
            -8.19083000e-04 * tc[1]
            +8.07010667e-07 * tc[2]
            -4.00710750e-10 * tc[3]
            +7.78139200e-14 * tc[4]
            +2.91476400e+04 * invT;
        /*species 11: CH3O */
        species[11] =
            +1.10620400e+00
            +3.60829750e-03 * tc[1]
            +1.77949067e-06 * tc[2]
            -1.84440900e-09 * tc[3]
            +4.15122200e-13 * tc[4]
            +9.78601100e+02 * invT;
        /*species 12: CH2O */
        species[12] =
            +6.52731000e-01
            +6.31572000e-03 * tc[1]
            -6.29389333e-06 * tc[2]
            +5.12507750e-09 * tc[3]
            -1.68264740e-12 * tc[4]
            -1.48654000e+04 * invT;
        /*species 13: HCO */
        species[13] =
            +1.89833000e+00
            +3.09957350e-03 * tc[1]
            -3.20769467e-06 * tc[2]
            +2.72456250e-09 * tc[3]
            -9.14977000e-13 * tc[4]
            +4.15992200e+03 * invT;
        /*species 14: CH2 */
        species[14] =
            +2.76223700e+00
            +5.79909500e-04 * tc[1]
            +8.29861667e-08 * tc[2]
            +2.20020900e-10 * tc[3]
            -1.46648700e-13 * tc[4]
            +4.53679100e+04 * invT;
        /*species 15: CH3 */
        species[15] =
            +1.43044300e+00
            +5.56205000e-03 * tc[1]
            -5.60073333e-06 * tc[2]
            +4.05457250e-09 * tc[3]
            -1.17299060e-12 * tc[4]
            +1.64237800e+04 * invT;
        /*species 16: CH4 */
        species[16] =
            -2.21258500e-01
            +8.73834000e-03 * tc[1]
            -9.27803000e-06 * tc[2]
            +7.62427000e-09 * tc[3]
            -2.44786200e-12 * tc[4]
            -9.82522900e+03 * invT;
        /*species 17: C2H3 */
        species[17] =
            +1.45927600e+00
            +3.68573800e-03 * tc[1]
            +7.03291000e-07 * tc[2]
            -3.30410500e-10 * tc[3]
            -2.36956800e-13 * tc[4]
            +3.33522500e+04 * invT;
        /*species 18: C2H4 */
        species[18] =
            -1.86148800e+00
            +1.39808150e-02 * tc[1]
            -1.12955900e-05 * tc[2]
            +6.96288000e-09 * tc[3]
            -1.94757580e-12 * tc[4]
            +5.57304600e+03 * invT;
        /*species 19: C2H5 */
        species[19] =
            +1.69070200e+00
            +4.35956650e-03 * tc[1]
            +1.47327967e-06 * tc[2]
            +2.33467575e-10 * tc[3]
            -7.85554600e-13 * tc[4]
            +1.28704000e+04 * invT;
        /*species 20: C3H4 */
        species[20] =
            +1.61307487e+00
            +6.06116855e-03 * tc[1]
            +6.18018000e-06 * tc[2]
            -8.63146187e-09 * tc[3]
            +3.06706778e-12 * tc[4]
            +2.15415642e+04 * invT;
        /*species 21: C3H5 */
        species[21] =
            +2.78794693e+00
            +4.74207167e-03 * tc[1]
            +8.07811227e-06 * tc[2]
            -9.14010025e-09 * tc[3]
            +2.97184712e-12 * tc[4]
            +1.86261218e+04 * invT;
        /*species 23: C3H7 */
        species[23] =
            +5.15518000e-02
            +1.29959900e-02 * tc[1]
            +7.93351333e-07 * tc[2]
            -4.90239225e-09 * tc[3]
            +1.87464940e-12 * tc[4]
            +1.06318630e+04 * invT;
        /*species 28: N2 */
        species[28] =
            +2.29867700e+00
            +7.04120000e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88971000e-13 * tc[4]
            -1.02090000e+03 * invT;
    } else {
        /*species 1: O2 */
        species[1] =
            +2.69757800e+00
            +3.06759850e-04 * tc[1]
            -4.19614000e-08 * tc[2]
            +4.43820250e-12 * tc[3]
            -2.27287000e-16 * tc[4]
            -1.23393000e+03 * invT;
        /*species 2: CO2 */
        species[2] =
            +3.45362300e+00
            +1.57008450e-03 * tc[1]
            -4.26137000e-07 * tc[2]
            +5.98499250e-11 * tc[3]
            -3.33806600e-15 * tc[4]
            -4.89669600e+04 * invT;
        /*species 3: H2O */
        species[3] =
            +1.67214600e+00
            +1.52814650e-03 * tc[1]
            -2.91008667e-07 * tc[2]
            +3.00249000e-11 * tc[3]
            -1.27832360e-15 * tc[4]
            -2.98992100e+04 * invT;
        /*species 4: CO */
        species[4] =
            +2.02507800e+00
            +7.21344500e-04 * tc[1]
            -1.87694267e-07 * tc[2]
            +2.54645250e-11 * tc[3]
            -1.38219040e-15 * tc[4]
            -1.42683500e+04 * invT;
        /*species 5: H2 */
        species[5] =
            +1.99142300e+00
            +3.50032200e-04 * tc[1]
            -1.87794300e-08 * tc[2]
            -2.30789450e-12 * tc[3]
            +3.16550400e-16 * tc[4]
            -8.35034000e+02 * invT;
        /*species 6: OH */
        species[6] =
            +1.88273000e+00
            +5.06987000e-04 * tc[1]
            -7.58959000e-08 * tc[2]
            +5.43671000e-12 * tc[3]
            -1.02526100e-16 * tc[4]
            +3.88688800e+03 * invT;
        /*species 7: H2O2 */
        species[7] =
            +3.57316700e+00
            +2.16806800e-03 * tc[1]
            -4.91563000e-07 * tc[2]
            +5.87226000e-11 * tc[3]
            -2.86330800e-15 * tc[4]
            -1.80069600e+04 * invT;
        /*species 8: HO2 */
        species[8] =
            +3.07219100e+00
            +1.06564800e-03 * tc[1]
            -1.76938167e-07 * tc[2]
            +1.52806725e-11 * tc[3]
            -5.68233000e-16 * tc[4]
            -1.57972700e+02 * invT;
        /*species 9: H */
        species[9] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716300e+04 * invT;
        /*species 10: O */
        species[10] =
            +1.54206000e+00
            -1.37753100e-05 * tc[1]
            -1.03426767e-09 * tc[2]
            +1.13776675e-12 * tc[3]
            -8.73610400e-17 * tc[4]
            +2.92308000e+04 * invT;
        /*species 11: CH3O */
        species[11] =
            +2.77080000e+00
            +3.93574850e-03 * tc[1]
            -8.85461333e-07 * tc[2]
            +9.86107750e-11 * tc[3]
            -4.22523200e-15 * tc[4]
            +1.27832500e+02 * invT;
        /*species 12: CH2O */
        species[12] =
            +1.99560600e+00
            +3.34066050e-03 * tc[1]
            -8.76318333e-07 * tc[2]
            +1.18428825e-10 * tc[3]
            -6.42503400e-15 * tc[4]
            -1.53203700e+04 * invT;
        /*species 13: HCO */
        species[13] =
            +2.55727100e+00
            +1.67278650e-03 * tc[1]
            -4.45002000e-07 * tc[2]
            +6.17643250e-11 * tc[3]
            -3.42770200e-15 * tc[4]
            +3.91632400e+03 * invT;
        /*species 14: CH2 */
        species[14] =
            +2.63640800e+00
            +9.66528500e-04 * tc[1]
            -5.62338667e-08 * tc[2]
            -2.52474750e-11 * tc[3]
            +3.61651200e-15 * tc[4]
            +4.53413400e+04 * invT;
        /*species 15: CH3 */
        species[15] =
            +1.84405200e+00
            +3.06898700e-03 * tc[1]
            -7.43448333e-07 * tc[2]
            +9.46290250e-11 * tc[3]
            -4.90431800e-15 * tc[4]
            +1.64378100e+04 * invT;
        /*species 16: CH4 */
        species[16] =
            +6.83479000e-01
            +5.11862000e-03 * tc[1]
            -1.29170967e-06 * tc[2]
            +1.69639625e-10 * tc[3]
            -9.00684600e-15 * tc[4]
            -1.00807900e+04 * invT;
        /*species 17: C2H3 */
        species[17] =
            +4.93346800e+00
            +2.00887300e-03 * tc[1]
            -1.32224667e-07 * tc[2]
            -3.60316750e-11 * tc[3]
            +4.75728800e-15 * tc[4]
            +3.18543500e+04 * invT;
        /*species 18: C2H4 */
        species[18] =
            +2.52841900e+00
            +5.74259000e-03 * tc[1]
            -1.47279500e-06 * tc[2]
            +1.96115025e-10 * tc[3]
            -1.05336960e-14 * tc[4]
            +4.42828900e+03 * invT;
        /*species 19: C2H5 */
        species[19] =
            +6.19048000e+00
            +3.24203850e-03 * tc[1]
            -2.14268833e-07 * tc[2]
            -5.86969750e-11 * tc[3]
            +7.76175400e-15 * tc[4]
            +1.06745500e+04 * invT;
        /*species 20: C3H4 */
        species[20] =
            +5.31694869e+00
            +5.56681310e-03 * tc[1]
            -1.32096339e-06 * tc[2]
            +1.58908444e-10 * tc[3]
            -7.57499770e-15 * tc[4]
            +2.01174617e+04 * invT;
        /*species 21: C3H5 */
        species[21] =
            +5.54761132e+00
            +6.65761230e-03 * tc[1]
            -1.59444367e-06 * tc[2]
            +1.92987453e-10 * tc[3]
            -9.23861616e-15 * tc[4]
            +1.72714707e+04 * invT;
        /*species 23: C3H7 */
        species[23] =
            +6.70269870e+00
            +8.02210150e-03 * tc[1]
            -1.76110733e-06 * tc[2]
            +1.90746475e-10 * tc[3]
            -7.87845680e-15 * tc[4]
            +8.29843360e+03 * invT;
        /*species 28: N2 */
        species[28] =
            +1.92664000e+00
            +7.43988500e-04 * tc[1]
            -1.89492033e-07 * tc[2]
            +2.52426000e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
    }

    /*species with midpoint at T=1388 kelvin */
    if (T < 1388) {
        /*species 22: C3H6 */
        species[22] =
            -6.05384556e-01
            +1.44553831e-02 * tc[1]
            -5.16289360e-06 * tc[2]
            +9.72035522e-10 * tc[3]
            -6.75780704e-14 * tc[4]
            +1.06688164e+03 * invT;
    } else {
        /*species 22: C3H6 */
        species[22] =
            +7.01595958e+00
            +6.85118170e-03 * tc[1]
            -1.55416578e-06 * tc[2]
            +1.80313601e-10 * tc[3]
            -8.34740252e-15 * tc[4]
            -1.87821271e+03 * invT;
    }

    /*species with midpoint at T=1390 kelvin */
    if (T < 1390) {
        /*species 25: C7H15O2 */
        species[25] =
            +1.37499334e+00
            +4.17325953e-02 * tc[1]
            -1.71299107e-05 * tc[2]
            +4.10544155e-09 * tc[3]
            -4.39010432e-13 * tc[4]
            -1.99237961e+04 * invT;
    } else {
        /*species 25: C7H15O2 */
        species[25] =
            +2.39023689e+01
            +1.75358460e-02 * tc[1]
            -4.01467687e-06 * tc[2]
            +4.68662055e-10 * tc[3]
            -2.17895582e-14 * tc[4]
            -2.82976050e+04 * invT;
    }

    /*species with midpoint at T=1391 kelvin */
    if (T < 1391) {
        /*species 0: NC7H16 */
        species[0] =
            -2.26836187e+00
            +4.27177910e-02 * tc[1]
            -1.75115595e-05 * tc[2]
            +4.07364302e-09 * tc[3]
            -4.04789850e-13 * tc[4]
            -2.56586565e+04 * invT;
    } else {
        /*species 0: NC7H16 */
        species[0] =
            +2.12148969e+01
            +1.73837875e-02 * tc[1]
            -3.94690430e-06 * tc[2]
            +4.58246195e-10 * tc[3]
            -2.12260532e-14 * tc[4]
            -3.42760081e+04 * invT;
    }

    /*species with midpoint at T=1396 kelvin */
    if (T < 1396) {
        /*species 26: C7KET12 */
        species[26] =
            -4.17566303e-01
            +5.06039345e-02 * tc[1]
            -2.55285332e-05 * tc[2]
            +7.51846515e-09 * tc[3]
            -9.65805584e-13 * tc[4]
            -4.68054419e+04 * invT;
    } else {
        /*species 26: C7KET12 */
        species[26] =
            +2.87472906e+01
            +1.53311147e-02 * tc[1]
            -3.51878633e-06 * tc[2]
            +4.11568357e-10 * tc[3]
            -1.91634335e-14 * tc[4]
            -5.66856828e+04 * invT;
    }
    //printf(" spec[%d]=%.12f \n", 26, species[26]);
    return;
}


/*compute the h/(RT) at the given temperature (Eq 20) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void speciesEnthalpy(double * restrict species, double * restrict tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=1382 kelvin */
    if (T < 1382) {
        /*species 24: C7H15-2 */
        species[24] =
            -3.79155767e-02
            +3.78363285e-02 * tc[1]
            -1.35824545e-05 * tc[2]
            +2.33169736e-09 * tc[3]
            -9.84721490e-14 * tc[4]
            -2.35605303e+03 * invT;
    } else {
        /*species 24: C7H15-2 */
        species[24] =
            +2.16368842e+01
            +1.61662402e-02 * tc[1]
            -3.64246023e-06 * tc[2]
            +4.20892650e-10 * tc[3]
            -1.94354818e-14 * tc[4]
            -1.05873616e+04 * invT;
    }

    /*species with midpoint at T=1383 kelvin */
    if (T < 1383) {
        /*species 27: C5H11CO */
        species[27] =
            +2.14479069e+00
            +3.08931781e-02 * tc[1]
            -1.24711563e-05 * tc[2]
            +2.83209488e-09 * tc[3]
            -2.73835396e-13 * tc[4]
            -1.43451172e+04 * invT;
    } else {
        /*species 27: C5H11CO */
        species[27] =
            +1.94783812e+01
            +1.25233014e-02 * tc[1]
            -2.84953782e-06 * tc[2]
            +3.31394860e-10 * tc[3]
            -1.53700659e-14 * tc[4]
            -2.07923937e+04 * invT;
    }

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 1: O2 */
        species[1] =
            +3.21293600e+00
            +5.63743000e-04 * tc[1]
            -1.91871667e-07 * tc[2]
            +3.28469250e-10 * tc[3]
            -1.75371080e-13 * tc[4]
            -1.00524900e+03 * invT;
        /*species 2: CO2 */
        species[2] =
            +2.27572500e+00
            +4.96103600e-03 * tc[1]
            -3.46970333e-06 * tc[2]
            +1.71667175e-09 * tc[3]
            -4.23456000e-13 * tc[4]
            -4.83731400e+04 * invT;
        /*species 3: H2O */
        species[3] =
            +3.38684200e+00
            +1.73749100e-03 * tc[1]
            -2.11823200e-06 * tc[2]
            +1.74214525e-09 * tc[3]
            -5.01317600e-13 * tc[4]
            -3.02081100e+04 * invT;
        /*species 4: CO */
        species[4] =
            +3.26245200e+00
            +7.55970500e-04 * tc[1]
            -1.29391833e-06 * tc[2]
            +1.39548600e-09 * tc[3]
            -4.94990200e-13 * tc[4]
            -1.43105400e+04 * invT;
        /*species 5: H2 */
        species[5] =
            +3.29812400e+00
            +4.12472100e-04 * tc[1]
            -2.71433833e-07 * tc[2]
            -2.36885850e-11 * tc[3]
            +8.26974400e-14 * tc[4]
            -1.01252100e+03 * invT;
        /*species 6: OH */
        species[6] =
            +3.63726600e+00
            +9.25455000e-05 * tc[1]
            -5.58721667e-07 * tc[2]
            +5.96800750e-10 * tc[3]
            -1.68628840e-13 * tc[4]
            +3.60678200e+03 * invT;
        /*species 7: H2O2 */
        species[7] =
            +3.38875400e+00
            +3.28461300e-03 * tc[1]
            -4.95004333e-08 * tc[2]
            -1.15645150e-09 * tc[3]
            +4.94303000e-13 * tc[4]
            -1.76631500e+04 * invT;
        /*species 8: HO2 */
        species[8] =
            +2.97996300e+00
            +2.49834850e-03 * tc[1]
            -1.26366567e-06 * tc[2]
            +5.88548000e-10 * tc[3]
            -1.61780480e-13 * tc[4]
            +1.76227400e+02 * invT;
        /*species 9: H */
        species[9] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716300e+04 * invT;
        /*species 10: O */
        species[10] =
            +2.94642900e+00
            -8.19083000e-04 * tc[1]
            +8.07010667e-07 * tc[2]
            -4.00710750e-10 * tc[3]
            +7.78139200e-14 * tc[4]
            +2.91476400e+04 * invT;
        /*species 11: CH3O */
        species[11] =
            +2.10620400e+00
            +3.60829750e-03 * tc[1]
            +1.77949067e-06 * tc[2]
            -1.84440900e-09 * tc[3]
            +4.15122200e-13 * tc[4]
            +9.78601100e+02 * invT;
        /*species 12: CH2O */
        species[12] =
            +1.65273100e+00
            +6.31572000e-03 * tc[1]
            -6.29389333e-06 * tc[2]
            +5.12507750e-09 * tc[3]
            -1.68264740e-12 * tc[4]
            -1.48654000e+04 * invT;
        /*species 13: HCO */
        species[13] =
            +2.89833000e+00
            +3.09957350e-03 * tc[1]
            -3.20769467e-06 * tc[2]
            +2.72456250e-09 * tc[3]
            -9.14977000e-13 * tc[4]
            +4.15992200e+03 * invT;
        /*species 14: CH2 */
        species[14] =
            +3.76223700e+00
            +5.79909500e-04 * tc[1]
            +8.29861667e-08 * tc[2]
            +2.20020900e-10 * tc[3]
            -1.46648700e-13 * tc[4]
            +4.53679100e+04 * invT;
        /*species 15: CH3 */
        species[15] =
            +2.43044300e+00
            +5.56205000e-03 * tc[1]
            -5.60073333e-06 * tc[2]
            +4.05457250e-09 * tc[3]
            -1.17299060e-12 * tc[4]
            +1.64237800e+04 * invT;
        /*species 16: CH4 */
        species[16] =
            +7.78741500e-01
            +8.73834000e-03 * tc[1]
            -9.27803000e-06 * tc[2]
            +7.62427000e-09 * tc[3]
            -2.44786200e-12 * tc[4]
            -9.82522900e+03 * invT;
        /*species 17: C2H3 */
        species[17] =
            +2.45927600e+00
            +3.68573800e-03 * tc[1]
            +7.03291000e-07 * tc[2]
            -3.30410500e-10 * tc[3]
            -2.36956800e-13 * tc[4]
            +3.33522500e+04 * invT;
        /*species 18: C2H4 */
        species[18] =
            -8.61488000e-01
            +1.39808150e-02 * tc[1]
            -1.12955900e-05 * tc[2]
            +6.96288000e-09 * tc[3]
            -1.94757580e-12 * tc[4]
            +5.57304600e+03 * invT;
        /*species 19: C2H5 */
        species[19] =
            +2.69070200e+00
            +4.35956650e-03 * tc[1]
            +1.47327967e-06 * tc[2]
            +2.33467575e-10 * tc[3]
            -7.85554600e-13 * tc[4]
            +1.28704000e+04 * invT;
        /*species 20: C3H4 */
        species[20] =
            +2.61307487e+00
            +6.06116855e-03 * tc[1]
            +6.18018000e-06 * tc[2]
            -8.63146187e-09 * tc[3]
            +3.06706778e-12 * tc[4]
            +2.15415642e+04 * invT;
        /*species 21: C3H5 */
        species[21] =
            +3.78794693e+00
            +4.74207167e-03 * tc[1]
            +8.07811227e-06 * tc[2]
            -9.14010025e-09 * tc[3]
            +2.97184712e-12 * tc[4]
            +1.86261218e+04 * invT;
        /*species 23: C3H7 */
        species[23] =
            +1.05155180e+00
            +1.29959900e-02 * tc[1]
            +7.93351333e-07 * tc[2]
            -4.90239225e-09 * tc[3]
            +1.87464940e-12 * tc[4]
            +1.06318630e+04 * invT;
        /*species 28: N2 */
        species[28] =
            +3.29867700e+00
            +7.04120000e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88971000e-13 * tc[4]
            -1.02090000e+03 * invT;
    } else {
        /*species 1: O2 */
        species[1] =
            +3.69757800e+00
            +3.06759850e-04 * tc[1]
            -4.19614000e-08 * tc[2]
            +4.43820250e-12 * tc[3]
            -2.27287000e-16 * tc[4]
            -1.23393000e+03 * invT;
        /*species 2: CO2 */
        species[2] =
            +4.45362300e+00
            +1.57008450e-03 * tc[1]
            -4.26137000e-07 * tc[2]
            +5.98499250e-11 * tc[3]
            -3.33806600e-15 * tc[4]
            -4.89669600e+04 * invT;
        /*species 3: H2O */
        species[3] =
            +2.67214600e+00
            +1.52814650e-03 * tc[1]
            -2.91008667e-07 * tc[2]
            +3.00249000e-11 * tc[3]
            -1.27832360e-15 * tc[4]
            -2.98992100e+04 * invT;
        /*species 4: CO */
        species[4] =
            +3.02507800e+00
            +7.21344500e-04 * tc[1]
            -1.87694267e-07 * tc[2]
            +2.54645250e-11 * tc[3]
            -1.38219040e-15 * tc[4]
            -1.42683500e+04 * invT;
        /*species 5: H2 */
        species[5] =
            +2.99142300e+00
            +3.50032200e-04 * tc[1]
            -1.87794300e-08 * tc[2]
            -2.30789450e-12 * tc[3]
            +3.16550400e-16 * tc[4]
            -8.35034000e+02 * invT;
        /*species 6: OH */
        species[6] =
            +2.88273000e+00
            +5.06987000e-04 * tc[1]
            -7.58959000e-08 * tc[2]
            +5.43671000e-12 * tc[3]
            -1.02526100e-16 * tc[4]
            +3.88688800e+03 * invT;
        /*species 7: H2O2 */
        species[7] =
            +4.57316700e+00
            +2.16806800e-03 * tc[1]
            -4.91563000e-07 * tc[2]
            +5.87226000e-11 * tc[3]
            -2.86330800e-15 * tc[4]
            -1.80069600e+04 * invT;
        /*species 8: HO2 */
        species[8] =
            +4.07219100e+00
            +1.06564800e-03 * tc[1]
            -1.76938167e-07 * tc[2]
            +1.52806725e-11 * tc[3]
            -5.68233000e-16 * tc[4]
            -1.57972700e+02 * invT;
        /*species 9: H */
        species[9] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716300e+04 * invT;
        /*species 10: O */
        species[10] =
            +2.54206000e+00
            -1.37753100e-05 * tc[1]
            -1.03426767e-09 * tc[2]
            +1.13776675e-12 * tc[3]
            -8.73610400e-17 * tc[4]
            +2.92308000e+04 * invT;
        /*species 11: CH3O */
        species[11] =
            +3.77080000e+00
            +3.93574850e-03 * tc[1]
            -8.85461333e-07 * tc[2]
            +9.86107750e-11 * tc[3]
            -4.22523200e-15 * tc[4]
            +1.27832500e+02 * invT;
        /*species 12: CH2O */
        species[12] =
            +2.99560600e+00
            +3.34066050e-03 * tc[1]
            -8.76318333e-07 * tc[2]
            +1.18428825e-10 * tc[3]
            -6.42503400e-15 * tc[4]
            -1.53203700e+04 * invT;
        /*species 13: HCO */
        species[13] =
            +3.55727100e+00
            +1.67278650e-03 * tc[1]
            -4.45002000e-07 * tc[2]
            +6.17643250e-11 * tc[3]
            -3.42770200e-15 * tc[4]
            +3.91632400e+03 * invT;
        /*species 14: CH2 */
        species[14] =
            +3.63640800e+00
            +9.66528500e-04 * tc[1]
            -5.62338667e-08 * tc[2]
            -2.52474750e-11 * tc[3]
            +3.61651200e-15 * tc[4]
            +4.53413400e+04 * invT;
        /*species 15: CH3 */
        species[15] =
            +2.84405200e+00
            +3.06898700e-03 * tc[1]
            -7.43448333e-07 * tc[2]
            +9.46290250e-11 * tc[3]
            -4.90431800e-15 * tc[4]
            +1.64378100e+04 * invT;
        /*species 16: CH4 */
        species[16] =
            +1.68347900e+00
            +5.11862000e-03 * tc[1]
            -1.29170967e-06 * tc[2]
            +1.69639625e-10 * tc[3]
            -9.00684600e-15 * tc[4]
            -1.00807900e+04 * invT;
        /*species 17: C2H3 */
        species[17] =
            +5.93346800e+00
            +2.00887300e-03 * tc[1]
            -1.32224667e-07 * tc[2]
            -3.60316750e-11 * tc[3]
            +4.75728800e-15 * tc[4]
            +3.18543500e+04 * invT;
        /*species 18: C2H4 */
        species[18] =
            +3.52841900e+00
            +5.74259000e-03 * tc[1]
            -1.47279500e-06 * tc[2]
            +1.96115025e-10 * tc[3]
            -1.05336960e-14 * tc[4]
            +4.42828900e+03 * invT;
        /*species 19: C2H5 */
        species[19] =
            +7.19048000e+00
            +3.24203850e-03 * tc[1]
            -2.14268833e-07 * tc[2]
            -5.86969750e-11 * tc[3]
            +7.76175400e-15 * tc[4]
            +1.06745500e+04 * invT;
        /*species 20: C3H4 */
        species[20] =
            +6.31694869e+00
            +5.56681310e-03 * tc[1]
            -1.32096339e-06 * tc[2]
            +1.58908444e-10 * tc[3]
            -7.57499770e-15 * tc[4]
            +2.01174617e+04 * invT;
        /*species 21: C3H5 */
        species[21] =
            +6.54761132e+00
            +6.65761230e-03 * tc[1]
            -1.59444367e-06 * tc[2]
            +1.92987453e-10 * tc[3]
            -9.23861616e-15 * tc[4]
            +1.72714707e+04 * invT;
        /*species 23: C3H7 */
        species[23] =
            +7.70269870e+00
            +8.02210150e-03 * tc[1]
            -1.76110733e-06 * tc[2]
            +1.90746475e-10 * tc[3]
            -7.87845680e-15 * tc[4]
            +8.29843360e+03 * invT;
        /*species 28: N2 */
        species[28] =
            +2.92664000e+00
            +7.43988500e-04 * tc[1]
            -1.89492033e-07 * tc[2]
            +2.52426000e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
    }

    /*species with midpoint at T=1388 kelvin */
    if (T < 1388) {
        /*species 22: C3H6 */
        species[22] =
            +3.94615444e-01
            +1.44553831e-02 * tc[1]
            -5.16289360e-06 * tc[2]
            +9.72035522e-10 * tc[3]
            -6.75780704e-14 * tc[4]
            +1.06688164e+03 * invT;
    } else {
        /*species 22: C3H6 */
        species[22] =
            +8.01595958e+00
            +6.85118170e-03 * tc[1]
            -1.55416578e-06 * tc[2]
            +1.80313601e-10 * tc[3]
            -8.34740252e-15 * tc[4]
            -1.87821271e+03 * invT;
    }

    /*species with midpoint at T=1390 kelvin */
    if (T < 1390) {
        /*species 25: C7H15O2 */
        species[25] =
            +2.37499334e+00
            +4.17325953e-02 * tc[1]
            -1.71299107e-05 * tc[2]
            +4.10544155e-09 * tc[3]
            -4.39010432e-13 * tc[4]
            -1.99237961e+04 * invT;
    } else {
        /*species 25: C7H15O2 */
        species[25] =
            +2.49023689e+01
            +1.75358460e-02 * tc[1]
            -4.01467687e-06 * tc[2]
            +4.68662055e-10 * tc[3]
            -2.17895582e-14 * tc[4]
            -2.82976050e+04 * invT;
    }

    /*species with midpoint at T=1391 kelvin */
    if (T < 1391) {
        /*species 0: NC7H16 */
        species[0] =
            -1.26836187e+00
            +4.27177910e-02 * tc[1]
            -1.75115595e-05 * tc[2]
            +4.07364302e-09 * tc[3]
            -4.04789850e-13 * tc[4]
            -2.56586565e+04 * invT;
    } else {
        /*species 0: NC7H16 */
        species[0] =
            +2.22148969e+01
            +1.73837875e-02 * tc[1]
            -3.94690430e-06 * tc[2]
            +4.58246195e-10 * tc[3]
            -2.12260532e-14 * tc[4]
            -3.42760081e+04 * invT;
    }

    /*species with midpoint at T=1396 kelvin */
    if (T < 1396) {
        /*species 26: C7KET12 */
        species[26] =
            +5.82433697e-01
            +5.06039345e-02 * tc[1]
            -2.55285332e-05 * tc[2]
            +7.51846515e-09 * tc[3]
            -9.65805584e-13 * tc[4]
            -4.68054419e+04 * invT;
    } else {
        /*species 26: C7KET12 */
        species[26] =
            +2.97472906e+01
            +1.53311147e-02 * tc[1]
            -3.51878633e-06 * tc[2]
            +4.11568357e-10 * tc[3]
            -1.91634335e-14 * tc[4]
            -5.66856828e+04 * invT;
    }
    return;
}


/*compute the S/R at the given temperature (Eq 21) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void speciesEntropy(double * restrict species, double * restrict tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1382 kelvin */
    if (T < 1382) {
        /*species 24: C7H15-2 */
        species[24] =
            -3.79155767e-02 * tc[0]
            +7.56726570e-02 * tc[1]
            -2.03736817e-05 * tc[2]
            +3.10892981e-09 * tc[3]
            -1.23090186e-13 * tc[4]
            +3.37321506e+01 ;
    } else {
        /*species 24: C7H15-2 */
        species[24] =
            +2.16368842e+01 * tc[0]
            +3.23324804e-02 * tc[1]
            -5.46369035e-06 * tc[2]
            +5.61190200e-10 * tc[3]
            -2.42943523e-14 * tc[4]
            -8.52209653e+01 ;
    }

    /*species with midpoint at T=1383 kelvin */
    if (T < 1383) {
        /*species 27: C5H11CO */
        species[27] =
            +2.14479069e+00 * tc[0]
            +6.17863563e-02 * tc[1]
            -1.87067345e-05 * tc[2]
            +3.77612650e-09 * tc[3]
            -3.42294245e-13 * tc[4]
            +2.23128045e+01 ;
    } else {
        /*species 27: C5H11CO */
        species[27] =
            +1.94783812e+01 * tc[0]
            +2.50466029e-02 * tc[1]
            -4.27430673e-06 * tc[2]
            +4.41859813e-10 * tc[3]
            -1.92125824e-14 * tc[4]
            -7.21995578e+01 ;
    }

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 1: O2 */
        species[1] =
            +3.21293600e+00 * tc[0]
            +1.12748600e-03 * tc[1]
            -2.87807500e-07 * tc[2]
            +4.37959000e-10 * tc[3]
            -2.19213850e-13 * tc[4]
            +6.03473800e+00 ;
        /*species 2: CO2 */
        species[2] =
            +2.27572500e+00 * tc[0]
            +9.92207200e-03 * tc[1]
            -5.20455500e-06 * tc[2]
            +2.28889567e-09 * tc[3]
            -5.29320000e-13 * tc[4]
            +1.01884900e+01 ;
        /*species 3: H2O */
        species[3] =
            +3.38684200e+00 * tc[0]
            +3.47498200e-03 * tc[1]
            -3.17734800e-06 * tc[2]
            +2.32286033e-09 * tc[3]
            -6.26647000e-13 * tc[4]
            +2.59023300e+00 ;
        /*species 4: CO */
        species[4] =
            +3.26245200e+00 * tc[0]
            +1.51194100e-03 * tc[1]
            -1.94087750e-06 * tc[2]
            +1.86064800e-09 * tc[3]
            -6.18737750e-13 * tc[4]
            +4.84889700e+00 ;
        /*species 5: H2 */
        species[5] =
            +3.29812400e+00 * tc[0]
            +8.24944200e-04 * tc[1]
            -4.07150750e-07 * tc[2]
            -3.15847800e-11 * tc[3]
            +1.03371800e-13 * tc[4]
            -3.29409400e+00 ;
        /*species 6: OH */
        species[6] =
            +3.63726600e+00 * tc[0]
            +1.85091000e-04 * tc[1]
            -8.38082500e-07 * tc[2]
            +7.95734333e-10 * tc[3]
            -2.10786050e-13 * tc[4]
            +1.35886000e+00 ;
        /*species 7: H2O2 */
        species[7] =
            +3.38875400e+00 * tc[0]
            +6.56922600e-03 * tc[1]
            -7.42506500e-08 * tc[2]
            -1.54193533e-09 * tc[3]
            +6.17878750e-13 * tc[4]
            +6.78536300e+00 ;
        /*species 8: HO2 */
        species[8] =
            +2.97996300e+00 * tc[0]
            +4.99669700e-03 * tc[1]
            -1.89549850e-06 * tc[2]
            +7.84730667e-10 * tc[3]
            -2.02225600e-13 * tc[4]
            +9.22272400e+00 ;
        /*species 9: H */
        species[9] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -4.60117600e-01 ;
        /*species 10: O */
        species[10] =
            +2.94642900e+00 * tc[0]
            -1.63816600e-03 * tc[1]
            +1.21051600e-06 * tc[2]
            -5.34281000e-10 * tc[3]
            +9.72674000e-14 * tc[4]
            +2.96399500e+00 ;
        /*species 11: CH3O */
        species[11] =
            +2.10620400e+00 * tc[0]
            +7.21659500e-03 * tc[1]
            +2.66923600e-06 * tc[2]
            -2.45921200e-09 * tc[3]
            +5.18902750e-13 * tc[4]
            +1.31521800e+01 ;
        /*species 12: CH2O */
        species[12] =
            +1.65273100e+00 * tc[0]
            +1.26314400e-02 * tc[1]
            -9.44084000e-06 * tc[2]
            +6.83343667e-09 * tc[3]
            -2.10330925e-12 * tc[4]
            +1.37848200e+01 ;
        /*species 13: HCO */
        species[13] =
            +2.89833000e+00 * tc[0]
            +6.19914700e-03 * tc[1]
            -4.81154200e-06 * tc[2]
            +3.63275000e-09 * tc[3]
            -1.14372125e-12 * tc[4]
            +8.98361400e+00 ;
        /*species 14: CH2 */
        species[14] =
            +3.76223700e+00 * tc[0]
            +1.15981900e-03 * tc[1]
            +1.24479250e-07 * tc[2]
            +2.93361200e-10 * tc[3]
            -1.83310875e-13 * tc[4]
            +1.71257800e+00 ;
        /*species 15: CH3 */
        species[15] =
            +2.43044300e+00 * tc[0]
            +1.11241000e-02 * tc[1]
            -8.40110000e-06 * tc[2]
            +5.40609667e-09 * tc[3]
            -1.46623825e-12 * tc[4]
            +6.78979400e+00 ;
        /*species 16: CH4 */
        species[16] =
            +7.78741500e-01 * tc[0]
            +1.74766800e-02 * tc[1]
            -1.39170450e-05 * tc[2]
            +1.01656933e-08 * tc[3]
            -3.05982750e-12 * tc[4]
            +1.37221900e+01 ;
        /*species 17: C2H3 */
        species[17] =
            +2.45927600e+00 * tc[0]
            +7.37147600e-03 * tc[1]
            +1.05493650e-06 * tc[2]
            -4.40547333e-10 * tc[3]
            -2.96196000e-13 * tc[4]
            +1.15562000e+01 ;
        /*species 18: C2H4 */
        species[18] =
            -8.61488000e-01 * tc[0]
            +2.79616300e-02 * tc[1]
            -1.69433850e-05 * tc[2]
            +9.28384000e-09 * tc[3]
            -2.43446975e-12 * tc[4]
            +2.42114900e+01 ;
        /*species 19: C2H5 */
        species[19] =
            +2.69070200e+00 * tc[0]
            +8.71913300e-03 * tc[1]
            +2.20991950e-06 * tc[2]
            +3.11290100e-10 * tc[3]
            -9.81943250e-13 * tc[4]
            +1.21382000e+01 ;
        /*species 20: C3H4 */
        species[20] =
            +2.61307487e+00 * tc[0]
            +1.21223371e-02 * tc[1]
            +9.27027000e-06 * tc[2]
            -1.15086158e-08 * tc[3]
            +3.83383472e-12 * tc[4]
            +1.02503319e+01 ;
        /*species 21: C3H5 */
        species[21] =
            +3.78794693e+00 * tc[0]
            +9.48414335e-03 * tc[1]
            +1.21171684e-05 * tc[2]
            -1.21868003e-08 * tc[3]
            +3.71480890e-12 * tc[4]
            +7.82822499e+00 ;
        /*species 23: C3H7 */
        species[23] =
            +1.05155180e+00 * tc[0]
            +2.59919800e-02 * tc[1]
            +1.19002700e-06 * tc[2]
            -6.53652300e-09 * tc[3]
            +2.34331175e-12 * tc[4]
            +2.11225590e+01 ;
        /*species 28: N2 */
        species[28] =
            +3.29867700e+00 * tc[0]
            +1.40824000e-03 * tc[1]
            -1.98161100e-06 * tc[2]
            +1.88050500e-09 * tc[3]
            -6.11213750e-13 * tc[4]
            +3.95037200e+00 ;
    } else {
        /*species 1: O2 */
        species[1] =
            +3.69757800e+00 * tc[0]
            +6.13519700e-04 * tc[1]
            -6.29421000e-08 * tc[2]
            +5.91760333e-12 * tc[3]
            -2.84108750e-16 * tc[4]
            +3.18916600e+00 ;
        /*species 2: CO2 */
        species[2] =
            +4.45362300e+00 * tc[0]
            +3.14016900e-03 * tc[1]
            -6.39205500e-07 * tc[2]
            +7.97999000e-11 * tc[3]
            -4.17258250e-15 * tc[4]
            -9.55395900e-01 ;
        /*species 3: H2O */
        species[3] =
            +2.67214600e+00 * tc[0]
            +3.05629300e-03 * tc[1]
            -4.36513000e-07 * tc[2]
            +4.00332000e-11 * tc[3]
            -1.59790450e-15 * tc[4]
            +6.86281700e+00 ;
        /*species 4: CO */
        species[4] =
            +3.02507800e+00 * tc[0]
            +1.44268900e-03 * tc[1]
            -2.81541400e-07 * tc[2]
            +3.39527000e-11 * tc[3]
            -1.72773800e-15 * tc[4]
            +6.10821800e+00 ;
        /*species 5: H2 */
        species[5] =
            +2.99142300e+00 * tc[0]
            +7.00064400e-04 * tc[1]
            -2.81691450e-08 * tc[2]
            -3.07719267e-12 * tc[3]
            +3.95688000e-16 * tc[4]
            -1.35511000e+00 ;
        /*species 6: OH */
        species[6] =
            +2.88273000e+00 * tc[0]
            +1.01397400e-03 * tc[1]
            -1.13843850e-07 * tc[2]
            +7.24894667e-12 * tc[3]
            -1.28157625e-16 * tc[4]
            +5.59571200e+00 ;
        /*species 7: H2O2 */
        species[7] =
            +4.57316700e+00 * tc[0]
            +4.33613600e-03 * tc[1]
            -7.37344500e-07 * tc[2]
            +7.82968000e-11 * tc[3]
            -3.57913500e-15 * tc[4]
            +5.01137000e-01 ;
        /*species 8: HO2 */
        species[8] =
            +4.07219100e+00 * tc[0]
            +2.13129600e-03 * tc[1]
            -2.65407250e-07 * tc[2]
            +2.03742300e-11 * tc[3]
            -7.10291250e-16 * tc[4]
            +3.47602900e+00 ;
        /*species 9: H */
        species[9] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -4.60117600e-01 ;
        /*species 10: O */
        species[10] =
            +2.54206000e+00 * tc[0]
            -2.75506200e-05 * tc[1]
            -1.55140150e-09 * tc[2]
            +1.51702233e-12 * tc[3]
            -1.09201300e-16 * tc[4]
            +4.92030800e+00 ;
        /*species 11: CH3O */
        species[11] =
            +3.77080000e+00 * tc[0]
            +7.87149700e-03 * tc[1]
            -1.32819200e-06 * tc[2]
            +1.31481033e-10 * tc[3]
            -5.28154000e-15 * tc[4]
            +2.92957500e+00 ;
        /*species 12: CH2O */
        species[12] =
            +2.99560600e+00 * tc[0]
            +6.68132100e-03 * tc[1]
            -1.31447750e-06 * tc[2]
            +1.57905100e-10 * tc[3]
            -8.03129250e-15 * tc[4]
            +6.91257200e+00 ;
        /*species 13: HCO */
        species[13] =
            +3.55727100e+00 * tc[0]
            +3.34557300e-03 * tc[1]
            -6.67503000e-07 * tc[2]
            +8.23524333e-11 * tc[3]
            -4.28462750e-15 * tc[4]
            +5.55229900e+00 ;
        /*species 14: CH2 */
        species[14] =
            +3.63640800e+00 * tc[0]
            +1.93305700e-03 * tc[1]
            -8.43508000e-08 * tc[2]
            -3.36633000e-11 * tc[3]
            +4.52064000e-15 * tc[4]
            +2.15656100e+00 ;
        /*species 15: CH3 */
        species[15] =
            +2.84405200e+00 * tc[0]
            +6.13797400e-03 * tc[1]
            -1.11517250e-06 * tc[2]
            +1.26172033e-10 * tc[3]
            -6.13039750e-15 * tc[4]
            +5.45269700e+00 ;
        /*species 16: CH4 */
        species[16] =
            +1.68347900e+00 * tc[0]
            +1.02372400e-02 * tc[1]
            -1.93756450e-06 * tc[2]
            +2.26186167e-10 * tc[3]
            -1.12585575e-14 * tc[4]
            +9.62339500e+00 ;
        /*species 17: C2H3 */
        species[17] =
            +5.93346800e+00 * tc[0]
            +4.01774600e-03 * tc[1]
            -1.98337000e-07 * tc[2]
            -4.80422333e-11 * tc[3]
            +5.94661000e-15 * tc[4]
            -8.53031300e+00 ;
        /*species 18: C2H4 */
        species[18] =
            +3.52841900e+00 * tc[0]
            +1.14851800e-02 * tc[1]
            -2.20919250e-06 * tc[2]
            +2.61486700e-10 * tc[3]
            -1.31671200e-14 * tc[4]
            +2.23038900e+00 ;
        /*species 19: C2H5 */
        species[19] =
            +7.19048000e+00 * tc[0]
            +6.48407700e-03 * tc[1]
            -3.21403250e-07 * tc[2]
            -7.82626333e-11 * tc[3]
            +9.70219250e-15 * tc[4]
            -1.47808900e+01 ;
        /*species 20: C3H4 */
        species[20] =
            +6.31694869e+00 * tc[0]
            +1.11336262e-02 * tc[1]
            -1.98144509e-06 * tc[2]
            +2.11877925e-10 * tc[3]
            -9.46874713e-15 * tc[4]
            -1.09718862e+01 ;
        /*species 21: C3H5 */
        species[21] =
            +6.54761132e+00 * tc[0]
            +1.33152246e-02 * tc[1]
            -2.39166550e-06 * tc[2]
            +2.57316605e-10 * tc[3]
            -1.15482702e-14 * tc[4]
            -9.27486841e+00 ;
        /*species 23: C3H7 */
        species[23] =
            +7.70269870e+00 * tc[0]
            +1.60442030e-02 * tc[1]
            -2.64166100e-06 * tc[2]
            +2.54328633e-10 * tc[3]
            -9.84807100e-15 * tc[4]
            -1.54801800e+01 ;
        /*species 28: N2 */
        species[28] =
            +2.92664000e+00 * tc[0]
            +1.48797700e-03 * tc[1]
            -2.84238050e-07 * tc[2]
            +3.36568000e-11 * tc[3]
            -1.68833775e-15 * tc[4]
            +5.98052800e+00 ;
    }

    /*species with midpoint at T=1388 kelvin */
    if (T < 1388) {
        /*species 22: C3H6 */
        species[22] =
            +3.94615444e-01 * tc[0]
            +2.89107662e-02 * tc[1]
            -7.74434040e-06 * tc[2]
            +1.29604736e-09 * tc[3]
            -8.44725880e-14 * tc[4]
            +2.19003736e+01 ;
    } else {
        /*species 22: C3H6 */
        species[22] =
            +8.01595958e+00 * tc[0]
            +1.37023634e-02 * tc[1]
            -2.33124867e-06 * tc[2]
            +2.40418134e-10 * tc[3]
            -1.04342532e-14 * tc[4]
            -2.00160668e+01 ;
    }

    /*species with midpoint at T=1390 kelvin */
    if (T < 1390) {
        /*species 25: C7H15O2 */
        species[25] =
            +2.37499334e+00 * tc[0]
            +8.34651906e-02 * tc[1]
            -2.56948660e-05 * tc[2]
            +5.47392207e-09 * tc[3]
            -5.48763040e-13 * tc[4]
            +2.53067342e+01 ;
    } else {
        /*species 25: C7H15O2 */
        species[25] =
            +2.49023689e+01 * tc[0]
            +3.50716920e-02 * tc[1]
            -6.02201530e-06 * tc[2]
            +6.24882740e-10 * tc[3]
            -2.72369478e-14 * tc[4]
            -9.73923542e+01 ;
    }

    /*species with midpoint at T=1391 kelvin */
    if (T < 1391) {
        /*species 0: NC7H16 */
        species[0] =
            -1.26836187e+00 * tc[0]
            +8.54355820e-02 * tc[1]
            -2.62673393e-05 * tc[2]
            +5.43152403e-09 * tc[3]
            -5.05987313e-13 * tc[4]
            +3.53732912e+01 ;
    } else {
        /*species 0: NC7H16 */
        species[0] =
            +2.22148969e+01 * tc[0]
            +3.47675750e-02 * tc[1]
            -5.92035645e-06 * tc[2]
            +6.10994927e-10 * tc[3]
            -2.65325665e-14 * tc[4]
            -9.23040196e+01 ;
    }

    /*species with midpoint at T=1396 kelvin */
    if (T < 1396) {
        /*species 26: C7KET12 */
        species[26] =
            +5.82433697e-01 * tc[0]
            +1.01207869e-01 * tc[1]
            -3.82927998e-05 * tc[2]
            +1.00246202e-08 * tc[3]
            -1.20725698e-12 * tc[4]
            +3.33331449e+01 ;
    } else {
        /*species 26: C7KET12 */
        species[26] =
            +2.97472906e+01 * tc[0]
            +3.06622294e-02 * tc[1]
            -5.27817950e-06 * tc[2]
            +5.48757810e-10 * tc[3]
            -2.39542919e-14 * tc[4]
            -1.22432490e+02 ;
    }
    return;
}


/*save molecular weights into array */
void molecularWeight(double * restrict wt)
{
    wt[0] = 100.205570; /*NC7H16 */
    wt[1] = 31.998800; /*O2 */
    wt[2] = 44.009950; /*CO2 */
    wt[3] = 18.015340; /*H2O */
    wt[4] = 28.010550; /*CO */
    wt[5] = 2.015940; /*H2 */
    wt[6] = 17.007370; /*OH */
    wt[7] = 34.014740; /*H2O2 */
    wt[8] = 33.006770; /*HO2 */
    wt[9] = 1.007970; /*H */
    wt[10] = 15.999400; /*O */
    wt[11] = 31.034460; /*CH3O */
    wt[12] = 30.026490; /*CH2O */
    wt[13] = 29.018520; /*HCO */
    wt[14] = 14.027090; /*CH2 */
    wt[15] = 15.035060; /*CH3 */
    wt[16] = 16.043030; /*CH4 */
    wt[17] = 27.046210; /*C2H3 */
    wt[18] = 28.054180; /*C2H4 */
    wt[19] = 29.062150; /*C2H5 */
    wt[20] = 40.065330; /*C3H4 */
    wt[21] = 41.073300; /*C3H5 */
    wt[22] = 42.081270; /*C3H6 */
    wt[23] = 43.089240; /*C3H7 */
    wt[24] = 99.197600; /*C7H15-2 */
    wt[25] = 131.196400; /*C7H15O2 */
    wt[26] = 146.187830; /*C7KET12 */
    wt[27] = 99.153970; /*C5H11CO */
    wt[28] = 28.013400; /*N2 */

    return;
}


/*save atomic weights into array */
void atomicWeight(double * restrict awt)
{
    awt[0] = 1.007970; /*H */
    awt[1] = 12.011150; /*C */
    awt[2] = 15.999400; /*O */
    awt[3] = 14.006700; /*N */

    return;
}

/* get temperature given internal energy in mass units and mass fracs */
void VGET_T_GIVEN_EY(int * restrict np, double * restrict e, double * restrict y, double * restrict T, int * ierr)
{
#ifdef CONVERGENCE
    const int maxiter = 5000;
    const double tol  = 1.e-12;
#else
    const int maxiter = 200;
    const double tol  = 1.e-6;
#endif
    double tmin[*np], tmax[*np]; 
    double emin[*np], emax[*np];
    double cv,t1,e1,dt;
    int iwrk;
    double rwrk;
    int i; /* loop counter */

    for (int k=0; k<(*np); k++) {
        tmin[k] = 90;
        tmax[k] = 4000;
    }

    VCKUBMS(np, tmin, y, emin);
    VCKUBMS(np, tmax, y, emax);
    for (int k=0; k<(*np); k++) {
        //printf(" IN VGET_T_GIVEN_EY  k= %d T= %.12f ", k, T[k]);
        //printf(" emin= %.12f \n", emin[k]);
        if (e[k] < emin[k]) {
            printf(" NOT GOOD I %d %.12f \n", k, emin[k]);
            /*Linear Extrapolation below tmin */
            CKCVBS(&tmin[k], &y[k*(29)], &iwrk, &rwrk, &cv);
            T[k] = tmin[k] - (emin[k]-e[k])/cv;
            *ierr = 1;
            return;
        }
        if (e[k] > emax[k]) {
            printf(" NOT GOOD II %d %.12f \n", k, emax[k]);
            /*Linear Extrapolation above tmax */
            CKCVBS(&tmax[k],&y[k*(29)],&iwrk,&rwrk,&cv);
            T[k] = tmax[k] - (emax[k]-e[k])/cv;
            *ierr = 1;
            return;
        }
        t1 = T[k];
        if (t1 < tmin[k] || t1 > tmax[k]) {
            t1 = tmin[k] + (tmax[k]-tmin[k])/(emax[k]-emin[k])*(e[k]-emin[k]);
        }
        for (i = 0; i < maxiter; ++i) {
            CKUBMS(&t1,&y[k*(29)],&iwrk,&rwrk,&e1);
            CKCVBS(&t1,&y[k*(29)],&iwrk,&rwrk,&cv);
            dt = (e[k] - e1) / cv;
            if (dt > 100.) { dt = 100.; }
            else if (dt < -100.) { dt = -100.; }
            else if (fabs(dt) < tol) break;
            else if (t1+dt == t1) break;
            t1 += dt;
        }
        T[k] = t1;
        *ierr = 0;
    }

    return;
}

/* get temperature given internal energy in mass units and mass fracs */
void GET_T_GIVEN_EY(double * restrict e, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict t, int * ierr)
{
#ifdef CONVERGENCE
    const int maxiter = 5000;
    const double tol  = 1.e-12;
#else
    const int maxiter = 200;
    const double tol  = 1.e-6;
#endif
    double ein  = *e;
    double tmin = 90;/*max lower bound for thermo def */
    double tmax = 4000;/*min upper bound for thermo def */
    double e1,emin,emax,cv,t1,dt;
    int i;/* loop counter */

    CKUBMS(&tmin, y, iwrk, rwrk, &emin);
    CKUBMS(&tmax, y, iwrk, rwrk, &emax);
    //printf(" IN GET_T_GIVEN_EY  T= %.12f ", *t);
    //printf(" emin= %.12f \n", emin);
    if (ein < emin) {
        printf(" NOT GOOD I %.12f ", emin);
        /*Linear Extrapolation below tmin */
        CKCVBS(&tmin, y, iwrk, rwrk, &cv);
        *t = tmin - (emin-ein)/cv;
        *ierr = 1;
        return;
    }
    if (ein > emax) {
        printf(" NOT GOOD II %.12f ", emax);
        /*Linear Extrapolation above tmax */
        CKCVBS(&tmax, y, iwrk, rwrk, &cv);
        *t = tmax - (emax-ein)/cv;
        *ierr = 1;
        return;
    }
    t1 = *t;
    //printf(" Lu dans ERC_nheptane %.12f \n", t1);
    if (t1 < tmin || t1 > tmax) {
        t1 = tmin + (tmax-tmin)/(emax-emin)*(ein-emin);
    }
    for (i = 0; i < maxiter; ++i) {
        //printf(" Loop dans ERC_nheptane i t1: %d %4.4f \n", i, t1);
        //printf(" Loop dans ERC_nheptane YN2: %4.4f \n", y[28]);
        CKUBMS(&t1,y,iwrk,rwrk,&e1);
        CKCVBS(&t1,y,iwrk,rwrk,&cv);
        dt = (ein - e1) / cv;
        //printf(" Loop dans ERC_nheptane cv: %4.8e \n", cv);
        //printf(" Loop dans ERC_nheptane e1 dt:%4.8e %4.4f \n", e1, dt);
        if (dt > 100.) { dt = 100.; }
        else if (dt < -100.) { dt = -100.; }
        else if (fabs(dt) < tol) break;
        else if (t1+dt == t1) break;
        t1 += dt;
        //printf(" Loop dans ERC_nheptane %.12f ", t1);
    }
    *t = t1;
    *ierr = 0;
    return;
}

/* get temperature given enthalpy in mass units and mass fracs */
void GET_T_GIVEN_HY(double * restrict h, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict t, int * ierr)
{
#ifdef CONVERGENCE
    const int maxiter = 5000;
    const double tol  = 1.e-12;
#else
    const int maxiter = 200;
    const double tol  = 1.e-6;
#endif
    double hin  = *h;
    double tmin = 90;/*max lower bound for thermo def */
    double tmax = 4000;/*min upper bound for thermo def */
    double h1,hmin,hmax,cp,t1,dt;
    int i;/* loop counter */
    CKHBMS(&tmin, y, iwrk, rwrk, &hmin);
    CKHBMS(&tmax, y, iwrk, rwrk, &hmax);
    if (hin < hmin) {
        /*Linear Extrapolation below tmin */
        CKCPBS(&tmin, y, iwrk, rwrk, &cp);
        *t = tmin - (hmin-hin)/cp;
        *ierr = 1;
        return;
    }
    if (hin > hmax) {
        /*Linear Extrapolation above tmax */
        CKCPBS(&tmax, y, iwrk, rwrk, &cp);
        *t = tmax - (hmax-hin)/cp;
        *ierr = 1;
        return;
    }
    t1 = *t;
    if (t1 < tmin || t1 > tmax) {
        t1 = tmin + (tmax-tmin)/(hmax-hmin)*(hin-hmin);
    }
    for (i = 0; i < maxiter; ++i) {
        CKHBMS(&t1,y,iwrk,rwrk,&h1);
        CKCPBS(&t1,y,iwrk,rwrk,&cp);
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


/*compute the critical parameters for each species */
void GET_CRITPARAMS(double * restrict Tci, double * restrict ai, double * restrict bi, double * restrict acentric_i)
{

    double   EPS[29];
    double   SIG[29];
    double    wt[29];
    double avogadro = 6.02214199e23;
    double boltzmann = 1.3806503e-16; //we work in CGS
    double Rcst = 83.144598; //in bar [CGS] !

    egtransetEPS(EPS);
    egtransetSIG(SIG);
    molecularWeight(wt);

    /*species 0: NC7H16 */
    Tci[0] = 1.316 * EPS[0] ; 
    ai[0] = (5.55 * pow(avogadro,2.0) * EPS[0]*boltzmann * pow(1e-8*SIG[0],3.0) ) / (pow(wt[0],2.0)); 
    bi[0] = 0.855 * avogadro * pow(1e-8*SIG[0],3.0) / (wt[0]); 
    acentric_i[0] = 0.0 ;

    /*species 1: O2 */
    /*Imported from NIST */
    Tci[1] = 154.581000 ; 
    ai[1] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[1],2.0) / (pow(31.998800,2.0) * 50.430466); 
    bi[1] = 0.08664 * Rcst * Tci[1] / (31.998800 * 50.430466); 
    acentric_i[1] = 0.022200 ;

    /*species 2: CO2 */
    /*Imported from NIST */
    Tci[2] = 304.120000 ; 
    ai[2] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[2],2.0) / (pow(44.009950,2.0) * 73.740000); 
    bi[2] = 0.08664 * Rcst * Tci[2] / (44.009950 * 73.740000); 
    acentric_i[2] = 0.225000 ;

    /*species 3: H2O */
    /*Imported from NIST */
    Tci[3] = 647.096000 ; 
    ai[3] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[3],2.0) / (pow(18.015340,2.0) * 220.640000); 
    bi[3] = 0.08664 * Rcst * Tci[3] / (18.015340 * 220.640000); 
    acentric_i[3] = 0.344300 ;

    /*species 4: CO */
    /*Imported from NIST */
    Tci[4] = 132.850000 ; 
    ai[4] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[4],2.0) / (pow(28.010000,2.0) * 34.940000); 
    bi[4] = 0.08664 * Rcst * Tci[4] / (28.010000 * 34.940000); 
    acentric_i[4] = 0.045000 ;

    /*species 5: H2 */
    /*Imported from NIST */
    Tci[5] = 33.145000 ; 
    ai[5] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[5],2.0) / (pow(2.015880,2.0) * 12.964000); 
    bi[5] = 0.08664 * Rcst * Tci[5] / (2.015880 * 12.964000); 
    acentric_i[5] = -0.219000 ;

    /*species 6: OH */
    Tci[6] = 1.316 * EPS[6] ; 
    ai[6] = (5.55 * pow(avogadro,2.0) * EPS[6]*boltzmann * pow(1e-8*SIG[6],3.0) ) / (pow(wt[6],2.0)); 
    bi[6] = 0.855 * avogadro * pow(1e-8*SIG[6],3.0) / (wt[6]); 
    acentric_i[6] = 0.0 ;

    /*species 7: H2O2 */
    Tci[7] = 1.316 * EPS[7] ; 
    ai[7] = (5.55 * pow(avogadro,2.0) * EPS[7]*boltzmann * pow(1e-8*SIG[7],3.0) ) / (pow(wt[7],2.0)); 
    bi[7] = 0.855 * avogadro * pow(1e-8*SIG[7],3.0) / (wt[7]); 
    acentric_i[7] = 0.0 ;

    /*species 8: HO2 */
    Tci[8] = 1.316 * EPS[8] ; 
    ai[8] = (5.55 * pow(avogadro,2.0) * EPS[8]*boltzmann * pow(1e-8*SIG[8],3.0) ) / (pow(wt[8],2.0)); 
    bi[8] = 0.855 * avogadro * pow(1e-8*SIG[8],3.0) / (wt[8]); 
    acentric_i[8] = 0.0 ;

    /*species 9: H */
    Tci[9] = 1.316 * EPS[9] ; 
    ai[9] = (5.55 * pow(avogadro,2.0) * EPS[9]*boltzmann * pow(1e-8*SIG[9],3.0) ) / (pow(wt[9],2.0)); 
    bi[9] = 0.855 * avogadro * pow(1e-8*SIG[9],3.0) / (wt[9]); 
    acentric_i[9] = 0.0 ;

    /*species 10: O */
    Tci[10] = 1.316 * EPS[10] ; 
    ai[10] = (5.55 * pow(avogadro,2.0) * EPS[10]*boltzmann * pow(1e-8*SIG[10],3.0) ) / (pow(wt[10],2.0)); 
    bi[10] = 0.855 * avogadro * pow(1e-8*SIG[10],3.0) / (wt[10]); 
    acentric_i[10] = 0.0 ;

    /*species 11: CH3O */
    Tci[11] = 1.316 * EPS[11] ; 
    ai[11] = (5.55 * pow(avogadro,2.0) * EPS[11]*boltzmann * pow(1e-8*SIG[11],3.0) ) / (pow(wt[11],2.0)); 
    bi[11] = 0.855 * avogadro * pow(1e-8*SIG[11],3.0) / (wt[11]); 
    acentric_i[11] = 0.0 ;

    /*species 12: CH2O */
    Tci[12] = 1.316 * EPS[12] ; 
    ai[12] = (5.55 * pow(avogadro,2.0) * EPS[12]*boltzmann * pow(1e-8*SIG[12],3.0) ) / (pow(wt[12],2.0)); 
    bi[12] = 0.855 * avogadro * pow(1e-8*SIG[12],3.0) / (wt[12]); 
    acentric_i[12] = 0.0 ;

    /*species 13: HCO */
    Tci[13] = 1.316 * EPS[13] ; 
    ai[13] = (5.55 * pow(avogadro,2.0) * EPS[13]*boltzmann * pow(1e-8*SIG[13],3.0) ) / (pow(wt[13],2.0)); 
    bi[13] = 0.855 * avogadro * pow(1e-8*SIG[13],3.0) / (wt[13]); 
    acentric_i[13] = 0.0 ;

    /*species 14: CH2 */
    Tci[14] = 1.316 * EPS[14] ; 
    ai[14] = (5.55 * pow(avogadro,2.0) * EPS[14]*boltzmann * pow(1e-8*SIG[14],3.0) ) / (pow(wt[14],2.0)); 
    bi[14] = 0.855 * avogadro * pow(1e-8*SIG[14],3.0) / (wt[14]); 
    acentric_i[14] = 0.0 ;

    /*species 15: CH3 */
    Tci[15] = 1.316 * EPS[15] ; 
    ai[15] = (5.55 * pow(avogadro,2.0) * EPS[15]*boltzmann * pow(1e-8*SIG[15],3.0) ) / (pow(wt[15],2.0)); 
    bi[15] = 0.855 * avogadro * pow(1e-8*SIG[15],3.0) / (wt[15]); 
    acentric_i[15] = 0.0 ;

    /*species 16: CH4 */
    /*Imported from NIST */
    Tci[16] = 190.560000 ; 
    ai[16] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[16],2.0) / (pow(16.043030,2.0) * 45.990000); 
    bi[16] = 0.08664 * Rcst * Tci[16] / (16.043030 * 45.990000); 
    acentric_i[16] = 0.011000 ;

    /*species 17: C2H3 */
    Tci[17] = 1.316 * EPS[17] ; 
    ai[17] = (5.55 * pow(avogadro,2.0) * EPS[17]*boltzmann * pow(1e-8*SIG[17],3.0) ) / (pow(wt[17],2.0)); 
    bi[17] = 0.855 * avogadro * pow(1e-8*SIG[17],3.0) / (wt[17]); 
    acentric_i[17] = 0.0 ;

    /*species 18: C2H4 */
    /*Imported from NIST */
    Tci[18] = 282.340000 ; 
    ai[18] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[18],2.0) / (pow(28.054000,2.0) * 50.410000); 
    bi[18] = 0.08664 * Rcst * Tci[18] / (28.054000 * 50.410000); 
    acentric_i[18] = 0.087000 ;

    /*species 19: C2H5 */
    Tci[19] = 1.316 * EPS[19] ; 
    ai[19] = (5.55 * pow(avogadro,2.0) * EPS[19]*boltzmann * pow(1e-8*SIG[19],3.0) ) / (pow(wt[19],2.0)); 
    bi[19] = 0.855 * avogadro * pow(1e-8*SIG[19],3.0) / (wt[19]); 
    acentric_i[19] = 0.0 ;

    /*species 20: C3H4 */
    Tci[20] = 1.316 * EPS[20] ; 
    ai[20] = (5.55 * pow(avogadro,2.0) * EPS[20]*boltzmann * pow(1e-8*SIG[20],3.0) ) / (pow(wt[20],2.0)); 
    bi[20] = 0.855 * avogadro * pow(1e-8*SIG[20],3.0) / (wt[20]); 
    acentric_i[20] = 0.0 ;

    /*species 21: C3H5 */
    Tci[21] = 1.316 * EPS[21] ; 
    ai[21] = (5.55 * pow(avogadro,2.0) * EPS[21]*boltzmann * pow(1e-8*SIG[21],3.0) ) / (pow(wt[21],2.0)); 
    bi[21] = 0.855 * avogadro * pow(1e-8*SIG[21],3.0) / (wt[21]); 
    acentric_i[21] = 0.0 ;

    /*species 22: C3H6 */
    Tci[22] = 1.316 * EPS[22] ; 
    ai[22] = (5.55 * pow(avogadro,2.0) * EPS[22]*boltzmann * pow(1e-8*SIG[22],3.0) ) / (pow(wt[22],2.0)); 
    bi[22] = 0.855 * avogadro * pow(1e-8*SIG[22],3.0) / (wt[22]); 
    acentric_i[22] = 0.0 ;

    /*species 23: C3H7 */
    Tci[23] = 1.316 * EPS[23] ; 
    ai[23] = (5.55 * pow(avogadro,2.0) * EPS[23]*boltzmann * pow(1e-8*SIG[23],3.0) ) / (pow(wt[23],2.0)); 
    bi[23] = 0.855 * avogadro * pow(1e-8*SIG[23],3.0) / (wt[23]); 
    acentric_i[23] = 0.0 ;

    /*species 24: C7H15-2 */
    Tci[24] = 1.316 * EPS[24] ; 
    ai[24] = (5.55 * pow(avogadro,2.0) * EPS[24]*boltzmann * pow(1e-8*SIG[24],3.0) ) / (pow(wt[24],2.0)); 
    bi[24] = 0.855 * avogadro * pow(1e-8*SIG[24],3.0) / (wt[24]); 
    acentric_i[24] = 0.0 ;

    /*species 25: C7H15O2 */
    Tci[25] = 1.316 * EPS[25] ; 
    ai[25] = (5.55 * pow(avogadro,2.0) * EPS[25]*boltzmann * pow(1e-8*SIG[25],3.0) ) / (pow(wt[25],2.0)); 
    bi[25] = 0.855 * avogadro * pow(1e-8*SIG[25],3.0) / (wt[25]); 
    acentric_i[25] = 0.0 ;

    /*species 26: C7KET12 */
    Tci[26] = 1.316 * EPS[26] ; 
    ai[26] = (5.55 * pow(avogadro,2.0) * EPS[26]*boltzmann * pow(1e-8*SIG[26],3.0) ) / (pow(wt[26],2.0)); 
    bi[26] = 0.855 * avogadro * pow(1e-8*SIG[26],3.0) / (wt[26]); 
    acentric_i[26] = 0.0 ;

    /*species 27: C5H11CO */
    Tci[27] = 1.316 * EPS[27] ; 
    ai[27] = (5.55 * pow(avogadro,2.0) * EPS[27]*boltzmann * pow(1e-8*SIG[27],3.0) ) / (pow(wt[27],2.0)); 
    bi[27] = 0.855 * avogadro * pow(1e-8*SIG[27],3.0) / (wt[27]); 
    acentric_i[27] = 0.0 ;

    /*species 28: N2 */
    /*Imported from NIST */
    Tci[28] = 126.192000 ; 
    ai[28] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[28],2.0) / (pow(28.013400,2.0) * 33.958000); 
    bi[28] = 0.08664 * Rcst * Tci[28] / (28.013400 * 33.958000); 
    acentric_i[28] = 0.037200 ;

    return;
}


#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetLENIMC EGTRANSETLENIMC
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetLENIMC egtransetlenimc
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetLENIMC egtransetlenimc_
#endif
void egtransetLENIMC(int* LENIMC ) {
    *LENIMC = 118;}


#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetLENRMC EGTRANSETLENRMC
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetLENRMC egtransetlenrmc
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetLENRMC egtransetlenrmc_
#endif
void egtransetLENRMC(int* LENRMC ) {
    *LENRMC = 16994;}


#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetNO EGTRANSETNO
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetNO egtransetno
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetNO egtransetno_
#endif
void egtransetNO(int* NO ) {
    *NO = 4;}


#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetKK EGTRANSETKK
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetKK egtransetkk
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetKK egtransetkk_
#endif
void egtransetKK(int* KK ) {
    *KK = 29;}


#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetNLITE EGTRANSETNLITE
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetNLITE egtransetnlite
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetNLITE egtransetnlite_
#endif
void egtransetNLITE(int* NLITE ) {
    *NLITE = 2;}


/*Patm in ergs/cm3 */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetPATM EGTRANSETPATM
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetPATM egtransetpatm
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetPATM egtransetpatm_
#endif
void egtransetPATM(double* PATM) {
    *PATM =   0.1013250000000000E+07;}


/*the molecular weights in g/mol */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetWT EGTRANSETWT
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetWT egtransetwt
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetWT egtransetwt_
#endif
void egtransetWT(double* WT ) {
    WT[0] = 1.00205570E+02;
    WT[1] = 3.19988000E+01;
    WT[2] = 4.40099500E+01;
    WT[3] = 1.80153400E+01;
    WT[4] = 2.80105500E+01;
    WT[5] = 2.01594000E+00;
    WT[6] = 1.70073700E+01;
    WT[7] = 3.40147400E+01;
    WT[8] = 3.30067700E+01;
    WT[9] = 1.00797000E+00;
    WT[10] = 1.59994000E+01;
    WT[11] = 3.10344600E+01;
    WT[12] = 3.00264900E+01;
    WT[13] = 2.90185200E+01;
    WT[14] = 1.40270900E+01;
    WT[15] = 1.50350600E+01;
    WT[16] = 1.60430300E+01;
    WT[17] = 2.70462100E+01;
    WT[18] = 2.80541800E+01;
    WT[19] = 2.90621500E+01;
    WT[20] = 4.00653300E+01;
    WT[21] = 4.10733000E+01;
    WT[22] = 4.20812700E+01;
    WT[23] = 4.30892400E+01;
    WT[24] = 9.91976000E+01;
    WT[25] = 1.31196400E+02;
    WT[26] = 1.46187830E+02;
    WT[27] = 9.91539700E+01;
    WT[28] = 2.80134000E+01;
};


/*the lennard-jones potential well depth eps/kb in K */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetEPS EGTRANSETEPS
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetEPS egtranseteps
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetEPS egtranseteps_
#endif
void egtransetEPS(double* EPS ) {
    EPS[14] = 1.44000000E+02;
    EPS[15] = 1.44000000E+02;
    EPS[16] = 1.41400000E+02;
    EPS[17] = 2.09000000E+02;
    EPS[18] = 2.80800000E+02;
    EPS[19] = 2.52300000E+02;
    EPS[20] = 3.24800000E+02;
    EPS[21] = 3.16000000E+02;
    EPS[22] = 2.66800000E+02;
    EPS[23] = 2.66800000E+02;
    EPS[24] = 5.64030000E+02;
    EPS[25] = 5.61000000E+02;
    EPS[26] = 5.81300000E+02;
    EPS[27] = 4.98600000E+02;
    EPS[0] = 5.64030000E+02;
    EPS[1] = 1.07400000E+02;
    EPS[28] = 9.75300000E+01;
    EPS[2] = 2.44000000E+02;
    EPS[3] = 5.72400000E+02;
    EPS[4] = 9.81000000E+01;
    EPS[5] = 3.80000000E+01;
    EPS[6] = 8.00000000E+01;
    EPS[7] = 1.07400000E+02;
    EPS[8] = 1.07400000E+02;
    EPS[9] = 1.45000000E+02;
    EPS[10] = 8.00000000E+01;
    EPS[11] = 4.17000000E+02;
    EPS[12] = 4.98000000E+02;
    EPS[13] = 4.98000000E+02;
};


/*the lennard-jones collision diameter in Angstroms */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetSIG EGTRANSETSIG
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetSIG egtransetsig
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetSIG egtransetsig_
#endif
void egtransetSIG(double* SIG ) {
    SIG[14] = 3.80000000E+00;
    SIG[15] = 3.80000000E+00;
    SIG[16] = 3.74600000E+00;
    SIG[17] = 4.10000000E+00;
    SIG[18] = 3.97100000E+00;
    SIG[19] = 4.30200000E+00;
    SIG[20] = 4.29000000E+00;
    SIG[21] = 4.22000000E+00;
    SIG[22] = 4.98200000E+00;
    SIG[23] = 4.98200000E+00;
    SIG[24] = 6.00400000E+00;
    SIG[25] = 6.31700000E+00;
    SIG[26] = 6.50600000E+00;
    SIG[27] = 6.00900000E+00;
    SIG[0] = 6.00400000E+00;
    SIG[1] = 3.45800000E+00;
    SIG[28] = 3.62100000E+00;
    SIG[2] = 3.76300000E+00;
    SIG[3] = 2.60500000E+00;
    SIG[4] = 3.65000000E+00;
    SIG[5] = 2.92000000E+00;
    SIG[6] = 2.75000000E+00;
    SIG[7] = 3.45800000E+00;
    SIG[8] = 3.45800000E+00;
    SIG[9] = 2.05000000E+00;
    SIG[10] = 2.75000000E+00;
    SIG[11] = 3.69000000E+00;
    SIG[12] = 3.59000000E+00;
    SIG[13] = 3.59000000E+00;
};


/*the dipole moment in Debye */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetDIP EGTRANSETDIP
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetDIP egtransetdip
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetDIP egtransetdip_
#endif
void egtransetDIP(double* DIP ) {
    DIP[14] = 0.00000000E+00;
    DIP[15] = 0.00000000E+00;
    DIP[16] = 0.00000000E+00;
    DIP[17] = 0.00000000E+00;
    DIP[18] = 0.00000000E+00;
    DIP[19] = 0.00000000E+00;
    DIP[20] = 0.00000000E+00;
    DIP[21] = 0.00000000E+00;
    DIP[22] = 0.00000000E+00;
    DIP[23] = 0.00000000E+00;
    DIP[24] = 0.00000000E+00;
    DIP[25] = 1.70000000E+00;
    DIP[26] = 2.00000000E+00;
    DIP[27] = 2.00000000E+00;
    DIP[0] = 0.00000000E+00;
    DIP[1] = 0.00000000E+00;
    DIP[28] = 0.00000000E+00;
    DIP[2] = 0.00000000E+00;
    DIP[3] = 1.84400000E+00;
    DIP[4] = 0.00000000E+00;
    DIP[5] = 0.00000000E+00;
    DIP[6] = 0.00000000E+00;
    DIP[7] = 0.00000000E+00;
    DIP[8] = 0.00000000E+00;
    DIP[9] = 0.00000000E+00;
    DIP[10] = 0.00000000E+00;
    DIP[11] = 1.70000000E+00;
    DIP[12] = 0.00000000E+00;
    DIP[13] = 0.00000000E+00;
};


/*the polarizability in cubic Angstroms */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetPOL EGTRANSETPOL
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetPOL egtransetpol
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetPOL egtransetpol_
#endif
void egtransetPOL(double* POL ) {
    POL[14] = 0.00000000E+00;
    POL[15] = 0.00000000E+00;
    POL[16] = 2.60000000E+00;
    POL[17] = 0.00000000E+00;
    POL[18] = 0.00000000E+00;
    POL[19] = 0.00000000E+00;
    POL[20] = 0.00000000E+00;
    POL[21] = 0.00000000E+00;
    POL[22] = 0.00000000E+00;
    POL[23] = 0.00000000E+00;
    POL[24] = 0.00000000E+00;
    POL[25] = 0.00000000E+00;
    POL[26] = 0.00000000E+00;
    POL[27] = 0.00000000E+00;
    POL[0] = 0.00000000E+00;
    POL[1] = 1.60000000E+00;
    POL[28] = 1.76000000E+00;
    POL[2] = 2.65000000E+00;
    POL[3] = 0.00000000E+00;
    POL[4] = 1.95000000E+00;
    POL[5] = 7.90000000E-01;
    POL[6] = 0.00000000E+00;
    POL[7] = 0.00000000E+00;
    POL[8] = 0.00000000E+00;
    POL[9] = 0.00000000E+00;
    POL[10] = 0.00000000E+00;
    POL[11] = 0.00000000E+00;
    POL[12] = 0.00000000E+00;
    POL[13] = 0.00000000E+00;
};


/*the rotational relaxation collision number at 298 K */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetZROT EGTRANSETZROT
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetZROT egtransetzrot
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetZROT egtransetzrot_
#endif
void egtransetZROT(double* ZROT ) {
    ZROT[14] = 0.00000000E+00;
    ZROT[15] = 0.00000000E+00;
    ZROT[16] = 1.30000000E+01;
    ZROT[17] = 1.00000000E+00;
    ZROT[18] = 1.50000000E+00;
    ZROT[19] = 1.50000000E+00;
    ZROT[20] = 1.00000000E+00;
    ZROT[21] = 1.00000000E+00;
    ZROT[22] = 1.00000000E+00;
    ZROT[23] = 1.00000000E+00;
    ZROT[24] = 1.00000000E+00;
    ZROT[25] = 1.00000000E+00;
    ZROT[26] = 1.00000000E+00;
    ZROT[27] = 1.00000000E+00;
    ZROT[0] = 1.00000000E+00;
    ZROT[1] = 3.80000000E+00;
    ZROT[28] = 4.00000000E+00;
    ZROT[2] = 2.10000000E+00;
    ZROT[3] = 4.00000000E+00;
    ZROT[4] = 1.80000000E+00;
    ZROT[5] = 2.80000000E+02;
    ZROT[6] = 0.00000000E+00;
    ZROT[7] = 3.80000000E+00;
    ZROT[8] = 1.00000000E+00;
    ZROT[9] = 0.00000000E+00;
    ZROT[10] = 0.00000000E+00;
    ZROT[11] = 2.00000000E+00;
    ZROT[12] = 2.00000000E+00;
    ZROT[13] = 0.00000000E+00;
};


/*0: monoatomic, 1: linear, 2: nonlinear */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetNLIN EGTRANSETNLIN
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetNLIN egtransetnlin
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetNLIN egtransetnlin_
#endif
void egtransetNLIN(int* NLIN) {
    NLIN[14] = 1;
    NLIN[15] = 1;
    NLIN[16] = 2;
    NLIN[17] = 2;
    NLIN[18] = 2;
    NLIN[19] = 2;
    NLIN[20] = 1;
    NLIN[21] = 2;
    NLIN[22] = 2;
    NLIN[23] = 2;
    NLIN[24] = 2;
    NLIN[25] = 2;
    NLIN[26] = 2;
    NLIN[27] = 2;
    NLIN[0] = 2;
    NLIN[1] = 1;
    NLIN[28] = 1;
    NLIN[2] = 1;
    NLIN[3] = 2;
    NLIN[4] = 1;
    NLIN[5] = 1;
    NLIN[6] = 1;
    NLIN[7] = 2;
    NLIN[8] = 2;
    NLIN[9] = 0;
    NLIN[10] = 0;
    NLIN[11] = 2;
    NLIN[12] = 2;
    NLIN[13] = 2;
};


/*Poly fits for the viscosities, dim NO*KK */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFETA EGTRANSETCOFETA
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFETA egtransetcofeta
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFETA egtransetcofeta_
#endif
void egtransetCOFETA(double* COFETA) {
    COFETA[0] = -2.18094173E+01;
    COFETA[1] = 3.34626976E+00;
    COFETA[2] = -2.62425362E-01;
    COFETA[3] = 8.03040291E-03;
    COFETA[4] = -1.60066324E+01;
    COFETA[5] = 2.16753735E+00;
    COFETA[6] = -1.97226850E-01;
    COFETA[7] = 8.50065468E-03;
    COFETA[8] = -2.27427934E+01;
    COFETA[9] = 4.58846966E+00;
    COFETA[10] = -4.93198047E-01;
    COFETA[11] = 2.05723081E-02;
    COFETA[12] = -1.48387789E+01;
    COFETA[13] = 5.22431742E-01;
    COFETA[14] = 1.42282751E-01;
    COFETA[15] = -1.05730679E-02;
    COFETA[16] = -1.55744720E+01;
    COFETA[17] = 1.93951724E+00;
    COFETA[18] = -1.68049819E-01;
    COFETA[19] = 7.25705251E-03;
    COFETA[20] = -1.35655453E+01;
    COFETA[21] = 8.80923449E-01;
    COFETA[22] = -3.20611225E-02;
    COFETA[23] = 1.44621531E-03;
    COFETA[24] = -1.39702115E+01;
    COFETA[25] = 1.44311220E+00;
    COFETA[26] = -1.02767741E-01;
    COFETA[27] = 4.39662983E-03;
    COFETA[28] = -1.59760846E+01;
    COFETA[29] = 2.16753735E+00;
    COFETA[30] = -1.97226850E-01;
    COFETA[31] = 8.50065468E-03;
    COFETA[32] = -1.59911252E+01;
    COFETA[33] = 2.16753735E+00;
    COFETA[34] = -1.97226850E-01;
    COFETA[35] = 8.50065468E-03;
    COFETA[36] = -1.89413383E+01;
    COFETA[37] = 3.00307192E+00;
    COFETA[38] = -3.02604779E-01;
    COFETA[39] = 1.29293647E-02;
    COFETA[40] = -1.40007593E+01;
    COFETA[41] = 1.44311220E+00;
    COFETA[42] = -1.02767741E-01;
    COFETA[43] = 4.39662983E-03;
    COFETA[44] = -2.24087190E+01;
    COFETA[45] = 3.94674149E+00;
    COFETA[46] = -3.62445949E-01;
    COFETA[47] = 1.31677945E-02;
    COFETA[48] = -2.30504537E+01;
    COFETA[49] = 4.13126403E+00;
    COFETA[50] = -3.78083741E-01;
    COFETA[51] = 1.35206873E-02;
    COFETA[52] = -2.30675266E+01;
    COFETA[53] = 4.13126403E+00;
    COFETA[54] = -3.78083741E-01;
    COFETA[55] = 1.35206873E-02;
    COFETA[56] = -1.88109795E+01;
    COFETA[57] = 2.98574028E+00;
    COFETA[58] = -3.00482823E-01;
    COFETA[59] = 1.28428419E-02;
    COFETA[60] = -1.87762823E+01;
    COFETA[61] = 2.98574028E+00;
    COFETA[62] = -3.00482823E-01;
    COFETA[63] = 1.28428419E-02;
    COFETA[64] = -1.85911466E+01;
    COFETA[65] = 2.94126130E+00;
    COFETA[66] = -2.95057136E-01;
    COFETA[67] = 1.26224736E-02;
    COFETA[68] = -2.17979543E+01;
    COFETA[69] = 4.11859840E+00;
    COFETA[70] = -4.38141332E-01;
    COFETA[71] = 1.84208675E-02;
    COFETA[72] = -2.40988020E+01;
    COFETA[73] = 4.91550523E+00;
    COFETA[74] = -5.28299456E-01;
    COFETA[75] = 2.18076243E-02;
    COFETA[76] = -2.34718924E+01;
    COFETA[77] = 4.67141821E+00;
    COFETA[78] = -5.02333224E-01;
    COFETA[79] = 2.09043008E-02;
    COFETA[80] = -2.46718735E+01;
    COFETA[81] = 5.05499609E+00;
    COFETA[82] = -5.36738729E-01;
    COFETA[83] = 2.18055223E-02;
    COFETA[84] = -2.45524344E+01;
    COFETA[85] = 5.04519391E+00;
    COFETA[86] = -5.37454833E-01;
    COFETA[87] = 2.19116017E-02;
    COFETA[88] = -2.40598301E+01;
    COFETA[89] = 4.83130731E+00;
    COFETA[90] = -5.20343991E-01;
    COFETA[91] = 2.15775746E-02;
    COFETA[92] = -2.40479948E+01;
    COFETA[93] = 4.83130731E+00;
    COFETA[94] = -5.20343991E-01;
    COFETA[95] = 2.15775746E-02;
    COFETA[96] = -2.18144723E+01;
    COFETA[97] = 3.34626976E+00;
    COFETA[98] = -2.62425362E-01;
    COFETA[99] = 8.03040291E-03;
    COFETA[100] = -2.16087336E+01;
    COFETA[101] = 3.28666380E+00;
    COFETA[102] = -2.55354726E-01;
    COFETA[103] = 7.75191236E-03;
    COFETA[104] = -2.09433490E+01;
    COFETA[105] = 2.98470238E+00;
    COFETA[106] = -2.12074394E-01;
    COFETA[107] = 5.73810426E-03;
    COFETA[108] = -2.29641946E+01;
    COFETA[109] = 3.92103559E+00;
    COFETA[110] = -3.50088705E-01;
    COFETA[111] = 1.22904697E-02;
    COFETA[112] = -1.55270326E+01;
    COFETA[113] = 1.92766908E+00;
    COFETA[114] = -1.66518287E-01;
    COFETA[115] = 7.19100649E-03;
};


/*Poly fits for the conductivities, dim NO*KK */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFLAM EGTRANSETCOFLAM
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFLAM egtransetcoflam
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFLAM egtransetcoflam_
#endif
void egtransetCOFLAM(double* COFLAM) {
    COFLAM[0] = -2.22813109E+01;
    COFLAM[1] = 9.12871940E+00;
    COFLAM[2] = -8.55768554E-01;
    COFLAM[3] = 2.79644618E-02;
    COFLAM[4] = -2.11868240E+00;
    COFLAM[5] = 2.98567915E+00;
    COFLAM[6] = -2.86878052E-01;
    COFLAM[7] = 1.23850362E-02;
    COFLAM[8] = -1.40124764E+01;
    COFLAM[9] = 7.06760500E+00;
    COFLAM[10] = -7.44393576E-01;
    COFLAM[11] = 2.89770539E-02;
    COFLAM[12] = 1.63387231E+01;
    COFLAM[13] = -5.91396878E+00;
    COFLAM[14] = 1.09104916E+00;
    COFLAM[15] = -5.51653874E-02;
    COFLAM[16] = 5.83479004E+00;
    COFLAM[17] = -4.84788151E-01;
    COFLAM[18] = 2.11374672E-01;
    COFLAM[19] = -1.13934767E-02;
    COFLAM[20] = 1.31845927E+01;
    COFLAM[21] = -2.24986526E+00;
    COFLAM[22] = 3.81760448E-01;
    COFLAM[23] = -1.57362849E-02;
    COFLAM[24] = 1.65430880E+01;
    COFLAM[25] = -4.32533954E+00;
    COFLAM[26] = 6.96440219E-01;
    COFLAM[27] = -3.14178275E-02;
    COFLAM[28] = 7.15383393E-01;
    COFLAM[29] = 1.39783521E+00;
    COFLAM[30] = 8.83903407E-03;
    COFLAM[31] = -4.07837021E-03;
    COFLAM[32] = 3.75044177E+00;
    COFLAM[33] = 1.50600087E-01;
    COFLAM[34] = 1.62975594E-01;
    COFLAM[35] = -1.03397983E-02;
    COFLAM[36] = 6.08572159E-01;
    COFLAM[37] = 3.00307192E+00;
    COFLAM[38] = -3.02604779E-01;
    COFLAM[39] = 1.29293647E-02;
    COFLAM[40] = 2.78453831E+00;
    COFLAM[41] = 1.44311220E+00;
    COFLAM[42] = -1.02767741E-01;
    COFLAM[43] = 4.39662983E-03;
    COFLAM[44] = -1.44483522E+01;
    COFLAM[45] = 6.16025334E+00;
    COFLAM[46] = -4.76986670E-01;
    COFLAM[47] = 1.17790622E-02;
    COFLAM[48] = -9.38972921E+00;
    COFLAM[49] = 4.18204320E+00;
    COFLAM[50] = -2.28240995E-01;
    COFLAM[51] = 1.31619204E-03;
    COFLAM[52] = -6.17124562E+00;
    COFLAM[53] = 3.31735519E+00;
    COFLAM[54] = -1.78565326E-01;
    COFLAM[55] = 1.63476390E-03;
    COFLAM[56] = -6.26629604E-01;
    COFLAM[57] = 2.16914616E+00;
    COFLAM[58] = -1.38473121E-01;
    COFLAM[59] = 4.57177535E-03;
    COFLAM[60] = 6.46405356E+00;
    COFLAM[61] = -1.29894806E+00;
    COFLAM[62] = 4.16834334E-01;
    COFLAM[63] = -2.38865075E-02;
    COFLAM[64] = -4.63640986E-01;
    COFLAM[65] = 1.08053308E+00;
    COFLAM[66] = 1.59664255E-01;
    COFLAM[67] = -1.48274907E-02;
    COFLAM[68] = -1.16573779E+01;
    COFLAM[69] = 5.56358342E+00;
    COFLAM[70] = -4.57331283E-01;
    COFLAM[71] = 1.31263584E-02;
    COFLAM[72] = -1.95439042E+01;
    COFLAM[73] = 8.57472524E+00;
    COFLAM[74] = -8.32075705E-01;
    COFLAM[75] = 2.87186161E-02;
    COFLAM[76] = -1.50716942E+01;
    COFLAM[77] = 6.61727994E+00;
    COFLAM[78] = -5.53654394E-01;
    COFLAM[79] = 1.57358185E-02;
    COFLAM[80] = -1.83879647E+01;
    COFLAM[81] = 8.18864543E+00;
    COFLAM[82] = -7.96677658E-01;
    COFLAM[83] = 2.77175780E-02;
    COFLAM[84] = -1.64475887E+01;
    COFLAM[85] = 7.38762685E+00;
    COFLAM[86] = -6.80516749E-01;
    COFLAM[87] = 2.22121073E-02;
    COFLAM[88] = -2.01926953E+01;
    COFLAM[89] = 8.72799504E+00;
    COFLAM[90] = -8.49917931E-01;
    COFLAM[91] = 2.93793219E-02;
    COFLAM[92] = -2.01496659E+01;
    COFLAM[93] = 8.77441819E+00;
    COFLAM[94] = -8.59800261E-01;
    COFLAM[95] = 2.99701447E-02;
    COFLAM[96] = -2.20007450E+01;
    COFLAM[97] = 9.00745770E+00;
    COFLAM[98] = -8.40910150E-01;
    COFLAM[99] = 2.73666973E-02;
    COFLAM[100] = -1.83835038E+01;
    COFLAM[101] = 7.59449970E+00;
    COFLAM[102] = -6.63087665E-01;
    COFLAM[103] = 1.98637263E-02;
    COFLAM[104] = -2.09694020E+01;
    COFLAM[105] = 8.74632836E+00;
    COFLAM[106] = -8.34121078E-01;
    COFLAM[107] = 2.80553332E-02;
    COFLAM[108] = -1.94899188E+01;
    COFLAM[109] = 8.15247795E+00;
    COFLAM[110] = -7.50931652E-01;
    COFLAM[111] = 2.42229661E-02;
    COFLAM[112] = 7.60997504E+00;
    COFLAM[113] = -1.18418698E+00;
    COFLAM[114] = 3.03558703E-01;
    COFLAM[115] = -1.54159597E-02;
};


/*Poly fits for the diffusion coefficients, dim NO*KK*KK */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFD EGTRANSETCOFD
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFD egtransetcofd
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFD egtransetcofd_
#endif
void egtransetCOFD(double* COFD) {
    COFD[0] = -2.18912228E+01;
    COFD[1] = 4.47348600E+00;
    COFD[2] = -2.88470870E-01;
    COFD[3] = 9.60875755E-03;
    COFD[4] = -2.09423904E+01;
    COFD[5] = 4.99879868E+00;
    COFD[6] = -4.12878607E-01;
    COFD[7] = 1.70210595E-02;
    COFD[8] = -2.28053513E+01;
    COFD[9] = 5.42628508E+00;
    COFD[10] = -4.45690571E-01;
    COFD[11] = 1.75981564E-02;
    COFD[12] = -2.05565862E+01;
    COFD[13] = 4.43728272E+00;
    COFD[14] = -2.83087542E-01;
    COFD[15] = 9.35217321E-03;
    COFD[16] = -2.06913646E+01;
    COFD[17] = 4.92442201E+00;
    COFD[18] = -4.05069196E-01;
    COFD[19] = 1.67525806E-02;
    COFD[20] = -1.59937162E+01;
    COFD[21] = 3.76265798E+00;
    COFD[22] = -2.66912858E-01;
    COFD[23] = 1.12628530E-02;
    COFD[24] = -1.94895485E+01;
    COFD[25] = 4.64703793E+00;
    COFD[26] = -3.72283709E-01;
    COFD[27] = 1.54551501E-02;
    COFD[28] = -2.09653714E+01;
    COFD[29] = 4.99879868E+00;
    COFD[30] = -4.12878607E-01;
    COFD[31] = 1.70210595E-02;
    COFD[32] = -2.09540999E+01;
    COFD[33] = 4.99879868E+00;
    COFD[34] = -4.12878607E-01;
    COFD[35] = 1.70210595E-02;
    COFD[36] = -1.99238711E+01;
    COFD[37] = 5.27803851E+00;
    COFD[38] = -4.42140176E-01;
    COFD[39] = 1.80237084E-02;
    COFD[40] = -1.94633191E+01;
    COFD[41] = 4.64703793E+00;
    COFD[42] = -3.72283709E-01;
    COFD[43] = 1.54551501E-02;
    COFD[44] = -2.22225767E+01;
    COFD[45] = 5.05132223E+00;
    COFD[46] = -3.76867047E-01;
    COFD[47] = 1.39005636E-02;
    COFD[48] = -2.16117599E+01;
    COFD[49] = 4.74703304E+00;
    COFD[50] = -3.29522725E-01;
    COFD[51] = 1.15765554E-02;
    COFD[52] = -2.15985720E+01;
    COFD[53] = 4.74703304E+00;
    COFD[54] = -3.29522725E-01;
    COFD[55] = 1.15765554E-02;
    COFD[56] = -2.15568608E+01;
    COFD[57] = 5.27353202E+00;
    COFD[58] = -4.41738406E-01;
    COFD[59] = 1.80132159E-02;
    COFD[60] = -2.15871654E+01;
    COFD[61] = 5.27353202E+00;
    COFD[62] = -4.41738406E-01;
    COFD[63] = 1.80132159E-02;
    COFD[64] = -2.15620449E+01;
    COFD[65] = 5.26193907E+00;
    COFD[66] = -4.40716330E-01;
    COFD[67] = 1.79871592E-02;
    COFD[68] = -2.25748277E+01;
    COFD[69] = 5.43256416E+00;
    COFD[70] = -4.51725369E-01;
    COFD[71] = 1.80507219E-02;
    COFD[72] = -2.27314400E+01;
    COFD[73] = 5.39729337E+00;
    COFD[74] = -4.37101137E-01;
    COFD[75] = 1.70521181E-02;
    COFD[76] = -2.27800509E+01;
    COFD[77] = 5.42267938E+00;
    COFD[78] = -4.44086802E-01;
    COFD[79] = 1.74875185E-02;
    COFD[80] = -2.27953025E+01;
    COFD[81] = 5.29006863E+00;
    COFD[82] = -4.17624653E-01;
    COFD[83] = 1.60089937E-02;
    COFD[84] = -2.28363010E+01;
    COFD[85] = 5.31907985E+00;
    COFD[86] = -4.22514116E-01;
    COFD[87] = 1.62610921E-02;
    COFD[88] = -2.30934035E+01;
    COFD[89] = 5.42271039E+00;
    COFD[90] = -4.42221513E-01;
    COFD[91] = 1.73390244E-02;
    COFD[92] = -2.31017092E+01;
    COFD[93] = 5.42271039E+00;
    COFD[94] = -4.42221513E-01;
    COFD[95] = 1.73390244E-02;
    COFD[96] = -2.18886889E+01;
    COFD[97] = 4.47348600E+00;
    COFD[98] = -2.88470870E-01;
    COFD[99] = 9.60875755E-03;
    COFD[100] = -2.20335821E+01;
    COFD[101] = 4.48693408E+00;
    COFD[102] = -2.90462167E-01;
    COFD[103] = 9.70336833E-03;
    COFD[104] = -2.18977186E+01;
    COFD[105] = 4.39616382E+00;
    COFD[106] = -2.76994914E-01;
    COFD[107] = 9.06224705E-03;
    COFD[108] = -2.24418859E+01;
    COFD[109] = 4.74495988E+00;
    COFD[110] = -3.29202526E-01;
    COFD[111] = 1.15609276E-02;
    COFD[112] = -2.06668820E+01;
    COFD[113] = 4.91831350E+00;
    COFD[114] = -4.04379413E-01;
    COFD[115] = 1.67265792E-02;
    COFD[116] = -2.09423904E+01;
    COFD[117] = 4.99879868E+00;
    COFD[118] = -4.12878607E-01;
    COFD[119] = 1.70210595E-02;
    COFD[120] = -1.47079646E+01;
    COFD[121] = 3.10657376E+00;
    COFD[122] = -1.85922460E-01;
    COFD[123] = 7.92680827E-03;
    COFD[124] = -1.73301684E+01;
    COFD[125] = 3.98940337E+00;
    COFD[126] = -2.94368046E-01;
    COFD[127] = 1.23728707E-02;
    COFD[128] = -2.04397606E+01;
    COFD[129] = 5.23624904E+00;
    COFD[130] = -4.38651682E-01;
    COFD[131] = 1.79454253E-02;
    COFD[132] = -1.44544231E+01;
    COFD[133] = 3.00511224E+00;
    COFD[134] = -1.73057456E-01;
    COFD[135] = 7.38278721E-03;
    COFD[136] = -1.15529950E+01;
    COFD[137] = 2.41720750E+00;
    COFD[138] = -1.00182262E-01;
    COFD[139] = 4.37376763E-03;
    COFD[140] = -1.34974186E+01;
    COFD[141] = 2.80526779E+00;
    COFD[142] = -1.47524446E-01;
    COFD[143] = 6.29376912E-03;
    COFD[144] = -1.47230052E+01;
    COFD[145] = 3.10657376E+00;
    COFD[146] = -1.85922460E-01;
    COFD[147] = 7.92680827E-03;
    COFD[148] = -1.47156580E+01;
    COFD[149] = 3.10657376E+00;
    COFD[150] = -1.85922460E-01;
    COFD[151] = 7.92680827E-03;
    COFD[152] = -1.37328153E+01;
    COFD[153] = 3.42674072E+00;
    COFD[154] = -2.25944624E-01;
    COFD[155] = 9.59602387E-03;
    COFD[156] = -1.34772622E+01;
    COFD[157] = 2.80526779E+00;
    COFD[158] = -1.47524446E-01;
    COFD[159] = 6.29376912E-03;
    COFD[160] = -1.94217512E+01;
    COFD[161] = 4.74958634E+00;
    COFD[162] = -3.84601412E-01;
    COFD[163] = 1.59505248E-02;
    COFD[164] = -1.97777929E+01;
    COFD[165] = 4.88737199E+00;
    COFD[166] = -4.00855012E-01;
    COFD[167] = 1.65924977E-02;
    COFD[168] = -1.97689123E+01;
    COFD[169] = 4.88737199E+00;
    COFD[170] = -4.00855012E-01;
    COFD[171] = 1.65924977E-02;
    COFD[172] = -1.54101814E+01;
    COFD[173] = 3.41763870E+00;
    COFD[174] = -2.24801739E-01;
    COFD[175] = 9.54815836E-03;
    COFD[176] = -1.54340467E+01;
    COFD[177] = 3.41763870E+00;
    COFD[178] = -2.24801739E-01;
    COFD[179] = 9.54815836E-03;
    COFD[180] = -1.53751619E+01;
    COFD[181] = 3.39339822E+00;
    COFD[182] = -2.21757211E-01;
    COFD[183] = 9.42063333E-03;
    COFD[184] = -1.67688618E+01;
    COFD[185] = 3.79802425E+00;
    COFD[186] = -2.70990986E-01;
    COFD[187] = 1.14191222E-02;
    COFD[188] = -1.78073254E+01;
    COFD[189] = 4.17648998E+00;
    COFD[190] = -3.16975886E-01;
    COFD[191] = 1.32846969E-02;
    COFD[192] = -1.75004523E+01;
    COFD[193] = 4.03367299E+00;
    COFD[194] = -2.99744218E-01;
    COFD[195] = 1.25908330E-02;
    COFD[196] = -1.84552630E+01;
    COFD[197] = 4.34181535E+00;
    COFD[198] = -3.36587321E-01;
    COFD[199] = 1.40618044E-02;
    COFD[200] = -1.83508803E+01;
    COFD[201] = 4.30984197E+00;
    COFD[202] = -3.32770986E-01;
    COFD[203] = 1.39094747E-02;
    COFD[204] = -1.79735387E+01;
    COFD[205] = 4.11071307E+00;
    COFD[206] = -3.09074203E-01;
    COFD[207] = 1.29678777E-02;
    COFD[208] = -1.79786166E+01;
    COFD[209] = 4.11071307E+00;
    COFD[210] = -3.09074203E-01;
    COFD[211] = 1.29678777E-02;
    COFD[212] = -2.09411622E+01;
    COFD[213] = 4.99879868E+00;
    COFD[214] = -4.12878607E-01;
    COFD[215] = 1.70210595E-02;
    COFD[216] = -2.10558774E+01;
    COFD[217] = 5.00399764E+00;
    COFD[218] = -4.13379526E-01;
    COFD[219] = 1.70361425E-02;
    COFD[220] = -2.12327710E+01;
    COFD[221] = 5.04534270E+00;
    COFD[222] = -4.17903950E-01;
    COFD[223] = 1.71998906E-02;
    COFD[224] = -2.06704012E+01;
    COFD[225] = 4.91327235E+00;
    COFD[226] = -4.03805612E-01;
    COFD[227] = 1.67047562E-02;
    COFD[228] = -1.44285949E+01;
    COFD[229] = 2.99858376E+00;
    COFD[230] = -1.72232643E-01;
    COFD[231] = 7.34804765E-03;
    COFD[232] = -2.28053513E+01;
    COFD[233] = 5.42628508E+00;
    COFD[234] = -4.45690571E-01;
    COFD[235] = 1.75981564E-02;
    COFD[236] = -1.73301684E+01;
    COFD[237] = 3.98940337E+00;
    COFD[238] = -2.94368046E-01;
    COFD[239] = 1.23728707E-02;
    COFD[240] = -2.03925966E+01;
    COFD[241] = 4.98567306E+00;
    COFD[242] = -4.11537609E-01;
    COFD[243] = 1.69766778E-02;
    COFD[244] = -2.15572149E+01;
    COFD[245] = 5.36854850E+00;
    COFD[246] = -4.31373380E-01;
    COFD[247] = 1.67325146E-02;
    COFD[248] = -1.70267180E+01;
    COFD[249] = 3.87660436E+00;
    COFD[250] = -2.80627494E-01;
    COFD[251] = 1.18138535E-02;
    COFD[252] = -1.29858564E+01;
    COFD[253] = 2.88754338E+00;
    COFD[254] = -1.58249632E-01;
    COFD[255] = 6.76021204E-03;
    COFD[256] = -1.60407063E+01;
    COFD[257] = 3.68166087E+00;
    COFD[258] = -2.57364559E-01;
    COFD[259] = 1.08881828E-02;
    COFD[260] = -1.73476278E+01;
    COFD[261] = 3.98940337E+00;
    COFD[262] = -2.94368046E-01;
    COFD[263] = 1.23728707E-02;
    COFD[264] = -1.73390885E+01;
    COFD[265] = 3.98940337E+00;
    COFD[266] = -2.94368046E-01;
    COFD[267] = 1.23728707E-02;
    COFD[268] = -1.64909711E+01;
    COFD[269] = 4.36229440E+00;
    COFD[270] = -3.39095464E-01;
    COFD[271] = 1.41645916E-02;
    COFD[272] = -1.60184872E+01;
    COFD[273] = 3.68166087E+00;
    COFD[274] = -2.57364559E-01;
    COFD[275] = 1.08881828E-02;
    COFD[276] = -2.17139201E+01;
    COFD[277] = 5.38758317E+00;
    COFD[278] = -4.48340933E-01;
    COFD[279] = 1.79944740E-02;
    COFD[280] = -2.19036727E+01;
    COFD[281] = 5.44375518E+00;
    COFD[282] = -4.52244836E-01;
    COFD[283] = 1.80389212E-02;
    COFD[284] = -2.18934539E+01;
    COFD[285] = 5.44375518E+00;
    COFD[286] = -4.52244836E-01;
    COFD[287] = 1.80389212E-02;
    COFD[288] = -1.81750360E+01;
    COFD[289] = 4.35070874E+00;
    COFD[290] = -3.37658377E-01;
    COFD[291] = 1.41049661E-02;
    COFD[292] = -1.82011239E+01;
    COFD[293] = 4.35070874E+00;
    COFD[294] = -3.37658377E-01;
    COFD[295] = 1.41049661E-02;
    COFD[296] = -1.81478704E+01;
    COFD[297] = 4.32863219E+00;
    COFD[298] = -3.35006470E-01;
    COFD[299] = 1.39983891E-02;
    COFD[300] = -1.98631073E+01;
    COFD[301] = 4.82865226E+00;
    COFD[302] = -3.94029099E-01;
    COFD[303] = 1.63270142E-02;
    COFD[304] = -2.08088198E+01;
    COFD[305] = 5.14475868E+00;
    COFD[306] = -4.29211829E-01;
    COFD[307] = 1.76287935E-02;
    COFD[308] = -2.04989678E+01;
    COFD[309] = 5.01033357E+00;
    COFD[310] = -4.14008048E-01;
    COFD[311] = 1.70559642E-02;
    COFD[312] = -2.13743701E+01;
    COFD[313] = 5.25833649E+00;
    COFD[314] = -4.40410969E-01;
    COFD[315] = 1.79800167E-02;
    COFD[316] = -2.13081407E+01;
    COFD[317] = 5.24403923E+00;
    COFD[318] = -4.39272407E-01;
    COFD[319] = 1.79576065E-02;
    COFD[320] = -2.09751443E+01;
    COFD[321] = 5.08056463E+00;
    COFD[322] = -4.21952793E-01;
    COFD[323] = 1.73552764E-02;
    COFD[324] = -2.09811595E+01;
    COFD[325] = 5.08056463E+00;
    COFD[326] = -4.21952793E-01;
    COFD[327] = 1.73552764E-02;
    COFD[328] = -2.28038033E+01;
    COFD[329] = 5.42628508E+00;
    COFD[330] = -4.45690571E-01;
    COFD[331] = 1.75981564E-02;
    COFD[332] = -2.29076429E+01;
    COFD[333] = 5.42460337E+00;
    COFD[334] = -4.45270089E-01;
    COFD[335] = 1.75725973E-02;
    COFD[336] = -2.29881956E+01;
    COFD[337] = 5.42338458E+00;
    COFD[338] = -4.43824576E-01;
    COFD[339] = 1.74633100E-02;
    COFD[340] = -2.27673631E+01;
    COFD[341] = 5.44729221E+00;
    COFD[342] = -4.52058774E-01;
    COFD[343] = 1.80072303E-02;
    COFD[344] = -1.69985939E+01;
    COFD[345] = 3.86929427E+00;
    COFD[346] = -2.79728651E-01;
    COFD[347] = 1.17769347E-02;
    COFD[348] = -2.05565862E+01;
    COFD[349] = 4.43728272E+00;
    COFD[350] = -2.83087542E-01;
    COFD[351] = 9.35217321E-03;
    COFD[352] = -2.04397606E+01;
    COFD[353] = 5.23624904E+00;
    COFD[354] = -4.38651682E-01;
    COFD[355] = 1.79454253E-02;
    COFD[356] = -2.15572149E+01;
    COFD[357] = 5.36854850E+00;
    COFD[358] = -4.31373380E-01;
    COFD[359] = 1.67325146E-02;
    COFD[360] = -1.58086292E+01;
    COFD[361] = 2.65424082E+00;
    COFD[362] = -1.10243743E-02;
    COFD[363] = -3.68457990E-03;
    COFD[364] = -2.02791813E+01;
    COFD[365] = 5.17643426E+00;
    COFD[366] = -4.32650593E-01;
    COFD[367] = 1.77521858E-02;
    COFD[368] = -1.60490005E+01;
    COFD[369] = 4.12894470E+00;
    COFD[370] = -3.11280304E-01;
    COFD[371] = 1.30569853E-02;
    COFD[372] = -1.83035172E+01;
    COFD[373] = 4.66435737E+00;
    COFD[374] = -3.74296088E-01;
    COFD[375] = 1.55332609E-02;
    COFD[376] = -1.97265566E+01;
    COFD[377] = 5.00923653E+00;
    COFD[378] = -4.13897836E-01;
    COFD[379] = 1.70524216E-02;
    COFD[380] = -1.97212974E+01;
    COFD[381] = 5.00923653E+00;
    COFD[382] = -4.13897836E-01;
    COFD[383] = 1.70524216E-02;
    COFD[384] = -1.88398861E+01;
    COFD[385] = 5.28766000E+00;
    COFD[386] = -4.43001094E-01;
    COFD[387] = 1.80463788E-02;
    COFD[388] = -1.82875708E+01;
    COFD[389] = 4.66435737E+00;
    COFD[390] = -3.74296088E-01;
    COFD[391] = 1.55332609E-02;
    COFD[392] = -1.93199869E+01;
    COFD[393] = 4.21681909E+00;
    COFD[394] = -2.50231374E-01;
    COFD[395] = 7.81267943E-03;
    COFD[396] = -2.03286949E+01;
    COFD[397] = 4.72161925E+00;
    COFD[398] = -3.25604803E-01;
    COFD[399] = 1.13856236E-02;
    COFD[400] = -2.03222242E+01;
    COFD[401] = 4.72161925E+00;
    COFD[402] = -3.25604803E-01;
    COFD[403] = 1.13856236E-02;
    COFD[404] = -2.05176849E+01;
    COFD[405] = 5.28313423E+00;
    COFD[406] = -4.42595579E-01;
    COFD[407] = 1.80356675E-02;
    COFD[408] = -2.05368957E+01;
    COFD[409] = 5.28313423E+00;
    COFD[410] = -4.42595579E-01;
    COFD[411] = 1.80356675E-02;
    COFD[412] = -2.09577528E+01;
    COFD[413] = 5.37851020E+00;
    COFD[414] = -4.48937706E-01;
    COFD[415] = 1.80857999E-02;
    COFD[416] = -2.14416154E+01;
    COFD[417] = 5.43806985E+00;
    COFD[418] = -4.51984607E-01;
    COFD[419] = 1.80451889E-02;
    COFD[420] = -2.15465990E+01;
    COFD[421] = 5.38978785E+00;
    COFD[422] = -4.35592956E-01;
    COFD[423] = 1.69676354E-02;
    COFD[424] = -2.16415067E+01;
    COFD[425] = 5.42347394E+00;
    COFD[426] = -4.43696048E-01;
    COFD[427] = 1.74526528E-02;
    COFD[428] = -2.15592524E+01;
    COFD[429] = 5.27728120E+00;
    COFD[430] = -4.15361581E-01;
    COFD[431] = 1.58893774E-02;
    COFD[432] = -2.15823819E+01;
    COFD[433] = 5.30395665E+00;
    COFD[434] = -4.19948713E-01;
    COFD[435] = 1.61283165E-02;
    COFD[436] = -2.19259317E+01;
    COFD[437] = 5.41758144E+00;
    COFD[438] = -4.41053871E-01;
    COFD[439] = 1.72705882E-02;
    COFD[440] = -2.19294503E+01;
    COFD[441] = 5.41758144E+00;
    COFD[442] = -4.41053871E-01;
    COFD[443] = 1.72705882E-02;
    COFD[444] = -2.05558125E+01;
    COFD[445] = 4.43728272E+00;
    COFD[446] = -2.83087542E-01;
    COFD[447] = 9.35217321E-03;
    COFD[448] = -2.01811055E+01;
    COFD[449] = 4.24520878E+00;
    COFD[450] = -2.57219171E-01;
    COFD[451] = 8.20894239E-03;
    COFD[452] = -1.99759660E+01;
    COFD[453] = 4.12684449E+00;
    COFD[454] = -2.39932282E-01;
    COFD[455] = 7.39544059E-03;
    COFD[456] = -2.04086296E+01;
    COFD[457] = 4.40679937E+00;
    COFD[458] = -2.81870764E-01;
    COFD[459] = 9.40471579E-03;
    COFD[460] = -2.01962832E+01;
    COFD[461] = 5.15333309E+00;
    COFD[462] = -4.30146489E-01;
    COFD[463] = 1.76625026E-02;
    COFD[464] = -2.06913646E+01;
    COFD[465] = 4.92442201E+00;
    COFD[466] = -4.05069196E-01;
    COFD[467] = 1.67525806E-02;
    COFD[468] = -1.44544231E+01;
    COFD[469] = 3.00511224E+00;
    COFD[470] = -1.73057456E-01;
    COFD[471] = 7.38278721E-03;
    COFD[472] = -1.70267180E+01;
    COFD[473] = 3.87660436E+00;
    COFD[474] = -2.80627494E-01;
    COFD[475] = 1.18138535E-02;
    COFD[476] = -2.02791813E+01;
    COFD[477] = 5.17643426E+00;
    COFD[478] = -4.32650593E-01;
    COFD[479] = 1.77521858E-02;
    COFD[480] = -1.42466030E+01;
    COFD[481] = 2.92213383E+00;
    COFD[482] = -1.62707676E-01;
    COFD[483] = 6.95211737E-03;
    COFD[484] = -1.15519070E+01;
    COFD[485] = 2.40345064E+00;
    COFD[486] = -9.91400939E-02;
    COFD[487] = 4.36235570E-03;
    COFD[488] = -1.32712433E+01;
    COFD[489] = 2.70501014E+00;
    COFD[490] = -1.34413475E-01;
    COFD[491] = 5.72152293E-03;
    COFD[492] = -1.44684500E+01;
    COFD[493] = 3.00511224E+00;
    COFD[494] = -1.73057456E-01;
    COFD[495] = 7.38278721E-03;
    COFD[496] = -1.44616016E+01;
    COFD[497] = 3.00511224E+00;
    COFD[498] = -1.73057456E-01;
    COFD[499] = 7.38278721E-03;
    COFD[500] = -1.34773721E+01;
    COFD[501] = 3.30798941E+00;
    COFD[502] = -2.10999668E-01;
    COFD[503] = 8.96870768E-03;
    COFD[504] = -1.32520179E+01;
    COFD[505] = 2.70501014E+00;
    COFD[506] = -1.34413475E-01;
    COFD[507] = 5.72152293E-03;
    COFD[508] = -1.90947006E+01;
    COFD[509] = 4.62819535E+00;
    COFD[510] = -3.70101251E-01;
    COFD[511] = 1.53707619E-02;
    COFD[512] = -1.94491481E+01;
    COFD[513] = 4.76713321E+00;
    COFD[514] = -3.86705929E-01;
    COFD[515] = 1.60350491E-02;
    COFD[516] = -1.94408353E+01;
    COFD[517] = 4.76713321E+00;
    COFD[518] = -3.86705929E-01;
    COFD[519] = 1.60350491E-02;
    COFD[520] = -1.51284988E+01;
    COFD[521] = 3.30256227E+00;
    COFD[522] = -2.10347339E-01;
    COFD[523] = 8.94258391E-03;
    COFD[524] = -1.51513486E+01;
    COFD[525] = 3.30256227E+00;
    COFD[526] = -2.10347339E-01;
    COFD[527] = 8.94258391E-03;
    COFD[528] = -1.51167688E+01;
    COFD[529] = 3.28831215E+00;
    COFD[530] = -2.08633772E-01;
    COFD[531] = 8.87391434E-03;
    COFD[532] = -1.65688462E+01;
    COFD[533] = 3.72373814E+00;
    COFD[534] = -2.62310192E-01;
    COFD[535] = 1.10815606E-02;
    COFD[536] = -1.74857137E+01;
    COFD[537] = 4.05605079E+00;
    COFD[538] = -3.02454805E-01;
    COFD[539] = 1.27004015E-02;
    COFD[540] = -1.71875552E+01;
    COFD[541] = 3.91643623E+00;
    COFD[542] = -2.85489387E-01;
    COFD[543] = 1.20121308E-02;
    COFD[544] = -1.81728784E+01;
    COFD[545] = 4.24023311E+00;
    COFD[546] = -3.24539148E-01;
    COFD[547] = 1.35842242E-02;
    COFD[548] = -1.80716029E+01;
    COFD[549] = 4.20964822E+00;
    COFD[550] = -3.20930845E-01;
    COFD[551] = 1.34421679E-02;
    COFD[552] = -1.76339606E+01;
    COFD[553] = 3.98770708E+00;
    COFD[554] = -2.94160646E-01;
    COFD[555] = 1.23644039E-02;
    COFD[556] = -1.76386567E+01;
    COFD[557] = 3.98770708E+00;
    COFD[558] = -2.94160646E-01;
    COFD[559] = 1.23644039E-02;
    COFD[560] = -2.06902559E+01;
    COFD[561] = 4.92442201E+00;
    COFD[562] = -4.05069196E-01;
    COFD[563] = 1.67525806E-02;
    COFD[564] = -2.08064719E+01;
    COFD[565] = 4.93178191E+00;
    COFD[566] = -4.05846844E-01;
    COFD[567] = 1.67795427E-02;
    COFD[568] = -2.09575447E+01;
    COFD[569] = 4.96410963E+00;
    COFD[570] = -4.09246419E-01;
    COFD[571] = 1.68966000E-02;
    COFD[572] = -2.03384381E+01;
    COFD[573] = 4.80405736E+00;
    COFD[574] = -3.91110136E-01;
    COFD[575] = 1.62109876E-02;
    COFD[576] = -1.42266828E+01;
    COFD[577] = 2.91778615E+00;
    COFD[578] = -1.62159133E-01;
    COFD[579] = 6.92895656E-03;
    COFD[580] = -1.59937162E+01;
    COFD[581] = 3.76265798E+00;
    COFD[582] = -2.66912858E-01;
    COFD[583] = 1.12628530E-02;
    COFD[584] = -1.15529950E+01;
    COFD[585] = 2.41720750E+00;
    COFD[586] = -1.00182262E-01;
    COFD[587] = 4.37376763E-03;
    COFD[588] = -1.29858564E+01;
    COFD[589] = 2.88754338E+00;
    COFD[590] = -1.58249632E-01;
    COFD[591] = 6.76021204E-03;
    COFD[592] = -1.60490005E+01;
    COFD[593] = 4.12894470E+00;
    COFD[594] = -3.11280304E-01;
    COFD[595] = 1.30569853E-02;
    COFD[596] = -1.15519070E+01;
    COFD[597] = 2.40345064E+00;
    COFD[598] = -9.91400939E-02;
    COFD[599] = 4.36235570E-03;
    COFD[600] = -9.61322634E+00;
    COFD[601] = 1.87738363E+00;
    COFD[602] = -2.92985344E-02;
    COFD[603] = 1.27857667E-03;
    COFD[604] = -1.09748446E+01;
    COFD[605] = 2.31040527E+00;
    COFD[606] = -8.74442763E-02;
    COFD[607] = 3.87525190E-03;
    COFD[608] = -1.15547543E+01;
    COFD[609] = 2.41720750E+00;
    COFD[610] = -1.00182262E-01;
    COFD[611] = 4.37376763E-03;
    COFD[612] = -1.15539008E+01;
    COFD[613] = 2.41720750E+00;
    COFD[614] = -1.00182262E-01;
    COFD[615] = 4.37376763E-03;
    COFD[616] = -1.05322133E+01;
    COFD[617] = 2.38287287E+00;
    COFD[618] = -9.25451677E-02;
    COFD[619] = 3.90477283E-03;
    COFD[620] = -1.09715176E+01;
    COFD[621] = 2.31040527E+00;
    COFD[622] = -8.74442763E-02;
    COFD[623] = 3.87525190E-03;
    COFD[624] = -1.48338675E+01;
    COFD[625] = 3.57131810E+00;
    COFD[626] = -2.43940431E-01;
    COFD[627] = 1.03434646E-02;
    COFD[628] = -1.50226829E+01;
    COFD[629] = 3.65330206E+00;
    COFD[630] = -2.53985209E-01;
    COFD[631] = 1.07540119E-02;
    COFD[632] = -1.50215914E+01;
    COFD[633] = 3.65330206E+00;
    COFD[634] = -2.53985209E-01;
    COFD[635] = 1.07540119E-02;
    COFD[636] = -1.16100827E+01;
    COFD[637] = 2.38028423E+00;
    COFD[638] = -9.22299290E-02;
    COFD[639] = 3.89194284E-03;
    COFD[640] = -1.16143127E+01;
    COFD[641] = 2.38028423E+00;
    COFD[642] = -9.22299290E-02;
    COFD[643] = 3.89194284E-03;
    COFD[644] = -1.15834311E+01;
    COFD[645] = 2.37417315E+00;
    COFD[646] = -9.15013712E-02;
    COFD[647] = 3.86299275E-03;
    COFD[648] = -1.26266190E+01;
    COFD[649] = 2.71852437E+00;
    COFD[650] = -1.36179939E-01;
    COFD[651] = 5.79856970E-03;
    COFD[652] = -1.33954202E+01;
    COFD[653] = 3.01958105E+00;
    COFD[654] = -1.74891189E-01;
    COFD[655] = 7.46027223E-03;
    COFD[656] = -1.32148679E+01;
    COFD[657] = 2.91941286E+00;
    COFD[658] = -1.62365908E-01;
    COFD[659] = 6.93774959E-03;
    COFD[660] = -1.39240907E+01;
    COFD[661] = 3.17763277E+00;
    COFD[662] = -1.94889823E-01;
    COFD[663] = 8.30436995E-03;
    COFD[664] = -1.38287706E+01;
    COFD[665] = 3.14964252E+00;
    COFD[666] = -1.91369421E-01;
    COFD[667] = 8.15663173E-03;
    COFD[668] = -1.35245705E+01;
    COFD[669] = 2.96195761E+00;
    COFD[670] = -1.67607052E-01;
    COFD[671] = 7.15329427E-03;
    COFD[672] = -1.35251055E+01;
    COFD[673] = 2.96195761E+00;
    COFD[674] = -1.67607052E-01;
    COFD[675] = 7.15329427E-03;
    COFD[676] = -1.59936160E+01;
    COFD[677] = 3.76265798E+00;
    COFD[678] = -2.66912858E-01;
    COFD[679] = 1.12628530E-02;
    COFD[680] = -1.60961940E+01;
    COFD[681] = 3.77271143E+00;
    COFD[682] = -2.68084668E-01;
    COFD[683] = 1.13083251E-02;
    COFD[684] = -1.62329249E+01;
    COFD[685] = 3.80480708E+00;
    COFD[686] = -2.71809953E-01;
    COFD[687] = 1.14521282E-02;
    COFD[688] = -1.57635543E+01;
    COFD[689] = 3.68244512E+00;
    COFD[690] = -2.57458142E-01;
    COFD[691] = 1.08919055E-02;
    COFD[692] = -1.15279021E+01;
    COFD[693] = 2.39753083E+00;
    COFD[694] = -9.83466917E-02;
    COFD[695] = 4.32694147E-03;
    COFD[696] = -1.94895485E+01;
    COFD[697] = 4.64703793E+00;
    COFD[698] = -3.72283709E-01;
    COFD[699] = 1.54551501E-02;
    COFD[700] = -1.34974186E+01;
    COFD[701] = 2.80526779E+00;
    COFD[702] = -1.47524446E-01;
    COFD[703] = 6.29376912E-03;
    COFD[704] = -1.60407063E+01;
    COFD[705] = 3.68166087E+00;
    COFD[706] = -2.57364559E-01;
    COFD[707] = 1.08881828E-02;
    COFD[708] = -1.83035172E+01;
    COFD[709] = 4.66435737E+00;
    COFD[710] = -3.74296088E-01;
    COFD[711] = 1.55332609E-02;
    COFD[712] = -1.32712433E+01;
    COFD[713] = 2.70501014E+00;
    COFD[714] = -1.34413475E-01;
    COFD[715] = 5.72152293E-03;
    COFD[716] = -1.09748446E+01;
    COFD[717] = 2.31040527E+00;
    COFD[718] = -8.74442763E-02;
    COFD[719] = 3.87525190E-03;
    COFD[720] = -1.22981977E+01;
    COFD[721] = 2.49197560E+00;
    COFD[722] = -1.06603364E-01;
    COFD[723] = 4.50971238E-03;
    COFD[724] = -1.35078099E+01;
    COFD[725] = 2.80526779E+00;
    COFD[726] = -1.47524446E-01;
    COFD[727] = 6.29376912E-03;
    COFD[728] = -1.35027460E+01;
    COFD[729] = 2.80526779E+00;
    COFD[730] = -1.47524446E-01;
    COFD[731] = 6.29376912E-03;
    COFD[732] = -1.25828466E+01;
    COFD[733] = 3.11263114E+00;
    COFD[734] = -1.86688345E-01;
    COFD[735] = 7.95911073E-03;
    COFD[736] = -1.22826905E+01;
    COFD[737] = 2.49197560E+00;
    COFD[738] = -1.06603364E-01;
    COFD[739] = 4.50971238E-03;
    COFD[740] = -1.76970903E+01;
    COFD[741] = 4.29054603E+00;
    COFD[742] = -3.30479889E-01;
    COFD[743] = 1.38185279E-02;
    COFD[744] = -1.83085512E+01;
    COFD[745] = 4.51789863E+00;
    COFD[746] = -3.57703548E-01;
    COFD[747] = 1.49085039E-02;
    COFD[748] = -1.83023101E+01;
    COFD[749] = 4.51789863E+00;
    COFD[750] = -3.57703548E-01;
    COFD[751] = 1.49085039E-02;
    COFD[752] = -1.42289701E+01;
    COFD[753] = 3.10519183E+00;
    COFD[754] = -1.85747616E-01;
    COFD[755] = 7.91942899E-03;
    COFD[756] = -1.42476859E+01;
    COFD[757] = 3.10519183E+00;
    COFD[758] = -1.85747616E-01;
    COFD[759] = 7.91942899E-03;
    COFD[760] = -1.41940170E+01;
    COFD[761] = 3.08518563E+00;
    COFD[762] = -1.83211692E-01;
    COFD[763] = 7.81220123E-03;
    COFD[764] = -1.55934124E+01;
    COFD[765] = 3.51456245E+00;
    COFD[766] = -2.36917200E-01;
    COFD[767] = 1.00534193E-02;
    COFD[768] = -1.63733545E+01;
    COFD[769] = 3.79864092E+00;
    COFD[770] = -2.71062595E-01;
    COFD[771] = 1.14218887E-02;
    COFD[772] = -1.62164844E+01;
    COFD[773] = 3.71034643E+00;
    COFD[774] = -2.60739292E-01;
    COFD[775] = 1.10202921E-02;
    COFD[776] = -1.70275170E+01;
    COFD[777] = 3.97787104E+00;
    COFD[778] = -2.92958143E-01;
    COFD[779] = 1.23153094E-02;
    COFD[780] = -1.69092218E+01;
    COFD[781] = 3.94140805E+00;
    COFD[782] = -2.88511662E-01;
    COFD[783] = 1.21342675E-02;
    COFD[784] = -1.66034546E+01;
    COFD[785] = 3.75923792E+00;
    COFD[786] = -2.66511261E-01;
    COFD[787] = 1.12471370E-02;
    COFD[788] = -1.66068325E+01;
    COFD[789] = 3.75923792E+00;
    COFD[790] = -2.66511261E-01;
    COFD[791] = 1.12471370E-02;
    COFD[792] = -1.94888119E+01;
    COFD[793] = 4.64703793E+00;
    COFD[794] = -3.72283709E-01;
    COFD[795] = 1.54551501E-02;
    COFD[796] = -1.95587320E+01;
    COFD[797] = 4.64075109E+00;
    COFD[798] = -3.71554667E-01;
    COFD[799] = 1.54269261E-02;
    COFD[800] = -1.97432870E+01;
    COFD[801] = 4.68838474E+00;
    COFD[802] = -3.77193241E-01;
    COFD[803] = 1.56502297E-02;
    COFD[804] = -1.91036850E+01;
    COFD[805] = 4.51902651E+00;
    COFD[806] = -3.57827404E-01;
    COFD[807] = 1.49129882E-02;
    COFD[808] = -1.32446869E+01;
    COFD[809] = 2.69833584E+00;
    COFD[810] = -1.33541002E-01;
    COFD[811] = 5.68346517E-03;
    COFD[812] = -2.09653714E+01;
    COFD[813] = 4.99879868E+00;
    COFD[814] = -4.12878607E-01;
    COFD[815] = 1.70210595E-02;
    COFD[816] = -1.47230052E+01;
    COFD[817] = 3.10657376E+00;
    COFD[818] = -1.85922460E-01;
    COFD[819] = 7.92680827E-03;
    COFD[820] = -1.73476278E+01;
    COFD[821] = 3.98940337E+00;
    COFD[822] = -2.94368046E-01;
    COFD[823] = 1.23728707E-02;
    COFD[824] = -1.97265566E+01;
    COFD[825] = 5.00923653E+00;
    COFD[826] = -4.13897836E-01;
    COFD[827] = 1.70524216E-02;
    COFD[828] = -1.44684500E+01;
    COFD[829] = 3.00511224E+00;
    COFD[830] = -1.73057456E-01;
    COFD[831] = 7.38278721E-03;
    COFD[832] = -1.15547543E+01;
    COFD[833] = 2.41720750E+00;
    COFD[834] = -1.00182262E-01;
    COFD[835] = 4.37376763E-03;
    COFD[836] = -1.35078099E+01;
    COFD[837] = 2.80526779E+00;
    COFD[838] = -1.47524446E-01;
    COFD[839] = 6.29376912E-03;
    COFD[840] = -1.47385123E+01;
    COFD[841] = 3.10657376E+00;
    COFD[842] = -1.85922460E-01;
    COFD[843] = 7.92680827E-03;
    COFD[844] = -1.47309355E+01;
    COFD[845] = 3.10657376E+00;
    COFD[846] = -1.85922460E-01;
    COFD[847] = 7.92680827E-03;
    COFD[848] = -1.37337211E+01;
    COFD[849] = 3.42674072E+00;
    COFD[850] = -2.25944624E-01;
    COFD[851] = 9.59602387E-03;
    COFD[852] = -1.34872388E+01;
    COFD[853] = 2.80526779E+00;
    COFD[854] = -1.47524446E-01;
    COFD[855] = 6.29376912E-03;
    COFD[856] = -1.91134554E+01;
    COFD[857] = 4.63832627E+00;
    COFD[858] = -3.71274069E-01;
    COFD[859] = 1.54160891E-02;
    COFD[860] = -1.97923482E+01;
    COFD[861] = 4.88737199E+00;
    COFD[862] = -4.00855012E-01;
    COFD[863] = 1.65924977E-02;
    COFD[864] = -1.97832077E+01;
    COFD[865] = 4.88737199E+00;
    COFD[866] = -4.00855012E-01;
    COFD[867] = 1.65924977E-02;
    COFD[868] = -1.54192952E+01;
    COFD[869] = 3.41763870E+00;
    COFD[870] = -2.24801739E-01;
    COFD[871] = 9.54815836E-03;
    COFD[872] = -1.54436104E+01;
    COFD[873] = 3.41763870E+00;
    COFD[874] = -2.24801739E-01;
    COFD[875] = 9.54815836E-03;
    COFD[876] = -1.53851569E+01;
    COFD[877] = 3.39339822E+00;
    COFD[878] = -2.21757211E-01;
    COFD[879] = 9.42063333E-03;
    COFD[880] = -1.67826233E+01;
    COFD[881] = 3.79802425E+00;
    COFD[882] = -2.70990986E-01;
    COFD[883] = 1.14191222E-02;
    COFD[884] = -1.78213641E+01;
    COFD[885] = 4.17648998E+00;
    COFD[886] = -3.16975886E-01;
    COFD[887] = 1.32846969E-02;
    COFD[888] = -1.75147592E+01;
    COFD[889] = 4.03367299E+00;
    COFD[890] = -2.99744218E-01;
    COFD[891] = 1.25908330E-02;
    COFD[892] = -1.84720157E+01;
    COFD[893] = 4.34181535E+00;
    COFD[894] = -3.36587321E-01;
    COFD[895] = 1.40618044E-02;
    COFD[896] = -1.83678207E+01;
    COFD[897] = 4.30984197E+00;
    COFD[898] = -3.32770986E-01;
    COFD[899] = 1.39094747E-02;
    COFD[900] = -1.79906619E+01;
    COFD[901] = 4.11071307E+00;
    COFD[902] = -3.09074203E-01;
    COFD[903] = 1.29678777E-02;
    COFD[904] = -1.79959176E+01;
    COFD[905] = 4.11071307E+00;
    COFD[906] = -3.09074203E-01;
    COFD[907] = 1.29678777E-02;
    COFD[908] = -2.09640855E+01;
    COFD[909] = 4.99879868E+00;
    COFD[910] = -4.12878607E-01;
    COFD[911] = 1.70210595E-02;
    COFD[912] = -2.10482827E+01;
    COFD[913] = 4.99490919E+00;
    COFD[914] = -4.12490091E-01;
    COFD[915] = 1.70086418E-02;
    COFD[916] = -2.11872013E+01;
    COFD[917] = 5.02181116E+00;
    COFD[918] = -4.15179449E-01;
    COFD[919] = 1.70945609E-02;
    COFD[920] = -2.06177961E+01;
    COFD[921] = 4.88872783E+00;
    COFD[922] = -4.01009699E-01;
    COFD[923] = 1.65983883E-02;
    COFD[924] = -1.44426226E+01;
    COFD[925] = 2.99858376E+00;
    COFD[926] = -1.72232643E-01;
    COFD[927] = 7.34804765E-03;
    COFD[928] = -2.09540999E+01;
    COFD[929] = 4.99879868E+00;
    COFD[930] = -4.12878607E-01;
    COFD[931] = 1.70210595E-02;
    COFD[932] = -1.47156580E+01;
    COFD[933] = 3.10657376E+00;
    COFD[934] = -1.85922460E-01;
    COFD[935] = 7.92680827E-03;
    COFD[936] = -1.73390885E+01;
    COFD[937] = 3.98940337E+00;
    COFD[938] = -2.94368046E-01;
    COFD[939] = 1.23728707E-02;
    COFD[940] = -1.97212974E+01;
    COFD[941] = 5.00923653E+00;
    COFD[942] = -4.13897836E-01;
    COFD[943] = 1.70524216E-02;
    COFD[944] = -1.44616016E+01;
    COFD[945] = 3.00511224E+00;
    COFD[946] = -1.73057456E-01;
    COFD[947] = 7.38278721E-03;
    COFD[948] = -1.15539008E+01;
    COFD[949] = 2.41720750E+00;
    COFD[950] = -1.00182262E-01;
    COFD[951] = 4.37376763E-03;
    COFD[952] = -1.35027460E+01;
    COFD[953] = 2.80526779E+00;
    COFD[954] = -1.47524446E-01;
    COFD[955] = 6.29376912E-03;
    COFD[956] = -1.47309355E+01;
    COFD[957] = 3.10657376E+00;
    COFD[958] = -1.85922460E-01;
    COFD[959] = 7.92680827E-03;
    COFD[960] = -1.47234717E+01;
    COFD[961] = 3.10657376E+00;
    COFD[962] = -1.85922460E-01;
    COFD[963] = 7.92680827E-03;
    COFD[964] = -1.37332818E+01;
    COFD[965] = 3.42674072E+00;
    COFD[966] = -2.25944624E-01;
    COFD[967] = 9.59602387E-03;
    COFD[968] = -1.34823780E+01;
    COFD[969] = 2.80526779E+00;
    COFD[970] = -1.47524446E-01;
    COFD[971] = 6.29376912E-03;
    COFD[972] = -1.91062232E+01;
    COFD[973] = 4.63832627E+00;
    COFD[974] = -3.71274069E-01;
    COFD[975] = 1.54160891E-02;
    COFD[976] = -1.97852399E+01;
    COFD[977] = 4.88737199E+00;
    COFD[978] = -4.00855012E-01;
    COFD[979] = 1.65924977E-02;
    COFD[980] = -1.97762272E+01;
    COFD[981] = 4.88737199E+00;
    COFD[982] = -4.00855012E-01;
    COFD[983] = 1.65924977E-02;
    COFD[984] = -1.54148567E+01;
    COFD[985] = 3.41763870E+00;
    COFD[986] = -2.24801739E-01;
    COFD[987] = 9.54815836E-03;
    COFD[988] = -1.54389517E+01;
    COFD[989] = 3.41763870E+00;
    COFD[990] = -2.24801739E-01;
    COFD[991] = 9.54815836E-03;
    COFD[992] = -1.53802870E+01;
    COFD[993] = 3.39339822E+00;
    COFD[994] = -2.21757211E-01;
    COFD[995] = 9.42063333E-03;
    COFD[996] = -1.67759053E+01;
    COFD[997] = 3.79802425E+00;
    COFD[998] = -2.70990986E-01;
    COFD[999] = 1.14191222E-02;
    COFD[1000] = -1.78145099E+01;
    COFD[1001] = 4.17648998E+00;
    COFD[1002] = -3.16975886E-01;
    COFD[1003] = 1.32846969E-02;
    COFD[1004] = -1.75077731E+01;
    COFD[1005] = 4.03367299E+00;
    COFD[1006] = -2.99744218E-01;
    COFD[1007] = 1.25908330E-02;
    COFD[1008] = -1.84638250E+01;
    COFD[1009] = 4.34181535E+00;
    COFD[1010] = -3.36587321E-01;
    COFD[1011] = 1.40618044E-02;
    COFD[1012] = -1.83595375E+01;
    COFD[1013] = 4.30984197E+00;
    COFD[1014] = -3.32770986E-01;
    COFD[1015] = 1.39094747E-02;
    COFD[1016] = -1.79822885E+01;
    COFD[1017] = 4.11071307E+00;
    COFD[1018] = -3.09074203E-01;
    COFD[1019] = 1.29678777E-02;
    COFD[1020] = -1.79874565E+01;
    COFD[1021] = 4.11071307E+00;
    COFD[1022] = -3.09074203E-01;
    COFD[1023] = 1.29678777E-02;
    COFD[1024] = -2.09528426E+01;
    COFD[1025] = 4.99879868E+00;
    COFD[1026] = -4.12878607E-01;
    COFD[1027] = 1.70210595E-02;
    COFD[1028] = -2.10363020E+01;
    COFD[1029] = 4.99490919E+00;
    COFD[1030] = -4.12490091E-01;
    COFD[1031] = 1.70086418E-02;
    COFD[1032] = -2.11749653E+01;
    COFD[1033] = 5.02181116E+00;
    COFD[1034] = -4.15179449E-01;
    COFD[1035] = 1.70945609E-02;
    COFD[1036] = -2.06065545E+01;
    COFD[1037] = 4.88872783E+00;
    COFD[1038] = -4.01009699E-01;
    COFD[1039] = 1.65983883E-02;
    COFD[1040] = -1.44357738E+01;
    COFD[1041] = 2.99858376E+00;
    COFD[1042] = -1.72232643E-01;
    COFD[1043] = 7.34804765E-03;
    COFD[1044] = -1.99238711E+01;
    COFD[1045] = 5.27803851E+00;
    COFD[1046] = -4.42140176E-01;
    COFD[1047] = 1.80237084E-02;
    COFD[1048] = -1.37328153E+01;
    COFD[1049] = 3.42674072E+00;
    COFD[1050] = -2.25944624E-01;
    COFD[1051] = 9.59602387E-03;
    COFD[1052] = -1.64909711E+01;
    COFD[1053] = 4.36229440E+00;
    COFD[1054] = -3.39095464E-01;
    COFD[1055] = 1.41645916E-02;
    COFD[1056] = -1.88398861E+01;
    COFD[1057] = 5.28766000E+00;
    COFD[1058] = -4.43001094E-01;
    COFD[1059] = 1.80463788E-02;
    COFD[1060] = -1.34773721E+01;
    COFD[1061] = 3.30798941E+00;
    COFD[1062] = -2.10999668E-01;
    COFD[1063] = 8.96870768E-03;
    COFD[1064] = -1.05322133E+01;
    COFD[1065] = 2.38287287E+00;
    COFD[1066] = -9.25451677E-02;
    COFD[1067] = 3.90477283E-03;
    COFD[1068] = -1.25828466E+01;
    COFD[1069] = 3.11263114E+00;
    COFD[1070] = -1.86688345E-01;
    COFD[1071] = 7.95911073E-03;
    COFD[1072] = -1.37337211E+01;
    COFD[1073] = 3.42674072E+00;
    COFD[1074] = -2.25944624E-01;
    COFD[1075] = 9.59602387E-03;
    COFD[1076] = -1.37332818E+01;
    COFD[1077] = 3.42674072E+00;
    COFD[1078] = -2.25944624E-01;
    COFD[1079] = 9.59602387E-03;
    COFD[1080] = -1.37059729E+01;
    COFD[1081] = 3.74572091E+00;
    COFD[1082] = -2.64900591E-01;
    COFD[1083] = 1.11831118E-02;
    COFD[1084] = -1.25810873E+01;
    COFD[1085] = 3.11263114E+00;
    COFD[1086] = -1.86688345E-01;
    COFD[1087] = 7.95911073E-03;
    COFD[1088] = -1.83321659E+01;
    COFD[1089] = 4.99755571E+00;
    COFD[1090] = -4.12758959E-01;
    COFD[1091] = 1.70174618E-02;
    COFD[1092] = -1.89125228E+01;
    COFD[1093] = 5.19819765E+00;
    COFD[1094] = -4.34963808E-01;
    COFD[1095] = 1.78329868E-02;
    COFD[1096] = -1.89119590E+01;
    COFD[1097] = 5.19819765E+00;
    COFD[1098] = -4.34963808E-01;
    COFD[1099] = 1.78329868E-02;
    COFD[1100] = -1.47105599E+01;
    COFD[1101] = 3.73941906E+00;
    COFD[1102] = -2.64151516E-01;
    COFD[1103] = 1.11534528E-02;
    COFD[1104] = -1.47128122E+01;
    COFD[1105] = 3.73941906E+00;
    COFD[1106] = -2.64151516E-01;
    COFD[1107] = 1.11534528E-02;
    COFD[1108] = -1.46504290E+01;
    COFD[1109] = 3.72374230E+00;
    COFD[1110] = -2.62310679E-01;
    COFD[1111] = 1.10815796E-02;
    COFD[1112] = -1.60801669E+01;
    COFD[1113] = 4.18258317E+00;
    COFD[1114] = -3.17706089E-01;
    COFD[1115] = 1.33139059E-02;
    COFD[1116] = -1.70609605E+01;
    COFD[1117] = 4.53774238E+00;
    COFD[1118] = -3.59864268E-01;
    COFD[1119] = 1.49858355E-02;
    COFD[1120] = -1.68129119E+01;
    COFD[1121] = 4.41680518E+00;
    COFD[1122] = -3.45824229E-01;
    COFD[1123] = 1.44423400E-02;
    COFD[1124] = -1.76744140E+01;
    COFD[1125] = 4.70907386E+00;
    COFD[1126] = -3.79700042E-01;
    COFD[1127] = 1.57519279E-02;
    COFD[1128] = -1.75284656E+01;
    COFD[1129] = 4.66504845E+00;
    COFD[1130] = -3.74375842E-01;
    COFD[1131] = 1.55363348E-02;
    COFD[1132] = -1.72372988E+01;
    COFD[1133] = 4.49320169E+00;
    COFD[1134] = -3.55062447E-01;
    COFD[1135] = 1.48161449E-02;
    COFD[1136] = -1.72375725E+01;
    COFD[1137] = 4.49320169E+00;
    COFD[1138] = -3.55062447E-01;
    COFD[1139] = 1.48161449E-02;
    COFD[1140] = -1.99238205E+01;
    COFD[1141] = 5.27803851E+00;
    COFD[1142] = -4.42140176E-01;
    COFD[1143] = 1.80237084E-02;
    COFD[1144] = -1.99886318E+01;
    COFD[1145] = 5.27453026E+00;
    COFD[1146] = -4.41827327E-01;
    COFD[1147] = 1.80155336E-02;
    COFD[1148] = -2.01190603E+01;
    COFD[1149] = 5.29848259E+00;
    COFD[1150] = -4.44011383E-01;
    COFD[1151] = 1.80753134E-02;
    COFD[1152] = -1.96413638E+01;
    COFD[1153] = 5.19930736E+00;
    COFD[1154] = -4.35079915E-01;
    COFD[1155] = 1.78369576E-02;
    COFD[1156] = -1.34540506E+01;
    COFD[1157] = 3.30341833E+00;
    COFD[1158] = -2.10450242E-01;
    COFD[1159] = 8.94670534E-03;
    COFD[1160] = -1.94633191E+01;
    COFD[1161] = 4.64703793E+00;
    COFD[1162] = -3.72283709E-01;
    COFD[1163] = 1.54551501E-02;
    COFD[1164] = -1.34772622E+01;
    COFD[1165] = 2.80526779E+00;
    COFD[1166] = -1.47524446E-01;
    COFD[1167] = 6.29376912E-03;
    COFD[1168] = -1.60184872E+01;
    COFD[1169] = 3.68166087E+00;
    COFD[1170] = -2.57364559E-01;
    COFD[1171] = 1.08881828E-02;
    COFD[1172] = -1.82875708E+01;
    COFD[1173] = 4.66435737E+00;
    COFD[1174] = -3.74296088E-01;
    COFD[1175] = 1.55332609E-02;
    COFD[1176] = -1.32520179E+01;
    COFD[1177] = 2.70501014E+00;
    COFD[1178] = -1.34413475E-01;
    COFD[1179] = 5.72152293E-03;
    COFD[1180] = -1.09715176E+01;
    COFD[1181] = 2.31040527E+00;
    COFD[1182] = -8.74442763E-02;
    COFD[1183] = 3.87525190E-03;
    COFD[1184] = -1.22826905E+01;
    COFD[1185] = 2.49197560E+00;
    COFD[1186] = -1.06603364E-01;
    COFD[1187] = 4.50971238E-03;
    COFD[1188] = -1.34872388E+01;
    COFD[1189] = 2.80526779E+00;
    COFD[1190] = -1.47524446E-01;
    COFD[1191] = 6.29376912E-03;
    COFD[1192] = -1.34823780E+01;
    COFD[1193] = 2.80526779E+00;
    COFD[1194] = -1.47524446E-01;
    COFD[1195] = 6.29376912E-03;
    COFD[1196] = -1.25810873E+01;
    COFD[1197] = 3.11263114E+00;
    COFD[1198] = -1.86688345E-01;
    COFD[1199] = 7.95911073E-03;
    COFD[1200] = -1.22676499E+01;
    COFD[1201] = 2.49197560E+00;
    COFD[1202] = -1.06603364E-01;
    COFD[1203] = 4.50971238E-03;
    COFD[1204] = -1.76771446E+01;
    COFD[1205] = 4.29054603E+00;
    COFD[1206] = -3.30479889E-01;
    COFD[1207] = 1.38185279E-02;
    COFD[1208] = -1.82888352E+01;
    COFD[1209] = 4.51789863E+00;
    COFD[1210] = -3.57703548E-01;
    COFD[1211] = 1.49085039E-02;
    COFD[1212] = -1.82828341E+01;
    COFD[1213] = 4.51789863E+00;
    COFD[1214] = -3.57703548E-01;
    COFD[1215] = 1.49085039E-02;
    COFD[1216] = -1.42149314E+01;
    COFD[1217] = 3.10519183E+00;
    COFD[1218] = -1.85747616E-01;
    COFD[1219] = 7.91942899E-03;
    COFD[1220] = -1.42331195E+01;
    COFD[1221] = 3.10519183E+00;
    COFD[1222] = -1.85747616E-01;
    COFD[1223] = 7.91942899E-03;
    COFD[1224] = -1.41789556E+01;
    COFD[1225] = 3.08518563E+00;
    COFD[1226] = -1.83211692E-01;
    COFD[1227] = 7.81220123E-03;
    COFD[1228] = -1.55744378E+01;
    COFD[1229] = 3.51456245E+00;
    COFD[1230] = -2.36917200E-01;
    COFD[1231] = 1.00534193E-02;
    COFD[1232] = -1.63541181E+01;
    COFD[1233] = 3.79864092E+00;
    COFD[1234] = -2.71062595E-01;
    COFD[1235] = 1.14218887E-02;
    COFD[1236] = -1.61969977E+01;
    COFD[1237] = 3.71034643E+00;
    COFD[1238] = -2.60739292E-01;
    COFD[1239] = 1.10202921E-02;
    COFD[1240] = -1.70058787E+01;
    COFD[1241] = 3.97787104E+00;
    COFD[1242] = -2.92958143E-01;
    COFD[1243] = 1.23153094E-02;
    COFD[1244] = -1.68874275E+01;
    COFD[1245] = 3.94140805E+00;
    COFD[1246] = -2.88511662E-01;
    COFD[1247] = 1.21342675E-02;
    COFD[1248] = -1.65815097E+01;
    COFD[1249] = 3.75923792E+00;
    COFD[1250] = -2.66511261E-01;
    COFD[1251] = 1.12471370E-02;
    COFD[1252] = -1.65847421E+01;
    COFD[1253] = 3.75923792E+00;
    COFD[1254] = -2.66511261E-01;
    COFD[1255] = 1.12471370E-02;
    COFD[1256] = -1.94626201E+01;
    COFD[1257] = 4.64703793E+00;
    COFD[1258] = -3.72283709E-01;
    COFD[1259] = 1.54551501E-02;
    COFD[1260] = -1.95315965E+01;
    COFD[1261] = 4.64075109E+00;
    COFD[1262] = -3.71554667E-01;
    COFD[1263] = 1.54269261E-02;
    COFD[1264] = -1.97158370E+01;
    COFD[1265] = 4.68838474E+00;
    COFD[1266] = -3.77193241E-01;
    COFD[1267] = 1.56502297E-02;
    COFD[1268] = -1.90774948E+01;
    COFD[1269] = 4.51902651E+00;
    COFD[1270] = -3.57827404E-01;
    COFD[1271] = 1.49129882E-02;
    COFD[1272] = -1.32254608E+01;
    COFD[1273] = 2.69833584E+00;
    COFD[1274] = -1.33541002E-01;
    COFD[1275] = 5.68346517E-03;
    COFD[1276] = -2.22225767E+01;
    COFD[1277] = 5.05132223E+00;
    COFD[1278] = -3.76867047E-01;
    COFD[1279] = 1.39005636E-02;
    COFD[1280] = -1.94217512E+01;
    COFD[1281] = 4.74958634E+00;
    COFD[1282] = -3.84601412E-01;
    COFD[1283] = 1.59505248E-02;
    COFD[1284] = -2.17139201E+01;
    COFD[1285] = 5.38758317E+00;
    COFD[1286] = -4.48340933E-01;
    COFD[1287] = 1.79944740E-02;
    COFD[1288] = -1.93199869E+01;
    COFD[1289] = 4.21681909E+00;
    COFD[1290] = -2.50231374E-01;
    COFD[1291] = 7.81267943E-03;
    COFD[1292] = -1.90947006E+01;
    COFD[1293] = 4.62819535E+00;
    COFD[1294] = -3.70101251E-01;
    COFD[1295] = 1.53707619E-02;
    COFD[1296] = -1.48338675E+01;
    COFD[1297] = 3.57131810E+00;
    COFD[1298] = -2.43940431E-01;
    COFD[1299] = 1.03434646E-02;
    COFD[1300] = -1.76970903E+01;
    COFD[1301] = 4.29054603E+00;
    COFD[1302] = -3.30479889E-01;
    COFD[1303] = 1.38185279E-02;
    COFD[1304] = -1.91134554E+01;
    COFD[1305] = 4.63832627E+00;
    COFD[1306] = -3.71274069E-01;
    COFD[1307] = 1.54160891E-02;
    COFD[1308] = -1.91062232E+01;
    COFD[1309] = 4.63832627E+00;
    COFD[1310] = -3.71274069E-01;
    COFD[1311] = 1.54160891E-02;
    COFD[1312] = -1.83321659E+01;
    COFD[1313] = 4.99755571E+00;
    COFD[1314] = -4.12758959E-01;
    COFD[1315] = 1.70174618E-02;
    COFD[1316] = -1.76771446E+01;
    COFD[1317] = 4.29054603E+00;
    COFD[1318] = -3.30479889E-01;
    COFD[1319] = 1.38185279E-02;
    COFD[1320] = -2.10593567E+01;
    COFD[1321] = 4.93529071E+00;
    COFD[1322] = -3.65604768E-01;
    COFD[1323] = 1.35569433E-02;
    COFD[1324] = -2.16448504E+01;
    COFD[1325] = 5.18393156E+00;
    COFD[1326] = -3.98839089E-01;
    COFD[1327] = 1.50174959E-02;
    COFD[1328] = -2.16361003E+01;
    COFD[1329] = 5.18393156E+00;
    COFD[1330] = -3.98839089E-01;
    COFD[1331] = 1.50174959E-02;
    COFD[1332] = -1.99931992E+01;
    COFD[1333] = 4.99234147E+00;
    COFD[1334] = -4.12227413E-01;
    COFD[1335] = 1.69999326E-02;
    COFD[1336] = -2.00168352E+01;
    COFD[1337] = 4.99234147E+00;
    COFD[1338] = -4.12227413E-01;
    COFD[1339] = 1.69999326E-02;
    COFD[1340] = -2.02417979E+01;
    COFD[1341] = 5.05995682E+00;
    COFD[1342] = -4.19583868E-01;
    COFD[1343] = 1.72643549E-02;
    COFD[1344] = -2.13316826E+01;
    COFD[1345] = 5.31504777E+00;
    COFD[1346] = -4.45287609E-01;
    COFD[1347] = 1.80982140E-02;
    COFD[1348] = -2.18246157E+01;
    COFD[1349] = 5.42978824E+00;
    COFD[1350] = -4.51570216E-01;
    COFD[1351] = 1.80515349E-02;
    COFD[1352] = -2.17177863E+01;
    COFD[1353] = 5.38391804E+00;
    COFD[1354] = -4.48822101E-01;
    COFD[1355] = 1.80505757E-02;
    COFD[1356] = -2.21211461E+01;
    COFD[1357] = 5.43170153E+00;
    COFD[1358] = -4.46963936E-01;
    COFD[1359] = 1.76740849E-02;
    COFD[1360] = -2.21010249E+01;
    COFD[1361] = 5.43812778E+00;
    COFD[1362] = -4.48750136E-01;
    COFD[1363] = 1.77861511E-02;
    COFD[1364] = -2.20454826E+01;
    COFD[1365] = 5.39682182E+00;
    COFD[1366] = -4.48776824E-01;
    COFD[1367] = 1.79851810E-02;
    COFD[1368] = -2.20504720E+01;
    COFD[1369] = 5.39682182E+00;
    COFD[1370] = -4.48776824E-01;
    COFD[1371] = 1.79851810E-02;
    COFD[1372] = -2.22213768E+01;
    COFD[1373] = 5.05132223E+00;
    COFD[1374] = -3.76867047E-01;
    COFD[1375] = 1.39005636E-02;
    COFD[1376] = -2.19824618E+01;
    COFD[1377] = 4.91387652E+00;
    COFD[1378] = -3.58234874E-01;
    COFD[1379] = 1.30693081E-02;
    COFD[1380] = -2.18763645E+01;
    COFD[1381] = 4.83720845E+00;
    COFD[1382] = -3.46627541E-01;
    COFD[1383] = 1.25098849E-02;
    COFD[1384] = -2.19972394E+01;
    COFD[1385] = 4.99433306E+00;
    COFD[1386] = -3.72400138E-01;
    COFD[1387] = 1.38127710E-02;
    COFD[1388] = -1.90471463E+01;
    COFD[1389] = 4.61496749E+00;
    COFD[1390] = -3.68575530E-01;
    COFD[1391] = 1.53120471E-02;
    COFD[1392] = -2.16117599E+01;
    COFD[1393] = 4.74703304E+00;
    COFD[1394] = -3.29522725E-01;
    COFD[1395] = 1.15765554E-02;
    COFD[1396] = -1.97777929E+01;
    COFD[1397] = 4.88737199E+00;
    COFD[1398] = -4.00855012E-01;
    COFD[1399] = 1.65924977E-02;
    COFD[1400] = -2.19036727E+01;
    COFD[1401] = 5.44375518E+00;
    COFD[1402] = -4.52244836E-01;
    COFD[1403] = 1.80389212E-02;
    COFD[1404] = -2.03286949E+01;
    COFD[1405] = 4.72161925E+00;
    COFD[1406] = -3.25604803E-01;
    COFD[1407] = 1.13856236E-02;
    COFD[1408] = -1.94491481E+01;
    COFD[1409] = 4.76713321E+00;
    COFD[1410] = -3.86705929E-01;
    COFD[1411] = 1.60350491E-02;
    COFD[1412] = -1.50226829E+01;
    COFD[1413] = 3.65330206E+00;
    COFD[1414] = -2.53985209E-01;
    COFD[1415] = 1.07540119E-02;
    COFD[1416] = -1.83085512E+01;
    COFD[1417] = 4.51789863E+00;
    COFD[1418] = -3.57703548E-01;
    COFD[1419] = 1.49085039E-02;
    COFD[1420] = -1.97923482E+01;
    COFD[1421] = 4.88737199E+00;
    COFD[1422] = -4.00855012E-01;
    COFD[1423] = 1.65924977E-02;
    COFD[1424] = -1.97852399E+01;
    COFD[1425] = 4.88737199E+00;
    COFD[1426] = -4.00855012E-01;
    COFD[1427] = 1.65924977E-02;
    COFD[1428] = -1.89125228E+01;
    COFD[1429] = 5.19819765E+00;
    COFD[1430] = -4.34963808E-01;
    COFD[1431] = 1.78329868E-02;
    COFD[1432] = -1.82888352E+01;
    COFD[1433] = 4.51789863E+00;
    COFD[1434] = -3.57703548E-01;
    COFD[1435] = 1.49085039E-02;
    COFD[1436] = -2.16448504E+01;
    COFD[1437] = 5.18393156E+00;
    COFD[1438] = -3.98839089E-01;
    COFD[1439] = 1.50174959E-02;
    COFD[1440] = -2.11948780E+01;
    COFD[1441] = 4.94957962E+00;
    COFD[1442] = -3.61249743E-01;
    COFD[1443] = 1.31396501E-02;
    COFD[1444] = -2.11862687E+01;
    COFD[1445] = 4.94957962E+00;
    COFD[1446] = -3.61249743E-01;
    COFD[1447] = 1.31396501E-02;
    COFD[1448] = -2.05736708E+01;
    COFD[1449] = 5.19175702E+00;
    COFD[1450] = -4.34287561E-01;
    COFD[1451] = 1.78097455E-02;
    COFD[1452] = -2.05970566E+01;
    COFD[1453] = 5.19175702E+00;
    COFD[1454] = -4.34287561E-01;
    COFD[1455] = 1.78097455E-02;
    COFD[1456] = -2.05462669E+01;
    COFD[1457] = 5.17356357E+00;
    COFD[1458] = -4.32338093E-01;
    COFD[1459] = 1.77409309E-02;
    COFD[1460] = -2.16008317E+01;
    COFD[1461] = 5.38202048E+00;
    COFD[1462] = -4.48902904E-01;
    COFD[1463] = 1.80661615E-02;
    COFD[1464] = -2.19239746E+01;
    COFD[1465] = 5.42213817E+00;
    COFD[1466] = -4.44597702E-01;
    COFD[1467] = 1.75307003E-02;
    COFD[1468] = -2.19859525E+01;
    COFD[1469] = 5.44658055E+00;
    COFD[1470] = -4.51507914E-01;
    COFD[1471] = 1.79664843E-02;
    COFD[1472] = -2.21293977E+01;
    COFD[1473] = 5.38638208E+00;
    COFD[1474] = -4.34919926E-01;
    COFD[1475] = 1.69301831E-02;
    COFD[1476] = -2.21276873E+01;
    COFD[1477] = 5.40064912E+00;
    COFD[1478] = -4.37771126E-01;
    COFD[1479] = 1.70895608E-02;
    COFD[1480] = -2.22599437E+01;
    COFD[1481] = 5.43626968E+00;
    COFD[1482] = -4.48222611E-01;
    COFD[1483] = 1.77528565E-02;
    COFD[1484] = -2.22648380E+01;
    COFD[1485] = 5.43626968E+00;
    COFD[1486] = -4.48222611E-01;
    COFD[1487] = 1.77528565E-02;
    COFD[1488] = -2.16105899E+01;
    COFD[1489] = 4.74703304E+00;
    COFD[1490] = -3.29522725E-01;
    COFD[1491] = 1.15765554E-02;
    COFD[1492] = -2.17219099E+01;
    COFD[1493] = 4.75637248E+00;
    COFD[1494] = -3.30964596E-01;
    COFD[1495] = 1.16469146E-02;
    COFD[1496] = -2.16433785E+01;
    COFD[1497] = 4.69172346E+00;
    COFD[1498] = -3.21053636E-01;
    COFD[1499] = 1.11654521E-02;
    COFD[1500] = -2.19851974E+01;
    COFD[1501] = 4.94717766E+00;
    COFD[1502] = -3.60880803E-01;
    COFD[1503] = 1.31216550E-02;
    COFD[1504] = -1.94158077E+01;
    COFD[1505] = 4.75821765E+00;
    COFD[1506] = -3.85636294E-01;
    COFD[1507] = 1.59920740E-02;
    COFD[1508] = -2.15985720E+01;
    COFD[1509] = 4.74703304E+00;
    COFD[1510] = -3.29522725E-01;
    COFD[1511] = 1.15765554E-02;
    COFD[1512] = -1.97689123E+01;
    COFD[1513] = 4.88737199E+00;
    COFD[1514] = -4.00855012E-01;
    COFD[1515] = 1.65924977E-02;
    COFD[1516] = -2.18934539E+01;
    COFD[1517] = 5.44375518E+00;
    COFD[1518] = -4.52244836E-01;
    COFD[1519] = 1.80389212E-02;
    COFD[1520] = -2.03222242E+01;
    COFD[1521] = 4.72161925E+00;
    COFD[1522] = -3.25604803E-01;
    COFD[1523] = 1.13856236E-02;
    COFD[1524] = -1.94408353E+01;
    COFD[1525] = 4.76713321E+00;
    COFD[1526] = -3.86705929E-01;
    COFD[1527] = 1.60350491E-02;
    COFD[1528] = -1.50215914E+01;
    COFD[1529] = 3.65330206E+00;
    COFD[1530] = -2.53985209E-01;
    COFD[1531] = 1.07540119E-02;
    COFD[1532] = -1.83023101E+01;
    COFD[1533] = 4.51789863E+00;
    COFD[1534] = -3.57703548E-01;
    COFD[1535] = 1.49085039E-02;
    COFD[1536] = -1.97832077E+01;
    COFD[1537] = 4.88737199E+00;
    COFD[1538] = -4.00855012E-01;
    COFD[1539] = 1.65924977E-02;
    COFD[1540] = -1.97762272E+01;
    COFD[1541] = 4.88737199E+00;
    COFD[1542] = -4.00855012E-01;
    COFD[1543] = 1.65924977E-02;
    COFD[1544] = -1.89119590E+01;
    COFD[1545] = 5.19819765E+00;
    COFD[1546] = -4.34963808E-01;
    COFD[1547] = 1.78329868E-02;
    COFD[1548] = -1.82828341E+01;
    COFD[1549] = 4.51789863E+00;
    COFD[1550] = -3.57703548E-01;
    COFD[1551] = 1.49085039E-02;
    COFD[1552] = -2.16361003E+01;
    COFD[1553] = 5.18393156E+00;
    COFD[1554] = -3.98839089E-01;
    COFD[1555] = 1.50174959E-02;
    COFD[1556] = -2.11862687E+01;
    COFD[1557] = 4.94957962E+00;
    COFD[1558] = -3.61249743E-01;
    COFD[1559] = 1.31396501E-02;
    COFD[1560] = -2.11778052E+01;
    COFD[1561] = 4.94957962E+00;
    COFD[1562] = -3.61249743E-01;
    COFD[1563] = 1.31396501E-02;
    COFD[1564] = -2.05681711E+01;
    COFD[1565] = 5.19175702E+00;
    COFD[1566] = -4.34287561E-01;
    COFD[1567] = 1.78097455E-02;
    COFD[1568] = -2.05912951E+01;
    COFD[1569] = 5.19175702E+00;
    COFD[1570] = -4.34287561E-01;
    COFD[1571] = 1.78097455E-02;
    COFD[1572] = -2.05402551E+01;
    COFD[1573] = 5.17356357E+00;
    COFD[1574] = -4.32338093E-01;
    COFD[1575] = 1.77409309E-02;
    COFD[1576] = -2.15926683E+01;
    COFD[1577] = 5.38202048E+00;
    COFD[1578] = -4.48902904E-01;
    COFD[1579] = 1.80661615E-02;
    COFD[1580] = -2.19156553E+01;
    COFD[1581] = 5.42213817E+00;
    COFD[1582] = -4.44597702E-01;
    COFD[1583] = 1.75307003E-02;
    COFD[1584] = -2.19774825E+01;
    COFD[1585] = 5.44658055E+00;
    COFD[1586] = -4.51507914E-01;
    COFD[1587] = 1.79664843E-02;
    COFD[1588] = -2.21195673E+01;
    COFD[1589] = 5.38638208E+00;
    COFD[1590] = -4.34919926E-01;
    COFD[1591] = 1.69301831E-02;
    COFD[1592] = -2.21177535E+01;
    COFD[1593] = 5.40064912E+00;
    COFD[1594] = -4.37771126E-01;
    COFD[1595] = 1.70895608E-02;
    COFD[1596] = -2.22499095E+01;
    COFD[1597] = 5.43626968E+00;
    COFD[1598] = -4.48222611E-01;
    COFD[1599] = 1.77528565E-02;
    COFD[1600] = -2.22547061E+01;
    COFD[1601] = 5.43626968E+00;
    COFD[1602] = -4.48222611E-01;
    COFD[1603] = 1.77528565E-02;
    COFD[1604] = -2.15974324E+01;
    COFD[1605] = 4.74703304E+00;
    COFD[1606] = -3.29522725E-01;
    COFD[1607] = 1.15765554E-02;
    COFD[1608] = -2.17079728E+01;
    COFD[1609] = 4.75637248E+00;
    COFD[1610] = -3.30964596E-01;
    COFD[1611] = 1.16469146E-02;
    COFD[1612] = -2.16291739E+01;
    COFD[1613] = 4.69172346E+00;
    COFD[1614] = -3.21053636E-01;
    COFD[1615] = 1.11654521E-02;
    COFD[1616] = -2.19720412E+01;
    COFD[1617] = 4.94717766E+00;
    COFD[1618] = -3.60880803E-01;
    COFD[1619] = 1.31216550E-02;
    COFD[1620] = -1.94074945E+01;
    COFD[1621] = 4.75821765E+00;
    COFD[1622] = -3.85636294E-01;
    COFD[1623] = 1.59920740E-02;
    COFD[1624] = -2.15568608E+01;
    COFD[1625] = 5.27353202E+00;
    COFD[1626] = -4.41738406E-01;
    COFD[1627] = 1.80132159E-02;
    COFD[1628] = -1.54101814E+01;
    COFD[1629] = 3.41763870E+00;
    COFD[1630] = -2.24801739E-01;
    COFD[1631] = 9.54815836E-03;
    COFD[1632] = -1.81750360E+01;
    COFD[1633] = 4.35070874E+00;
    COFD[1634] = -3.37658377E-01;
    COFD[1635] = 1.41049661E-02;
    COFD[1636] = -2.05176849E+01;
    COFD[1637] = 5.28313423E+00;
    COFD[1638] = -4.42595579E-01;
    COFD[1639] = 1.80356675E-02;
    COFD[1640] = -1.51284988E+01;
    COFD[1641] = 3.30256227E+00;
    COFD[1642] = -2.10347339E-01;
    COFD[1643] = 8.94258391E-03;
    COFD[1644] = -1.16100827E+01;
    COFD[1645] = 2.38028423E+00;
    COFD[1646] = -9.22299290E-02;
    COFD[1647] = 3.89194284E-03;
    COFD[1648] = -1.42289701E+01;
    COFD[1649] = 3.10519183E+00;
    COFD[1650] = -1.85747616E-01;
    COFD[1651] = 7.91942899E-03;
    COFD[1652] = -1.54192952E+01;
    COFD[1653] = 3.41763870E+00;
    COFD[1654] = -2.24801739E-01;
    COFD[1655] = 9.54815836E-03;
    COFD[1656] = -1.54148567E+01;
    COFD[1657] = 3.41763870E+00;
    COFD[1658] = -2.24801739E-01;
    COFD[1659] = 9.54815836E-03;
    COFD[1660] = -1.47105599E+01;
    COFD[1661] = 3.73941906E+00;
    COFD[1662] = -2.64151516E-01;
    COFD[1663] = 1.11534528E-02;
    COFD[1664] = -1.42149314E+01;
    COFD[1665] = 3.10519183E+00;
    COFD[1666] = -1.85747616E-01;
    COFD[1667] = 7.91942899E-03;
    COFD[1668] = -1.99931992E+01;
    COFD[1669] = 4.99234147E+00;
    COFD[1670] = -4.12227413E-01;
    COFD[1671] = 1.69999326E-02;
    COFD[1672] = -2.05736708E+01;
    COFD[1673] = 5.19175702E+00;
    COFD[1674] = -4.34287561E-01;
    COFD[1675] = 1.78097455E-02;
    COFD[1676] = -2.05681711E+01;
    COFD[1677] = 5.19175702E+00;
    COFD[1678] = -4.34287561E-01;
    COFD[1679] = 1.78097455E-02;
    COFD[1680] = -1.62211837E+01;
    COFD[1681] = 3.73344904E+00;
    COFD[1682] = -2.63450213E-01;
    COFD[1683] = 1.11260598E-02;
    COFD[1684] = -1.62382314E+01;
    COFD[1685] = 3.73344904E+00;
    COFD[1686] = -2.63450213E-01;
    COFD[1687] = 1.11260598E-02;
    COFD[1688] = -1.61936541E+01;
    COFD[1689] = 3.71781439E+00;
    COFD[1690] = -2.61615179E-01;
    COFD[1691] = 1.10544512E-02;
    COFD[1692] = -1.76822324E+01;
    COFD[1693] = 4.17392086E+00;
    COFD[1694] = -3.16668286E-01;
    COFD[1695] = 1.32724062E-02;
    COFD[1696] = -1.86827277E+01;
    COFD[1697] = 4.53129263E+00;
    COFD[1698] = -3.59157354E-01;
    COFD[1699] = 1.49603162E-02;
    COFD[1700] = -1.84057376E+01;
    COFD[1701] = 4.40577261E+00;
    COFD[1702] = -3.44465891E-01;
    COFD[1703] = 1.43864134E-02;
    COFD[1704] = -1.93089486E+01;
    COFD[1705] = 4.69775512E+00;
    COFD[1706] = -3.78328786E-01;
    COFD[1707] = 1.56963078E-02;
    COFD[1708] = -1.91790889E+01;
    COFD[1709] = 4.65710841E+00;
    COFD[1710] = -3.73456882E-01;
    COFD[1711] = 1.55008064E-02;
    COFD[1712] = -1.88470808E+01;
    COFD[1713] = 4.48692836E+00;
    COFD[1714] = -3.54362855E-01;
    COFD[1715] = 1.47902773E-02;
    COFD[1716] = -1.88500135E+01;
    COFD[1717] = 4.48692836E+00;
    COFD[1718] = -3.54362855E-01;
    COFD[1719] = 1.47902773E-02;
    COFD[1720] = -2.15562373E+01;
    COFD[1721] = 5.27353202E+00;
    COFD[1722] = -4.41738406E-01;
    COFD[1723] = 1.80132159E-02;
    COFD[1724] = -2.16217667E+01;
    COFD[1725] = 5.27002901E+00;
    COFD[1726] = -4.41426690E-01;
    COFD[1727] = 1.80051108E-02;
    COFD[1728] = -2.17475829E+01;
    COFD[1729] = 5.29328698E+00;
    COFD[1730] = -4.43508195E-01;
    COFD[1731] = 1.80599224E-02;
    COFD[1732] = -2.12691260E+01;
    COFD[1733] = 5.19287816E+00;
    COFD[1734] = -4.34405339E-01;
    COFD[1735] = 1.78137969E-02;
    COFD[1736] = -1.51077105E+01;
    COFD[1737] = 3.29803842E+00;
    COFD[1738] = -2.09804181E-01;
    COFD[1739] = 8.92085522E-03;
    COFD[1740] = -2.15871654E+01;
    COFD[1741] = 5.27353202E+00;
    COFD[1742] = -4.41738406E-01;
    COFD[1743] = 1.80132159E-02;
    COFD[1744] = -1.54340467E+01;
    COFD[1745] = 3.41763870E+00;
    COFD[1746] = -2.24801739E-01;
    COFD[1747] = 9.54815836E-03;
    COFD[1748] = -1.82011239E+01;
    COFD[1749] = 4.35070874E+00;
    COFD[1750] = -3.37658377E-01;
    COFD[1751] = 1.41049661E-02;
    COFD[1752] = -2.05368957E+01;
    COFD[1753] = 5.28313423E+00;
    COFD[1754] = -4.42595579E-01;
    COFD[1755] = 1.80356675E-02;
    COFD[1756] = -1.51513486E+01;
    COFD[1757] = 3.30256227E+00;
    COFD[1758] = -2.10347339E-01;
    COFD[1759] = 8.94258391E-03;
    COFD[1760] = -1.16143127E+01;
    COFD[1761] = 2.38028423E+00;
    COFD[1762] = -9.22299290E-02;
    COFD[1763] = 3.89194284E-03;
    COFD[1764] = -1.42476859E+01;
    COFD[1765] = 3.10519183E+00;
    COFD[1766] = -1.85747616E-01;
    COFD[1767] = 7.91942899E-03;
    COFD[1768] = -1.54436104E+01;
    COFD[1769] = 3.41763870E+00;
    COFD[1770] = -2.24801739E-01;
    COFD[1771] = 9.54815836E-03;
    COFD[1772] = -1.54389517E+01;
    COFD[1773] = 3.41763870E+00;
    COFD[1774] = -2.24801739E-01;
    COFD[1775] = 9.54815836E-03;
    COFD[1776] = -1.47128122E+01;
    COFD[1777] = 3.73941906E+00;
    COFD[1778] = -2.64151516E-01;
    COFD[1779] = 1.11534528E-02;
    COFD[1780] = -1.42331195E+01;
    COFD[1781] = 3.10519183E+00;
    COFD[1782] = -1.85747616E-01;
    COFD[1783] = 7.91942899E-03;
    COFD[1784] = -2.00168352E+01;
    COFD[1785] = 4.99234147E+00;
    COFD[1786] = -4.12227413E-01;
    COFD[1787] = 1.69999326E-02;
    COFD[1788] = -2.05970566E+01;
    COFD[1789] = 5.19175702E+00;
    COFD[1790] = -4.34287561E-01;
    COFD[1791] = 1.78097455E-02;
    COFD[1792] = -2.05912951E+01;
    COFD[1793] = 5.19175702E+00;
    COFD[1794] = -4.34287561E-01;
    COFD[1795] = 1.78097455E-02;
    COFD[1796] = -1.62382314E+01;
    COFD[1797] = 3.73344904E+00;
    COFD[1798] = -2.63450213E-01;
    COFD[1799] = 1.11260598E-02;
    COFD[1800] = -1.62558809E+01;
    COFD[1801] = 3.73344904E+00;
    COFD[1802] = -2.63450213E-01;
    COFD[1803] = 1.11260598E-02;
    COFD[1804] = -1.62118658E+01;
    COFD[1805] = 3.71781439E+00;
    COFD[1806] = -2.61615179E-01;
    COFD[1807] = 1.10544512E-02;
    COFD[1808] = -1.77048073E+01;
    COFD[1809] = 4.17392086E+00;
    COFD[1810] = -3.16668286E-01;
    COFD[1811] = 1.32724062E-02;
    COFD[1812] = -1.87055896E+01;
    COFD[1813] = 4.53129263E+00;
    COFD[1814] = -3.59157354E-01;
    COFD[1815] = 1.49603162E-02;
    COFD[1816] = -1.84288731E+01;
    COFD[1817] = 4.40577261E+00;
    COFD[1818] = -3.44465891E-01;
    COFD[1819] = 1.43864134E-02;
    COFD[1820] = -1.93344144E+01;
    COFD[1821] = 4.69775512E+00;
    COFD[1822] = -3.78328786E-01;
    COFD[1823] = 1.56963078E-02;
    COFD[1824] = -1.92047221E+01;
    COFD[1825] = 4.65710841E+00;
    COFD[1826] = -3.73456882E-01;
    COFD[1827] = 1.55008064E-02;
    COFD[1828] = -1.88728753E+01;
    COFD[1829] = 4.48692836E+00;
    COFD[1830] = -3.54362855E-01;
    COFD[1831] = 1.47902773E-02;
    COFD[1832] = -1.88759638E+01;
    COFD[1833] = 4.48692836E+00;
    COFD[1834] = -3.54362855E-01;
    COFD[1835] = 1.47902773E-02;
    COFD[1836] = -2.15865030E+01;
    COFD[1837] = 5.27353202E+00;
    COFD[1838] = -4.41738406E-01;
    COFD[1839] = 1.80132159E-02;
    COFD[1840] = -2.16530054E+01;
    COFD[1841] = 5.27002901E+00;
    COFD[1842] = -4.41426690E-01;
    COFD[1843] = 1.80051108E-02;
    COFD[1844] = -2.17791442E+01;
    COFD[1845] = 5.29328698E+00;
    COFD[1846] = -4.43508195E-01;
    COFD[1847] = 1.80599224E-02;
    COFD[1848] = -2.12993900E+01;
    COFD[1849] = 5.19287816E+00;
    COFD[1850] = -4.34405339E-01;
    COFD[1851] = 1.78137969E-02;
    COFD[1852] = -1.51305610E+01;
    COFD[1853] = 3.29803842E+00;
    COFD[1854] = -2.09804181E-01;
    COFD[1855] = 8.92085522E-03;
    COFD[1856] = -2.15620449E+01;
    COFD[1857] = 5.26193907E+00;
    COFD[1858] = -4.40716330E-01;
    COFD[1859] = 1.79871592E-02;
    COFD[1860] = -1.53751619E+01;
    COFD[1861] = 3.39339822E+00;
    COFD[1862] = -2.21757211E-01;
    COFD[1863] = 9.42063333E-03;
    COFD[1864] = -1.81478704E+01;
    COFD[1865] = 4.32863219E+00;
    COFD[1866] = -3.35006470E-01;
    COFD[1867] = 1.39983891E-02;
    COFD[1868] = -2.09577528E+01;
    COFD[1869] = 5.37851020E+00;
    COFD[1870] = -4.48937706E-01;
    COFD[1871] = 1.80857999E-02;
    COFD[1872] = -1.51167688E+01;
    COFD[1873] = 3.28831215E+00;
    COFD[1874] = -2.08633772E-01;
    COFD[1875] = 8.87391434E-03;
    COFD[1876] = -1.15834311E+01;
    COFD[1877] = 2.37417315E+00;
    COFD[1878] = -9.15013712E-02;
    COFD[1879] = 3.86299275E-03;
    COFD[1880] = -1.41940170E+01;
    COFD[1881] = 3.08518563E+00;
    COFD[1882] = -1.83211692E-01;
    COFD[1883] = 7.81220123E-03;
    COFD[1884] = -1.53851569E+01;
    COFD[1885] = 3.39339822E+00;
    COFD[1886] = -2.21757211E-01;
    COFD[1887] = 9.42063333E-03;
    COFD[1888] = -1.53802870E+01;
    COFD[1889] = 3.39339822E+00;
    COFD[1890] = -2.21757211E-01;
    COFD[1891] = 9.42063333E-03;
    COFD[1892] = -1.46504290E+01;
    COFD[1893] = 3.72374230E+00;
    COFD[1894] = -2.62310679E-01;
    COFD[1895] = 1.10815796E-02;
    COFD[1896] = -1.41789556E+01;
    COFD[1897] = 3.08518563E+00;
    COFD[1898] = -1.83211692E-01;
    COFD[1899] = 7.81220123E-03;
    COFD[1900] = -2.02417979E+01;
    COFD[1901] = 5.05995682E+00;
    COFD[1902] = -4.19583868E-01;
    COFD[1903] = 1.72643549E-02;
    COFD[1904] = -2.05462669E+01;
    COFD[1905] = 5.17356357E+00;
    COFD[1906] = -4.32338093E-01;
    COFD[1907] = 1.77409309E-02;
    COFD[1908] = -2.05402551E+01;
    COFD[1909] = 5.17356357E+00;
    COFD[1910] = -4.32338093E-01;
    COFD[1911] = 1.77409309E-02;
    COFD[1912] = -1.61936541E+01;
    COFD[1913] = 3.71781439E+00;
    COFD[1914] = -2.61615179E-01;
    COFD[1915] = 1.10544512E-02;
    COFD[1916] = -1.62118658E+01;
    COFD[1917] = 3.71781439E+00;
    COFD[1918] = -2.61615179E-01;
    COFD[1919] = 1.10544512E-02;
    COFD[1920] = -1.61688299E+01;
    COFD[1921] = 3.70240968E+00;
    COFD[1922] = -2.59811278E-01;
    COFD[1923] = 1.09842260E-02;
    COFD[1924] = -1.76468178E+01;
    COFD[1925] = 4.15115465E+00;
    COFD[1926] = -3.13942179E-01;
    COFD[1927] = 1.31634405E-02;
    COFD[1928] = -1.86611512E+01;
    COFD[1929] = 4.51482946E+00;
    COFD[1930] = -3.57370004E-01;
    COFD[1931] = 1.48965915E-02;
    COFD[1932] = -1.83550733E+01;
    COFD[1933] = 4.37621420E+00;
    COFD[1934] = -3.40819809E-01;
    COFD[1935] = 1.42360254E-02;
    COFD[1936] = -1.92602241E+01;
    COFD[1937] = 4.66778124E+00;
    COFD[1938] = -3.74693124E-01;
    COFD[1939] = 1.55486448E-02;
    COFD[1940] = -1.91516731E+01;
    COFD[1941] = 4.63565403E+00;
    COFD[1942] = -3.70964470E-01;
    COFD[1943] = 1.54041133E-02;
    COFD[1944] = -1.88205264E+01;
    COFD[1945] = 4.46474713E+00;
    COFD[1946] = -3.51718970E-01;
    COFD[1947] = 1.46847779E-02;
    COFD[1948] = -1.88237651E+01;
    COFD[1949] = 4.46474713E+00;
    COFD[1950] = -3.51718970E-01;
    COFD[1951] = 1.46847779E-02;
    COFD[1952] = -2.15613443E+01;
    COFD[1953] = 5.26193907E+00;
    COFD[1954] = -4.40716330E-01;
    COFD[1955] = 1.79871592E-02;
    COFD[1956] = -2.16639108E+01;
    COFD[1957] = 5.26764804E+00;
    COFD[1958] = -4.41215097E-01;
    COFD[1959] = 1.79996259E-02;
    COFD[1960] = -2.17991940E+01;
    COFD[1961] = 5.29320616E+00;
    COFD[1962] = -4.43500540E-01;
    COFD[1963] = 1.80596976E-02;
    COFD[1964] = -2.13381291E+01;
    COFD[1965] = 5.19899042E+00;
    COFD[1966] = -4.35046771E-01;
    COFD[1967] = 1.78358248E-02;
    COFD[1968] = -1.50951030E+01;
    COFD[1969] = 3.28342843E+00;
    COFD[1970] = -2.08038899E-01;
    COFD[1971] = 8.84974116E-03;
    COFD[1972] = -2.25748277E+01;
    COFD[1973] = 5.43256416E+00;
    COFD[1974] = -4.51725369E-01;
    COFD[1975] = 1.80507219E-02;
    COFD[1976] = -1.67688618E+01;
    COFD[1977] = 3.79802425E+00;
    COFD[1978] = -2.70990986E-01;
    COFD[1979] = 1.14191222E-02;
    COFD[1980] = -1.98631073E+01;
    COFD[1981] = 4.82865226E+00;
    COFD[1982] = -3.94029099E-01;
    COFD[1983] = 1.63270142E-02;
    COFD[1984] = -2.14416154E+01;
    COFD[1985] = 5.43806985E+00;
    COFD[1986] = -4.51984607E-01;
    COFD[1987] = 1.80451889E-02;
    COFD[1988] = -1.65688462E+01;
    COFD[1989] = 3.72373814E+00;
    COFD[1990] = -2.62310192E-01;
    COFD[1991] = 1.10815606E-02;
    COFD[1992] = -1.26266190E+01;
    COFD[1993] = 2.71852437E+00;
    COFD[1994] = -1.36179939E-01;
    COFD[1995] = 5.79856970E-03;
    COFD[1996] = -1.55934124E+01;
    COFD[1997] = 3.51456245E+00;
    COFD[1998] = -2.36917200E-01;
    COFD[1999] = 1.00534193E-02;
    COFD[2000] = -1.67826233E+01;
    COFD[2001] = 3.79802425E+00;
    COFD[2002] = -2.70990986E-01;
    COFD[2003] = 1.14191222E-02;
    COFD[2004] = -1.67759053E+01;
    COFD[2005] = 3.79802425E+00;
    COFD[2006] = -2.70990986E-01;
    COFD[2007] = 1.14191222E-02;
    COFD[2008] = -1.60801669E+01;
    COFD[2009] = 4.18258317E+00;
    COFD[2010] = -3.17706089E-01;
    COFD[2011] = 1.33139059E-02;
    COFD[2012] = -1.55744378E+01;
    COFD[2013] = 3.51456245E+00;
    COFD[2014] = -2.36917200E-01;
    COFD[2015] = 1.00534193E-02;
    COFD[2016] = -2.13316826E+01;
    COFD[2017] = 5.31504777E+00;
    COFD[2018] = -4.45287609E-01;
    COFD[2019] = 1.80982140E-02;
    COFD[2020] = -2.16008317E+01;
    COFD[2021] = 5.38202048E+00;
    COFD[2022] = -4.48902904E-01;
    COFD[2023] = 1.80661615E-02;
    COFD[2024] = -2.15926683E+01;
    COFD[2025] = 5.38202048E+00;
    COFD[2026] = -4.48902904E-01;
    COFD[2027] = 1.80661615E-02;
    COFD[2028] = -1.76822324E+01;
    COFD[2029] = 4.17392086E+00;
    COFD[2030] = -3.16668286E-01;
    COFD[2031] = 1.32724062E-02;
    COFD[2032] = -1.77048073E+01;
    COFD[2033] = 4.17392086E+00;
    COFD[2034] = -3.16668286E-01;
    COFD[2035] = 1.32724062E-02;
    COFD[2036] = -1.76468178E+01;
    COFD[2037] = 4.15115465E+00;
    COFD[2038] = -3.13942179E-01;
    COFD[2039] = 1.31634405E-02;
    COFD[2040] = -1.92128578E+01;
    COFD[2041] = 4.60990601E+00;
    COFD[2042] = -3.67990757E-01;
    COFD[2043] = 1.52895090E-02;
    COFD[2044] = -2.02603120E+01;
    COFD[2045] = 4.97397948E+00;
    COFD[2046] = -4.10308789E-01;
    COFD[2047] = 1.69343395E-02;
    COFD[2048] = -2.00303008E+01;
    COFD[2049] = 4.87104087E+00;
    COFD[2050] = -3.98982435E-01;
    COFD[2051] = 1.65207609E-02;
    COFD[2052] = -2.09109892E+01;
    COFD[2053] = 5.13479452E+00;
    COFD[2054] = -4.28125931E-01;
    COFD[2055] = 1.75896529E-02;
    COFD[2056] = -2.07944331E+01;
    COFD[2057] = 5.09987004E+00;
    COFD[2058] = -4.24162283E-01;
    COFD[2059] = 1.74396530E-02;
    COFD[2060] = -2.04475834E+01;
    COFD[2061] = 4.93091032E+00;
    COFD[2062] = -4.05754875E-01;
    COFD[2063] = 1.67763601E-02;
    COFD[2064] = -2.04521806E+01;
    COFD[2065] = 4.93091032E+00;
    COFD[2066] = -4.05754875E-01;
    COFD[2067] = 1.67763601E-02;
    COFD[2068] = -2.25737490E+01;
    COFD[2069] = 5.43256416E+00;
    COFD[2070] = -4.51725369E-01;
    COFD[2071] = 1.80507219E-02;
    COFD[2072] = -2.26516120E+01;
    COFD[2073] = 5.43033150E+00;
    COFD[2074] = -4.51600069E-01;
    COFD[2075] = 1.80513342E-02;
    COFD[2076] = -2.27589369E+01;
    COFD[2077] = 5.44369424E+00;
    COFD[2078] = -4.52241964E-01;
    COFD[2079] = 1.80389808E-02;
    COFD[2080] = -2.23500495E+01;
    COFD[2081] = 5.38231713E+00;
    COFD[2082] = -4.48905748E-01;
    COFD[2083] = 1.80649697E-02;
    COFD[2084] = -1.65468024E+01;
    COFD[2085] = 3.71875916E+00;
    COFD[2086] = -2.61726361E-01;
    COFD[2087] = 1.10588035E-02;
    COFD[2088] = -2.27314400E+01;
    COFD[2089] = 5.39729337E+00;
    COFD[2090] = -4.37101137E-01;
    COFD[2091] = 1.70521181E-02;
    COFD[2092] = -1.78073254E+01;
    COFD[2093] = 4.17648998E+00;
    COFD[2094] = -3.16975886E-01;
    COFD[2095] = 1.32846969E-02;
    COFD[2096] = -2.08088198E+01;
    COFD[2097] = 5.14475868E+00;
    COFD[2098] = -4.29211829E-01;
    COFD[2099] = 1.76287935E-02;
    COFD[2100] = -2.15465990E+01;
    COFD[2101] = 5.38978785E+00;
    COFD[2102] = -4.35592956E-01;
    COFD[2103] = 1.69676354E-02;
    COFD[2104] = -1.74857137E+01;
    COFD[2105] = 4.05605079E+00;
    COFD[2106] = -3.02454805E-01;
    COFD[2107] = 1.27004015E-02;
    COFD[2108] = -1.33954202E+01;
    COFD[2109] = 3.01958105E+00;
    COFD[2110] = -1.74891189E-01;
    COFD[2111] = 7.46027223E-03;
    COFD[2112] = -1.63733545E+01;
    COFD[2113] = 3.79864092E+00;
    COFD[2114] = -2.71062595E-01;
    COFD[2115] = 1.14218887E-02;
    COFD[2116] = -1.78213641E+01;
    COFD[2117] = 4.17648998E+00;
    COFD[2118] = -3.16975886E-01;
    COFD[2119] = 1.32846969E-02;
    COFD[2120] = -1.78145099E+01;
    COFD[2121] = 4.17648998E+00;
    COFD[2122] = -3.16975886E-01;
    COFD[2123] = 1.32846969E-02;
    COFD[2124] = -1.70609605E+01;
    COFD[2125] = 4.53774238E+00;
    COFD[2126] = -3.59864268E-01;
    COFD[2127] = 1.49858355E-02;
    COFD[2128] = -1.63541181E+01;
    COFD[2129] = 3.79864092E+00;
    COFD[2130] = -2.71062595E-01;
    COFD[2131] = 1.14218887E-02;
    COFD[2132] = -2.18246157E+01;
    COFD[2133] = 5.42978824E+00;
    COFD[2134] = -4.51570216E-01;
    COFD[2135] = 1.80515349E-02;
    COFD[2136] = -2.19239746E+01;
    COFD[2137] = 5.42213817E+00;
    COFD[2138] = -4.44597702E-01;
    COFD[2139] = 1.75307003E-02;
    COFD[2140] = -2.19156553E+01;
    COFD[2141] = 5.42213817E+00;
    COFD[2142] = -4.44597702E-01;
    COFD[2143] = 1.75307003E-02;
    COFD[2144] = -1.86827277E+01;
    COFD[2145] = 4.53129263E+00;
    COFD[2146] = -3.59157354E-01;
    COFD[2147] = 1.49603162E-02;
    COFD[2148] = -1.87055896E+01;
    COFD[2149] = 4.53129263E+00;
    COFD[2150] = -3.59157354E-01;
    COFD[2151] = 1.49603162E-02;
    COFD[2152] = -1.86611512E+01;
    COFD[2153] = 4.51482946E+00;
    COFD[2154] = -3.57370004E-01;
    COFD[2155] = 1.48965915E-02;
    COFD[2156] = -2.02603120E+01;
    COFD[2157] = 4.97397948E+00;
    COFD[2158] = -4.10308789E-01;
    COFD[2159] = 1.69343395E-02;
    COFD[2160] = -2.11347359E+01;
    COFD[2161] = 5.25546346E+00;
    COFD[2162] = -4.40168690E-01;
    COFD[2163] = 1.79744174E-02;
    COFD[2164] = -2.09622452E+01;
    COFD[2165] = 5.17994790E+00;
    COFD[2166] = -4.33032053E-01;
    COFD[2167] = 1.77658771E-02;
    COFD[2168] = -2.15980732E+01;
    COFD[2169] = 5.33463542E+00;
    COFD[2170] = -4.46573990E-01;
    COFD[2171] = 1.81071379E-02;
    COFD[2172] = -2.15309986E+01;
    COFD[2173] = 5.32069306E+00;
    COFD[2174] = -4.45511865E-01;
    COFD[2175] = 1.80889760E-02;
    COFD[2176] = -2.13630095E+01;
    COFD[2177] = 5.22906675E+00;
    COFD[2178] = -4.38093497E-01;
    COFD[2179] = 1.79353664E-02;
    COFD[2180] = -2.13677101E+01;
    COFD[2181] = 5.22906675E+00;
    COFD[2182] = -4.38093497E-01;
    COFD[2183] = 1.79353664E-02;
    COFD[2184] = -2.27303300E+01;
    COFD[2185] = 5.39729337E+00;
    COFD[2186] = -4.37101137E-01;
    COFD[2187] = 1.70521181E-02;
    COFD[2188] = -2.28218940E+01;
    COFD[2189] = 5.40010147E+00;
    COFD[2190] = -4.37662496E-01;
    COFD[2191] = 1.70835057E-02;
    COFD[2192] = -2.28533080E+01;
    COFD[2193] = 5.38163297E+00;
    COFD[2194] = -4.33974802E-01;
    COFD[2195] = 1.68774596E-02;
    COFD[2196] = -2.26848257E+01;
    COFD[2197] = 5.42187576E+00;
    COFD[2198] = -4.44522343E-01;
    COFD[2199] = 1.75259406E-02;
    COFD[2200] = -1.74553113E+01;
    COFD[2201] = 4.04791428E+00;
    COFD[2202] = -3.01466078E-01;
    COFD[2203] = 1.26603027E-02;
    COFD[2204] = -2.27800509E+01;
    COFD[2205] = 5.42267938E+00;
    COFD[2206] = -4.44086802E-01;
    COFD[2207] = 1.74875185E-02;
    COFD[2208] = -1.75004523E+01;
    COFD[2209] = 4.03367299E+00;
    COFD[2210] = -2.99744218E-01;
    COFD[2211] = 1.25908330E-02;
    COFD[2212] = -2.04989678E+01;
    COFD[2213] = 5.01033357E+00;
    COFD[2214] = -4.14008048E-01;
    COFD[2215] = 1.70559642E-02;
    COFD[2216] = -2.16415067E+01;
    COFD[2217] = 5.42347394E+00;
    COFD[2218] = -4.43696048E-01;
    COFD[2219] = 1.74526528E-02;
    COFD[2220] = -1.71875552E+01;
    COFD[2221] = 3.91643623E+00;
    COFD[2222] = -2.85489387E-01;
    COFD[2223] = 1.20121308E-02;
    COFD[2224] = -1.32148679E+01;
    COFD[2225] = 2.91941286E+00;
    COFD[2226] = -1.62365908E-01;
    COFD[2227] = 6.93774959E-03;
    COFD[2228] = -1.62164844E+01;
    COFD[2229] = 3.71034643E+00;
    COFD[2230] = -2.60739292E-01;
    COFD[2231] = 1.10202921E-02;
    COFD[2232] = -1.75147592E+01;
    COFD[2233] = 4.03367299E+00;
    COFD[2234] = -2.99744218E-01;
    COFD[2235] = 1.25908330E-02;
    COFD[2236] = -1.75077731E+01;
    COFD[2237] = 4.03367299E+00;
    COFD[2238] = -2.99744218E-01;
    COFD[2239] = 1.25908330E-02;
    COFD[2240] = -1.68129119E+01;
    COFD[2241] = 4.41680518E+00;
    COFD[2242] = -3.45824229E-01;
    COFD[2243] = 1.44423400E-02;
    COFD[2244] = -1.61969977E+01;
    COFD[2245] = 3.71034643E+00;
    COFD[2246] = -2.60739292E-01;
    COFD[2247] = 1.10202921E-02;
    COFD[2248] = -2.17177863E+01;
    COFD[2249] = 5.38391804E+00;
    COFD[2250] = -4.48822101E-01;
    COFD[2251] = 1.80505757E-02;
    COFD[2252] = -2.19859525E+01;
    COFD[2253] = 5.44658055E+00;
    COFD[2254] = -4.51507914E-01;
    COFD[2255] = 1.79664843E-02;
    COFD[2256] = -2.19774825E+01;
    COFD[2257] = 5.44658055E+00;
    COFD[2258] = -4.51507914E-01;
    COFD[2259] = 1.79664843E-02;
    COFD[2260] = -1.84057376E+01;
    COFD[2261] = 4.40577261E+00;
    COFD[2262] = -3.44465891E-01;
    COFD[2263] = 1.43864134E-02;
    COFD[2264] = -1.84288731E+01;
    COFD[2265] = 4.40577261E+00;
    COFD[2266] = -3.44465891E-01;
    COFD[2267] = 1.43864134E-02;
    COFD[2268] = -1.83550733E+01;
    COFD[2269] = 4.37621420E+00;
    COFD[2270] = -3.40819809E-01;
    COFD[2271] = 1.42360254E-02;
    COFD[2272] = -2.00303008E+01;
    COFD[2273] = 4.87104087E+00;
    COFD[2274] = -3.98982435E-01;
    COFD[2275] = 1.65207609E-02;
    COFD[2276] = -2.09622452E+01;
    COFD[2277] = 5.17994790E+00;
    COFD[2278] = -4.33032053E-01;
    COFD[2279] = 1.77658771E-02;
    COFD[2280] = -2.06551644E+01;
    COFD[2281] = 5.04973552E+00;
    COFD[2282] = -4.18409723E-01;
    COFD[2283] = 1.72193358E-02;
    COFD[2284] = -2.14709052E+01;
    COFD[2285] = 5.27933403E+00;
    COFD[2286] = -4.42255843E-01;
    COFD[2287] = 1.80267390E-02;
    COFD[2288] = -2.13958881E+01;
    COFD[2289] = 5.26174490E+00;
    COFD[2290] = -4.40699830E-01;
    COFD[2291] = 1.79867710E-02;
    COFD[2292] = -2.11150115E+01;
    COFD[2293] = 5.12506158E+00;
    COFD[2294] = -4.27045545E-01;
    COFD[2295] = 1.75498196E-02;
    COFD[2296] = -2.11198124E+01;
    COFD[2297] = 5.12506158E+00;
    COFD[2298] = -4.27045545E-01;
    COFD[2299] = 1.75498196E-02;
    COFD[2300] = -2.27789100E+01;
    COFD[2301] = 5.42267938E+00;
    COFD[2302] = -4.44086802E-01;
    COFD[2303] = 1.74875185E-02;
    COFD[2304] = -2.28608979E+01;
    COFD[2305] = 5.42202839E+00;
    COFD[2306] = -4.44180027E-01;
    COFD[2307] = 1.74980028E-02;
    COFD[2308] = -2.29391484E+01;
    COFD[2309] = 5.42308167E+00;
    COFD[2310] = -4.43120534E-01;
    COFD[2311] = 1.74086740E-02;
    COFD[2312] = -2.27315839E+01;
    COFD[2313] = 5.44654266E+00;
    COFD[2314] = -4.51461923E-01;
    COFD[2315] = 1.79629833E-02;
    COFD[2316] = -1.71615857E+01;
    COFD[2317] = 3.90982668E+00;
    COFD[2318] = -2.84689179E-01;
    COFD[2319] = 1.19797702E-02;
    COFD[2320] = -2.27953025E+01;
    COFD[2321] = 5.29006863E+00;
    COFD[2322] = -4.17624653E-01;
    COFD[2323] = 1.60089937E-02;
    COFD[2324] = -1.84552630E+01;
    COFD[2325] = 4.34181535E+00;
    COFD[2326] = -3.36587321E-01;
    COFD[2327] = 1.40618044E-02;
    COFD[2328] = -2.13743701E+01;
    COFD[2329] = 5.25833649E+00;
    COFD[2330] = -4.40410969E-01;
    COFD[2331] = 1.79800167E-02;
    COFD[2332] = -2.15592524E+01;
    COFD[2333] = 5.27728120E+00;
    COFD[2334] = -4.15361581E-01;
    COFD[2335] = 1.58893774E-02;
    COFD[2336] = -1.81728784E+01;
    COFD[2337] = 4.24023311E+00;
    COFD[2338] = -3.24539148E-01;
    COFD[2339] = 1.35842242E-02;
    COFD[2340] = -1.39240907E+01;
    COFD[2341] = 3.17763277E+00;
    COFD[2342] = -1.94889823E-01;
    COFD[2343] = 8.30436995E-03;
    COFD[2344] = -1.70275170E+01;
    COFD[2345] = 3.97787104E+00;
    COFD[2346] = -2.92958143E-01;
    COFD[2347] = 1.23153094E-02;
    COFD[2348] = -1.84720157E+01;
    COFD[2349] = 4.34181535E+00;
    COFD[2350] = -3.36587321E-01;
    COFD[2351] = 1.40618044E-02;
    COFD[2352] = -1.84638250E+01;
    COFD[2353] = 4.34181535E+00;
    COFD[2354] = -3.36587321E-01;
    COFD[2355] = 1.40618044E-02;
    COFD[2356] = -1.76744140E+01;
    COFD[2357] = 4.70907386E+00;
    COFD[2358] = -3.79700042E-01;
    COFD[2359] = 1.57519279E-02;
    COFD[2360] = -1.70058787E+01;
    COFD[2361] = 3.97787104E+00;
    COFD[2362] = -2.92958143E-01;
    COFD[2363] = 1.23153094E-02;
    COFD[2364] = -2.21211461E+01;
    COFD[2365] = 5.43170153E+00;
    COFD[2366] = -4.46963936E-01;
    COFD[2367] = 1.76740849E-02;
    COFD[2368] = -2.21293977E+01;
    COFD[2369] = 5.38638208E+00;
    COFD[2370] = -4.34919926E-01;
    COFD[2371] = 1.69301831E-02;
    COFD[2372] = -2.21195673E+01;
    COFD[2373] = 5.38638208E+00;
    COFD[2374] = -4.34919926E-01;
    COFD[2375] = 1.69301831E-02;
    COFD[2376] = -1.93089486E+01;
    COFD[2377] = 4.69775512E+00;
    COFD[2378] = -3.78328786E-01;
    COFD[2379] = 1.56963078E-02;
    COFD[2380] = -1.93344144E+01;
    COFD[2381] = 4.69775512E+00;
    COFD[2382] = -3.78328786E-01;
    COFD[2383] = 1.56963078E-02;
    COFD[2384] = -1.92602241E+01;
    COFD[2385] = 4.66778124E+00;
    COFD[2386] = -3.74693124E-01;
    COFD[2387] = 1.55486448E-02;
    COFD[2388] = -2.09109892E+01;
    COFD[2389] = 5.13479452E+00;
    COFD[2390] = -4.28125931E-01;
    COFD[2391] = 1.75896529E-02;
    COFD[2392] = -2.15980732E+01;
    COFD[2393] = 5.33463542E+00;
    COFD[2394] = -4.46573990E-01;
    COFD[2395] = 1.81071379E-02;
    COFD[2396] = -2.14709052E+01;
    COFD[2397] = 5.27933403E+00;
    COFD[2398] = -4.42255843E-01;
    COFD[2399] = 1.80267390E-02;
    COFD[2400] = -2.20068020E+01;
    COFD[2401] = 5.38409914E+00;
    COFD[2402] = -4.48760682E-01;
    COFD[2403] = 1.80447792E-02;
    COFD[2404] = -2.19620149E+01;
    COFD[2405] = 5.37894961E+00;
    COFD[2406] = -4.48923484E-01;
    COFD[2407] = 1.80825469E-02;
    COFD[2408] = -2.18467821E+01;
    COFD[2409] = 5.31331054E+00;
    COFD[2410] = -4.45223645E-01;
    COFD[2411] = 1.81014539E-02;
    COFD[2412] = -2.18525195E+01;
    COFD[2413] = 5.31331054E+00;
    COFD[2414] = -4.45223645E-01;
    COFD[2415] = 1.81014539E-02;
    COFD[2416] = -2.27938535E+01;
    COFD[2417] = 5.29006863E+00;
    COFD[2418] = -4.17624653E-01;
    COFD[2419] = 1.60089937E-02;
    COFD[2420] = -2.29016678E+01;
    COFD[2421] = 5.29681598E+00;
    COFD[2422] = -4.18727614E-01;
    COFD[2423] = 1.60648947E-02;
    COFD[2424] = -2.29044907E+01;
    COFD[2425] = 5.26490219E+00;
    COFD[2426] = -4.13127922E-01;
    COFD[2427] = 1.57702008E-02;
    COFD[2428] = -2.29185460E+01;
    COFD[2429] = 5.38576397E+00;
    COFD[2430] = -4.34796260E-01;
    COFD[2431] = 1.69232691E-02;
    COFD[2432] = -1.81473610E+01;
    COFD[2433] = 4.23397947E+00;
    COFD[2434] = -3.23805200E-01;
    COFD[2435] = 1.35554854E-02;
    COFD[2436] = -2.28363010E+01;
    COFD[2437] = 5.31907985E+00;
    COFD[2438] = -4.22514116E-01;
    COFD[2439] = 1.62610921E-02;
    COFD[2440] = -1.83508803E+01;
    COFD[2441] = 4.30984197E+00;
    COFD[2442] = -3.32770986E-01;
    COFD[2443] = 1.39094747E-02;
    COFD[2444] = -2.13081407E+01;
    COFD[2445] = 5.24403923E+00;
    COFD[2446] = -4.39272407E-01;
    COFD[2447] = 1.79576065E-02;
    COFD[2448] = -2.15823819E+01;
    COFD[2449] = 5.30395665E+00;
    COFD[2450] = -4.19948713E-01;
    COFD[2451] = 1.61283165E-02;
    COFD[2452] = -1.80716029E+01;
    COFD[2453] = 4.20964822E+00;
    COFD[2454] = -3.20930845E-01;
    COFD[2455] = 1.34421679E-02;
    COFD[2456] = -1.38287706E+01;
    COFD[2457] = 3.14964252E+00;
    COFD[2458] = -1.91369421E-01;
    COFD[2459] = 8.15663173E-03;
    COFD[2460] = -1.69092218E+01;
    COFD[2461] = 3.94140805E+00;
    COFD[2462] = -2.88511662E-01;
    COFD[2463] = 1.21342675E-02;
    COFD[2464] = -1.83678207E+01;
    COFD[2465] = 4.30984197E+00;
    COFD[2466] = -3.32770986E-01;
    COFD[2467] = 1.39094747E-02;
    COFD[2468] = -1.83595375E+01;
    COFD[2469] = 4.30984197E+00;
    COFD[2470] = -3.32770986E-01;
    COFD[2471] = 1.39094747E-02;
    COFD[2472] = -1.75284656E+01;
    COFD[2473] = 4.66504845E+00;
    COFD[2474] = -3.74375842E-01;
    COFD[2475] = 1.55363348E-02;
    COFD[2476] = -1.68874275E+01;
    COFD[2477] = 3.94140805E+00;
    COFD[2478] = -2.88511662E-01;
    COFD[2479] = 1.21342675E-02;
    COFD[2480] = -2.21010249E+01;
    COFD[2481] = 5.43812778E+00;
    COFD[2482] = -4.48750136E-01;
    COFD[2483] = 1.77861511E-02;
    COFD[2484] = -2.21276873E+01;
    COFD[2485] = 5.40064912E+00;
    COFD[2486] = -4.37771126E-01;
    COFD[2487] = 1.70895608E-02;
    COFD[2488] = -2.21177535E+01;
    COFD[2489] = 5.40064912E+00;
    COFD[2490] = -4.37771126E-01;
    COFD[2491] = 1.70895608E-02;
    COFD[2492] = -1.91790889E+01;
    COFD[2493] = 4.65710841E+00;
    COFD[2494] = -3.73456882E-01;
    COFD[2495] = 1.55008064E-02;
    COFD[2496] = -1.92047221E+01;
    COFD[2497] = 4.65710841E+00;
    COFD[2498] = -3.73456882E-01;
    COFD[2499] = 1.55008064E-02;
    COFD[2500] = -1.91516731E+01;
    COFD[2501] = 4.63565403E+00;
    COFD[2502] = -3.70964470E-01;
    COFD[2503] = 1.54041133E-02;
    COFD[2504] = -2.07944331E+01;
    COFD[2505] = 5.09987004E+00;
    COFD[2506] = -4.24162283E-01;
    COFD[2507] = 1.74396530E-02;
    COFD[2508] = -2.15309986E+01;
    COFD[2509] = 5.32069306E+00;
    COFD[2510] = -4.45511865E-01;
    COFD[2511] = 1.80889760E-02;
    COFD[2512] = -2.13958881E+01;
    COFD[2513] = 5.26174490E+00;
    COFD[2514] = -4.40699830E-01;
    COFD[2515] = 1.79867710E-02;
    COFD[2516] = -2.19620149E+01;
    COFD[2517] = 5.37894961E+00;
    COFD[2518] = -4.48923484E-01;
    COFD[2519] = 1.80825469E-02;
    COFD[2520] = -2.19186788E+01;
    COFD[2521] = 5.37442568E+00;
    COFD[2522] = -4.49173971E-01;
    COFD[2523] = 1.81244325E-02;
    COFD[2524] = -2.17818407E+01;
    COFD[2525] = 5.29865944E+00;
    COFD[2526] = -4.44028510E-01;
    COFD[2527] = 1.80758390E-02;
    COFD[2528] = -2.17876516E+01;
    COFD[2529] = 5.29865944E+00;
    COFD[2530] = -4.44028510E-01;
    COFD[2531] = 1.80758390E-02;
    COFD[2532] = -2.28348262E+01;
    COFD[2533] = 5.31907985E+00;
    COFD[2534] = -4.22514116E-01;
    COFD[2535] = 1.62610921E-02;
    COFD[2536] = -2.29409328E+01;
    COFD[2537] = 5.32468208E+00;
    COFD[2538] = -4.23462699E-01;
    COFD[2539] = 1.63101476E-02;
    COFD[2540] = -2.29281623E+01;
    COFD[2541] = 5.28652757E+00;
    COFD[2542] = -4.17049438E-01;
    COFD[2543] = 1.59799372E-02;
    COFD[2544] = -2.29245190E+01;
    COFD[2545] = 5.40000327E+00;
    COFD[2546] = -4.37643030E-01;
    COFD[2547] = 1.70824209E-02;
    COFD[2548] = -1.80453507E+01;
    COFD[2549] = 4.20312679E+00;
    COFD[2550] = -3.20161664E-01;
    COFD[2551] = 1.34118992E-02;
    COFD[2552] = -2.30934035E+01;
    COFD[2553] = 5.42271039E+00;
    COFD[2554] = -4.42221513E-01;
    COFD[2555] = 1.73390244E-02;
    COFD[2556] = -1.79735387E+01;
    COFD[2557] = 4.11071307E+00;
    COFD[2558] = -3.09074203E-01;
    COFD[2559] = 1.29678777E-02;
    COFD[2560] = -2.09751443E+01;
    COFD[2561] = 5.08056463E+00;
    COFD[2562] = -4.21952793E-01;
    COFD[2563] = 1.73552764E-02;
    COFD[2564] = -2.19259317E+01;
    COFD[2565] = 5.41758144E+00;
    COFD[2566] = -4.41053871E-01;
    COFD[2567] = 1.72705882E-02;
    COFD[2568] = -1.76339606E+01;
    COFD[2569] = 3.98770708E+00;
    COFD[2570] = -2.94160646E-01;
    COFD[2571] = 1.23644039E-02;
    COFD[2572] = -1.35245705E+01;
    COFD[2573] = 2.96195761E+00;
    COFD[2574] = -1.67607052E-01;
    COFD[2575] = 7.15329427E-03;
    COFD[2576] = -1.66034546E+01;
    COFD[2577] = 3.75923792E+00;
    COFD[2578] = -2.66511261E-01;
    COFD[2579] = 1.12471370E-02;
    COFD[2580] = -1.79906619E+01;
    COFD[2581] = 4.11071307E+00;
    COFD[2582] = -3.09074203E-01;
    COFD[2583] = 1.29678777E-02;
    COFD[2584] = -1.79822885E+01;
    COFD[2585] = 4.11071307E+00;
    COFD[2586] = -3.09074203E-01;
    COFD[2587] = 1.29678777E-02;
    COFD[2588] = -1.72372988E+01;
    COFD[2589] = 4.49320169E+00;
    COFD[2590] = -3.55062447E-01;
    COFD[2591] = 1.48161449E-02;
    COFD[2592] = -1.65815097E+01;
    COFD[2593] = 3.75923792E+00;
    COFD[2594] = -2.66511261E-01;
    COFD[2595] = 1.12471370E-02;
    COFD[2596] = -2.20454826E+01;
    COFD[2597] = 5.39682182E+00;
    COFD[2598] = -4.48776824E-01;
    COFD[2599] = 1.79851810E-02;
    COFD[2600] = -2.22599437E+01;
    COFD[2601] = 5.43626968E+00;
    COFD[2602] = -4.48222611E-01;
    COFD[2603] = 1.77528565E-02;
    COFD[2604] = -2.22499095E+01;
    COFD[2605] = 5.43626968E+00;
    COFD[2606] = -4.48222611E-01;
    COFD[2607] = 1.77528565E-02;
    COFD[2608] = -1.88470808E+01;
    COFD[2609] = 4.48692836E+00;
    COFD[2610] = -3.54362855E-01;
    COFD[2611] = 1.47902773E-02;
    COFD[2612] = -1.88728753E+01;
    COFD[2613] = 4.48692836E+00;
    COFD[2614] = -3.54362855E-01;
    COFD[2615] = 1.47902773E-02;
    COFD[2616] = -1.88205264E+01;
    COFD[2617] = 4.46474713E+00;
    COFD[2618] = -3.51718970E-01;
    COFD[2619] = 1.46847779E-02;
    COFD[2620] = -2.04475834E+01;
    COFD[2621] = 4.93091032E+00;
    COFD[2622] = -4.05754875E-01;
    COFD[2623] = 1.67763601E-02;
    COFD[2624] = -2.13630095E+01;
    COFD[2625] = 5.22906675E+00;
    COFD[2626] = -4.38093497E-01;
    COFD[2627] = 1.79353664E-02;
    COFD[2628] = -2.11150115E+01;
    COFD[2629] = 5.12506158E+00;
    COFD[2630] = -4.27045545E-01;
    COFD[2631] = 1.75498196E-02;
    COFD[2632] = -2.18467821E+01;
    COFD[2633] = 5.31331054E+00;
    COFD[2634] = -4.45223645E-01;
    COFD[2635] = 1.81014539E-02;
    COFD[2636] = -2.17818407E+01;
    COFD[2637] = 5.29865944E+00;
    COFD[2638] = -4.44028510E-01;
    COFD[2639] = 1.80758390E-02;
    COFD[2640] = -2.15435537E+01;
    COFD[2641] = 5.18484753E+00;
    COFD[2642] = -4.33562002E-01;
    COFD[2643] = 1.77848067E-02;
    COFD[2644] = -2.15494363E+01;
    COFD[2645] = 5.18484753E+00;
    COFD[2646] = -4.33562002E-01;
    COFD[2647] = 1.77848067E-02;
    COFD[2648] = -2.30919032E+01;
    COFD[2649] = 5.42271039E+00;
    COFD[2650] = -4.42221513E-01;
    COFD[2651] = 1.73390244E-02;
    COFD[2652] = -2.31810931E+01;
    COFD[2653] = 5.42262618E+00;
    COFD[2654] = -4.42380851E-01;
    COFD[2655] = 1.73519942E-02;
    COFD[2656] = -2.32239159E+01;
    COFD[2657] = 5.40902423E+00;
    COFD[2658] = -4.39389466E-01;
    COFD[2659] = 1.71787145E-02;
    COFD[2660] = -2.30186448E+01;
    COFD[2661] = 5.43600467E+00;
    COFD[2662] = -4.48146589E-01;
    COFD[2663] = 1.77480458E-02;
    COFD[2664] = -1.76051897E+01;
    COFD[2665] = 3.97979335E+00;
    COFD[2666] = -2.93193093E-01;
    COFD[2667] = 1.23248995E-02;
    COFD[2668] = -2.31017092E+01;
    COFD[2669] = 5.42271039E+00;
    COFD[2670] = -4.42221513E-01;
    COFD[2671] = 1.73390244E-02;
    COFD[2672] = -1.79786166E+01;
    COFD[2673] = 4.11071307E+00;
    COFD[2674] = -3.09074203E-01;
    COFD[2675] = 1.29678777E-02;
    COFD[2676] = -2.09811595E+01;
    COFD[2677] = 5.08056463E+00;
    COFD[2678] = -4.21952793E-01;
    COFD[2679] = 1.73552764E-02;
    COFD[2680] = -2.19294503E+01;
    COFD[2681] = 5.41758144E+00;
    COFD[2682] = -4.41053871E-01;
    COFD[2683] = 1.72705882E-02;
    COFD[2684] = -1.76386567E+01;
    COFD[2685] = 3.98770708E+00;
    COFD[2686] = -2.94160646E-01;
    COFD[2687] = 1.23644039E-02;
    COFD[2688] = -1.35251055E+01;
    COFD[2689] = 2.96195761E+00;
    COFD[2690] = -1.67607052E-01;
    COFD[2691] = 7.15329427E-03;
    COFD[2692] = -1.66068325E+01;
    COFD[2693] = 3.75923792E+00;
    COFD[2694] = -2.66511261E-01;
    COFD[2695] = 1.12471370E-02;
    COFD[2696] = -1.79959176E+01;
    COFD[2697] = 4.11071307E+00;
    COFD[2698] = -3.09074203E-01;
    COFD[2699] = 1.29678777E-02;
    COFD[2700] = -1.79874565E+01;
    COFD[2701] = 4.11071307E+00;
    COFD[2702] = -3.09074203E-01;
    COFD[2703] = 1.29678777E-02;
    COFD[2704] = -1.72375725E+01;
    COFD[2705] = 4.49320169E+00;
    COFD[2706] = -3.55062447E-01;
    COFD[2707] = 1.48161449E-02;
    COFD[2708] = -1.65847421E+01;
    COFD[2709] = 3.75923792E+00;
    COFD[2710] = -2.66511261E-01;
    COFD[2711] = 1.12471370E-02;
    COFD[2712] = -2.20504720E+01;
    COFD[2713] = 5.39682182E+00;
    COFD[2714] = -4.48776824E-01;
    COFD[2715] = 1.79851810E-02;
    COFD[2716] = -2.22648380E+01;
    COFD[2717] = 5.43626968E+00;
    COFD[2718] = -4.48222611E-01;
    COFD[2719] = 1.77528565E-02;
    COFD[2720] = -2.22547061E+01;
    COFD[2721] = 5.43626968E+00;
    COFD[2722] = -4.48222611E-01;
    COFD[2723] = 1.77528565E-02;
    COFD[2724] = -1.88500135E+01;
    COFD[2725] = 4.48692836E+00;
    COFD[2726] = -3.54362855E-01;
    COFD[2727] = 1.47902773E-02;
    COFD[2728] = -1.88759638E+01;
    COFD[2729] = 4.48692836E+00;
    COFD[2730] = -3.54362855E-01;
    COFD[2731] = 1.47902773E-02;
    COFD[2732] = -1.88237651E+01;
    COFD[2733] = 4.46474713E+00;
    COFD[2734] = -3.51718970E-01;
    COFD[2735] = 1.46847779E-02;
    COFD[2736] = -2.04521806E+01;
    COFD[2737] = 4.93091032E+00;
    COFD[2738] = -4.05754875E-01;
    COFD[2739] = 1.67763601E-02;
    COFD[2740] = -2.13677101E+01;
    COFD[2741] = 5.22906675E+00;
    COFD[2742] = -4.38093497E-01;
    COFD[2743] = 1.79353664E-02;
    COFD[2744] = -2.11198124E+01;
    COFD[2745] = 5.12506158E+00;
    COFD[2746] = -4.27045545E-01;
    COFD[2747] = 1.75498196E-02;
    COFD[2748] = -2.18525195E+01;
    COFD[2749] = 5.31331054E+00;
    COFD[2750] = -4.45223645E-01;
    COFD[2751] = 1.81014539E-02;
    COFD[2752] = -2.17876516E+01;
    COFD[2753] = 5.29865944E+00;
    COFD[2754] = -4.44028510E-01;
    COFD[2755] = 1.80758390E-02;
    COFD[2756] = -2.15494363E+01;
    COFD[2757] = 5.18484753E+00;
    COFD[2758] = -4.33562002E-01;
    COFD[2759] = 1.77848067E-02;
    COFD[2760] = -2.15553889E+01;
    COFD[2761] = 5.18484753E+00;
    COFD[2762] = -4.33562002E-01;
    COFD[2763] = 1.77848067E-02;
    COFD[2764] = -2.31001838E+01;
    COFD[2765] = 5.42271039E+00;
    COFD[2766] = -4.42221513E-01;
    COFD[2767] = 1.73390244E-02;
    COFD[2768] = -2.31900283E+01;
    COFD[2769] = 5.42262618E+00;
    COFD[2770] = -4.42380851E-01;
    COFD[2771] = 1.73519942E-02;
    COFD[2772] = -2.32330814E+01;
    COFD[2773] = 5.40902423E+00;
    COFD[2774] = -4.39389466E-01;
    COFD[2775] = 1.71787145E-02;
    COFD[2776] = -2.30269244E+01;
    COFD[2777] = 5.43600467E+00;
    COFD[2778] = -4.48146589E-01;
    COFD[2779] = 1.77480458E-02;
    COFD[2780] = -1.76098861E+01;
    COFD[2781] = 3.97979335E+00;
    COFD[2782] = -2.93193093E-01;
    COFD[2783] = 1.23248995E-02;
    COFD[2784] = -2.18886889E+01;
    COFD[2785] = 4.47348600E+00;
    COFD[2786] = -2.88470870E-01;
    COFD[2787] = 9.60875755E-03;
    COFD[2788] = -2.09411622E+01;
    COFD[2789] = 4.99879868E+00;
    COFD[2790] = -4.12878607E-01;
    COFD[2791] = 1.70210595E-02;
    COFD[2792] = -2.28038033E+01;
    COFD[2793] = 5.42628508E+00;
    COFD[2794] = -4.45690571E-01;
    COFD[2795] = 1.75981564E-02;
    COFD[2796] = -2.05558125E+01;
    COFD[2797] = 4.43728272E+00;
    COFD[2798] = -2.83087542E-01;
    COFD[2799] = 9.35217321E-03;
    COFD[2800] = -2.06902559E+01;
    COFD[2801] = 4.92442201E+00;
    COFD[2802] = -4.05069196E-01;
    COFD[2803] = 1.67525806E-02;
    COFD[2804] = -1.59936160E+01;
    COFD[2805] = 3.76265798E+00;
    COFD[2806] = -2.66912858E-01;
    COFD[2807] = 1.12628530E-02;
    COFD[2808] = -1.94888119E+01;
    COFD[2809] = 4.64703793E+00;
    COFD[2810] = -3.72283709E-01;
    COFD[2811] = 1.54551501E-02;
    COFD[2812] = -2.09640855E+01;
    COFD[2813] = 4.99879868E+00;
    COFD[2814] = -4.12878607E-01;
    COFD[2815] = 1.70210595E-02;
    COFD[2816] = -2.09528426E+01;
    COFD[2817] = 4.99879868E+00;
    COFD[2818] = -4.12878607E-01;
    COFD[2819] = 1.70210595E-02;
    COFD[2820] = -1.99238205E+01;
    COFD[2821] = 5.27803851E+00;
    COFD[2822] = -4.42140176E-01;
    COFD[2823] = 1.80237084E-02;
    COFD[2824] = -1.94626201E+01;
    COFD[2825] = 4.64703793E+00;
    COFD[2826] = -3.72283709E-01;
    COFD[2827] = 1.54551501E-02;
    COFD[2828] = -2.22213768E+01;
    COFD[2829] = 5.05132223E+00;
    COFD[2830] = -3.76867047E-01;
    COFD[2831] = 1.39005636E-02;
    COFD[2832] = -2.16105899E+01;
    COFD[2833] = 4.74703304E+00;
    COFD[2834] = -3.29522725E-01;
    COFD[2835] = 1.15765554E-02;
    COFD[2836] = -2.15974324E+01;
    COFD[2837] = 4.74703304E+00;
    COFD[2838] = -3.29522725E-01;
    COFD[2839] = 1.15765554E-02;
    COFD[2840] = -2.15562373E+01;
    COFD[2841] = 5.27353202E+00;
    COFD[2842] = -4.41738406E-01;
    COFD[2843] = 1.80132159E-02;
    COFD[2844] = -2.15865030E+01;
    COFD[2845] = 5.27353202E+00;
    COFD[2846] = -4.41738406E-01;
    COFD[2847] = 1.80132159E-02;
    COFD[2848] = -2.15613443E+01;
    COFD[2849] = 5.26193907E+00;
    COFD[2850] = -4.40716330E-01;
    COFD[2851] = 1.79871592E-02;
    COFD[2852] = -2.25737490E+01;
    COFD[2853] = 5.43256416E+00;
    COFD[2854] = -4.51725369E-01;
    COFD[2855] = 1.80507219E-02;
    COFD[2856] = -2.27303300E+01;
    COFD[2857] = 5.39729337E+00;
    COFD[2858] = -4.37101137E-01;
    COFD[2859] = 1.70521181E-02;
    COFD[2860] = -2.27789100E+01;
    COFD[2861] = 5.42267938E+00;
    COFD[2862] = -4.44086802E-01;
    COFD[2863] = 1.74875185E-02;
    COFD[2864] = -2.27938535E+01;
    COFD[2865] = 5.29006863E+00;
    COFD[2866] = -4.17624653E-01;
    COFD[2867] = 1.60089937E-02;
    COFD[2868] = -2.28348262E+01;
    COFD[2869] = 5.31907985E+00;
    COFD[2870] = -4.22514116E-01;
    COFD[2871] = 1.62610921E-02;
    COFD[2872] = -2.30919032E+01;
    COFD[2873] = 5.42271039E+00;
    COFD[2874] = -4.42221513E-01;
    COFD[2875] = 1.73390244E-02;
    COFD[2876] = -2.31001838E+01;
    COFD[2877] = 5.42271039E+00;
    COFD[2878] = -4.42221513E-01;
    COFD[2879] = 1.73390244E-02;
    COFD[2880] = -2.18861678E+01;
    COFD[2881] = 4.47348600E+00;
    COFD[2882] = -2.88470870E-01;
    COFD[2883] = 9.60875755E-03;
    COFD[2884] = -2.20307099E+01;
    COFD[2885] = 4.48693408E+00;
    COFD[2886] = -2.90462167E-01;
    COFD[2887] = 9.70336833E-03;
    COFD[2888] = -2.18947133E+01;
    COFD[2889] = 4.39616382E+00;
    COFD[2890] = -2.76994914E-01;
    COFD[2891] = 9.06224705E-03;
    COFD[2892] = -2.24393654E+01;
    COFD[2893] = 4.74495988E+00;
    COFD[2894] = -3.29202526E-01;
    COFD[2895] = 1.15609276E-02;
    COFD[2896] = -2.06657732E+01;
    COFD[2897] = 4.91831350E+00;
    COFD[2898] = -4.04379413E-01;
    COFD[2899] = 1.67265792E-02;
    COFD[2900] = -2.20335821E+01;
    COFD[2901] = 4.48693408E+00;
    COFD[2902] = -2.90462167E-01;
    COFD[2903] = 9.70336833E-03;
    COFD[2904] = -2.10558774E+01;
    COFD[2905] = 5.00399764E+00;
    COFD[2906] = -4.13379526E-01;
    COFD[2907] = 1.70361425E-02;
    COFD[2908] = -2.29076429E+01;
    COFD[2909] = 5.42460337E+00;
    COFD[2910] = -4.45270089E-01;
    COFD[2911] = 1.75725973E-02;
    COFD[2912] = -2.01811055E+01;
    COFD[2913] = 4.24520878E+00;
    COFD[2914] = -2.57219171E-01;
    COFD[2915] = 8.20894239E-03;
    COFD[2916] = -2.08064719E+01;
    COFD[2917] = 4.93178191E+00;
    COFD[2918] = -4.05846844E-01;
    COFD[2919] = 1.67795427E-02;
    COFD[2920] = -1.60961940E+01;
    COFD[2921] = 3.77271143E+00;
    COFD[2922] = -2.68084668E-01;
    COFD[2923] = 1.13083251E-02;
    COFD[2924] = -1.95587320E+01;
    COFD[2925] = 4.64075109E+00;
    COFD[2926] = -3.71554667E-01;
    COFD[2927] = 1.54269261E-02;
    COFD[2928] = -2.10482827E+01;
    COFD[2929] = 4.99490919E+00;
    COFD[2930] = -4.12490091E-01;
    COFD[2931] = 1.70086418E-02;
    COFD[2932] = -2.10363020E+01;
    COFD[2933] = 4.99490919E+00;
    COFD[2934] = -4.12490091E-01;
    COFD[2935] = 1.70086418E-02;
    COFD[2936] = -1.99886318E+01;
    COFD[2937] = 5.27453026E+00;
    COFD[2938] = -4.41827327E-01;
    COFD[2939] = 1.80155336E-02;
    COFD[2940] = -1.95315965E+01;
    COFD[2941] = 4.64075109E+00;
    COFD[2942] = -3.71554667E-01;
    COFD[2943] = 1.54269261E-02;
    COFD[2944] = -2.19824618E+01;
    COFD[2945] = 4.91387652E+00;
    COFD[2946] = -3.58234874E-01;
    COFD[2947] = 1.30693081E-02;
    COFD[2948] = -2.17219099E+01;
    COFD[2949] = 4.75637248E+00;
    COFD[2950] = -3.30964596E-01;
    COFD[2951] = 1.16469146E-02;
    COFD[2952] = -2.17079728E+01;
    COFD[2953] = 4.75637248E+00;
    COFD[2954] = -3.30964596E-01;
    COFD[2955] = 1.16469146E-02;
    COFD[2956] = -2.16217667E+01;
    COFD[2957] = 5.27002901E+00;
    COFD[2958] = -4.41426690E-01;
    COFD[2959] = 1.80051108E-02;
    COFD[2960] = -2.16530054E+01;
    COFD[2961] = 5.27002901E+00;
    COFD[2962] = -4.41426690E-01;
    COFD[2963] = 1.80051108E-02;
    COFD[2964] = -2.16639108E+01;
    COFD[2965] = 5.26764804E+00;
    COFD[2966] = -4.41215097E-01;
    COFD[2967] = 1.79996259E-02;
    COFD[2968] = -2.26516120E+01;
    COFD[2969] = 5.43033150E+00;
    COFD[2970] = -4.51600069E-01;
    COFD[2971] = 1.80513342E-02;
    COFD[2972] = -2.28218940E+01;
    COFD[2973] = 5.40010147E+00;
    COFD[2974] = -4.37662496E-01;
    COFD[2975] = 1.70835057E-02;
    COFD[2976] = -2.28608979E+01;
    COFD[2977] = 5.42202839E+00;
    COFD[2978] = -4.44180027E-01;
    COFD[2979] = 1.74980028E-02;
    COFD[2980] = -2.29016678E+01;
    COFD[2981] = 5.29681598E+00;
    COFD[2982] = -4.18727614E-01;
    COFD[2983] = 1.60648947E-02;
    COFD[2984] = -2.29409328E+01;
    COFD[2985] = 5.32468208E+00;
    COFD[2986] = -4.23462699E-01;
    COFD[2987] = 1.63101476E-02;
    COFD[2988] = -2.31810931E+01;
    COFD[2989] = 5.42262618E+00;
    COFD[2990] = -4.42380851E-01;
    COFD[2991] = 1.73519942E-02;
    COFD[2992] = -2.31900283E+01;
    COFD[2993] = 5.42262618E+00;
    COFD[2994] = -4.42380851E-01;
    COFD[2995] = 1.73519942E-02;
    COFD[2996] = -2.20307099E+01;
    COFD[2997] = 4.48693408E+00;
    COFD[2998] = -2.90462167E-01;
    COFD[2999] = 9.70336833E-03;
    COFD[3000] = -2.20424944E+01;
    COFD[3001] = 4.44475081E+00;
    COFD[3002] = -2.85218358E-01;
    COFD[3003] = 9.48659980E-03;
    COFD[3004] = -2.18910162E+01;
    COFD[3005] = 4.34657950E+00;
    COFD[3006] = -2.70778845E-01;
    COFD[3007] = 8.80340499E-03;
    COFD[3008] = -2.23779324E+01;
    COFD[3009] = 4.67548322E+00;
    COFD[3010] = -3.20205650E-01;
    COFD[3011] = 1.11747902E-02;
    COFD[3012] = -2.07822044E+01;
    COFD[3013] = 4.92605576E+00;
    COFD[3014] = -4.05242624E-01;
    COFD[3015] = 1.67586331E-02;
    COFD[3016] = -2.18977186E+01;
    COFD[3017] = 4.39616382E+00;
    COFD[3018] = -2.76994914E-01;
    COFD[3019] = 9.06224705E-03;
    COFD[3020] = -2.12327710E+01;
    COFD[3021] = 5.04534270E+00;
    COFD[3022] = -4.17903950E-01;
    COFD[3023] = 1.71998906E-02;
    COFD[3024] = -2.29881956E+01;
    COFD[3025] = 5.42338458E+00;
    COFD[3026] = -4.43824576E-01;
    COFD[3027] = 1.74633100E-02;
    COFD[3028] = -1.99759660E+01;
    COFD[3029] = 4.12684449E+00;
    COFD[3030] = -2.39932282E-01;
    COFD[3031] = 7.39544059E-03;
    COFD[3032] = -2.09575447E+01;
    COFD[3033] = 4.96410963E+00;
    COFD[3034] = -4.09246419E-01;
    COFD[3035] = 1.68966000E-02;
    COFD[3036] = -1.62329249E+01;
    COFD[3037] = 3.80480708E+00;
    COFD[3038] = -2.71809953E-01;
    COFD[3039] = 1.14521282E-02;
    COFD[3040] = -1.97432870E+01;
    COFD[3041] = 4.68838474E+00;
    COFD[3042] = -3.77193241E-01;
    COFD[3043] = 1.56502297E-02;
    COFD[3044] = -2.11872013E+01;
    COFD[3045] = 5.02181116E+00;
    COFD[3046] = -4.15179449E-01;
    COFD[3047] = 1.70945609E-02;
    COFD[3048] = -2.11749653E+01;
    COFD[3049] = 5.02181116E+00;
    COFD[3050] = -4.15179449E-01;
    COFD[3051] = 1.70945609E-02;
    COFD[3052] = -2.01190603E+01;
    COFD[3053] = 5.29848259E+00;
    COFD[3054] = -4.44011383E-01;
    COFD[3055] = 1.80753134E-02;
    COFD[3056] = -1.97158370E+01;
    COFD[3057] = 4.68838474E+00;
    COFD[3058] = -3.77193241E-01;
    COFD[3059] = 1.56502297E-02;
    COFD[3060] = -2.18763645E+01;
    COFD[3061] = 4.83720845E+00;
    COFD[3062] = -3.46627541E-01;
    COFD[3063] = 1.25098849E-02;
    COFD[3064] = -2.16433785E+01;
    COFD[3065] = 4.69172346E+00;
    COFD[3066] = -3.21053636E-01;
    COFD[3067] = 1.11654521E-02;
    COFD[3068] = -2.16291739E+01;
    COFD[3069] = 4.69172346E+00;
    COFD[3070] = -3.21053636E-01;
    COFD[3071] = 1.11654521E-02;
    COFD[3072] = -2.17475829E+01;
    COFD[3073] = 5.29328698E+00;
    COFD[3074] = -4.43508195E-01;
    COFD[3075] = 1.80599224E-02;
    COFD[3076] = -2.17791442E+01;
    COFD[3077] = 5.29328698E+00;
    COFD[3078] = -4.43508195E-01;
    COFD[3079] = 1.80599224E-02;
    COFD[3080] = -2.17991940E+01;
    COFD[3081] = 5.29320616E+00;
    COFD[3082] = -4.43500540E-01;
    COFD[3083] = 1.80596976E-02;
    COFD[3084] = -2.27589369E+01;
    COFD[3085] = 5.44369424E+00;
    COFD[3086] = -4.52241964E-01;
    COFD[3087] = 1.80389808E-02;
    COFD[3088] = -2.28533080E+01;
    COFD[3089] = 5.38163297E+00;
    COFD[3090] = -4.33974802E-01;
    COFD[3091] = 1.68774596E-02;
    COFD[3092] = -2.29391484E+01;
    COFD[3093] = 5.42308167E+00;
    COFD[3094] = -4.43120534E-01;
    COFD[3095] = 1.74086740E-02;
    COFD[3096] = -2.29044907E+01;
    COFD[3097] = 5.26490219E+00;
    COFD[3098] = -4.13127922E-01;
    COFD[3099] = 1.57702008E-02;
    COFD[3100] = -2.29281623E+01;
    COFD[3101] = 5.28652757E+00;
    COFD[3102] = -4.17049438E-01;
    COFD[3103] = 1.59799372E-02;
    COFD[3104] = -2.32239159E+01;
    COFD[3105] = 5.40902423E+00;
    COFD[3106] = -4.39389466E-01;
    COFD[3107] = 1.71787145E-02;
    COFD[3108] = -2.32330814E+01;
    COFD[3109] = 5.40902423E+00;
    COFD[3110] = -4.39389466E-01;
    COFD[3111] = 1.71787145E-02;
    COFD[3112] = -2.18947133E+01;
    COFD[3113] = 4.39616382E+00;
    COFD[3114] = -2.76994914E-01;
    COFD[3115] = 9.06224705E-03;
    COFD[3116] = -2.18910162E+01;
    COFD[3117] = 4.34657950E+00;
    COFD[3118] = -2.70778845E-01;
    COFD[3119] = 8.80340499E-03;
    COFD[3120] = -2.17238327E+01;
    COFD[3121] = 4.24131806E+00;
    COFD[3122] = -2.55356304E-01;
    COFD[3123] = 8.07539912E-03;
    COFD[3124] = -2.22828865E+01;
    COFD[3125] = 4.60177838E+00;
    COFD[3126] = -3.09086396E-01;
    COFD[3127] = 1.06404928E-02;
    COFD[3128] = -2.09328632E+01;
    COFD[3129] = 4.95825742E+00;
    COFD[3130] = -4.08633441E-01;
    COFD[3131] = 1.68756216E-02;
    COFD[3132] = -2.24418859E+01;
    COFD[3133] = 4.74495988E+00;
    COFD[3134] = -3.29202526E-01;
    COFD[3135] = 1.15609276E-02;
    COFD[3136] = -2.06704012E+01;
    COFD[3137] = 4.91327235E+00;
    COFD[3138] = -4.03805612E-01;
    COFD[3139] = 1.67047562E-02;
    COFD[3140] = -2.27673631E+01;
    COFD[3141] = 5.44729221E+00;
    COFD[3142] = -4.52058774E-01;
    COFD[3143] = 1.80072303E-02;
    COFD[3144] = -2.04086296E+01;
    COFD[3145] = 4.40679937E+00;
    COFD[3146] = -2.81870764E-01;
    COFD[3147] = 9.40471579E-03;
    COFD[3148] = -2.03384381E+01;
    COFD[3149] = 4.80405736E+00;
    COFD[3150] = -3.91110136E-01;
    COFD[3151] = 1.62109876E-02;
    COFD[3152] = -1.57635543E+01;
    COFD[3153] = 3.68244512E+00;
    COFD[3154] = -2.57458142E-01;
    COFD[3155] = 1.08919055E-02;
    COFD[3156] = -1.91036850E+01;
    COFD[3157] = 4.51902651E+00;
    COFD[3158] = -3.57827404E-01;
    COFD[3159] = 1.49129882E-02;
    COFD[3160] = -2.06177961E+01;
    COFD[3161] = 4.88872783E+00;
    COFD[3162] = -4.01009699E-01;
    COFD[3163] = 1.65983883E-02;
    COFD[3164] = -2.06065545E+01;
    COFD[3165] = 4.88872783E+00;
    COFD[3166] = -4.01009699E-01;
    COFD[3167] = 1.65983883E-02;
    COFD[3168] = -1.96413638E+01;
    COFD[3169] = 5.19930736E+00;
    COFD[3170] = -4.35079915E-01;
    COFD[3171] = 1.78369576E-02;
    COFD[3172] = -1.90774948E+01;
    COFD[3173] = 4.51902651E+00;
    COFD[3174] = -3.57827404E-01;
    COFD[3175] = 1.49129882E-02;
    COFD[3176] = -2.19972394E+01;
    COFD[3177] = 4.99433306E+00;
    COFD[3178] = -3.72400138E-01;
    COFD[3179] = 1.38127710E-02;
    COFD[3180] = -2.19851974E+01;
    COFD[3181] = 4.94717766E+00;
    COFD[3182] = -3.60880803E-01;
    COFD[3183] = 1.31216550E-02;
    COFD[3184] = -2.19720412E+01;
    COFD[3185] = 4.94717766E+00;
    COFD[3186] = -3.60880803E-01;
    COFD[3187] = 1.31216550E-02;
    COFD[3188] = -2.12691260E+01;
    COFD[3189] = 5.19287816E+00;
    COFD[3190] = -4.34405339E-01;
    COFD[3191] = 1.78137969E-02;
    COFD[3192] = -2.12993900E+01;
    COFD[3193] = 5.19287816E+00;
    COFD[3194] = -4.34405339E-01;
    COFD[3195] = 1.78137969E-02;
    COFD[3196] = -2.13381291E+01;
    COFD[3197] = 5.19899042E+00;
    COFD[3198] = -4.35046771E-01;
    COFD[3199] = 1.78358248E-02;
    COFD[3200] = -2.23500495E+01;
    COFD[3201] = 5.38231713E+00;
    COFD[3202] = -4.48905748E-01;
    COFD[3203] = 1.80649697E-02;
    COFD[3204] = -2.26848257E+01;
    COFD[3205] = 5.42187576E+00;
    COFD[3206] = -4.44522343E-01;
    COFD[3207] = 1.75259406E-02;
    COFD[3208] = -2.27315839E+01;
    COFD[3209] = 5.44654266E+00;
    COFD[3210] = -4.51461923E-01;
    COFD[3211] = 1.79629833E-02;
    COFD[3212] = -2.29185460E+01;
    COFD[3213] = 5.38576397E+00;
    COFD[3214] = -4.34796260E-01;
    COFD[3215] = 1.69232691E-02;
    COFD[3216] = -2.29245190E+01;
    COFD[3217] = 5.40000327E+00;
    COFD[3218] = -4.37643030E-01;
    COFD[3219] = 1.70824209E-02;
    COFD[3220] = -2.30186448E+01;
    COFD[3221] = 5.43600467E+00;
    COFD[3222] = -4.48146589E-01;
    COFD[3223] = 1.77480458E-02;
    COFD[3224] = -2.30269244E+01;
    COFD[3225] = 5.43600467E+00;
    COFD[3226] = -4.48146589E-01;
    COFD[3227] = 1.77480458E-02;
    COFD[3228] = -2.24393654E+01;
    COFD[3229] = 4.74495988E+00;
    COFD[3230] = -3.29202526E-01;
    COFD[3231] = 1.15609276E-02;
    COFD[3232] = -2.23779324E+01;
    COFD[3233] = 4.67548322E+00;
    COFD[3234] = -3.20205650E-01;
    COFD[3235] = 1.11747902E-02;
    COFD[3236] = -2.22828865E+01;
    COFD[3237] = 4.60177838E+00;
    COFD[3238] = -3.09086396E-01;
    COFD[3239] = 1.06404928E-02;
    COFD[3240] = -2.25450838E+01;
    COFD[3241] = 4.83444539E+00;
    COFD[3242] = -3.45616828E-01;
    COFD[3243] = 1.24414117E-02;
    COFD[3244] = -2.03017794E+01;
    COFD[3245] = 4.79345500E+00;
    COFD[3246] = -3.89853526E-01;
    COFD[3247] = 1.61611187E-02;
    COFD[3248] = -2.06668820E+01;
    COFD[3249] = 4.91831350E+00;
    COFD[3250] = -4.04379413E-01;
    COFD[3251] = 1.67265792E-02;
    COFD[3252] = -1.44285949E+01;
    COFD[3253] = 2.99858376E+00;
    COFD[3254] = -1.72232643E-01;
    COFD[3255] = 7.34804765E-03;
    COFD[3256] = -1.69985939E+01;
    COFD[3257] = 3.86929427E+00;
    COFD[3258] = -2.79728651E-01;
    COFD[3259] = 1.17769347E-02;
    COFD[3260] = -2.01962832E+01;
    COFD[3261] = 5.15333309E+00;
    COFD[3262] = -4.30146489E-01;
    COFD[3263] = 1.76625026E-02;
    COFD[3264] = -1.42266828E+01;
    COFD[3265] = 2.91778615E+00;
    COFD[3266] = -1.62159133E-01;
    COFD[3267] = 6.92895656E-03;
    COFD[3268] = -1.15279021E+01;
    COFD[3269] = 2.39753083E+00;
    COFD[3270] = -9.83466917E-02;
    COFD[3271] = 4.32694147E-03;
    COFD[3272] = -1.32446869E+01;
    COFD[3273] = 2.69833584E+00;
    COFD[3274] = -1.33541002E-01;
    COFD[3275] = 5.68346517E-03;
    COFD[3276] = -1.44426226E+01;
    COFD[3277] = 2.99858376E+00;
    COFD[3278] = -1.72232643E-01;
    COFD[3279] = 7.34804765E-03;
    COFD[3280] = -1.44357738E+01;
    COFD[3281] = 2.99858376E+00;
    COFD[3282] = -1.72232643E-01;
    COFD[3283] = 7.34804765E-03;
    COFD[3284] = -1.34540506E+01;
    COFD[3285] = 3.30341833E+00;
    COFD[3286] = -2.10450242E-01;
    COFD[3287] = 8.94670534E-03;
    COFD[3288] = -1.32254608E+01;
    COFD[3289] = 2.69833584E+00;
    COFD[3290] = -1.33541002E-01;
    COFD[3291] = 5.68346517E-03;
    COFD[3292] = -1.90471463E+01;
    COFD[3293] = 4.61496749E+00;
    COFD[3294] = -3.68575530E-01;
    COFD[3295] = 1.53120471E-02;
    COFD[3296] = -1.94158077E+01;
    COFD[3297] = 4.75821765E+00;
    COFD[3298] = -3.85636294E-01;
    COFD[3299] = 1.59920740E-02;
    COFD[3300] = -1.94074945E+01;
    COFD[3301] = 4.75821765E+00;
    COFD[3302] = -3.85636294E-01;
    COFD[3303] = 1.59920740E-02;
    COFD[3304] = -1.51077105E+01;
    COFD[3305] = 3.29803842E+00;
    COFD[3306] = -2.09804181E-01;
    COFD[3307] = 8.92085522E-03;
    COFD[3308] = -1.51305610E+01;
    COFD[3309] = 3.29803842E+00;
    COFD[3310] = -2.09804181E-01;
    COFD[3311] = 8.92085522E-03;
    COFD[3312] = -1.50951030E+01;
    COFD[3313] = 3.28342843E+00;
    COFD[3314] = -2.08038899E-01;
    COFD[3315] = 8.84974116E-03;
    COFD[3316] = -1.65468024E+01;
    COFD[3317] = 3.71875916E+00;
    COFD[3318] = -2.61726361E-01;
    COFD[3319] = 1.10588035E-02;
    COFD[3320] = -1.74553113E+01;
    COFD[3321] = 4.04791428E+00;
    COFD[3322] = -3.01466078E-01;
    COFD[3323] = 1.26603027E-02;
    COFD[3324] = -1.71615857E+01;
    COFD[3325] = 3.90982668E+00;
    COFD[3326] = -2.84689179E-01;
    COFD[3327] = 1.19797702E-02;
    COFD[3328] = -1.81473610E+01;
    COFD[3329] = 4.23397947E+00;
    COFD[3330] = -3.23805200E-01;
    COFD[3331] = 1.35554854E-02;
    COFD[3332] = -1.80453507E+01;
    COFD[3333] = 4.20312679E+00;
    COFD[3334] = -3.20161664E-01;
    COFD[3335] = 1.34118992E-02;
    COFD[3336] = -1.76051897E+01;
    COFD[3337] = 3.97979335E+00;
    COFD[3338] = -2.93193093E-01;
    COFD[3339] = 1.23248995E-02;
    COFD[3340] = -1.76098861E+01;
    COFD[3341] = 3.97979335E+00;
    COFD[3342] = -2.93193093E-01;
    COFD[3343] = 1.23248995E-02;
    COFD[3344] = -2.06657732E+01;
    COFD[3345] = 4.91831350E+00;
    COFD[3346] = -4.04379413E-01;
    COFD[3347] = 1.67265792E-02;
    COFD[3348] = -2.07822044E+01;
    COFD[3349] = 4.92605576E+00;
    COFD[3350] = -4.05242624E-01;
    COFD[3351] = 1.67586331E-02;
    COFD[3352] = -2.09328632E+01;
    COFD[3353] = 4.95825742E+00;
    COFD[3354] = -4.08633441E-01;
    COFD[3355] = 1.68756216E-02;
    COFD[3356] = -2.03017794E+01;
    COFD[3357] = 4.79345500E+00;
    COFD[3358] = -3.89853526E-01;
    COFD[3359] = 1.61611187E-02;
    COFD[3360] = -1.42056656E+01;
    COFD[3361] = 2.91297621E+00;
    COFD[3362] = -1.61544771E-01;
    COFD[3363] = 6.90271324E-03;
};


/*List of specs with small weight, dim NLITE */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetKTDIF EGTRANSETKTDIF
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetKTDIF egtransetktdif
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetKTDIF egtransetktdif_
#endif
void egtransetKTDIF(int* KTDIF) {
    KTDIF[0] = 6;
    KTDIF[1] = 10;
};


/*Poly fits for thermal diff ratios, dim NO*NLITE*KK */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFTD EGTRANSETCOFTD
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFTD egtransetcoftd
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFTD egtransetcoftd_
#endif
void egtransetCOFTD(double* COFTD) {
    COFTD[0] = 2.42391280E-01;
    COFTD[1] = 3.48237221E-04;
    COFTD[2] = -1.24474741E-07;
    COFTD[3] = 1.34873240E-11;
    COFTD[4] = 4.64678698E-01;
    COFTD[5] = 1.14030680E-05;
    COFTD[6] = 6.58534950E-10;
    COFTD[7] = -2.03732653E-13;
    COFTD[8] = 3.87094237E-01;
    COFTD[9] = 1.44275136E-04;
    COFTD[10] = -5.38923894E-08;
    COFTD[11] = 6.22591981E-12;
    COFTD[12] = 1.43936166E-01;
    COFTD[13] = 3.48305588E-04;
    COFTD[14] = -1.21364880E-07;
    COFTD[15] = 1.28856345E-11;
    COFTD[16] = 4.62058474E-01;
    COFTD[17] = 2.48465300E-06;
    COFTD[18] = 4.84734452E-09;
    COFTD[19] = -7.72670148E-13;
    COFTD[20] = 0.00000000E+00;
    COFTD[21] = 0.00000000E+00;
    COFTD[22] = 0.00000000E+00;
    COFTD[23] = 0.00000000E+00;
    COFTD[24] = 4.29427457E-01;
    COFTD[25] = -1.14523829E-05;
    COFTD[26] = 1.15901946E-08;
    COFTD[27] = -1.76574544E-12;
    COFTD[28] = 4.68174867E-01;
    COFTD[29] = 1.14888629E-05;
    COFTD[30] = 6.63489663E-10;
    COFTD[31] = -2.05265506E-13;
    COFTD[32] = 4.66477093E-01;
    COFTD[33] = 1.14472000E-05;
    COFTD[34] = 6.61083606E-10;
    COFTD[35] = -2.04521137E-13;
    COFTD[36] = -1.66147177E-01;
    COFTD[37] = -1.85683391E-05;
    COFTD[38] = 6.10633845E-09;
    COFTD[39] = -7.25879508E-13;
    COFTD[40] = 4.22965568E-01;
    COFTD[41] = -1.12800511E-05;
    COFTD[42] = 1.14157891E-08;
    COFTD[43] = -1.73917506E-12;
    COFTD[44] = 2.61165595E-01;
    COFTD[45] = 2.74720736E-04;
    COFTD[46] = -9.97038966E-08;
    COFTD[47] = 1.09552287E-11;
    COFTD[48] = 2.45552422E-01;
    COFTD[49] = 2.89729034E-04;
    COFTD[50] = -1.04573532E-07;
    COFTD[51] = 1.14298265E-11;
    COFTD[52] = 2.44404444E-01;
    COFTD[53] = 2.88374527E-04;
    COFTD[54] = -1.04084642E-07;
    COFTD[55] = 1.13763911E-11;
    COFTD[56] = 3.73760279E-01;
    COFTD[57] = 4.08535711E-05;
    COFTD[58] = -1.33521153E-08;
    COFTD[59] = 1.58645972E-12;
    COFTD[60] = 3.81177034E-01;
    COFTD[61] = 4.16642536E-05;
    COFTD[62] = -1.36170695E-08;
    COFTD[63] = 1.61794081E-12;
    COFTD[64] = 3.89344085E-01;
    COFTD[65] = 4.00827368E-05;
    COFTD[66] = -1.28679798E-08;
    COFTD[67] = 1.52645237E-12;
    COFTD[68] = 3.87048767E-01;
    COFTD[69] = 1.07274020E-04;
    COFTD[70] = -3.98201643E-08;
    COFTD[71] = 4.65925878E-12;
    COFTD[72] = 3.45838836E-01;
    COFTD[73] = 1.64815034E-04;
    COFTD[74] = -6.14814629E-08;
    COFTD[75] = 7.02051905E-12;
    COFTD[76] = 3.64215419E-01;
    COFTD[77] = 1.44167839E-04;
    COFTD[78] = -5.38598920E-08;
    COFTD[79] = 6.20486961E-12;
    COFTD[80] = 3.36086947E-01;
    COFTD[81] = 2.03523790E-04;
    COFTD[82] = -7.55091594E-08;
    COFTD[83] = 8.52162103E-12;
    COFTD[84] = 3.41784477E-01;
    COFTD[85] = 1.98002827E-04;
    COFTD[86] = -7.35531775E-08;
    COFTD[87] = 8.31904346E-12;
    COFTD[88] = 3.71298285E-01;
    COFTD[89] = 1.62124329E-04;
    COFTD[90] = -6.05409586E-08;
    COFTD[91] = 6.94230156E-12;
    COFTD[92] = 3.72133279E-01;
    COFTD[93] = 1.62488922E-04;
    COFTD[94] = -6.06771060E-08;
    COFTD[95] = 6.95791374E-12;
    COFTD[96] = 2.42292158E-01;
    COFTD[97] = 3.48094816E-04;
    COFTD[98] = -1.24423839E-07;
    COFTD[99] = 1.34818086E-11;
    COFTD[100] = 2.41861537E-01;
    COFTD[101] = 3.54602014E-04;
    COFTD[102] = -1.26621487E-07;
    COFTD[103] = 1.37080380E-11;
    COFTD[104] = 2.33549496E-01;
    COFTD[105] = 3.65335951E-04;
    COFTD[106] = -1.30033054E-07;
    COFTD[107] = 1.40389934E-11;
    COFTD[108] = 2.62756149E-01;
    COFTD[109] = 3.25893055E-04;
    COFTD[110] = -1.17339768E-07;
    COFTD[111] = 1.27963772E-11;
    COFTD[112] = 4.62401641E-01;
    COFTD[113] = 1.96551308E-06;
    COFTD[114] = 5.10373738E-09;
    COFTD[115] = -8.08460072E-13;
    COFTD[116] = -5.87585797E-02;
    COFTD[117] = 6.03873322E-04;
    COFTD[118] = -1.89767502E-07;
    COFTD[119] = 1.88165158E-11;
    COFTD[120] = 3.04089430E-01;
    COFTD[121] = 2.65271032E-04;
    COFTD[122] = -9.70946444E-08;
    COFTD[123] = 1.07635556E-11;
    COFTD[124] = 1.26448444E-01;
    COFTD[125] = 4.59133576E-04;
    COFTD[126] = -1.57182261E-07;
    COFTD[127] = 1.64763301E-11;
    COFTD[128] = -5.62840660E-02;
    COFTD[129] = 5.52051117E-04;
    COFTD[130] = -1.73184146E-07;
    COFTD[131] = 1.71538798E-11;
    COFTD[132] = 3.19162811E-01;
    COFTD[133] = 2.41892210E-04;
    COFTD[134] = -8.90518568E-08;
    COFTD[135] = 9.93959623E-12;
    COFTD[136] = 1.66147177E-01;
    COFTD[137] = 1.85683391E-05;
    COFTD[138] = -6.10633845E-09;
    COFTD[139] = 7.25879508E-13;
    COFTD[140] = 3.40789020E-01;
    COFTD[141] = 1.86618018E-04;
    COFTD[142] = -6.94241287E-08;
    COFTD[143] = 7.87383759E-12;
    COFTD[144] = 3.05228036E-01;
    COFTD[145] = 2.66264290E-04;
    COFTD[146] = -9.74581973E-08;
    COFTD[147] = 1.08038577E-11;
    COFTD[148] = 3.04675603E-01;
    COFTD[149] = 2.65782378E-04;
    COFTD[150] = -9.72818075E-08;
    COFTD[151] = 1.07843038E-11;
    COFTD[152] = 0.00000000E+00;
    COFTD[153] = 0.00000000E+00;
    COFTD[154] = 0.00000000E+00;
    COFTD[155] = 0.00000000E+00;
    COFTD[156] = 3.38244124E-01;
    COFTD[157] = 1.85224419E-04;
    COFTD[158] = -6.89056931E-08;
    COFTD[159] = 7.81503847E-12;
    COFTD[160] = 5.03120986E-03;
    COFTD[161] = 5.43289883E-04;
    COFTD[162] = -1.76510875E-07;
    COFTD[163] = 1.78710494E-11;
    COFTD[164] = -3.16410304E-02;
    COFTD[165] = 5.63965226E-04;
    COFTD[166] = -1.79753402E-07;
    COFTD[167] = 1.79812649E-11;
    COFTD[168] = -3.15672407E-02;
    COFTD[169] = 5.62650008E-04;
    COFTD[170] = -1.79334201E-07;
    COFTD[171] = 1.79393310E-11;
    COFTD[172] = 2.23754436E-01;
    COFTD[173] = 3.08293806E-04;
    COFTD[174] = -1.10427107E-07;
    COFTD[175] = 1.19868185E-11;
    COFTD[176] = 2.25931288E-01;
    COFTD[177] = 3.11293120E-04;
    COFTD[178] = -1.11501425E-07;
    COFTD[179] = 1.21034353E-11;
    COFTD[180] = 2.31556629E-01;
    COFTD[181] = 3.09931883E-04;
    COFTD[182] = -1.11173779E-07;
    COFTD[183] = 1.20831276E-11;
    COFTD[184] = 1.57705710E-01;
    COFTD[185] = 4.13918237E-04;
    COFTD[186] = -1.43688545E-07;
    COFTD[187] = 1.52134514E-11;
    COFTD[188] = 9.15485422E-02;
    COFTD[189] = 4.74751339E-04;
    COFTD[190] = -1.60420888E-07;
    COFTD[191] = 1.66648306E-11;
    COFTD[192] = 1.15932628E-01;
    COFTD[193] = 4.55157930E-04;
    COFTD[194] = -1.55343432E-07;
    COFTD[195] = 1.62485973E-11;
    COFTD[196] = 6.03774138E-02;
    COFTD[197] = 5.11860295E-04;
    COFTD[198] = -1.70557374E-07;
    COFTD[199] = 1.75527891E-11;
    COFTD[200] = 6.66780752E-02;
    COFTD[201] = 5.07632941E-04;
    COFTD[202] = -1.69602339E-07;
    COFTD[203] = 1.74854347E-11;
    COFTD[204] = 1.05547641E-01;
    COFTD[205] = 4.76273236E-04;
    COFTD[206] = -1.61709766E-07;
    COFTD[207] = 1.68539304E-11;
    COFTD[208] = 1.05666055E-01;
    COFTD[209] = 4.76807568E-04;
    COFTD[210] = -1.61891188E-07;
    COFTD[211] = 1.68728388E-11;
    COFTD[212] = -5.87465680E-02;
    COFTD[213] = 6.03749875E-04;
    COFTD[214] = -1.89728709E-07;
    COFTD[215] = 1.88126693E-11;
    COFTD[216] = -5.79524543E-02;
    COFTD[217] = 6.06256451E-04;
    COFTD[218] = -1.90635813E-07;
    COFTD[219] = 1.89099601E-11;
    COFTD[220] = -6.51750002E-02;
    COFTD[221] = 6.10372090E-04;
    COFTD[222] = -1.91133211E-07;
    COFTD[223] = 1.89104839E-11;
    COFTD[224] = -3.34113480E-02;
    COFTD[225] = 5.91141763E-04;
    COFTD[226] = -1.88390203E-07;
    COFTD[227] = 1.88436542E-11;
    COFTD[228] = 3.20290850E-01;
    COFTD[229] = 2.40548948E-04;
    COFTD[230] = -8.85886011E-08;
    COFTD[231] = 9.89229207E-12;
};

/* End of file  */




#if 0




\\
\\
\\  This is the mechanism file
\\
\\
!Mechanism imported from chem.bin                                                
ELEMENTS
 H C O N
END
SPECIES
NC7H16   O2       CO2      H2O
CO       H2       OH       H2O2     HO2
H        O        CH3O     CH2O     HCO
CH2      CH3      CH4      C2H3     C2H4
C2H5     C3H4     C3H5     C3H6     C3H7
C7H15-2  C7H15O2  C7KET12  C5H11CO  N2
END

TRANS ALL
NC7H16             2   564.030     6.004     0.000     0.000     1.000          
O2                 1   107.400     3.458     0.000     1.600     3.800          
N2                 1    97.530     3.621     0.000     1.760     4.000          
H2                 1    38.000     2.920     0.000     0.790   280.000          
H2O                2   572.400     2.605     1.844     0.000     4.000          
H2O2               2   107.400     3.458     0.000     0.000     3.800          
HO2                2   107.400     3.458     0.000     0.000     1.000          
CO2                1   244.000     3.763     0.000     2.650     2.100          
CO                 1    98.100     3.650     0.000     1.950     1.800          
OH                 1    80.000     2.750     0.000     0.000     0.000          
H                  0   145.000     2.050     0.000     0.000     0.000          
O                  0    80.000     2.750     0.000     0.000     0.000          
CH3O               2   417.000     3.690     1.700     0.000     2.000          
CH2O               2   498.000     3.590     0.000     0.000     2.000          
HCO                2   498.000     3.590     0.000     0.000     0.000          
CH2                1   144.000     3.800     0.000     0.000     0.000          
CH3                1   144.000     3.800     0.000     0.000     0.000          
CH4                2   141.400     3.746     0.000     2.600    13.000          
C2H3               2   209.000     4.100     0.000     0.000     1.000          
C2H4               2   280.800     3.971     0.000     0.000     1.500          
C2H5               2   252.300     4.302     0.000     0.000     1.500          
C3H4               1   324.800     4.290     0.000     0.000     1.000          
C3H5               2   316.000     4.220     0.000     0.000     1.000          
C3H6               2   266.800     4.982     0.000     0.000     1.000          
C3H7               2   266.800     4.982     0.000     0.000     1.000          
C7H15-2            2   564.030     6.004     0.000     0.000     1.000          
C7H15O2            2   561.000     6.317     1.700     0.000     1.000          
C7KET12            2   581.300     6.506     2.000     0.000     1.000          
C5H11CO            2   498.600     6.009     2.000     0.000     1.000          
END

REACTIONS
NC7H16+H<=>H2+C7H15-2                  4.380000E+07       2.00      4.759484E+03
NC7H16+OH<=>H2O+C7H15-2                9.700000E+09       1.30      1.689817E+03
NC7H16+HO2<=>H2O2+C7H15-2              1.650000E+13       0.00      1.694816E+04
NC7H16+O2<=>HO2+C7H15-2                2.000000E+15       0.00      4.737486E+04
O2+C7H15-2<=>C7H15O2                   1.560000E+12       0.00      0.000000E+00
O2+C7H15O2<=>OH+C7KET12                4.500000E+14       0.00      1.823073E+04
C7KET12<=>OH+CH2O+C5H11CO              9.530000E+14       0.00      4.109554E+04
C5H11CO<=>CO+C2H4+C3H7                 9.840000E+15       0.00      4.019564E+04
C7H15-2<=>C2H4+C2H5+C3H6               7.045000E+14       0.00      3.459625E+04
C3H7<=>CH3+C2H4                        9.600000E+13       0.00      3.094664E+04
C3H7<=>H+C3H6                          1.250000E+14       0.00      3.689600E+04
CH3+C3H6<=>CH4+C3H5                    9.000000E+12       0.00      8.479080E+03
O2+C3H5<=>HO2+C3H4                     6.000000E+11       0.00      9.998915E+03
OH+C3H4<=>CH2O+C2H3                    1.000000E+12       0.00      0.000000E+00
OH+C3H4<=>HCO+C2H4                     1.000000E+12       0.00      0.000000E+00
HO2+CH3<=>OH+CH3O                      5.000000E+13       0.00      0.000000E+00
OH+CH3<=>H2O+CH2                       7.500000E+06       2.00      4.999458E+03
OH+CH2<=>H+CH2O                        2.500000E+13       0.00      0.000000E+00
O2+CH2<=>OH+HCO                        4.300000E+10       0.00     -4.999458E+02
O2+CH2<=>CO2+H2                        6.900000E+11       0.00      4.999458E+02
O2+CH2<=>H2O+CO                        2.000000E+10       0.00     -9.998915E+02
O2+CH2<=>O+CH2O                        5.000000E+13       0.00      8.999024E+03
O2+CH2<=>CO2+ 2H                       1.600000E+12       0.00      9.998915E+02
O2+CH2<=>CO+OH+H                       8.600000E+10       0.00     -4.999458E+02
CO+CH3O<=>CO2+CH3                      1.570000E+14       0.00      1.179872E+04
CO+OH<=>CO2+H                          8.987000E+07       1.38      5.232309E+03
OH+O<=>O2+H                            4.000000E+14      -0.50      0.000000E+00
HO2+H<=> 2OH                           1.700000E+14       0.00      8.749051E+02
 2OH<=>H2O+O                           6.000000E+08       1.30      0.000000E+00
O2+H+M<=>HO2+M                         3.600000E+17      -0.72      0.000000E+00
        CO2 /   5.00/H2O /  21.00/CO /   2.00/H2 /   3.30/
H2O2+M<=> 2OH+M                        1.000000E+16       0.00      4.549506E+04
        CO2 /   5.00/H2O /  21.00/CO /   2.00/H2 /   3.30/
H2+OH<=>H2O+H                          1.170000E+09       1.30      3.625607E+03
 2HO2<=>O2+H2O2                        3.000000E+12       0.00      0.000000E+00
OH+CH2O<=>H2O+HCO                      5.563000E+10       1.09     -7.650870E+01
HO2+CH2O<=>H2O2+HCO                    3.000000E+12       0.00      7.999132E+03
O2+HCO<=>CO+HO2                        3.300000E+13      -0.40      0.000000E+00
HCO+M<=>CO+H+M                         1.591000E+18       0.95      5.670618E+04
        
CH3O+CH3<=>CH2O+CH4                    4.300000E+14       0.00      0.000000E+00
OH+C2H4<=>CH2O+CH3                     6.000000E+13       0.00      9.598959E+02
OH+C2H4<=>H2O+C2H3                     8.020000E+13       0.00      5.954354E+03
O2+C2H3<=>CH2O+HCO                     4.000000E+12       0.00     -2.499729E+02
HCO+C2H3<=>CO+C2H4                     6.034000E+13       0.00      0.000000E+00
O2+C2H5<=>HO2+C2H4                     2.000000E+10       0.00     -2.199761E+03
O2+CH4<=>HO2+CH3                       7.900000E+13       0.00      5.599393E+04
OH+HO2<=>O2+H2O                        7.500000E+12       0.00      0.000000E+00
O2+CH3<=>OH+CH2O                       3.800000E+11       0.00      8.999024E+03
H+CH4<=>H2+CH3                         6.600000E+08       1.60      1.083882E+04
OH+CH4<=>H2O+CH3                       1.600000E+06       2.10      2.459733E+03
O+CH4<=>OH+CH3                         1.020000E+09       1.50      8.603067E+03
HO2+CH4<=>H2O2+CH3                     9.000000E+11       0.00      1.869797E+04
CH2+CH4<=> 2CH3                        4.000000E+12       0.00     -5.699382E+02
C3H6<=>CH3+C2H3                        3.150000E+15       0.00      8.549072E+04
END

\\
\\
\\  This is the therm file
\\
\\
!Mechanism imported from chem.bin                                                
THERMO
    300.000  1000.000  5000.000
NC7H16                  H  16C   7          G   300.000  5000.000  1391.000    1
 2.22148969E+01 3.47675750E-02-1.18407129E-05 1.83298478E-09-1.06130266E-13    2
-3.42760081E+04-9.23040196E+01-1.26836187E+00 8.54355820E-02-5.25346786E-05    3
 1.62945721E-08-2.02394925E-12-2.56586565E+04 3.53732912E+01                   4
O2                      O   2               G   300.000  5000.000  1000.000    1
 3.69757800E+00 6.13519700E-04-1.25884200E-07 1.77528100E-11-1.13643500E-15    2
-1.23393000E+03 3.18916600E+00 3.21293600E+00 1.12748600E-03-5.75615000E-07    3
 1.31387700E-09-8.76855400E-13-1.00524900E+03 6.03473800E+00                   4
N2                      N   2               G   300.000  5000.000  1000.000    1
 2.92664000E+00 1.48797700E-03-5.68476100E-07 1.00970400E-10-6.75335100E-15    2
-9.22797700E+02 5.98052800E+00 3.29867700E+00 1.40824000E-03-3.96322200E-06    3
 5.64151500E-09-2.44485500E-12-1.02090000E+03 3.95037200E+00                   4
CO2                     C   1O   2          G   300.000  5000.000  1000.000    1
 4.45362300E+00 3.14016900E-03-1.27841100E-06 2.39399700E-10-1.66903300E-14    2
-4.89669600E+04-9.55395900E-01 2.27572500E+00 9.92207200E-03-1.04091100E-05    3
 6.86668700E-09-2.11728000E-12-4.83731400E+04 1.01884900E+01                   4
H2O                     H   2O   1          G   300.000  5000.000  1000.000    1
 2.67214600E+00 3.05629300E-03-8.73026000E-07 1.20099600E-10-6.39161800E-15    2
-2.98992100E+04 6.86281700E+00 3.38684200E+00 3.47498200E-03-6.35469600E-06    3
 6.96858100E-09-2.50658800E-12-3.02081100E+04 2.59023300E+00                   4
CO                      C   1O   1          G   300.000  5000.000  1000.000    1
 3.02507800E+00 1.44268900E-03-5.63082800E-07 1.01858100E-10-6.91095200E-15    2
-1.42683500E+04 6.10821800E+00 3.26245200E+00 1.51194100E-03-3.88175500E-06    3
 5.58194400E-09-2.47495100E-12-1.43105400E+04 4.84889700E+00                   4
H2                      H   2               G   300.000  5000.000  1000.000    1
 2.99142300E+00 7.00064400E-04-5.63382900E-08-9.23157800E-12 1.58275200E-15    2
-8.35034000E+02-1.35511000E+00 3.29812400E+00 8.24944200E-04-8.14301500E-07    3
-9.47543400E-11 4.13487200E-13-1.01252100E+03-3.29409400E+00                   4
OH                      H   1O   1          G   300.000  5000.000  1000.000    1
 2.88273000E+00 1.01397400E-03-2.27687700E-07 2.17468400E-11-5.12630500E-16    2
 3.88688800E+03 5.59571200E+00 3.63726600E+00 1.85091000E-04-1.67616500E-06    3
 2.38720300E-09-8.43144200E-13 3.60678200E+03 1.35886000E+00                   4
H2O2                    H   2O   2          G   300.000  5000.000  1000.000    1
 4.57316700E+00 4.33613600E-03-1.47468900E-06 2.34890400E-10-1.43165400E-14    2
-1.80069600E+04 5.01137000E-01 3.38875400E+00 6.56922600E-03-1.48501300E-07    3
-4.62580600E-09 2.47151500E-12-1.76631500E+04 6.78536300E+00                   4
HO2                     H   1O   2          G   300.000  5000.000  1000.000    1
 4.07219100E+00 2.13129600E-03-5.30814500E-07 6.11226900E-11-2.84116500E-15    2
-1.57972700E+02 3.47602900E+00 2.97996300E+00 4.99669700E-03-3.79099700E-06    3
 2.35419200E-09-8.08902400E-13 1.76227400E+02 9.22272400E+00                   4
H                       H   1               G   300.000  5000.000  1000.000    1
 2.50000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
 2.54716300E+04-4.60117600E-01 2.50000000E+00 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00 2.54716300E+04-4.60117600E-01                   4
O                       O   1               G   300.000  5000.000  1000.000    1
 2.54206000E+00-2.75506200E-05-3.10280300E-09 4.55106700E-12-4.36805200E-16    2
 2.92308000E+04 4.92030800E+00 2.94642900E+00-1.63816600E-03 2.42103200E-06    3
-1.60284300E-09 3.89069600E-13 2.91476400E+04 2.96399500E+00                   4
CH3O                    H   3C   1O   1     G   300.000  5000.000  1000.000    1
 3.77080000E+00 7.87149700E-03-2.65638400E-06 3.94443100E-10-2.11261600E-14    2
 1.27832500E+02 2.92957500E+00 2.10620400E+00 7.21659500E-03 5.33847200E-06    3
-7.37763600E-09 2.07561100E-12 9.78601100E+02 1.31521800E+01                   4
CH2O                    H   2C   1O   1     G   300.000  5000.000  1000.000    1
 2.99560600E+00 6.68132100E-03-2.62895500E-06 4.73715300E-10-3.21251700E-14    2
-1.53203700E+04 6.91257200E+00 1.65273100E+00 1.26314400E-02-1.88816800E-05    3
 2.05003100E-08-8.41323700E-12-1.48654000E+04 1.37848200E+01                   4
HCO                     H   1C   1O   1     G   300.000  5000.000  1000.000    1
 3.55727100E+00 3.34557300E-03-1.33500600E-06 2.47057300E-10-1.71385100E-14    2
 3.91632400E+03 5.55229900E+00 2.89833000E+00 6.19914700E-03-9.62308400E-06    3
 1.08982500E-08-4.57488500E-12 4.15992200E+03 8.98361400E+00                   4
CH2                     H   2C   1          G   300.000  5000.000  1000.000    1
 3.63640800E+00 1.93305700E-03-1.68701600E-07-1.00989900E-10 1.80825600E-14    2
 4.53413400E+04 2.15656100E+00 3.76223700E+00 1.15981900E-03 2.48958500E-07    3
 8.80083600E-10-7.33243500E-13 4.53679100E+04 1.71257800E+00                   4
CH3                     H   3C   1          G   300.000  5000.000  1000.000    1
 2.84405200E+00 6.13797400E-03-2.23034500E-06 3.78516100E-10-2.45215900E-14    2
 1.64378100E+04 5.45269700E+00 2.43044300E+00 1.11241000E-02-1.68022000E-05    3
 1.62182900E-08-5.86495300E-12 1.64237800E+04 6.78979400E+00                   4
CH4                     H   4C   1          G   300.000  5000.000  1000.000    1
 1.68347900E+00 1.02372400E-02-3.87512900E-06 6.78558500E-10-4.50342300E-14    2
-1.00807900E+04 9.62339500E+00 7.78741500E-01 1.74766800E-02-2.78340900E-05    3
 3.04970800E-08-1.22393100E-11-9.82522900E+03 1.37221900E+01                   4
C2H3                    H   3C   2          G   300.000  5000.000  1000.000    1
 5.93346800E+00 4.01774600E-03-3.96674000E-07-1.44126700E-10 2.37864400E-14    2
 3.18543500E+04-8.53031300E+00 2.45927600E+00 7.37147600E-03 2.10987300E-06    3
-1.32164200E-09-1.18478400E-12 3.33522500E+04 1.15562000E+01                   4
C2H4                    H   4C   2          G   300.000  5000.000  1000.000    1
 3.52841900E+00 1.14851800E-02-4.41838500E-06 7.84460100E-10-5.26684800E-14    2
 4.42828900E+03 2.23038900E+00-8.61488000E-01 2.79616300E-02-3.38867700E-05    3
 2.78515200E-08-9.73787900E-12 5.57304600E+03 2.42114900E+01                   4
C2H5                    H   5C   2          G   300.000  5000.000  1000.000    1
 7.19048000E+00 6.48407700E-03-6.42806500E-07-2.34787900E-10 3.88087700E-14    2
 1.06745500E+04-1.47808900E+01 2.69070200E+00 8.71913300E-03 4.41983900E-06    3
 9.33870300E-10-3.92777300E-12 1.28704000E+04 1.21382000E+01                   4
C3H4                    H   4C   3          G   300.000  5000.000  1000.000    1
 6.31694869E+00 1.11336262E-02-3.96289018E-06 6.35633775E-10-3.78749885E-14    2
 2.01174617E+04-1.09718862E+01 2.61307487E+00 1.21223371E-02 1.85405400E-05    3
-3.45258475E-08 1.53353389E-11 2.15415642E+04 1.02503319E+01                   4
C3H5                    H   5C   3          G   300.000  5000.000  1000.000    1
 6.54761132E+00 1.33152246E-02-4.78333100E-06 7.71949814E-10-4.61930808E-14    2
 1.72714707E+04-9.27486841E+00 3.78794693E+00 9.48414335E-03 2.42343368E-05    3
-3.65604010E-08 1.48592356E-11 1.86261218E+04 7.82822499E+00                   4
C3H6                    H   6C   3          G   300.000  5000.000  1388.000    1
 8.01595958E+00 1.37023634E-02-4.66249733E-06 7.21254402E-10-4.17370126E-14    2
-1.87821271E+03-2.00160668E+01 3.94615444E-01 2.89107662E-02-1.54886808E-05    3
 3.88814209E-09-3.37890352E-13 1.06688164E+03 2.19003736E+01                   4
C3H7                    H   7C   3          G   300.000  5000.000  1000.000    1
 7.70269870E+00 1.60442030E-02-5.28332200E-06 7.62985900E-10-3.93922840E-14    2
 8.29843360E+03-1.54801800E+01 1.05155180E+00 2.59919800E-02 2.38005400E-06    3
-1.96095690E-08 9.37324700E-12 1.06318630E+04 2.11225590E+01                   4
C7H15-2                 H  15C   7          G   300.000  5000.000  1382.000    1
 2.16368842E+01 3.23324804E-02-1.09273807E-05 1.68357060E-09-9.71774091E-14    2
-1.05873616E+04-8.52209653E+01-3.79155767E-02 7.56726570E-02-4.07473634E-05    3
 9.32678943E-09-4.92360745E-13-2.35605303E+03 3.37321506E+01                   4
C7H15O2                 H  15C   7O   2     G   300.000  5000.000  1390.000    1
 2.49023689E+01 3.50716920E-02-1.20440306E-05 1.87464822E-09-1.08947791E-13    2
-2.82976050E+04-9.73923542E+01 2.37499334E+00 8.34651906E-02-5.13897320E-05    3
 1.64217662E-08-2.19505216E-12-1.99237961E+04 2.53067342E+01                   4
C7KET12                 H  14C   7O   3     G   300.000  5000.000  1396.000    1
 2.97472906E+01 3.06622294E-02-1.05563590E-05 1.64627343E-09-9.58171675E-14    2
-5.66856828E+04-1.22432490E+02 5.82433697E-01 1.01207869E-01-7.65855996E-05    3
 3.00738606E-08-4.82902792E-12-4.68054419E+04 3.33331449E+01                   4
C5H11CO                 H  11C   6O   1     G   300.000  5000.000  1383.000    1
 1.94783812E+01 2.50466029E-02-8.54861346E-06 1.32557944E-09-7.68503296E-14    2
-2.07923937E+04-7.21995578E+01 2.14479069E+00 6.17863563E-02-3.74134690E-05    3
 1.13283795E-08-1.36917698E-12-1.43451172E+04 2.23128045E+01                   4
END

\\
\\
\\  This is the tran file
\\
\\

#endif
