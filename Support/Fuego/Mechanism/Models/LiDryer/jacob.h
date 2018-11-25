
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void eval_jacob (const double, const double, double *, const double * restrict, double * restrict);
void eval_rxn_rates (const double, const double * restrict, double * restrict, double * restrict);
void get_rxn_pres_mod (const double, const double, const double * restrict, double * restrict);
void eval_spec_rates (const double * restrict, const double * restrict, const double * restrict, double * restrict, double * restrict);

void eval_conc (const double, const double, const double * restrict, double * restrict, double * restrict, double * restrict, double * restrict);
//void eval_conc_rho (const double, const double, const double * restrict, double * restrict, double * restrict, double * restrict, double * restrict);
void eval_h (const double, double * restrict);
//void eval_u (const double, double * restrict);
//void eval_cv (const double, double * restrict);
void eval_cp (const double, double * restrict);
