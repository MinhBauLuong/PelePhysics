#ifndef _main_F_H_
#define _main_F_H_

#include <AMReX_BLFort.H>

extern "C"
{
    void extern_init
    (const int* name, const int* namlen);

    void extern_close();

    void get_num_spec(int* nspec);

    void initialize_data_byhand
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(massfraction));

    void initialize_data
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(massfraction),
     const BL_FORT_FAB_ARG_3D(rY_src_ext));

    void read_data_from_txt
    (const int* name, const int* namlen);

    void react_state
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(rYo),
     const BL_FORT_FAB_ARG_3D(rYn),
     const BL_FORT_FAB_ARG_3D(rY_src_ext),
     const BL_FORT_IFAB_ARG_3D(mask),
     const BL_FORT_FAB_ARG_3D(cost),
     const Real* time, const int* ndt, const Real* dt);
}

#endif // _main_F_H_
