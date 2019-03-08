# AMREX_HOME defines the directory in which we will find all the BoxLib code
# If you set AMREX_HOME as an environment variable, this line will be ignored

DEBUG      = FALSE
USE_MPI    = FALSE
USE_OMP    = FALSE

USE_CUDA   = TRUE
ifeq ($(USE_CUDA), TRUE)
    USE_CUDA_CVODE = TRUE
    USE_CVODE_CUSOLVER = FALSE
    COMP           = gcc
else
    USE_CUDA_CVODE = FALSE
    USE_CVODE_CUSOLVER = FALSE
    COMP           = gcc
endif

DIM        = 2
PRECISION  = DOUBLE
EBASE     = main

FUEGO_GAS  = TRUE

#######################
#include ./Make.cvodeFiniterface 
include Make.CVODE

USE_DVODE=FALSE
ifeq ($(USE_DVODE), TRUE)
    DEFINES  += -DUSE_DVODE
endif

USE_SUNDIALS_3x4x=TRUE
#CVODE_LIB_DIR = /home/amfelden/CODES/PeleLM/Source/ChemDriver/sundials-3.2.1/instdir/lib/ 
CVODE_LIB_DIR = /home/amfelden/CODES/PeleLM/Source/ChemDriver/cvode-4.0.2/instdir/lib/

USE_KLU = FALSE
ifeq ($(USE_KLU), TRUE)
    DEFINES  += -DUSE_KLU
endif

USE_FORTRAN_CVODE=FALSE
ifeq ($(USE_FORTRAN_CVODE), TRUE)
    DEFINES  += -DUSE_FORTRAN_CVODE
    USE_PYJAC = FALSE
    ifeq ($(USE_PYJAC), TRUE)
        DEFINES  += -DUSE_PYJAC
    endif
else
    USE_PYJAC = FALSE
endif

USE_SDC_FORTRAN = FALSE
ifeq ($(USE_SDC_FORTRAN), TRUE)
    DEFINES  += -DUSE_SDC_F
endif

#######################

ifeq ($(FUEGO_GAS), TRUE)
  Eos_dir         = Fuego
  #Chemistry_Model = ERC_nheptane
  #Chemistry_Model = LiDryer
  Chemistry_Model = dodecane_wang
  Reactions_dir   = Fuego
  Transport_dir   = Simple
else
  Eos_dir       = GammaLaw
  Reactions_dir = Null
  Transport_dir = Constant
endif

Bpack   := ./Make.package
Blocs   := .

include $(PELE_PHYSICS_HOME)/Testing/Exec/Make.PelePhysics
include $(AMREX_HOME)/Src/Amr/Make.package 
