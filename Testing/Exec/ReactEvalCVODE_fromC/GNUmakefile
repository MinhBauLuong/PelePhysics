# AMREX_HOME defines the directory in which we will find all the BoxLib code
# If you set AMREX_HOME as an environment variable, this line will be ignored

DEBUG      = TRUE
USE_MPI    = FALSE
USE_OMP    = FALSE
USE_CUDA   = FALSE
DEBUG       = TRUE
COMP       = gnu
##FCOMP	   = pgfortran
DIM        = 2
PRECISION  = DOUBLE
EBASE     = main

FUEGO_GAS  = TRUE

#USE_CVODE = FALSE
#USE_SUNDIALS3=TRUE
USE_SUNDIALS_3x4x=TRUE
USE_ARKODE=FALSE
USE_CUDA_CVODE = FALSE
CVODE_LIB_DIR ?= ${CVODE_LIB}


ifeq ($(FUEGO_GAS), TRUE)
  Eos_dir         = Fuego
  #Chemistry_Model = ERC_nheptane
  #Chemistry_Model = drm19
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
