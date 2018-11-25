PRECISION  = DOUBLE

DEBUG      = FALSE

DIM        = 2

USE_MPI    = FALSE
USE_OMP    = FALSE

USE_CUDA           = FALSE
USE_CVODE_CUSOLVER = FALSE

USE_SUNDIALS_3x4x=TRUE
#CVODE_LIB_DIR ?= ${CVODE_LIB}

ifeq ($(USE_CUDA), TRUE)
    USE_CUDA_CVODE = TRUE
    COMP           = pgi
else
    USE_CUDA_CVODE = FALSE
    COMP           = gnu
endif

EBASE     = main

FUEGO_GAS  = TRUE

ifeq ($(FUEGO_GAS), TRUE)
  Eos_dir         = Fuego
  #Chemistry_Model = ERC_nheptane
  Chemistry_Model = LiDryer
  #Chemistry_Model = dodecane_wang
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
include Make.CVODE
