# AMREX_HOME defines the directory in which we will find all the BoxLib code
# If you set AMREX_HOME as an environment variable, this line will be ignored

DEBUG      = FALSE
USE_MPI    = FALSE
USE_OMP    = FALSE

USE_CUDA   = FALSE

COMP       = gcc
FCOMP	   = gfortran
DIM        = 2
PRECISION  = DOUBLE
EBASE     = main

FUEGO_GAS  = TRUE

#USE_CVODE = FALSE
#USE_SUNDIALS3=TRUE
USE_SUNDIALS_3x4x=TRUE
USE_DVODE=TRUE

USE_CUDA_CVODE = FALSE

USE_SDC_FORTRAN = FALSE
ifeq ($(USE_SDC_FORTRAN), TRUE)
    DEFINES  += -DUSE_SDC_F
endif

#USE_PYJAC = FALSE
#ifeq ($(USE_PYJAC), TRUE)
#    #CPPFLAGS += -DUSE_PYJAC 
#    DEFINES  += -DUSE_PYJAC 
#endif

CVODE_LIB_DIR = /home/amfelden/CODES/PeleLM/Source/ChemDriver/sundials-3.2.1/instdir/lib/ 


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
