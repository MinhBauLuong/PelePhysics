# AMREX_HOME defines the directory in which we will find all the BoxLib code
# If you set AMREX_HOME as an environment variable, this line will be ignored

DEBUG      = FALSE

DIM        = 2

COMP       = gcc
FCOMP	   = mpif90

USE_MPI    = FALSE
USE_OMP    = FALSE

FUEGO_GAS  = TRUE

USE_CVODE = FALSE
#USE_PETSC = TRUE
USE_SUNDIALS3=TRUE
USE_ARKODE=FALSE

PRECISION  = DOUBLE
EBASE     = main

ifeq ($(FUEGO_GAS), TRUE)
  Eos_dir         = Fuego
  #Chemistry_Model = chem-H
  #Chemistry_Model = LiDryer
  #Chemistry_Model = BurkeDryer_mod
  #Chemistry_Model = BurkeDryer
  #Chemistry_Model = air
  #Chemistry_Model = grimech12
  #Chemistry_Model = grimech30
  #Chemistry_Model = grimech30-noArN
  #Chemistry_Model = drm19
  #Chemistry_Model = dodecane_wang
  Chemistry_Model = ERC_nheptane
  #Chemistry_Model = LLNL_GasolineSurr_RedI
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
