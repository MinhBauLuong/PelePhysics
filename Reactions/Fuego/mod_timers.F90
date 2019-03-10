MODULE mod_timers

   IMPLICIT NONE

   TYPE timer_t
      CHARACTER(LEN=256) :: name
      REAL :: start_time, stop_time
      REAL :: total_time   
      REAL :: fraction
      CONTAINS
         PROCEDURE :: init => timer_init
         PROCEDURE :: print => timer_print
         PROCEDURE :: add => timer_add
         PROCEDURE :: start => timer_start
         PROCEDURE :: stop => timer_stop
         PROCEDURE :: reset => timer_reset
         PROCEDURE :: norm => timer_norm
   END TYPE timer_t

   TYPE(timer_t) :: t_total
   TYPE(timer_t) :: t_init
   TYPE(timer_t) :: t_bechem
   TYPE(timer_t) :: t_ck
   TYPE(timer_t) :: t_eos
   TYPE(timer_t) :: t_readData
   TYPE(timer_t) :: t_ReInit
   TYPE(timer_t) :: t_CVODE
   TYPE(timer_t) :: t_AJac
   TYPE(timer_t) :: t_ckJac


   CONTAINS

SUBROUTINE timer_init(this, newname)      
   CLASS(timer_t), INTENT(INOUT) :: this
   CHARACTER(LEN=*), INTENT(IN) :: newname
   this%name = TRIM(newname)
   this%total_time = 0.0d0
   this%fraction = 0.0d0
END SUBROUTINE timer_init

SUBROUTINE timer_print(this)
   CLASS(timer_t), INTENT(IN) :: this
   WRITE(*,'(A5,A30,ES10.4,A,F6.2,A)') " --> ",ADJUSTR(TRIM(this%name)//" time (s): "), this%total_time, " (", this%fraction, "%)"
END SUBROUTINE timer_print

SUBROUTINE timer_add(this, delta_t)
   CLASS(timer_t), INTENT(INOUT) :: this
   REAL :: delta_t
   this%total_time = this%total_time + delta_t
END SUBROUTINE timer_add

SUBROUTINE timer_start(this)
   CLASS(timer_t), INTENT(INOUT) :: this
   CALL cpu_time(this%start_time)
END SUBROUTINE timer_start

SUBROUTINE timer_stop(this)
   CLASS(timer_t), INTENT(INOUT) :: this
   CALL cpu_time(this%stop_time)
   this%total_time = this%total_time + ( this%stop_time - this%start_time )
   this%start_time = 0.0d0
END SUBROUTINE timer_stop

SUBROUTINE timer_reset(this)
   CLASS(timer_t), INTENT(INOUT) :: this
   this%total_time = 0.0d0
   this%start_time = 0.0d0
   this%fraction = 0.0d0
END SUBROUTINE timer_reset

SUBROUTINE timer_norm(this, that)
   CLASS(timer_t), INTENT(INOUT) :: this
   TYPE(timer_t), INTENT(IN) :: that
   this%fraction = this%total_time / that%total_time * 100.0d0
END SUBROUTINE timer_norm

END MODULE mod_timers
