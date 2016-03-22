module sat_f
   use iso_c_binding
   implicit none
   interface
      function sat_f_setdim(dim)
         use, intrinsic :: iso_c_binding, only: C_PTR
         integer :: dim
         type(C_PTR) :: sat_f_setdim
      end function

      subroutine sat_f_setvec(sat, n, num, buff, id)
         use, intrinsic :: iso_c_binding, only: C_PTR
         type(C_PTR) :: sat
         integer, dimension(:) :: n
         integer :: num, buff
         character(len=*) :: id
      end subroutine

      subroutine sat_f_destroy(sat)
         use, intrinsic :: iso_c_binding, only: C_PTR
         type(C_PTR) :: sat
      end subroutine

      subroutine sat_f_choose(sat, para)
         use, intrinsic :: iso_c_binding, only: C_PTR
         type(C_PTR) :: sat
         integer :: para
      end subroutine

      subroutine sat_f_update(sat, para, val)
         use, intrinsic :: iso_c_binding, only: C_PTR
         type(C_PTR) :: sat
         integer :: para
         double precision :: val
      end subroutine

      function sat_f_terminated(sat)
         use, intrinsic :: iso_c_binding, only: C_PTR
         type(C_PTR) :: sat
         integer :: sat_f_terminated
      end function 
   end interface
end module
