program main
   use sat_f
   use iso_c_binding
   implicit none
   type(C_PTR) :: sat
   integer, dimension(:) :: n(3)

   n(1) = 10
   n(2) = 10
   n(3) = 10
   sat = sat_f_setdim(3)
   call sat_f_setvec(sat, n, 16, 100, "./lhd_16_1.csv" // char(0))

   call sat_f_destroy(sat)
end program
