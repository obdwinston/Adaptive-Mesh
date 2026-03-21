module config_mod
   implicit none

   ! === GRID ===

   integer, parameter :: nx_coarse = 60
   integer, parameter :: ny_coarse = 40
   real(8), parameter :: Lx = 30.0d0
   real(8), parameter :: Ly = 20.0d0

   ! === FLOW ===

   real(8), parameter :: Re = 40.0d0
   real(8), parameter :: U_inf = 1.0d0
   real(8), parameter :: cfl = 0.15d0

   ! === BODY ===

   real(8), parameter :: cx_body = 5.0d0
   real(8), parameter :: cy_body = 10.0d0
   real(8), parameter :: side = 1.0d0
   integer, parameter :: body_buffer = 2
   character(len=64), parameter :: body_file = 'body.dat'

   ! === PRESSURE SOLVER ===

   real(8), parameter :: omega = 1.7d0
   real(8), parameter :: pressure_tol = 1.0d-5
   integer, parameter :: pressure_maxiter = 200000
   integer, parameter :: n_pressure_cycles = 3

   ! === SAMR ===

   integer, parameter :: n_fine_levels = 2
   integer, parameter :: max_fine_levels = 4
   integer, parameter :: fine_ratio(max_fine_levels) = (/ 2, 2, 2, 2 /)
   integer, parameter :: fine_nesting = 0
   integer, parameter :: fine_buffer = 2
   character(len=16), parameter :: flag_method = 'vorticity'  ! 'vorticity' or 'gradient'
   real(8), parameter :: flag_threshold(max_fine_levels) = (/ 0.5d0, 1.0d0, 2.0d0, 4.0d0 /)

   integer, parameter :: regrid_interval = 20
   real(8), parameter :: regrid_x_lim(2) = (/ 0.0d0, Lx - 5.0d0 /)
   real(8), parameter :: regrid_y_lim(2) = (/ 0.0d0, Ly /)

   ! === OUTPUT ===

   integer, parameter :: nt = 20000
   integer, parameter :: plot_every = 40
   character(len=64), parameter :: output_dir = 'output'

   ! === INTERNAL CONSTANTS ===

   real(8), parameter :: geom_tol = 1.0d-10
   integer, parameter :: max_patches = 100

end module config_mod
