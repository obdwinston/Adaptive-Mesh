program solver
   use config_mod
   use grid_mod
   use flow_mod
   use body_mod
   use io_mod
   implicit none

   real(8) :: nu, dx_c, dy_c, dt, t, max_div, l2_div
   integer :: global_step, lev, lev_idx
   integer :: n_total_patches
   integer(8) :: clock_start, clock_end, clock_rate
   real(8) :: elapsed

   type(Hierarchy) :: hier
   type(Body) :: bd

   ! VTK output tracking
   integer, parameter :: max_vtk_steps = nt / plot_every + 10
   integer :: vtk_steps(max_vtk_steps)
   real(8) :: vtk_times(max_vtk_steps)
   integer :: vtk_npatches(max_vtk_steps)
   integer :: n_vtk_written

   call execute_command_line('mkdir -p ' // trim(output_dir))

   ! === DERIVED QUANTITIES ===

   nu = U_inf * side / Re
   dx_c = Lx / dble(nx_coarse)
   dy_c = Ly / dble(ny_coarse)

   ! === BODY ===

   call read_body(bd)

   ! === HIERARCHY ===

   hier%n_levels = 1 + n_fine_levels
   allocate(hier%levels(hier%n_levels))
   hier%domain = (/ 0.0d0, 0.0d0, Lx, Ly /)
   hier%nu = nu

   ! coarse level
   hier%levels(1)%refinement_ratio = 1
   hier%levels(1)%n_patches = 1
   allocate(hier%levels(1)%patches(1))
   call initialise_patch(hier%levels(1)%patches(1), 0.0d0, 0.0d0, &
      nx_coarse, ny_coarse, dx_c, dy_c, 0, nu)

   ! fine levels (initially empty)
   do lev_idx = 1, n_fine_levels
      hier%levels(lev_idx + 1)%refinement_ratio = fine_ratio(lev_idx)
      hier%levels(lev_idx + 1)%n_patches = 0
      allocate(hier%levels(lev_idx + 1)%patches(0))
   end do

   ! === INITIAL CONDITIONS ===

   call initialise_flow(hier%levels(1)%patches(1))
   call apply_force_patch(hier%levels(1)%patches(1), bd, .false.)
   call apply_bc_velocity_patch(hier%levels(1)%patches(1), hier%domain, .false.)
   call apply_bc_pressure_patch(hier%levels(1)%patches(1), hier%domain)

   call regrid_hierarchy(hier)
   call apply_force_hierarchy(hier, bd, .false.)

   call compute_timestep_hierarchy(hier, dt)

   n_total_patches = 0
   do lev = 1, hier%n_levels
      n_total_patches = n_total_patches + hier%levels(lev)%n_patches
   end do

   write(*, '(A,I0,A,I0,A,F8.4,A,F8.4)') &
      'Grid: ', nx_coarse, ' x ', ny_coarse, ', dx = ', dx_c, ', dy = ', dy_c
   write(*, '(A,F6.1,A,F10.6,A,F10.6)') &
      'Re = ', Re, ', nu = ', nu, ', dt = ', dt
   write(*, '(A,I0,A,F5.1,A,F5.1,A,F5.1)') &
      'Body: ', bd%n_pts, ' pts at (', cx_body, ', ', cy_body, &
      '), side = ', side
   write(*, '(A,I0,A,I0)') &
      'Fine levels: ', n_fine_levels, ', patches: ', n_total_patches

   ! === TIME INTEGRATION ===

   global_step = 0
   t = 0.0d0
   n_vtk_written = 0

   call system_clock(clock_start, clock_rate)

   do while (global_step < nt)

      ! 1. periodic regrid
      if (global_step > 0 .and. mod(global_step, regrid_interval) == 0) &
         call regrid_hierarchy(hier)

      ! 2. compute timestep
      call compute_timestep_hierarchy(hier, dt)

      ! 3. fill ghost layers and apply velocity BCs
      call fill_ghost_hierarchy(hier)
      call apply_bc_velocity_hierarchy(hier, .false.)

      ! 4. predict velocity, apply body forcing and predictor BCs
      call predict_velocity_hierarchy(hier, dt)
      call apply_force_hierarchy(hier, bd, .true.)
      call apply_bc_velocity_hierarchy(hier, .true.)

      ! 5. solve pressure
      call solve_pressure_hierarchy(hier, dt)

      ! 6. correct velocity
      call correct_velocity_hierarchy(hier, dt)

      t = t + dt
      global_step = global_step + 1

      ! diagnostics and VTK output
      if (mod(global_step, plot_every) == 0) then
         call compute_divergence_hierarchy(hier, max_div, l2_div)
         n_total_patches = 0
         do lev = 1, hier%n_levels
            n_total_patches = n_total_patches + hier%levels(lev)%n_patches
         end do

         write(*, '(A,I6,A,F10.4,A,I3,A,ES10.2,A,ES10.2)') &
            'step ', global_step, ', t = ', t, &
            ', patches = ', n_total_patches, &
            ', max_div = ', max_div, ', l2_div = ', l2_div

         call write_output_hierarchy(hier, global_step)

         n_vtk_written = n_vtk_written + 1
         vtk_steps(n_vtk_written) = global_step
         vtk_times(n_vtk_written) = t
         vtk_npatches(n_vtk_written) = n_total_patches
      end if

      ! 7. zero body interior for next step
      call apply_force_hierarchy(hier, bd, .false.)
   end do

   call system_clock(clock_end)
   elapsed = dble(clock_end - clock_start) / dble(clock_rate)
   write(*, '(A,F10.2,A)') 'Elapsed: ', elapsed, ' s'

   ! === PVD OUTPUT ===

   if (n_vtk_written > 0) then
      call write_output_collection(vtk_steps, vtk_times, n_vtk_written, vtk_npatches)
      write(*, '(A,I0,A)') 'Wrote PVD collection with ', n_vtk_written, ' timesteps.'
   end if

   write(*, '(A)') 'Done.'

end program solver
