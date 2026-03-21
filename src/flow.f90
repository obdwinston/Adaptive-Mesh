module flow_mod
   use config_mod
   use grid_mod
   implicit none
   private

   public :: apply_bc_velocity_patch, apply_bc_velocity_hierarchy
   public :: apply_bc_pressure_patch, apply_bc_pressure_hierarchy
   public :: predict_velocity_patch, predict_velocity_hierarchy
   public :: solve_pressure_patch, solve_pressure_hierarchy
   public :: correct_velocity_patch, correct_velocity_hierarchy
   public :: flag_gradient, flag_vorticity, regrid_level, regrid_hierarchy
   public :: initialise_flow
   public :: compute_divergence_patch, compute_divergence_hierarchy
   public :: compute_timestep_patch, compute_timestep_hierarchy

contains

   ! === BOUNDARY CONDITIONS ===

   subroutine apply_bc_velocity_patch(pa, domain, predictor)
      ! Apply velocity boundary conditions on domain-boundary edges.
      !
      ! Inlet:  Dirichlet u = U_inf, v = 0.
      ! Outlet: Neumann du/dx = 0, dv/dx = 0.
      ! Walls:  slip (du/dy = 0, v = 0).
      !
      ! When predictor is true, operates on u_star/v_star instead of u/v.
      !
      ! Args:
      !     pa: Patch whose boundary velocities are set.
      !     domain: Domain extents (x_lo, y_lo, x_hi, y_hi).
      !     predictor: If true, operate on u_star/v_star instead of u/v.
      !
      ! Vars:
      !     at_l, at_r, at_b, at_t: Flags indicating which domain edges the
      !         patch touches.

      type(Patch), intent(inout), target :: pa
      real(8), intent(in) :: domain(4)
      logical, intent(in) :: predictor

      real(8), pointer :: u(:,:), v(:,:)
      logical :: at_l, at_r, at_b, at_t
      integer :: nx, ny, j, i

      if (predictor) then
         u => pa%u_star; v => pa%v_star
      else
         u => pa%u; v => pa%v
      end if

      nx = pa%nx; ny = pa%ny

      at_l = abs(pa%x0 - domain(1)) < geom_tol
      at_r = abs(patch_x_end(pa) - domain(3)) < geom_tol
      at_b = abs(pa%y0 - domain(2)) < geom_tol
      at_t = abs(patch_y_end(pa) - domain(4)) < geom_tol

      ! inlet
      if (at_l) then
         do j = 1, ny
            u(j, 0) = U_inf
         end do
         do j = 0, ny
            v(j, 0) = -v(j, 1)
         end do
      end if

      ! outlet
      if (at_r) then
         do j = 1, ny
            u(j, nx) = u(j, nx - 1)
         end do
         do j = 0, ny
            v(j, nx + 1) = v(j, nx)
         end do
      end if

      ! bottom wall
      if (at_b) then
         do i = 0, nx
            u(0, i) = u(1, i)
         end do
         do i = 1, nx
            v(0, i) = 0.0d0
         end do
      end if

      ! top wall
      if (at_t) then
         do i = 0, nx
            u(ny + 1, i) = u(ny, i)
         end do
         do i = 1, nx
            v(ny, i) = 0.0d0
         end do
      end if
   end subroutine apply_bc_velocity_patch

   subroutine apply_bc_velocity_hierarchy(hier, predictor)
      ! Apply velocity boundary conditions on all patches in the hierarchy.
      !
      ! Args:
      !     hier: Grid hierarchy whose patches receive boundary conditions.
      !     predictor: If true, operate on u_star/v_star instead of u/v.

      type(Hierarchy), intent(inout) :: hier
      logical, intent(in) :: predictor

      integer :: lev, pidx

      do lev = 1, hier%n_levels
         do pidx = 1, hier%levels(lev)%n_patches
            call apply_bc_velocity_patch(hier%levels(lev)%patches(pidx), &
               hier%domain, predictor)
         end do
      end do
   end subroutine apply_bc_velocity_hierarchy

   subroutine apply_bc_pressure_patch(pa, domain)
      ! Apply pressure boundary conditions on domain-boundary edges.
      !
      ! Inlet:  Neumann dp/dx = 0.
      ! Outlet: Dirichlet p = 0.
      ! Walls:  Neumann dp/dy = 0.
      !
      ! Args:
      !     pa: Patch whose boundary pressures are set.
      !     domain: Domain extents (x_lo, y_lo, x_hi, y_hi).
      !
      ! Vars:
      !     at_l, at_r, at_b, at_t: Flags indicating which domain edges the
      !         patch touches.

      type(Patch), intent(inout) :: pa
      real(8), intent(in) :: domain(4)

      logical :: at_l, at_r, at_b, at_t
      integer :: nx, ny, i, j

      nx = pa%nx; ny = pa%ny

      at_l = abs(pa%x0 - domain(1)) < geom_tol
      at_r = abs(patch_x_end(pa) - domain(3)) < geom_tol
      at_b = abs(pa%y0 - domain(2)) < geom_tol
      at_t = abs(patch_y_end(pa) - domain(4)) < geom_tol

      ! inlet Neumann
      if (at_l) then
         do j = 1, ny
            pa%p(j, 0) = pa%p(j, 1)
         end do
      end if

      ! outlet Dirichlet
      if (at_r) then
         do j = 1, ny
            pa%p(j, nx + 1) = -pa%p(j, nx)
         end do
      end if

      ! bottom Neumann
      if (at_b) then
         do i = 1, nx
            pa%p(0, i) = pa%p(1, i)
         end do
      end if

      ! top Neumann
      if (at_t) then
         do i = 1, nx
            pa%p(ny + 1, i) = pa%p(ny, i)
         end do
      end if
   end subroutine apply_bc_pressure_patch

   subroutine apply_bc_pressure_hierarchy(hier)
      ! Apply pressure boundary conditions on all patches in the hierarchy.
      !
      ! Args:
      !     hier: Grid hierarchy whose patches receive boundary conditions.

      type(Hierarchy), intent(inout) :: hier

      integer :: lev, pidx

      do lev = 1, hier%n_levels
         do pidx = 1, hier%levels(lev)%n_patches
            call apply_bc_pressure_patch(hier%levels(lev)%patches(pidx), &
               hier%domain)
         end do
      end do
   end subroutine apply_bc_pressure_hierarchy

   ! === PROJECTION METHOD ===

   subroutine predict_velocity_patch(pa, dt)
      ! Advance velocity one step using explicit Euler for advection and diffusion.
      !
      ! Args:
      !     pa: Patch to advance.
      !     dt: Timestep.
      !
      ! Vars:
      !     u_e, u_w, u_n, u_s: East/west/north/south face velocities for u-advection.
      !     v_e, v_w, v_n, v_s: East/west/north/south face velocities for v-advection.
      !     adv_x, adv_y: Advection terms in x and y.
      !     diff_x, diff_y: Diffusion terms in x and y.

      type(Patch), intent(inout) :: pa
      real(8), intent(in) :: dt

      integer :: i, j, nx, ny
      real(8) :: dx, dy, nu, dx2, dy2
      real(8) :: u_e, u_w, u_n, u_s
      real(8) :: v_e, v_w, v_n, v_s
      real(8) :: adv_x, adv_y, diff_x, diff_y

      nx = pa%nx; ny = pa%ny
      dx = pa%dx; dy = pa%dy; nu = pa%nu
      dx2 = dx * dx; dy2 = dy * dy

      pa%u_star = pa%u
      pa%v_star = pa%v

      ! u-momentum interior
      do j = 1, ny
         do i = 1, nx - 1
            u_e = 0.5d0 * (pa%u(j, i) + pa%u(j, i + 1))
            u_w = 0.5d0 * (pa%u(j, i - 1) + pa%u(j, i))
            adv_x = (u_e * u_e - u_w * u_w) / dx

            u_n = 0.5d0 * (pa%u(j, i) + pa%u(j + 1, i))
            v_n = 0.5d0 * (pa%v(j, i) + pa%v(j, i + 1))
            u_s = 0.5d0 * (pa%u(j - 1, i) + pa%u(j, i))
            v_s = 0.5d0 * (pa%v(j - 1, i) + pa%v(j - 1, i + 1))
            adv_y = (u_n * v_n - u_s * v_s) / dy

            diff_x = (pa%u(j, i + 1) - 2.0d0 * pa%u(j, i) + pa%u(j, i - 1)) / dx2
            diff_y = (pa%u(j + 1, i) - 2.0d0 * pa%u(j, i) + pa%u(j - 1, i)) / dy2

            pa%u_star(j, i) = pa%u(j, i) + dt * (-adv_x - adv_y + nu * (diff_x + diff_y))
         end do
      end do

      ! v-momentum interior
      do j = 1, ny - 1
         do i = 1, nx
            u_e = 0.5d0 * (pa%u(j, i) + pa%u(j + 1, i))
            v_e = 0.5d0 * (pa%v(j, i) + pa%v(j, i + 1))
            u_w = 0.5d0 * (pa%u(j, i - 1) + pa%u(j + 1, i - 1))
            v_w = 0.5d0 * (pa%v(j, i - 1) + pa%v(j, i))
            adv_x = (u_e * v_e - u_w * v_w) / dx

            v_n = 0.5d0 * (pa%v(j, i) + pa%v(j + 1, i))
            v_s = 0.5d0 * (pa%v(j - 1, i) + pa%v(j, i))
            adv_y = (v_n * v_n - v_s * v_s) / dy

            diff_x = (pa%v(j, i + 1) - 2.0d0 * pa%v(j, i) + pa%v(j, i - 1)) / dx2
            diff_y = (pa%v(j + 1, i) - 2.0d0 * pa%v(j, i) + pa%v(j - 1, i)) / dy2

            pa%v_star(j, i) = pa%v(j, i) + dt * (-adv_x - adv_y + nu * (diff_x + diff_y))
         end do
      end do
   end subroutine predict_velocity_patch

   subroutine predict_velocity_hierarchy(hier, dt)
      ! Advance velocity one step on all patches.
      !
      ! Args:
      !     hier: Grid hierarchy whose patches are advanced.
      !     dt: Timestep.

      type(Hierarchy), intent(inout) :: hier
      real(8), intent(in) :: dt

      integer :: lev, pidx

      do lev = 1, hier%n_levels
         do pidx = 1, hier%levels(lev)%n_patches
            call predict_velocity_patch(hier%levels(lev)%patches(pidx), dt)
         end do
      end do
   end subroutine predict_velocity_hierarchy

   subroutine solve_pressure_patch(pa, domain, dt)
      ! Solve the pressure Poisson equation on one patch using SOR.
      !
      ! Args:
      !     pa: Patch on which the pressure Poisson equation is solved.
      !     domain: Domain extents (x_lo, y_lo, x_hi, y_hi).
      !     dt: Timestep used to form the divergence RHS.
      !
      ! Vars:
      !     coeff: Precomputed Laplacian denominator for SOR update.
      !     p_nb: Sum of neighbour pressure contributions.
      !     p_new: Updated pressure after one SOR relaxation.

      type(Patch), intent(inout) :: pa
      real(8), intent(in) :: domain(4), dt

      integer :: i, j, k, nx, ny
      real(8) :: dx, dy, dx2, dy2, coeff, err
      real(8) :: p_nb, p_new, p_diff

      nx = pa%nx; ny = pa%ny
      dx = pa%dx; dy = pa%dy
      dx2 = dx * dx; dy2 = dy * dy
      coeff = 1.0d0 / (2.0d0 / dx2 + 2.0d0 / dy2)

      ! compute RHS: divergence of predictor velocity
      do j = 1, ny
         do i = 1, nx
            pa%rhs(j, i) = (1.0d0 / dt) * ( &
               (pa%u_star(j, i) - pa%u_star(j, i - 1)) / dx + &
               (pa%v_star(j, i) - pa%v_star(j - 1, i)) / dy)
         end do
      end do

      do k = 1, pressure_maxiter
         err = 0.0d0

         do j = 1, ny
            do i = 1, nx
               p_nb = (pa%p(j, i + 1) + pa%p(j, i - 1)) / dx2 + &
                  (pa%p(j + 1, i) + pa%p(j - 1, i)) / dy2
               p_new = (1.0d0 - omega) * pa%p(j, i) + &
                  omega * coeff * (p_nb - pa%rhs(j, i))
               p_diff = abs(p_new - pa%p(j, i))
               if (p_diff > err) err = p_diff
               pa%p(j, i) = p_new
            end do
         end do
         call apply_bc_pressure_patch(pa, domain)

         if (err < pressure_tol) exit
      end do
   end subroutine solve_pressure_patch

   subroutine solve_pressure_hierarchy(hier, dt)
      ! Solve pressure across all levels using composite V-cycles.
      !
      ! Args:
      !     hier: Grid hierarchy containing all levels and patches.
      !     dt: Timestep used to form the divergence RHS.

      type(Hierarchy), intent(inout) :: hier
      real(8), intent(in) :: dt

      integer :: cyc, lev, fp_idx, pp_idx

      do cyc = 1, n_pressure_cycles

         ! solve on coarse level
         do pp_idx = 1, hier%levels(1)%n_patches
            call solve_pressure_patch(hier%levels(1)%patches(pp_idx), &
               hier%domain, dt)
         end do

         ! propagate to each fine level
         do lev = 2, hier%n_levels

            ! update fine ghost pressure from coarse
            do fp_idx = 1, hier%levels(lev)%n_patches
               do pp_idx = 1, hier%levels(lev - 1)%n_patches
                  call fill_ghost_pressure(hier%levels(lev)%patches(fp_idx), &
                     hier%levels(lev - 1)%patches(pp_idx))
               end do
               call apply_bc_pressure_patch(hier%levels(lev)%patches(fp_idx), &
                  hier%domain)
            end do

            ! solve on fine level
            do fp_idx = 1, hier%levels(lev)%n_patches
               call solve_pressure_patch(hier%levels(lev)%patches(fp_idx), &
                  hier%domain, dt)
            end do

            ! restrict fine pressure back to coarse
            do fp_idx = 1, hier%levels(lev)%n_patches
               do pp_idx = 1, hier%levels(lev - 1)%n_patches
                  call restrict_pressure(hier%levels(lev)%patches(fp_idx), &
                     hier%levels(lev - 1)%patches(pp_idx))
               end do
            end do
            do pp_idx = 1, hier%levels(lev - 1)%n_patches
               call apply_bc_pressure_patch(hier%levels(lev - 1)%patches(pp_idx), &
                  hier%domain)
            end do
         end do
      end do
   end subroutine solve_pressure_hierarchy

   subroutine correct_velocity_patch(pa, dt)
      ! Correct velocity on one patch: u = u_star - dt * grad(p).
      !
      ! Args:
      !     pa: Patch whose velocity is corrected in place.
      !     dt: Timestep.

      type(Patch), intent(inout) :: pa
      real(8), intent(in) :: dt

      integer :: i, j, nx, ny
      real(8) :: dx, dy

      nx = pa%nx; ny = pa%ny; dx = pa%dx; dy = pa%dy

      do j = 1, ny
         do i = 1, nx - 1
            pa%u(j, i) = pa%u_star(j, i) - dt * (pa%p(j, i + 1) - pa%p(j, i)) / dx
         end do
      end do

      do j = 1, ny - 1
         do i = 1, nx
            pa%v(j, i) = pa%v_star(j, i) - dt * (pa%p(j + 1, i) - pa%p(j, i)) / dy
         end do
      end do
   end subroutine correct_velocity_patch

   subroutine correct_velocity_hierarchy(hier, dt)
      ! Apply composite velocity correction across all levels, finest first,
      ! then restrict corrections onto coarser levels.
      !
      ! Args:
      !     hier: Grid hierarchy containing all levels and patches.
      !     dt: Timestep.

      type(Hierarchy), intent(inout) :: hier
      real(8), intent(in) :: dt

      integer :: lev, fp_idx, pp_idx

      do lev = hier%n_levels, 1, -1
         do fp_idx = 1, hier%levels(lev)%n_patches
            call correct_velocity_patch(hier%levels(lev)%patches(fp_idx), dt)
         end do
      end do

      do lev = hier%n_levels, 2, -1
         do fp_idx = 1, hier%levels(lev)%n_patches
            do pp_idx = 1, hier%levels(lev - 1)%n_patches
               call restrict_velocity(hier%levels(lev)%patches(fp_idx), &
                  hier%levels(lev - 1)%patches(pp_idx))
            end do
         end do
      end do
   end subroutine correct_velocity_hierarchy

   ! === GRID ADAPTATION ===

   subroutine flag_gradient(pa, threshold, nesting, flagged)
      ! Flag cells where the speed gradient magnitude exceeds a threshold.
      !
      ! Args:
      !     pa: Patch to evaluate.
      !     threshold: Gradient magnitude above which a cell is flagged.
      !     nesting: Number of cells to inset from the patch boundary
      !         (avoids flagging near coarse-fine interfaces).
      !     flagged: Output boolean mask (pa%ny x pa%nx).
      !
      ! Vars:
      !     speed: Cell-centred speed magnitude.
      !     S: Speed padded with one ghost layer for central differencing.
      !     dSdx, dSdy: Central-difference speed gradients.

      type(Patch), intent(in) :: pa
      real(8), intent(in) :: threshold
      integer, intent(in) :: nesting
      logical, intent(out) :: flagged(pa%ny, pa%nx)

      integer :: i, j, nx, ny, n
      real(8) :: u_cc, v_cc, dSdx, dSdy, grad_mag, xc, yc
      real(8), allocatable :: speed(:,:), S(:,:)

      nx = pa%nx; ny = pa%ny; n = nesting
      allocate(speed(ny, nx), S(0:ny+1, 0:nx+1))

      ! compute cell-centred speed
      do j = 1, ny
         do i = 1, nx
            u_cc = 0.5d0 * (pa%u(j, i - 1) + pa%u(j, i))
            v_cc = 0.5d0 * (pa%v(j - 1, i) + pa%v(j, i))
            speed(j, i) = sqrt(u_cc * u_cc + v_cc * v_cc)
         end do
      end do

      ! embed in padded array for central differences
      S = 0.0d0
      S(1:ny, 1:nx) = speed
      S(0, 1:nx) = speed(1, :)
      S(ny+1, 1:nx) = speed(ny, :)
      do j = 0, ny + 1
         S(j, 0) = S(j, 1)
         S(j, nx+1) = S(j, nx)
      end do

      flagged = .false.

      do j = 1 + n, ny - n
         yc = pa%y0 + (j - 0.5d0) * pa%dy
         if (yc < regrid_y_lim(1) .or. yc > regrid_y_lim(2)) cycle
         do i = 1 + n, nx - n
            xc = pa%x0 + (i - 0.5d0) * pa%dx
            if (xc < regrid_x_lim(1) .or. xc > regrid_x_lim(2)) cycle
            dSdx = (S(j, i + 1) - S(j, i - 1)) / (2.0d0 * pa%dx)
            dSdy = (S(j + 1, i) - S(j - 1, i)) / (2.0d0 * pa%dy)
            grad_mag = sqrt(dSdx * dSdx + dSdy * dSdy)
            if (grad_mag > threshold) flagged(j, i) = .true.
         end do
      end do

      deallocate(speed, S)
   end subroutine flag_gradient


   subroutine flag_vorticity(pa, threshold, nesting, flagged)
      ! Flag cells where the vorticity magnitude exceeds a threshold.
      !
      ! Args:
      !     pa: Patch to evaluate.
      !     threshold: Vorticity magnitude above which a cell is flagged.
      !     nesting: Number of cells to inset from the patch boundary
      !         (avoids flagging near coarse-fine interfaces).
      !     flagged: Output boolean mask (pa%ny x pa%nx).

      type(Patch), intent(in) :: pa
      real(8), intent(in) :: threshold
      integer, intent(in) :: nesting
      logical, intent(out) :: flagged(pa%ny, pa%nx)

      integer :: i, j, nx, ny, n
      real(8) :: dvdx, dudy, omega_z, xc, yc

      nx = pa%nx; ny = pa%ny; n = nesting

      flagged = .false.

      do j = 1 + n, ny - n
         yc = pa%y0 + (j - 0.5d0) * pa%dy
         if (yc < regrid_y_lim(1) .or. yc > regrid_y_lim(2)) cycle
         do i = 1 + n, nx - n
            xc = pa%x0 + (i - 0.5d0) * pa%dx
            if (xc < regrid_x_lim(1) .or. xc > regrid_x_lim(2)) cycle

            ! dv/dx at cell centre: difference of v on vertical faces
            dvdx = (pa%v(j, i) - pa%v(j, i - 1)) / pa%dx
            ! du/dy at cell centre: difference of u on horizontal faces
            dudy = (pa%u(j, i) - pa%u(j - 1, i)) / pa%dy

            omega_z = dvdx - dudy
            if (abs(omega_z) > threshold) flagged(j, i) = .true.
         end do
      end do

   end subroutine flag_vorticity


   subroutine regrid_level(hier, target_level, threshold, buffer, nesting)
      ! Rebuild patches at the target level based on velocity gradients and
      ! body coverage on the parent level, preserving data where possible.
      !
      ! Args:
      !     hier: Grid hierarchy to modify.
      !     target_level: Level index to rebuild (must be >= 2).
      !     threshold: Gradient magnitude threshold for flagging cells.
      !     buffer: Number of coarse cells to pad around each cluster.
      !     nesting: Nesting width passed to flag_gradient.
      !
      ! Vars:
      !     old_fine: Snapshot of existing patches before regridding.
      !     new_fine: Newly created patches.
      !     boxes: Bounding-box descriptors from cluster_flagged.
      !     bx_lo, bx_hi, by_lo, by_hi: Body bounding box with buffer.

      type(Hierarchy), intent(inout) :: hier
      integer, intent(in) :: target_level, buffer, nesting
      real(8), intent(in) :: threshold

      integer :: r, n_old, n_new, n_total
      integer :: par_idx, old_idx, box_idx, n_boxes, k, i, j
      real(8) :: boxes(6, max_patches)
      type(Patch), allocatable :: old_fine(:), new_fine(:)
      type(Patch) :: new_patch
      logical, allocatable :: flagged(:,:)
      real(8) :: xc, yc, bx_lo, bx_hi, by_lo, by_hi

      r = hier%levels(target_level)%refinement_ratio

      ! snapshot existing patches
      n_old = hier%levels(target_level)%n_patches
      if (n_old > 0) then
         allocate(old_fine(n_old))
         do k = 1, n_old
            old_fine(k) = hier%levels(target_level)%patches(k)
         end do
      end if

      ! accumulate new patches
      allocate(new_fine(max_patches))
      n_new = 0

      do par_idx = 1, hier%levels(target_level - 1)%n_patches
         associate(parent => hier%levels(target_level - 1)%patches(par_idx))
            allocate(flagged(parent%ny, parent%nx))
            if (trim(flag_method) == 'vorticity') then
               call flag_vorticity(parent, threshold, nesting, flagged)
            else
               call flag_gradient(parent, threshold, nesting, flagged)
            end if

            ! flag body region
            bx_lo = cx_body - 0.5d0 * side - body_buffer * parent%dx
            bx_hi = cx_body + 0.5d0 * side + body_buffer * parent%dx
            by_lo = cy_body - 0.5d0 * side - body_buffer * parent%dy
            by_hi = cy_body + 0.5d0 * side + body_buffer * parent%dy
            do j = 1, parent%ny
               yc = parent%y0 + (j - 0.5d0) * parent%dy
               do i = 1, parent%nx
                  xc = parent%x0 + (i - 0.5d0) * parent%dx
                  if (xc >= bx_lo .and. xc <= bx_hi .and. &
                     yc >= by_lo .and. yc <= by_hi) then
                     flagged(j, i) = .true.
                  end if
               end do
            end do

            call cluster_flagged(parent, flagged, parent%ny, parent%nx, r, buffer, &
               boxes, n_boxes)
            deallocate(flagged)

            do box_idx = 1, n_boxes
               call initialise_patch(new_patch, &
                  boxes(1, box_idx), boxes(2, box_idx), &
                  nint(boxes(3, box_idx)), nint(boxes(4, box_idx)), &
                  boxes(5, box_idx), boxes(6, box_idx), &
                  target_level - 1, hier%nu)
               call fill_ghost_patch(new_patch, parent)
               call prolong_patch(parent, new_patch)

               do old_idx = 1, n_old
                  call copy_patch(old_fine(old_idx), new_patch)
               end do

               n_new = n_new + 1
               new_fine(n_new) = new_patch
            end do
         end associate
      end do

      ! reassemble level
      n_total = n_new
      if (allocated(hier%levels(target_level)%patches)) &
         deallocate(hier%levels(target_level)%patches)
      hier%levels(target_level)%n_patches = n_total
      allocate(hier%levels(target_level)%patches(max(n_total, 0)))
      do k = 1, n_new
         hier%levels(target_level)%patches(k) = new_fine(k)
      end do

      ! apply BCs on new patches
      do box_idx = 1, n_new
         call apply_bc_velocity_patch( &
            hier%levels(target_level)%patches(box_idx), &
            hier%domain, .false.)
         call apply_bc_pressure_patch( &
            hier%levels(target_level)%patches(box_idx), &
            hier%domain)
      end do

      if (allocated(old_fine)) deallocate(old_fine)
      deallocate(new_fine)
   end subroutine regrid_level

   subroutine regrid_hierarchy(hier)
      ! Rebuild adaptive patches at all fine levels.
      !
      ! Args:
      !     hier: Grid hierarchy to regrid.

      type(Hierarchy), intent(inout) :: hier

      integer :: lev_idx

      do lev_idx = 1, n_fine_levels
         call regrid_level(hier, lev_idx + 1, flag_threshold(lev_idx), &
            fine_buffer, fine_nesting)
      end do
   end subroutine regrid_hierarchy

   ! === SOLVER UTILITIES ===

   subroutine initialise_flow(pa)
      ! Set uniform flow and zero pressure as initial conditions.
      !
      ! Args:
      !     pa: Patch to initialise.

      type(Patch), intent(inout) :: pa

      pa%p = 0.0d0
      pa%v = 0.0d0
      pa%u = U_inf
   end subroutine initialise_flow

   subroutine compute_divergence_patch(pa, max_div, sum_sq, n_cells)
      ! Compute divergence statistics on a single patch.
      !
      ! Args:
      !     pa: Patch to evaluate.
      !     max_div: Maximum absolute divergence found (output).
      !     sum_sq: Sum of squared divergence values (output).
      !     n_cells: Number of cells evaluated (output).

      type(Patch), intent(in) :: pa
      real(8), intent(out) :: max_div, sum_sq
      integer, intent(out) :: n_cells

      integer :: i, j
      real(8) :: div_val

      max_div = 0.0d0
      sum_sq = 0.0d0
      n_cells = pa%nx * pa%ny
      do j = 1, pa%ny
         do i = 1, pa%nx
            div_val = (pa%u(j, i) - pa%u(j, i - 1)) / pa%dx + &
               (pa%v(j, i) - pa%v(j - 1, i)) / pa%dy
            if (abs(div_val) > max_div) max_div = abs(div_val)
            sum_sq = sum_sq + div_val * div_val
         end do
      end do
   end subroutine compute_divergence_patch

   subroutine compute_divergence_hierarchy(hier, max_div, l2_div)
      ! Compute divergence statistics across all patches.
      !
      ! Args:
      !     hier: Grid hierarchy to evaluate.
      !     max_div: Maximum absolute divergence found (output).
      !     l2_div: RMS (L2-norm) divergence across all cells (output).

      type(Hierarchy), intent(in) :: hier
      real(8), intent(out) :: max_div, l2_div

      integer :: lev, pidx, n_cells, total_cells
      real(8) :: patch_div, patch_sq, total_sq

      max_div = 0.0d0
      total_sq = 0.0d0
      total_cells = 0
      do lev = 1, hier%n_levels
         do pidx = 1, hier%levels(lev)%n_patches
            call compute_divergence_patch(hier%levels(lev)%patches(pidx), &
               patch_div, patch_sq, n_cells)
            if (patch_div > max_div) max_div = patch_div
            total_sq = total_sq + patch_sq
            total_cells = total_cells + n_cells
         end do
      end do
      l2_div = sqrt(total_sq / dble(max(total_cells, 1)))
   end subroutine compute_divergence_hierarchy

   subroutine compute_timestep_patch(pa, dt_patch)
      ! Compute the minimum stable timestep for a single patch.
      !
      ! Args:
      !     pa: Patch to evaluate.
      !     dt_patch: Minimum stable timestep for this patch (output).
      !
      ! Vars:
      !     u_max, v_max: Maximum absolute velocity components on the patch.
      !     dt_diff: Diffusive stability limit.
      !     dt_adv: Advective CFL limit.

      type(Patch), intent(in) :: pa
      real(8), intent(out) :: dt_patch

      integer :: i, j
      real(8) :: u_max, v_max, dt_diff, dt_adv
      real(8) :: abs_val

      u_max = 1.0d-10
      do i = lbound(pa%u, 2), ubound(pa%u, 2)
         do j = lbound(pa%u, 1), ubound(pa%u, 1)
            abs_val = abs(pa%u(j, i))
            if (abs_val > u_max) u_max = abs_val
         end do
      end do

      v_max = 1.0d-10
      do i = lbound(pa%v, 2), ubound(pa%v, 2)
         do j = lbound(pa%v, 1), ubound(pa%v, 1)
            abs_val = abs(pa%v(j, i))
            if (abs_val > v_max) v_max = abs_val
         end do
      end do

      ! diffusive stability limit
      if (pa%nu > 0.0d0) then
         dt_diff = 0.25d0 / (pa%nu * (1.0d0 / (pa%dx * pa%dx) + &
            1.0d0 / (pa%dy * pa%dy)))
      else
         dt_diff = 1.0d30
      end if

      ! advective CFL limit
      dt_adv = min(pa%dx / u_max, pa%dy / v_max)

      dt_patch = cfl * min(dt_diff, dt_adv)
   end subroutine compute_timestep_patch

   subroutine compute_timestep_hierarchy(hier, dt)
      ! Compute the minimum stable timestep across all patches.
      !
      ! Args:
      !     hier: Grid hierarchy to evaluate.
      !     dt: Minimum stable timestep (output).
      !
      ! Vars:
      !     u_max, v_max: Maximum absolute velocity components on a patch.
      !     dt_diff: Diffusive stability limit.
      !     dt_adv: Advective CFL limit.
      !     dt_patch: Per-patch minimum of dt_diff and dt_adv scaled by CFL.

      type(Hierarchy), intent(in) :: hier
      real(8), intent(out) :: dt

      integer :: lev, pidx
      real(8) :: dt_patch

      dt = 1.0d30

      do lev = 1, hier%n_levels
         do pidx = 1, hier%levels(lev)%n_patches
            call compute_timestep_patch(hier%levels(lev)%patches(pidx), dt_patch)
            if (dt_patch < dt) dt = dt_patch
         end do
      end do
   end subroutine compute_timestep_hierarchy

end module flow_mod
