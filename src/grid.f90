module grid_mod
   use config_mod
   implicit none
   private

   ! === DATA STRUCTURES ===

   type, public :: Patch
      real(8) :: x0, y0        ! lower-left corner
      integer :: nx, ny        ! number of pressure cells
      real(8) :: dx, dy        ! cell spacing
      integer :: level         ! refinement level (0 = coarsest)
      real(8) :: nu            ! kinematic viscosity

      ! Solution arrays allocated with 0-based lower bound.
      ! u(0:ny+1, 0:nx),  v(0:ny, 0:nx+1),  p(0:ny+1, 0:nx+1)
      ! u_star same shape as u, v_star same shape as v, rhs(1:ny, 1:nx)
      real(8), allocatable :: u(:,:), v(:,:), p(:,:)
      real(8), allocatable :: u_star(:,:), v_star(:,:), rhs(:,:)
   end type Patch

   type, public :: Level
      integer :: refinement_ratio
      integer :: n_patches
      type(Patch), allocatable :: patches(:)
   end type Level

   type, public :: Hierarchy
      integer :: n_levels
      type(Level), allocatable :: levels(:)
      real(8) :: domain(4)     ! (x_lo, y_lo, x_hi, y_hi)
      real(8) :: nu
   end type Hierarchy

   public :: initialise_patch
   public :: patch_x_end, patch_y_end
   public :: interpolate_value
   public :: fill_ghost_field
   public :: fill_ghost_pressure, fill_ghost_velocity
   public :: fill_ghost_patch, fill_ghost_hierarchy
   public :: restrict_pressure, restrict_velocity
   public :: prolong_pressure, prolong_velocity, prolong_patch
   public :: copy_pressure, copy_velocity, copy_patch
   public :: connect_components, cluster_flagged

contains

   ! === PATCH LIFECYCLE ===

   subroutine initialise_patch(pa, x0, y0, nx, ny, dx, dy, lev, nu)
      ! Allocate all arrays on a patch with 0-based indexing.
      !
      ! Args:
      !     pa: Patch to initialise.
      !     x0: X-coordinate of the lower-left corner.
      !     y0: Y-coordinate of the lower-left corner.
      !     nx: Number of pressure cells in x.
      !     ny: Number of pressure cells in y.
      !     dx: Cell spacing in x.
      !     dy: Cell spacing in y.
      !     lev: Refinement level (0 = coarsest).
      !     nu: Kinematic viscosity.

      type(Patch), intent(inout) :: pa
      real(8), intent(in) :: x0, y0, dx, dy, nu
      integer, intent(in) :: nx, ny, lev

      pa%x0 = x0; pa%y0 = y0
      pa%nx = nx; pa%ny = ny
      pa%dx = dx; pa%dy = dy
      pa%level = lev; pa%nu = nu

      if (allocated(pa%u)) deallocate(pa%u)
      if (allocated(pa%v)) deallocate(pa%v)
      if (allocated(pa%p)) deallocate(pa%p)
      if (allocated(pa%u_star)) deallocate(pa%u_star)
      if (allocated(pa%v_star)) deallocate(pa%v_star)
      if (allocated(pa%rhs)) deallocate(pa%rhs)

      allocate(pa%u(0:ny+1, 0:nx)); pa%u = 0.0d0
      allocate(pa%v(0:ny, 0:nx+1)); pa%v = 0.0d0
      allocate(pa%p(0:ny+1, 0:nx+1)); pa%p = 0.0d0
      allocate(pa%u_star(0:ny+1, 0:nx)); pa%u_star = 0.0d0
      allocate(pa%v_star(0:ny, 0:nx+1)); pa%v_star = 0.0d0
      allocate(pa%rhs(1:ny, 1:nx)); pa%rhs = 0.0d0
   end subroutine initialise_patch

   real(8) function patch_x_end(pa)
      ! Return the right edge coordinate of a patch.
      !
      ! Args:
      !     pa: Patch whose right edge is queried.

      type(Patch), intent(in) :: pa
      patch_x_end = pa%x0 + pa%nx * pa%dx
   end function patch_x_end

   real(8) function patch_y_end(pa)
      ! Return the top edge coordinate of a patch.
      !
      ! Args:
      !     pa: Patch whose top edge is queried.

      type(Patch), intent(in) :: pa
      patch_y_end = pa%y0 + pa%ny * pa%dy
   end function patch_y_end

   ! === FIELD INTERPOLATION ===

   real(8) function interpolate_value(data, nj, ni, ox, oy, dx, dy, x, y)
      ! Bilinear interpolation of data(0:nj-1, 0:ni-1) at position (x, y).
      !
      ! Args:
      !     data: 2-D field array with 0-based indexing.
      !     nj: Number of rows in data.
      !     ni: Number of columns in data.
      !     ox: X-origin of the data grid.
      !     oy: Y-origin of the data grid.
      !     dx: Grid spacing in x.
      !     dy: Grid spacing in y.
      !     x: X-coordinate of the interpolation point.
      !     y: Y-coordinate of the interpolation point.
      !
      ! Vars:
      !     fi, fj: Continuous index positions in x and y.
      !     tx, ty: Interpolation weights clamped to [0, 1].

      integer, intent(in) :: nj, ni
      real(8), intent(in) :: data(0:nj-1, 0:ni-1)
      real(8), intent(in) :: ox, oy, dx, dy, x, y

      real(8) :: fi, fj, tx, ty
      integer :: i0, j0

      fi = (x - ox) / dx
      fj = (y - oy) / dy
      i0 = max(0, min(int(floor(fi)), ni - 2))
      j0 = max(0, min(int(floor(fj)), nj - 2))
      tx = max(0.0d0, min(1.0d0, fi - dble(i0)))
      ty = max(0.0d0, min(1.0d0, fj - dble(j0)))

      interpolate_value = (1.0d0 - tx) * (1.0d0 - ty) * data(j0, i0) &
         + tx * (1.0d0 - ty) * data(j0, i0+1) &
         + (1.0d0 - tx) * ty * data(j0+1, i0) &
         + tx * ty * data(j0+1, i0+1)
   end function interpolate_value

   ! === LEVEL TRANSFERS ===

   subroutine fill_ghost_field(fine_data, nj_f, ni_f, f_ox, f_oy, f_dx, f_dy, &
      coarse_data, nj_c, ni_c, c_ox, c_oy, c_dx, c_dy)
      ! Fill the outermost ghost ring of fine_data by interpolating from coarse_data.
      !
      ! Args:
      !     fine_data: Fine-grid field array whose ghost ring is written.
      !     nj_f: Number of rows in fine_data.
      !     ni_f: Number of columns in fine_data.
      !     f_ox: X-origin of the fine data grid.
      !     f_oy: Y-origin of the fine data grid.
      !     f_dx: Fine grid spacing in x.
      !     f_dy: Fine grid spacing in y.
      !     coarse_data: Coarse-grid field array used as interpolation source.
      !     nj_c: Number of rows in coarse_data.
      !     ni_c: Number of columns in coarse_data.
      !     c_ox: X-origin of the coarse data grid.
      !     c_oy: Y-origin of the coarse data grid.
      !     c_dx: Coarse grid spacing in x.
      !     c_dy: Coarse grid spacing in y.

      integer, intent(in) :: nj_f, ni_f, nj_c, ni_c
      real(8), intent(inout) :: fine_data(0:nj_f-1, 0:ni_f-1)
      real(8), intent(in) :: coarse_data(0:nj_c-1, 0:ni_c-1)
      real(8), intent(in) :: f_ox, f_oy, f_dx, f_dy
      real(8), intent(in) :: c_ox, c_oy, c_dx, c_dy

      integer :: i, j
      real(8) :: xp, yp, y_top, x_right

      ! skip if fine ghost positions fall outside coarse extent
      if (f_ox < c_ox - geom_tol) return
      if (f_ox + (ni_f - 1) * f_dx > c_ox + (ni_c - 1) * c_dx + geom_tol) return
      if (f_oy < c_oy - geom_tol) return
      if (f_oy + (nj_f - 1) * f_dy > c_oy + (nj_c - 1) * c_dy + geom_tol) return

      ! bottom ghost row
      do i = 0, ni_f - 1
         xp = f_ox + i * f_dx
         fine_data(0, i) = interpolate_value(coarse_data, nj_c, ni_c, &
            c_ox, c_oy, c_dx, c_dy, xp, f_oy)
      end do

      ! top ghost row
      y_top = f_oy + (nj_f - 1) * f_dy
      do i = 0, ni_f - 1
         xp = f_ox + i * f_dx
         fine_data(nj_f - 1, i) = interpolate_value(coarse_data, nj_c, ni_c, &
            c_ox, c_oy, c_dx, c_dy, xp, y_top)
      end do

      ! left ghost column (interior rows only)
      do j = 1, nj_f - 2
         yp = f_oy + j * f_dy
         fine_data(j, 0) = interpolate_value(coarse_data, nj_c, ni_c, &
            c_ox, c_oy, c_dx, c_dy, f_ox, yp)
      end do

      ! right ghost column (interior rows only)
      x_right = f_ox + (ni_f - 1) * f_dx
      do j = 1, nj_f - 2
         yp = f_oy + j * f_dy
         fine_data(j, ni_f - 1) = interpolate_value(coarse_data, nj_c, ni_c, &
            c_ox, c_oy, c_dx, c_dy, x_right, yp)
      end do
   end subroutine fill_ghost_field

   subroutine fill_ghost_pressure(fine, coarse)
      ! Fill pressure ghost layer on a fine patch from a coarse patch.
      !
      ! Args:
      !     fine: Fine patch whose pressure ghosts are written.
      !     coarse: Coarse patch used as interpolation source.

      type(Patch), intent(inout) :: fine
      type(Patch), intent(in) :: coarse

      real(8) :: f_ox, f_oy, c_ox, c_oy
      integer :: nj_f, ni_f, nj_c, ni_c

      f_ox = fine%x0 - 0.5d0 * fine%dx
      f_oy = fine%y0 - 0.5d0 * fine%dy
      c_ox = coarse%x0 - 0.5d0 * coarse%dx
      c_oy = coarse%y0 - 0.5d0 * coarse%dy
      nj_f = fine%ny + 2; ni_f = fine%nx + 2
      nj_c = coarse%ny + 2; ni_c = coarse%nx + 2
      call fill_ghost_field(fine%p, nj_f, ni_f, f_ox, f_oy, fine%dx, fine%dy, &
         coarse%p, nj_c, ni_c, c_ox, c_oy, coarse%dx, coarse%dy)
   end subroutine fill_ghost_pressure

   subroutine fill_ghost_velocity(fine, coarse)
      ! Fill velocity ghost layers on a fine patch from a coarse patch.
      !
      ! Args:
      !     fine: Fine patch whose velocity ghosts are written.
      !     coarse: Coarse patch used as interpolation source.

      type(Patch), intent(inout) :: fine
      type(Patch), intent(in) :: coarse

      real(8) :: f_ox, f_oy, c_ox, c_oy
      integer :: nj_f, ni_f, nj_c, ni_c

      ! u ghost layer
      f_ox = fine%x0; f_oy = fine%y0 - 0.5d0 * fine%dy
      c_ox = coarse%x0; c_oy = coarse%y0 - 0.5d0 * coarse%dy
      nj_f = fine%ny + 2; ni_f = fine%nx + 1
      nj_c = coarse%ny + 2; ni_c = coarse%nx + 1
      call fill_ghost_field(fine%u, nj_f, ni_f, f_ox, f_oy, fine%dx, fine%dy, &
         coarse%u, nj_c, ni_c, c_ox, c_oy, coarse%dx, coarse%dy)

      ! v ghost layer
      f_ox = fine%x0 - 0.5d0 * fine%dx; f_oy = fine%y0
      c_ox = coarse%x0 - 0.5d0 * coarse%dx; c_oy = coarse%y0
      nj_f = fine%ny + 1; ni_f = fine%nx + 2
      nj_c = coarse%ny + 1; ni_c = coarse%nx + 2
      call fill_ghost_field(fine%v, nj_f, ni_f, f_ox, f_oy, fine%dx, fine%dy, &
         coarse%v, nj_c, ni_c, c_ox, c_oy, coarse%dx, coarse%dy)
   end subroutine fill_ghost_velocity

   subroutine fill_ghost_patch(fine, coarse)
      ! Fill all ghost layers on a fine patch from a coarse patch.
      !
      ! Args:
      !     fine: Fine patch whose ghost layers are written.
      !     coarse: Coarse patch used as interpolation source.

      type(Patch), intent(inout) :: fine
      type(Patch), intent(in) :: coarse

      call fill_ghost_pressure(fine, coarse)
      call fill_ghost_velocity(fine, coarse)
   end subroutine fill_ghost_patch

   subroutine fill_ghost_hierarchy(hier)
      ! Fill ghost layers on all fine patches from their parent level.
      !
      ! Args:
      !     hier: Grid hierarchy whose fine-level ghosts are filled.

      type(Hierarchy), intent(inout) :: hier

      integer :: lev, fp_idx, pp_idx

      do lev = 2, hier%n_levels
         do fp_idx = 1, hier%levels(lev)%n_patches
            do pp_idx = 1, hier%levels(lev - 1)%n_patches
               call fill_ghost_patch(hier%levels(lev)%patches(fp_idx), &
                  hier%levels(lev - 1)%patches(pp_idx))
            end do
         end do
      end do
   end subroutine fill_ghost_hierarchy

   subroutine restrict_pressure(fine, coarse)
      ! Restrict pressure from fine to coarse using r x r cell averaging.
      !
      ! Args:
      !     fine: Fine patch supplying high-resolution pressure.
      !     coarse: Coarse patch whose pressure cells are overwritten.
      !
      ! Vars:
      !     r: Refinement ratio derived from grid spacings.
      !     i_f, j_f: Fine-grid cell indices matching each coarse cell.
      !     fine_x_end, fine_y_end: Spatial extent of the fine patch.

      type(Patch), intent(in) :: fine
      type(Patch), intent(inout) :: coarse

      integer :: i_c, j_c, i_f, j_f, r, ii, jj
      real(8) :: dx_c, dy_c, dx_f, dy_f, xc, yc, avg
      real(8) :: fine_x_end, fine_y_end

      dx_c = coarse%dx; dy_c = coarse%dy
      dx_f = fine%dx; dy_f = fine%dy
      r = nint(dx_c / dx_f)
      fine_x_end = patch_x_end(fine)
      fine_y_end = patch_y_end(fine)

      do i_c = 1, coarse%nx
         xc = coarse%x0 + (i_c - 0.5d0) * dx_c
         if (xc < fine%x0 - geom_tol .or. xc > fine_x_end + geom_tol) cycle
         do j_c = 1, coarse%ny
            yc = coarse%y0 + (j_c - 0.5d0) * dy_c
            if (yc < fine%y0 - geom_tol .or. yc > fine_y_end + geom_tol) cycle
            i_f = nint((xc - dx_c * 0.5d0 - fine%x0) / dx_f) + 1
            j_f = nint((yc - dy_c * 0.5d0 - fine%y0) / dy_f) + 1
            avg = 0.0d0
            do jj = 0, r - 1
               do ii = 0, r - 1
                  avg = avg + fine%p(j_f + jj, i_f + ii)
               end do
            end do
            coarse%p(j_c, i_c) = avg / dble(r * r)
         end do
      end do
   end subroutine restrict_pressure

   subroutine restrict_velocity(fine, coarse)
      ! Restrict velocity from fine to coarse. Uses r-in-y averaging for u
      ! and r-in-x averaging for v at coincident face locations.
      !
      ! Args:
      !     fine: Fine patch supplying high-resolution velocity.
      !     coarse: Coarse patch whose velocity faces are overwritten.
      !
      ! Vars:
      !     r: Refinement ratio derived from grid spacings.
      !     i_f, j_f: Fine-grid face indices matching each coarse face.
      !     fine_x_end, fine_y_end: Spatial extent of the fine patch.

      type(Patch), intent(in) :: fine
      type(Patch), intent(inout) :: coarse

      integer :: i_c, j_c, i_f, j_f, r, kk
      real(8) :: dx_c, dy_c, dx_f, dy_f, xc, yc, avg
      real(8) :: fine_x_end, fine_y_end

      dx_c = coarse%dx; dy_c = coarse%dy
      dx_f = fine%dx; dy_f = fine%dy
      r = nint(dx_c / dx_f)
      fine_x_end = patch_x_end(fine)
      fine_y_end = patch_y_end(fine)

      ! restrict u: average r fine values in y at coincident x-faces
      do i_c = 1, coarse%nx - 1
         xc = coarse%x0 + i_c * dx_c
         if (xc < fine%x0 + geom_tol .or. xc > fine_x_end - geom_tol) cycle
         i_f = nint((xc - fine%x0) / dx_f)
         do j_c = 1, coarse%ny
            yc = coarse%y0 + (j_c - 0.5d0) * dy_c
            if (yc < fine%y0 - geom_tol .or. yc > fine_y_end + geom_tol) cycle
            j_f = nint((coarse%y0 + (j_c - 1) * dy_c - fine%y0) / dy_f) + 1
            avg = 0.0d0
            do kk = 0, r - 1
               avg = avg + fine%u(j_f + kk, i_f)
            end do
            coarse%u(j_c, i_c) = avg / dble(r)
         end do
      end do

      ! restrict v: average r fine values in x at coincident y-faces
      do j_c = 1, coarse%ny - 1
         yc = coarse%y0 + j_c * dy_c
         if (yc < fine%y0 + geom_tol .or. yc > fine_y_end - geom_tol) cycle
         j_f = nint((yc - fine%y0) / dy_f)
         do i_c = 1, coarse%nx
            xc = coarse%x0 + (i_c - 0.5d0) * dx_c
            if (xc < fine%x0 - geom_tol .or. xc > fine_x_end + geom_tol) cycle
            i_f = nint((coarse%x0 + (i_c - 1) * dx_c - fine%x0) / dx_f) + 1
            avg = 0.0d0
            do kk = 0, r - 1
               avg = avg + fine%v(j_f, i_f + kk)
            end do
            coarse%v(j_c, i_c) = avg / dble(r)
         end do
      end do
   end subroutine restrict_velocity

   subroutine prolong_pressure(coarse, fine)
      ! Prolong pressure from coarse to fine interior via bilinear interpolation.
      !
      ! Args:
      !     coarse: Coarse patch used as interpolation source.
      !     fine: Fine patch whose pressure interior is written.

      type(Patch), intent(in) :: coarse
      type(Patch), intent(inout) :: fine

      integer :: i, j
      real(8) :: xp, yp
      real(8) :: c_pox, c_poy
      integer :: nj_cp, ni_cp

      c_pox = coarse%x0 - 0.5d0 * coarse%dx
      c_poy = coarse%y0 - 0.5d0 * coarse%dy
      nj_cp = coarse%ny + 2; ni_cp = coarse%nx + 2

      do i = 1, fine%nx
         xp = fine%x0 + (i - 0.5d0) * fine%dx
         do j = 1, fine%ny
            yp = fine%y0 + (j - 0.5d0) * fine%dy
            fine%p(j, i) = interpolate_value(coarse%p, nj_cp, ni_cp, &
               c_pox, c_poy, coarse%dx, coarse%dy, xp, yp)
         end do
      end do
   end subroutine prolong_pressure

   subroutine prolong_velocity(coarse, fine)
      ! Prolong velocity from coarse to fine interior via bilinear interpolation.
      !
      ! Args:
      !     coarse: Coarse patch used as interpolation source.
      !     fine: Fine patch whose velocity interior is written.

      type(Patch), intent(in) :: coarse
      type(Patch), intent(inout) :: fine

      integer :: i, j
      real(8) :: xp, yp
      real(8) :: c_uox, c_uoy, c_vox, c_voy
      integer :: nj_cu, ni_cu, nj_cv, ni_cv

      c_uox = coarse%x0
      c_uoy = coarse%y0 - 0.5d0 * coarse%dy
      nj_cu = coarse%ny + 2; ni_cu = coarse%nx + 1

      c_vox = coarse%x0 - 0.5d0 * coarse%dx
      c_voy = coarse%y0
      nj_cv = coarse%ny + 1; ni_cv = coarse%nx + 2

      ! u interior
      do i = 1, fine%nx - 1
         xp = fine%x0 + i * fine%dx
         do j = 1, fine%ny
            yp = fine%y0 + (j - 0.5d0) * fine%dy
            fine%u(j, i) = interpolate_value(coarse%u, nj_cu, ni_cu, &
               c_uox, c_uoy, coarse%dx, coarse%dy, xp, yp)
         end do
      end do

      ! v interior
      do i = 1, fine%nx
         xp = fine%x0 + (i - 0.5d0) * fine%dx
         do j = 1, fine%ny - 1
            yp = fine%y0 + j * fine%dy
            fine%v(j, i) = interpolate_value(coarse%v, nj_cv, ni_cv, &
               c_vox, c_voy, coarse%dx, coarse%dy, xp, yp)
         end do
      end do
   end subroutine prolong_velocity

   subroutine prolong_patch(coarse, fine)
      ! Prolong all fields from coarse to fine interior.
      !
      ! Args:
      !     coarse: Coarse patch used as interpolation source.
      !     fine: Fine patch whose interior is written.

      type(Patch), intent(in) :: coarse
      type(Patch), intent(inout) :: fine

      call prolong_pressure(coarse, fine)
      call prolong_velocity(coarse, fine)
   end subroutine prolong_patch

   subroutine copy_pressure(old_p, new_p)
      ! Copy pressure from an old patch to a new patch where they overlap.
      !
      ! Args:
      !     old_p: Source patch (must share dx/dy with new_p).
      !     new_p: Destination patch whose overlapping pressure cells are overwritten.

      type(Patch), intent(in) :: old_p
      type(Patch), intent(inout) :: new_p

      real(8) :: dx, dy
      real(8) :: o_pox, o_poy, n_pox, n_poy
      real(8) :: x_lo, x_hi, y_lo, y_hi
      real(8) :: old_x_end, old_y_end, new_x_end, new_y_end
      integer :: io_s, io_e, jo_s, jo_e
      integer :: in_s, in_e, jn_s, jn_e
      integer :: ii, jj

      dx = old_p%dx; dy = old_p%dy
      old_x_end = patch_x_end(old_p)
      old_y_end = patch_y_end(old_p)
      new_x_end = patch_x_end(new_p)
      new_y_end = patch_y_end(new_p)

      o_pox = old_p%x0 - 0.5d0 * dx; o_poy = old_p%y0 - 0.5d0 * dy
      n_pox = new_p%x0 - 0.5d0 * dx; n_poy = new_p%y0 - 0.5d0 * dy
      x_lo = max(old_p%x0 + 0.5d0*dx, new_p%x0 + 0.5d0*dx)
      x_hi = min(old_x_end - 0.5d0*dx, new_x_end - 0.5d0*dx)
      y_lo = max(old_p%y0 + 0.5d0*dy, new_p%y0 + 0.5d0*dy)
      y_hi = min(old_y_end - 0.5d0*dy, new_y_end - 0.5d0*dy)
      if (x_lo <= x_hi + geom_tol .and. y_lo <= y_hi + geom_tol) then
         io_s = nint((x_lo - o_pox) / dx); io_e = nint((x_hi - o_pox) / dx)
         jo_s = nint((y_lo - o_poy) / dy); jo_e = nint((y_hi - o_poy) / dy)
         in_s = nint((x_lo - n_pox) / dx); in_e = nint((x_hi - n_pox) / dx)
         jn_s = nint((y_lo - n_poy) / dy); jn_e = nint((y_hi - n_poy) / dy)
         do ii = 0, io_e - io_s
            do jj = 0, jo_e - jo_s
               new_p%p(jn_s + jj, in_s + ii) = old_p%p(jo_s + jj, io_s + ii)
            end do
         end do
      end if
   end subroutine copy_pressure

   subroutine copy_velocity(old_p, new_p)
      ! Copy velocity from an old patch to a new patch where they overlap.
      !
      ! Args:
      !     old_p: Source patch (must share dx/dy with new_p).
      !     new_p: Destination patch whose overlapping velocity faces are overwritten.

      type(Patch), intent(in) :: old_p
      type(Patch), intent(inout) :: new_p

      real(8) :: dx, dy
      real(8) :: o_ox, o_oy, n_ox, n_oy
      real(8) :: x_lo, x_hi, y_lo, y_hi
      real(8) :: old_x_end, old_y_end, new_x_end, new_y_end
      integer :: io_s, io_e, jo_s, jo_e
      integer :: in_s, in_e, jn_s, jn_e
      integer :: ii, jj

      dx = old_p%dx; dy = old_p%dy
      old_x_end = patch_x_end(old_p)
      old_y_end = patch_y_end(old_p)
      new_x_end = patch_x_end(new_p)
      new_y_end = patch_y_end(new_p)

      ! u velocity
      o_ox = old_p%x0; o_oy = old_p%y0 - 0.5d0 * dy
      n_ox = new_p%x0; n_oy = new_p%y0 - 0.5d0 * dy
      x_lo = max(old_p%x0, new_p%x0)
      x_hi = min(old_x_end, new_x_end)
      y_lo = max(old_p%y0 + 0.5d0*dy, new_p%y0 + 0.5d0*dy)
      y_hi = min(old_y_end - 0.5d0*dy, new_y_end - 0.5d0*dy)
      if (x_lo <= x_hi + geom_tol .and. y_lo <= y_hi + geom_tol) then
         io_s = nint((x_lo - o_ox) / dx); io_e = nint((x_hi - o_ox) / dx)
         jo_s = nint((y_lo - o_oy) / dy); jo_e = nint((y_hi - o_oy) / dy)
         in_s = nint((x_lo - n_ox) / dx); in_e = nint((x_hi - n_ox) / dx)
         jn_s = nint((y_lo - n_oy) / dy); jn_e = nint((y_hi - n_oy) / dy)
         do ii = 0, io_e - io_s
            do jj = 0, jo_e - jo_s
               new_p%u(jn_s + jj, in_s + ii) = old_p%u(jo_s + jj, io_s + ii)
            end do
         end do
      end if

      ! v velocity
      o_ox = old_p%x0 - 0.5d0 * dx; o_oy = old_p%y0
      n_ox = new_p%x0 - 0.5d0 * dx; n_oy = new_p%y0
      x_lo = max(old_p%x0 + 0.5d0*dx, new_p%x0 + 0.5d0*dx)
      x_hi = min(old_x_end - 0.5d0*dx, new_x_end - 0.5d0*dx)
      y_lo = max(old_p%y0, new_p%y0)
      y_hi = min(old_y_end, new_y_end)
      if (x_lo <= x_hi + geom_tol .and. y_lo <= y_hi + geom_tol) then
         io_s = nint((x_lo - o_ox) / dx); io_e = nint((x_hi - o_ox) / dx)
         jo_s = nint((y_lo - o_oy) / dy); jo_e = nint((y_hi - o_oy) / dy)
         in_s = nint((x_lo - n_ox) / dx); in_e = nint((x_hi - n_ox) / dx)
         jn_s = nint((y_lo - n_oy) / dy); jn_e = nint((y_hi - n_oy) / dy)
         do ii = 0, io_e - io_s
            do jj = 0, jo_e - jo_s
               new_p%v(jn_s + jj, in_s + ii) = old_p%v(jo_s + jj, io_s + ii)
            end do
         end do
      end if
   end subroutine copy_velocity

   subroutine copy_patch(old_p, new_p)
      ! Copy all fields from an old patch to a new patch where they overlap.
      !
      ! Args:
      !     old_p: Source patch (must share dx/dy with new_p).
      !     new_p: Destination patch whose overlapping data is overwritten.

      type(Patch), intent(in) :: old_p
      type(Patch), intent(inout) :: new_p

      call copy_pressure(old_p, new_p)
      call copy_velocity(old_p, new_p)
   end subroutine copy_patch

   ! === CELL CLUSTERING ===

   subroutine connect_components(mask, ny, nx, labels, count)
      ! Label connected components in a boolean mask using BFS with 4-connectivity.
      !
      ! Args:
      !     mask: Boolean flag array (ny x nx).
      !     ny: Number of rows in mask.
      !     nx: Number of columns in mask.
      !     labels: Output label array; each cell receives its component id.
      !     count: Total number of connected components found.
      !
      ! Vars:
      !     qj, qi: BFS queue arrays for row and column indices.
      !     dj, di: Direction offsets for 4-connectivity neighbours.

      integer, intent(in) :: ny, nx
      logical, intent(in) :: mask(ny, nx)
      integer, intent(out) :: labels(ny, nx)
      integer, intent(out) :: count

      integer :: i, j, ci, cj, ni2, nj2, d, head, tail
      integer, allocatable :: qj(:), qi(:)
      integer, parameter :: dj(4) = (/ -1, 1, 0, 0 /)
      integer, parameter :: di(4) = (/ 0, 0, -1, 1 /)

      labels = 0
      count = 0
      allocate(qj(nx * ny), qi(nx * ny))

      do j = 1, ny
         do i = 1, nx
            if (mask(j, i) .and. labels(j, i) == 0) then
               count = count + 1
               head = 1; tail = 1
               qj(1) = j; qi(1) = i
               labels(j, i) = count
               do while (head <= tail)
                  cj = qj(head); ci = qi(head)
                  head = head + 1
                  do d = 1, 4
                     nj2 = cj + dj(d); ni2 = ci + di(d)
                     if (nj2 >= 1 .and. nj2 <= ny .and. &
                        ni2 >= 1 .and. ni2 <= nx) then
                        if (mask(nj2, ni2) .and. labels(nj2, ni2) == 0) then
                           labels(nj2, ni2) = count
                           tail = tail + 1
                           qj(tail) = nj2; qi(tail) = ni2
                        end if
                     end if
                  end do
               end do
            end if
         end do
      end do

      deallocate(qj, qi)
   end subroutine connect_components

   subroutine cluster_flagged(coarse, flagged, ny_f, nx_f, ratio, buffer, &
      boxes, n_boxes)
      ! Convert a boolean flag mask into fine-patch bounding boxes.
      !
      ! Args:
      !     coarse: Parent coarse patch used for coordinate mapping.
      !     flagged: Boolean flag mask (ny_f x nx_f).
      !     ny_f: Number of rows in the flag mask.
      !     nx_f: Number of columns in the flag mask.
      !     ratio: Refinement ratio (coarse-to-fine cells).
      !     buffer: Number of coarse cells to pad around each cluster.
      !     boxes: Output array; each column holds (x0, y0, nx, ny, dx, dy).
      !     n_boxes: Number of bounding boxes produced.
      !
      ! Vars:
      !     labels: Connected-component label array from BFS.
      !     bb: Bounding-box extents in 0-based coarse-cell indices.
      !     dx_fine, dy_fine: Fine-grid cell spacings.

      type(Patch), intent(in) :: coarse
      integer, intent(in) :: ny_f, nx_f, ratio, buffer
      logical, intent(in) :: flagged(ny_f, nx_f)
      real(8), intent(out) :: boxes(6, max_patches)
      integer, intent(out) :: n_boxes

      integer :: labels(ny_f, nx_f), count
      integer :: bb(4, max_patches)
      integer :: lbl, i, j, a, b
      integer :: i_min, i_max, j_min, j_max
      integer :: a_imin, a_jmin, a_imax, a_jmax, b_imin, b_jmin, b_imax, b_jmax
      logical :: any_flagged, did_merge
      real(8) :: dx_fine, dy_fine

      ! check if any cell is flagged
      any_flagged = .false.
      do j = 1, ny_f
         do i = 1, nx_f
            if (flagged(j, i)) then
               any_flagged = .true.
               exit
            end if
         end do
         if (any_flagged) exit
      end do
      if (.not. any_flagged) then
         n_boxes = 0
         return
      end if

      call connect_components(flagged, ny_f, nx_f, labels, count)

      ! build one bounding box per component (0-based coarse-cell indices)
      n_boxes = 0
      do lbl = 1, count
         i_min = nx_f; i_max = 0
         j_min = ny_f; j_max = 0
         do j = 1, ny_f
            do i = 1, nx_f
               if (labels(j, i) == lbl) then
                  if (i - 1 < i_min) i_min = i - 1
                  if (i - 1 > i_max) i_max = i - 1
                  if (j - 1 < j_min) j_min = j - 1
                  if (j - 1 > j_max) j_max = j - 1
               end if
            end do
         end do

         ! expand by buffer, clip to patch extent
         i_min = max(i_min - buffer, 0)
         i_max = min(i_max + buffer, coarse%nx - 1)
         j_min = max(j_min - buffer, 0)
         j_max = min(j_max + buffer, coarse%ny - 1)
         n_boxes = n_boxes + 1
         bb(1, n_boxes) = i_min
         bb(2, n_boxes) = j_min
         bb(3, n_boxes) = i_max
         bb(4, n_boxes) = j_max
      end do

      ! merge overlapping boxes
      did_merge = .true.
      do while (did_merge)
         did_merge = .false.
         outer: do a = 1, n_boxes
            do b = a + 1, n_boxes
               a_imin = bb(1, a); a_jmin = bb(2, a); a_imax = bb(3, a); a_jmax = bb(4, a)
               b_imin = bb(1, b); b_jmin = bb(2, b); b_imax = bb(3, b); b_jmax = bb(4, b)
               if (a_imin <= b_imax .and. b_imin <= a_imax .and. &
                  a_jmin <= b_jmax .and. b_jmin <= a_jmax) then
                  bb(1, a) = min(a_imin, b_imin)
                  bb(2, a) = min(a_jmin, b_jmin)
                  bb(3, a) = max(a_imax, b_imax)
                  bb(4, a) = max(a_jmax, b_jmax)
                  do i = b, n_boxes - 1
                     bb(:, i) = bb(:, i + 1)
                  end do
                  n_boxes = n_boxes - 1
                  did_merge = .true.
                  exit outer
               end if
            end do
         end do outer
      end do

      ! convert to fine-grid dimensions
      dx_fine = coarse%dx / dble(ratio)
      dy_fine = coarse%dy / dble(ratio)
      do a = 1, n_boxes
         boxes(1, a) = coarse%x0 + bb(1, a) * coarse%dx
         boxes(2, a) = coarse%y0 + bb(2, a) * coarse%dy
         boxes(3, a) = dble((bb(3, a) - bb(1, a) + 1) * ratio)
         boxes(4, a) = dble((bb(4, a) - bb(2, a) + 1) * ratio)
         boxes(5, a) = dx_fine
         boxes(6, a) = dy_fine
      end do
   end subroutine cluster_flagged

end module grid_mod
