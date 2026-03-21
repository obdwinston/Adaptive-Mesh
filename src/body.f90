module body_mod
   use config_mod
   use grid_mod
   implicit none
   private

   type, public :: Body
      integer :: n_pts                      ! number of polygon points
      real(8), allocatable :: xb(:), yb(:)  ! polygon coordinates (1:n_pts)
      real(8) :: cx, cy                     ! centroid after translation
      real(8) :: half_w, half_h             ! bounding box half-extents
   end type Body

   public :: read_body, apply_force_patch, apply_force_hierarchy

contains

   ! === BODY GEOMETRY ===

   subroutine read_body(b)
      ! Read polygon coordinates from a text file and scale/translate them.
      !
      ! Args:
      !     b: Body structure to populate with scaled/translated polygon.
      !
      ! Vars:
      !     D: Characteristic length (max of width, height) before scaling.
      !     sc: Scale factor to normalize body to config side length.

      type(Body), intent(inout) :: b

      integer :: iu, ios, n, k
      real(8) :: xt, yt
      real(8) :: x_min, x_max, y_min, y_max, w, h, D, sc

      open(newunit=iu, file=body_file, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         write(*, '(A,A)') 'ERROR: cannot open body file ', trim(body_file)
         stop 1
      end if

      ! count lines
      n = 0
      do
         read(iu, *, iostat=ios) xt, yt
         if (ios /= 0) exit
         n = n + 1
      end do
      rewind(iu)

      b%n_pts = n
      if (allocated(b%xb)) deallocate(b%xb)
      if (allocated(b%yb)) deallocate(b%yb)
      allocate(b%xb(n), b%yb(n))

      do k = 1, n
         read(iu, *) b%xb(k), b%yb(k)
      end do
      close(iu)

      ! bounding box of raw coordinates
      x_min = minval(b%xb); x_max = maxval(b%xb)
      y_min = minval(b%yb); y_max = maxval(b%yb)
      w = x_max - x_min; h = y_max - y_min
      D = max(w, h)

      ! scale so characteristic length = side
      sc = side / D
      b%xb = b%xb * sc
      b%yb = b%yb * sc

      ! translate centroid of bounding box to (cx_body, cy_body)
      x_min = minval(b%xb); x_max = maxval(b%xb)
      y_min = minval(b%yb); y_max = maxval(b%yb)
      b%xb = b%xb - 0.5d0 * (x_min + x_max) + cx_body
      b%yb = b%yb - 0.5d0 * (y_min + y_max) + cy_body

      b%cx = cx_body; b%cy = cy_body
      b%half_w = 0.5d0 * (maxval(b%xb) - minval(b%xb))
      b%half_h = 0.5d0 * (maxval(b%yb) - minval(b%yb))
   end subroutine read_body

   logical function in_polygon(x, y, b)
      ! Ray-casting point-in-polygon test.
      !
      ! Args:
      !     x: Query x-coordinate.
      !     y: Query y-coordinate.
      !     b: Body polygon to test against.

      real(8), intent(in) :: x, y
      type(Body), intent(in) :: b

      integer :: i, j, n
      logical :: inside

      n = b%n_pts
      inside = .false.
      j = n
      do i = 1, n
         if (((b%yb(i) > y) .neqv. (b%yb(j) > y)) .and. &
            (x < (b%xb(j) - b%xb(i)) * (y - b%yb(i)) / &
            (b%yb(j) - b%yb(i)) + b%xb(i))) then
            inside = .not. inside
         end if
         j = i
      end do
      in_polygon = inside
   end function in_polygon

   ! === BODY FORCING ===

   subroutine apply_force_patch(pa, bd, predictor)
      ! Zero velocity at all u-face and v-face positions inside a solid body.
      !
      ! Args:
      !     pa: Patch to apply forcing on.
      !     bd: Body polygon defining the solid region.
      !     predictor: If true, operates on u_star/v_star instead of u/v.
      !
      ! Vars:
      !     hw, hh: Half-width and half-height of body bounding box.
      !     xu, yu: Physical position of each u-face.
      !     xv, yv: Physical position of each v-face.

      type(Patch), intent(inout), target :: pa
      type(Body), intent(in) :: bd
      logical, intent(in) :: predictor

      real(8), pointer :: u(:,:), v(:,:)
      integer :: i, j
      real(8) :: hw, hh, xu, yu, xv, yv

      if (predictor) then
         u => pa%u_star; v => pa%v_star
      else
         u => pa%u; v => pa%v
      end if

      hw = bd%half_w; hh = bd%half_h

      ! zero u-faces inside the body
      do i = 0, pa%nx
         xu = pa%x0 + i * pa%dx
         do j = 0, pa%ny + 1
            yu = pa%y0 + (j - 0.5d0) * pa%dy
            if (xu < bd%cx - hw - geom_tol .or. xu > bd%cx + hw + geom_tol .or. &
               yu < bd%cy - hh - geom_tol .or. yu > bd%cy + hh + geom_tol) cycle
            if (in_polygon(xu, yu, bd)) u(j, i) = 0.0d0
         end do
      end do

      ! zero v-faces inside the body
      do i = 0, pa%nx + 1
         xv = pa%x0 + (i - 0.5d0) * pa%dx
         do j = 0, pa%ny
            yv = pa%y0 + j * pa%dy
            if (xv < bd%cx - hw - geom_tol .or. xv > bd%cx + hw + geom_tol .or. &
               yv < bd%cy - hh - geom_tol .or. yv > bd%cy + hh + geom_tol) cycle
            if (in_polygon(xv, yv, bd)) v(j, i) = 0.0d0
         end do
      end do
   end subroutine apply_force_patch

   subroutine apply_force_hierarchy(hier, bd, predictor)
      ! Zero velocity inside a solid body on all patches in the hierarchy.
      !
      ! Args:
      !     hier: Grid hierarchy whose patches receive body forcing.
      !     bd: Body polygon defining the solid region.
      !     predictor: If true, operates on u_star/v_star instead of u/v.

      type(Hierarchy), intent(inout) :: hier
      type(Body), intent(in) :: bd
      logical, intent(in) :: predictor

      integer :: lev, pidx

      do lev = 1, hier%n_levels
         do pidx = 1, hier%levels(lev)%n_patches
            call apply_force_patch(hier%levels(lev)%patches(pidx), bd, predictor)
         end do
      end do
   end subroutine apply_force_hierarchy

end module body_mod
