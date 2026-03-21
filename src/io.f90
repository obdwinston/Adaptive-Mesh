module io_mod
   use config_mod
   use grid_mod
   implicit none
   private

   public :: write_output_patch, write_output_hierarchy, write_output_collection

contains

   ! === PATCH OUTPUT ===

   subroutine write_output_patch(pa, patch_id, step)
      ! Write a single patch to a VTK XML StructuredGrid (.vts) file containing
      ! cell-centre coordinates and point data (u, v, speed, pressure, vorticity).
      !
      ! Args:
      !     pa: Patch to write.
      !     patch_id: Global patch identifier used in the filename.
      !     step: Timestep number used in the filename.
      !
      ! Vars:
      !     u_cc, v_cc: Cell-centred velocity components interpolated from faces.
      !     dvdx, dudy: Velocity derivatives used for vorticity computation.

      type(Patch), intent(in) :: pa
      integer, intent(in) :: patch_id, step

      character(len=256) :: filename
      integer :: iu, i, j, nx, ny
      real(8) :: xc, yc, dvdx, dudy
      real(8), allocatable :: u_cc(:,:), v_cc(:,:)

      nx = pa%nx; ny = pa%ny

      allocate(u_cc(ny, nx), v_cc(ny, nx))
      do j = 1, ny
         do i = 1, nx
            u_cc(j, i) = 0.5d0 * (pa%u(j, i - 1) + pa%u(j, i))
            v_cc(j, i) = 0.5d0 * (pa%v(j - 1, i) + pa%v(j, i))
         end do
      end do

      write(filename, '(A,"/step_",I5.5,"_patch_",I2.2,".vts")') &
         trim(output_dir), step, patch_id

      open(newunit=iu, file=trim(filename), status='replace', action='write')

      write(iu, '(A)') '<?xml version="1.0"?>'
      write(iu, '(A,I0,A,I0,A)') &
         '<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">'
      write(iu, '(A,I0,A,I0,A)') &
         '  <StructuredGrid WholeExtent="0 ', nx-1, ' 0 ', ny-1, ' 0 0">'
      write(iu, '(A,I0,A,I0,A)') &
         '    <Piece Extent="0 ', nx-1, ' 0 ', ny-1, ' 0 0">'

      ! field data (metadata)
      write(iu, '(A)') '      <FieldData>'
      write(iu, '(A)') '        <DataArray type="Int32" Name="level" NumberOfTuples="1" format="ascii">'
      write(iu, '(I0)') pa%level
      write(iu, '(A)') '        </DataArray>'
      write(iu, '(A)') '        <DataArray type="Float64" Name="dx" NumberOfTuples="1" format="ascii">'
      write(iu, '(ES18.8E3)') pa%dx
      write(iu, '(A)') '        </DataArray>'
      write(iu, '(A)') '        <DataArray type="Float64" Name="dy" NumberOfTuples="1" format="ascii">'
      write(iu, '(ES18.8E3)') pa%dy
      write(iu, '(A)') '        </DataArray>'
      write(iu, '(A)') '      </FieldData>'

      ! points (cell centres)
      write(iu, '(A)') '      <Points>'
      write(iu, '(A)') '        <DataArray type="Float64" ' // &
         'NumberOfComponents="3" format="ascii">'
      do j = 1, ny
         do i = 1, nx
            xc = pa%x0 + (i - 0.5d0) * pa%dx
            yc = pa%y0 + (j - 0.5d0) * pa%dy
            write(iu, '(3ES18.8E3)') xc, yc, 0.0d0
         end do
      end do
      write(iu, '(A)') '        </DataArray>'
      write(iu, '(A)') '      </Points>'

      write(iu, '(A)') '      <PointData>'

      write(iu, '(A)') '        <DataArray type="Float64" Name="u_cc" format="ascii">'
      do j = 1, ny
         do i = 1, nx
            write(iu, '(ES18.8E3)') u_cc(j, i)
         end do
      end do
      write(iu, '(A)') '        </DataArray>'

      write(iu, '(A)') '        <DataArray type="Float64" Name="v_cc" format="ascii">'
      do j = 1, ny
         do i = 1, nx
            write(iu, '(ES18.8E3)') v_cc(j, i)
         end do
      end do
      write(iu, '(A)') '        </DataArray>'

      write(iu, '(A)') '        <DataArray type="Float64" Name="speed" format="ascii">'
      do j = 1, ny
         do i = 1, nx
            write(iu, '(ES18.8E3)') sqrt(u_cc(j, i)**2 + v_cc(j, i)**2)
         end do
      end do
      write(iu, '(A)') '        </DataArray>'

      write(iu, '(A)') '        <DataArray type="Float64" Name="pressure" format="ascii">'
      do j = 1, ny
         do i = 1, nx
            write(iu, '(ES18.8E3)') pa%p(j, i)
         end do
      end do
      write(iu, '(A)') '        </DataArray>'

      ! vorticity (dv/dx - du/dy)
      write(iu, '(A)') '        <DataArray type="Float64" Name="vorticity" format="ascii">'
      do j = 1, ny
         do i = 1, nx
            if (i > 1 .and. i < nx) then
               dvdx = (v_cc(j, i + 1) - v_cc(j, i - 1)) / (2.0d0 * pa%dx)
            else if (i == 1) then
               dvdx = (v_cc(j, i + 1) - v_cc(j, i)) / pa%dx
            else
               dvdx = (v_cc(j, i) - v_cc(j, i - 1)) / pa%dx
            end if

            if (j > 1 .and. j < ny) then
               dudy = (u_cc(j + 1, i) - u_cc(j - 1, i)) / (2.0d0 * pa%dy)
            else if (j == 1) then
               dudy = (u_cc(j + 1, i) - u_cc(j, i)) / pa%dy
            else
               dudy = (u_cc(j, i) - u_cc(j - 1, i)) / pa%dy
            end if

            write(iu, '(ES18.8E3)') dvdx - dudy
         end do
      end do
      write(iu, '(A)') '        </DataArray>'

      write(iu, '(A)') '      </PointData>'
      write(iu, '(A)') '    </Piece>'
      write(iu, '(A)') '  </StructuredGrid>'
      write(iu, '(A)') '</VTKFile>'

      close(iu)
      deallocate(u_cc, v_cc)
   end subroutine write_output_patch

   ! === HIERARCHY OUTPUT ===

   subroutine write_output_hierarchy(hier, step)
      ! Write one .vts file for every patch in the hierarchy.
      !
      ! Args:
      !     hier: Grid hierarchy to write.
      !     step: Timestep number used in filenames.

      type(Hierarchy), intent(in) :: hier
      integer, intent(in) :: step

      integer :: lev, pidx, global_id

      global_id = 0
      do lev = 1, hier%n_levels
         do pidx = 1, hier%levels(lev)%n_patches
            global_id = global_id + 1
            call write_output_patch(hier%levels(lev)%patches(pidx), &
               global_id, step)
         end do
      end do
   end subroutine write_output_hierarchy

   subroutine write_output_collection(steps, times, n_steps, n_patches_at_step)
      ! Write a ParaView collection (.pvd) file referencing all written timesteps.
      !
      ! Args:
      !     steps: Array of timestep indices that were written.
      !     times: Array of simulation times corresponding to each step.
      !     n_steps: Number of entries in steps and times.
      !     n_patches_at_step: Number of patches written at each step.

      integer, intent(in) :: n_steps
      integer, intent(in) :: steps(n_steps)
      real(8), intent(in) :: times(n_steps)
      integer, intent(in) :: n_patches_at_step(n_steps)

      character(len=256) :: pvd_file, vts_file
      integer :: iu, s, pid

      write(pvd_file, '(A,"/solution.pvd")') trim(output_dir)
      open(newunit=iu, file=trim(pvd_file), status='replace', action='write')

      write(iu, '(A)') '<?xml version="1.0"?>'
      write(iu, '(A)') '<VTKFile type="Collection" version="0.1">'
      write(iu, '(A)') '  <Collection>'

      do s = 1, n_steps
         do pid = 1, n_patches_at_step(s)
            write(vts_file, '("step_",I5.5,"_patch_",I2.2,".vts")') steps(s), pid
            write(iu, '(A,ES14.6E3,A,A,A)') &
               '    <DataSet timestep="', times(s), &
               '" file="', trim(vts_file), '"/>'
         end do
      end do

      write(iu, '(A)') '  </Collection>'
      write(iu, '(A)') '</VTKFile>'
      close(iu)
   end subroutine write_output_collection

end module io_mod
