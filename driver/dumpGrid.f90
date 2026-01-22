subroutine dump_grid_pre_connectivity(nnodes, coords, ncells, conn, cellType, iblank)
  use mpi
  implicit none

  integer, intent(in) :: nnodes, ncells
  real(8), intent(in) :: coords(3, nnodes)
  integer, intent(in) :: conn(:, :)
  integer, intent(in) :: cellType(ncells)
  integer, intent(in) :: iblank(ncells)

  integer :: i, c, v
  integer :: rank, ierr, unit
  character(len=128) :: fname

  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

  write(fname, '(A,I0,A)') 'grid_pre_holecut_rank_', rank, '.dat'
  unit = 100 + rank

  open(unit=unit, file=fname, status='replace', action='write')

  ! ---- Nodes ----
  write(unit,*) 'NODES ', nnodes
  do i = 1, nnodes
     write(unit,'(I8,3ES20.12)') i, coords(1,i), coords(2,i), coords(3,i)
  end do

  ! ---- Cells ----
  write(unit,*) 'CELLS ', ncells
  do c = 1, ncells
     write(unit,'(I8,100I8)') c, (conn(v,c), v=1,size(conn,1))
  end do

  ! ---- Cell types ----
  write(unit,*) 'CELLTYPES'
  do c = 1, ncells
     write(unit,'(I8,I8)') c, cellType(c)
  end do

  ! ---- Iblank ----
  write(unit,*) 'IBLANK'
  do c = 1, ncells
     write(unit,'(I8,I8)') c, iblank(c)
  end do

  close(unit)

end subroutine dump_grid_pre_connectivity
