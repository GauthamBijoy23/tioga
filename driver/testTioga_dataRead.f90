program testTioga
  !
  use gridtype
  !
  use run_tioga, only : tioga_init_conn, tioga_solutions, tioga_fin, gr	
  !
  implicit none
  !
  include 'mpif.h'
  !
!  type(grid), target, allocatable :: gr(:)
!  type(grid), pointer :: g
!  integer :: myid,numprocs,ierr,nsave,runtime
!  integer :: itime,ntimesteps,iter,nsubiter
!  real*8 :: t0
!  integer :: blockid
!  logical :: iclip,exists
!  integer :: ntypes, t
!  integer :: nv1,nv2,nv3,ios1,ios2,rbase,ibase
!  real*8 :: t1,t2
!  integer :: i,n,m,ib,j,nblocks
!  real*8 :: xt(3),rnorm
!  integer :: dcount,fcount
!  integer, allocatable :: receptorInfo(:),inode(:)
!  real*8, allocatable :: frac(:)
!  character*6 :: integer_string
!  character*64 :: fname
!  character(len=20) :: t_str
  
call tioga_init_conn('/home/gautham_linux/tioga/case/input')
call tioga_solutions(gr(1)%q, gr(1)%nvar,gr(1)%nv)
call tioga_fin()


end program testTioga
