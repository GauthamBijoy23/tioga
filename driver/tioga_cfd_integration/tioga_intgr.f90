module tioga_overall

  use gridtype
  implicit none
  include 'mpif.h'

  type(grid), target, allocatable :: gr(:)
  type(grid), pointer :: g
  integer :: myid,numprocs,ierr,nsave
  integer :: itime,ntimesteps,iter,nsubiter
  real*8 :: t0
  integer :: blockid
  logical :: iclip,exists
  integer :: ntypes
  integer :: nv1,nv2,nv3
  real*8 :: t1,t2
  integer :: i,n,m,ib,j,nblocks
  real*8 :: xt(3),rnorm
  integer :: dcount,fcount
  integer, allocatable :: receptorInfo(:),inode(:)
  real*8, allocatable :: frac(:)
  character*6 :: integer_string
  character*64 :: fname

contains

subroutine initialize_mpi()
  
  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world,myid,ierr)
  call mpi_comm_size(mpi_comm_world,numprocs,ierr)
  if (myid==0) write(6,*) '# tioga test on ',numprocs,' processes'

end subroutine
  
subroutine read_grids()
 
  nblocks=0
  exists=.true.
  do while(exists)
   write(integer_string,"(I6)") 100000+100*myid+nblocks
   fname='input/flow'//trim(adjustl(integer_string(2:)))//'.dat'
   !write(6,*) fname
   inquire(FILE=fname,EXIST=exists) 
   if (exists) nblocks=nblocks+1
  enddo
  !write(6,*) myid,nblocks
  allocate(gr(nblocks))
  do ib=1,nblocks
     write(integer_string,"(I6)") 100000+100*myid+ib-1
     fname='input/flow'//trim(adjustl(integer_string(2:)))//'.dat'
     call readGrid_from_file(fname,gr(ib))
  enddo
  !
  if (myid==0) write(6,*) '# tioga test : finished reading grids '

end subroutine

!===================================================================
subroutine initialize_tioga()
!===================================================================
  
  call tioga_init_f90(mpi_comm_world)
  call mpi_barrier(mpi_comm_world,ierr)
    
  ntypes=1
  nv1=6
  nv2=8
  nv3=4     !NUMBER OF VERTICES PER TET CELL
  do ib=1,nblocks
   g=>gr(ib)
   if (g%n6 > 0)  then
    call tioga_registergrid_data_mb(ib,g%bodytag(1),g%nv,g%x,g%iblank,g%nwbc,g%nobc,g%wbcnode,g%obcnode,&
       ntypes,nv1,g%n6,g%ndc6)
   else if (g%n8 > 0) then
    call tioga_registergrid_data_mb(ib,g%bodytag(1),g%nv,g%x,g%iblank,g%nwbc,g%nobc,g%wbcnode,g%obcnode,&
       ntypes,nv2,g%n8,g%ndc8)
   !-----------------------------------------------------------------------------------------------------
   !   REGISTERING TET FILE GRID DATA
   !
   else if (g%n4 > 0) then
    call tioga_registergrid_data_mb(ib,g%bodytag(1),g%nv,g%x,g%iblank,g%nwbc,g%nobc,g%wbcnode,g%obcnode,&
       ntypes,nv3,g%n4,g%ndc4)
   !
   !------------------------------------------------------------------------------------------------------
   endif
  enddo
end subroutine

subroutine tioga_connectivity()
  call tioga_preprocess_grids                  !< preprocess the grids (call again if dynamic) 
  call mpi_barrier(mpi_comm_world,ierr)

  call cpu_time(t1)         
  call tioga_performconnectivity               !< determine iblanking and interpolation patterns
  call tioga_reduce_fringes()
  call cpu_time(t2)

  call mpi_barrier(mpi_comm_world,ierr)
  if (myid==0) write(6,*) 'connectivity time=',t2-t1
end subroutine 
subroutine tioga_registerandupdate()
  do ib=1,nblocks
    g=>gr(ib)
    call tioga_registersolution(g%bodytag(1),g%q)
  enddo
  call tioga_dataupdate_mb(gr(1)%nvar,'row')    !< update the q-variables (can be called anywhere)
                                             !< nvar = number of field variables per node
                                             !< if fields are different arrays, you can also 
                                             !< call this multiple times for each field
  call mpi_barrier(mpi_comm_world,ierr)
    call cpu_time(t2)
  if (myid==0) write(6,*) 'data update time=',t2-t1
end subroutine


subroutine tioga_finalise()

  call mpi_barrier(mpi_comm_world,ierr)
  call tioga_delete
  call mpi_finalize(ierr)
end subroutine 
end module tioga_overall

program solver

implicit none
use module tioga_overall
integer :: ntimesteps
ntimesteps = 100
call initialize_mpi()
call read_grids()
call initialize_tioga()
  do timestep = 1, ntimesteps
      
     call tioga_connectivity()
     ! call cfd_solve_step()

     call tioga_registerandupdate()

  enddo

call tioga_finalise()
end program solver
