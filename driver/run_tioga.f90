!
! Tioga is a library for overset grid assembly on parallel distributed systems
! Copyright (C) 2015 Jay Sitaraman
!
!>
!> Test code to integrate TIOGA into CFD solver (Based on testTioga_dataRead.f90)
!> 
!> Gautham Bijoy
!>
!> 17/02/2026
!>
module run_tioga
  !
  use gridtype
  use iso_c_binding, only: c_int
  !
  implicit none
  !
  include 'mpif.h'
  !
  type(grid), target, allocatable :: gr(:)
  type(grid), pointer :: g
  integer :: myid,numprocs,ierr,nsave,runtime
  integer :: itime,ntimesteps,iter,nsubiter
  real*8 :: t0
  integer :: blockid
  logical :: iclip,exists
  integer :: ntypes, t
  integer :: nv1,nv2,nv3,ios1,rbase,ibase
  real*8 :: t1,t2
  integer :: i,n,m,ib,j,nblocks,d
  real*8 :: xt(3),rnorm
  integer :: dcount,fcount
  integer, allocatable :: receptorInfo(:),inode(:)
  real*8, allocatable :: frac(:)
  character*6 :: integer_string
  character*64 :: fname  
  character(len=20) :: t_str
  integer :: nc,nc_t                                ! Total number of cells - allocating iblank cell array
  integer, allocatable :: iblankcells(:), ibc(:)
  integer, allocatable :: iblanknodes(:)
  real(8), allocatable :: cellvals(:,:)
  integer, allocatable, save :: donor_receptorInfo(:), donor_inode(:), donor_dcount(:), donor_roff(:), donor_foff(:), node_off(:)
real*8,  allocatable, save :: donor_frac(:)

!=============================================================
!Variables for higher order Cell-Based Interpolation (mod10f)
!=============================================================

integer(c_int), bind(C, name='tioga_cb_meshtag') :: cb_tag
integer, save :: cb_neq = 0
real(8), allocatable, save :: qrow(:)
integer, save :: is_CellBased = 1     ! 0 = node-based (IDW), 1 = cell-based Interpolation  (mod10f)
  
contains 

!===============================================================================
!Subroutines for Cell Based interpolation
!===============================================================================

integer function cb_blk()
  integer :: ib
  cb_blk = 0
  do ib = 1, nblocks
    if (gr(ib)%bodytag(1) == cb_tag) then
      cb_blk = ib
      return
    end if
  end do
  write(*,*) 'cb_blk: no block matches tag', cb_tag
  stop 1
end function

subroutine cb_nodes_per_cell(cell, npts)
  integer :: cell, npts
  npts = 1
end subroutine

subroutine cb_receptor_nodes(cell, npts, xyz)
  integer :: cell, npts, ib, m
  real(8) :: xyz(3)
  ib = cb_blk()
  xyz = 0.0d0
  do m = 1, 4
    xyz = xyz + 0.25d0*gr(ib)%x(3*gr(ib)%ndc4(m,cell)-2:3*gr(ib)%ndc4(m,cell))
  end do
end subroutine

subroutine cb_inclusion(cell, xp, passflag, rst)
  integer :: cell, passflag, ib, m
  real(8) :: xp(3), rst(3), v(3,4), a(3,3), b(3), det, l(4)
  ib = cb_blk()
  do m = 1, 4
    v(:,m) = gr(ib)%x(3*gr(ib)%ndc4(m,cell)-2:3*gr(ib)%ndc4(m,cell))
  end do
  a(:,1) = v(:,1)-v(:,4); a(:,2) = v(:,2)-v(:,4); a(:,3) = v(:,3)-v(:,4)
  b = xp - v(:,4)
  det  = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) - a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1)) + a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
  l(1) = (b(1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) - a(1,2)*(b(2)*a(3,3)-a(2,3)*b(3)) + a(1,3)*(b(2)*a(3,2)-a(2,2)*b(3)))/det
  l(2) = (a(1,1)*(b(2)*a(3,3)-a(2,3)*b(3)) - b(1)*(a(2,1)*a(3,3)-a(2,3)*a(3,1)) + a(1,3)*(a(2,1)*b(3)-b(2)*a(3,1)))/det
  l(3) = (a(1,1)*(a(2,2)*b(3)-b(2)*a(3,2)) - a(1,2)*(a(2,1)*b(3)-b(2)*a(3,1)) + b(1)*(a(2,1)*a(3,2)-a(2,2)*a(3,1)))/det
  l(4) = 1.0d0 - l(1) - l(2) - l(3)
  passflag = 0
  if (all(l >= -1.0d-8) .and. all(l <= 1.0d0+1.0d-8)) passflag = 1
  rst = 0.0d0
end subroutine

subroutine cb_donor_frac(cell, xp, nweights, inode, frac, rst, ndim)
  integer :: cell, nweights, inode(*), ndim
  real(8) :: xp(3), frac(*), rst(3)
  nweights = 1
  inode(1) = (cell-1)*cb_neq + 1
  frac(1)  = 1.0d0
  rst = 0.0d0
end subroutine

subroutine cb_to_modal(cell, npts_in, qtmp, npts, index_out, qout)
  integer :: cell, npts_in, npts, index_out
  real(8) :: qtmp(*), qout(*)
  npts = 1
  index_out = (cell-1)*cb_neq + 1
  qout(1:cb_neq) = qtmp(1:cb_neq)
end subroutine

subroutine cb_write_donormap()
  integer :: ib, d, dc, fc, rb, e, u, nfr, rcell, dcell
  integer :: coff(nblocks+1), fstart(nblocks+1)
  integer, allocatable :: ri(:), ino(:), fcell(:)
  real(8), allocatable :: fr(:)
  coff(1) = 0
  do ib = 1, nblocks
    coff(ib+1) = coff(ib) + gr(ib)%ncells
  end do
  allocate(fcell(nc_t)); nfr = 0; fstart(1) = 1
  do ib = 1, nblocks
    do e = coff(ib)+1, coff(ib+1)
      if (iblankcells(e) == -1) then
        nfr = nfr + 1; fcell(nfr) = e
      end if
    end do
    fstart(ib+1) = nfr + 1
  end do
  open(newunit=u, file='checking_nodes/donorinfo_cell.dat', status='replace')
  write(u,*) '# rank rec_pt rec_ib nw | donor_offset weight | rec_cell don_cell (global, 1-based)'
  do ib = 1, nblocks
    call tioga_getdonorcount2(gr(ib)%bodytag(1), dc, fc)
    write(u,*) 'Block Num=', gr(ib)%bodytag(1), 'Donor Count=', dc, 'Fcount=', fc
    if (dc == 0) cycle
    if (fc /= dc) stop 'donormap: expected nweights = 1'
    allocate(ri(4*dc), ino(fc), fr(fc))
    call tioga_getdonorinfo2(gr(ib)%bodytag(1), ri, ino, fr)
    do d = 1, dc
      rb    = ri(4*d-1) + 1
      rcell = fcell(fstart(rb) + ri(4*d-2))
      dcell = coff(ib) + (ino(d)-1)/cb_neq + 1
      write(u,*) ri(4*d-3:4*d), ino(d), fr(d), rcell, dcell
    end do
    deallocate(ri, ino, fr)
  end do
  close(u)
  deallocate(fcell)
end subroutine
!-----------------------CALLBACKS OVER-------------------------------

subroutine tioga_cell_conn(neq_in)
  integer, intent(in) :: neq_in
  if (is_CellBased == 1) then
    cb_neq = neq_in
    call tioga_performconnectivity_highorder()
    call cb_write_donormap()    !Writing cell donor map
  end if
end subroutine

subroutine tioga_solutions_cell(cu, neq)
  real(8), intent(inout) :: cu(:,:)
  integer, intent(in) :: neq
  integer :: ib, ie, k, off, nnanf
  character(len=64) :: finame
  real(8), allocatable :: uq_pre(:)
  integer :: nok, nunch, nbad, nnan, cl(5), cg(5)
  logical :: chg
  if (.not. allocated(qrow)) allocate(qrow(nc_t*neq))
  do ie = 1, nc_t
    do k = 1, neq
      qrow((ie-1)*neq+k) = cu(ie,k)
    end do
  end do
  nnanf = 0
  write(finame,'(A,I4.4,A)') 'checking_nodes/qcell_preinterp_', myid, '.dat'  !Qcell Pre interp (mod10h)
  open(90, file=trim(finame), status='replace')
  do ie = 1, nc_t
    write(90,*) ie, iblankcells(ie), (qrow((ie-1)*neq+k), k=1,neq)
  end do
  close(90)                                                                   !END Qcell Pre interp (mod10h)
  allocate(uq_pre(nc_t*neq)); uq_pre = qrow(1:nc_t*neq)
  off = 0
  do ib = 1, nblocks
    call tioga_registersolution(gr(ib)%bodytag(1), qrow(off*neq+1:))
    off = off + gr(ib)%ncells
  end do
  call tioga_dataupdate_mb(neq, 'row')
  write(finame,'(A,I4.4,A)') 'checking_nodes/qcell_postinterp_', myid, '.dat'  !Qcell Post interp (mod10h)
  open(99, file=trim(finame), status='replace')
  do ie = 1, nc_t
    write(99,*) ie, iblankcells(ie), (qrow((ie-1)*neq+k), k=1,neq)
  end do
  close(99)                                                                    !End Post interp write (mod10h)

  nok = 0; nunch = 0; nbad = 0; nnan = 0                                       !Checks for nan,nonfringe change,etc. (mod10x)
  do ie = 1, nc_t
    chg = any(qrow((ie-1)*neq+1:ie*neq) /= uq_pre((ie-1)*neq+1:ie*neq))
    if (iblankcells(ie) == -1) then
      if (chg) then
        nok = nok + 1
      else
        nunch = nunch + 1
      end if
      if (any(qrow((ie-1)*neq+1:ie*neq) /= qrow((ie-1)*neq+1:ie*neq))) nnan = nnan + 1
    else
      if (chg) nbad = nbad + 1
      if (iblankcells(ie) == 1) then
        if (any(qrow((ie-1)*neq+1:ie*neq) /= qrow((ie-1)*neq+1:ie*neq))) nnanf = nnanf + 1
      end if
    end if
  end do
  deallocate(uq_pre)
  cl = (/nok, nunch, nbad, nnan, nnanf/)
  call mpi_allreduce(cl, cg, 5, mpi_integer, mpi_sum, mpi_comm_world, ierr)
  if (myid == 0) print*, 'fringe updated / unchanged / nonfringe changed / fringe NaN / field NaN:', cg !End checks (mod10x)
  do ie = 1, nc_t
    if (iblankcells(ie) == -1) cu(ie,1:neq) = qrow((ie-1)*neq+1:ie*neq)
  end do
end subroutine


!===================================================================================================================

subroutine tioga_init_conn(input_dir)

character(len=*), intent(in) :: input_dir
integer :: offset2,tn,ib,i                           !offset counter for writing iblank nodes (mod8)
integer :: jj, nd                                    !Var for writing donor node xyz in getdonorinfo output (mod10h)
real(8) :: sd(4)                                     !Same as above (mod10h)
integer :: ndon_nf                                   !Checking if all donors are field (mod10i)
  !
  ! initialize mpi
  !
!  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world,myid,ierr)
  call mpi_comm_size(mpi_comm_world,numprocs,ierr)
  !
  ! read grid based on that generated by the
  ! strand/Cart generator
  !
  if (myid==0) write(6,*) '# tioga test on ',numprocs,' processes'
 
  nblocks=0
  exists=.true.
  do while(exists)
   write(integer_string,"(I6)") 100000+100*myid+nblocks
   fname=trim(input_dir)//'/flow'//trim(adjustl(integer_string(2:)))//'.dat'  ! Hardcoded inp address -> Input dir as env var
   !write(6,*) fname
   inquire(FILE=fname,EXIST=exists) 
   if (exists) nblocks=nblocks+1
  enddo
  !write(6,*) myid,nblocks
  allocate(gr(nblocks))
  do ib=1,nblocks
     write(integer_string,"(I6)") 100000+100*myid+ib-1
     fname=trim(input_dir)//'/flow'//trim(adjustl(integer_string(2:)))//'.dat' ! Hardcoded inp address -> Input dir as env var
     call readGrid_from_file(fname,gr(ib))
  enddo
  !
  if (myid==0) write(6,*) '# tioga test : finished reading grids '
  !
  ! initialize tioga
  !
  call tioga_init_f90(mpi_comm_world)
  call mpi_barrier(mpi_comm_world,ierr)
  !  
  ntypes=1
  nv1=6
  nv2=8
  nv3=4     !NUMBER OF VERTICES PER TET CELL

ndon_nf = 0          !Chcek if all donors are field (mod10i)

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

  call tioga_preprocess_grids                  !< preprocess the grids (call again if dynamic) 
  call mpi_barrier(mpi_comm_world,ierr)
  write(t_str, '(I0)') t
  call cpu_time(t1)         

!Guarded call based on is_CellBased flag (mod10f)-------------------------  
if (is_CellBased == 1) call tioga_set_highorder_callback(cb_nodes_per_cell, cb_receptor_nodes, cb_inclusion, cb_donor_frac, cb_to_modal)
call tioga_performconnectivity
if (is_CellBased == 0) call tioga_reduce_fringes()
!-------------------------------------------------------------------------

!  call tioga_performconnectivity               !< determine iblanking and interpolation patterns
!  call tioga_reduce_fringes()
  call cpu_time(t2)

tn = 0
do ib = 1, nblocks
  g => gr(ib)
  tn = tn + g%nv
end do
write(*,*)"tn = ",tn
if (.not. allocated(iblanknodes)) allocate(iblanknodes(tn))

offset2 = 0
do ib = 1, nblocks
  g => gr(ib)
  do i = 1, g%nv
    iblanknodes(offset2 + i) = g%iblank(i)
  end do
  offset2 = offset2 + g%nv
end do

  call mpi_barrier(mpi_comm_world,ierr)
  if (myid==0) write(6,*) 'connectivity time=',t2-t1
 
!  call cpu_time(t1)
  do ib=1,nblocks
  g=>gr(ib)
  m=1
  do i=1,g%nv
     xt(:)=g%x(3*i-2:3*i)
          do n=1,g%nvar
       g%q(m)=(xt(1)+xt(2)+xt(3))
!       write(98,*) g%q(m)
       m=m+1
     enddo
  enddo
  enddo
!close(98)
  if (1==0) then
  do ib=1,nblocks
   call tioga_getdonorcount(ib,dcount,fcount)
   allocate(receptorInfo(4*dcount))
   allocate(inode(fcount),frac(fcount))
   !write(6,*) dcount,fcount
   !
   ! use this API if you need to interpolate the field variables yourself
   !
   call tioga_getdonorinfo(gr(ib)%bodytag(1),receptorInfo,inode,frac,dcount)
   !
   !> dcount = number of donors
   !> receptorInfo = {receptorProcess Id, receptor Index, receptor Block id, number of fractions}*dcount
   !> inode = indices for each receptor one group after the other
   !> frac  = weights for each receptor one group after the other
   !>
   m=0
   do i=1,dcount
    write(1000+myid,"(20(1x,I6))") ib,(inode(m+j),j=1,receptorInfo(4*i)+1),receptorInfo(4*i-3:4*i)
    m=m+(receptorInfo(4*i)+1)
   enddo 
   call flush(1000+myid)
   deallocate(receptorInfo,inode,frac)
  end do
  endif
!--------------------Get donor map-----------------------------



 !---------- Writing Iblank Cell ---------------(mod5)  
nc_t = 0                             !Initialising total no. cells = 0, looping meshblocks to get total no.
do ib=1,nblocks
   g=>gr(ib)
   if (g%n6 > 0) then
      nc_t = nc_t + g%n6
   else if (g%n8 > 0) then
      nc_t = nc_t + g%n8
   else if (g%n4 > 0) then
      nc_t = nc_t + g%n4
   endif
enddo
allocate(iblankcells(nc_t))      

m = 1                                  
do ib=1,nblocks
   g=>gr(ib)
   if (allocated(ibc)) deallocate(ibc)
   if (g%n6 > 0) then
      allocate(ibc(g%n6))
      nc = g%n6
   else if (g%n8 > 0) then
      allocate(ibc(g%n8))
      nc = g%n8
   else if (g%n4 > 0) then
      allocate(ibc(g%n4))
      nc = g%n4
   endif
   call tioga_getiblank_cell(g%bodytag(1), ibc, nc)

   do i = 1,nc              
      iblankcells(m) = ibc(i)   
      m = m + 1
   enddo

   deallocate(ibc)
enddo

!--------------------- writing donor info ---------------------------------------
if(is_CellBased==0)then
open(unit=23, file='checking_nodes/donorinfo.dat', status='replace', &
     action='write', iostat=ios1)
if (ios1 /= 0) then
write(*,*) 'ERROR opening donorinfo.dat'
endif
do blockid=1,nblocks
    call tioga_getdonorcount(blockid, dcount, fcount)
    allocate(receptorInfo(4*dcount))
    allocate(inode(fcount))
    allocate(frac(fcount))
    write(6,*)"dcount:",dcount,"fcount:",fcount
    call tioga_getdonorinfo(blockid,receptorInfo,inode,frac,dcount)
    write(23,*)"--------------------------------------------------"
    write(23,*)"Block Num=",blockid,"Donor Count = ",dcount,"Fcount=",fcount
    write(23,*)"--------------------------------------------------"
!    do d = 1, dcount
!    rbase = 4*(d-1)
!    ibase = 5*(d-1)
!     write(23,'(4(I8,1X),2X,5(I8,1X),2X,5(E16.8,1X))') &
!        receptorInfo(rbase+1), receptorInfo(rbase+2), &
!        receptorInfo(rbase+3), receptorInfo(rbase+4), &
!        inode(ibase+1), inode(ibase+2), inode(ibase+3), &
!        inode(ibase+4), inode(ibase+5), &
!        frac(ibase+1), frac(ibase+2), frac(ibase+3), &
!        frac(ibase+4), frac(ibase+5)
!    enddo
    do d = 1, dcount                            !This loop also calculates xyz of donor nodes into output (mod10h)
      rbase = 4*(d-1)
      ibase = 5*(d-1)
      do jj = 1, 4
        nd = inode(ibase+jj) + 1

        if (gr(blockid)%iblank(nd) /= 1) then      !Checking if all donors are field (mod10i)
          ndon_nf = ndon_nf + 1
          write(*,*) 'non-field donor: block', blockid, 'node', nd-1, 'iblank', gr(blockid)%iblank(nd)
          endif                                    !End mod10i

        sd(jj) = gr(blockid)%x(3*nd-2) + gr(blockid)%x(3*nd-1) + gr(blockid)%x(3*nd)
      enddo
      write(23,'(4(I8,1X),2X,5(I8,1X),2X,5(F16.8,1X),2X,4(F16.8,1X))') &
         receptorInfo(rbase+1:rbase+4), inode(ibase+1:ibase+5), &
         frac(ibase+1:ibase+5), sd
    enddo
    deallocate(receptorInfo,inode,frac)
enddo
write(*,*) 'non-field donor nodes found:', ndon_nf
close(23)
endif
!--------------------------------------------------------------------------------

if (.not. allocated(donor_dcount)) allocate(donor_dcount(nblocks), donor_roff(nblocks+1), donor_foff(nblocks+1), node_off(nblocks+1))
node_off(1) = 0
donor_roff(1) = 0
donor_foff(1) = 0
do ib=1,nblocks
   g=>gr(ib)
   node_off(ib+1) = node_off(ib) + g%nv
   call tioga_getdonorcount(g%bodytag(1),dcount,fcount)   ! fixed: bodytag, not ib
   donor_dcount(ib)  = dcount
   donor_roff(ib+1)  = donor_roff(ib) + 4*dcount
   donor_foff(ib+1)  = donor_foff(ib) + fcount
enddo
allocate(donor_receptorInfo(donor_roff(nblocks+1)))
allocate(donor_inode(donor_foff(nblocks+1)))
allocate(donor_frac(donor_foff(nblocks+1)))
do ib=1,nblocks
   g=>gr(ib)
   call tioga_getdonorinfo(g%bodytag(1), &
        donor_receptorInfo(donor_roff(ib)+1:donor_roff(ib+1)), &
        donor_inode(donor_foff(ib)+1:donor_foff(ib+1)), &
        donor_frac(donor_foff(ib)+1:donor_foff(ib+1)), &
        donor_dcount(ib))
enddo




end subroutine tioga_init_conn

!----------------------------------------------------------------------------------------------------------

subroutine tioga_solutions(q_node, neq, nodes)

real(8), allocatable :: uq_pre(:)                         !Var. declarations for pre and post interp check (mod10e)
integer :: uk_ok, uk_fr_unch, uk_nf_chg, uk_l(3), uk_g(3)
logical :: uk_chg
integer :: lt_ib, lt_i, lt_gn                             !Linear test Vars. 1 (mod10h)
real(8) :: lt_s                                           !Linear test Vars. 2 (mod10h)
character(len=64) :: uk_fn                                ! End var dec for mod10e

  real(8), intent(inout) :: q_node(:)
  integer, intent(in) :: neq,nodes
  integer :: k,mnode,m,i  
integer :: offset, coff           !Offset counters for loops (mod10)
character(len=64) :: finame
t=0

  do lt_ib = 1, nblocks                                          !Linear test q setting (mod10h)
    do lt_i = 1, gr(lt_ib)%nv
      lt_gn = node_off(lt_ib) + lt_i
      if (iblanknodes(lt_gn) /= 1) cycle
      lt_s = gr(lt_ib)%x(3*lt_i-2) + gr(lt_ib)%x(3*lt_i-1) + gr(lt_ib)%x(3*lt_i)
      do k = 1, neq
        q_node((lt_gn-1)*neq+k) = lt_s + dble(k)
      end do
    end do
  end do                                                         ! Linear test q set end (mod10h)

open(unit=88, file='checking_nodes/qvals_preinterp.dat', status='replace')  !Qvals pre interp (mod10g)
do i = 1, nodes
   WRITE(88,*) i, (q_node((i-1)*neq+k), k=1,neq), merge(1,0,iblanknodes(i)==-1)   !Adding flag for receptors (mod10h)
enddo
close(88)                                                                   !End Qvals pre interp

offset = 0
do ib=1,nblocks
   g=>gr(ib)
   call tioga_registersolution(g%bodytag(1), q_node(offset*neq+1:))
   offset = offset + g%nv
enddo
allocate(uq_pre(nodes*neq)); uq_pre = q_node(1:nodes*neq)               !mod10e
call cpu_time(t1)
call tioga_dataupdate_mb(neq, 'row')

open(unit=89, file='checking_nodes/qvals_postinterp.dat', status='replace')  !Qvals post interp (mod10g)
do i = 1, nodes
   WRITE(89,*) i, (q_node((i-1)*neq+k), k=1,neq), merge(1,0,iblanknodes(i)==-1)   !Adding flag for receptors (mod10h)
enddo
close(89)                                                                   !End Qvals post interp

uk_ok = 0; uk_fr_unch = 0; uk_nf_chg = 0                                 !Actual pre and post interp check (mod10e)
write(uk_fn,'(A,I4.4,A)') 'checking_nodes/update_check_', myid, '.dat'
open(unit=96, file=trim(uk_fn), status='replace')
write(96,'(A)') '# node iblank issue'
do i = 1, nodes
  uk_chg = any(q_node((i-1)*neq+1:i*neq) /= uq_pre((i-1)*neq+1:i*neq))
  if (iblanknodes(i) == -1) then
    if (uk_chg) then
      uk_ok = uk_ok + 1
    else
      uk_fr_unch = uk_fr_unch + 1
      write(96,*) i, iblanknodes(i), 'fringe_unchanged'
    end if
  else if (uk_chg) then
    uk_nf_chg = uk_nf_chg + 1
    write(96,*) i, iblanknodes(i), 'nonfringe_changed'
  end if
end do
close(96)
deallocate(uq_pre)
uk_l = (/uk_ok, uk_fr_unch, uk_nf_chg/)
call mpi_allreduce(uk_l, uk_g, 3, mpi_integer, mpi_sum, mpi_comm_world, ierr)
if (myid == 0) print*, 'fringe updated / fringe unchanged / nonfringe changed:', uk_g !END mod10e

!  do ib=1,nblocks
!    g=>gr(ib)
!    call tioga_registersolution(g%bodytag(1),q_node)
!  enddo
!  call tioga_dataupdate_mb(neq,'row')    !< update the q-variables (can be called anywhere)
                                             !< nvar = number of field variables per node
                                             !< if fields are different arrays, you can also 
                                             !< call this multiple times for each field
  call mpi_barrier(mpi_comm_world,ierr)
    call cpu_time(t2)
  if (myid==0) write(6,*) 'data update time=',t2-t1


if (.not. allocated(cellvals)) allocate(cellvals(nc_t,neq))
cellvals = 0.0d0
coff = 0
do ib = 1, nblocks
  g => gr(ib)
  do i = 1, g%ncells
    if (iblankcells(coff+i) /= -1) cycle
    do k = 1, neq
      do m = 1, 4
        mnode = node_off(ib) + g%ndc4(m,i)
        cellvals(coff+i,k) = cellvals(coff+i,k) + q_node((mnode-1)*neq+k)
      enddo
      cellvals(coff+i,k) = cellvals(coff+i,k) / 4.0d0
    enddo
  enddo
  coff = coff + g%ncells
enddo

  !
  ! compute error in interpolation
  ! should be machine-zero for the linear
  ! problem here
  ! 
!  do ib=1,nblocks
!   g=>gr(ib)
!  m=1
!   do i=1,g%nv
!    xt(:)=g%x(3*i-2:3*i)
!     do n=1,g%nvar
!       rnorm=rnorm+(g%q(m)-(xt(1)+xt(2)+xt(3)))**2
!       m=m+1
!     enddo
!   enddo
!  enddo
!  if (myid==0) write(6,"(A36)") '-- Interpolation error statistics --'
!  if (myid==0) write(6,"(A15)") 'ProcId    Error'
!  call flush(6)
!  call mpi_barrier(mpi_comm_world,ierr)     ! global synchronisation point (all process must reach here b4 cont.)
!  write(6,"(I4,3x,E15.7)") myid,sqrt(rnorm/g%nv/g%nvar)
  call mpi_barrier(mpi_comm_world,ierr)
!  call tioga_writeoutputfiles(g%nvar,'row',t) !< write output files, if need be
!  write(*,*) 'Moving block:', g%bodytag(1)
!  call movegrid(g)
end subroutine tioga_solutions

!---------------------------------------------------------------------------------------------------------------

subroutine tioga_fin()
  call mpi_barrier(mpi_comm_world,ierr)
  call tioga_delete
!  call mpi_finalize(ierr)
end subroutine tioga_fin
  !
end module run_tioga
