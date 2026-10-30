program calcqvals
implicit none
integer :: p, ios, blockno, dcount, fcount, i, junk, rec1, r_blk, dcell,dcount_tot,x2
integer :: m,ntot,etot,c,rmax,rmin,n1,n2
real(8) :: w1,w2,w3,w4,junkr,s1,s2,s3,s4
integer :: nn,rflag
real(8) :: q1,q2,q3,q4,q5,q6,q7,q8,q9,q10
real(8) :: vals1(10), vals2(10)
character(len=256) :: line
real(8),allocatable :: replacement(:,:)
logical,allocatable :: has_replacement(:)

open(unit=1,file='qvals_preinterp.dat',status='old',action='read')
open(unit=2,file='donorinfo.dat',status='old',action='read')
open(unit=5,file='qvals_postintcheck.dat',status='replace')
open(unit=7,file='receptorsolns.dat',status='replace')

ntot = 146280
etot = 805194
c=0
rmax = -1
rmin = 999999
!------------------------Read section 1 in donorinfo ---------------------
read(2,*)
read(2,'(A)')line

p = index(line, "Block Num=")
read(line(p+len('Block Num='):),*,iostat=ios) blockno
if (ios/=0) stop "Error Reading Block Number"

p = index(line, "Donor Count =")
read(line(p+len('Donor Count ='):),*,iostat=ios) dcount
if (ios/=0) stop "Error Reading Donor Count"

p = index(line, "Fcount=")
read(line(p+len('Fcount='):),*,iostat=ios) fcount
if (ios/=0) stop "Error Reading F Count"

write(*,*)"Block  = ",blockno
write(*,*)"Donors = ",dcount
write(*,*)"Fcount = ",fcount
read(2,*)

do i=1,dcount
  read(2,*)junk,rec1,r_blk,junk,junk,junk,junk,junk,dcell,w1,w2,w3,w4,junkr,s1,s2,s3,s4
  if(r_blk ==1) then
     rec1 = rec1+48000
  endif
  q1=((s1+1)*w1)+((s2+1)*w2)+((s3+1)*w3)+((s4+1)*w4)
  q2=((s1+2)*w1)+((s2+2)*w2)+((s3+2)*w3)+((s4+2)*w4)
  q3=((s1+3)*w1)+((s2+3)*w2)+((s3+3)*w3)+((s4+3)*w4)
  q4=((s1+4)*w1)+((s2+4)*w2)+((s3+4)*w3)+((s4+4)*w4)
  q5=((s1+5)*w1)+((s2+5)*w2)+((s3+5)*w3)+((s4+5)*w4)
  q6=((s1+6)*w1)+((s2+6)*w2)+((s3+6)*w3)+((s4+6)*w4)
  q7=((s1+7)*w1)+((s2+7)*w2)+((s3+7)*w3)+((s4+7)*w4)
  q8=((s1+8)*w1)+((s2+8)*w2)+((s3+8)*w3)+((s4+8)*w4)
  q9=((s1+9)*w1)+((s2+9)*w2)+((s3+9)*w3)+((s4+9)*w4)
  q10=((s1+10)*w1)+((s2+10)*w2)+((s3+10)*w3)+((s4+10)*w4)
  write(7,*)(rec1+1),q1,q2,q3,q4,q5,q6,q7,q8,q9,q10
  if((rec1+1)>rmax)rmax=(rec1+1)
  if((rec1+1)<rmin)rmin=(rec1+1)
enddo 

dcount_tot=0
dcount_tot = dcount_tot+dcount
write(*,*)"Total Dcount = ",dcount_tot

!----------------Read second section------------------------
read(2,*)
read(2,'(A)')line

p = index(line, "Block Num=")
read(line(p+len('Block Num='):),*,iostat=ios) blockno
if (ios/=0) stop "Error Reading Block Number"

p = index(line, "Donor Count =")
read(line(p+len('Donor Count ='):),*,iostat=ios) dcount
if (ios/=0) stop "Error Reading Donor Count"

p = index(line, "Fcount=")
read(line(p+len('Fcount='):),*,iostat=ios) fcount
if (ios/=0) stop "Error Reading F Count"

write(*,*)"Block  = ",blockno
write(*,*)"Donors = ",dcount
write(*,*)"Fcount = ",fcount
read(2,*)

do i=1,dcount
  read(2,*)junk,rec1,r_blk,junk,junk,junk,junk,junk,dcell,w1,w2,w3,w4,junkr,s1,s2,s3,s4
  if(r_blk ==1) then
     rec1 = rec1+48000
  endif
  q1=((s1+1)*w1)+((s2+1)*w2)+((s3+1)*w3)+((s4+1)*w4)
  q2=((s1+2)*w1)+((s2+2)*w2)+((s3+2)*w3)+((s4+2)*w4)
  q3=((s1+3)*w1)+((s2+3)*w2)+((s3+3)*w3)+((s4+3)*w4)
  q4=((s1+4)*w1)+((s2+4)*w2)+((s3+4)*w3)+((s4+4)*w4)
  q5=((s1+5)*w1)+((s2+5)*w2)+((s3+5)*w3)+((s4+5)*w4)
  q6=((s1+6)*w1)+((s2+6)*w2)+((s3+6)*w3)+((s4+6)*w4)
  q7=((s1+7)*w1)+((s2+7)*w2)+((s3+7)*w3)+((s4+7)*w4)
  q8=((s1+8)*w1)+((s2+8)*w2)+((s3+8)*w3)+((s4+8)*w4)
  q9=((s1+9)*w1)+((s2+9)*w2)+((s3+9)*w3)+((s4+9)*w4)
  q10=((s1+10)*w1)+((s2+10)*w2)+((s3+10)*w3)+((s4+10)*w4)
  write(7,*)(rec1+1),q1,q2,q3,q4,q5,q6,q7,q8,q9,q10
  if((rec1+1)>rmax)rmax=(rec1+1)
  if((rec1+1)<rmin)rmin=(rec1+1)
enddo

close(7)

dcount_tot = dcount_tot+dcount
write(*,*)"Total Dcount = ",dcount_tot
write(*,*)"Max rec node = ",rmax
write(*,*)"Min rec node = ",rmin
allocate(replacement(10,ntot),has_replacement(ntot))
has_replacement = .false.
!-----------------------Reading qvals post interp------------------------------------
open(unit=7,file='receptorsolns.dat',status='old')
do
  read(7,*,iostat=ios)n2,vals2
  if(ios/=0)exit
  replacement(:,n2)=vals2
  has_replacement(n2)=.true.
enddo

c=0
x2=0
do
  read(1,*,iostat=ios)nn,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,rflag
  if(ios/=0)exit
  if(rflag/=1)then
    write(5,*)nn,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,rflag
  elseif(rflag ==1) then
    write(5,*) nn,replacement(:,nn)
    c=c+1
  else
    x2=x2+1
  endif
enddo
write(*,*)"Replaced Values =",c
if(x2/=0)write(*,*)"UHOHHHHH....."
close(1);close(2);close(5);close(7)
end program
