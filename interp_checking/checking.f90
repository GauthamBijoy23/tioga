program checkingqvals
implicit none

integer :: ios,rflag,nn1,nn2,c,d
real(8) :: vals1(10),vals2(10)

open(unit=1,file='qvals_postinterp.dat',status='old')
open(unit=2,file='qvals_postintcheck.dat',status='old')

c=0
d=0
do
  read(1,*,iostat=ios)nn1,vals1,rflag
  if(ios/=0)exit
  read(2,*)nn2,vals2
  if(rflag/=0)then
    d=d+1
    if(any(abs(vals1-vals2) .gt. 1.0d-6))then
      c=c+1
    endif
  endif
enddo

write(*,*)"Loop executed",d,"times"

if(c/=0)then
  write(*,*)"Diffs exist"
else
  write(*,*)"ALL GOOD"
endif

close(1);close(2)
end program
