FILES NEEDED :

1. donorinfo.dat
2.qvals_preinterp.dat
3.qvals_postinterp.dat

CHECK FLOW:

1. calcqvals2.f90 - obj: cqv2 - calculates what qvals should be post interp - outputs receptorsoln and qvals_postintcheck.dat

2. checking.f90 - obj: check - checks the calculated qvals post interp with the qvals pre interp 
