!Routine for fast downsampling of an array of doubles.


subroutine DownSample(a,newlen,alen,out,err)
  implicit none

  integer alen, newlen
!f2py intent(in) alen, newlen
  
  real*8 a(alen)
!f2py intent(in) a

  real*8 out(newlen),err(newlen)
!f2py intent(out) out, err


  integer nsteps,i, k , count

  real*8 suma,suma2,fnsteps,laststep

  nsteps = alen/newlen
  fnsteps = real(nsteps)

  suma = 0
  suma2= 0
  k = 1
  laststep = fnsteps + real(alen - newlen*nsteps)
  count = 0

  do i=1,alen

     if (i .eq. alen-laststep+1) then
        nsteps  = int(laststep)
        fnsteps = laststep
     end if

     !Calc mean and stderr of samples in step
     suma  = suma  + a(i)
     suma2 = suma2 + a(i)*a(i)
     
     if (mod(i-count,nsteps) .eq. 0) then             
        out(k) = suma /fnsteps
        err(k) = suma2/fnsteps - suma*suma/fnsteps/fnsteps
        !If last err only has one element
        ! Use error of previous box.
        if (err(k) .eq. 0 .and. k .gt. 1) then
           err(k) = err(k-1)
        end if
        
        !reset suma and suma2
        suma = 0
        suma2= 0
        
        !increase output index k
        k = k + 1
        count = i
     end if

        
  end do

end subroutine DownSample
  
