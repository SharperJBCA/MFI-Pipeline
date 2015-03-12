!!Quick routine to wrap the SLALIB planet position function.
!!It is inaccurate but fine for QUIJOTE accuracy. It downsamples
!! the data by 100x!

subroutine plpos(mjd,lon,lat,planet,gal,ra,dec)
  implicit none
  
  integer  planet, gal
!f2py intent(in) planet, gal

  real*8 mjd, lon, lat
!f2py intent(in) mjd, lon, lat

  real*8 ra, dec
!f2py intent(out) ra, dec

  real*8 diam, l, b , pi
  
  pi = 4. *atan(1.)

  call sla_rdplan(mjd,planet,lon,lat,ra,dec,diam)

  if (gal .eq. 1) then
     call sla_eqgal(ra,dec,l,b)
     ra = l
     dec = b
  end if
  
end subroutine plpos



subroutine planet(mjd,NP,PV)
  implicit none
  
  integer  NP
!f2py intent(in) NP

  real*8 mjd
!f2py intent(in) mjd

  real*8 PV(6)
!f2py intent(out) PV

  integer JSTAT

  call sla_planet(mjd,NP,PV,JSTAT)
  
end subroutine planet
