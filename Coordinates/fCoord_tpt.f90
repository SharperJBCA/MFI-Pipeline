!! COORDINATES_TPOINT.F90
!
! COORDINATES_TPOINT: A SET OF ROUTINES FOR CONVERTING BETWEEN DIFFERENT COORDINATE SYSTEMS.
!  EACH FUNCTION INCLUDES TPOINT VARIABLE INPUTS FOR AZEL TELESCOPE MOUNTS AND SKY ANGLES FOR NONE
!  FOCAL PLANE GEOMETRY.
!

subroutine get_e2h(ra, dec, mjd, lng, lat, az,el, n)
  implicit none

  integer n
!f2py intent(in) n
  real*8 lat, lng
!f2py intent(in) lat, lng
  real*8 ra(n), dec(n), mjd(n)
!f2py intent(in) ra, dec, mjd
  real*8 az(n), el(n)
!f2py intent(out) az, el

  real*8, external :: sla_gmst
  real*8  ha,lst
  integer  i

  do i=1,n
     lst = sla_gmst(mjd(i))
     !lst = lst - lng
     ha  = lst - ra(i) - lng
     call sla_de2h(ha,dec(i),lat,az(i),el(i))
  end do

end subroutine get_e2h

subroutine get_h2e_tpoint(az, el, xpos, ypos, Pf, Px, Py, Pc, Pn, Pa, Pb, lat,lng, mjd,gal, ra,dec,p, n)
  implicit none

  integer n, gal
!f2py intent(in) n, gal
  
  real*8 mjd(n)
  real*8 az(n), el(n), lat,lng
!f2py intent(in) az, el ,mjd, lat,lng

  real*8 xpos,ypos,Pf, Px, Py, Pc, Pn, Pa, Pb
!f2py intent(in) xpos,ypos,Pf, Px, Py, Pc, Pn, Pa, Pb

  real*8 ra(n), dec(n), p(n)
!f2py intent(out) ra, dec, p

  !Internal arrays:
  integer i

  real*8 azc,elc,l,b
  real*8 ep0,ha,lst
  real*8, external :: sla_gmst

  !Apply corrections:
  do i=1,n
     call ffunc_h2e(az(i),el(i),xpos,ypos,Pf,Px,Py,Pc,Pn,Pa,Pb,azc,elc)


     !Convert to RA/DEC
     call sla_dh2e(azc,elc,lat,ha,dec(i))
     lst   = sla_gmst(mjd(i))
     lst   = lst - lng
     ra(i) = lst - ha
     
     !Precess:
     ep0 = 2000 + (mjd(i)-51545D0)/365.25D0
     call sla_preces('FK5',ep0,2000D0, ra(i),dec(i))
     
     !Parallactic Angle:
     p(i) = atan2( sin(ha), cos(dec(i))*tan(lat) - sin(dec(i))*cos(ha) )
     
     if (gal == 1) then
        call sla_eqgal(ra(i),dec(i),l,b)
        ra(i) = l
        dec(i) = b
     end if

  end do 
  
end subroutine get_h2e_tpoint

subroutine get_horizon_tpoint(az, el, xpos, ypos, Pf, Px, Py, Pc, Pn, Pa, Pb,iaz,iel, n)
  implicit none

  integer n
!f2py intent(in) n
  
  real*8 az(n), el(n)
!f2py intent(in) az, el

  real*8 xpos,ypos,Pf, Px, Py, Pc, Pn, Pa, Pb
!f2py intent(in) xpos,ypos,Pf, Px, Py, Pc, Pn, Pa, Pb

  real*8 iaz(n), iel(n)
!f2py intent(out) iaz,iel

  !Internal arrays:
  integer i

  real*8 azc,elc

  !Apply corrections:
  do i=1,n
     call ffunc_h2e(az(i),el(i),xpos,ypos,Pf,Px,Py,Pc,Pn,Pa,Pb,iaz(i),iel(i))
  end do 
  
end subroutine get_horizon_tpoint


!!!! DENIS' POINTING MODEL !!!!

subroutine ffunc_h2e(az, el, xpos, ypos, Pf, Px, Py, Pc, Pn, Pa, Pb,azc,elc)
  implicit none

  real*8 az,el,xpos,ypos,Pf,Px,Py,Pc,Pn,Pa,Pb

  real*8 top,bot,daz,del
  real*8 azc,elc
  real*8 x,y,z,w
  real*8 x0,y0,z0
  real*8 xa,ya,za

  real*8 dx0,dy0,dz0,dxn,dyn,dzn

  integer j
  ! Declare local constant Pi
  real*8, PARAMETER :: Pi = 3.141592653589793238462643383279502884197169399
  integer, PARAMETER :: maxiter = 100

  !Horn focal plane position:
  IF (abs(ypos) > 0) THEN
     top = xpos*sin(el)*cos(el) * (1.0/tan(el) + tan(el))
     bot = cos(el) - ypos*sin(el)
     daz = atan2(top,bot)
     azc = az + daz 
  ELSE
     azc = az
  ENDIF

  IF (abs(xpos) > 0) THEN
     top = sin(daz)/xpos - (cos(el)*cos(daz))
     bot = sin(el)
     del = atan(top/bot)
     elc = del
  ELSE
     elc = el
  ENDIF


  IF ((abs(xpos) > 0) .AND. (abs(ypos) > 0)) THEN 

     !Encoder offsets:
     azc = azc - Pa
     elc = elc + Pb

     !Convert to cartisian
     x = -cos(elc) * cos(azc)
     y =  cos(elc) * sin(azc)
     z =  sin(elc)
     
     !Non-perpendicularities:
     x0 = x
     y0 = y
     z0 = z
  
     do j=1,maxiter
        w = (Pc + Pn*z) / sqrt(x*x + y*y)
        dx0 = - w*(y - w*x)
        dy0 =   w*(x + w*y)
        dz0 =   0.0
        
        xa = x0 + dx0
        ya = y0 + dy0
        za = z0 + dz0
        
        w   =     (Pc+ Pn*za) / sqrt(xa*xa + ya*ya)
        dxn = - w*(ya - w*xa)
        dyn =   w*(xa + w*ya)
        dzn =   0.0
        
        x = xa
        y = ya
        z = za
        
        IF (abs(dx0 - dxn) < 0.00001 .AND. abs(dy0 - dyn) < 0.00001 .AND. abs(dz0 - dzn) < 0.00001) THEN
           EXIT
        END IF
     end do
  
     x = x0 + dxn
     y = y0 + dyn
     z = z0 + dzn
  
     !Roll-axis misalignment:
     x0 = x
     y0 = y
     z0 = z
     
     do j=1,maxiter
        dx0 = - Px*z
        dy0 = - Py*z
        dz0 =   Px*x + Py*y
        
        xa = x0 + dx0
        ya = y0 + dy0
        za = z0 + dz0
        
        dxn = - Px*za
        dyn = - Py*za
        dzn =   Px*xa + Py*ya
        
        x = xa
        y = ya
        z = za
     
        IF (abs(dx0 - dxn) < 0.00001 .AND. abs(dy0 - dyn) < 0.00001 .AND. abs(dz0 - dzn) < 0.00001) THEN
           EXIT
        END IF
     end do
  
     x = x0 + dxn
     y = y0 + dyn
     z = z0 + dzn
     
     !Hookes Law:

     x0 = x
     y0 = y
     z0 = z
  
     do j=1,maxiter
        dx0 =   Pf*z*x
        dy0 =   Pf*z*y
        dz0 = - Pf*(x*x + y*y)
        
        xa = x0 + dx0
        ya = y0 + dy0
        za = z0 + dz0
     
        dxn =   Pf*za*xa
        dyn =   Pf*za*ya 
        dzn = - Pf*(xa*xa + ya*ya) 
     
        x = xa
        y = ya
        z = za
     
        IF (abs(dx0 - dxn) < 0.00001 .AND. abs(dy0 - dyn) < 0.00001 .AND. abs(dz0 - dzn) < 0.00001) THEN
           EXIT
        END IF
     end do
  
     x = x0 + dxn
     y = y0 + dyn
     z = z0 + dzn
  
     !Convert back to AZ,EL
     elc = asin( z/sqrt(x*x + y*y + z*z) ) 
     azc = pi - atan2(y,x)



  ELSE
     azc = az
     elc = el
  END IF 


end subroutine ffunc_h2e
