module system_mod

use global_mod, only: rcut,dim,V0,Nmax,dr
use bessel_mod, only: Bessk

implicit none

contains

!-----------------------------------------------------------------------

  function LogPsi(opt,Rm,rij)

    logical,save     :: FirstCall = .True.
    real (kind=8)    :: Rm
    real (kind=8)    :: rij,zij
    real (kind=8)    :: K0,K1
    real (kind=8)    :: dfdr,d2fdr2
    real (kind=8)    :: LogPsi
    integer (kind=4) :: opt
    
    real (kind=8),save :: Lmax,C1,C2,C3
    
    if (FirstCall) then

       Lmax = 2.d0*rcut

       C3 = 0.5d0*sqrt(4.d0*Rm)*(Lmax-Rm)**2/(Lmax*(Lmax-2.d0*Rm))*&
            & Bessk(1,sqrt(4.d0/Rm))/Bessk(0,sqrt(4.d0/Rm))
       C2 = exp(4.d0*C3/Lmax)
       C1 = C2*exp(-C3*(1.d0/Rm+1.d0/(Lmax-Rm)))/Bessk(0,sqrt(4.d0/Rm))
       
       FirstCall = .False.
       
    end if

    zij = sqrt(4.d0/rij)
    
    K0 = Bessk(0,zij)
    K1 = Bessk(1,zij)

    if (opt==0) then

       if (rij<=Rm) then
      
          LogPsi = log(C1*K0)
    
       else
   
          LogPsi = log(C2)-C3*(1.d0/rij+1.d0/(Lmax-rij))
   
       end if

    else if (opt==1) then
       
       if (rij<=Rm) then
      
          dfdr   = zij*K1/(2.d0*rij) 
          LogPsi = dfdr/K0

       else

          LogPsi = C3*(1.d0/rij**2-1.d0/(Lmax-rij)**2)
      
       end if
 
    else if (opt==2) then

       if (rij<=Rm) then

          dfdr   = zij*K1/(2.d0*rij)
          d2fdr2 = K0/rij**3-dfdr/rij
          LogPsi = (d2fdr2*K0-dfdr*dfdr)/(K0*K0)
     
       else

          LogPsi = -2.d0*C3*(1.d0/rij**3+1.d0/(Lmax-rij)**3)

       end if

    else 

       print *, 'The parameter opt in the function LogPsi crash!!!'
       stop

    end if
    
    return
  end function LogPsi

!-----------------------------------------------------------------------

  function Potential(xij,rij)

    implicit none
    
    real (kind=8) :: Potential,rij
    
    real (kind=8),dimension (dim) :: xij

    Potential = (1.d0-3.d0*V0*(xij(1)/rij)**2)/rij**3
    
  end function Potential

!-----------------------------------------------------------------------  

  function Force(k,xij,rij)

    implicit none 
    
    real (kind=8)    :: rij,Force
    real (kind=8)    :: costheta
    real (kind=8)    :: dVdr
    integer (kind=4) :: k

    real (kind=8),dimension (dim) :: xij
    
    costheta = xij(1)/rij
    dVdr     = -3.d0/rij**4
    Force    = dVdr*(1.d0-5.d0*V0*costheta**2)*xij(k)/rij

    if (k==1) Force = Force+2.d0*V0*dVdr*xij(k)/rij 

!!$    rij_i = 1.d0/rij
!!$
!!$    costheta = xij(1)*rij_i
!!$    dVdr     = -3.d0*rij_i**4
!!$    Force    = dVdr*(1.d0-5.d0*V0*costheta**2)*xij(k)*rij_i
!!$    
!!$    if (k==1) Force = Force+2.d0*V0*dVdr*xij(k)*rij_i
    
    return 
  end function Force

!-----------------------------------------------------------------------

end module system_mod
