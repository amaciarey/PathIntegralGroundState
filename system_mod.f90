module system_mod

use global_mod, only: rcut,dim,V0,Nmax,dr
use bessel_mod, only: Bessk

implicit none

contains

!-----------------------------------------------------------------------

  function LogPsi(opt,Rm,rij)

    real (kind=8)    :: Rm,rij
    real (kind=8)    :: LogPsi
    integer (kind=4) :: opt
    
    if (opt==0) then
       
       LogPsi = -0.5d0*(Rm/rij)**5
            
    else if (opt==1) then
       
       LogPsi = 2.5d0*(Rm/rij)**5/rij
       
    else if (opt==2) then

       LogPsi = -15.d0*(Rm/rij)**5/rij**2

    else 

       print *, 'The parameter opt in the function LogPsi crash!!!'
       stop
       
    end if

    return
  end function LogPsi

!-----------------------------------------------------------------------

  function Potential(xij,rij)

    implicit none

    logical, save :: FirstCall = .true.
    real (kind=8) :: Potential,Hx,rij
    real (kind=8) :: dij,dij2,dij4,dij6

    real (kind=8), save            :: A,alpha,C6,C8,C10,D,V0
    real (kind=8), dimension (dim) :: xij

    if (FirstCall) then

       !Parameters of the Aziz potential

       V0    = 5.81817d0
       A     = 0.54485046d6
       alpha = 13.353384d0
       C6    = 1.3732412d0
       C8    = 0.4253785d0
       C10   = 0.1781d0
       D     = 1.241314d0

       FirstCall = .false.

    end if

    dij  = rij*2.556d0/2.9673d0
    dij2 = dij*dij
    dij4 = dij2*dij2
    dij6 = dij4*dij2

    if (dij<=D) then
       Hx = exp(-(D/dij-1.d0)**2)
    else
       Hx = 1.d0
    end if

    Potential = V0*(A*exp(-alpha*dij)-(C6+C8/dij2+C10/dij4)*Hx/dij6)

  end function Potential

!-----------------------------------------------------------------------

!!$  function Potential(xij,rij)
!!$
!!$    implicit none
!!$
!!$    real (kind=8) :: Potential,V0,rij,rij6
!!$
!!$    real (kind=8), dimension(dim) :: xij
!!$
!!$    V0   = 22.0228d0
!!$    rij6 = rij**6
!!$
!!$    Potential = V0*(1.d0/rij6-1.d0)/rij6
!!$
!!$  end function Potential

!-----------------------------------------------------------------------

   function Force(k,xij,rij)

    implicit none 
    
    real (kind=8)    :: rij,Force
    real (kind=8)    :: dVdr,V0
    integer (kind=4) :: k

    real (kind=8),dimension (dim) :: xij
    
    V0 = 22.0228d0
    
    dVdr  = -6.d0*V0*(2.d0/rij**7-1.d0)/rij**6
    Force = dVdr*xij(k)/rij

    return 
  end function Force

!-----------------------------------------------------------------------  

end module system_mod

