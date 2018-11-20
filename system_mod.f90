module system_mod

use global_mod, only: rcut,dim,Nmax,dr
use bessel_mod, only: Bessk

implicit none

real (kind=8) :: Rm
real (kind=8), dimension (:), allocatable :: a_ho

contains

!-----------------------------------------------------------------------

  subroutine ReadSystemParameters(trap)
    
    implicit none

    logical :: trap

    namelist /jastrow/ Rm
    namelist /extpot/ a_ho
  
    if (trap) then
       allocate (a_ho(dim))
       read (5,nml=extpot)
       rewind (5)
    end if
  
    read (5,nml=jastrow)
    rewind (5)

    return
  end subroutine ReadSystemParameters

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

    !LogPsi = 0.d0

    return
  end function LogPsi

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

!!$  function Potential(xij,rij)
!!$
!!$    !Aziz I HFDHE2 Potential
!!$
!!$    implicit none
!!$
!!$    logical, save :: FirstCall = .true.
!!$    real (kind=8) :: Potential,Hx,rij
!!$    real (kind=8) :: dij,dij2,dij4,dij6
!!$
!!$    real (kind=8), save            :: A,alpha,C6,C8,C10,D,E_0,V0,rm
!!$    real (kind=8), dimension (dim) :: xij
!!$
!!$    if (FirstCall) then
!!$
!!$       !Parameters of the Aziz potential
!!$
!!$       E_0   = 10.8d0
!!$       rm    = 2.9673d0
!!$       A     = 0.54485046d6
!!$       alpha = 13.353384d0
!!$       C6    = 1.3732412d0
!!$       C8    = 0.4253785d0
!!$       C10   = 0.1781d0
!!$       D     = 1.241314d0
!!$
!!$       V0    = E_0/1.85505153154686d0

!!$       FirstCall = .false.
!!$
!!$    end if
!!$
!!$    dij  = rij*2.556d0/rm
!!$    dij2 = dij*dij
!!$    dij4 = dij2*dij2
!!$    dij6 = dij4*dij2
!!$
!!$    if (dij<=D) then
!!$       Hx = exp(-(D/dij-1.d0)**2)
!!$    else
!!$       Hx = 1.d0
!!$    end if
!!$
!!$    Potential = V0*(A*exp(-alpha*dij)-(C6+C8/dij2+C10/dij4)*Hx/dij6)
!!$
!!$  end function Potential

!-----------------------------------------------------------------------

  function Potential(xij,rij)

    !Aziz II HFD-B(HE) Potential

    implicit none

    logical, save :: FirstCall = .true.
    real (kind=8) :: Potential,Hx,rij
    real (kind=8) :: dij,dij2,dij4,dij6

    real (kind=8), save            :: A,alpha,beta,C6,C8,C10,D,E_0,V0,rm
    real (kind=8), dimension (dim) :: xij

    if (FirstCall) then

       !Parameters of the Aziz potential

       E_0   = 10.948d0
       rm    = 2.963d0
       A     = 1.8443101d5
       alpha = 10.43329537d0
       beta  = -2.27965105d0
       C6    = 1.36745214d0
       C8    = 0.42123807d0
       C10   = 0.17473318d0
       D     = 1.4826d0

       V0    = E_0/1.85505153154686d0

       FirstCall = .false.

    end if

    dij  = rij*2.556d0/rm
    dij2 = dij*dij
    dij4 = dij2*dij2
    dij6 = dij4*dij2

    if (dij<=D) then
       Hx = exp(-(D/dij-1.d0)**2)
    else
       Hx = 1.d0
    end if

    Potential = V0*(A*exp(-alpha*dij+beta*dij2)-(C6+C8/dij2+C10/dij4)*Hx/dij6)

  end function Potential

!-----------------------------------------------------------------------  

  function Force(k,xij,rij)

    implicit none 

    real (kind=8)    :: rij,Force
    integer (kind=4) :: k
    
    real (kind=8), dimension (dim) :: xij
    
!!$    real (kind=8)    :: rij,Force
!!$    real (kind=8)    :: costheta
!!$    real (kind=8)    :: dVdr
!!$    integer (kind=4) :: k
!!$
!!$    real (kind=8),dimension (dim) :: xij
!!$    
!!$    costheta = xij(1)/rij
!!$    dVdr     = -3.d0/rij**4
!!$    Force    = dVdr*(1.d0-5.d0*V0*costheta**2)*xij(k)/rij
!!$
!!$    if (k==1) Force = Force+2.d0*V0*dVdr*xij(k)/rij 

    return 
  end function Force

!-----------------------------------------------------------------------

   function TrapPsi(opt,a_osc,x)

    implicit none

    real (kind=8)    :: TrapPsi,a_osc,x
    integer (kind=4) :: opt

    if (opt==0) then
       TrapPsi = -0.5d0*(x/a_osc)**2
    else if (opt==1) then
       TrapPsi = -(x/a_osc**2)
    else if (opt==2) then
       TrapPsi = -1.d0/a_osc**2
    else
       print *, 'The parameter opt in the function TrapPsi crash!!!'
       stop
    end if

    !TrapPsi = 0.d0

    return
  end function TrapPsi

!-----------------------------------------------------------------------

  function TrapPot(opt,a_osc,x)

    implicit none

    real (kind=8)    :: TrapPot,a_osc,x
    integer (kind=4) :: opt

    if (opt==0) then
       TrapPot = 0.5d0*x**2/a_osc**4
    else if (opt==1) then
       TrapPot = x/a_osc**4
    end if

    return
  end function TrapPot

!-----------------------------------------------------------------------  

end module system_mod
