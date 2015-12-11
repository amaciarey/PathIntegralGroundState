module global_mod

implicit none

logical          :: wf_table,v_table
real (kind=8)    :: pi
real (kind=8)    :: V0
real (kind=8)    :: Rm
real (kind=8)    :: rbin,dr
real (kind=8)    :: rcut,rcut2
real (kind=8)    :: CWorm
integer (kind=4) :: dim,Np,Nbin,Nb,Nmax,Npw

real (kind=8),dimension (:),allocatable :: Lbox,LboxHalf,qbin

contains

!-----------------------------------------------------------------------
  
  function GreenFunction(opt,ib,dt,Pot,F2)

    implicit none

    real (kind=8)    :: GreenFunction
    real (kind=8)    :: dt,Pot,F2
    real (kind=8)    :: Ve,dVe
    real (kind=8)    :: Vc,dVc
    integer (kind=4) :: ib,opt

    GreenFunction = 0.d0

    if (opt==0) then

       Ve = Pot
       Vc = Pot+dt**2*F2/6.d0

       if (ib==0) then
          GreenFunction = dt*Ve/3.d0
       else if (ib==2*Nb) then
          GreenFunction = dt*Ve/3.d0
       else
          if (mod(ib,2)==0) then
             GreenFunction = 2.d0*dt*Ve/3.d0
          else
             GreenFunction = 4.d0*dt*Vc/3.d0
          end if
       end if

    else if (opt==1) then
       
       dVe = Pot
       dVc = Pot+dt**2*F2/2.d0

       if (ib==0) then
          GreenFunction = dVe/3.d0
       else if (ib==2*Nb) then
          GreenFunction = dVe/3.d0
       else
          if (mod(ib,2)==0) then
             GreenFunction = 2.d0*dVe/3.d0
          else
             GreenFunction = 4.d0*dVc/3.d0
          end if
       end if

    end if
       
    return
  end function GreenFunction

!-----------------------------------------------------------------------

end module global_mod
