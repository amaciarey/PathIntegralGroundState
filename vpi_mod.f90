module vpi_mod

use random_mod
use global_mod
use system_mod
use pbc_mod

implicit none

contains

!-----------------------------------------------------------------------

  function GreenFunction(opt,ib,dt,Pot,F2)

    implicit none

    real (kind=8)    :: GreenFunction
    real (kind=8)    :: dt,Pot,F2
    integer (kind=4) :: ib,opt

    GreenFunction = 0.d0

    if (opt==0) then

       if (mod(ib,2)==0) then
          GreenFunction = dt*2.d0*(Pot+a_1*dt*dt*F2/6.d0)/3.d0
       else
          GreenFunction = dt*4.d0*(Pot+(1.d0-a_1)*dt*dt*F2/12.d0)/3.d0
       end if

    else if (opt==1) then
       
       if (mod(ib,2)==0) then
          GreenFunction = 2.d0*(Pot+0.5d0*a_1*dt*dt*F2)/3.d0
       else
          GreenFunction = 4.d0*(Pot+0.25d0*(1.d0-a_1)*dt*dt*F2)/3.d0
       end if

    end if
       
    return
  end function GreenFunction

!-----------------------------------------------------------------------

  function Force(k,xij,rij)

    implicit none 
    
    real (kind=8)    :: rij,Force
    real (kind=8)    :: costheta
    real (kind=8)    :: dVdr
    integer (kind=4) :: k

    real (kind=8),dimension (dim) :: xij
    
    costheta = xij(1)/rij

    dVdr = -3.d0/rij**4

    Force = dVdr*(1.d0-5.d0*V0*costheta**2)*xij(k)/rij

    if (k==1) Force = Force+2.d0*V0*dVdr*xij(k)/rij 
    
    return 
  end function Force

!!$  function Force(k,xij,rij)
!!$
!!$    implicit none 
!!$    
!!$    real (kind=8)    :: rij,Force
!!$    real (kind=8)    :: dVdr
!!$    integer (kind=4) :: k
!!$
!!$    real (kind=8),dimension (dim) :: xij
!!$    
!!$    dVdr  = 0.5d0*(Potential(abs(rij+dr))-Potential(abs(rij-dr)))/dr
!!$    Force = dVdr*xij(k)/rij
!!$
!!$    return 
!!$  end function Force

!-----------------------------------------------------------------------

  function OBDMGuess(r)

    implicit none

    real (kind=8) :: r,OBDMGuess

    OBDMGuess = (1.d0-N0)*exp(-AK*r*r)+N0

    return
  end function OBDMGuess

!-----------------------------------------------------------------------

  subroutine ReadParameters(resume,crystal,diagonal,wf_table,sampling,&
       & density,alpha,dt,a_1,t_0,delta_cm,Rm,Ak,N0,dim,Np,Nb,seed,&
       & Lstag,Nlev,Nstag,Nmax,Nobdm,Nblock,Nstep,Nbin,Nk)
    
    implicit none

    logical           :: resume, crystal, diagonal, wf_table
    character (len=3) :: sampling
    real (kind=8)     :: density, alpha, dt, a_1, t_0, delta_cm
    real (kind=8)     :: Rm, Ak, N0
    integer (kind=4)  :: dim, Np, Nb, seed, Lstag, Nstag, Nmax
    integer (kind=4)  :: Nobdm, Nblock, Nstep, Nbin, Nk, Nlev

    open (unit=1, file='vpi.in', status='old')

    read (1,*)
    read (1,*)
    read (1,*) resume
    read (1,*)
    read (1,*)
    read (1,*) 
    read (1,*) dim
    read (1,*) Np
    read (1,*) density
    read (1,*) alpha
    read (1,*) crystal
    read (1,*)
    read (1,*)
    read (1,*)
    read (1,*) diagonal
    read (1,*) dt
    read (1,*) Nb
    read (1,*) seed
    read (1,*) a_1
    read (1,*) t_0
    read (1,*) delta_cm
    read (1,*) sampling
    if (sampling=="sta") then
       read (1,*) Lstag
       read (1,*) 
    else
       if (sampling=="bis") then
          read (1,*) 
          read (1,*) Nlev
       else
          print *, ' '
          print *, 'ERROR!!!'
          print *, ' '
          print *, 'The sampling method must be one of these two options:'
          print *, ' '
          print *, ' 1. sta = For staging method'
          print *, ' 2. bis = For bisection method'
          print *, ' '
          print *, 'Set the parameter "sampling" to one of the valid options'
          print *, 'and try again'
          print *, ' '
          stop
       end if
    end if
    read (1,*) Nstag
    read (1,*)
    read (1,*)
    read (1,*) 
    read (1,*) Nmax
    read (1,*) Rm
    read (1,*) wf_table
    read (1,*)
    read (1,*)
    read (1,*) 
    read (1,*) AK
    read (1,*) N0
    read (1,*) Nobdm
    read (1,*) Npw
    read (1,*)
    read (1,*)
    read (1,*) 
    read (1,*) Nblock
    read (1,*) Nstep
    read (1,*) Nbin
    read (1,*) Nk
    
    close (unit=1)

    return
  end subroutine ReadParameters

!-----------------------------------------------------------------------

  subroutine JastrowTable(rmax,Rm,WF)

    implicit none
    
    real (kind=8)    :: Rm,rmax
    real (kind=8)    :: r
    integer (kind=4) :: i
    
    real (kind=8),dimension (0:Nmax+1) :: WF
    
    dr = rmax/real(Nmax-1)

    open (unit=1,file='jastrow.out')

    do i=1,Nmax
       
       r     = (i-1)*dr
       WF(i) = LogPsi(0,Rm,r)
       write (1,'(20g20.10e3)') r,exp(WF(i))

    end do
    
    close (unit=1)

    WF(0)      = WF(2)
    WF(Nmax+1) = WF(Nmax)

    return
  end subroutine JastrowTable
  
!-----------------------------------------------------------------------

  subroutine init(seed,Path,xend,crystal,resume)

    implicit none
    
    logical          :: crystal,resume
    integer (kind=4) :: seed
    integer (kind=4) :: j,k,ip,ib
    real (kind=8),dimension(dim,Np,0:2*Nb) :: Path
    real (kind=8),dimension(dim,2)         :: xend
    real (kind=8),dimension(dim,Np)        :: R
    
    if (resume) then
       
       open (unit=2,file='checkpoint.dat',status='old')

       do ip=1,Np

          do ib=0,2*Nb

             read (2,*) (Path(k,ip,ib),k=1,dim)

          end do
          
       end do

       do j=1,2
          
          read (2,*) (xend(k,j),k=1,dim)

       end do

       call mtgetf('rand_state','u')

    else

       call sgrnd(seed)

       if (crystal) then

          open (unit=2,file='config_ini.in',status='old')
       
          read (2,*)
          read (2,*)
          read (2,*)

          do ip=1,Np
             read (2,*) (R(k,ip),k=1,dim)
          end do

       else

          do ip=1,Np
             do k=1,dim
                R(k,ip) = Lbox(k)*grnd()
             end do
          end do
    
       end if
       
       do ib=0,2*Nb
          do ip=1,Np
             do k=1,dim
                Path(k,ip,ib) = R(k,ip)
             end do
          end do
       end do

       do j=1,2
          do k=1,dim
             xend(k,j) = Path(k,Np,Nb)+1.d-5*(2.d0*grnd()-1.d0)
          end do
       end do

    end if

    return
  end subroutine init

!-----------------------------------------------------------------------

  subroutine CheckPoint(Path,xend)
    
    implicit none

    integer (kind=4) :: j,k,ip,ib

    real (kind=8),dimension(dim,Np,0:2*Nb) :: Path
    real (kind=8),dimension(dim,2)         :: xend

    open (unit=3,file='checkpoint.dat')

    do ip=1,Np
         
       do ib=0,2*Nb

          write (3,*) (Path(k,ip,ib),k=1,dim)
            
       end do

    end do

    do j=1,2
       write (3,*) (xend(k,j),k=1,dim)
    end do

    close (unit=3)

    call mtsavef('rand_state','u')

    return
  end subroutine CheckPoint

!-----------------------------------------------------------------------

  subroutine TranslateChain(delta,LogWF,dt,ip,Path,accepted)

    implicit none

    logical          :: accept
    real (kind=8)    :: delta,dt
    real (kind=8)    :: DeltaS,SumDeltaS
    integer (kind=4) :: ip,ib,k,accepted

    real (kind=8),dimension (dim)           :: xold,xnew,dx
    real (kind=8),dimension (0:Nmax+1)      :: LogWF
    real (kind=8),dimension (dim,Np,0:2*Nb) :: Path
    real (kind=8),dimension (dim,0:2*Nb)    :: NewChain

    do k=1,dim
       dx(k) = delta*(2.d0*grnd()-1.d0)
    end do
    
    SumDeltaS = 0.d0
    
    do ib=0,2*Nb
       
       do k=1,dim
          
          xold(k) = Path(k,ip,ib)
          xnew(k) = xold(k)+dx(k)
          
          call BoundaryConditions(k,xnew(k))

          NewChain(k,ib) = xnew(k)
          
       end do
       
       call UpdateAction(LogWF,Path,ip,ib,xnew,xold,dt,DeltaS)
       
       SumDeltaS = SumDeltaS+DeltaS
       
    end do
    
    !Metropolis question
    
    if (exp(-SumDeltaS)>=1.d0) then
       accept = .True.
    else
       if (exp(-SumDeltaS)>=grnd()) then
          accept = .True.
       else
          accept = .False.
       end if
    end if
    
    if (accept) then
       
       accepted = accepted+1         
       
       do ib=0,2*Nb
          do k=1,dim
             Path(k,ip,ib) = NewChain(k,ib)
          end do
       end do
       
    end if
    
    return
  end subroutine TranslateChain

!-----------------------------------------------------------------------

  subroutine TranslateHalfChain(half,delta,LogWF,dt,ip,Path,xend,accepted)

    implicit none

    logical          :: accept
    real (kind=8)    :: delta,dt
    real (kind=8)    :: DeltaS,SumDeltaS
    real (kind=8)    :: rij2,rij
    real (kind=8)    :: Snew,Sold
    integer (kind=4) :: ip,ib,k,accepted
    integer (kind=4) :: ibi,ibf
    integer (kind=4) :: half
    
    real (kind=8),dimension (dim)           :: xij
    real (kind=8),dimension (dim)           :: xold,xnew,dx
    real (kind=8),dimension (dim,2)         :: xend
    real (kind=8),dimension (0:Nmax+1)      :: LogWF
    real (kind=8),dimension (dim,Np,0:2*Nb) :: Path
    real (kind=8),dimension (dim,0:2*Nb)    :: OldChain

    Snew = 0.d0
    Sold = 0.d0

    do k=1,dim
       Path(k,ip,Nb) = xend(k,half)
    end do

    do k=1,dim
       dx(k) = delta*(2.d0*grnd()-1.d0)
    end do
    
    SumDeltaS = 0.d0
    
    if (half==1) then
       ibi = 0
       ibf = Nb
    else
       ibi = Nb
       ibf = 2*Nb
    end if

    do ib=ibi,ibf
       do k=1,dim
          OldChain(k,ib) = Path(k,ip,ib)
       end do
    end do

    do ib=ibi,ibf
       
       do k=1,dim
          
          xold(k) = Path(k,ip,ib)
          xnew(k) = xold(k)+dx(k)
          
          call BoundaryConditions(k,xnew(k))

          Path(k,ip,ib) = xnew(k)

       end do
       
       call UpdateAction(LogWF,Path,ip,ib,xnew,xold,dt,DeltaS)

       SumDeltaS = SumDeltaS+DeltaS

    end do
           
!!$    if (half==1) then
!!$       
!!$       do k=1,dim
!!$          xij(k) = Path(k,ip,Nb)-xend(k,2)
!!$       end do
!!$
!!$       call MinimumImage(xij,rij2)
!!$
!!$       if (rij2<=rcut2) then
!!$          rij  = sqrt(rij2)
!!$          Snew = log(OBDMGuess(rij)*rij**(dim-1))
!!$       end if
!!$    
!!$       do k=1,dim
!!$          xij(k) = OldChain(k,Nb)-xend(k,2)
!!$       end do       
!!$
!!$       call MinimumImage(xij,rij2)
!!$
!!$       if (rij2<=rcut2) then
!!$          rij  = sqrt(rij2)
!!$          Sold = log(OBDMGuess(rij)*rij**(dim-1))
!!$       end if
!!$      
!!$       SumDeltaS = SumDeltaS+(Snew-Sold)
!!$    
!!$    end if
!!$
!!$    if (half==2) then
!!$
!!$       do k=1,dim
!!$          xij(k) = Path(k,ip,Nb)-xend(k,1)
!!$       end do
!!$
!!$       call MinimumImage(xij,rij2)
!!$
!!$       if (rij2<=rcut2) then
!!$          rij  = sqrt(rij2)
!!$          Snew = log(OBDMGuess(rij)*rij**(dim-1))
!!$       end if
!!$
!!$       do k=1,dim
!!$          xij(k) = OldChain(k,Nb)-xend(k,1)
!!$       end do
!!$
!!$       call MinimumImage(xij,rij2)
!!$
!!$       if (rij2<=rcut2) then
!!$          rij  = sqrt(rij2)
!!$          Sold = log(OBDMGuess(rij)*rij**(dim-1))
!!$       end if
!!$
!!$       SumDeltaS = SumDeltaS+(Snew-Sold)
!!$
!!$    end if
    
    if (exp(-SumDeltaS)>=1.d0) then
       accept = .True.
    else
       if (exp(-SumDeltaS)>=grnd()) then
          accept = .True.
       else
          accept = .False.
       end if
    end if
    
    if (accept) then
      
       accepted = accepted+1
              
       do k=1,dim
          xend(k,half) = Path(k,ip,Nb)
       end do

    else

       do ib=ibi,ibf
          do k=1,dim
             Path(k,ip,ib) = OldChain(k,ib)
          end do
       end do
       
    end if
    
    return
  end subroutine TranslateHalfChain

!-----------------------------------------------------------------------

  subroutine BeadSampling(LogWF,dt,ip,Path,accepted)

    implicit none

    logical          :: accept
    real (kind=8)    :: dt
    real (kind=8)    :: gauss1,gauss2
    real (kind=8)    :: DeltaS,sigma
    integer (kind=4) :: ip,ib,k,accepted

    real (kind=8),dimension (dim)           :: xnew,xold
    real (kind=8),dimension (dim)           :: xmid,xprev,xnext
    real (kind=8),dimension (0:Nmax+1)      :: LogWF
    real (kind=8),dimension (dim,Np,0:2*Nb) :: Path

    do ib=0,2*Nb

        !Free particle sampling

        do k=1,dim

           xold(k) = Path(k,ip,ib)               
           
           call rangauss(1.d0,0.d0,gauss1,gauss2)
               
           if (ib==0) then
              xnext(k) = xold(k)-Path(k,ip,ib+1)
              if (xnext(k)<-LboxHalf(k)) xnext(k) = xnext(k)+Lbox(k)
              if (xnext(k)> LboxHalf(k)) xnext(k) = xnext(k)-Lbox(k)
              xnext(k) = xold(k)-xnext(k)

              xmid(k) = xnext(k)
              sigma   = sqrt(dt)
           else if (ib==2*Nb) then
              xprev(k) = Path(k,ip,ib-1)-xold(k)
              if (xprev(k)<-LboxHalf(k)) xprev(k) = xprev(k)+Lbox(k)
              if (xprev(k)> LboxHalf(k)) xprev(k) = xprev(k)-Lbox(k)
              xprev(k) = xold(k)+xprev(k)
              
              xmid(k) = xprev(k)
              sigma   = sqrt(dt)
           else
              xprev(k) = Path(k,ip,ib-1)-xold(k)
              if (xprev(k)<-LboxHalf(k)) xprev(k) = xprev(k)+Lbox(k)
              if (xprev(k)> LboxHalf(k)) xprev(k) = xprev(k)-Lbox(k)
              xprev(k) = xold(k)+xprev(k)
              
              xnext(k) = xold(k)-Path(k,ip,ib+1)
              if (xnext(k)<-LboxHalf(k)) xnext(k) = xnext(k)+Lbox(k)
              if (xnext(k)> LboxHalf(k)) xnext(k) = xnext(k)-Lbox(k)
              xnext(k) = xold(k)-xnext(k)
              
              xmid(k)  = 0.5d0*(xprev(k)+xnext(k))
              sigma    = sqrt(0.5d0*dt)
           end if

           xnew(k) = xmid(k)+sigma*gauss1
           
           !Periodic boundary conditions
           
           call BoundaryConditions(k,xnew(k))

        end do
        
        call UpdateAction(LogWF,Path,ip,ib,xnew,xold,dt,DeltaS)
        
        !Metropolis question
        
        if (exp(-DeltaS)>=1.d0) then
           accept = .True.
        else
           if (exp(-DeltaS)>=grnd()) then
              accept = .True.
           else
              accept = .False.
           end if
        end if
        
        if (accept) then
           
           accepted = accepted+1
           
           do k=1,dim
              Path(k,ip,ib) = xnew(k)
           end do
           
        end if
        
     end do
     
    return
  end subroutine BeadSampling

!-----------------------------------------------------------------------

  subroutine Staging(LogWF,dt,Lstag,ip,Path,accepted)

    implicit none 

    logical          :: accept
    real (kind=8)    :: dt,sigma
    real (kind=8)    :: gauss1,gauss2
    real (kind=8)    :: DeltaS,SumDeltaS
    integer (kind=4) :: ip,ib,k,accepted
    integer (kind=4) :: Lstag,j
    integer (kind=4) :: ii,ie

    real (kind=8),dimension (dim)           :: xnew,xold
    real (kind=8),dimension (dim)           :: xmid,xprev,xnext
    real (kind=8),dimension (0:Nmax+1)      :: LogWF
    real (kind=8),dimension (dim,Np,0:2*Nb) :: Path
    real (kind=8),dimension (dim,0:2*Nb)    :: OldChain
    
    ii = int((2*Nb-Lstag+1)*grnd())
    ie = ii+Lstag

    do ib=ii,ie
       do k=1,dim
          OldChain(k,ib) = Path(k,ip,ib)
       end do
    end do
   
    SumDeltaS = 0.d0

    do j=1,Lstag-1

       do k=1,dim

          xold(k) = Path(k,ip,ii+j)               
           
          call rangauss(1.d0,0.d0,gauss1,gauss2)
               
          xprev(k) = Path(k,ip,ii+j-1)-xold(k)
          if (xprev(k)<-LboxHalf(k)) xprev(k) = xprev(k)+Lbox(k)
          if (xprev(k)> LboxHalf(k)) xprev(k) = xprev(k)-Lbox(k)
          xprev(k) = xold(k)+xprev(k)
          
          xnext(k) = xold(k)-Path(k,ip,ii+Lstag)
          if (xnext(k)<-LboxHalf(k)) xnext(k) = xnext(k)+Lbox(k)
          if (xnext(k)> LboxHalf(k)) xnext(k) = xnext(k)-Lbox(k)
          xnext(k) = xold(k)-xnext(k)
          
          xmid(k)  = (xnext(k)+xprev(k)*(Lstag-j))/real(Lstag-j+1)
          sigma    = sqrt((real(Lstag-j)/real(Lstag-j+1))*dt)
          xnew(k)  = xmid(k)+sigma*gauss1
          
          !Periodic boundary conditions
          
          call BoundaryConditions(k,xnew(k))

          Path(k,ip,ii+j) = xnew(k)
          
       end do
       
       call UpdateAction(LogWF,Path,ip,ii+j,xnew,xold,dt,DeltaS)       

       SumDeltaS = SumDeltaS+DeltaS
       
    end do

    !Metropolis question
    
    if (exp(-SumDeltaS)>=1.d0) then
       accept = .True.
    else
       if (exp(-SumDeltaS)>=grnd()) then
          accept = .True.
       else
          accept = .False.
       end if
    end if
    
    if (accept) then
       
       accepted = accepted+1         
    
    else

       do ib=ii,ie
          do k=1,dim
             Path(k,ip,ib) = OldChain(k,ib)
          end do
       end do
       
    end if

    return
  end subroutine Staging

!-----------------------------------------------------------------------

  subroutine StagingHalfChain(half,LogWF,dt,Lstag,ip,Path,xend,accepted)

    implicit none 

    logical          :: accept
    real (kind=8)    :: dt,sigma
    real (kind=8)    :: gauss1,gauss2
    real (kind=8)    :: DeltaS,SumDeltaS
    integer (kind=4) :: ip,ib,k,accepted
    integer (kind=4) :: Lstag,j
    integer (kind=4) :: ii,ie
    integer (kind=4) :: half

    real (kind=8),dimension (dim)           :: xnew,xold
    real (kind=8),dimension (dim)           :: xmid,xprev,xnext
    real (kind=8),dimension (dim,2)         :: xend
    real (kind=8),dimension (0:Nmax+1)      :: LogWF
    real (kind=8),dimension (dim,Np,0:2*Nb) :: Path
    real (kind=8),dimension (dim,0:2*Nb)    :: OldChain
    
    do k=1,dim
       Path(k,ip,Nb) = xend(k,half)
    end do
    
    if (half==1) then
       ii = int((Nb-Lstag+1)*grnd())
       ie = ii+Lstag
    else
       ii = int((Nb-Lstag+1)*grnd())+Nb
       ie = ii+Lstag
    end if
     
    do ib=ii,ie
       do k=1,dim
          OldChain(k,ib) = Path(k,ip,ib)
       end do
    end do
   
    SumDeltaS = 0.d0

    do j=1,Lstag-1

       do k=1,dim

          xold(k) = Path(k,ip,ii+j)
           
          call rangauss(1.d0,0.d0,gauss1,gauss2)
               
          xprev(k) = Path(k,ip,ii+j-1)-xold(k)
          if (xprev(k)<-LboxHalf(k)) xprev(k) = xprev(k)+Lbox(k)
          if (xprev(k)> LboxHalf(k)) xprev(k) = xprev(k)-Lbox(k)
          xprev(k) = xold(k)+xprev(k)
          
          xnext(k) = xold(k)-Path(k,ip,ii+Lstag)
          if (xnext(k)<-LboxHalf(k)) xnext(k) = xnext(k)+Lbox(k)
          if (xnext(k)> LboxHalf(k)) xnext(k) = xnext(k)-Lbox(k)
          xnext(k) = xold(k)-xnext(k)
          
          xmid(k) = (xnext(k)+xprev(k)*(Lstag-j))/real(Lstag-j+1)
          sigma   = sqrt((real(Lstag-j)/real(Lstag-j+1))*dt)
          xnew(k) = xmid(k)+sigma*gauss1
          
          !Periodic boundary conditions
          
          call BoundaryConditions(k,xnew(k))

          Path(k,ip,ii+j) = xnew(k)
          
       end do
       
       call UpdateAction(LogWF,Path,ip,ii+j,xnew,xold,dt,DeltaS) 

       SumDeltaS = SumDeltaS+DeltaS
       
    end do

    !Metropolis question
    
    if (exp(-SumDeltaS)>=1.d0) then
       accept = .True.
    else
       if (exp(-SumDeltaS)>=grnd()) then
          accept = .True.
       else
          accept = .False.
       end if
    end if
    
    if (accept) then
       
       accepted = accepted+1  

       do k=1,dim
          xend(k,half) = Path(k,ip,Nb)
       end do

    else

       !Restore the original chain

       do ib=ii,ie
          do k=1,dim
             Path(k,ip,ib) = OldChain(k,ib)
          end do
       end do

    end if
    
    return
  end subroutine StagingHalfChain

!-----------------------------------------------------------------------

  subroutine MoveHead(LogWF,dt,Lmax,ip,Path,accepted)

    implicit none 

    logical          :: accept
    real (kind=8)    :: dt,sigma
    real (kind=8)    :: gauss1,gauss2
    real (kind=8)    :: DeltaS,SumDeltaS
    integer (kind=4) :: ip,ib,k,accepted
    integer (kind=4) :: j
    integer (kind=4) :: ii,ie
    integer (kind=4) :: Lmax,Ls

    real (kind=8),dimension (dim)           :: xnew,xold
    real (kind=8),dimension (dim)           :: xmid,xprev,xnext
    real (kind=8),dimension (0:Nmax+1)      :: LogWF
    real (kind=8),dimension (dim,Np,0:2*Nb) :: Path
    real (kind=8),dimension (dim,0:2*Nb)    :: OldChain
    
    Ls = int((Lmax-1)*grnd())+2

    ii = 0
    ie = ii+Ls

    !Save the original positions of the piece of the chain
    !that will be displaced
      
    do ib=ii,ie
       do k=1,dim
          OldChain(k,ib) = Path(k,ip,ib)
       end do
    end do

    SumDeltaS = 0.d0

    !Make an initial guess for the position of first bead

    do k=1,dim

       xold(k) = Path(k,ip,ii)               
           
       call rangauss(1.d0,0.d0,gauss1,gauss2)
               
       xnext(k) = xold(k)-Path(k,ip,ie)
       if (xnext(k)<-LboxHalf(k)) xnext(k) = xnext(k)+Lbox(k)
       if (xnext(k)> LboxHalf(k)) xnext(k) = xnext(k)-Lbox(k)
       xnext(k) = xold(k)-xnext(k)
       
       xmid(k)  = xnext(k)
       sigma    = sqrt(real(Ls)*dt)
       xnew(k)  = xmid(k)+sigma*gauss1
       
       !Periodic boundary conditions
       
       call BoundaryConditions(k,xnew(k))
       
       Path(k,ip,ii) = xnew(k)
       
    end do

    call UpdateAction(LogWF,Path,ip,ii,xnew,xold,dt,DeltaS)       

    SumDeltaS = SumDeltaS+DeltaS

    !Reconstruction of the whole chain piece using Staging

    do j=1,Ls-1

       do k=1,dim

          xold(k) = Path(k,ip,ii+j)               
           
          call rangauss(1.d0,0.d0,gauss1,gauss2)
               
          xprev(k) = Path(k,ip,ii+j-1)-xold(k)
          if (xprev(k)<-LboxHalf(k)) xprev(k) = xprev(k)+Lbox(k)
          if (xprev(k)> LboxHalf(k)) xprev(k) = xprev(k)-Lbox(k)
          xprev(k) = xold(k)+xprev(k)
          
          xnext(k) = xold(k)-Path(k,ip,ii+Ls)
          if (xnext(k)<-LboxHalf(k)) xnext(k) = xnext(k)+Lbox(k)
          if (xnext(k)> LboxHalf(k)) xnext(k) = xnext(k)-Lbox(k)
          xnext(k) = xold(k)-xnext(k)
          
          xmid(k)  = (xnext(k)+xprev(k)*(Ls-j))/real(Ls-j+1)
          sigma    = sqrt((real(Ls-j)/real(Ls-j+1))*dt)
          xnew(k)  = xmid(k)+sigma*gauss1
          
          !Periodic boundary conditions
          
          call BoundaryConditions(k,xnew(k))

          Path(k,ip,ii+j) = xnew(k)
          
       end do

       call UpdateAction(LogWF,Path,ip,ii+j,xnew,xold,dt,DeltaS)       

       SumDeltaS = SumDeltaS+DeltaS
       
    end do

    !Metropolis question
    
    if (exp(-SumDeltaS)>=1.d0) then
       accept = .True.
    else
       if (exp(-SumDeltaS)>=grnd()) then
          accept = .True.
       else
          accept = .False.
       end if
    end if
    
    if (accept) then
       
       accepted = accepted+1         
    
    else

       do ib=ii,ie
          do k=1,dim
             Path(k,ip,ib) = OldChain(k,ib)
          end do
       end do
       
    end if

    return
  end subroutine MoveHead

!-----------------------------------------------------------------------

  subroutine MoveHeadHalfChain(half,LogWF,dt,Lmax,ip,Path,xend,accepted)

    implicit none 

    logical          :: accept
    real (kind=8)    :: dt,sigma
    real (kind=8)    :: gauss1,gauss2
    real (kind=8)    :: DeltaS,SumDeltaS
    real (kind=8)    :: rij,rij2
    real (kind=8)    :: Snew,Sold
    integer (kind=4) :: ip,ib,k,accepted
    integer (kind=4) :: j
    integer (kind=4) :: ii,ie
    integer (kind=4) :: Lmax,Ls
    integer (kind=4) :: half

    real (kind=8),dimension (dim)           :: xij
    real (kind=8),dimension (dim)           :: xnew,xold
    real (kind=8),dimension (dim)           :: xmid,xprev,xnext
    real (kind=8),dimension (dim,2)         :: xend
    real (kind=8),dimension (0:Nmax+1)      :: LogWF
    real (kind=8),dimension (dim,Np,0:2*Nb) :: Path
    real (kind=8),dimension (dim,0:2*Nb)    :: OldChain
        
    Ls = int((Lmax-1)*grnd())+2
    
    Sold = 0.d0
    Snew = 0.d0

    do k=1,dim
       Path(k,ip,Nb) = xend(k,half)
    end do

    if (half==1) then
       ii = 0
       ie = ii+Ls
    else
       ii = Nb
       ie = ii+Ls
    end if
          
    do ib=ii,ie
       do k=1,dim
          OldChain(k,ib) = Path(k,ip,ib)
       end do
    end do

    SumDeltaS = 0.d0

    !Make an initial guess for the position of first bead

    do k=1,dim

       xold(k) = Path(k,ip,ii)              
           
       call rangauss(1.d0,0.d0,gauss1,gauss2)
               
       xnext(k) = xold(k)-Path(k,ip,ie)
       if (xnext(k)<-LboxHalf(k)) xnext(k) = xnext(k)+Lbox(k)
       if (xnext(k)> LboxHalf(k)) xnext(k) = xnext(k)-Lbox(k)
       xnext(k) = xold(k)-xnext(k)
       
       xmid(k) = xnext(k)
       sigma   = sqrt(real(Ls)*dt)
       xnew(k) = xmid(k)+sigma*gauss1
       
       !Periodic boundary conditions
       
       call BoundaryConditions(k,xnew(k))
       
       Path(k,ip,ii) = xnew(k)
       
    end do

    call UpdateAction(LogWF,Path,ip,ii,xnew,xold,dt,DeltaS)       

    SumDeltaS = SumDeltaS+DeltaS

    !Reconstruction of the whole chain piece using Staging

    do j=1,Ls-1

       do k=1,dim

          xold(k) = Path(k,ip,ii+j)               
           
          call rangauss(1.d0,0.d0,gauss1,gauss2)
               
          xprev(k) = Path(k,ip,ii+j-1)-xold(k)
          if (xprev(k)<-LboxHalf(k)) xprev(k) = xprev(k)+Lbox(k)
          if (xprev(k)> LboxHalf(k)) xprev(k) = xprev(k)-Lbox(k)
          xprev(k) = xold(k)+xprev(k)
          
          xnext(k) = xold(k)-Path(k,ip,ii+Ls)
          if (xnext(k)<-LboxHalf(k)) xnext(k) = xnext(k)+Lbox(k)
          if (xnext(k)> LboxHalf(k)) xnext(k) = xnext(k)-Lbox(k)
          xnext(k) = xold(k)-xnext(k)
          
          xmid(k)  = (xnext(k)+xprev(k)*(Ls-j))/real(Ls-j+1)
          sigma    = sqrt((real(Ls-j)/real(Ls-j+1))*dt)
          xnew(k)  = xmid(k)+sigma*gauss1
          
          !Periodic boundary conditions
          
          call BoundaryConditions(k,xnew(k))

          Path(k,ip,ii+j) = xnew(k)
          
       end do

       call UpdateAction(LogWF,Path,ip,ii+j,xnew,xold,dt,DeltaS)

       SumDeltaS = SumDeltaS+DeltaS
       
    end do

!!$    if (half==2) then
!!$       
!!$       do k=1,dim
!!$          xij(k) = Path(k,ip,Nb)-xend(k,1)
!!$       end do
!!$
!!$       call MinimumImage(xij,rij2)
!!$
!!$       if (rij2<=rcut2) then
!!$          rij  = sqrt(rij2)
!!$          Snew = log(OBDMGuess(rij)*rij**(dim-1))
!!$       end if
!!$          
!!$       do k=1,dim
!!$          xij(k) = OldChain(k,Nb)-xend(k,1)
!!$       end do
!!$
!!$       call MinimumImage(xij,rij2)
!!$
!!$       if (rij2<=rcut2) then
!!$          rij  = sqrt(rij2)
!!$          Sold = log(OBDMGuess(rij)*rij**(dim-1))
!!$       end if
!!$
!!$       SumDeltaS = SumDeltaS+(Snew-Sold)
!!$
!!$    end if

    !Metropolis question
    
    if (exp(-SumDeltaS)>=1.d0) then
       accept = .True.
    else
       if (exp(-SumDeltaS)>=grnd()) then
          accept = .True.
       else
          accept = .False.
       end if
    end if
    
    if (accept) then
       
       accepted = accepted+1  
              
       do k=1,dim
          xend(k,half) = Path(k,ip,Nb)
       end do

    else

       !Restore the original chain

       do ib=ii,ie
          do k=1,dim
             Path(k,ip,ib) = OldChain(k,ib)
          end do
       end do

    end if

    return
  end subroutine MoveHeadHalfChain

!-----------------------------------------------------------------------

  subroutine MoveTail(LogWF,dt,Lmax,ip,Path,accepted)

    implicit none 
    
    logical          :: accept
    real (kind=8)    :: dt,sigma
    real (kind=8)    :: gauss1,gauss2
    real (kind=8)    :: DeltaS,SumDeltaS
    integer (kind=4) :: ip,ib,k,accepted
    integer (kind=4) :: j
    integer (kind=4) :: ii,ie
    integer (kind=4) :: Lmax,Ls
    
    real (kind=8),dimension (dim)           :: xnew,xold
    real (kind=8),dimension (dim)           :: xmid,xprev,xnext
    real (kind=8),dimension (0:Nmax+1)      :: LogWF
    real (kind=8),dimension (dim,Np,0:2*Nb) :: Path
    real (kind=8),dimension (dim,0:2*Nb)    :: OldChain
    
    Ls = int((Lmax-1)*grnd())+2
    
    ii = 2*Nb-Ls
    ie = 2*Nb

    do ib=ii,ie
       do k=1,dim
          OldChain(k,ib) = Path(k,ip,ib)
       end do
    end do

    SumDeltaS = 0.d0
    
    !Make an initial guess for the position of last bead
    
    do k=1,dim
       
       xold(k) = Path(k,ip,ie)               
       
       call rangauss(1.d0,0.d0,gauss1,gauss2)
       
       xprev(k) = Path(k,ip,ii)-xold(k)
       if (xprev(k)<-LboxHalf(k)) xprev(k) = xprev(k)+Lbox(k)
       if (xprev(k)> LboxHalf(k)) xprev(k) = xprev(k)-Lbox(k)
       xprev(k) = xold(k)+xprev(k)
       
       xmid(k)  = xprev(k)
       sigma    = sqrt(real(Ls)*dt)
       
       xnew(k) = xmid(k)+sigma*gauss1
       
       !Periodic boundary conditions
       
       call BoundaryConditions(k,xnew(k))
       
       Path(k,ip,ie) = xnew(k)
       
    end do
    
    call UpdateAction(LogWF,Path,ip,ie,xnew,xold,dt,DeltaS)       
    
    SumDeltaS = SumDeltaS+DeltaS
    
    !Reconstruction of the whole chain piece using Staging

    do j=1,Ls-1
       
       do k=1,dim
          
          xold(k) = Path(k,ip,ii+j)               
          
          call rangauss(1.d0,0.d0,gauss1,gauss2)
          
          xprev(k) = Path(k,ip,ii+j-1)-xold(k)
          if (xprev(k)<-LboxHalf(k)) xprev(k) = xprev(k)+Lbox(k)
          if (xprev(k)> LboxHalf(k)) xprev(k) = xprev(k)-Lbox(k)
          xprev(k) = xold(k)+xprev(k)
          
          xnext(k) = xold(k)-Path(k,ip,ie)
          if (xnext(k)<-LboxHalf(k)) xnext(k) = xnext(k)+Lbox(k)
          if (xnext(k)> LboxHalf(k)) xnext(k) = xnext(k)-Lbox(k)
          xnext(k) = xold(k)-xnext(k)
          
          xmid(k)  = (xnext(k)+xprev(k)*(Ls-j))/real(Ls-j+1)
          sigma    = sqrt((real(Ls-j)/real(Ls-j+1))*dt)
          
          xnew(k) = xmid(k)+sigma*gauss1
          
          !Periodic boundary conditions
          
          call BoundaryConditions(k,xnew(k))
          
          Path(k,ip,ii+j) = xnew(k)
          
       end do
       
       call UpdateAction(LogWF,Path,ip,ii+j,xnew,xold,dt,DeltaS)       
       
       SumDeltaS = SumDeltaS+DeltaS
       
    end do
    
    !Metropolis question
    
    if (exp(-SumDeltaS)>=1.d0) then
       accept = .True.
    else
       if (exp(-SumDeltaS)>=grnd()) then
          accept = .True.
       else
          accept = .False.
       end if
    end if
    
    if (accept) then
       
       accepted = accepted+1         
       
    else
       
       do ib=ii,ie
          do k=1,dim
             Path(k,ip,ib) = OldChain(k,ib)
          end do
       end do
       
    end if
      
    return
  end subroutine MoveTail

!-----------------------------------------------------------------------

  subroutine MoveTailHalfChain(half,LogWF,dt,Lmax,ip,Path,xend,accepted)

    implicit none 
    
    logical          :: accept
    real (kind=8)    :: dt,sigma
    real (kind=8)    :: gauss1,gauss2
    real (kind=8)    :: DeltaS,SumDeltaS
    real (kind=8)    :: rij,rij2
    real (kind=8)    :: Snew,Sold
    integer (kind=4) :: ip,ib,k,accepted
    integer (kind=4) :: j
    integer (kind=4) :: ii,ie
    integer (kind=4) :: Lmax,Ls
    integer (kind=4) :: half
    
    real (kind=8),dimension (dim)           :: xij
    real (kind=8),dimension (dim)           :: xnew,xold
    real (kind=8),dimension (dim)           :: xmid,xprev,xnext
    real (kind=8),dimension (dim,2)         :: xend
    real (kind=8),dimension (0:Nmax+1)      :: LogWF
    real (kind=8),dimension (dim,Np,0:2*Nb) :: Path
    real (kind=8),dimension (dim,0:2*Nb)    :: OldChain
        
    Ls = int((Lmax-1)*grnd())+2

    Sold = 0.d0
    Snew = 0.d0
    
    do k=1,dim
       Path(k,ip,Nb) = xend(k,half)
    end do

    if (half==1) then
       ii = Nb-Ls
       ie = Nb
    else
       ii = 2*Nb-Ls
       ie = 2*Nb
    end if

    do ib=ii,ie
       do k=1,dim
          OldChain(k,ib) = Path(k,ip,ib)
       end do
    end do
    
    SumDeltaS = 0.d0
    
    !Make an initial guess for the position of last bead
    
    do k=1,dim
       
       xold(k) = Path(k,ip,ii+Ls)               
       
       call rangauss(1.d0,0.d0,gauss1,gauss2)
       
       xprev(k) = Path(k,ip,ii)-xold(k)
       if (xprev(k)<-LboxHalf(k)) xprev(k) = xprev(k)+Lbox(k)
       if (xprev(k)> LboxHalf(k)) xprev(k) = xprev(k)-Lbox(k)
       xprev(k) = xold(k)+xprev(k)
       
       xmid(k) = xprev(k)
       sigma   = sqrt(real(Ls)*dt)
       xnew(k) = xmid(k)+sigma*gauss1
       
       !Periodic boundary conditions
       
       call BoundaryConditions(k,xnew(k))
       
       Path(k,ip,ii+Ls) = xnew(k)
       
    end do
    
    call UpdateAction(LogWF,Path,ip,ii+Ls,xnew,xold,dt,DeltaS) 
    
    SumDeltaS = SumDeltaS+DeltaS

    !Reconstruction of the whole chain piece using Staging

    do j=1,Ls-1
       
       do k=1,dim
          
          xold(k) = Path(k,ip,ii+j)               
          
          call rangauss(1.d0,0.d0,gauss1,gauss2)
          
          xprev(k) = Path(k,ip,ii+j-1)-xold(k)
          if (xprev(k)<-LboxHalf(k)) xprev(k) = xprev(k)+Lbox(k)
          if (xprev(k)> LboxHalf(k)) xprev(k) = xprev(k)-Lbox(k)
          xprev(k) = xold(k)+xprev(k)
          
          xnext(k) = xold(k)-Path(k,ip,ii+Ls)
          if (xnext(k)<-LboxHalf(k)) xnext(k) = xnext(k)+Lbox(k)
          if (xnext(k)> LboxHalf(k)) xnext(k) = xnext(k)-Lbox(k)
          xnext(k) = xold(k)-xnext(k)
          
          xmid(k)  = (xnext(k)+xprev(k)*(Ls-j))/real(Ls-j+1)
          sigma    = sqrt((real(Ls-j)/real(Ls-j+1))*dt)
          
          xnew(k) = xmid(k)+sigma*gauss1
          
          !Periodic boundary conditions
          
          call BoundaryConditions(k,xnew(k))
          
          Path(k,ip,ii+j) = xnew(k)
          
       end do
       
       call UpdateAction(LogWF,Path,ip,ii+j,xnew,xold,dt,DeltaS) 

       SumDeltaS = SumDeltaS+DeltaS
       
    end do
      
!!$    if (half==1) then
!!$       
!!$       do k=1,dim
!!$          xij(k) = Path(k,ip,Nb)-xend(k,2)
!!$       end do
!!$
!!$       call MinimumImage(xij,rij2)
!!$
!!$       if (rij2<=rcut2) then
!!$          rij  = sqrt(rij2)
!!$          Snew = log(OBDMGuess(rij)*rij**(dim-1))
!!$       end if
!!$          
!!$       do k=1,dim
!!$          xij(k) = OldChain(k,Nb)-xend(k,2)
!!$       end do
!!$
!!$       call MinimumImage(xij,rij2)
!!$
!!$       if (rij2<=rcut2) then
!!$          rij  = sqrt(rij2)
!!$          Sold = log(OBDMGuess(rij)*rij**(dim-1))
!!$       end if
!!$
!!$       SumDeltaS = SumDeltaS+(Snew-Sold)
!!$    
!!$    end if
    
    !Metropolis question
    
    if (exp(-SumDeltaS)>=1.d0) then
       accept = .True.
    else
       if (exp(-SumDeltaS)>=grnd()) then
          accept = .True.
       else
          accept = .False.
       end if
    end if
    
    if (accept) then
             
       accepted = accepted+1         

       do k=1,dim
          xend(k,half) = Path(k,ip,Nb)
       end do
       
    else
    
       !Restore the original chain

       do ib=ii,ie
          do k=1,dim
             Path(k,ip,ib) = OldChain(k,ib)
          end do
       end do

    end if

    return
  end subroutine MoveTailHalfChain

!-----------------------------------------------------------------------

  subroutine Bisection(LogWF,dt,level,ip,Path,accepted)

    implicit none

    logical          :: accept
    real (kind=8)    :: dt,dt_bis
    real (kind=8)    :: gauss1,gauss2
    real (kind=8)    :: DeltaS,sigma
    real (kind=8)    :: SumDeltaS,TotDeltaS
    real (kind=8)    :: LevelDeltaS,PrevDeltaS
    integer (kind=4) :: j
    integer (kind=4) :: ii,ie
    integer (kind=4) :: ip,ib,k,accepted
    integer (kind=4) :: ilev,Nlev,level,delta_ib
    integer (kind=4) :: iprev,inext,icurr
        
    real (kind=8),dimension (dim)           :: xnew,xold
    real (kind=8),dimension (dim)           :: xmid,xprev,xnext
    real (kind=8),dimension (0:Nmax+1)      :: LogWF
    real (kind=8),dimension (dim,Np,0:2*Nb) :: Path
    real (kind=8),dimension (dim,0:2*Nb)    :: OldChain
    
    Nlev = level

    !Pick a random bead that is the starting point of the displaced 
    !piece of chain
    
    ii = int((2*Nb-2**(Nlev)+1)*grnd())
    ie = ii+2**Nlev

    !Save the original chain

    do ib=ii,ie
       do k=1,dim
          OldChain(k,ib) = Path(k,ip,ib)
       end do
    end do
    
    PrevDeltaS = 0.d0
    SumDeltaS  = 0.d0

    do ilev=1,Nlev

       delta_ib  = 2**(Nlev-ilev+1)
       dt_bis    = 0.5d0*real(delta_ib)*dt
       sigma     = sqrt(0.5d0*dt_bis)
          
       LevelDeltaS = 0.d0
    
       do j=1,2**(ilev-1)

          iprev = ii+(j-1)*delta_ib
          inext = ii+j*delta_ib
          icurr = (iprev+inext)/2

          !Free particle sampling

          do k=1,dim

             xold(k) = Path(k,ip,icurr)               

             call rangauss(1.d0,0.d0,gauss1,gauss2)

             xprev(k) = Path(k,ip,iprev)-xold(k)
             if (xprev(k)<-LboxHalf(k)) xprev(k) = xprev(k)+Lbox(k)
             if (xprev(k)> LboxHalf(k)) xprev(k) = xprev(k)-Lbox(k)
             xprev(k) = xold(k)+xprev(k)
             
             xnext(k) = xold(k)-Path(k,ip,inext)
             if (xnext(k)<-LboxHalf(k)) xnext(k) = xnext(k)+Lbox(k)
             if (xnext(k)> LboxHalf(k)) xnext(k) = xnext(k)-Lbox(k)
             xnext(k) = xold(k)-xnext(k)
             
             xmid(k) = 0.5d0*(xprev(k)+xnext(k))
             xnew(k) = xmid(k)+sigma*gauss1

             !Periodic boundary conditions

             call BoundaryConditions(k,xnew(k))

             Path(k,ip,icurr) = xnew(k)
             
          end do

          call UpdateAction(LogWF,Path,ip,icurr,xnew,xold,dt,DeltaS)
          
          LevelDeltaS = LevelDeltaS+2**(Nlev-ilev)*DeltaS
          SumDeltaS   = SumDeltaS+DeltaS

       end do

       !Evaluation of the action corresponding to the current bisection
       !level. It is important to note that the "exact" action is only
       !evaluated at the last level, in the rest of the levels of sampling
       !we use an approximate action that is S_k = 2**(Nlev-k)*S_1. 
              
       if (ilev==Nlev) then
          TotDeltaS = SumDeltaS-PrevDeltaS
       else
          TotDeltaS = LevelDeltaS-PrevDeltaS
       end if
       
       !Metropolis question

       if (exp(-TotDeltaS)>=1.d0) then
          accept     = .True.
          PrevDeltaS = LevelDeltaS
       else
          if (exp(-TotDeltaS)>=grnd()) then
             accept     = .True.
             PrevDeltaS = LevelDeltaS
          else
             accept = .False.
             exit             
          end if
       end if

    end do

    if (accept) then
    
       accepted = accepted+1
    
    else
    
       do ib=ii,ie
          do k=1,dim
             Path(k,ip,ib) = OldChain(k,ib)
          end do
       end do
    
    end if
   
    return
  end subroutine Bisection

!-----------------------------------------------------------------------

  subroutine BisectionHalf(half,LogWF,dt,level,ip,Path,xend,accepted)

    implicit none

    logical          :: accept
    real (kind=8)    :: dt,dt_bis
    real (kind=8)    :: gauss1,gauss2
    real (kind=8)    :: DeltaS,sigma
    real (kind=8)    :: SumDeltaS,TotDeltaS
    real (kind=8)    :: LevelDeltaS,PrevDeltaS
    integer (kind=4) :: j
    integer (kind=4) :: ii,ie
    integer (kind=4) :: ip,ib,k,accepted
    integer (kind=4) :: ilev,Nlev,level,delta_ib
    integer (kind=4) :: iprev,inext,icurr
    integer (kind=4) :: half
        
    real (kind=8),dimension (dim,2)         :: xend
    real (kind=8),dimension (dim)           :: xnew,xold
    real (kind=8),dimension (dim)           :: xmid,xprev,xnext
    real (kind=8),dimension (0:Nmax+1)      :: LogWF
    real (kind=8),dimension (dim,Np,0:2*Nb) :: Path
    real (kind=8),dimension (dim,0:2*Nb)    :: OldChain
    
    Nlev = level

    do k=1,dim
       Path(k,ip,Nb) = xend(k,half)
    end do

    !Pick a random bead that is the starting point of the displaced 
    !piece of chain
    
    if (half==1) then
       ii = int((Nb-2**(Nlev)+1)*grnd())
       ie = ii+2**Nlev
    else
       ii = int((Nb-2**(Nlev)+1)*grnd())+Nb
       ie = ii+2**Nlev
    end if

    !Save the original chain

    do ib=ii,ie
       do k=1,dim
          OldChain(k,ib) = Path(k,ip,ib)
       end do
    end do
    
    PrevDeltaS = 0.d0
    SumDeltaS  = 0.d0

    do ilev=1,Nlev

       delta_ib  = 2**(Nlev-ilev+1)
       dt_bis    = 0.5d0*real(delta_ib)*dt
       sigma     = sqrt(0.5d0*dt_bis)
          
       LevelDeltaS = 0.d0
    
       do j=1,2**(ilev-1)

          iprev = ii+(j-1)*delta_ib
          inext = ii+j*delta_ib
          icurr = (iprev+inext)/2

          !Free particle sampling

          do k=1,dim

             xold(k) = Path(k,ip,icurr)               

             call rangauss(1.d0,0.d0,gauss1,gauss2)

             xprev(k) = Path(k,ip,iprev)-xold(k)
             if (xprev(k)<-LboxHalf(k)) xprev(k) = xprev(k)+Lbox(k)
             if (xprev(k)> LboxHalf(k)) xprev(k) = xprev(k)-Lbox(k)
             xprev(k) = xold(k)+xprev(k)
             
             xnext(k) = xold(k)-Path(k,ip,inext)
             if (xnext(k)<-LboxHalf(k)) xnext(k) = xnext(k)+Lbox(k)
             if (xnext(k)> LboxHalf(k)) xnext(k) = xnext(k)-Lbox(k)
             xnext(k) = xold(k)-xnext(k)
             
             xmid(k) = 0.5d0*(xprev(k)+xnext(k))
             xnew(k) = xmid(k)+sigma*gauss1

             !Periodic boundary conditions

             call BoundaryConditions(k,xnew(k))

             Path(k,ip,icurr) = xnew(k)
             
          end do

          call UpdateAction(LogWF,Path,ip,icurr,xnew,xold,dt,DeltaS)
          
          LevelDeltaS = LevelDeltaS+2**(Nlev-ilev)*DeltaS
          SumDeltaS   = SumDeltaS+DeltaS

       end do

       !Evaluation of the action corresponding to the current bisection
       !level. It is important to note that the "exact" action is only
       !evaluated at the last level, in the rest of the levels of sampling
       !we use an approximate action that is S_k = 2**(Nlev-k)*S_1. 
              
       if (ilev==Nlev) then
          TotDeltaS = SumDeltaS-PrevDeltaS
       else
          TotDeltaS = LevelDeltaS-PrevDeltaS
       end if
       
       !Metropolis question

       if (exp(-TotDeltaS)>=1.d0) then
          accept     = .True.
          PrevDeltaS = LevelDeltaS
       else
          if (exp(-TotDeltaS)>=grnd()) then
             accept     = .True.
             PrevDeltaS = LevelDeltaS
          else
             accept = .False.
             exit             
          end if
       end if

    end do

    if (accept) then
    
       accepted = accepted+1

       do k=1,dim
          xend(k,half) = Path(k,ip,Nb)
       end do
    
    else
    
       do ib=ii,ie
          do k=1,dim
             Path(k,ip,ib) = OldChain(k,ib)
          end do
       end do
    
    end if
   
    return
  end subroutine BisectionHalf

!-----------------------------------------------------------------------

  subroutine MoveHeadBisection(LogWF,dt,level,ip,Path,accepted)

    implicit none

    logical          :: accept
    real (kind=8)    :: dt,dt_bis
    real (kind=8)    :: gauss1,gauss2
    real (kind=8)    :: DeltaS,sigma
    real (kind=8)    :: SumDeltaS,TotDeltaS
    real (kind=8)    :: LevelDeltaS,PrevDeltaS
    integer (kind=4) :: j
    integer (kind=4) :: ii,ie
    integer (kind=4) :: ip,ib,k,accepted
    integer (kind=4) :: ilev,Nlev,level,delta_ib
    integer (kind=4) :: iprev,inext,icurr
        
    real (kind=8),dimension (dim)           :: xnew,xold
    real (kind=8),dimension (dim)           :: xmid,xprev,xnext
    real (kind=8),dimension (0:Nmax+1)      :: LogWF
    real (kind=8),dimension (dim,Np,0:2*Nb) :: Path
    real (kind=8),dimension (dim,0:2*Nb)    :: OldChain
    
    Nlev = int((level-1)*grnd())+2

    ii = 0
    ie = ii+2**Nlev

    !Save the original positions of the piece of the chain
    !that will be displaced

    do ib=ii,ie
       do k=1,dim
          OldChain(k,ib) = Path(k,ip,ib)
       end do
    end do
    
    !Make an initial guess for the position of first bead

    SumDeltaS = 0.d0

    do k=1,dim

       xold(k) = Path(k,ip,ii)               
           
       call rangauss(1.d0,0.d0,gauss1,gauss2)
               
       xnext(k) = xold(k)-Path(k,ip,ie)
       if (xnext(k)<-LboxHalf(k)) xnext(k) = xnext(k)+Lbox(k)
       if (xnext(k)> LboxHalf(k)) xnext(k) = xnext(k)-Lbox(k)
       xnext(k) = xold(k)-xnext(k)
       
       xmid(k)  = xnext(k)
       sigma    = sqrt(2**Nlev*dt)
       
       xnew(k) = xmid(k)+sigma*gauss1
       
       !Periodic boundary conditions
       
       call BoundaryConditions(k,xnew(k))
       
       Path(k,ip,ii) = xnew(k)
       
    end do

    call UpdateAction(LogWF,Path,ip,ii,xnew,xold,dt,DeltaS)    

    SumDeltaS  = SumDeltaS+DeltaS
    PrevDeltaS = 0.d0
    
    do ilev=1,Nlev

       delta_ib = 2**(Nlev-ilev+1)
       dt_bis   = 0.5d0*real(delta_ib)*dt
       sigma    = sqrt(0.5d0*dt_bis)
          
       LevelDeltaS = 0.d0
    
       do j=1,2**(ilev-1)

          iprev = ii+(j-1)*delta_ib
          inext = ii+j*delta_ib
          icurr = (iprev+inext)/2

          !Free particle sampling

          do k=1,dim

             xold(k) = Path(k,ip,icurr)               

             call rangauss(1.d0,0.d0,gauss1,gauss2)

             xprev(k) = Path(k,ip,iprev)-xold(k)
             if (xprev(k)<-LboxHalf(k)) xprev(k) = xprev(k)+Lbox(k)
             if (xprev(k)> LboxHalf(k)) xprev(k) = xprev(k)-Lbox(k)
             xprev(k) = xold(k)+xprev(k)
             
             xnext(k) = xold(k)-Path(k,ip,inext)
             if (xnext(k)<-LboxHalf(k)) xnext(k) = xnext(k)+Lbox(k)
             if (xnext(k)> LboxHalf(k)) xnext(k) = xnext(k)-Lbox(k)
             xnext(k) = xold(k)-xnext(k)
             
             xmid(k) = 0.5d0*(xprev(k)+xnext(k))
             xnew(k) = xmid(k)+sigma*gauss1

             !Periodic boundary conditions

             call BoundaryConditions(k,xnew(k))

             Path(k,ip,icurr) = xnew(k)
             
          end do

          call UpdateAction(LogWF,Path,ip,icurr,xnew,xold,dt,DeltaS)
          
          LevelDeltaS = LevelDeltaS+2**(Nlev-ilev)*DeltaS
          SumDeltaS   = SumDeltaS+DeltaS

       end do

       !Evaluation of the action corresponding to the current bisection
       !level. It is important to note that the "exact" action is only
       !evaluated at the last level, in the rest of the levels of sampling
       !we use an approximate action that is S_k = 2**(Nlev-k)*S_1. 
              
       if (ilev==Nlev) then
          TotDeltaS = SumDeltaS-PrevDeltaS
       else
          TotDeltaS = LevelDeltaS-PrevDeltaS
       end if
       
       !Metropolis question

       if (exp(-TotDeltaS)>=1.d0) then
          accept     = .True.
          PrevDeltaS = LevelDeltaS
       else
          if (exp(-TotDeltaS)>=grnd()) then
             accept     = .True.
             PrevDeltaS = LevelDeltaS
          else
             accept = .False.
             exit             
          end if
       end if

    end do

    if (accept) then
    
       accepted = accepted+1
    
    else
    
       do ib=ii,ie
          do k=1,dim
             Path(k,ip,ib) = OldChain(k,ib)
          end do
       end do
    
    end if
   
    return
  end subroutine MoveHeadBisection

!-----------------------------------------------------------------------

  subroutine MoveHeadHalfBisection(half,LogWF,dt,level,ip,Path,xend,accepted)

    implicit none

    logical          :: accept
    real (kind=8)    :: dt,dt_bis
    real (kind=8)    :: gauss1,gauss2
    real (kind=8)    :: DeltaS,sigma
    real (kind=8)    :: SumDeltaS,TotDeltaS
    real (kind=8)    :: LevelDeltaS,PrevDeltaS
    real (kind=8)    :: rij,rij2
    real (kind=8)    :: Snew,Sold
    integer (kind=4) :: j
    integer (kind=4) :: ii,ie
    integer (kind=4) :: ip,ib,k,accepted
    integer (kind=4) :: ilev,Nlev,level,delta_ib
    integer (kind=4) :: iprev,inext,icurr
    integer (kind=4) :: half
        
    real (kind=8),dimension (dim,2)         :: xend
    real (kind=8),dimension (dim)           :: xij
    real (kind=8),dimension (dim)           :: xnew,xold
    real (kind=8),dimension (dim)           :: xmid,xprev,xnext
    real (kind=8),dimension (0:Nmax+1)      :: LogWF
    real (kind=8),dimension (dim,Np,0:2*Nb) :: Path
    real (kind=8),dimension (dim,0:2*Nb)    :: OldChain
    
    Nlev = level

    Sold = 0.d0
    Snew = 0.d0

    do k=1,dim
       Path(k,ip,Nb) = xend(k,half)
    end do

    !Pick a random bead that is the starting point of the displaced 
    !piece of chain
    
    if (half==1) then
       ii = 0
       ie = ii+2**Nlev
    else
       ii = Nb
       ie = ii+2**Nlev
    end if

    !Save the original chain

    do ib=ii,ie
       do k=1,dim
          OldChain(k,ib) = Path(k,ip,ib)
       end do
    end do
    
    !Make an initial guess for the position of first bead

    SumDeltaS = 0.d0

    do k=1,dim

       xold(k) = Path(k,ip,ii)              
           
       call rangauss(1.d0,0.d0,gauss1,gauss2)
               
       xnext(k) = xold(k)-Path(k,ip,ie)
       if (xnext(k)<-LboxHalf(k)) xnext(k) = xnext(k)+Lbox(k)
       if (xnext(k)> LboxHalf(k)) xnext(k) = xnext(k)-Lbox(k)
       xnext(k) = xold(k)-xnext(k)
       
       xmid(k) = xnext(k)
       sigma   = sqrt(real(2**Nlev)*dt)
       xnew(k) = xmid(k)+sigma*gauss1
       
       !Periodic boundary conditions
       
       call BoundaryConditions(k,xnew(k))
       
       Path(k,ip,ii) = xnew(k)
       
    end do

    call UpdateAction(LogWF,Path,ip,ii,xnew,xold,dt,DeltaS)   

    SumDeltaS  = SumDeltaS+DeltaS
    PrevDeltaS = 0.d0
    
    do ilev=1,Nlev

       delta_ib  = 2**(Nlev-ilev+1)
       dt_bis    = 0.5d0*real(delta_ib)*dt
       sigma     = sqrt(0.5d0*dt_bis)
          
       LevelDeltaS = 0.d0
    
       do j=1,2**(ilev-1)

          iprev = ii+(j-1)*delta_ib
          inext = ii+j*delta_ib
          icurr = (iprev+inext)/2

          !Free particle sampling

          do k=1,dim

             xold(k) = Path(k,ip,icurr)               

             call rangauss(1.d0,0.d0,gauss1,gauss2)

             xprev(k) = Path(k,ip,iprev)-xold(k)
             if (xprev(k)<-LboxHalf(k)) xprev(k) = xprev(k)+Lbox(k)
             if (xprev(k)> LboxHalf(k)) xprev(k) = xprev(k)-Lbox(k)
             xprev(k) = xold(k)+xprev(k)
             
             xnext(k) = xold(k)-Path(k,ip,inext)
             if (xnext(k)<-LboxHalf(k)) xnext(k) = xnext(k)+Lbox(k)
             if (xnext(k)> LboxHalf(k)) xnext(k) = xnext(k)-Lbox(k)
             xnext(k) = xold(k)-xnext(k)
             
             xmid(k) = 0.5d0*(xprev(k)+xnext(k))
             xnew(k) = xmid(k)+sigma*gauss1

             !Periodic boundary conditions

             call BoundaryConditions(k,xnew(k))

             Path(k,ip,icurr) = xnew(k)
             
          end do

          call UpdateAction(LogWF,Path,ip,icurr,xnew,xold,dt,DeltaS)
          
          LevelDeltaS = LevelDeltaS+2**(Nlev-ilev)*DeltaS
          SumDeltaS   = SumDeltaS+DeltaS

       end do

       !Evaluation of the action corresponding to the current bisection
       !level. It is important to note that the "exact" action is only
       !evaluated at the last level, in the rest of the levels of sampling
       !we use an approximate action that is S_k = 2**(Nlev-k)*S_1. 
              
       if (ilev==Nlev) then
          
          if (half==2) then
             
             do k=1,dim
                xij(k) = Path(k,ip,Nb)-xend(k,1)
             end do
             
             call MinimumImage(xij,rij2)
             
             if (rij2<=rcut2) then
                rij  = sqrt(rij2)
                Snew = log(OBDMGuess(rij)*rij**(dim-1))
             end if
             
             do k=1,dim
                xij(k) = OldChain(k,Nb)-xend(k,1)
             end do
             
             call MinimumImage(xij,rij2)
             
             if (rij2<=rcut2) then
                rij  = sqrt(rij2)
                Sold = log(OBDMGuess(rij)*rij**(dim-1))
             end if
             
             SumDeltaS = SumDeltaS+(Snew-Sold)
    
          end if
          
          TotDeltaS = SumDeltaS-PrevDeltaS
       else
          TotDeltaS = LevelDeltaS-PrevDeltaS
       end if
       
       !Metropolis question

       if (exp(-TotDeltaS)>=1.d0) then
          accept     = .True.
          PrevDeltaS = LevelDeltaS
       else
          if (exp(-TotDeltaS)>=grnd()) then
             accept     = .True.
             PrevDeltaS = LevelDeltaS
          else
             accept = .False.
             exit             
          end if
       end if

    end do

    if (accept) then
    
       accepted = accepted+1

       do k=1,dim
          xend(k,half) = Path(k,ip,Nb)
       end do
       
    else
    
       do ib=ii,ie
          do k=1,dim
             Path(k,ip,ib) = OldChain(k,ib)
          end do
       end do
    
    end if
   
    return
  end subroutine MoveHeadHalfBisection

!-----------------------------------------------------------------------

  subroutine MoveTailBisection(LogWF,dt,level,ip,Path,accepted)

    implicit none

    logical          :: accept
    real (kind=8)    :: dt,dt_bis
    real (kind=8)    :: gauss1,gauss2
    real (kind=8)    :: DeltaS,sigma
    real (kind=8)    :: SumDeltaS,TotDeltaS
    real (kind=8)    :: LevelDeltaS,PrevDeltaS
    integer (kind=4) :: j
    integer (kind=4) :: ii,ie
    integer (kind=4) :: ip,ib,k,accepted
    integer (kind=4) :: ilev,Nlev,level,delta_ib
    integer (kind=4) :: iprev,inext,icurr
        
    real (kind=8),dimension (dim)           :: xnew,xold
    real (kind=8),dimension (dim)           :: xmid,xprev,xnext
    real (kind=8),dimension (0:Nmax+1)      :: LogWF
    real (kind=8),dimension (dim,Np,0:2*Nb) :: Path
    real (kind=8),dimension (dim,0:2*Nb)    :: OldChain
    
    Nlev = int((level-1)*grnd())+2

    !Pick a random bead that is the starting point of the displaced 
    !piece of chain

    ii = 2*Nb-2**Nlev
    ie = 2*Nb

    !Save the original chain

    do ib=ii,ie
       do k=1,dim
          OldChain(k,ib) = Path(k,ip,ib)
       end do
    end do
    
    !Make an initial guess for the position of first bead

    SumDeltaS = 0.d0

    do k=1,dim
       
       xold(k) = Path(k,ip,ie)               
       
       call rangauss(1.d0,0.d0,gauss1,gauss2)
       
       xprev(k) = Path(k,ip,ii)-xold(k)
       if (xprev(k)<-LboxHalf(k)) xprev(k) = xprev(k)+Lbox(k)
       if (xprev(k)> LboxHalf(k)) xprev(k) = xprev(k)-Lbox(k)
       xprev(k) = xold(k)+xprev(k)
       
       xmid(k)  = xprev(k)
       sigma    = sqrt(2**Nlev*dt)
       
       xnew(k) = xmid(k)+sigma*gauss1
       
       !Periodic boundary conditions
       
       call BoundaryConditions(k,xnew(k))
       
       Path(k,ip,ie) = xnew(k)
       
    end do

    call UpdateAction(LogWF,Path,ip,ie,xnew,xold,dt,DeltaS)       

    SumDeltaS  = SumDeltaS+DeltaS
    PrevDeltaS = 0.d0
    
    do ilev=1,Nlev

       delta_ib  = 2**(Nlev-ilev+1)
       dt_bis    = 0.5d0*real(delta_ib)*dt
       sigma     = sqrt(0.5d0*dt_bis)
          
       LevelDeltaS = 0.d0
    
       do j=1,2**(ilev-1)

          iprev = ii+(j-1)*delta_ib
          inext = ii+j*delta_ib
          icurr = (iprev+inext)/2

          !Free particle sampling

          do k=1,dim

             xold(k) = Path(k,ip,icurr)               

             call rangauss(1.d0,0.d0,gauss1,gauss2)

             xprev(k) = Path(k,ip,iprev)-xold(k)
             if (xprev(k)<-LboxHalf(k)) xprev(k) = xprev(k)+Lbox(k)
             if (xprev(k)> LboxHalf(k)) xprev(k) = xprev(k)-Lbox(k)
             xprev(k) = xold(k)+xprev(k)
             
             xnext(k) = xold(k)-Path(k,ip,inext)
             if (xnext(k)<-LboxHalf(k)) xnext(k) = xnext(k)+Lbox(k)
             if (xnext(k)> LboxHalf(k)) xnext(k) = xnext(k)-Lbox(k)
             xnext(k) = xold(k)-xnext(k)
             
             xmid(k) = 0.5d0*(xprev(k)+xnext(k))
             xnew(k) = xmid(k)+sigma*gauss1

             !Periodic boundary conditions

             call BoundaryConditions(k,xnew(k))

             Path(k,ip,icurr) = xnew(k)
             
          end do

          call UpdateAction(LogWF,Path,ip,icurr,xnew,xold,dt,DeltaS)
          
          LevelDeltaS = LevelDeltaS+2**(Nlev-ilev)*DeltaS
          SumDeltaS   = SumDeltaS+DeltaS

       end do

       !Evaluation of the action corresponding to the current bisection
       !level. It is important to note that the "exact" action is only
       !evaluated at the last level, in the rest of the levels of sampling
       !we use an approximate action that is S_k = 2**(Nlev-k)*S_1. 
              
       if (ilev==Nlev) then
          TotDeltaS = SumDeltaS-PrevDeltaS
       else
          TotDeltaS = LevelDeltaS-PrevDeltaS
       end if
       
       !Metropolis question

       if (exp(-TotDeltaS)>=1.d0) then
          accept     = .True.
          PrevDeltaS = LevelDeltaS
       else
          if (exp(-TotDeltaS)>=grnd()) then
             accept     = .True.
             PrevDeltaS = LevelDeltaS
          else
             accept = .False.
             exit             
          end if
       end if

    end do

    if (accept) then
    
       accepted = accepted+1
    
    else
    
       do ib=ii,ie
          do k=1,dim
             Path(k,ip,ib) = OldChain(k,ib)
          end do
       end do
    
    end if
   
    return
  end subroutine MoveTailBisection

 !-----------------------------------------------------------------------

  subroutine MoveTailHalfBisection(half,LogWF,dt,level,ip,Path,xend,accepted)

    implicit none

    logical          :: accept
    real (kind=8)    :: dt,dt_bis
    real (kind=8)    :: gauss1,gauss2
    real (kind=8)    :: DeltaS,sigma
    real (kind=8)    :: SumDeltaS,TotDeltaS
    real (kind=8)    :: LevelDeltaS,PrevDeltaS
    real (kind=8)    :: rij,rij2
    real (kind=8)    :: Snew,Sold
    integer (kind=4) :: j
    integer (kind=4) :: ii,ie
    integer (kind=4) :: ip,ib,k,accepted
    integer (kind=4) :: ilev,Nlev,level,delta_ib
    integer (kind=4) :: iprev,inext,icurr
    integer (kind=4) :: half
        
    real (kind=8),dimension (dim,2)         :: xend
    real (kind=8),dimension (dim)           :: xij
    real (kind=8),dimension (dim)           :: xnew,xold
    real (kind=8),dimension (dim)           :: xmid,xprev,xnext
    real (kind=8),dimension (0:Nmax+1)      :: LogWF
    real (kind=8),dimension (dim,Np,0:2*Nb) :: Path
    real (kind=8),dimension (dim,0:2*Nb)    :: OldChain
    
    Nlev = level
    
    Sold = 0.d0
    Snew = 0.d0

    do k=1,dim
       Path(k,ip,Nb) = xend(k,half)
    end do

    !Pick a random bead that is the starting point of the displaced 
    !piece of chain
    
    if (half==1) then
       ii = Nb-2**Nlev
       ie = Nb
    else
       ii = 2*Nb-2**Nlev
       ie = 2*Nb
    end if

    !Save the original chain

    do ib=ii,ie
       do k=1,dim
          OldChain(k,ib) = Path(k,ip,ib)
       end do
    end do
    
    !Make an initial guess for the position of first bead

    SumDeltaS = 0.d0

    do k=1,dim
       
       xold(k) = Path(k,ip,ie)               
       
       call rangauss(1.d0,0.d0,gauss1,gauss2)
       
       xprev(k) = Path(k,ip,ii)-xold(k)
       if (xprev(k)<-LboxHalf(k)) xprev(k) = xprev(k)+Lbox(k)
       if (xprev(k)> LboxHalf(k)) xprev(k) = xprev(k)-Lbox(k)
       xprev(k) = xold(k)+xprev(k)
       
       xmid(k)  = xprev(k)
       sigma    = sqrt(2**Nlev*dt)
       
       xnew(k) = xmid(k)+sigma*gauss1
       
       !Periodic boundary conditions
       
       call BoundaryConditions(k,xnew(k))
       
       Path(k,ip,ie) = xnew(k)
       
    end do

    call UpdateAction(LogWF,Path,ip,ie,xnew,xold,dt,DeltaS)    

    SumDeltaS  = SumDeltaS+DeltaS
    PrevDeltaS = 0.d0
    
    do ilev=1,Nlev

       delta_ib  = 2**(Nlev-ilev+1)
       dt_bis    = 0.5d0*real(delta_ib)*dt
       sigma     = sqrt(0.5d0*dt_bis)
          
       LevelDeltaS = 0.d0
    
       do j=1,2**(ilev-1)

          iprev = ii+(j-1)*delta_ib
          inext = ii+j*delta_ib
          icurr = (iprev+inext)/2

          !Free particle sampling

          do k=1,dim

             xold(k) = Path(k,ip,icurr)               

             call rangauss(1.d0,0.d0,gauss1,gauss2)

             xprev(k) = Path(k,ip,iprev)-xold(k)
             if (xprev(k)<-LboxHalf(k)) xprev(k) = xprev(k)+Lbox(k)
             if (xprev(k)> LboxHalf(k)) xprev(k) = xprev(k)-Lbox(k)
             xprev(k) = xold(k)+xprev(k)
             
             xnext(k) = xold(k)-Path(k,ip,inext)
             if (xnext(k)<-LboxHalf(k)) xnext(k) = xnext(k)+Lbox(k)
             if (xnext(k)> LboxHalf(k)) xnext(k) = xnext(k)-Lbox(k)
             xnext(k) = xold(k)-xnext(k)
             
             xmid(k) = 0.5d0*(xprev(k)+xnext(k))
             xnew(k) = xmid(k)+sigma*gauss1

             !Periodic boundary conditions

             call BoundaryConditions(k,xnew(k))

             Path(k,ip,icurr) = xnew(k)
             
          end do

          call UpdateAction(LogWF,Path,ip,icurr,xnew,xold,dt,DeltaS)
          
          LevelDeltaS = LevelDeltaS+2**(Nlev-ilev)*DeltaS
          SumDeltaS   = SumDeltaS+DeltaS

       end do

       !Evaluation of the action corresponding to the current bisection
       !level. It is important to note that the "exact" action is only
       !evaluated at the last level, in the rest of the levels of sampling
       !we use an approximate action that is S_k = 2**(Nlev-k)*S_1. 
              
       if (ilev==Nlev) then
          
          if (half==1) then
             
             do k=1,dim
                xij(k) = Path(k,ip,Nb)-xend(k,2)
             end do
             
             call MinimumImage(xij,rij2)
             
             if (rij2<=rcut2) then
                rij  = sqrt(rij2)
                Snew = log(OBDMGuess(rij)*rij**(dim-1))
             end if
             
             do k=1,dim
                xij(k) = OldChain(k,Nb)-xend(k,2)
             end do
             
             call MinimumImage(xij,rij2)
             
             if (rij2<=rcut2) then
                rij  = sqrt(rij2)
                Sold = log(OBDMGuess(rij)*rij**(dim-1))
             end if
             
             SumDeltaS = SumDeltaS+(Snew-Sold)
    
          end if
          
          TotDeltaS = SumDeltaS-PrevDeltaS
       else
          TotDeltaS = LevelDeltaS-PrevDeltaS
       end if
       
       !Metropolis question

       if (exp(-TotDeltaS)>=1.d0) then
          accept     = .True.
          PrevDeltaS = LevelDeltaS
       else
          if (exp(-TotDeltaS)>=grnd()) then
             accept     = .True.
             PrevDeltaS = LevelDeltaS
          else
             accept = .False.
             exit             
          end if
       end if

    end do

    if (accept) then
    
       accepted = accepted+1

       do k=1,dim
          xend(k,half) = Path(k,ip,Nb)
       end do
       
    else
    
       do ib=ii,ie
          do k=1,dim
             Path(k,ip,ib) = OldChain(k,ib)
          end do
       end do
    
    end if
   
    return
  end subroutine MoveTailHalfBisection

!-----------------------------------------------------------------------

  subroutine OpenChain(LogWF,density,dt,Lmax,ip,Path,xend,isopen,accepted)
    
    implicit none

    logical          :: isopen,accept
    real (kind=8)    :: density
    real (kind=8)    :: dt,sigma
    real (kind=8)    :: gauss1,gauss2
    real (kind=8)    :: DeltaS,SumDeltaS,DeltaK
    real (kind=8)    :: rij2
    integer (kind=4) :: half
    integer (kind=4) :: ip,ib,k,accepted
    integer (kind=4) :: j
    integer (kind=4) :: ii,ie
    integer (kind=4) :: Lmax,Ls

    real (kind=8),dimension (dim,2)         :: xend
    real (kind=8),dimension (dim)           :: xnew,xold,xij
    real (kind=8),dimension (dim)           :: xmid,xprev,xnext
    real (kind=8),dimension (0:Nmax+1)      :: LogWF
    real (kind=8),dimension (dim,Np,0:2*Nb) :: Path
    real (kind=8),dimension (dim,0:2*Nb)    :: OldChain

    Ls   = int((Lmax-1)*grnd())+2
    half = int(grnd()*2)+1
    
    SumDeltaS = -log(CWorm*density)

    if (half==1) then

       ii = Nb-Ls
       ie = Nb
       
       !Save the original positions of the piece of the chain
       !that will be displaced
      
       do ib=ii,ie
          do k=1,dim
             OldChain(k,ib) = Path(k,ip,ib)
          end do
       end do
       
       !Make an initial guess for the position of central bead
    
       do k=1,dim
       
          xold(k) = Path(k,ip,ie)               
       
          call rangauss(1.d0,0.d0,gauss1,gauss2)
          
          xprev(k) = Path(k,ip,ii)-xold(k)
          if (xprev(k)<-LboxHalf(k)) xprev(k) = xprev(k)+Lbox(k)
          if (xprev(k)> LboxHalf(k)) xprev(k) = xprev(k)-Lbox(k)
          xprev(k) = xold(k)+xprev(k)
          
          xmid(k)  = xprev(k)
          sigma    = sqrt(real(Ls)*dt)
          
          xnew(k) = xmid(k)+sigma*gauss1
          
          !Periodic boundary conditions
          
          call BoundaryConditions(k,xnew(k))
          
          Path(k,ip,ie) = xnew(k)
          
       end do
       
       !call UpdateAction(LogWF,Path,ip,ie,xnew,xold,dt,DeltaS)       
    
       !SumDeltaS = SumDeltaS+DeltaS

    else

       ii = Nb
       ie = Nb+Ls

       !Save the original positions of the piece of the chain
       !that will be displaced
      
       do ib=ii,ie
          do k=1,dim
             OldChain(k,ib) = Path(k,ip,ib)
          end do
       end do
       
       SumDeltaS = 0.d0
       
       !Make an initial guess for the position of central bead
       
       do k=1,dim
          
          xold(k) = Path(k,ip,ii)               
          
          call rangauss(1.d0,0.d0,gauss1,gauss2)
          
          xnext(k) = xold(k)-Path(k,ip,ie)
          if (xnext(k)<-LboxHalf(k)) xnext(k) = xnext(k)+Lbox(k)
          if (xnext(k)> LboxHalf(k)) xnext(k) = xnext(k)-Lbox(k)
          xnext(k) = xold(k)-xnext(k)
          
          xmid(k)  = xnext(k)
          sigma    = sqrt(real(Ls)*dt)
          xnew(k)  = xmid(k)+sigma*gauss1
          
          !Periodic boundary conditions
          
          call BoundaryConditions(k,xnew(k))
          
          Path(k,ip,ii) = xnew(k)
          
       end do
       
       !call UpdateAction(LogWF,Path,ip,ii,xnew,xold,dt,DeltaS)       
       
       !SumDeltaS = SumDeltaS+DeltaS
              
    end if

    !Reconstruction of the whole chain piece using Staging

    do j=1,Ls-1
       
       do k=1,dim
          
          xold(k) = Path(k,ip,ii+j)               
          
          call rangauss(1.d0,0.d0,gauss1,gauss2)
          
          xprev(k) = Path(k,ip,ii+j-1)-xold(k)
          if (xprev(k)<-LboxHalf(k)) xprev(k) = xprev(k)+Lbox(k)
          if (xprev(k)> LboxHalf(k)) xprev(k) = xprev(k)-Lbox(k)
          xprev(k) = xold(k)+xprev(k)
          
          xnext(k) = xold(k)-Path(k,ip,ie)
          if (xnext(k)<-LboxHalf(k)) xnext(k) = xnext(k)+Lbox(k)
          if (xnext(k)> LboxHalf(k)) xnext(k) = xnext(k)-Lbox(k)
          xnext(k) = xold(k)-xnext(k)
          
          xmid(k)  = (xnext(k)+xprev(k)*(Ls-j))/real(Ls-j+1)
          sigma    = sqrt((real(Ls-j)/real(Ls-j+1))*dt)
          
          xnew(k) = xmid(k)+sigma*gauss1
          
          !Periodic boundary conditions
          
          call BoundaryConditions(k,xnew(k))
          
          Path(k,ip,ii+j) = xnew(k)
          
       end do
       
       call UpdateAction(LogWF,Path,ip,ii+j,xnew,xold,dt,DeltaS)       
       
       SumDeltaS = SumDeltaS+DeltaS
       
    end do

    !Evaluation of the change in the kinetic action due to the fact
    !that the link is broken

    do k=1,dim
       xij(k) = Path(k,ip,ii)-Path(k,ip,ie)
    end do
    
    call MinimumImage(xij,rij2)
    
    DeltaK = -0.5d0*rij2/(real(Ls)*dt)-0.5d0*real(dim)*log(2.d0*pi*real(Ls)*dt)
     
    !Metropolis question
    
    if (exp(-SumDeltaS-DeltaK)>=1.d0) then
       accept = .True.
    else
       if (exp(-SumDeltaS-DeltaK)>=grnd()) then
          accept = .True.
       else
          accept = .False.
       end if
    end if

    if (accept) then
    
       isopen   = .true.
       accepted = accepted+1
       
       if (half==1) then
          do k=1,dim
             xend(k,1) = Path(k,ip,Nb)
             xend(k,2) = OldChain(k,Nb)
          end do
       else
          do k=1,dim
             xend(k,1) = OldChain(k,Nb)
             xend(k,2) = Path(k,ip,Nb)
          end do
       end if
    
    else
       
       do ib=ii,ie
          do k=1,dim
             Path(k,ip,ib) = OldChain(k,ib)
          end do
       end do

       do k=1,dim
          xend(k,1) = Path(k,ip,ib)
          xend(k,2) = xend(k,1)
       end do

    end if

    return
  end subroutine OpenChain

!-----------------------------------------------------------------------

  subroutine CloseChain(LogWF,density,dt,Lmax,ip,Path,xend,isopen,accepted)

    implicit none

    logical          :: isopen,accept
    real (kind=8)    :: density
    real (kind=8)    :: dt,sigma
    real (kind=8)    :: gauss1,gauss2
    real (kind=8)    :: DeltaS,SumDeltaS,DeltaK
    real (kind=8)    :: rij2
    integer (kind=4) :: half
    integer (kind=4) :: ip,ib,k,accepted
    integer (kind=4) :: j
    integer (kind=4) :: ii,ie
    integer (kind=4) :: Lmax,Ls

    real (kind=8),dimension (dim,2)         :: xend
    real (kind=8),dimension (dim)           :: xnew,xold,xij
    real (kind=8),dimension (dim)           :: xmid,xprev,xnext
    real (kind=8),dimension (0:Nmax+1)      :: LogWF
    real (kind=8),dimension (dim,Np,0:2*Nb) :: Path
    real (kind=8),dimension (dim,0:2*Nb)    :: OldChain

    Ls = int((Lmax-1)*grnd())+2

    SumDeltaS = log(CWorm*density)

    half = int(grnd()*2)+1

    if (half==1) then

       ii = Nb-Ls
       ie = Nb
       
       !Save the original positions of the piece of the chain
       !that will be displaced
      
       do ib=ii,ie
          do k=1,dim
             OldChain(k,ib) = Path(k,ip,ib)
          end do
       end do

       do k=1,dim
          Path(k,ip,ie) = xend(k,2)
       end do

    else
       
       ii = Nb
       ie = Nb+Ls
       
       !Save the original positions of the piece of the chain
       !that will be displaced
      
       do ib=ii,ie
          do k=1,dim
             OldChain(k,ib) = Path(k,ip,ib)
          end do
       end do

       do k=1,dim
          Path(k,ip,ii) = xend(k,1)
       end do

    end if

    !Reconstruction of the whole chain piece using Staging

    do j=1,Ls-1
       
       do k=1,dim
          
          xold(k) = Path(k,ip,ii+j)               
          
          call rangauss(1.d0,0.d0,gauss1,gauss2)
          
          xprev(k) = Path(k,ip,ii+j-1)-xold(k)
          if (xprev(k)<-LboxHalf(k)) xprev(k) = xprev(k)+Lbox(k)
          if (xprev(k)> LboxHalf(k)) xprev(k) = xprev(k)-Lbox(k)
          xprev(k) = xold(k)+xprev(k)
          
          xnext(k) = xold(k)-Path(k,ip,ie)
          if (xnext(k)<-LboxHalf(k)) xnext(k) = xnext(k)+Lbox(k)
          if (xnext(k)> LboxHalf(k)) xnext(k) = xnext(k)-Lbox(k)
          xnext(k) = xold(k)-xnext(k)
          
          xmid(k)  = (xnext(k)+xprev(k)*(Ls-j))/real(Ls-j+1)
          sigma    = sqrt((real(Ls-j)/real(Ls-j+1))*dt)
          
          xnew(k) = xmid(k)+sigma*gauss1
          
          !Periodic boundary conditions
          
          call BoundaryConditions(k,xnew(k))
          
          Path(k,ip,ii+j) = xnew(k)
          
       end do
       
       call UpdateAction(LogWF,Path,ip,ii+j,xnew,xold,dt,DeltaS)       
       
       SumDeltaS = SumDeltaS+DeltaS
       
    end do

    !Evaluation of the change in the kinetic action due to the fact
    !that the link is broken

    do k=1,dim
       xij(k) = Path(k,ip,ii)-Path(k,ip,ie)
    end do
    
    call MinimumImage(xij,rij2)
    
    DeltaK = -0.5d0*rij2/(real(Ls)*dt)-0.5d0*real(dim)*log(2.d0*pi*real(Ls)*dt)
    
    !Metropolis question
    
    if (exp(-SumDeltaS+DeltaK)>=1.d0) then
       accept = .True.
    else
       if (exp(-SumDeltaS+DeltaK)>=grnd()) then
          accept = .True.
       else
          accept = .False.
       end if
    end if
       
    if (accept) then
    
       isopen   = .false.
       accepted = accepted+1
       
       if (half==1) then
          do k=1,dim
             xend(k,1) = Path(k,ip,Nb)
             xend(k,2) = xend(k,1)
          end do
       else
          do k=1,dim
             xend(k,1) = Path(k,ip,Nb)
             xend(k,2) = xend(k,1)
          end do
       end if
    
    else
       
       do ib=ii,ie
          do k=1,dim
             Path(k,ip,ib) = OldChain(k,ib)
          end do
       end do

    end if
    

    return
  end subroutine CloseChain

!-----------------------------------------------------------------------

  subroutine UpdateAction(LogWF,Path,ip,ib,xnew,xold,dt,DeltaS)

    implicit none

    real (kind=8)    :: dt,DeltaS
    real (kind=8)    :: DeltaLogPsi
    real (kind=8)    :: DeltaPot,DeltaF2
    integer (kind=4) :: ib
    integer (kind=4) :: k
    integer (kind=4) :: ip,jp
        
    real (kind=8),dimension(0:Nmax+1)      :: LogWF
    real (kind=8),dimension(dim,Np,0:2*Nb) :: Path
    real (kind=8),dimension(dim,Np)        :: R
    real (kind=8),dimension(dim)           :: xnew,xold
    
    !Evaluating the difference of the potential energy between the new
    !and old configuration of the displaced bead
    
    do jp=1,Np
       do k=1,dim
          R(k,jp) = Path(k,jp,ib)
       end do
    end do
    
    call UpdatePot(ip,R,xnew,xold,DeltaPot,DeltaF2)
    
    DeltaLogPsi = 0.d0

    if (ib==0) then
       call UpdateWf(LogWF,ip,R,xnew,xold,DeltaLogPsi)
    else if (ib==2*Nb) then
       call UpdateWf(LogWF,ip,R,xnew,xold,DeltaLogPsi)
    end if

    DeltaS = -DeltaLogPsi+GreenFunction(0,ib,dt,DeltaPot,DeltaF2)
    
    return
  end subroutine UpdateAction

!-----------------------------------------------------------------------

  subroutine UpdateWf(LogWF,ip,R,xnew,xold,DeltaPsi)

    implicit none

    real (kind=8)    :: DeltaPsi
    real (kind=8)    :: urold,urnew
    real (kind=8)    :: PsiOld,PsiNew
    real (kind=8)    :: rijold2,rijnew2,rijold,rijnew
    real (kind=8)    :: Interpolate
    
    integer (kind=4) :: ip,jp
    integer (kind=4) :: k
    
    real (kind=8),dimension (0:Nmax+1) :: LogWF
    real (kind=8),dimension (dim,Np)   :: R
    real (kind=8),dimension (dim)      :: xnew,xold,xijold,xijnew
    
    PsiOld = 0.d0
    PsiNew = 0.d0
    
    do jp=1,Np
       
       if (jp/=ip) then
          
          !Distance between each pair of dipoles considering only the 
          !nearest image of each particle
          
          rijold2 = 0.d0
          rijnew2 = 0.d0
          
          do k=1,dim
             
             xijold(k) = xold(k)-R(k,jp)
             xijnew(k) = xnew(k)-R(k,jp)
             
          end do

          call MinimumImage(xijnew,rijnew2)
          call MinimumImage(xijold,rijold2)
          
          !Evaluation of the difference between the new and the old 
          !wave functions.
          !We evaluate the logarithm of the quotient of probabilities.
          
          if (rijold2<=rcut2) then
             
             rijold = sqrt(rijold2)

             if (wf_table) then
                urold = Interpolate(0,Nmax,dr,LogWF,rijold)
             else
                urold = LogPsi(0,Rm,rijold)
             end if
             
             PsiOld = PsiOld+urold
             
          end if
          
          if (rijnew2<=rcut2) then
             
             rijnew = sqrt(rijnew2)

             if (wf_table) then
                urnew = Interpolate(0,Nmax,dr,LogWF,rijnew)
             else
                urnew = LogPsi(0,Rm,rijnew)
             end if
                
             PsiNew = PsiNew+urnew
             
          end if
          
       end if
       
    end do
    
    DeltaPsi = PsiNew-PsiOld
    
  end subroutine UpdateWf
  
!-----------------------------------------------------------------------

  subroutine UpdatePot(ip,R,xnew,xold,DeltaPot,DeltaF2)

    implicit none

    real (kind=8)    :: DeltaPot,DeltaF2
    real (kind=8)    :: PotOld,PotNew
    real (kind=8)    :: Fnew2,Fold2
    real (kind=8)    :: rijnew2,rijold2
    real (kind=8)    :: rijnew,rijold
    integer (kind=4) :: ip,jp
    integer (kind=4) :: k
        
    real (kind=8),dimension(dim,Np) :: R
    real (kind=8),dimension(dim)    :: xnew,xold
    real (kind=8),dimension(dim)    :: xijnew,xijold
    real (kind=8),dimension(dim)    :: Fnew,Fold
    
    PotNew = 0.d0
    PotOld = 0.d0
    
    Fnew   = 0.d0
    Fold   = 0.d0
    
    do jp=1,Np
       
       if (jp/=ip) then
          
          rijnew2 = 0.d0
          rijold2 = 0.d0
          
          do k=1,dim
             
             xijnew(k) = xnew(k)-R(k,jp)
             xijold(k) = xold(k)-R(k,jp)
             
          end do

          call MinimumImage(xijold,rijold2)
          call MinimumImage(xijnew,rijnew2)
          
          if (rijnew2<=rcut2) then
             
             rijnew = sqrt(rijnew2)
             PotNew = PotNew+Potential(xijnew,rijnew)
             
             do k=1,dim
                
                Fnew(k) = Fnew(k)+Force(k,xijnew,rijnew)
                
             end do
             
          end if
          
          if (rijold2<=rcut2) then
             
             rijold = sqrt(rijold2)
             PotOld = PotOld+Potential(xijold,rijold)
             
             do k=1,dim
                
                Fold(k) = Fold(k)+Force(k,xijold,rijold)
                
             end do
             
          end if
          
       end if
       
    end do
    
    Fnew2 = 0.d0
    Fold2 = 0.d0
    
    do k=1,dim
       
       Fnew2 = Fnew2+Fnew(k)*Fnew(k)
       Fold2 = Fold2+Fold(k)*Fold(k)
       
    end do
    
    DeltaPot = PotNew-PotOld
    DeltaF2  = Fnew2-Fold2
    
    return 
  end subroutine UpdatePot
  
!-----------------------------------------------------------------------

end module vpi_mod
