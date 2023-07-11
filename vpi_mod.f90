module vpi_mod

use random_mod
use global_mod
use system_mod
use pbc_mod

implicit none

contains

!-----------------------------------------------------------------------

  subroutine ReadParameters(resume,crystal,wf_table,v_table,swapping,trap,&
       & sampling,density,dt,delta_cm,Rm,dim,Np,Nb,seed,&
       & CMFreq,Lstag,Nlev,Nstag,Nmax,Nobdm,Nblock,Nstep,Nbin,Nk)
    
    implicit none

    logical           :: resume, crystal, trap, swapping
    logical           :: wf_table, v_table
    character (len=3) :: sampling
    real (kind=8)     :: density, dt, delta_cm, Rm
    integer (kind=4)  :: ios
    integer (kind=4)  :: dim, Np, Nb, seed, Lstag, Nstag, Nmax
    integer (kind=4)  :: Nobdm, Nblock, Nstep, Nbin, Nk, Nlev, CMFreq

    namelist /system/ dim, Np, density, crystal, trap
    namelist /samp/ resume,dt,Nb,seed,delta_cm,CMFreq,sampling,&
         Lstag,Nlev,Nstag,Nblock,Nstep,Nbin,Nk
    namelist /obdm/ swapping, CWorm, Nobdm, Npw
    namelist /wavefun/ Nmax, wf_table, v_table


    ! Set default values to some variables

    ! system nml

    crystal = .false.
    trap    = .false.
    
    ! samp nml

    resume = .false.
    seed   = 1982
    Lstag  = 2
    Nlev   = 1
    
    ! obdm nml

    swapping = .false.
    CWorm    = 0.d0
    Nobdm    = 0
    Npw      = 0

    ! wavefun nml

    Nmax     = 10000
    wf_table = .false.
    v_table  = .false.

    ! Read from standard input

    read (5,nml=system,iostat=ios)
    write (*,nml=system)
    rewind (5)
    read (5,nml=samp,iostat=ios)
    write (*,nml=samp)
    rewind (5)
    read (5,nml=obdm,iostat=ios)
    write (*,nml=obdm)
    rewind (5)
    read (5,nml=wavefun,iostat=ios)
    write (*,nml=wavefun)
    rewind (5)

    call ReadSystemParameters(trap)

    return
  end subroutine ReadParameters

!-----------------------------------------------------------------------

  subroutine JastrowTable(rmax,WF)

    implicit none
    
    real (kind=8)    :: rmax
    real (kind=8)    :: r
    integer (kind=4) :: i
    
    real (kind=8),dimension (0:Nmax+1) :: WF
    
    dr = rmax/real(Nmax-1)

    open (unit=1,file='jastrow.out')

    do i=1,Nmax
       
       r     = (i-1)*dr
       WF(i) = LogPsi(0,Rm,r)
       write (1,'(20g20.10e3)') r,exp(WF(i)),WF(i)

    end do
    
    close (unit=1)

    WF(0)      = WF(2)
    WF(Nmax+1) = WF(Nmax)

    return
  end subroutine JastrowTable

!-----------------------------------------------------------------------

  subroutine PotentialTable(rmax,VTable)

    implicit none
    
    real (kind=8)    :: rmax
    real (kind=8)    :: r
    integer (kind=4) :: i
    
    real (kind=8),dimension (0:Nmax+1) :: VTable
    real (kind=8),dimension (dim)      :: x
    
    dr = rmax/real(Nmax-1)

    open (unit=1,file='potential.out')

    do i=1,Nmax
       
       r         = (i-1)*dr
       VTable(i) = Potential(x,r)
       write (1,'(20g20.10e3)') r,VTable(i)

    end do
    
    close (unit=1)

    VTable(0)      = VTable(2)
    VTable(Nmax+1) = VTable(Nmax)

    return
  end subroutine PotentialTable
  
!-----------------------------------------------------------------------

  subroutine init(trap,seed,Path,xend,crystal,resume,isopen,iworm)

    implicit none
    
    logical          :: trap,crystal,resume,isopen
    real (kind=8)    :: gauss1,gauss2
    integer (kind=4) :: seed,iworm
    integer (kind=4) :: j,k,ip,ib
    
    real (kind=8),dimension(dim,Np,0:2*Nb) :: Path
    real (kind=8),dimension(dim,2)         :: xend
    real (kind=8),dimension(dim,Np)        :: R
    
    if (resume) then
       
       open (unit=2,file='checkpoint.dat',status='old')

       read (2,*) trap
       read (2,*) isopen
       read (2,*) iworm

       do ip=1,Np

          do ib=0,2*Nb
             read (2,*) (Path(k,ip,ib),k=1,dim)
          end do
          
       end do

       read (2,*)
       read (2,*)

       do j=1,2
          read (2,*) (xend(k,j),k=1,dim)
       end do

       call mtgetf('rand_state','u')

    else

       call sgrnd(seed)

       if (trap) then

!!$          do ip=1,Np
!!$
!!$             do k=1,dim
!!$                call rangauss(1.d0,0.d0,gauss1,gauss2)
!!$                if (k==1) then
!!$                   R(k,ip) = a_ho(k)*gauss1
!!$                else if (k==2) then
!!$                   R(k,ip) = a_ho(k)*gauss1
!!$                else if (k==3) then
!!$                   R(k,ip) = a_ho(k)*gauss1
!!$                end if
!!$             end do
!!$
!!$          end do

          do ip=1,Np
             
             do k=1,dim
                R(k,ip) = 2.d0*a_ho(k)*(grnd()-0.5d0)
             end do

          end do

       else

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
                   R(k,ip) = Lbox(k)*(grnd()-0.5d0)
                end do
             end do
    
          end if

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
             xend(k,j) = Path(k,Np,Nb)
          end do
       end do

    end if

    return
  end subroutine init

!-----------------------------------------------------------------------

  subroutine CheckPoint(trap,Path,xend,isopen,iworm)
    
    implicit none

    logical          :: isopen,trap
    integer (kind=4) :: j,k,ip,ib,iworm

    real (kind=8),dimension(dim,Np,0:2*Nb) :: Path
    real (kind=8),dimension(dim,2)         :: xend

    open (unit=3,file='checkpoint.dat')

    if (trap) then
       write (3,*) ".True."
    else
       write (3,*) ".False."
    end if

    if (isopen) then
       write (3,*) ".True."
    else
       write (3,*) ".False."
    end if

    write (3,*) iworm

    do ip=1,Np
       do ib=0,2*Nb

          write (3,*) (Path(k,ip,ib),k=1,dim)
            
       end do
    end do

    write (3,*)
    write (3,*)

    do j=1,2
       write (3,*) (xend(k,j),k=1,dim)
    end do

    close (unit=3)

    call mtsavef('rand_state','u')

    return
  end subroutine CheckPoint

!-----------------------------------------------------------------------

  subroutine TranslateChain(trap,delta,LogWF,VTable,dt,ip,Path,accepted)

    implicit none

    logical          :: accept,trap
    real (kind=8)    :: delta,dt
    real (kind=8)    :: DeltaS,SumDeltaS
    integer (kind=4) :: ip,ib,k,accepted

    real (kind=8),dimension (dim)           :: xold,xnew,dx
    real (kind=8),dimension (0:Nmax+1)      :: LogWF,VTable
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
          
          if (trap .eqv. .false.) then
             call BoundaryConditions(k,xnew(k))
          end if

          NewChain(k,ib) = xnew(k)
          
       end do
       
       call UpdateAction(trap,LogWF,VTable,Path,ip,ib,xnew,xold,dt,DeltaS)
       
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

  subroutine TranslateHalfChain(trap,half,delta,LogWF,VTable,dt,ip,Path,xend,accepted)

    implicit none

    logical          :: accept,trap
    real (kind=8)    :: delta,dt
    real (kind=8)    :: DeltaS,SumDeltaS
    integer (kind=4) :: ip,ib,k,accepted
    integer (kind=4) :: ibi,ibf
    integer (kind=4) :: half
    
    real (kind=8),dimension (dim)           :: xold,xnew,dx
    real (kind=8),dimension (dim,2)         :: xend
    real (kind=8),dimension (0:Nmax+1)      :: LogWF,VTable
    real (kind=8),dimension (dim,Np,0:2*Nb) :: Path
    real (kind=8),dimension (dim,0:2*Nb)    :: OldChain

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
          
          if (trap .eqv. .false.) then
             call BoundaryConditions(k,xnew(k))
          end if

          Path(k,ip,ib) = xnew(k)

       end do
       
       call UpdateAction(trap,LogWF,VTable,Path,ip,ib,xnew,xold,dt,DeltaS)

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

       do ib=ibi,ibf
          do k=1,dim
             Path(k,ip,ib) = OldChain(k,ib)
          end do
       end do
       
    end if
    
    return
  end subroutine TranslateHalfChain

!-----------------------------------------------------------------------

  subroutine Staging(trap,LogWF,VTable,dt,Lstag,ip,Path,accepted)

    implicit none 

    logical          :: accept,trap
    real (kind=8)    :: dt,sigma
    real (kind=8)    :: gauss1,gauss2
    real (kind=8)    :: DeltaS,SumDeltaS
    integer (kind=4) :: ip,ib,k,accepted
    integer (kind=4) :: Lstag,j
    integer (kind=4) :: ii,ie

    real (kind=8),dimension (dim)           :: xnew,xold
    real (kind=8),dimension (dim)           :: xmid,xprev,xnext
    real (kind=8),dimension (0:Nmax+1)      :: LogWF,VTable
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
          if (trap .eqv. .false.) then
             if (xprev(k)<-LboxHalf(k)) xprev(k) = xprev(k)+Lbox(k)
             if (xprev(k)> LboxHalf(k)) xprev(k) = xprev(k)-Lbox(k)
          end if
          xprev(k) = xold(k)+xprev(k)
          
          xnext(k) = xold(k)-Path(k,ip,ii+Lstag)
          if (trap .eqv. .false.) then
             if (xnext(k)<-LboxHalf(k)) xnext(k) = xnext(k)+Lbox(k)
             if (xnext(k)> LboxHalf(k)) xnext(k) = xnext(k)-Lbox(k)
          end if
          xnext(k) = xold(k)-xnext(k)
          
          sigma   = sqrt((real(Lstag-j)/real(Lstag-j+1))*dt)
          xmid(k) = (xnext(k)+xprev(k)*(Lstag-j))/real(Lstag-j+1)
          xnew(k) = xmid(k)+sigma*gauss1
          
          !Periodic boundary conditions
          
          if (trap .eqv. .false.) then
             call BoundaryConditions(k,xnew(k))
          end if

          Path(k,ip,ii+j) = xnew(k)
          
       end do
       
       call UpdateAction(trap,LogWF,VTable,Path,ip,ii+j,xnew,xold,dt,DeltaS)       

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

  subroutine MoveHead(trap,LogWF,VTable,dt,Lmax,ip,Path,accepted)

    implicit none 

    logical          :: accept,trap
    real (kind=8)    :: dt,sigma
    real (kind=8)    :: gauss1,gauss2
    real (kind=8)    :: DeltaS,SumDeltaS
    integer (kind=4) :: ip,ib,k,accepted
    integer (kind=4) :: j
    integer (kind=4) :: ii,ie
    integer (kind=4) :: Lmax,Ls

    real (kind=8),dimension (dim)           :: xnew,xold
    real (kind=8),dimension (dim)           :: xmid,xprev,xnext
    real (kind=8),dimension (0:Nmax+1)      :: LogWF,VTable
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
       if (trap .eqv. .false.) then
          if (xnext(k)<-LboxHalf(k)) xnext(k) = xnext(k)+Lbox(k)
          if (xnext(k)> LboxHalf(k)) xnext(k) = xnext(k)-Lbox(k)
       end if
       xnext(k) = xold(k)-xnext(k)
       
       sigma   = sqrt(real(Ls)*dt)
       xmid(k) = xnext(k)
       xnew(k) = xmid(k)+sigma*gauss1
       
       !Periodic boundary conditions
       if (trap .eqv. .false.) then
          call BoundaryConditions(k,xnew(k))
       end if

       Path(k,ip,ii) = xnew(k)
       
    end do

    call UpdateAction(trap,LogWF,VTable,Path,ip,ii,xnew,xold,dt,DeltaS)       

    SumDeltaS = SumDeltaS+DeltaS

    !Reconstruction of the whole chain piece using Staging

    do j=1,Ls-1

       do k=1,dim

          xold(k) = Path(k,ip,ii+j)               
           
          call rangauss(1.d0,0.d0,gauss1,gauss2)
               
          xprev(k) = Path(k,ip,ii+j-1)-xold(k)
          if (trap .eqv. .false.) then
             if (xprev(k)<-LboxHalf(k)) xprev(k) = xprev(k)+Lbox(k)
             if (xprev(k)> LboxHalf(k)) xprev(k) = xprev(k)-Lbox(k)
          end if
          xprev(k) = xold(k)+xprev(k)
          
          xnext(k) = xold(k)-Path(k,ip,ii+Ls)
          if (trap .eqv. .false.) then
             if (xnext(k)<-LboxHalf(k)) xnext(k) = xnext(k)+Lbox(k)
             if (xnext(k)> LboxHalf(k)) xnext(k) = xnext(k)-Lbox(k)
          end if
          xnext(k) = xold(k)-xnext(k)
       
          sigma   = sqrt((real(Ls-j)/real(Ls-j+1))*dt)
          xmid(k) = (xnext(k)+xprev(k)*(Ls-j))/real(Ls-j+1)
          xnew(k) = xmid(k)+sigma*gauss1
          
          !Periodic boundary conditions
          
          if (trap .eqv. .false.) then
             call BoundaryConditions(k,xnew(k))
          end if

          Path(k,ip,ii+j) = xnew(k)
          
       end do

       call UpdateAction(trap,LogWF,VTable,Path,ip,ii+j,xnew,xold,dt,DeltaS)       

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

  subroutine MoveTail(trap,LogWF,VTable,dt,Lmax,ip,Path,accepted)

    implicit none 
    
    logical          :: accept,trap
    real (kind=8)    :: dt,sigma
    real (kind=8)    :: gauss1,gauss2
    real (kind=8)    :: DeltaS,SumDeltaS
    integer (kind=4) :: ip,ib,k,accepted
    integer (kind=4) :: j
    integer (kind=4) :: ii,ie
    integer (kind=4) :: Lmax,Ls
    
    real (kind=8),dimension (dim)           :: xnew,xold
    real (kind=8),dimension (dim)           :: xmid,xprev,xnext
    real (kind=8),dimension (0:Nmax+1)      :: LogWF,VTable
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
       if (trap .eqv. .false.) then
          if (xprev(k)<-LboxHalf(k)) xprev(k) = xprev(k)+Lbox(k)
          if (xprev(k)> LboxHalf(k)) xprev(k) = xprev(k)-Lbox(k)
       end if
       xprev(k) = xold(k)+xprev(k)
       
       sigma   = sqrt(real(Ls)*dt)
       xmid(k) = xprev(k)
       xnew(k) = xmid(k)+sigma*gauss1
       
       !Periodic boundary conditions
       
       if (trap .eqv. .false.) then
          call BoundaryConditions(k,xnew(k))
       end if

       Path(k,ip,ie) = xnew(k)
       
    end do
    
    call UpdateAction(trap,LogWF,VTable,Path,ip,ie,xnew,xold,dt,DeltaS)       
    
    SumDeltaS = SumDeltaS+DeltaS
    
    !Reconstruction of the whole chain piece using Staging

    do j=1,Ls-1
       
       do k=1,dim
          
          xold(k) = Path(k,ip,ii+j)               
          
          call rangauss(1.d0,0.d0,gauss1,gauss2)
          
          xprev(k) = Path(k,ip,ii+j-1)-xold(k)
          if (trap .eqv. .false.) then
             if (xprev(k)<-LboxHalf(k)) xprev(k) = xprev(k)+Lbox(k)
             if (xprev(k)> LboxHalf(k)) xprev(k) = xprev(k)-Lbox(k)
          end if
          xprev(k) = xold(k)+xprev(k)
          
          xnext(k) = xold(k)-Path(k,ip,ie)
          if (trap .eqv. .false.) then
             if (xnext(k)<-LboxHalf(k)) xnext(k) = xnext(k)+Lbox(k)
             if (xnext(k)> LboxHalf(k)) xnext(k) = xnext(k)-Lbox(k)
          end if
          xnext(k) = xold(k)-xnext(k)
          
          sigma   = sqrt((real(Ls-j)/real(Ls-j+1))*dt)
          xmid(k) = (xnext(k)+xprev(k)*(Ls-j))/real(Ls-j+1)
          xnew(k) = xmid(k)+sigma*gauss1
          
          !Periodic boundary conditions
          
          if (trap .eqv. .false.) then
             call BoundaryConditions(k,xnew(k))
          end if
          
          Path(k,ip,ii+j) = xnew(k)
          
       end do
       
       call UpdateAction(trap,LogWF,VTable,Path,ip,ii+j,xnew,xold,dt,DeltaS)       
       
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

  subroutine Bisection(trap,LogWF,VTable,dt,level,ip,Path,accepted)

    implicit none

    logical          :: accept,trap
    real (kind=8)    :: dt,dt_bis
    real (kind=8)    :: gauss1,gauss2
    real (kind=8)    :: DeltaS,sigma
    real (kind=8)    :: LevelDeltaS
    integer (kind=4) :: j
    integer (kind=4) :: ii,ie
    integer (kind=4) :: ip,ib,k,accepted
    integer (kind=4) :: ilev,Nlev,level,delta_ib
    integer (kind=4) :: iprev,inext,icurr
        
    real (kind=8),dimension (dim)           :: xnew,xold
    real (kind=8),dimension (dim)           :: xmid,xprev,xnext
    real (kind=8),dimension (0:Nmax+1)      :: LogWF,VTable
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
    
    !Start the multilevel sampling with bisection updates

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
             if (trap .eqv. .false.) then
                if (xprev(k)<-LboxHalf(k)) xprev(k) = xprev(k)+Lbox(k)
                if (xprev(k)> LboxHalf(k)) xprev(k) = xprev(k)-Lbox(k)
             end if
             xprev(k) = xold(k)+xprev(k)
             
             xnext(k) = xold(k)-Path(k,ip,inext)
             if (trap .eqv. .false.) then
                if (xnext(k)<-LboxHalf(k)) xnext(k) = xnext(k)+Lbox(k)
                if (xnext(k)> LboxHalf(k)) xnext(k) = xnext(k)-Lbox(k)
             end if
             xnext(k) = xold(k)-xnext(k)
             
             xmid(k) = 0.5d0*(xprev(k)+xnext(k))
             xnew(k) = xmid(k)+sigma*gauss1

             !Periodic boundary conditions

             if (trap .eqv. .false.) then
                call BoundaryConditions(k,xnew(k))
             end if

             Path(k,ip,icurr) = xnew(k)
             
          end do

          call UpdateAction(trap,LogWF,VTable,Path,ip,icurr,xnew,xold,dt,DeltaS)
          
          LevelDeltaS = LevelDeltaS+DeltaS
          
       end do

       !Metropolis question for the current level

       if (exp(-LevelDeltaS)>=1.d0) then
          accept = .True.
       else
          if (exp(-LevelDeltaS)>=grnd()) then
             accept = .True.
          else
             accept = .False.
             exit             
          end if
       end if

    end do

!!$    if (exp(-SumDeltaS)>=1.d0) then
!!$       accept = .True.
!!$    else
!!$       if (exp(-SumDeltaS)>=grnd()) then
!!$          accept = .True.
!!$       else
!!$          accept = .False.
!!$       end if
!!$    end if
    
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

  subroutine MoveHeadBisection(trap,LogWF,VTable,dt,level,ip,Path,accepted)

    implicit none

    logical          :: accept,trap,continue_bisecting
    real (kind=8)    :: dt,dt_bis
    real (kind=8)    :: gauss1,gauss2
    real (kind=8)    :: DeltaS,sigma
    real (kind=8)    :: LevelDeltaS
    integer (kind=4) :: j
    integer (kind=4) :: ii,ie
    integer (kind=4) :: ip,ib,k,accepted
    integer (kind=4) :: ilev,Nlev,level,delta_ib
    integer (kind=4) :: iprev,inext,icurr
        
    real (kind=8),dimension (dim)           :: xnew,xold
    real (kind=8),dimension (dim)           :: xmid,xprev,xnext
    real (kind=8),dimension (0:Nmax+1)      :: LogWF,VTable
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

    do k=1,dim

       xold(k) = Path(k,ip,ii)               
           
       call rangauss(1.d0,0.d0,gauss1,gauss2)
               
       xnext(k) = xold(k)-Path(k,ip,ie)
       if (trap .eqv. .false.) then
          if (xnext(k)<-LboxHalf(k)) xnext(k) = xnext(k)+Lbox(k)
          if (xnext(k)> LboxHalf(k)) xnext(k) = xnext(k)-Lbox(k)
       end if
       xnext(k) = xold(k)-xnext(k)
       
       sigma   = sqrt(2**Nlev*dt)
       xmid(k) = xnext(k)
       xnew(k) = xmid(k)+sigma*gauss1
       
       !Periodic boundary conditions
       
       if (trap .eqv. .false.) then
          call BoundaryConditions(k,xnew(k))
       end if
       
       Path(k,ip,ii) = xnew(k)
       
    end do

    call UpdateAction(trap,LogWF,VTable,Path,ip,ii,xnew,xold,dt,DeltaS)    

    if (exp(-DeltaS)>=1.d0) then
       continue_bisecting = .true.
    else
       if (exp(-DeltaS)>=grnd()) then
          continue_bisecting = .true.
       else
          continue_bisecting = .false.
       end if
    end if
    
    !If the guess for the first bead is accepted we can begin the funny 
    !multilevel part of the sampling
    
    if (continue_bisecting) then    

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
                if (trap .eqv. .false.) then
                   if (xprev(k)<-LboxHalf(k)) xprev(k) = xprev(k)+Lbox(k)
                   if (xprev(k)> LboxHalf(k)) xprev(k) = xprev(k)-Lbox(k)
                end if
                xprev(k) = xold(k)+xprev(k)
                
                xnext(k) = xold(k)-Path(k,ip,inext)
                if (trap .eqv. .false.) then
                   if (xnext(k)<-LboxHalf(k)) xnext(k) = xnext(k)+Lbox(k)
                   if (xnext(k)> LboxHalf(k)) xnext(k) = xnext(k)-Lbox(k)
                end if
                xnext(k) = xold(k)-xnext(k)
                
                xmid(k) = 0.5d0*(xprev(k)+xnext(k))
                xnew(k) = xmid(k)+sigma*gauss1
                
                !Periodic boundary conditions
                
                if (trap .eqv. .false.) then
                   call BoundaryConditions(k,xnew(k))
                end if
                
                Path(k,ip,icurr) = xnew(k)
                
             end do

             call UpdateAction(trap,LogWF,VTable,Path,ip,icurr,xnew,xold,dt,DeltaS)
             
             LevelDeltaS = LevelDeltaS+DeltaS
             
          end do

          !Metropolis question for the current level
          
          if (exp(-LevelDeltaS)>=1.d0) then
             accept = .True.
          else
             if (exp(-LevelDeltaS)>=grnd()) then
                accept = .True.
             else
                accept = .False.
                exit             
             end if
          end if
          
       end do

    else
       
       accept = .false.
       
    end if

!!$    if (exp(-SumDeltaS)>=1.d0) then
!!$       accept     = .True.
!!$    else
!!$       if (exp(-SumDeltaS)>=grnd()) then
!!$          accept = .True.
!!$       else
!!$          accept = .False.
!!$       end if
!!$    end if

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

  subroutine MoveTailBisection(trap,LogWF,VTable,dt,level,ip,Path,accepted)

    implicit none

    logical          :: accept,trap,continue_bisecting
    real (kind=8)    :: dt,dt_bis
    real (kind=8)    :: gauss1,gauss2
    real (kind=8)    :: DeltaS,sigma
    real (kind=8)    :: LevelDeltaS
    integer (kind=4) :: j
    integer (kind=4) :: ii,ie
    integer (kind=4) :: ip,ib,k,accepted
    integer (kind=4) :: ilev,Nlev,level,delta_ib
    integer (kind=4) :: iprev,inext,icurr
        
    real (kind=8),dimension (dim)           :: xnew,xold
    real (kind=8),dimension (dim)           :: xmid,xprev,xnext
    real (kind=8),dimension (0:Nmax+1)      :: LogWF,VTable
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
    
    !Make an initial guess for the position of last bead

    do k=1,dim
       
       xold(k) = Path(k,ip,ie)               
       
       call rangauss(1.d0,0.d0,gauss1,gauss2)
       
       xprev(k) = Path(k,ip,ii)-xold(k)
       if (trap .eqv. .false.) then
          if (xprev(k)<-LboxHalf(k)) xprev(k) = xprev(k)+Lbox(k)
          if (xprev(k)> LboxHalf(k)) xprev(k) = xprev(k)-Lbox(k)
       end if
       xprev(k) = xold(k)+xprev(k)
       
       sigma   = sqrt(2**Nlev*dt)
       xmid(k) = xprev(k)
       xnew(k) = xmid(k)+sigma*gauss1
       
       !Periodic boundary conditions
       
       if (trap .eqv. .false.) then
          call BoundaryConditions(k,xnew(k))
       end if
       
       Path(k,ip,ie) = xnew(k)
       
    end do

    call UpdateAction(trap,LogWF,VTable,Path,ip,ie,xnew,xold,dt,DeltaS)       

    if (exp(-DeltaS)>=1.d0) then
       continue_bisecting = .true.
    else
       if (exp(-DeltaS)>=grnd()) then
          continue_bisecting = .true.
       else
          continue_bisecting = .false.
       end if
    end if
        
    !If the guess for the first bead is accepted we can begin the funny 
    !multilevel part of the sampling

    if (continue_bisecting) then

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
                if (trap .eqv. .false.) then
                   if (xprev(k)<-LboxHalf(k)) xprev(k) = xprev(k)+Lbox(k)
                   if (xprev(k)> LboxHalf(k)) xprev(k) = xprev(k)-Lbox(k)
                end if
                xprev(k) = xold(k)+xprev(k)
                
                xnext(k) = xold(k)-Path(k,ip,inext)
                if (trap .eqv. .false.) then
                   if (xnext(k)<-LboxHalf(k)) xnext(k) = xnext(k)+Lbox(k)
                   if (xnext(k)> LboxHalf(k)) xnext(k) = xnext(k)-Lbox(k)
                end if
                xnext(k) = xold(k)-xnext(k)
                
                xmid(k) = 0.5d0*(xprev(k)+xnext(k))
                xnew(k) = xmid(k)+sigma*gauss1
                
                !Periodic boundary conditions

                if (trap .eqv. .false.) then
                   call BoundaryConditions(k,xnew(k))
                end if

                Path(k,ip,icurr) = xnew(k)
             
             end do

             call UpdateAction(trap,LogWF,VTable,Path,ip,icurr,xnew,xold,dt,DeltaS)
          
             LevelDeltaS = LevelDeltaS+DeltaS
         
          end do

          !Metropolis question for the current level

          if (exp(-LevelDeltaS)>=1.d0) then
             accept = .True.
          else
             if (exp(-LevelDeltaS)>=grnd()) then
                accept = .True.
             else
                accept = .False.
                exit             
             end if
          end if

       end do

    else

       accept = .false.

    end if

!!$    if (exp(-SumDeltaS)>=1.d0) then
!!$       accept = .True.
!!$    else
!!$       if (exp(-SumDeltaS)>=grnd()) then
!!$          accept = .True.
!!$       else
!!$          accept = .False.
!!$       end if
!!$    end if

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

  subroutine StagingHalfChain(trap,half,LogWF,VTable,dt,Lstag,ip,Path,xend,accepted)

    implicit none 

    logical          :: accept,trap
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
    real (kind=8),dimension (0:Nmax+1)      :: LogWF,VTable
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
          if (trap .eqv. .false.) then
             if (xprev(k)<-LboxHalf(k)) xprev(k) = xprev(k)+Lbox(k)
             if (xprev(k)> LboxHalf(k)) xprev(k) = xprev(k)-Lbox(k)
          end if
          xprev(k) = xold(k)+xprev(k)
          
          xnext(k) = xold(k)-Path(k,ip,ii+Lstag)
          if (trap .eqv. .false.) then
             if (xnext(k)<-LboxHalf(k)) xnext(k) = xnext(k)+Lbox(k)
             if (xnext(k)> LboxHalf(k)) xnext(k) = xnext(k)-Lbox(k)
          end if
          xnext(k) = xold(k)-xnext(k)
          
          sigma   = sqrt((real(Lstag-j)/real(Lstag-j+1))*dt)
          xmid(k) = (xnext(k)+xprev(k)*(Lstag-j))/real(Lstag-j+1)
          xnew(k) = xmid(k)+sigma*gauss1
          
          !Periodic boundary conditions
          
          if (trap .eqv. .false.) then
             call BoundaryConditions(k,xnew(k))
          end if

          Path(k,ip,ii+j) = xnew(k)
          
       end do
       
       call UpdateAction(trap,LogWF,VTable,Path,ip,ii+j,xnew,xold,dt,DeltaS) 

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

  subroutine MoveHeadHalfChain(trap,half,LogWF,VTable,dt,Lmax,ip,Path,xend,accepted)

    implicit none 

    logical          :: accept,trap
    real (kind=8)    :: dt,sigma
    real (kind=8)    :: gauss1,gauss2
    real (kind=8)    :: DeltaS,SumDeltaS
    real (kind=8)    :: Snew,Sold
    integer (kind=4) :: ip,ib,k,accepted
    integer (kind=4) :: j
    integer (kind=4) :: ii,ie
    integer (kind=4) :: Lmax,Ls
    integer (kind=4) :: half

    real (kind=8),dimension (dim)           :: xnew,xold
    real (kind=8),dimension (dim)           :: xmid,xprev,xnext
    real (kind=8),dimension (dim,2)         :: xend
    real (kind=8),dimension (0:Nmax+1)      :: LogWF,VTable
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
       if (trap .eqv. .false.) then
          if (xnext(k)<-LboxHalf(k)) xnext(k) = xnext(k)+Lbox(k)
          if (xnext(k)> LboxHalf(k)) xnext(k) = xnext(k)-Lbox(k)
       end if
       xnext(k) = xold(k)-xnext(k)
       
       sigma   = sqrt(real(Ls)*dt)
       xmid(k) = xnext(k)
       xnew(k) = xmid(k)+sigma*gauss1
       
       !Periodic boundary conditions
       
       if (trap .eqv. .false.) then
          call BoundaryConditions(k,xnew(k))
       end if
       
       Path(k,ip,ii) = xnew(k)
       
    end do

    call UpdateAction(trap,LogWF,VTable,Path,ip,ii,xnew,xold,dt,DeltaS)       

    if (half==1) then
       SumDeltaS = SumDeltaS+DeltaS
    else
       SumDeltaS = SumDeltaS+0.5d0*DeltaS
    end if

    !Reconstruction of the whole chain piece using Staging

    do j=1,Ls-1

       do k=1,dim

          xold(k) = Path(k,ip,ii+j)               
           
          call rangauss(1.d0,0.d0,gauss1,gauss2)
               
          xprev(k) = Path(k,ip,ii+j-1)-xold(k)
          if (trap .eqv. .false.) then
             if (xprev(k)<-LboxHalf(k)) xprev(k) = xprev(k)+Lbox(k)
             if (xprev(k)> LboxHalf(k)) xprev(k) = xprev(k)-Lbox(k)
          end if
          xprev(k) = xold(k)+xprev(k)
          
          xnext(k) = xold(k)-Path(k,ip,ii+Ls)
          if (trap .eqv. .false.) then
             if (xnext(k)<-LboxHalf(k)) xnext(k) = xnext(k)+Lbox(k)
             if (xnext(k)> LboxHalf(k)) xnext(k) = xnext(k)-Lbox(k)
          end if
          xnext(k) = xold(k)-xnext(k)
       
          sigma   = sqrt((real(Ls-j)/real(Ls-j+1))*dt)
          xmid(k) = (xnext(k)+xprev(k)*(Ls-j))/real(Ls-j+1)
          xnew(k) = xmid(k)+sigma*gauss1
          
          !Periodic boundary conditions
          
          if (trap .eqv. .false.) then
             call BoundaryConditions(k,xnew(k))
          end if

          Path(k,ip,ii+j) = xnew(k)
          
       end do

       call UpdateAction(trap,LogWF,VTable,Path,ip,ii+j,xnew,xold,dt,DeltaS)

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
  end subroutine MoveHeadHalfChain

!-----------------------------------------------------------------------

  subroutine MoveTailHalfChain(trap,half,LogWF,VTable,dt,Lmax,ip,Path,xend,accepted)

    implicit none 
    
    logical          :: accept,trap
    real (kind=8)    :: dt,sigma
    real (kind=8)    :: gauss1,gauss2
    real (kind=8)    :: DeltaS,SumDeltaS
    integer (kind=4) :: ip,ib,k,accepted
    integer (kind=4) :: j
    integer (kind=4) :: ii,ie
    integer (kind=4) :: Lmax,Ls
    integer (kind=4) :: half
    
    real (kind=8),dimension (dim)           :: xnew,xold
    real (kind=8),dimension (dim)           :: xmid,xprev,xnext
    real (kind=8),dimension (dim,2)         :: xend
    real (kind=8),dimension (0:Nmax+1)      :: LogWF,VTable
    real (kind=8),dimension (dim,Np,0:2*Nb) :: Path
    real (kind=8),dimension (dim,0:2*Nb)    :: OldChain
        
    Ls = int((Lmax-1)*grnd())+2

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
       if (trap .eqv. .false.) then
          if (xprev(k)<-LboxHalf(k)) xprev(k) = xprev(k)+Lbox(k)
          if (xprev(k)> LboxHalf(k)) xprev(k) = xprev(k)-Lbox(k)
       end if
       xprev(k) = xold(k)+xprev(k)
       
       sigma   = sqrt(real(Ls)*dt)
       xmid(k) = xprev(k)
       xnew(k) = xmid(k)+sigma*gauss1
       
       !Periodic boundary conditions
       
       if (trap .eqv. .false.) then
          call BoundaryConditions(k,xnew(k))
       end if
       
       Path(k,ip,ii+Ls) = xnew(k)
       
    end do
    
    call UpdateAction(trap,LogWF,VTable,Path,ip,ii+Ls,xnew,xold,dt,DeltaS) 
    
    if (half==1) then
       SumDeltaS = SumDeltaS+0.5d0*DeltaS
    else
       SumDeltaS = SumDeltaS+DeltaS
    end if

    !Reconstruction of the whole chain piece using Staging

    do j=1,Ls-1
       
       do k=1,dim
          
          xold(k) = Path(k,ip,ii+j)               
          
          call rangauss(1.d0,0.d0,gauss1,gauss2)
          
          xprev(k) = Path(k,ip,ii+j-1)-xold(k)
          if (trap .eqv. .false.) then
             if (xprev(k)<-LboxHalf(k)) xprev(k) = xprev(k)+Lbox(k)
             if (xprev(k)> LboxHalf(k)) xprev(k) = xprev(k)-Lbox(k)
          end if
          xprev(k) = xold(k)+xprev(k)
          
          xnext(k) = xold(k)-Path(k,ip,ii+Ls)
          if (trap .eqv. .false.) then
             if (xnext(k)<-LboxHalf(k)) xnext(k) = xnext(k)+Lbox(k)
             if (xnext(k)> LboxHalf(k)) xnext(k) = xnext(k)-Lbox(k)
          end if
          xnext(k) = xold(k)-xnext(k)
          
          sigma   = sqrt((real(Ls-j)/real(Ls-j+1))*dt)
          xmid(k) = (xnext(k)+xprev(k)*(Ls-j))/real(Ls-j+1)
          xnew(k) = xmid(k)+sigma*gauss1
          
          !Periodic boundary conditions
          
          if (trap .eqv. .false.) then
             call BoundaryConditions(k,xnew(k))
          end if
          
          Path(k,ip,ii+j) = xnew(k)
          
       end do
       
       call UpdateAction(trap,LogWF,VTable,Path,ip,ii+j,xnew,xold,dt,DeltaS) 

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
  end subroutine MoveTailHalfChain

!-----------------------------------------------------------------------

  subroutine OpenChain(trap,LogWF,VTable,density,dt,Lmax,ip,Path,xend,&
                      &isopen,accepted,new_perm_cycle)
    
    implicit none

    logical          :: new_perm_cycle
    logical          :: isopen,accept,trap
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
    real (kind=8),dimension (0:Nmax+1)      :: LogWF,VTable
    real (kind=8),dimension (dim,Np,0:2*Nb) :: Path
    real (kind=8),dimension (dim,0:2*Nb)    :: OldChain

    Ls   = 2*int(((Lmax-2)/2)*grnd())+2
    half = int(grnd()*2)+1
    
    SumDeltaS = -log(CWorm*density)

    if (half==1) then

       ii = Nb-Ls
       ie = Nb

       !Evaluation of the change in the kinetic action due to the fact
       !that the link is broken
       
       do k=1,dim
          xij(k) = Path(k,ip,ii)-Path(k,ip,ie)
       end do
       
       if (trap) then
          rij2 = 0.d0
          do k=1,dim
             rij2 = rij2+xij(k)*xij(k)
          end do
       else
          call MinimumImage(xij,rij2)
       end if

       DeltaK = -0.5d0*rij2/(real(Ls)*dt)-&
              & 0.5d0*real(dim)*log(2.d0*pi*real(Ls)*dt)
              
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
          if (trap .eqv. .false.) then
             if (xprev(k)<-LboxHalf(k)) xprev(k) = xprev(k)+Lbox(k)
             if (xprev(k)> LboxHalf(k)) xprev(k) = xprev(k)-Lbox(k)
          end if
          xprev(k) = xold(k)+xprev(k)
          
          sigma   = sqrt(real(Ls)*dt)
          xmid(k) = xprev(k)
          xnew(k) = xmid(k)+sigma*gauss1
          
          !Periodic boundary conditions
          
          if (trap .eqv. .false.) then
             call BoundaryConditions(k,xnew(k))
          end if

          Path(k,ip,ie) = xnew(k)
          
       end do
       
       call UpdateAction(trap,LogWF,VTable,Path,ip,ie,xnew,xold,dt,DeltaS)       
    
    else

       ii = Nb
       ie = Nb+Ls

       !Evaluation of the change in the kinetic action due to the fact
       !that the link is broken
       
       do k=1,dim
          xij(k) = Path(k,ip,ii)-Path(k,ip,ie)
       end do
       
       if (trap) then
          rij2 = 0.d0
          do k=1,dim
             rij2 = rij2+xij(k)*xij(k)
          end do
       else
          call MinimumImage(xij,rij2)
       end if

       DeltaK = -0.5d0*rij2/(real(Ls)*dt)-&
              & 0.5d0*real(dim)*log(2.d0*pi*real(Ls)*dt)
       
       !Save the original positions of the piece of the chain
       !that will be displaced
      
       do ib=ii,ie
          do k=1,dim
             OldChain(k,ib) = Path(k,ip,ib)
          end do
       end do
       
       !Make an initial guess for the position of central bead
       
       do k=1,dim
          
          xold(k) = Path(k,ip,ii)               
          
          call rangauss(1.d0,0.d0,gauss1,gauss2)
          
          xnext(k) = xold(k)-Path(k,ip,ie)
          if (trap .eqv. .false.) then
             if (xnext(k)<-LboxHalf(k)) xnext(k) = xnext(k)+Lbox(k)
             if (xnext(k)> LboxHalf(k)) xnext(k) = xnext(k)-Lbox(k)
          end if
          xnext(k) = xold(k)-xnext(k)
          
          sigma   = sqrt(real(Ls)*dt)
          xmid(k) = xnext(k)
          xnew(k) = xmid(k)+sigma*gauss1
          
          !Periodic boundary conditions
          
          if (trap .eqv. .false.) then
             call BoundaryConditions(k,xnew(k))
          end if

          Path(k,ip,ii) = xnew(k)
          
       end do
       
       call UpdateAction(trap,LogWF,VTable,Path,ip,ii,xnew,xold,dt,DeltaS)       
       
    end if
    
    SumDeltaS = SumDeltaS+0.5d0*DeltaS
    
    !Reconstruction of the whole chain piece using Staging

    do j=1,Ls-1
       
       do k=1,dim
          
          xold(k) = Path(k,ip,ii+j)               
          
          call rangauss(1.d0,0.d0,gauss1,gauss2)
          
          xprev(k) = Path(k,ip,ii+j-1)-xold(k)
          if (trap .eqv. .false.) then
             if (xprev(k)<-LboxHalf(k)) xprev(k) = xprev(k)+Lbox(k)
             if (xprev(k)> LboxHalf(k)) xprev(k) = xprev(k)-Lbox(k)
          end if
          xprev(k) = xold(k)+xprev(k)
          
          xnext(k) = xold(k)-Path(k,ip,ie)
          if (trap .eqv. .false.) then
             if (xnext(k)<-LboxHalf(k)) xnext(k) = xnext(k)+Lbox(k)
             if (xnext(k)> LboxHalf(k)) xnext(k) = xnext(k)-Lbox(k)
          end if
          xnext(k) = xold(k)-xnext(k)
          
          sigma   = sqrt((real(Ls-j)/real(Ls-j+1))*dt)
          xmid(k) = (xnext(k)+xprev(k)*(Ls-j))/real(Ls-j+1)
          xnew(k) = xmid(k)+sigma*gauss1
          
          !Periodic boundary conditions
          
          if (trap .eqv. .false.) then
             call BoundaryConditions(k,xnew(k))
          end if

          Path(k,ip,ii+j) = xnew(k)
          
       end do
       
       call UpdateAction(trap,LogWF,VTable,Path,ip,ii+j,xnew,xold,dt,DeltaS)       
       
       SumDeltaS = SumDeltaS+DeltaS
       
    end do

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

       new_perm_cycle = .true.
    
    else
       
       do ib=ii,ie
          do k=1,dim
             Path(k,ip,ib) = OldChain(k,ib)
          end do
       end do

       do k=1,dim
          xend(k,1) = Path(k,ip,Nb)
          xend(k,2) = xend(k,1)
       end do

       new_perm_cycle = .false.

    end if

    return
  end subroutine OpenChain

!-----------------------------------------------------------------------

  subroutine CloseChain(trap,LogWF,VTable,density,dt,Lmax,ip,Path,xend,&
                       &isopen,accepted,end_perm_cycle)

    implicit none

    logical          :: end_perm_cycle
    logical          :: isopen,accept,trap
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
    real (kind=8),dimension (0:Nmax+1)      :: LogWF,VTable
    real (kind=8),dimension (dim,Np,0:2*Nb) :: Path
    real (kind=8),dimension (dim,0:2*Nb)    :: OldChain

    Ls   = 2*int(((Lmax-2)/2)*grnd())+2
    half = int(grnd()*2)+1
    
    SumDeltaS = log(CWorm*density)

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
          xold(k)       = OldChain(k,ie)
          xnew(k)       = Path(k,ip,ie)          
       end do

       call UpdateAction(trap,LogWF,VTable,Path,ip,ie,xnew,xold,dt,DeltaS)       
    
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
          xold(k)       = OldChain(k,ii)
          xnew(k)       = Path(k,ip,ii)
       end do

       call UpdateAction(trap,LogWF,VTable,Path,ip,ii,xnew,xold,dt,DeltaS)       
    
    end if

    SumDeltaS = SumDeltaS+0.5d0*DeltaS
    
    !Reconstruction of the whole chain piece using Staging

    do j=1,Ls-1
       
       do k=1,dim
          
          xold(k) = Path(k,ip,ii+j)               
          
          call rangauss(1.d0,0.d0,gauss1,gauss2)
          
          xprev(k) = Path(k,ip,ii+j-1)-xold(k)
          if (trap .eqv. .false.) then
             if (xprev(k)<-LboxHalf(k)) xprev(k) = xprev(k)+Lbox(k)
             if (xprev(k)> LboxHalf(k)) xprev(k) = xprev(k)-Lbox(k)
          end if
          xprev(k) = xold(k)+xprev(k)
          
          xnext(k) = xold(k)-Path(k,ip,ie)
          if (trap .eqv. .false.) then
             if (xnext(k)<-LboxHalf(k)) xnext(k) = xnext(k)+Lbox(k)
             if (xnext(k)> LboxHalf(k)) xnext(k) = xnext(k)-Lbox(k)
          end if
          xnext(k) = xold(k)-xnext(k)
          
          sigma   = sqrt((real(Ls-j)/real(Ls-j+1))*dt)
          xmid(k) = (xnext(k)+xprev(k)*(Ls-j))/real(Ls-j+1)
          xnew(k) = xmid(k)+sigma*gauss1
          
          !Periodic boundary conditions
          
          if (trap .eqv. .false.) then
             call BoundaryConditions(k,xnew(k))
          end if

          Path(k,ip,ii+j) = xnew(k)
          
       end do
       
       call UpdateAction(trap,LogWF,VTable,Path,ip,ii+j,xnew,xold,dt,DeltaS)       
       
       SumDeltaS = SumDeltaS+DeltaS
       
    end do

    !Evaluation of the change in the kinetic action due to the fact
    !that the link is broken
    
    do k=1,dim
       xij(k) = Path(k,ip,ii)-Path(k,ip,ie)
    end do
    
    if (trap) then
       rij2 = 0.d0
       do k=1,dim
          rij2 = rij2+xij(k)*xij(k)
       end do
    else
       call MinimumImage(xij,rij2)
    end if

    DeltaK = -0.5d0*rij2/(real(Ls)*dt)-&
           & 0.5d0*real(dim)*log(2.d0*pi*real(Ls)*dt)
    
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
    
       end_perm_cycle = .true.

    else
       
       do ib=ii,ie
          do k=1,dim
             Path(k,ip,ib) = OldChain(k,ib)
          end do
       end do

       end_perm_cycle = .false.

    end if
    

    return
  end subroutine CloseChain

!-----------------------------------------------------------------------

  subroutine Swap(trap,LogWF,VTable,dt,Lmax,iw,Path,xend,accepted,ipar,&
                 &swap_accepted)
    implicit none

    logical          :: accept,trap
    logical          :: swap_accepted
    real (kind=8)    :: dt,sigma
    real (kind=8)    :: gauss1,gauss2
    real (kind=8)    :: DeltaS,SumDeltaS
    real (kind=8)    :: rij2
    real (kind=8)    :: Sk,Sw
    real (kind=8)    :: uran,sum
    integer (kind=4) :: iit
    integer (kind=4) :: ip,ib,k,accepted
    integer (kind=4) :: j,ik,iw,ipar
    integer (kind=4) :: ii,ie
    integer (kind=4) :: Lmax,Ls

    real (kind=8),dimension (Np)            :: Pp
    real (kind=8),dimension (dim,2)         :: xend
    real (kind=8),dimension (dim)           :: xnew,xold,xij
    real (kind=8),dimension (dim)           :: xmid,xprev,xnext
    real (kind=8),dimension (0:Nmax+1)      :: LogWF,VTable
    real (kind=8),dimension (dim,Np,0:2*Nb) :: Path
    real (kind=8),dimension (dim,0:2*Nb)    :: OldChain,OldWorm 
    integer (kind=4),dimension (Np)         :: Particles_in_perm_cycle
    integer (kind=4),dimension (Np)         :: Perm_histogram
    
    swap_accepted = .false.
    
    Ls = 2*int(((Lmax-2)/2)*grnd())+2
    
    ii = Nb-Ls
    ie = Nb

    !The first step is to select the particle that will become the 
    !partner of the Worm in the swap update

    Sw = 0.d0
    Pp = 0.d0

    do ip=1,Np

       do k=1,dim
          xij(k) = Path(k,ip,ii)-xend(k,2)
       end do

       if (trap) then
          rij2 = 0.d0
          do k=1,dim
             rij2 = rij2+xij(k)*xij(k)
          end do
       else
          call MinimumImage(xij,rij2)
       end if

       Pp(ip) = exp(-0.5d0*rij2/(real(Ls)*dt))
       Sw     = Sw+Pp(ip)

    end do

    !Selecting a random particle according to the probabilities that
    !we have defined above

    uran = grnd()
    ip   = 0
    sum  = 0.d0

    do
       ip  = ip+1
       sum = sum+Pp(ip)/Sw 
       if (uran <= sum) then
          ik = ip
          exit
       end if
    end do

    !Check if the chosen partner is different from the Worm itself, 
    !otherwise the update is automatically rejected.

    if (ik /= iw) then
   
       Sk = 0.d0

       do ip=1,Np
          
          do k=1,dim
             xij(k) = Path(k,ip,ii)-Path(k,ik,ie)
          end do
          
          if (trap) then
             rij2 = 0.d0
             do k=1,dim
                rij2 = rij2+xij(k)*xij(k)
             end do
          else
             call MinimumImage(xij,rij2)
          end if

          Sk = Sk+exp(-0.5d0*rij2/(real(Ls)*dt))

       end do
       
       if (grnd() <= Sw/Sk) then

          !Saving configuration of the partner of the Worm and of the
          !Worm itself

          do ib=0,2*Nb
             do k=1,dim
                OldChain(k,ib) = Path(k,ik,ib)
                OldWorm(k,ib)  = Path(k,iw,ib)
             end do
          end do

          !Set the new position of the bead Nb of chain ik to the tail of the 
          !Worm (xend(k,2))

          do k=1,dim
             Path(k,ik,ie) = xend(k,2)
          end do

          !Reconstruction of the whole chain piece using Staging

          SumDeltaS = 0.d0

          do j=1,Ls-1
       
             do k=1,dim
          
                xold(k) = Path(k,ik,ii+j)               
          
                call rangauss(1.d0,0.d0,gauss1,gauss2)
          
                xprev(k) = Path(k,ik,ii+j-1)-xold(k)
                if (trap .eqv. .false.) then
                   if (xprev(k)<-LboxHalf(k)) xprev(k) = xprev(k)+Lbox(k)
                   if (xprev(k)> LboxHalf(k)) xprev(k) = xprev(k)-Lbox(k)
                end if
                xprev(k) = xold(k)+xprev(k)
                
                xnext(k) = xold(k)-Path(k,ik,ie)
                if (trap .eqv. .false.) then
                   if (xnext(k)<-LboxHalf(k)) xnext(k) = xnext(k)+Lbox(k)
                   if (xnext(k)> LboxHalf(k)) xnext(k) = xnext(k)-Lbox(k)
                end if
                xnext(k) = xold(k)-xnext(k)
                
                sigma   = sqrt((real(Ls-j)/real(Ls-j+1))*dt)
                xmid(k) = (xnext(k)+xprev(k)*(Ls-j))/real(Ls-j+1)
                xnew(k) = xmid(k)+sigma*gauss1
                
                !Periodic boundary conditions
                
                if (trap .eqv. .false.) then
                   call BoundaryConditions(k,xnew(k))
                end if

                Path(k,ik,ii+j) = xnew(k)
                
             end do
       
             call UpdateAction(trap,LogWF,VTable,Path,ik,ii+j,xnew,xold,dt,DeltaS)  

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

             do ib=Nb,2*Nb
                do k=1,dim
                   Path(k,iw,ib) = Path(k,ik,ib)
                   Path(k,ik,ib) = OldWorm(k,ib)
                end do
             end do
             
             do k=1,dim
                xend(k,2)     = OldChain(k,Nb)
                Path(k,iw,Nb) = xend(k,2)
             end do

             swap_accepted = .true.
             ipar = ik
             
          else
             
             do ib=0,2*Nb
                do k=1,dim
                   Path(k,ik,ib) = OldChain(k,ib)
                   Path(k,iw,ib) = OldWorm(k,ib)
                end do
             end do
             
             swap_accepted = .false.
             
          end if

       end if

    end if

    return
  end subroutine Swap

!-----------------------------------------------------------------------

  subroutine UpdateAction(trap,LogWF,VTable,Path,ip,ib,xnew,xold,dt,DeltaS)

    implicit none

    logical          :: trap
    real (kind=8)    :: dt,DeltaS
    real (kind=8)    :: DeltaLogPsi
    real (kind=8)    :: DeltaPot,DeltaF2
    integer (kind=4) :: ib
    integer (kind=4) :: ip
        
    real (kind=8),dimension(0:Nmax+1)      :: LogWF,VTable
    real (kind=8),dimension(dim,Np,0:2*Nb) :: Path
    real (kind=8),dimension(dim)           :: xnew,xold
    
    !Evaluating the difference of the potential energy between the new
    !and old configuration of the displaced bead
    
    if (mod(ib,2)==0) then
       call UpdatePot(trap,VTable,ip,Path(:,:,ib),xnew,xold,DeltaPot)
       DeltaF2 = 0.d0
    else
       call UpdatePot(trap,VTable,ip,Path(:,:,ib),xnew,xold,DeltaPot,DeltaF2)
    end if

    !Contribution of the wave function at the begining and the end of the 
    !chain

    if (ib==0) then
       call UpdateWf(trap,LogWF,ip,Path(:,:,ib),xnew,xold,DeltaLogPsi)
    else if (ib==2*Nb) then
       call UpdateWf(trap,LogWF,ip,Path(:,:,ib),xnew,xold,DeltaLogPsi)
    else
       DeltaLogPsi = 0.d0
    end if

    DeltaS = -DeltaLogPsi+GreenFunction(0,ib,dt,DeltaPot,DeltaF2)
    
    return
  end subroutine UpdateAction

!-----------------------------------------------------------------------

  subroutine UpdateWf(trap,LogWF,ip,R,xnew,xold,DeltaPsi)

    implicit none

    logical          :: trap
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

    if (trap) then
       do k=1,dim
          PsiOld = PsiOld+TrapPsi(0,a_ho(k),xold(k))
          PsiNew = PsiNew+TrapPsi(0,a_ho(k),xnew(k))
       end do
    end if
    
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

          if (trap) then
             rijold2 = 0.d0
             rijnew2 = 0.d0
             do k=1,dim
                rijold2 = rijold2+xijold(k)*xijold(k)
                rijnew2 = rijnew2+xijnew(k)*xijnew(k)
             end do
          else
             call MinimumImage(xijnew,rijnew2)
             call MinimumImage(xijold,rijold2)
          end if
       
          !Evaluation of the difference between the new and the old 
          !wave functions.
          !We evaluate the logarithm of the quotient of probabilities.

          if (trap) then
             
             rijold = sqrt(rijold2)

             if (wf_table) then
                urold = Interpolate(0,Nmax,dr,LogWF,rijold)
             else
                urold = LogPsi(0,Rm,rijold)
             end if
             
             PsiOld = PsiOld+urold

             rijnew = sqrt(rijnew2)

             if (wf_table) then
                urnew = Interpolate(0,Nmax,dr,LogWF,rijnew)
             else
                urnew = LogPsi(0,Rm,rijnew)
             end if
             
             PsiNew = PsiNew+urnew

          else
          
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
          
       end if
       
    end do
    
    DeltaPsi = PsiNew-PsiOld

    return    
  end subroutine UpdateWf
  
!-----------------------------------------------------------------------

  subroutine UpdatePot(trap,VTable,ip,R,xnew,xold,DeltaPot,DeltaF2)

    implicit none

    logical          :: trap
    real (kind=8)    :: DeltaPot
    real (kind=8)    :: PotOld,PotNew
    real (kind=8)    :: Fnew2,Fold2
    real (kind=8)    :: rijnew2,rijold2
    real (kind=8)    :: rijnew,rijold
    real (kind=8)    :: Interpolate
    integer (kind=4) :: ip,jp
    integer (kind=4) :: k

    real (kind=8), optional :: DeltaF2
    
    real (kind=8),dimension(0:Nmax+1) :: VTable
    real (kind=8),dimension(dim,Np)   :: R
    real (kind=8),dimension(dim)      :: xnew,xold
    real (kind=8),dimension(dim)      :: xijnew,xijold
    real (kind=8),dimension(dim)      :: Fnew,Fold
    
    PotNew = 0.d0
    PotOld = 0.d0

    Fnew   = 0.d0
    Fold   = 0.d0

    if (trap) then
       do k=1,dim
          PotNew  = PotNew+TrapPot(0,a_ho(k),xnew(k))
          PotOld  = PotOld+TrapPot(0,a_ho(k),xold(k))
          Fold(k) = TrapPot(1,a_ho(k),xold(k))
          Fnew(k) = TrapPot(1,a_ho(k),xnew(k))
       end do
    end if
    
    do jp=1,Np
       
       if (jp/=ip) then
          
          rijnew2 = 0.d0
          rijold2 = 0.d0
          
          do k=1,dim
             
             xijnew(k) = xnew(k)-R(k,jp)
             xijold(k) = xold(k)-R(k,jp)
             
          end do

          if (trap) then
             rijold2 = 0.d0
             rijnew2 = 0.d0
             do k=1,dim
                rijold2 = rijold2+xijold(k)*xijold(k)
                rijnew2 = rijnew2+xijnew(k)*xijnew(k)
             end do
          else
             call MinimumImage(xijold,rijold2)
             call MinimumImage(xijnew,rijnew2)
          end if

          if (trap) then
             
             rijnew = sqrt(rijnew2)

             if (v_table) then
                PotNew = PotNew+Interpolate(0,Nmax,dr,VTable,rijnew)
             else
                PotNew = PotNew+Potential(xijnew,rijnew)
             end if
             
             if (present(DeltaF2)) then
                if (v_table) then
                   do k=1,dim
                      Fnew(k) = Fnew(k)+Interpolate(1,Nmax,dr,VTable,rijnew)*xijnew(k)/rijnew
                   end do
                else
                   do k=1,dim
                      Fnew(k) = Fnew(k)+Force(k,xijnew,rijnew)
                   end do
                end if
             end if
             
             if (rijold2<=rcut2) then
                
                rijold = sqrt(rijold2)
                
                if (v_table) then
                   PotOld = PotOld+Interpolate(0,Nmax,dr,VTable,rijold)
                else
                   PotOld = PotOld+Potential(xijold,rijold)
                end if
                
                if (present(DeltaF2)) then
                   if (v_table) then
                      do k=1,dim
                         Fold(k) = Fold(k)+Interpolate(1,Nmax,dr,VTable,rijold)*xijold(k)/rijold
                      end do
                   else
                      do k=1,dim
                         Fold(k) = Fold(k)+Force(k,xijold,rijold)
                      end do
                   end if
                end if
             
             end if

          else

             if (rijnew2<=rcut2) then
                
                rijnew = sqrt(rijnew2)
                
                if (v_table) then
                   PotNew = PotNew+Interpolate(0,Nmax,dr,VTable,rijnew)
                else
                   PotNew = PotNew+Potential(xijnew,rijnew)
                end if
                
                if (present(DeltaF2)) then
                   if (v_table) then
                      do k=1,dim
                         Fnew(k) = Fnew(k)+Interpolate(1,Nmax,dr,VTable,rijnew)*xijnew(k)/rijnew
                      end do
                   else
                      do k=1,dim
                         Fnew(k) = Fnew(k)+Force(k,xijnew,rijnew)
                      end do
                   end if
                end if
                
             end if
             
             if (rijold2<=rcut2) then
                
                rijold = sqrt(rijold2)
                
                if (v_table) then
                   PotOld = PotOld+Interpolate(0,Nmax,dr,VTable,rijold)
                else
                   PotOld = PotOld+Potential(xijold,rijold)
                end if
                
                if (present(DeltaF2)) then
                   if (v_table) then
                      do k=1,dim
                         Fold(k) = Fold(k)+Interpolate(1,Nmax,dr,VTable,rijold)*xijold(k)/rijold
                      end do
                   else
                      do k=1,dim
                         Fold(k) = Fold(k)+Force(k,xijold,rijold)
                      end do
                   end if
                end if
                
             end if
             
          end if
          
       end if
       
    end do
    
    Fnew2 = 0.d0
    Fold2 = 0.d0

    if (present(DeltaF2))then
       do k=1,dim
          
          Fnew2 = Fnew2+Fnew(k)*Fnew(k)
          Fold2 = Fold2+Fold(k)*Fold(k)
          
       end do
       DeltaF2  = Fnew2-Fold2
    end if

    DeltaPot = PotNew-PotOld
   
    return 
  end subroutine UpdatePot
  
!-----------------------------------------------------------------------

end module vpi_mod
