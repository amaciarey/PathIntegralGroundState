program vpi

use random_mod
use global_mod
use system_mod
use sample_mod
use vpi_mod

implicit none 

logical          :: crystal,resume,isopen
real (kind=8)    :: dt,delta_cm
real (kind=8)    :: density,alpha
real (kind=8)    :: E1,E2
real (kind=8)    :: E,Kin,Pot
real (kind=8)    :: Et,Kt
real (kind=8)    :: BlockAvE,BlockAvK,BlockAvV
real (kind=8)    :: BlockAvE2,BlockAvK2,BlockAvV2
real (kind=8)    :: BlockVarE,BlockVarK,BlockVarV
real (kind=8)    :: BlockAvEt,BlockAvKt,BlockAvVt
real (kind=8)    :: BlockAvEt2,BlockAvKt2,BlockAvVt2
real (kind=8)    :: BlockVarEt,BlockVarKt,BlockVarVt
real (kind=8)    :: AvE,AvK,AvV
real (kind=8)    :: AvE2,AvK2,AvV2
real (kind=8)    :: VarE,VarK,VarV
real (kind=8)    :: AvEt,AvKt,AvVt
real (kind=8)    :: AvEt2,AvKt2,AvVt2
real (kind=8)    :: VarEt,VarKt,VarVt
real (kind=8)    :: numz_block
real (kind=8)    :: end,begin
real (kind=8)    :: stag_move,stag_half
real (kind=8)    :: attempted,attemp_half
integer (kind=4) :: seed
integer (kind=4) :: Lstag,Nlev
integer (kind=4) :: k,j
integer (kind=4) :: ip
integer (kind=4) :: Nstep,istep
integer (kind=4) :: iblock,Nblock
integer (kind=4) :: istag,Nstag
integer (kind=4) :: ngr
integer (kind=4) :: Nobdm,iobdm
integer (kind=4) :: Nk
integer (kind=4) :: acc_bd,acc_cm,acc_head,acc_tail
integer (kind=4) :: acc_bd_half,acc_cm_half,acc_head_half,acc_tail_half
integer (kind=4) :: acc_open,try_open
integer (kind=4) :: acc_close,try_close
integer (kind=4) :: iworm,iupdate
integer (kind=4) :: idiag,idiag_block,idiag_aux
integer (kind=4) :: obdm_bl,diag_bl

character (len=3) :: sampling

real (kind=8),dimension(:),allocatable     :: LogWF
real (kind=8),dimension(:,:),allocatable   :: xend
real (kind=8),dimension(:,:,:),allocatable :: Path
real (kind=8),dimension(:,:),allocatable   :: nrho,AvNr,AvNr2,VarNr
real (kind=8),dimension(:,:),allocatable   :: Sk,AvSk,AvSk2,VarSk
real (kind=8),dimension(:),allocatable     :: gr,AvGr,AvGr2,VarGr

!Reading input parameters

call ReadParameters(resume,crystal,wf_table,sampling,&
     & density,alpha,dt,a_1,t_0,delta_cm,Rm,dim,Np,Nb,seed,&
     & Lstag,Nlev,Nstag,Nmax,Nobdm,Nblock,Nstep,Nbin,Nk)

pi = acos(-1.d0)
V0 = (sin(alpha))**2

allocate (Lbox(dim),LboxHalf(dim),qbin(dim))

if (crystal) then

   open (unit=2,file='config_ini.in',status='old')

   read (2,*) Np
   read (2,*) (Lbox(k),k=1,dim)
   read (2,*) density

   close (unit=2)
   
else

   do k=1,dim
      Lbox(k) = (real(Np)/density)**(1.d0/real(dim))
   end do

end if  

do k=1,dim
   LboxHalf(k) = 0.5d0*Lbox(k)
   qbin(k)     = 2.d0*pi/Lbox(k)
end do

rcut     = minval(LboxHalf)
rcut2    = rcut*rcut
rbin     = rcut/real(Nbin)
delta_cm = delta_cm/density**(1.d0/real(dim))

!Definition of the parameters of the propagator

t_1 = 0.5d0-t_0
u_0 = (1.d0-1.d0/(1.d0-2.d0*t_0)+1.d0/(6.d0*(1.d0-2.d0*t_0)**3))/12.d0
v_1 = 1.d0/(6.d0*(1.d0-2.d0*t_0)**2)
v_2 = 1.d0-2.d0*v_1

!Generate an initial Path

allocate (Path(dim,Np,0:2*Nb))
allocate (LogWF(0:Nmax+1))
allocate (xend(dim,2))

call init(seed,Path,xend,crystal,resume)
call JastrowTable(rcut,Rm,LogWF)

!Definition of used formats

103 format (x,a,x,i5)
104 format (x,a,x,G13.6e2)
105 format (x,a,x,3G13.6e2)

!Printing the simulation parameters

print *,    ''
print *,    '=============================================================='
print *,    '       VPI Monte Carlo for homogeneous 2D dipoles             '
print *,    '=============================================================='
print *,    ''
print *,    ' '
if (sampling=="sta") then
   print *, '# The Monte Carlo sampling will be performed using STAGING'
   print *, '  algorithm'
else
   print *, '# The Monte Carlo sampling will be performed using BISECTION'
   print *, '  algorithm'
end if
print *,    ' '
print *,    '# Simulation parameters:'
print *,    ''
print 103,  '  > Dimensions          :',dim
print 103,  '  > Number of particles :',Np
print 104,  '  > Density             :',density
print 104,  '  > Polarization        :',alpha
print 105,  '  > Size of the box     :',Lbox
print 103,  '  > Number of beads     :',Nb
print 104,  '  > Time step           :',dt
print 103,  '  > Number of blocks    :',Nblock
print 103,  '  > MC steps per block  :',Nstep
print *,    ''

!Initializing accumulators

allocate (gr(Nbin),AvGr(Nbin),AvGr2(Nbin),VarGr(Nbin))
allocate (nrho(0:Npw,Nbin),AvNr(0:Npw,Nbin),AvNr2(0:Npw,Nbin),VarNr(0:Npw,Nbin))
allocate (Sk(dim,Nk),AvSk(dim,Nk),AvSk2(dim,Nk),VarSk(dim,Nk))

open (unit=2,file='e_vpi.out')

AvE  = 0.d0
AvK  = 0.d0
AvV  = 0.d0

AvEt  = 0.d0
AvKt  = 0.d0
AvVt  = 0.d0

AvE2 = 0.d0
AvK2 = 0.d0
AvV2 = 0.d0

AvEt2 = 0.d0
AvKt2 = 0.d0
AvVt2 = 0.d0

AvGr  = 0.d0
AvGr2 = 0.d0
VarGr = 0.d0

AvSk  = 0.d0
AvSk2 = 0.d0
VarSk = 0.d0

AvNr  = 0.d0
AvNr2 = 0.d0
VarNr = 0.d0
 
nrho = 0.d0

isopen    = .false.
iworm     = 0
idiag     = 0
idiag_aux = 0
obdm_bl   = 0
diag_bl   = 0

!Begin the main Monte Carlo loop

do iblock=1,Nblock
   
   call cpu_time(begin)

   !Initializing block accumulators
   
   try_open  = 0
   acc_open  = 0

   try_close = 0
   acc_close = 0

   attempted = 0
   acc_cm    = 0

   stag_move = 0
   acc_bd    = 0
   acc_head  = 0
   acc_tail  = 0
   
   attemp_half = 0
   acc_cm_half = 0
   
   stag_half     = 0
   acc_bd_half   = 0
   acc_head_half = 0
   acc_tail_half = 0

   idiag_block   = 0 

   ngr = 0
   gr  = 0.d0
   Sk  = 0.d0
   
   BlockAvE = 0.d0
   BlockAvK = 0.d0
   BlockAvV = 0.d0

   BlockAvEt = 0.d0
   BlockAvKt = 0.d0
   BlockAvVt = 0.d0

   BlockAvE2 = 0.d0
   BlockAvK2 = 0.d0
   BlockAvV2 = 0.d0
   
   BlockAvEt2 = 0.d0
   BlockAvKt2 = 0.d0
   BlockAvVt2 = 0.d0
   
   do istep=1,Nstep

      !Open and Close updates... let's see what happens...

      iupdate = int(grnd()*2)

      if (isopen) then
         if (iupdate == 0) then
            call CloseChain(LogWF,density,dt,Lstag,iworm,Path,xend,isopen,acc_close)
            try_close = try_close+1
         end if
      else
         if (iupdate == 1) then
            iworm = int(grnd()*Np)+1
            call OpenChain(LogWF,density,dt,Lstag,iworm,Path,xend,isopen,acc_open)
            try_open = try_open+1
         end if
      end if

      if (isopen) then

         do ip=1,Np

            if (ip/=iworm) then

               attempted = attempted+1
            
               call TranslateChain(delta_cm,LogWF,dt,ip,Path,acc_cm)
               
               do istag=1,Nstag
                  stag_move = stag_move+1
                  call MoveHead(LogWF,dt,Lstag,ip,Path,acc_head)
                  call MoveTail(LogWF,dt,Lstag,ip,Path,acc_tail)
                  call Staging(LogWF,dt,Lstag,ip,Path,acc_bd)
               end do

            else

               do iobdm=1,Nobdm
               
                  do j=1,2

                     attemp_half = attemp_half+1
               
                     call TranslateHalfChain(j,delta_cm,LogWF,dt,ip,Path,xend,acc_cm_half)
                     
                     do istag=1,Nstag

                        stag_half = stag_half+1
                               
                        call MoveHeadHalfChain(j,LogWF,dt,Lstag,ip,Path,xend,acc_head_half)
                        call MoveTailHalfChain(j,LogWF,dt,Lstag,ip,Path,xend,acc_tail_half)
                        call StagingHalfChain(j,LogWF,dt,Lstag,ip,Path,xend,acc_bd_half)
                        
                     end do

                  end do

                  call OBDM(xend,nrho)

               end do

            end if
            
         end do

      else
         
         idiag       = idiag+1
         idiag_aux   = idiag_aux+1
         idiag_block = idiag_block+1
      
         do ip=1,Np

            attempted = attempted+1

            call TranslateChain(delta_cm,LogWF,dt,ip,Path,acc_cm)
            
            do istag=1,Nstag

               stag_move = stag_move+1

               call MoveHead(LogWF,dt,Lstag,ip,Path,acc_head)
               call MoveTail(LogWF,dt,Lstag,ip,Path,acc_tail)
               call Staging(LogWF,dt,Lstag,ip,Path,acc_bd)
            
            end do
               
         end do

         !Energy calculation using mixed estimator
      
         call LocalEnergy(LogWF,Path(:,:,0),Rm,E1,Kin,Pot)
         call LocalEnergy(LogWF,Path(:,:,2*Nb),Rm,E2,Kin,Pot)

         E = 0.5d0*(E1+E2)

         !call PotentialEnergy(Path(:,:,Nb),Pot)
      
         !Energy evaluation using thermodynamic estimator

         call ThermEnergy(Path,dt,Et,Kt,Pot)

         Kin = E-Pot         
            
         !Accumulating energy averages     
      
         call Accumulate(E,Kin,Pot,BlockAvE,BlockAvK,BlockAvV)
         call Accumulate(Et,Kt,Pot,BlockAvEt,BlockAvKt,BlockAvVt)

         call Accumulate(E**2,Kin**2,Pot**2,BlockAvE2,BlockAvK2,BlockAvV2)
         call Accumulate(Et**2,Kt**2,Pot**2,BlockAvEt2,BlockAvKt2,BlockAvVt2)

         !Structural quantities
      
         ngr = ngr+1
      
         call PairCorrelation(Path(:,:,Nb),gr)
         call StructureFactor(Nk,Path(:,:,Nb),Sk)
         
      end if

   end do

   if (idiag_block/=0) then

      diag_bl = diag_bl+1
      
      call NormalizeGr(density,ngr,gr)
      call NormalizeSk(Nk,ngr,Sk)
      
      call AccumGr(gr,AvGr,AvGr2)
      call AccumSk(Nk,Sk,AvSk,AvSk2)

      !Normalizing averages and evaluating variances per block
   
      call NormalizeAv(idiag_block,BlockAvE,BlockAvK,BlockAvV)
      call NormalizeAv(idiag_block,BlockAvE2,BlockAvK2,BlockAvV2)

      BlockVarE = Var(idiag_block,BlockAvE,BlockAvE2)
      BlockVarK = Var(idiag_block,BlockAvK,BlockAvK2)
      BlockVarV = Var(idiag_block,BlockAvV,BlockAvV2)

      call NormalizeAv(idiag_block,BlockAvEt,BlockAvKt,BlockAvVt)
      call NormalizeAv(idiag_block,BlockAvEt2,BlockAvKt2,BlockAvVt2)

      BlockVarEt = Var(idiag_block,BlockAvEt,BlockAvEt2)
      BlockVarKt = Var(idiag_block,BlockAvKt,BlockAvKt2)
      BlockVarVt = Var(idiag_block,BlockAvVt,BlockAvVt2)

      !Accumulating global averages
  
      call Accumulate(BlockAvE,BlockAvK,BlockAvV,AvE,AvK,AvV)
      call Accumulate(BlockAvE**2,BlockAvK**2,BlockAvV**2,AvE2,AvK2,AvV2)

      call Accumulate(BlockAvEt,BlockAvKt,BlockAvVt,AvEt,AvKt,AvVt)
      call Accumulate(BlockAvEt**2,BlockAvKt**2,BlockAvVt**2,AvEt2,AvKt2,AvVt2)

      !Outputs of the block

      write (2,'(5g20.10e3)') real(iblock),BlockAvE/Np,BlockAvK/Np,BlockAvV/Np
      write (99,'(5g20.10e3)') real(iblock),BlockAvEt/Np,BlockAvKt/Np,BlockAvVt/Np

   end if

   if (idiag_aux/Nstep>=1) then
      
      obdm_bl    = obdm_bl+1
      numz_block = real(idiag_aux)
      
      call NormalizeNr(density,numz_block,Nobdm,nrho)
      call AccumNr(nrho,AvNr,AvNr2)
      
      !Restarting counters
   
      idiag_aux = 0
      nrho      = 0.d0
   
   end if
     
   if (mod(iblock,10)==0) then

      call CheckPoint(Path,xend)

   end if

   call cpu_time(end)

101 format (x,a,x,f6.2,x,a)
102 format (a,x,G16.8e2,x,a,x,G16.8e2)

   print *,   '-----------------------------------------------------------'
   print *,   'BLOCK NUMBER :',iblock
   print *,   ' '
   print *,   '# Block results:'
   print *,   ' '
   print 102, '  > <E>  =',BlockAvE/Np,'+/-',BlockVarE/Np
   print 102, '  > <Ec> =',BlockAvK/Np,'+/-',BlockVarK/Np
   print 102, '  > <Ep> =',BlockAvV/Np,'+/-',BlockVarV/Np
   print *,   ' '
   print 102, '  > <Et> =',BlockAvEt/Np,'+/-',BlockVarEt/Np
   print 102, '  > <Kt> =',BlockAvKt/Np,'+/-',BlockVarKt/Np
   print 102, '  > <Vt> =',BlockAvVt/Np,'+/-',BlockVarVt/Np
   print *,   ''
   print *,   '# Acceptance of diagonal movements:'
   print *,   ' '
   print 101, '> CM movements      =',100*real(acc_cm)/attempted,'%'
   print 101, '> Staging movements =',100*real(acc_bd)/stag_move,'%'
   print 101, '> Head movements    =',100*real(acc_head)/stag_move,'%'
   print 101, '> Tail movements    =',100*real(acc_tail)/stag_move,'%'
   print *,   ' '
   print *,   '# Acceptance of off-diagonal movements:'
   print *,   ' '
   print 101, '> CM movements      =',100*real(acc_cm_half)/attemp_half,'%'
   print 101, '> Staging movements =',100*real(acc_bd_half)/stag_half,'%'
   print 101, '> Head movements    =',100*real(acc_head_half)/stag_half,'%'
   print 101, '> Tail movements    =',100*real(acc_tail_half)/stag_half,'%'
   print *,   ' '
   print *,   '# Acceptance open/close updates:'
   print *,   ' '
   print 101, '> Diagonal conf.    =',100.d0*real(idiag_block)/real(Nstep),'%'
   print 101, '> Open acc          =',100.d0*real(acc_open)/real(try_open),'%'
   print 101, '> Close try         =',100.d0*real(acc_close)/real(try_close),'%'
   print *,   ' '
   print 101, '# Time per block    =',end-begin,'seconds'
  
end do

close (unit=2)

deallocate (LogWF)

!Normalizing global averages and evaluating final variances

call NormalizeAv(diag_bl,AvE,AvK,AvV)
call NormalizeAv(diag_bl,AvE2,AvK2,AvV2)

VarE = Var(diag_bl,AvE,AvE2)
VarK = Var(diag_bl,AvK,AvK2)
VarV = Var(diag_bl,AvV,AvV2)

call NormalizeAv(diag_bl,AvEt,AvKt,AvVt)
call NormalizeAv(diag_bl,AvEt2,AvKt2,AvVt2)

VarEt = Var(diag_bl,AvEt,AvEt2)
VarKt = Var(diag_bl,AvKt,AvKt2)
VarVt = Var(diag_bl,AvVt,AvVt2)

print *, '=============================================================='
print *, 'FINAL RESULTS:'
print *, ''
print *, '# Final averages:'
print *, ''
print 102, '  > <E>  =',AvE/Np,'+/-',VarE/Np
print 102, '  > <Ec> =',AvK/Np,'+/-',VarK/Np
print 102, '  > <Ep> =',AvV/Np,'+/-',VarV/Np
print *, ''
print 102, '  > <Et> =',AvEt/Np,'+/-',VarEt/Np
print 102, '  > <Kt> =',AvKt/Np,'+/-',VarKt/Np
print 102, '  > <Vt> =',AvVt/Np,'+/-',VarVt/Np
print *, ''
print *, '=============================================================='
print *, ''

call NormAvGr(diag_bl,AvGr,AvGr2,VarGr)
call NormAvSk(diag_bl,Nk,AvSk,AvSk2,VarSk)
call NormAvNr(obdm_bl,AvNr,AvNr2,VarNr)

deallocate (Path)
deallocate (nrho,AvNr,AvNr2,VarNr)
deallocate (gr,AvGr,AvGr2,VarGr)
deallocate (Sk,AvSk,AvSk2,VarSk)

end program vpi

!-----------------------------------------------------------------------
