program vpi

use random_mod
use global_mod
use system_mod
use sample_mod
use vpi_mod

implicit none 

logical          :: crystal,resume
logical          :: isopen,swapping
logical          :: trap
real (kind=8)    :: r8_gamma
real (kind=8)    :: dt,delta_cm
real (kind=8)    :: density
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
real (kind=8)    :: try_stag,try_stag_half
real (kind=8)    :: try_cm,try_cm_half
integer (kind=4) :: seed
integer (kind=4) :: Lstag,Nlev
integer (kind=4) :: k,j
integer (kind=4) :: ip
integer (kind=4) :: Nstep,istep
integer (kind=4) :: iblock,Nblock
integer (kind=4) :: istag,Nstag
integer (kind=4) :: CMFreq
integer (kind=4) :: ngr
integer (kind=4) :: Nobdm,iobdm
integer (kind=4) :: Nk
integer (kind=4) :: acc_bd,acc_cm,acc_head,acc_tail
integer (kind=4) :: acc_bd_half,acc_cm_half,acc_head_half,acc_tail_half
integer (kind=4) :: acc_open,try_open
integer (kind=4) :: acc_close,try_close
integer (kind=4) :: acc_swap,try_swap
integer (kind=4) :: iworm,iupdate
integer (kind=4) :: idiag,idiag_block,idiag_aux
integer (kind=4) :: obdm_bl,diag_bl

character (len=3) :: sampling

real (kind=8),dimension(:),allocatable     :: LogWF,VTable
real (kind=8),dimension(:,:),allocatable   :: xend
real (kind=8),dimension(:,:,:),allocatable :: Path
real (kind=8),dimension(:,:),allocatable   :: nrho,AvNr,AvNr2,VarNr
real (kind=8),dimension(:,:),allocatable   :: Sk,AvSk,AvSk2,VarSk
real (kind=8),dimension(:),allocatable     :: gr,AvGr,AvGr2,VarGr
real (kind=8),dimension(:,:),allocatable   :: dens


logical                                    :: new_perm_cycle=.false.
logical                                    :: end_perm_cycle=.false.
logical                                    :: swap_accepted
integer (kind=4)                           :: iperm,ik
integer (kind=4),dimension(:),allocatable  :: Particles_in_perm_cycle
integer (kind=4),dimension(:),allocatable  :: Perm_histogram

!Reading input parameters

call ReadParameters(resume,crystal,wf_table,v_table,swapping,trap,&
     & sampling,density,dt,delta_cm,Rm,dim,Np,Nb,seed,&
     & CMFreq,Lstag,Nlev,Nstag,Nmax,Nobdm,Nblock,Nstep,Nbin,Nk)

pi = acos(-1.d0)

if (trap) then
   
   rcut = 1.d0
   
   do k=1,dim
      rcut=3.d0*rcut*a_ho(k)
   end do

   density  = real(Np)/(pi**(0.5d0*dim)*rcut/r8_gamma(0.5d0*dim+1.d0)) 
   rcut     = rcut**(1.d0/real(dim))
   rcut     = 10.d0*rcut
   delta_cm = delta_cm*minval(a_ho)
   
else

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
   delta_cm = delta_cm/density**(1.d0/real(dim))

end if

rcut2    = rcut*rcut
rbin     = rcut/real(Nbin)
isopen   = .false.
iworm    = 0

!Generate an initial Path

allocate (Path(dim,Np,0:2*Nb))
!allocate (LogWF(0:Nmax+1),VTable(0:Nmax+1))
allocate (xend(dim,2))

if (swapping) then
   allocate (Particles_in_perm_cycle(Np))
   allocate (Perm_histogram(Np))
   Perm_histogram(:) = 0
end if

call init(trap,seed,Path,xend,crystal,resume,isopen,iworm)

if (wf_table) then
   allocate (LogWF(0:Nmax+1))
   call JastrowTable(rcut,LogWF)
end if
if (v_table) then
   allocate (VTable(0:Nmax+1))
   call PotentialTable(rcut,VTable)
end if

!Printing the simulation parameters

103 format (x,a,x,i5)
104 format (x,a,x,G13.6e2)
105 format (x,a,x,3G13.6e2)

print *,    ''
print *,    '=============================================================='
print *,    '                      VPI Monte Carlo                         '
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
if (swapping) then
   print *, '# The Monte Carlo sampling will use swap updates'
else
   print *, '# The Monte Carlo sampling will not use swap updates'
end if
print *,    ' '
print *,    '# Simulation parameters:'
print *,    ''
print 103,  '  > Dimensions          :',dim
print 103,  '  > Number of particles :',Np
if (trap) then
   print 105,  '  > Trapping length     :',a_ho
else 
   print 104,  '  > Density             :',density
   print 105,  '  > Size of the box     :',Lbox
end if
print 103,  '  > Number of beads     :',Nb
print 104,  '  > Time step           :',dt
print 103,  '  > Number of blocks    :',Nblock
print 103,  '  > MC steps per block  :',Nstep
print *,    ''

!Initializing accumulators

allocate (gr(Nbin),AvGr(Nbin),AvGr2(Nbin),VarGr(Nbin),dens(Nbin,Nbin))
allocate (nrho(0:Npw,Nbin),AvNr(0:Npw,Nbin),AvNr2(0:Npw,Nbin),VarNr(0:Npw,Nbin))
allocate (Sk(dim,Nk),AvSk(dim,Nk),AvSk2(dim,Nk),VarSk(dim,Nk))

open (unit=1,file='e_vpi.out')
open (unit=2,file='et_vpi.out')

AvE  = 0.d0
AvK  = 0.d0
AvV  = 0.d0

AvEt = 0.d0
AvKt = 0.d0
AvVt = 0.d0

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

idiag     = 0
idiag_aux = 0
obdm_bl   = 0
diag_bl   = 0

dens = 0.d0

!Begin the main Monte Carlo loop

do iblock=1,Nblock
   
   !Initializing block accumulators

   call cpu_time(begin)
   
   try_open  = 0
   acc_open  = 0

   try_close = 0
   acc_close = 0

   try_cm = 0
   acc_cm = 0

   try_stag = 0
   acc_bd   = 0
   acc_head = 0
   acc_tail = 0
   
   try_cm_half = 0
   acc_cm_half = 0
   
   try_stag_half = 0
   acc_bd_half   = 0
   acc_head_half = 0
   acc_tail_half = 0

   try_swap = 0
   acc_swap = 0

   idiag_block = 0 

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

      !Open and Close updates to determine if the system is in a diagonal
      !or in an off-diagonal configuration

      iupdate = int(grnd()*2)

      if (isopen) then
         if (iupdate == 0) then
            call CloseChain(trap,LogWF,VTable,density,dt,Lstag,iworm,Path,xend,&
                           &isopen,acc_close,end_perm_cycle)
            try_close = try_close+1
            call PermutationSampling(new_perm_cycle,end_perm_cycle,iperm,&
                                    &Particles_in_perm_cycle,Perm_histogram,&
                                    &isopen,iworm)
         end if
      else
         if (iupdate == 1) then
            iworm = int(grnd()*Np)+1
            call OpenChain(trap,LogWF,VTable,density,dt,Lstag,iworm,Path,xend,&
                          &isopen,acc_open,new_perm_cycle)
            try_open = try_open+1
            call PermutationSampling(new_perm_cycle,end_perm_cycle,iperm,&
                                    &Particles_in_perm_cycle,Perm_histogram,&
                                    &isopen,iworm)
         end if
      end if

      !Sampling of the paths 

      if (isopen) then

         if (mod(istep,CMFreq)==0) then

            do ip=1,Np

               if (ip/=iworm) then
               
                  try_cm = try_cm+1
                  call TranslateChain(trap,delta_cm,LogWF,VTable,dt,ip,Path,acc_cm)
               
               end if

            end do

         end if

         do istag=1,Nstag

            do ip=1,Np

               if (ip/=iworm) then
                  
                  try_stag = try_stag+1

                  if (sampling=="sta") then
                     call MoveHead(trap,LogWF,VTable,dt,Lstag,ip,Path,acc_head)
                     call MoveTail(trap,LogWF,VTable,dt,Lstag,ip,Path,acc_tail)
                     call Staging(trap,LogWF,VTable,dt,Lstag,ip,Path,acc_bd)
                  else
                     call MoveHeadBisection(trap,LogWF,VTable,dt,Nlev,ip,Path,acc_head)
                     call MoveTailBisection(trap,LogWF,VTable,dt,Nlev,ip,Path,acc_tail)
                     call Bisection(trap,LogWF,VTable,dt,Nlev,ip,Path,acc_bd)
                  end if

               end if

            end do

         end do

         ! Worm movements

         do iobdm=1,Nobdm

            ip = iworm

            do j=1,2
               try_cm_half = try_cm_half+1
               call TranslateHalfChain(trap,j,delta_cm,LogWF,VTable,dt,ip,Path,xend,acc_cm_half)
            end do

            do j=1,2
                  
               try_stag_half = try_stag_half+1
                        
               call MoveHeadHalfChain(trap,j,LogWF,VTable,dt,Lstag,ip,Path,xend,acc_head_half)
               call MoveTailHalfChain(trap,j,LogWF,VTable,dt,Lstag,ip,Path,xend,acc_tail_half)
               call StagingHalfChain(trap,j,LogWF,VTable,dt,Lstag,ip,Path,xend,acc_bd_half)
                        
            end do

            if (swapping) then
                        
               try_swap = try_swap+1
               
               call Swap(trap,LogWF,VTable,dt,Lstag,ip,Path,xend,acc_swap,ik,swap_accepted)
               call PermutationSampling(new_perm_cycle,end_perm_cycle,iperm,&
                                       &Particles_in_perm_cycle,Perm_histogram,&
                                       &isopen,iworm,ik,swap_accepted)
                        
            end if

            if (trap .eqv. .false.) then
               call OBDM(xend,nrho)
            end if

         end do

      else
         
         idiag       = idiag+1
         idiag_aux   = idiag_aux+1
         idiag_block = idiag_block+1

         if (mod(istep,CMFreq)==0) then
            
            do ip=1,Np
               try_cm = try_cm+1
               call TranslateChain(trap,delta_cm,LogWF,VTable,dt,ip,Path,acc_cm)
            end do

         end if
         
         do istag=1,Nstag
            
            do ip=1,Np
               
               try_stag = try_stag+1

               if (sampling=="sta") then
                  call MoveHead(trap,LogWF,VTable,dt,Lstag,ip,Path,acc_head)
                  call MoveTail(trap,LogWF,VTable,dt,Lstag,ip,Path,acc_tail)
                  call Staging(trap,LogWF,VTable,dt,Lstag,ip,Path,acc_bd)
               else
                  call MoveHeadBisection(trap,LogWF,VTable,dt,Nlev,ip,Path,acc_head)
                  call MoveTailBisection(trap,LogWF,VTable,dt,Nlev,ip,Path,acc_tail)
                  call Bisection(trap,LogWF,VTable,dt,Nlev,ip,Path,acc_bd)
               end if

            end do

         end do

         !Energy calculation using mixed estimator
      
         call LocalEnergy(trap,LogWF,VTable,Path(:,:,0),E1,Kin,Pot)
         call LocalEnergy(trap,LogWF,VTable,Path(:,:,2*Nb),E2,Kin,Pot)

         E = 0.5d0*(E1+E2)

         !Energy evaluation using thermodynamic estimator

         call ThermEnergy(trap,VTable,Path,dt,Et,Kt,Pot)

         Kin = E-Pot         
            
         !Accumulating energy averages     
      
         call Accumulate(E,Kin,Pot,BlockAvE,BlockAvK,BlockAvV)
         call Accumulate(Et,Kt,Pot,BlockAvEt,BlockAvKt,BlockAvVt)

         call Accumulate(E**2,Kin**2,Pot**2,BlockAvE2,BlockAvK2,BlockAvV2)
         call Accumulate(Et**2,Kt**2,Pot**2,BlockAvEt2,BlockAvKt2,BlockAvVt2)

         !Structural quantities
      
         ngr = ngr+1
      
         if (trap .eqv. .false.) then
            call PairCorrelation(Path(:,:,Nb),gr)
            call StructureFactor(Nk,Path(:,:,Nb),Sk)
         end if
         
         !call DensityProfile(Path(:,:,Nb),dens)
         
      end if

   end do

   if (idiag_block/=0) then

      !Normalize final results of the block
      
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

      !Update all the accumulators and observables

      diag_bl = diag_bl+1
      
      call Accumulate(BlockAvE,BlockAvK,BlockAvV,AvE,AvK,AvV)
      call Accumulate(BlockAvE**2,BlockAvK**2,BlockAvV**2,AvE2,AvK2,AvV2)

      call Accumulate(BlockAvEt,BlockAvKt,BlockAvVt,AvEt,AvKt,AvVt)
      call Accumulate(BlockAvEt**2,BlockAvKt**2,BlockAvVt**2,AvEt2,AvKt2,AvVt2)
      
      if (trap .eqv. .false.) then
      
         call NormalizeGr(density,ngr,gr)
         call NormalizeSk(Nk,ngr,Sk)

         call AccumGr(gr,AvGr,AvGr2)
         call AccumSk(Nk,Sk,AvSk,AvSk2)      

      end if

      !Outputs of the block

      write (1,'(5g20.10e3)') real(iblock),BlockAvE/Np,BlockAvK/Np,BlockAvV/Np
      write (2,'(5g20.10e3)') real(iblock),BlockAvEt/Np,BlockAvKt/Np,BlockAvVt/Np

   end if

   if (idiag_aux/Nstep>=1) then
      
      obdm_bl    = obdm_bl+1
      numz_block = real(idiag_aux)
      
      if (trap .eqv. .false.) then

         call NormalizeNr(density,numz_block,Nobdm,nrho)
         call AccumNr(nrho,AvNr,AvNr2)

      end if
      
      !Restarting counters
   
      idiag_aux = 0
      nrho      = 0.d0
   
   end if
     
   if (mod(iblock,1)==0) then

      call CheckPoint(trap,Path,xend,isopen,iworm)

   end if

   call cpu_time(end)

101 format (x,a,x,f7.2,x,a)
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
   print 101, '> CM movements      =',100*real(acc_cm)/try_cm,'%'
   print 101, '> Staging movements =',100*real(acc_bd)/try_stag,'%'
   print 101, '> Head movements    =',100*real(acc_head)/try_stag,'%'
   print 101, '> Tail movements    =',100*real(acc_tail)/try_stag,'%'
   print *,   ' '
   print *,   '# Acceptance of off-diagonal movements:'
   print *,   ' '
   print 101, '> CM movements      =',100*real(acc_cm_half)/try_cm_half,'%'
   print 101, '> Staging movements =',100*real(acc_bd_half)/try_stag_half,'%'
   print 101, '> Head movements    =',100*real(acc_head_half)/try_stag_half,'%'
   print 101, '> Tail movements    =',100*real(acc_tail_half)/try_stag_half,'%'
   print *,   ' '
   print *,   '# Acceptance open/close updates:'
   print *,   ' '
   print 101, '> Diagonal conf.    =',100.d0*real(idiag_block)/real(Nstep),'%'
   print 101, '> Open acc          =',100.d0*real(acc_open)/real(try_open),'%'
   print 101, '> Close acc         =',100.d0*real(acc_close)/real(try_close),'%'
   print 101, '> Swap acc          =',100.d0*real(acc_swap)/real(try_swap),'%'
   print 101, ' '
   print 101, '# Time per block    =',end-begin,'seconds'

end do

do ip=1,Np
   write (99,*) ip,Perm_histogram(ip)
end do

close (unit=1)
close (unit=2)

if (wf_table) then
   deallocate (LogWF)
end if
if (v_table) then
   deallocate (VTable)
end if

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

if (trap .eqv. .false.) then

   call NormAvGr(diag_bl,AvGr,AvGr2,VarGr)
   call NormAvSk(diag_bl,Nk,AvSk,AvSk2,VarSk)
   call NormAvNr(obdm_bl,AvNr,AvNr2,VarNr)

end if

if (swapping) then
   deallocate (Particles_in_perm_cycle,Perm_histogram)
end if

deallocate (Path)
deallocate (nrho,AvNr,AvNr2,VarNr)
deallocate (gr,AvGr,AvGr2,VarGr)
deallocate (Sk,AvSk,AvSk2,VarSk)

end program vpi

!-----------------------------------------------------------------------
