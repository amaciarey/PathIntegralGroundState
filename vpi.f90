program vpi

use random_mod
use global_mod
use system_mod
use sample_mod
use vpi_mod

implicit none 

logical          :: crystal,diagonal,resume
real (kind=8)    :: dt
real (kind=8)    :: density,alpha
real (kind=8)    :: delta_cm
real (kind=8)    :: E1,E2
real (kind=8)    :: E,Kin,Pot
real (kind=8)    :: BlockAvE,BlockAvK,BlockAvV
real (kind=8)    :: BlockAvE2,BlockAvK2,BlockAvV2
real (kind=8)    :: BlockVarE,BlockVarK,BlockVarV
real (kind=8)    :: AvE,AvK,AvV
real (kind=8)    :: AvE2,AvK2,AvV2
real (kind=8)    :: VarE,VarK,VarV
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

character (len=3) :: sampling

real (kind=8),dimension(:),allocatable     :: LogWF
real (kind=8),dimension(:,:),allocatable   :: xend
real (kind=8),dimension(:,:,:),allocatable :: Path
real (kind=8),dimension(:,:),allocatable   :: nrho,AvNr,AvNr2,VarNr
real (kind=8),dimension(:,:),allocatable   :: Sk,AvSk,AvSk2,VarSk
real (kind=8),dimension(:),allocatable     :: gr,AvGr,AvGr2,VarGr

!Reading input parameters

call ReadParameters(resume,crystal,diagonal,wf_table,sampling,&
     & density,alpha,dt,a_1,t_0,delta_cm,Rm,Ak,N0,dim,Np,Nb,seed,&
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

rcut        = minval(LboxHalf)
rcut2       = rcut*rcut
rbin        = rcut/real(Nbin)
delta_cm    = delta_cm/density**(1.d0/real(dim))
attempted   = real(Nstep*Np)
stag_move   = real(attempted*Nstag)
attemp_half = 2.d0*real(Nstep*Nobdm)
stag_half   = 2.d0*real(Nstep*Nobdm*Nstag)

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

print *, ''
print *, '=============================================================='
print *, '       VPI Monte Carlo for homogeneous 2D dipoles             '
print *, '=============================================================='
print *, ''
if (diagonal) then
   print *, '# The simulation will consider DIAGONAL configurations only'
else 
   print *, '# The simulation will consider OFF-DIAGONAL configurations'
end if
print *, ' '
if (sampling=="sta") then
   print *, '# The Monte Carlo sampling will be performed using STAGING'
   print *, '  algorithm'
else
   print *, '# The Monte Carlo sampling will be performed using BISECTION'
   print *, '  algorithm'
end if
print *, ' '
print *, '# Simulation parameters:'
print *, ''
print 103, '  > Dimensions          :',dim
print 103, '  > Number of particles :',Np
print 104, '  > Density             :',density
print 104, '  > Polarization        :',alpha
print 105, '  > Size of the box     :',Lbox
print 103, '  > Number of beads     :',Nb
print 104, '  > Time step           :',dt
print 103, '  > Number of blocks    :',Nblock
print 103, '  > MC steps per block  :',Nstep
print *, ''

!Initializing accumulators

allocate (gr(Nbin),AvGr(Nbin),AvGr2(Nbin),VarGr(Nbin))
allocate (nrho(0:Npw,Nbin),AvNr(0:Npw,Nbin),AvNr2(0:Npw,Nbin),VarNr(0:Npw,Nbin))
allocate (Sk(dim,Nk),AvSk(dim,Nk),AvSk2(dim,Nk),VarSk(dim,Nk))

open (unit=2,file='e_vpi.out')

AvE  = 0.d0
AvK  = 0.d0
AvV  = 0.d0

AvE2 = 0.d0
Avk2 = 0.d0
AvV2 = 0.d0

ngr = 0

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

!Begin the main Monte Carlo loop

do iblock=1,Nblock
   
   call cpu_time(begin)

   BlockAvE = 0.d0
   BlockAvK = 0.d0
   BlockAvV = 0.d0

   BlockAvE2 = 0.d0
   BlockAvK2 = 0.d0
   BlockAvV2 = 0.d0

   acc_cm   = 0
   acc_bd   = 0
   acc_head = 0
   acc_tail = 0

   acc_cm_half   = 0
   acc_bd_half   = 0
   acc_head_half = 0
   acc_tail_half = 0

   ngr  = 0
   gr   = 0.d0
   Sk   = 0.d0
   nrho = 0.d0
   
   do istep=1,Nstep

      if (diagonal) then

         do ip=1,Np
            
            if (mod(istep,1)==0) then
               call TranslateChain(delta_cm,LogWF,dt,ip,Path,acc_cm)
            end if

            do istag=1,Nstag
               
               if (sampling=="sta") then                     
                  call MoveHead(LogWF,dt,Lstag,ip,Path,acc_head)
                  call MoveTail(LogWF,dt,Lstag,ip,Path,acc_tail)
                  call Staging(LogWF,dt,Lstag,ip,Path,acc_bd)
               else
                  call Bisection(LogWF,dt,Nlev,ip,Path,acc_bd)
                  call MoveHeadBisection(LogWF,dt,Nlev,ip,Path,acc_head)
                  call MoveTailBisection(LogWF,dt,Nlev,ip,Path,acc_tail)
               end if

            end do

         end do

      else
          
         do ip=1,Np-1

            if (mod(istep,1)==0) then
               call TranslateChain(delta_cm,LogWF,dt,ip,Path,acc_cm)
            end if

            do istag=1,Nstag            
               
               if (sampling=="sta") then
                  call MoveHead(LogWF,dt,Lstag,ip,Path,acc_head)
                  call MoveTail(LogWF,dt,Lstag,ip,Path,acc_tail)
                  call Staging(LogWF,dt,Lstag,ip,Path,acc_bd)
               else
                  call Bisection(LogWF,dt,Nlev,ip,Path,acc_bd)
                  call MoveHeadBisection(LogWF,dt,Nlev,ip,Path,acc_head)
                  call MoveTailBisection(LogWF,dt,Nlev,ip,Path,acc_tail)
               end if

            end do

         end do

         !Move the last chain that represents the off-diagonal configuration of the 
         !system

         do iobdm=1,Nobdm
               
            do j=1,2
               
               if (mod(istep,1)==0) then
                  call TranslateHalfChain(j,delta_cm,LogWF,dt,Np,Path,xend,acc_cm_half)
               end if

               do istag=1,Nstag
                               
                  if (sampling=="sta") then
                     call MoveHeadHalfChain(j,LogWF,dt,Lstag,Np,Path,xend,acc_head_half)
                     call MoveTailHalfChain(j,LogWF,dt,Lstag,Np,Path,xend,acc_tail_half)
                     call StagingHalfChain(j,LogWF,dt,Lstag,Np,Path,xend,acc_bd_half)
                  else
                     call MoveHeadHalfBisection(j,LogWF,dt,Nlev,Np,Path,xend,acc_head_half)
                     call MoveTailHalfBisection(j,LogWF,dt,Nlev,Np,Path,xend,acc_tail_half)
                     call BisectionHalf(j,LogWF,dt,Nlev,Np,Path,xend,acc_bd_half)
                  end if

               end do

            end do

            call OBDM(xend,nrho)

         end do

      end if

      !Energy calculation using mixed estimator
      
      call LocalEnergy(LogWF,Path(:,:,0),Rm,E1,Kin,Pot)
      call LocalEnergy(LogWF,Path(:,:,2*Nb),Rm,E2,Kin,Pot)

      E = 0.5d0*(E1+E2)

      call PotentialEnergy(Path(:,:,Nb),Pot)

      Kin = E-Pot
            
      !Accumulating energy averages     
      
      call Accumulate(E,Kin,Pot,BlockAvE,BlockAvK,BlockAvV)
      call Accumulate(E**2,Kin**2,Pot**2,BlockAvE2,BlockAvK2,BlockAvV2)

      !Structural quantities
      
      ngr = ngr+1
      
      call PairCorrelation(Path(:,:,Nb),gr)
      call StructureFactor(Nk,Path(:,:,Nb),Sk)

   end do

   !Normalizing averages and evaluating variances per block
   
   call NormalizeAv(Nstep,BlockAvE,BlockAvK,BlockAvV)
   call NormalizeAv(Nstep,BlockAvE2,BlockAvK2,BlockAvV2)

   BlockVarE = Var(Nstep,BlockAvE,BlockAvE2)
   BlockVarK = Var(Nstep,BlockAvK,BlockAvK2)
   BlockVarV = Var(Nstep,BlockAvV,BlockAvV2)
   
   !Accumulating global averages
  
   call Accumulate(BlockAvE,BlockAvK,BlockAvV,AvE,AvK,AvV)
   call Accumulate(BlockAvE**2,BlockAvK**2,BlockAvV**2,AvE2,AvK2,AvV2)

   !Outputs of the block

   write (2,'(5g20.10e3)') real(iblock),BlockAvE/Np,BlockAvK/Np,BlockAvV/Np

   if (mod(iblock,10)==0) then

      call CheckPoint(Path,xend)

   end if

   call Normalize(density,Nk,ngr,gr,Sk,nrho)
   call AccumGr(gr,AvGr,AvGr2)
   call AccumSk(Nk,Sk,AvSk,AvSk2)
   if (diagonal .eqv. .false.) then
      call AccumNr(nrho,AvNr,AvNr2)
   end if

   call cpu_time(end)

101 format (x,a,x,f5.2,x,a)
102 format (a,x,G16.8e2,x,a,x,G16.8e2)

   print *, '-----------------------------------------------------------'
   print *, 'BLOCK NUMBER :',iblock
   print *, ' '
   print *, '# Block results:'
   print *, ' '
   print 102, '  > <E>  =',BlockAvE/Np,'+/-',BlockVarE/Np
   print 102, '  > <Ec> =',BlockAvK/Np,'+/-',BlockVarK/Np
   print 102, '  > <Ep> =',BlockAvV/Np,'+/-',BlockVarV/Np
   print *, ''
   print *, '# Acceptance of diagonal movements:'
   print *, ' '
   print 101, '> CM movements      =',100*real(acc_cm)/attempted,'%'
   print 101, '> Staging movements =',100*real(acc_bd)/stag_move,'%'
   print 101, '> Head movements    =',100*real(acc_head)/stag_move,'%'
   print 101, '> Tail movements    =',100*real(acc_tail)/stag_move,'%'
   if (diagonal .eqv. .false.) then
      print *, ' '
      print *, '# Acceptance of off-diagonal movements:'
      print *, ' '
      print 101, '> CM movements      =',100*real(acc_cm_half)/attemp_half,'%'
      print 101, '> Staging movements =',100*real(acc_bd_half)/stag_half,'%'
      print 101, '> Head movements    =',100*real(acc_head_half)/stag_half,'%'
      print 101, '> Tail movements    =',100*real(acc_tail_half)/stag_half,'%'
   end if
   print *, ' '
   print 101, '# Time per block    =',end-begin,'seconds'
  
end do

close (unit=2)

deallocate (LogWF)

!Normalizing global averages and evaluating final variances

call NormalizeAv(Nblock,AvE,AvK,AvV)
call NormalizeAv(Nblock,AvE2,AvK2,AvV2)

VarE = Var(Nblock,AvE,AvE2)
VarK = Var(Nblock,AvK,AvK2)
VarV = Var(Nblock,AvV,AvV2)

print *, '=============================================================='
print *, 'FINAL RESULTS:'
print *, ''
print *, '# Final averages:'
print *, ''
print 102, '  > <E>  =',AvE/Np,'+/-',VarE/Np
print 102, '  > <Ec> =',AvK/Np,'+/-',VarK/Np
print 102, '  > <Ep> =',AvV/Np,'+/-',VarV/Np
print *, ''
print *, '=============================================================='
print *, ''

call NormAvGr(Nblock,AvGr,AvGr2,VarGr)
call NormAvSk(Nblock,Nk,AvSk,AvSk2,VarSk)
call NormAvNr(Nblock,AvNr,AvNr2,VarNr)

deallocate (Path)
deallocate (nrho)
deallocate (gr,AvGr,AvGr2,VarGr)
deallocate (Sk,AvSk,AvSk2,VarSk)

end program vpi

!-----------------------------------------------------------------------
