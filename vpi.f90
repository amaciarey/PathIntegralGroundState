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
real (kind=8)    :: E,Kin,Pot
real (kind=8)    :: BlockAvE,BlockAvK,BlockAvV
real (kind=8)    :: BlockAvE2,BlockAvK2,BlockAvV2
real (kind=8)    :: BlockVarE,BlockVarK,BlockVarV
real (kind=8)    :: AvE,AvK,AvV
real (kind=8)    :: AvE2,AvK2,AvV2
real (kind=8)    :: VarE,VarK,VarV
real (kind=8)    :: end,begin
integer (kind=4) :: seed
integer (kind=4) :: Lstag
integer (kind=4) :: k,j
integer (kind=4) :: ip
integer (kind=4) :: attempted
integer (kind=4) :: Nstep,istep
integer (kind=4) :: iblock,Nblock
integer (kind=4) :: ibis,Nbis
integer (kind=4) :: ngr
integer (kind=4) :: Nobdm,iobdm
integer (kind=4) :: Nk
integer (kind=4) :: acc_bd,acc_cm,acc_head,acc_tail

real (kind=8),dimension(:),allocatable     :: LogWF
real (kind=8),dimension(:,:),allocatable   :: xend
real (kind=8),dimension(:,:,:),allocatable :: Path
real (kind=8),dimension(:,:),allocatable   :: nrho,AvNr,AvNr2,VarNr
real (kind=8),dimension(:,:),allocatable   :: Sk,AvSk,AvSk2,VarSk
real (kind=8),dimension(:),allocatable     :: gr,AvGr,AvGr2,VarGr

!Reading input parameters

call ReadParameters(resume,crystal,diagonal,wf_table,density,&
       & alpha,dt,a_1,t_0,delta_cm,Rm,Ak,N0,dim,Np,Nb,seed,Lstag,Nmax,&
       & Nobdm,Nblock,Nstep,Nbin,Nk)

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

rcut      = minval(LboxHalf)
rcut2     = rcut*rcut
rbin      = rcut/real(Nbin)
delta_cm  = delta_cm/density**(1.d0/real(dim))
attempted = Nstep*Np
Nbis      = int(real(2*Nb+1)/real(Lstag)+1)

!Definition of the parameters of the propagator

t_1 = 0.5d0-t_0
u_0 = (1.d0-1.d0/(1.d0-2.d0*t_0)+1.d0/(6.d0*(1.d0-2.d0*t_0)**3))/12.d0
v_1 = 1.d0/(6.d0*(1.d0-2.d0*t_0)**2)
v_2 = 1.d0-2.d0*v_1

!Generate an initial Path

allocate (Path(dim,Np,0:2*Nb))
allocate (LogWF(0:Nmax+1))
allocate (xend(dim,2))

call sgrnd(seed)
call init(Path,xend,crystal,resume)
call Jastrow_Table(rcut,Rm,LogWF)

!Printing the simulation parameters

print *, ''
print *, '=============================================================='
print *, '       VPI Monte Carlo for homogeneous 2D dipoles             '
print *, '=============================================================='
print *, ''
print *, '# Simulation parameters:'
print *, ''
print *, '  > Dimensions            :',dim
print *, '  > Number of particles   :',Np
print *, '  > Density               :',density
print *, '  > Polarization          :',alpha
print *, '  > Size of the box       :',Lbox
print *, '  > Number of beads       :',Nb
print *, '  > Time step             :',dt
print *, '  > Number of blocks      :',Nblock
print *, '  > MC steps per block    :',Nstep
print *, ''

!Begin of the Monte Carlo sampling

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

   ngr  = 0
   gr   = 0.d0
   Sk   = 0.d0
   nrho = 0.d0
   
   do istep=1,Nstep

      if (diagonal) then
      
         do ip=1,Np
            
            call TranslateChain(delta_cm,LogWF,dt,ip,Path,acc_cm)
                  
            do ibis=1,Nbis            
                     
               call MoveHead(LogWF,dt,Lstag,ip,Path,acc_head)
               call MoveTail(LogWF,dt,Lstag,ip,Path,acc_tail)
               call Staging(LogWF,dt,Lstag,ip,Path,acc_bd)
               
            end do

         end do

      else
          
         do ip=1,Np-1

            call TranslateChain(delta_cm,LogWF,dt,ip,Path,acc_cm)
                  
            do ibis=1,Nbis            
               
               call MoveHead(LogWF,dt,Lstag,ip,Path,acc_head)
               call MoveTail(LogWF,dt,Lstag,ip,Path,acc_tail)
               call Staging(LogWF,dt,Lstag,ip,Path,acc_bd)
               
            end do

         end do

         !Move the last chain that is broken in two pieces

         do iobdm=1,Nobdm
               
            do j=1,2

               call TranslateHalfChain(j,delta_cm,LogWF,dt,Np,Path,xend)
                  
               do ibis=1,Nbis
                  
                  call MoveHeadHalfChain(j,LogWF,dt,Lstag,Np,Path,xend)
                  call MoveTailHalfChain(j,LogWF,dt,Lstag,Np,Path,xend)
                  call StagingHalfChain(j,LogWF,dt,Lstag,Np,Path,xend)
                  
               end do

            end do

            call OBDM(xend,nrho)

         end do

      end if

      !Energy calculation using mixed estimator
      
      call LocalEnergy(LogWF,Path(:,:,0),Rm,E,Kin,Pot)
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

   print *, '-----------------------------------------------------------'
   print *, 'BLOCK NUMBER :',iblock
   print *, ''
   print *, '# Block results:'
   print *, ''
   print *, '  > <E>      =',BlockAvE/Np,'+/-',BlockVarE/Np
   print *, '  > <Ec>     =',BlockAvK/Np,'+/-',BlockVarK/Np
   print *, '  > <Ep>     =',BlockAvV/Np,'+/-',BlockVarV/Np
   print *, ''
   print *, '# Acceptance CM moves      =',real(acc_cm)/real(attempted)
   print *, '# Acceptance staging moves =',real(acc_bd)/real(attempted*Nbis)
   print *, '# Acceptance head moves    =',real(acc_head)/real(attempted*Nbis)
   print *, '# Acceptance tail moves    =',real(acc_tail)/real(attempted*Nbis)
   print *, '# Time per block           =',end-begin
  
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
print *, '  > <E>  =',AvE/Np,'+/-',VarE/Np
print *, '  > <Ec> =',AvK/Np,'+/-',VarK/Np
print *, '  > <Ep> =',AvV/Np,'+/-',VarV/Np
print *, ''
print *, '=============================================================='

call NormAvGr(Nblock,AvGr,AvGr2,VarGr)
call NormAvSk(Nblock,Nk,AvSk,AvSk2,VarSk)
call NormAvNr(Nblock,AvNr,AvNr2,VarNr)

deallocate (Path)
deallocate (nrho)
deallocate (gr,AvGr,AvGr2,VarGr)
deallocate (Sk,AvSk,AvSk2,VarSk)

end program vpi

!-----------------------------------------------------------------------
