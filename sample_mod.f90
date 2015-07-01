module sample_mod

use vpi_mod 
use global_mod
use system_mod
use pbc_mod

implicit none

contains

!-----------------------------------------------------------------------

  subroutine PotentialEnergy(R,Pot,F2)

    implicit none
    
    real (kind=8)    :: Pot
    real (kind=8)    :: rij2,rij
    integer (kind=4) :: ip,jp
    integer (kind=4) :: k

    real (kind=8),optional :: F2
       
    real (kind=8),dimension(dim,Np) :: R,F
    real (kind=8),dimension(dim)    :: xij
    
    Pot = 0.d0
        
    do ip=1,Np
       do k=1,dim
          F(k,ip) = 0.d0
       end do
    end do
    
    do ip=1,Np-1
       do jp=ip+1,Np
          
          rij2 = 0.d0
          
          do k=1,dim
             
             xij(k) = R(k,ip)-R(k,jp)
             
          end do

          call MinimumImage(xij,rij2)
          
          if (rij2<=rcut2) then
             
             rij = sqrt(rij2)
             
             Pot = Pot+Potential(xij,rij)
             
             if (present(F2)) then
                do k=1,dim
                   
                   F(k,ip) = F(k,ip)+Force(k,xij,rij)
                   F(k,jp) = F(k,jp)-Force(k,xij,rij)

                end do
             end if
             
          end if
          
       end do
    end do
    
    if (present(F2)) then

       F2 = 0.d0
       
       do ip=1,Np
          do k=1,dim
             F2 = F2+F(k,ip)*F(k,ip)
          end do
       end do

    end if
    
    return 
  end subroutine PotentialEnergy
  
!-----------------------------------------------------------------------

  subroutine LocalEnergy(LogWF,R,Rm,E,Kin,Pot)

    implicit none
    
    real (kind=8)    :: E,Kin,Pot,Rm
    real (kind=8)    :: rij,rij2
    real (kind=8)    :: dudr,d2udr2
    real (kind=8)    :: LapLogPsi
    real (kind=8)    :: Interpolate
    
    integer (kind=4) :: i,j
    integer (kind=4) :: k
        
    real (kind=8),dimension (0:Nmax+1) :: LogWF
    real (kind=8),dimension (dim,Np)   :: R,F
    real (kind=8),dimension (dim)      :: xij
    
    Kin       = 0.d0
    Pot       = 0.d0
    LapLogPsi = 0.d0
    
    do i=1,Np
       do k=1,dim
          F(k,i) = 0.d0
       end do
    end do
    
    do i=1,Np-1
       do j=i+1,Np
          
          !Calculation of the distance between each pair of dipoles
          !choosing only the nearest image of each particle
          
          do k=1,dim
          
             xij(k) = R(k,i)-R(k,j)
          
          end do

          call MinimumImage(xij,rij2)
          
          !Calculation of the local kinetic energy
          !We use in this routine the following prescription:
          
          ! Kin = -Laplacian(Psi)/(2*Psi)
          ! T   = -Laplacian(log(Psi))/4
          ! F   = Grad(Psi)/(sqrt(2)*Psi)
          ! Kin = 2*T-F*F
          
          !We use linear interpolation between the nearest points of the
          !tabulated wave function, as in the previous subroutine we use:
          !
          !    x_{i-1} < z < x_i
          !
          !    psi(z) = psi(x_i)*(z-x_{i-1})/dr+psi(x_{i-1})*(x_i-z)/dr   
          !
          !The evaluation of derivatives is also numeric by interpolating 
          !in the next and previous intervals.
          
          if (rij2<=rcut2) then
             
             rij = sqrt(rij2)

             if (wf_table) then
                
                dudr   = Interpolate(1,Nmax,dr,LogWF,rij)
                d2udr2 = Interpolate(2,Nmax,dr,LogWF,rij)
                
             else
             
                dudr   = LogPsi(1,Rm,rij)
                d2udr2 = LogPsi(2,Rm,rij)
             
             end if

             LapLogPsi = LapLogPsi+((real(dim)-1)*dudr/rij+d2udr2)        
             
             do k=1,dim
                F(k,i) = F(k,i)+dudr*xij(k)/rij
                F(k,j) = F(k,j)-dudr*xij(k)/rij
             end do
             
             !Calculation of the local potential energy 
             
             Pot = Pot+Potential(xij,rij)
             
          end if
          
       end do
    end do
    
    !Kinetic energy and drift force calculation
    
    Kin = 2.d0*LapLogPsi
    
    do i=1,Np
       do k=1,dim
          Kin = Kin+F(k,i)*F(k,i)
       end do
    end do
    
    !Computing energy
    
    Kin = -0.5d0*Kin
    E   = Kin+Pot
    
    return
  end subroutine LocalEnergy

!-----------------------------------------------------------------------

  subroutine ThermEnergy(Path,dt,E,Ec,Ep)

    implicit none

    real (kind=8)    :: E,Ec,Ep
    real (kind=8)    :: dt
    real (kind=8)    :: Pot,F2
    real (kind=8)    :: rij2
    integer (kind=4) :: ib
    integer (kind=4) :: ip
    integer (kind=4) :: k
        
    real (kind=8),dimension(dim,Np,0:2*Nb) :: Path
    real (kind=8),dimension(dim)           :: xij
        
    Ep = 0.d0
    Ec = 0.d0
    E  = 0.d0
    
    do ib=0,2*Nb-1
       
       if (mod(ib,2)==0) then
          call PotentialEnergy(Path(:,:,ib),Pot)
          F2 = 0.d0
       else
          call PotentialEnergy(Path(:,:,ib),Pot,F2)
       end if

       if (ib==Nb) then
          Ep = Pot
       end if
       
       E = E+GreenFunction(1,ib,dt,Pot,F2)
              
       do ip=1,Np
          
          rij2 = 0.d0
          
          do k=1,dim
             
             xij(k) = Path(k,ip,ib)-Path(k,ip,ib+1)
             
          end do
          
          call MinimumImage(xij,rij2)

          if (rij2<=rcut2) E = E-0.5d0*rij2/(dt*dt)
          
       end do
       
    end do

    call PotentialEnergy(Path(:,:,2*Nb),Pot)
    
    E  = E+GreenFunction(1,2*Nb,dt,Pot,0.d0)
    E  = 0.5d0*(E/real(Nb)+real(dim*Np)/dt)
    Ec = E-Ep
    
    return
  end subroutine ThermEnergy
  
!-----------------------------------------------------------------------

  subroutine PairCorrelation(R,gr)

    implicit none
    
    real (kind=8)    :: rij2,rij
    integer (kind=4) :: ip,jp
    integer (kind=4) :: k
    integer (kind=4) :: ibin
    
    real (kind=8),dimension (dim)    :: xij
    real (kind=8),dimension (dim,Np) :: R
    real (kind=8),dimension (Nbin)   :: gr
    
    do ip=1,Np-1
       do jp=ip+1,Np
          
          rij2 = 0.d0
          
          do k=1,dim
             
             xij(k) = R(k,ip)-R(k,jp)
             
          end do
          
          call MinimumImage(xij,rij2)

          if (rij2<=rcut2) then
             
             rij = sqrt(rij2)
             
             ibin     = int(rij/rbin)+1
             gr(ibin) = gr(ibin)+2.d0
             
          end if
          
       end do
    end do
    
    return
  end subroutine PairCorrelation

!-----------------------------------------------------------------------

  subroutine StructureFactor(Nk,R,Sk)

    implicit none

    real (kind=8)    :: qr
    real (kind=8)    :: SumCos,SumSin
    integer (kind=4) :: iq,Nk
    integer (kind=4) :: k
    integer (kind=4) :: ip
  
    real (kind=8),dimension(dim)    :: x
    real (kind=8),dimension(dim,Np) :: R
    real (kind=8),dimension(dim,Nk) :: Sk

    SumCos = 0.d0
    SumSin = 0.d0

    do iq=1,Nk
   
       do k=1,dim

          SumCos = 0.d0
          SumSin = 0.d0
          
          do ip=1,Np
             
             x(k) = R(k,ip)
             qr   = real(iq)*qbin(k)*x(k)
             
             SumCos = SumCos+cos(qr)
             SumSin = SumSin+sin(qr)
             
          end do
          
          Sk(k,iq) = Sk(k,iq)+(SumCos*SumCos+SumSin*SumSin)
          
       end do
       
    end do

    return
  end subroutine StructureFactor

!-----------------------------------------------------------------------

  subroutine OBDM(xend,nrho)

    implicit none

    real (kind=8)    :: rij2,rij
    real (kind=8)    :: costheta,sintheta
    real (kind=8)    :: cos2mtheta,sin2mtheta
    complex (kind=8) :: exptheta,exp2theta,exp2mtheta
    integer (kind=4) :: ibin,k,m

    real (kind=8), dimension (dim,2)      :: xend
    real (kind=8), dimension (dim)        :: xij
    real (kind=8), dimension (0:Npw,Nbin) :: nrho
    
    do k=1,dim
       xij(k) = xend(k,1)-xend(k,2)
    end do

    call MinimumImage(xij,rij2)

    if (rij2<=rcut2) then
       
       rij  = sqrt(rij2)
       ibin = int(rij/rbin)+1

       sintheta  = xij(2)/rij
       costheta  = xij(1)/rij
       exptheta  = cmplx(costheta,sintheta,8)
       exp2theta = exptheta*exptheta

       exp2mtheta = 1.d0

       do m=0,Npw

          cos2mtheta = real(exp2mtheta)
          sin2mtheta = aimag(exp2mtheta)

          nrho(m,ibin) = nrho(m,ibin)+cos2mtheta

          exp2mtheta = exp2mtheta*exp2theta

       end do

    end if

    return
  end subroutine OBDM
  
!-----------------------------------------------------------------------

  subroutine NormalizeGr(density,ngr,gr)

    implicit none

    real (kind=8)    :: r8_gamma
    real (kind=8)    :: density
    real (kind=8)    :: nid,r,norm
    real (kind=8)    :: k_n
    integer (kind=4) :: ibin
    integer (kind=4) :: ngr
    
    real (kind=8),dimension (Nbin) :: gr
        
    k_n  = pi**(0.5d0*dim)/r8_gamma(0.5d0*dim+1.d0)
    norm = real(Np)*real(ngr)

    do ibin=1,Nbin
       r   = (real(ibin)-0.5d0)*rbin
       nid = density*k_n*((r+0.5d0*rbin)**dim-(r-0.5d0*rbin)**dim)
       gr(ibin) = gr(ibin)/(nid*norm)
    end do

    return
  end subroutine NormalizeGr

!-----------------------------------------------------------------------

  subroutine NormalizeSk(Nk,ngr,Sk)

    implicit none
    
    real (kind=8)    :: norm
    integer (kind=4) :: k,ik,Nk
    integer (kind=4) :: ngr

    real (kind=8), dimension (dim,Nk) :: Sk

    norm = real(Np)*real(ngr)

    do ik=1,Nk
       do k=1,dim
          Sk(k,ik) = Sk(k,ik)/norm
       end do
    end do

    return
  end subroutine NormalizeSk

!-----------------------------------------------------------------------

  subroutine NormalizeNr(density,zconf,Nobdm,nrho)

    implicit none

    real (kind=8)    :: r8_gamma
    real (kind=8)    :: density
    real (kind=8)    :: k_n,nid
    real (kind=8)    :: zconf
    real (kind=8)    :: r

    integer (kind=4) :: ibin,m
    integer (kind=4) :: Nobdm

    real (kind=8), dimension (0:Npw,Nbin) :: nrho
    
    k_n = pi**(0.5d0*dim)/r8_gamma(0.5d0*dim+1.d0)

    do ibin=1,Nbin
       r   = (real(ibin)-0.5d0)*rbin
       nid = density*k_n*((r+0.5d0*rbin)**dim-(r-0.5d0*rbin)**dim)
       do m=0,Npw
          nrho(m,ibin) = nrho(m,ibin)/(CWorm*nid*zconf*real(Nobdm))
       end do
       !write (98,'(20g20.10e3)') r,(nrho(m,ibin),m=0,Npw)
    end do
    !write (98,*) 
    !write (98,*)
    
    return
  end subroutine NormalizeNr

!-----------------------------------------------------------------------

  subroutine AccumGr(gr,AvGr,AvGr2)

    implicit none
    
    integer (kind=4) :: j

    real (kind=8),dimension(Nbin) :: gr,AvGr,AvGr2

    do j=1,Nbin
       AvGr(j)  = AvGr(j)+gr(j)
       AvGr2(j) = AvGr2(j)+gr(j)*gr(j)
    end do

    return
  end subroutine AccumGr

!-----------------------------------------------------------------------

  subroutine AccumSk(Nk,Sk,AvSk,AvSk2)

    implicit none
    
    integer (kind=4) :: j,k,Nk

    real (kind=8),dimension(dim,Nk) :: Sk,AvSk,AvSk2

    do j=1,Nk
       do k=1,dim
          AvSk(k,j)  = AvSk(k,j)+Sk(k,j)
          AvSk2(k,j) = AvSk2(k,j)+Sk(k,j)*Sk(k,j)
       end do
    end do

    return
  end subroutine AccumSk

!-----------------------------------------------------------------------

 subroutine AccumNr(nrho,AvNr,AvNr2)

    implicit none
    
    integer (kind=4) :: j,m

    real (kind=8),dimension(0:Npw,Nbin) :: nrho,AvNr,AvNr2

    do j=1,Nbin
       do m=0,Npw
          AvNr(m,j)  = AvNr(m,j)+nrho(m,j)
          AvNr2(m,j) = AvNr2(m,j)+nrho(m,j)*nrho(m,j)
       end do
    end do

    return
  end subroutine AccumNr

!-----------------------------------------------------------------------

  subroutine NormAvGr(Nitem,AvGr,AvGr2,VarGr)

    implicit none

    real (kind=8)    :: r
    integer (kind=4) :: Nitem,j

    real (kind=8),dimension(Nbin) :: AvGr,AvGr2,VarGr
    
    open (unit=11,file='gr_vpi.out')

    do j=1,Nbin
       r        = (real(j)-0.5d0)*rbin
       AvGr(j)  = AvGr(j)/real(Nitem)
       AvGr2(j) = AvGr2(j)/real(Nitem)
       VarGr(j) = Var(Nitem,AvGr(j),AvGr2(j))
       write (11,'(20g20.10e3)') r,AvGr(j),VarGr(j)
    end do
    
    close (unit=11)

    return
  end subroutine NormAvGr

!-----------------------------------------------------------------------

  subroutine NormAvSk(Nitem,Nk,AvSk,AvSk2,VarSk)

    implicit none

    integer (kind=4) :: Nitem,Nk,j,k

    real (kind=8),dimension(dim,Nk) :: AvSk,AvSk2,VarSk
    
    open (unit=11,file='sk_vpi.out')

    do j=1,Nk
       do k=1,dim
          AvSk(k,j)  = AvSk(k,j)/real(Nitem)
          AvSk2(k,j) = AvSk2(k,j)/real(Nitem)
          VarSk(k,j) = Var(Nitem,AvSk(k,j),AvSk2(k,j))
       end do
       write (11,'(20g20.10e3)') (j*qbin(k),AvSk(k,j),VarSk(k,j),k=1,dim)
    end do

    close (unit=11)

    return
  end subroutine NormAvSk

!-----------------------------------------------------------------------

  subroutine NormAvNr(Nitem,AvNr,AvNr2,VarNr)

    implicit none

    real (kind=8)    :: r
    integer (kind=4) :: Nitem,j,m

    real (kind=8),dimension(0:Npw,Nbin) :: AvNr,AvNr2,VarNr
    
    open (unit=11,file='nr_vpi.out')

    do j=1,Nbin
       r = (real(j)-0.5d0)*rbin
       do m=0,Npw
          AvNr(m,j)  = AvNr(m,j)/real(Nitem)
          AvNr2(m,j) = AvNr2(m,j)/real(Nitem)
          VarNr(m,j) = Var(Nitem,AvNr(m,j),AvNr2(m,j))
       end do
       write (11,'(20g20.10e3)') r,(AvNr(m,j),VarNr(m,j),m=0,Npw)
    end do
    
    close (unit=11)

    return
  end subroutine NormAvNr

!-----------------------------------------------------------------------
!
! This subroutine simply accumulates the value of the energies at each 
! step/block in order to obtain average values.
!
!-----------------------------------------------------------------------

  subroutine Accumulate(E,Ec,Ep,SumE,SumEc,SumEp)

    implicit none

    real (kind=8) :: E,SumE
    real (kind=8) :: Ec,SumEc
    real (kind=8) :: Ep,SumEp
    
    SumE  = SumE+E
    SumEc = SumEc+Ec
    SumEp = SumEp+Ep
    
    return
  end subroutine Accumulate

!-----------------------------------------------------------------------
!
! This subroutine evaluates the average of a magnitude along a step, 
! block or the whole simulation.
!
!-----------------------------------------------------------------------

  subroutine NormalizeAv(Nitem,SumE,SumEc,SumEp)

    implicit none

    real (kind=8)    :: SumE,SumEc,SumEp
    integer (kind=4) :: Nitem

    SumE  = SumE/real(Nitem)
    SumEc = SumEc/real(Nitem)
    SumEp = SumEp/real(Nitem)

    return
  end subroutine NormalizeAv

!-----------------------------------------------------------------------
!
! This function simply evaluates the variance of a variable.
!
!-----------------------------------------------------------------------

  function Var(Nitem,Sum,Sum2)

    implicit none
    
    real (kind=8)    :: Var
    real (kind=8)    :: Sum,Sum2
    integer (kind=4) :: Nitem

    Var = sqrt((Sum2-Sum*Sum)/real(Nitem))

    return
  end function Var

!-----------------------------------------------------------------------

end module sample_mod
