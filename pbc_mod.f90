module pbc_mod

use global_mod, only: dim,Lbox,LboxHalf

implicit none

contains

!-----------------------------------------------------------------------

  subroutine BoundaryConditions(k,xij)

    implicit none

    real (kind=8)    :: xij
    integer (kind=4) :: k
    
    !if (xij>Lbox(k)) xij = xij-Lbox(k)
    !if (xij<0.d0)    xij = xij+Lbox(k)

    !if (xij> LboxHalf(k)) xij = xij-Lbox(k)
    !if (xij<-LboxHalf(k)) xij = xij+Lbox(k)

    if (abs(xij) > LboxHalf(k)) xij = xij-sign(Lbox(k),xij)

    return
  end subroutine BoundaryConditions

!-----------------------------------------------------------------------

  subroutine MinimumImage(xij,rij2)

    implicit none

    real (kind=8)    :: rij2
    integer (kind=4) :: k

    real (kind=8),dimension (dim) :: xij

    rij2 = 0.d0
    
    do k=1,dim
       
       !if (xij(k)> LboxHalf(k)) xij(k) = xij(k)-Lbox(k)
       !if (xij(k)<-LboxHalf(k)) xij(k) = xij(k)+Lbox(k)
       
       if (abs(xij(k)) > LboxHalf(k)) then
          xij(k) = xij(k)-sign(Lbox(k),xij(k))
       end if
       
       rij2 = rij2+xij(k)*xij(k)
       
    end do

    return
  end subroutine MinimumImage

!-----------------------------------------------------------------------

  subroutine MinimumImageDistance(k,xij)

    implicit none
    
    real (kind=8)    :: xij
    integer (kind=4) :: k
    
    !if (xij < -LboxHalf(k)) xij = xij+Lbox(k)
    !if (xij >  LboxHalf(k)) xij = xij-Lbox(k)

    if (abs(xij) > LboxHalf(k)) xij = xij-sign(Lbox(k),xij)

    return
  end subroutine MinimumImageDistance

!-----------------------------------------------------------------------

end module pbc_mod
