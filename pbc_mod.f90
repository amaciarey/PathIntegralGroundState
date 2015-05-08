module pbc_mod

use global_mod, only: dim,Lbox,LboxHalf

implicit none

contains

!-----------------------------------------------------------------------

  subroutine BoundaryConditions(k,xij)

    implicit none

    real (kind=8)    :: xij
    integer (kind=4) :: k
    
    if (xij>Lbox(k)) xij = xij-Lbox(k)
    if (xij<0.d0)    xij = xij+Lbox(k)

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
       
       if (xij(k)> LboxHalf(k)) xij(k) = xij(k)-Lbox(k)
       if (xij(k)<-LboxHalf(k)) xij(k) = xij(k)+Lbox(k)
       
       rij2 = rij2+xij(k)*xij(k)
       
    end do
    
    return
  end subroutine MinimumImage

!-----------------------------------------------------------------------

end module pbc_mod
