module global_mod

implicit none

logical          :: table = .True.
real (kind=8)    :: pi
real (kind=8)    :: V0
real (kind=8)    :: Rm
real (kind=8)    :: N0,AK
real (kind=8)    :: a_1,t_0,t_1,u_0,v_1,v_2
real (kind=8)    :: rbin,dr
real (kind=8)    :: rcut,rcut2
integer (kind=4) :: dim,Np,Nbin,Nb,Nmax,Npw

real (kind=8),dimension (:),allocatable :: Lbox,LboxHalf,qbin

end module global_mod
