module global_mod

implicit none

logical          :: wf_table
real (kind=8)    :: pi
real (kind=8)    :: V0
real (kind=8)    :: Rm
real (kind=8)    :: rbin,dr
real (kind=8)    :: rcut,rcut2
real (kind=8)    :: CWorm
integer (kind=4) :: dim,Np,Nbin,Nb,Nmax,Npw

real (kind=8),dimension (:),allocatable :: Lbox,LboxHalf,qbin

end module global_mod
