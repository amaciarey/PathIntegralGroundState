function Interpolate(opt,N,dx,F,x)

implicit none

real (kind=8)    :: Interpolate,dx,x
real (kind=8)    :: aux1,aux2
real (kind=8)    :: Fbefore,Fafter,Fcurr
integer (kind=4) :: opt
integer (kind=4) :: N,ix

real (kind=8),dimension (0:N+1) :: F

ix   = int(x/dx)+1
aux1 = x-(ix-1)*dx
aux2 = dx-aux1 

Interpolate = 0.d0

if (opt==0) then
   
   Interpolate = (aux1*F(ix)+aux2*F(ix-1))/dx

else if (opt==1) then

   Fbefore = (aux1*F(ix-1)+aux2*F(ix-2))/dx
   Fafter  = (aux1*F(ix+1)+aux2*F(ix))/dx

   Interpolate = 0.5d0*(Fafter-Fbefore)/dx

else if (opt==2) then

   Fbefore = (aux1*F(ix-1)+aux2*F(ix-2))/dx
   Fcurr   = (aux1*F(ix)+aux2*F(ix-1))/dx
   Fafter  = (aux1*F(ix+1)+aux2*F(ix))/dx

   Interpolate = (Fafter-2.d0*Fcurr+Fbefore)/(dx*dx)

else 
   
   print *, 'Parameter opt in function interpolate crash!!!!!'

end if

return
end function Interpolate
