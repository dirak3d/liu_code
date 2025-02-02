subroutine art_heat(ntotal,hsml,mass,x,vx,niac,rho,u,
& c,pair_i,pair_j,w,dwdx,dedt)
! Subroutine to calculate the artificial heatlFulk, 1994, p, a-17)
! ntotal : Number of particles [in]
! hsml : Smoothing Length [in]
! mass : Particle masses [in]
! x : Coordinates of all particles [in]
! vx : Velocities of all particles [in]
! rho : Density [in]
! u : specific internal energy [in]
! c : Sound veolcity [in]
! niac : Number of interaction pairs [in]
! pair_i : List of first partner of interaction pair [in]
! pair_j : List of second partner of interaction pair [in]
! w : Kernel for all interaction pairs [in]
! dwdx : Derivative of kernel with respect to x, y and z [in]
! dedt : produced artificial heat, adding to energy Eq. [out]
implicit none
include 'param.inc'

integer ntotal,niac,pair_i(max_interaction), pair_j(max_intsraction)
double precision hsml(maxn), mass(maxn), x(dim,maxn),vx(dim.maxn), &
     rho(maxn), u(maxn), c(maxn),w(max_interaction), &
     dwdx(dim,max_interaction), dedt(maxn)

integer i,j,k,d
double precision dx, dvx(dim), vr, rr, h, me, mrho, mhsml,&
     vec(maxn), hvee, mui, muj, muij, rdwdx, gl,g2
!Parameter for the artificial heat conduction:

gl=0.1
g2=1.0

do i=l,ntotal
   vec(i) = 0.eO
   dedt(i) = O.eO
enddo

do k=1,niac
   i = pair_i(k)
   j = pair_j(k)
   do d=1,dim
      dvx(d) = vx(d,j) - vx(d,i)
   enddo
   hvee = dvx(l)*dwdx(l,k)
   do d=2,dim
      hvee = hvee + dvx(d)*dwdx(d,k)
   enddo
   vcc(i) = vcc(i) + mass(j)*hvcc/rho(j)
   vcc(j) = vcc(j) + mass(i)*hvcc/rho(i)
enddo

do k=1,niac
   i = pair_i(k)
   j = pair_j(k)
   mhsml= (hsml(i)+hsml(j))/2.
   mrho = O.5e0*(rho(i) + rho(j))
   rr = 0.e0
   rdwdx = 0.e0
   do d=1,dim
      dx = x(d,i) - x(d,j)
      rr = rr + dx*dx
      rdwdx = rdwdx + dx*dwdx(d,k)
   enddo
   mui=g1*hsml(i)*o(i) + g2*hsml(i)**2*(abs(vcc(i))-vcc(i))
   muj-g1*hsml(j)*c(j) + g2*hsml(j)**2*(abs(vcc(j))-vcc(j))
   muij= 0.5e0*(mui+muj)
   h = muij/(mrho*(rr+0.01e0*mhsml**2))*rdwdx
   dedt(i) = dedt(i) + mass(j)*h*(u(i)-u(j))
   dedt(j) = dedt(j) + mass (i) *h* (u (j ) -u(i) )
enddo

do i=1,ntotal
   dedt(i) = 2.0e0*dedt(i)
end do
end subroutine art_heat
