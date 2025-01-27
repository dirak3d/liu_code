subroutine sum_density(ntotal,hsml,mass,niac,pair_i,pair_j , w, itype,rho)
  ! Subroutine to calculate the density with SPH summation algorithm.
  ! ntotal : Number of particles [in]
  ! hsml : Smoothing Length [in]
  ! mass : Particle masses [in]
  ! niac : Number of interaction pairs [in]
  ! pair_i : List of first partner of interaction pair [in]
  ! pair_j : List of second partner of interaction pair [in]
  ! w : Kernel for all interaction pairs [in]
  ! itype : type of particles [in]
  ! x : Coordinates of all particles [in]
  ! rho : Density [out]
  implicit none
  include 'param.inc'
  integer ntotal, niac, pair_i(max^interaction), pair_j(max_interaction), itype(maxn)
  double precision hsml(maxn),mass(maxn), w(max_interaction), rho (maxn)
  integer i, j, k, d
  double precision selfdens, hv(dim), r, wi(maxn)
  ! wi(maxn) integration'of the kernel itself
  do d=l,dim
     hv(d) = O.eO
  enddo
  ! Self density of each particle: Wii (Kernel for distance 0)
  ! and take contribution of particle itself:
  r=0.0e0
  ! Firstly calculate the integration of the kernel over the space
  do i=l,ntbtal
     call" "Kernel (r.hv.hsml (i) , selfdens, hv)
     wi(i)=selfdens*mass(i)/rho(i)
  enddo
  
  do k=l,niac
     i = pair_i(k)
     j = pair_j(k)
     wi(i) = wi(i) + mass (j) /rho (j)*w(k)
     wi(j) = wi(j) + mass (i)/rho(i) *w(k)
  enddo
  
  ! Secondly calculate the rho integration over the space
  do i=l,ntotal
     call kernel(r,hv,hsml(i),selfdens,hv)
     rho(i) = selfdens*mass(i)
  enddo
  ! Calculate SPH sum for rho:
  do k=l,niac
     i = pair i(kl
     j = pair~j(k)
     rho(i) - rho(i) + mass(j)*w(k)
     rho(j) - rho(j) + mass(i)*w(k)
  enddo
  c Thirdly, calculate the normalized rho, rho=sum(rho)/sum(w)
  if (nor_density) then
     do i=l, ntotal
        rho(i)=rho(i)/wi(i)
     enddo
  endif
end subroutine sum_density

subroutine con_density(ntotal,mass,niac,pair_i,pair_j, dwdx,vx,itype,x,rho, drhodt)
  ! Subroutine to calculate 'the density with SPH continuiity approach.
  ! ntotal : Number of particles [in]
  ! mass : Particle masses [in]
  ! niac : Number of interaction pairs [in]
  ! pair_i : List of first partner of interaction pair [in]
  ! pair_j : List of second partner of interaction pair [in]
  ! dwdx : derivation of Kernel for all interaction pairs [in]
  ! vx : Velocities of all particles [in]
  ! itype : type of particles [in]
  ! x : Coordinates of all particles [in]
  ! rho : Density [in]
  ! drhodt : Density change rate of each particle [out]
  implicit none
  include 'param.inc'
  integer ntotal,niac,pair_i(max_interaction), pair_j(max_interaction), itype(maxn)
  double precision mass(maxn), dwdx(dim, max_interaction), vx(dim,maxn), x(dim,maxn), rho(maxn), drhodt(maxn)
  integer i,j,k,d
  double precision vcc, dvx(dim)
  do i = 1, ntotal
     drhodt(i) = 0.0e0
  enddo
  do k=1,niac
     i = pair_i(k)
     j = pair_j(k)
     do d=1,dim
        dvx(d) = vx(d,i) - vx(d,j)
     enddo
     vcc = dvx(1)*dwdx(1,k)
     do d=2, dim
        vcc = vcc + dvx(d)*dwdx(d,k)
     enddo
     drhodt(i) = drhodt(i) + mass(j)*vcc
     drhodt(j) = drhodt(j) + mass(i)*vcc
  enddo
end subroutine con_density
