subroutine direc_find(itimestep, ntotal, hsml, niac, pair_i, pair_j, w, dwdx, countiac)
  implicit none
  include 'param.inc'

  integer itimestep, ntotal, niac, pair_i(max_interaction), pair_j(max_interaction), countiac(maxn)
  double precision hsml(maxn), x(dim,maxn), w(max_interaction), dwdx(dim,max_interaction)
  integer i, j, d, sumiac, maxiac, miniac, noiac, maxp, minp, scale_k
  double precision dxiac(dim), driac, r, mhsml, tdwdx(dim)

  if(skf.eq. 1) then
     scale_k = 2
  else if (skf.eq. 2) then
     scale_k = 3
  else if (skf.eq. 3) then
     scale_k = 3
  endif


  do i=1,ntotal
     countiac(i) = 0
  enddo
  
  niac = 0

  do i=1,ntotal-1
     do j=1+1, ntotal
        dxiac(1) = x(1,i) -x(1,j)
        driac = dxiac(1)*dxiac(1)
        do d=2,dim
           dxiac(d) = x(d,i) - x(d,j)
           driac = driac + dxiac(d)*dxiac(d)
        enddo
        mhsml = (hsml(i) + hsml(j))/2.0
        if(sqrt(driac).lt.scale_k*mhsml) then
           if(niac.lt.max_interaction) then
        !Neigboring pair list, and totalinteraction nubmer and
        !the interaction nubmer for each particle
              
              niac = niac + 1
              pair_i(niac) = i
              pair_j(niac) = j
              r = sqrt(driac)
              countiac(i) = countiac(i) + 1
              countiac(j) = countiac(j) + 1
              !Kernel and derviation of kernel

              call kernel(r,dxiac,mhsml, w(niac), tdwdx)
              do d=1,dim
                 dwdx(d,niac) = tdwdx(d)
              enddo
           else
              print *, ">>> ERROR <<<: Too many interactions"
              stop
           endif
        endif
     enddo
  enddo

  !Estadísticas para la interacción

  sumiac = 0
  maxiac = 0
  miniac = 1000
  noiac = 0

  do i=1,ntotal
     sumiac = sumiac + countiac(i)
     if (countiac(i).gt. maxiac) then
        maxiac = countiac(i)
        maxp =  i
     endif
     if (countiac(i).lt. miniac) then
        miniac = countiac(i)
        minp = i
     endif
     if (countiac(i).eq.0) noiac = noiac + 1 
  enddo

  if (mod(itimestep, print_step).eq.0) then
     if(int_stat) then
        print *, ">> Statistics: interactions per particles:"
        print *, "... Particle: ", maxp, "maximal interactions: ", maxiac
        print *, "... Particle: ", minp, "minimal interacitons: ", miniac
        print *, "... Average: ", real(sumiac)/real(ntotal)
        print *, "... Total pairs : ", niac
        print *, "... Particles with no interactions", noiac
        
     endif
  endif
  
end subroutine direc_find
