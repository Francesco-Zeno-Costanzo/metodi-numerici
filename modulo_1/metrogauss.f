	program campgauss
	implicit none
	real*8 :: aver, s2, d, start, x, y, q, q_try, z, acc, ran2, Gauss
	integer(16) :: N, i
	
	call ranstart
	N=2**24 + 10000 !con l'intenzione di scartare i primi 10000 per termalizzazione
	aver=5.		!media della distribuzione da campionare
	s2=1.		!varianza della distribuzione da campionare
	d=0.1		!passo dell'algoritmo
	start=0.0	!punto di partenza
	
	acc=0
	q = start
	
	open(1, file='metrogauss.txt', status="unknown")
	write(1, *) N		!grandezza del campione generato
	write(1, *) N-10000	!grandezza del campione da analizzare
	
	do i=1, N		!campionamento tramite metropolis
		x = ran2()
		y = ran2()
		
		q_try = q + d*(2.0*x - 1.0)
		
		z=Gauss(q_try, aver, s2)/Gauss(q, aver, s2)
		
		if(y<z) then
			q = q_try
	  		acc = acc+1                     
		endif
        	write(1, *) i, q, acc                  
	enddo
	
	call ranfinish
	end program campgauss
	
	real*8 function Gauss(x, m, s)		!distribuzione da campionare
		implicit none
		real*8 :: x, m, s
		Gauss = exp((-(x-m)**2.)/(2.0*s))
		return
	end function Gauss
	
C=============================================================================
C Generazione di numeri casuali
C=============================================================================	
      function ran2()
      implicit real*4 (a-h,o-z)
      implicit integer*4 (i-n)
      integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      real*8 ran2,am,eps,rnmx
      parameter(im1=2147483563,im2=2147483399,am=1./im1,imm1=im1-1,
     &          ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,
     &          ir2=3791,ntab=32,ndiv=1+imm1/ntab,eps=1.2e-7,
     &          rnmx=1.-eps)
      integer idum2,j,k,iv,iy
      common /dasav/ idum,idum2,iv(ntab),iy
c      save iv,iy,idum2
c      data idum2/123456789/, iv/NTAB*0/, iy/0/

      if(idum.le.0) then
         idum=max0(-idum,1)
         idum2=idum
         do j=ntab+8,1,-1
            k=idum/iq1
            idum=ia1*(idum-k*iq1)-k*ir1
            if(idum.lt.0) idum=idum+im1
            if(j.le.ntab) iv(j)=idum
         enddo
         iy=iv(1)
      endif
      k=idum/iq1
      idum=ia1*(idum-k*iq1)-k*ir1
      if(idum.lt.0) idum=idum+im1
      k=idum2/iq2
      idum2=ia2*(idum2-k*iq2)-k*ir2
      if(idum2.lt.0) idum2=idum2+im2
      j=1+iy/ndiv
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1) iy=iy+imm1
      ran2=min(am*iy,rnmx)

      return
      end

c=============================================================================
      subroutine ranstart
      implicit real*4 (a-h,o-z)
      implicit integer*4 (i-n)
      common /dasav/ idum,idum2,iv(32),iy

      open(unit=23, file='randomseed', status='unknown')
      read(23,*) idum
      read(23,*,end=117) idum2
      do i=1,32
         read(23,*) iv(i)
      enddo
      read(23,*) iy
      close(23)
      goto 118                          !!takes account of the first start
 117  if(idum.ge.0) idum = -idum -1     !!
      close(23)
 118  continue                          !!

      return
      end

c=============================================================================
      subroutine ranfinish
      implicit real*4 (a-h,o-z)
      implicit integer*4 (i-n)
      common /dasav/ idum,idum2,iv(32),iy

      open(unit=23, file='randomseed', status='unknown')
      write(23,*) idum
      write(23,*) idum2
      do i=1,32
         write(23,*) iv(i)
      enddo
      write(23,*) iy
      close(23)

      return
      end
      
	
