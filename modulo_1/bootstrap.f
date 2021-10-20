	program bootstrap
	implicit none
	real*8, dimension(:), allocatable ::  u, s, g, z
	real*8 :: ran2, a, media, sig, b
	integer(16) :: i, N, j, l, M, k
	integer(16) :: nb, Dd
	integer :: R, q
	
	call ranstart
	open(unit=1, file="metrogausserr.txt", status="unknown")
	open(unit=0, file="metrogauss.txt", status="old", action="read")
	
	read(0, *) N			!grandezza totale del campione
	read(0, *) M			!grandezza del campione analizzato
	R=150				!numero dei resampling da fare
	
	allocate(z(M), g(M), u(R), s(N))
				
	do i = 1, N			!prima e terza colonna non interessano
		read(0, *) a,  s(i), a	!per l'analisi, quindi sono messi nella
	enddo				!variabile sempre sovrascritta e mai più usata
	
	
	do i=1, M			!elimino la parte in cui la catena 
		g(i)=s((N - M) + i)	!potrebbe non essere termalizzata
	enddo
	
	do q=1, 20
		Dd = 2**q		!dimensioni dei blocchi
		nb=M/Dd			!grandezza dei blocchi
		do l =1, R		!ciclo dei ricampionamenti
			do i = 1, nb
				b = ran2()
				j = int(b*M +1)		!scelgo un sito a caso
				do k = 1, Dd		!inserisco i successivi Dd 
							!elementi in un array temporaneo
	
					z((i-1)*Dd+k) = g(mod(j+k-2,M)+1)	
				enddo
			enddo
			u(l) = sum(z)/float(M)		!caclolo media dei ricampionamenti
		enddo
		
		media = sum(u)/float(R)			!calcolo media della media ed errore
		sig=0					
		
		do i=1, R
			sig = sig + (u(i) -media)**2
		enddo
		
		sig = sqrt(sig/float(R - 1))
		write(1,*) sig, Dd
	enddo
		
	call ranfinish
	
	end program bootstrap
	
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
