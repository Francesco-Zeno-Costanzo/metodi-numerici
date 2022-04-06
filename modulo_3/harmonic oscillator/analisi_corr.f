	program analisi_correlazione
	
	real*8, dimension(:, :), allocatable :: Y
	real*8, dimension(:), allocatable :: aver_y, daver_y
	real*8, dimension(:), allocatable :: y2
	real*8 :: daver_y2
	integer :: R
	common N, R
	
	call cpu_time(start)
	call ranstart
	
	open(1, file='datiy2c.dat', status='old')
	open(2, file='y2toy2c.dat', status='old')
	open(3, file='datiy2cplot.dat', status='unknown')
	
	R = 100
	
	read(1, *) N
	read(1, *) nret
	read(1, *) K
	
	do l = 1, nret
		
		allocate(Y(N, K), aver_y(K), daver_y(K))
		allocate(y2(N))
		
		do i = 1, N				!ci sono N curve per ogni eta
			do j = 1, K			!ogni curva è lunga k.
				read(1, *) Y(i,j)	!elemento della curva
			enddo
		enddo
		
		!Y(:,i) contiene quindi tutti i primi punti delle n curve
		!la cui media darà il primo della curva finale
		
		!calcolo della media di y^2 per calcolare il correlatore connesso
		aver_y2 = 0
		do i = 1, N		
			read(2, *) y2(i)
		enddo
		aver_y2 = sum(y2)/float(N)
		call errore(y2, daver_y2)
		
		!calcolo correlatore connesso
		do i = 1, K
			aver_y(i) = sum(Y(:,i))/float(N) - aver_y2**2
			call errore(Y(:,i), daver_y(i))
		enddo
		daver_y = sqrt(daver_y**2 + (2*aver_y2*daver_y2)**2)
		write(3,*) aver_y, daver_y
		
		deallocate(Y, aver_y, daver_y, y2)
	enddo

	
	call ranfinish
	
	call cpu_time(finish)
	print '("tempo di esecuzione = ", f8.4," secondi.")', finish-start
	
	end program analisi_correlazione
	
C=============================================================================

	subroutine errore(x, dx)
	common N, R
	real*8, dimension(:), allocatable :: z, a
	real*8, dimension(N) :: x
	integer(16) :: nb, Dd, i, j, l
	real*8 :: media_x, dx
	integer :: g, R
	
	dx = 0
	
	allocate(z(N), a(R))

	!il calcolo dell'errore avviene tramite il binned bootstrap poichè
	!a priori non si può sapere qual è il tempo di decorellazione migliore
	!da inserire nella simulazione, e per non sprecare tempo macchina,
	!bisogna tenere conto di una possibile correlazione frai i dati
	
    	Dd = 2**14			!dimensione dei blocchi	
	nb=N/Dd				!numero dei blocchi
	do l = 1, R			!ciclo sui ricampionamenti
		do i = 1, nb
		
			j = int(ran2()*N +1) !scelgo sito a caso
			do g = 1, Dd		 	
				z((i-1)*Dd+g) = x(mod(j+g-2,N)+1) !ricampiono a blocchi	
			enddo
			
		enddo
		
		!calcolo della media e della varianza dei ricampionamenti
		
		a(l) = sum(z)/float(N)
			
	enddo
	
	!calcolo la media degli estimatori
	media_x  = sum(a)/float(R)
		
	do i = 1, R
		!calcolo scarto quadratico
		dx  = dx  + (a(i) - media_x )**2

	enddo
	
	dx = sqrt(dx/float(R - 1))		!prendo l'errore sul campione
						!perchè  sono stai effettuati
						!R ricampionamenti
						!quindi divido solo per R-1
						!e non R(R-1)
	return
	end
	

	
C=============================================================================
      function ran2()
      implicit real*4 (a-h,o-z)
      implicit integer*4 (i-n)
      integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      real*4 ran2,am,eps,rnmx
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
