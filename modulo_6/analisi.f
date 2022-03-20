      program analisi
      
      real*8, dimension(:), allocatable :: fm2, ds2, dt2, dens
      real*8 :: dfm2, dds2, ddt2, ddens
      integer :: R, P
      integer, dimension(:), allocatable :: cp
      character(len=30) file_name, formatstring
      common N, R, M
	
      call cpu_time(start)
      call ranstart
C===================================================================
C Il codice è implementato per leggere P file da una cartella 'dati'
C labellati da un numero e queste informazini devono essere scritte
C nel file analisi_parm.txt
C===================================================================    
      open(unit=0, file='analisi_parm.txt', status="unknown")
      read(0,*) P
      allocate(cp(P))
      do j = 1, P
          read(0, *) cp(j)
      enddo
      
      open(unit=20, file="misure.dat", status="unknown")
      open(unit=21, file="densm.dat", status="unknown")
      
      R = 200
      
      do j = 1, P
          
          if (cp(j) < 10) then
              formatstring = "(A5,i1,A4)"
          else if ((cp(j)>=10).AND.(cp(j) < 100)) then
              formatstring = "(A5,i2,A4)"
          else if (cp(j) >= 100) then
              formatstring = "(A5,i3,A4)"
          endif
          
          write(file_name, formatstring) "dati/",cp(j),".dat"
          
          open(1, file=file_name, status="unknown")
          
          read(1, *) N
          read(1, *) M
          read(1, *) nx, nt
      
     
          allocate(fm2(N), ds2(N), dt2(N), dens(N-M)) 
      
          do i = 1, N
              read(1, *) fm2(i), ds2(i), dt2(i)				
          enddo
          
	  dens = fm2(M:) + ds2(M:) - dt2(M:)
	  
          a = sum(fm2(M:))/float(N-M)
          b = sum(ds2(M:))/float(N-M)
          c = sum(dt2(M:))/float(N-M)
          d = sum(dens)/float(N-M)
          
          call errore(fm2(M:), dfm2)
          call errore(ds2(M:), dds2)
          call errore(dt2(M:), ddt2)
          call errore(dens, ddens)      
      
          write(20,*) nt, a, b, c, dfm2, dds2, ddt2
          write(21,*) nt, d, ddens
          close(1)
          deallocate(fm2, ds2, dt2, dens)
          
      enddo
         
      call ranfinish
      
      call cpu_time(finish)
      print'("tempo di esecuzione= ", f16.8," secondi.")',finish-start
      	
      end program analisi
      
C=============================================================================

      subroutine errore(x, dx)
      common N, R, M
      real*8, dimension(:), allocatable :: z, a
      real*8, dimension(N-M) :: x
      integer(16) :: nb, Dd, i, j, l
      real*8 :: media_x, dx
      integer :: g, R
	
      dx = 0.0
	
      allocate(z(N-M), a(R))

      !il calcolo dell'errore avviene tramite il binned bootstrap poichè
      !a priori non si può sapere qual è il tempo di decorellazione migliore
      !da inserire nella simulazione, e per non sprecare tempo macchina,
      !bisogna tenere conto di una possibile correlazione frai i dati
	
      Dd = 2**16			!dimensione dei blocchi	
      nb=(N-M)/Dd			!numero dei blocchi
      do l = 1, R			!ciclo sui ricampionamenti
          do i = 1, nb
		
	      j = int(ran2()*(N-M) +1) !scelgo sito a caso
	      do g = 1, Dd		 	
	          z((i-1)*Dd+g) = x(mod(j+g-2,(N-M))+1) !ricampiono a blocchi	
	      enddo
			
	  enddo	
	  !calcolo della media e della varianza dei ricampionamenti		
	  a(l) = sum(z)/float(N-M)
			
      enddo
	
      !calcolo la media degli estimatori
      media_x  = sum(a)/float(R)
		
      do i = 1, R
          !calcolo scarto quadratico
          dx  = dx  + (a(i) - media_x )**2

      enddo
	
      dx = sqrt(dx/float(R - 1))	!prendo l'errore sul campione
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
