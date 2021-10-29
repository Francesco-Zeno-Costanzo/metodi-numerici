	program ising
	parameter (nlatt=50, nvol=nlatt**2)
	integer campo
	common/reticolo/campo(nlatt, nlatt)
	common/vario/bext
	
	call cpu_time(start)
	call ranstart

	open(1, file='input.txt',status='old')		!file coi parametri
	open(2, file='dati50.dat',status='unknown')	!file coi risultati
	
	read(1,*) misure         !numero di misure
        read(1,*) i_dec          !updating fra una misura e l'altra
        read(1,*) bext      	 !valore del campo esterno
        read(1,*) bmin 		 !temperatura inversa minima
        read(1,*) bmax		 !temperatura inversa missima
        read(1,*) npassi	 !numero di temperature
       
	write(2, *) misure	 !scrivo quantità che seriviranno
	write(2, *) nvol	 !nell'analisi
	write(2, *) npassi
	
	call bordi()		 !condizioni di bordo continuo
	call init()		 !inizializzo la matrice
	
	!la scrittura avviene su un unico file che poi verra letto a blocchi
	!ogni blocco corrisponde ad una temperatura diversa
	!ci sono npassi blocchi, ciascuno lungo misure
	
	do k=1, npassi		 !ciclo sulle temperature
	
		beta = bmin + (k)*(bmax-bmin)/float(npassi)
		
		do i=1, misure 	!ciclo sulle misure a fissa temperatura
		
			do j=1, i_dec
				call metropolis(beta)   !decorrela la matrice
			enddo
			
			call magnetizzazione(mag)	!misurazione delle osservabili
			call energia(ene)
			write(2,*) abs(mag), ene	!salvataggio risultati
		enddo
	enddo
	
	call ranfinish
	
	call cpu_time(finish)
	print '("tempo di esecuzione= ", f16.8," secondi.")', finish-start
	
	end program ising

C============================================================================

	subroutine bordi() 		
	   
	parameter (nlatt = 50,nvol = nlatt**2)
      	common/passi/nl(nlatt),ns(nlatt)
      	
      	do i =1, nlatt
      		nl(i) = i + 1		!codizioni per ogni sito
      		ns(i) = i - 1 
      	enddo 
      	nl(nlatt) = 1			!condizioni periodiche per i siti esterni
      	ns(1) = nlatt
      	
      	return
      	end
 
C============================================================================
      	
      	subroutine init()
      	
      	parameter (nlatt = 50,nvol = nlatt**2)
      	common/reticolo/campo(nlatt, nlatt)
      	
      	do i=1, nlatt
      		do j=1, nlatt
      			x=ran2()
      			if(x<0.5) then
      				campo(i,j) = 1
      			else
      				campo(i,j) = -1
      			endif
      		enddo
      	enddo
      	return
      	end

C============================================================================

	subroutine metropolis(beta)
	
	parameter (nlatt = 50,nvol = nlatt**2)
	common/reticolo/campo(nlatt, nlatt)
	common/passi/nl(nlatt),ns(nlatt)
	common/vario/bext
	
	do iv = 1, nvol				!ciclo su tutti i siti
		
		i = int(nlatt*ran2()+1.0)	!scelta random di un sito
		j = int(nlatt*ran2()+1.0)
		
		ip = nl(i)			!calcolo dei primi vicini
		im = ns(i)
		jp = nl(j)
		jm = ns(j)
		
		F = campo(i, jp) + campo(i, jm) + 
     & 		    campo(ip, j) + campo(im, j)	!sommo i primi vicini
     		F = beta*(F + bext)		!aggiungo eventuale campo esterno
     		
     		ispin = campo(i, j)
     		
     		p =  exp(-2.0*ispin*F)		!probabilità di accettare la mossa
     		
     		x = ran2()			!numero casuale per il test
     		
     		if(x < p) then			!test di accettanza
     			campo(i, j) = -ispin	!se F è negativo il test è
     		endif				!passato di Default
	
	enddo
	return
	end
	

C============================================================================

	subroutine magnetizzazione(mag)
	
	parameter (nlatt = 50,nvol = nlatt**2)
	common/reticolo/campo(nlatt, nlatt)
	
	mag = 0					!inzializzo la variabile
	do i=1, nlatt				!ciclo su tutto il reticolo
		do j=1, nlatt			!e sommo ogni sito
			mag = mag + campo(i, j)
		enddo
	enddo
	return
	end

C============================================================================

	subroutine energia(ene)
	
	parameter (nlatt = 50,nvol = nlatt**2)
	common/reticolo/campo(nlatt, nlatt)
	common/passi/nl(nlatt),ns(nlatt)
	common/vario/bext
	
	ene=0					!inzializzo la variabile
	do i = 1, nlatt				!ciclo su tutto il reticolo
		do j =1, nlatt
			ip = nl(i)		!calcolo primi vicini
			im = ns(i)		
			jp = nl(j)
			jm = ns(j)
			
			F = campo(i, jp) + campo(i, jm) + 
     & 		    	    campo(ip, j) + campo(im, j)
     			ene = ene - 0.5*F*campo(i, j)	!0.5 per non sovracontare
     			ene = ene - bext*campo(i, j)	!eventuale campo esterno
		enddo
	enddo
	
	return
	end


C============================================================================


c============================================================================
c  RANDOM NUMBER GENERATOR: standard ran2 from numerical recipes
c============================================================================
      function ran2()
      implicit real*4 (a-h,o-z)
      implicit integer*4 (i-n)
      integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      real ran2,am,eps,rnmx
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
c=============================================================================

