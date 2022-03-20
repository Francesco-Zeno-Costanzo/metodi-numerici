      program campo2d
      	
      include "par2d.f"
C===============================================================
C file che permette di cambiare più facilmente i parametri della
C simulazione facendolo una sola volta piuttosto che diverse
C volte all'interno del codice 
C===============================================================

      call cpu_time(start)
      call ranstart
      
      open(1,file='input.txt',status='old')
      open(2,file='dati/14mo.dat',status='unknown')

      read(1,*) N
      read(1,*) M
      read(1,*) i_dec      
      read(1,*) massa
      
      misure = N + M
      
      write(2,*) misure 
      write(2,*) M
      write(2,*) nx, nt
      
      
      massa2 = massa**2
      massa2p4 = massa**2 + 4.0 !ricorrente
      
      call bordi()
      call init()

      do iter = 1, misure
          do j = 1, i_dec
              call heatbath()
              call overrelax()
              call overrelax()
              call overrelax()
              call overrelax()
          enddo


          call misurazioni(fm2, ds2, dt2)
          write(2,*) fm2, ds2, dt2

      enddo 
     
      call ranfinish
      
      call cpu_time(finish)
      print'("tempo di esecuzione =", f16.8," secondi.")',finish-start

      end program campo2d

C============================================================================

      subroutine bordi() 		
      
      include "par2d.f"
      
      integer, dimension(2) :: c
      c = (/ nx, nt /)
      
      do j = 1, 2
          do i = 1, c(j)
              nl(i, j) = i + 1
              ns(i, j) = i - 1
          enddo
          nl(c(j), j) = 1
          ns(1, j) = c(j)
      enddo
      
      return
      end
      
c=================================================================

      subroutine init()
      
      include "par2d.f"
      	
      do i_x = 1, nx
          do i_t = 1, nt
              x = ran2()
              campo(i_x, i_t) = 1.0 - 2.0*x
          enddo
      enddo
      
      return
      end

C==================================================================
      subroutine heatbath()
      
      include "par2d.f"
      
      pi = 4.0*ATAN(1.0) ! pigreco
      
      do i_x = 1, nx
          do i_t = 1, nt
              f = 0.0 
              phi = campo(i_x, i_t)
              f = campo(nl(i_x, 1), i_t)+
     &            campo(ns(i_x, 1), i_t)+
     &            campo(i_x, nl(i_t, 2))+
     &            campo(i_x, ns(i_t, 2))                


              sigma2 = 1.0/massa2p4
              aver = f*sigma2

              x = sqrt(sigma2)*sqrt(-2.0*log(1-ran2())) 
              y = x*cos(2.0*pi*ran2()) + aver
              campo(i_x, i_t) = y
          enddo
      enddo


      return
      end



C==================================================================

      subroutine overrelax()
      
      include "par2d.f"
      
      do i_x = 1, nx
          do i_t = 1, nt
              f = 0.0 
              phi = campo(i_x, i_t)
              f = campo(nl(i_x, 1), i_t)+
     &            campo(ns(i_x, 1), i_t)+
     &            campo(i_x, nl(i_t, 2))+
     &            campo(i_x, ns(i_t, 2))             

              aver = f/massa2p4
                 
              campo(i_x, i_t) = 2.0*aver - phi
              
          enddo
      enddo


      return
      end
      
c=================================================================

      subroutine misurazioni(fm2, ds2, dt2)
      
      include "par2d.f"

      fm2 = 0.0
      ds2 = 0.0
      dt2 = 0.0
      do i_x = 1, nx
          do i_t = 1, nt
              phi = campo(i_x, i_t)
              f_s = campo(nl(i_x,1), i_t)
              f_t = campo(i_x, nl(i_t,2))

              fm2 = fm2 + massa2*phi**2 
              ds2 = ds2 - 2.0*phi*f_s  + 2.0*phi**2
              dt2 = dt2 - 2.0*phi*f_t  + 2.0*phi**2
          enddo
      enddo

      fm2 = fm2/float(nvol)
      ds2 = ds2/float(nvol)
      dt2 = dt2/float(nvol)

      return
      end

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


