	program blocking
	implicit none
	real*8, dimension(:), allocatable :: g, t
	integer(16) :: i, j, N, M
	real*8 :: media, a, q
	integer :: P
	
	open(unit=1, file="metrogausserrdb.txt", status="unknown")
	open(unit=0, file="metrogauss.txt", status="old", action="read")
	
	read(0, *) N			!grandezza totale del campione
	read(0, *) M			!grandezza del campione analizzato
	
	P=20
	allocate(g(N), t(M))	
	media=0
	
	do i = 1, N			!prima e terza colonna non interessano
		read(0, *) q,  g(i), q	!per l'analisi, quindi sono messi nella
	enddo				!variabile sempre sovrascritta e mai più usata
	
	do i=1, M			!elimino la parte in cui la catena 
		t(i) = g((N - M) + i)	!potrebbe non essere termalizzata
	enddo
	
	media=sum(t)/float(M) 
	print*, media
	
	do j=1, P
		a = 0.
		do i=1, M/2**j				!i primi 2**j elementi vengo sempre
			t(i) = (t(2*i- 1) + t(2*i))/2.	!sovascritti per poter iterare
			a = a + (t(i) - media)**2	!poichè solo a quelli si accede
		end do
		a=sqrt(a/((M/2.0**j)*(M/2.0**j - 1)))	
		write (1, *) a, 2**j
	end do
	
	deallocate(g, t)
	
	close(0)
	close(1)
	
	end program blocking
