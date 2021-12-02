	program correlazione
	implicit none
	real*8, dimension(:), allocatable :: step, g, acc, P
	real*8 :: media, acceptance, c, t
	integer(16) :: i, j, N, M, k 
	
	open(unit=1, file="metrogausscorr.txt", status="unknown")
	open(unit=0, file="metrogauss.txt", status="old", action="read")
	
	read(0, *) N			!grandezza totale del campione
	read(0, *) M			!grandezza del campione analizzato
	
	media = 0
	k = 5000
	t = 0
	
	allocate(step(N), g(N), acc(N), P(M)) 
	
	do i = 1, N
		read(0, *) step(i), g(i), acc(i)				
	enddo
	
	acceptance=acc(N)/float(N)

	do i=1, M				!elimino la parte in cui la catena 
		P(i) = g((N - M) + i)		!potrebbe non essere termalizzata
	enddo
	
	media = sum(P)/float(M)
	
	do j = 1, k			!calcolo della autocorrelazione
		c=0			!va come O(N^2); ci vuole un po'
		do i = 1, M-j
			c = c + (P(i) - media)*(P(i + j) - media)
		end do
		c = c/float(M - j)
		t = t + c		!caclolo tempo di correlazione integrato
		
		write(1, *) c
	end do
	
	deallocate(step, g, acc, P)
	close(0)
	close(1)
	print *, media, acceptance, t
	
	end program correlazione
