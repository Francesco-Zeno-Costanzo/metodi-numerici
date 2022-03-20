      
      parameter(nx=160, nt=14, nvol=nx*nt)
      parameter(nlatt = nx)
      
      real massa, massa2, massa2p4
      real campo
      
      common/varie/massa,massa2,massa2p4 
      common/reticolo/campo(nx, nt) 
      common/passi/nl(nlatt, 2),ns(nlatt, 2)

C===============================================================
C file che permette di cambiare più facilmente i parametri della
C simulazione facendolo una sola volta qui piuttosto che diverse
C volte all'interno del codice della simulazione
C===============================================================
