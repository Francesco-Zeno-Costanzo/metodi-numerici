      
      parameter (nlatt=20, nvol=nlatt**2)
	integer campo
	common/reticolo/campo(nlatt, nlatt)
	common/passi/nl(nlatt),ns(nlatt)
	common/vario/bext

C================================================================
C file che permette di cambiare più facilmente i parametri della
C simulazione facendolo una sola volta qui piuttosto che diverse
C volte all'interno del codice della simulazione; fare un ciclo
C sui reticoli allocando e deallocando il campo vorrebbe dire una
C simulazione molto lunga con il rischio di non poter verificare
C che tutto sia funzionando bene e quindi potenziale perdita di
C tempo
C================================================================
