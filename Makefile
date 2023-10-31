all:
	f2py -c --opt=-O3 -m --quiet rffw rffw.f90 hrftn.f sacsubf.f > rffw.log
	pip install -e .
clean:
	rm drive.exe
	rm Recfun*
