all: deepgp$(SHLIB_EXT)

deepgp.so:
	${R_HOME}/bin/R CMD SHLIB -o deepgp$(SHLIB_EXT) *.c

clean:
	rm -rf *.o *.so
