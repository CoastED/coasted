input: clean
	cd src; gmsh -2 sphere.geo; ../../../bin/gmsh2triangle -s sphere.msh

clean:
	rm -f  fluidity.* *.ele *.face *.node *.vtu *.poly *.edge *.s *.d.1 *.stat \
	matrixdump matrixdump.info
	rm -f src/*.msh
