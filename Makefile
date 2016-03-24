
MEX=$(HOME)/bin/mex

all: mexLoadAgilentTraces mexLoadAgilentParams

mexLoadAgilentTraces:
	mex -g -v  -largeArrayDims mexLoadAgilentTraces.c mexCommons.c varian.c

mexLoadAgilentParams:
	mex -g -v  -largeArrayDims mexLoadAgilentTraces.c mexCommons.c varian.c

clean:
	rm *.o *.mexmaci64
