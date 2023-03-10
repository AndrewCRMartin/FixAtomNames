OFILES = fixlabels.o FixAtomLabels.o
LIBS   = -lbiop -lgen -lm -lxml2
LIBDIR = $(HOME)/lib
INCDIR = $(HOME)/include
COPT   = -O3  -I $(INCDIR)
LOPT   = -L $(LIBDIR)

fixlabels : $(OFILES) 
	cc $(LOPT) -o $@ $(OFILES) $(LIBS)

.c.o : 
	cc $(COPT) -c -o $@ $<
