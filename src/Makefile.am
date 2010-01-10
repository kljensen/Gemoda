
# Define the program gemodas
bin_PROGRAMS = gemoda-s gemoda-p gemoda-r

# Make everything in this subdirectory first!
SUBDIRS = FastaSeqIO

# the sources that make gemodas are listed here
# fastaSeqIO is not listed because we link it in as a convenience library
gemoda_s_SOURCES = gemoda-s.c align.c words.c bitSet.c matrices.c \
		   newConv.c convll.c patStats.c
gemoda_p_SOURCES = gemoda-p.c bitSet.c newConv.c convll.c \
		   protAlign.c pdb.c patStats.c
gemoda_r_SOURCES = gemoda-r.c realIo.c realCompare.c bitSet.c \
		   newConv.c convll.c protAlign.c patStats.c

#  we require these header files
EXTRA_DIST = bitSet.h convll.h  spat.h matrixmap.h matrices.h \
	     pdb.h protAlign.h realIo.h FastaSeqIO/fastaSeqIO.h \
	     realCompare.h matdata.h patStats.h

INCLUDES = -I$(top_builddir) -I$(top_srcdir)

gemoda_s_LDADD = $(MATHLIB) FastaSeqIO/libfastaSeqIO.la
gemoda_p_LDADD = $(SCILIBS)
gemoda_r_LDADD = $(SCILIBS) FastaSeqIO/libfastaSeqIO.la
