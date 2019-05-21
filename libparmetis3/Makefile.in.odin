
# Which compiler to use
CC = mpicc


# What optimization level to use
OPTFLAGS = -O3 

# Include directories for the compiler
INCDIR = 

# What options to be used by the compiler
COPTIONS = 

# Which loader to use
LD = mpicc

# In which directories to look for any additional libraries
LIBDIR =

# What additional libraries to link the programs with (eg., -lmpi)
#XTRALIBS = -lefence
#XTRALIBS = -ldmalloc

# What archiving to use
AR = ar rv

# What to use for indexing the archive
#RANLIB = ranlib
RANLIB = ar -ts

VERNUM = 
