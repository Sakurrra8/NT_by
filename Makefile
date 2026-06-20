PRGM		=edit
FC	        =ifort
CC			=g++ -fopenmp $(STATIC_LIBS)
STATIC_LIBS ?=
MPICC		=mpicc
FFLAGS		=-fast
CFLAGS 		=
CXXSTD ?= -std=c++17
BUILD ?= release
ARCH_FLAGS ?=

ifeq ($(BUILD),debug)
CPPFLAGS	=$(CXXSTD) -O0 -g
else ifeq ($(BUILD),profile)
CPPFLAGS	=$(CXXSTD) -O3 -g -DNDEBUG -fno-omit-frame-pointer $(ARCH_FLAGS)
LDFLAGS 	+= -fno-omit-frame-pointer
else
CPPFLAGS	=$(CXXSTD) -O3 -DNDEBUG $(ARCH_FLAGS)
endif

MCFLAGS		=-O2
LDFLAGS 	+=
NETCDFDIR	=
INCS		=
HDF5_DIR ?= /data/leuven/379/vsc37950/nt/opt/hdf5
INCLUDE		=-I ./include -I$(HDF5_DIR)/include -fopenmp
LIBS		=-lm -L$(HDF5_DIR)/lib -lhdf5_cpp -lhdf5
FCOMPILE 	= $(FC) $(FFLAGS)
CCOMPILE 	= $(CC) $(CPPFLAGS) $(CFLAGS)
LINKCC 		= $(CC) $(LDFLAGS)
LINKMPICC	= $(CC) $(LDFLAGS)

SRC_DIR = ./src
OBJ_DIR = ./obj
INC_DIR = ./include
BIN_DIR = ./bin

SRCS := $(wildcard $(SRC_DIR)/*.cpp)
OBJS := $(patsubst $(SRC_DIR)/%.cpp,%.o,$(SRCS))

VPATH = $(SRC_DIR)
vpath %.h $(INC_DIR)
vpath %.o $(OBJ_DIR)

all : $(BIN_DIR)/$(PRGM)

$(BIN_DIR)/$(PRGM):$(OBJS)
	$(LINKCC) $(addprefix $(OBJ_DIR)/, $(OBJS)) $(LIBS) -o $@
	
%.o:%.cpp Makefile $(wildcard $(INC_DIR)/*.h)
	$(CCOMPILE) -c $(INCS) $(INCLUDE) $< -o $(OBJ_DIR)/$@

%.o: %.f90
	$(FCOMPILE) -c $(INCS) $< -o $@

clean:
	rm -f $(OBJ_DIR)/*.o $(BIN_DIR)/$(PRGM) 
clean-data:
	rm -f *~ *.txt *.dat *.nc errfile outfile
