PRGM		=edit
FC	        =ifort
CC			=g++ -fopenmp -static-libstdc++ -static-libgcc
MPICC		=mpicc
FFLAGS		=-fast
CFLAGS 		=
CPPFLAGS	=-std=c++17 -g 
MCFLAGS		=-O2
LDFLAGS 	=
NETCDFDIR	=
INCS		=
INCLUDE		=-I ./include -fopenmp
LIBS		=-lm -L/data/leuven/379/vsc37950/nt/opt/hdf5/lib -lhdf5_cpp -lhdf5
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
	
%.o:%.cpp
	$(CCOMPILE) -c $(INCS) $(INCLUDE) $< -o $(OBJ_DIR)/$@

%.o: %.f90
	$(FCOMPILE) -c $(INCS) $< -o $@

clean:
	rm -f $(OBJ_DIR)/*.o $(BIN_DIR)/$(PRGM) 
clean-data:
	rm -f *~ *.txt *.dat *.nc errfile outfile

