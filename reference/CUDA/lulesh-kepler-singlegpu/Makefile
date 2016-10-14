NVCC		= nvcc
FLAGS		= -arch=sm_35
DFLAGS	= -G -g -lineinfo
RFLAGS 	= -O3 -DNDEBUG 

#SILO_INCLUDES := /usr/local/silo-4.8/include
#SILO_LIBS := /usr/local/silo-4.8/lib

LINKFLAGS = 
#LINKFLAGS += -L$(SILO_LIBS) -lsilo

#INC_SILO:= -I$(SILO_INCLUDES)

all: release 

debug: LINKFLAGS += -G -g

release: 	FLAGS += $(RFLAGS)
debug: 		FLAGS += $(DFLAGS)

release: lulesh
debug: lulesh

lulesh: allocator.o lulesh.o
	$(NVCC) $(LINKFLAGS) allocator.o lulesh.o -o lulesh

allocator.o: allocator.cu vector.h
	$(NVCC) $(FLAGS) allocator.cu -I ./ -c -o allocator.o

lulesh.o: lulesh.cu util.h vector.h texture_objAPI.h allocator.h
	$(NVCC) $(FLAGS) lulesh.cu -I ./  $(INC_SILO) -c -o lulesh.o

clean: 
	rm -rf allocator.o  lulesh.o lulesh




