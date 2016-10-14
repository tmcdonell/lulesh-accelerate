# This is the Loci base directory.
LOCI_BASE ?= /home/lush/Work/Loci/OBJ
# List your object files here
OBJS  =  main.o utils.o elem.o node.o timestep.o varsFileInputs.o grid_io.o 
# List the name of your compiled program here
TARGET = hydro

#############################################################################
# No changes needed below this line

include $(LOCI_BASE)/Loci.conf

#COPTLESS = -O0 -g
#COPT=$(COPTLESS)

default:
	$(MAKE) $(TARGET)

$(TARGET): $(OBJS)
	$(LD) -o $(TARGET) $(OBJS) $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)


clean:
	rm -fr $(OBJS) $(TARGET) 

# Junk files that are created while editing and running cases
JUNK = $(wildcard *~) $(wildcard crash_dump.*)  core debug output $(wildcard .schedule*)
# ".cc" files created from .loci files
LOCI_INTERMEDIATE_FILES = $(subst .loci,.cc, $(wildcard *.loci) )

distclean:
	rm -fr $(OBJS) $(TARGET) $(JUNK) $(LOCI_INTERMEDIATE_FILES) $(DEPEND_FILES) 

DEPEND_FILES=$(subst .o,.d,$(OBJS))

#include automatically generated dependencies
-include $(DEPEND_FILES)
