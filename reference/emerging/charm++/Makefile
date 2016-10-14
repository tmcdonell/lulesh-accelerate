CHARMBASE = $(HOME)/charm
OPTS = -O3
CHARMC = $(CHARMBASE)/bin/charmc $(OPTS)

OUT := lulesh
C_CI := $(wildcard *.ci)
C_SRCS := $(wildcard *.cc)
C_OBJS := $(C_SRCS:.cc=.o)
DECL := $(C_CI:.ci=.decl.h)
DEF := $(C_CI:.ci=.def.h)
DECL_DEF := $(DECL) $(DEF)

.PHONY: all projections loadbalance clean distclean

all: $(OUT)

$(OUT): $(C_SRCS) $(DECL_DEF)
	$(CHARMC) -language charm++ $(C_SRCS) -lm -o $(OUT)

projections: $(C_SRCS) $(DECL_DEF)
	$(CHARMC) -tracemode projections -tracemode summary -module CommonLBs -language charm++ $(C_SRCS) -lm -o $(OUT)

loadbalance: $(C_SRCS) $(DECL_DEF)
	$(CHARMC) -language charm++ -module CommonLBs $(C_SRCS) -lm -o $(OUT)

$(DECL_DEF): $(C_CI)
	$(CHARMC) $(C_CI)

clean:
	rm -f *.decl.h *.def.h *.o charmrun $(OUT)

distclean: clean
