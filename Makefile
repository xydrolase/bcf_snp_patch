EXE_NAME=bcf_snp_patch

CPP=g++

CPPFLAGS=-O3 -std=c++11 -I$(HOME)/boost/include -I$(HOME)/sources/htslib
LDFLAGS=-L$(HOME)/boost/lib -L$(HOME)/lib -lhts -lboost_program_options

all: $(EXE_NAME)

SRCS=$(wildcard *.cpp)
OBJS=$(SRCS:.cpp=.o)

$(EXE_NAME): $(OBJS)
	$(CPP) -o $(EXE_NAME) $^ $(LDFLAGS)

.PHONY: clean

clean:
	rm $(EXE_NAME) *.o
