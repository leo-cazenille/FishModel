VERBOSE = 0

# Verbose commands
AT_0 := @
AT_1 := 
AT = $(AT_$(VERBOSE))
ECHO_0 := @echo
ECHO_1 := @true
ECHO = $(ECHO_$(VERBOSE))

SRC_CC = $(wildcard *.cpp **/*.cpp ../src/random.cpp)
SRC_CC := $(filter-out pythonInterface.cpp, $(SRC_CC))
SRC_CC := $(filter-out testModel.cpp, $(SRC_CC))
OBJ = $(SRC_CC:.cpp=.o)
LDFLAGS = $(shell pkg-config --libs opencv) -fopenmp -fPIC
#CCFLAGS = -O3 -std=c++11 -fopenmp $(shell pkg-config --cflags opencv) -I../src/ #-DENABLE_OPENMP
CCFLAGS = -O0 -g -std=c++11 -fopenmp $(shell pkg-config --cflags opencv) -I../src/ `python-config --includes` -fPIC
TARGETEXEC := model

CPPC = g++
RM = rm -f

all: executable pythonlib

executable: $(TARGETEXEC)


$(TARGETEXEC): $(OBJ) testModel.o
	$(ECHO) "[LD]	" $@
	$(AT)$(CPPC) -o $@ $^ $(LDFLAGS)


%.o: %.cpp
	$(ECHO) "[CPP]	" $@
	$(AT)$(CPPC) -o $@ -c $< $(CCFLAGS) $(INCLUDE)

pythonlib: $(OBJ) pythonInterface.o
	$(ECHO) "[SHARED]" $@
	$(AT)$(CPPC) -shared -o $(TARGETEXEC).so $^ `python-config --ldflags` $(LDFLAGS) 

clean:
	$(AT)$(RM) $(OBJ) pythonInterface.o testModel.o

mrproper: clean
	$(AT)$(RM) $(TARGETEXEC) $(TARGETEXEC).so


