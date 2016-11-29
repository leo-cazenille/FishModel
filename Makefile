VERBOSE = 0

# Verbose commands
AT_0 := @
AT_1 := 
AT = $(AT_$(VERBOSE))
ECHO_0 := @echo
ECHO_1 := @true
ECHO = $(ECHO_$(VERBOSE))

SRC_CC = $(wildcard *.cpp **/*.cpp ../src/random.cpp)
OBJ = $(SRC_CC:.cpp=.o)
LDFLAGS = $(shell pkg-config --libs opencv) -fopenmp
#CCFLAGS = -O3 -std=c++11 -fopenmp $(shell pkg-config --cflags opencv) -I../src/ #-DENABLE_OPENMP
CCFLAGS = -O0 -g -std=c++11 -fopenmp $(shell pkg-config --cflags opencv) -I../src/
TARGETEXEC := model

CPPC = g++
RM = rm -f


executable: $(TARGETEXEC)


$(TARGETEXEC): $(OBJ)
	$(ECHO) "[LD]	" $@
	$(CPPC) -o $@ $^ $(LDFLAGS)


%.o: %.cpp
	$(ECHO) "[CPP]	" $@
	$(AT)$(CPPC) -o $@ -c $< $(CCFLAGS) $(INCLUDE)


clean:
	$(AT)$(RM) $(OBJ)

mrproper: clean
	$(AT)$(RM) $(TARGETEXEC)


