target = mytest
CC = gcc
CPP = g++
FROTRAN = ifort
CFLAGS = -O2 -I./include
LIBS = 
LINKER = $(CPP)
LINKFLAG = -I.

# source files
source   := $(wildcard ./src/*.cpp)
srcfiles := $(notdir $(source))

# object files
objects	 := $(patsubst %.cpp, %.o, $(source))

all: $(target)

$(target): $(objects)
	@echo "linking the program"
	$(LINKER) $(LINKFLAG) $(objects) -o$(target) 
	@echo "linking is done!"


# compile rules

.cpp.o:
	$(CPP) $(CFLAGS) -c $< -o $@
.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	@rm *.o 
