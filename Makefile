HOST=$(shell hostname)
LIBPATH=/usr/local/lib
DIR=$(HOME)/basketball
HEADER=$(DIR)/include
LIBDIR=$(DIR)/lib
CC=g++
LIBS=-larmadillo -lm
FLAGS=-I$(HEADER) -pthread -std=c++11
RELEASE=-O3 -DNDEBUG
DEBUG=-DDEBUG -g -Wall -pedantic -Wextra #-Weffc++	
SHARED_OBJ=$(LIBDIR)/libbasket.so
OBJDIR=$(DIR)/obj
SRCS=$(wildcard $(DIR)/src/*.cpp $(DIR)/src/*.c)
OBJS=$(addsuffix .o,$(addprefix $(OBJDIR)/,$(basename $(notdir $(SRCS)))))
#$(info $$OBJS is [${OBJS}])

# for compiling everything, release and debug modes
release: FLAGS += $(RELEASE)
release: all
debug: FLAGS += $(DEBUG)
debug: all


all: $(SHARED_OBJ)

$(SHARED_OBJ) : $(OBJS)
	$(CC) -shared $(FLAGS) -o $@ $^ -L$(LIBPATH) $(LIBS)

$(OBJDIR)/player.o : $(DIR)/src/player.cpp $(HEADER)/player.hpp $(HEADER)/optim.h $(HEADER)/kinematics.h
	$(CC) -c -fPIC $(FLAGS) -o $@ $<
	
$(OBJDIR)/kalman.o : $(DIR)/src/kalman.cpp $(HEADER)/kalman.h
	$(CC) -c -fPIC $(FLAGS) -o $@ $<
	
$(OBJDIR)/extkalman.o : $(DIR)/src/extkalman.cpp $(HEADER)/kalman.h
	$(CC) -c -fPIC $(FLAGS) -o $@ $<
	
$(OBJDIR)/sl_basketball_interface.o : $(DIR)/src/sl_basketball_interface.cpp $(HEADER)/kalman.h $(HEADER)/optim.h $(HEADER)/kinematics.h \
							$(HEADER)/constants.h $(HEADER)/sl_basketball_interface.h
	$(CC) -c -fPIC $(FLAGS) -o $@ $<
	
$(OBJDIR)/optim.o : $(DIR)/src/optim.cpp $(HEADER)/kinematics.h $(HEADER)/optim.h $(HEADER)/utils.h
	$(CC) -c -fPIC $(FLAGS) -o $@ $<
	
$(OBJDIR)/kinematics.o : $(DIR)/src/kinematics.cpp $(HEADER)/kinematics.h
	$(CC) -c -fPIC $(FLAGS) -o $@ $<
	
$(OBJDIR)/utils.o : $(DIR)/src/utils.c $(HEADER)/utils.h
	$(CC) -c -fPIC $(FLAGS) -o $@ $<
	
clean:
	rm -rf obj/*.o lib/*.so
	
##### ALL TESTS ARE INCLUDED HERE
test:
	$(CC) $(FLAGS) $(DEBUG) test/test_optim.cpp -o unit_tests \
	               $(SHARED_OBJ) -L$(LIBPATH) $(LIBS) $(LIBPATH)/libboost_unit_test_framework.a -lnlopt
				

.PHONY: all test clean
