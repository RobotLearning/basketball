HOST=$(shell hostname)
LIBPATH=/usr/local/lib
DIR=$(HOME)/basketball
HEADER=$(DIR)/include
LIBDIR=$(DIR)/lib
CC=g++
LIBS=-larmadillo -lm -lboost_program_options
FLAGS=-I$(HEADER) -pthread -std=c++11
DEBUG=-DDEBUG -g -Wall -pedantic -Wextra -Weffc++
SHARED_OBJ=$(LIBDIR)/libbasket.so
OBJDIR=$(DIR)/obj
SRCS=$(wildcard $(DIR)/src/*.cpp $(DIR)/src/*.c)
OBJS=$(addsuffix .o,$(addprefix $(OBJDIR)/,$(basename $(notdir $(SRCS)))))
#$(info $$OBJS is [${OBJS}])

all: $(SHARED_OBJ)

$(SHARED_OBJ) : $(OBJS)
	$(CC) -shared $(FLAGS) -o $@ $^ -L$(LIBPATH) $(LIBS)

$(OBJDIR)/player.o : $(DIR)/src/player.cpp $(DIR)/include/player.hpp $(DIR)/include/optim.h $(DIR)/include/kinematics.h
	$(CC) -c -fPIC $(FLAGS) -o $@ $<
	
$(OBJDIR)/utils.o : $(DIR)/src/utils.c $(DIR)/include/utils.h
	gcc -c -fPIC -I$(HEADER) -pthread -o $@ $<	

$(OBJDIR)/optim.o : $(DIR)/src/optim.cpp $(DIR)/include/kinematics.h $(DIR)/include/optim.h $(DIR)/include/utils.h
	$(CC) -c -fPIC $(FLAGS) -o $@ $<
	
$(OBJDIR)/kalman.o : $(DIR)/src/kalman.cpp $(DIR)/include/kalman.h
	$(CC) -c -fPIC $(FLAGS) -o $@ $<
	
$(OBJDIR)/extkalman.o : $(DIR)/src/extkalman.cpp $(DIR)/include/kalman.h
	$(CC) -c -fPIC $(FLAGS) -o $@ $<
	
$(OBJDIR)/sl_interface.o : $(DIR)/src/sl_interface.cpp $(DIR)/include/kalman.h $(DIR)/include/optim.h $(DIR)/include/kinematics.h
	$(CC) -c -fPIC $(FLAGS) -o $@ $<
	
$(OBJDIR)/kinematics.o : $(DIR)/src/kinematics.cpp $(DIR)/include/kinematics.h
	$(CC) -c -fPIC $(FLAGS) -o $@ $<
	
clean:
	rm -rf obj/*.o lib/*.so
	
##### ALL TESTS ARE INCLUDED HERE
test:
	$(CC) $(FLAGS) test/test_optim.cpp -o unit_tests \
	               $(SHARED_OBJ) -L$(LIBPATH) $(LIBS) $(LIBPATH)/libboost_unit_test_framework.a -lnlopt
				

.PHONY: all test clean
