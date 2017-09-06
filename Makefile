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
OBJ_DIR=$(DIR)/obj
SRCS=$(wildcard $(DIR)/src/*.cpp)
OBJS=$(addsuffix .o,$(addprefix $(OBJ_DIR)/,$(basename) $(notdir $(SRC))))

all: $(SHARED_OBJ)

$(SHARED_OBJ) : $(OBJS)
	$(CC) -shared $(FLAGS) -o $@ $^ -L$(LIBPATH) $(LIBS)

#ADD HEADER PREREQUISITES!
$(OBJ_DIR)/%.o : $(DIR)/src/%.cpp
	$(CC) -c -fPIC $(FLAGS) -o $@ $<

clean:
	rm -rf obj/*.o lib/*.so

.PHONY: all test clean
