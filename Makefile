DIR_INC := ./inc
DIR_SRC := ./src
DIR_OBJ := ./obj

PREFIX ?= /usr/local
BINDIR ?= $(PREFIX)/bin
INCLUDE_DIRS ?=
LIBRARY_DIRS ?=

SRC := $(wildcard ${DIR_SRC}/*.cpp)
OBJ := $(patsubst %.cpp,${DIR_OBJ}/%.o,$(notdir ${SRC}))

TARGET := seq2sat

BIN_TARGET := ${TARGET}

CXX ?= g++
CXXFLAGS := -std=c++11 -g -O3 -I${DIR_INC} $(foreach includedir,$(INCLUDE_DIRS),-I$(includedir)) ${CXXFLAGS}
LIBS := -lz -lpthread
LD_FLAGS := $(foreach librarydir,$(LIBRARY_DIRS),-L$(librarydir)) $(LIBS) $(LD_FLAGS)

# Default target
all: ${BIN_TARGET}

${BIN_TARGET}:${OBJ}
	$(CXX) $(OBJ) -o $@ $(LD_FLAGS)

# Debug target
debug: CXXFLAGS += $(DEBUG_CXXFLAGS)
debug: ${BIN_TARGET}

${DIR_OBJ}/%.o:${DIR_SRC}/%.cpp make_obj_dir
	$(CXX) -c $< -o $@ $(CXXFLAGS)

make_obj_dir:
	@if test ! -d $(DIR_OBJ) ; \
	then \
		mkdir $(DIR_OBJ) ; \
	fi
	
.PHONY:clean
clean:
	@if test -d $(DIR_OBJ) ; \
	then \
		find $(DIR_OBJ) -name *.o -delete; \
	fi
	@if test -e $(TARGET) ; \
	then \
		rm $(TARGET) ; \
	fi

install:
	install $(TARGET) $(BINDIR)/$(TARGET)
	@echo "Installed."