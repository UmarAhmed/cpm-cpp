
# List directories
SRC_DIR := src
OBJ_DIR := obj
BIN_DIR := bin

# Executable name
EXE := $(BIN_DIR)/main

# List source files
SRC := $(wildcard $(SRC_DIR)/*.cc)

# List object files
OBJ := $(SRC:$(SRC_DIR)/%.cc=$(OBJ_DIR)/%.o)


# Flags
CC := g++
CPPFLAGS := -I include -MMD -MP
CXXFLAGS := -std=c++14
LIB := -larmadillo -O2
PLOT := -DWITHOUT_NUMPY -I/usr/include/python3.8/ -lpython3.8


.PHONY: all clean

all: $(EXE)


$(EXE): $(OBJ) | $(BIN_DIR)
	$(CC) $^ $(PLOT) $(LIB) -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc | $(OBJ_DIR)
	$(CC) $(CPPFLAGS) $(CXXFLAGS) $(PLOT) $(LIB) -c $< -o $@

$(BIN_DIR) $(OBJ_DIR):
	mkdir -p $@

clean:
	@$(RM) -rv $(BIN_DIR) $(OBJ_DIR)

-include $(OBJ:.o=.d)
