# compiler to use and flags
CC := gcc
CC_FLAGS := -Wall -pedantic -std=c11 -ggdb

GNUPLOT = gnuplot
GNUPLOT_FLAGS = -pedantic

LIBS := -lm

# location of header files
INCLUDE := include
# location of binary files
BIN := bin
# location of .c files
SRC := src
# location of .gnu files
PLOT := plot

# We create a list of all the sources by looking for all the .cpp and .c files
SOURCES := $(wildcard $(SRC)/*.c)

# We create a list of object files by replacing the .cpp or .c extension with .o in the list of sources
OBJECTS := $(patsubst $(SRC)/%.c, $(BIN)/%.o, $(filter %.c, $(SOURCES))) 

# We need to tell the compiler where to find the headers
HEADERS := $(wildcard $(INCLUDE)/*.h)

#  .PHONY target specifies that all and clean are not real files, but are just targets that don't produce output files.
.PHONY: clean

$(BIN)/pendulum: $(BIN)/pendulum.o $(BIN)/fields.o $(BIN)/rk78.o
	$(CC) $(CC_FLAGS) $^ -o $@ $(LIBS)

$(BIN)/flow_sample: $(BIN)/flow_sample.o $(BIN)/flow.o $(BIN)/fields.o $(BIN)/rk78.o
	$(CC) $(CC_FLAGS) $^ -o $@ $(LIBS)

$(BIN)/lorenz: $(BIN)/lorenz.o $(BIN)/flow.o $(BIN)/fields.o $(BIN)/rk78.o
	$(CC) $(CC_FLAGS) $^ -o $@ $(LIBS)

$(BIN)/rtbps: $(BIN)/rtbps.o $(BIN)/flow.o $(BIN)/fields.o $(BIN)/rk78.o $(BIN)/rtbps.o
	$(CC) $(CC_FLAGS) $^ -o $@ $(LIBS)

$(BIN)/ex2: $(BIN)/flow.o $(BIN)/fields.o $(BIN)/rk78.o $(BIN)/ex2.o
	$(CC) $(CC_FLAGS) $^ -o $@ $(LIBS)

$(BIN)/ex43: $(BIN)/flow.o $(BIN)/rk78.o $(BIN)/ex43.o $(BIN)/fields.o $(BIN)/gauss.o $(BIN)/cmani.o
	$(CC) $(CC_FLAGS) $^ -o $@ $(LIBS)

$(BIN)/cmani_rtbp: $(BIN)/flow.o $(BIN)/rk78.o $(BIN)/cmani_rtbp.o $(BIN)/fields.o $(BIN)/gauss.o $(BIN)/cmani.o
	$(CC) $(CC_FLAGS) $^ -o $@ $(LIBS)

# We compile the .c files
$(BIN)/%.o: $(SRC)/%.c
	$(CC) $(CC_FLAGS) -I$(INCLUDE) -c $< -o $@ $(LIBS)

clean:
	rm -f $(BIN)/*