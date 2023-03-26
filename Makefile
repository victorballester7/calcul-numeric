CC = gcc                       					   # compiler to use
CC_FLAGS = -Wall -pedantic -std=c99	-ggdb  # flags

GNUPLOT = gnuplot
GNUPLOT_FLAGS = -pedantic

# location of binary files
BIN := bin
# location of .c files
SRC := src

# location of .gnu files
PLOT := plot

# location of header files
INCLUDE := include

LIBRARIES := -lm
EXECUTABLE := main 

H_MIN=0.01
H_MAX=0.05
TOL=0.000001

SIGMA=10
RHO=28
BETA=2.6
TF=50
NT=5000

COMPILE=$(CC) $(CC_FLAGS) -I$(INCLUDE) $^ -o $(BIN)/$@ $(LIBRARIES)
EXECUTE=./$(BIN)/$^
# Rules are of the form:
# target: prerequisites
# 	 recipe

# This defines all the targets that are not files.
.PHONY: plotPendulum runPendulum runFlowSample plotLorenz runLorenz rf_pendulum  rf_pendulum flowSample lorenz_int clean 

# default target to be executed is the first encountered. To change it, uncomment the following line and replace all with another target.
# .DEFAULT_GOAL := all
# all: $(BIN)/$(EXECUTABLE) 

# # call it with 'make run' because it is not the default one.
# run: clean all
# 	@./$(BIN)/$(EXECUTABLE)

plotPendulum: runPendulum
	@$(GNUPLOT) $(GNUPLOT_FLAGS) $(PLOT)/pendulum.gnu

runPendulum: rf_pendulum
	@$(EXECUTE) $(H_MIN) $(H_MAX) $(TOL)

runFlowSample: flowSample
	@$(EXECUTE)

plotLorenz: runLorenz
	@$(GNUPLOT) $(GNUPLOT_FLAGS) $(PLOT)/lorenz.gnu

runLorenz: lorenz_int
	@$(EXECUTE) $(SIGMA) $(RHO) $(BETA) $(TF) $(NT)


# $^ evaluates to $(SRC)/rf_pendulum.c $(SRC)/pendulum.c $(SRC)/rk78.c
# S@ evaluates to rf_pendulum
# The option -I is will tell the compiler where to find the header file.
rf_pendulum: $(SRC)/rf_pendulum.c $(SRC)/fields.c $(SRC)/rk78.c
	@$(COMPILE)

flowSample: $(SRC)/flow_sample.c $(SRC)/flow.c $(SRC)/fields.c $(SRC)/rk78.c
	@$(COMPILE)

lorenz_int: $(SRC)/lorenz_int.c $(SRC)/flow.c $(SRC)/fields.c $(SRC)/rk78.c
	@$(COMPILE)

rtbps_int: $(SRC)/rtbps_int.c $(SRC)/flow.c $(SRC)/fields.c $(SRC)/rk78.c $(SRC)/rtbps.c
	@$(COMPILE)

ex2: $(SRC)/flow.c $(SRC)/fields.c $(SRC)/rk78.c $(SRC)/ex2.c
	@$(COMPILE)

clean:
	@-rm -f $(BIN)/*