set encoding utf8; 
# set terminal qt size 1750, 850;

# set the colors for each block
set linetype 1 linecolor rgb "red"
set linetype 2 linecolor rgb "blue"
set linetype 3 linecolor rgb "green"
set linetype 4 linecolor rgb "orange"
set linetype 5 linecolor rgb "purple"
set linetype 6 linecolor rgb "gray"
set linetype 7 linecolor rgb "black"

set view equal xyz
set xlabel 'x'; 
set ylabel 'y'; 
set zlabel 'z';
set key outside bottom center; 

# define the data file and block separator
datafile = 'data/output_lorenz.txt'

# count the number of blocks in the file
stats datafile using 2:3 nooutput
nblocks = STATS_blocks

splot for[i=0:nblocks-1] datafile index i using 2:3:4 title columnheader(1) with lines linecolor i+1