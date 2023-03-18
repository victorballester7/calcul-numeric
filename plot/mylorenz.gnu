set encoding utf8; 
set terminal qt size 1750, 850;

# set the colors for each block
set linetype 1 linecolor rgb "red"
set linetype 2 linecolor rgb "blue"
set linetype 3 linecolor rgb "green"
set linetype 4 linecolor rgb "orange"
set linetype 5 linecolor rgb "purple"
set linetype 6 linecolor rgb "gray"
set linetype 7 linecolor rgb "black"

set view projection xz;
set xlabel 'z'; 
set ylabel 'y'; 
set zlabel 'x';
# set key outside bottom center; 

# set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb "#black" behind
# set style fill transparent solid 0.5

# set title ""
# unset xlabel
# unset ylabel
# unset border
# unset xtics
# unset ztics
# unset key

# define the data file and block separator
datafile = 'data/output_lorenz.txt'

# count the number of blocks in the file
stats datafile using 2:3 nooutput
nblocks = STATS_blocks

splot for[i=0:nblocks-1] datafile index i using 4:3:2 with lines linecolor i+1