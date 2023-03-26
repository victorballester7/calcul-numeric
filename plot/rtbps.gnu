set view equal xyz
set xlabel 'x' ; set ylabel 'y' ; set zlabel 'z'
set xrange [-1.2:-0.8]

splot 'data/output_rtbps.txt' u 2:3:4 w l,'data/L1.txt' w p,'data/trr.txt' w p,'data/lln.txt' w p