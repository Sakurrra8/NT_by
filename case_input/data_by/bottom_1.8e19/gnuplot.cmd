set title noenhanced
set key noenhanced
set nolog x
set nolog y
set format y "%g"
set zero 1e-50
set timestamp
plot \
"gnuplot.data" using 1:2 title "" with linespoints
pause 3600
