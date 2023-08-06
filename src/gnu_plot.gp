print "Enter the max value of L"

max_value = int(system("read max_value_input && echo $max_value_input"))

xtics_values = ""
do for [i=0:10] {
    xtic_value = i * (max_value / 10)
    xtics_values = xtics_values.sprintf("\"%d\" %d, ", xtic_value, xtic_value)
}
xtics_values = substr(xtics_values, 1, strlen(xtics_values) - 2)
#set terminal eps
#set output "hyperraf.eps"
#set key below width -2 vertical maxrows 1
set multiplot layout 1,2 #title "RDF"

set style line 1 lt 1 lw 1.5 lc rgb "black" dashtype 2
set style line 2 lt 3 lw 3.5 lc  rgb "black"

xtics_values = substr(xtics_values, 1, strlen(xtics_values) - 2)

#set yrange[0:1]
set ylabel 'CPU Time (sec)'
set ytics nomirror
set xlabel "{/:Bold L}" font ",11" offset 0,0.5,0
set xtics nomirror
set xtics scale 0.0

set title "{/:Bold Mathematica}" font ",10"
plot "~/cpu_test_mtime.dat"  u 2 w lines ls 1 t "MHypergeometric", "~/cpu_use_mtime.dat"  u 2 w lines ls 2 t "MHyperRAF"
set title "{/:Bold Julia}" font ",10"
plot "~/cpu_test_jtime.dat"  u 2 w lines ls 1 t "JHypergeometric", "~/cpu_use_jtime.dat"  u 2 w lines ls 2 t "JHyperRAF"

unset multiplot
