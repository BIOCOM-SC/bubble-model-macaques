# Fig. 1C: one simulated animal (seed 16) at the fit-A optimum, day 50. Panel only (2 x 2 in); the
# three-panel Fig. 1 adds it next to panels A-B. Data: ../data/lesions_day50/seed16.txt, from
#   mpirun -n 1 bubble-model -noplot -dlobe -N 1 -seed 16 -snapshot 50   (params: fit-A optimum)
# PDF (cairo) instead of EPS: the lobes are drawn with transparency.
# Another snapshot: gnuplot -e "F='file'; OUTF='out.pdf'" lung_lesions.plt
load "common.gp"
set terminal pdfcairo enhanced font "Helvetica,12" size 2,2
if (!exists("OUTF")) OUTF = O."lung_lesions.pdf"
set output OUTF
if (!exists("F")) F = D."lesions_day50/seed16.txt"
LOBES = "../../lung/macaque/lobes/"
# colours of L1 ... L7 as in Fig. 1
array LCOL[7] = ["#90ff70", "#ff7000", "#ffee00", "#0fffee", "#0090ff", "#000090", "#000fff"]

# lesions -> arrays R, X, Y, Z
stats F using 1 nooutput
NLES = STATS_records
array R[NLES]; array X[NLES]; array Y[NLES]; array Z[NLES]
stats F using (i = int($0) + 1, R[i] = $1, X[i] = $2, Y[i] = $3, Z[i] = $4) nooutput

unset tics; unset border; unset key; unset colorbox
set xyplane relative 0
set xrange [-80:85]; set yrange [-105:65]; set zrange [-120:90]
set view 60, 30, 2.8
set pm3d depthorder
set parametric
set urange [-pi/2:pi/2]        # latitude
set vrange [0:2*pi]            # longitude
set samples 40; set isosamples 20

splot for [k=1:7] LOBES.sprintf("L%d.txt", k) skip 1 with polygons \
          fs transparent solid 0.15 noborder fc rgb LCOL[k], \
      for [i=1:NLES] X[i] + R[i]*cos(u)*cos(v), Y[i] + R[i]*cos(u)*sin(v), Z[i] + R[i]*sin(u) \
          with polygons fs solid noborder fc "black"
