gnuplot << EOF
set term png enhanced size 1920,1080 font "Times, 48"
set output "~/Dropbox/plot.png"
set logscale y 10
plot "mle" u 1:2 w lp, "" u 1:3 w lp
EOF
