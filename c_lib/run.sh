#!\bin\bash
exec="./simu"
max_test=10000
output="/Dropbox/research/MATLAB/syn/est_channel/fisher_data/q16_8ary.csv"
$exec $max_test 1 > ~/$output
$exec $max_test 3 >> ~/$output
$exec $max_test 6 >> ~/$output
$exec $max_test 12 >> ~/$output
$exec $max_test 24 >> ~/$output
$exec $max_test 48 >> ~/$output
$exec $max_test 96 >> ~/$output

