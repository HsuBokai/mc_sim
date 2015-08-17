#!\bin\bash
exec="./simu"
max_test=100000
output="/Dropbox/research/MATLAB/syn/est_channel/fisher_data/q16_skew_rate.csv"
$exec $max_test 2 16 > ~/$output
#n1=24
k=3
for c in `seq 1 6`
do
#	n1=`expr $n1 \* 2`
	$exec $max_test $k 16 >> ~/$output
	k=`expr $k \* 2`
	echo $n1
done

