#!\bin\bash
exec="./simu"
max_test=100000
output="/Dropbox/research/MATLAB/syn/est_channel/fisher_data/pitman_scale.csv"
n1=12
$exec $max_test 1 $n1 > ~/$output
for c in `seq 1 7`
do
	n1=`expr $n1 \* 2`
	$exec $max_test 1 $n1 >> ~/$output
	echo $n1
done

