use Data::Dumper;
while($line = <>){

	@array = split ",", $line;
	@sorted = sort {abs($a) <=> abs($b)} @array;
	print Dumper(\@sorted);
	#1 is because of the blank and -3 is to ignore the outliers
	print "Lowest and Highest: $sorted[1],$sorted[-3]\n";
	foreach $i (0..$#sorted){
		$temp = abs($sorted[$i]);
		$sum += $temp;
	}
	$total = $sum / $#array;
	print "Total: $total\n"
}
