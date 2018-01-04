use Data::Dumper;

sub ReadSuperMatrix {
	
	
	$count = 0;
	$SuperName = ""; 
    local *SuperName = $_[0];
    %FastaHash = (); $name = "";
    $seq = "";
    #print "$SuperName\n";
    open(Supermatrix, "$SuperName");
	while($line = <Supermatrix>){
		
		chomp $line;
		if($line =~ /^>/){
			if($count != 0){
				
				$FastaHash{$name} = $seq;
			}
			$name = ($line =~ m/>(.*)/)[0];
			$seq = "";
		}else{
			$seq .= $line;
		}
		$count++;
	}
	$FastaHash{$name} = $seq;
	return %FastaHash;
	
	
}

#Get the likelihood for a given tree
sub TreeOne {
	
	@temp_array = ();
	@temp_array = @_;
	$likelihood_one = 0;
	$count = 0;
	@array = ();	
	#print "$temp_array[0]\n";
	open(file1, "$temp_array[0]");
	while($line = <file1>){
	
		if($count == $temp_array[1]){
		
			@array = split " ", $line;
			foreach $i (1..$#array){
				$likelihood_one += $array[$i];
			}
		}
		$count++;
	}
	$count = 0;
	return $likelihood_one;
}

#Get all the details from the Parts file
sub GetParts {
	
	@temp_array = ();
	@temp_array = @_;
	$PartFile = $temp_array[0];
	$IsItATest = $temp_array[1];
	$CountOfGenes = 0;
	open(Parts, "$PartFile")||die "No Parts file\nperl SSLL_INDV_TREE.pl help";
	while($line = <Parts>){

		chomp $line;
		$location = ($line =~ m/.*? = (.*)/)[0];
		#This can be used to specify only a few genes
		#CHANGE AFTER TESTING!!!!!!!!!!!!!!!!
		if($IsItATest eq "True"){
			if($CountOfGenes < 2){
				
				push @loc, $location;
		
			}else{}
		}else{
			push @loc, $location;
		}
		$CountOfGenes++;		
	}
	return @loc;
}

#Function to Condense the Sites to Genes
sub CondenseToGene {
	@gene_array = @_;
	$genes = "$gene_array[2]";
	$concatenated = "";
	$count = 0;
	open(Model,"$gene_array[0]")||die "No model\n";
	while($line = <Model>){

        chomp $line;
        $len = ($line =~ m/.*?=(.*)/)[0];
        ($start,$stop) = split "-", $len, 2;
        $start -= 1;
        $stop -= 1;
		push @stop_array, $stop;
		push @start_array, $start;
	}
	open(sites, "$gene_array[1]")||die "No sites\n";
	while($line = <sites>){
		
		chomp $line;
		@array = split " ", $line;

		$name = $array[0];
		splice @array, 0, 1;
		#print Dumper(\@array);
		if($count != 0){
			#print "$name\t";
			$concatenated .= "$name\t";
			foreach $i (0..$#stop_array){
				$start = $start_array[$i];
				$stop = $stop_array[$i];
				foreach $j ($start..$stop){
					
						$value += $array[$j];
						if($array[$j] eq ""){
							#print "ERROR Condensing to genes Kill, Kill, Kill!!!!!!!!!!!\n";
							#print StatsOut "ERROR Condensing to Gene File is not being read in correctly!!!!!!!!!!!\n";
						}
						#print "At Site $j: $array[$j] $value\n";
						
				}
			#print "$value ";
			$concatenated .= "$value ";
			$value = 0;
			}
			#print "\n";
			$concatenated .= "\n";
		}else{
			
				#print "$line\n";
				($trees,$sit) = split " ", $line;
				$concatenated .= "  $trees  $gene_array[2]\n";
		}
		$count++;
	}
	return $concatenated;
}


sub GetBipartitions {
	
	@array_of_conflict = (); @sorted_conflict = ();
	@array_of_clade = (); @sorted_clade = ();
	%Conf_node = (); $count_con = 1;
	
	open(log_file, "bp.log") || die "No log File";
	while($line = <log_file>){
		
		chomp $line;
		@array_of_conflict = split ",", $conflicting_node;
		@sorted_conflict = sort @array_of_conflict;
		if ($line =~ /CLADE/){
			
			$clade = ($line =~ m/CLADE: (.*?) \| .*/)[0];
			@array_of_clade = split " ", $clade;
			@sorted_clade = sort @array_of_clade;
		}
		if(@sorted_conflict ~~ @sorted_clade){
			#This will contain the clade of interest
			if($line =~ /CLADE/){
				
				print StatsOut "The concordant relationship is a.k.a 0: $conflicting_node\n"; 
				$Conf_node{"0"} = $conflicting_node;
				#push @Questionable, $conflicting_node;
			}
			if($line =~ /COUNT/ && $line !~ /ICA/){
				
				$conflict = ($line =~ m/\t (.*?) \| .*/)[0];
				$conflict =~ s/ /,/g;
				print StatsOut "The conflicting relationships are number $count_con: $conflict\n";
				$Conf_node{"$count_con"} = "$conflict";
				#push @Questionable, $conflict;
				$count_con++;
			}			
		}	
	}
	return %Conf_node;
}

#Get the bipartitions that each tree is associated with
sub TreeBiparts {
	
	%hash = @_;
	for $keys (sort {lc $a cmp lc $b} keys %hash){
		push @Questionable, $hash{$keys};
	}
	$tree_count = 0;
	open(trees, "Unique.tre")||die "No tree file\nperl SSLL_INDV_TREE.pl help";
	while($line = <trees>){
		
		chomp $line;
		open(tempout, ">temptesttre");
		print tempout "$line\n";
		system("$pxbp -t temptesttre > temp.log");
		open(temp, "temp.log");
		while($temp = <temp>){
		
			chomp $line;
			$clade = ($temp =~ m/CLADE: (.*?) \t.*/)[0];
			@temp = split " ", $clade;
			@sort_temp = sort @temp;
			foreach $i (0..$#Questionable){
				
				@other_temp = split ",", $Questionable[$i];
				@sort_other_temp = sort @other_temp;
				if(@sort_other_temp ~~ @sort_temp){
					
					print StatsOut "Tree $tree_count: Matches Conflict $i\n";
					#TREE HASH will have the tree along with the conflicts it supports
					$TREE_HASH{$tree_count} .= "$i,";
					
				}
			}
		}
		$tree_count++;
		
	}
	return %TREE_HASH;
	system("rm temp.log temptesttre");
}
#Get Stats for the Supermatrix Gene Comparison
sub SuperMatrixGeneStats {
	
	@array = @_;
	@Sort = ();
	@Trees = ();
	@all_array = ();
	@GeneInfoArray = ();
	@lines = split "\n", $array[0];
	$count = 0;
	foreach $i (0..$#lines){
		$line = $lines[$i];
		if($count != 0){
			@array = split " ", $line;
			push @all_array, [@array];
		}else{
			
			($trees,$sites) = split " ", $line;
		}
		$count++;
	}
	$count -= 2;
	#sites is the number of genes
	foreach $i (1..$sites){
	
		foreach $j (0..$count){
			
			#print "$all_array[$j][$i]\t";
			@numbs = split ",", $TREE_HASH{$j};
			#$BEST_HASH now has each genes value bipartition => Likelihood:Tree#, Likelihood:Tree#
			foreach $k (0..$#numbs){
					
				$BEST_HASH{$numbs[$k]} .= "$all_array[$j][$i]:$j,";
			}
		}	
		#print Dumper(\%BEST_HASH);
		#Get the best tree for each bipartition
		
		for $keys (sort {lc $a cmp lc $b} keys %BEST_HASH){
		
			@Elegant = split ",", $BEST_HASH{$keys};
			#Get the best likelihood for the bipartition
			#Of note this only considers on tree to be able
			#To contribute top to the bipartition, maybe good
			#To look into for lowering parameters
			@Sort = sort {$b <=> $a} @Elegant;
			#print Dumper(\@Sort);
			#@Sort = sort @Elegant; !!!!!!!!!!
			push @Trees, "$Sort[0]:$keys";

		}
		
		#Sorted by best Likelihood and in the format Likelihood:TreeItComesFrom:BipartitionItSupports
		@SortTrees = sort {$b <=> $a} @Trees;
		#@SortTrees = sort @Trees; !!!!!!!!!!!!!!!!
		#Get the top likelihood bipartition while keeping in mind multiple biparitions may have it
		push @GeneInfoArray, [@SortTrees];
		#print Dumper(\@SortTrees);
		#print Dumper(\%BEST_HASH);
		%HELD_HASH = ();
		%BEST_HASH = ();
		@Trees = ();
		
	}
	return @GeneInfoArray;
}
#Get the frequency
sub FrequencyCalcs {
	
	%FrequencyHash = (); 
	local *FrequencyHash = $_[0];
    %ConflictRelationship = (); 
    local *ConflictRelationship = $_[1];
	for $keys (sort {lc $a cmp lc $b} keys %FrequencyHash){
		
		print StatsOut "Relationship $keys\t$ConflictRelationship{$keys}:\t$FrequencyHash{$keys}\n";
	}
    #print Dumper(\%FrequencyHash);
}
#Get the information for each Likelihood calculation
sub LikelihoodCalcs{
	
	$Tree1Like = 0; 
	local *Tree1Like = $_[0];
    %ConflictRelationship = (); 
    local *ConflictRelationship = $_[1];
    %TreeTotalLikes = (); 
    local *TreeTotalLikes = $_[2];
    $Tree2Like = 0;
    local *Tree2Like = $_[3];
	print StatsOut "The Likelihood of Your first tree: $Tree1Like\n";
	print StatsOut "The Likelihood of Your second tree: $Tree2Like\n";
	@array = ();
	@sorted_array = ();
	push @array, $Tree1Like;
	push @array, $Tree2Like;
	for $keys (sort {lc $a cmp lc $b} keys %TreeTotalLikes){
		print StatsOut "Relationship $keys\t$ConflictRelationship{$keys}:\t$TreeTotalLikes{$keys}\n";
		push @array, $TreeTotalLikes{$keys};
	}
	@sorted_array = sort {$b <=> $a} @array;
	print StatsOut "The Best Likelihood is: $sorted_array[0]\n";
}

sub MGWB_BL_LikelihoodCalcs{
	    
	%TreeTotalLikes = (); local *TreeTotalLikes = $_[0];
	%TreeRels = (); local *TreeRels = $_[1];
	@array = (); @sorted_array = ();
	for $keys (sort {lc $a cmp lc $b} keys %TreeTotalLikes){
		print StatsOut "Relationship $keys\t$TreeRels{$keys}:\t$TreeTotalLikes{$keys}\n";
		push @array, $TreeTotalLikes{$keys};
	}
	@sorted_array = sort {$b <=> $a} @array;
	print StatsOut "The Best Likelihood is: $sorted_array[0]\n";
	
}

#Get the total parameters for each analysis Number of trees used * (2*Taxa - 3) + (6*Number of Genes)
sub ParameterStatsForMGWB{
	

    %ConflictRelationship = (); 
    local *ConflictRelationship = $_[0];
	%ParametersUsed = (); 
	local *ParametersUsed = $_[1];
	$TotalBranches = 0;
	local *TotalBranches = $_[2];
	$TotalGenes = 0;
	local *TotalGenes = $_[3];
	$params_of_tree = (2*$TotalBranches - 3) + ($TotalGenes*6);
	print StatsOut "The Parameters in your first tree: $params_of_tree\n";
	print StatsOut "The Parameters in your second tree: $params_of_tree\n";
	for $keys (sort {lc $a cmp lc $b} keys %ConflictRelationship){
		#$params_of_tree = ($ParametersUsed{$keys} * (2*$TotalBranches - 3)) + ($TotalGenes*6);
		$params_of_tree = ($ParametersUsed{$keys} * (2*$TotalBranches - 3)) + ($TotalGenes*6);
		print StatsOut "Relationship $keys\t$ConflictRelationship{$keys}:\t$params_of_tree\n";
	}
	
}

sub MSWBparameters{


    %ConflictRelationship = (); 
    local *ConflictRelationship = $_[0];
	%ParametersUsed = (); 
	local *ParametersUsed = $_[1];
	$TotalBranches = 0;
	local *TotalBranches = $_[2];
	$TotalGenes = 0;
	local *TotalGenes = $_[3];
	$loud = 0;
	local *loud = $_[4];
	%HashofParameters = ();
	$params_of_tree = 0; $params_of_bipart = 0;
	$params_of_tree = (2*$TotalBranches - 3) + ($TotalGenes*6);
	if($loud eq "True"){
		print StatsOut "The Parameters in your first tree: $params_of_tree\n";
		print StatsOut "The Parameters in your second tree: $params_of_tree\n";
	}
	for $keys (sort {lc $a cmp lc $b} keys %ConflictRelationship){
		$params_of_bipart = (($ParametersUsed{$keys} * 6) + (((2*$TotalBranches)-3) * $ParametersUsed{$keys}));
		if($loud eq "True"){
			print StatsOut "Relationship $keys\t$ConflictRelationship{$keys}:\t$params_of_bipart\n";
		}
		$HashofParameters{$keys} = $params_of_bipart;
	}
	return %HashofParameters;
}

sub MSWB_BLparameters{
	
    %ConflictRelationship = (); 
    local *ConflictRelationship = $_[0];
	%ParametersUsed = (); 
	local *ParametersUsed = $_[1];
	$TotalBranches = 0;
	local *TotalBranches = $_[2];
	$TotalGenes = 0;
	local *TotalGenes = $_[3];	
	$loud = 0;
	local *loud = $_[4];
	for $keys (sort {lc $a cmp lc $b} keys %ConflictRelationship){
		$params_of_bipart = (($ParametersUsed{$keys} * 6) + (((2*$TotalBranches)-3) * $ParametersUsed{$keys}));
		if($loud eq "True"){
			print StatsOut "Relationship $keys\t$ConflictRelationship{$keys}:\t$params_of_bipart\n";
		}
		$HashofParameters{$keys} = $params_of_bipart;	
	}
	return %HashofParameters;
}

#Get the AIC values
sub AICcalcs{
	
	%ConflictRelationship = (); 
    local *ConflictRelationship = $_[0];
	%ParametersUsed = (); 
	local *ParametersUsed = $_[1];
	$TotalBranches = 0;
	local *TotalBranches = $_[2];
	$TotalGenes = 0;
	local *TotalGenes = $_[3];
	$Tree1Like = 0; 
	local *Tree1Like = $_[4];
	%TreeTotalLikes = (); 
    local *TreeTotalLikes = $_[5];
    $Tree2Like = 0;
    local *Tree2Like = $_[6];
	$params_of_tree = (2*$TotalBranches - 3) + ($TotalGenes*6);
	$aic = (-2*$Tree1Like) + (2*$params_of_tree);
	@array = ();
	@sorted_array = ();
	push @array, $aic;
	print StatsOut "The AIC of your first tree: $aic\n";
	$aic = (-2*$Tree2Like) + (2*$params_of_tree);
	print StatsOut "The AIC of your second tree: $aic\n";
	push @array, $aic;
	for $keys (sort {lc $a cmp lc $b} keys %ConflictRelationship){
		$params_of_tree = ($ParametersUsed{$keys} * (2*$TotalBranches - 3)) + ($TotalGenes*6);
		#print "$params_of_tree\n";
		$aic = (2*$params_of_tree) - (2*$TreeTotalLikes{$keys});
		print StatsOut "Relationship $keys\t$ConflictRelationship{$keys}:\t$aic\n";
		push @array, $aic;
	}
	@sorted_array = sort {$a <=> $b} @array;
	print StatsOut "The Best AIC is: $sorted_array[0]\n";
	return $sorted_array[0];
}
#Get the difference from best AIC
sub DeltaAIC {
	
	%ConflictRelationship = (); 
    local *ConflictRelationship = $_[0];
	%ParametersUsed = (); 
	local *ParametersUsed = $_[1];
	$TotalBranches = 0;
	local *TotalBranches = $_[2];
	$TotalGenes = 0;
	local *TotalGenes = $_[3];
	$Tree1Like = 0; 
	local *Tree1Like = $_[4];
	%TreeTotalLikes = (); 
    local *TreeTotalLikes = $_[5];
    $top_AIC = 0;
    local *top_AIC = $_[6];
    $Tree2Like = 0;
    local *Tree2Like = $_[7];
    $TotalDelta = 0;
    $aic = (-2*$Tree1Like) + (2*$params_of_tree);
    $change = $aic - $top_AIC;
    $TotalDelta += exp(-0.5 * $change);
    print StatsOut "The Delta AIC for Tree 1 is: $change\n";
    $aic = (-2*$Tree2Like) + (2*$params_of_tree);
    $change = $aic - $top_AIC;
    $TotalDelta += exp(-0.5 * $change);
    print StatsOut "The Delta AIC for Tree 2 is: $change\n";
	for $keys (sort {lc $a cmp lc $b} keys %ConflictRelationship){
		$params_of_tree = ($ParametersUsed{$keys} * (2*$TotalBranches - 3)) + ($TotalGenes*6);
		$aic = (2*$params_of_tree) - (2*$TreeTotalLikes{$keys});
		$change = $aic - $top_AIC;
		print StatsOut "Relationship $keys\t$ConflictRelationship{$keys}:\t$change\n";
		$TotalDelta += exp(-0.5 * $change);	
		
	}
	
	$aic = (-2*$Tree1Like) + (2*$params_of_tree);
	$change = $aic - $top_AIC;
	$weight = (exp(-0.5 * $change) / $TotalDelta);
	print StatsOut "######################## AIC Weight####################\n";
	print StatsOut "The AIC weight for Tree 1 is: $weight\n";
	$aic = (-2*$Tree2Like) + (2*$params_of_tree);
	$change = $aic - $top_AIC;
	$weight = (exp(-0.5 * $change) / $TotalDelta);
	print StatsOut "The AIC weight for Tree 2 is: $weight\n";
	for $keys (sort {lc $a cmp lc $b} keys %ConflictRelationship){
		$params_of_tree = ($ParametersUsed{$keys} * (2*$TotalBranches - 3)) + ($TotalGenes*6);
		$aic = (2*$params_of_tree) - (2*$TreeTotalLikes{$keys});
		$change = $aic - $top_AIC;
		$weight = (exp(-0.5 * $change) / $TotalDelta);
		print StatsOut "Relationship $keys\t$ConflictRelationship{$keys}:\t$weight\n";

	}
}
#Get bipartition support by site
sub SiteBipartSupport{
	
	$gene_len = 0; local *gene_len = $_[0];
	%TreesAndBiparts = (); local *TreesAndBiparts = $_[1];
	$super = 0; local *super = $_[2];
	$count = 0; @AllArray = ();
	$q = 0; %BipartLikes = ();
	%SiteBipart = ();
	($start, $stop) = split "-", $gene_len;
	#print "$start\n";
	#print Dumper(\%TreesAndBiparts);
	%BipartLikes = ();
	print "MSWB Calc for gene of length: $gene_len\n";
	open(SuperMat, "$super");
	while($line = <SuperMat>){
		
		if($count != 0){
			$q = $count - 1;
			@array = split " ", $line;
			foreach $i ($start..$stop){
				#print "Tree $count: $array[$i]\t";
				#print "Bipartition: $TreesAndBiparts{$count}\n";
				@AllArray = split ",", $TreesAndBiparts{$q};
				
				#Get Likelihood associated with each bipartition
				foreach $j (0..$#AllArray){
					#print "Site $i Tree $q Has a Likelihood: $array[$i] Corresponding to Bipartition: $AllArray[$j]\n";
					#Order is Likelihood Tree and then bipartition
					$BipartLikes{$i} .= "$array[$i]:$q:$AllArray[$j],";
				}
			}
			#print Dumper(\%BipartLikes);
			#print "\n";
		}
		$count++;
	}
	#print Dumper(\%BipartLikes);
	return %BipartLikes;
}
#get the AIC score of the MSWBs
sub MSWBAic{
	
	%AllParameters = ();
	local *AllParameters = $_[0];
	%ConflictRelationship = (); 
    local *ConflictRelationship = $_[1];
	$Tree1Like = 0; 
	local *Tree1Like = $_[2];
	%TreeTotalLikes = (); 
    local *TreeTotalLikes = $_[3];
    $Tree2Like = 0;
    local *Tree2Like = $_[4];
    $taxa_used = 0;
    local *taxa_used = $_[5];
    $aic_of_tree1 = 0; $aic_of_tree2 = 0;
    $aic_of_bipart = 0; $param = 0;
    @all_vals = (); @only_mswb = ();
    @sorted_only_mswb = (); @sorted_array = ();
    
    $param = (((2 * $taxa_used) - 3) + 6);
    $aic_of_tree1 = (-2 * $Tree1Like) + (2 * $param);
    $aic_of_tree2 = (-2 * $Tree2Like) + (2 * $param);
    print StatsOut "The AIC score of the first tree: $aic_of_tree1\n";
	print StatsOut "The AIC of the second tree: $aic_of_tree2\n";
	push @all_vals, $aic_of_tree1;
	push @all_vals, $aic_of_tree2;

    foreach my $keys (sort {lc $a cmp lc $b} keys %ConflictRelationship){
		
		$aic_of_bipart = (-2 * $TreeTotalLikes{$keys}) + (2 * $AllParameters{$keys});
		print StatsOut "Relationship $keys\t$ConflictRelationship{$keys}:\t$aic_of_bipart\n";
		push @all_vals, $aic_of_bipart;
		push @only_mswb, $aic_of_bipart;
	}
	@sorted_array = sort {$a <=> $b} @all_vals;
	@sorted_only_mswb = sort  {$a <=> $b} @only_mswb;
	print StatsOut "The best AIC is: $sorted_array[0]\n";
	return ($sorted_array[0], $sorted_only_mswb[0]);
	
}
#MSWBdeltaAIC(\%TotalParams, \%ConflictHash, \%GeneLikelihood, \$TreeOneLikelihood, \$TreeTwoLikelihood, \$best_aic_of_all, \$count_of_seqs);
sub MSWBdeltaAIC {
	
	
	%AllParameters = (); local *AllParameters = $_[0];
	%ConflictRelationship = (); local *ConflictRelationship = $_[1];
    %TreeTotalLikes = (); local *TreeTotalLikes = $_[2];
	$Tree1Like = 0; local *Tree1Like = $_[3];
    $Tree2Like = 0; local *Tree2Like = $_[4];
    $best_aic = 0; local *best_aic = $_[5];
    $taxa_used = 0; local *taxa_used = $_[6];
	$aic_of_tree1 = 0; $aic_of_tree2 = 0;
    $aic_of_bipart = 0; $param = 0;
    $delta_aic = 0; $TotalDelta = 0;
    @all_vals = (); @only_mswb = ();
    @sorted_only_mswb = (); @sorted_array = ();
    @VerboseAIC = (); @SortedVerbose = ();
    $param = (((2 * $taxa_used) - 3) + 6);
    
    $aic_of_tree1 = (-2 * $Tree1Like) + (2 * $param);
    push @VerboseAIC, "$aic_of_tree1:-2";
    $delta_aic = $aic_of_tree1 - $best_aic;
    $TotalDelta += exp(-0.5 * $delta_aic);
    print StatsOut "Delta AIC of Tree 1: $delta_aic\n";
    $aic_of_tree2 = (-2 * $Tree2Like) + (2 * $param);
    push @VerboseAIC, "$aic_of_tree2:-1";
    $delta_aic = $aic_of_tree2 - $best_aic;
    $TotalDelta += exp(-0.5 * $delta_aic);
	print StatsOut "Delta AIC of Tree 2: $delta_aic\n";
    $aic_of_tree2 = (-2 * $Tree2Like) + (2 * $taxa_used);
	foreach my $keys (sort {lc $a cmp lc $b} keys %ConflictRelationship){
		
		$aic_of_bipart = (-2 * $TreeTotalLikes{$keys}) + (2 * $AllParameters{$keys});
		push @VerboseAIC, "$aic_of_bipart:$keys";
		$delta_aic = $aic_of_bipart - $best_aic;
		$TotalDelta += exp(-0.5 * $delta_aic);
		print StatsOut "Relationship $keys\t$ConflictRelationship{$keys}:\t$delta_aic\n";
	}

	$aic_of_tree1 = 0; $aic_of_tree2 = 0;
    $aic_of_bipart = 0; $param = 0;
    $delta_aic = 0; $weight = 0;
    $param = (((2 * $taxa_used) - 3) + 6);
	print StatsOut "#########################AIC Weight#####################\n";
    
    $aic_of_tree1 = (-2 * $Tree1Like) + (2 * $param);
    $delta_aic = $aic_of_tree1 - $best_aic;
	$weight = (exp(-0.5 * $delta_aic) / $TotalDelta);
	print StatsOut "The Weight of Tree 1: $weight\n";
	
	$aic_of_tree2 = (-2 * $Tree2Like) + (2 * $param);
	$delta_aic = $aic_of_tree2 - $best_aic;
	$weight = (exp(-0.5 * $delta_aic) / $TotalDelta);
	print StatsOut "The Weight of Tree 2: $weight\n";
	
	foreach my $keys (sort {lc $a cmp lc $b} keys %ConflictRelationship){
		
		$aic_of_bipart = (-2 * $TreeTotalLikes{$keys}) + (2 * $AllParameters{$keys});
		$delta_aic = $aic_of_bipart - $best_aic;
		$weight = (exp(-0.5 * $delta_aic) / $TotalDelta);
		print StatsOut "Relationship $keys\t$ConflictRelationship{$keys}:\t$weight\n";
	}
	@SortedVerbose = sort {$a <=> $b} @VerboseAIC;
	return @SortedVerbose;
}




open(Configure, "$ARGV[0]")||die "Please See Configure File, Goodbye\n";
while($line = <Configure>){
	
	if($line =~ /^pxrmt:/){
		$pxrmt = ($line =~ /.*?: (.*)/)[0];	
	}elsif($line =~ /^pxbp:/){
		$pxbp = ($line =~ /.*?: (.*)/)[0];	
	}elsif($line =~ /^raxml:/){
		$raxml = ($line =~ /.*?: (.*)/)[0];
	}elsif($line =~ /^Species:/){
		$conflicting_node = ($line =~ /.*?: (.*)/)[0];
	}elsif($line =~ /^outfile:/){
		$outfile = ($line =~ /.*?: (.*)/)[0];
	}elsif($line =~ /^Supermatrix:/){
		$SuperMatrix = ($line =~ /.*?: (.*)/)[0]; 
	}elsif($line =~ /^PartsFile:/){
		$PartFile = ($line =~ /.*?: (.*)/)[0]; 
	}elsif($line =~ /^Set:/){
		$TreeFile = ($line =~ /.*?: (.*)/)[0]; 
	}elsif($line =~ /^Threads:/){
		$threads = ($line =~ /.*?: (.*)/)[0]; 
	}elsif($line =~ /^Test:/){
		$IsItATest = ($line =~ /.*?: (.*)/)[0];
	}elsif($line =~ /^Verbose:/){
		$verbose  = ($line =~ /.*?: (.*)/)[0];	
	}elsif($line =~ /^RUN_MGWB:/){
		$MGWB = ($line =~ /.*?: (.*)/)[0];	
	}elsif($line =~ /^RUN_MSWB:/){
		$MSWB = ($line =~ /.*?: (.*)/)[0];	
	}elsif($line =~ /^MGWB_BL:/){
		$MGWB_BL = ($line =~ /.*?: (.*)/)[0];	
	}elsif($line =~ /^MSWB_BL:/){
		$MSWB_BL = ($line =~ /.*?: (.*)/)[0];
	}elsif($line =~ /^secret:/){
		$secret = ($line =~ /.*?: (.*)/)[0];
	}
}
open(StatsOut, ">$outfile")||die "In program give a name of the outfile\n";

######################################################################
#Designed to factor branch lengths into the SSLL calculation
#If taxa are missing genes this will remove them
#This should all be in a config file in the future
#Really rough early code!!!!!!!!!!!

#Location of pxrmt
#$pxrmt = "pxrmt";
#Location of RAxML
#$raxml = "raxmlHPC";
#Name of outfile
#$outfile = "Stats";
#Location of pxbp
#$pxbp = "pxbp";
#Species in your clade of interest, so in the first tree of your tree set
#What are the species that are in that relationship? Put them in the quotes
#As comma separated
#$conflicting_node = "alligator,caiman,phrynops,caretta,chelonoidis_nigra,emys_orbicularis";
#$conflicting_node = "Gallus,Taeniopygia,alligator,caiman,phrynops,caretta,chelonoidis_nigra,emys_orbicularis";
#$conflicting_node = "Homo,Monodelphis,Ornithorhynchus";
#open(StatsOut, ">$outfile")||die "In program give a name of the outfile\n";

#This is the order run it by doing perl EdgeTester.pl ExampleConcat.fa ExampleGenes.model ExampleAllGeneTrees.tre 
#$SuperMatrix = $ARGV[0];
#$PartFile = $ARGV[1];
#$TreeFile = $ARGV[2];
#Set the number of threads
#$threads = 3;
#Set this to True for a quick 2 gene analysis, anything else for a full analysis
#$IsItATest = "False";
#$IsItATest = "True";
#This prints out a ton of the stats, good complement to a test run
#$verbose = "True";
#Put to true for the MGWB calcs
#$MGWB = "Tre";
#Put to true to run the MSWB calcs
#$MSWB = "Tre";
#Put True for the MGWB with branch lengths estimated separately
#$MGWB_BL = "Tre";
#Put True for the MSWB with branch lengths estimated separately
#$MSWB_BL = "True";
#######################################################################
print StatsOut "################################################################\n";
print StatsOut "If any of this is wrong the analysis won't work, double check!!!!\n";
print StatsOut "You have set pxrmt to be in the path: $pxrmt\n";
print StatsOut "You have set pxbp to be in the path: $pxbp\n";
print StatsOut "You have set raxml to be in the path: $raxml\n";
print StatsOut "The species in the relationship you are looking at are: $conflicting_node\n";
print StatsOut "Your supermatrix is file called: $SuperMatrix\n";
print StatsOut "Your partition file is called: $PartFile\n";
print StatsOut "Your tree set file is (newick hopefully?): $TreeFile\n";
print StatsOut "Running a test?: $IsItATest\n";
print StatsOut "Running it in verbose?: $verbose\n";
print StatsOut "Conducting MGWB?: $MGWB\n";
print StatsOut "Conducting MSWB?: $MSWB\n";
print StatsOut "Conducting MGWB_BL?: $MGWB_BL\n";
print StatsOut "Conducting MSWB_BL?: $MSWB_BL\n";
print StatsOut "################################################################\n";


#Creat a file called Unique.tre with the Unique trees
#and a file called bp.log
system("pxbp -t $TreeFile -u | grep \"\(\" > Unique.tre");
system("pxbp -t Unique.tre -e -v > bp.log");
#CountOfGenes is the number of genes
#location is the length of the gene
#loc is an array with the length of genes
@loc = GetParts($PartFile, $IsItATest);
$CountOfGenes = ($#loc+1);

#Get a Hash with Conflict
(%ConflictHash) = GetBipartitions();

#Hash with trees and there bipartitions
%TreeHash = TreeBiparts(%ConflictHash);

#Shared by MGWB and MSWB
if($MGWB eq "True" or $MSWB eq "True"){
	
	print "Running RAxML Supermatrix Analysis\n";
	print "This can take a while...\n";
	if($secret eq "True"){
	}else{
		system("$raxml -f g -T $threads -s $SuperMatrix -m GTRGAMMA -z Unique.tre -q $PartFile  -n SUPER_EX_SSLL | grep \"Tree\" | grep \":\"");
	}	
	#system("$raxml -f g -T 3 -s $SuperMatrix -m GTRGAMMA -z Unique.tre -q $PartFile  -n EX_SSLL");
	system("cp RAxML_perSiteLLs.SUPER_EX_SSLL SuperMatrixSiteLLs.SSLL");
	$count_of_seqs = 0;
	open(file, "$SuperMatrix")||die "No Supermatrix";
	while($line = <file>){
		if($line =~ /^>/){
			$count_of_seqs++;
		}
	}
	
}
if($MGWB_BL eq "True" or $MSWB_BL eq "True"){
	
	print StatsOut "################################################################\n";
	print StatsOut "#YOU ARE RUNNING WITH BRANCH LENGTHS VARYING\n";
	print StatsOut "#HERE ARE THE TAXA AND MISSING ONES FROM EACH GENE:\n";
	print "Estimating Gene Tree Likelihoods\n";
	%FastaHash = ();
	$count = 0; %TotalTaxaHash = ();
	%FastaHash = ReadSuperMatrix(\$SuperMatrix);
	open(MGWB_BL_OUT, ">SuperMatrixSiteLLsBL.SSLL");
	foreach $i (0..$#loc){
		#This block creates a temporary Fasta for the gene
		open(out, ">temp.fa");
		#print "$loc[$i]\n";
		($start,$stop) = split "-", $loc[$i], 2;
		$dif = $stop - $start + 1;
		$to_remove = "";
		$seq_count = 0;
		for $keys (sort keys %FastaHash){
		
			$seq = substr $FastaHash{$keys}, ($start-1), $dif;
			if($AMINO eq "TRUE"){
				$missing = $seq =~ tr/\-X/\-X/;
			}else{
				$missing = $seq =~ tr/N\-X/N\-X/;
			}
			#print "Diff = $dif\tMissing = $missing\n";
			if($missing != $dif){
				print out ">$keys\n$seq\n";
				#All the seqs going into the likelihood calc, helps with
				#Number of parameters
				$seq_count++;
			}else{
				#print out ">$keys\n$seq\n"; #Temp
				$to_remove .= "$keys,"
			}
		}
		$TotalTaxaHash{$i} = $seq_count;
		print StatsOut "Gene_$i: Has $seq_count in it\n";
		#Anything that was missing from the gene gets eliminated here
		print StatsOut "These were removed from likelihood calc due to all missing data: $to_remove\n";
		if($to_remove ne ""){
			system("$pxrmt -t Unique.tre -n $to_remove > TempTree.tre");
		}else{
			system("cp Unique.tre TempTree.tre");
		}
		#print "$raxml\n";
		#Deletes RAxML Files
		system("rm RAxML_*");
		print "(☞ﾟヮﾟ)☞ Processing Gene $i\n";
		system("$raxml -f g -T $threads -s temp.fa -m GTRGAMMA -z TempTree.tre -n EX_SSLL | grep \"Tree\" | grep \":\"");
		$sitecount = 0;
		open(RAXML, "RAxML_perSiteLLs.EX_SSLL");
		while($line = <RAXML>){
			
			chomp $line;
			if($sitecount != 0){
				($name, $sites) = split " ", $line, 2;
				#if(exists $HASH{$name}){
					
				#	$HASH{$name} = "$sites";
				#}else{
					#print "Here $name: $sites\n";
					$HASH{$name} .= "$sites";
				#}
			}
			$sitecount++;
			
		}
		
	}
	#print Dumper(\%HASH);
	$test = "";
	$hash_size = keys %HASH;
	print MGWB_BL_OUT "  $hash_size  $stop\n";
	foreach $keys (sort {lc $a cmp lc $b} keys %HASH){
		print MGWB_BL_OUT "$keys $HASH{$keys}\n";
		
	}
	
	
}

############
#MGWB Method
############
if($MGWB eq "True"){
	
	print StatsOut "##########################################\n";
	print StatsOut "MGWB Analysis (Calculated Across Matrix):\n";
	
	#Condense the data to genes (I think this happens twice? but code works so play with later?)
	$CondensedGenes = CondenseToGene($PartFile, "SuperMatrixSiteLLs.SSLL",$CountOfGenes);
	#Get the stats for the first two trees, make sure its not getting paramaterized to death
	$TreeOneLikelihood = TreeOne("SuperMatrixSiteLLs.SSLL", 1);
	$TreeTwoLikelihood = TreeOne("SuperMatrixSiteLLs.SSLL", 2);
	#print "$TreeTwoLikelihood\n";
	
	#Condensed Genes Contains the condensed down to genes
	#$CondensedGenes = CondenseToGene($PartFile, "SuperMatrixSiteLLs.SSLL",$CountOfGenes);
	
	#GeneArray will contain a sorted array of best likelihood with the trees
	@GeneArray = SuperMatrixGeneStats($CondensedGenes);
	
	#Summarize the results
	$top_likelihood = 0; %TopTree = ();
	%Frequency = (); %TreeLikelihoods = ();
	$TotalDelta = 0; @array = ();
	
	#Loop all the genes in the analysis
	foreach $i (0..$#GeneArray){
		$gene = $i + 1;
		@array = @{$GeneArray[$i]};
		#Start looping each gene
		foreach $j (0..$#array){
			@StatsArray = split ":", $array[$j];
			$TreeLikelihoods{"$StatsArray[2]"} += $StatsArray[0];
			
			if($verbose eq "True"){
				print StatsOut "Gene: $gene\tLikelihood:\t$StatsArray[0]\tFrom Tree: $StatsArray[1]\tSupports Bipartion: $StatsArray[2]\n";
			}
			if($j == 0){
			
				$top_likelihood = $StatsArray[0];
				$top_stats .= "$gene\t$StatsArray[0]\t$StatsArray[2]\t$ConflictHash{$StatsArray[2]}\t$StatsArray[1]\n";
				$Frequency{"$StatsArray[2]"}++;
			}else{
				#Get counts if multiple hits for the best likelihood
				if($top_likelihood == $StatsArray[0]){
					$top_stats .= "$gene\t$StatsArray[0]\t$StatsArray[2]\t$ConflictHash{$StatsArray[2]}\t$StatsArray[1]\n";
					$Frequency{"$StatsArray[2]"}++;

				}
			}
			$TopTree{"$StatsArray[2]"} .= "$StatsArray[1],";
		

		}
	
	}
	
	%ParamHash = (); @array = ();
	%ParamsPerBipart = ();
	
	#Creates something called ParamsPerBipart with total different trees in it
	foreach $keys (sort {lc $a cmp lc $b} keys %TopTree){
		@array = split ",", $TopTree{$keys};
		foreach $i (0..$#array){
			if(exists $ParamHash{$array[$i]}){
				
			}else{
				$ParamHash{$array[$i]}++;
			}
		}
		$count = keys %ParamHash;
		$ParamsPerBipart{"$keys"} = $count;
		%ParamHash = ();
	}

	print StatsOut "#########################MGWB INDV GENE STATS############################\n";
	print StatsOut "Gene\tLikelihood\tBipartition\tRelationship\tTree\n";
	print StatsOut $top_stats;
	
	#Summarize everything
	
	#Frequency of bipartition
	print StatsOut "#####################MGWB Frequency Stats##############\n";
	FrequencyCalcs(\%Frequency, \%ConflictHash);
	
	#Likelihood
	print StatsOut "#####################MGWB Likelihood Stats##############\n";
	LikelihoodCalcs(\$TreeOneLikelihood, \%ConflictHash, \%TreeLikelihoods, \$TreeTwoLikelihood);
	
	#Parameters
	print StatsOut "#####################MGWB Parameters Used###############\n";
	ParameterStatsForMGWB(\%ConflictHash, \%ParamsPerBipart, \$count_of_seqs, \$CountOfGenes);
	
	#AIC
	print StatsOut "#####################MGWB AIC ##########################\n";
	$best_aic = AICcalcs(\%ConflictHash, \%ParamsPerBipart, \$count_of_seqs, \$CountOfGenes, \$TreeOneLikelihood, \%TreeLikelihoods, \$TreeTwoLikelihood);
	
	#Delta AIC and AIC Weights
	print StatsOut "#####################MGWB Delta AIC#####################\n";
	DeltaAIC(\%ConflictHash, \%ParamsPerBipart, \$count_of_seqs, \$CountOfGenes, \$TreeOneLikelihood, \%TreeLikelihoods, \$best_aic, \$TreeTwoLikelihood);
	
	print "MGWB Analysis finished :) , however, this is the quickest though :(\n";
}

############
#MSWB Method
############
if($MSWB eq "True"){
	
	print StatsOut "##########################################\n";
	print StatsOut "MSWB Analysis (Calculated Across Matrix):\n";
	
	if($verbose eq "True"){
		
		print StatsOut "#####################################################################\n";
		print StatsOut "This is verbose mode the calculations for the tree without any form of\n";
		print StatsOut "MSWB calculation will also be shown and compared for delta AIC this\n";
		print StatsOut "is not factored in in the non verbose MSWB analysis!\n";
		print StatsOut "#####################################################################\n";
		
	}else{
		
		print StatsOut "This is non verbose mode the AIC scores of the genes\n";
		print StatsOut "are only for trees treated using the MSWB analysis\n";
		
	}
	
	%FrequencyOfSites = (); %TreesUsed = ();
	%AICHash = (); %VerbAICFreq = ();
	%NonVerbAICFreq = (); %ParamsOfAllGenes = ();
	%TotalBipartLikelihood = ();
	$TotalTreeOne = 0; $TotalTreeTwo = 0;
	#Condensed Genes Contains the condensed down to genes
	$CondensedGenes = CondenseToGene($PartFile, "SuperMatrixSiteLLs.SSLL",$CountOfGenes);
	
	if($IsItATest eq "True"){
		
		$analysis = 1;
	}else{
		
		$analysis = $#loc;
	}
	
	#Get something saying site X is from bipart Y and supports Tree Z
	foreach $i (0..$analysis){
		
		@gene_array = (); $q = $i + 1;
		%BipartLikes = (); @array = ();
		@sorted_array = (); @NewSorted = ();
		$top_likely = 0; $gene_len = 0;
		@split_array = ();
		#Get the basics of the Supermatrix
		@gene_array = split "\n", $CondensedGenes;
		@split_array = split " ", $gene_array[1];
		$TreeOneLikelihood = $split_array[$q];
		@split_array = split " ", $gene_array[2];
		$TreeTwoLikelihood = $split_array[$q];
		#print "$TreeOneLikelihood\t$TreeTwoLikelihood\n";

		#Get a HASH with Site => TopLikelihood:Tree:Bipartition
		%BipartLikes = SiteBipartSupport(\$loc[$i],\%TreeHash, \"SuperMatrixSiteLLs.SSLL");		
		($start, $end) = split "-", $loc[$i];
		
		#Get Frequency For each Bipartition as Best (Checks out)
		
		foreach $keys ($start..$end){
			
			@sorted_array = ();
			@array = split ",", $BipartLikes{$keys};
			@sorted_array = sort {$b <=> $a} @array;

			foreach $j (0..$#sorted_array){
				@NewSorted = split ":", @sorted_array[$j];
				if($j == 0){
					$FrequencyOfSites{$NewSorted[2]}++;
					$top_likely = $NewSorted[0];

				}else{
					if($top_likely == $NewSorted[0]){
						$FrequencyOfSites{$NewSorted[2]}++;
					}
				}
			}
		}

		#Get gene by gene stats
		%TotalTreesUsed = (); %GeneLikelihood = ();
		foreach $keys ($start..$end){
			
			@sorted_array = (); @array = ();
			@temp_array = (); %HASH = ();
			@array = split ",", $BipartLikes{$keys};
			
			#Get a Hash with bipartition => Likelihood:Tree
			foreach $i (0..$#array){
				@temp_array = split ":", $array[$i];
				$HASH{$temp_array[2]} .= "$temp_array[0]:$temp_array[1],";
			}
			
			%BipartAndLike = (); %BipartAndTree = ();
			#Get the top Likelihood for each site
			for $keys2 (sort keys %HASH){
				@temp_array = split ",", $HASH{$keys2};
				@sorted_array = sort {$b <=> $a} @temp_array;
				@temp_array = split ":", $sorted_array[0];
				$BipartAndLike{$keys2} = $temp_array[0];
				$BipartAndTree{$keys2} = $temp_array[1];
				$GeneLikelihood{$keys2} += $temp_array[0];
				$TotalTreesUsed{$keys2} .= "$temp_array[1],";
				$TreesUsed{$keys2} .= "$temp_array[1],";
				
			}
			#Add this data together for genes

			#@sorted_array = sort {$b <=> $a} @array;
			
		}
		#print Dumper(\%TotalTreesUsed);
		#print Dumper(\%GeneLikelihood);
		
		#Condense the tree to the number of parameters
		%ParamHash = (); %ParamsPerBipart = ();
		foreach $keys (sort {lc $a cmp lc $b} keys %TotalTreesUsed){
			@array = split ",", $TotalTreesUsed{$keys};
			foreach $k (0..$#array){
				if(exists $ParamHash{$array[$k]}){
				
				}else{
					$ParamHash{$array[$k]}++;
				}
			}
			$count = keys %ParamHash;
			$ParamsPerBipart{"$keys"} = $count;
			%ParamHash = ();
		}
		#If it's verbose give all stats for all genes
		%TotalParams = ();
		if($verbose eq "True"){
			
			$best_aic_of_all = 0;
			$z = $i + 1; @VerbAIC = ();
			print StatsOut "################Gene $z Parameters############\n";
			%TotalParams = MSWBparameters(\%ConflictHash, \%ParamsPerBipart, \$count_of_seqs, \1, \$verbose);
						
			print StatsOut "################Gene $z Likelihood############\n";
			LikelihoodCalcs(\$TreeOneLikelihood, \%ConflictHash, \%GeneLikelihood, \$TreeTwoLikelihood);
			$TotalTreeOne += $TreeOneLikelihood;
			$TotalTreeTwo += $TreeTwoLikelihood;
			#!!!!!!!!!!Check count of genes probably change to 1
			#Really review this!!!!!!!!!!!!!
			print StatsOut "#####################Gene $z AIC##########################\n";
			($best_aic_of_all, $best_aic_of_mswb) = MSWBAic(\%TotalParams, \%ConflictHash, \$TreeOneLikelihood, \%GeneLikelihood, \$TreeTwoLikelihood, \$count_of_seqs);
			
			print StatsOut "#####################Gene $z Delta AIC####################\n";
			@VerbAIC = MSWBdeltaAIC(\%TotalParams, \%ConflictHash, \%GeneLikelihood, \$TreeOneLikelihood, \$TreeTwoLikelihood, \$best_aic_of_all, \$count_of_seqs);
			
			$hold = 0; @only_some = ();
			foreach $i (0..$#VerbAIC){
				
				@array = split ":", $VerbAIC[$i];
				
				if ($i == 0){
					
					$hold = $array[0];
					$VerbAICFreq{$array[1]}++;
				}else{
					
					if($hold == $array[0]){
						
						$VerbAICFreq{$array[1]}++;
					}
				}
				
			}
		
		
		}
		#Contains Likelihoods
		#print Dumper(\%GeneLikelihood);
		#Contains Parameters
		#print Dumper(\%ParamsPerBipart);
		
		#Get all the parameters used
		%TotalParams = MSWBparameters(\%ConflictHash, \%ParamsPerBipart, \$count_of_seqs, \1);
		@GeneAIC = (); $aic = 0; @SortedGeneAIC = ();
		foreach $keys (sort {lc $a cmp lc $b} keys %ParamsPerBipart){
			
			$ParamsOfAllGenes{$keys} += ($ParamsPerBipart{$keys}*6);
			
		}
		#print Dumper(\%ParamsOfAllGenes);
		foreach $keys (sort {lc $a cmp lc $b} keys %GeneLikelihood){
			$TotalBipartLikelihood{$keys} += $GeneLikelihood{$keys};
		}
		foreach $keys (sort {lc $a cmp lc $b} keys %GeneLikelihood){
			
			$aic = (-2 * $GeneLikelihood{$keys}) + (2 * $TotalParams{$keys});
			push @GeneAIC, "$aic:$keys";
		}
		#print Dumper(\@GeneAIC);
		@SortedGeneAIC = sort {$a <=> $b} @GeneAIC;
		@array = ();
		foreach $i (0..$#SortedGeneAIC){
			@array = split ":", $SortedGeneAIC[$i];
			if($i == 0){
				$AICHash{$array[1]}++;
				$hold = $array[0];
			}else{
				if($hold == $array[0]){
				
					$AICHash{$array[1]}++;
				}
			}
		}
	}
	#print Dumper(\%AICHash);
	#print Dumper(\%TotalBipartLikelihood);
	
	#TreesUsed Contains Total trees used from the analysis
	#print Dumper(\%TreesUsed);
	if($verbose eq "True"){
		
		
		print StatsOut "################################################\n";
		print StatsOut "This is output unique to verbose mode!\n";
		print StatsOut "#############GENE STATS Verbose#################\n";
		foreach $keys (sort {lc $a cmp lc $b} keys %VerbAICFreq){
		
			if($keys == -1){
				print StatsOut "Number of Genes supporting Tree 1: $VerbAICFreq{$keys}\n";
			}elsif($keys == -2){
				print StatsOut "Number of Genes supporting Tree 2: $VerbAICFreq{$keys}\n";
			}else{
				print StatsOut "Relationship $keys\t$ConflictHash{$keys}:\t$VerbAICFreq{$keys}\n";
			}
		}
		#print Dumper(\%ParamsOfAllGenes);
		#print Dumper(\%ParamsPerBipart);
		#print Dumper(\%ParamsOfAllGenes);
		
		
		#ParamsPerBipart is the number of trees used per bipartition
		#ParamsOfAllGenes is the total parameters from genes
		
		print StatsOut "##########Parameters MSWB Verbose#################\n";
		$params = (($count_of_seqs * 2) - 3) + (($analysis + 1)* 6);
		print StatsOut "Total Parameters of Tree 1: $params\n";
		print StatsOut "Total Parameters of Tree 2: $params\n";
		foreach $keys (sort {lc $a cmp lc $b} keys %ParamsPerBipart){
			$parameters = ($ParamsPerBipart{$keys} * (($count_of_seqs * 2) - 3)) + ($ParamsOfAllGenes{$keys});
			print StatsOut "Relationship $keys\t$ConflictHash{$keys}:\t$parameters\n";
		}
		print StatsOut "#########Likelihood MSWB Verbose#################\n";
		@likelihood_array = (); @best_likelihood = ();
		print StatsOut "The likelihood of Tree 1: $TotalTreeOne\n";
		print StatsOut "The likelihood of Tree 2: $TotalTreeTwo\n";
		push @likelihood_array, $TotalTreeOne;
		push @likelihood_array, $TotalTreeTwo;
		foreach $keys (sort {lc $a cmp lc $b} keys %TotalBipartLikelihood){
	
			print StatsOut "Relationship $keys\t$ConflictHash{$keys}:\t$TotalBipartLikelihood{$keys}\n";
			push @likelihood_array, $TotalBipartLikelihood{$keys};
		}
		@best_likelihood = sort {$b <=> $a} @likelihood_array;
		print StatsOut "Best Likelihood: $best_likelihood[0]\n";
		$aic = 0; $best_aic = 0;
		@aic_array = (); @sort_aic = ();
		print StatsOut "################AIC MSWB Verbose#################\n";
		$aic = (-2 * $TotalTreeOne) + (2 * $params);
		print StatsOut "The AIC of Tree 1: $aic\n";
		push @aic_array, $aic;
		$aic = (-2 * $TotalTreeTwo) + (2 * $params);
		print StatsOut "The AIC of Tree 2: $aic\n";
		push @aic_array, $aic;
		foreach $keys (sort {lc $a cmp lc $b} keys %TotalBipartLikelihood){
			$parameters = ($ParamsPerBipart{$keys} * (($count_of_seqs * 2) - 3)) + ($ParamsOfAllGenes{$keys});
			$aic = (-2 * $TotalBipartLikelihood{$keys}) + (2 * $parameters);
			print StatsOut "Relationship $keys\t$ConflictHash{$keys}:\t$aic\n";
			$AICs{$keys} = $aic;
			push @aic_array, $aic;
		}
		@sort_aic = sort {$a <=> $b} @aic_array;
		print StatsOut "The Best AIC is: $sort_aic[0]\n";
		$best_aic = $sort_aic[0];
		$delta_aic = 0; $TotalDelta = 0;
		%Deltas = (); $tree1_delta = 0;
		$tree2_delta = 0;
		print StatsOut "############MSWB Delta AIC Scores Verbose############\n";
		$aic = (-2 * $TotalTreeOne) + (2 * $params);
		$delta_aic = $aic - $best_aic;
		$tree1_delta  = $delta_aic;
		$TotalDelta += exp(-0.5 * $delta_aic);
		print StatsOut "The Delta AIC of Tree 1: $delta_aic\n";
		$delta_aic = ($aic - $best_aic);
		$tree2_delta = $delta_aic;
		$TotalDelta += exp(-0.5 * $delta_aic);
		$aic = (-2 * $TotalTreeTwo) + (2 * $params);
		print StatsOut "The Delta AIC of Tree 2: $delta_aic\n";
		foreach $keys (sort {lc $a cmp lc $b} keys %AICs){
			$delta_aic = $AICs{$keys} - $best_aic;
			print StatsOut "Relationship $keys\t$ConflictHash{$keys}:\t$delta_aic\n";
			$Deltas{$keys} = $delta_aic;
			$TotalDelta += exp(-0.5 * $delta_aic);
		}
		$weight = 0;
		print StatsOut "###############MSWB AIC Weights Verbose################\n";
		$weight = ((exp(-0.5 *$tree1_delta)) / $TotalDelta);
		print StatsOut "The Weight of Tree 1: $weight\n";
		$weight = ((exp(-0.5 *$tree2_delta)) / $TotalDelta);
		print StatsOut "The Weight of Tree 2: $weight\n";
		foreach $keys (sort {lc $a cmp lc $b} keys %Deltas){
		
			$weight = (exp(-0.5 * $Deltas{$keys}) / $TotalDelta);
			print StatsOut "Relationship $keys\t$ConflictHash{$keys}:\t$weight\n";
		}
	}
	
	
	print StatsOut "################################################\n";
	print StatsOut "This is printed regardless of mode, the AIC scores\n";
	print StatsOut "Frequency's etc...Are only using MSWB stats\n";
	print StatsOut "################################################\n";

	print StatsOut "##################MSWB Frequency Gene Stats##############\n";
	foreach $keys (sort {lc $a cmp lc $b} keys %AICHash){
		print StatsOut "Relationship $keys\t$ConflictHash{$keys}:\t$AICHash{$keys}\n";
	}
	#Frequency of Bipartition
	print StatsOut "#####################MSWB Frequency Site Stats##############\n";
	FrequencyCalcs(\%FrequencyOfSites, \%ConflictHash);

	#Total Parameters
	print StatsOut "#####################MSWB Parameters##############\n";
	foreach $keys (sort {lc $a cmp lc $b} keys %ParamsOfAllGenes){
		$parameters = ($ParamsPerBipart{$keys} * (($count_of_seqs * 2) - 3)) + ($ParamsOfAllGenes{$keys});
		print StatsOut "Relationship $keys\t$ConflictHash{$keys}:\t$parameters\n";
	}
	#Total Likelihood
	@likelihood_array = (); @best_likelihood = ();
	print StatsOut "#####################MSWB Likelihoods##############\n";
	foreach $keys (sort {lc $a cmp lc $b} keys %TotalBipartLikelihood){
	
		print StatsOut "Relationship $keys\t$ConflictHash{$keys}:\t$TotalBipartLikelihood{$keys}\n";
		push @likelihood_array, $TotalBipartLikelihood{$keys};
	}
	@best_likelihood = sort {$b <=> $a} @likelihood_array;
	print StatsOut "Best Likelihood: $best_likelihood[0]\n";
	
	#Best AIC
	$aic = 0; @aic_array = (); %AICs = ();
	@sorted_aic = (); $best_aic = 0; @NonVerbAIC = ();
	@sorted_nonverb = ();
	print StatsOut "#####################MSWB AIC Scores##############\n";
	foreach $keys (sort {lc $a cmp lc $b} keys %TotalBipartLikelihood){
		$parameters = ($ParamsPerBipart{$keys} * (($count_of_seqs * 2) - 3)) + ($ParamsOfAllGenes{$keys});
		$aic = (-2 * $TotalBipartLikelihood{$keys}) + (2 * $parameters);
		print StatsOut "Relationship $keys\t$ConflictHash{$keys}:\t$aic\n";
		$AICs{$keys} = $aic;
		push @aic_array, $aic;
		
	}
	@sorted_aic = sort {$a <=> $b} @aic_array;
	$best_aic = $sorted_aic[0];
	$delta_aic = 0; $TotalDelta = 0;
	%Deltas = ();
	print StatsOut "Best AIC: $best_aic\n";
	print StatsOut "#################MSWB Delta AIC Scores#############\n";
	foreach $keys (sort {lc $a cmp lc $b} keys %AICs){
		$delta_aic = $AICs{$keys} - $best_aic;
		print StatsOut "Relationship $keys\t$ConflictHash{$keys}:\t$delta_aic\n";
		$Deltas{$keys} = $delta_aic;
		$TotalDelta += exp(-0.5 * $delta_aic);
	}
	$weight = 0;
	print StatsOut "###################MSWB AIC Weights################\n";
	foreach $keys (sort {lc $a cmp lc $b} keys %Deltas){
		
		$weight = (exp(-0.5 * $Deltas{$keys}) / $TotalDelta);
		print StatsOut "Relationship $keys\t$ConflictHash{$keys}:\t$weight\n";
	}

}


###############
#MGWB_BL Method
###############

if($MGWB_BL eq "True"){
	
	print StatsOut "################################################################\n";
	print StatsOut "#Doing the MGWB Branch Lengths included Analysis\n";
	if($verbose eq "True"){
		
		print StatsOut "#Running with verbose, lots of info will be printed\n";
	}
	$CondensedGenes = "";
	system("cp SuperMatrixSiteLLsBL.SSLL Holder.txt");
	$CondensedGenes = CondenseToGene($PartFile, "Holder.txt",$CountOfGenes);
	#print "$CondensedGenes";
	@GeneArray = ();
	@GeneArray = SuperMatrixGeneStats($CondensedGenes);
	#Likelihood:Tree:Bipartition it supports
	#print Dumper(\@GeneArray);
	#Summarize the results
	$top_likelihood = 0; %TopTree = ();
	%Frequency = (); %TreeLikelihoods = ();
	$TotalDelta = 0; @array = ();
	@StatsArray = ();%TaxaParams = ();
	$TaxaParams = 0; $top_stats = "";
	#Loop all the genes in the analysis
	foreach $i (0..$#GeneArray){
		$gene = $i + 1;
		@array = @{$GeneArray[$i]};
		$branches = $TotalTaxaHash{$i};
		$TaxaParams += (2*$branches - 3);
		#Start looping each gene
		foreach $j (0..$#array){
			@StatsArray = split ":", $array[$j];
			$TreeLikelihoods{"$StatsArray[2]"} += $StatsArray[0];
			#$TaxaParams{"$StatsArray[2]"} += (2*$branches - 3);
			
			if($verbose eq "True"){
				print StatsOut "Gene: $gene\tLikelihood: $StatsArray[0]\tFrom Tree: $StatsArray[1]\tSupports Bipartion: $StatsArray[2]\tTotal taxa in tree: $branches\n";
			}
			if($j == 0){
			
				$top_likelihood = $StatsArray[0];
				$top_stats .= "$gene\t$StatsArray[0]\t$StatsArray[2]\t$ConflictHash{$StatsArray[2]}\t$StatsArray[1]\t$branches\n";
				$Frequency{"$StatsArray[2]"}++;
			}else{
				#Get counts if multiple hits for the best likelihood
				if($top_likelihood == $StatsArray[0]){
					$top_stats .= "$gene\t$StatsArray[0]\t$StatsArray[2]\t$ConflictHash{$StatsArray[2]}\t$StatsArray[1]\t$branches\n";
					$Frequency{"$StatsArray[2]"}++;

				}
			}
			$TopTree{"$StatsArray[2]"} .= "$StatsArray[1],";
		

		}
	
	}
	#%TaxaParams should have the same number for all its the sum of the branches estimated
	#print Dumper(\%Frequency);
	$parameters = 0;
	print StatsOut "#############################################################\n";
	print StatsOut "#Summary Stats\n";
	print StatsOut "Gene\tLikelihood\tBipartition\tRelationship\tTree\tTaxaInTree\n";
	print StatsOut "$top_stats";
	#Frequency of bipartition
	print StatsOut "##################MGWB With BL Frequency Stats###############\n";
	FrequencyCalcs(\%Frequency, \%ConflictHash);
	
	print StatsOut "##################MGWB With BL Likelihoods###################\n";
	MGWB_BL_LikelihoodCalcs(\%TreeLikelihoods, \%ConflictHash);
	
	print StatsOut "##################MGWB With BL Parameters Used###############\n";
	$parameters = $TaxaParams + (6 * $CountOfGenes);
	print StatsOut "Total Parameters (same for all): $parameters\n";
	
	$aic = 0; @array = (); @sorted_array = (); $best_aic = 0;
	print StatsOut "##################MGWB With BL AIC Scores###############\n";
	for $keys (sort {lc $a cmp lc $b} keys %TreeLikelihoods){
		
		$aic = (-2 * $TreeLikelihoods{$keys}) + (2 * $parameters);
		print StatsOut "Relationship $keys\t$ConflictHash{$keys}:\t$aic\n";
		push @array, $aic;
	}
	@sorted_array = sort {$a <=> $b} @array;
	$best_aic = $sorted_array[0];
	print StatsOut "Best AIC: $sorted_array[0]\n";
	
	$aic = 0; $delta_aic = 0; $TotalDelta = 0;
	print StatsOut "##############MGWB With BL Delta AIC Scores###############\n";
	for $keys (sort {lc $a cmp lc $b} keys %TreeLikelihoods){
		$aic = (-2 * $TreeLikelihoods{$keys}) + (2 * $parameters);
		$delta_aic = $aic - $best_aic;
		$TotalDelta += exp(-0.5 * $delta_aic);
		print StatsOut "Relationship $keys\t$ConflictHash{$keys}:\t$delta_aic\n";
		
	}
	$weight = 0;
	print StatsOut "##############MGWB With BL AIC Weights###################\n";
		for $keys (sort {lc $a cmp lc $b} keys %TreeLikelihoods){
		$aic = (-2 * $TreeLikelihoods{$keys}) + (2 * $parameters);
		$delta_aic = $aic - $best_aic;
		$weight = (exp(-0.5 * $delta_aic) / $TotalDelta);
		print StatsOut "Relationship $keys\t$ConflictHash{$keys}:\t$weight\n";
		
	}
}

###############
#MSWB_BL Method
###############

if($MSWB_BL eq "True"){
	
	print StatsOut "################################################################\n";
	print StatsOut "#Doing the MSWB Branch Lengths included Analysis\n";
	if($verbose eq "True"){
		
		print StatsOut "#Running with verbose, lots of info will be printed\n";
	}
	$CondensedGenes = "";
	system("cp SuperMatrixSiteLLsBL.SSLL Holder.txt");
	$CondensedGenes = CondenseToGene($PartFile, "Holder.txt",$CountOfGenes);
	if($IsItATest eq "True"){
		print StatsOut "#Running Test Mode\n";
		$analysis = 1;
	}else{
		$analysis = $#loc;
	}
	
	%FrequencyOfSites = (); %TotalTreeLikelihood = ();
	%AllTheParams = (); %AllAicFreqs = ();
	foreach $i (0..$analysis){
		
		#print Dumper(\%TreeHash);
		
		#BipartLikes has Likelihood:Tree:Bipartition
		%BipartLikes = SiteBipartSupport(\$loc[$i],\%TreeHash, \"Holder.txt");
		($start, $end) = split "-", $loc[$i];
		foreach $keys ($start..$end){
			
			
			@sorted_array = (); @array = ();
			@array = split ",", $BipartLikes{$keys};
			@sorted_array = sort {$b <=> $a} @array;

			
			#Will get frequency of every cite supporting a given tree
			foreach $j (0..$#sorted_array){
				@NewSorted = split ":", @sorted_array[$j];
				if($j == 0){
					$FrequencyOfSites{$NewSorted[2]}++;
					$top_likely = $NewSorted[0];

				}else{
					if($top_likely == $NewSorted[0]){
						$FrequencyOfSites{$NewSorted[2]}++;
					}
				}
			}
			
			
			#print Dumper(\%FrequencyOfSites);
				
		}
		%TotalTreesUsed = (); %GeneLikelihood = ();
		foreach $keys ($start..$end){
			
			@sorted_array = (); @array = ();
			@temp_array = (); %HASH = ();
			@array = split ",", $BipartLikes{$keys};
			
			#Get a Hash with bipartition => Likelihood:Tree
			foreach $i (0..$#array){
				@temp_array = split ":", $array[$i];
				$HASH{$temp_array[2]} .= "$temp_array[0]:$temp_array[1],";
			}
			
			%BipartAndLike = (); %BipartAndTree = ();
			#Get the top Likelihood for each site
			for $keys2 (sort keys %HASH){
				@temp_array = split ",", $HASH{$keys2};
				@sorted_array = sort {$b <=> $a} @temp_array;
				@temp_array = split ":", $sorted_array[0];
				$BipartAndLike{$keys2} = $temp_array[0];
				$BipartAndTree{$keys2} = $temp_array[1];
				$GeneLikelihood{$keys2} += $temp_array[0];
				$TotalTreesUsed{$keys2} .= "$temp_array[1],";
				$TreesUsed{$keys2} .= "$temp_array[1],";
				$TotalTreeLikelihood{$keys2} += $temp_array[0];
				
			}
			#Add this data together for genes

			#@sorted_array = sort {$b <=> $a} @array;
			
		}

		#print Dumper(\%BipartAndLike);
		%ParamHash = (); %ParamsPerBipart = ();
		foreach $keys (sort {lc $a cmp lc $b} keys %TotalTreesUsed){
			@array = split ",", $TotalTreesUsed{$keys};
			foreach $k (0..$#array){
				if(exists $ParamHash{$array[$k]}){
				
				}else{
					$ParamHash{$array[$k]}++;
				}
			}
			$count = keys %ParamHash;
			$ParamsPerBipart{"$keys"} = $count;
			%ParamHash = ();
		}
		#print Dumper(\%ParamsPerBipart);
		
		#Do everything one gene at a time
		%TotalParams = ();
		if($verbose eq "True"){
			
			$best_aic_of_all = 0;
			$z = $i + 1; @VerbAIC = ();
			print StatsOut "################Gene $z Parameters############\n";
			print $counts_of_seqs;
			%TotalParams = MSWB_BLparameters(\%ConflictHash, \%ParamsPerBipart, \$TotalTaxaHash{$i}, \1, \$verbose);
			
						
			print StatsOut "################Gene $z Likelihood############\n";
			foreach $keys (sort {lc $a cmp lc $b} keys %GeneLikelihood){
				
				print StatsOut "Relationship $keys $ConflictHash{$keys} $GeneLikelihood{$keys}\n";
				
			}
			#!!!!!!!!!!Check count of genes probably change to 1
			#Really review this!!!!!!!!!!!!!
			$aic = 0; $best_aic = 0; @array = (); @sorted_array = ();
			print StatsOut "#####################Gene $z AIC##########################\n";

			foreach $keys (sort {lc $a cmp lc $b} keys %GeneLikelihood){
					
					$aic = (-2 * $GeneLikelihood{$keys}) + (2 * $TotalParams{$keys});
					print StatsOut "Relationship $keys\t$ConflictHash{$keys}:\t$aic\n";
					push @array, $aic;
			}
			@sorted_array = sort {$a <=> $b} @array;
			$best_aic = $sorted_array[0];
			print StatsOut "Best AIC is: $best_aic\n";
			$delta_aic = 0; $TotalDelta = 0; $aic_of_bipart = 0;
			print StatsOut "#####################Gene $z Delta AIC####################\n";
			foreach my $keys (sort {lc $a cmp lc $b} keys %GeneLikelihood){
		
				$aic_of_bipart = (-2 * $GeneLikelihood{$keys}) + (2 * $TotalParams{$keys});
				$delta_aic = $aic_of_bipart - $best_aic;
				$TotalDelta += exp(-0.5 * $delta_aic);
				print StatsOut "Relationship $keys\t$ConflictHash{$keys}:\t$delta_aic\n";
			}
			$delta_aic = 0; $aic_of_bipart = 0;
			print StatsOut "#####################Gene $z AIC weight###################\n";
			foreach my $keys (sort {lc $a cmp lc $b} keys %GeneLikelihood){
				$aic_of_bipart = (-2 * $GeneLikelihood{$keys}) + (2 * $TotalParams{$keys});
				$delta_aic = $aic_of_bipart - $best_aic;
				$weight = (exp(-0.5 * $delta_aic) / $TotalDelta);
				print StatsOut "Relationship $keys\t$ConflictHash{$keys}:\t$weight\n";
					
			}
		
		}
		%TotalParams = MSWB_BLparameters(\%ConflictHash, \%ParamsPerBipart, \$TotalTaxaHash{$i}, \1, \"False");
		#Get the best gene a.k.a gene with best AIC frequencies %AllAicFreqs
		@AICmax = (); @sorted_aic_max = ();
		foreach $keys (sort {lc $a cmp lc $b} keys %GeneLikelihood){
					
			$aic = (-2 * $GeneLikelihood{$keys}) + (2 * $TotalParams{$keys});
			push @AICmax, "$aic:$keys";
		}
		@sorted_aic_max = sort {$a <=> $b} @AICmax;
		#print Dumper(\@AICmax);
		#print Dumper(\@sorted_aic_max);
		@array = (); $best = 0;
		foreach $keys (0..$#sorted_aic_max){
			@array = split ":", $sorted_aic_max[$keys];
			if($keys == 0){
				
				$AllAicFreqs{$array[1]}++;
				$best = $array[0];
				#print "$best\n";
			}else{
				if($best == $array[0]){
					$AllAicFreqs{$array[1]}++;
				}	
			}
			
				
		}
		
		
		#print Dumper(\%FrequencyOfSites);
		#print Dumper(\%TotalParams);
		foreach my $keys (sort {lc $a cmp lc $b} keys %TotalParams){
			
			$AllTheParams{$keys} += $TotalParams{$keys};
		}
	}

	#print Dumper(\%AllAicFreqs);
	print StatsOut "############################################################\n";
	print StatsOut "#Summary of everything\n";
	#Frequency of Bipartition
	print StatsOut "##################MSWB BL Frequency Site Stats##############\n";
	FrequencyCalcs(\%FrequencyOfSites, \%ConflictHash);
	print StatsOut "##################MSWB BL Frequency of Best AIC#############\n";
	foreach my $keys (sort {lc $a cmp lc $b} keys %AllAicFreqs){
		print StatsOut "Relationship $keys\t$ConflictHash{$keys}:\t$AllAicFreqs{$keys}\n";
	}
	print StatsOut "##################MSWB BL Parameters Used###################\n";
	foreach my $keys (sort {lc $a cmp lc $b} keys %AllTheParams){
		print StatsOut "Relationship $keys\t$ConflictHash{$keys}:\t$AllTheParams{$keys}\n";
	}
	print StatsOut "##################MSWB BL Likelihood Stats##################\n";
	foreach my $keys (sort {lc $a cmp lc $b} keys %TotalTreeLikelihood){
		print StatsOut "Relationship $keys\t$ConflictHash{$keys}:\t$TotalTreeLikelihood{$keys}\n";
	}
	$aic = 0; $best_aic = 0; @array = (); @sorted_array = ();
	print StatsOut "##################MSWB BL AIC Scores########################\n";
	foreach my $keys (sort {lc $a cmp lc $b} keys %TotalTreeLikelihood){
		$aic = (-2 * $TotalTreeLikelihood{$keys}) + (2 * $TotalParams{$keys});
		print StatsOut "Relationship $keys\t$ConflictHash{$keys}:\t$aic\n";
		push @array, $aic;
	}
	@sorted_array = sort {$a <=> $b} @array;
	$best_aic = $sorted_array[0];
	print StatsOut "Best AIC is: $best_aic\n";
	$TotalDelta = 0;
	print StatsOut "##################MSWB BL Delta AIC########################\n";
	foreach my $keys (sort {lc $a cmp lc $b} keys %TotalTreeLikelihood){
		
		$aic_of_bipart = (-2 * $TotalTreeLikelihood{$keys}) + (2 * $TotalParams{$keys});
		$delta_aic = $aic_of_bipart - $best_aic;
		$TotalDelta += exp(-0.5 * $delta_aic);
		print StatsOut "Relationship $keys\t$ConflictHash{$keys}:\t$delta_aic\n";
	}
	$weight = 0;
	print StatsOut "#####################MSWB BL AIC weight###################\n";
	foreach my $keys (sort {lc $a cmp lc $b} keys %TotalTreeLikelihood){
		$aic_of_bipart = (-2 * $TotalTreeLikelihood{$keys}) + (2 * $TotalParams{$keys});
		$delta_aic = $aic_of_bipart - $best_aic;
		$weight = (exp(-0.5 * $delta_aic) / $TotalDelta);
		print StatsOut "Relationship $keys\t$ConflictHash{$keys}:\t$weight\n";				
	}
}

print "Files Made during process have been moved to a folder called: $SuperMatrix\_OutfilesAndStuff\n";
print "Data Has been written to: $outfile\n";

system("mkdir $SuperMatrix\_OutfilesAndStuff");
system("mv Holder.txt bp.log phyx.logfile RAxML_*EX_SSLL SuperMatrixSiteLLsBL.SSLL SuperMatrixSiteLLs.SSLL temp.fa temp.fa.reduced temp.log temptesttre TempTree.tre Unique.tre $SuperMatrix\_OutfilesAndStuff");
