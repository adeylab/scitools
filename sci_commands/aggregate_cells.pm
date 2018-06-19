package sci_commands::aggregate_cells;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("aggregate_cells");

sub aggregate_cells {

@ARGV = @_;
# Defaults
$range_default = "1-2";
$oversample = 1;
$xdim = 1;
$ydim = 2;
$aggN = 15;
$minN = 5;
$theme = "Clean";

getopts("O:D:N:A:a:C:c:x:y:R:XK:T:n:v:", \%opt);

$die2 = "
scitools dims-aggregate [options] [input dims/lambda/values]
   or    aggregate

Aggregates cells via k-means or similar. Produces a merged dims for the
new aggregate cells, and a merge annotation file (can merge matrix
using this new annotation file)

Options:
   -O   [STR]   Output prefix (default is [input].aggregate.annot)
   -D   [RNG]   Range of dimensions to aggregate from (range format)
                (e.g. 1-5,6,8; def = $range_default; 1 if lambda input)
   -N   [INT]   Target number of cells to aggregate (def = $aggN)
   -n   [INT]   Min number of cells in cluster (def = $minN)
   -v   [FLT]   Oversampling value (ie mean assignments for each cell)
                Note: this will create an annot file with comma separated
                annotation memberships if >1. (def = $oversample)
   -K   [INT]   Number of clusters (overrides -N to be the min/aggregate)
                This will utilize standard k-means to aggregate cells
                For lambda clustering it will evenly space them
   -A   [STR]   Annotation file (will only merge within annot)
   -a   [STR]   Comma separated list of annotations to include
                  (requires -A to be specified)

Plotting Options (def is input dims file coordinates - only for dims files):
   -C   [STR]   Color coding file (annot (tab) #hexColor)
   -c   [STR]   Color coding string
                  Annot=#hexColor,Annot2=#hexColor
   -x   [INT]   X-dimension to plot (def = $xdim)
   -y   [INT]   Y-dimension to plot (def = $ydim)
   -T   [STR]   Theme: (def = $theme)
                  Clean = no axis lines (for tSNE)
                  Reg = regular, include axes and grid (for PCA)
   -R   [STR]   Rscript call (def = $Rscript)
   -X           Retain intermediate files (def = delete)

";

#die "ERROR: This function is under development and cannot be used in its current form!\n";

if (!defined $ARGV[0]) {die $die2};
if (defined $opt{'a'} && !defined $opt{'A'}) {die "\nMust provide an annotaiton file (-A) if specifying annotations to plot (-a)!\n$die2"};
if (defined $opt{'C'} && defined $opt{'c'}) {die "\nSpecify either a color string (-c) or a color coding file (-C), not both!\n$die2"};
if (!defined $opt{'D'}) {$opt{'D'} = $range_default};
if (defined $opt{'N'}) {$aggN = $opt{'N'}};
if (defined $opt{'n'}) {$minN = $opt{'n'}};
if (defined $opt{'T'}) {$theme = $opt{'T'}};
if (defined $opt{'v'}) {$oversample = $opt{'v'}};

if ($ARGV[0] =~ /\.(lambda|values)$/) {
	print STDERR "INFO: A lambda or values file was detected - will aggregate in one dimension\n";
	if (!defined $opt{'D'}) {$opt{'D'} = 1};
	read_values($ARGV[0]);
	if (defined $opt{'A'}) {
		read_annot($opt{'A'});
	} else {
		foreach $cellID (keys %CELLID_value) {
			$CELLID_annot{$cellID} = "Agg";
		}
		$ANNOT_include{"Agg"} = 1;
	}
} else {
	read_dims($ARGV[0]);
	if (defined $opt{'A'}) {
		read_annot($opt{'A'});
	} else {
		foreach $cellID (keys %CELLID_DIMS) {
			$CELLID_annot{$cellID} = "Agg";
		}
		$ANNOT_include{"Agg"} = 1;
	}
}

if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
$opt{'O'} =~ s/\.(dims|lambda|values)$//;

if (defined $opt{'R'}) {$Rscript = $opt{'R'}};

if (defined $opt{'a'}) {
	@ANNOT_LIST = split(/,/, $opt{'a'});
	foreach $annot (@ANNOT_LIST) {
		$ANNOT_include{$annot} = 1;
	}
} else {
	foreach $annot (keys %ANNOT_count) {
		$ANNOT_include{$annot} = 1;
	}
}

if (defined $opt{'C'}) {read_color_file($opt{'C'})};
if (defined $opt{'c'}) {read_color_string($opt{'c'})};

if (defined $opt{'x'}) {$xdim = $opt{'x'}};
if (defined $opt{'y'}) {$ydim = $opt{'y'}};
$xpos = $xdim-1; $ypos = $ydim-1;

read_ranges($opt{'D'});

if ($ARGV[0] =~ /\.(lambda|values)$/) {
	# figure out even spaced vs. even cell N
	# if -K, then will be evenly spaced on lambda - uniform for all annotations
	# elsif -N, then will be even number of cells - varable for annotations
	if (defined $opt{'K'}) { # find centroids that will be constant across all annotations
		%CLUST_global_centroid = ();
		$clust_span = $value_range/$opt{'K'};
		$center = $clust_span/2;
		$clustID = 0;
		$CLUST_global_centroid{$clustID} = $center;
		for ($clustID = 1; $clustID < $opt{'K'}; $clustID++) {
			$center += $clust_span;
			$CLUST_global_centroid{$clustID} = $center;
		}
	}
}

open CNT, ">$opt{'O'}.centroids.dims";
open OUT, ">$opt{'O'}.aggregate.annot";

# loop through annotations individually
foreach $annot (keys %ANNOT_include) {
if ($ARGV[0] =~ /\.(lambda|values)$/) {
	# lambda file
	# setup new CLUST center hash for the annot - all can be annot-only, no need for global to be kept
	%CLUST_center = (); $assignment_count = 0; %CLUST_assignments = (); %CELLID_initial_cluster = (); %CLUST_cellIDs = ();
	# setup cluster centroids based on K or N
	if (defined $opt{'K'}) { # copy over the clusters and add in the annot
		foreach $clustID (keys %CLUST_global_centroid) {
			$CLUST_center{$annot."_".$clustID} = $CLUST_global_centroid{$clustID};
		}
		# now do initial assignments based on closest center
		foreach $cellID (keys %CELLID_value) {
			if ($CELLID_annot{$cellID} eq $annot) {
				$minDist = 1e9; $winAsn = "NA";
				foreach $clustID (keys %CLUST_center) {
					$dist = abs($CELLID_value{$cellID} - $CLUST_center{$clustID});
					if ($dist<$minDist) {
						$minDist = $dist;
						$winAsn = $clustID;
					}
				}
				$CELLID_initial_cluster{$cellID} = $winAsn;
				push @{$CLUST_cellIDs{$winAsn}}, $cellID;
				$CLUST_assignments{$winAsn}++;
				$assignment_count++;
			}
		}
		# now check for reaching min n, remove other clusters and do orphan assignment
		%CELLID_orphan = ();
		foreach $clustID (keys %CLUST_center) {
			if ($CLUST_assignments{$clustID} < $minN) {
				print STDERR "INFO: Cluster $clustID has $CLUST_assignments{$clustID} cells assigned which is < $minN, and will be removed and have cells re-assigned.\n";
				delete $CLUST_center{$clustID};
				foreach $cellID (@{$CLUST_cellIDs{$clustID}}) {
					$CELLID_orphan{$cellID} = 1;
				}
				delete $CLUST_cellIDs{$clustID};
			}
		}
		foreach $cellID (keys %CELLID_orphan) {
			$minDist = 1e9; $winAsn = "NA";
			foreach $clustID (keys %CLUST_center) {
				$dist = abs($CELLID_value{$cellID} - $CLUST_center{$clustID});
				if ($dist<$minDist) {
					$minDist = $dist;
					$winAsn = $clustID;
				}
			}
			$CELLID_initial_cluster{$cellID} = $winAsn;
			push @{$CLUST_cellIDs{$winAsn}}, $cellID;
			$CLUST_assignments{$winAsn}++;
		}
	} else {
		# go through cells and calculate the center for each grouping - also do initial assignment
		$clustNum = 0; $clust_sum = 0; $clust_memberCT = 0; $clustID = $annot."_".$clustNum; $annot_cellCT_in_lambda = 0;
		foreach $cellID (keys %CELLID_value) {
			if ($CELLID_annot{$cellID} eq $annot) {
				$annot_cellCT_in_lambda++;
			}
		}
		$cluster_ct = int($annot_cellCT_in_lambda/$aggN);
		if ($cluster_ct < 2) {
			$cluster_ct = 2;
			print STDERR "WARNING: The N value specified, $aggN, results in cluster numbers less than 2 (cells in $annot = $annot_cellCT_in_lambda), setting to 2\n";
		}
		$target_assignments = int($annot_cellCT_in_lambda/$cluster_ct);
		foreach $cellID (sort {$CELLID_value{$a}<=>$CELLID_value{$b}} keys %CELLID_value) {
			if ($CELLID_annot{$cellID} eq $annot) {
				$clust_memberCT++;
				$clust_sum += $CELLID_value{$cellID};
				$CELLID_initial_cluster{$cellID} = $clustID;
				push @{$CLUST_cellIDs{$clustID}}, $cellID;
				$assignment_count++;
				if ($clust_memberCT >= $target_assignments) {
					$center = $clust_sum/$clust_memberCT;
					$CLUST_center{$clustID} = $center;
					$clustNum++; $clust_sum = 0; $clust_memberCT = 0; $clustID = $annot."_".$clustNum; 
				}
			}
		}
		# last one
		$center = $clust_sum/$clust_memberCT;
		$CLUST_center{$clustID} = $center;
	}
	
	# now check if there is oversampling
	if ($oversample>1) {
		%CELLID_cluster = ();
		if (defined $opt{'K'}) { # oversample for each cluster
			foreach $clustID (keys %CLUST_center) {
				$target_assignments = int(($CLUST_assignments{$clustID}*$oversample)+1);
				%CELLID_dist = ();
				foreach $cellID (keys %CELLID_value) {
					if ($CELLID_annot{$cellID} eq $annot && $CELLID_initial_cluster{$cellID} ne $clustID) {
						$CELLID_dist{$cellID} = abs($CELLID_value{$cellID} - $CLUST_center{$clustID});
					} elsif ($CELLID_annot{$cellID} eq $annot && $CELLID_initial_cluster{$cellID} eq $clustID) {
						$CELLID_cluster{$cellID} .= "$clustID,";
					}
				}
				foreach $cellID (sort {$CELLID_dist{$a}<=>$CELLID_dist{$b}} keys %CELLID_dist) {
					if ($CLUST_assignments{$clustID}<$target_assignments) {
						$CELLID_cluster{$cellID} .= "$clustID,";
						$CLUST_assignments{$clustID}++;
						push @{$CLUST_cellIDs{$clustID}}, $cellID;
					}
				}
			}
			foreach $clustID (keys %CLUST_center) { # no center changes, so print centroids
				print CNT "$clustID\t$CLUST_center{$clustID}\n";
			}
		} else {
			$added_assignments = ($annot_cellCT_in_lambda*($oversample-1));
			$target_assignments += int(($added_assignments/$clustNum)+1);
			foreach $clustID (keys %CLUST_center) {
				%CELLID_dist = ();
				foreach $cellID (keys %CELLID_value) {
					if ($CELLID_annot{$cellID} eq $annot && $CELLID_initial_cluster{$cellID} ne $clustID) {
						$CELLID_dist{$cellID} = abs($CELLID_value{$cellID} - $CLUST_center{$clustID});
					} elsif ($CELLID_annot{$cellID} eq $annot && $CELLID_initial_cluster{$cellID} eq $clustID) {
						$CELLID_cluster{$cellID} .= "$clustID,";
					}
				}
				foreach $cellID (sort {$CELLID_dist{$a}<=>$CELLID_dist{$b}} keys %CELLID_dist) {
					if ($CLUST_assignments{$clustID}<$target_assignments) {
						$CELLID_cluster{$cellID} .= "$clustID,";
						$CLUST_assignments{$clustID}++;
						push @{$CLUST_cellIDs{$clustID}}, $cellID;
					}
				}
			}
			# remake centroid & print
			foreach $clustID (keys %CLUST_center) {
				$clust_sum = 0; $clust_memberCT = 0;
				foreach $cellID (@{$CLUST_cellIDs{$clustID}}) {
					$clust_sum+=$CELLID_value{$cellID};
					$clust_memberCT++;
				}
				$CLUST_center{$clustID} = $clust_sum/$clust_memberCT;
				print CNT "$clustID\t$CLUST_center{$clustID}\n";
			}
		}
		# print out annot
		foreach $cellID (keys %CELLID_cluster) {
			$cluster_set = $CELLID_cluster{$cellID};
			$cluster_set =~ s/,$//;
			print OUT "$cellID\t$cluster_set\n";
		}
	} else { # no oversampling
		foreach $clustID (keys %CLUST_center) {
			print CNT "$clustID\t$CLUST_center{$clustID}\n";
			foreach $cellID (@{$CLUST_cellIDs{$clustID}}) {
				print OUT "$cellID\t$clustID\n";
			}
		}
	}
	
} else {
	# Dimensions file
	
	%CLUST_include = ();
	%CELLID_initial_cluster = ();
	%CLUST_cellCT = ();
	$ANNOT_count{$annot} = 0;
	
	foreach $cellID (keys %CELLID_DIMS) {
		if ($CELLID_annot{$cellID} eq $annot) {$ANNOT_count{$annot}++};
	}
	
	if (defined $opt{'K'}) {
		$K = $opt{'K'};
		$aggN = int($ANNOT_count{$annot}/$K);
	} else {
		$K = int(($ANNOT_count{$annot}/$aggN)+0.5);
		if ($K<1) {$K=1};
	}

	open AO, ">$opt{'O'}.$annot.agg_dims";
	foreach $cellID (keys %CELLID_DIMS) {
		if ($CELLID_annot{$cellID} eq $annot) {
			print AO "$cellID";
			for ($dimID = 0; $dimID < @RANGE_VALUES; $dimID++) {
				print AO "\t$CELLID_DIMS{$cellID}[$RANGE_VALUES[$dimID]]";
			}
			print AO "\n";
		}
	} close AO;
	open R, ">$opt{'O'}.$annot.agg_dims.r";
	if ($ARGV[0] =~ /\.(lambda|values)$/) {
	print R "
D<-read.table(\"$opt{'O'}.$annot.agg_dims\",row.names=1)";
	} else {
	print R "
D<-read.table(\"$opt{'O'}.$annot.agg_dims\",row.names=1)[,$range_R_set]";
	}
	print R"
K<-kmeans(D,$K)
ANN<-as.matrix(K\$cluster)
ANN[,1]<-sub(\"\^\", \"K\", ANN[,1])
write.table(ANN,file=\"$opt{'O'}.$annot.agg_dims.annot\",col.names=FALSE,row.names=TRUE,quote=FALSE,sep=\"\\t\")
"; close R;
	system("$Rscript $opt{'O'}.$annot.agg_dims.r 2>/dev/null");

	# find passing clusters (needed for separation)
	open IN, "$opt{'O'}.$annot.agg_dims.annot";
	while ($l = <IN>) {
		chomp $l;
		($cellID,$clust_raw) = split(/\t/, $l);
		$cluster = "$annot\_$clust_raw";
		$CLUST_cellCT{$cluster}++;
	} close IN;
	foreach $cluster (keys %CLUST_cellCT) {
		if ($CLUST_cellCT{$cluster}>=$minN) {
			$CLUST_include{$cluster} = 1;
		}
	}

	# build assignments
	open IN, "$opt{'O'}.$annot.agg_dims.annot";
	while ($l = <IN>) {
		chomp $l;
		($cellID,$clust_raw) = split(/\t/, $l);
		$cluster = "$annot\_$clust_raw";
		if (defined $CLUST_include{$cluster}) {
			$CELLID_initial_cluster{$cellID} = "$cluster";
		}
		push @{$CLUSTER_CELLS{$cluster}}, $cellID;
	} close IN;
	
	# now compute centroids
	foreach $cluster (keys %CLUST_include) {
		@DIM_SUMS = ();
		@{$CLUSTER_centroidDims{$cluster}} = ();
		foreach $cellID (@{$CLUSTER_CELLS{$cluster}}) {
			for ($dimID = 0; $dimID < @RANGE_VALUES; $dimID++) {
				$DIM_SUMS[$dimID] += $CELLID_DIMS{$cellID}[$RANGE_VALUES[$dimID]];
			}
		}
		for ($dimID = 0; $dimID < @DIM_SUMS; $dimID++) {
			$mean_dim = $DIM_SUMS[$dimID]/$CLUST_cellCT{$cluster};
			$CLUSTER_centroidDims{$cluster}[$dimID] = $mean_dim;
		}
	}
	
	# now assign cells based on the centroids & factor in oversampling
	if ($oversample>1) {
		foreach $cluster (keys %CLUST_include) {
			# calculate cell to centroid distances
			%CELLID_distance = ();
			foreach $cellID (keys %CELLID_initial_cluster) {
				$dist_sum = 0;
				for ($dimID = 0; $dimID < @{$CLUSTER_centroidDims{$cluster}}; $dimID++) {
					$dist_sum += ($CLUSTER_centroidDims{$cluster}[$dimID] - $CELLID_DIMS{$cellID}[$RANGE_VALUES[$dimID]])**2;
				}
				$CELLID_distance{$cellID} = sqrt($dist_sum);
			}
			
			# now calculate target number of cells in cluster
			if ($CLUST_cellCT{$cluster}>=($oversample*$aggN)) {
				$CLUST_memberCT{$cluster} = $CLUST_cellCT{$cluster};
			} else {
				$CLUST_memberCT{$cluster} = int($oversample*$aggN);
			}
			
			# now assign & recompute centroids
			$added = 0; @DIM_SUMS = ();
			foreach $cellID (sort {$CELLID_distance{$a}<=>$CELLID_distance{$b}} keys %CELLID_distance) {
				if ($added < $CLUST_memberCT{$cluster}) {
					for ($dimID = 0; $dimID < @RANGE_VALUES; $dimID++) {
						$DIM_SUMS[$dimID] += $CELLID_DIMS{$cellID}[$RANGE_VALUES[$dimID]];
					}
					$CELLID_cluster{$cellID} .= "$cluster,";
					$added++;
				}
			}
			for ($dimID = 0; $dimID < @DIM_SUMS; $dimID++) {
				$mean_dim = $DIM_SUMS[$dimID]/$CLUST_memberCT{$cluster};
				$CLUSTER_centroidDims{$cluster}[$dimID] = $mean_dim;
			}
		}
	} else {
		foreach $cellID (keys %CELLID_initial_cluster) {
			$CELLID_cluster{$cellID} = $CELLID_initial_cluster{$cellID};
		}
	}
	
	foreach $cluster (keys %CLUST_include) {
		print CNT "$cluster";
		for ($dimID = 0; $dimID < @{$CLUSTER_centroidDims{$cluster}}; $dimID++) {
			print CNT "\t$CLUSTER_centroidDims{$cluster}[$dimID]";
		} print CNT "\n";
	}
	foreach $cellID (keys %CELLID_initial_cluster) {
		$CELLID_cluster{$cellID} =~ s/,$//;
		if ($CELLID_cluster{$cellID} ne "") {
			print OUT "$cellID\t$CELLID_cluster{$cellID}\n";
		}
	}
	
	if (!defined $opt{'X'}) {
		system("rm -f $opt{'O'}.$annot.agg_dims.annot $opt{'O'}.$annot.agg_dims.r $opt{'O'}.$annot.agg_dims");
	}
}
}
close CNT; close OUT;

# Plot the projections of cells to their centroids
if ($ARGV[0] =~ /\.(lambda|values)$/) {
	# labmda data plots?
	# NOT CURRENTLY SUPPORTED
} else {
# Dimensions data & plot
open PD, ">$opt{'O'}.aggregate_data";
foreach $cellID (keys %CELLID_cluster) {
	@CLUSTER_SET = split(/,/, $CELLID_cluster{$cellID});
	foreach $cluster (@CLUSTER_SET) {
		print PD "$cellID:$cluster\t$CELLID_annot{$cellID}\t$CELLID_DIMS{$cellID}[$xdim]\t$CLUSTER_centroidDims{$cluster}[$xpos]\t$CELLID_DIMS{$cellID}[$ydim]\t$CLUSTER_centroidDims{$cluster}[$ypos]\n";
	}
} close PD;

open R, ">$opt{'O'}.aggregate_data.r";
print R "
library(ggplot2)
CELLS<-read.table(\"$opt{'O'}.aggregate_data\",row.names=1)
colnames(CELLS)<-c(\"annot\",\"xdim\",\"centx\",\"ydim\",\"centy\")
PLT<-ggplot() +";

if (!defined $opt{'c'} && !defined $opt{'C'} && !defined $opt{'A'}) {
	print R "
	geom_segment(aes(x=CELLS\$xdim,xend=CELLS\$centx,y=CELLS\$ydim,yend=CELLS\$centy),size=0.5,color=\"lightsteelblue4\",alpha=0.25) +
	geom_point(aes(CELLS\$xdim,CELLS\$ydim),color=\"lightsteelblue4\",size=0.5) +
	geom_point(aes(CELLS\$centx,CELLS\$centy),color=\"gray20\",size=1.75) +
	geom_point(aes(CELLS\$centx,CELLS\$centy),color=\"lightsteelblue4\",size=1) +";
} elsif (!defined $opt{'V'}) {
	print R "
	geom_segment(aes(x=CELLS\$xdim,xend=CELLS\$centx,y=CELLS\$ydim,yend=CELLS\$centy,color=CELLS\$annot),alpha=0.25,size=0.5) +
	geom_point(aes(CELLS\$xdim,CELLS\$ydim,color=CELLS\$annot),size=0.5) +
	geom_point(aes(CELLS\$centx,CELLS\$centy),color=\"gray20\",size=1.75) +
	geom_point(aes(CELLS\$centx,CELLS\$centy,color=CELLS\$annot),size=1) +
	guides(colour = guide_legend(override.aes = list(size=4))) +";
}

if ($color_mapping !~ /none/i) {
	print R "
	scale_colour_manual(values = c($color_mapping)) +";
}

if ($theme =~ /Clean/i) {
	print R "
	theme_bw() +
	theme(panel.border=element_blank(),
		  panel.grid=element_blank(),
		  axis.line=element_blank(),
		  axis.ticks=element_blank(),
		  legend.background=element_blank(),
		  legend.title=element_blank(),
		  panel.background=element_blank(),
		  axis.text=element_blank(),
		  axis.title.x=element_blank(),
		  axis.title.y=element_blank(),
		  plot.background=element_blank(),
		  plot.margin=unit(c(0,0,0,0),\"pt\"))\n";
} else {
	print R "
	theme_bw()\n";
}

print R "
ggsave(plot=PLT,filename=\"$opt{'O'}.aggregate.png\",width=5,height=4,dpi=900)
ggsave(plot=PLT,filename=\"$opt{'O'}.aggregate.pdf\",width=6,height=5)
";
close R;

system("$Rscript $opt{'O'}.aggregate_data.r");

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.aggregate_data $opt{'O'}.aggregate_data.r");
}
}

}
1;
