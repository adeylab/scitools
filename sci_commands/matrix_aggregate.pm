package sci_commands::matrix_aggregate;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_aggregate");

sub matrix_aggregate {

@ARGV = @_;
use Getopt::Std; %opt = ();
getopts("O:r:c:a:bm", \%opt);

$die2 = "
scitools matrix-aggregate [options] [counts matrix] [annotation file]
   or    aggregate-matrix

Note: if a non counts matrix is provided, it will still sum the values.
Support for comma-separated annotations (i.e. cell aggregation with
oversampling) has been added.

Options:
   -O   [STR]   Output prefix (default is [input annot].matrix)
   -c   [STR]   File with list of CellIDs to retain
   -r   [STR]   File with list of RowIDs to retain
   -a   [STR]   Comma separated list of annotations to include
   -b           Treat each cell as a binary value (resulting matrix is
                  the number of cells at the peak. def = total counts)
   -m           Instead aggregate by the annotations and print a mean
                  and standard deviation for the aggregates
                  (overrides -b)
   
";

#name output 
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[1]; $opt{'O'} =~ s/\.matrix$//};


#read in annotation, scitools approach
if (!defined $ARGV[1]) {die $die2}
else {read_annot($ARGV[1])};

#define on your own which annot to include, scitools approach
if (defined $opt{'a'}) {
	@ANNOT_LIST = split(/,/, $opt{'a'});
	foreach $annot (@ANNOT_LIST) {
		$ANNOT_include{$annot} = 1;
	}
}

#define cells to keep
if (defined $opt{'c'}) {
	open IN, "$opt{'c'}";
	while ($cellID = <IN>) {
		chomp $cellID;
		$CELLID_list_include{$cellID} = 1;
	} close IN;
}

#define feature to keep
if (defined $opt{'r'}) {
	open IN, "$opt{'r'}";
	while ($rowID = <IN>) {
		chomp $rowID;
		$ROWID_list_include{$rowID} = 1;
	} close IN;
}

#read in matrix scitools standard
read_matrix($ARGV[0]);

#define your matrix based on which cells and features are allowed to be kept
#same style as CELL_FEATURE_VALUE
%ANNOT_FEATURE_aggregate_value = ();
#annot that we aggreagate by
@MATRIX_aggregate_COLNAMES = ();

if (defined $opt{'m'}) {
	%ANNOT_FEATURE_list = ();
	%ANNOT_FEATURE_aggregate_stdev = ();
}

#then the cells and aggregate signal
foreach $cellID (keys %CELLID_FEATURE_value) {
	#define cells to include based -c and -a options
	#start out with all Cells being included
	$includeCell = 1;
	#remove cells that are not in include list provided
	if (defined $opt{'c'} && !defined $CELLID_list_include{$cellID}) {$includeCell = 0};
	
	# handle multi-annot
	@ANNOTS = split(/,/, $CELLID_annot{$cellID});
	foreach $annot (@ANNOTS) {
		#remove cells that are not in annot list provided
		if (defined $opt{'a'} && !defined $ANNOT_include{$annot}) {$includeCell = 0};
		if ($includeCell>0) {
			foreach $rowID (keys %{$CELLID_FEATURE_value{$cellID}}) {
				$includeRow = 1;
				#define which rows to include based on -r option
				if (defined $opt{'r'} && !defined $ROWID_list_include{$rowID}) {$includeRow = 0};
				if ($includeRow>0 ) {
					if (defined $opt{'m'}) {
						push @{$ANNOT_FEATURE_list{$annot}{$rowID}}, $CELLID_FEATURE_value{$cellID}{$rowID};
					} elsif (defined $opt{'b'}) {
						if ($CELLID_FEATURE_value{$cellID}{$rowID}>0) {
							if (defined $ANNOT_FEATURE_aggregate_value{$annot}{$rowID}) {
								$ANNOT_FEATURE_aggregate_value{$annot}{$rowID}++;
							}
						} else {
							$ANNOT_FEATURE_aggregate_value{$annot}{$rowID}=1;
						}
					} else {
						if (defined $ANNOT_FEATURE_aggregate_value{$annot}{$rowID}) {
							$ANNOT_FEATURE_aggregate_value{$annot}{$rowID}+=$CELLID_FEATURE_value{$cellID}{$rowID};
						} else {
							$ANNOT_FEATURE_aggregate_value{$annot}{$rowID}=$CELLID_FEATURE_value{$cellID}{$rowID};
						}
					}
				}  
			}
			$included_cells++;
			$CELLID_include{$cellID} = 1;
			#define annotations to include
			$ANNOT_include{$annot} = 1;  
		}
	}
}

if (defined $opt{'m'}) {
	foreach $annot (keys %ANNOT_FEATURE_list) {
		foreach $rowID (keys %{$ANNOT_FEATURE_list{$annot}}) {
			$element_count = 0; $element_sum = 0;
			foreach $value (@{$ANNOT_FEATURE_list{$annot}{$rowID}}) {
				$element_count++; $element_sum += $value;
			}
			if ($element_count>0) {
				$element_mean = $element_sum/$element_count;
				$ANNOT_FEATURE_aggregate_value{$annot}{$rowID} = $element_mean;
				$element_stdev_sum = 0;
				foreach $value (@{$ANNOT_FEATURE_list{$annot}{$rowID}}) {
					$element_stdev_sum += ($element_mean-$value)**2;
				}
				$ANNOT_FEATURE_aggregate_stdev{$annot}{$rowID} = sqrt($element_stdev_sum/$element_count);
			} else {
				$ANNOT_FEATURE_aggregate_value{$annot}{$rowID} = 0;
				$ANNOT_FEATURE_aggregate_stdev{$annot}{$rowID} = 0;
			}
		}
	}
	%ANNOT_FEATURE_list = ();
}

#create a matrix of included annot
foreach $annot (keys %ANNOT_include) {
	push(@MATRIX_aggregate_COLNAMES,$annot);
}

#then look which rows are allowed to be kept
foreach $rowID (keys %MATRIX_feature_nonZero) {
	if (!defined $opt{'r'} || defined $ROWID_list_include{$rowID}) {
		$ROWID_count{$rowID} = 0;
		foreach $cellID (keys %CELLID_include) {
			if ($CELLID_FEATURE_value{$cellID}{$rowID}>0) {
				$ROWID_count{$rowID}++;
			}
		}
		$ROWID_include{$rowID} = 1;
		$included_rows++;
	}       
}

#do write out depending on b option
if (defined $opt{'m'}) {
	open OUT, ">$opt{'O'}.mean.aggregate.matrix";
	open SDV, ">$opt{'O'}.stdev.aggregate.matrix";
	open LOG, ">$opt{'O'}.mean.aggregate.matrix.log"
} elsif (defined $opt{'b'}) {
	open OUT, ">$opt{'O'}.binary.aggregate.matrix";
	open LOG, ">$opt{'O'}.binary.aggregate.matrix.log";
} else {
	open OUT, ">$opt{'O'}.aggregate.matrix";
	open LOG, ">$opt{'O'}.aggregate.matrix.log";
}

#first define header
$out_header = join("\t",@MATRIX_aggregate_COLNAMES);

print OUT "$out_header\n";
if (defined $opt{'m'}) {print SDV "$out_header\n"};

for ($rowNum = 0; $rowNum < @MATRIX_ROWNAMES; $rowNum++) {
	$rowID = $MATRIX_ROWNAMES[$rowNum];
	if (defined $ROWID_include{$rowID}) {
		print OUT "$rowID"; if (defined $opt{'m'}) {print SDV "$rowID"};
		for ($aggrNum = 0; $aggrNum < @MATRIX_aggregate_COLNAMES; $aggrNum++) {
			print OUT "\t$ANNOT_FEATURE_aggregate_value{$MATRIX_aggregate_COLNAMES[$aggrNum]}{$rowID}";
			if (defined $opt{'m'}) {print SDV "\t$ANNOT_FEATURE_aggregate_stdev{$MATRIX_aggregate_COLNAMES[$aggrNum]}{$rowID}"};
		} print OUT "\n"; if (defined $opt{'m'}) {print SDV "\n"};
	}
} close OUT;
if (defined $opt{'m'}) {close SDV};

$ts = localtime(time);
print LOG "$ts scitools annot-aggr
Matrix file = $ARGV[1]
Annotation file = $ARGV[0]

Options:
";
foreach $option (keys %opt) {
	print LOG "   $option   $opt{$option}\n";
}
print LOG "Total cells in input: $matrix_colNum
Total cells retained in output: $included_cells
Total rows in input: $matrix_rowNum
Total rows retained in output: $included_rows\n";
close LOG;

exit;

}
1;
