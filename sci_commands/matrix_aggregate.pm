package sci_commands::matrix_aggregate;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_aggregate");

sub matrix_aggregate {

@ARGV = @_;
use Getopt::Std; %opt = ();
getopts("O:a:o:tk:zL:W:", \%opt);

# defaults
$operations = "add";
$constant = "rownum";

$die2 = "
scitools matrix-aggregate [options] [matrix or sparseMatrix values file] [annotation file]
   or    aggregate-matrix

Options:
   -O   [STR]   Output prefix (default is [input annot].matrix)
   -a   [STR]   Comma separated list of annotations to include, will also order output
   -o   [STR]   Operation(s), comma separated list of the follwing (def = $operations)
                  Each will be a seprate 
                  add/a       adds values of merged cells
                  mean/m      mean of merged cells
                  stdev/s     standard deviation
                  median/d    median of cells
                  count/c     number of nonzero cells at the site
                  fraction/f  fraction of nonzero cells
   -t           Factor in term-frequency prior to performing operation
                  (value of site / sum of all sites for the cell) * constant
   -k   [STR]   Constant to multiply for term-frequency can be a number or 'rownum'
                  def = $constant
   -L   [STR]  Columns file for sparseMatrix (will try to auto-detect)
   -W   [STR]  Rows file for sparseMatrix (will try to auto-detect)
   -z           Gzip output
   
";

if (!defined $ARGV[1]) {die $die2};

if (!defined $opt{'O'}) {
	$prefix = $ARGV[0];
	$prefix =~ s/\.gz$//;
	$prefix =~ s/\.matrix$//;
	$prefix =~ s/\.sparseMatrix$//;
	$afile = $ARGV[1];
	$afile =~ s/\.annot$//;
	$opt{'O'} = $prefix.".".$afile;
}

# get sparseMatrix components
if ($ARGV[0] =~ /sparseMatrix/i) {  die "ERROR: SParsematrix support in progress!\n"; ######## SPARSE DIE ########
	$sparse = 1;
	if (defined $opt{'L'}) {
		$col_file = $opt{'L'};
	} else {
		if (-e "$prefix.sparseMatrix.cols") {
			$col_file = "$prefix.sparseMatrix.cols";
		} elsif (-e "$prefix.sparseMatrix.cols.gz") {
			$col_file = "$prefix.sparseMatrix.cols.gz";
		} else {
			die "ERROR: Cannot detect cols file (e.g. $prefix.sparseMatrix.cols), please provide as -L\n";
		}
	}
	if (defined $opt{'W'}) {
		$row_file = $opt{'W'};
	} else {
		if (-e "$prefix.sparseMatrix.rows") {
			$row_file = "$prefix.sparseMatrix.rows";
		} elsif (-e "$prefix.sparseMatrix.rows.gz") {
			$row_file = "$prefix.sparseMatrix.rows.gz";
		} else {
			die "ERROR: Cannot detect rows file (e.g. $prefix.sparseMatrix.rows), please provide as -W\n";
		}
	}
	open IN, "$col_file"; @COLS = (); push @COLS, "cols";
	while ($cellID = <IN>) {
		chomp $cellID;
		push @COLS, $cellID;
	} close IN;
	open IN, "$row_file"; @ROWS = (); push @ROWS, "rows";
	while ($rowID = <IN>) {
		chomp $rowID;
		push @ROWS, $rowID;
	} close IN;
} else {$sparse = 0};

# load operations
if (defined $opt{'o'}) {$operations = $opt{'o'}};
@OPS = split(/,/, $operations);
foreach $operation (@OPS) {
	if ($operation =~ /(^a$|sum|add)/) {$OPS{'a'} = 1}
	elsif ($operation =~ /(^m$|mean|average)/) {$OPS{'m'} = 1}
	elsif ($operation =~ /(^d$|median)/) {$OPS{'d'} = 1}
	elsif ($operation =~ /(^s$|stdev)/) {$OPS{'s'} = 1}
	elsif ($operation =~ /(^c$|count)/) {$OPS{'c'} = 1}
	elsif ($operation =~ /(^f$|fraction)/) {$OPS{'f'} = 1}
	else {die "ERROR: Cannot determine which operation '$operation' represents!\n"};
}

# get annot include
if (defined $opt{'a'}) {
	@ANNOT_INCLUDE = split(/,/, $opt{'a'});
	foreach $annot (@ANNOT_INCLUDE) {
		$ANNOT_include{$annot} = 1;
	}
} else {
	@ANNOT_INCLUDE = ();
}

# read in annotations for cellIDs
open IN, "$ARGV[1]";
while ($l = <IN>) {
	chomp $l;
	($cellID,$annot_list) = split(/\t/, $l);
	@{$CELLID_annots{$cellID}} = split(/,/, $annot_list);
	foreach $annot (@{$CELLID_annots{$cellID}}) {
		if (!defined $opt{'a'} || defined $ANNOT_include{$annot}) {
			$ANNOT_CellIDs{$annot}{$cellID} = 1;
			if (!defined $ANNOT_list{$annot} && !defined $opt{'a'}) {
				push @ANNOT_INCLUDE, $annot;
				$ANNOT_list{$annot} = 1;
			}
		}
	}
} close IN;

# cellID sums for term-frequency
if (defined $opt{'t'}) {
	$rowCount = 0;
	if ($sparse > 0.5) {
		if ($ARGV[0] =~ /\.gz/) {
			open IN, "$zcat $ARGV[0] |";
		} else {
			open IN, "$ARGV[0]";
		}
		while ($l = <IN>) {
			chomp $l;
			($rowNum,$cellNum,$value) = split(/\t/, $l);
			$CELLID_sum{$COLS[$cellNum]} += $value;
			if ($rowNum>$rowCount){$rowCount = $rowNum};
		} close IN;
	} else {
		if ($ARGV[0] =~ /\.gz/) {
			open IN, "$zcat $ARGV[0] |";
		} else {
			open IN, "$ARGV[0]";
		}
		$h = <IN>; chomp $h; @COLS = split(/\t/, $h);
		while ($l = <IN>) {
			$rowCount++;
			@P = split(/\t/, $l);
			$rowID = shift(@P);
			for ($i = 0; $i < @COLS; $i++) {
				$CELLID_sum{$COLS[$i]} += $P[$i];
			}
		} close IN;
	}
	if (defined $opt{'k'}) {
		if ($opt{'k'} =~ /row/) {
			$constant = $rowCount;
		} else {
			$constant = $opt{'k'};
		}
	} else {
		$constant = $rowCount;
	}
}

# perform operations row-by-row
if ($sparse>0.5) {
	# open OUTs
	if (defined $OPS{'a'}) {if (defined $opt{'z'}) {open ADD, "| $gzip > $opt{'O'}.sparseMatrix.sum.gz"} else {open ADD, ">$opt{'O'}.sparseMatrix.sum"}};
	if (defined $OPS{'m'}) {if (defined $opt{'z'}) {open MEAN, "| $gzip > $opt{'O'}.sparseMatrix.mean.gz"} else {open MEAN, ">$opt{'O'}.sparseMatrix.matrix"}};
	if (defined $OPS{'d'}) {if (defined $opt{'z'}) {open MEDIAN, "| $gzip > $opt{'O'}.sparseMatrix.median.gz"} else {open MEDIAN, ">$opt{'O'}.sparseMatrix.median"}};
	if (defined $OPS{'c'}) {if (defined $opt{'z'}) {open COUNT, "| $gzip > $opt{'O'}.sparseMatrix.count.gz"} else {open COUNT, ">$opt{'O'}.sparseMatrix.count"}};
	if (defined $OPS{'s'}) {if (defined $opt{'z'}) {open STDEV, "| $gzip > $opt{'O'}.sparseMatrix.stdev.gz"} else {open STDEV, ">$opt{'O'}.sparseMatrix.stdev"}};
	if (defined $opt{'z'}) {system("cp $row_file $opt{'O'}.sparseMatrix.rows")} else {system("cp $row_file $opt{'O'}.sparseMatrix.rows.gz")}
	
	###### PICK UP SPARSE HERE ######
	
} else {
	# make column header
	$col_names = join("\t", @ANNOT_INCLUDE);
	# open OUTs
	if (defined $OPS{'a'}) {if (defined $opt{'z'}) {open ADD, "| $gzip > $opt{'O'}.sum.matrix.gz"} else {open ADD, ">$opt{'O'}.sum.matrix"}; print ADD "$col_names\n"};
	if (defined $OPS{'m'}) {if (defined $opt{'z'}) {open MEAN, "| $gzip > $opt{'O'}.mean.matrix.gz"} else {open MEAN, ">$opt{'O'}.mean.matrix"}; print MEAN "$col_names\n"};
	if (defined $OPS{'d'}) {if (defined $opt{'z'}) {open MEDIAN, "| $gzip > $opt{'O'}.median.matrix.gz"} else {open MEDIAN, ">$opt{'O'}.median.matrix"}; print MEDIAN "$col_names\n"};
	if (defined $OPS{'c'}) {if (defined $opt{'z'}) {open COUNT, "| $gzip > $opt{'O'}.count.matrix.gz"} else {open COUNT, ">$opt{'O'}.count.matrix"}; print COUNT "$col_names\n"};
	if (defined $OPS{'s'}) {if (defined $opt{'z'}) {open STDEV, "| $gzip > $opt{'O'}.stdev.matrix.gz"} else {open STDEV, ">$opt{'O'}.stdev.matrix"}; print STDEV "$col_names\n"};
	if (defined $OPS{'f'}) {if (defined $opt{'z'}) {open FRAC, "| $gzip > $opt{'O'}.fraction.matrix.gz"} else {open FRAC, ">$opt{'O'}.fraction.matrix"}; print FRAC "$col_names\n"};

	if ($ARGV[0] =~ /\.gz/) {
		open IN, "$zcat $ARGV[0] |";
	} else {
		open IN, "$ARGV[0]";
	}
	$h = <IN>; chomp $h; @COLS = split(/\t/, $h);
	
	# get cell counts for each annot
	for ($i = 0; $i < @COLS; $i++) {
		foreach $annot (keys %ANNOT_CellIDs) {
			if (defined $ANNOT_CellIDs{$annot}{$COLS[$i]}) {
				$ANNOT_cellCount{$annot}++;
			}
		}
	}
	
	while ($l = <IN>) {
		@P = split(/\t/, $l);
		$rowID = shift(@P);
		
		# load in values
		%ANNOT_values = ();
		for ($i = 0; $i < @P; $i++) {
			foreach $annot (keys %ANNOT_CellIDs) {
				if (defined $ANNOT_CellIDs{$annot}{$COLS[$i]}) {
					if (!defined $opt{'t'}) {
						push @{$ANNOT_values{$annot}}, $P[$i];
					} else {
						$value = ($P[$i]/$CELLID_sum{$COLS[$i]})*$constant;
						push @{$ANNOT_values{$annot}}, $value;
					}
				}
			}
		}
		
		if (defined $OPS{'c'} || defined $OPS{'f'}) {
			if (defined $OPS{'c'}) {print COUNT "$rowID"};
			if (defined $OPS{'f'}) {print FRAC "$rowID"};
			for ($i = 0; $i < @ANNOT_INCLUDE; $i++) {
				$annot = $ANNOT_INCLUDE[$i]; $nonZero = 0;
				for ($j = 0; $j < @{$ANNOT_values{$annot}}; $j++) {
					if ($ANNOT_values{$annot}[$j] != 0) {$nonZero++};
				}
				if (defined $OPS{'c'}) {print COUNT "\t$nonZero"};
				if (defined $OPS{'f'}) {
					$frac = ($nonZero/$ANNOT_cellCount{$annot});
					print FRAC "\t$frac";
				}
			}
			if (defined $OPS{'c'}) {print COUNT "\n"};
			if (defined $OPS{'f'}) {print FRAC "\n"};
		}
		if (defined $OPS{'d'}) {
			print MEDIAN "$rowID";
			for ($i = 0; $i < @ANNOT_INCLUDE; $i++) {
				$annot = $ANNOT_INCLUDE[$i];
				@SORTED = sort {$a<=>$b} @{$ANNOT_values{$annot}};
				$median = $SORTED[int((@SORTED-1)/2)];
				print MEDIAN "\t$median";
			}
			print MEDIAN "\n";
		}
		if (defined $OPS{'a'} || defined $OPS{'m'} || defined $OPS{'s'}) {
			if (defined $OPS{'a'}) {print ADD "$rowID"; %ANNOT_sum = ()};
			if (defined $OPS{'m'}) {print MEAN "$rowID"; %ANNOT_mean = ()};
			if (defined $OPS{'s'}) {print STDEV "$rowID"; %ANNOT_stdev_sum = ()};
			for ($i = 0; $i < @ANNOT_INCLUDE; $i++) {
				$annot = $ANNOT_INCLUDE[$i];
				$ANNOT_sum{$annot} = 0;
				for ($j = 0; $j < @{$ANNOT_values{$annot}}; $j++) {
					$ANNOT_sum{$annot} += $ANNOT_values{$annot}[$j];
				}
				if (defined $OPS{'a'}) {print ADD "\t$ANNOT_sum{$annot}"};
				if (defined $OPS{'m'} || defined $OPS{'s'}) {
					$ANNOT_mean{$annot} = $ANNOT_sum{$annot}/$ANNOT_cellCount{$annot};
				}
				if (defined $OPS{'m'}) {print MEAN "\t$ANNOT_mean{$annot}"};
				if (defined $OPS{'s'}) {
					for ($j = 0; $j < @{$ANNOT_values{$annot}}; $j++) {
						$ANNOT_stdev_sum{$annot} += ($ANNOT_values{$annot}[$j] - $ANNOT_mean{$annot})**2;
					}
					$stdev = sqrt($ANNOT_stdev_sum{$annot}/$ANNOT_cellCount{$annot});
					print STDEV "\t$stdev";
				}
			}
			if (defined $OPS{'a'}) {print ADD "\n"};
			if (defined $OPS{'m'}) {print MEAN "\n"};
			if (defined $OPS{'s'}) {print STDEV "\n"};
		}
		
	} close IN;
}

exit;
}
1;
