package sci_commands::matrix_filter;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_filter");

sub matrix_filter {

@ARGV = @_;

# Defaults
$colMin = 1;
$rowMin = 1;

getopts("O:C:R:c:r:A:a:B:vb:f:V:G:zM:", \%opt);

$die2 = "
scitools matrix-filter [options] [matrix OR sparseMatrix values/tfidf file]
   or    filter-matrix
   
Options:
   -O   [STR]   Output prefix (default is [input prefix].filt.matrix)
   -C   [INT]   Number of nonZero sites per column to retain (def = $colMin)
   -R   [INT]   Number of nonZero sites per row to retain (def = $rowMin)
   -c   [STR]   List of CellIDs to retain
   -r   [STR]   List of RowIDs to retain
   -B   [BED]   If chr_start_end features, will pull only those overlapping bed file
   -v           If -B is toggeled will exclude those peaks
   -f   [STR]   If -B has annotations, just consider those with the name -f
   -A   [STR]   Annotation file
   -a   [STR]   Comma separated list of annotations to include (requires -A)
   -V   [STR]   Values file to filter on (e.g. FRIP)
   -G   [MIN,MAX] Min and max values to include from values file (required if -V)
   -M   [STR]   Comma separated list of chromosomes to exclude (eg chrX,chr3)
                 Only works if row names are chr_start_end formatted.
   -z           Gzip output (def = no)
   -b   [STR]   Bedtools call (for -B filtering; def = $bedtools)

Note: -C and -R filters are applied after all other filtering.

";

if (!defined $ARGV[0]) {die $die2};
if (defined $opt{'C'}) {$colMin = $opt{'C'}};
if (defined $opt{'R'}) {$rowMin = $opt{'R'}};
if (defined $opt{'V'} && !defined $opt{'G'}) {die "ERROR: if a values fiel si provided (-V), a range must be provided as well (-R)\n"};
if ($ARGV[0] =~ /sparse/i) {
	$sparse = 1;
	$sparse_prefix = $ARGV[0]; $sparse_prefix =~ s/\.gz$//; $sparse_prefix =~ s/\.(values|tfidf)$//;
	$rowID_file = "na";
	if (-e "$sparse_prefix.rows.gz") {$rowID_file = "$sparse_prefix.rows.gz"};
	if (-e "$sparse_prefix.rows") {$rowID_file = "$sparse_prefix.rows"};
	if ($rowID_file eq "na") {die "ERROR: SParse matrix provided ($ARGV[0]) but the rows file ($sparse_prefix.rows(.gz)) cannot be found!\n"};
	$colID_file = "na";
	if (-e "$sparse_prefix.cols.gz") {$colID_file = "$sparse_prefix.cols.gz"};
	if (-e "$sparse_prefix.cols") {$colID_file = "$sparse_prefix.cols"};
	if ($colID_file eq "na") {die "ERROR: SParse matrix provided ($ARGV[0]) but the cols file ($sparse_prefix.cols(.gz)) cannot be found!\n"};
} else {
	$sparse = 0;
}
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
$opt{'O'} =~ s/\.gz$//i;
$opt{'O'} =~ s/\.(matrix|values|tfidf)$//i;
$opt{'O'} =~ s/\.sparseMatrix$//i;
$opt{'O'} .= ".filt";

if (defined $opt{'a'} && !defined $opt{'A'}) {die "\nMust provide an annotaiton file (-A) if specifying annotations to filter (-a)!\n$die2"};

if (defined $opt{'A'}) {read_annot($opt{'A'})};

if (defined $opt{'a'}) {
	@ANNOT_LIST = split(/,/, $opt{'a'});
	foreach $annot (@ANNOT_LIST) {
		$ANNOT_include{$annot} = 1;
	}
}
if (defined $opt{'M'}) {
	@CHR_EXCLUDE = split(/,/, $opt{'M'});
	foreach $chr (@CHR_EXCLUDE) {
		$CHR_exclude{$chr} = 1;
	}
}

if (defined $opt{'c'}) {
	open IN, "$opt{'c'}";
	while ($cellID = <IN>) {
		chomp $cellID;
		$CELLID_list_include{$cellID} = 1;
	} close IN;
}

if (defined $opt{'V'}) {
	($minV,$maxV) = split(/,/, $opt{'G'});
	open IN, "$opt{'V'}";
	while ($l = <IN>) {
		chomp $l;
		($cellID,$value) = split(/\t/, $l);
		if ($value >= $minV && $value <= $maxV) {
			$CELLID_values_include{$cellID} = 1;
		}
	} close IN;
}

if (defined $opt{'r'}) {
	open IN, "$opt{'r'}";
	while ($rowID = <IN>) {
		chomp $rowID;
		$ROWID_list_include{$rowID} = 1;
	} close IN;
}

if (defined $opt{'B'}) {
	if ($sparse == 0) {
		if ($ARGV[0] =~ /\.gz$/) {
			open MATRIX, "$zcat $ARGV[0] |"; $null = <MATRIX>;
		} else {
			open MATRIX, "$ARGV[0]"; $null = <MATRIX>;
		}
	} else {
		if ($rowID_file =~ /\.gz$/) {
			open MATRIX, "$zcat $rowID_file |";
		} else {
			open MATRIX, "$rowID_file";
		}
	}
	open BED, "| $bedtools sort -i - > $opt{'O'}.matrix_rows.bed";
	while ($l = <MATRIX>) {
		chomp $l; @P = split(/\t/, $l);
		($chr,$start,$end) = split(/_/, $P[0]);
		print BED "$chr\t$start\t$end\n";
	}
	close MATRIX; close BED;
	if (defined $opt{'f'}) {
		$opt_f_features = 0;
		open IN, "$opt{'B'}";
		open OUT, "| $bedtools sort -i - > $opt{'O'}.selected_features.bed";
		while ($l = <IN>) {
			chomp $l;
			@P = split(/\t/, $l);
			if ($P[3] eq "$opt{'f'}") {
				$opt_f_features++;
				print OUT "$P[0]\t$P[1]\t$P[2]\n";
			}
		} close IN; close OUT;
		if ($opt_f_features<1) {
			die "ERROR: Feature $opt{'f'} was specified for bed file $opt{'B'} and no features were found!\n";
		}
		$opt{'B'} = "$opt{'O'}.selected_features.bed";
	}
	if (!defined $opt{'v'}) {
		open INT, "$bedtools intersect -a $opt{'O'}.matrix_rows.bed -b $opt{'B'} -wa -wb |";
	} else {
		open INT, "$bedtools intersect -v -a $opt{'O'}.matrix_rows.bed -b $opt{'B'} -wa -wb |";
	}
	while ($l = <INT>) {
		chomp $l;
		@P = split(/\t/, $l);
		$rowID = "$P[0]_$P[1]_$P[2]";
		$ROWID_bed_include{$rowID} = 1;
	} close INT;
	system("rm -f $opt{'O'}.matrix_rows.bed");
	if (defined $opt{'f'}) {
		system("rm -f $opt{'O'}.selected_features.bed");
	}
} elsif (defined $opt{'M'}) {
	$opt{'B'} = "toggled via -M $opt{'M'}"; # trigger to use same filtering mechanism
	if ($sparse == 0) {
		if ($ARGV[0] =~ /\.gz$/) {
			open MATRIX, "$zcat $ARGV[0] |"; $null = <MATRIX>;
		} else {
			open MATRIX, "$ARGV[0]"; $null = <MATRIX>;
		}
	} else {
		if ($rowID_file =~ /\.gz$/) {
			open MATRIX, "$zcat $rowID_file |";
		} else {
			open MATRIX, "$rowID_file";
		}
	}
	while ($l = <MATRIX>) {
		chomp $l; @P = split(/\t/, $l);
		($chr,$start,$end) = split(/_/, $P[0]);
		$rowID = "$chr\_$start\_$end";
		if (!defined $CHR_exclude{$chr}) {
			$ROWID_bed_include{$rowID} = 1;
		}
	}
}

if ($sparse < 0.5) {
	if ($ARGV[0] =~ /\.gz$/) {
		open MATRIX, "$zcat $ARGV[0] |";
	} else {
		open MATRIX, "$ARGV[0]";
	}
	$h = <MATRIX>; chomp $h; @MATRIX_COLNAMES = split(/\t/, $h);
	$matrix_colNum = @MATRIX_COLNAMES;
	for ($cellNum = 0; $cellNum < @MATRIX_COLNAMES; $cellNum++) {
		$cellID = $MATRIX_COLNAMES[$cellNum];
		$COLNAME_nonzero{$cellID} = 0;
	}
	$matrix_rowNum = 0;
} else {
	if ($rowID_file =~ /\.gz$/) {
		open ROWS, "$zcat $rowID_file |";
	} else {
		open ROWS, "$rowID_file";
	}
	$rowNum = 0;
	while ($rowID = <ROWS>) {
		chomp $rowID;
		$rowNum++;
		$ROWnum_ROWname{$rowNum} = $rowID;
	} close ROWS;
	$matrix_rowNum = $rowNum;
	
	if ($colID_file =~ /\.gz$/) {
		open COLS, "$zcat $colID_file |";
	} else {
		open COLS, "$colID_file";
	}
	$colNum = 0;
	while ($colID = <COLS>) {
		chomp $colID;
		$colNum++;
		$COLnum_COLname{$colNum} = $colID;
	} close COLS;
	$matrix_colNum = $colNum;
	
	if ($ARGV[0] =~ /\.gz$/) {
		open MATRIX, "$zcat $ARGV[0] |";
	} else {
		open MATRIX, "$ARGV[0]";
	}
}

while ($l = <MATRIX>) {
	chomp $l;
	if ($sparse < 0.5) {
		$matrix_rowNum++;
		@P = split(/\t/, $l);
		$rowID = shift(@P);
		$ROWNAME_nonzero{$rowID} = 0;
		if ((!defined $opt{'r'} && !defined $opt{'B'}) ||
			(defined $opt{'r'} && defined $ROWID_list_include{$rowID}) ||
			(defined $opt{'B'} && defined $ROWID_bed_include{$rowID})) {
			for ($cellNum = 0; $cellNum < @P; $cellNum++) {
				$cellID = $MATRIX_COLNAMES[$cellNum];
				if ((!defined $opt{'c'} || defined $CELLID_list_include{$cellID}) &&
					(!defined $opt{'a'} || defined $ANNOT_include{$CELLID_annot{$cellID}}) &&
					(!defined $opt{'V'} || defined $CELLID_values_include{$cellID})) {
					if (abs($P[$cellNum]) > 0) {
						$COLNAME_nonzero{$cellID}++;
						$ROWNAME_nonzero{$rowID}++;
					}
				}
			}
		}
	} else {
		($rowNum,$colNum,$value) = split(/\s/, $l);
		$rowID = $ROWnum_ROWname{$rowNum};
		$cellID = $COLnum_COLname{$colNum};
		if (!defined $COLNAME_nonzero{$cellID}) {$COLNAME_nonzero{$cellID} = 0};
		if (!defined $ROWNAME_nonzero{$rowID}) {$ROWNAME_nonzero{$rowID} = 0};

		if ((!defined $opt{'r'} && !defined $opt{'B'}) ||
		    (defined $opt{'r'} && defined $ROWID_list_include{$rowID}) ||
		    (defined $opt{'B'} && defined $ROWID_bed_include{$rowID})) {
		    if ((!defined $opt{'c'} || defined $CELLID_list_include{$cellID}) &&
			(!defined $opt{'a'} || defined $ANNOT_include{$CELLID_annot{$cellID}}) &&
			(!defined $opt{'V'} || defined $CELLID_values_include{$cellID})) {
			        if (abs($value)>0) {
					$COLNAME_nonzero{$cellID}++;
					$ROWNAME_nonzero{$rowID}++;
				}
		    }
		}
	}
} close MATRIX;

if ($sparse < 0.5) {
	if (defined $opt{'z'}) {
		open OUT, "| $gzip > $opt{'O'}.matrix.gz";
	} else {
		open OUT, ">$opt{'O'}.matrix";
	}
	$out_header = ""; $included_cells = 0; $included_rows = 0;
	for ($cellNum = 0; $cellNum < @MATRIX_COLNAMES; $cellNum++) {
		$cellID = $MATRIX_COLNAMES[$cellNum];
		if ($COLNAME_nonzero{$cellID}>=$colMin) {
			$out_header .= "$cellID\t";
			$included_cells++;
		}
	} $out_header =~ s/\t$//;
	print OUT "$out_header\n";
	if ($ARGV[0] =~ /\.gz$/) {
		open MATRIX, "$zcat $ARGV[0] |"; $null = <MATRIX>;
	} else {
		open MATRIX, "$ARGV[0]"; $null = <MATRIX>;
	}
	while ($l = <MATRIX>) {
		chomp $l;
		@P = split(/\t/, $l);
		$rowID = shift(@P);
		if (($ROWNAME_nonzero{$rowID}>=$rowMin) &&
			((!defined $opt{'r'} && !defined $opt{'B'}) ||
			 (defined $opt{'r'} && defined $ROWID_list_include{$rowID}) ||
			 (defined $opt{'B'} && defined $ROWID_bed_include{$rowID}))) {
				$included_rows++;
				print OUT "$rowID";
				for ($cellNum = 0; $cellNum < @P; $cellNum++) {
					$cellID = $MATRIX_COLNAMES[$cellNum];
					if ($COLNAME_nonzero{$cellID}>=$colMin) {
						print OUT "\t$P[$cellNum]";
					}
				}
				print OUT "\n";
		}
	} close MATRIX; close OUT;
} else {
	if (defined $opt{'z'}) {
		open OUT, "| $gzip > $opt{'O'}.sparseMatrix.rows.gz";
	} else {
		open OUT, ">$opt{'O'}.sparseMatrix.rows";
	}
	if ($rowID_file =~ /\.gz$/) {
		open ROWS, "$zcat $rowID_file |";
	} else {
		open ROWS, "$rowID_file";
	}
	%oldRowNum_newRowNum = ();
	$newRowNum = 1; $rowNum = 1;
	while ($rowID = <ROWS>) {
		chomp $rowID;
		if (($ROWNAME_nonzero{$rowID}>=$rowMin) &&
			((!defined $opt{'r'} && !defined $opt{'B'}) ||
			 (defined $opt{'r'} && defined $ROWID_list_include{$rowID}) ||
			 (defined $opt{'B'} && defined $ROWID_bed_include{$rowID}))) {
				$included_rows++;
				$oldRowNum_newRowNum{$rowNum} = $newRowNum;
				$newRowNum++;
				print OUT "$rowID\n";
		}
		$rowNum++;
	}
	close OUT; close ROWS;
	
	if (defined $opt{'z'}) {
		open OUT, "| $gzip > $opt{'O'}.sparseMatrix.cols.gz";
	} else {
		open OUT, ">$opt{'O'}.sparseMatrix.cols";
	}
	if ($colID_file =~ /\.gz$/) {
		open COLS, "$zcat $colID_file |";
	} else {
		open COLS, "$colID_file";
	}
	%oldColNum_newColNum = ();
	$colNum = 1; $newColNum = 1;
	while ($cellID = <COLS>) {
		chomp $cellID;
		if ($COLNAME_nonzero{$cellID}>=$colMin) {
			$included_cells++;
			$oldColNum_newColNum{$colNum} = $newColNum;
			$newColNum++;
			print OUT "$cellID\n";
		}
		$colNum++;
	}
	close OUT; close COLS;
	
	if (defined $opt{'z'}) {
		open OUT, "| $gzip > $opt{'O'}.sparseMatrix.values.gz";
	} else {
		open OUT, ">$opt{'O'}.sparseMatrix.values";
	}
	if ($ARGV[0] =~ /\.gz$/) {
		open MATRIX, "$zcat $ARGV[0] |";
	} else {
		open MATRIX, "$ARGV[0]";
	}
	while ($l = <MATRIX>) {
		chomp $l;
		($rowNum,$colNum,$value) = split(/\s/, $l);
		if (defined $oldRowNum_newRowNum{$rowNum} && defined $oldColNum_newColNum{$colNum}) {
			print OUT "$oldRowNum_newRowNum{$rowNum}\t$oldColNum_newColNum{$colNum}\t$value\n";
		}
	} close MATRIX; close OUT;
}

open LOG, ">$opt{'O'}.log";
$ts = localtime(time);
print LOG "$ts scitools atac-filter
Matrix file = $ARGV[0]
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

}
1;
