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

getopts("O:C:R:c:r:A:a:", \%opt);

$die2 = "
scitools matrix-filter [options] [counts matrix]
   or    filter-matrix
   
Options:
   -O   [STR]   Output prefix (default is [input].filt_[colMin]_[rowMin].matrix)
   -C   [INT]   Number of nonZero sites per column to retain (def = $colMin)
   -R   [INT]   Number of nonZero sites per row to retain (def = $rowMin)
   -c   [STR]   List of CellIDs to retain
   -r   [STR]   List of RowIDs to retain
   -A   [STR]   Annotation file
   -a   [STR]   Comma separated list of annotations to include (requires -A)

Note: -C and -R filters are applied after all other filtering.

";

if (!defined $ARGV[0]) {die $die2};
if (defined $opt{'C'}) {$colMin = $opt{'C'}};
if (defined $opt{'R'}) {$rowMin = $opt{'R'}};
if (!defined $opt{'O'}) {
	$opt{'O'} = $ARGV[0];
	$opt{'O'} =~ s/\.matrix$//;
	$opt{'O'} .= ".filt_$colMin\_$rowMin";
}
if (defined $opt{'a'} && !defined $opt{'A'}) {die "\nMust provide an annotaiton file (-A) if specifying annotations to filter (-a)!\n$die2"};

if (defined $opt{'A'}) {read_annot($opt{'A'})};

if (defined $opt{'a'}) {
	@ANNOT_LIST = split(/,/, $opt{'a'});
	foreach $annot (@ANNOT_LIST) {
		$ANNOT_include{$annot} = 1;
	}
}

if (defined $opt{'c'}) {
	open IN, "$opt{'c'}";
	while ($cellID = <IN>) {
		chomp $cellID;
		$CELLID_list_include{$cellID} = 1;
	} close IN;
}

if (defined $opt{'r'}) {
	open IN, "$opt{'r'}";
	while ($rowID = <IN>) {
		chomp $rowID;
		$ROWID_list_include{$rowID} = 1;
	} close IN;
}

open MATRIX, "$ARGV[0]";
$h = <MATRIX>; chomp $h; @MATRIX_COLNAMES = split(/\t/, $h);
$matrix_colNum = @MATRIX_COLNAMES;
for ($cellNum = 0; $cellNum < @MATRIX_COLNAMES; $cellNum++) {
	$cellID = $MATRIX_COLNAMES[$cellNum];
	$COLNAME_nonzero{$cellID} = 0;
}
$matrix_rowNum = 0;
while ($l = <MATRIX>) {
	$matrix_rowNum++;
	chomp $l;
	@P = split(/\t/, $l);
	$rowID = shift(@P);
	$ROWNAME_nonzero{$rowID} = 0;
	if (!defined $opt{'r'} || defined $ROWID_list_include{$rowID}) {
		for ($cellNum = 0; $cellNum < @P; $cellNum++) {
			$cellID = $MATRIX_COLNAMES[$cellNum];
			if ((!defined $opt{'c'} || defined $CELLID_list_include{$cellID}) &&
				(!defined $opt{'a'} || defined $ANNOT_include{$CELLID_annot{$cellID}})) {
				if (abs($P[$cellNum]) > 0) {
					$COLNAME_nonzero{$cellID}++;
					$ROWNAME_nonzero{$rowID}++;
				}
			}
		}
	}
} close MATRIX;

open OUT, ">$opt{'O'}.matrix";

$out_header = ""; $included_cells = 0; $included_rows = 0;
for ($cellNum = 0; $cellNum < @MATRIX_COLNAMES; $cellNum++) {
	$cellID = $MATRIX_COLNAMES[$cellNum];
	if ($COLNAME_nonzero{$cellID}>=$colMin) {
		$out_header .= "$cellID\t";
		$included_cells++;
	}
} $out_header =~ s/\t$//;
print OUT "$out_header\n";

open MATRIX, "$ARGV[0]"; $null = <MATRIX>;
while ($l = <MATRIX>) {
	chomp $l;
	@P = split(/\t/, $l);
	$rowID = shift(@P);
	if ($ROWNAME_nonzero{$rowID}>=$rowMin) {
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
