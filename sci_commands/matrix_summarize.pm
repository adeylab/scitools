package sci_commands::matrix_summarize;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_summarize");

sub matrix_summarize {

@ARGV = @_;

getopts("O:R:XA:a:", \%opt);

$die2 = "
scitools matrix-summarize [options] [matrix of any kind]
   or    summarize-matrix
   
Generates summaries on the matrix file and plots some of the properties.

Options:
   -O   [STR]   Output prefix (default is [input_matrix].suffix)
   -R   [STR]   Rscript call (def = $Rscript)
   -X           Retain intermediate files (def = delete)

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};

# initialize stats
$ncol = 0; $nrow = 0;
$maxVal = "na"; $minVal = "na";
$sum = 0;
$nonzero = 0; $nrow_nonzero = 0; $ncol_nonzero = 0;
$mean = 0; $stdev = 0;

open STATS, ">$opt{'O'}.stats";
open ROW_DATA, ">$opt{'O'}.row.stats";
print ROW_DATA "mean\tstdev\tmin\tmax\tsum\tnonzero\n";
open COL_DATA, ">$opt{'O'}.col.stats";
print COL_DATA "mean\tstdev\tmin\tmax\tsum\tnonzero\n";
open COL_DIST, ">$opt{'O'}.col.vals";

open IN, "$ARGV[0]";
$h = <IN>; chomp $h; @H = split(/\t/, $h);
$ncol = @H;

for ($i = 0; $i < @H; $i++) {
	@{$COL_values[$i]} = ();
}

while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	$rowID = shift(@P);
	$row_sum = 0; $row_mean = 0; $row_stdev_sum = 0; $row_stdev = 0; $row_nonzero = 0; $row_min = $P[0]; $row_max = $P[0];
	for ($i = 0; $i < @H; $i++) {
		$sum+=$P[0];
		if ($P[$i] != 0) {$nonzero++; $row_nonzero++};
		if ($minVal eq "na" || $P[$i] < $minVal) {$minVal = $P[$i]};
		if ($maxVal eq "na" || $P[$i] > $maxVal) {$maxVal = $P[$i]};
		$row_sum += $P[$i];
		if ($P[$i] < $row_min) {$row_min = $P[$i]};
		if ($P[$i] > $row_max) {$row_max = $P[$i]};
		push @{$COL_values[$i]}, $P[$i];
	}
	$row_mean = $row_sum/$ncol;
	for ($i = 0; $i < @H; $i++) {
		$row_stdev_sum += ($P[$i] - $row_mean)**2;
	}
	$row_stdev = sqrt($row_stdev_sum/$ncol);
	print ROW_DATA "$rowID\t$row_mean\t$row_stdev\t$row_min\t$row_max\t$row_sum\t$row_nonzero\n";
	if ($row_nonzero>0) {$nrow_nonzero++};
	$nrow++;
} close IN; close ROW_DATA;

$mean = $sum/($ncol*$nrow);
$stdevSum = 0;

for ($i = 0; $i < @H; $i++) {
	$col_sum = 0; $col_mean = 0; $col_stdev_sum = 0; $col_stdev = 0; $col_nonzero = 0; $col_min = $COL_values[$i][0]; $col_max = $COL_values[$i][0];
	for ($j = 0; $j < @{$COL_values[$i]}; $j++) {
		$val = $COL_values[$i][$j];
		$col_sum += $val;
		if ($val > $col_max) {$col_max = $val};
		if ($val < $col_min) {$col_min = $val};
		if ($val != 0) {$col_nonzero++};
	}
	$col_mean = $col_sum/@{$COL_values[$i]};
	for ($j = 0; $j < @{$COL_values[$i]}; $j++) {
		$val = $COL_values[$i][$j];
		$col_stdev_sum += ($val - $col_mean)**2;
		$stdev_sum += ($val - $mean)**2;
		print COL_DIST "$H[$i]\t$val\n";
	}
	$col_stdev = sqrt($col_stdev_sum/@{$COL_values[$i]});
	print COL_DATA "$H[$i]\t$col_mean\t$col_stdev\t$col_min\t$col_max\t$col_sum\t$col_nonzero\n";
	if ($col_nonzero>0) {$ncol_nonzero++};
} close COL_DATA; close COL_DIST;

$stdev = sqrt($stdev_sum/($ncol*$nrow));

print STATS "Matrix        $ARGV[0]
Columns       $ncol
Rows          $nrow
Mean          $mean
Stdev         $stdev
Min           $minVal
Max           $maxVal
Nonzero_rows  $nrow_nonzero
Nonzero_cols  $ncol_nonzero
Nonzero_vals  $nonzero
Matrix_sum    $sum
"; close STATS;

}
1;
