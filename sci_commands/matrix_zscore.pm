package sci_commands::matrix_zscore;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_zscore");

sub matrix_zscore {

@ARGV = @_;
getopts("O:CRG", \%opt);

$die2 = "
scitools matrix-zscore [options] [matrix]
   or    zscore-matrix
   
Performs z-scoring on the matrix

Options:
   -O   [STR]   Output prefix (default is [input].zscore.matrix)
   -R           Z-score rows
   -C           Z-score columns
   -G           Z-score globally

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix$//};
if ((defined $opt{'R'} && defined $opt{'C'}) ||
    (defined $opt{'R'} && defined $opt{'G'}) ||
	(defined $opt{'C'} && defined $opt{'G'}) ||
	(!defined $opt{'R'} && !defined $opt{'C'} && !defined $opt{'G'})) {die "ERROR: Please perform row, column, OR global z-scoring.\n"};

if (defined $opt{'R'}) {
	open OUT, ">$opt{'O'}.R_zscore.matrix";
} elsif (defined $opt{'C'}) {
	open OUT, ">$opt{'O'}.C_zscore.matrix";
	@SUMS = (); @COUNTS = (); @MEANS = (); @STDEV_SUM = (); @STDEV = ();
} elsif (defined $opt{'G'}) {
	open OUT, ">$opt{'O'}.G_zscore.matrix";
	$sum = 0; $count = 0; $stdev_sum = 0; $stdev = 0;
}

open IN, "$ARGV[0]";
$h = <IN>; chomp $h; @H = split(/\t/, $h);
print OUT "$h\n";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	$siteID = shift(@P);
	if (defined $opt{'R'}) {
		$sum = 0;
		for ($i = 0; $i < @P; $i++) {
			$sum += $P[$i];
		}
		$mean = $sum/(@P);
		$stdev_sum = 0;
		for ($i = 0; $i < @P; $i++) {
			$stdev_sum += ($mean-$P[$i])**2;
		}
		$stdev = sqrt($stdev_sum/(@P));
		if ($stdev != 0) {
			print OUT "$siteID";
			for ($i = 0; $i < @P; $i++) {
				$zscore = ($P[$i]-$mean)/$stdev;
				print OUT "\t$zscore";
			}
			print OUT "\n";
		} else {
			print STDERR "WARNING: row $siteID has a stdev of 0.\n";
		}
	} elsif (defined $opt{'C'}) {
		for ($i = 0; $i < @P; $i++) {
			$SUMS[$i] += $P[$i];
			$COUNTS[$i]++;
		}
	} elsif (defined $opt{'G'}) {
		for ($i = 0; $i < @P; $i++) {
			$sum += $P[$i]; $count++;
		}
	}
} close IN;

if (defined $opt{'R'}) {close OUT; exit};

if (defined $opt{'C'}) {
	for ($i = 0; $i < @SUMS; $i++) {
		$MEANS[$i] = $SUMS[$i]/$COUNTS[$i];
	}
} else {
	$mean = $sum/$count;
}

open IN, "$ARGV[0]";
$h = <IN>;
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	$siteID = shift(@P);
	if (defined $opt{'C'}) {
		for ($i = 0; $i < @P; $i++) {
			$STDEV_SUM[$i] += ($MEANS[$i]-$P[$i])**2;
		}
	} elsif (defined $opt{'G'}) {
		for ($i = 0; $i < @P; $i++) {
			$stdev_sum += ($mean-$P[$i])**2;
		}
	}
} close IN;

if (defined $opt{'C'}) {
	for ($i = 0; $i < @SUMS; $i++) {
		$STDEV[$i] = sqrt($STDEV_SUM[$i]/$COUNTS[$i]);
	}
} else {
	$stdev = sqrt($stdev_sum/$count);
}

open IN, "$ARGV[0]";
$h = <IN>;
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	$siteID = shift(@P);
	print OUT "$siteID";
	if (defined $opt{'C'}) {
		for ($i = 0; $i < @P; $i++) {
			$zscore = ($P[$i]-$MEANS[$i])/$STDEV[$i];
			print OUT "\t$zscore";
		}
	} elsif (defined $opt{'G'}) {
		for ($i = 0; $i < @P; $i++) {
			$zscore = ($P[$i]-$mean)/$stdev;
			print OUT "\t$zscore";
		}
	}
	print OUT "\n";
} close IN;

}
1;
