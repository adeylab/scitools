package sci_commands::matrix_binarize;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_binarize");

sub matrix_binarize {

@ARGV = @_;

getopts("O:T:", \%opt);
$thresh = 1;

$die2 = "
scitools matrix-binarize [options] [binarized matrix]
   or    binarize-matrix

Binarizes a matrix.

Options:
   -O   [STR]   Output prefix (default is [input prefix].binarized.matrix)
   -T   [VAL]   Threshold for a 1 vs 0 (def = $thresh)
   -A           Use absolute value (negatives become -1)

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix$//};
if (defined $opt{'T'}) {$thresh = $opt{'T'}};


open IN, "$ARGV[0]";
open OUT, ">$opt{'O'}.binarized.matrix";
$h = <IN>; print OUT "$h";
while ($l = <IN>) {
	chomp $l;
	@P = split (/\t/, $l);
	$rowID = shift(@P);
	print OUT "$rowID";
	for ($i = 0; $i < @P; $i++) {
		if (!defined $opt{'A'}) {
			if ($P[$i] >= $thresh) {
				print OUT "\t1";
			} else {
				print OUT "\t0";
			}
		} else {
			if (abs($P[$i]) >= $thresh) {
				if ($P[$i] > 0) {
					print OUT "\t1";
				} else {
					print OUT "\t-1";
				}
			} else {
				print OUT "\t0";
			}
		}
	}
	print OUT "\n";
} close IN; close OUT;

}
1;
