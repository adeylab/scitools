package sci_commands::matrix_naomit;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_naomit");

sub matrix_naomit {

@ARGV = @_;

getopts("O:", \%opt);

$die2 = "
scitools matrix-naomit [options] [matrix]
   or    naomit-matrix
   
Note: Print out number of NA's in matrix then remove NA's from matrix

Options:
   -O   [STR]   Output prefix (default is [input].naomit.matrix)

";

if (!defined $ARGV[0]) {die $die2};

if (!defined $opt{'O'}) {
	$opt{'O'} = $ARGV[0];
	$opt{'O'} =~ s/\.matrix$//;
}
#count how many
$matrix_NA=0;
open IN, $ARGV[0];
open OUT, ">$opt{'O'}.naomit.matrix";
$head = <IN>; print OUT "$head";
while ($line = <IN>) {
	if ($line=~ m/NA/) {
		$matrix_NA++;
	} else {
		print OUT "$line";
	}
} close IN; close OUT;

open LOG, ">$opt{'O'}.log";
$ts = localtime(time);
print LOG "$ts scitools naomit
Matrix file = $ARGV[0]
Total rows with NA: $matrix_NA. \n
";
close LOG;

}
1;
