package sci_commands::matrix_tfidf;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_tfidf");

sub matrix_tfidf {

@ARGV = @_;

getopts("O:", \%opt);

$die2 = "
scitools matrix-tfidf [options] [counts matrix]
   or    tfidf

Options:
   -O   [STR]   Output prefix (default is [input].tfidf)

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix$//};

read_matrix_stats($ARGV[0]);

open IN, "$ARGV[0]";
open OUT, ">$opt{'O'}.tfidf";
$h = <IN>; print OUT "$h";
while ($l = <IN>) {
	chomp $l;
	@P = split (/\t/, $l);
	$rowID = shift(@P);
	print OUT "$rowID";
	for ($i = 0; $i < @P; $i++) {
		$tf = ($P[$i]/$MATRIX_CellID_signal{$MATRIX_COLNAMES[$i]});
		$idf = (log(1+($matrix_colNum/($MATRIX_feature_signal{$rowID}+1))));
		$score = ($tf*$idf);
		print OUT "\t$score";
	}
	print OUT "\n";
} close IN; close OUT;

}
1;
