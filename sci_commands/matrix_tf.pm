package sci_commands::matrix_tf;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_tf");

sub matrix_tf {

@ARGV = @_;

getopts("O:N:", \%opt);

$die2 = "
scitools matrix-tf [options] [counts matrix]
   or    tf

Term frequency normalization of matrix

Options:
   -O   [STR]   Output prefix (default is [input].tf)
   -N   [VAL]   Normalization constant (def = row number)

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix$//};

read_matrix_stats($ARGV[0]);
if (!defined $opt{'N'}) {$opt{'N'} = $matrix_rowNum};

open IN, "$ARGV[0]";
open OUT, ">$opt{'O'}.tf";
$h = <IN>; print OUT "$h";
while ($l = <IN>) {
	chomp $l;
	@P = split (/\t/, $l);
	$rowID = shift(@P);
	print OUT "$rowID";
	for ($i = 0; $i < @P; $i++) {
		$tf = ($P[$i]/$MATRIX_CellID_signal{$MATRIX_COLNAMES[$i]});
		$score = $tf*$opt{'N'};
		print OUT "\t$score";
	}
	print OUT "\n";
} close IN; close OUT;

}
1;
