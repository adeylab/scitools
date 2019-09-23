package sci_commands::matrix_sparse;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_sparse");

sub matrix_sparse {

@ARGV = @_;
getopts("O:z", \%opt);

$die2 = "

scitools matrix-sparse (options) [matrix file]
  or     make-sparse

Options:
   -O   [STR]   Output prefix (def = matrix prefix)
   -z           Gzip output

";

if (!defined $ARGV[0]) {die $die2};

if (!defined $opt{'O'}) {
	$opt{'O'} = $ARGV[0];
	$opt{'O'} =~ s/\.matrix//i;
}

open IN, $ARGV[0];
$header = <IN>; chomp $header;
@CELLIDs = split(/\s/, $header);
if (defined $opt{'z'}) {
	open OUT, "| $gzip > $opt{'O'}.sparseMatrix.cols.gz";
} else {
	open OUT, ">$opt{'O'}.sparseMatrix.cols";
}
for ($i = 0; $i < @CELLIDs; $i++) {
	print OUT "$CELLIDs[$i]\n";
} close OUT;
$rowNum = 0;
if (defined $opt{'z'}) {
	open ROW, "| $gzip > $opt{'O'}.sparseMatrix.rows.gz";
} else {
	open ROW, ">$opt{'O'}.sparseMatrix.rows";
}
if (defined $opt{'z'}) {
	open MTX, "| $gzip > $opt{'O'}.sparseMatrix.values.gz";
} else {
	open MTX, ">$opt{'O'}.sparseMatrix.values";
}
#print MTX "rowID\tcellID\tvalue\n";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	$rowNum++;
	$rowID = shift(@P);
	print ROW "$rowID\n";
	for ($i = 0; $i < @P; $i++) {
		if ($P[$i] != 0) {
			print MTX "$rowNum\t".($i+1)."\t$P[$i]\n";
		}
	}
} close ROW; close MTX;

}
1;
