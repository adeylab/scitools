package sci_commands::matrix_pull;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_pull");

sub matrix_pull {

@ARGV = @_;

getopts("Fr", \%opt);

$die2 = "
scitools matrix-pull [options] [matrix] [column,row] [column,row] etc... OR [file]
   or    pull-matrix

Pull a value or list of values from a matrix. Outputs to STDOUT, tab delimited:
column (tab) row (tab) value

Can specify 'all' for all rows / columns

If a file is specified, it should be tab-delimited: column (tab) row

Options:
   -F   File is the input (will try to auto detect)
   -r   Reverse output order (will be row (tab) col (tab) value)

";

if (!defined $ARGV[1]) {die $die2};

%ROW_COL = ();
if (defined $opt{'F'} || $ARGV[1] !~ /,/) {
	open IN, "$ARGV[1]";
	while ($l = <IN>) {
		chomp $l;
		($col,$row) = split(/\t/, $l);
		if ($col =~ /all/i) {$col = "all"};
		if ($row =~ /all/i) {$row = "all"};
		$ROW_COL{$row}{$col} = 1;
	} close IN;
} else {
	for ($i = 1; $i < @ARGV; $i++) {
		($col,$row) = split(/,/, $ARGV[$i]);
		if ($col =~ /all/i) {$col = "all"};
		if ($row =~ /all/i) {$row = "all"};
		$ROW_COL{$row}{$col} = 1;
	}
}

open IN, "$ARGV[0]";
$h = <IN>; chomp $h; @H = split(/\t/, $h);
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	$row = shift(@P);
	if (defined $ROW_COL{$row} || defined $ROW_COL{'all'}) {
		for ($i = 0; $i < @H; $i++) {
			if (defined $ROW_COL{$row}{$H[$i]} ||
				defined $ROW_COL{$row}{'all'} ||
				defined $ROW_COL{'all'}{'all'} ||
				defined $ROW_COL{'all'}{$H[$i]}) {
				if (!defined $opt{'r'}) {
					print "$H[$i]\t$row\t$P[$i]\n";
				} else {
					print "$row\t$H[$i]\t$P[$i]\n";
				}
			}
		}
	}
} close IN;

}
1;