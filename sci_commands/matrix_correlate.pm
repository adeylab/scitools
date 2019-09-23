package sci_commands::matrix_correlate;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_correlate");

sub matrix_correlate {

@ARGV = @_;

# Defaults
$method = "pearson";
getopts("O:M:R", \%opt);

$die2 = "
scitools matrix-correlate [options] [input matrix]

Will correlate a matrix and produce a correlation matrix.

   -O   [STR]   Output prefix (default is input prefix)
   -M   [STR]   Method (pearson or spearman, def = $method)
   -R           Correlate rows (def = columns)

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix$//};
if (defined $opt{'M'}) {
	if ($opt{'M'} =~ /^p/i) {$method = "pearson"}
	elsif ($opt{'M'} =~ /^s/i) {$method = "spearman"}
	else {die "ERROR: -M must be p/pearson or s/spearman!\n$die2"}
}

open R, ">$opt{'O'}.$method.r";
if (!defined $opt{'R'}) {
	print R "IN<-read.table(\"$ARGV[0]\")\n";
} else {
	print R "IN<-t(read.table(\"$ARGV[0]\"))\n";
}
print R "COR<-cor(IN,method=\"$method\")
write.table(COR,file=\"$opt{'O'}.$method.matrix\",quote=FALSE,col.names=TRUE,row.names=TRUE,sep=\"\\t\")\n";
close R;

system("Rscript $opt{'O'}.$method.r && rm -f $opt{'O'}.$method.r");

}
1;