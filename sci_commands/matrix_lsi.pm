package sci_commands::matrix_lsi;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_lsi");

sub matrix_lsi {

@ARGV = @_;
# Defaults
$range_default = "1-15";

getopts("O:D:d:XR:", \%opt);

$die2 = "
scitools matrix-lsi [options] [tfidf matrix]
   or    lsi

Options:
   -O   [STR]   Output prefix (def = [input].tfidf.LSI.matrix)
   -D   [RNG]   Range of dimensions to include (range format)
                (e.g. 1-5,6,8; def = $range_default)
   -X           Retain intermediate files (def = delete)
   -R   [STR]   Rscript call (def = $Rscript)

Note: Requires svd R package

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
if (!defined $opt{'D'}) {$opt{'D'} = $range_default};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};
read_ranges($opt{'D'});

open PARAMS, ">$opt{'O'}.LSI.params";
print PARAMS "scitools matrix-lsi ON $ARGV[0]
Dimensions: $range_R_set\n";
close PARAMS;

open R, ">$opt{'O'}.LSI.SVD.r";
print R "
library(svd)
IN<-read.table(\"$ARGV[0]\")
SVD<-svd(as.matrix(IN))
ds<-diag(SVD\$d[$range_R_set])
us<-as.matrix(SVD\$u[,$range_R_set])
vs<-as.matrix(SVD\$v[,$range_R_set])
write.table(ds,file=\"$opt{'O'}.LSI.SVD.d\",sep=\"\\t\",row.names=FALSE,col.names=FALSE)
write.table(us,file=\"$opt{'O'}.LSI.SVD.u\",sep=\"\\t\",row.names=FALSE,col.names=FALSE)
write.table(vs,file=\"$opt{'O'}.LSI.SVD.v\",sep=\"\\t\",row.names=FALSE,col.names=FALSE)
";
close R;

system("$Rscript $opt{'O'}.LSI.SVD.r");

open DIAG, "$opt{'O'}.LSI.SVD.d";
$rowNum = 0;
while ($l = <DIAG>) {
	chomp $l;
	@P = split(/\t/, $l);
	$diagSum+=$P[$rowNum];
	$rowNum++;
} close DIAG;

open DIAG, "$opt{'O'}.LSI.SVD.d";
open OUT, ">$opt{'O'}.LSI.SVD.d.norm";

$rowNum = 0;
while ($l = <DIAG>) {
	chomp $l;
	@P = split(/\t/, $l);
	$P[$rowNum] = ($P[$rowNum]/$diagSum);
	$out = join("\t", @P);
	print OUT "$out\n";
	$rowNum++;
} close DIAG; close OUT;

open R, ">$opt{'O'}.LSI.SVD.reconstruct.r";
print R "
D<-as.matrix(read.table(\"$opt{'O'}.LSI.SVD.d.norm\"))
U<-as.matrix(read.table(\"$opt{'O'}.LSI.SVD.u\"))
V<-as.matrix(read.table(\"$opt{'O'}.LSI.SVD.v\"))
LSI<-t(scale(V %*% D %*% t(U)))
write.table(LSI,file=\"$opt{'O'}.LSI.SVD.reconstructed\",sep=\"\\t\",row.names=FALSE,col.names=FALSE)
";
close R;

system("$Rscript $opt{'O'}.LSI.SVD.reconstruct.r");

open IN, "$ARGV[0]";
open LSI, "$opt{'O'}.LSI.SVD.reconstructed";
open OUT, ">$opt{'O'}.LSI.matrix";
$h = <IN>; chomp $h;
print OUT "$h\n";
while ($l = <IN>) {
	$v = <LSI>;
	chomp $l; chomp $v;
	@P = split(/\t/, $l);
	$siteID = shift(@P);
	print OUT "$siteID\t$v\n";
} close IN; close OUT; close LSI;

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.LSI.SVD.r $opt{'O'}.LSI.SVD.d $opt{'O'}.LSI.SVD.d.norm $opt{'O'}.LSI.SVD.u $opt{'O'}.LSI.SVD.v $opt{'O'}.LSI.SVD.reconstructed $opt{'O'}.LSI.SVD.reconstruct.r");
}

}
1;
