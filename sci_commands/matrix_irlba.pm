package sci_commands::matrix_irlba;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_irlba");

sub matrix_irlba {

$dims = 50;

@ARGV = @_;
getopts("O:XR:D:C:", \%opt);

$die2 = "
scitools matrix-irlba [options] [input matrix or sparseMatrix values/tfidf file]
   or    irlba

Produces a dims file with irlba coordinates.
If sparseMatrix is provided, there must be an associated .cols file.

Options:
   -O   [STR]   Output prefix (default is [input].irlba.dims)
   -D   [INT]   Max dimensions to compute (def = $dims)
   -C   [STR]   Cols file (for sparseMatrix - will try to auto-detect)
   -X           Retain intermediate files (def = delete)
   -R   [STR]   Rscript call (def = $Rscript)

Note: Requires irlba R package

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {
	$prefix = $ARGV[0];
	$prefix =~ s/\.gz$//;
	$prefix =~ s/\.(tfidf|matrix)$//;
	$prefix =~ s/\.sparseMatrix$//;
	$opt{'O'} = $prefix;
}
if ($ARGV[0] =~ /sparseMatrix/i) {
	$sparse = 1;
	if (defined $opt{'C'}) {
		$col_file = $opt{'C'};
	} else {
		if (-e "$prefix.sparseMatrix.cols") {
			$col_file = "$prefix.sparseMatrix.cols";
		} elsif (-e "$prefix.sparseMatrix.cols.gz") {
			$col_file = "$prefix.sparseMatrix.cols.gz";
		} else {
			die "ERROR: Cannot detect cols file (e.g. $prefix.sparseMatrix.cols), please provide as -C\n";
		}
	}
} else {$sparse = 0};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};

open R, ">$opt{'O'}.irlba.r";
print R "
library(Matrix)
library(irlba)
library(ggplot2)
IN<-as.matrix(read.table(\"$ARGV[0]\"))";

if ($sparse > 0.5) {
	print R "
IN<-sparseMatrix(i=IN[,1],j=IN[,2],x=IN[,3])
COLS<-read.table(\"$col_file\")
colnames(IN)<-COLS\$V1";
}

print R "
I<-irlba(IN, $dims)
rownames(I\$v)<-colnames(IN)

var<-as.data.frame(I\$d,ncol=1)
colnames(var)<-c(\"d\")
var\$var_names<-paste0(\"DIM_\",row.names(var))
var\$var_names<-factor(var\$var_names,levels=var\$var_names)

var\$variance<-prop.table(var\$d^2)
ggplot()+geom_histogram(aes(x=var\$var_names,y=var\$variance),stat=\"identity\")+xlab(\"PC\")+ylab(\"Variance Explained\")+theme_minimal()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(\"$opt{'O'}.irlba_variance.png\")

df<-as.data.frame(I\$v)
PCA_colnames<-c()
for (i in 1:$dims) {
  p<-paste(\"DIM\",i, sep = \"\")
  PCA_colnames <- append(PCA_colnames, p)
}
colnames(df)<-PCA_colnames
dft<-t(df)

write.table(dft, file=\"$opt{'O'}.irlba_$dims.matrix\", quote=FALSE, sep=\"\\t\")

";
close R;

system("$Rscript $opt{'O'}.irlba.r");

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.irlba.r");
}

}
1;
