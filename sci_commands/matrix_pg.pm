package sci_commands::matrix_pg;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_pg");

sub matrix_pg {

@ARGV = @_;


getopts("O:XR:k:", \%opt);

$die2 = "
scitools matrix-pg [options] [input matrix or dims] 

This scripts uses Louvain method to find clusters of cells based on provided dims or matrix

Options:
   -O   [STR]   Output prefix (default is [input].pg.annot)
   -k	[INT]	Number of nearest-neighbor cells to use for clustering. (Default:50)
   -X           Remove intermediate files (def = Keep if unflagged)
   -R   [STR]   R script call (def = $Rscript)

";


if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'k'}) {$opt{'k'}=50};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
$opt{'O'} =~ s/\.matrix$//;
$opt{'O'} =~ s/\.dims$//;
if (defined $opt{'P'}) {$Pscript = $opt{'P'}};
read_matrix($ARGV[0]);

open OUT, ">$opt{'O'}.k$opt{'k'}.phenograph.R";
print OUT"
library(Rphenograph)
";


if ($ARGV[0] =~ /\.dims/) {
print OUT "
dat_matrix<-as.data.frame(read.table(\"$ARGV[0]\",header=T))
";
} else {
print OUT "
dat_matrix<-as.data.frame(t(read.table(\"$ARGV[0]\",header=T)))
";
}


#detect communities by given k value.
print OUT "
cellID<-row.names(dat_matrix)

dat_pg<-Rphenograph(dat_matrix,k=$opt{'k'})
pg_clusterID<-paste(\"pg_\",membership(dat_pg[[2]]),sep=\"\")
dat_out<-do.call(rbind, Map(data.frame, A=cellID, B=pg_clusterID))

write.table(dat_out,file=\"$opt{'O'}.k$opt{'k'}.pg.annot\",col.names=F,row.names=F,quote=F,sep=\"\\t\")
";
close OUT;
system("$Rscript $opt{'O'}.k$opt{'k'}.phenograph.R");

if (defined $opt{'X'}) {
	system("rm -f $opt{'O'}.k$opt{'k'}.phenograph.R");
}

}
1;