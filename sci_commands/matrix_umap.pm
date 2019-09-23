package sci_commands::matrix_umap;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_umap");

sub matrix_umap {

@ARGV = @_;
# Defaults
$dims = 2;
$neigh = 30;
$mdist=0.1;
$metric="euclidean";

getopts("O:n:d:m:D:XP:pR:", \%opt);

$die2 = "
scitools umap [options] [matrix/dims]
UMAP is an alternative dimensionality reduction technique to [tsne|pca|swne]. 
There are two implementations available, both R and python versions. 
R takes longer to run and is less memory efficient
but does not require a conda environment set up.
R implementation also has added features.

To run the python implementation:
  1. Enter into the command line: \"conda activate py3\"
  2. Run scitools umap with the -p flag.

Options:
   -O   [STR]   Output prefix (default is matrix file prefix)
   -D   [INT]   Dimensions to embed UMAP in (def = $dims)
   -n   [INT]   number of neighbors to use for UMAP high dim embedding (def = $neigh)
   -d   [INT]   min distance for mapping (def = $mdist)
   -m   [STR]   metric to use (def = $metric)
   -X           Retain intermediate files (def = delete)
   -p           Run the Python version of UMAP (faster but requires conda environment)
   -P   [STR]   python script call (def = $Pscript)
   -R   [STR]   R script call (def = $Rscript)

   
Note: Requires python, numpy, and umap to be installed and callable
      This command is a wrapper for executing the python code.

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
$opt{'O'} =~ s/\.matrix$//; $opt{'O'} =~ s/\.dims$//;
if (defined $opt{'D'}) {$dims = $opt{'D'}};
if (defined $opt{'n'}) {$neigh = $opt{'n'}};
if (defined $opt{'d'}) {$mdist = $opt{'d'}};
if (defined $opt{'m'}) {$metric = $opt{'m'}};
if (defined $opt{'P'}) {$Pscript = $opt{'P'}};
read_matrix($ARGV[0]);

if (defined $opt{'p'}){


open OUT, ">$opt{'O'}.UMAP.py";
print OUT"
import numpy
import umap
";
if ($ARGV[0] =~ /\.dims/) {
print OUT "
data_matrix=numpy.loadtxt(\"$ARGV[0]\",skiprows=1,usecols=range(1,$matrix_colNum))
#data_matrix=ori_matrix.transpose()\n";
} else {
print OUT "data_matrix=numpy.loadtxt(\"$ARGV[0]\",skiprows=1,usecols=range(1,$matrix_colNum))\n";
}

print OUT "

fit = umap.UMAP(
       n_neighbors=$neigh,
       min_dist=$mdist,
       n_components=$dims,
       metric=\'$metric\'
   )
data=data_matrix.T
u = fit.fit_transform(data)
numpy.savetxt(\"$opt{'O'}.temp.UMAP.dims\",u,delimiter=\"\\t\")
";
close OUT;
system("$Pscript $opt{'O'}.UMAP.py");


open OUT, ">$opt{'O'}.UMAP.dims";
$counter=0;
open IN, "$opt{'O'}.temp.UMAP.dims"; 
while($l=<IN>)
{
$l =~ s/"//g;   
print OUT $MATRIX_COLNAMES[$counter]."\t".$l;
$counter++;
}
close(IN);
close OUT;

} else {

# R version (default)

open OUT, ">$opt{'O'}.UMAP.R";
print OUT "
library(umap)
";
if ($ARGV[0] =~ /\.dims$/) {
	print OUT "dat<-t(read.table(\"$ARGV[0]\",row.names=1,header=F))";
} else {
	print OUT "dat<-read.table(\"$ARGV[0]\",row.names=1,header=T)";
}
print OUT "

custom.config = umap.defaults
custom.config\$n_components = $dims

umap_dims<-umap(as.data.frame(t(dat)),method=c(\"naive\"), config = custom.config)
saveRDS(umap_dims,file=\"$opt{'O'}.UMAP.rds\")
write.table(as.matrix(umap_dims\$layout),file=\"$opt{'O'}.UMAP.dims\",col.names=F,row.names=T,sep=\"\\t\",quote=F)";

close OUT;
system("$Rscript $opt{'O'}.UMAP.R");

}

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.temp.UMAP.dims $opt{'O'}.UMAP.py $opt{'O'}.UMAP.R");
}

}
1;
