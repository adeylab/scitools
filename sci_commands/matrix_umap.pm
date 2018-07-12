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

getopts("O:n:d:m:D:XP:", \%opt);

$die2 = "
scitools umap [options] [matrix]

Options:
   -O   [STR]   Output prefix (default is matrix file prefix)
   -D   [INT]   Dimensions to embed UMAP in (def = $dims)
   -n   [INT]   number of neighbors to use for UMAP high dim embedding (def = $neigh)
   -d   [INT]   min distance for mapping (def = $mdist)
   -m   [STR]   metric to use (def = $metric)
   -X           Retain intermediate files (def = delete)
   -P   [STR]   python script call (def = $Pscript)
   
Note: Requires python, numpy, and umap to be installed and callable
      This command is a wrapper for executing the python code.

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
$opt{'O'} =~ s/\.matrix$//;
if (defined $opt{'D'}) {$dims = $opt{'D'}};
if (defined $opt{'n'}) {$neigh = $opt{'n'}};
if (defined $opt{'d'}) {$mdist = $opt{'d'}};
if (defined $opt{'m'}) {$metric = $opt{'m'}};
if (defined $opt{'P'}) {$Pscript = $opt{'P'}};
read_matrix($ARGV[0]);

open OUT, ">$opt{'O'}.UMAP.py";
print OUT"
import numpy
import umap
data_matrix=numpy.loadtxt(\"$ARGV[0]\",skiprows=1,usecols=range(1,$matrix_colNum))

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



if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.temp.UMAP.dims $opt{'O'}.UMAP.py");
}

}
1;
