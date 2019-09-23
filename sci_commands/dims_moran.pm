package sci_commands::dims_moran;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
use File::Basename;

@EXPORT = ("dims_moran");

sub dims_moran {

@ARGV = @_;
getopts("O:Xn:k:c:R", \%opt);

$die2 = "
scitools dims_moran [options] [Scores Matrix] [dims file]

Will take in a dims file, perform k-nearst neighbor clustering and perform 
Moran's I test to look for spectral biases in the plot 
(viz enrichment of TF/CCAN accessibility in different cell groupings)
Output is a p-value ordered matrix with each feature as a row.
P-value output is generated through monte-carlo simulation.

[Scores Matrix]     File generated through the atac_chromvar command. 
					By default named: deviation_scores.matrix
[dims file] 		Dimension dims file (2D Only) generated through scitools matrix_[tsne|umap|swne]

Options:
   -O   [STR]   Output prefix (default is [input].moran.matrix)
   -X           Retain intermediate files (def = delete)
   -n 	[INT] 	Number of Monte-Carlo simulations for empirical t-test on Moran's I value. (Def: 99)
   -k 	[INT] 	Number of cells to define a neighborhood in the 2D dims file. (Def: 25)
   -c 	[INT] 	Number of cores for monte-carlo simulation (Def: 10)
   -R   [STR]   Rscript call (def = $Rscript)
";

if (!defined $ARGV[0]) {die $die2};
if (!defined $ARGV[1]) {die $die2};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};
if (!defined $opt{'n'}) {$opt{'n'}=99}
if (!defined $opt{'k'}) {$opt{'k'}=25}
if (!defined $opt{'c'}) {$opt{'c'}=10}
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.dims$//};

open R, ">$opt{'O'}.moran_test.r";
print R "
library(spdep)
library(RANN)
library(parallel)

dims<-read.table(\"$ARGV[1]\",header=F,row.names=1)
colnames(dims)<-c(\"x\",\"y\")
scores_mat<-read.table(\"$ARGV[0]\")
warning(\"Read in files.\")

#generate k-nearest neighbor neiborhoods
knn<-knearneigh(as.matrix(dims), k=$opt{'k'}, longlat = FALSE, RANN=TRUE)
warning(\"Generated kNN.\")

#read in dims to nb matrix
w<-knn2nb(knn,row.names=row.names(knn\$x))

#generate spatial weight matrix
wm<-nb2mat(w,style=\'W\',zero.policy=FALSE)

#SUBET spatial weight matrix
warning(\"Using $opt{'c'} cores for Moran's Test.\")

m_test_list<-mclapply(1:nrow(scores_mat),FUN=function(x){
out<-tryCatch({
	dat_temp<-t(na.omit(t(scores_mat[x,])))
wm_temp<-wm[row.names(wm) \%in\% colnames(dat_temp),row.names(wm) \%in\% colnames(dat_temp)]
dat_temp<-dat_temp[colnames(dat_temp) \%in\% row.names(wm_temp)]
rwm <- mat2listw(wm_temp, style=\'W\')
#Use monte-carlo approach for significance
m_test<-moran.mc(dat_temp, rwm, nsim=$opt{'n'},zero.policy=TRUE,na.action=na.omit)
out<-c(row.names(scores_mat)[x],\"PASS\",m_test\$statistic, m_test\$parameter, m_test\$p.value)
}, error=function(e) {out<-c(row.names(scores_mat)[x],\"FAIL\",\"NA\",\"NA\",\"NA\")})
return(out)
},mc.cores=$opt{'c'})

warning(\"Moran's Test complete. Formatting and writing out data table.\")

m_test_list_2<-as.data.frame(do.call(rbind,lapply(m_test_list, function(x) unlist(x))))
colnames(m_test_list_2)<-c(\"feature\",\"PASS/FAIL\",\"I_value\",\"observed_rank\",\"montecarlo_pvalue\")
m_test_list_2<-m_test_list_2[order(m_test_list_2\$montecarlo_pvalue),]
write.table(m_test_list_2,\"$opt{'O'}.moran.matrix\",sep=\"\\t\",col.names=T,row.names=F,quote=F)
";
close R;
system("$Rscript $opt{'O'}.moran_test.r");
if (!defined $opt{'X'}) {
        system("rm -f $opt{'O'}.moran_test.r");
}
}
1;



