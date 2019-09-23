package sci_commands::matrix_anova;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_anova");

sub matrix_anova {

@ARGV = @_;


getopts("O:a:F:R:X", \%opt);

$die2 = "
scitools plot-dims [options] [matrix file] [annot]

Options:
   -O   [STR]   Output prefix (default is dims file 1 prefix)
   -a 	[STR]		Comma separated list of annoations to include in anova
   -F   [NUM]	fdr.level for correction def:(0.05)
   -R   [STR]   Rscript call (def = $Rscript)
   -X           Do not delete intermediate files (def = delete)

Note: Requires R package qvalue

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $ARGV[1]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};
if (!defined $opt{'F'}) {$fdr = 0.05};
read_matrix($ARGV[0]);
read_annot($ARGV[1]);

if (defined $opt{'a'}) {
	@ANNOT_LIST = split(/,/, $opt{'a'});
	foreach $annot (@ANNOT_LIST) {
		$ANNOT_include{$annot} = 1;
	}
}

open ANOVA, ">$opt{'O'}.Annova_sum.txt";
print ANOVA "Feature_name\tF\tpval\n"; 
close ANOVA;
	
	foreach $feature (keys %MATRIX_feature_nonZero) {
		$feature_polished = $feature;
		$feature_polished =~ s/(:|;|'|\.|\(|\)|\{|\}|\[|\]|\s|\|)/_/;
		open VALS, ">$opt{'O'}.$feature_polished.values";
		foreach $cellID (keys %CELLID_FEATURE_value) {
			$annot=$CELLID_annot{$cellID};
			if (defined $CELLID_FEATURE_value{$cellID}{$feature}) {
				if (defined $opt{'a'})
				{
					if (defined $ANNOT_include{$annot})
				{
				print VALS "$cellID\t$CELLID_FEATURE_value{$cellID}{$feature}\t$annot\n";
				}
				}
				else
				{
				print VALS "$cellID\t$CELLID_FEATURE_value{$cellID}{$feature}\t$annot\n";	
				}
			}
		} close VALS;



open R, ">$opt{'O'}.$feature_polished.anova.r";
print R "
library(ggplot2)
IN<-read.table(\"$opt{'O'}.$feature_polished.values\")

feature<-\"$feature_polished\"
fit <- aov(V2 ~ V3, data=IN)
F<-summary(fit)[[1]][[\"F value\"]][1]
P<-summary(fit)[[1]][[\"Pr(>F)\"]][1]

output<-c(feature,F,P)
write(output,file=\"$opt{'O'}.Annova_sum.txt\",append=TRUE,ncolumns=3,sep = \"\\t\")

#individual comparisons with TukeyHSD
comp<-TukeyHSD(fit)
write.table(as.matrix(comp\$V3),file=\"$opt{'O'}.$feature_polished.TukeyHSD.txt\",quote=FALSE,sep=\"\\t\",row.names=T,col.names=T)";

close R;

system("$Rscript $opt{'O'}.$feature_polished.anova.r");



if (!defined $opt{'X'}) {system("rm -f $opt{'O'}.$feature_polished.values $opt{'O'}.$feature_polished.anova.r")};





}

open R, ">$opt{'O'}.final.anova.r";
print R "
library(qvalue)
IN<-read.delim(\"$opt{'O'}.Annova_sum.txt\")
qobj <- qvalue(IN\$pval, fdr.level=$fdr, pi0.method=\"bootstrap\", adj=1.2)
OUT<-data.frame(Feature_name=IN\$Feature_name,F=IN\$F,pval=IN\$pval,qval=qobj\$qvalues,lfdr=qobj\$lfdr,sign=qobj\$significant)


write.table(as.matrix(OUT),file=\"./$opt{'O'}.annova.qval.txt\",quote=FALSE,sep=\"\\t\",row.names=T,col.names=T)";
close R;
system("$Rscript $opt{'O'}.final.anova.r");
if (!defined $opt{'X'}) {system("rm -f $opt{'O'}.final.anova.r")};


}
1;
