package sci_commands::plot_complexity;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("plot_complexity");

sub plot_complexity {

@ARGV = @_;
getopts("O:A:a:C:c:R:M", \%opt);

$die2 = "
scitools plot-complexity [options] [complexity file(s) can be comma separated]

Options:
   -O   [STR]   Output prefix (default is complexity file 1 prefix)
   -M           Run mixed model to determine read count cutoff for cells (def = no)
   -A   [STR]   Annotation file (to color code points)
   -a   [STR]   Comma separated list of annoations to include in plot
                (requires -A to be specified)
   -C   [STR]   Color coding file (annot (tab) #hexColor)
   -c   [STR]   Color coding string
                Annot=#hexColor,Annot2=#hexColor
   -R   [STR]   Rscript call (def = $Rscript)

Note: Requires ggplot2 R package

";

if (!defined $ARGV[0]) {die $die2};
if (defined $opt{'a'} && !defined $opt{'A'}) {die "\nMust provide an annotaiton file (-A) if specifying annotations to plot (-a)!\n$die2"};
if (defined $opt{'C'} && defined $opt{'c'}) {die "\nSpecify either a color string (-c) or a color coding file (-C), not both!\n$die2"};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};
if (defined $opt{'A'}) {read_annot($opt{'A'})};
if (defined $opt{'a'}) {
	@ANNOT_LIST = split(/,/, $opt{'a'});
	foreach $annot (@ANNOT_LIST) {
		$ANNOT_include{$annot} = 1;
	}
}
if (defined $opt{'C'}) {read_color_file($opt{'C'})};
if (defined $opt{'c'}) {read_color_string($opt{'c'})};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
$opt{'O'} =~ s/\.txt$//;

read_complexity($ARGV[0]);

open OUT, ">$opt{'O'}.plot.txt";
foreach $cellID (keys %CELLID_complexity) {
	if (defined $opt{'a'}) {
		$annot = $CELLID_annot{$cellID};
		if (defined $ANNOT_include{$annot} && defined $CELLID_annot{$cellID}) {
			print OUT "$cellID\t$CELLID_annot{$cellID}\t$CELLID_uniq_reads{$cellID}\t$CELLID_complexity{$cellID}\n";
		}
	} elsif (defined $opt{'A'} && defined $CELLID_annot{$cellID}) {
		print OUT "$cellID\t$CELLID_annot{$cellID}\t$CELLID_uniq_reads{$cellID}\t$CELLID_complexity{$cellID}\n";
	} else {
		print OUT "$cellID\tCell\t$CELLID_uniq_reads{$cellID}\t$CELLID_complexity{$cellID}\n";
	}
} close OUT;

open R, ">$opt{'O'}.plot.r";
print R "
library(ggplot2)
IN<-read.table(\"$opt{'O'}.plot.txt\")
PLT<-ggplot(data=subset(IN,V4<100&V4>0)) + theme_bw() +
   geom_point(aes(V4,log10(V3),color=V2),size=1,alpha=0.3) +
   geom_density2d(aes(V4,log10(V3),color=V2),size=0.3) +
";
if (defined $opt{'C'} || defined $opt{'c'}) {
	print R "   scale_colour_manual(values = c($color_mapping)) +
	guides(colour = guide_legend(override.aes = list(size=4))) +";
}
print R "
   scale_x_continuous(limits=c(0,100)) +
   scale_y_continuous(limits=c(0,6)) +
   xlab(\"Complexity\") +
   ylab(\"log10 Unique Reads\") +";
if (defined $opt{'A'}) {
print R "
	theme(legend.background=element_blank(),legend.title=element_blank())";
} else {
print R "
	theme(legend.position=\"none\")";
}
print R "
ggsave(plot=PLT,filename=\"$opt{'O'}.png\",width=7,height=6)
ggsave(plot=PLT,filename=\"$opt{'O'}.pdf\",width=7,height=6)
";

if (defined $opt{'M'}) {
	print R "
IN_sub=subset(IN,V4<100&V4>0)

#take unique aligned 
IN_sub\$unique_aligned<-IN_sub\$V3 
x1 <- IN_sub\$unique_aligned[IN_sub\$unique_aligned != 0]

# trimodal fit
km <- kmeans(log10(x1),centers=3)
clustr <- as.factor(km\$cluster)
data_1<-data.frame(\"val\"=log10(x1),\"km\"=clustr)
highest_cluster<-which.max(km\$centers)
second<-km\$center[-highest_cluster,1]
second_highest_cluster<-as.numeric(names(which.max(second)))
cluster_top<-data_1\$val[which(data_1\$km==highest_cluster)]
border1=mean(max(data_1\$val[which(data_1\$km==second_highest_cluster)]),min(data_1\$val[which(data_1\$km==highest_cluster)]))

#in case the groups seperate we do normal fit
sink(\"$opt{'O'}.log\")
fit <- fitdistr(cluster_top, \"normal\")
para <- fit\$estimate
threshold_final_1=para[1]-1.96*para[2]

#in case 3 groups do not sepearte well do a mixed model instead of normal
MMS<-normalmixEM(log10(cluster_top),maxrestarts=1000,maxit = 10000)
threshold_final_2 =10**(MMS\$mu[2]-1.96*MMS\$sigma[2])

p1<-ggplot(data_1, aes(x=data_1\$val)) + geom_histogram(aes(fill=data_1\$km,y=..count../sum(..count..)),color=\"grey50\")+ylab(\"Density\")+stat_density(geom=\"line\",color=\"red\")+geom_vline(xintercept = threshold_final_1)+xlab(\"Log10 Unique Reads\")+theme_bw()+theme(legend.background=element_blank(),text = element_text(size=10),legend.title=element_blank())
ggsave(p1,filename = \"$opt{'O'}.dist.threshold_normal.pdf\")
ggsave(p1,filename = \"$opt{'O'}.dist.threshold_normal.png\")

p2<-ggplot(data_1, aes(x=data_1\$val)) + geom_histogram(aes(fill=data_1\$km,y=..count../sum(..count..)),color=\"grey50\")+ylab(\"Density\")+stat_density(geom=\"line\",color=\"red\")+geom_vline(xintercept = threshold_final_2)+xlab(\"Log10 Unique Reads\")+theme_bw()+theme(legend.background=element_blank(),text = element_text(size=10),legend.title=element_blank())
ggsave(p2,filename = \"$opt{'O'}.dist.threshold_mixed.pdf\")
ggsave(p2,filename = \"$opt{'O'}.dist.threshold_mixed.png\")

print(paste(\"the number of total cells: \",length(IN_sub\$unique_aligned)))
print(paste(\"threshold 1: \",10**threshold_final_1))
print(paste(\"the number of cells above this threshold 1: \",length(which(IN_sub\$unique_aligned>round(10**threshold_final_1)))))

print(paste(\"threshold 2: \",10**threshold_final_2))
print(paste(\"the number of cells above this threshold 2: \",length(which(IN_sub\$unique_aligned>round(10**threshold_final_2)))))
sink()
";
}

print R "
PLT<-ggplot(data=subset(IN,V4<100&V4>0)) + theme_bw() +
	geom_histogram(aes(log10(V3),fill=V2)) +";
if (defined $opt{'C'} || defined $opt{'c'}) {
	print R "   scale_colour_manual(values = c($color_mapping)) +";
}
print R "
	xlab(\"log10 Unique Reads\") +
	ylab(\"Counts\") +
	scale_x_continuous(limits=c(0,6)) +";
if (defined $opt{'A'}) {
print R "
	theme(legend.background=element_blank(),legend.title=element_blank())";
} else {
print R "
	theme(legend.position=\"none\")";
}
print R "
ggsave(plot=PLT,filename=\"$opt{'O'}.hist.png\",width=7,height=6)
ggsave(plot=PLT,filename=\"$opt{'O'}.hist.pdf\",width=7,height=6)
";

close R;

system("$Rscript $opt{'O'}.plot.r");

}
1;
