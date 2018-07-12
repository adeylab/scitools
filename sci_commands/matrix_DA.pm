package sci_commands::matrix_DA;


use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_DA");

sub matrix_DA {

@ARGV = @_;
use Getopt::Std; %opt = ();
getopts("O:A:I:", \%opt);

$die2 = "
scitools matrix-da [options] [counts matrix] [aggr annotation file]
   or    da-matrix

This script will perform DA analysis on aggregate matrix. Aggregate annotation file that is output by the aggregate_cells

Options:
   -O   [STR]   Output prefix (default is [input annot].matrix)
   -A   [STR]   provided annotation file which consist of aggragate_centroid_name\tgroupthat you want to compare
		e.g.: IND1_1	IND1
		      IND1_2	IND1
		      IND2_3	IND2
		Warning: This will be used for comparisons instead of agg annot file
   -I	[STR]   If defined script compares an individual group to all others combined as opposed to comparing group by group (default) 	

   
";

#name output and create folder 
if (!defined $ARGV[1]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix$//};


$name_out = "DA_plots";
	
if (-e "$opt{'O'}.$name_out") {
	die "\nFATAL: $opt{'O'}.$name_out directory already exists! Exiting!\n$die2";
}

system("mkdir $opt{'O'}.$name_out");
system("mkdir $opt{'O'}.$name_out/plots");
	
#read in annotation, scitools approach
if (defined $opt{'A'}) 
	{
	read_annot($opt{'A'});
	for my $CELLID (sort keys %CELLID_annot)
	{
		push(@{$ANNOT_AGGID{$CELLID_annot{$CELLID}}},$CELLID);
	}
	
	}
elsif (!defined $opt{'A'} && defined $ARGV[1]) {
	read_annot($ARGV[1]);
	for my $aggannot (sort keys %ANNOT_count)
	{
		@annotagg = split(/_/, $aggannot);
		push(@{$ANNOT_AGGID{$annotagg[0]}},$aggannot);
	}
	} else {die $die2}



#read in matrix for basic stats, scitools standard
read_matrix_stats($ARGV[0]);


#create contrast annot

#contrast of individual groups against all other groups together

if (defined $opt{'I'})
{
print "Doing all vs ind comparison\n";
for my $group1 (sort keys %ANNOT_AGGID)
	{
               	$contrast="$group1\_vs_all_as_ref";
               	$contrast_hash{$contrast}++;
                open OUT, "> $opt{'O'}.$name_out/Diff_acc_$contrast.annot"; 
                for my $group2 (sort keys %ANNOT_AGGID)
		{
			if ($group1 ne $group2)
			{
				for my $AGGID (@{$ANNOT_AGGID{$group1}})
				{
               			print OUT $AGGID. "\t" .$group1. "\n";
				}
                        
				for my $AGGID (@{$ANNOT_AGGID{$group2}})
				{
               			print OUT $AGGID. "\t" .$contrast."\n";
				}
			}	  
		}
	 close(OUT);
         }
        
 } else {
print "Doing ind vs ind comparison\n";
#contrast of individual groups against other individual groups
for my $group1 (sort keys %ANNOT_AGGID)
	{
		for my $group2 (sort keys %ANNOT_AGGID)
		{
		$contrast="$group1\_vs_$group2\_as_ref";
               	if (($group1 ne $group2) && (!defined $contrast_hash{$contrast}))
               	{     	
               	 
                 $contrast_hash{$contrast}++;
                        open OUT, "> $opt{'O'}.$name_out/Diff_acc_$contrast.annot";   
                        for my $AGGID (@{$ANNOT_AGGID{$group1}})
               		{
               			print OUT $AGGID. "\t" .$group1. "\n";
               		}
                        
               		for my $AGGID (@{$ANNOT_AGGID{$group2}})
               		{
               			print OUT $AGGID. "\t" .$contrast."\n";
               		}
               	}
               	close(OUT);
		}
            

	}
}
for $contrast (sort keys %contrast_hash)
	{
print "Analyzing : ". $contrast."\n";        
open R, ">$opt{'O'}.$name_out/Diff_acc_$contrast.r";
print R "
#library(\"reshape2\")
library(\"ggplot2\")
library(\"DESeq2\")
library(\"BiocParallel\")
#library(calibrate)
library(qvalue)
register(MulticoreParam(30)) 
counts_mat<-as.matrix(read.delim(\"$ARGV[0]\"))
counts_mat<-na.omit(counts_mat)
coldata<-read.table(file = \"$opt{'O'}.$name_out/Diff_acc_$contrast.annot\",sep = \"\\t\",row.names = 1)
colnames(coldata)<-c(\"condition\")
counts_mat <- counts_mat[, rownames(coldata)]
all(rownames(coldata) == colnames(counts_mat))
dds <- DESeqDataSetFromMatrix(countData = counts_mat,
                              colData = coldata,
                              design = ~ condition)
#what you compare against
dds\$condition <- relevel(dds\$condition, ref = \"$contrast\")
dds <- DESeq(dds,parallel = TRUE)
res <- results(dds)
write.table(as.matrix(res),file = \"$opt{'O'}.$name_out/Differential_acc_$contrast.txt\", col.names = TRUE, row.names = TRUE, sep = \"\t\", quote = FALSE)

res <- lfcShrink(dds, coef=2)
write.table(as.matrix(res),file = \"$opt{'O'}.$name_out/Differential_acc_$contrast\_shrunk.txt\", col.names = TRUE, row.names = TRUE, sep = \"\t\", quote = FALSE)


shrunk_corr<-data.frame(log2FoldChange=res\$log2FoldChange,pvalue=res\$pvalue,padj=res\$padj)
row.names(shrunk_corr)<-row.names(res)
qval<-qvalue(shrunk_corr\$pvalue,fdr.level=0.01)
shrunk_corr\$qval<-qval\$qvalues
#write.table(as.matrix(shrunk_corr),file = \"Diff_acc_shrunk_$contrast\_combined.txt\", col.names = TRUE, row.names = TRUE, sep = \"\t\", quote = FALSE)
df<-shrunk_corr


output<-data.frame(\"annotation\"=row.names(b),\"pval\"=df\$pvalue,\"pval_adjust\"=df\$padj,\"log2fold\"= df\$log2FoldChange,\"qval\"= df\$qval)

##Highlight genes that have an absolute fold change > 2 and a p-value < Bonferroni cut-off
output\$threshold = as.factor(output\$log2fold > 1 & output\$qval < 0.05)


write.table(as.matrix(output),file = \"$opt{'O'}.$name_out/Differential_acc_$contrast\_shrunk_q001.txt\", col.names = TRUE, row.names = FALSE, sep = \"\t\", quote = FALSE)

##Construct the plot object
g <- ggplot(data=output, aes(x=log2fold, y=-log10(pval), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  xlab(\"log2 fold change\") + ylab(\"-log10 p-value\")+theme_bw()
ggsave(plot = g,filename = \"$opt{'O'}.$name_out/plots/Differential_acc_$contrast\_shrunk_qval_001_threshold_plotpval.png\")
ggsave(plot = g,filename = \"$opt{'O'}.$name_out/plots/Differential_acc_$contrast\_shrunk_qval_001_threshold_plotpval.pdf\")

##Construct the plot object
g <- ggplot(data=output, aes(x=log2fold, y=qval, colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  xlab(\"log2 fold change\") + ylab(\"q-value\")+theme_bw()
ggsave(plot = g,filename = \"$opt{'O'}.$name_out/plots/Differential_acc_$contrast\_shrunk_qval_001_threshold_plotqval.png\")
ggsave(plot = g,filename = \"$opt{'O'}.$name_out/plots/Differential_acc_$contrast\_shrunk_qval_001_threshold_plotqval.pdf\")



qval<-qvalue(shrunk_corr\$pvalue,fdr.level=0.2)
shrunk_corr\$qval<-qval\$qvalues
res<-shrunk_corr
df<-as.data.frame(res)
output<-data.frame(\"annotation\"=row.names(b),\"pval\"=df\$pvalue,\"pval_adjust\"=df\$padj,\"log2fold\"= df\$log2FoldChange,\"qval\"= df\$qval)

output\$threshold = as.factor(output\$log2fold > 1 & output\$qval < 0.05)
write.table(as.matrix(output),file = \"$opt{'O'}.$name_out/Differential_acc_$contrast\_shrunk_q02.txt\", col.names = TRUE, row.names = FALSE, sep = \"\t\", quote = FALSE)

         
";

close(R);
system("Rscript $opt{'O'}.$name_out/Diff_acc_$contrast.r > $opt{'O'}.$name_out/Diff_acc_$contrast.stdout 2> $opt{'O'}.$name_out/Diff_acc_$contrast.stderr");	
system("rm -f $opt{'O'}.$name_out/Diff_acc_$contrast.r $opt{'O'}.$name_out/Diff_acc_$contrast.stdout $opt{'O'}.$name_out/Diff_acc_$contrast.stderr $opt{'O'}.$name_out/Differential_acc_$contrast.txt $opt{'O'}.$name_out/Differential_acc_$contrast\_shrunk.txt");	
}

#from here we are looking at peaks that are specifially differentially accessible only in that contrast

#compare 0.01 FDR contrast peaks to 0.2 FDR of other groups 
# it is only worth it to do this in the one vs all comparisons
if (defined $opt{'I'})
{
print "Doing all vs ind comparison FDR 0.2 peaks removal\n";
for my $contrast1 (sort keys %contrast_hash)
	{
	my ($group1a,$group1b) = split(/_vs_/, $contrast1);
	
        my %signf_peaks=();
        open IN, "$opt{'O'}.$name_out/Differential_acc_$contrast\_shrunk_q001.txt"; 
        my $firstline = <IN>;
                        while (my $l = <IN>) 
                        {
                          chomp $l;
                          my @P = split(/\t/, $l);
                          if ($P[5] eq "TRUE")
                          {
                          $P[0]=~ s/\s//g;
                          $signf_peaks{$P[0]}=$l;
                          }

                         }
                          close(IN);

	for my $contrast2 (sort keys %contrast_hash)
               {
                my ($group2a,$group2b) = split(/_vs_/, $contrast2);           
                  if (($group1a eq $group2a) && ($group1b ne $group2b))
                        {
                  open IN2, "$opt{'O'}.$name_out/Differential_acc_$contrast\_shrunk_q02.txt"; 
                  $firstline = <IN2>;
                  while (my $l = <IN2>) 
                  {
                  chomp $l;
                  my @P = split(/\t/, $l);
                  $P[0]=~ s/\s//g;	
                      if ($P[5] eq "TRUE")
                      {
                      	if(exists $signf_peaks{$P[0]})
                      	{
                             #remove peaks significant in other pops    
                        delete $signf_peaks{$P[0]};
                        }
                      }
                  }
            	 close(IN2);
                      }
                }
        
        
        open (OUT,">","$opt{'O'}.$name_out/Diff_acc_shrunk_$contrast1\_filtered_final.txt");
        open (OUT2,">","$opt{'O'}.$name_out/Diff_acc_shrunk_$contrast1\_filtered_final_just_good_peaks.txt");
        open (OUT3,">","$opt{'O'}.$name_out/Diff_acc_shrunk_$contrast1\_not_filtered_final_just_good_peaks.txt");
	open IN3, "$opt{'O'}.$name_out/Differential_acc_$contrast\_shrunk_q001.txt"; 
	$firstline= <IN3>;
	chomp $firstline;
	$firstline =~ s/annotation\.//g;
	print OUT $firstline. "\tFinal_filter_pass\n";
        while (my $l = <IN3>) 
	{
          chomp $l;

          my @P = split(/\t/, $l);
          $P[0]=~ s/\s//g;	
          if (exists $signf_peaks{$P[0]})
            {
              print OUT $l. "\tTRUE\n";
              print OUT2 $l. "\n";
            }
            else
            {
              print OUT $l. "\tFALSE\n";
            }
	
          if ($P[5] eq "TRUE")
            {
              print OUT3 $l. "\n";
            }
	}
	close(IN3);
	close(OUT);
	close(OUT2);
        close(OUT3);

        open R, ">$opt{'O'}.$name_out/Diff_acc_$contrast1.r";
        print R "
        library(\"ggplot2\")
        a<-read.delim(\"./Diff_acc_shrunk_$contrast1\_filtered_final.txt\")
        
        ##Construct the plot object
        g <- ggplot(data=a, aes(x=log2fold, y=-log10(pval), colour=Final_filter_pass)) +
        geom_point(alpha=0.4, size=1.75) +
        xlab(\"log2 fold change\") + ylab(\"-log10 p-value\")+theme_bw()
        ggsave(plot = g,filename = \"$opt{'O'}.$name_out/plots/Differential_access_$contrast1\_q001_q02_removed.png\")
        ggsave(plot = g,filename = \"$opt{'O'}.$name_out/plots/Differential_access_$contrast1\_q001_q02_removed.pdf\")
        ";
        close(R);
        system("Rscript $opt{'O'}.$name_out/Diff_acc_$contrast1.r > $opt{'O'}.$name_out/Diff_acc_$contrast1.stdout 2> $opt{'O'}.$name_out/Diff_acc_$contrast1.stderr");	
        system("rm -f $opt{'O'}.$name_out/Diff_acc_$contrast1.r $opt{'O'}.$name_out/Diff_acc_$contrast1.stdout $opt{'O'}.$name_out/Diff_acc_$contrast1.stderr");	
        }
}

}
1;
