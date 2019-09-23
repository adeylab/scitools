package sci_commands::matrix_DA;


use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_DA");

sub matrix_DA {

@ARGV = @_;
use Getopt::Std; %opt = ();
#190425 RM Correction: Some option flags listed by unused.
getopts("O:A:IT:Xn:", \%opt);

$die2 = "
scitools matrix-da [options] [aggregate matrix] [aggregate annotation file]
   or    da-matrix

This script will perform differential accessibility analysis on:
[aggregate matrix]:           an aggregate matrix generated through scitools aggregate-matrix
[aggregate annotation file]:  an aggregate annotation file that is output by scitools aggregate_cells

With the statistical test that is specified by option -T.

Options:
   -O   [STR]     Output prefix (default is [input annot].matrix)
   -A   [STR]     provided annotation file which consist of aggregate_centroid_name\tgroupthat you want to compare
                           e.g.: IND1_1   IND1
                                 IND1_2   IND1
                                 IND2_3   IND2
                           Warning: This will be used for comparisons instead of [aggregate annotation file]
   -T   [STR]     Currently accepts [Wald|LRT]
                  Type of test performed: negative binomial \"Wald\" or Likelihood ratio test (\"LRT\"). binomialff test will be added later default: \"Wald\" 
   -I   [FLAG]    If defined script compares an individual group to all others combined.
                  (Default: If not flagged performs all group by group comparisons.)  
   -n   [INT]     Number of cores to use for comparisons. (Default=1)
   -X   [FLAG]    Retain intermediate files (def = delete)
 
";

#name output and create folder 
if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix$//};
if (!defined $opt{'T'}) {$opt{'T'} = "Wald"};
if (!defined $opt{'n'}) {$opt{'n'} = 1};

#now creates separate directory for Wald, binomialff and LRT

if ($opt{'T'} eq "Wald")
  {
   $name_out = "DA_Wald"; 
  }
elsif ($opt{'T'} eq "LRT")
  {
    $name_out = "DA_LRT";
  }
elsif ($opt{'T'} eq "binomialff")
  {
    $name_out = "DA_binomialff";
  }
  else {print "This is not a method defined"; die}

#RM: Just overwrite folder contents if they already exist
#if (-e "$opt{'O'}.$name_out") {
#	die "\nFATAL: $opt{'O'}.$name_out directory already exists! Exiting!\n$die2";
#}

system("mkdir $opt{'O'}.$name_out");
system("mkdir $opt{'O'}.$name_out/plots");
	
#read in annotation, scitools approach
if (defined $opt{'A'}) {
   read_annot($opt{'A'});
   for my $CELLID (sort keys %CELLID_annot){
   push(@{$ANNOT_AGGID{$CELLID_annot{$CELLID}}},$CELLID);
   }
  } elsif (!defined $opt{'A'} && defined $ARGV[1]) {
   read_annot($ARGV[1]);
   for my $aggannot (sort keys %ANNOT_count){
      @annotagg = split(/_/, $aggannot);
	  pop(@annotagg);
	  $aggannotID = join("_", @annotagg);
      push(@{$ANNOT_AGGID{$annotagg[0]}},$aggannot);
   }
   } else {die $die2}



#read in matrix for basic stats, scitools standard
read_matrix_stats($ARGV[0]);


#create contrast annot
#contrast of individual groups against all other groups together

if (defined $opt{'I'}){
  print "Doing all vs ind comparison\n";
  for my $group1 (sort keys %ANNOT_AGGID){
      $contrast="$group1\_vs_all_as_ref";
      $contrast_hash{$contrast}++;
    open OUT, "> $opt{'O'}.$name_out/Diff_acc_$contrast.annot"; 
    for my $group2 (sort keys %ANNOT_AGGID){
      if ($group1 ne $group2){
         for my $AGGID (@{$ANNOT_AGGID{$group1}}){
          print OUT $AGGID. "\t" .$group1. "\n";
           }
        for my $AGGID (@{$ANNOT_AGGID{$group2}}){
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
  library(ggplot2)
  library(DESeq2)
  library(BiocParallel)
  library(qvalue)
  register(MulticoreParam($opt{'n'})) 
  counts_mat<-as.matrix(read.delim(\"$ARGV[0]\"))
  counts_mat<-na.omit(counts_mat)
  coldata<-read.table(file = \"$opt{'O'}.$name_out/Diff_acc_$contrast.annot\",sep = \"\\t\")
  colnames(coldata)<-c(\"aggID\",\"condition\")

  coldata<-unique(coldata[1:2])
  rownames(coldata)<-coldata\$aggID

  counts_mat <- counts_mat[, rownames(coldata)]
  all(rownames(coldata) == colnames(counts_mat))

  dds <- DESeqDataSetFromMatrix(countData = counts_mat,
                              colData = coldata,
                              design = ~ condition)
  #what you compare against
  dds\$condition <- relevel(dds\$condition, ref = \"$contrast\")
  ";
  if ($opt{'T'} eq "Wald")
  {
  print "This is Wald\n";
  print R "
  dds <- DESeq(dds,parallel = TRUE)
  res <- results(dds)
  write.table(as.matrix(res),file = \"$opt{'O'}.$name_out/Differential_acc_$contrast\_wald.txt\", col.names = TRUE, row.names = TRUE, sep = \"\\t\", quote = FALSE)

  res <- lfcShrink(dds, coef=2)
  write.table(as.matrix(res),file = \"$opt{'O'}.$name_out/Differential_acc_$contrast\_shrunk_wald.txt\", col.names = TRUE, row.names = TRUE, sep = \"\\t\", quote = FALSE)

  shrunk_corr<-data.frame(log2FoldChange=res\$log2FoldChange,pvalue=res\$pvalue,padj=res\$padj)
  row.names(shrunk_corr)<-row.names(res)
  qval<-qvalue(shrunk_corr\$pvalue,fdr.level=0.01)
  shrunk_corr\$qval<-qval\$qvalues
  # write.table(as.matrix(shrunk_corr),file = \"Diff_acc_shrunk_$contrast\_combined_wald.txt\", col.names = TRUE, row.names = TRUE, sep = \"\\t\", quote = FALSE)
  df<-shrunk_corr


  output<-data.frame(\"annotation\"=row.names(df),\"pval\"=df\$pvalue,\"pval_adjust\"=df\$padj,\"log2fold\"= df\$log2FoldChange,\"qval\"= df\$qval)

  ##Highlight genes that have an absolute fold change > 2 and a p-value < Bonferroni cut-off
  output\$threshold = as.factor(output\$log2fold > 1 & output\$qval < 0.05)

  write.table(as.matrix(output),file = \"$opt{'O'}.$name_out/Differential_acc_$contrast\_shrunk_q001_wald.txt\", col.names = TRUE, row.names = FALSE, sep = \"\\t\", quote = FALSE)

  ##Construct the plot object
  g <- ggplot(data=output, aes(x=log2fold, y=-log10(output\$pval), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  xlab(\"log2 fold change\") + ylab(\"-log10 p-value\")+theme_bw()
  ggsave(plot = g,filename = \"$opt{'O'}.$name_out/plots/Differential_acc_$contrast\_shrunk_qval_001_threshold_plotpval_wald.png\")
  ggsave(plot = g,filename = \"$opt{'O'}.$name_out/plots/Differential_acc_$contrast\_shrunk_qval_001_threshold_plotpval_wald.pdf\")

  ##Construct the plot object
  g <- ggplot(data=output, aes(x=log2fold, y=output\$qval, colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  xlab(\"log2 fold change\") + ylab(\"q-value\")+theme_bw()
  ggsave(plot = g,filename = \"$opt{'O'}.$name_out/plots/Differential_acc_$contrast\_shrunk_qval_001_threshold_plotqval_wald.png\")
  ggsave(plot = g,filename = \"$opt{'O'}.$name_out/plots/Differential_acc_$contrast\_shrunk_qval_001_threshold_plotqval_wald.pdf\")



  qval<-qvalue(shrunk_corr\$pvalue,fdr.level=0.2)
  shrunk_corr\$qval<-qval\$qvalues
  res<-shrunk_corr
  df<-as.data.frame(res)
  output<-data.frame(\"annotation\"=row.names(df),\"pval\"=df\$pvalue,\"pval_adjust\"=df\$padj,\"log2fold\"= df\$log2FoldChange,\"qval\"= df\$qval)

  output\$threshold = as.factor(output\$log2fold > 1 & output\$qval < 0.05)
  write.table(as.matrix(output),file = \"$opt{'O'}.$name_out/Differential_acc_$contrast\_shrunk_q02_wald.txt\", col.names = TRUE, row.names = FALSE, sep = \"\\t\", quote = FALSE)

  ";
  close(R);
  }
  elsif ($opt{'T'} eq "LRT")
  {
  print "This is LRT\n";
  print R "
  # an alternate analysis: likelihood ratio test
  ddsLRT <- DESeq(dds, test=\"LRT\", reduced= ~ 1)
  res <- results(ddsLRT)
  write.table(as.matrix(res),file = \"$opt{'O'}.$name_out/Differential_acc_$contrast\_LRT.txt\", col.names = TRUE, row.names = TRUE, sep = \"\\t\", quote = FALSE)

  df-data.frame(log2FoldChange=res\$log2FoldChange,pvalue=res\$pvalue,padj=res\$padj)
  row.names(df)<-row.names(res)
  qval<-qvalue(shrunk_corr\$pvalue,fdr.level=0.01)
  df\$qval<-qval\$qvalues
  output<-data.frame(\"annotation\"=row.names(df),\"pval\"=df\$pvalue,\"pval_adjust\"=df\$padj,\"log2fold\"= df\$log2FoldChange,\"qval\"= df\$qval)

  ##Highlight genes that have an absolute fold change > 2 and a p-value < Bonferroni cut-off
  output\$threshold = as.factor(output\$log2fold > 1 & output\$qval < 0.05)


  write.table(as.matrix(output),file = \"$opt{'O'}.$name_out/Differential_acc_$contrast\_q001_LRT.txt\", col.names = TRUE, row.names = FALSE, sep = \"\\t\", quote = FALSE)

  ##Construct the plot object
  g <- ggplot(data=output, aes(x=log2fold, y=-log10(output\$pval), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  xlab(\"log2 fold change\") + ylab(\"-log10 p-value\")+theme_bw()
  ggsave(plot = g,filename = \"$opt{'O'}.$name_out/plots/Differential_acc_$contrast\_shrunk_qval_001_threshold_plotpval_LRT.png\")
  ggsave(plot = g,filename = \"$opt{'O'}.$name_out/plots/Differential_acc_$contrast\_shrunk_qval_001_threshold_plotpval_LRT.pdf\")

  ##Construct the plot object
  g <- ggplot(data=output, aes(x=log2fold, y=output\$qval, colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  xlab(\"log2 fold change\") + ylab(\"q-value\")+theme_bw()
  ggsave(plot = g,filename = \"$opt{'O'}.$name_out/plots/Differential_acc_$contrast\_shrunk_qval_001_threshold_plotqval_LRT.png\")
  ggsave(plot = g,filename = \"$opt{'O'}.$name_out/plots/Differential_acc_$contrast\_shrunk_qval_001_threshold_plotqval_LRT.pdf\")



  qval<-qvalue(shrunk_corr\$pvalue,fdr.level=0.2)
  shrunk_corr\$qval<-qval\$qvalues
  res<-shrunk_corr
  df<-as.data.frame(res)
  output<-data.frame(\"annotation\"=row.names(df),\"pval\"=df\$pvalue,\"pval_adjust\"=df\$padj,\"log2fold\"= df\$log2FoldChange,\"qval\"= df\$qval)

  output\$threshold = as.factor(output\$log2fold > 1 & output\$qval < 0.05)
  write.table(as.matrix(output),file = \"$opt{'O'}.$name_out/Differential_acc_$contrast\_q02_LRT.txt\", col.names = TRUE, row.names = FALSE, sep = \"\\t\", quote = FALSE)
  ";
  close(R);
  }
  elsif ($opt{'T'} eq "binomialff")
  {
  print "This is binomialff\n";
  system("scitools matrix-make-cds -O binomialfftemp $ARGV[0] $opt{'O'}.$name_out/Diff_acc_$contrast.annot $opt{'D'}");
  print R "
  suppressWarnings(library(cicero))
  message(\"Loading Cicero/monocle\")
  # reading in matrix, annotation, and dimension data

  cds_cell_data <- read.delim(\"binomialfftemp.cds_files/cds_cell_data.txt\")
  cds_dims_data <- read.delim(\"binomialfftemp.cds_files/cds_dims_data.txt\",header=F)
  cds_site_data <- read.delim(\"binomialfftemp.cds_files/cds_site_data.txt\")
  cds_counts_matrix <- read.table(\"binomialfftemp.cds_files/cds_counts_matrix.txt\")
  message(\"Read in CDS Files.\")

  # generating cell data set
  feature_data <- new(\"AnnotatedDataFrame\", data = cds_site_data)
  sample_data <- new(\"AnnotatedDataFrame\", data = cds_cell_data)

  if(ncol(cds_dims_data)==4){
  colnames(cds_dims_data)<-c(\"cellID\",\"Dimension_1\",\"Dimension_2\",\"Dimension_3\")
  row.names(cds_dims_data)<-cds_dims_data\$cellID
  dimension_reduction<-cds_dims_data[2:ncol(cds_dims_data)]
  }else{
  colnames(cds_dims_data)<-c(\"cellID\",\"Dimension_1\",\"Dimension_2\")
  row.names(cds_dims_data)<-cds_dims_data\$cellID
  dimension_reduction<-cds_dims_data[2:ncol(cds_dims_data)]
  }

  message(\"Setting up CDS matrix, binarized for DA.\")


  cds <- suppressWarnings(newCellDataSet(as.matrix(cds_counts_matrix), phenoData = sample_data, featureData = feature_data,expressionFamily = negbinomial.size(),lowerDetectionLimit = 0))

  cds\@expressionFamily <- binomialff()
  cds\@expressionFamily\@vfamily <- \"binomialff\"

  pData(cds)\$temp <- NULL
  fData(cds)\$chr <- as.numeric(as.character(fData(cds)\$chr))
  fData(cds)\$bp1 <- as.numeric(as.character(fData(cds)\$bp1))
  fData(cds)\$bp2 <- as.numeric(as.character(fData(cds)\$bp2))
  cds <- cds[order(fData(cds)\$chr, fData(cds)\$bp1),]

  set.seed(2017)
  pData(cds)\$cells <- NULL
  cds <- detectGenes(cds)
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)

  #now do the DA
  diff_DA <- differentialGeneTest(agg_cds,fullModelFormulaStr=\"~timepoint + num_genes_expressed\")
  qval001<-qvalue(diff_DA\$pvalue,fdr.level=0.01)
  qval02<-qvalue(diff_DA\$pvalue,fdr.level=0.2)
  #will find out what log2FoldChange is
  output<-data.frame(\"annotation\"=row.names(diff_DA),\"pval\"=diff_DA\$pvalue,\"pval_adjust\"=qval001,\"log2fold\"= qval001,\"qval\"= qval001)
  ##Highlight genes that have an absolute fold change > 2 and a p-value < corrected
  output\$threshold = as.factor(output\$qval < 0.05)
  write.table(as.matrix(output),file = \"$opt{'O'}.$name_out/Differential_acc_$contrast\_q001_binomialff.txt\", col.names = TRUE, row.names = FALSE, sep = \"\\t\", quote = FALSE)
  output<-data.frame(\"annotation\"=row.names(diff_DA),\"pval\"=diff_DA\$pvalue,\"pval_adjust\"=qval02,\"log2fold\"= qval02,\"qval\"= qval02)
  output\$threshold = as.factor(output\$qval < 0.05)
  write.table(as.matrix(output),file = \"$opt{'O'}.$name_out/Differential_acc_$contrast\_q02_binomialff.txt\", col.names = TRUE, row.names = FALSE, sep = \"\\t\", quote = FALSE
  ";
  close(R);
  } else { print "This is not a method defined"; die};

  system("Rscript $opt{'O'}.$name_out/Diff_acc_$contrast.r > $opt{'O'}.$name_out/Diff_acc_$contrast.stdout 2> $opt{'O'}.$name_out/Diff_acc_$contrast.stderr");
if (!defined $opt{'X'}) {
  print "Removing txt and r files\n";
  system("rm -f $opt{'O'}.$name_out/Diff_acc_$contrast.r $opt{'O'}.$name_out/Diff_acc_$contrast.stdout $opt{'O'}.$name_out/Diff_acc_$contrast.stderr"); 
}
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
        if($opt{'T'} eq "Wald"){
        open IN, "$opt{'O'}.$name_out/Differential_acc_$contrast\_shrunk_q001_wald.txt";
        } 
        elsif ($opt{'T'} eq "LRT") {
          open IN, "$opt{'O'}.$name_out/Differential_acc_$contrast\_q001_LRT.txt";
        }
        elsif ($opt{'T'} eq "binomialff") {
          open IN, "$opt{'O'}.$name_out/Differential_acc_$contrast\_q001_binomialff.txt";
        } 
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
                        if($opt{'T'} eq "Wald"){
                          open IN2, "$opt{'O'}.$name_out/Differential_acc_$contrast\_shrunk_q02_wald.txt";
                        } 
                        elsif ($opt{'T'} eq "LRT") {
                         open IN2, "$opt{'O'}.$name_out/Differential_acc_$contrast\_q02_LRT.txt";
                        }
                        elsif ($opt{'T'} eq "binomialff") {
                        open IN2, "$opt{'O'}.$name_out/Differential_acc_$contrast\_q02_binomialff.txt";
                        } 
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
	       if($opt{'T'} eq "Wald"){
        open IN3, "$opt{'O'}.$name_out/Differential_acc_$contrast\_shrunk_q001_wald.txt";
        } 
        elsif ($opt{'T'} eq "LRT") {
          open IN3, "$opt{'O'}.$name_out/Differential_acc_$contrast\_q001_LRT.txt";
        }
        elsif ($opt{'T'} eq "binomialff") {
          open IN3, "$opt{'O'}.$name_out/Differential_acc_$contrast\_q001_binomialff.txt";
        } 
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

        open R, ">$opt{'O'}.$name_out/Diff_acc_$contrast1\_plot.r";
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
        system("Rscript $opt{'O'}.$name_out/Diff_acc_$contrast1\_plot.r > $opt{'O'}.$name_out/Diff_acc_$contrast1\_plot.stdout 2> $opt{'O'}.$name_out/Diff_acc_$contrast1\_plot.stderr");	
		if (!defined $opt{'X'}) {
			system("rm -f $opt{'O'}.$name_out/Diff_acc_$contrast1\_plot.r $opt{'O'}.$name_out/Diff_acc_$contrast1\.r $opt{'O'}.$name_out/Diff_acc_$contrast1.stdout $opt{'O'}.$name_out/Diff_acc_$contrast1.stderr");
		}
        }
}

}
1;
