package sci_commands::cicero_plot;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("cicero_plot");

sub cicero_plot {

@ARGV = @_;
# Defaults

getopts("O:R:XG:f:c:", \%opt);

$die2 = "
scitools cicero-plot [options] [directory containing cds files][bed or CCAN file or list of windows or list of genes]

Prior to running, ensure you have ran scitools matrix-makecds and cds_cicero. 

ARGV0 = cicero_analysis directory
  should have cicero.CCANS.txt and cicero.output.txt
ARGV1 = list of windows (bed)
  or:  $2== file
  or: chrN_start_end,chrN_start_end,etc...
     OR list of gene symbols (file with list)
  or: symbol1,symbol2,etc...
      will add $flank kbp in either direction
      of transcript length
  
Options:
   -O    [STR]   Output Directory (default is [current working directory]/cicero_output)
   -R    [STR]   Rscript call (def = $Rscript)
   -G    [STR]   Genome to be used. Must be of list: hg19,mm10 hg38 has to be made (Default: hg38)
   -f    [INT]   flacking regions around plotted gene regions provided. (Default: 500000)
   -c    [NUM]   correlation cutoff  (Default: 0.15)
   -X            Retain intermediate files (Default = delete)

Note: currently mm10 and hg19, and hg38
                  
";


if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = "$ARGV[0]/cicero_plots"};
if (-e "$opt{'O'}") {
   print STDERR "\n$opt{'O'} directory found - will output here.\n";
} else {
   print STDERR "\n$opt{'O'} directory not found, creating it.\n";
   system("mkdir $opt{'O'}");
}
$ARGV[0] =~ s/\/$//;
$opt{'O'} =~ s/\/$//;
if (!defined $opt{'G'}) {$opt{'G'} = "hg38"};
if (!defined $opt{'f'}) {$flank = "500000"}else {$flank = $opt{'f'}};
if (!defined $opt{'c'}) {$corr_cutoff = 0.15} else {$corr_cutoff = $opt{'c'}};
if ($opt{'G'} eq "mm10") {$gene_annot = "/home/groups/oroaklab/refs/mm10/mm10.ensembl.cicero.annotations.txt"}elsif($opt{'G'} eq "hg19"){$gene_annot = "/home/groups/oroaklab/refs/hg19/hg19.ensembl.cicero.annotations.txt"}elsif($opt{'G'} eq "hg38"){$gene_annot = "/home/groups/oroaklab/refs/hg38/hg38.ensembl.cicero.annotations.txt"} else {die $die2};


$bed = -1;
$ccan_file = -1;
if (-e "$ARGV[1]") {
   #first determine if it is file
   print STDERR "$ARGV[1] determined to be a file.\n";
   open IN, "$ARGV[1]";
   while ($l = <IN>) {
      chomp $l;
      if ($l =~ /\t/) {
         if ($bed < 0) {print STDERR "$ARGV[1] determined to be a bed or ccan file.\n"}
         $bed = 1;
         @P = split(/\t/, $l);
         if (@P == 2) {
            if ($ccan_file < 0) {print STDERR "$ARGV[1] determined to be a ccan file.\n"}
            $ccan_file = 1;
            ($chr,$start,$end) = split(/_/, $P[0]);
            if (defined $CCAN_chr{$P[1]}) {
               if ($start < $CCAN_start{$P[1]}) {$CCAN_start{$P[1]} = $start}
               if ($end > $CCAN_end{$P[1]}) {$CCAN_end{$P[1]} = $end}
            } else {
               $CCAN_chr{$P[1]} = $chr;
               $CCAN_start{$P[1]} = $start;
               $CCAN_end{$P[1]} = $end;
            }
         } elsif (@P > 2) {
            $chr = $P[0]; $start = $P[1]; $end = $P[2];
            $winID = "$chr\_$start\_$end";
            $WINID_chr{$winID} = $chr;
            $WINID_start{$winID} = $start;
            $WINID_end{$winID} = $end;
         }
      } else {
         $bed = 0;
         if ($bed < 0) {print STDERR "$ARGV[1] determined to be a gene symbol list file.\n"}
         $GENE_chr{$l} = "null";
      }
   } close IN;
#not reading in file but is given a list 
} else {
   print STDERR "$ARGV[1] determined to be a list.\n";
   @WINS = split(/,/, $ARGV[1]);
   if ($WINS[0] =~ /^chr/) {
      $bed = 1;
      print STDERR "$ARGV[1] determined to be bed positions.\n";
      foreach $winID (@WINS) {
         #define windows and swich bed 
         ($chr,$start,$end) = split(/_/, $winID);
         $WINID_chr{$winID} = $chr;
         $WINID_start{$winID} = $start;
         $WINID_end{$winID} = $end;
      }
   } else {
      #in case of gene symbols save gene names
      $bed = 0;
      print STDERR "$ARGV[1] determined to be gene symbols.\n";
      foreach $winID (@WINS) {
         $GENE_chr{$winID} = "null";
      }
   }
}

#when known that ccan file was used define windows based on that
if ($ccan_file > 0) {
   foreach $ccan (keys %CCAN_chr) {
      $winID = "$CCAN_chr{$ccan}_$CCAN_start{$ccan}_$CCAN_end{$ccan}";
      $WINID_chr{$winID} = $CCAN_chr{$ccan};
      $WINID_start{$winID} = $CCAN_start{$ccan};
      $WINID_end{$winID} = $CCAN_end{$ccan};
      $WINID_CCAN{$winID}=$ccan;

   }
}

#if not CCAN or bed file then look up gene coordinates
if ($bed == 0) {
   open IN, "$gene_annot";
   $null = <IN>;
   while ($l = <IN>) {
      chomp $l;
      @P = split(/\t/, $l);
      if (defined $GENE_chr{$P[8]}) {
         if (!defined $GENE_start{$P[8]}) {$GENE_start{$P[8]} = $P[2];$GENE_end{$P[8]} = $P[3]} else {
         if ($P[2] < $GENE_start{$P[8]}) {$GENE_start{$P[8]} = $P[2]};
         if ($P[3] > $GENE_end{$P[8]}) {$GENE_end{$P[8]} = $P[3]};
         }
         $GENE_chr{$P[8]} = $P[1];
         $included++;
      }
   } close IN;
   
   if ($included < 1) {die "\nFATAL ERROR: No gene symbols in provided list were found in $gene_annot file.\n"};
   
   #now define windows based on the provided gene list
   foreach $gene (keys %GENE_chr) {
      if ($GENE_chr{$gene} eq "null") {
         print STDERR "WARNING: Cannot find $gene in $gene_annot file!\n";
      } else {
         $tempstart=$GENE_start{$gene} - $flank;
         if ($tempstart < 1) {$tempstart = 1};
          $tempend=$GENE_end{$gene} + $flank;
         $winID = "$GENE_chr{$gene}\_$tempstart\_$tempend";
         $WINID_chr{$winID}=$GENE_chr{$gene};
         $WINID_start{$winID} = $tempstart;
         $WINID_end{$winID} = $tempend;
         $WINID_gene{$winID}=$gene;
      }
   }
}

# make sure hits are present within plot windows

open BED, ">$opt{'O'}/site_windows.tmp";
foreach $winID (sort keys %WINID_chr) {
   $chr = $WINID_chr{$winID};
   $start = $WINID_start{$winID};
   $end = $WINID_end{$winID};
   if (($chr !~ /(Y|M|L|Un)/))
   {
      if ($start !~ /J/)
      {
   print BED "$chr\t$start\t$end\n";   
      }  
   }
} close BED;

system("sort -k 1,1 -k2,2n $opt{'O'}/site_windows.tmp > $opt{'O'}/site_windows.bed");

open COR, ">$opt{'O'}/corr_windows.tmp";
#get cons windows that are above the threshold
open IN, "$ARGV[0]/cicero.output.txt";
   $h = <IN>;
   while ($l = <IN>) {
      chomp $l;
      @P = split(/\t/, $l);
      if ($P[2] > $corr_cutoff) {
         ($chr1,$start1,$end1) = split(/_/, $P[0]);
         ($chr2,$start2,$end2) = split(/_/, $P[1]);
         if ($start1<$start2) {$start = $start1} else {$start = $start2};
         if ($end1>$end2) {$end = $end1} else {$end = $end2};
         if ($chr !~ /(Y|M|L|Un)/)
         {
            if ($start !~ /J/)
            {
            print COR "$chr\t$start\t$end\n";
            }
         }
      }
   } close IN;
close COR;

system("sort -k 1,1 -k2,2n $opt{'O'}/corr_windows.tmp > $opt{'O'}/corr_windows.bed");

open INT, "intersectBed -a $opt{'O'}/corr_windows.bed -b $opt{'O'}/site_windows.bed -wa -wb |";
while ($l = <INT>) {
   chomp $l;
   @P = split(/\t/, $l);
   $winID = "$P[3]_$P[4]_$P[5]";
   $INCLUDE{$winID}++;
} close INT;

$total = 0; $passing = 0;
foreach $winID (keys %WINID_chr) {
   $total++;
   if ((defined $INCLUDE{$winID}) ) {
      $passing++;
   }
}

print STDERR "Of $total window views, $passing have at least one cicero link meeting $corr_cutoff corr cutoff.\n";

if ($passing > 0) {

$timeID = time;
$ts = localtime(time);

open R, ">$ARGV[0]/plot_scriptID_$timeID.r";
print R "
# $ts
# Cicero plots using:
# $ARGV[0]
# $ARGV[1]
# $gene_annot

library(monocle)
library(cicero)

gene_annotations <- as.data.frame(read.table(\"$gene_annot\"))
  gene_annotations\$chromosome <- as.factor(gene_annotations\$chromosome)
  gene_annotations\$start <- as.integer(gene_annotations\$start)
  gene_annotations\$end <- as.integer(gene_annotations\$end)
  gene_annotations\$strand <- as.character(gene_annotations\$strand)
  gene_annotations\$feature <- as.character(gene_annotations\$feature)
  gene_annotations\$gene <- as.character(gene_annotations\$gene)
  gene_annotations\$transcript <- as.character(gene_annotations\$transcript)
  gene_annotations\$symbol <- as.character(gene_annotations\$symbol)
# Loading cons file for chromosome $chr.
cons <- read.table(\"$ARGV[0]/cicero.output.txt\",header=T)\n
# Plotting window $winID



";

foreach $winID (sort keys %INCLUDE) {
   $chr = $WINID_chr{$winID};
   $start = $WINID_start{$winID};
   $end = $WINID_end{$winID};
      if (defined $WINID_gene{$winID})
      {
      print R "
      pdf(\"$opt{'O'}/$WINID_gene{$winID}_$winID\_cutoff\_$corr_cutoff.pdf\",width=24,height=12)
      try(plot_connections(cons,\"$chr\", $start, $end, gene_model = gene_annotations, coaccess_cutoff = .15, connection_width = .5, collapseTranscripts = \"longest\" ))
      dev.off()\n";


         } elsif (defined $WINID_CCAN{$winID}) {
      print R "
      pdf(\"$opt{'O'}/CCAN_$WINID_CCAN{$winID}_$winID\_cutoff\_$corr_cutoff.pdf\",width=24,height=12)
      try(plot_connections(cons,\"$chr\", $start, $end, gene_model = gene_annotations, coaccess_cutoff = .15, connection_width = .5, collapseTranscripts = \"longest\" ))
      dev.off()\n";
         } else {

      print R "
      pdf(\"$opt{'O'}/$winID\_cutoff\_$corr_cutoff.pdf\",width=24,height=12)
      try(plot_connections(cons,\"$chr\", $start, $end, gene_model = gene_annotations, coaccess_cutoff = .15, connection_width = .5, collapseTranscripts = \"longest\" ))
      dev.off()\n";
      }
      }

close R;
system("Rscript $ARGV[0]/plot_scriptID_$timeID.r");
if (!defined $opt{'X'}) {
    system("rm -f $ARGV[0]/plot_scriptID_$timeID.r");
    system("rm -f $opt{'O'}/*.tmp");
}
}
}
1;
