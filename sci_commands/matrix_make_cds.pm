package sci_commands::matrix_make_cds;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_make_cds");

sub matrix_make_cds {

@ARGV = @_;
# Defaults

getopts("O:R:G:X", \%opt);

$die2 = "
scitools matrix_make_cds (options) [input matrix] [annotation file] [dims file]
   or    matrix_makecds
         matrix2cds
		 
 Generates CDS format files necessary for Monocle3 and Cicero calls. Places formatted files in a subdirectory.
 Generates 4 txt files to be used in scitools atac-monocle3 and scitools atac-cicero calls for quick processing.
 Alternatively, the files are ready to be loaded into monocle or cicero in R, external to scitools. This enables
   use of additional cicero or monocle commands that are not included in the scitools wrapper.
 For more informaiton on monocle and cicero we recommend visiting the trapnell lab github
 
	1. cds_site_data.txt       :   Feature set file (called peaks used in counts matrix)
	2. cds_cell_data.txt       :   Restructured Annotation file, with timepoints included.
	3. cds_dims_data.txt       :   Restructured dims file
	4. cds_counts_matrix.txt   :   Restructured counts matrix file.

	[input matrix]      =   filtered counts matrix
	[annotation file]   =   <.annot> formatted file, generated through scitools make-annot command
	[dims file]         =   <.dims> formated file, generated through scitools matrix-[umap|swne|pca|irlba|tsne] command

Options:
   -O   [STR]   Output Directory (default is [current working directory]/cds_files)
   -G   [bed]   Bed of promters with the name as a gene (will include gene column)
   -R   [STR]   Rscript call (def = $Rscript)
   -b   [STR]   Bedtools call (def = $bedtools)
   -X           Retain intermediate files (for -G)
                  
";


if (!defined $ARGV[0]) {die $die2};
if (!defined $ARGV[1]) {die $die2};
if (!defined $ARGV[2]) {die $die2};

if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix$//};
$opt{'O'} =~ s/\.cds_files$//;

system("mkdir $opt{'O'}.cds_files");

if (defined $opt{'G'}) {
	open IN, "$ARGV[0]";
	open OUT, ">$opt{'O'}.cds_files/cds_sites.bed";
	$null = <IN>;
	while ($l = <IN>) {
		chomp $l;
		@P = split(/\t/, $l);
		$site = shift(@P);
		($chr,$bp1,$bp2) = split(/[_]/, $site);
		$siteName = "$chr\_$bp1\_$bp2";
		$SITEID_gene{$siteName} = "NA";
		print OUT "$chr\t$bp1\t$bp2\t$siteName\n";
	} close IN; close OUT;
	open IN, "$bedtools intersect -a $opt{'O'}.cds_files/cds_sites.bed -b $opt{'G'} -wa -wb |";
	if (defined $opt{'X'}) {open OUT, ">$opt{'O'}.cds_files/cds_site_gene_links.txt"};
	while ($l = <IN>) {
		chomp $l;
		@P = split(/\t/, $l);
		if ($SITEID_gene{$P[3]} eq "NA") {
			$SITEID_gene{$P[3]} = $P[7];
		}
		if (defined $opt{'X'}) {
			print OUT "$P[3]\t$P[7]\n";
		}
	} close IN;
	if (!defined $opt{'X'}) {
		system("rm -f $opt{'O'}.cds_files/cds_sites.bed");
	} else {
		close OUT;
	}
}

open SITE_DATA, ">$opt{'O'}.cds_files/cds_site_data.txt";
open CELL_DATA, ">$opt{'O'}.cds_files/cds_cell_data.txt";
open DIMS_DATA, ">$opt{'O'}.cds_files/cds_dims_data.txt";
open COUNTS, ">$opt{'O'}.cds_files/cds_counts_matrix.txt";

read_annot($ARGV[1]);
read_dims($ARGV[2]);

#make sure annot matches matrix
open IN, "$ARGV[0]";
$h = <IN>; chomp $h; @H = split(/\t/, $h);
$h_out = "";
for ($i = 0; $i < @H; $i++) {
	if ((defined $CELLID_annot{$H[$i]}) && defined ($CELLID_DIMS{$H[$i]})) {
		$h_out .= "$H[$i]\t";
	}
}
#finish printing out header by removing last \t
$h_out =~ s/\t$//;
print COUNTS "$h_out\n";
#print out header

if (defined $opt{'G'}) {
	print SITE_DATA "site_name\tchr\tbp1\tbp2\tnum_cells_expressed\tsite_length\tgene\n";
} else {
	print SITE_DATA "site_name\tchr\tbp1\tbp2\tnum_cells_expressed\tsite_length\n";
}
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	$site = shift(@P);
	($chr,$bp1,$bp2) = split(/[_]/, $site);
	$siteName = "$chr\_$bp1\_$bp2";
	$siteOut = "$siteName";
	$siteOut2 = "$siteName";
	$chrID = $chr; $chrID =~ s/chr//;
	my $SITENAME_maxSignal=0;
	for ($i = 0; $i < @P; $i++) {
		if ((defined $CELLID_annot{$H[$i]}) && defined ($CELLID_DIMS{$H[$i]})) {
			$SITENAME_totalSignal{$siteName}+=$P[$i];
			
			if ($P[$i] > $SITENAME_maxSignal{$siteName}){
				$SITENAME_maxSignal{$siteName}=$P[$i];
				#print $siteName."\t$P[$i]\t".$SITENAME_maxSignal{$siteName}."\n";
			}
			
			if ($P[$i]>0) {
				$SITENAME_expressed{$siteName}++;
			}
			
			
			
			if ($P[$i]>0) {
				$CELLID_expressed{$H[$i]}++;
			}
			$siteOut .= "\t$P[$i]";
		}
	}
	print COUNTS "$siteOut\n";
	
	
	for ($i = 0; $i < @P; $i++) 
		{
			if ((defined $CELLID_annot{$H[$i]}) && defined ($CELLID_DIMS{$H[$i]})) 
			{
		
			if ($P[$i]>0) {
				$tempval = $P[$i]/$SITENAME_maxSignal{$siteName};
				#print $siteName."\t$P[$i]\t$SITENAME_maxSignal{$siteName}\t".$tempval."\n";
				$siteOut2 .= "\t$tempval";
			}
			else
			{
				$siteOut2 .= "\t0";
			}
			}
		
		}
	
	if (defined $opt{'G'}) {
		if (!defined $SITEID_gene{$siteName}) {
			$SITEID_gene{$siteName} = "NA";
			print STDERR "WARNING: $siteName was not found when identifying gene tags! - setting to NA\n";
		}
		print SITE_DATA "$siteName\t$siteName\t$chrID\t$bp1\t$bp2\t$SITENAME_expressed{$siteName}\t".($bp2-$bp1)."\t$SITEID_gene{$siteName}\n";
	} else {
		print SITE_DATA "$siteName\t$siteName\t$chrID\t$bp1\t$bp2\t$SITENAME_expressed{$siteName}\t".($bp2-$bp1)."\n";
	}
} close IN;
close COUNTS; close SITE_DATA; 

print CELL_DATA "cells\ttimepoint\tnum_genes_expressed\n";
#might need to modify this to be more general than "timepoints"
####Removed DIMS_DATA HEADER
#print DIMS_DATA "dimension_1\tdimension_2\n";

for ($i = 0; $i < @H; $i++) {
	if ((defined $CELLID_annot{$H[$i]}) && defined ($CELLID_DIMS{$H[$i]})) {
		print DIMS_DATA  join("\t",@{$CELLID_DIMS{$H[$i]}})."\n";
		print CELL_DATA "$H[$i]\t$H[$i]\t$CELLID_annot{$H[$i]}\t$CELLID_expressed{$H[$i]}\n";
	}
} close CELL_DATA; close DIMS_DATA;
}
1;
