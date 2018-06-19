package sci_commands::plot_reads;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("plot_reads");

sub plot_reads {

@ARGV = @_;
# Defaults
$height = 6;
$width = 6;
$pt_size = 0.5;
$flanking_size = 100000;
$gene_scale_factor = 1;
$gene_text_size = 1.5;

getopts("O:A:a:C:c:R:Xs:Dh:w:rp:B:G:S:f:t:V:v:", \%opt);

$die2 = "
scitools plot-reads [options] [rmdup sci bam file] [chrN:start-end] [region 2] ...

Options:
   -O   [STR]   Output prefix (default is input prefix)
   -D           Output to a directory (will add to it if it exists)
   -B   [BED]   Bed file of peaks (optional)
   -G   {STR]   Gene info (refGene.txt formats)
                Shortcut eg: hg38, hg19, mm10 (see scitools -h for more details)
   -S   [INT]   Flanking size if gene names are specified
                (reuires -G, def = $flanking_size)
   -A   [STR]   Annotation file (to color code reads)
   -a   [STR]   Comma separated list of annotations to include in plot
                  (requires -A, will plot annots in the specified order)
   -C   [STR]   Color coding file (annot (tab) #hexColor)
   -c   [STR]   Color coding string
                  Annot=#hexColor,Annot2=#hexColor
   -h   [IN]    Height (inches, def = $height)
   -w   [IN]    Width (inches, def = $width)
   -p   [FLT]   Point size (def = $pt_size)
   -f   [FLT]   Gene plot spacing factor (def = $gene_scale_factor)
                  (for -G, larger values = more vertical spread)
   -t   [FLT]   Gene name text size (def = $gene_text_size)
   -r           Order cells by start of first read (def = randomize)
   -V   [STR]   Values or lambda file for cell ordering
                  (overrides -r, only cells in file will plot)
   -v           Reverse order of the provided values/lambda file
   -R   [STR]   Rscript call (def = $Rscript)
   -s   [STR]   Samtools call (def = $samtools)
   -X           Do not delete intermediate files (def = delete)

Note: Requires ggplot2 R package

";

if (!defined $ARGV[1]) {die $die2};
if (defined $opt{'a'} && !defined $opt{'A'}) {die "\nMust provide an annotaiton file (-A) if specifying annotations to plot (-a)!\n$die2"};
if (defined $opt{'C'} && defined $opt{'c'}) {die "\nSpecify either a color string (-c) or a color coding file (-C), not both!\n$die2"};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};
if (defined $opt{'s'}) {$samtools = $opt{'s'}};
if (defined $opt{'h'}) {$height = $opt{'h'}};
if (defined $opt{'w'}) {$width = $opt{'w'}};
if (defined $opt{'S'}) {$flanking_size = $opt{'S'}};
if (defined $opt{'f'}) {$gene_scale_factor = $opt{'f'}};
if (defined $opt{'t'}) {$gene_text_size = $opt{'t'}};
if (defined $opt{'A'}) {
	read_annot($opt{'A'});
	foreach $annot (keys %ANNOT_count) {
		$ANNOT_ids{$annot} = 0;
	}
} else {
	$ANNOT_ids{"Cell"} = 0;
}
if (defined $opt{'a'}) {
	@ANNOT_LIST = split(/,/, $opt{'a'});
	foreach $annot (@ANNOT_LIST) {
		$ANNOT_include{$annot} = 1;
	}
}
if (defined $opt{'C'}) {read_color_file($opt{'C'})};
if (defined $opt{'c'}) {read_color_string($opt{'c'})};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
$opt{'O'} =~ s/\.bam$//;

if (defined $opt{'G'}) {
	if (defined $REF{$opt{'G'}}) {
		$ref_file = $REF{$opt{'G'}};
		$opt{'G'} = $ref_file;
		$opt{'G'} =~ s/\.fa$/\.refGene.txt/;
	}
	read_refgene($opt{'G'});
}

if (defined $opt{'D'}) {
	if (-e "$opt{'O'}.read_plots") {
		print STDERR "INFO: Output directory exists - adding new regions to $opt{'O'}.read_plots.\n";
	} else {
		system("mkdir $opt{'O'}.read_plots");
	}
	$opt{'O'} .= ".read_plots/";
} else {$opt{'O'} .= "read_plot\."};

$bam = shift(@ARGV);
if (-e "$bam.bai") {
#	print STDERR "\n\nINFO: bai file found.\n";
} else {
	print STDERR "\n\nINFO: No bai file found - generating.\n";
	system("$samtools index $bam");
}

# LOOP THROUGH SPECIFIED REGIONS
foreach $region (@ARGV) {
print STDERR "INFO: Processing $region\n";
if ($region =~ /.+:.+-.+/) { # region
	($chr,$start,$end) = split(/[:-]/, $region);
	$region_name = "$chr\_$start\_$end";
} else { # gene?
	if (!defined $opt{'G'}) {
		print STDERR "WARNING: Cannot interpret coordinates $region - is it a gene name? If so, specify -G\n";
		$region = "skip";
	} else {
		if (defined $GENENAME_coords{$region}) {
			$region_name = $region; $region_name =~ s/[-:,\.\+\|\(\)]/_/g;
			$region = $GENENAME_coords{$region};
			($chr,$start,$end) = split(/[:-]/, $region);
			$region = "$chr:".($start-$flanking_size)."-".($end+$flanking_size);
		} elsif (defined $GENEID_coords{$region}) {
			$region_name = $region; $region_name =~ s/[-:,\.\+\|\(\)]/_/g;
			$region = $GENEID_coords{$region};
			($chr,$start,$end) = split(/[:-]/, $region);
			$region = "$chr:".($start-$flanking_size)."-".($end+$flanking_size);
		} else {
			print STDERR "WARNING: Cannot interpret coordinates $region - if it is a gene name or accession, it is not in $opt{'G'}\n";
			$region = "skip";
		}
	}
}

if ($region ne "skip") {
($chr,$start,$end) = split(/[:-]/, $region);
open OUT, ">$opt{'O'}$region_name.data";
open IN, "$samtools view $bam $region |";
$reads_in_region = 0; $cell_ct = 0; %CELLS_present = ();
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	$cellID = $P[0]; $cellID =~ s/:.+$//; 
	if (!defined $CELLS_present{$cellID}) {
		$CELLS_present{$cellID} = 1;
		$cell_ct++;
	}
	$annot = "null";
	if (!defined $opt{'A'}) {
		$annot = "Cell"; $CELLID_annot{$cellID} = "Cell";
	} elsif (defined $CELLID_annot{$cellID}) {
		if ((defined $opt{'a'} && defined $ANNOT_include{$CELLID_annot{$cellID}}) || !defined $opt{'a'}) {
			$annot = $CELLID_annot{$cellID};
		}
	}
	if ($annot ne "null") {
		$posn = int($P[3]+(length($P[9])/2));
		if (defined $CELLID_id{$cellID}) {
			$id = $CELLID_id{$cellID};
		} else {
			$id = $ANNOT_ids{$annot}; $ANNOT_ids{$annot}++;
			$CELLID_id{$cellID} = $id;
			if (!defined $opt{'r'}) {$ID_randVal{$cellID} = rand(1e20)};
		}
		print OUT "$annot\t$id\t$posn\n";
		$reads_in_region++;
	}
} close IN; close OUT;

if ($reads_in_region>0) {

# order annot if -a specified, otherwise order however
@ANNOT_ORDER = ();
if (!defined $opt{'a'}) {
	foreach $annot (sort keys %ANNOT_ids) {
		push @ANNOT_ORDER, $annot;
	}
} else {
	for ($i = 0; $i < @ANNOT_LIST; $i++) {
		if (defined $ANNOT_ids{$ANNOT_LIST[$i]}) {
			$annot = $ANNOT_LIST[$i];
			push @ANNOT_ORDER, $annot;
		}
	}
}

# DO reordering stuff here
$newID = 0;
if (!defined $opt{'V'}) {
	for ($i = 0; $i < @ANNOT_ORDER; $i++) {
		$annot = $ANNOT_ORDER[$i];
		if (!defined $opt{'r'}) {
			foreach $cellID (sort {$ID_randVal{$a}<=>$ID_randVal{$b}} keys %ID_randVal) {
				if ($CELLID_annot{$cellID} eq $annot) {
					$ANNOT_CELLID_newID{$annot}{$CELLID_id{$cellID}} = $newID;
					$newID++;
				}
			}
		} else {
			foreach $cellID (sort {$CELLID_id{$a}<=>$CELLID_id{$b}} keys %CELLID_id) {
				if ($CELLID_annot{$cellID} eq $annot) {
					$ANNOT_CELLID_newID{$annot}{$CELLID_id{$cellID}} = $newID;
					$newID++;
				}
			}
		}
	}
} else {
	read_values($opt{'V'});
	if (!defined $opt{'v'}) {
		foreach $cellID (sort {$CELLID_value{$a}<=>$CELLID_value{$b}} keys %CELLID_value) {
			if (defined $CELLID_id{$cellID} && defined $CELLID_annot{$cellID}) {
				$ANNOT_CELLID_newID{$CELLID_annot{$cellID}}{$CELLID_id{$cellID}} = $newID;
				$newID++;
			}
		}
	} else {
		foreach $cellID (sort {$CELLID_value{$b}<=>$CELLID_value{$a}} keys %CELLID_value) {
			if (defined $CELLID_id{$cellID} && defined $CELLID_annot{$cellID}) {
				$ANNOT_CELLID_newID{$CELLID_annot{$cellID}}{$CELLID_id{$cellID}} = $newID;
				$newID++;
			}
		}
	}
}

open IN, "$opt{'O'}$region_name.data";
open OUT, ">$opt{'O'}$region_name.data.full";
while ($l = <IN>) {
	chomp $l;
	($annot,$id,$posn) = split(/\t/, $l);
	if (defined $ANNOT_CELLID_newID{$annot}{$id}) {
		$newID = $ANNOT_CELLID_newID{$annot}{$id};
		print OUT "READ\t$annot\t$newID\t$posn\n";
	}
} close IN;

if (defined $opt{'B'}) {
	$peaks_in_region = 0;
	open IN, "$opt{'B'}";
	while ($l = <IN>) {
		chomp $l;
		@P = split(/\t/, $l);
		if ($P[0] eq $chr && $P[1] > $start && $P[2] < $end) {
			$posn = int(($P[1]+$P[2])/2);
			print OUT "PEAK\tpeak\t0\t$posn\n";
			$peaks_in_region++;
		}
	}
	if ($peaks_in_region==0) {
		print STDERR "WARNING: There are no peaks in $region_name, proceeding without peak plotting.\n";
	}
}

close OUT;
if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}$region_name.data && mv $opt{'O'}$region_name.data.full $opt{'O'}$region_name.data");
} else {
	system("mv $opt{'O'}$region_name.data $opt{'O'}$region_name.data.initial && mv $opt{'O'}$region_name.data.full $opt{'O'}$region_name.data");
}

if (defined $opt{'G'}) {
	$scale_increment = int(($cell_ct/80)*$gene_scale_factor);
	$gene_num = -1*$scale_increment;
	$genes_in_region = 0;
	open IN, "$opt{'G'}";
	open OUT, ">$opt{'O'}$region_name.gene_data";
	while ($l = <IN>) {
		chomp $l;
		if ($l !~ /^#/) {
			@P = split(/\t/, $l);
			if ($P[2] eq $chr && !defined $GENES{$P[12]}) {
				if ($P[4]>=$start && $P[5]<=$end) { # internal to window
					$gene_start = $P[4]; $gene_end = $P[5];
				} elsif ($P[4]<$start && $P[5]>=$end) { # starts before and ends after
					$gene_start = $start; $gene_end = $end;
				} elsif ($P[4]<$start && $P[5]>$start) { # starts before, ends internal
					$gene_start = $start; $gene_end = $P[5];
				} elsif ($P[4]>$start && $P[4]<$end) { # starts internal, ends after
					$gene_start = $P[4]; $gene_end = $end;
				} else { #external
					$gene_start = "null";
				}
				if ($gene_start ne "null") { # use it
					$genes_in_region++;
					$gene_num-=$scale_increment;
					$GENES{$P[12]} = 1;
					if ($P[3] eq "+") {
						print OUT "T\t$gene_start\t$gene_end\t".($gene_start-5000)."\t$P[12]\t$gene_num\n";
					} else {
						print OUT "T\t$gene_end\t$gene_start\t".($gene_start-5000)."\t$P[12]\t$gene_num\n";
					}
					@E_STARTS = split(/,/, $P[9]);
					@E_ENDS = split(/,/, $P[10]);
					for ($exon = 0; $exon < @E_STARTS; $exon++) {
						$e_start = $E_STARTS[$exon]; $e_end = $E_ENDS[$exon];
						if ($e_start>=$start && $e_end<=$end) { # internal to window
							# keep values
						} elsif ($e_start<$start && $e_end>=$end) { # starts before and ends after
							$e_start = $start; $e_end = $end;
						} elsif ($e_start<$start && $e_end>$start) { # starts before, ends internal
							$e_start = $start;
						} elsif ($e_start>$start && $e_start<$end) { # starts internal, ends after
							$e_end = $end;
						} else { #external
							$e_start = "null";
						}
						if ($e_start ne "null") {
							print OUT "E\t$e_start\t$e_end\t$exon\t$P[12]\t$gene_num\n";
						}
					}
				}
			}
		}
	} close IN; close OUT;
}


open R, ">$opt{'O'}$region_name.r";
print R "# plotting region $region
library(ggplot2)
IN <- subset(read.table(\"$opt{'O'}$region_name.data\"),V1==\"READ\")
colnames(IN) <- c(\"type\",\"annot\",\"y\",\"x\")";

if (defined $opt{'B'} && $peaks_in_region>0) {
	print R "
PEAKS <- subset(read.table(\"$opt{'O'}$region_name.data\"),V1==\"PEAK\")
colnames(PEAKS) <- c(\"type\",\"annot\",\"y\",\"x\")";
}

if (defined $opt{'G'} && $genes_in_region>0) {
	print R "
TRANS <- subset(read.table(\"$opt{'O'}$region_name.gene_data\"),V1==\"T\")
colnames(TRANS) <- c(\"type\",\"xstart\",\"xend\",\"namepos\",\"name\",\"y\")
EXONS <- subset(read.table(\"$opt{'O'}$region_name.gene_data\"),V1==\"E\")
colnames(EXONS) <- c(\"type\",\"start\",\"end\",\"exon\",\"name\",\"y\")";
}

print R "
PLT <- ggplot() + theme_bw() +";

if (defined $opt{'B'} && $peaks_in_region>0) {
	print R "
	geom_vline(aes(xintercept=PEAKS\$x),size=0.75,color=\"gray95\") +";
}

if (defined $opt{'G'} && $genes_in_region>0) {
	print R "
	geom_segment(aes(x=TRANS\$xstart,xend=TRANS\$xend,y=TRANS\$y,yend=TRANS\$y),size=0.3,color=\"darkblue\",arrow=arrow(length=unit(3,\"points\"))) +
	geom_segment(aes(x=EXONS\$start,xend=EXONS\$end,y=EXONS\$y,yend=EXONS\$y),size=1,color=\"darkblue\") +
	geom_text(aes(x=TRANS\$namepos,y=TRANS\$y,label=TRANS\$name),size=$gene_text_size,color=\"black\",hjust=1) +";
}

if (!defined $opt{'c'} && !defined $opt{'C'} && !defined $opt{'A'}) {
	print R "
	geom_point(aes(IN\$x,IN\$y),color=\"lightsteelblue4\",size=$pt_size,shape=15) +";
} else {
	print R "
	geom_point(aes(IN\$x,IN\$y,color=IN\$annot),size=$pt_size,shape=15) +"
}

if ($color_mapping !~ /none/i) {
	print R "
	scale_colour_manual(values = c($color_mapping)) +";
}

print R "
	scale_y_reverse(expand=c(0.02,0.02)) +
	scale_x_continuous(expand=c(0,0)) +
	xlab(\"$region_name\") +
	ylab(\"Cells\") +
	guides(colour = guide_legend(override.aes = list(size=4))) +
	theme(strip.background=element_blank(),";
if (defined $opt{'B'} && $peaks_in_region>0) {
print R "
		panel.grid.major.x=element_blank(),
		panel.grid.minor.x=element_blank(),";
}
if (defined $opt{'A'}) {
	print R "
		axis.ticks.y=element_blank(),
		axis.text.y=element_blank(),
		legend.title=element_blank(),
		panel.background=element_blank(),
		plot.background=element_blank())";
} else {
	print R "
		axis.ticks.y=element_blank(),
		axis.text.y=element_blank(),
		legend.position=\"none\",
		panel.background=element_blank(),
		plot.background=element_blank())";
}
print R "
ggsave(PLT,filename=\"$opt{'O'}$region_name.png\",width=$width,height=$height,dpi=900)
ggsave(PLT,filename=\"$opt{'O'}$region_name.pdf\",width=$width,height=$height)
";

system("$Rscript $opt{'O'}$region_name.r");

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}$region_name.r $opt{'O'}$region_name.*data*");
}

} else {
	print STDERR "WARNING: $region_name has no reads in the provided bam file. Skipping!\n";
}
}
}

}
1;
