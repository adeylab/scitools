package sci_commands::plot_dims;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("plot_dims");

sub plot_dims {

@ARGV = @_;
# Defaults
$xdim = 1;
$ydim = 2;
$minV = -10;
$maxV = 10;
$theme = "Clean";
$gradient_def = "BuG90Rd";
$ptSize = 1;
$alpha = 1;
$binary_thresh = 4;
$binary_fail_color = "gray75";
$binary_pass_color = "red3";

getopts("O:A:a:C:c:R:x:y:T:V:M:XS:s:G:p:f:Bb:k:", \%opt);

$die2 = "
scitools plot-dims [options] [dimensions file(s), comma sep]

Options - general:
   -O   [STR]   Output prefix (default is dims file 1 prefix)
   -x   [INT]   X-dimension to plot (def = $xdim)
   -y   [INT]   Y-dimension to plot (def = $ydim)
   -T   [STR]   Theme: (def = $theme)
                  Clean = no axis lines (for tSNE)
                  Reg = regular, include axes and grid (for PCA)
   -p   [FLT]   Point size (def = $ptSize)
   -f   [FLT]   Alpha for plotting points (def = $alpha)

Plotting by annotations:
   -A   [STR]   Annotation file (to color code points)
   -a   [STR]   Comma separated list of annotations to include in plot
                  (requires -A to be specified)
   -C   [STR]   Color coding file (annot (tab) #hexColor)
   -c   [STR]   Color coding string
                  Annot=#hexColor,Annot2=#hexColor  

Plotting by values:				  
   -V   [STR]   Values file. tab-delimited, cellID (tab) value
                  Will plot as teh color of points (overrides
                  annotation colors)
   -S   [MIN,MAX]   Min and max values for scaling (if -V specified)
                  (def = $minV,$maxV)
   -G   [GRD]   Color gradient (def = $gradient_def)
                  For all available gradients, run 'scitools gradients'
   -B           Binary on/off for values plotting (Overrides -S & -G)
   -b   [FLT]   Binary threshold (def = $binary_thresh)
   -k   [NEG,POS]   Colors for failing and passing cells by binary threshold
                  (def = $binary_fail_color,$binary_pass_color)

To plot multiple values specified in a matrix file:
   -M   [STR]   Matrix file (e.g. deviation z-scores from chromVAR)
                  Will create a folder as -O option and produce
                  a plot for each row of the matrix. Overrides -V.
                  Only includes rows with at least one non-zero value.

Other options:
   -R   [STR]   Rscript call (def = $Rscript)
   -s   [STR]   scitools call (def = $scitools)
   -X           Do not delete intermediate files (def = delete)

Note: Requires ggplot2 R package

";

if (!defined $ARGV[0]) {die $die2};
if (defined $opt{'a'} && !defined $opt{'A'}) {die "\nMust provide an annotaiton file (-A) if specifying annotations to plot (-a)!\n$die2"};
if (defined $opt{'C'} && defined $opt{'c'}) {die "\nSpecify either a color string (-c) or a color coding file (-C), not both!\n$die2"};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
$opt{'O'} =~ s/\.txt$//; $opt{'O'} =~ s/\.dims$//;
if (!defined $opt{'G'}) {$opt{'G'} = $gradient_def};
$gradient_function = get_gradient($opt{'G'});
if (defined $opt{'p'}) {$ptSize = $opt{'p'}};
if (defined $opt{'f'}) {$alpha = $opt{'f'}};

read_dims($ARGV[0]);

if (defined $opt{'M'}) {
	print STDERR "SCITOOLS: Matrix file plotting detected! Will plot a separate value-based plot for each row entry of the matrix.\n";
	
	$matrix_out = "matrix_plots";
	
	if (-e "$opt{'O'}.$matrix_out") {
		die "\nFATAL: $opt{'O'}.$matrix_out directory already exists! Exiting!\n$die2";
	}
	system("mkdir $opt{'O'}.$matrix_out");
	
	read_matrix($opt{'M'});
	
	$common_opts = "";
	if (defined $opt{'A'}) {$common_opts .= "-A $opt{'A'} "};
	if (defined $opt{'a'}) {$common_opts .= "-a $opt{'a'} "};
	if (defined $opt{'x'}) {$common_opts .= "-x $opt{'x'} "};
	if (defined $opt{'y'}) {$common_opts .= "-y $opt{'y'} "};
	if (defined $opt{'T'}) {$common_opts .= "-T $opt{'T'} "};
	if (defined $opt{'S'}) {$common_opts .= "-S $opt{'S'} "};
	if (defined $opt{'R'}) {$common_opts .= "-R $opt{'R'} "};
	if (defined $opt{'p'}) {$common_opts .= "-p $opt{'p'} "};
	if (defined $opt{'f'}) {$common_opts .= "-f $opt{'f'} "};
	if (defined $opt{'k'}) {$common_opts .= "-k $opt{'k'} "};
	if (defined $opt{'b'}) {$common_opts .= "-b $opt{'b'} "};
	if (defined $opt{'X'}) {$common_opts .= "-X "};
	if (defined $opt{'B'}) {$common_opts .= "-B "};
	$common_opts =~ s/\s$//;
	
	foreach $feature (keys %MATRIX_feature_nonZero) {
		$feature_polished = $feature;
		$feature_polished =~ s/(:|;|'|\.|\(|\)|\{|\}|\[|\]|\s|\|)/_/;
		open VALS, ">$opt{'O'}.$matrix_out/$feature_polished.values";
		foreach $cellID (keys %CELLID_FEATURE_value) {
			if (defined $CELLID_FEATURE_value{$cellID}{$feature}) {
				print VALS "$cellID\t$CELLID_FEATURE_value{$cellID}{$feature}\n";
			}
		} close VALS;
		system("$scitools plot-dims $common_opts -O $opt{'O'}.$matrix_out/$feature_polished -V $opt{'O'}.$matrix_out/$feature_polished.values -G $opt{'G'} $ARGV[0]");
		if (!defined $opt{'X'}) {system("rm -f $opt{'O'}.$matrix_out/$feature_polished.values")};
	}
	
	exit;
}

if (defined $opt{'R'}) {$Rscript = $opt{'R'}};
if (defined $opt{'s'}) {$scitools = $opt{'s'}};
if (defined $opt{'A'}) {read_annot($opt{'A'})};
if (defined $opt{'a'}) {
	@ANNOT_LIST = split(/,/, $opt{'a'});
	foreach $annot (@ANNOT_LIST) {
		$ANNOT_include{$annot} = 1;
	}
}
if (defined $opt{'C'}) {read_color_file($opt{'C'})};
if (defined $opt{'c'}) {read_color_string($opt{'c'})};
if (defined $opt{'V'}) {read_values($opt{'V'})};
if (defined $opt{'x'}) {$xdim = $opt{'x'}};
if (defined $opt{'y'}) {$ydim = $opt{'y'}};
if (defined $opt{'S'}) {($minV,$maxV) = split(/,/, $opt{'S'})};
if (defined $opt{'b'}) {$binary_thresh = $opt{'b'}};
if (defined $opt{'k'}) {($binary_fail_color,$binary_pass_color) = split(/,/, $opt{'k'})};
if ($binary_fail_color !~ /^#/) {$binary_fail_color = "\"".$binary_fail_color."\""};
if ($binary_pass_color !~ /^#/) {$binary_pass_color = "\"".$binary_pass_color."\""};

open DATA, ">$opt{'O'}.plot.txt";
foreach $cellID (keys %CELLID_DIMS) {

	if (defined $opt{'A'}) {
		if (defined $opt{'a'}) {
			if (defined $CELLID_annot{$cellID}) {
				if (defined $ANNOT_include{$CELLID_annot{$cellID}}) {
					$annot = $CELLID_annot{$cellID}
				} else {$annot = "Exclude"};
			} else {$annot = "Exclude"};
		} else {
			if (defined $CELLID_annot{$cellID}) {
				$annot = $CELLID_annot{$cellID};
			} else {$annot = "Exclude"};
		}
	} else {
		$annot = "Cell";
	}
	
	if ($annot !~ /Exclude/i) {
		if (!defined $opt{'V'}) { # qualitative annotations
			print DATA "$cellID\t$annot\t$CELLID_DIMS{$cellID}[$xdim]\t$CELLID_DIMS{$cellID}[$ydim]\n";
		} else { # value annotations
			if (defined $CELLID_value{$cellID}) {
				if (defined $opt{'B'}) { # binary threshold
					if ($CELLID_value{$cellID}>=$binary_thresh) {
						$annot = "PASS";
					} else {
						$annot = "FAIL";
					}
				}
				print DATA "$cellID\t$CELLID_value{$cellID}\t$CELLID_DIMS{$cellID}[$xdim]\t$CELLID_DIMS{$cellID}[$ydim]\t$annot\n";
			} else {
				print STDERR "WARNING: $cellID does not have values specified in $opt{'V'}, skipping.\n";
			}
		}
	}
	
} close DATA;

open R, ">$opt{'O'}.plot.r";

print R "
library(ggplot2)
IN<-read.table(\"$opt{'O'}.plot.txt\")
$gradient_function";

if (defined $opt{'B'} && defined $opt{'V'}) {
print R "
fail<-subset(IN,V5==\"FAIL\")
pass<-subset(IN,V5==\"PASS\")";
}

print R "
PLT<-ggplot() +";

if (!defined $opt{'c'} && !defined $opt{'C'} && !defined $opt{'A'} && !defined $opt{'V'}) { # no special mode specified
	print R "
	geom_point(aes(IN\$V3,IN\$V4),color=\"lightsteelblue4\",size=$ptSize,alpha=$alpha,shape=16) +";
} elsif (!defined $opt{'V'}) { # annotation specified
	print R "
	geom_point(aes(IN\$V3,IN\$V4,color=IN\$V2),size=$ptSize,alpha=$alpha,shape=16) +
	guides(colour = guide_legend(override.aes = list(size=4,alpha=1))) +";
} else { # values file specified
	if (!defined $opt{'B'}) {
	print R "
	geom_point(aes(IN\$V3,IN\$V4,color=IN\$V2),size=$ptSize,alpha=$alpha,shape=16) +";
	} else {
	print R "
	geom_point(aes(fail\$V3,fail\$V4),color=$binary_fail_color,size=$ptSize,alpha=$alpha,shape=16) +
	geom_point(aes(pass\$V3,pass\$V4),color=$binary_pass_color,size=$ptSize,alpha=$alpha,shape=16) +";
	}
}

if ($color_mapping !~ /none/i && !defined $opt{'V'}) {
	print R "
	scale_colour_manual(values = c($color_mapping)) +";
} elsif (defined $opt{'V'} && !defined $opt{'B'}) {
	print R "
	scale_colour_gradientn(colours=gradient_funct(21),limits=c($minV,$maxV)) +";
}

if ($theme =~ /Clean/i) {
	if (defined $opt{'A'}) {
	print R "
	theme_bw() +
	theme(panel.border=element_blank(),
		  panel.grid=element_blank(),
		  axis.line=element_blank(),
		  axis.ticks=element_blank(),
		  legend.background=element_blank(),
		  legend.title=element_blank(),
		  panel.background=element_blank(),
		  axis.text=element_blank(),
		  axis.title.x=element_blank(),
		  axis.title.y=element_blank(),
		  plot.background=element_blank(),
		  plot.margin=unit(c(0,0,0,0),\"pt\"))\n";
	} else {
	print R "
	theme_bw() +
	theme(panel.border=element_blank(),
		  panel.grid=element_blank(),
		  axis.line=element_blank(),
		  axis.ticks=element_blank(),
		  legend.position=\"none\",
		  panel.background=element_blank(),
		  axis.text=element_blank(),
		  axis.title.x=element_blank(),
		  axis.title.y=element_blank(),
		  plot.background=element_blank(),
		  plot.margin=unit(c(0,0,0,0),\"pt\"))\n";
	}
} else {
	if (defined $opt{'A'}) {
	print R "
	theme_bw()\n";
	} else {
	print R "
	theme_bw(legend.position=\"none\")\n";
	}
}

print R "
ggsave(plot=PLT,filename=\"$opt{'O'}.plot.png\",width=5,height=4,dpi=900)
ggsave(plot=PLT,filename=\"$opt{'O'}.plot.pdf\",width=6,height=5)
";

if (defined $opt{'A'} && defined $opt{'V'}) {
print R "
#Violin plot over opt A 

Violin<-ggplot() + theme_bw() +
	geom_violin(aes(IN\$V5,IN\$V2,fill=IN\$V5),alpha=0.5,color=\"gray30\",size=0.5) +
	geom_jitter(aes(IN\$V5,IN\$V2),color=\"gray30\",size=0.15) +";
if ($color_mapping !~ /none/i) {
	print R "
	scale_fill_manual(values = c($color_mapping)) +";
}
print R "
	guides(fill=FALSE,colour=FALSE) +
	xlab(\"Annot\") +
	ylab(\"Feature value\") 
ggsave(plot=Violin,filename=\"$opt{'O'}.violin.png\",width=7,height=3,dpi=900)
ggsave(plot=Violin,filename=\"$opt{'O'}.violin.pdf\",width=7,height=3)";
}

close R;

system("$Rscript $opt{'O'}.plot.r");

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.plot.r $opt{'O'}.plot.txt");
}

}
1;
