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
$panel_neg_color = "#e0e0e0";
$panel_pos_color = "#08115a";
$binary_fail_color = "gray75";
$binary_pass_color = "red3";
$panel_pass_color = "black";
$width = 5;
$height = 4;
$frames = 100;
$pause = 10;
$gif_height = 600;
$gif_width = 600;

getopts("O:A:a:C:c:R:x:y:T:V:M:XS:s:G:p:f:Bb:k:w:h:m:Wr:N:F:j:k:zu:d", \%opt);

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
   -w   [FLT]   Plot width (inches, def = $width)
   -h   [FLT]   Plot height (inches, def = $height)
   -z           Z-score the dimensions (def = no)

Plotting by annotations:
   -A   [STR]   Annotation file (to color code points)
   -a   [STR]   Comma separated list of annotations to include in plot
                  (requires -A to be specified)
   -C   [STR]   Color coding file (annot (tab) #hexColor)
   -c   [STR]   Color coding string
                  Annot=#hexColor,Annot2=#hexColor
   -W           Wrap panels, one panel per annotation, each with all
                  other annotations grayed out.
				Note: width and height (-w and -h) will be per-panel
   -k   [NEG,POS] Colors for false and true annotation assignments.
                  (def = $panel_neg_color,$panel_pos_color)
   -r   [INT]   Number of rows for panel wrapping (def = $panel_nrow)

Plotting by values:				  
   -V   [STR]   Values file. tab-delimited, cellID (tab) value
                  Will plot as the color of points (overrides
                  annotation colors)
   -S   [MIN,MAX] Min and max values for scaling (if -V specified)
                  (def = $minV,$maxV)
   -G   [GRD]   Color gradient (def = $gradient_def)
                  For all available gradients, run 'scitools gradients'
   -B           Binary on/off for values plotting (Overrides -S & -G)
   -b   [FLT]   Binary threshold (def = $binary_thresh)
   -k   [NEG,POS] Colors for failing and passing cells by binary threshold
                  (def = $binary_fail_color,$binary_pass_color)

To plot multiple values specified in a matrix file:
   -M   [STR]   Matrix file (e.g. deviation z-scores from chromVAR)
                  Will create a folder as -O option and produce
                  a plot for each row of the matrix. Overrides -V.
                  Only includes rows with at least one non-zero value.
				  
To plot Monocle output branch points:
   -m 	[STR] 	Monocle branchpoints.txt file to be plotted underneath 
				the point plots. Generated through cds_monocle command call.
				
Animated plots:
   -N   [MODE]  Toggle animation (requires gganimate & gifski libraries)
                  Height and width values will be interpreted as points (-h,-w)
                  Def height = $gif_height, Def width = $gif_width
                  Modes:
                    Transition/T   Transition to a second dims file that
                                   is provided as an additional argument
                                   after the first dims file. -z recommended.
   -j   [INT]   For transition mode, x-dimension (def = $xdim)
   -k   [INT]   For transition mode, y-dimension (def = $ydim)
   -F   [INT]   Number of frames (def = $frames)
   -u   [INT]   Additional frames to pause for at start and end (def = $pause)
   -d           Rewind plot (def = no)
   
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
if (defined $opt{'h'}) {$height = $opt{'h'}};
if (defined $opt{'w'}) {$width = $opt{'w'}};
if (defined $opt{'M'} && defined $opt{'W'}) {die "\nERROR: Panel wrapping and matrix plotting are not currently compatible.\n"};
if (defined $opt{'W'} && !defined $opt{'A'}) {die "\nERROR: An annotaiton file must be provided for panel wrapping.\n"};
if (defined $opt{'W'} && defined $opt{'N'}) {die "\nERROR: Panel wrapping and animation plotting are not currently compatible.\n"};
if (defined $opt{'M'} && defined $opt{'N'}) {die "\nERROR: Matrix plotting and animation plotting are not currently compatible.\n"};
if ($opt{'N'} =~ /^T/i && !defined $ARGV[1]) {die "\nERROR: Must provide a second dims file for transition animation plotting.\n"};

if (defined $opt{'M'} && defined $opt{'k'}) {
	($panel_neg_color,$panel_pos_color) = split(/,/, $opt{'k'});
}
if (defined $opt{'r'}) {$panel_nrow = $opt{'r'}};

read_dims($ARGV[0]);

if (defined $opt{'M'} && !defined $opt{'N'}) { # matrix plotting
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
	if (defined $opt{'w'}) {$common_opts .= "-w $opt{'w'} "};
	if (defined $opt{'h'}) {$common_opts .= "-h $opt{'h'} "};
	if (defined $opt{'m'}) {$common_opts .= "-m $opt{'m'} "};
	if (defined $opt{'z'}) {$common_opts .= "-z "};
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
if (defined $opt{'W'}) {
	if (defined $opt{'k'}) {($binary_fail_color,$panel_pass_color) = split(/,/, $opt{'k'})};
	if ($binary_fail_color !~ /^#/) {$binary_fail_color = "\"".$binary_fail_color."\""};
	if ($panel_pass_color !~ /^#/) {$panel_pass_color = "\"".$panel_pass_color."\""};
} else {
	if (defined $opt{'k'}) {($binary_fail_color,$binary_pass_color) = split(/,/, $opt{'k'})};
	if ($binary_fail_color !~ /^#/) {$binary_fail_color = "\"".$binary_fail_color."\""};
	if ($binary_pass_color !~ /^#/) {$binary_pass_color = "\"".$binary_pass_color."\""};
}

# read 2nd dims for animation & set up the transition frames
if ($opt{'N'} =~ /^T/i) {

	if (defined $opt{'h'}) {$gif_height = $opt{'h'}};
	if (defined $opt{'w'}) {$gif_width = $opt{'w'}};
	
	if (defined $opt{'u'}) {$pause = $opt{'u'}};

	if (defined $opt{'j'}) {$xdim2 = $opt{'j'}} else {$xdim2 = $xdim};
	if (defined $opt{'k'}) {$ydim2 = $opt{'k'}} else {$ydim2 = $ydim};
	
	if (!defined $ARGV[1]) {die $die};
	open IN, "$ARGV[1]";
	while ($l = <IN>) {
		chomp $l;
		@P = split(/\t/, $l);
		$cellID = $P[0];
		@{$CELLID_DIMS_2{$cellID}} = @P;
	} close IN;
	
	%CELLID_DIM1_FRAME_value = ();
	%CELLID_DIM2_FRAME_value = ();
	
}

$cell_ct = 0;
foreach $cellID (keys %CELLID_DIMS) {
	
	$exclude_cell = 0;
	
	if (defined $opt{'A'}) {
		if (defined $opt{'a'}) {
			if (defined $CELLID_annot{$cellID}) {
				if (defined $ANNOT_include{$CELLID_annot{$cellID}}) {
					$annot = $CELLID_annot{$cellID}
				} else {$exclude_cell = 1};
			} else {$exclude_cell = 1};
		} else {
			if (defined $CELLID_annot{$cellID}) {
				$annot = $CELLID_annot{$cellID};
			} else {$exclude_cell = 1};
		}
	} else {
		$annot = "Cell";
		$CELLID_annot{$cellID} = $annot;
	}
	
	if ($opt{'N'} =~ /^T/i) {
		if (!defined $CELLID_DIMS_2{$cellID}[$xdim2] || !defined $CELLID_DIMS_2{$cellID}[$ydim2]) {
			$exclude_cell = 1;
		}
	}
	
	if ($exclude_cell < 0.5) {
		$CELLID_include{$cellID} = 1;
		$cell_ct++;
	}
	
}

if (defined $opt{'z'}) { # ZSCORE!
	$dim1_sum = 0; $dim2_sum = 0;
	$dim1_mean = 0; $dim2_mean = 0;
	$dim1_stdev_sum = 0; $dim2_stdev_sum = 0;
	$dim1_stdev = 0; $dim2_stdev = 0;
	if ($opt{'N'} =~ /^T/i) {
		$dim1_sum2 = 0; $dim2_sum2 = 0;
		$dim1_mean2 = 0; $dim2_mean2 = 0;
		$dim1_stdev_sum2 = 0; $dim2_stdev_sum2 = 0;
		$dim1_stdev2 = 0; $dim2_stdev2 = 0;
	}
	
	# mean
	foreach $cellID (keys %CELLID_include) {
		$dim1_sum += $CELLID_DIMS{$cellID}[$xdim];
		$dim2_sum += $CELLID_DIMS{$cellID}[$ydim];
		if ($opt{'N'} =~ /^T/i) {
			$dim1_sum2 += $CELLID_DIMS_2{$cellID}[$xdim2];
			$dim2_sum2 += $CELLID_DIMS_2{$cellID}[$ydim2];
		}
	}
	$dim1_mean = $dim1_sum/$cell_ct;
	$dim2_mean = $dim2_sum/$cell_ct;
	if ($opt{'N'} =~ /^T/i) {
		$dim1_mean2 = $dim1_sum2/$cell_ct;
		$dim2_mean2 = $dim2_sum2/$cell_ct;
	}
	
	# stdev
	foreach $cellID (keys %CELLID_include) {
		$dim1_stdev_sum += ($CELLID_DIMS{$cellID}[$xdim]-$dim1_mean)**2;
		$dim2_stdev_sum += ($CELLID_DIMS{$cellID}[$ydim]-$dim2_mean)**2;
		if ($opt{'N'} =~ /^T/i) {
			$dim1_stdev_sum2 += ($CELLID_DIMS_2{$cellID}[$xdim2]-$dim1_mean2)**2;
			$dim2_stdev_sum2 += ($CELLID_DIMS_2{$cellID}[$ydim2]-$dim2_mean2)**2;
		}
	}
	$dim1_stdev = sqrt($dim1_stdev_sum/$cell_ct);
	$dim2_stdev = sqrt($dim2_stdev_sum/$cell_ct);
	if ($opt{'N'} =~ /^T/i) {
		$dim1_stdev2 = sqrt($dim1_stdev_sum2/$cell_ct);
		$dim2_stdev2 = sqrt($dim2_stdev_sum2/$cell_ct);
	}
	
	# replace values in dims object
	foreach $cellID (keys %CELLID_include) {
		$value = ($CELLID_DIMS{$cellID}[$xdim]-$dim1_mean)/$dim1_stdev;
		$CELLID_DIMS{$cellID}[$xdim] = $value;
		$value = ($CELLID_DIMS{$cellID}[$ydim]-$dim2_mean)/$dim2_stdev;
		$CELLID_DIMS{$cellID}[$ydim] = $value;
		if ($opt{'N'} =~ /^T/i) {
			$value = ($CELLID_DIMS_2{$cellID}[$xdim2]-$dim1_mean2)/$dim1_stdev2;
			$CELLID_DIMS_2{$cellID}[$xdim2] = $value;
			$value = ($CELLID_DIMS_2{$cellID}[$ydim2]-$dim2_mean2)/$dim2_stdev2;
			$CELLID_DIMS_2{$cellID}[$ydim2] = $value;
		}
	}
}

open DATA, ">$opt{'O'}.plot.txt";
foreach $cellID (keys %CELLID_include) {
	$annot = $CELLID_annot{$cellID};
	$ANNOT_include{$annot} = 1;
	
	if ($opt{'N'} =~ /^T/) {
	
		$xstart = $CELLID_DIMS{$cellID}[$xdim];
		$xend = $CELLID_DIMS_2{$cellID}[$xdim2];
		$ystart = $CELLID_DIMS{$cellID}[$ydim];
		$yend = $CELLID_DIMS_2{$cellID}[$ydim2];
		
		$slope = ($yend-$ystart)/($xend-$xstart);
		$y_int = $ystart - ($xstart*$slope);
		
		# x progression
		if ($xend > $xstart) {$prog_multiplier = 1} else {$prog_multiplier = -1};
		$xdist = abs($xend-$xstart);
		$xincrement = $xdist/$frames;
		for ($frame = 0; $frame < $frames; $frame++) {
			$xpos = $xstart+($prog_multiplier*$xincrement*$frame);
			$ypos = ($slope*$xpos) + $y_int;
			
			if (!defined $opt{'V'}) { # qualitative annotations
				print DATA "$cellID\t$frame\t$annot\t$xpos\t$ypos\n";
			} else { # value annotations
				######################## Value plotting for animations in progress
				die "ERROR: -V not currently compatible with animation plotting.\n";
				########################
			}
		}
		
		
	} else {
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

if (defined $opt{'W'}) { # panel plotting

print R "
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)

IN<-read.table(\"$opt{'O'}.plot.txt\")\n";

$panel_levels = "c("; $panel_values = "";

foreach $annot (keys %ANNOT_include) {
	$panel_levels .= "\"$annot\",";
	$panel_values .= "\"$panel_pos_color\",";
	push @level_list, $annot; 
}
$panel_levels =~ s/,$//; $panel_levels .= ")";

print R "IN\$V2<-factor(IN\$V2, levels = $panel_levels)

PLT<-ggplot() +
	geom_point(aes(IN\$V3,IN\$V4,color=IN\$V2),size=$ptSize,alpha=$alpha,shape=16) +
	guides(colour = guide_legend(override.aes = list(size=4,alpha=1))) +";

if ($theme =~ /Clean/i) {
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
	theme_bw()\n";
}

$panel_values =~ s/,$//;
print R "
ALL<-PLT + scale_color_manual(values=c($panel_values))\n";
$grid_list = "ALL";
@value_list = split(/,/, $panel_values);
for ($plotID = 0; $plotID < @value_list; $plotID++) {
	$panel_values = "";
	for ($valID = 0; $valID < @value_list; $valID++) {
		if ($valID == $plotID) {
			$panel_values .= "\"$panel_pos_color\",";
		} else {
			$panel_values .= "\"$panel_neg_color\",";
		}
	}
	$panel_values =~ s/,$//;
	# reorder so select annot is plotted on top
#	print R "SELECT<-subset(IN, V2==\"$level_list[$plotID]\")
#OTHER<-subset(IN, V2==\"$level_list[$plotID]\")
#IN<-rbind(OTHER,SELECT)\n";
	# plot
	print R "PLT_$plotID<-PLT + scale_color_manual(values=c($panel_values))\n";
	$grid_list .= ",PLT_$plotID";
}

if (!defined $opt{'r'}) {
	$panel_nrow = int(sqrt(@TYPES));
}

$ncol_factor = (@value_list+1)/$panel_nrow;
$panel_ncol = int($ncol_factor);
if ($ncol_factor>$panel_ncol) {$panel_ncol++};
$grid_width = $width*$panel_ncol;
$grid_height = $height*$panel_nrow;

print R "
PLT_grid<-grid.arrange($grid_list,nrow=$panel_nrow)

ggsave(plot=PLT_grid,filename=\"$opt{'O'}.plot.png\",width=$grid_width,height=$grid_height,dpi=900)
ggsave(plot=PLT_grid,filename=\"$opt{'O'}.plot.pdf\",width=$grid_width,height=$grid_height)
";

}

if (!defined $opt{'W'}) { # Not panel plotting

print R "
library(ggplot2)
IN<-read.table(\"$opt{'O'}.plot.txt\")
$gradient_function";

if (defined $opt{'N'}) {
print R "
library(gganimate)
library(gifski)";
if (defined $opt{'V'}) {
print R "
colnames(IN)<-c(\"cellid\",\"frame\",\"annot\",\"xdim\",\"ydim\",\"value\")";
} else {
print R "
colnames(IN)<-c(\"cellid\",\"frame\",\"annot\",\"xdim\",\"ydim\")";
}
} else {
if (defined $opt{'V'}) {
print R "
colnames(IN)<-c(\"cellid\",\"value\",\"xdim\",\"ydim\",\"annot\")";
} else {
print R "
colnames(IN)<-c(\"cellid\",\"annot\",\"xdim\",\"ydim\")";
}
}


if (defined $opt{'B'} && defined $opt{'V'}) {
print R "
fail<-subset(IN,annot==\"FAIL\")
pass<-subset(IN,annot==\"PASS\")";
}

if (defined $opt{'m'}){
print R "
branch<-read.table(\"$opt{'m'}\",header=T)";
}
print R "
PLT<-ggplot(data=IN) +";
if (defined $opt{'m'}){
print R "
geom_segment(aes(x=branch\$source_prin_graph_dim_1,y=branch\$source_prin_graph_dim_2,xend=branch\$target_prin_graph_dim_1,yend=branch\$target_prin_graph_dim_2)) +";
}

if (!defined $opt{'c'} && !defined $opt{'C'} && !defined $opt{'A'} && !defined $opt{'V'}) { # no special mode specified
	print R "
	geom_point(aes(xdim,ydim),color=\"lightsteelblue4\",size=$ptSize,alpha=$alpha,shape=16) +";
} elsif (!defined $opt{'V'}) { # annotation specified
	print R "
	geom_point(aes(xdim,ydim,color=annot),size=$ptSize,alpha=$alpha,shape=16) +
	guides(colour = guide_legend(override.aes = list(size=4,alpha=1))) +";
} else { # values file specified
	if (!defined $opt{'B'}) {
	print R "
	geom_point(aes(xdim,ydim,color=value),size=$ptSize,alpha=$alpha,shape=16) +";
	} else { # NEED TO VERIFY THIS WITH ADDING THE data=IN IN PLOT INITIALIZATION!
	print R "
	geom_point(aes(fail\$xdim,fail\$ydim),color=$binary_fail_color,size=$ptSize,alpha=$alpha,shape=16) +
	geom_point(aes(pass\$xdim,pass\$ydim),color=$binary_pass_color,size=$ptSize,alpha=$alpha,shape=16) +";
	}
}

if ($color_mapping !~ /none/i && !defined $opt{'V'}) {
	print R "
	scale_colour_manual(values = c($color_mapping)) +";
} elsif (defined $opt{'V'} && !defined $opt{'B'}) {
	print R "
	scale_colour_gradientn(colours=gradient_funct(21),limits=c($minV,$maxV)) +";
}

if (defined $opt{'N'}) {
	print R "
	transition_time(frame) +";
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

if (defined $opt{'N'}) {
	if (defined $opt{'d'}) {
print R "
ANIM<-animate(PLT,start_pause = $pause,end_pause = $pause,rewind = TRUE,height=$gif_height,width=$gif_width)
anim_save(\"$opt{'O'}.animation.gif\",ANIM)
";
	} else {
print R "
ANIM<-animate(PLT,start_pause = $pause,end_pause = $pause,height=$gif_height,width=$gif_width)
anim_save(\"$opt{'O'}.animation.gif\",ANIM)
";
	}
} else {
print R "
ggsave(plot=PLT,filename=\"$opt{'O'}.plot.png\",width=$width,height=$height,dpi=900)
ggsave(plot=PLT,filename=\"$opt{'O'}.plot.pdf\",width=$width,height=$height)
";
}

if (defined $opt{'A'} && defined $opt{'V'} && !defined $opt{'N'}) {
print R "
#Violin plot over opt A 

Violin<-ggplot() + theme_bw() +
	geom_violin(aes(IN\$annot,IN\$value,fill=IN\$annot),alpha=0.5,color=\"gray30\",size=0.5) +
	geom_jitter(aes(IN\$annot,IN\$value),color=\"gray30\",size=0.15) +";
if ($color_mapping !~ /none/i) {
	print R "
	scale_fill_manual(values = c($color_mapping)) +";
}
print R "

	guides(fill=FALSE,colour=FALSE) +
	xlab(\"Annot\") +
	ylab(\"Feature value\") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot=Violin,filename=\"$opt{'O'}.violin.png\",width=7,height=3,dpi=900)
ggsave(plot=Violin,filename=\"$opt{'O'}.violin.pdf\",width=7,height=3)

#add anova here 

folder<-strsplit(\"$opt{'O'}\", \"\/\")[[1]]
feature<-strsplit(\"$opt{'O'}\", \"\/\")[[1]]

fit <- aov(value ~ annot, data=IN)
F<-summary(fit)[[1]][[\"F value\"]][1]
P<-summary(fit)[[1]][[\"Pr(>F)\"]][1]

output<-c(feature[2],F,P)
write(output,file=paste0(folder[1],\"/Annova_sum.txt\"),append=TRUE,ncolumns=3,sep = \"\\t\")

#individual comparisons with TukeyHSD
comp<-TukeyHSD(fit)
write.table(as.matrix(comp\$annot),file=\"./$opt{'O'}.TukeyHSD.txt\",quote=FALSE,sep=\"\\t\",row.names=T,col.names=T)";

}

close R;

}

system("$Rscript $opt{'O'}.plot.r");

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.plot.r $opt{'O'}.plot.txt");
}

}
1;
