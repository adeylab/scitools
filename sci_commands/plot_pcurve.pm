package sci_commands::plot_pcurve;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("plot_pcurve");

sub plot_pcurve {

@ARGV = @_;
# Defaults
$xdim = 1;
$ydim = 2;
$alpha = 1;
$ptSize = 1;
$theme = "Clean";
$gradient_def = "YlOrRd";

getopts("O:A:a:C:c:R:x:y:T:Xo:p:l:G:M:f:s:", \%opt);

$die2 = "
scitools plot-dims [options] [pcurve prefix]

Will search for [prefix].orig.dims
                [prefix].proj.dims
                [prefix].lambda (or [prefix].centered.lambda)
                                (or [prefix].pruned.lambda)

Options:
   -O   [STR]   Output prefix (default is input prefix)
   -o   [STR]   Original dims file (def = auto find [prefix].orig.dims)
   -p   [STR]   Pcurve dims file (def = auto find [prefix].proj.dims)
   -l   [STR]   Lambda file (def = auto find [prefix].lambda)
                (if [prefix].centered.lambda is found, it will it)
   -M   [STR]   Matrix file - will plot each value set in a matrix
                along the lambda values of the pcurve
   -f   [FLT]   Alpha of points (for -M plotting; def = $alpha)
   -s   [FLT]   Point size (def = $ptSize)
   -A   [STR]   Annotation file (to color code points)
   -a   [STR]   Comma separated list of annoations to include in plot
                  (requires -A to be specified)
   -x   [INT]   X-dimension to plot (def = $xdim)
   -y   [INT]   Y-dimension to plot (def = $ydim)
   -T   [STR]   Theme: (def = $theme)
                  Clean = no axis lines (for tSNE)
                  Reg = regular, include axes and grid (for PCA)
   -C   [STR]   Color coding file (annot (tab) #hexColor)
   -c   [STR]   Color coding string
                  Annot=#hexColor,Annot2=#hexColor
   -G   [GRD]   Color gradient for labmda (def = $gradient_def)
                  For all available gradients, run 'scitools gradient'
   -R   [STR]   Rscript call (def = $Rscript)
   -X           Do not delete intermediate files (def = delete)

Note: Requires ggplot2 R package

If -O -o -p and -l are all defined [pcurve prefix] does not need to be
specified to execute this command.

";

if (!defined $ARGV[0]) { 
    if(!defined $opt{'O'} &&
	   !defined $opt{'o'} &&
	   !defined $opt{'p'} &&
	   !defined $opt{'l'}) {die $die2};
}
if (defined $opt{'a'} && !defined $opt{'A'}) {die "\nMust provide an annotaiton file (-A) if specifying annotations to plot (-a)!\n$die2"};
if (defined $opt{'C'} && defined $opt{'c'}) {die "\nSpecify either a color string (-c) or a color coding file (-C), not both!\n$die2"};

if (defined $opt{'R'}) {$Rscript = $opt{'R'}};
if (defined $opt{'s'}) {$ptSize = $opt{'s'}};
if (defined $opt{'f'}) {$alpha = $opt{'f'}};
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
$opt{'O'} =~ s/\.txt$//; $opt{'O'} =~ s/\.dims$//;
if (defined $opt{'x'}) {$xdim = $opt{'x'}};
if (defined $opt{'y'}) {$ydim = $opt{'y'}};
if (!defined $opt{'G'}) {$opt{'G'} = $gradient_def};
$gradient_function = get_gradient($opt{'G'});

if (!defined $opt{'o'}) {
	if (-e "$ARGV[0].orig.dims") {
		$opt{'o'} = "$ARGV[0].orig.dims";
	} else {
		die "\n\nERROR: Cannot find $ARGV[0].orig.dims, check your prefix ($ARGV[0])\n";
	}
}
if (!defined $opt{'p'}) {
	if (-e "$ARGV[0].proj.dims") {
		$opt{'p'} = "$ARGV[0].proj.dims";
	} else {
		die "\n\nERROR: Cannot find $ARGV[0].proj.dims, check your prefix ($ARGV[0])\n";
	}
}
if (!defined $opt{'l'}) {
	if (-e "$ARGV[0].centered.lambda") {
		$opt{'l'} = "$ARGV[0].centered.lambda";
	} elsif (-e "$ARGV[0].lambda") {
		$opt{'l'} = "$ARGV[0].lambda";
	} else {
		die "\n\nERROR: Cannot find $ARGV[0].centered.lambda OR $ARGV[0].lambda, check your prefix ($ARGV[0])\n";
	}
}

read_dims($opt{'o'});
read_pcurve_dims($opt{'p'});
read_values($opt{'l'});


if (defined $opt{'M'}) {
	print STDERR "SCITOOLS: Matrix file plotting detected! Will plot a separate value-based plot for each row entry of the matrix.\n";
	
	$matrix_out = "matrix_plots";
	
	if (-e "$opt{'O'}.$matrix_out\pcurve") {
		die "\nFATAL: $opt{'O'}.$matrix_out\pcurve directory already exists! Exiting!\n$die2";
	}
	system("mkdir $opt{'O'}.$matrix_out\pcurve");
	
	read_matrix($opt{'M'});
	
	foreach $feature (keys %MATRIX_feature_nonZero) {
		$feature_polished = $feature;
		$feature_polished =~ s/(:|;|'|\.|\(|\)|\{|\}|\[|\]|\s|\|)/_/;
		open VALS, ">$opt{'O'}.$matrix_out\pcurve/$feature_polished.values";

		foreach $cellID (keys %CELLID_FEATURE_value) {
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
					} else {$annot = "Cell"};
				}
			} else {
				$annot = "Cell";
			}
			if (defined $CELLID_FEATURE_value{$cellID}{$feature} && $CELLID_value{$cellID} && ($annot !~ /Exclude/i)) {
				print VALS "$cellID\t$CELLID_value{$cellID}\t$CELLID_FEATURE_value{$cellID}{$feature}\t$annot\n";
			}
		} close VALS;
		
		open R, ">$opt{'O'}.plot.r";

# lambda vals plot
	print R "
library(ggplot2)
IN<-read.table(\"$opt{'O'}.$matrix_out\pcurve/$feature_polished.values\",header=FALSE)
# Make the lambda and matrix val plot
PointPlot<-ggplot() + ";

# plot the original dims w/ respective color
if (!defined $opt{'c'} && !defined $opt{'C'} && !defined $opt{'A'}) {
	print R "
	geom_point(aes(x=IN\$V2,y=IN\$V3),color=\"lightsteelblue4\",size=$ptSize,alpha=$alpha,shape=16) +";
} else {
	print R "
	geom_point(aes(x=IN\$V2,y=IN\$V3,color=IN\$V4),size=$ptSize,alpha=$alpha,shape=16) +
	guides(colour = guide_legend(override.aes = list(size=4))) +";
}

if ($color_mapping !~ /none/i) {
	print R "
	scale_colour_manual(values = c($color_mapping)) +";
}

print R "
geom_smooth(aes(x=IN\$V2,y=IN\$V3)) +
theme_bw() + xlab(\"Lambda\") + ylab(\"Feature value\") +";
if (defined $opt{'A'}) {
print R "
	theme(legend.background=element_blank(),
		legend.title=element_blank(),
		panel.background=element_blank(),
		plot.background=element_blank())";
} else {
print R "
	theme(legend.position=\"none\",
		panel.background=element_blank(),
		plot.background=element_blank())";
}
print R "

ggsave(plot=PointPlot,filename=\"$opt{'O'}.$matrix_out\pcurve/$feature_polished.lambda.png\",width=5,height=4,dpi=900)
ggsave(plot=PointPlot,filename=\"$opt{'O'}.$matrix_out\pcurve/$feature_polished.lambda.plot.pdf\",width=5,height=4)

a=cor(IN\$V2,IN\$V3,method =\"pearson\")
b=cor(IN\$V2,IN\$V3,method =\"spearman\")
output<-c(\"$feature_polished\",a)
output2<-c(\"$feature_polished\",b)
write(output,file=\"$opt{'O'}.$matrix_out\pcurve/corr.pearson.txt\",append=TRUE,ncolumns=2,sep = \"\\t\")
write(output2,file=\"$opt{'O'}.$matrix_out\pcurve/corr.spearman.txt\",append=TRUE,ncolumns=2,sep = \"\\t\")

";
	close R;
	
		system("$Rscript $opt{'O'}.plot.r");
		if (!defined $opt{'X'}) {system("rm -f $opt{'O'}.$matrix_out\pcurve/$feature_polished.values $opt{'O'}.plot.r")};
	}
	
	
	open R, ">$opt{'O'}.plot.r";
	print R "
	library(ggplot2)
	IN<-read.table(\"$opt{'O'}.$matrix_out\pcurve/corr.pearson.txt\",header=FALSE)
	# Make the corr distribution plot
	
	HISTPlot<-ggplot(IN,aes(x=V2)) + geom_histogram()+ 
	theme_bw() + xlab(\"Correlation\") + ylab(\"Counts\") +";
if (defined $opt{'A'}) {
print R "
	theme(legend.background=element_blank(),
		legend.title=element_blank(),
		panel.background=element_blank(),
		plot.background=element_blank())";
} else {
print R "
	theme(legend.position=\"none\",
		panel.background=element_blank(),
		plot.background=element_blank())";
}
print R "

	ggsave(plot=HISTPlot,filename=\"$opt{'O'}.$matrix_out\pcurve/Hist.corr.Pearson.lambda.png\",width=5,height=4,dpi=900)
	ggsave(plot=HISTPlot,filename=\"$opt{'O'}.$matrix_out\pcurve/Hist.corr.Pearson.lambda.plot.pdf\",width=5,height=4)
	
		IN<-read.table(\"$opt{'O'}.$matrix_out\pcurve/corr.spearman.txt\",header=FALSE)
	# Make the corr distribution plot
	
	HISTPlot<-ggplot(IN,aes(x=V2)) + geom_histogram()+ 
	theme_bw() + xlab(\"Correlation\") + ylab(\"Counts\") +";
if (defined $opt{'A'}) {
print R "
	theme(legend.background=element_blank(),
		legend.title=element_blank(),
		panel.background=element_blank(),
		plot.background=element_blank())";
} else {
print R "
	theme(legend.position=\"none\",
		panel.background=element_blank(),
		plot.background=element_blank())";
}
print R "

	ggsave(plot=HISTPlot,filename=\"$opt{'O'}.$matrix_out\pcurve/Hist.corr.Spearman.lambda.png\",width=5,height=4,dpi=900)
	ggsave(plot=HISTPlot,filename=\"$opt{'O'}.$matrix_out\pcurve/Hist.corr.Spearman.lambda.plot.pdf\",width=5,height=4)
	";
	close(R);
	system("$Rscript $opt{'O'}.plot.r");
	if (!defined $opt{'X'}) {system("rm -f $opt{'O'}.plot.r")};
	
	exit;
}


open DATA, ">$opt{'O'}.plot.txt";
print DATA "cellID\tannot\tlambda\todim1\todim2\tpdim1\tpdim2\n";
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
			} else {
				$annot = "Cell";
			}
		}
	} else {
		$annot = "Cell";
	}
	
	if ($annot !~ /Exclude/i && defined $CELLID_PCURVE_DIMS{$cellID}[$xdim]) {
		print DATA "$cellID\t$annot\t$CELLID_value{$cellID}\t$CELLID_DIMS{$cellID}[$xdim]\t$CELLID_DIMS{$cellID}[$ydim]\t$CELLID_PCURVE_DIMS{$cellID}[$xdim]\t$CELLID_PCURVE_DIMS{$cellID}[$ydim]\n";
	}
	
} close DATA;

open R, ">$opt{'O'}.plot.r";

# DIMS PLOT
print R "
library(ggplot2)
IN<-read.table(\"$opt{'O'}.plot.txt\",header=TRUE)

# Make the odim and pdim plot
PointPlot<-ggplot() +";

# plot the original dims w/ respective color
if (!defined $opt{'c'} && !defined $opt{'C'} && !defined $opt{'A'}) {
	print R "
	geom_point(aes(IN\$odim1,IN\$odim2),color=\"lightsteelblue4\",size=$ptSize,alpha=$alpha,shape=16) +";
} else {
	print R "
	geom_point(aes(IN\$odim1,IN\$odim2,color=IN\$annot),size=$ptSize,alpha=$alpha,shape=16) +
	guides(colour = guide_legend(override.aes = list(size=4))) +";
}

if ($color_mapping !~ /none/i) {
	print R "
	scale_colour_manual(values = c($color_mapping)) +";
}

# plot the projected points and color by lambda
print R "
	geom_point(aes(IN\$pdim1,IN\$pdim2),color=\"gray30\",size=$ptSize,alpha=$alpha,shape=16) +";

# set theme
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

# plot them
print R "
ggsave(plot=PointPlot,filename=\"$opt{'O'}.plot.png\",width=5,height=4,dpi=900)
ggsave(plot=PointPlot,filename=\"$opt{'O'}.plot.pdf\",width=5,height=4)";

# DIMS PLOT - "CENTIPEDE"
print R "
# Make the odim and pdim \"centipede\" style plot
$gradient_function
CentPlot<-ggplot() +";

# plot the original dims w/ respective color
print R "geom_segment(aes(x=IN\$odim1,xend=IN\$pdim1,y=IN\$odim2,yend=IN\$pdim2),size=0.2,color=\"lightsteelblue4\") +
	geom_point(aes(IN\$odim1,IN\$odim2),color=\"lightsteelblue4\",size=0.5,shape=16) +
	geom_point(aes(IN\$pdim1,IN\$pdim2,color=IN\$lambda),size=0.75,shape=16) +
	scale_color_gradientn(colours=gradient_funct(15)) +
	labs(color=\"Lambda\") +";

# set theme
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

# plot them
print R "
ggsave(plot=CentPlot,filename=\"$opt{'O'}.vector_plot.png\",width=5,height=4,dpi=900)
ggsave(plot=CentPlot,filename=\"$opt{'O'}.vector_plot.pdf\",width=5,height=4)";

# DENSITY LAMBDA PLOT
print R "
# Plot the density over lambda

DensityPlot<-ggplot() +";
if (!defined $opt{'c'} && !defined $opt{'C'} && !defined $opt{'A'}) {
	print R "
	geom_density(aes(IN\$lambda),color=\"lightsteelblue4\",size=1) +";
} else {
	print R "
	geom_density(aes(IN\$lambda,color=IN\$annot),size=1) +
	guides(colour = guide_legend(override.aes = list(size=4))) +";
}
if ($color_mapping !~ /none/i) {
	print R "
	scale_colour_manual(values = c($color_mapping)) +";
}
print R "
	theme_bw() + xlab(\"Lambda\") +	ylab(\"Density\") +";
if (defined $opt{'A'}) {
print R "
	theme(legend.background=element_blank(),
		legend.title=element_blank(),
		panel.background=element_blank(),
		plot.background=element_blank())";
} else {
print R "
	theme(legend.position=\"none\",
		panel.background=element_blank(),
		plot.background=element_blank())";
}
print R "

ggsave(plot=DensityPlot,filename=\"$opt{'O'}.lambda.png\",width=5,height=3,dpi=900)
ggsave(plot=DensityPlot,filename=\"$opt{'O'}.lambda.pdf\",width=5,height=3)";

# VIOLIN LAMBDA PLOT
if (defined $opt{'A'}) {
print R "
# Horizontal violin plot over lambda

Violin<-ggplot() + theme_bw() +
	geom_violin(aes(IN\$annot,IN\$lambda,fill=IN\$annot),alpha=0.5,color=\"gray30\",size=0.5) +
	geom_jitter(aes(IN\$annot,IN\$lambda),color=\"gray30\",size=0.15) +";
if ($color_mapping !~ /none/i) {
	print R "
	scale_fill_manual(values = c($color_mapping)) +";
}
print R "
	scale_x_discrete(limits = rev(levels(IN\$V2))) +
	guides(fill=FALSE,colour=FALSE) +
	xlab(\"Annot\") +
	ylab(\"Lambda\") +
	coord_flip()
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
