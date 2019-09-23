package sci_commands::plot_project_dims;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("plot_project_dims");

sub plot_project_dims {

@ARGV = @_;
# Defaults

$theme = "Reg";
$gradient_def = "BuG90Rd";

$panel_pass_color = "black";
$width = 5;
$height = 4;

getopts("O:A:a:C:c:R:T:V:M:Xs:w:h:", \%opt);

$die2 = "
scitools plot-project-dims [options] [dimensions file]

Project cell density onto individual dimensions, Useful for linear dim reductions


Options - general:
   -O   [STR]   Output prefix (default is dims file 1 prefix)
   -T   [STR]   Theme: (def = $theme)
                  Clean = no axis lines (for tSNE)
                  Reg = regular, include axes and grid (for PCA)
   -w   [FLT]   Plot width (inches, def = $width)
   -h   [FLT]   Plot height (inches, def = $height)

Plotting by annotations:
   -A   [STR]   Annotation file (to color code points)
   -a   [STR]   Comma separated list of annotations to include in plot
                  (requires -A to be specified)
   -C   [STR]   Color coding file (annot (tab) #hexColor)
   -c   [STR]   Color coding string
                  Annot=#hexColor,Annot2=#hexColor

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
if (defined $opt{'h'}) {$height = $opt{'h'}};
if (defined $opt{'w'}) {$width = $opt{'w'}};
$opt{'O'} =~ s/\.txt$//; $opt{'O'} =~ s/\.dims$//;



$filename = $ARGV[0];
open($fh, '<:encoding(UTF-8)', $filename)
  or die "Could not open file '$filename' $!";
 
 
while ($row = <$fh>) {
  chomp $row;
  @row_array = split(/\t/, $row);
  $cellID=shift @row_array;
  for ($i = 0; $i < @row_array; $i++) 
  {$CELLID_DIMS{$cellID}[$i]=@row_array[$i];}
  $counter=scalar @row_array;
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

for (my $i=0; $i <= $counter; $i++) {
open DATA, ">$opt{'O'}.dim$i.plot.txt";
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

			print DATA "$cellID\t$annot\t$CELLID_DIMS{$cellID}[$i]\n";
		}
	}
close DATA;	
 

open R, ">$opt{'O'}.plot.r";

print R "
library(ggplot2)
IN<-read.table(\"$opt{'O'}.dim$i.plot.txt\")";


if (!defined $opt{'c'} && !defined $opt{'C'} && !defined $opt{'A'}) { # no special mode specified
	print R "
	PLT<-ggplot(IN,aes(V3,color=\"lightsteelblue4\") + geom_density() +"; 
} else { # color
	print R "
	PLT<-ggplot(IN,aes(V3,color=V2)) + geom_density() +"; 
}

if ($color_mapping !~ /none/i) {
	print R "
	scale_colour_manual(values = c($color_mapping)) +";
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
ggsave(plot=PLT,filename=\"$opt{'O'}.dim$i.plot.png\",width=$width,height=$height,dpi=900)
ggsave(plot=PLT,filename=\"$opt{'O'}.dim$i.plot.pdf\",width=$width,height=$height)
";

close R;


system("$Rscript $opt{'O'}.plot.r");


if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.plot.r $opt{'O'}.dim$i.plot.txt");
}
}
}
1;
