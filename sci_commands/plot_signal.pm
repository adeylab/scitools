package sci_commands::plot_signal;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("plot_signal");

sub plot_signal {

@ARGV = @_;
# Defaults
$gradient_def = "WtPu";
$min = 0;
$max = 5;
$width_spec = "each=1";
$height=6;

getopts("O:r:G:R:Xh:w:V:H:", \%opt);

$die2 = "
scitools plot-signal [options] [(zscored) signal file]
   or    signal-plot

Options:
   -O   [STR]   Output prefix (default is signal file prefix)
   -r   [m,M]   Range: Min,Max values to impose for plotting (def = $min,$max)
   -G   [GRD]   Color gradient (def = $gradient_def)
                  For all available gradients, run 'scitools gradient'
   -h   [INS]   Height of plot (inches, def = $height)
   -w   [STR]   Width of plot (in, each=[for each annot] or all=[total width]
                  if a single number will assume total width def = $width_spec)
   -V   [STR]   Vertical line deliminators file (.vlines.txt)
   -H   [STR]   Horizontal line deliminators file (.hlines.txt)
   -R   [STR]   Rscript call (def = $Rscript)
   -X           Do not delete intermediate files (def = delete)

Note: Requires ggplot2 R package

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
if (!defined $opt{'G'}) {$opt{'G'} = $gradient_def};
$gradient_function = get_gradient($opt{'G'});
if (defined $opt{'r'}) {($min,$max) = split(/,/, $opt{'r'})};
if (defined $opt{'h'}) {$height = $opt{'h'}};
if (defined $opt{'w'}) {$width_spec = $opt{'w'}};
if ($width_spec =~ /=/) {
	($width_type,$width_value) = split(/=/, $width_spec);
} else {
	$width_type = "all"; $width_value = $width_spec;
}

open IN, "$ARGV[0]";
open OUT, ">$opt{'O'}.plot";
$h = <IN>; chomp $h; @H = split(/\t/, $h);
$win_pos = 0;
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	$winID = shift(@P);
	for ($i = 0; $i < @P; $i++) {
		if ($P[$i] < $min) {$P[$i] = $min};
		if ($P[$i] > $max) {$P[$i] = $max};
		print OUT "$win_pos\t$i\t$P[$i]\n";
	}
	$win_pos++;
} close OUT;

# determine plot width
%ANNOT_list = (); $annot_count = 0;
for ($i = 0; $i < @H; $i++) {
	($annot,$subWin) = split(/_window_/, $H[$i]);
	if (!defined $ANNOT_list{$annot}) {
		$ANNOT_list{$annot} = 1;
		$annot_count++;
	}
}
if ($width_type =~ /all/i) {
	$width = $width_value;
} elsif ($annot_count <= 1) {
	$width = $width_value+0.5;
} else {
	$width = ($width_value*$annot_count)+0.5;
}

open R, ">$opt{'O'}.plot.r";
print R "
library(ggplot2)
$gradient_function
IN<-read.table(\"$opt{'O'}.plot\")";

if (defined $opt{'V'}) {
	print R "
VL<-read.table(\"$opt{'V'}\")";
}

if (defined $opt{'H'}) {
	print R "
HL<-read.table(\"$opt{'H'}\")";
}

print R "
colnames(IN)<-c(\"row\",\"col\",\"val\")
PLT<-ggplot() + theme_bw() +
	geom_tile(aes(IN\$col,IN\$row,fill=IN\$val)) +";

if (defined $opt{'V'}) {
	print R "
	geom_vline(aes(xintercept = VL\$V1),linetype=\"dotted\",color=\"black\",size=0.25) +";
}

if (defined $opt{'H'}) {
	print R "
	geom_hline(aes(yintercept = HL\$V1),linetype=\"dashed\",color=\"black\",size=0.25) +";
}

print R "
	theme(axis.title.x=element_blank(),
		axis.text.x=element_blank(),
		axis.ticks.x=element_blank(),
		axis.title.y=element_blank(),
		axis.text.y=element_blank(),
		axis.ticks.y=element_blank(),
		panel.grid.major=element_blank(),
		panel.grid.minor=element_blank(),
		axis.line=element_blank(),
		panel.background=element_blank()) +
	labs(fill=\"Signal\") +
	scale_y_reverse(expand=c(0,0)) +
	scale_x_continuous(expand=c(0,0)) +
	scale_fill_gradientn(colours=gradient_funct(99))
ggsave(plot=PLT,filename=\"$opt{'O'}.plot.png\",width=$width,height=$height,dpi=900)
ggsave(plot=PLT,filename=\"$opt{'O'}.plot.pdf\",width=$width,height=$height)
";
close R;

system("$Rscript $opt{'O'}.plot.r");

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.plot.r $opt{'O'}.plot");
}

}
1;
