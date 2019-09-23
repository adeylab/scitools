package sci_commands::plot_factors;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("plot_factors");

sub plot_factors {

@ARGV = @_;
getopts("O:XR:", \%opt);

$die2 = "
scitools plot-factors [options] [factors.txt file]

plot-factors  plots  k-range vs error txt

Options:
   -O   [STR]   Output prefix (default is [input].plot.png and [input].plot.pdf)
   -X           Retain intermediate files (def = delete)
   -R   [STR]   Rscript call (def = $Rscript)

Note: Requires ggplot package 
";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix$//};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};

open R, ">$opt{'O'}.plot.factors.r";
print R "
#read packages
library(ggplot2)
#input
k_input<-read.table(file=\"$ARGV[0]\",header=TRUE)

PLT<-ggplot(data=k_input,aes(x=k,y=recon.err)) + theme_bw() +
	geom_point(shape=16) + geom_line() +
	xlab(\"Number of factors\") + ylab(\"Reconstruction error\") +
	theme(strip.background=element_rect(fill=\"transparent\")) +
	theme(strip.background=element_blank(),
		panel.grid=element_blank(),
		panel.background=element_blank(),
		plot.background=element_blank())
ggsave(plot=PLT,filename=\"$opt{'O'}.plot.png\",width=5,height=4,dpi=900)
ggsave(plot=PLT,filename=\"$opt{'O'}.plot.pdf\",width=5,height=4);
";

close R;

system("$Rscript $opt{'O'}.plot.factors.r");

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.plot.factors.r");
}

}
1;
