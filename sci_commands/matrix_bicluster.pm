package sci_commands::matrix_bicluster;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_bicluster");

sub matrix_bicluster {

# Defaults
($width,$height) = (12,12);
$res = 600;
$gradient_def = "BuRd";
$imageType = "png";

@ARGV = @_;
getopts("O:d:t:s:A:a:C:c:r:N:n:G:L:R:VvPXpS:", \%opt);

$die2 = "
scitools matrix-bicluster [options] [input matrix]
   or    bicluster-matrix
   or    bicluster

Biclusters matrix and plots with heatmap2 and complexheatmap

Options:
   -O   [STR]   Output (default is matrix file prefix)
   -d   [INT,INT] width,height (inches, def = $width,$height)
   -t   [pdf/png] output image type
   -s   [INT]   Resolution for png (def = $res)
   -A   [STR]   Annotation file (for column coloring)
   -C   [STR]   Color coding file (annot (tab) #hexColor)
   -c   [STR]   Color coding string
                  Annot=#hexColor,Annot2=#hexColor
   -r   [STR]   Row annotation file
   -N   [STR]   Color coding file for rows
   -n   [STR]   Color coding string for rows
   -G   [GRD]   Color gradient (def = $gradient_def)
                  For all available gradients, run 'scitools gradient'
   -R   [STR]   Rscript call (def = $Rscript)
   -V           Do not cluster columns
   -v           Do not cluster rows
   -S 	      Based on silhuette plot do row clustering with optimal number of clusters
   -P 			Do not show column names
   -p 			Do not show row names
   -L 			List of names that are part of the row annotation (eg tf name) it will plot them in the heatmap eg Rfx,Jun will mark all the tfs with this name
   -X           Do not delete intermediate files (def = delete)
   
";

if (!defined $ARGV[0]) {die $die2};
if (defined $opt{'C'} && defined $opt{'c'}) {die "\nSpecify either a color string (-c) or a color coding file (-C), not both!\n$die2"};
if (defined $opt{'N'} && defined $opt{'n'}) {die "\nSpecify either a color string (-n) or a color coding file (-N), not both!\n$die2"};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
$opt{'O'} =~ s/\.matrix$//;

if (!defined $opt{'G'}) {$opt{'G'} = $gradient_def};
$gradient_function = get_gradient($opt{'G'});

if (defined $opt{'d'}) {($width,$height) = split(/,/, $opt{'d'})};
if (defined $opt{'t'}) {if ($opt{'t'} =~ /(pdf|png)/) {$imageType = $opt{'t'}} else {die "\nERROR: Must specift pdf OR png for option -t.\n"}};
if (defined $opt{'s'}) {$res = $opt{'s'}};

if (defined $opt{'R'}) {$Rscript = $opt{'R'}};

if (defined $opt{'r'}) {read_annot($opt{'r'})};
%ROWID_annot = %CELLID_annot; %CELLID_annot = ();

if (defined $opt{'A'}) {read_annot($opt{'A'})};
if (defined $opt{'a'}) {
	@ANNOT_LIST = split(/,/, $opt{'a'});
	foreach $annot (@ANNOT_LIST) {
		$ANNOT_include{$annot} = 1;
	}
}

if (defined $opt{'N'}) {read_color_file($opt{'N'})};
if (defined $opt{'n'}) {read_color_string($opt{'n'})};
%ROW_ANNOT_color = %ANNOT_color; %ANNOT_color = ();
$color_mapping_row=$color_mapping;

if (defined $opt{'C'}) {read_color_file($opt{'C'})};
if (defined $opt{'c'}) {read_color_string($opt{'c'})};
# make col or row vector
open MATRIX, "$ARGV[0]";
$h = <MATRIX>; chomp $h; @H = split(/\t/, $h);
if (defined $opt{'A'} && (defined $opt{'c'} || defined $opt{'C'})) {
	open COLC, ">$opt{'O'}.col_colors.list";
	for ($colID = 0; $colID < @H; $colID++) {
		if (defined $ANNOT_color{$CELLID_annot{$H[$colID]}}) {
			print COLC "$ANNOT_color{$CELLID_annot{$H[$colID]}}\n";
		} else {
			print COLC "gray50\n";
		}
	}
	close COLC;
}
if (defined $opt{'r'} && (defined $opt{'n'} || defined $opt{'N'})) {
	open ROWC, ">$opt{'O'}.row_colors.list";
	while ($l = <MATRIX>) {
		chomp $l;
		@P = split(/\t/, $l);
		$rowID = shift(@P);
		if (defined $ROW_ANNOT_color{$ROWID_annot{$rowID}}) {
			print ROWC "$ROW_ANNOT_color{$ROWID_annot{$rowID}}\n";
		} else {
			print ROWC "gray50\n";
		}
	}
	close ROWC;
}
close MATRIX;

# make r script
open R, ">$opt{'O'}.heatmap2.r";

print R "
library(gplots)
library(grid)
#second type of heatmap
library(ComplexHeatmap)
library(circlize)
library(cluster)
set.seed(123)


$gradient_function
IN <- read.table(\"$ARGV[0]\");
IN_ch<-read.delim(\"$ARGV[0]\")\n";


if ($imageType =~ /pdf/i) {
	print R "$imageType(\"$opt{'O'}.silhouette.$imageType\",width=$width,height=$height)\n";
} else {
	print R "$imageType(\"$opt{'O'}.silhouette.$imageType\",width=$width,height=$height,units=\"in\",res=$res)\n";
}
if (defined $opt{'S'}){
print R "
k.max <- 15
sil <- rep(0, k.max)
for(i in 2:k.max){
pr = pam(IN_ch,k=i)
sil[i] <-mean(silhouette(pr)[, \"sil_width\"])
}

# Plot the  average silhouette width
plot(1:k.max, sil, type = \"b\", pch = 19,
     frame = FALSE, xlab = \"Number of clusters k\")
abline(v = which.max(sil), lty = 2)
dev.off()
kmaval=which.max(sil)
pa=pam(IN_ch,k=kmaval)
write.table(pa\$clustering,rownames=T,colnames=T,sep=\"t\",quote=F)
\n";


}

#complexheatmap implementation


#complex heatmap with col annot
if (defined $opt{'A'}) {
print R "annot_ch<-read.table(\"$opt{'A'}\",header=F) 
#order annot file in matrix order
annot_col<-annot_ch[match(colnames(IN_ch), annot_ch\$V1),]
\n"};


if (defined $opt{'A'} && (defined $opt{'c'} || defined $opt{'C'})) {
print R "color_ch<-list(Type=c($color_mapping))
ha_col<-HeatmapAnnotation(Type = annot_col\$V2,col=color_ch)\n"
} elsif(defined $opt{'A'} && (!defined $opt{'c'} || !defined $opt{'C'})){print R "ha_col<-HeatmapAnnotation(Type = annot_col\$V2)\n"};


#now with rows

if (defined $opt{'r'}) {
print R "
annot_ch_row<-read.table(\"$opt{'r'}\",header=F) 
#order annot file in matrix order
annot_row<-annot_ch_row[,match(rownames(IN_ch), annot_ch_row\$V1)]
\n"};


if (defined $opt{'r'} && (defined $opt{'n'} || defined $opt{'N'})) {
print R "
color_ch_row<-list(Type=c($color_mapping_row))
ha_row<-HeatmapAnnotation(Type = annot_ch_row\$V2,col=color_ch_row,which = \"row\")\n"
};


if (defined $opt{'L'})
{

@input=split(',',$opt{'L'});
@qinput=map{"'$_'"} @input;
$string="c(".join(',',@qinput).")";
print R "
names <- $string
aloc<-sapply(names,function(x) grep(x,rownames(IN_ch)))
aloc<-unlist(aloc)
alabel<-rownames(IN_ch)[aloc]
labels<-alabel
subset<-aloc
ha_row_script<-rowAnnotation(link = row_anno_link(at = subset, labels = labels),width = unit(1, \"cm\") + max_text_width(labels))\n";

}







if ($imageType =~ /pdf/i) {
	print R "$imageType(\"$opt{'O'}.complexheatmap.$imageType\",width=$width,height=$height)\n";
} else {
	print R "$imageType(\"$opt{'O'}.complexheatmap.$imageType\",width=$width,height=$height,units=\"in\",res=$res)\n";
}


print R "H <- Heatmap(IN_ch";
if (defined $opt{'S'}) {
	print R " ,split = paste0(\"pam\", pa\$clustering)";
}

if (defined $opt{'A'}) {
	print R " ,bottom_annotation=ha_col";
}

if (defined $opt{'v'}) {
	print R ",cluster_rows=FALSE";
}
if (defined $opt{'V'}) {
	print R ",cluster_columns=FALSE";
}

if (defined $opt{'p'}) {
	print R ",show_row_names=F";
}
if (defined $opt{'P'}) {
	print R ",show_column_names=F";
}

print R ")";


if (defined $opt{'r'} && (defined $opt{'n'} || defined $opt{'N'})) {
	print R "
	H=H+ha_row";
}

if (defined $opt{'L'})
{
	print R "
	H=H+ha_row_script";	
}



print R "\nH\n dev.off()\n";






#heatmap2 implementation

if (defined $opt{'r'} && (defined $opt{'n'} || defined $opt{'N'})) {print R "ROW <- read.delim(\"$opt{'O'}.row_colors.list\",header= F)\n"};
if (defined $opt{'A'} && (defined $opt{'c'} || defined $opt{'C'})) {print R "COL <- read.delim(\"$opt{'O'}.col_colors.list\",header= F)\n"};


if ($imageType =~ /pdf/i) {
	print R "$imageType(\"$opt{'O'}.heatmap2.$imageType\",width=$width,height=$height)\n";
} else {
	print R "$imageType(\"$opt{'O'}.heatmap2.$imageType\",width=$width,height=$height,units=\"in\",res=$res)\n";
}

print R "HM2 <- heatmap.2(as.matrix(IN),
	trace=\"none\",
	#col=gradient_funct(99),
	tracecol=\"black\"";

if (defined $opt{'r'} && (defined $opt{'n'} || defined $opt{'N'}) && defined $opt{'A'} && (defined $opt{'c'} || defined $opt{'C'})) {
	print R ",
	ColSideColors=as.vector(COL\$V1),
	RowSideColors=as.vector(ROW\$V1)";
} elsif (defined $opt{'r'} && (defined $opt{'n'} || defined $opt{'N'})) {
	print R ",
	RowSideColors=as.vector(ROW\$V1)";
} elsif (defined $opt{'A'} && (defined $opt{'c'} || defined $opt{'C'})) {
	print R ",
	ColSideColors=as.vector(COL\$V1)";
}

if (defined $opt{'v'}) {
	print R ",
	Rowv=FALSE";
}
if (defined $opt{'V'}) {
	print R ",
	Colv=FALSE";
}

print R ")\ndev.off()\n";





close R;

system("Rscript $opt{'O'}.heatmap2.r");

if (!defined $opt{'X'}) {
	if (defined $opt{'r'} && (defined $opt{'n'} || defined $opt{'N'})) {system("rm -f $opt{'O'}.row_colors.list")};
	if (defined $opt{'A'} && (defined $opt{'c'} || defined $opt{'C'})) {system("rm -f $opt{'O'}.col_colors.list")};
	system("rm -f $opt{'O'}.heatmap2.r");
}

}
1;