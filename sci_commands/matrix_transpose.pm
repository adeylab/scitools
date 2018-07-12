package sci_commands::matrix_transpose;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_transpose");

sub matrix_transpose {

@ARGV = @_;
# Defaults

getopts("O:X:", \%opt);

$die2 = "
scitools matrix-transpose [options] [matrix]
   or    transpose

Options:
   -O   [STR]   Output prefix (def = [input].transposed.matrix)
   -X           Retain intermediate files (def = delete)
";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix$//};

#this is so much faster than the R so use this
$command="awk \'
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = \$i
    }
}
NF>p { p = NF }
END { 
    
     str=a[1,1]
        for(i=2; i<=NR; i++){
            str=str\"\\t\"a[i,1];
        }
     print str
    
    
    for(j=2; j<=p; j++) {
        str=\"row\"j-1\"\\t\"a[1,j]
        for(i=2; i<=NR; i++){
            str=str\"\\t\"a[i,j];
        }
        print str
    }
}\' $ARGV[0] >$opt{'O'}.transposed.matrix";

system($command);


open R, ">$opt{'O'}.transpose.r";
print R "

IN<-read.table(\"$ARGV[0]\",header=F)
n <- IN[,1]
tIN <- as.data.frame(t(IN[,-1]))
colnames(tIN)<-n
#rownames(tIN)<-1:length(rownames(tIN))
write.table(tIN,file=\"$opt{'O'}.transpose.matrix\",sep=\"\\t\",row.names=TRUE,col.names=TRUE,quote=FALSE)
";
#system("Rscript $opt{'O'}.transpose.r");

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.transpose.r");
}

}
1;
