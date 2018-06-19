package sci_commands::matrix_merge;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_merge");

sub matrix_merge {

@ARGV = @_;
$joiner = "_";

getopts("O:U", \%opt);

$die2 = "
scitools matrix-merge [options] [matrix1] [matrix2] (optional matrix N etc...)
   or    merge-matrix

Will merge matrices if names of cells do not match (others can be further developed).
Features will be the intersect across the matrices.

Options:
   -O   [STR]   Output file name / prefix (def = matrix1 prefix w/ merge)
   -U           Print the union of all features (those not present in a
                 matrix will be set to 0)

";

if (!defined $ARGV[1]) {die $die2};

if (!defined $opt{'O'}) {
$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.counts.matrix$//; $opt{'O'} =~ s/\.binary.matrix$//; $opt{'O'} =~ s/\.matrix$//;
};
#
if ($ARGV[0] =~ /\.matrix$/) 
{
if ($ARGV[0] =~ /\.counts.matrix$/) {$opt{'O'} .= ".merged.counts";} elsif ($ARGV[0] =~ /\.binary.matrix$/) {$opt{'O'} .= ".merged.binary";} else {print "Warning: Type of matrix not in name\n"; $opt{'O'} .= ".merged";};
} else {print "ERROR: This is not a matrix, rename\n"; die};

#create a merged matrix with CELLID_FEATURE_value structure
print "reading matrix $ARGV[0]\n";
read_matrix($ARGV[0]);
%CELLID_FEATURE_merged_value=%CELLID_FEATURE_value;
#create an included rowid hash 
foreach $cellID (keys %CELLID_FEATURE_value) 
{
	foreach $rowID (keys %{$CELLID_FEATURE_value{$cellID}}) 
	{	
	$All_included_rowid{$rowID}=1;
	}
}
#create an all cells Array to preserve order. 
@MATRIX_ALL_COLNAMES=@MATRIX_COLNAMES;

for ($matrixID = 1; $matrixID < @ARGV; $matrixID++) 
{
	print "reading matrix $ARGV[$matrixID]\n";
	read_matrix($ARGV[$matrixID]);
	
#quick check for cellID overlaps (faster and can change to warning later)
	foreach $cellID (keys %CELLID_FEATURE_value) 
	{
		if (defined $CELLID_FEATURE_merged_value{$cellID})
		{
		print "ERROR: Overlapping cellIDs, renaming not yet developed\nThis is the cell name $cellID\n"; die
		}
	}
	#if no overlap continue
	foreach $cellID (keys %CELLID_FEATURE_value) 
	{
		if (!defined $CELLID_FEATURE_merged_value{$cellID})
		{
			foreach $rowID (keys %{$CELLID_FEATURE_value{$cellID}}) 
			{
				if (defined $All_included_rowid{$rowID})
				{# add values to the merged matrix
				$CELLID_FEATURE_merged_value{$cellID}{$rowID}=$CELLID_FEATURE_value{$cellID}{$rowID};
				}
				else
				{# if not defined, define and add values
				$CELLID_FEATURE_merged_value{$cellID}{$rowID}=$CELLID_FEATURE_value{$cellID}{$rowID};
				$All_included_rowid{$rowID}=1;
				}
			}	
		}
		else {print "ERROR: Overlapping cellIDs, renaming not yet developed\nThis is the cell name $cellID\n"; die};
		#can add part that handles overlapping IDs
		#need to consider two cases: 1. cellIDs are same because different capacity run on same cells, 2. CellIDs are same but different cells, need to rename
	}
	#add new cellnames
	push(@MATRIX_ALL_COLNAMES,@MATRIX_COLNAMES);
	
}

open OUT, ">$opt{'O'}.matrix";



print OUT join("\t",@MATRIX_ALL_COLNAMES)."\n";
foreach $rowID (sort keys %All_included_rowid) 
	{
	#UNION vs INTERSECT
	$in_all = 1;
	#OUTPUT hash
	@OUTPUT_line=();
	push(@OUTPUT_line,$rowID);
		foreach $cellID (@MATRIX_ALL_COLNAMES) 
		{
			if (defined $CELLID_FEATURE_merged_value{$cellID}{$rowID})
			{#if defined input value
			push(@OUTPUT_line,$CELLID_FEATURE_merged_value{$cellID}{$rowID});
			}
			else
			{#if not defined input  0
			push(@OUTPUT_line,0);
			$in_all = -1;
			}
		}
		if (defined $opt{'U'} || $in_all>0)
		{
			print OUT join("\t",@OUTPUT_line)."\n";
		}
	}
close OUT;

}
1;
