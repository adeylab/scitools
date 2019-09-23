scitools Version: 0.2.0
Adey Lab (www.adeylab.org, www.github.com/adeylab/scitools)

WARNING: scitools is still in a development stage - a number of
commands are in active development.

v0.2.0 update: This version contains a number of tools that were
internally used and moved here for public release. Note: many tools
require addiitonal dependencies not listed in the dependencies
function. This is being developed now and in progress. This only
applies to new functions, most of which have dependencies that are
self-explanatory (e.g. cistopic requires the R cistopic library),
and others have the dependencies in the individual function's
help text when called with no argumants.

v0.2.1 will have an updated dependencies and new document for the
recommended workflow of analysis.

SETUP: To set up scitools, download the repository to the
location you want the executable to be. Keep all index files
and folders within the scitools directory. We recommend adding
the scitools execulable to your path. We also recommend editing the
scitools.cfg file to add in command-line callable executables used
by scitools, defaults assume they are all in your path. Once it is
set up, run 'scitools dependencies' which will check for external
software that scitools relies on.

DESCRIPTION: scitools is a set of scripts designed for working with
single-cell combinatorial indexing data. It includes tools to go from
fastq files (after bcl2fastq) to a processed dataset. It includes a
number of external tools and R packages that are called by scitools.
If scitools is used, be sure to cite those tools!

scitools commands are in the form of [class]-[operation]. Most can
be specified in reverse order, or if the operation is unique to the
class of files or the class type can be determined by the files in
the arguments then just the operation name can be specified.
Run 'scitools list' for a list of commands, or
'scitools list [keyword]' to list commands that include the keyword,
e.g. 'scitools list fastq'

DEPENDENCIES: (command-line callable, can be specified in options)
Run scitools dependencies to check software and R requirements.
The default executables can be altered in the scitools.cfg file.
A user may also have their own personal scitools.cfg file in their
home directory to override the use of the config file that is present
with the scitools executable.

Executables:
   gzip         For gzipped fastq files. Default: gzip & zcat
                (these are hardoded at the beginning of the scitools
                 code and not available as options)
   samtools     Bam-related commands. Default: samtools
   bedtools     Bed-related commands. Default: bedtools
   Rscript      Numerous R operations. Default: Rscript
   bwa          For alignment only. Default: bwa
   macs2        For atac-callpeak only. Default: macs2
   scitools     Can call itself. Default: scitools
   (bismark)    For sci-MET
   (bowtie2)    For sci-MET / used by bismark

R packages:
   ggplot2         For plotting commands
   svd             For Latent Semantic Indexing (LSI)
   Rtsne           For tSNE visualization
   methods         For PCA
   dbscan          For density-based clustering
   irlba           For fast SVD (recommended)

File Types used by scitools:

   fastq    Standard fastq files, input can be gzipped or not, output fastq
               files will always be gzipped.
   
   bam      Bam files. For scitools, it is expected the barcode sequence
               or cellID is in the read name ([barcode]:[read_number])
               Sam files are not supported to encourage (force) space saving
			   
   annot    Annotation file: tab delimited, 2 columns
               Col 1 = cellID (barcode), could also be feature name
               Col 2 = annotation information (e.g. experimental condition)
            values file A special case of an annotaiton file where the
               annotation is a continuous variable and not discrete.
			   
   matrix   Matrix file, tab delimited
               Row 1 = CellIDs, has 1 less column than all other rows
               Rows = field 1 is name, then 1 field for each cellID
			   
   sparseMatrix   Sparse formatting for a matrix - triplet format with an
               associated cols and rows file that go with a values file.
			   
   dims     Dimensions file, tab delimited
               Col 1 = cellID
               Col 2-N = dimensions for each cell
            range    Range format is not a file type but a way to input a
               range of values (e.g. which dimensions to use from a dims file)
               This is a comma separated list and can include dashes for a
               range of values (e.g. 1,3-6,8-10,13)
			   
   colors   A file that has the annotaiton then a tab and an adssociated
               color for convenience during plotting.
