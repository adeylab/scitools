package sci_commands::met_NOMeextract;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("met_metextract");

sub met_NOMeextract {

@ARGV = @_;

$gzip = "gzip"; #DEFAULT=gzip

# To enable "hg38", "hg19", and "mm10" shorcut usage, ensure all files are present
%GENOMES = (
   "hg19" => "/home/groups/oroaklab/refs/hg19/bismark",
   "hg38" => "/home/groups/oroaklab/refs/hg38/bismark",
   "mm10" => "/home/groups/oroaklab/refs/mm10/bismark"
);

use Getopt::Std; %opt = ();
use Cwd;
use Data::Dumper;

getopts("b:O:T:o:", \%opt);
$die2 = "
scitools met_NOMeextract [options] [any_C File]

sciNOMe generation of filtered GC, CG, and CH sorted bed files from gzipped any_C_context file.
[any_C File]                            =               File generated from scitools bam-methylationextraction call.

1. Takes a Bismark any_C file of methylation extraction, subsets to HCG [H= A,T,C] HCH and GCH
2. Then sorts the output and saves (for faster subsequent runs) as <prefix>.[CG/CH/GC].sorted.bed.
                NOTE: Temporary files generated through the sort function will be written to the current directory by default.
                Change the directory of temporary file output with the -T option.

                CG/GC/CH sorted Bed Files are modified bed files:
                <CHROM> <START> <END> <methylated C count> <non-methylated C count> <Percent Methylation at C> <Barcode>

Options:
   -O   [STR]   Output Prefix.
                                (Default: all text within last \"/\" and first \"\.\"\ of [Input File])
   -T   [STR]   Temporary Directory for sort files. (Default = Current Directory)
   -o   [STR]   Output Directory
                                (Default: Current working directory)
   -b   [STR]   REQUIRED: Bismark Reference Genome, full path to folder. Can also supply hg19, hg38, or mm10 as character arguments.
";

if (!defined $GENOMES{$opt{'b'}}) {die "\n\nERROR: Genome $ARGV[0] is not a proper genome reference! Exiting!\n"} else {$genome = $GENOMES{$opt{'b'}}};

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'T'}) {$opt{'T'} = getcwd()};
if (!defined $opt{'o'}) {$opt{'o'} = getcwd()};
if (!defined $opt{'O'}) {$opt{'O'}=$ARGV[1]; @O = split(/\./,$opt{'O'}); $opt{'O'}=$O[0]};

sub extract_chromosome_name {
    ## Bowtie extracts the first string after the inition > in the FASTA file, so we are doing this as well
    my $fasta_header = shift;
    if ($fasta_header =~ s/^>//){
    my ($chromosome_name) = split (/\s+/,$fasta_header);
    return $chromosome_name;
      }
      else{
    die "The specified chromosome ($fasta_header) file doesn't seem to be in FASTA format as required!\n";
      }
  }

sub read_genome_into_memory{
  $genome_folder = $genome;
  ## reading in and storing the specified genome in the %chromosomes hash
  chdir ($genome_folder) or die "Can't move to $genome_folder: $!";
  warn "Now reading in and storing sequence information of the genome specified in: $genome_folder\n\n";

  my @chromosome_filenames =  <*.fa>;

  ### if there aren't any genomic files with the extension .fa we will look for files with the extension .fa.gz
  unless (@chromosome_filenames){
      warn "Couldn't find files ending in .fa, trying .fa.gz instead\n";
      @chromosome_filenames =  <*.fa.gz>;
  }

  ### if there aren't any genomic files with the extension .fa or .fq.gz we will look for files with the extension .fasta
  unless (@chromosome_filenames){
      @chromosome_filenames =  <*.fasta>;
  }
  ### if there aren't any genomic files with the extension .fa or .fa.gz or .fasta we will look for files with the extension .fasta.gz
  unless (@chromosome_filenames){
      @chromosome_filenames =  <*.fasta.gz>;
  }

  unless (@chromosome_filenames){
    die "The specified genome folder $genome_folder does not contain any sequence files in FastA format (with .fa, .fa.gz, .fasta or .fasta.gz file extensions)\n";
  }

  foreach my $chromosome_filename (@chromosome_filenames){

    # skipping the tophat entire mouse genome fasta file
    next if ($chromosome_filename eq 'Mus_musculus.NCBIM37.fa');

    if( $chromosome_filename =~ /\.gz$/){
  open (CHR_IN,"gunzip -c $chromosome_filename |") or die "Failed to read from sequence file $chromosome_filename $!\n";
    }
    else{
  open (CHR_IN,$chromosome_filename) or die "Failed to read from sequence file $chromosome_filename $!\n";
    }

    ### first line needs to be a fastA header
    my $first_line = <CHR_IN>;
    chomp $first_line;
    $first_line =~ s/\r//; # removing /r carriage returns

    ### Extracting chromosome name from the FastA header
    my $chromosome_name = extract_chromosome_name($first_line);

    my $sequence;
    while (<CHR_IN>){
      chomp;
      $_ =~ s/\r//; # removing /r carriage returns

      if ($_ =~ /^>/){
  ### storing the previous chromosome in the %chromosomes hash, only relevant for Multi-Fasta-Files (MFA)
  if (exists $chromosomes{$chromosome_name}){
    warn "chr $chromosome_name (",length $sequence ," bp)\n";
    die "Exiting because chromosome name already exists. Please make sure all chromosomes have a unique name!\n";
  }
  else {
    if (length($sequence) == 0){
      warn "Chromosome $chromosome_name in the multi-fasta file $chromosome_filename did not contain any sequence information!\n";
    }
    warn "chr $chromosome_name (",length $sequence ," bp)\n";
    $chromosomes{$chromosome_name} = $sequence;
    $processed{$chromosome_name} = 0; # processed chromosomes will be set to 1 later to allow a record of which chromosome has been processed
  }
  ### resetting the sequence variable
  $sequence = '';
  ### setting new chromosome name
  $chromosome_name = extract_chromosome_name($_);
      }
      else{
  $sequence .= uc$_;
      }
    }

    if (exists $chromosomes{$chromosome_name}){
      warn "chr $chromosome_name (",length $sequence ," bp)\t";
      die "Exiting because chromosome name already exists. Please make sure all chromosomes have a unique name.\n";
    }
    else{
      if (length($sequence) == 0){
  warn "Chromosome $chromosome_name in the file $chromosome_filename did not contain any sequence information!\n";
      }
      warn "chr $chromosome_name (",length $sequence ," bp)\n";
      $chromosomes{$chromosome_name} = $sequence;
      $processed{$chromosome_name} = 0; # processed chromosomes will be set to 1 later to allow a record of which chromosome has been processed
    }
  }
  warn "\n";
  chdir $parent_dir or die "Failed to move to directory $parent_dir\n";
}

sub handle_filehandles{
    my $my_chr = shift;
    my $CG_OUT = "$opt{'o'}/$opt{'O'}.HCG.NOMe.bed";
    my $CH_OUT = "$opt{'o'}/$opt{'O'}.HCH.NOMe.bed";
    my $GC_OUT = "$opt{'o'}/$opt{'O'}.GCH.NOMe.bed";

    open (CG,">", "$CG_OUT");
    open (CH,">", "$CH_OUT");
    open (GC,">", "$GC_OUT");

    warn ">>> Writing genome-wide HCH bed file to: $CH_OUT <<<\n\n";
    warn ">>> Writing genome-wide HCG bed file to: $CG_OUT <<<\n\n";
    warn ">>> Writing genome-wide GCH bed file to: $GC_OUT <<<\n\n";
}

sub generate_genome_wide_cytosine_report {
#intial setup for function variables
    warn  "="x78,"\n";
    warn "Methylation information will now be written into a genome-wide cytosine report\n";
    warn  "="x78,"\n\n";
    sleep (2);

    if ($ARGV[0] =~ /gz$/){
    open (IN,"$gzip -dc $opt{'o'}/$ARGV[0] |") or die "Failed to read from gzipped file $ARGV[0]: $!\n";
    } else {
    open (IN, "$opt{'o'}/$ARGV[0]") or die "Failed to read from file $ARGV[0]: $!\n";
    }

    handle_filehandles();

    while (<IN>){
#reading in a line and setting variables
        chomp;
        #modified to read in anyC file instead
        my ($readID,$metcall,$chrom,$start,$methID,$readstart,$readend,$readdirection) = (split /\t/);
        my @R = split(/:/,$readID);
        $cellID=$R[0];
        if ($metcall =~ /\+/){$meth = 1;$nonmeth = 0;
        } elsif ($metcall =~ /\-/){$nonmeth = 1;$meth = 0;
        } else {next;}

        my $tri_nt;
        my $context;
 
#gather chromosome context for working chrom

        if ($readdirection =~ /\+/){    # C on forward strand
        $tri_nt = substr ($chromosomes{$chrom},($start-2),3);   # positions are 0-based!
        } elsif ($readdirection =~ /\-/){ # C on reverse strand
            if ($start-3 < 0) {next; # VB 28 09 2017 ; also, it would be more efficient to put this out of the loop...
          } else{
          $tri_nt = substr ($chromosomes{$chrom},($start-2),3);   # positions are 0-based!
          }
          $tri_nt = reverse $tri_nt;
          $tri_nt =~ tr/ACTG/TGAC/;
        }
        next if (length$tri_nt < 3); # trinucleotide sequence could not be extracted
#Assigning cytosine context for each site
        ### determining cytosine context
        if ($tri_nt =~ /^[ATC]CG$/){$context = 'HCG';
        print CG join("\t",$chrom,$start,$start,$meth,$nonmeth,$cellID),"\n";
        } elsif ($tri_nt =~ /^GC[ATC]$/){$context = 'GCH';
        print GC join("\t",$chrom,$start,$start,$meth,$nonmeth,$cellID),"\n";
        } elsif ($tri_nt =~ /^[ATC]C[ATC]$/){$context = 'HCH';
        print CH join("\t",$chrom,$start,$start,$meth,$nonmeth,$cellID),"\n";
        } else{next;}
      }
      close CG;
      close GC;
      close CH;
      }

read_genome_into_memory();
generate_genome_wide_cytosine_report();

$sort_cg="sort -T $opt{'T'} -k1,1 -k2,2n -k6,6 $CG_OUT > $opt{'o'}/$opt{'O'}.HCG.NOMe.sorted.bed";
$gzip_cg="gzip $opt{'o'}/$opt{'O'}.HCG.NOMe.sorted.bed";
$rm_cg="rm -f $CG_OUT";

$sort_ch="sort -T $opt{'T'} -k1,1 -k2,2n -k6,6 $CH_OUT > $opt{'o'}/$opt{'O'}.HCH.NOMe.sorted.bed";
$gzip_ch="gzip $opt{'o'}/$opt{'O'}.HCH.NOMe.sorted.bed";
$rm_ch="rm -f $CH_OUT";

$sort_gc="sort -T $opt{'T'} -k1,1 -k2,2n -k6,6 $GC_OUT > $opt{'o'}/$opt{'O'}.GCH.NOMe.sorted.bed";
$gzip_gc="gzip $opt{'o'}/$opt{'O'}.GCH.NOMe.sorted.bed";
$rm_gc="rm -f $GC_OUT";

system($sort_cg);
system($sort_gc);
system($sort_ch);

system($gzip_cg);
system($gzip_gc);
system($gzip_ch);

system($rm_cg);
system($rm_gc);
system($rm_ch);

}
1;