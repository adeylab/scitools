package sci_commands::cds_cistopicmatrix;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("cds_cistopicmatrix");

sub cds_cistopicmatrix {

@ARGV = @_;
# Defaults

getopts("O:R:T:X", \%opt);

$die2 = "
scitools cds-cistopicmatrix [options] [cistopic rds]

Prior to running, ensure you have run cistopic 

Options:
   -O    [STR]   Output Directory (default is [current working directory]/cicero_output)
   -R    [STR]   Rscript call (def = $Rscript)
   -T    [STR]   Type of cistopic matrix output. Chose method between \"Z-score\", \"Probability\", \"Imputed\", imputed is a large matrix with imputed peak cell values, defined as Z-score
   -X           Retain intermediate files (Default = delete)
                  
";


if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
$opt{'O'}=~s/.rds//;
if (!defined $opt{'T'}) {$opt{'T'} = "Z-score"};


open R, ">cds_matrix.r";

print R "
modelMatSelection <- function(
  object,
  target,
  method,
  all.regions=FALSE
){
  # Check info
  if (length(object\@selected.model) < 1){
    stop(\'Please, run selectModel() first.\')
  }
  
  if (target == \'cell\'){
    if (method == 'Z-score'){
      modelMat <- scale(object\@selected.model\$document_expects, center=TRUE, scale=TRUE)
    }
    
    else if (method == \'Probability\'){
      alpha <- object\@calc.params[[\'runModels\']]\$alpha/length(object\@selected.model\$topic_sums)
      modelMat <- apply(object\@selected.model\$document_sums, 2, function(x) {(x + alpha)/sum(x + alpha)})
    }
    else{
      stop(\'Incorrect method selected. Chose method between \"Z-score\" and \"Probability\".\')
    }
    colnames(modelMat) <- object\@cell.names
    rownames(modelMat) <- paste0(\'Topic\', 1:nrow(modelMat))
  }
  
  else if (target == \'region\'){
    if (!all.regions){
      if (length(object\@binarized.cisTopics) < 1){
        stop(\'Please, use binarizecisTopics() first for defining the high confidence regions for dimensionality reduction!\')
      }
      else {
        regions <- unique(unlist(lapply(object\@binarized.cisTopics, rownames)))
      }
    }
    
    topic.mat <- object\@selected.model\$topics
    
    if (method == \'NormTop\'){
      normalizedTopics <- topic.mat/(rowSums(topic.mat) + 1e-05)
      modelMat <- apply(normalizedTopics, 2, function(x) x * (log(x + 1e-05) - sum(log(x + 1e-05))/length(x)))
    }
    
    else if (method == \'Z-score\'){
      modelMat <- scale(object\@selected.model\$topics, center=TRUE, scale=TRUE)
    }
    
    else if (method == \'Probability\'){
      beta <- object\@calc.params[[\'runModels\']]\$beta
      topic.mat <- object\@selected.model\$topics
      modelMat <-  (topic.mat + beta)/rowSums(topic.mat + beta)
    }
    
    else{
      stop(\'Incorrect method selected. Chose \"NormTop\", \"Z-score\" and \"Probability\".\')
    }
    
    colnames(modelMat) <- object\@region.names
    rownames(modelMat) <- paste0(\'Topic\', 1:nrow(modelMat))
    
    if (!all.regions){
      modelMat <- modelMat[,regions]
    }
  }
  
  else{
    stop(\'Please, provide target=\"cell\" or \"region\".\')
  }
  
  return(modelMat)
}





cellwritecistopicmatrix <-function (object, method = \"Z-score\", ...)
{
    if (!\"cisTopic\" %in% installed.packages()) {
        stop(\"Please, install cistopic: \")
    }
    else {
        require(cisTopic)
    }
    if (length(object\@selected.model) < 1) {
        stop(\"Please, run selectModel() first.\")
    }
	
	if (method == \'Imputed\') {
	pred.matrix <- as.data.frame(predictiveDistribution(object))
    colnames(pred.matrix) <- object\@cell.names
	write.table(pred.matrix,file=paste0(\"$opt{'O'}\",\"_\",\"$opt{'T'}\",\".matrix\"),col.names=T,row.names=T,sep=\"\\t\",quote=F)
    }
	else {
    topic.mat <- modelMatSelection(object, \"cell\", method)
    rownames(topic.mat) <- paste(\"Topic\", seq(1, nrow(topic.mat)))
    colnames(topic.mat) <- object\@cell.names
	write.table(topic.mat,file=paste0(\"$opt{'O'}\",\"_\",\"$opt{'T'}\",\"topic.matrix\"),col.names=T,row.names=T,sep=\"\\t\",quote=F)
	}
}

cistopicObject<-readRDS(\"$ARGV[0]\")
cellwritecistopicmatrix(cistopicObject,method=\"$opt{'T'}\")
";




close R;
system("$Rscript cds_matrix.r");


if (!defined $opt{'X'}) {
    system("rm -f cds_matrix.r");
}
}
1;
