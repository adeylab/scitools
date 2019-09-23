package sci_commands::plot_plotly;
use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("plot_plotly");
sub plot_plotly {
@ARGV = @_;
# Defaults
getopts("O:R:XA:V:p:B:C:M:m:", \%opt);
$die2 = "
scitools plot_plotly [options] [2D or 3D dims file]
   or    plotly
This function behaves similarly to plot-dims. Will take in a 2D or 3D [umap|tsne|swne|pca] dims file.
Will generate a single html widget with interactive graph for all files supplied to plot on the dims points.
Currently accepted files:
Annotation files:                            Will generate a legend with each annotation.
Monocle monocle3_cells.txt:                  Will plot pseudotime value over cells as continuous color scale.
Value files:                                 Will plot values over cells as continuous color scale.
ChromVar deviation_scores.matrix             Will plot scores over cells as continuous color scale.
                                             Limited to 10 specified transcription factors.
Monocle monocle3_segmentvectors.txt output:  Will plot cds_monocle cell trajectories.

Options:
   -O   [STR]   Output Directory/Prefix (default is all [2D or 3D dims file].plotly)
   -R   [STR]   Rscript call (def = $Rscript)
   -A   [STR]   Comma-separated list of .annot files. Will generate a legend with each annotation.
   -V   [STR]   Comma-separated list of .val files. Will generate a legend item for each value.
   -p   [STR]   scitools cds_monocle generated monocle3_cells.txt file. Will generate a legend item for pseudotime.
   -M   [STR]   scitools chromvar generated deviation_scores.matrix file. Will read in transcription factor deviation scores file.
                  If specified, -m must also be specified.
   -m   [STR]   comma separated list of transcription factors. CAN ONLY BE LENGTH 10 OR LESS.   
                  If specified, -M must also be specified.
   -B   [STR]   scitools cds_monocle generated monocle3_segmentvectors.txt file. Will plot a cell trajectory line plot.
   -C   [STR]   scitools_rmdup generated complexity plot. Will plot read depth and unique percentage.
   -X           Retain intermediate files (Default = delete)
                  
";
if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0];  $opt{'O'} =~ s/\.dims$//;};
open R, ">$opt{'O'}.plotly.r";
print R "
library(plotly)
library(htmlwidgets)
library(RColorBrewer)
dims<-read.table(\"$ARGV[0]\",header=F)
if (ncol(dims)==2){
   colnames(dims)<-c(\"cellID\",\"x\",\"y\")
} else {
   colnames(dims)<-c(\"cellID\",\"x\",\"y\",\"z\")
}
write(paste(\"Read in dims file:\",\"$ARGV[0]\",sep=\"\\t\"),stderr())
write(paste(\"Dimensions Detected:\",as.character(ncol(dims)-1),sep=\"\\t\"),stderr())

";
$alternative_color_files=0;

#Read in given annotation files
if (defined $opt{'A'}) {
   #Rewrite comma separated list for R list format
   @annot_list = split(/,/,$opt{'A'});
   foreach (@annot_list) {$_ = "'$_'";}
   $opt{'A'}=join(', ', @annot_list );
   $alternative_color_files=1;
   print R "
   for (i in c($opt{'A'})){
   # read in file base name for listing annotation
   i_name<-basename(i)
   #strip the \".annot\" text
   i_name<-unlist(strsplit(i_name,\"\.annot\")[1])
   #read in the file
   working_df<-read.table(i,header=F)
   #set column names
   colnames(working_df)<-c(\"cellID\",paste(i_name))
   #generate a new data frame for each file read in
   assign(paste0(\"col_annot_\",i_name),working_df)

   write(paste(\"Read in annot file:\",i,sep=\"\\t\"),stderr())
   }";
};

#Read in given value files
if (defined $opt{'V'}) {
   #Rewrite comma separated list for R list format
   @val_list = split(/,/,$opt{'V'});
   foreach (@val_list) {$_ = "'$_'";}
   $opt{'V'}=join(', ', @val_list );
   $alternative_color_files=1;
   print R "
   for (i in c($opt{'V'})){
   # read in file base name for listing annotation
   i_name<-basename(i)
   #strip the \".val\" text
   i_name<-unlist(strsplit(i_name,\"\.val\")[1])
   #read in the file
   working_df<-read.table(i,header=F)
   #set column names
   colnames(working_df)<-c(\"cellID\",paste(i_name))
   #generate a new data frame for each file read in
   assign(paste0(\"col_val_\",i_name),working_df)
   
   write(paste(\"Read in value file:\",i,sep=\"\\t\"),stderr())
   }";
};

#Read in given monocle.cells.txt file for pseudotime
if (defined $opt{'p'}) {

   $alternative_color_files=1;
   print R "
   if(!grepl(\"monocle3_cells.txt\$\", \"$opt{'p'}\")) {
      stop(\"Error: ensure your -p option is file named monocle3_cells.txt\")
   }
   col_m_file<-read.table{\"$opt{'p'}\",header=T}
   col_m_file\$cellID<-row.names(col_m_file)
   col_m_file<-col_m_file[,c(18,17)]
   write(paste(\"Read in pseudotime file:\",\"$opt{'p'}\",sep=\"\\t\"),stderr())

   ";
};

#Read in given complexity.txt
if (defined $opt{'C'}) {
   $alternative_color_files=1;
   print R "
   if(!grepl(\"complexity.txt\$\", \"$opt{'C'}\")) {
      stop(\"Error: ensure your -C option is a complexity file.\")
   }
   col_c_file<-read.table(\"$opt{'C'}\",header=F, row.names=1)
   colnames(col_c_file)<-c(\"cellID\",\"total_reads\",\"uniq_reads\",\"percent_uniq\")
   write(paste(\"Read in complexity file:\",\"$opt{'C'}\",sep=\"\\t\"),stderr())

   ";
};

#Read in deviation_scores.txt and subset to only the list of TF given.
if (defined $opt{'M'} && defined $opt{'m'}) {
   $alternative_color_files=1;
   @tf_list = split(/,/,$opt{'m'});
   foreach (@tf_list) {$_ = "'$_'";}
   $opt{'m'}=join(', ', @tf_list );
   $alternative_color_files=1;
   print R "
   if(!grepl(\"deviation_scores.matrix\$\", \"$opt{'M'}\")) {
      stop(\"Error: ensure your -M option is a deviation_scores.matrix file\")
   }
   if(length(c($opt{'m'}))>10){
      stop(\"Error: Limit your transcription factor list to 10.\")
   }
   M_file<-read.table(\"$opt{'M'}\",header=T)

   tf_list<-c($opt{'m'})
   out<- c()
   for (x in tf_list){
   check_each <- rownames(M_file)[grepl(x,rownames(M_file),ignore.case=TRUE)]
   out<-c(out,check_each)
   }
   write(\"Found transcription factor matches:\",stderr())
   write(paste(out,sep=\"\\n\"),stderr())


   col_M_file<-subset(M_file,rownames(M_file) \%in\% out)
   col_M_file<-data.frame(t(col_M_file))
   col_M_file\$cellID<-rownames(col_M_file)
   write(paste(\"Read in chromvar file:\",\"$opt{'M'}\",sep=\"\\t\"),stderr())

   ";
};

#Read in given monocle3_segmentvectors.txt file for branchpoints
if (defined $opt{'B'}) {
   $alternative_color_files=1;
   print R "
   if(!grepl(\"monocle3_segmentvectors.txt\$\", \"$opt{'B'}\")) {
      stop(\"Error: ensure your -B option is file named monocle3_segmentvectors.txt\")
   }
   b_file<-read.table(\"$opt{'B'}\",header=T)
   b_file\$line_segment<-as.factor(b_file\$line_segment)
   write(paste(\"Read in monocle branchpoint file:\",\"$opt{'B'}\",sep=\"\\t\"),stderr())

   ";
};


#Loop in color files into a single data frame
if ($alternative_color_files==1) {
print R "
#list all files for coloring
color_frames<-ls(pattern=\"col_\")
#merge first file with dims
for (i in color_frames) {
   dims<-merge(dims,get(i),by=\"cellID\")
}
";
};

print R "
#make a plotly loop to go through files.
#Generate 3D Plot
if (\"z\" %in% colnames(dims)) {
   #initiate 3d plot
p<-plot_ly(type=\"scatter3d\")
p<- p %>% add_trace(p,mode=\"markers\",data=dims,x = ~x, y = ~y, z = ~z,color=~dims[,5],size=10,opacity=0.7,hovertext=~cellID)
for (i in 6:ncol(dims)) {
   if (sapply(dims[i],is.factor)) {
      p<- p %>% add_trace(mode=\"markers\",data=dims,x = ~x, y = ~y, z = ~z,color=~dims[,i],size=10,opacity=0.7,hovertext=~cellID,visible=\"legendonly\")
      } else {
      p<- p %>% add_trace(mode=\"markers\",data=dims,x = ~x, y = ~y, z = ~z,color=~dims[,i],size=10,opacity=0.7,hovertext=~cellID,visible=\"legendonly\")
      }
   }
if (exists(\"b_file\")) {
   p<- p %>% add_trace(p,data=b_file,mode=\"lines\",split=~line_segment,x=~prin_graph_dim_1,y=~prin_graph_dim_2,z=~prin_graph_dim_3,color=I(\"black\"),opacity=0.7,legendgroup=\"Trajectory\")
}

p<-p %>% layout(scene = list(xaxis = list(title = \'X\'),yaxis = list(title = \'Y\'), zaxis = list(title = \'Z\')))

} else {
   #initiate 2D plot
p<-plot_ly(type=\"scatter\")
p<- p %>% add_trace(p,mode=\"markers\",data=dims,x = ~x, y = ~y, color=~dims[,5],size=10,opacity=0.7,hovertext=~cellID)
for (i in 6:ncol(dims)) {
   if (sapply(dims[i],is.factor)) {
      p<- p %>% add_trace(mode=\"markers\",data=dims,x = ~x, y = ~y, color=~dims[,i],size=10,opacity=0.7,hovertext=~cellID,visible=\"legendonly\")
      } else {
      p<- p %>% add_trace(mode=\"markers\",data=dims,x = ~x, y = ~y, color=~dims[,i],size=10,opacity=0.7,hovertext=~cellID,visible=\"legendonly\")
      }
   }
if (exists(\"b_file\")) {
   p<- p %>% add_trace(p,data=b_file,mode=\"lines\",split=~line_segment,x=~prin_graph_dim_1,y=~prin_graph_dim_2,color=I(\"black\"),opacity=0.7,legendgroup=\"Trajectory\")
}

p<-p %>% layout(scene = list(xaxis = list(title = \'X\'),yaxis = list(title = \'Y\')))

}


htmlwidgets::saveWidget(widget=p,\"$opt{'O'}.html\", selfcontained = FALSE)

";
close R;
system("$Rscript $opt{'O'}.plotly.r");
if (!defined $opt{'X'}) {
    system("rm -f $opt{'O'}.plotly.r");
}
}
1;
