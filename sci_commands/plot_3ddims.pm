package sci_commands::plot_3ddims;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
use File::Basename;

@EXPORT = ("plot_3ddims");

sub plot_3ddims {

@ARGV = @_;
getopts("O:XR:A:T:t:a:p:o:ud:", \%opt);

$die2 = "
scitools plot-3ddims [options] [3D dims file]

Generate rotating 3D Plots as animated GIFs or webGL output.
[3D dims file]      can be either a *CDS.rds file, output from scitools cds_monocle function call with the -D 3 option
                    or a 3-dimension dims file generated through scitools matrix_[tsne|umap|swne] call with -D 3 option

Note: To have branch points plotted, you must use the CDS.rds file for [3D dims file].

Options:
   -o   [STR]   Output prefix (default is [input].3dplot.gif)
   -O   [STR]   Output directory (default is [3D dims file] directory)
   -A   [STR]   Annotation file to color points by
   -T   [STR]   ChromVar Deviation Scores Matrix (Required for -t option)
   -t   [STR]   Comma separated list of Transcription Factors to Output Z scored-colored points
   -X           Retain intermediate files (def = delete)
   -u           Will generate unannotated, uncolored gif.
   -a   [NUM]   0 to 1 range for point alpha values (Default = 0.6)
   -d   [NUM]   If defined -T and -t, can set color cut off for top -d percent accessible 
                cells rather than a color range.
   -p   [INT]   Point size (Default=3)
   -R   [STR]   Rscript call (def = $Rscript)
";

if (!defined $ARGV[0]) {die $die2};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};
if (!defined $opt{'a'}) {$opt{'a'}=0.6}
if (!defined $opt{'p'}) {$opt{'p'}=3}
if (!defined $opt{'O'}) {$opt{'O'}= dirname($ARGV[0])}

if ($ARGV[0] =~ /rds$/){
    if (!defined $opt{'o'}) {$opt{'o'} = $ARGV[0]; $opt{'o'} =~ s/\.rds$//; $opt{'o'}=basename($opt{'o'})};

    open R, ">$opt{'O'}/$opt{'o'}.3Dplot.r";
    print R "
        #Defined 3D Trajectory function
        library(monocle)
        library(rgl)
        library(RColorBrewer)

        dir=\"$opt{'O'}\"
        if(dir==\".\"){dir=getwd()}

        cds<-readRDS(file=\"$ARGV[0]\")
        cell_alpha=$opt{'a'}
        cell_size=$opt{'p'}
        backbone_segment_color=\"#000000\"
        gene_short_name <- NA
        sample_name <- NA
        sample_state <- pData(cds)\$State
        data_dim_1 <- NA
        data_dim_2 <- NA
        lib_info_with_pseudo <- pData(cds)
        if (is.null(cds\@dim_reduce_type)) {
            stop(\"Error: dimensionality not yet reduced. Please call reduceDimension() before calling this function.\")}

        reduced_dim_coords <- reducedDimK(cds)
        ica_space_df <- data.frame(Matrix::t(reduced_dim_coords[1:3,]))
        colnames(ica_space_df) <- c(\"prin_graph_dim_1\", \"prin_graph_dim_2\",\"prin_graph_dim_3\")
        ica_space_df\$sample_name <- row.names(ica_space_df)
        ica_space_df\$sample_state <- row.names(ica_space_df)

        dp_mst <- minSpanningTree(cds)
        if (is.null(dp_mst)) {
            stop(\"Error: Learned Graph not yet called. Please call learnGraph() before calling this function\")}
        edge_list <- as.data.frame(get.edgelist(dp_mst))
        colnames(edge_list) <- c(\"source\", \"target\")
        edge_df <- merge(ica_space_df, edge_list, by.x = \"sample_name\",by.y = \"source\", all = TRUE)
        edge_df <- plyr::rename(edge_df, c(prin_graph_dim_1 = \"source_prin_graph_dim_1\",prin_graph_dim_2 = \"source_prin_graph_dim_2\", prin_graph_dim_3 = \"source_prin_graph_dim_3\"))
        edge_df <- merge(edge_df, ica_space_df[, c(\"sample_name\",\"prin_graph_dim_1\", \"prin_graph_dim_2\", \"prin_graph_dim_3\")],by.x = \"target\", by.y = \"sample_name\", all = TRUE)
        edge_df <- plyr::rename(edge_df, c(prin_graph_dim_1 = \"target_prin_graph_dim_1\",prin_graph_dim_2 = \"target_prin_graph_dim_2\", prin_graph_dim_3 = \"target_prin_graph_dim_3\"))

        S_matrix <- reducedDimS(cds)
        data_df <- data.frame(t(S_matrix[1:3, ]))
        colnames(data_df) <- c(\"data_dim_1\", \"data_dim_2\", \"data_dim_3\")
        data_df\$sample_name <- row.names(data_df)
        lib_info_with_pseudo\$sample_name <- row.names(lib_info_with_pseudo)
        data_df <- merge(data_df, lib_info_with_pseudo, by = \"sample_name\")
        point_colors_df <- data.frame(sample_name = data_df\$sample_name,point_colors = \"darkgray\")

        point_colors_df\$point_alpha = cell_alpha
        point_colors_df\$point_alpha[is.na(point_colors_df\$point_colors)] = 0
        point_colors_df = merge(point_colors_df, data_df,by=\"sample_name\")
        lib_info_with_pseudo\$sample_name<-rownames(lib_info_with_pseudo)
        pdat<-lib_info_with_pseudo[,c(1,ncol(lib_info_with_pseudo))]
        point_colors_df = merge(point_colors_df,pdat,by=\"sample_name\")

    ";
    if (defined $opt{'u'}) {
    print R "
        #Generate Uncolored Rotating Gif

        open3d(windowRect = c(0, 50, 800, 800))
        segments3d(x=as.vector(t(edge_df[, c(3,7)])),y=as.vector(t(edge_df[, c(4,8)])),z=as.vector(t(edge_df[, c(5,9)])), lwd = 2, col = backbone_segment_color,line_antialias = TRUE)
        points3d(point_colors_df[, c(\"data_dim_1\", \"data_dim_2\", \"data_dim_3\")], size = cell_size, col = point_colors_df\$point_colors, alpha = point_colors_df\$point_alpha,point_antialias = TRUE)
        movie3d(spin3d(axis=c(0,0,1),rpm=5),duration=30,movie=\"Unannotated\",dir=dir,convert=TRUE)

    ";
    }
    if (defined $opt{'A'}){
    print R "
            # Generate Annotation Plot
            annot<-read.table(file=\"$opt{'A'}\",header=F)
            anno_name_list<-unlist(strsplit(\"$opt{'A'}\",\"/\"))
            anno_name_list2<-unlist(strsplit(anno_name_list[length(anno_name_list)],\"[.]\"))
            anno_name<-anno_name_list2[length(anno_name_list2)-1]
            colnames(annot)<-c(\"sample_name\",\"anno\")
            point_colors_df<-merge(point_colors_df,annot,by.x=\"sample_name\")
        colors<-colorRampPalette(brewer.pal(length(unique(point_colors_df\$anno)),\"Set2\"))(length(unique(point_colors_df\$anno)))
        point_colors_df<-transform(point_colors_df,colsplit=as.numeric(factor(point_colors_df\$anno)))
        point_colors_df\$point_colors<-colors[point_colors_df\$colsplit]

        open3d(windowRect = c(0, 50, 800, 800))
        segments3d(x=as.vector(t(edge_df[, c(3,7)])),y=as.vector(t(edge_df[, c(4,8)])),z=as.vector(t(edge_df[, c(5,9)])), lwd = 2, col = backbone_segment_color,line_antialias = TRUE)
        point_colors_df\$point_alpha = cell_alpha
        point_colors_df\$point_alpha[is.na(point_colors_df\$point_colors)] = 0
        points3d(point_colors_df[, c(\"data_dim_1\", \"data_dim_2\", \"data_dim_3\")], size = cell_size, col = point_colors_df\$point_colors, alpha = point_colors_df\$point_alpha,point_antialias = TRUE)
        legend3d(\"bottomright\",legend=unique(point_colors_df\$anno),pch = 16,col=unique(point_colors_df\$point_colors),cex=1, inset=c(0.02))
            movie3d(spin3d(axis=c(0,0,1),rpm=5),duration=30,movie=anno_name,dir=dir,convert=TRUE)
    ";
    }

    if (defined $opt{'T'} && defined $opt{'t'}){
        @tf_list = split(/,/,$opt{'t'});
        foreach (@tf_list) {$_ = "'$_'";}
        $opt{'t'}=join(', ', @tf_list );
        print R "

        #Generate TF Accessibility Colored Rotating Gif
        dev<-read.table(\"$opt{'T'}\",header=T)

        tf_list<-c($opt{'t'})
        out<- c()
        for (x in tf_list){
        check_each <- rownames(dev)[grepl(x,rownames(dev),ignore.case=TRUE)]
        out<-c(out,check_each)
        }

        goi<-subset(dev,rownames(dev) %in% out)
        goi<-data.frame(t(goi))
        goi\$cellID<-rownames(goi)
        point_colors_df<-merge(point_colors_df,goi,by.x=\"sample_name\",by.y=\"cellID\")

        colors<-colorRampPalette(brewer.pal(6,\"RdBu\"))(100)
        ";
        if (!defined $opt{'d'}){
            print R "
            for (i in grep(\"ENS\",colnames(point_colors_df))){
            title=colnames(point_colors_df)[i]
            point_colors_df\$colsplit<-as.numeric(cut(scale(point_colors_df[i]),100))
            point_colors_df\$point_colors<-colors[point_colors_df\$colsplit]

            open3d(windowRect = c(200, 200, 1024, 1024))
            segments3d(x=as.vector(t(edge_df[, c(3,7)])),y=as.vector(t(edge_df[, c(4,8)])),z=as.vector(t(edge_df[, c(5,9)])), lwd = 2, col = backbone_segment_color,line_antialias = TRUE)
            point_colors_df\$point_alpha = 0.8
            point_colors_df\$point_alpha[is.na(point_colors_df\$point_colors)] = 0
            points3d(point_colors_df[, c(\"data_dim_1\", \"data_dim_2\", \"data_dim_3\")], size = cell_size, col = point_colors_df\$point_colors, alpha = point_colors_df\$point_alpha, point_antialias = TRUE)
            legend3d(\"bottomright\",legend=colnames(point_colors_df)[i],inset=c(0.02),cex=0.8)
            movie3d(spin3d(axis=c(0,0,1),rpm=5),duration=30,movie=paste(colnames(point_colors_df)[i]),dir=dir,convert=TRUE)
        }
        ";
         } else {
        print R "

        for (i in grep(\"ENS\",colnames(point_colors_df))){
        title=colnames(point_colors_df)[i]
        bot_d_val=min(tail(sort(as.numeric(unlist(point_colors_df[i]))),round((nrow(point_colors_df[i])/100)*$opt{'d'})))
        levels(point_colors_df\$point_colors)<-c(\"darkgrey\",\"red\")
        point_colors_df[point_colors_df[i]>=bot_d_val,]\$point_colors<-\"red\"

        open3d(windowRect = c(200, 200, 1024, 1024))
        segments3d(x=as.vector(t(edge_df[, c(3,7)])),y=as.vector(t(edge_df[, c(4,8)])),z=as.vector(t(edge_df[, c(5,9)])), lwd = 2, col = backbone_segment_color,line_antialias = TRUE)
        point_colors_df\$point_alpha = 0.8
        point_colors_df\$point_alpha[is.na(point_colors_df\$point_colors)] = 0
        points3d(point_colors_df[, c(\"data_dim_1\", \"data_dim_2\", \"data_dim_3\")], size = cell_size, col = point_colors_df\$point_colors, alpha = point_colors_df\$point_alpha, point_antialias = TRUE)
        legend3d(\"bottomright\",legend=colnames(point_colors_df)[i],inset=c(0.02),cex=0.8)
        movie3d(spin3d(axis=c(0,0,1),rpm=5),duration=30,movie=paste(colnames(point_colors_df)[i]),dir=dir,convert=TRUE)

        }
        ";}

    } }else {

    if (!defined $opt{'o'}) {$opt{'o'} = $ARGV[0]; $opt{'o'} =~ s/\.dims$//; $opt{'o'}=basename($opt{'o'})};
    open R, ">$opt{'O'}/$opt{'o'}.3Dplot.r";
    print R "
        #Defined 3D Trajectory function
        library(rgl)
        library(RColorBrewer)

        dir=\"$opt{'O'}\"
        if(dir==\".\"){dir=getwd()}

        data_df<-read.table(file=\"$ARGV[0]\",header=F,row.names=1)
        cell_alpha=$opt{'a'}
        cell_size=$opt{'p'}
        backbone_segment_color=\"#000000\"

        colnames(data_df) <- c(\"data_dim_1\", \"data_dim_2\", \"data_dim_3\")
        data_df\$sample_name <- row.names(data_df)
        point_colors_df <- data.frame(sample_name = data_df\$sample_name,point_colors = \"darkgray\")

        point_colors_df\$point_alpha = cell_alpha
        point_colors_df\$point_alpha[is.na(point_colors_df\$point_colors)] = 0
        point_colors_df = merge(point_colors_df, data_df,by=\"sample_name\")
    ";

    if (defined $opt{'u'}) {
    print R "
        #Generate Uncolored Rotating Gif

        open3d(windowRect = c(0, 50, 800, 800))
        points3d(point_colors_df[, c(\"data_dim_1\", \"data_dim_2\", \"data_dim_3\")], size = cell_size, col = point_colors_df\$point_colors, alpha = point_colors_df\$point_alpha,point_antialias = TRUE)
        movie3d(spin3d(axis=c(0,0,1),rpm=5),duration=30,movie=\"Unannotated\",dir=dir,convert=TRUE)

    ";
    }
    if (defined $opt{'A'}){
    print R "
            # Generate Annotation Plot
            annot<-read.table(file=\"$opt{'A'}\",header=F)
            anno_name_list<-unlist(strsplit(\"$opt{'A'}\",\"/\"))
            anno_name_list2<-unlist(strsplit(anno_name_list[length(anno_name_list)],\"[.]\"))
            anno_name<-anno_name_list2[length(anno_name_list2)-1]
            colnames(annot)<-c(\"sample_name\",\"anno\")
            point_colors_df<-merge(point_colors_df,annot,by.x=\"sample_name\")
        colors<-colorRampPalette(brewer.pal(length(unique(point_colors_df\$anno)),\"Set2\"))(length(unique(point_colors_df\$anno)))
        point_colors_df<-transform(point_colors_df,colsplit=as.numeric(factor(point_colors_df\$anno)))
        point_colors_df\$point_colors<-colors[point_colors_df\$colsplit]

        open3d(windowRect = c(0, 50, 800, 800))
        point_colors_df\$point_alpha = cell_alpha
        point_colors_df\$point_alpha[is.na(point_colors_df\$point_colors)] = 0
        points3d(point_colors_df[, c(\"data_dim_1\", \"data_dim_2\", \"data_dim_3\")], size = cell_size, col = point_colors_df\$point_colors, alpha = point_colors_df\$point_alpha,point_antialias = TRUE)
        legend3d(\"bottomright\",legend=unique(point_colors_df\$anno),pch = 16,col=unique(point_colors_df\$point_colors),cex=1, inset=c(0.02))
            movie3d(spin3d(axis=c(0,0,1),rpm=5),duration=30,movie=anno_name,dir=dir,convert=TRUE)
    ";
    }

    if (defined $opt{'T'} && defined $opt{'t'}){
    @tf_list = split(/,/,$opt{'t'});
    foreach (@tf_list) {$_ = "'$_'";}
    $opt{'t'}=join(', ', @tf_list );
    print R "

        #Generate TF Accessibility Colored Rotating Gif
        dev<-read.table(\"$opt{'T'}\",header=T)

        tf_list<-c($opt{'t'})
        out<- c()
        for (x in tf_list){
        check_each <- rownames(dev)[grepl(x,rownames(dev),ignore.case=TRUE)]
        out<-c(out,check_each)
        }

        goi<-subset(dev,rownames(dev) %in% out)
        goi<-data.frame(t(goi))
        goi\$cellID<-rownames(goi)
        point_colors_df<-merge(point_colors_df,goi,by.x=\"sample_name\",by.y=\"cellID\")

        colors<-colorRampPalette(brewer.pal(6,\"RdBu\"))(100)

        for (i in grep(\"ENS\",colnames(point_colors_df))){
        title=colnames(point_colors_df)[i]
        point_colors_df\$colsplit<-as.numeric(cut(scale(point_colors_df[i]),100))
        point_colors_df\$point_colors<-colors[point_colors_df\$colsplit]

        open3d(windowRect = c(200, 200, 1024, 1024))
        point_colors_df\$point_alpha = 0.8
        point_colors_df\$point_alpha[is.na(point_colors_df\$point_colors)] = 0
        points3d(point_colors_df[, c(\"data_dim_1\", \"data_dim_2\", \"data_dim_3\")], size = cell_size, col = point_colors_df\$point_colors, alpha = point_colors_df\$point_alpha, point_antialias = TRUE)
        legend3d(\"bottomright\",legend=colnames(point_colors_df)[i],inset=c(0.02),cex=0.8)
        movie3d(spin3d(axis=c(0,0,1),rpm=5),duration=30,movie=paste(colnames(point_colors_df)[i]),dir=dir,convert=TRUE)
    }

    ";
    }





close R;
system("$Rscript $opt{'O'}/$opt{'o'}.3Dplot.r");
if (!defined $opt{'X'}) {
        system("rm -f $opt{'O'}/$opt{'o'}.3Dplot.r");
}
}
}
1;
