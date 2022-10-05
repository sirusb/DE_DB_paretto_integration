#' Install any missing package
#'
#' @return
#' @export
#'
#' @examples
freeze <- function() {
    package_mat <- installed.packages()[,c("Package","Version","LibPath")]
    requirement_lst <- apply(package_mat, 1, function(pp) {
    paste(paste0(pp["Package"],"==",pp["Version"]),pp["LibPath"],sep="\t")
    })

    rootDir = dirname(funr::sys.script())
    requirements <- unique(unname(sort(unlist(requirement_lst))))
    write.table(requirements,glue::glue("{rootDir}/R.freeze.txt"),quote=F,row.names=F,col.names=F)
}


installMissingPkg <- function(...){

  requiredPackages<-unlist(list(...))

  if(!all(requiredPackages %in% .packages(all.available = TRUE))){
    missingPackages = setdiff(requiredPackages, .packages(all.available = TRUE))
    pak::pkg_install(missingPackages)
    missingPackages = setdiff(requiredPackages, .packages(all.available = TRUE))
    BiocManager::install(missingPackages,update=FALSE)
    missingPackages = setdiff(requiredPackages, .packages(all.available = TRUE))
    for(pkg in missingPackages){
      if( !is.element(pkg, .packages(all.available = TRUE)) ) {
        message(paste0(pkg," missing, installing it ..."))
        BiocManager::install(pkg,update=FALSE)
      }
    }    
  }

}

using<-function(..., character.only = FALSE) {

    if(character.only){
      libs <- unlist(list(...))
    }else{
      libs <- unlist(list(as.character(substitute(...) ) ))
    }
  
    req<-unlist( lapply(libs,function(p) suppressPackageStartupMessages(require(p,character.only=TRUE)) ) )
    need<-libs[req==FALSE]
    if(length(need)>0){ 
        installMissingPkg(need)
        lapply(need,require,character.only=TRUE)
    }
}


checkPackages <- function(){
  if(!suppressPackageStartupMessages(require("pak",character.only = TRUE))){
    install.packages("pak", repos = "https://r-lib.github.io/p/pak/dev/")
  }

  requiredPackages <- c("cli","logger","openxlsx","fs","funr","readxl","readr",
                        "glue","tidyverse", "argparse",
                        "ggplot2","ggsci","formattable",
                        'logger','rPref','ggsci',"plotly",         
                        'GenomicRanges',"GenomicFeatures"                        
  )
  
  suppressPackageStartupMessages(using(requiredPackages, character.only=TRUE))
  #freeze()  
}

prepareFolders <- function(outDir){

  logger::log_info("Preparing the output directory")

  fs::dir_create(outDir)
  fs::dir_create(glue("{outDir}/Objects"))
  fs::dir_create(glue("{outDir}/Figures"))
  fs::dir_create(glue("{outDir}/Summary"))

}

parseParameters <- function(){

  suppressPackageStartupMessages(using(logger))

  Supported_versions = c("hg19","hg38")

  # Parse paramters
  parser <- argparse::ArgumentParser(
    description = "Main script for the Differential binding and Differential expression Pareto integrative analysis."
  );


  parser$add_argument(
    "-b",
    "--db",
    dest='dbResults',
    help = "Path to the differential binding results. Should contain the following two columns: 'Region' (coordinates of the peak in the format chr:start-end), 'log2FC', 'FDR' ",
    type = "character"
  );


  parser$add_argument(
    "-e",
    "--de",
    dest='deResults',
    help = "Path to the differential gene expression tsv file. Should containt the following columns: 'gene', 'log2FC', 'FDR'. Make sure to clean the lowly expressed genes first.",
    type= "character"
  );
  

  parser$add_argument(
    "-g",
    "--genome",
    dest='genomeVersion',
    help = paste0("Version of the reference genome. Accepts the same values as the one of the following:",paste(Supported_versions,collapse = ', ')), # ugly code just in case glue is not installed.
    type = "character"
  );

  parser$add_argument(
    "-t",
    "--top",
    dest='ntop',
    help = "How many top correlated peak-gene relationship you want to find. Defaut: 1000",
    type = "integer",
    default=1000
  );

  parser$add_argument(
    "-f",
    "--rnaFC",
    dest='rnaFC',
    help = "RNA-seq fold-change threshold to highlight the most significant DE genes with. Default:2",
    type = "double",
    default=2.0
  );

  parser$add_argument(
    "-r",
    "--regionFC",
    dest='regionFC',
    help = "Region fold-change threshold to highlight the most significant DB regions with. Default:2.",
    type = "double",
    default=2.0
  );

  parser$add_argument(
    "-q",
    "--qval",
    dest='qvalThr',
    help = "qvalue (FDR) threshold to apply to select the signicant genes or regions. Default: 1, no FDR filtering applied.",
    type = "double",
    default=1.0
  );
  

  parser$add_argument(
    "-o",
    "--out",
    dest='outDir',
    help = "Path to the output directory",
    type = "character",
    default=getwd()
  );

 parser$add_argument(
    "-p",
    "--prefix",
    dest='prefix',
    help = "Analysis prefix",
    type = "character",
    default="pareto"
  );

  # Parse arguments
  opt <- parser$parse_args();

  for(n in names(opt)){
    if(is.null(opt[[n]])){
      cli::cli_blockquote("RNA-seq and Differential Binding integrative analysis")
      parser$print_help()
      cli::cli_alert_danger("Couldn't parse correctly the provided arguments.")
      stop()
    }
  }

  ## check that the files exist 
  if(!file.exists(opt$dbResults)){
    cli::cli_alert_danger(glue::glue("Couldn't find {opt$dbResults}."))
    stop()
  }

  if(!file.exists(opt$deResults)){
    cli::cli_alert_danger(glue::glue("Couldn't find {opt$deResults}."))
    stop()
  }
  
  ## check if the genome version is supported 
  if(!opt$genomeVersion %in% Supported_versions ){
    cli::cli_alert_danger(glue("Only the following genomes are supported {glue_collapse(Supported_versions, sep=', ')}"))
    stop()
  }
  

  prepareFolders(opt$outDir)

  config = list(
    DBPath = opt$dbResults,
    DEResults  = opt$deResults,    
    Genome = opt$genomeVersion,
    ntop = opt$ntop,
    rnaFC = opt$rnaFC,
    regionFC = opt$regionFC, 
    FDR = opt$qvalThr,
    assay = "ChIP-seq",
    prefix = opt$prefix,
    outDir = opt$outDir
  )

  readr::write_rds(config, file=glue("{config$outDir}/Objects/{config$prefix}_config.rds"))
  print(config)
  return(config)
}

#' Identify and fix duplicated genes by combining their geneSymbol with geneID
#'
#' @param mat data.frame with columns eneSymbol" and "geneID" present
#' @param returnMat Logical Wheter to return the whole matrix with the row.names containing the fixed gene names, 
#' or just trturn the fixed gene names as avector.
#'
#' @return
#' A
#' @export
#'
#' @examples
fixDuplicatedGenes <- function(mat, returnMat=TRUE){

  if(!all(c("geneSymbol","geneID") %in% colnames(mat))){
    mat = as.data.frame(mat)
    rownames(mat) = mat$gene
    return(mat)
  }

  mat = as.data.frame(mat)
  genes = mat$geneSymbol
  ind_nonUniqGenes <- duplicated(genes)

  genes[ind_nonUniqGenes] <- paste(genes[duplicated(genes)],"_",mat$geneID[ind_nonUniqGenes], sep="")

  if(returnMat){
    rownames(mat) <- genes
    return(mat)
  }else{
    return(genes)
  }
}


plotCorrelatedInteractions <- function(config, DEG_pareto_res, correlationType= "correlated"){

  if(! correlationType %in% c("correlated","anticorrelated","all")){
    cli::cli_abort('Please make sure that the correlationType is one of the following values: "correlated","anticorrelated","all" ')    
  }

  dout = glue("{config$outDir}/Figures")

  corr_plots = DEG_pareto_res %>%
    imap(.f = function(lst,cmp){
      df = lst[['Combined_zScores']]
      colnames(df)[2] = "Combined_zscore"

      df$Text = glue::glue("\nGene: {df$gene}\nRegion: {df$Region}\n Distance To TSS: {prettyNum(df$distanceToTSS,big.mark = ',')}bp")
      
      col = grep("DB_fc", colnames(df),value = TRUE)
      col = grep('zscore',col,invert = TRUE,value = TRUE)

      

      if(correlationType == "correlated"){
        upRelation = subset(df, Relationship == "Both Up")
        downRelation = subset(df, Relationship == "Both Down")  
      }

      if(correlationType == "anticorrelated"){
        upRelation = subset(df, Relationship == "RNA_Up CR_Down")
        downRelation = subset(df, Relationship == "RNA_Down CR_Up")
      }

      if(correlationType == "all"){

        upRelation = subset(df, Relationship == "Both Up")
        downRelation = subset(df, Relationship == "Both Down")  

        upAntiRelation = subset(df, Relationship == "RNA_Up CR_Down" )
        downAntiRelation = subset(df, Relationship == "RNA_Down CR_Up")
      }
      
      
      top_up = upRelation %>%
        arrange(-abs(Combined_zscore)) %>%
        head(n=10)
      
      top_down = downRelation %>%
        arrange(-abs(Combined_zscore)) %>%
        head(n=10)


      if(correlationType == "all"){
        top_anti_up = upAntiRelation %>%
            arrange(-abs(Combined_zscore)) %>%
            head(n=10)
      
        top_anti_down = downAntiRelation %>%
          arrange(-abs(Combined_zscore)) %>%
          head(n=10)      
      }

      
      ylab = gsub("DB_fc_","",col)
      ylab = gsub("_vs_"," vs ", ylab)
      ylab = paste0(config$assay," log2FC\n ( ",ylab," )")
      
      
      geom_x_lbl = "Correlated with gene up-regulation"
      geom_y_lbl = "Correlated with gene down-regulation"

      if(correlationType == "anticorrelated"){

        geom_x_lbl = "Anti-correlated with gene up-regulation"
        geom_y_lbl = "Anti-correlated with gene down-regulation"

      }

      if(correlationType == "all"){
        geom_x_lbl = "Correlated with gene up-regulation"
        geom_y_lbl = "Correlated with gene down-regulation"      

        geom_x_lbl2 = "Anti-correlated with gene up-regulation"
        geom_y_lbl2 = "Anti-correlated with gene down-regulation"
      }


      p = ggplot(df, aes_string(x = "DE_FC", y= col)) +
            geom_point(size=1.2, color="grey80",mapping = aes(text=Text)) +
            geom_point(data = upRelation, size=1.2, ,mapping = aes(text=Text, color=geom_x_lbl)) +
            geom_point(data = downRelation,size=1.2,mapping = aes(text=Text, color=geom_y_lbl)) 

      if(correlationType == "all"){
        p = p +
            geom_point(data = upAntiRelation, size=1.2, ,mapping = aes(text=Text, color=geom_x_lbl2)) +
            geom_point(data = downAntiRelation,size=1.2,mapping = aes(text=Text, color=geom_y_lbl2)) 
      }


    breaks = c(geom_x_lbl, geom_y_lbl)
    cols = c("#c63254","#005073")
    if(correlationType == "all"){
        breaks = c(geom_x_lbl, geom_y_lbl, geom_x_lbl2, geom_y_lbl2)    
        cols = c("#c63254","#005073","orange","purple")
    }

    p = p + geom_hline(yintercept = c(-1,1) * log2(config$rnaFC),linetype='dashed') +
            geom_vline(xintercept = c(-1,1) * log2(config$regionFC),linetype='dashed') +
            scale_color_manual(name = "Relationship", 
                                      breaks = breaks, 
                                      values = cols,
                                      labels = breaks) +
            theme_bw() 

    if(correlationType != "all"){
      p =  p + ggtitle("Peak-gene pairs showing associated gene-peak interactions",
                      subtitle = glue::glue("(Total={nrow(df)} pairs, {geom_x_lbl}={nrow(upRelation)}, {geom_y_lbl}={nrow(downRelation)})"))
    }else{
      p =  p + ggtitle("Peak-gene pairs showing associated gene-peak interactions",
                      subtitle = glue::glue("(Total={nrow(df)} pairs,\n{geom_x_lbl}={nrow(upRelation)}, {geom_y_lbl}={nrow(downRelation)})\n{geom_x_lbl2}={nrow(upAntiRelation)}, {geom_y_lbl2}={nrow(downAntiRelation)})"))
    }

            
      p = p + theme(plot.title = element_text(hjust = 0.5,size = 12, color="black"), 
                    plot.subtitle = element_text(hjust = 0.5,size = 10, color="black"),
                    axis.text = element_text(color="black")) +
              xlab(glue("RNA-seq log2FC\n( {gsub('_vs_|__VS__',' vs ',config$prefix)} )")) +
              ylab(ylab)
            
      p = p + geom_point(data = top_up, color="black",mapping = aes(test=Text))
      p = p + geom_point(data = top_down, color="black",mapping = aes(test=Text))    
      
      p = p + ggrepel::geom_text_repel(data = top_up, aes_string(x = "DE_FC", y= col, label="gene"), color="black",size=4)
      p = p + ggrepel::geom_text_repel(data = top_down, aes_string(x = "DE_FC", y= col, label="gene"), color="black", size=4)

      if(correlationType == "all"){

        p = p + geom_point(data = top_anti_up, color="black",mapping = aes(test=Text))
        p = p + geom_point(data = top_anti_down, color="black",mapping = aes(test=Text))      

        p = p + ggrepel::geom_text_repel(data = top_anti_up, aes_string(x = "DE_FC", y= col, label="gene"), color="black",size=4)
        p = p + ggrepel::geom_text_repel(data = top_anti_down, aes_string(x = "DE_FC", y= col, label="gene"), color="black", size=4)      
      }


      fout = glue::glue("{dout}/{cmp}_{correlationType}_peaks.pdf")
      ggsave(filename = fout,width = 12,height = 8,plot = p)
      
      fout = glue::glue("{dout}/{cmp}_{correlationType}_peaks.png")
      ggsave(filename = fout,width = 12,height = 8,plot = p)
          
      ## Create interactive plot
      plt = ggplotly(p)
      fout = glue::glue("{dout}/{cmp}_{correlationType}_peaks.html")
      htmlwidgets::saveWidget(plt, fout)    
        
      return(p)
    })

  logger::log_info("plotCorrelatedInteractions - Saving {correlationType} plots")
  fout=glue("{config$outDir}/Objects/{config$prefix}_{correlationType}_plots.rds")
  readr::write_rds(corr_plots, fout)

  invisible(corr_plots)
}



annotateParetoResults <- function(config, DEG_pareto_res){

  logger::log_info("annotateParetoResults - Filtering correlated interactions ... ")

  DEG_pareto_res_filtered_correlated = DEG_pareto_res %>%
  imap(.f = function(lst, cmp){        
    allGenesRes = lst[['Combined_zScores']] %>% rownames_to_column(var = "Interaction") 
            
    colnames(allGenesRes)[3] ='Combined_zscore'

    #allGenesRes2[['Combined_zscore']] = allGenesRes[,2]                
    res = subset(allGenesRes, Combined_zscore>0) # Same FC direction

    res2 = subset(res, abs(DE_FC)>= log2( config$rnaFC ) & Gene_DE_FDR  < config$FDR ) 

    pos_db_fc = grep("DB_fc_", colnames(res2))
    colnames(res2)[pos_db_fc] = "DB_FC"

    res2 = subset(res2, abs(DB_FC) >= log2( config$regionFC) &  Region_FDR < config$FDR)

    if(nrow(res2) == 0){
        return(res2)
    }

    max_cor  = split(res, res$gene) %>%
      imap(.f= function(df,gene){        
        reg_score= exp(-( (4 * df$distanceToTSS/100e3) ) )
        df$RegScore =df$Combined_zscore * reg_score
        
        df = subset(df, RegScore >= median(df$RegScore))
        return(df)
      }) %>%
      do.call(what = "rbind")
    
    
    res2_bothUp = subset(res2, DE_FC > 0)
    res2_bothDown = subset(res2, DE_FC < 0) 
    
    max_cor$Relationship = 'NONE'
  
    max_cor$Relationship[max_cor$Interaction %in% res2_bothUp$Interaction ] = "Both Up"
    max_cor$Relationship[max_cor$Interaction %in% res2_bothDown$Interaction ] = "Both Down"
    
    return(max_cor)
  })


  logger::log_info("annotateParetoResults - Filtering anti-correlated interactions ... ")

  DEG_pareto_res_filtered_anticorrelated = DEG_pareto_res %>%
    imap(.f = function(lst, cmp){
      allGenesRes = lst[['Combined_zScores']] %>% rownames_to_column(var = "Interaction") 

      colnames(allGenesRes)[3] ='Combined_zscore'                      
      
      res = subset(allGenesRes, Combined_zscore < 0) # different FC direction
      res2 = subset(res,  abs(DE_FC)>= log2( config$rnaFC ) & Gene_DE_FDR  < config$FDR ) 

      pos_db_fc = grep("DB_fc_", colnames(res2))
      colnames(res2)[pos_db_fc] = "DB_FC"
      
      res2 = subset(res2, abs(DB_FC) >= log2( config$regionFC) &  Region_FDR < config$FDR)

      if(nrow(res2) == 0){
        return(res2)
      }
      
      max_cor  = split(res, res$gene) %>%
        imap(.f= function(df,gene){          
          reg_score= exp(-( (4 * df$distanceToTSS/100e3) ) )
          df$RegScore = -df$Combined_zscore * reg_score
          df = subset(df, RegScore >= median(df$RegScore))
          return(df)
        }) %>%
        do.call(what = "rbind")
      
      
      res2_UpDown = subset(res2, DE_FC > 0)
      res2_DownUp = subset(res2, DE_FC < 0) 
      
      max_cor$Relationship = 'NONE'
    
      max_cor$Relationship[max_cor$Interaction %in% res2_UpDown$Interaction] = "RNA_Up CR_Down"
      max_cor$Relationship[max_cor$Interaction %in% res2_DownUp$Interaction] = "RNA_Down CR_Up"      
      return(max_cor)
    })

  ## Aggregate results 
  for(cmp in names(DEG_pareto_res)){

    DEG_pareto_res[[cmp]][["Combined_zScores"]]$Relationship = "NONE"
    if(nrow(DEG_pareto_res_filtered_anticorrelated[[cmp]]) > 0){    
      DEG_pareto_res[[cmp]][["Combined_zScores"]][ DEG_pareto_res_filtered_anticorrelated[[cmp]]$Interaction, ]$Relationship = DEG_pareto_res_filtered_anticorrelated[[cmp]]$Relationship
    }

    if(nrow(DEG_pareto_res_filtered_correlated[[cmp]]) > 0 ){
      DEG_pareto_res[[cmp]][["Combined_zScores"]][ DEG_pareto_res_filtered_correlated[[cmp]]$Interaction, ]$Relationship = DEG_pareto_res_filtered_correlated[[cmp]]$Relationship
    }    

  }
  

  fout=glue("{config$outDir}/Objects/{config$prefix}_correlated_anticorrelated_res.rds")
  logger::log_info("annotateParetoResults - Saving results to {fout}")

  readr::write_rds(DEG_pareto_res, fout)
  logger::log_info("{symbol$tick} annotateParetoResults - DONE!")
  return(DEG_pareto_res)
}


#' Calculate the Pareto correlation.
#'
#' @param DE_res The differential expression table.
#' @param chip_fc A matrix that contains the log2FC of the peaks
#' @param gns GRanges object, contains the coordinates of the genes in the colulm 'fixedName'
#' @param custom_cond character vector, contains the constrains to apply, By default the peaks that
#' are correlated and anti-correlated with all signals are found.
#' @param nlevels integer, the number of the top levels to return (similat to top PCs).
#' @return
#' A list that contains the top peak-gene intractions in 'correlated', 'anti-correlate' and 'custom' settings.
#' The combined z-scores for all the peak-gene pairs are returned too.
#' @export
#'
#' @examples
paretoCorr <- function(DE_res, chip_fc, DB_df, gns, custom_cond=NULL, nlevels=10){

  DE_res$zscore = DE_res$log2FC/stats::sd(DE_res$log2FC)

  chip_zscore = apply(chip_fc,2,function(fc) fc/stats::sd(fc))


  pos <- match( rownames(DE_res), gns$fixed_names)

  smp_genes = gns[pos,]
  start(smp_genes) = start(smp_genes) - 5e4
  end(smp_genes) = end(smp_genes) + 5e4

  repro_peaks.GR = GRanges(rownames(chip_fc))

  smp_ovp = findOverlaps(smp_genes, repro_peaks.GR, ignore.strand=TRUE) %>% as.data.frame()

  smp_ovp_lst = split(smp_ovp, smp_ovp$queryHits)

  logger::log_info("Calculating peak-genes combined scores")

  pb = progress_estimated(length(smp_ovp_lst))

  peaksCorrelation = smp_ovp_lst %>%
    imap(.f = function(df,id){
      peaks = repro_peaks.GR[df$subjectHits]

      expr_fc = DE_res$zscore[as.numeric(id)]
      regions = Signac::GRangesToString(peaks,sep = c(":","-"))

      gene.gr = gns[gns$fixed_names == rownames(DE_res)[as.numeric(id)] ]

      d = distance(peaks,promoters(gene.gr,1,1))

      if(length(df$subjectHits)==1){
        chip_z = matrix(chip_zscore[df$subjectHits, ], nrow=1)
        colnames(chip_z) = colnames(chip_zscore)
        rownames(chip_z) = rownames(chip_zscore)[df$subjectHits]
      }else{
        chip_z =  chip_zscore[df$subjectHits, ]
      }

      cors <- (expr_fc * chip_z) %>%  data.frame() %>%  rownames_to_column(var = "Region") %>%  mutate(gene =  rownames(DE_res)[as.numeric(id)])
      #cors = cors[,-2]
      #colnames(cors)[2] = gsub("_1$","",colnames(cors)[2])

      if(length(df$subjectHits)==1){
        chip_fc2 = matrix(chip_fc[df$subjectHits, ], nrow=1)
        colnames(chip_fc2) = colnames(chip_fc)
        rownames(chip_fc2) = rownames(chip_fc)[df$subjectHits]
      }else{
        chip_fc2 = chip_fc[df$subjectHits, ]
      }


      colnames(chip_fc2) = glue::glue("DB_fc_{colnames(chip_fc2)}")
      colnames(chip_z) = glue::glue("DB.zscore.{colnames(chip_z)}")

      DB_FDRs <- DB_df %>% dplyr::select(starts_with('FDR'))
      colnames(DB_FDRs) = glue("Region_{colnames(DB_FDRs)}")
      ## add more info to the results
      rpos = match(cors$Region, DB_df$Region)
      
      cors$DE_FC = DE_res$log2FC[as.numeric(id)]
      cors$Gene_DE_FDR= DE_res$FDR[as.numeric(id)]
      cors$DE_FC.zscore = expr_fc
      cors = cbind(cors, DB_df[pos,])
      #cors$Region_FDR = DB_df$FDR[rpos]  
      #cors$Region_Regulation =  DB_df$Regulation[rpos]    
      #cors$FeatureAssignment =  DB_df$FeatureAssignment[rpos]    
      cors = cbind(cors, chip_fc2)
      cors = cbind(cors, chip_z)
      cors$distanceToTSS = d

      pb$tick()$print()
      return( cors)
    })

  peaksCorrelation =  peaksCorrelation %>%
    do.call(what = "rbind")
    
  ###
  conditions <- list(correlated = rep("high",ncol(chip_fc)),
                     anticorrelated = rep("low",ncol(chip_fc))
  )

  if(!is.null(custom_cond)){
    if(length(custom_cond) != ncol(chip_fc)){
      cli::cli_abort("Please provide a condition list that have the same size are the signals.")
    }
    conditions[['custom']] = custom_cond
  }


  res = peaksCorrelation

  ret <- list()

  for(cnd in names(conditions)){
    logger::log_info(glue::glue("Identifying {cnd} gene-peaks ..."))
    ## Find correlated regions
    p <- rPref::empty()
    for(i in 1:length(conditions[[cnd]])){
      mrk = colnames(res)[2:(ncol(chip_fc)+1)][i]
      if(conditions[[cnd]][i]=="high"){
        p = p * rPref::high_(mrk)
      }else{
        p = p * rPref::low_(mrk)
      }
    }

    rownames(res) = glue::glue("{res$Region}__{res$gene}")
    df_final = res[,2:(ncol(chip_fc)+1)]
    psel_res <- rPref::psel(df = df_final, pref = p, top = nlevels)
    names(psel_res)[names(psel_res) == '.level'] <- 'front'
    psel_res = cbind(psel_res, res[rownames(psel_res),])

    tokeep = !duplicated(colnames(psel_res))
    psel_res = psel_res[,tokeep]
    colnames(psel_res) = gsub("DB_fc_zscore.","DB_fc.",colnames(psel_res))

    ## if we require only one dataset to be high or low,
    ## we only keep the positisitive or negative z-scroes
    if(length(unique(conditions[[cnd]])) >1){
      tbl = table(conditions[[cnd]])
      touse = names(tbl)[tbl==1]
      p.touse = which(conditions[[cnd]] == touse)

      # if(length(touse) ==1 & touse == "low"){
      #   pos = which(psel_res[[p.touse]] < 0)
      # }else{
      #   pos = which(psel_res[[p.touse]] > 0)
      # }
      # psel_res = psel_res[pos,]
    }

    ret[[cnd]] = psel_res
  }

  tokeep = grep("_2$", colnames(res), invert = TRUE)

  res = res[,tokeep]  
  colnames(res) = gsub("_1$","", colnames(res))
  colnames(res)[2] = glue("Combined.zscore.{colnames(res)[2]}")
  ret[['Combined_zScores']] = res
  return(ret)
}


getGenesCoordinated <- function(config){

  info_file = switch(config$Genome,
				"hg19" = glue("{config$sourceDir}/geneInfo/gencode.v19.annotation.gtf.gene"),
				"hg38" = glue("{config$sourceDir}/geneInfo/gencode.v31.primary_assembly.annotation.gtf.gene")
				)
  logger::log_info("Loading gene info from : {info_file}")
  genome_genes = readr::read_tsv(info_file, show_col_types = FALSE)

  colnames(genome_genes) = gsub("Ensembl_ID","geneID", colnames(genome_genes))
  colnames(genome_genes) = gsub("Gene_Symbol","geneSymbol", colnames(genome_genes))

  colnames(genome_genes) = gsub("Chr","seqnames", colnames(genome_genes))
  colnames(genome_genes) = gsub("Start","start", colnames(genome_genes))
  colnames(genome_genes) = gsub("End","End", colnames(genome_genes))

  genome_genes$fixed_names = fixDuplicatedGenes(genome_genes, returnMat = FALSE)
  genome_gene.gr = GRanges(genome_genes)
  
  ## We remove level-3 genes for the moment
  genome_gene.gr = subset( genome_gene.gr, Annotation_Level != 3)
  return(genome_gene.gr )
}

doParetoAnalysis <- function(config){

  logger::log_info(glue("Reading DEG data at : {config$DEResults}"))
  DEG_res = readr::read_tsv(config$DEResults, show_col_types = FALSE)   
    

  logger::log_info(glue("Reading differential binding data at: {config$DBPath}"))
  cmpName = strsplit(x =  basename(config$DBPath),split = '\\.') %>% unlist()
  cmpName = cmpName[1]
  DB_dfs = list() 
  DB_dfs[[cmpName]] = readr::read_tsv(config$DBPath, show_col_types = FALSE)   

  reproPeaks = GRanges(DB_dfs[[cmpName]]$Region)

  logger::log_info("Calculating z-scores")
  
  #Peaks_FC <- matrix(0, nrow=nrow(DB_dfs[[cmpName]]), ncol = length(DB_dfs))
  Peaks_FC = DB_dfs[[cmpName]] %>% dplyr::select(starts_with("Log2FC")) %>% data.matrix()
  colnames(Peaks_FC) = gsub("Log2FC.","",colnames(Peaks_FC),  ignore.case = TRUE)
  rownames(Peaks_FC) = DB_dfs[[cmpName]]$Region

  #for(cmp in names(DB_dfs)){ 
  #  Peaks_FC[,cmp] = DB_dfs[[cmp]]$log2FC
  #}

  colnames(Peaks_FC) = gsub("__VS__","_vs_",colnames(Peaks_FC))


  DEG_pareto_res <- list()

  #for(cmp in colnames(Peaks_FC)){
    
  cli::cli_h1(glue::glue("Processing {cmpName}"))
    
  DEG = DEG_res
  
  DEG = fixDuplicatedGenes(DEG, returnMat = TRUE)

  colnames(DEG) = gsub("log2FC","logFC", colnames(DEG),ignore.case = TRUE)   
  colnames(DEG) = gsub("geneSymbol","gene", colnames(DEG),ignore.case = TRUE) 
  
  
  #fc_mat = matrix(0, ncol=2, nrow=nrow(Peaks_FC))
  #rownames(fc_mat) = rownames(Peaks_FC)

  #fc_mat[,1] = Peaks_FC[,cmp]
  #fc_mat[,2] = Peaks_FC[,cmp]
  #colnames(fc_mat) = glue("{cmp}_{1:2}")

  
  custom_cond = ifelse(grepl(cmp, colnames(Peaks_FC)),"high",'low')

  logger::log_info("Loading {config$Genome} gene coordinates ...")
  genes = getGenesCoordinated(config)
  logger::log_info("{symbol$tick} [DONE] Gene coordinates loaded")

  pos = which(rownames(DEG) %in% genes$fixed_names)
  
  if(length(pos) == 0){      
    cli::cli_abort("ERROR [{format(Sys.time(), '%Y-%m-%d %H:%M:%S')}] doParetoAnalysis - Make sure the genome version you are using is consistant, couldn't match names.")
  }

  DEG = DEG[pos,]
  
  logger::log_info("Running integrative analysis for {cmpName}")
  parRes = paretoCorr(DE_res = DEG, 
            chip_fc= Peaks_FC,
            DB_df = DB_dfs[[cmpName]],
            gns= genes,
            custom_cond = custom_cond,
            nlevels = config$ntop
  )

  DEG_pareto_res[[cmpName]] = parRes
  logger::log_info("{symbol$tick} {cmpName} Done!")
  #}

  ## saving results
  fout = glue("{config$outDir}/Objects/{config$prefix}_DEG_pareto_res.rds")

  logger::log_info("doParetoAnalysis - Saving initial resuts to : {fout}")  
  readr::write_rds(DEG_pareto_res, fout)


  ## For each gene select the most correlated peaks


  logger::log_info("doParetoAnalysis - DONE.")  
  return(DEG_pareto_res)
}


generatePlots <- function(config, paretoRes){

  for(correlationType in c("correlated","anticorrelated","all")){
    logger::log_info("generatePlots - Generating {correlationType} plots ...")
    plotCorrelatedInteractions(config = config, 
                               DEG_pareto_res =  paretoRes, 
                               correlationType = correlationType
                              )
  }

  logger::log_info("generatePlots - DONE.")
}








SaveResults <- function(config, results){

  cmp = names(results)[1]  
  fout = glue("{config$outDir}/Summary/{cmp}_correlation_results.xlsx")

  logger::log_info("Saving results to {fout}")

  res = results[[1]][['Combined_zScores']] %>% rownames_to_column(var = "Interaction") 
  excel_cols = c("#C5D9F1","#EBF1DE","#E4DFEC","#DAEEF3","#FDE9D9","#F2F2F2","#ACB9CA")
  
  wb <- createWorkbook()
  addWorksheet(wb = wb,sheetName = "DEG_DB_association")

  cols1 = c("Interaction","Region","FeatureAssignment","Region_Regulation","gene","distanceToTSS","Relationship")
  othercols = setdiff(colnames(res), cols1)
  res = res[,c(cols1, othercols)]

  writeDataTable(wb,
                  sheet = 1,
                  x = res,
                  tableStyle ="TableStyleLight15", 
                  bandedCols = F,
                  bandedRows = F)

  col_groups = list(1:6,7,8, 9:11, 12:14)
  
  for(i in 1:length(col_groups)){
    
    for(j in col_groups[[i]]){
      hs1 <- createStyle(fgFill = excel_cols[i],fontColour = "black",textRotation = 60,textDecoration = "Bold")
      addStyle(wb, 1, style = hs1, rows = 1, cols=j)
    }
  }      
  saveWorkbook(wb, fout, overwrite = TRUE)

  logger::log_info("{symbol$tick} SaveResults - DONE.")
}


