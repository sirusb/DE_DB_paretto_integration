#!/bin/env Rscript
#########################################################################
# Copyright (c) 2021 ~ Nadhir Djekidel && St Jude
#
# This source code is released for free distribution under the terms of the
# Attribution-NonCommercial-ShareAlike 4.0 International (cc)
#
#*Author:       Nadhir Djekidel < djek.nad [at] gmail DOT com >
# File Name: compMappingReport.R
# Description:
#########################################################################

prepareEnv <- function(){
    
    x <- "BiocManager"
    if (!requireNamespace(x, quietly = TRUE)) { install.packages(x) }

    ## Check if funr is installed (it is essential to get the source file path).
    if(!'funr' %in% .packages(all.available = TRUE)){
    message("The essential package 'funr' is missing, installing it ...")
    BiocManager::install('funr')
    }

    suppressPackageStartupMessages(require(funr))
    suppressPackageStartupMessages(require(logger))

    rootDir = dirname(funr::sys.script())
    sourceDir = paste0(rootDir,"/R")
    tmp =lapply( list.files(sourceDir,pattern = ".R",full.names = TRUE) , source)


    ## check if all the packages are installed
    checkPackages()

    ## Parse parameters    
    log_info("Checking and loading packages ...")
    config  = parseParameters() 
    config$sourceDir = rootDir
    return(config)   
}

main <- function(config){

    ## Run analysis
    log_info("Running Pareto analysis")
    DE_pareto_res <- doParetoAnalysis(config)

    ## Annotate the results based on the selected thresholds
    corredRes = annotateParetoResults(config, DE_pareto_res)

    ## Plot correlations
    generatePlots(config, corredRes)    

    ## Save results
    SaveResults(config, corredRes)
    
    logger::log_info("Integrative analysis done. Thanks.")
    cli::cli_h1("")
}

tryCatch(
    expr = {
        config = prepareEnv()
        main(config)
    },
    error= function(e){
        SendNotification(config, success=FALSE)
        print(rlang::trace_back())
    }
)

