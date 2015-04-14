################################################################################
##########Load necessary Packages
################################################################################
#sorry, but I'm using a lot of packages :/
library(knitr)
library(Hmisc)
library(GenomicRanges)
library(ggplot2)
library(reshape2)
library(stringr)
 # library(chipseq)
library(parallel)
library(rtracklayer)
library(Rsamtools)
library(BSgenome.Dmelanogaster.UCSC.dm3)
library(testthat)
library(pryr)
library(plyr)
library(data.table)
library(dplyr)
library(magrittr)
library(tidyr)




chrs.keep<-seqlevels(Dmelanogaster)[!grepl('U|M', seqlevels(Dmelanogaster))]
si<-seqinfo(Dmelanogaster)[chrs.keep]
bigchrs = chrs.keep[1:6]


################################################################################
##########Table Manipulations
################################################################################




#' This function takes in a table, a bare column name, and a string indicating
#' the seperator. It then puts each seperated value on it's own row, with
#' All other columns duplicated
#' @param x - a table or Granges object
#' @param colname - a column name, bare
#' @param split - a character to split the column by
#' @return a table where the concatenated column has been split into individual rows
#' splitByCol(tbl,columntosplit,seperatorstring)

splitByCol =   function(x,colname,split){


    isgr = is(x,'GenomicRanges')
    if(isgr) x=GR2DT(x)

    colname = substitute(colname)%>%as.character%>%paste0(collapse='')

    stopifnot(colname %in% colnames(x))


    split = strsplit(x[[colname]],split)
    
    x = x[rep(seq_along(split),times = vapply(split,length,666)),]
    
    x[[colname]] = unlist(split)
    
    if(isgr) x%<>%DT2GR
    
    x

  }

#' Convert a GRanges object to a data table
GR2DT = function(gr){
  if(is.data.table(gr)) {
    warning('already a data table')
    return(gr)
  }
  #all columns must be vectors
  # stopifnot( mcols(gr)%>%vapply(is.vector,FUN.VAL=TRUE)%>%all)

  dt = as.data.frame(gr,row.names=NULL,optional=FALSE)%>%as.data.table
  dt$strand= as.character(dt$strand)
  # setkey(dt,seqnames,strand,start,end)
  dt

  expect_true( 
    all(colnames( mcols(gr)) %in% colnames(dt) ),
    "no columns lost"
  )
  dt
}
gr2dt=GR2DT

#' Convert a data table to a GRanges object
#' By default it checks that the GRanges object
#' has the chromosomes in si
DT2GR = function(dt,seqinf=si,checksi=TRUE){

  if(is(dt,'GenomicRanges')) {
    warning('already a GRanges Object')
    return(dt)
  }


  stopifnot(w('seqnames start')%in%colnames(dt))
  stopifnot(w('width')%in%colnames(dt)|w('end')%in%colnames(dt))
  if(checksi){stopifnot(dt[['seqnames']] %in% seqlevels(seqinf))
  }else{seqinf=NULL}
  
  hasend=FALSE
  haswidth=FALSE

  if(! is.null(dt[['end']])){
    stopifnot (dt[['end']] %>% is_greater_than(0) %>%all)
    hasend=TRUE
  }
  if(! is.null(dt[['width']])){
    stopifnot (dt[['width']] %>% is_greater_than(0) %>%all)
    haswidth=TRUE
  }
  
  stopifnot(dt[['start']] %>% is.numeric)
  stopifnot(hasend|haswidth )
  
  if(haswidth & ! hasend ){
    dt[['end']]  = dt[['start']]+dt[['width']]-1 
  } 
  if(hasend ){

  } 

  #strand
  if(! 'strand' %in% colnames(dt)){
    dt[['strand']] ='*'
  }

  stopifnot(dt[['strand']] %in% c('+','-','*'))
  



  mdatcols = colnames(dt) %>% setdiff(w('seqnames start width strand end')) 
  #create basic granges object
  if(checksi){
    gr=GRanges(dt[['seqnames']],IRanges(dt[['start']],dt[['end']]),strand=dt[['strand']],seqinfo=seqinf)
  }else{    gr=GRanges(dt[['seqnames']],IRanges(dt[['start']],dt[['end']]),strand=dt[['strand']])}

  #add in the metadata if we need to
  if(length(mdatcols)){
    if(is.data.table(dt)){ mcols(gr) = dt[,mdatcols,with=FALSE]%>%as("DataFrame")
    }else{ mcols(gr) = dt[,mdatcols]%>%as("DataFrame")}
  }

  expect_true( 
    all(colnames(dt) %in% c(colnames(mcols(gr)),'seqnames','start','end','width' ,'strand')),
    "no columns lost"
  )

  gr
}



