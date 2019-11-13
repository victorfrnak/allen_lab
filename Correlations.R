

############################################################################################################################################
##This script is going to look at all the time-compound correlations
#for all three of the genotypes FN301952,"Will82, FN300012 Mo7, Mo0, Mo4
###in a single line using the custon corstars() function
#################################################################################################################################

# Libraries

#install.packages("tidyverse")
library(tidyverse)

#install.packages("hrbrthemes")
library(hrbrthemes)

#install.packages("viridis")
library(viridis)

#install.packages("plotly")
library(plotly)

#install.packages("d3heatmap")
library(d3heatmap)
​


#install.packages("heatmaply")
library(heatmaply)


#custom function corstars (Note that I didn't write, it's available from: 
#http://www.sthda.com/english/wiki/elegant-correlation-table-using-xtable-r-package

# x is a matrix containing the data
# method : correlation method. "pearson"" or "spearman"" is supported
# removeTriangle : remove upper or lower triangle
# results :  if "html" or "latex"
  # the results will be displayed in html or latex format
corstars <-function(x, method=c("pearson", "spearman"), removeTriangle=c("upper", "lower"),
                     result=c("none", "html", "latex")){
    #Compute correlation matrix
    require(Hmisc)
    x <- as.matrix(x)
    correlation_matrix<-rcorr(x, type=method[1])
    R <- correlation_matrix$r # Matrix of correlation coeficients
    p <- correlation_matrix$P # Matrix of p-value 
    
    ## Define notions for significance levels; spacing is important.
    mystars <- ifelse(p < .0001, "****", ifelse(p < .001, "*** ", ifelse(p < .01, "**  ", ifelse(p < .05, "*   ", "    "))))
    
    ## trunctuate the correlation matrix to two decimal
    R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[,-1]
    
    ## build a new matrix that includes the correlations with their apropriate stars
    Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(x))
    diag(Rnew) <- paste(diag(R), " ", sep="")
    rownames(Rnew) <- colnames(x)
    colnames(Rnew) <- paste(colnames(x), "", sep="")
    
    ## remove upper triangle of correlation matrix
    if(removeTriangle[1]=="upper"){
      Rnew <- as.matrix(Rnew)
      Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
      Rnew <- as.data.frame(Rnew)
    }
    
    ## remove lower triangle of correlation matrix
    else if(removeTriangle[1]=="lower"){
      Rnew <- as.matrix(Rnew)
      Rnew[lower.tri(Rnew, diag = TRUE)] <- ""
      Rnew <- as.data.frame(Rnew)
    }
    
    ## remove last column and return the correlation matrix
    Rnew <- cbind(Rnew[1:length(Rnew)-1])
    if (result[1]=="none") return(Rnew)
    else{
      if(result[1]=="html") print(xtable(Rnew), type="html")
      else print(xtable(Rnew), type="latex") 
    }
} 


allThreeGenotypes = read.table(file = "data_for_correlations.csv", sep = ",", header = TRUE, comment.char = "",quote = "",stringsAsFactors = FALSE, na.strings = "NA")


##note well that variable names are given using the alternative (and shorter) names
##i.e.  MO7 = FN30195, MO0 = Will82, MO4 = FN300012



#FN301952/MO7
measurementDataFN301952  = allThreeGenotypes[grepl("FN301952", allThreeGenotypes$X),] 
rownames(measurementDataFN301952) = measurementDataFN301952$X
measurementDataFN301952$X = NULL

#Will82/MO0
measurementDataWill82  = allThreeGenotypes[grepl("Will82", allThreeGenotypes$X),] 
rownames(measurementDataWill82) = measurementDataWill82$X
measurementDataWill82$X = NULL

#MO4/FN300012
measurementFN300012  = allThreeGenotypes[grepl("FN300012", allThreeGenotypes$X),] 
rownames(measurementFN300012) = measurementFN300012$X
measurementFN300012$X = NULL


#output all of the correlation matrices into their
#respective files
write.table(corstars(measurementFN300012), file = "MO4_FN300012_corstars.txt", sep = "\t")
write.table(corstars(measurementDataWill82 ), file = "MO0_Will82_corstars.txt", sep = "\t")
write.table(corstars(measurementDataFN301952), file = "MO7_FN_301952_corstars.txt", sep = "\t")

