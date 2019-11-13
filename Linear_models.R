

############################################################################################################################################################################################################################
###
##This script is going to display linear models (A ~ B) between all possible pairs of A ~ B for the three genotypes
## across timepoints.  For each A ~ B model, this script will also store the p-values for the interaction term of each Time-X vs Time-Y for
#all times X and Y in a single accession
#The input file is provided.  Note that some of the columns in the file are not used for this analysis
###
############################################################################################################################################################################################################################
 

#install the libraries
library(RColorBrewer)
display.brewer.all()

#read in the data
exampleTableSubset = read.table(file = "exampleTableSubset_Biomass_LinMod.csv", sep = ",", header = TRUE, comment.char = "",quote = "",stringsAsFactors = FALSE, na.strings = "NA")


#read in the data
exampleTableSubset = read.table(file = "mo0_mo4_mo7_biomass_linmod.csv", sep = ",", header = TRUE, comment.char = "",quote = "",stringsAsFactors = FALSE, na.strings = "NA")

#get the column info for dummy variables (Genotype and Stage)
genotypeDummy = exampleTableSubset$Genotype
stageDummy = exampleTableSubset$Stage

#remove all extraneous columns
columnsOfexampleTable = colnames(exampleTableSubset)
columnsOfexampleTable = columnsOfexampleTable[columnsOfexampleTable != "Time"]
columnsOfexampleTable = columnsOfexampleTable[columnsOfexampleTable != "Rep"]
columnsOfexampleTable = columnsOfexampleTable[columnsOfexampleTable != "Replicate"]
columnsOfexampleTable = columnsOfexampleTable[columnsOfexampleTable != "Genotype"]
columnsOfexampleTable = columnsOfexampleTable[columnsOfexampleTable != "Stage"]
columnsOfexampleTable = columnsOfexampleTable[columnsOfexampleTable != "Sample"]
columnsOfexampleTable = columnsOfexampleTable[columnsOfexampleTable != "Protein"]

pValsAllInfo = data.frame("Name" = "", "p-value-slop-difference-R5" = 0, "p-value-slop-difference-R6" = 0, "p-value-slop-difference-R7" = 0, "p-value-slop-difference-R7.5" = 0, "p-value-slop-difference-R8" = 0, stringsAsFactors = FALSE)

#have the table ith all the p-values
pValsAllComparisons = data.frame("Compound 1" = 0, "Compound 2" = 0, "Time 1" = 0, "Time 2" = 0, "Line" = "", "p-value" = 0, "baseline", "relative to baseline change in slope is",stringsAsFactors = FALSE)



####we are going to loop through metabolite-metabolite comparison
##


pdf("LinearModels_SeedDevelopment.pdf")




for(a in 1:length(columnsOfexampleTable))
{
  #all palette available from RColorBrewer
  
  
  for(b in 1:length(columnsOfexampleTable))
  {
  
    #for(c in 1:length(columnsOfexampleTable))
    #{
    
      #won't actually need C ... unless we go back to ratios
      #so could be helpful
      if((a != b) && (b != a))
      {
           
        compoundA = columnsOfexampleTable[a]
        compoundB = columnsOfexampleTable[b]
     
        vecA = exampleTableSubset[, colnames(exampleTableSubset) == compoundA]
        vecB = exampleTableSubset[,colnames(exampleTableSubset) == compoundB ]
    
  
          print(compoundA)
  
	    #make a model for all the A,B pairs at each timepoint
        theSubsetII <- data.frame(ACol=numeric(length(vecA)), B=numeric(length(vecB)),  theIndicatorI = numeric(length(vecA)), theIndicatorII = numeric(length(vecA)))
  
        theSubsetII$theIndicatorI = stageDummy
        theSubsetII$theIndicatorI = as.factor(theSubsetII$theIndicatorI)
  
        theSubsetII$theIndicatorII = genotypeDummy
        theSubsetII$theIndicatorII = as.factor(theSubsetII$theIndicatorII)
  
        theSubsetII$ACol = vecA
        theSubsetII$B = vecB


      

       #go through the SubsetII
  		#which is 
      for(j in 1:length(unique(theSubsetII$theIndicatorII)))
      {
        whichLine = unique(theSubsetII$theIndicatorII)[j]
    
        justLineSubset = theSubsetII[theSubsetII$theIndicatorII == whichLine,]
            
        
         #select the colors that will be used
      	cols<-brewer.pal(n=5,name="Set1")
  		
  		#cols contain the names of four different colors
  		#create a color vector corresponding to levels in the T1 variable in dat
  
        cols_t1<-cols[justLineSubset$theIndicatorI]
        
        #code for the plotting
        plot(ACol ~ B,justLineSubset,col=cols_t1,pch = 16, main=paste(compoundA, " vs ", compoundB, " in ", whichLine), ylab=compoundA, xlab=compoundB)
    
    
       #code for the predictive model
        m<-lm(ACol~theIndicatorI + B*theIndicatorI, justLineSubset)
       
       
       	#code for the predictions
        new_X<-expand.grid(B=seq(0,max(justLineSubset$B),length=10),theIndicatorI=unique(justLineSubset$theIndicatorI))
    
        pred<-predict(m,new_X)
    
        xs<-seq(0,max(justLineSubset$B),length=10)
    
    
    
        #add trend lines
        lines(xs,pred[1:10],col=cols[1],lty=1,lwd=3)
        lines(xs,pred[11:20],col=cols[2],lty=1,lwd=3)
        lines(xs,pred[21:30],col=cols[3],lty=1,lwd=3)
        lines(xs,pred[31:40],col=cols[4],lty=1,lwd=3)
        lines(xs,pred[41:50],col=cols[5],lty=1,lwd=3)
    
    
    		#add legend
        legend(x="topright",legend=rep(unique(justLineSubset$theIndicatorI),times=1),col=rep(cols,times=1),pch=rep(c(16,18),each=5),lwd=1,lty=rep(c(1,2),each=5),bty="n",ncol=2,cex=0.7,pt.cex=0.7)
        
  
    
        #to get the p-value, look for the last four and theIndicatorIIMO7:ACol
        summary(m)
        
        #name
        modelName = paste("correlation between", compoundB, " and ", compoundA, " in ", whichLine)
             
        #the pvalues for each model
        #to store
        vecOfPvalsAdd = c(summary(m)$coefficients[,4][6], summary(m)$coefficients[,4][7],summary(m)$coefficients[,4][8], summary(m)$coefficients[,4][9],summary(m)$coefficients[,4][10])
        
        pValsAllInfo = rbind(pValsAllInfo,c(modelName,vecOfPvalsAdd)) 
        
        
        #in the same Line and , Number: compare modelT1 vs T2 
        for (c in 1:length(unique(justLineSubset$theIndicatorI)))
        {
          
          cat1 = unique(justLineSubset$theIndicatorI)[c]
          print(cat1)
          for (d in 1:length(unique(justLineSubset$theIndicatorI)))
          {
            cat2 = unique(justLineSubset$theIndicatorI)[d]
            
            if((c < d))
            {
              
              #if compound A and B are different, then subset the table
              #by just these times and make a model
              print(cat2)
              part1 = justLineSubset[justLineSubset$theIndicatorI == cat1,]
              part2 = justLineSubset[justLineSubset$theIndicatorI == cat2,]
              
              
              ##bind together the two subsetted tables
              fullTable = rbind(part1, part2)
              
              #just do two at a time
              mJustTwoAtATime<-lm(ACol~theIndicatorI + B*theIndicatorI, fullTable )
              
              #extract the p-value
              myPvalue = summary(mJustTwoAtATime)$coefficients[,4][4]
              
              #get the dummy variable name
              baseline = as.character(unique(fullTable$theIndicatorI)[1])
         
              
              #get the coefficient 
              theCoefficient = summary(mJustTwoAtATime)$coefficients[,1][4]
              
       
              print(baseline)
              print("is the baseline")
              
              #include all of the compound names, coefficients and p-values
              toAddToTheTable = c(compoundA, compoundB, as.character(cat1), as.character(cat2), as.character(whichLine), myPvalue, baseline, theCoefficient)
              
              #now, we add the data for the model to the table
              pValsAllComparisons = rbind(pValsAllComparisons, toAddToTheTable)
            }
          }
        }
       }
      }    
    }
}
dev.off()

#output the table 
write.table(pValsAllComparisons, file = "interaction_pvalues_Table.txt")

