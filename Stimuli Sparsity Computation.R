
#  .o oOOOOOOOo                                            OOOo
#   Ob.OOOOOOOo  OOOo.      oOOo.                      .adOOOOOOO
#   OboO"""""""""""".OOo. .oOOOOOo.    OOOo.oOOOOOo.."""""""""'OO
#   OOP.oOOOOOOOOOOO "POOOOOOOOOOOo.   `"OOOOOOOOOP,OOOOOOOOOOOB'
#   `O'OOOO'     `OOOOo"OOOOOOOOOOO` .adOOOOOOOOO"oOOO'    `OOOOo
#   .OOOO'            `OOOOOOOOOOOOOOOOOOOOOOOOOO'            `OO
#   OOOOO                 '"OOOOOOOOOOOOOOOO"`                oOO
#  oOOOOOba.                .adOOOOOOOOOOba               .adOOOOo.
# oOOOOOOOOOOOOOba.    .adOOOOOOOOOO@^OOOOOOOba.     .adOOOOOOOOOOOO
#OOOOOOOOOOOOOOOOO.OOOOOOOOOOOOOO"`  '"OOOOOOOOOOOOO.OOOOOOOOOOOOOO
#OOOO"       "YOoOOOOMOIONODOO"`  .   '"OOROAOPOEOOOoOY"     "OOO"
#   Y           'OOOOOOOOOOOOOO: .oOOo. :OOOOOOOOOOO?'         :`
#   :            .oO%OOOOOOOOOOo.OOOOOO.oOOOOOOOOOOOO?         .
#   .            oOOP"%OOOOOOOOoOOOOOOO?oOOOOO?OOOO"OOo
#                '%o  OOOO"%OOOO%"%OOOOO"OOOOOO"OOO':
#                     `$"  `OOOO' `O"Y ' `OOOO'  o             .
#   .                  .     OP"          : o     .
#                            :
#
#
#-Stimuli Sparsity Calculations
#-RYAN BOYD 
#-UNIVERSITY OF TEXAS AT AUSTIN
#
#----Built for R 3.3.2, 2017-06-12
#----Version 1.00
#
#----Requires the MASS and smacof packages
#----uncomment the line below to install both packages if
#----they are not already installed
#
#install.packages(c('MASS', 'smacof'), dependencies=T)
#
#----This particular script also uses the 'gtools' package for
#----data restructuring. This is not necessary to run the
#----function itself, but it is necessary to run this version
#----of the script using the data it was designed for
#
#install.packages('gtools')






# 
# 
#  _____        __ _              _____                 _ _            _____      _        ______                _   _             
# |  __ \      / _(_)            |  __ \               (_) |          / ____|    | |      |  ____|              | | (_)            
# | |  | | ___| |_ _ _ __   ___  | |  | | ___ _ __  ___ _| |_ _   _  | |     __ _| | ___  | |__ _   _ _ __   ___| |_ _  ___  _ __  
# | |  | |/ _ \  _| | '_ \ / _ \ | |  | |/ _ \ '_ \/ __| | __| | | | | |    / _` | |/ __| |  __| | | | '_ \ / __| __| |/ _ \| '_ \ 
# | |__| |  __/ | | | | | |  __/ | |__| |  __/ | | \__ \ | |_| |_| | | |___| (_| | | (__  | |  | |_| | | | | (__| |_| | (_) | | | |
# |_____/ \___|_| |_|_| |_|\___| |_____/ \___|_| |_|___/_|\__|\__, |  \_____\__,_|_|\___| |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|
#                                                              __/ |                                                               
#                                                             |___/                                                                
# 

#this is the function that does all of the heavy lifting
#this function will do some basic matrix completion (replacement of missing values with median)
#and uses the smacof algorithm for multidimensional scaling
#it will also automatically create a symmetric matrix from the original dataset,
#even if the original data is not presented in an entirely symmetric format
#it will also generate a plot of the first 2 dimensions of the MDS result, with the
#MDS stress value appended to the top
#
#this function can also do any combination of diplaying the plot in an R window
#and saving a copy of the plot to disk. one option does not rely upon the other,
#so it is possible to save the plot without displaying it and vice versa
#
#as a note, "datatype" should be "ratio" by default for this script. While the
#input data that you are providing might be interval data, the similarity matrix
#will convert this to ratio data in that it averages out the stimulus ratings
#across participants
CompSparsity = function(df, Dimensions, stim1, stim2, stim1prop,
                       showPlot=F, MDSmaxiter=500, datatype='ratio',
                       stim2prop, stim_ratings, OutputFilename=NA,
                       OutputPlotFilename=NAj, AvgMatrixOutputFilename=NA){

  
    library(MASS)
    library(smacof)
  

    #get the row/col names, which we'll use later for both
    #naming and positioning information
    
    cols = unique(as.character(df[,stim1]))
    rows = unique(as.character(df[,stim2]))
    
    #we want a nice, symmetric
    cols = unique(c(cols, rows))
    rows = cols
    
    stim_properties = data.frame(c(as.character(df[,stim1]),
                                   as.character(df[,stim2])),
                                 c(as.character(df[,stim1prop]),
                                   as.character(df[,stim2prop])))
    stim_properties = unique(stim_properties)
    names(stim_properties) = c('stim', 'stim_prop')
    

    
    
    
    #initialize our matrix
    StimulusMatrix = matrix(0, nrow = length(rows), ncol = length(cols))
    rownames(StimulusMatrix) = rows
    colnames(StimulusMatrix) = cols

    StimulusMatrix_Counts = matrix(0, nrow = length(rows), ncol = length(cols))
    
    cat("\n\n")
    cat("Computing Mean Stimulus Rating Matrix...")
    cat("\n")
    pb = txtProgressBar(min=0, max=nrow(df),
                        initial=0, char='*', width=45,
                        style=3, label="") 
    
    
    
    attach(df)
    for(i in 1:nrow(df)){
      #add co-occurrence weights to our overall matrix
      matrix_col = match(df[i, stim1], cols)
      matrix_row = match(df[i, stim2], rows)
      
      #note that "stim_ratings" for the data for which this whole function was made are "stimulus similarity" ratings
      #depending on the nature of your stimulus ratings, you may need to reverse-code or change the nature of your
      #data prior to running this entire function
      StimulusMatrix[matrix_row, matrix_col] = StimulusMatrix[matrix_row, matrix_col] + df[i, stim_ratings]
      StimulusMatrix_Counts[matrix_row, matrix_col] = StimulusMatrix_Counts[matrix_row, matrix_col] + 1
      
      setTxtProgressBar(pb,i)
      
    }
    detach(df)
    
    #add in the transposed values so that it's a congruent matrix
    StimulusMatrix = StimulusMatrix + t(StimulusMatrix)
    
    StimulusMatrix_Counts = StimulusMatrix_Counts + t(StimulusMatrix_Counts)
    StimulusMatrix_Counts[StimulusMatrix_Counts == 0] = NA
    #get the averages
    StimulusMatrix = StimulusMatrix / StimulusMatrix_Counts
    
    
    #if there are missing values for stimuli pairs, we replace them
    #with the median similarity from the entire matrix
    #this is because we cannot perform an MDS with missing values, and
    #we don't want to just set them to 0, as this could seriously deteriorate
    #the results from the MDS procedure
    StimulusMatrix[is.na(StimulusMatrix)] = median(StimulusMatrix, na.rm=T)
    #the diagonal should always be 0
    diag(StimulusMatrix) = 0
    
    rm(StimulusMatrix_Counts)
    
    
    
    #perform scaling
    #mds_result = isoMDS(d=as.dist(dist(StimulusMatrix)), k=Dimensions, maxit=500)
    mds_result = smacofSym(delta=dist(StimulusMatrix),
                           ndim = Dimensions, itmax=MDSmaxiter,
                           type=datatype)
    stress = mds_result$stress
    mds_result = as.data.frame(mds_result$conf)
    
    

    #plot this stuff
    x <- mds_result[,1]
    y <- mds_result[,2]
    
    if(showPlot==T){
      plot(x, y, xlab="Dim 1", ylab="Dim 2",
           main=paste("Metric MDS Locations (Stress = ", round(stress, digits=5), ")", sep=""),	type="n")
      text(x, y, labels = row.names(StimulusMatrix), cex=.7)
    }
    
    if(!is.na(OutputPlotFilename)){
      png(filename=as.character(OutputPlotFilename))
      plot(x, y, xlab="Dim 1", ylab="Dim 2", 
           main=paste("Metric MDS Locations (Stress = ", round(stress, digits=5), ")", sep=""),	type="n")
      text(x, y, labels = row.names(StimulusMatrix), cex=.7)
      dev.off()
    }
    
    rm(x, y)
    

    

    #mds_result$stress = stress
    #compute each row's distance from each other row
    mds_result$Sparsity = NA
    
    #make sure that we set up the Stimulus properties
    mds_result$stim = row.names(mds_result)
    mds_result = merge(mds_result, stim_properties, by='stim', all=T)



    cat("\n\n")
    cat("Calculating Sparsities...")
    cat("\n")
    pb = txtProgressBar(min=0, max=nrow(mds_result), initial=0,
                        char='*', width=45, style=3, label="") 
    #goes through and calculates sparsity    
    for(i in 1:nrow(mds_result)){
      euc_dist = 0
      
      
      Condition_Matches = 0
      for(j in 1:nrow(mds_result)){
        
          #we only want to calculate the euclidean distance to other items with
          #the same Stimulus property -- e.g., only within the same valence, etc.
          if((i != j) && (mds_result$stim_prop[i] == mds_result$stim_prop[j])){
            euc_dist = euc_dist + as.numeric(dist(rbind(mds_result[i, 2:(Dimensions+1)],
                                                        mds_result[j, 2:(Dimensions+1)])))
            #track the number in case the property distribution is imbalanced -- make no assumptions
            Condition_Matches = Condition_Matches + 1
          }
        
      }

      euc_dist = euc_dist / Condition_Matches
      mds_result$Sparsity[i] = euc_dist 
      setTxtProgressBar(pb,i)

    }
    
    #sort by stimulus property
    mds_result = mds_result[with(mds_result, order(stim_prop)), ]
    
    
    
    cat("\n\n")
    cat("Returning result...\n")

    if(!is.na(OutputFilename)){
      write.csv(mds_result, as.character(OutputFilename),
                row.names=F, na="")
    }
    
    if(!is.na(AvgMatrixOutputFilename)){
      write.csv(StimulusMatrix, as.character(AvgMatrixOutputFilename),
                row.names=T, na="")
    }
    
    return(mds_result)

}












# 
#  _____                _       __  _____                                 _____        _        
# |  __ \              | |     / / |  __ \                               |  __ \      | |       
# | |__) |___  __ _  __| |    / /  | |__) | __ ___ _ __   __ _ _ __ ___  | |  | | __ _| |_ __ _ 
# |  _  // _ \/ _` |/ _` |   / /   |  ___/ '__/ _ \ '_ \ / _` | '__/ _ \ | |  | |/ _` | __/ _` |
# | | \ \  __/ (_| | (_| |  / /    | |   | | |  __/ |_) | (_| | | |  __/ | |__| | (_| | || (_| |
# |_|  \_\___|\__,_|\__,_| /_/     |_|   |_|  \___| .__/ \__,_|_|  \___| |_____/ \__,_|\__\__,_|

#                                                 | |                                           
#                                                 |_|                                           
                                                                                               

#set up the working directory and read in your dataset
setwd('X:/Directory')
df = read.table('sample data - simulated.txt', sep='\t', header=T)



#sets up a couple of lists so that we can separate out the conditions
#this is dataset-specific and won't apply to all datasets
selfconditions = c('self_pos_neg', 'self_neg_pos', 'self_pos_pos', 'self_neg_neg')
otherconditions = c('other_pos_neg', 'other_neg_pos', 'other_pos_pos', 'other_neg_neg')



#set up the items to be able to identify which factor they belong to
#for example, "sp" is "self-positive" and "on" is "other-negative"
#again, this is dataset specific
#however, you will need to pass your dataframe PLUS a total of 5 columns
#to the function
#------2 columns with stimuli labels/names
#------2 columns with stimuli property labels (e.g., "postive" and "negative")
#------1 column with the participant similarity ratings of the stimuli pairs


#for this dataset, we have the stimuli and the trial conditions, but not actual columns
#with the stimuli properties. As such, we make 2 new variables in the dataset and assign
#stimulus properties as a function of the trial conditions, which have the stimulus property
#information embedded within them
df$stimprop2 = NA
df$stimprop3 = NA

df$stimprop2 = ifelse(df$trialcode == 'self_pos_neg', 'self',  as.character(df$stimprop2))
df$stimprop3 = ifelse(df$trialcode == 'self_pos_neg', 'self',  as.character(df$stimprop3))
df$stimprop2 = ifelse(df$trialcode == 'self_neg_pos', 'self',  as.character(df$stimprop2))
df$stimprop3 = ifelse(df$trialcode == 'self_neg_pos', 'self',  as.character(df$stimprop3))
df$stimprop2 = ifelse(df$trialcode == 'self_pos_pos', 'self',  as.character(df$stimprop2))
df$stimprop3 = ifelse(df$trialcode == 'self_pos_pos', 'self',  as.character(df$stimprop3))
df$stimprop2 = ifelse(df$trialcode == 'self_neg_neg', 'self',  as.character(df$stimprop2))
df$stimprop3 = ifelse(df$trialcode == 'self_neg_neg', 'self',  as.character(df$stimprop3))

df$stimprop2 = ifelse(df$trialcode == 'other_pos_neg', 'other',  as.character(df$stimprop2))
df$stimprop3 = ifelse(df$trialcode == 'other_pos_neg', 'other',  as.character(df$stimprop3))
df$stimprop2 = ifelse(df$trialcode == 'other_neg_pos', 'other',  as.character(df$stimprop2))
df$stimprop3 = ifelse(df$trialcode == 'other_neg_pos', 'other',  as.character(df$stimprop3))
df$stimprop2 = ifelse(df$trialcode == 'other_pos_pos', 'other',  as.character(df$stimprop2))
df$stimprop3 = ifelse(df$trialcode == 'other_pos_pos', 'other',  as.character(df$stimprop3))
df$stimprop2 = ifelse(df$trialcode == 'other_neg_neg', 'other',  as.character(df$stimprop2))
df$stimprop3 = ifelse(df$trialcode == 'other_neg_neg', 'other',  as.character(df$stimprop3))



pos_stimuli = c("SINCERE", "HONEST", "UNDERSTANDING", "LOYAL", "TRUTHFUL",
                "TRUSTWORTHY", "INTELLIGENT", "DEPENDABLE", "HAPPY", "CLEAN",
                "WISE", "CONSIDERATE", "GOOD-NATURED", "RELIABLE", "MATURE",
                "WARM", "EARNEST", "KIND", "FRIENDLY", "INTERESTING")

neg_stimuli = c("CRUEL", "DISRESPECTFUL", "OBNOXIOUS", "DISHONORABLE",
                "MALICIOUS", "DECEITFUL", "ANNOYING", "DISLIKABLE", "HOSTILE",
                "INSULTING", "GREEDY", "SPITEFUL", "CONCEITED", "RUDE", "THOUGHTLESS",
                "HEARTLESS", "INSOLENT", "VULGAR", "NARROW-MINDED", "SELFISH")


df$stimprop2 = ifelse(as.character(df$stimulusitem2) %in% pos_stimuli,
                      paste(as.character(df$stimprop2), "_pos", sep=''),
                      paste(as.character(df$stimprop2), "_neg", sep=""))

df$stimprop3 = ifelse(as.character(df$stimulusitem3) %in% pos_stimuli,
                      paste(as.character(df$stimprop3), "_pos", sep=''),
                      paste(as.character(df$stimprop3), "_neg", sep=""))




#  _____                                    _                     
# |  __ \                 /\               | |                    
# | |__) |   _ _ __      /  \   _ __   __ _| |_   _ ___  ___  ___ 
# |  _  / | | | '_ \    / /\ \ | '_ \ / _` | | | | / __|/ _ \/ __|
# | | \ \ |_| | | | |  / ____ \| | | | (_| | | |_| \__ \  __/\__ \
# |_|  \_\__,_|_| |_| /_/    \_\_| |_|\__,_|_|\__, |___/\___||___/
#                                              __/ |              
#                                             |___/               
                                           



#now we can actually run the code
#this first line just does it for data in the "self" condition
Self_Sparsity = CompSparsity(df = subset(df, trialcode %in% selfconditions),
                           Dimensions=3, stim1="stimulusitem2", stim2="stimulusitem3",
                           stim1prop="stimprop2", stim2prop="stimprop3", stim_ratings="response",
                           OutputFilename="Self_Sparsity.csv", OutputPlotFilename="Self-Sparsity Plot.png",
                           AvgMatrixOutputFilename="Self_Rating_Avg_Matrix.csv",
                           showPlot=T, MDSmaxiter=500, datatype='ratio')


Other_Sparsity = CompSparsity(df = subset(df, trialcode %in% otherconditions),
                           Dimensions=3, stim1="stimulusitem2", stim2="stimulusitem3",
                           stim1prop="stimprop2", stim2prop="stimprop3", stim_ratings="response",
                           OutputFilename="Other_Sparsity.csv", OutputPlotFilename="Other-Sparsity Plot.png",
                           AvgMatrixOutputFilename="Other_Rating_Avg_Matrix.csv",
                           showPlot=T, MDSmaxiter=500, datatype='ratio')



#this is also specific to this dataset -- 
#combine the self and other information
#...then save it to disk
library(gtools)
Sparsity_Combined = smartbind(Self_Sparsity, Other_Sparsity)

write.csv(Sparsity_Combined, "All Sparsity Information.csv", row.names=F, na='')