#SirMixaPlot version 8.3 - Less is more
#----------------------------------------------------------------------------------------------------------

#This program is designed to take a csv (comma separated values) file and generate scatterplots based on columns.

# Changes from version 8.2:
#   -Removed from redundant functions. Mostly the gate() and negPopCorrect since they no longer serve a purpose
#   -Streamlined some functions to remove unnecessary global variables
#   -Changed normalizer to use LMean_ANC_edu instead of the now defunct ALIMean_ANC_edu

# Changes from version 8.1:
#   -Changed the default EdU neg setting on normalizer()
#   -Added a work-around for concater() should rbinding fail due to missing columns
#   -imaGen() now also creates a concatanated csv of all the ROIs in the WholeCell folder with the _ROI_ tag
#   -Added runPCA() function to take a dataframe and generate a PCA biplot; still in beta testing

#Changes from SirMixaPlot version 8.0:
#   -Added some comments and typos
#   -Added the reMap() function for re-creating a facsimile of the original image

#Changes from SirMixaPlot version 7.0:
#   -Takes data in from the YggData.imj macro
#   -imaGen() is now much bigger and generates combined CSVs for nuclear vs whole cell comparisons

#Things still that need to be coded:

#Program opens required packages
#If these packages aren't install, use the '0_setup.R' function to install them.

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(RecordLinkage))
#suppressPackageStartupMessages(library(ggbiplot))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(plotly))



#------------------------------------------------------------------------------------------------------------

# Functions for your pleasure

#----------------------------------------------------------------
#sirmixaplot() loads the csv file (filo) and will call joiner() to start graphing if run is TRUE. This may lead to negPopCorrect if the appropriate names aren't present
#To re-run the negPopCorrect, simply call adjust as TRUE
sirmixaplot <- function(filo, 
                        run=TRUE, 
                        adjust=TRUE){
  #The script, m'lord
  
  #This is the beginning of the readline prompts for the program. User inputs desired output conditionals. Also, certain columns are added for safety...
  thingee <<- filo
  cat(paste("Opening file: ", thingee, sep=""))
  cat("\n")
  
  #This checks if the file has the '_cells', which signals it has been previously processed. Otherwise it extracts the appropriate rows and creates a cells file
  if (grepl("_cells.csv", thingee, fixed = TRUE)){
    cells <<- read_csv(file=thingee)
  } else {
    cells <<- read_csv(file=thingee)
    cat("Creating a new file, sir or madam.")
    cat("\n")
    thingee <<- paste0(substr(thingee, 1, nchar(thingee)-4), "_cells.csv")
    cat(paste(thingee, "created."))
    if ("IntDen_ANC_dna" %in% names(cells) & "IntDen_ANC_anchor" %in% names(cells)){
      log2Name <- readline(prompt = "Both 'DNA' and 'anchor' found. Which should be used for log2_dna? (ANCHOR/dna)? ")
      if (log2Name == "dna"){
        cells$log2_dna <<- log(cells$IntDen_ANC_dna, 2)
      } else {
        cells$log2_dna <<- log(cells$IntDen_ANC_anchor, 2)
      }
    } else if ("IntDen_ANC_dna" %in% names(cells)){
      cells$log2_dna <<- log(cells$IntDen_ANC_dna, 2)
    }
    
    for (i in 1:ncol(cells)){
      if (grepl("Mean_", names(cells)[i])){
        cells[paste0("L", names(cells)[i])] <<- log(cells[,i], 10)
      }
    }
    #writes the file
    write_csv(cells, thingee)
  }
  cat("\n")
  
  #Checks that some size correctiong has occured and begins making the graph
  if(run==TRUE){
    joiner()
    #writes the file
    write_csv(cells, thingee)
  }
  if(!"dna_norm" %in% names(cells) & TRUE %in% grepl("ANC_edu", names(cells))){
    normalizer()
    #writes the file
    write_csv(cells, thingee)
  }
  
}

#--------------------------------------------------------------------------------------
#joiner() connects the sirmixaplot() with grapho() and can be used to graph new plots without running through sirmixaplot(). Default dataset is cells but any dataset can be used
joiner <- function(df = cells){
  #newAxis is reset for the changeAxis() function
  newAxis <<- FALSE
  
  #cats the dataset columns
  cat("\n")
  cat("Data categories detected:")
  cat("\n")
  print(sort(names(df)))
  cat("\n")
  
  #This sets the x and y variables
  px<<-readline(prompt = "px - What parameter is the x-axis: ")
  xname<<-readline(prompt = "What is the name of the x axis: ")
  cat("\n")
  py<<-readline(prompt = "py - What parameter is the y-axis: ")
  yname<<-readline(prompt = "What is the name of the y axis: ")
  cat("Working...")
  cat("\n")
  
  #calls the graphing function
  grapho(df = cells,
         X = px, 
         Y = py, 
         Xn = xname,
         Yn = yname)
}

#grapho() actually makes the scatter plots
grapho <- function(df = cells, 
                   X = px, 
                   Y = py, 
                   Xn = xname,
                   Yn = yname){
  
  #backup dataframe is generated
  interim <<- df

  cat("Data frame generated")
  cat("\n")
  
  #Initial scatter plot is generated
  qq<<-ggplot(data=df, aes_string(X, Y))+
    theme_linedraw()+
    theme(plot.title = element_text(size=40, face="bold"))+
    theme(axis.title = element_text(size=30, face="bold"))+
    geom_point(size=2)+
    scale_x_continuous(name=Xn)+
    scale_y_continuous(name=Yn)+
    theme(axis.text = element_text(face='bold', size=16), axis.ticks = element_blank())
  
  if(!exists("newAxis")){
    newAxis <<- FALSE
  }
  if(newAxis == TRUE){
    changeAxis(store[1], store[2], store[3], store[4])
  }
  
  #Program finishes up
  cat("Plot generated")
  cat("\n")
  print(qq)
}

#----------------------------------------------------------------------------------
#This function creates quadrants and calculates the percentage in each quadrant
#   Geom_line function can be altered to change the aesthetics of the quadrant lines

modthequad <- function(xi = 1.5,
                        yi = 1.5,
                        df = cells,
                        saveit = F, 
                        pause = T, 
                        graphit = T){
  if (pause == T){
    xi<-readline(prompt = "What is the x-intercept: ")
    yi<-readline(prompt = "What is the y-intercept: ")
  }
  quad <- c("Q1", "Q2", "Q3", "Q4")
  todos <- nrow(df)
  qt1 <<- nrow(subset(df, df[px] < xi & df[py] > yi))
  qt2 <<- nrow(subset(df, df[px] > xi & df[py] > yi))
  qt3 <<- nrow(subset(df, df[px] < xi & df[py] < yi))
  qt4 <<- nrow(subset(df, df[px] > xi & df[py] < yi))
  PercentTotal <- c(format(round(qt1/todos*100, 2), nsmall = 2), format(round(qt2/todos*100, 2), nsmall = 2), format(round(qt3/todos*100, 2), nsmall = 2), format(round(qt4/todos*100, 2), nsmall = 2))
  quads <- data.frame(quad, PercentTotal)
  xmin <- min(df[px])
  xmax <- max(df[px])
  ymin <- min(df[py])
  ymax <- max(df[py])
  qqmod <<- qq
  qqmod <<- qqmod+
    geom_hline(yintercept = as.numeric(yi), color = "red", linetype = "dashed", linewidth = 2)+
    geom_vline(xintercept = as.numeric(xi), color = "red", linetype = "dashed", linewidth = 2)+
    annotate("text", x = xmin*1.05, y = ymax*0.95, label = paste0("Q1: ", as.character(format(round((qt1/todos*100), 2), nsmall = 2)), "%"))+
    annotate("text", x = xmax*0.95, y = ymax*0.95, label = paste0("Q2: ", as.character(format(round((qt2/todos*100), 2), nsmall = 2)), "%"))+
    annotate("text", x = xmin*1.05, y = ymin*1.05, label = paste0("Q3: ", as.character(format(round((qt3/todos*100), 2), nsmall = 2)), "%"))+
    annotate("text", x = xmax*0.95, y = ymin*1.05, label = paste0("Q4: ", as.character(format(round((qt4/todos*100), 2), nsmall = 2)), "%"))
  print(quads)
  if(graphit==T){
    print(qqmod+geom_density_2d())
    if(saveit==T){
      ggsave(paste0(px,"_by_",py,"_mod.pdf"))
    }
    cat("\n")
    cat("call 'qqmod' to see plot with quadrants.")
    cat("\n")
  }
}

#---------------------------------------------------------------------------
#This function will remove certain points from the cells dataframe.
cull <- function(){
  print(sort(names(cells)))
  para <- readline(prompt = "Which parameter should be culled: ")
  bORs <- readline(prompt = "Data with values (b)igger, (s)maller, or (=) to the target: ")
  vale <- readline(prompt = "What is the target value: ")

  if (bORs=="b"){
    #print("culling bigger")
    cells <<- cells[unlist(cells[,para]) < as.numeric(vale),]
  } else if (bORs=="s"){
    #print("culling smaller")
    cells <<- cells[unlist(cells[,para]) > as.numeric(vale),]
  } else if (bORs=="="){
    #print("culling equal")
    cells <<- cells <<- cells[unlist(cells[,para]) == as.numeric(vale),]
  } else{
    grapho(cells)
  }
  grapho(cells)
  cat("\n")
}

#----------------------------------------------------------------
#This function will automatically pseudo-color density of a cell cycle profile and provides the basic functions for density pseudo-coloring
get_density <- function(px, 
                        py,
                        ...){
  theme_set(theme_bw(base_size = 16))
  dens <<- MASS::kde2d(px, py, ...)
  ix <<- findInterval(px, dens$x)
  iy <<- findInterval(py, dens$y)
  ii <<- cbind(ix, iy)
  return(dens$z[ii])
}

color_density<-function(){
  cells$density <<- get_density(cells[px][,1], cells[py][,1], n=100)
  qq <<- ggplot(data = cells, aes_string(px, py)) + geom_point(aes(color = density)) + scale_color_viridis()+xlab(xname)+ylab(yname)
  
  if(newAxis == TRUE){
    changeAxis(store[1], store[2], store[3], store[4])
  }
  print(qq)
}


#----------------------------------------------------------------
#This function asks you to define a G1 peak and re-assigns that value to  '1'
normalizer <- function(){
  
  #create a column holder (just in case...)
  cells$edu<<-"Negative"
  cells$ploidy<<-"2N"
  if(T %in% is.na(cells$LMean_ANC_edu)){
    print(paste0("Detected ", nrow(cells[is.na(cells$LMean_ANC_edu),]), " NA values in LMean_ANC_edu (", round(100*nrow(cells[is.na(cells$LMean_ANC_edu),])/nrow(cells), 2),"% of tota cells)."))
    remove_nas <- readline(prompt = "Should these be removed? (Y/n) ")
    if(remove_nas != "n"){
      print("Removing rows with NA values")
      cells <<- cells[!is.na(cells$LMean_ANC_edu),]
    }
  }
  if(T %in% is.na(cells$log2_dna)){
    print(paste0("Detected ", nrow(cells[is.na(cells$log2_dna),]), " NA values in log2_dna (", round(100*nrow(cells[is.na(cells$log2_dna),])/nrow(cells), 2),"% of tota cells)."))
    remove_nas <- readline(prompt = "Should these be removed? (Y/n) ")
    if(remove_nas != "n"){
      print("Removing rows with NA values")
      cells <<- cells[!is.na(cells$log2_dna),]
    }
  }
  cat("Determining the G1 population now.")
  cat("\n")
  
  #This part calls an EdU negative popuation and creates the normalized values
  bigBin <- which.max(density(cells[!is.na(cells$LMean_ANC_edu) & !is.infinite(cells$LMean_ANC_edu),]$LMean_ANC_edu)$y)
  edu_neg <- density(cells[!is.na(cells$LMean_ANC_edu) & !is.infinite(cells$LMean_ANC_edu),]$LMean_ANC_edu)$x[bigBin]+0.1
  edu_hist <- ggplot(cells, aes(LMean_ANC_edu))+geom_density()+geom_vline(xintercept = edu_neg)+xlab("Log EdU")
  print(edu_hist)
  g1_good <- readline(prompt = paste0("Is ", format(round(edu_neg, 2), nsmall = 2), " representative of the EdU cutoff? (Y/n) "))
  if (g1_good == "n"){
    edu_neg = as.numeric(readline(prompt = "What value should EdU be cutoff as negative? "))
    edu_hist <- ggplot(cells, aes(LMean_ANC_edu))+geom_density()+geom_vline(xintercept = edu_neg)+xlab("Log EdU")
    print(edu_hist)
  }
  cat(paste0("Using ", format(round(edu_neg, 2), nsmall = 2), " as the EdU cutoff."))
  cat("\n")

  cells[cells$LMean_ANC_edu > edu_neg & !is.na(cells$LMean_ANC_edu),]$edu <<- "Positive"
  
  eduNegCells <<- subset(cells, LMean_ANC_edu < edu_neg)
  #print(head(eduNegCells))
  cells$edu_norm <<- (cells$LMean_ANC_edu+1)-mean(eduNegCells$LMean_ANC_edu)
  
  #This part calls the 2N peak that is EdU-negative   
  bigBin <- which.max(density(eduNegCells$log2_dna)$y)
  diploid <- density(eduNegCells$log2_dna)$x[bigBin]
  dna_hist <- ggplot(eduNegCells, aes(log2_dna))+geom_density()+geom_vline(xintercept = diploid, color = "green", size = 1.5)+xlab("DNA content")
  print(dna_hist)
  g1_good <- readline(prompt = paste0("Is ", format(round(diploid, 2), nsmall = 2), " representative of the 2N peak? (Y/n/manual) "))
  if (g1_good == "n"){
    diploid <- as.numeric(readline(prompt = "What is the value of the 2N peak? "))
    dna_hist <- ggplot(eduNegCells, aes(log2_dna))+geom_density()+geom_vline(xintercept = diploid, color = "green", size = 1.5)+xlab("DNA content")
    print(dna_hist)
  }
  
  
  if (g1_good == "manual"){
    diploid_old <- diploid
    diploid <- readline(prompt = "What is the value of the 2N peak? ")
    if(diploid == ""){
      diploid <- diploid_old
    } else {
      diploid <- as.numeric(diploid)
    }
    dna_hist <- ggplot(eduNegCells, aes(log2_dna))+geom_density()+geom_vline(xintercept = diploid, color = "green", size = 1.5)+xlab("DNA content")
    print(dna_hist)
    peak_2n <- as.numeric(readline(prompt = "What should be the 2N - 4N trough value: "))

    dna_hist <- dna_hist+geom_vline(xintercept = peak_2n, color = "blue", size = 1.5)
    print(dna_hist)
    peak_4n <- as.numeric(readline(prompt = "What should be the 4N - >4N trough value: "))
    dna_hist <- dna_hist+geom_vline(xintercept = peak_4n, color = "red", size = 1.5)
    print(dna_hist)
    
    cells$dna_norm <<- (cells$log2_dna+1)-diploid
    cells[cells$log2_dna > peak_2n,]$ploidy <<- "4N"
    cells[cells$log2_dna > peak_4n,]$ploidy <<- ">4N"
    
  }else {
    cat(paste0("Using ", format(round(diploid, 2), nsmall = 2), " as the 2N peak value."))
    cat("\n")
    #This part bins each point's ploidy        
    cells$dna_norm <<- (cells$log2_dna+1)-diploid
    cells[cells$dna_norm > 1.6,]$ploidy <<- "4N"
    cells[cells$dna_norm > 2.6,]$ploidy <<- ">4N"
  }
  
  
  ckF <- readline(prompt = paste0("Is the filename ", thingee, "? (Y/n) "))
  if (ckF != "n"){
    write_csv(cells, thingee)
  } else{
    thingee <<-readline(prompt = 'What is the name of the file: ')
    write_csv(cells, thingee)
  }
  px <<- "dna_norm"
  py <<- "edu_norm"
  xname <<- "DNA content (Log 2)"
  yname <<- "EdU content (Log 10)"
  grapho(cells)
}


#-----------------------------------------------------------------
# So this is the compensating function in case (god forbid) you have spectral overlap. Make sure to run this having graphed the offending colors against eachother
# IMPORTANT - The skewed color must be graphed on the y axis
# IMPORTANT - This should be obsolete if purgo.R has been run correctly
compensate <- function(df = cells){
  # First it generates a datframe to hold point. This is kept open to eventually morph this into a shape fitting gate
  #gridIron <<- data.frame(X=c(1,1), Y=c(1,1))
  cat("Thank you for choosing coompensate(). To begin, lets assign a first point. Be sure the intended line follows the overlap line:")
  # Then you pick a point on the line of the overlapping population
  pointOne_x <- as.numeric(readline(prompt = "What is the x-value of the first point? "))
  pointOne_y <- as.numeric(readline(prompt = "What is the y-value of the first point? "))
  print(pointOne_x)
  print(pointOne_y)
  pointOne <- data.frame(X=pointOne_x, Y=pointOne_y)
  print(pointOne)
  cat("\n")
  print(qq+geom_point(data = pointOne, aes(x=X, y=Y, color = "red", size = 16)))
  cat(paste0("Great, I have added that point. Now lets move on to the second point"))

  # Then you pick a second point
  pointOne_x <- as.numeric(readline(prompt = "What is the x-value of the second point? "))
  pointOne_y <- as.numeric(readline(prompt = "What is the y-value of the second point? "))
  pointTwo <- data.frame(X=pointOne_x, Y=pointOne_y)
  pointOne <- rbind(pointOne, pointTwo)

  print(qq+geom_point(data = pointOne, aes(x=X, y=Y, color = "red", size = 16)))

  yORn <- readline(prompt = "Does this look right (Y/n)? ")
  if (yORn == "n"){
    cat("Sorry, lets try that again.")
    cat("\n")
    compensate()
  } else {
    cat("Excellent, generating the linear equation now:")
    cat("\n")
    rise <- pointOne[2,2]-pointOne[1,2]
    run <- pointOne[2,1]-pointOne[1,1]
    slop <- rise/run
    yint <- pointOne[1,2]-(pointOne[1,1]*slop)
    print(qq+geom_point(data = pointOne, aes(x=X, y=Y, color = "red", size = 16))+
            geom_abline(slope = slop, intercept = yint, size = 2, color = "blue"))
    cat(paste0("Slope has been determined to be: ", slop))
    cat("\n")
    cat(paste0("Y-intercept has been determined to be: ", yint))
    yORn <- readline(prompt = "Does this look right (Y/n)? ")
    if (yORn == "n"){
      cat("Sorry, lets try that again.")
      cat("\n")
      compensate()
    } else {
      cat("Applying compensation.")
      check <- df
      check[py] <- check[py]-(slop*check[px]+yint)
      cat("Regraphing:")
      grapho(check)
      yORn <- readline(prompt = "Does this look right (Y/n)? ")
      if (yORn == "n"){
        cat("Sorry, lets try that again.")
        cat("\n")
        compensate()
      } else {
        cells <<- check
        cat(paste0("Cells saved with compensated ", py))
        yORn <- readline(prompt = "Would you like to recalculate the means (y/n)? ")
        if (yORn == "y"){
          for (i in 1:ncol(cells)){
            if (grepl(py, names(cells)[i])){
              cells[paste0("I", names(cells)[i])] <<- cells[,i]*cells$Area
            }
          }
          for (i in 1:ncol(cells)){
            if (grepl(py, names(cells)[i])){
              cells[paste0("L", names(cells)[i])] <<- log(cells[,i]+1, 10)
            }
          }
          cat("Means recalculated. Thank you for your patience and patronage.")
          write_csv(cells, thingee)

        }
      }
    }
  }
}

#------------------------------------------------------

# This function take the nuclear locations and sizes and re-maps the original image
#   df is the dataframe, defaulted to 'cells'
reMap <- function(df = cells, 
                  #   highlight is a vector containing the Number(s) to be highlighted in remapped image
                  highlight = F, 
                  #   X and Y refer to the beginning of the X and Y position names in the dataframe, defaulted to the nuclear geometric position
                  Y = "X_ANC_", 
                  X = "Y_ANC_", 
                  #   S is the size of the ROI, defaulted to nucelar size
                  S = "Area_ANC_", 
                  #   Xn and Yn are the axis names and don't really need to be changed
                  Xn = "X position", 
                  Yn = "Y position",
                  #   colz is a vector of colors to use for the highlights variable. Be sure to have enough colors for the number of populations in the highlights variable
                  colz = cb_black) {
  #These set the variable names using regex to anchor tot he front of the variable name
  X <- names(df)[str_detect(names(df), paste0("^",X))]
  Y <- names(df)[str_detect(names(df), paste0("^",Y))]
  S <- names(df)[str_detect(names(df), paste0("^",S))]
  #This creates a ggplot object for graphing
  test <- ggplot(data = df, aes(y = -unlist(df[X]), x = unlist(df[Y]), size = unlist(df[S])))
  #If the highlight variable is a vector and is not F, its graphed here
  if (is.vector(highlight) & highlight != F) {
    #print(1)
    check <- df
    check$highlight <- check$Number %in% highlight
    qq <<- test+geom_point(aes(color = check$highlight))+
      ylab(Yn)+
      xlab(Xn)+
      scale_color_manual(values = colz)+
      theme_classic()+
      theme(axis.text = element_text(face='bold', size=16), axis.ticks = element_blank())+
      theme(plot.title = element_text(size=40, face="bold"))+
      theme(axis.title = element_text(size=30, face="bold"))
  } 
  # If highlight is a not a vector (such as variable called) and isn't F, then its graphed
  else if (length(highlight) != 1 & highlight != F){
    #print(2)
    qq <<- test+geom_point(aes(color = highlight))+
      ylab(Yn)+
      xlab(Xn)+
      theme_classic()+
      theme(axis.text = element_text(face='bold', size=16), axis.ticks = element_blank())+
      theme(plot.title = element_text(size=40, face="bold"))+
      theme(axis.title = element_text(size=30, face="bold"))+
      scale_color_manual(values = colz)
  }
  # Finally, if highlight IS FALSE, then its graphed without a color option
  else {
    #print(3)
    qq <<- test+geom_point()+
      ylab(Yn)+
      xlab(Xn)+
      theme_classic()+
      theme(axis.text = element_text(face='bold', size=16), axis.ticks = element_blank())+
      theme(plot.title = element_text(size=40, face="bold"))+
      theme(axis.title = element_text(size=30, face="bold"))
  }
  #The plot is graphed
  print(qq)
}

#-----------------------------------------------------------------
# Runs either a PCA of t-SNE
runCA <- function(df=cells,
                   plotType = "tsne",
                   saveFile = "figures/pcaPlot.png",
                   subz = F, 
                   groupz = T){  
  cellist <- df
  
  # Get rid of NAs and Infs
  cellist[is.na(cellist)] <- 0
  cellist[cellist == Inf] <- 0
  cellist[cellist == -Inf] <- 0 

  if (groupz ==  T){
    elly <- T
    print(names(cellist))
    grope <- readline(prompt = paste0("What should the ", plotType, " be grouped by from the 'cellist' dataframe: "))
    groupList <- unlist(df[grope])
  }
  # get rid of non-numerics
  num <- unlist(lapply(df, is.numeric))
  cellist_nums <- df[ , num]
  
  # Get only the nuclear columns & remove NA
  if (subz != F){
    cellist_nums <- cellist_nums[str_detect(names(cellist_nums), subz) | str_detect(names(cellist_nums), "Number")]
    cellist <- cellist[str_detect(names(cellist), subz) | str_detect(names(cellist), "Number")]
  }
  
  #Get rid of non-unique rows
  cellist <- unique(cellist)
  cellist_nums <- unique(cellist_nums)
  
  if (plotType != "pca" & plotType != "tsne"){
    plotType <- readline(prompt = "Are you trying to make a tsne or pca plot: ")
  }
  
  # Make the pca object
  if(plotType == "pca"){
    # Get rid of columns with no variance
    x <- foo(cellist_nums)
    cellist_nums <- cellist_nums[,-x]
    cellist_nums[is.na(cellist_nums)] <- 0
    
    cell.pca <- prcomp(cellist_nums, center = TRUE, scale = TRUE)
    percentage <- round(cell.pca$sdev / sum(cell.pca$sdev) * 100, 2)

    if (groupz == T){
      pca_out <- data.frame(x=cell.pca$x[,1], y = cell.pca$x[,2], groupBy = cellist[,grope])
      pp <<- ggplot(data = pca_out, aes(x = x, y = y))+
        geom_point(aes(color = groupBy), size = 3, alpha = 0.5)+
        geom_density_2d(aes(color = groupBy))+
        xlab(paste0("PC1 (", percentage[1], " %)"))+
        ylab(paste0("PC2 (", percentage[2], " %)"))+
        theme_classic2()
    } else {
      pca_out <- data.frame(x=cell.pca$x[,1], y = cell.pca$x[,2])
      pp <<- ggplot(data = pca_out, aes(x = x, y = y))+
        geom_point(size = 3, alpha = 0.5)+
        geom_density_2d()+
        xlab(paste0("PC1 (", percentage[1], " %)"))+
        ylab(paste0("PC2 (", percentage[2], " %)"))+
        theme_classic2()
    }
  } else{
    cellist_m <- as.matrix(cellist_nums)
    cellist_m[is.na(cellist_m)] <- 0
    tsne_out <- Rtsne::Rtsne(cellist_m, perplexity = 100)
    if (groupz == T){
      tsne_plot <- data.frame(x = tsne_out$Y[,1], y = tsne_out$Y[,2], groupBy = cellist[,grope])
      pp <<- ggplot(data = tsne_plot, aes(x=x, y=y))+geom_point(aes(color = groupBy), size = 3, alpha = 0.5)+geom_density_2d(aes(color = groupBy))+theme_classic2()
    } else{
      tsne_plot <- data.frame(x = tsne_out$Y[,1], y = tsne_out$Y[,2])
      pp <<- ggplot(data = tsne_plot, aes(x=x, y=y))+geom_point(size = 3, alpha = 0.5)+geom_density_2d()+theme_classic2()
    }
  }
  print(pp)

  #Report how much of the data is shown
  howManyLeft <- round((nrow(cellist_nums)*ncol(cellist_nums))/(nrow(df)*ncol(df))*100, digits = 1)
  cat(paste0("Approximately ", howManyLeft, "% of the data is represented after trimming."))
  cat("\n")
  
  
  ggsave(saveFile, dpi = 300, height = 4, width = 4)
}

foo <- function(dat) {
  out <- lapply(dat, function(x) length(unique(x)))
  want <- which(!out > 1)
  unlist(want)
}

#-----------------------------------------------------------------
#nearest neighbor function for determining the nearest cell from one population to any cell of another population
nearestN <- function(home=cells, target=cells, label="CellToCell"){
  cells[paste0("nearest_",label)] <<- "NotHome"
  cells[paste0("nearest_",label,"_distance")] <<- 0
  for (i in 1:nrow(cells)){
    if(cells$Number[i] %in% unique(home$Number)){
      pool <- subset(target, Number != cells$Number[i])
      pool$distance <- sqrt((pool$X_ANC_dna-cells$X_ANC_dna[i])**2+(pool$Y_ANC_dna-cells$Y_ANC_dna[i])**2)
      cells[paste0("nearest_",label)][i,] <<- subset(pool, distance == min(pool$distance))$Number[1]
      cells[paste0("nearest_",label, "_distance")][i,] <<- subset(pool, distance == min(pool$distance))$distance[1]
    }
  }
}

#-----------------------------------------------------------------
# This function uses the general principle of nearestN but assumes you're clustering cells based on a gated variable.
# The function sorts all other cells based on their distance to each cell, then counts and records the nearby cells that are within a selected category until it encounters a nearby cell that is within the edge category
# It then stops can records the size of the cluster centered around that cell, and the numberIDs of those cells within the cluster
# After characterizing the clusters of each cell, it then compares the various clusters
clusterCount <- function(df = cells, clustPop="tag", center="Positive", edge="Negative", clusterCutoff=0.5){
  cells["clusterSize"] <<- 0
  cells["clusterMates"] <<- "None"
  cells["clusterID"] <<- "None"
  clusterID <- 1
  home <- cells[cells[clustPop]==center,]
  for (i in 1:nrow(cells)){
    if(cells$Number[i] %in% unique(home$Number)){
      pool <- subset(cells, Number != cells$Number[i])
      pool$distance <- sqrt((pool$X_ANC_dna-cells$X_ANC_dna[i])**2+(pool$Y_ANC_dna-cells$Y_ANC_dna[i])**2)
      pool <- pool[order(pool$distance),]
      clusterSize <- 1
      clusterMates <- c(cells$Number[i])
      for (j in 1:nrow(pool)){
        if(pool[clustPop][j,] != edge){
          clusterSize <- clusterSize+1
          clusterMates <- append(clusterMates, as.character(pool[j,]$Number))
        } else {
          break
        }
      }
      clusterMates <- sort(clusterMates)
      cells["clusterSize"][i,] <<- clusterSize
      cells["clusterMates"][i,] <<- paste(unlist(clusterMates), collapse = "_")
    }
  }
  for (i in 1:nrow(cells)){
    if(cells$Number[i] %in% unique(home$Number & cells$clusterID[i] == "None")){
      pool <- subset(home, Number != cells$Number[i])
      pool$lScore <- 
      pCluster <- subset(pool, lScore > clusterCutoff)
      cells$clusterID[cells$Number %in% pCluster$Number,] <- clusterID
      clusterID <- clusterID + 1
    }
  }
}

#-----------------------------------------------------------------
#-----------------------------------------------------------------
#Color palettes:
Color_blind<- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cb_black <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#The cellCycle_colors palette are the hexdecimal codes for: 
#                        G1,         S,       G2/M,    mitosis, & premature mitosis respecitvely, as used in Justice et al, 2019
cellCycle_colors <- c("#D4D4D4", "#98C84C", "#23B8CC", "#F16B1A", "#E5001C")



#-----------------------------------------------------------------------------------------------------------
