#!/usr/bin/Rscript

#The imaGen() script takes any number of csv files that are within your working directory and combines them. Take caution, however, as this script is designed for a particular imagej csv outpt
#It is important that the jetData uses the Set Measurements with the following settings active:
#   Area    Standard Deviation    Min & Max gray value    Center of Mass    Mean gray value   Perimenter    Display label
# Additionally, for the best results, decimal places should be set to 9

# If you have used the YggData macro, this will all be automatically set for you

#This is the joinR function used to assign nuclear number to each ROI that may or may not be localized to the nucleus
joinR <- function(y, x = "dna.csv"){
  # Files are opened
  tagger <- read.csv(x)
  fraction <- read.csv(y)
  # The number column is named, if not present
  if (!"Number" %in% names(tagger)){
    names(tagger)[1] <- "Number"
    write.csv(tagger, file = x, row.names = FALSE)
  }
  if (!"Number" %in% names(fraction)){
    names(fraction)[1] <- "Number"
    write.csv(fraction, file = y, row.names = FALSE)
  }
  #Default cellid is given
  fraction$cell <- "unknown"
  #Default distance is given
  tagger$distance <- 0
  fraction$NucDist <- 0
  # The closest nucleus is calculated for each ROI
  for (i in 1:nrow(fraction)){
    tagger$distance <- sqrt((tagger$X-fraction$X[i])^2+(tagger$Y-fraction$Y[i])^2)
    closest <- min(tagger$distance)
    fraction$NucDist[i] <- closest
    fraction$cell[i] <- unique(subset(tagger, distance == closest)$Number)
    #If two or more nuclei are picked, it lets you know
    if (length(closest) > 1){
      print("PING! More than one nucleus at minimum distance. Something is wrong (probably)")
    }
  }
  write.csv(fraction, file = y, row.names = FALSE)
}

imaGen <- function(directory="./", colorz = T, label_ids = T){
  # The directory is set and a list of CSV's is generated
  directoryN <- paste0(directory, "/Nuclear/")
  directoryC <- "WholeCell/"
  setwd(directoryN)
  filez <- list.files(pattern = ".csv")
  # The first file is opened to serve as a template
  cells <- read.table(filez[1], sep = ",", header = TRUE)
  
  # If the csv name is not simply the target name (colorz == F), then you will be asked to name a target for each image
  # REMEMBER: One of the targets must be 'dna'
  if(colorz == F){
    cat("REMEMBER: One of these color's must be labeled 'dna'")
    cat("\n")
    cat(paste0("This file is ", filez[1]))
    cat("\n")
    colo <- readline(prompt= paste0("What color is in ", filez[1], ": "))
  } else{
    # Otherwise, the target is taken from the csv name
    colo <- substr(filez[1], 1, nchar(filez[1])-4)
  }
  
  cells$Area <- cells$Area*100
  names(cells) <- paste0(names(cells), "_NUC_", colo)
  # The first column with designated as the cell number which is dictated by the nucleus and is consistent across all downstream applications
  names(cells)[1] <- "Number"
  # The label is removed as this isn't required after 'colo' labeling
  cells<-cells[,-2]
  # Default area is quite small, so its boosted to make the numbers more "real"
  
  # The rest of the csv's in the Nuclear directory are opened and cbinded to the first set
  for (i in filez[2:length(filez)]){
    # Already bound csv's (_all tagged) are ignored
    if (!grepl("_all.csv", i)){
      interim <- read.table(i, sep = ",", header = TRUE)
      if(colorz == F){
        cat("REMEMBER: One of these color's must be labeled 'dna'")
        cat("/n")
        cat(paste0("This file is ", i))
        cat("/n")
        colo <- readline(prompt= paste0("What color is in ", i, ": "))
      } else{
        colo <- substr(i, 1, nchar(i)-4)
      }
      names(interim) <- paste0(names(interim), "_NUC_", colo)
      cells <- cbind(cells, interim[,4:8], interim[,11:12], interim[,15:16])
    }
  }
  filnam <- yL
  write.csv(cells, file = paste0(filnam,"_NUC_all.csv"), row.names = FALSE)
  xoo <- list.files(pattern = "dna")
  setwd("../")
  
#Now that the nuclear data is collected, the ROI data will be parsed and added
  
  # First, you navigate to the wholeCell folder
  setwd(directoryC)
  # Then you get the list of files
  filez <- list.files()
  # Then you designate which tagger file to use. If colorz == T, the script assumes your file is called dna.csv and uses the same csv as found in the Nuclear folder
  if(colorz==F){
    pathy <- getwd()
    cat(paste0("Current directory is ", pathy))
    cat("\n")
    tagger <- readline(prompt = "What is the path to the dna file: ")
  } else{
    tagger <- paste0("../Nuclear/", xoo)
  }
  # Every ROI is designated a nucleus number
  for (j in filez){
    if (!grepl("_all.csv", j)){
      joinR(j, xoo)
    }
  }
  # The tagger file is opened
  cells <- read.csv(tagger)
  # Default area is quite small, so its boosted to make the numbers more "real"
  cells$Area <- cells$Area*100
  # The dna values are labeled as such
  names(cells) <- paste0(names(cells), "_WC_dna")
  # The first column with designated as the cell number which is dictated by the nucleus and is consistent across all downstream applications
  names(cells)[1] <- "Number"
  # The label is removed as this isn't required after 'colo' labeling
  cells<-cells[,-2]
  #Each ROI file is summarized and added to the nucleus data
  for (i in filez){
    if (!grepl("_all.csv", i)){
      if(colorz == F){
        cat("REMEMBER: One of these color's must be labeled 'dna'")
        cat("\n")
        cat(paste0("This file is ", i))
        cat("\n")
        colo <- readline(prompt= paste0("What color is in ", i, ": "))
      } else{
        colo <- substr(i, 1, nchar(i)-4)
      }
    }
    # Temporary dataset is made for easy column labeling
    interim <- cells
    # The ROI files is opened
    if (!grepl("_all.csv", i)){
      fraction <- read.csv(i)
    
#      write.csv(fraction, file = i, row.names = F)
    #Each nucleus number is analyzed for ROI's designated to it
      for (j in 1:nrow(interim)){
        interim$ROI_Num[j] <- 0
        interim$ROI_Area[j] <- 0
        interim$ROI_Mean[j] <- 0
        interim$ROI_Stdev[j] <- 0
        interim$ROI_Mode[j] <- 0
        interim$ROI_Perimeter[j] <- 0
        interim$ROI_IntDens[j] <- 0
        interim$ROI_IntTotal[j] <- 0
        interim$ROI_NucDist[j] <- 0
        if(interim$Number[j] %in% unique(fraction$cell)){
        # The ROIs that are assigned to a particular nucleus number are subsetted
          hitz <- subset(fraction, cell == interim$Number[j])
        
        # Individual ROI nuclear distances are calculated
        #The number of ROIS assigned to that cell
          interim$ROI_Num[j] <- nrow(hitz)
        
        #The average average of all the ROIs
          interim$ROI_Area[j] <- mean(hitz$Area)*100
        
        #The average Mean intensity of all the ROIs
          interim$ROI_Mean[j] <- mean(hitz$Mean)
        
        #The average standard deviation of all the ROIs
          interim$ROI_Stdev[j] <- mean(hitz$StdDev)
        
        #The average mode of pixel intensities of all the ROIs
          interim$ROI_Mode[j] <- mean(hitz$Mode)
        
        #The average perimeter of all the ROIs
          interim$ROI_Perimeter[j] <- mean(hitz$Perim.)
        
        #The average integrated mean/density of the ROIs
          interim$ROI_IntDens[j] <- mean(hitz$RawIntDen)
        
        #The total of the integrated means/densities of all the ROIs
          interim$ROI_IntTotal[j] <- sum(hitz$RawIntDen)
        
        #The average distance from the nucleus of all the ROIs
          interim$ROI_NucDist[j] <- mean(hitz$NucDist)
        } 
      }
    # A target tag is added for later identification
      names(interim) <- paste0(names(interim), "_", colo)
    # The new variables are added to their cell ID
      cells <- cbind(cells, interim[,(ncol(cells)+1):ncol(interim)])
    }
  }
  
  # This ends the WC_all collection
  # Now begins generating a final ROI file for funnies
  interim <- read.csv(filez[1])
  colo <- substr(filez[1], 1, nchar(i)-2)
  interim$roi <- colo
  for (i in filez[2:length(filez)]){
    if(!grepl("_all.csv", i)){
      x <- read.csv(i)
      colo <- substr(i, 1, nchar(i)-4)
      for (j in names(interim)){
        if (!j %in% names(x)){
          x[j] <- "NA"
        }
      }
      for (j in names(x)){
        if (!j %in% names(interim)){
          interim[j] <- "NA"
        }
      }
      x$roi <- colo
      interim <- rbind(interim, x)
    }
  }
  
  #Now things get saved
  write.csv(interim, file = paste0("../../",filnam, "_ROI_all.csv"), row.names = F)
  write.csv(cells, file = paste0(filnam, "_WC_all.csv"), row.names = F)
  setwd("../")
  
  #Now that we have a nuclear dataset and a WC dataset, we might as well put it all together...
  nucka <- paste0("Nuclear/", filnam, "_NUC_all.csv")
  nuke <- read.csv(nucka)
  wucka <- paste0("WholeCell/", filnam, "_WC_all.csv")
  wuke <- read.csv(wucka)
  cells <- cbind(nuke, wuke[,22:ncol(wuke)])
  
  #And we'll add some metadata and rename the NumberIDs since this is no longer assumed to be a sitched image
  cells$image <- yL
  cells$Number <- paste0(yL,"_",cells$Number)
  schema <- read.csv("../../../../schema.csv")
  if (xL %in% unique(schema$notebook_id)){
    hit <- subset(schema, notebook_id == xL)
  } else {
    hit <- subset(schema, name_id == xL)
  }
  for (ab in 6:ncol(hit)){
    cells[names(hit)[ab]] <- hit[ab]
  }
  #And finally save the whole thing
  write.csv(cells, file = paste0("../",filnam, "_WN_all.csv"), row.names = F)
}

#setwd("../../files/")
xList <- list.files()
for (xL in xList){
  if (!grepl("ijm", xL)){
    setwd(xL)
    yList <- list.files()
    for(yL in yList){
      setwd(paste0(yL, "/PNGS/"))
      imaGen()
      setwd("../../")
    }
    setwd("../")
  }
}
