#!/usr/bin/Rscript

#The imaGen() script takes any number of csv files that are within your working directory and combines them. Take caution, however, as this script is designed for a particular imagej csv outpt
#It is important that the jetData uses the Set Measurements with the following settings active:
#   Area    Standard Deviation    Min & Max gray value    Center of Mass    Mean gray value   Perimenter    Display label
# Additionally, for the best results, decimal places should be set to 9

# If you have used the YggData macro, this will all be automatically set for you

#This is the joinR function used to assign nuclear number to each ROI that may or may not be localized to the nucleus
joinR <- function(y, x = "dna.csv"){
  # Files are opened
  tagman <- read.csv(x)
  fraction <- read.csv(y)
  # The number column is named, if not present
  if (!"Number" %in% names(tagman)){
    names(tagman)[1] <- "Number"
    write.csv(tagman, file = x, row.names = FALSE)
  }
  if (!"Number" %in% names(fraction)){
    names(fraction)[1] <- "Number"
    write.csv(fraction, file = y, row.names = FALSE)
  }

  # Now we skip fractions with empty lists
  if(nrow(fraction) > 0){
    #Default cellid is given
    fraction$cell <- "unknown"
    #Default distance is given
    tagman$distance <- 0
    fraction$NucDist <- 0
    
    # The closest nucleus is calculated for each ROI
    for (i in 1:nrow(fraction)){
      tagman$distance <- sqrt((tagman$X-fraction$X[i])^2+(tagman$Y-fraction$Y[i])^2)
      closest <- min(tagman$distance)
      fraction$NucDist[i] <- closest
      fraction$cell[i] <- unique(subset(tagman, distance == closest)$Number)
      #If two or more nuclei are picked, it lets you know
      if (length(closest) > 1){
        print("PING! More than one nucleus at minimum distance. Something is wrong (probably)")
      }
    }
    write.csv(fraction, file = y, row.names = FALSE)
    print(paste0("Saved as ", y))
  } else {
    print(paste0("No data found for file: ", y))
  }
}

imaGen <- function(directory="./", 
                   colorz = T, 
                   label_ids = T,
                   peri = T,
                   wc = T){
  # The directory is set and a list of CSV's is generated
  directoryN <- paste0(directory, "/Anchor_extraction/")
  directoryC <- paste0(directory, "/NonAnchor_extraction/")
  directoryP <- paste0(directory, "/PeriAnchor_extraction/")
  
  # Goes to the anchor directory
  setwd(directoryN)
  # Generates the list of data files
  filez <- list.files(pattern = ".csv")
  filez <- filez[!grepl("all.csv", filez)]
  
  # If there is an explicit anchor file, it uses that, otherwise it expects dna or asks the user to define it
  if ("anchor.csv" %in% filez){
    anchorName <- "anchor"
  } else if ("dna.csv" %in% filez){
    anchorName <- "dna"
  } else {
    print(filez)
    anchorName <- readline(prompt = "Anchor file not detected. Which file should be the anchor?")
  }
  
  # The first file is opened to serve as a template
  # it doesn't matter which file, since they all share the same order based on the ancher-dependent extraction
  cells <- read.csv(filez[1])
  
  # If the csv name is not simply the target name (colorz == F), 
  # then you will be asked to name a target for each image
  
  # REMEMBER: One of the targets must be 'dna'
  if(colorz == F){
    cat("REMEMBER: One of these color's must be match the given anchor name.")
    cat("\n")
    cat(paste0("This file is ", filez[1]))
    cat("\n")
    colo <- readline(prompt= paste0("What color is in ", filez[1], ": "))
  } else{
    # Otherwise, the target is taken from the csv name
    colo <- substr(filez[1], 1, nchar(filez[1])-4)
  }
  
  # The area is artificially increased to make small numbers bigger
  cells$Area <- cells$Area*100
  
  # The template's columns are appended with the anchor marker and the 'color
  names(cells) <- paste0(names(cells), "_ANC_", colo)
  
  # The first column with designated as the cell number which is dictated by the nucleus 
  # and is consistent across all downstream applications
  names(cells)[1] <- "Number"
  
  # The label is removed as this isn't required after 'colo' labeling
  cells <- cells[,-2]
  
  # The rest of the csv's in the Nuclear directory are opened and cbinded to the first set
  for (i in filez[2:length(filez)]){
    # Already bound csv's (_all tagged) are ignored
    if (!grepl("_all.csv", i)){
      interim <- read.csv(i)
      if(colorz == F){
        cat("REMEMBER: One of these color's must match the given anchor name")
        cat("/n")
        cat(paste0("This file is ", i))
        cat("/n")
        colo <- readline(prompt= paste0("What color is in ", i, ": "))
      } else{
        colo <- substr(i, 1, nchar(i)-4)
      }
      names(interim) <- paste0(names(interim), "_ANC_", colo)
      cells <- cbind(cells, interim[,4:8], interim[,11:12], interim[,15:16])
    }
  }
  # The anchor dependent file is saved
  filnam <- yL
  xoo <- paste0(filnam,"_ANC_all.csv")
  write.csv(cells, file = paste0(filnam,"_ANC_all.csv"), row.names = FALSE)
  setwd("../")
  
#Now that the nuclear data is collected, the ROI data will be parsed and added
  # First, you navigate to the NonAnchor file folder
  if(wc==T){
    cat("Parsing anchor-indepedent data...\n")
    setwd(directoryC)
    
  # Then you get the list of files
    filez <- list.files()
    filez <- filez[!grepl("all.csv", filez)]
    
    # Every ROI is designated an anchor number and the file is saved
    for (j in filez){
      joinR(j, paste0("../Anchor_extraction/", anchorName, ".csv"))
    }
    
    # The anchor-extraction file is opened
    cells <- read.csv(paste0("../Anchor_extraction/", xoo))

    #Each ROI file is summarized and added to the anchor data
    for (i in filez){
      # First let smake sure they have labels (they do)
      if (!grepl("_all.csv", i)){
        if(colorz == F){
          cat("REMEMBER: One of these color's must match the given anchor'")
          cat("\n")
          cat(paste0("This file is ", i))
          cat("\n")
          colo <- readline(prompt= paste0("What color is in ", i, ": "))
        } else{
          colo <- substr(i, 1, nchar(i)-4)
        }
      }
      
      # Temporary dataset is made for easy column labeling and given summary columns
      interim <- cells
      interim$ROI_Num <- 0
      interim$ROI_Area <- 0
      interim$ROI_Mean <- 0
      interim$ROI_Stdev <- 0
      interim$ROI_Mode <- 0
      interim$ROI_Perimeter <- 0
      interim$ROI_IntDen <- 0
      interim$ROI_IntTotal <- 0
      interim$ROI_NucDist <- 0
      
      # The ROI files is opened as 'fraction'
      if (!grepl("_all.csv", i)){
        fraction <- read.csv(i)
        
        # if there are ROI values, it continues
        if (nrow(fraction > 0)){
          #Each nucleus number is analyzed for ROI's designated to it
          for (j in 1:nrow(interim)){
            
            if(interim$Number[j] %in% unique(fraction$cell)){
              # The ROIs that are assigned to a particular nucleus number are subsetted
              hitz <- subset(fraction, cell == interim$Number[j])
              if (colo != anchorName){
                cat(paste0("found ", nrow(hitz), " ", colo, " ROIS for anchor number ", interim$Number[j], "\n"))
              }
              
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
              interim$ROI_IntDen[j] <- mean(hitz$RawIntDen)
              
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
    }
    write.csv(cells, file = paste0(filnam, "_WC_all.csv"), row.names = F)
    
    # Now begins generating a final ROI file for funnsies
    interim <- read.csv(filez[1])
    colo <- substr(filez[1], 1, nchar(i)-2)
    interim$roi <- colo

    if(length(filez) > 1){
      for (i in filez[2:length(filez)]){
        if(!grepl("_all.csv", i)){
          x <- read.csv(i)
          if (nrow(x > 0)){
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
      }
    }

    #Now things get saved
    write.csv(interim, file = paste0("../../",filnam, "_ROI_all.csv"), row.names = F)
    write.csv(cells, file = paste0("../../", filnam, "_WN_all.csv"), row.names = F)
    setwd("../")
  }
  #This ends wc collection------------------------------------
  
  #start of peri collection
  if(peri==T){
    setwd(directoryP)
    # Then you get the list of files
    filez <- list.files()
    filez <- filez[!grepl("_all.csv", filez)]
    
    # If the ROI-independent files have been called, that's what we want to pull
    if (wc){
      nucleoid <- paste0("../../", filnam, "_WN_all.csv")
    } else {
      nucleoid <- paste0("../Anchor_extraction/", filnam, "_ANC_all.csv")
    }
    print(nucleoid)

    # Every ROI is designated a nucleus number
    for (j in filez){
      joinR(j, paste0("../Anchor_extraction/", anchorName, ".csv"))
    }
    print("Peris joined")
    
    # The tagger file is opened
    cells <- read.csv(nucleoid)

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
          interim$PERI_Area[j] <- 0
          interim$PERI_Mean[j] <- 0
          interim$PERI_Stdev[j] <- 0
          interim$PERI_Mode[j] <- 0
          interim$PERI_Perimeter[j] <- 0
          interim$PERI_IntDen[j] <- 0
          if(interim$Number[j] %in% unique(fraction$cell)){
            # The PERI's that are assigned to a particular nucleus number are subsetted
            hitz <- subset(fraction, cell == interim$Number[j])
            
            # Individual ROI nuclear distances are calculated
            #The number of ROIS assigned to that cell
            
            #The average average of all the ROIs
            interim$PERI_Area[j] <- mean(hitz$Area)*100
            
            #The average Mean intensity of all the ROIs
            interim$PERI_Mean[j] <- mean(hitz$Mean)
            
            #The average standard deviation of all the ROIs
            interim$PERI_Stdev[j] <- mean(hitz$StdDev)
            
            #The average mode of pixel intensities of all the ROIs
            interim$PERI_Mode[j] <- mean(hitz$Mode)
            
            #The average perimeter of all the ROIs
            interim$PERI_Perimeter[j] <- mean(hitz$Perim.)
            
            #The average integrated mean/density of the ROIs
            interim$PERI_IntDen[j] <- mean(hitz$RawIntDen)
          } 
        }
        # A target tag is added for later identification
        names(interim) <- paste0(names(interim), "_", colo)
        # The new variables are added to their cell ID
        cells <- cbind(cells, interim[,(ncol(cells)+1):ncol(interim)])
      }
    }
    # Now begins generating a final ROI file for funnies
    interim <- read.csv(filez[1])
    colo <- substr(filez[1], 1, nchar(i)-2)
    interim$peri <- colo
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
        x$peri <- colo
        interim <- rbind(interim, x)
      }
    }
    
    #Now things get saved
    write.csv(interim, file = paste0("../../",filnam, "_PERI_all.csv"), row.names = F)
    write.csv(cells, file = paste0(filnam, "_PN_all.csv"), row.names = F)
    write.csv(cells, file = paste0("../../", filnam, "_WN_all.csv"), row.names = F)
    setwd("../")
  }
  #-------------------------------------------------------------------------
  # This ends the PERI_all collection

  # #Now that we have a nuclear dataset and a WC dataset, we might as well put it all together...
  # nucka <- paste0("Anchor_extraction/", filnam, "_ANC_all.csv")
  # nuke <- read.csv(nucka)
  # if(wc==T & peri==T){
  #   print("Appending peri- and non-anchor extractions...")
  #   wucka <- paste0("NonAnchor_extraction/", filnam, "_WC_all.csv")
  #   wuke <- read.csv(wucka)
  #   pucka <- paste0("PeriAnchor_extraction/", filnam, "_PN_all.csv")
  #   puke <- read.csv(pucka)
  #   cells <- cbind(nuke, puke[,22:ncol(puke)], wuke[,22:ncol(wuke)])
  # }
  # if(wc == T){
  #   print("Appending non-anchor extractions...")
  #   wucka <- paste0("NonAnchor_extraction/", filnam, "_WC_all.csv")
  #   wuke <- read.csv(wucka)
  #   cells <- cbind(nuke, wuke[,22:ncol(wuke)])
  # }
  # if(peri==T){
  #   print("Appending peri-anchor extractions...")
  #   pucka <- paste0("PeriAnchor_extraction/", filnam, "_PN_all.csv")
  #   puke <- read.csv(pucka)
  #   cells <- cbind(nuke, puke[,22:ncol(puke)])
  # }
  if (wc | peri){
    cells <- read.csv(paste0("../", filnam, "_WN_all.csv"))
  } else {
    cells <- read.csv(paste0("Anchor_extraction/", filnam,"_ANC_all.csv"))
  }
  
  if(peri==T){
    thinkTank <- names(cells)[grepl("PERI_Area", names(cells))]
    for (z in thinkTank){
      zim <- strsplit(z, "Area")[[1]][2]
      tic <- names(cells)[grepl("Area_ANC", names(cells))][1]
      cells[paste0("PERI_SubArea", zim)] <- cells[paste0("PERI_Area", zim)]-cells[tic]
      tic <- paste0("IntDen_ANC", zim)
      cells[paste0("PERI_SubIntDen", zim)] <- cells[paste0("PERI_IntDen", zim)]-cells[tic]
      cells[paste0("PERI_SubMean", zim)] <- (100*cells[paste0("PERI_SubIntDen", zim)])/cells[paste0("PERI_SubArea", zim)]
    }
  }
  
  #And we'll add some metadata and rename the NumberIDs since this is no longer assumed to be a sitched image
  cells$image <- yL
  cells$Number <- paste0(yL,"_",cells$Number)
  schema <- read.csv("../../../../schema.csv")
  if (xL %in% unique(schema$notebook_id)){
    hit <- subset(schema, notebook_id == xL)
  } else {
    hit <- subset(schema, name_id == xL)
  }
  if (nrow(hit)>1){
    print("More than one scheme match found.")
    print("The first match will be used but it is recommended that your ammend the schema file to avoid other issues")
    hit <- hit[1,]
  }
  for (ab in 7:ncol(hit)){
    cells[names(hit)[ab]] <- hit[ab]
    if (peri==T){
      #print("adding peri meta")
      peris <- read.csv(paste0("../",filnam, "_PERI_all.csv"))
      #print(head(peris))
      peris[names(hit)[ab]] <- hit[ab]
      #print(head(peris))
      write.csv(peris, file = paste0("../",filnam, "_PERI_all.csv"), row.names = F)
    }
    if(wc==T){
      rois <- read.csv(paste0("../",filnam, "_ROI_all.csv"))
      rois[names(hit)[ab]] <- hit[ab]
      write.csv(rois, file = paste0("../",filnam, "_ROI_all.csv"), row.names = F)
    }
  }
  #And finally save the whole thing
  write.csv(cells, file = paste0("../",filnam, "_WN_all.csv"), row.names = F)
}

setwd("files")

periGo <- FALSE
wcGo <- FALSE

xList <- list.files()
for (xL in xList){
  if (!grepl("ijm", xL)){
    setwd(xL)
    print(paste0("Running ImaGen on the ", xL," folder:"))
    yList <- list.files()
    for(yL in yList){
      #print(yL)
      if(!grepl("_all.csv", yL) & yL != "Thumbs.db"){
        setwd(paste0(yL, "/PNGS/"))
        checkList <- list.files()
        if("PeriAnchor_extraction" %in% checkList){
          periGo <- TRUE
        }
        if("NonAnchor_extraction" %in% checkList){
          wcGo <- TRUE
        }
        if (length(list.files(path = "Anchor_extraction/", pattern = "ANC_all.csv")) == 0){
          imaGen(peri = periGo,
                 wc = wcGo)
          print(paste0("Combining data from image ", yL))
        } else {
          imaGen(peri = periGo,
                 wc = wcGo)
          print(paste0("Combined data already detected for: ", yL))
        }
        setwd("../../")
      }
    }
    setwd("../")
  }
}

setwd("../")