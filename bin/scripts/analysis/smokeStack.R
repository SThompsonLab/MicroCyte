library(ggpubr)
library(viridis)

smokeStack <- function(df = cells, 
                       normBy = "zscore",
                       setX = F,
                       setY = F,
                       bins = 20,
                       binSize = F,
                       minBinNumber = 100,
                       lims = F,
                       facet = F,
                       facet_within =T,
                       fixedCoord = T,
                       mid_set=F,
                       saveFile = "figures/coloredBar.svg",
                       saveData = F,
                       xdim = 5,
                       ydim = 0.5,
                       viridisPallet = "inferno",
                       call_icellate = F, 
                       fuse_icellate = T,
                       number_icellate = 10){
  # Asks for a variable to bin on if one isn't given
  if (setX == F){
    print(names(df))
    xName <- readline(prompt = "Which variable should be mapped on the x-axis? ")
  } else {
    xName <- setX
  }
  
  # Makes sure the number of bins is a numeric value
  if(binSize != F){
    if (is.numeric(binSize)){
      print(paste("Using bin size at", binSize))
    } else{
      print("Non-numeric bin size detected. Ignoring bin size")
      binSize <- F
    }
  }
  
  # Creates the dataframe to store the bin variables
  magicSet <- data.frame("Bin" = c(1:bins))
  
  # Asks for categorical variables to map on the y-axis if a set isn't supplied
  if (setY == F){
    varNumb <- as.numeric(readline(prompt = paste0("How many variables should be mapped across ", xName, "? ")))
    varMat <- data.frame("varNumb" = c(1:varNumb), "varName" = "HOLD")
    
    for (i in 1:varNumb){
      magicSet[,1+i] <- readline(prompt = paste0("What is variable number ", i, "? "))
    }
  } else {
    for (i in 1:length(setY)){
      magicSet[,1+i] <- setY[i]
    }
  }
  for(i in 2:length(names(magicSet))){
    names(magicSet)[i] <- unique(magicSet[,i])
    magicSet[,i] <- 0
  }

  # Now we bin the original data set based on the supplied x-axis variable
  # First create a fungible df
  interim <- df
  # Set its bin variable
  interim$bin <- 0
  interim["Name"] <-  df[,xName]
  # Find the bin range and set the first bin minimum value to the global minimum value of the binning variable
  xNameBinRange <- (max(df[xName])-min(df[xName]))/bins
  binSet <- c(min(df[xName]))
  
  #If a lim double is supplied, the data is trimmed accordingly
  if(is.double(lims)){
    interim <- interim[interim$Name >= lims[1] & interim$Name <= lims[2],]
  }
  
  # Now find the minimum bin value for each bin
  for (i in 1:bins-1){
    binSet <- append(binSet, rev(binSet)[1]+xNameBinRange)
  }
  # Create the bin variable
  binID <- 0
  # Iterate through the ranges to bin cells into their respective bins
  for (i in 2:length(binSet)){
    binID <- binID+1
    if(nrow(interim[interim$Name >= binSet[i-1] & interim$Name <= binSet[i],])>0){
      interim[interim$Name >= binSet[i-1] & interim$Name <= binSet[i],]$bin <- binID
    }
  }
  
  if(facet != F){
    for (silliness in unique(unlist(interim[facet]))){
      interim_magic <- magicSet
      
      new_interim <- interim[interim[facet] == silliness,]

      for (i in names(magicSet)[2:length(magicSet)]){
        if (normBy == "zscore"){

          if(facet_within){
            new_interim$var <- (unlist(new_interim[i])-mean(unlist(interim[i])))/sd(unlist(interim[i]))
          }else {
            new_interim$var <- (unlist(new_interim[i])-mean(unlist(new_interim[i])))/sd(unlist(new_interim[i]))
          }
          legendName <- "Z-Score"
        } else if (normBy == "lq"){
          legendName <- "Lower quartile"
        } else if(normBy == "uq"){
          legendName <- "Upper quartile"
        } else if (normBy == "mean"){
          if(facet_within){
            new_interim$var <- 100*unlist(new_interim[i])/mean(unlist(interim[i]))
          } else {
            new_interim$var <- 100*unlist(new_interim[i])/mean(unlist(new_interim[i]))
          }
          legendName <- "% of Mean"
        } else if (normBy == "iqr"){
          legendName <- "Interquartile range"
        } else {
          legendName <- "Intensity"
        }
        for (j in 1:nrow(interim_magic)){
          if (normBy == "zscore"){
            interim_magic[i][j,] <- median(new_interim[new_interim$bin == j,]$var)
            legendName <- "Z-Score"
          } else if (normBy == "lq"){
            interim_magic[i][j,] <- quantile(unlist(new_interim[new_interim$bin == j,][i]))[2]
          } else if(normBy == "uq"){
            legendName <- "Upper quartile"
            interim_magic[i][j,] <- quantile(unlist(new_interim[new_interim$bin == j,][i]))[4]
          } else if (normBy == "mean"){
            interim_magic[i][j,] <- median(new_interim[new_interim$bin == j,]$var)
          } else if (normBy == "iqr"){
            legendName <- "Interquartile range"
            interim_magic[i][j,] <- quantile(unlist(new_interim[new_interim$bin==j,][i]))[4]-quantile(unlist(new_interim[new_interim$bin==j,][i]))[2]
          } else {
            interim_magic[i][j,] <- median(unlist(new_interim[new_interim$bin == j,][i]))
            legendName <- "Intensity"
          }
          if(nrow(new_interim[new_interim$bin == j,]) < minBinNumber){
            interim_magic[i][j,] <- NA
          }
          if(call_icellate){
            icellate(targetCells = new_interim[new_interim$bin == j,],
                     folderName = paste0(facet, "-", silliness, "_Bin_", j),
                     fuse = fuse_icellate,
                     samplingNumber = number_icellate,
                     lineAnalyses = F,
                     heatIntensity = T, 
                     saveData = "icellates/bin_data.csv")
          }
        }
      }
      if(!exists("dumbness")){
        dumbness <- interim_magic
        dumbness[facet] <- silliness
      } else {
        interim_magic[facet] <- silliness
        dumbness <- rbind(dumbness, interim_magic)
      }
    }
    magicSet <- dumbness
  } else {
    for (i in names(magicSet)[2:length(magicSet)]){
      if (normBy == "zscore"){
        interim$var <- (unlist(interim[i])-mean(unlist(interim[i])))/sd(unlist(interim[i]))
        legendName <- "Z-Score"
      } else if (normBy == "lq"){
        legendName <- "Lower quartile"
      } else if(normBy == "uq"){
        legendName <- "Upper quartile"
      } else if (normBy == "mean"){
        interim$var <- 100*unlist(interim[i])/mean(unlist(interim[i]))
        legendName <- "% of Mean"
      } else if (normBy == "iqr"){
        legendName <- "Interquartile range"
      } else {
        legendName <- "Intensity"
      }
      for (j in 1:nrow(magicSet)){
        if (normBy == "zscore"){
          magicSet[i][j,] <- median(interim[interim$bin == j,]$var)
          legendName <- "Z-Score"
        } else if (normBy == "lq"){
          magicSet[i][j,] <- quantile(unlist(interim[interim$bin == j,][i]))[2]
        } else if(normBy == "uq"){
          legendName <- "Upper quartile"
          magicSet[i][j,] <- quantile(unlist(interim[interim$bin == j,][i]))[4]
        } else if (normBy == "mean"){
          magicSet[i][j,] <- median(interim[interim$bin == j,]$var)
        } else if (normBy == "iqr"){
          legendName <- "Interquartile range"
          magicSet[i][j,] <- quantile(unlist(interim[interim$bin==j,][i]))[4]-quantile(unlist(interim[interim$bin==j,][i]))[2]
        } else {
          magicSet[i][j,] <- median(unlist(interim[interim$bin == j,][i]))
          legendName <- "Intensity"
        }
        if(nrow(interim[interim$bin == j,]) < minBinNumber){
          magicSet[i][j,] <- NA
        }
        
        if(call_icellate){
          print(nrow(interim[interim$bin == j,]))
          icellate(targetCells = interim[interim$bin == j,],
                   folderName = paste0("Bin_", j),
                   fuse = fuse_icellate,
                   samplingNumber = number_icellate,
                   lineAnalyses = F, 
                   heatIntensity = T, 
                   saveData = "icellates/bin_data.csv")
        }
      }
    }
  }
  if (facet != F){
    magicSet <- magicSet %>% pivot_longer(cols = c(2:(ncol(magicSet)-1)), names_to = "Target", values_to = "Value")
    magicSet$Target <- paste0(facet, ": ", unlist(magicSet[facet]), " - ", magicSet$Target)
  }else {
    magicSet <- magicSet %>% pivot_longer(cols = c(2:ncol(magicSet)), names_to = "Target", values_to = "Value")
  }

  cBin <- c(1:bins)
  cBin <- cBin[c(T, F)]
  
  binSet <- binSet[2:length(binSet)]
  binSet <- round(binSet - (xNameBinRange/2),1)
  binSet <- binSet[c(T, F)]
  
  draft <- ggplot(data = magicSet, aes(x = Bin, y = Target, fill = Value))+
    geom_tile()+
    theme_classic2()+
    scale_fill_viridis(option = viridisPallet)+
    xlab(xName)+
    ylab("")+
    guides(fill=guide_legend(title = legendName))+
    theme(legend.position = "bottom")+
    scale_x_continuous(breaks = cBin, labels = binSet)
  if(fixedCoord==T){
    draft  <- draft+coord_fixed()
  }
  if (mid_set != F & is.numeric(mid_set)){
    draft <- draft+scale_fill_gradient2(low = "blue", midpoint = mid_set, mid = "black", high = "red")
  }
  
  print(draft)
  if(saveFile != F & is.character(saveFile)){
    ggsave(saveFile)
  }
  if(saveData != F & is.character(saveData)){
    write.csv(magicSet, saveData, row.names = F)
  }
}

