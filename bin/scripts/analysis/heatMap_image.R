suppressPackageStartupMessages(library(imager))
suppressPackageStartupMessages(library(magick))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))

heatMap_image <- function(image_location,
                          heatColor = "inferno",
                          intensity_binning = F
                          ){
  if(!grepl(pattern = ".png", x = image_location) & grepl(pattern = ".tif", x = image_location)){
    print("TIF detected. Creating PNG.")
    tick <- image_read(image_location)
    image_location <- str_replace(image_location, "tif", "png")
    image_write(tick, path = image_location, format = "png")
  }
  
  if (intensity_binning){
    oldImage <- load.image(image_location)
    # plot(oldImage)
    # pause <- readline(prompt = "Pausing...1")
    thick <- as.data.frame(oldImage)
    if (4 %in% unique(thick$cc)){
      print("Alpha channel detected. Removing and re-saving.")
      thick <- thick[thick$cc != 4,]
      newPlot <- suppressWarnings(as.cimg(thick))
      newPlot <- cimg2magick(newPlot, rotate = T)
      newPlot <- image_flop(newPlot)
      image_write(newPlot, path = image_location, format = "png")
      oldImage <- load.image(image_location)
      thick <- as.data.frame(oldImage)
      print("3-Channel PNG saved.")
    }
    check_gray <- as.data.frame(grayscale(oldImage))
    
    check_gray$R <- 0
    check_gray$G <- 0
    check_gray$B <- 0
    # 
    # print(head(check_gray))
    # pause <- readline(prompt = "Pausing...2")
    
    #check_gray$value <- log(check_gray$value, 10)
    print("Creating binned pallet hexes.")
    fake_plot <- ggplot(data = check_gray, aes(fill=value, x=x, y=y))+
      geom_point()+
      scale_fill_viridis_b(option = heatColor)
    fake_plot <- ggplot_build(fake_plot)
    fake_plot <- as.data.frame(fake_plot$data)
    # print(head(fake_plot))
    # pause <- readline(prompt = "Pausing...3")
    check_gray$hex <- fake_plot$fill
    plaything <<- check_gray
    #print(head(check_gray))
    for (hexcode in unique(check_gray$hex)){
      #print(hexcode)
      ick <- check_gray[check_gray$hex == hexcode,]
      #print(head(ick))
      if (!exists("hexData")){
        hexData <- data.frame("Image" = image_location,
                              "HexCode"= hexcode,
                              "Mean_X" = median(check_gray[check_gray$hex == hexcode,]$x),
                              "Mean_Y" = median(check_gray[check_gray$hex == hexcode,]$y),
                              "Mean" = mean(check_gray[check_gray$hex == hexcode,]$value),
                              "Median" = median(check_gray[check_gray$hex == hexcode,]$value),
                              "Range" = max(check_gray[check_gray$hex == hexcode,]$value)-min(check_gray[check_gray$hex == hexcode,]$value),
                              "Total" = sum(check_gray[check_gray$hex == hexcode,]$value),
                              "Area_px" = nrow(check_gray[check_gray$hex == hexcode,])
        )
        #print(hexData)
      } else {
        interimHex <- data.frame("Image" = image_location,
                                 "HexCode"= hexcode,
                                 "Mean_X" = median(check_gray[check_gray$hex == hexcode,]$x),
                                 "Mean_Y" = median(check_gray[check_gray$hex == hexcode,]$y),
                                 "Mean" = mean(check_gray[check_gray$hex == hexcode,]$value),
                                 "Median" = median(check_gray[check_gray$hex == hexcode,]$value),
                                 "Range" = max(check_gray[check_gray$hex == hexcode,]$value)-min(check_gray[check_gray$hex == hexcode,]$value),
                                 "Total" = sum(check_gray[check_gray$hex == hexcode,]$value),
                                 "Area_px" = nrow(check_gray[check_gray$hex == hexcode,])
        )
        #print(interimHex)
        hexData <- rbind(hexData, interimHex)
      }
    }
    write_csv(hexData, str_replace(image_location, ".png", "_heatBins.csv"))
    # print(head(check_gray))
    # pause <- readline(prompt = "Pausing...4")
    # print(unique(check_gray$hex))
    
    for (viro in unique(check_gray$hex)){
      rgb_code <- col2rgb(viro)
      check_gray[check_gray$hex == viro,]$R <- rgb_code[1]
      check_gray[check_gray$hex == viro,]$G <- rgb_code[2]
      check_gray[check_gray$hex == viro,]$B <- rgb_code[3]
    }
    # print(head(check_gray))
    # pause <- readline(prompt = "Pausing...5")
    
    check_gray$R <- check_gray$R/255
    check_gray$G <- check_gray$G/255
    check_gray$B <- check_gray$B/255
    
    #  print(head(thick))
    if (!"cc" %in% names(thick)){
      print("Greyscale image detected.")
      thick$cc <- 1
      thick2 <- thick
      thick3 <- thick
      thick2$cc <- 2
      thick3$cc <- 3
      thick <- rbind(thick, thick2, thick3)
    }
    print("Applying new pallet")
    thick[thick$cc == 1,]$value <- check_gray$R
    print("R set")
    thick[thick$cc == 2,]$value <- check_gray$G
    print("G set")
    thick[thick$cc == 3,]$value <- check_gray$B
    print("B set")
    
    newPlot <- suppressWarnings(as.cimg(thick))
    newPlot <- cimg2magick(newPlot, rotate = T)
    newPlot <- image_flop(newPlot)
    image_location <- str_replace(image_location, ".png", "_binnedHeatMapped.png")
    image_write(newPlot, path = image_location, format = "png")
    print("Complete")
  } else {
    oldImage <- load.image(image_location)
    # plot(oldImage)
    # pause <- readline(prompt = "Pausing...1")
    thick <- as.data.frame(oldImage)
    if (4 %in% unique(thick$cc)){
      print("Alpha channel detected. Removing and re-saving.")
      thick <- thick[thick$cc != 4,]
      newPlot <- suppressWarnings(as.cimg(thick))
      newPlot <- cimg2magick(newPlot, rotate = T)
      newPlot <- image_flop(newPlot)
      image_write(newPlot, path = image_location, format = "png")
      oldImage <- load.image(image_location)
      thick <- as.data.frame(oldImage)
      print("3-Channel PNG saved.")
    }
    check_gray <- as.data.frame(grayscale(oldImage))
    
    check_gray$R <- 0
    check_gray$G <- 0
    check_gray$B <- 0
    # 
    # print(head(check_gray))
    # pause <- readline(prompt = "Pausing...2")
    
    print("Creating pallet hexes.")
    fake_plot <- ggplot(data = check_gray, aes(fill=value, x=x, y=y))+
      geom_point()+
      scale_fill_viridis_c(option = heatColor)
    fake_plot <- ggplot_build(fake_plot)
    fake_plot <- as.data.frame(fake_plot$data)
    # print(head(fake_plot))
    # pause <- readline(prompt = "Pausing...3")
    check_gray$hex <- fake_plot$fill
    # print(head(check_gray))
    # pause <- readline(prompt = "Pausing...4")
    # print(unique(check_gray$hex))
    
    for (viro in unique(check_gray$hex)){
      rgb_code <- col2rgb(viro)
      check_gray[check_gray$hex == viro,]$R <- rgb_code[1]
      check_gray[check_gray$hex == viro,]$G <- rgb_code[2]
      check_gray[check_gray$hex == viro,]$B <- rgb_code[3]
    }
    # print(head(check_gray))
    # pause <- readline(prompt = "Pausing...5")
    
    check_gray$R <- check_gray$R/255
    check_gray$G <- check_gray$G/255
    check_gray$B <- check_gray$B/255
    
    #  print(head(thick))
    if (!"cc" %in% names(thick)){
      print("Greyscale image detected.")
      thick$cc <- 1
      thick2 <- thick
      thick3 <- thick
      thick2$cc <- 2
      thick3$cc <- 3
      thick <- rbind(thick, thick2, thick3)
    }
    print("Applying new pallet")
    thick[thick$cc == 1,]$value <- check_gray$R
    print("R set")
    thick[thick$cc == 2,]$value <- check_gray$G
    print("G set")
    thick[thick$cc == 3,]$value <- check_gray$B
    print("B set")
    
    newPlot <- suppressWarnings(as.cimg(thick))
    newPlot <- cimg2magick(newPlot, rotate = T)
    newPlot <- image_flop(newPlot)
    image_location <- str_replace(image_location, ".png", paste0("_",heatColor,"_HeatMapped.png"))
    image_write(newPlot, path = image_location, format = "png")
    print("Complete")
  }
}

heatmap_it_all <- function(heatmap_color = "inferno", 
                           binned_colors = F){
  setwd("files/")
  folder_set <- list.files()
  for (i in folder_set){
    if (!grepl(".csv", i)){
      setwd(i)
      print(paste0("Working in file: ", i))
      image_set <- list.files()
      for (j in image_set){
        if (!grepl(".csv", j) & !grepl("ratio", j) & !grepl(heatmap_color, j)){
          setwd(paste0(j, "/PNGS"))
          print(paste0("Working in image: ", j))
          final_set <- list.files(pattern = ".png")
          for (k in final_set){
            if (!grepl("HeatMapped", k)){
              print(paste0("Converting: ", k))
              suppressWarnings(heatMap_image(k, heatColor = heatmap_color, intensity_binning = binned_colors))
            }
          }
          setwd("../../")
        }
      }
      setwd("../")
    }
  }
  setwd("../")
}