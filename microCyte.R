defaultW <- getOption("warn")
options(warn = -1)
partial_amount <- 1000
#--------------------------------------
# First, Lets check if all the packages are installed
recPackages <- c("ggplot2",
                 "ggpubr",
                 "ggtern",
                 "tidyverse",
                 "MASS",
                 "viridis",
                 "Rtsne",
                 "stringr",
                 "imager",
                 "magick",
                 "RecordLinkage",
                 "hash",
                 "plotly",
                 "Rtsne")


installedPackages <- installed.packages()
for (i in recPackages){
  if(!i %in% installedPackages){
    install.packages(i)
  }
}
#----------------------------------------
scripts <- data.frame(
  "Number" = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
  "Script_Name" = c("Setup",
                    "dirGen",
                    "purgo",
                    "reName",
                    "imaGen",
                    "concato",
                    "unite",
                    "full run",
                    "reset names",
                    "analyze")
)
print("Welcome to MicroCyte. Be sure to have edited your schema file before beginning. Here are your options:")
print(scripts)

while (TRUE) {
  option <- readline(prompt = paste0("What would you like to run (0-", max(scripts$Number), ", or 'end'): "))
  if (option == "end"){
    print("Thank you for using MicroCyte.")
    break
  } else if (option == "0"){
    
# Runs set up script
    print("Running setup now...")
    source("bin/scripts/scripts/0_setup.R")
    print("Set up complete")
  } else if (option == "1"){
    
# Runs Dirgen only
    print("Generating directories...")
    source("bin/scripts/scripts/1_dirGen.R")
    print("Directories complete")
  } else if (option == "2"){
    
# Runs Purgo and asks if you want to rename
    print("Beginning the purge...")
    source("bin/scripts/scripts/2_purgo.R")
    print("Purge complete")
    contun <- readline(prompt = "Would you like to now run reName (Y/n)? ")
    if (contun != "n"){
      print("Renaming images per the schema file...")
      source("bin/scripts/scripts/3_reName.R")
      reName()
      print("Images renamed. Please run YggData to continue.")
    }
  } else if (option == "3"){
    
# Runs ReName only
    print("Renaming images per the schema file...")
    source("bin/scripts/scripts/3_reName.R")
    reName()
    print("Images renamed. Please run YggData to continue.")
  } else if (option == "4"){
    
# Runs Imagen and starts a finish chain
    print("Running Imagen now...")
    source("bin/scripts/scripts/5_imaGen.R")
    print("Datasets combined.")
    contun <- readline(prompt = "Would you like to now run concato (Y/n)? ")
    if (contun != "n"){
      print("Concatenating image sets...")
      source("bin/scripts/scripts/6_concato.R")
      normal <- readline(prompt = "Should the numeric values be normalized among image sets (Y/n)? ")
      if (normal != "n"){
        concato(normalization = F)
      } else {
        concato()
      }
      print("Image sets combined")
      contun <- readline(prompt = "Would you like to now run Unite (Y/n)? ")
      if (contun != "n"){
        print("Uniting sample datasets...")
        source("bin/scripts/scripts/7_unite.R")
        
        full <- readline(prompt = "Should the entire dataset be used (FULL), or equivalent cells per sample (partial): ")
        if (full == "partial"){
          full <- F
          partial_amount <- as.numeric(readline(prompt = "How many cells should be pulled: "))
        } else {
          full <- T
        }
        
        patt <- readline(prompt = "Should samples be united on a variable (VAR), or altogether (all): ")
        if (patt != "all"){
          patt_2 <- readline(prompt = "Is the variable a part of the (name), or a (SCHEMA) variable: ")
          if (patt_2 == "name"){
            patt <- readline(prompt = "What is the exact name you wish to unite: ")
          }else {
            schema <- read.csv("schema.csv")
            print(names(schema))
            schemeName <- readline(prompt = "Which shema variable would you like to unite by: ")
            schemeName <- unique(schema[schemeName])
          }
        } else {
          patt <- F
        }
        
        if(patt != F){
          if (patt_2 == "name"){
            unite(full = full, pattern = patt, sampleSize = partial_amount)
          } else {
            for (schemaNames in schemeName){
              schemaNames <- as.character(schemaNames)
              for (nameo in schemaNames){
                unite(full = full, pattern = nameo, sampleSize = partial_amount)
              }
            }
          }
          
        } else {
          unite(full = full, sampleSize = partial_amount)
        }
        
        print("Sample datasets combined")
      }
    }
  } else if (option == "5"){
    
#Runs concato and finish chain
    print("Concatenating image sets...")
    source("bin/scripts/scripts/6_concato.R")
    normal <- readline(prompt = "Should the numeric values be normalized among image sets (Y/n)? ")
    if (normal != "n"){
      concato(normalization = F)
    } else {
      concato()
    }
    print("Image sets combined")
    contun <- readline(prompt = "Would you like to now run Unite (Y/n)? ")
    if (contun != "n"){
      print("Uniting sample datasets...")
      source("bin/scripts/scripts/7_unite.R")
      
      full <- readline(prompt = "Should the entire dataset be used (FULL), or equivalent cells per sample (partial): ")
      if (full == "partial"){
        full <- F
        partial_amount <- as.numeric(readline(prompt = "How many cells should be pulled: "))
      } else {
        full <- T
      }
      
      patt <- readline(prompt = "Should samples be united on a variable (VAR), or altogether (all): ")
      if (patt != "all"){
        patt_2 <- readline(prompt = "Is the variable a part of the (name), or a (SCHEMA) variable: ")
        if (patt_2 == "name"){
          patt <- readline(prompt = "What is the exact name you wish to unite: ")
        }else {
          schema <- read.csv("schema.csv")
          print(names(schema))
          schemeName <- readline(prompt = "Which shema variable would you like to unite by: ")
          schemeName <- unique(schema[schemeName])
        }
      } else {
        patt <- F
      }
      
      if(patt != F){
        if (patt_2 == "name"){
          unite(full = full, pattern = patt, sampleSize = partial_amount)
        } else {
          for (schemaNames in schemeName){
            schemaNames <- as.character(schemaNames)
            for (nameo in schemaNames){
              unite(full = full, pattern = nameo, sampleSize = partial_amount)
            }
          }
        }
        
      } else {
        unite(full = full, sampleSize = partial_amount)
      }
      
      
      print("Sample datasets combined")
    }
  } else if (option == "6"){
    
# Runs unite only
    print("Uniting sample datasets...")
    source("bin/scripts/scripts/7_unite.R")
    
    full <- readline(prompt = "Should the entire dataset be used (FULL), or equivalent cells per sample (partial): ")
    if (full == "partial"){
      full <- F
      partial_amount <- as.numeric(readline(prompt = "How many cells should be pulled: "))
    } else {
      full <- T
    }
    
    patt <- readline(prompt = "Should samples be united on a variable (VAR), or altogether (all): ")
    if (patt != "all"){
      patt_2 <- readline(prompt = "Is the variable a part of the (name), or a (SCHEMA) variable: ")
      if (patt_2 == "name"){
        patt <- readline(prompt = "What is the exact name you wish to unite: ")
      }else {
        schema <- read.csv("schema.csv")
        print(names(schema))
        schemeName <- readline(prompt = "Which shema variable would you like to unite by: ")
        schemeName <- unique(schema[schemeName])
      }
    } else {
      patt <- F
    }
    
    if(patt != F){
      if (patt_2 == "name"){
        unite(full = full, pattern = patt, sampleSize = partial_amount)
      } else {
        for (schemaNames in schemeName){
          schemaNames <- as.character(schemaNames)
          for (nameo in schemaNames){
            unite(full = full, pattern = nameo, sampleSize = partial_amount)
          }
        }
      }
      
    } else {
      unite(full = full, sampleSize = partial_amount)
    }
    
    print("Sample datasets combined")
  } else if (option == "7"){
    
#Runs everything
    print("Running the entire process")
    
    normal <- readline(prompt = "Should the numeric values be normalized among image sets (Y/n)? ")
    full <- readline(prompt = "Should the entire dataset be used (FULL), or equivalent cells per sample (partial): ")
    if (full == "partial"){
      full <- F
      partial_amount <- as.numeric(readline(prompt = "How many cells should be pulled: "))
    } else {
      full <- T
    }
    
    patt <- readline(prompt = "Should samples be united on a variable (VAR), or altogether (all): ")
    if (patt != "all"){
      patt_2 <- readline(prompt = "Is the variable a part of the (name), or a (SCHEMA) variable: ")
      if (patt_2 == "name"){
        patt <- readline(prompt = "What is the exact name you wish to unite: ")
      }else {
        schema <- read.csv("schema.csv")
        print(names(schema))
        schemeName <- readline(prompt = "Which shema variable would you like to unite by: ")
        schemeName <- unique(schema[schemeName])
      }
    } else {
      patt <- F
    }
    
    dirOverride <- readline(prompt = "Should directories be generated (Y/n)? ")
    if (dirOverride != "n"){
      print("Generating directories...")
      source("bin/scripts/scripts/1_dirGen.R")
      print("Directories complete")
      
      print("Please take your images and continue when ready...")
      pause <- readline(prompt = "Continue:")
    }
    
    purgeOverride <- readline(prompt = "Have images been taken (Y/n)? ")
    if (purgeOverride != "n"){
      print("Beginning the purge...")
      source("bin/scripts/scripts/2_purgo.R")
      print("Purge complete")
    
      print("Renaming images per the schema file...")
      source("bin/scripts/scripts/3_reName.R")
      reName()
      print("Images renamed. Please run YggData and continue when finished...")
      pause <- readline(prompt = "Continue:")
    }
    
    print("Running Imagen now...")
    source("bin/scripts/scripts/5_imaGen.R")
    print("Datasets combined.")
    
    print("Concatenating image sets...")
    source("bin/scripts/scripts/6_concato.R")
    if (normal != "n"){
      concato(normalization = F)
    } else {
      concato()
    }
    print("Image sets combined")
    
    print("Uniting sample datasets...")
    source("bin/scripts/scripts/7_unite.R")
    if(patt != F){
      if (patt_2 == "name"){
        unite(full = full, pattern = patt, sampleSize = partial_amount)
      } else {
        for (schemaNames in schemeName){
          schemaNames <- as.character(schemaNames)
          for (nameo in schemaNames){
            unite(full = full, pattern = nameo, sampleSize = partial_amount)
          }
        }
      }
      
    } else {
      unite(full = full, sampleSize = partial_amount)
    }
    print("Sample datasets combined")
    
  } else if (option == "8"){
    source("bin/scripts/scripts/3_reName.R")
    reset_names()
  } else if (option == "9"){
    scriptList <- list.files("bin/scripts/analysis/", pattern = ".R")
    print("Sourcing the analysis scripts...")
    cat("\n")
    if (length(scriptList) == 0){
      print("No scripts detected. Make sure you are in the right directory and try again.")
    } else {
      for(scriptin in scriptList){
        cat(paste0("Sourcing script: ", scriptin))
        cat("\n")
        source(paste0("bin/scripts/analysis/", scriptin))
      }
    }
    cat("\n")
    print("Good luck with your analysis!")
    break
  } else {
    print("Sorry, I didn't catch that. Please try again, or type 'end' to finish.")
    print(scripts)
  }
}

options(warn = defaultW)