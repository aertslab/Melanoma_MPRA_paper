get_enhancer_counts<-function(enhancer.barcode.list, 
                              barcode.counts, 
                              column.name) {
  enhancer.counts<-lapply(X = enhancer.barcode.list, FUN = function(bcs.str) {
    barcodes<-unlist(strsplit(bcs.str, ","))
    mask<-barcode.counts$V1 %in% barcodes
    if(sum(mask) > 0) {
      return (barcode.counts[[column.name]][mask])
    }
  })
  names(enhancer.counts)<-p5P.enhancer.barcodes$V1
  return(enhancer.counts)
}


ProcessCHEQseq <- function(env, 
                           line, 
                           library, #"long" or "500bp"
                           plasmid.counts.file.path, 
                           cDNA.counts.file.path, 
                           filter.method = "min.counts", # or min.barcodes
                           min.counts = 5,
                           min.barcodes = 5,
                           barcode.matching = FALSE,
                           aggregate.barcodes.by.enhancer = FALSE) {
  library(dplyr)
  # Define library type
  proliferative.enhancers.long<-c("CDH1_24-I","CDH1_49-I","SOX10_15-3U","SOX10_-3-D","SOX10_-34-D","SOX10_-35-D","SOX10_-56-D","IRF4_4-I","GPM6B_P","KIT_114-D","MLANA_5-I","SGCD_15-I","SGCD_26-I","MITF_P","TYR_-9-D","TYR_-1-D","TYR_P","RRAGD_P")
  proliferative.enhancers.500bp <- c("CDH1_24-I", "CDH1_49-I", "SOX10_15-3UA", "SOX10_15-3UB", "SOX10_-34-DA", "SOX10_-34-DB", "SOX10_-56-DA", "SOX10_-56-DB", "IRF4_4-I", "KIT_114-DA", "KIT_114-DB", "MLANA_5-I", "SGCD_15-IA", "SGCD_15-IB", "SGCD_26-I", "MITF_P", "TYR_-9-D", "TYR_-1-D")  
  if(library == "long"){
    warning("Processing long enhancer library data... \n")
    proliferative.enhancers <- proliferative.enhancers.long
  } else if(library == "500bp") {
    warning("Processing 500bp enhancer library data... \n")
    proliferative.enhancers <- proliferative.enhancers.500bp
  } else {
    stop("Incorrect library name. Must be 'long' or '500bp'")
  }

  # Create environment and load raw data
  env[[line]] <- list()
  env[[line]]$plasmid.raw <- read.table(file = plasmid.counts.file.path, header = F, sep = "\t", quote = '', stringsAsFactors = F) %>%
    magrittr::set_colnames(c("enhancer", "barcode", "counts"))
  env[[line]]$cDNA.raw  <- read.table(file = cDNA.counts.file.path, header = F, sep = "\t", quote = '', stringsAsFactors = F) %>%
    magrittr::set_colnames(c("enhancer", "barcode", "counts"))
  
  # Rename incorrect enhancer name
  if(library == "500bp"){
    env[[line]]$plasmid.raw[env[[line]]$plasmid.raw$enhancer == "ALPK_2","enhancer"] <- "ALPK2_2"
    env[[line]]$cDNA.raw[env[[line]]$cDNA.raw$enhancer == "ALPK_2","enhancer"] <- "ALPK2_2"
  }
  
  # Rename enhancers with new names
  new.names <- read.table(file = "/staging/leuven/stg_00002/lcb/lcb_projects/CH1/5p_intron/analysis/region_new_names.tsv", header = T, sep = "\t")
  new.names.500 <- read.table(file = "/staging/leuven/stg_00002/lcb/lcb_projects/CH1/5p_intron/analysis/region_new_names_500.tsv", header = T, sep = "\t")
  if (library == "long") {
    for (i in new.names$old_names){
      env[[line]]$plasmid.raw[env[[line]]$plasmid.raw$enhancer == i,"enhancer"] <- as.character(new.names[new.names$old_names == i, "new_names"])
      env[[line]]$cDNA.raw[env[[line]]$cDNA.raw$enhancer == i,"enhancer"] <- as.character(new.names[new.names$old_names == i, "new_names"])
    }
  } else {
    for (i in new.names.500$old_names){
      env[[line]]$plasmid.raw[env[[line]]$plasmid.raw$enhancer == i,"enhancer"] <- as.character(new.names.500[new.names.500$old_names == i, "new_names"])
      env[[line]]$cDNA.raw[env[[line]]$cDNA.raw$enhancer == i,"enhancer"] <- as.character(new.names.500[new.names.500$old_names == i, "new_names"])
    }
  }
  
  if(barcode.matching) {
    # Combine matching barcodes
    warning("Filtering BC matching in cDNA and Plasmid... \n")
    env[[line]]$merged <- inner_join(env[[line]]$plasmid.raw, 
                                     env[[line]]$cDNA.raw, 
                                     by = "barcode", 
                                     suffix = c(".plasmid", ".cDNA"))
    env[[line]]$merged <- env[[line]]$merged[,-which(names(env[[line]]$merged) %in% c("enhancer.cDNA"))] %>% 
      rename(enhancer = enhancer.plasmid)
    env[[line]]$merged$enhancer <- as.factor(env[[line]]$merged$enhancer)
  } else {
    # Aggregate all barcodes by enhancer before CPM and input normalization 
    warning("Aggregating all BC counts by enhancer... \n")
    env[[line]]$plasmid.raw <- env[[line]]$plasmid.raw %>% 
      group_by(enhancer) %>%
      summarise(sum(counts),
                n()) %>%
      magrittr::set_colnames(c("enhancer", "counts", "num.barcodes"))
    env[[line]]$cDNA.raw <- env[[line]]$cDNA.raw %>% 
      group_by(enhancer) %>%
      summarise(sum(counts),
                n()) %>%
      magrittr::set_colnames(c("enhancer", "counts", "num.barcodes"))
    env[[line]]$merged <- inner_join(env[[line]]$plasmid.raw, 
                                                  env[[line]]$cDNA.raw, 
                                                  by = "enhancer", 
                                                  suffix = c(".plasmid", ".cDNA"))
    env[[line]]$merged$enhancer <- as.factor(env[[line]]$merged$enhancer)
  }
  
  # Add phenotype info
  env[[line]]$merged$phenotype <- "Invasive"
  env[[line]]$merged$phenotype[env[[line]]$merged$enhancer %in% proliferative.enhancers] <- "Proliferative"
  env[[line]]$merged$phenotype[env[[line]]$merged$enhancer == "NEG_CTRL"] <- "Control"
  
  if(filter.method == "min.counts") {
    warning("Apply min.counts filtering... \n")
    # Filter the enhancer-barcode based on the counts
    env[[line]]$merged <- env[[line]]$merged %>%
      filter(counts.plasmid >= min.counts & counts.cDNA >= min.counts)
  } else if (filter.method == "min.barcodes" & barcode.matching) {
    stop("min.barcodes filtering on barcode matched enhancer is not available")
      } else if (filter.method == "min.barcodes") {
    warning("Apply min.barcodes filtering...\n")
    # Filter on the number of barcodes / enhancers
    env[[line]]$merged <- env[[line]]$merged %>%
      filter(num.barcodes.plasmid >= min.barcodes & num.barcodes.cDNA >= min.barcodes)
     } else {
    stop("This filtering method is not implemented. Available methods: min.counts, min.barcodes")
  }

  warning("Input normalization... \n")
  library(preprocessCore)
  library(affy)
  if(!barcode.matching) {
    ## Normalize aggregated values
    # CPM Input normalization
    env[[line]]$merged.cpm.input.norm <- data.frame("enhancer"=env[[line]]$merged$enhancer, 
                                                                 stringsAsFactors = F)
    y <- edgeR::cpm(y = env[[line]]$merged[,c("counts.plasmid","counts.cDNA")])
    row.names(x = y) <- env[[line]]$merged$enhancer
    env[[line]]$merged.cpm.input.norm[["CPMNorm.counts.plasmid"]] <- y[,1]
    env[[line]]$merged.cpm.input.norm[["CPMNorm.counts.cDNA"]] <- y[,2]
    env[[line]]$merged.cpm.input.norm[["CPMNorm.FC"]] <- y[,2]/y[,1]
  } else if(aggregate.barcodes.by.enhancer) {
    # Aggregate all barcodes by enhancer before CPM and input normalization
    env[[line]]$merged <- env[[line]]$merged %>%
      group_by(enhancer) %>%
      summarise(sum(counts.plasmid),
                sum(counts.cDNA)) %>%
      magrittr::set_colnames(c("enhancer", "counts.plasmid", "counts.cDNA"))
    # CPM input normalization
    env[[line]]$merged.cpm.input.norm <- data.frame("enhancer"=env[[line]]$merged$enhancer, 
                                                                 stringsAsFactors = F)
    y <- edgeR::cpm(y = env[[line]]$merged[,c(2,3)])
    row.names(x = y) <- env[[line]]$merged$enhancer
    env[[line]]$merged.cpm.input.norm[["CPMNorm.counts.plasmid"]] <- y[,1]
    env[[line]]$merged.cpm.input.norm[["CPMNorm.counts.cDNA"]] <- y[,2]
    env[[line]]$merged.cpm.input.norm[["CPMNorm.FC"]] <- y[,2]/y[,1]
  } else {
    ## Normalize BC matched values
    # CPM Input normalization
    env[[line]]$merged.cpm.input.norm <- data.frame("enhancer"=env[[line]]$merged$enhancer, 
                                                                 "barcode"=env[[line]]$merged$barcode, 
                                                                 stringsAsFactors = F)
    y <- edgeR::cpm(y = env[[line]]$merged[,c("counts.plasmid","counts.cDNA")])
    row.names(x = y) <- env[[line]]$merged$barcode
    env[[line]]$merged.cpm.input.norm[["CPMNorm.counts.plasmid"]] <- y[,1]
    env[[line]]$merged.cpm.input.norm[["CPMNorm.counts.cDNA"]] <- y[,2]
    env[[line]]$merged.cpm.input.norm[["CPMNorm.FC"]] <- y[,2]/y[,1]
  }
  
  # Add phenotype info
  env[[line]]$merged.cpm.input.norm$phenotype <- "Invasive"
  env[[line]]$merged.cpm.input.norm$phenotype[env[[line]]$merged.cpm.input.norm$enhancer %in% proliferative.enhancers] <- "Proliferative"
  env[[line]]$merged.cpm.input.norm$phenotype[env[[line]]$merged.cpm.input.norm$enhancer == "NEG_CTRL"] <- "Control"
  
  if(!barcode.matching && "NEG_CTRL" %in% env[[line]]$merged.cpm.input.norm$enhancer) {
    warning("Basal activity normalization... \n")
    env[[line]]$merged.cpm.input.basal.norm <- data.frame("enhancer"=env[[line]]$merged$enhancer, 
                                                          stringsAsFactors = F)
    env[[line]]$merged.cpm.input.basal.norm[["BasalNorm.FC"]] <- env[[line]]$merged.cpm.input.norm[["CPMNorm.FC"]]/env[[line]]$merged.cpm.input.norm[env[[line]]$merged.cpm.input.norm$enhancer == "NEG_CTRL","CPMNorm.FC"]
    env[[line]]$merged.cpm.input.basal.norm <- env[[line]]$merged.cpm.input.basal.norm[env[[line]]$merged.cpm.input.basal.norm$enhancer != "NEG_CTRL",]
    env[[line]]$merged.cpm.input.basal.norm$enhancer <- droplevels(env[[line]]$merged.cpm.input.basal.norm$enhancer)
    env[[line]]$meta.data$activity.thresh <- env[[line]]$merged.cpm.input.norm[env[[line]]$merged.cpm.input.norm$enhancer == "NEG_CTRL","CPMNorm.FC"]
    
    # Add phenotype info
    env[[line]]$merged.cpm.input.basal.norm$phenotype <- "Invasive"
    env[[line]]$merged.cpm.input.basal.norm$phenotype[env[[line]]$merged.cpm.input.basal.norm$enhancer%in%proliferative.enhancers] <- "Proliferative"
    
  } else {

    # Define activity threshold (Mean FC of the opposite phenotype enhancers * 1.5)
    invasive.lines <- c("MM029", "MM047", "MM099")
    warning("Calculating activity threshold... \n")  
    if(line %in% invasive.lines){
      x <- env[[line]]$merged.cpm.input.norm[env[[line]]$merged.cpm.input.norm$phenotype == "Proliferative",] %>% arrange(CPMNorm.FC)
      env[[line]]$meta.data$activity.thresh <- mean(x[1:round(x = nrow(x), digits = 0),"CPMNorm.FC"])*1.5
    } else {
      x <- env[[line]]$merged.cpm.input.norm[env[[line]]$merged.cpm.input.norm$phenotype == "Invasive",] %>% arrange(CPMNorm.FC)
      env[[line]]$meta.data$activity.thresh <- mean(x[1:round(x = nrow(x), digits = 0),"CPMNorm.FC"])*1.5
    }
  }
  
  # Add total read count per cell line to meta.data
  env[[line]]$meta.data$sum.counts.plasmid <- sum(env[[line]]$plasmid.raw$counts)
  env[[line]]$meta.data$sum.counts.cDNA <- sum(env[[line]]$cDNA.raw$counts)
  return (env)
}


ProcessSTARRseq <- function(env, project, 
                            line, 
                            plasmid.counts.file.path, 
                            cDNA.counts.file.path, 
                            out.dir) {
  library(dplyr)
  warning("Processing ", line, " STARRseq data...\n")
  # Load data
  env[[line]] <- new.env()
  env[[line]]$plasmid <- read.table(file = plasmid.counts.file.path, header = F, sep = "\t", stringsAsFactors = F)
  colnames(env[[line]]$plasmid) <- c("enhancer", "BC", "counts")
  env[[line]]$cDNA <- read.table(file = cDNA.counts.file.path, header = F, sep = "\t", stringsAsFactors = F)
  colnames(env[[line]]$cDNA) <- c("enhancer", "BC", "counts")
  
  # Rename incorrect enhancer name
  env[[line]]$plasmid[env[[line]]$plasmid$enhancer == "ALPK_2","enhancer"] <- "ALPK2_2"
  env[[line]]$cDNA[env[[line]]$cDNA$enhancer == "ALPK_2","enhancer"] <- "ALPK2_2"
  
  # Rename enhancers with new names
  new.names.500 <- read.table(file = "/staging/leuven/stg_00002/lcb/lcb_projects/CH1/5p_intron/analysis/region_new_names_500.tsv", header = T, sep = "\t")
  for (i in new.names.500$old_names) {
    env[[line]]$plasmid[env[[line]]$plasmid$enhancer == i,"enhancer"] <- as.character(new.names.500[new.names.500$old_names == i, "new_names"])
    env[[line]]$cDNA[env[[line]]$cDNA$enhancer == i,"enhancer"] <- as.character(new.names.500[new.names.500$old_names == i, "new_names"])
  }
  
  # Merge plasmid and cDNA counts
  warning("Merging plasmid and cDNA counts...\n")
  env[[line]]$merged <- inner_join(env[[line]]$plasmid,
                                   env[[line]]$cDNA,
                                   by = "enhancer",
                                   suffix = c(".plasmid", ".cDNA"))
  
  # Normalize
  warning("Normalizing...\n")
  library(preprocessCore)
  library(affy)
  ## CPM
  y <- edgeR::cpm(y = env[[line]]$merged[,c(3,5)])
  env[[line]]$merged$CPMNorm.counts.plasmid <- y[,1]
  env[[line]]$merged$CPMNorm.counts.cDNA <- y[,2]
  
  # Input normalization
  env[[line]]$merged.input.norm <- data.frame("enhancer"= env[[line]]$merged$enhancer)
  env[[line]]$merged.input.norm[["CPMNorm.FC"]] <- env[[line]]$merged$CPMNorm.counts.cDNA/env[[line]]$merged$CPMNorm.counts.plasmid
  proliferative.enhancers<-c("CDH1_24-I", "CDH1_49-I", "SOX10_15-3UA", "SOX10_15-3UB", "SOX10_-34-DA", "SOX10_-34-DB", "SOX10_-56-DA", "SOX10_-56-DB", "IRF4_4-I", "KIT_114-DA", "KIT_114-DB", "MLANA_5-I", "SGCD_15-IA", "SGCD_15-IB", "SGCD_26-I", "MITF_P", "TYR_-9-D", "TYR_-1-D")
  env[[line]]$merged.input.norm$phenotype <- "Invasive"
  env[[line]]$merged.input.norm$phenotype[env[[line]]$merged.input.norm$enhancer %in% proliferative.enhancers] <- "Proliferative"
  env[[line]]$merged.input.norm$phenotype[env[[line]]$merged.input.norm$enhancer == "NEG_CTRL"] <- "Control"
  # Save
  warning("Saving...\n")
  write.table(x = env[[line]]$merged.input.norm, file = file.path(out.dir, paste0(project, "__",line,"__cpm_input_normalized.tsv")), quote = F, sep = "\t", row.names = F, col.names = T)
  return (env)
}