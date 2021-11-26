process.OLS.plasmid <- function(env, 
                                plasmid.counts.file.path, 
                                min.plasmid.bcs = 5,
                                out.dir) {
  warning("Processing Plasmid input data with filtering set to ", min.plasmid.bcs," barcodes.\n")
  env[["Plasmid"]] <- new.env()
  env[["Plasmid"]]$plasmid.raw <- read.table(file = plasmid.counts.file.path, header = F, sep = "\t", stringsAsFactors = F)
  
  # Aggegate the plasmid counts by enhancer tile
  env[["Plasmid"]]$plasmid.raw %>%
    group_by(V1) %>% 
    dplyr::summarise(
      counts = sum(V3),
      nb_barcodes = length(V1)
    ) %>%
    filter(nb_barcodes >= min.plasmid.bcs) %>%
    magrittr::set_colnames(c("synthetic.region", "counts", "nb_barcodes")) -> env[["Plasmid"]]$plasmid.aggregated
  warning("Number of tiles remaining after filtering: ", nrow(env[["Plasmid"]]$plasmid.aggregated),"\n")
  
  # CPM Normalization
  warning("CPM normalizing...")
  library(preprocessCore)
  library(affy)
  env[["Plasmid"]]$plasmid.aggregated$CPM.counts <- edgeR::cpm(y = env[["Plasmid"]]$plasmid.aggregated[,"counts"])
  warning("[OK]\n")
  
  # Assign tile class
  env[["Plasmid"]]$plasmid.aggregated$class <- "unassigned"
  env[["Plasmid"]]$plasmid.aggregated[grep(pattern = "shuffled", x = env[["Plasmid"]]$plasmid.aggregated$synthetic.region),"class"] <- "shuffled"
  env[["Plasmid"]]$plasmid.aggregated[grep(pattern = "==wt$", x = env[["Plasmid"]]$plasmid.aggregated$synthetic.region), "class"] <- "wt"
  env[["Plasmid"]]$plasmid.aggregated[grep(pattern = "==mut___", x = env[["Plasmid"]]$plasmid.aggregated$synthetic.region), "class"] <- "mut"
  # Adding metadata
  env[["Plasmid"]]$metadata <- list()
  env[["Plasmid"]]$metadata$BC.plasmid.filtering <- min.plasmid.bcs
  
  warning("Saving...")
  saveRDS(object = env$Plasmid, file = file.path(out.dir, paste0("plasmid_", environmentName(env), ".RDS")))
  warning("[OK]")
  
  hist(x = env[["Plasmid"]]$plasmid.aggregated$nb_barcodes, breaks = 200)
  abline(v = min.plasmid.bcs, col = "red")
  return (env)
}


process.OLS.cDNA <- function(env, 
                             line, 
                             cDNA.counts.file.path, 
                             min.cDNA.bcs = 5, 
                             out.dir = NULL) {
  if (!("Plasmid" %in% names(env))){
    stop("Process first the plasmid data with process.OLS.plasmid()\n")
  } else {
    
    warning("Processing ", line, " with filtering set to ", min.cDNA.bcs," barcodes for cDNA.\n")
    env[[line]] <- new.env()
    env[[line]]$cDNA.raw <- read.table(file = cDNA.counts.file.path, header = F, sep = "\t", stringsAsFactors = F)
    
    # Aggregate the cDNA counts by enhancer tile
    library(dplyr)
    env[[line]]$cDNA.raw %>% 
      group_by(V1) %>% 
      dplyr::summarise(
        counts = sum(V3),
        nb_barcodes = length(V1)
      ) %>%
      filter(nb_barcodes >= min.cDNA.bcs) %>%
      magrittr::set_colnames(c("synthetic.region", "counts", "nb_barcodes")) -> env[[line]]$cDNA.aggregated
    
    # CPM Normalization
    warning("CPM normalizing...")
    library(preprocessCore)
    library(affy)
    env[[line]]$cDNA.aggregated$CPM.counts <- edgeR::cpm(y = env[[line]]$cDNA.aggregated[,"counts"])
    warning("[OK]\n")
    
    # Merge plasmid and cDNA counts
    warning("Merging plasmid and cDNA counts...")
    env[[line]]$merged <- inner_join(env[["Plasmid"]]$plasmid.aggregated, 
                                     env[[line]]$cDNA.aggregated, 
                                     by = "synthetic.region", 
                                     suffix = c(".plasmid", ".cDNA"))
    warning("[OK]\n")
    warning("Number of tiles remaining after combining cDNA and Plasmid: ", nrow(env[[line]]$merged),"\n")
    warning("cDNA reads:    Median = ", median(x = env[[line]]$merged$counts.cDNA),"  Mean = ", mean(x = env[[line]]$merged$counts.cDNA), "\n")
    warning("Plasmid reads: Median = ", median(x = env[[line]]$merged$counts.plasmid),"  Mean = ", mean(x = env[[line]]$merged$counts.plasmid), "\n")
 
    env <- process.OLS.part2(env = env, 
                             line = line,
                             min.cDNA.bcs = min.cDNA.bcs)
    
    env[[line]]$metadata$BC.cDNA.filtering <- min.cDNA.bcs
    env[[line]]$metadata$BC.plasmid.filtering <- env[["Plasmid"]]$metadata$BC.plasmid.filtering
    
    # Save
    if (!is.null(out.dir)){
      write.table(x = env[[line]]$merged.input.basal.norm[,c(1,3)],
                  file = file.path(out.dir, paste0(line,"__filter_",min.cDNA.bcs,"-",env[["Plasmid"]]$metadata$BC.plasmid.filtering,"-pooled__cpm_basal_normalized.tsv")),
                  quote = F, sep = "\t", row.names = F, col.names = T)
    }
  }
  return (env)
} 


process.OLS.cDNA.plasmid <- function(env, 
                                     line, 
                                     cDNA.counts.file.path, 
                                     plasmid.counts.file.path, 
                                     min.plasmid.bcs = 5,
                                     min.cDNA.bcs = 5,
                                     tiling, #"A" or "B"
                                     keep.shuffled = T,
                                     out.dir = NULL) {
  warning("Processing ", line, " with filtering set to ", min.cDNA.bcs," barcodes for cDNA and ", min.plasmid.bcs," barcodes for plasmid.\n")
  env[[line]] <- new.env()
  env[[line]]$plasmid.raw <- read.table(file = plasmid.counts.file.path, header = F, sep = "\t", stringsAsFactors = F)
  env[[line]]$cDNA.raw <- read.table(file = cDNA.counts.file.path, header = F, sep = "\t", stringsAsFactors = F)

  # Aggregate counts by enhancer tile
  library(dplyr)
  env[[line]]$plasmid.raw %>%
    group_by(V1) %>% 
    dplyr::summarise(
      counts = sum(V3),
      nb_barcodes = length(V1)
    ) %>%
    filter(nb_barcodes >= min.plasmid.bcs) %>%
    magrittr::set_colnames(c("synthetic.region", "counts", "nb_barcodes")) -> env[[line]]$plasmid.aggregated
  env[[line]]$cDNA.raw %>% 
    group_by(V1) %>% 
    dplyr::summarise(
      counts = sum(V3),
      nb_barcodes = length(V1)
    ) %>%
    filter(nb_barcodes >= min.cDNA.bcs) %>%
    magrittr::set_colnames(c("synthetic.region", "counts", "nb_barcodes")) -> env[[line]]$cDNA.aggregated

  # CPM Normalization
  warning("CPM normalizing...")
  library(preprocessCore)
  library(affy)
  env[[line]]$plasmid.aggregated$CPM.counts <- edgeR::cpm(y = env[[line]]$plasmid.aggregated[,"counts"])
  env[[line]]$cDNA.aggregated$CPM.counts <- edgeR::cpm(y = env[[line]]$cDNA.aggregated[,"counts"])
  warning("[OK]\n")

  # Merge plasmid and cDNA counts
  warning("Merging plasmid and cDNA counts...")
  env[[line]]$merged <- inner_join(env[[line]]$plasmid.aggregated, 
                                   env[[line]]$cDNA.aggregated, 
                                   by = "synthetic.region", 
                                   suffix = c(".plasmid", ".cDNA"))
  
  # Assign tile class
  env[[line]]$merged$class <- "unassigned"
  env[[line]]$merged[grep(pattern = "shuffled", x = env[[line]]$merged$synthetic.region),"class"] <- "shuffled"
  env[[line]]$merged[grep(pattern = "==wt$", x = env[[line]]$merged$synthetic.region), "class"] <- "wt"
  env[[line]]$merged[grep(pattern = "==mut___", x = env[[line]]$merged$synthetic.region), "class"] <- "mut"
  
  # Assign phenotype
  library(stringr)
  proliferative<-c("CDH1_1","CDH1_2","SOX10_1","SOX10_2","SOX10_3","SOX10_4","SOX10_5","IRF4_1","GPM6B_2","KIT_1","MLANA_1","SGCD_2","SGCD_3","MITF_1","TYR_1","TYR_2","TYR_3","RRAGD_1")
  enh <- str_split(env[[line]]$merged$synthetic.region, '@@',simplify=T)[,2] %>% { str_split( . ,  '==',simplify=T)[,1] }
  env[[line]]$merged$phenotype <- "Invasive"
  env[[line]]$merged[str_split(env[[line]]$merged$synthetic.region, '@@',simplify=T)[,2] %>% { str_split( . ,  '==',simplify=T)[,1] } %in% proliferative, "phenotype"] <- "Proliferative"
  env[[line]]$merged[grep(pattern = "shuffled", x = env[[line]]$merged$synthetic.region),"phenotype"] <- "Shuffled"

  
  warning("[OK]\n")
  warning("Number of tiles remaining after combining cDNA and Plasmid: ", nrow(env[[line]]$merged),"\n")
  warning("cDNA reads:    Median = ", median(x = env[[line]]$merged$counts.cDNA),"  Mean = ", mean(x = env[[line]]$merged$counts.cDNA), "\n")
  warning("Plasmid reads: Median = ", median(x = env[[line]]$merged$counts.plasmid),"  Mean = ", mean(x = env[[line]]$merged$counts.plasmid), "\n")

  env <- process.OLS.part2(env = env, 
                           line = line,
                           min.cDNA.bcs = min.cDNA.bcs,
                           tiling = tiling,
                           keep.shuffled = keep.shuffled)
  
  env[[line]]$metadata$BC.cDNA.filtering <- min.cDNA.bcs
  env[[line]]$metadata$BC.plasmid.filtering <- min.plasmid.bcs
  
  # Save
  if (!is.null(out.dir)){
    if (keep.shuffled == T){
    write.table(x = env[[line]]$merged.input.basal.norm[, c("synthetic.region", "log.CPM.Input.BasalNorm")],
                file = file.path(out.dir, paste0(line,"_",environmentName(env),"__filter_",min.cDNA.bcs,"-",min.plasmid.bcs,"__cpm_basal_normalized.tsv")),
                quote = F, sep = "\t", row.names = F, col.names = T)
    }else{
      write.table(x = env[[line]]$merged.input.basal.norm[, c("synthetic.region", "log.CPM.Input.BasalNorm")],
                  file = file.path(out.dir, paste0(line,"_",environmentName(env),"__filter_",min.cDNA.bcs,"-",min.plasmid.bcs,"_without_shuffle__cpm_basal_normalized.tsv")),
                  quote = F, sep = "\t", row.names = F, col.names = T)
    }
  }
  return (env)
} 


process.OLS.part2 <- function(env, 
                              line,
                              min.cDNA.bcs,
                              tiling,
                              keep.shuffled) {
  warning("Normalizing activity...")
  # Input normalization
  env[[line]]$merged.input.norm <- data.frame("synthetic.region" = env[[line]]$merged$synthetic.region,
                                              "class" = env[[line]]$merged$class,
                                              "phenotype" = env[[line]]$merged$phenotype,
                                              stringsAsFactors = F)
  env[[line]]$merged.input.norm[["CPM.InputNorm"]] <- unlist(env[[line]]$merged[,"CPM.counts.cDNA"]/env[[line]]$merged[,"CPM.counts.plasmid"])
  
  # Basal expression normalization
  shuffled <- env[[line]]$merged.input.norm[grep(pattern = "shuffled", x = env[[line]]$merged.input.norm$synthetic.region),]
  basal.ex <- median(shuffled[,"CPM.InputNorm"])
  basal.mean <- mean(shuffled[,"CPM.InputNorm"])
  env[[line]]$merged.input.basal.norm <- data.frame("synthetic.region" = env[[line]]$merged.input.norm$synthetic.region,
                                                    "class" = env[[line]]$merged.input.norm$class,
                                                    "phenotype" = env[[line]]$merged$phenotype,
                                                    stringsAsFactors = F)
  env[[line]]$merged.input.basal.norm$CPM.Input.BasalNorm <- env[[line]]$merged.input.norm[,"CPM.InputNorm"]/basal.ex
  
  # Remove shuffled tiles
  if (keep.shuffled == F) {
  env[[line]]$merged.input.basal.norm <- env[[line]]$merged.input.basal.norm[!(env[[line]]$merged.input.basal.norm$class == "shuffled"),]
  }
  env[[line]]$merged.input.basal.norm$log.CPM.Input.BasalNorm <- log2(env[[line]]$merged.input.basal.norm$CPM.Input.BasalNorm)
  
  # Adding metadata
  env[[line]]$metadata <- list()
  env[[line]]$metadata$Basal.expression.median <- basal.ex
  env[[line]]$metadata$Basal.expression.mean <- basal.mean

  warning("[OK]\n")
  
  # Adding data to combined environment
  warning("Combining data with other samples...")
  ## Create combined dataset if it is not there yet
  if (!("Combined" %in% names(env))){
    ## With basal normalized dataframe
    env[["Combined"]] <- new.env()
    env[["Combined"]]$Input.normalised <- data.frame("synthetic.region" = unlist(read.table(file = paste0('/staging/leuven/stg_00002/lcb/dmauduit_HPC/CHEQseq_synth/Final_libraries/OLS-',tiling,'_tiles_name.txt'),
                                                                                     stringsAsFactors = F)),
                                                     stringsAsFactors = F)
    
    # Assign tile class
    env[["Combined"]]$Input.normalised$class <- "unassigned"
    env[["Combined"]]$Input.normalised[grep(pattern = "shuffled", x = env[["Combined"]]$Input.normalised$synthetic.region),"class"] <- "shuffled"
    env[["Combined"]]$Input.normalised[grep(pattern = "==wt$", x = env[["Combined"]]$Input.normalised$synthetic.region), "class"] <- "wt"
    env[["Combined"]]$Input.normalised[grep(pattern = "==mut___", x = env[["Combined"]]$Input.normalised$synthetic.region), "class"] <- "mut"

    # Assign phenotype
    proliferative<-c("CDH1_1","CDH1_2","SOX10_1","SOX10_2","SOX10_3","SOX10_4","SOX10_5","IRF4_1","GPM6B_2","KIT_1","MLANA_1","SGCD_2","SGCD_3","MITF_1","TYR_1","TYR_2","TYR_3","RRAGD_1")
    enh <- str_split(env[["Combined"]]$Input.normalised$synthetic.region, '@@',simplify=T)[,2] %>% { str_split( . ,  '==',simplify=T)[,1] }
    env[["Combined"]]$Input.normalised$phenotype <- "Invasive"
    env[["Combined"]]$Input.normalised[str_split(env[["Combined"]]$Input.normalised$synthetic.region, '@@',simplify=T)[,2] %>% { str_split( . ,  '==',simplify=T)[,1] } %in% proliferative, "phenotype"] <- "Proliferative"
    env[["Combined"]]$Input.normalised[grep(pattern = "shuffled", x = env[["Combined"]]$Input.normalised$synthetic.region),"phenotype"] <- "Shuffled"
    env[["Combined"]]$Basal.normalised <- env[["Combined"]]$Input.normalised 
  } 
  ## Remove data if the sample is already present in the combined dataset
  if (paste0("CPM.InputNorm.",line) %in% names(env[["Combined"]]$Input.normalised)){

    env[["Combined"]]$Input.normalised <- env[["Combined"]]$Input.normalised[, !names(env[["Combined"]]$Input.normalised) %in% paste0("CPM.InputNorm.",line)]
    env[["Combined"]]$Basal.normalised <- env[["Combined"]]$Basal.normalised[, !names(env[["Combined"]]$Basal.normalised) %in% paste0("CPM.InputNorm.",line)]
  }
  ## Add sample data to combined dataframe
  env[["Combined"]]$Input.normalised <- full_join(env[["Combined"]]$Input.normalised, 
                                                  env[[line]]$merged.input.norm[, c("synthetic.region", "CPM.InputNorm")], 
                                                  by = "synthetic.region")
  names(env[["Combined"]]$Input.normalised)[length(names(env[["Combined"]]$Input.normalised))] <- paste0("CPM.InputNorm.",line)
  env[["Combined"]]$Basal.normalised <- full_join(env[["Combined"]]$Basal.normalised, 
                                                  env[[line]]$merged.input.basal.norm[, c("synthetic.region", "CPM.Input.BasalNorm")], 
                                                  by = "synthetic.region")
  names(env[["Combined"]]$Basal.normalised)[length(names(env[["Combined"]]$Basal.normalised))] <- paste0("CPM.InputNorm.",line)
  
  
  warning("[OK]\n")
  
  # QC
  library(ggplot2)
  library(ggpubr)
  library(reshape2)
  hist(x = env[[line]]$merged$nb_barcodes.plasmid, breaks = 200)
  abline(v = env[[line]]$metadata$BC.plasmid.filtering, col = "red")
  hist(x = env[[line]]$merged$nb_barcodes.cDNA, breaks = 200)
  abline(v = min.cDNA.bcs, col = "red")
  plot3 <- ggplot(melt(env[[line]]$merged[,!names(env[[line]]$merged) %in% c("counts.plasmid","class","counts.cDNA")]), aes(x = log2(value), fill = variable)) +
    geom_density(alpha = .3) +
    ggtitle(label = paste0("CHEQ-seq ", line," - Plasmid / cDNA barcode distribution"), subtitle = "Log CPM Normalized")
  plot4 <- ggplot(env[[line]]$merged, aes(counts.plasmid, counts.cDNA)) +
    geom_point() +
    ggtitle(label = paste0(environmentName(env)," - ", line))
  plot5 <- ggplot(env[[line]]$merged, aes(log(counts.plasmid), log(counts.cDNA))) +
    geom_point() +
    ggtitle(label = paste0(environmentName(env)," - ", line," (log scale)"))
  p <- ggarrange(plot3, ggarrange(plot4, plot5,
                                  labels = c("B", "C"),
                                  ncol = 2),
                 labels = c("A"),
                 nrow = 2)
  print(p)
  
  # MA plot
  ## No normalization
  matrix <- as.matrix(env[[line]]$merged[,c("CPM.counts.plasmid","CPM.counts.cDNA")])
  ma.plot(A = rowMeans(log2(matrix)), 
          M = log2(matrix[,2])-log2(matrix[,1]), 
          cex = 1, pch = 16)
  title(main = paste0("CHEQ-seq Synthetic Enhancer - ",line), sub = "CPM Normalized MA plot")
  ## Loes normalization
  matrix_norm <- normalize.loess(matrix)
  ma.plot(A = rowMeans(log2(matrix_norm)), 
          M = log2(matrix_norm[,2])-log2(matrix_norm[,1]), 
          cex = 1, pch = 16)
  title(main = paste0("CHEQ-seq Synthetic Enhancer - ",line), sub = "Loes CPM Normalized MA plot")
  
  plot6 <- ggplot(data = env[[line]]$merged.input.norm, aes(x = class, y = log2(CPM.InputNorm))) +
    geom_boxplot() +
    # geom_hline(yintercept = log2(basal.ex), linetype="dashed", color = "red", size = 1) +
    ggtitle(label = paste0(environmentName(env), " - ", line)) +
    theme(panel.background = element_rect(fill = "white", colour = "white", size = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(size=16, face = "bold"),
          legend.key = element_rect(fill = "white"),
          legend.text = element_text(face = "bold", size = 12),
          legend.title = element_text(size = 14, face = "bold"),
          axis.text.x = element_text(angle = 60, hjust = 1),
          axis.text = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 14, face = "bold")) +
    labs(x = "Class", y = "logFC")
  print(plot6)
  
  #   labs(color = "Phenotype") +
  #   scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  return (env)
}


reg.subset <- function(tile.start, 
                       tile.end, 
                       OLS, #"A", "B"
                       MM.lines
                       ) {  
  library(stringr)
  library(dplyr)
  #Define variables
  MM.lines <- MM.lines
  chr <- str_split(tile.start, ":", simplify=T)[,1]
  start <- as.numeric(str_split(str_split(tile.start, ":", simplify=T)[,2],"-",simplify=T)[,1])
  range <- as.numeric(str_split(str_split(tile.end, ":", simplify=T)[,2],"-",simplify=T)[,1])-start
  
  #Create DF
  tmp <- data.frame("MM.line" = MM.lines, "OLS" = paste0("Tiling ",OLS), "tile.name.wt" = NA, "tile.name.mut" = NA, "wt.value" = NA, "mut.value" = NA, "wt.padj" = NA, "mut.padj" = NA)
  
  #Find max value within the defined region for wt and mut for each line
  for (i in MM.lines){
    #Subset wt and mut tiles
    if (OLS == "A"){
      wt.df <- CSE.OLSA[[i]][["merged.input.basal.norm"]][CSE.OLSA[[i]][["merged.input.basal.norm"]]$class == "wt",]
      mut.df <- CSE.OLSA[[i]][["merged.input.basal.norm"]][CSE.OLSA[[i]][["merged.input.basal.norm"]]$class == "mut",]
    }else if (OLS == "B"){
      wt.df <- CSE.OLSB[[i]][["merged.input.basal.norm"]][CSE.OLSB[[i]][["merged.input.basal.norm"]]$class == "wt",]
      mut.df <- CSE.OLSB[[i]][["merged.input.basal.norm"]][CSE.OLSB[[i]][["merged.input.basal.norm"]]$class == "mut",]
    }else{
      stop()
    }
    print(nrow(wt.df))
  #Subset by chromosome number 
  wt.df <- wt.df[grep(pattern = chr, x = wt.df$synthetic.region),]
  mut.df <- mut.df[grep(pattern = chr, x = mut.df$synthetic.region),]
  print(nrow(wt.df))

  #Extract start coordinate for each tile
  wt.df$chr.filt.wt <- as.numeric(str_split(string = str_split(string = wt.df[,"synthetic.region"],
                                                               pattern = ":", simplify=T)[,2],
                                            pattern = "-", simplify=T)[,1])
  mut.df$chr.filt.mut <- as.numeric(str_split(str_split(mut.df[,"synthetic.region"], ":", simplify=T)[,2],"-",simplify=T)[,1])

  #Substract tile start coordinate to start of region
  wt.df$substraction.start.wt <- wt.df$chr.filt.wt-start
  mut.df$substraction.start.mut <- mut.df$chr.filt.mut-start

  #Filter tiles within the range of the region and record max enhancer expression
  tmp[tmp$MM.line == i, "wt.value"] <- max(wt.df[wt.df$substraction.start.wt >= 0 & wt.df$substraction.start.wt <= range, "CPM.Input.BasalNorm"])
  tmp[tmp$MM.line == i, "wt.padj"] <- wt.df[wt.df$CPM.Input.BasalNorm == tmp[tmp$MM.line == i, "wt.value"], "padj"][1]
  tmp[tmp$MM.line == i, "tile.name.wt"] <- wt.df[wt.df$CPM.Input.BasalNorm == tmp[tmp$MM.line == i, "wt.value"], "synthetic.region"][1]
  tmp[tmp$MM.line == i, "mut.value"] <- max(mut.df[mut.df$substraction.start.mut >= 0 & mut.df$substraction.start.mut <= range, "CPM.Input.BasalNorm"])
  tmp[tmp$MM.line == i, "mut.padj"] <- mut.df[mut.df$CPM.Input.BasalNorm == tmp[tmp$MM.line == i, "mut.value"], "padj"][1]
  tmp[tmp$MM.line == i, "tile.name.mut"] <- mut.df[mut.df$CPM.Input.BasalNorm == tmp[tmp$MM.line == i, "mut.value"], "synthetic.region"][1]
  
  #Remove name suffix
  # tmp$MM.line <- str_split(string = tmp[,"MM.line"], pattern = "_", simplify = T)[,1]
  }
  return(tmp)
}