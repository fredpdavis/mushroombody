################################################################################
#
# analyzeKenyonTapinSeq_ms.R - analyze TAPIN-seq profiles of Kenyon cell types
#
# Fred P. Davis, fredpdavis@gmail.com. 2018
#
################################################################################


main <- function( data,
                  sampleFn = NA,
                  outDir = NA,
                  makeFigs = FALSE,
                  returnData = FALSE) {
# Purpose: Run all code to load data and generate figures and tables

   if (missing(data)){
      data <- list()
      data$specs <- setSpecs(sampleFn, outDir)
   }

   if (! "expr" %in% names(data)) {
      data$expr  <- loadExpr( data$specs,
                              printTxExprMat    = TRUE,
                              printTxCleanMat   = TRUE,
                              printGeneExprMat  = TRUE,
                              printGeneCleanMat = TRUE)
      data$specs$sampleLists <- makeSampleLists(data)
      if (returnData) return(data)
   }

   if (! "subtypeMarkers" %in% names(data)) {
      data$subtypeMarkers <- defineSubtypeMarkers.onemat(data)
      if (returnData) return(data)
   }

   if (makeFigs) {
      makeFigs(data = data)
   }

   return(data)

}


loadData <- function(exprOnly = FALSE) {

   startup()
   data <- list()
   data$specs             <- setSpecs()
   data$expr              <- loadExpr(specs = data$specs,
                                      printTxExprMat    = TRUE,
                                      printTxCleanMat   = TRUE,
                                      printGeneExprMat  = TRUE,
                                      printGeneCleanMat = TRUE)
   if (exprOnly) {return(data)}

   data$specs$sampleLists <- makeSampleLists(data)

   return(data)

}


startup <- function() {

   library(extrafont)
   library(pheatmap)
   library(calibrate)
   library(magrittr)
   library(dplyr)
   library(tidyr)
   library(limma)
   library(sleuth)
   library(digest)
   library(org.Mm.eg.db)
   library(png)
   library(RColorBrewer)
   library(statmod)

   font_import("/gpfs/gsfs7/users/davisfp/R/fonts")
   loadfonts()

}




makeFigs <- function(data, figList="all") {

   print("Making all tables")

   if (length(intersect(c("all","dataTable"), figList)) > 0) {
      writeDataTables(data, figName = "msDataTable")
   }

   if (length(intersect(c("all","tableS1"), figList)) > 0) {
      writeLibraryStatsTable(data,
                             tabName = "msTableS1")
   }

   if (length(intersect(c("all","geneTables"), figList)) > 0) {
      writeSubtypeMarkers(data, figName1 = "msTable2", figName2 = "msTableS2")
   }

   print("Making all figures")
   if (length(intersect(c("all","2B"), figList)) > 0) {
      compareReplicates(   data, mode="bioReps")
   }

   if (length(intersect(c("all","2C"), figList)) > 0) {
      plotHeatmap.MBsamples(data,
                            figName = "msFig2C",
                            tpmOverlay = TRUE,
                            mode = "driver",
                            pdfHeight = 2.5,
                            geneList = c("repo","elav","rut","ey","trio","Fas2","sNPF"))
   }

   if (length(intersect(c("all","3"), figList)) > 0) {
      plotHeatmap.MBsamples(data,
         figName = "msFig3A",
         mode = "driver",
         pdfHeight=4,
         geneList = unique(unlist(data$subtypeMarkers$diffGenes$cellType)))

      plotHeatmap.MBsamples(data,
         figName = "msFig3B",
         mode = "driver",
         showGeneNames = FALSE,
         pdfHeight=4,
         geneList = unique(unlist(data$subtypeMarkers$diffGenes$lobe)))

      plotHeatmap.MBsamples(data,
         figName = "msFig3B_LARGE",
         mode = "driver",
         pdfHeight = 10,
         geneList = unique(unlist(data$subtypeMarkers$diffGenes$lobe)))
   }

   if (length(intersect(c("all","4"), figList)) > 0) {
      barplotTPM(dat, geneList=c("Fas2","trio","sNPF","Gad1"),
                 figName = "msFig4")
   }

   if (length(intersect(c("all","5"), figList)) > 0) {
      tx <- plotHeatmap.MBsamples(data, figName = "msFig5A",
                                  clusterRows = TRUE,
                                  mode = "driver",
                                  geneListFn=
            paste0(data$specs$baseDir,"/data/misc/Gene.list.Mx20180212.txt"),
                                  geneListFn.pdfHeights= list("NT" = 2.5,
                                                              "NT.receptor" = 5,
                                                              "NP" = 4,
                                                              "NP.receptor" = 5,
                                                              "gap.junction" = 1.5),
                                  geneListFn.FBgnCol = 1,
                                  geneListFn.geneCol = 2,
                                  geneListFn.groupCol = 3)
   }

   return(1);

}

writeLibraryStatsTable <- function(data, tabName) {

   outTab <- data.frame(sample = character(),
                        cellType = character(),
                        num.reads = numeric(),
                        num.reads.kallisto.pseudoaligned = numeric(),
                        num.reads.star.uniq.aligned = numeric(),
                        num.detected.genes = numeric(),
                        stringsAsFactors=FALSE)

# read in all STAR results first; then more complicated kallisto parsing,
#   since have to read file to figure out which fastqSampleNum

   for (i in 1:nrow(data$specs$sampleInfo)) {
      sample_name <- data$specs$sampleInfo$sample_name[i]
      fastqSampleNum <- data$specs$sampleInfo$fastqSampleNum[i]

      starLogFn <- paste0(dat$specs$baseDir,
                          "/results/RNAseq/star_align.BDGP6.91/",
                          fastqSampleNum,
                          "/Log.final.out")

      numreads.total.star <-system(paste0("grep 'input reads' ", starLogFn),intern=TRUE)
      numreads.total.star <- gsub("[^0-9]","",  numreads.total.star)

      numreads.uniqalign <-system(paste0("grep 'reads number' ", starLogFn),
                                  intern=TRUE)
      numreads.uniqalign <- gsub("[^0-9]","",  numreads.uniqalign)

      numDetGenes <- sum(dat$expr$geneExpr[,paste0("tpm.",sample_name)] >= 1)

      outTab[sample_name, ] <- c(
                                 sample_name,
                                 data$specs$sampleInfo$cellType[i],
         prettyNum(numreads.total.star,big.mark=",",scientific=FALSE),
         NA,
         prettyNum(numreads.uniqalign,big.mark=",",scientific=FALSE),
         prettyNum(numDetGenes,big.mark=",",scientific=FALSE))
   }

   kallistoLogDir <- paste0(dat$specs$baseDir,
                            "/run/cleanrun.20180206/slurm_out") 

   kallistoLogFns <- Sys.glob(paste0(kallistoLogDir,"/process_rnaseq_samples.60986934.*.err"))

   for (kallistoLogFn in kallistoLogFns) {
      curFastqSampleNum <-system(paste0("grep 'seqtk.*_R1_001' ",kallistoLogFn,
                                        " | sed 's/.*S//' | sed 's/_.*//'"), intern=TRUE)
      if (!(curFastqSampleNum %in% data$specs$sampleInfo$fastqSampleNum)) next;

      curSample <- data$specs$sampleInfo$sample_name[
                     data$specs$sampleInfo$fastqSampleNum == curFastqSampleNum]
      numAlnReads.kallisto <-system(paste0("grep pseudoaligned ",
                                           kallistoLogFn),
                                    intern=TRUE)
      numreads.total <- gsub(" reads,.*","",numAlnReads.kallisto)
      numreads.total <- gsub(".*processed ","",numreads.total)

      numreads.pa <- gsub(".* reads, ","",numAlnReads.kallisto)
      numreads.pa <- gsub(" reads.*", "",numreads.pa)

      print(paste0("sample ", curSample," fastqsamplenum:",curFastqSampleNum," has numreads.total=",numreads.total))

      if (outTab$num.reads[outTab$sample == curSample] == numreads.total) {
         outTab$num.reads.kallisto.pseudoaligned[outTab$sample == curSample] <- numreads.pa
      } else {
         print(paste0(" -- OH FUCK READ NUMBER FUCK UP!: from stAR=",outTab$num.reads[outTab$sample==curSample]," vs from kallisto=",numreads.total))
      }

   }


   outFn <- paste0(dat$specs$outDir,"/",tabName,"_library_stats.tsv")
   write.table(outTab, file=outFn,
               col.names=TRUE,
               row.names=FALSE,
               quote=FALSE,
               sep="\t")

}


plotHeatmap.MBsamples <- function(data,
                                  mode = "driver", # "sample"
                                  pdfHeight,
                                  pdfWidth = 2.5,
                                  geneList,
                                  geneListFn,
                                  geneListFn.pdfHeights,
                                  geneListFn.FBgnCol,
                                  geneListFn.geneCol,
                                  geneListFn.groupCol,
                                  clusterCols = FALSE,
                                  clusterRows = FALSE,
                                  figName = "heatmap",
                                  tpmOverlay = FALSE,
                                  showGeneNames = TRUE) {

   hmapSamples <- c()
   hmapSampleLabels <- c()
   cellGroups <- setCellGroups()
   tx.cells <- unlist(cellGroups$cells)
   tx.labels  <- unlist(cellGroups$labels.simple)
   for (i in 1:length(tx.cells)) {
      driver <- tx.cells[i]
      label <- tx.labels[i]

      if (mode == "sample") {
         curSamples <- dat$specs$sampleInfo$cellType[dat$specs$sampleInfo$driver == driver]
      } else if (mode == "driver") {
         curSamples <- unique(dat$specs$sampleInfo$cellType[
            dat$specs$sampleInfo$driver == driver])
      }

      hmapSamples <- c(hmapSamples, curSamples)
      hmapSampleLabels <- c(hmapSampleLabels,rep(label, length(curSamples)))
   }

   tx.colGaps <- NA
   for (i in 1:length(cellGroups$cells)) {

      if (mode == "sample") {
         nSamples <- length(dat$specs$sampleInfo$sample_name[
               dat$specs$sampleInfo$driver %in% cellGroups$cells[[i]]])
      } else if (mode == "driver") {
         nSamples <- length(unique(dat$specs$sampleInfo$driver[
               dat$specs$sampleInfo$driver %in% cellGroups$cells[[i]]]))
      }

      if (is.na(tx.colGaps)) {
         tx.colGaps <- nSamples 
      } else {
         tx.colGaps <- c(tx.colGaps, max(tx.colGaps) + nSamples)
      }
   }

   if (showGeneNames) {
      show_rownames <- TRUE
   } else {
      show_rownames <- FALSE
   }

   if (!missing(geneListFn)) {
      x <- read.table(geneListFn, sep="\t", header=TRUE, stringsAsFactors=FALSE)
      geneList <- list()

      for (i in 1:nrow(x)) {
         
         if (!missing(geneListFn.FBgnCol)) {
            curGene <- unique(dat$specs$transcriptInfo$gene_name[
                              dat$specs$transcriptInfo$gene_id == x[i,geneListFn.FBgnCol]])
            if (length(curGene) < 1) {print(paste0("didn't find gene=",
                                             x[i,geneListFn.FBgnCol]))}
         } else {
            curGene <- x[i,geneListFn.geneCol]
            if (length(curGene) < 1) {print(paste0("didn't find gene=",
                                             x[i,geneListFn.geneCol]))}
         }
         if (length(curGene) < 1) {next;}

         curGroups <- unlist(strsplit(x[i,geneListFn.groupCol], ";"))
         print(paste0("cur Group = ",curGroups))
         for (curGroup in curGroups) {
            if (!curGroup %in% names(geneList)) {geneList[[curGroup]] <- c()}
            geneList[[ curGroup  ]] <- c(geneList[[curGroup]], curGene)
         }
      }

   } else {
      geneList <- list(genes = geneList)
   }

   if (mode == "sample") {
      tpmMat <- data$expr$geneExpr
   } else if (mode == "driver") {
      tpmMat <- data$expr$geneExpr.bycell
   }

   for (curListName in names(geneList)) {
         geneList[[curListName]] <- unique(geneList[[curListName]])

      if (missing(pdfHeight)) {
         pdfHeight <- 3
         if (length(geneList[[curListName]]) >= 50) {pdfHeight = 6}
      }

      if (!missing(geneListFn.pdfHeights) &
          curListName %in% names(geneListFn.pdfHeights)) {
         pdfHeight <- geneListFn.pdfHeights[[curListName]]
      }

      plotHeatmap(      data,
                        figName= paste0(figName,".",curListName),
                        pdfWidth = pdfWidth,
                        pdfHeight= pdfHeight,
                        tpmOverlay = tpmOverlay,
                        tpmMat = tpmMat,
                        plotSamples = hmapSamples,
                        sampleLabels = hmapSampleLabels,
                        show_rownames = show_rownames,
                        colGaps = tx.colGaps,
                        geneList = geneList[[curListName]],
                        clusterCols = clusterCols,
                        clusterRows = clusterRows,
                        minMaxExprThresh = 0,
                        minFc = 0)
   }
   return(1)
}

setSpecs <- function(sampleFn, outDir){
# Purpose: centralized routine to specify parameters

   specs <- list()

# SETUP PARAMETERS


# - Paths

   specs$baseDir <- "/data/davisfp/projects/dubnau.mb"
   specs$timeStamp<-format(Sys.time(), "%Y%m.%H%M")
   if (is.na(outDir)) {
      specs$outDir <- paste0(specs$baseDir,"/analysis/20180824.ms_figs")
   } else {
      specs$outDir <- outDir
   }

   if (!file.exists(specs$outDir)) dir.create(specs$outDir, recursive = TRUE)

# THRESHOLDS
   specs$thresh$exprGenes.minTPM <- 10
   specs$thresh$deGenes.FC <- 1.5
   specs$thresh$deGenes.FC.qnl <- 4
   specs$thresh$sleuth.qval <- 0.05
   specs$thresh$voom.qval <- 0.05

# EXPR data files
   specs$kallistoBaseDir <- paste0(specs$baseDir,
      "/results/RNAseq/kallisto.BDGP6.91/")

   specs$sleuthBaseDir <- paste0(specs$baseDir, "/results/RNAseq/sleuth.BDGP6.91/")


# - Debug

   specs$debug                 <- TRUE
   specs$debugSubset           <- FALSE  # boolean to only process some genes
   if (specs$debug)            print("WARNING! RUNNING IN DEBUG MODE")


# LOAD GLOBALLY ACCESSIBLE INFORMATION

   if (is.na(sampleFn)) {
      specs$sampleFn <- paste0(specs$baseDir,
         "/meta/kenyoncell_ms_rnaseq_samples.txt")
   } else {
      specs$sampleFn <- sampleFn
   }
   specs$figSamples <- list(
      skip.samples      = c( ))

   print("LOADING sample info")
   specs$sampleInfo <- loadSamples(specs)

   if (is.na(sampleFn)) {
   specs$sampleInfo$sample_name <- gsub("air_","rep", specs$sampleInfo$sample_name)
   }

# - Define 'celltypes' that are combination, control, or single drivers
   specs$allDrivers     <- unique(specs$sampleInfo$driver)

   specs$controlDrivers <- 
      specs$sampleInfo$driver[specs$sampleInfo$driver == "mock"]

# Load file: Transcript info
   specs$transcriptInfoFn <- paste0(specs$baseDir,
      "/data/misc_files.BDGP6.91/BDGP6.91.ERCC.INTACT.transcript_info.txt")
   specs$transcriptInfo <- read.table( specs$transcriptInfoFn,
                                       quote     = "",
                                       header    = TRUE,
                                       sep       = "\t",
                                       as.is     = TRUE )

# Load file: ERCC concentrations
   specs$erccConcFn <- paste0(specs$baseDir,
                              "/data/external/ERCC/cms_095046.txt")

# Load files: AnimalTFDB
   specs$animalTFDB.fn <- list(
      chromremodel = paste0(specs$baseDir,"/data/external/AnimalTFDB/",
          "Drosophila_melanogaster_chromatin_remodeling_factors_gene_list.txt"),
      cofactor = paste0(specs$baseDir,"/data/external/AnimalTFDB/",
          "Drosophila_melanogaster_transcription_co-factors_gene_list.txt"),
      tf  = paste0(specs$baseDir,"/data/external/AnimalTFDB/",
          "Drosophila_melanogaster_transcription_factors_gene_list.txt")
   )

   return(specs)
}


loadSamples <- function( specs ){
# Purpose: Load information about optic lobe RNA-seq samples

   print("Load sample information")
   sampleInfo <- read.table( specs$sampleFn,
                             header       = TRUE,
                             sep          = "\t",
                             as.is        = TRUE,
                             comment.char = '#')

   sampleInfo$driver <- sapply(strsplit(sampleInfo$sampleName,"_"),"[",1)

   sampleInfo$sample_name <- sampleInfo$sampleName
   sampleInfo$sampleName <- NULL

   sampleInfo$yield.nuclei <- sampleInfo$yield.nuclei / 1E3
   sampleInfo$yield.cdna <- sampleInfo$ND_ng_uL

   print("set yield.cdna = ND_ng_uL")


   return(sampleInfo)
}


loadExpr <- function( specs,
                      printTxExprMat  = FALSE,
                      printTxCleanMat = FALSE,
                      printGeneExprMat  = FALSE,
                      printGeneCleanMat  = FALSE,
                      plotTxLenDistr  = FALSE,
                      loadExons       = FALSE) {

# Purpose: load KALLISTO tables of transcript abundance,
#          generates merged tables of transcript and gene estimates,
#          generates tables of exon coverage counts
#          plots rRNA content vs yield

   library(plyr)

# Load BDGP transcript information

   print("Loading transcript info")
   transcriptInfo <- specs$transcriptInfo

# Build gene_id <-> gene_name map
   geneInfo <- unique(transcriptInfo[,c("gene_id","gene_name")])

   nSamples <- nrow(specs$sampleInfo)
   print(paste("Will process ",nSamples, "samples"))

   txExprMat <- NULL
   exonMat <- NULL
   intronMat <- NULL
   for (i in 1:nSamples){
      print(paste("Reading expression from sample # ", i,
                  ": ", specs$sampleInfo$sample_name[i], " ", date()))

#target_id       length  eff_length      est_counts      tpm
#FBtr0082757     1844    1595    9       0.718348

      abundFn <- paste0( specs$kallistoBaseDir,
                            specs$sampleInfo$fastqSampleNum[i],
                            "/abundance.tsv" )

      curAbund <- read.table(abundFn,
                             header     = TRUE,
                             sep        = "\t",
                             colClasses = c("character", "numeric",
                                            "NULL", "numeric", "numeric"))

      colnames(curAbund) <- c( "transcript_id",
                               "length",
                               paste0("est_counts.",
                                      specs$sampleInfo$sample_name[i]),
                               paste0("tpm.", specs$sampleInfo$sample_name[i]))

      if (is.null(txExprMat)) {
         txExprMat <- merge( transcriptInfo,
                             curAbund,
                             all.y      = TRUE,
                             by         = "transcript_id" )
      } else {
         txExprMat$length <- NULL
         txExprMat <- merge( txExprMat,
                             curAbund,
                             by = c("transcript_id"))
      }

#      2L      7528    8116    FBgn0031208     11
#      2L      8192    9484    FBgn0031208     5

      exonFn<-paste0(specs$baseDir, "/results/RNAseq/star_align/",
                specs$sampleInfo$sample_name[i], "/",
                specs$sampleInfo$sample_name[i],".exon_coverage_pergene.txt")

      if (file.exists(exonFn)) {
         if (specs$debug) print(paste("-> reading exon coverage: ", date()))

         curExon <- read.table(exonFn,
                               header     = FALSE,
                               sep        = "\t",
                               colClasses = c("character", "numeric",
                                              "numeric", "character", "numeric"))

         colnames(curExon) <- c("chr", "start", "end", "gene_id",
                                paste0("exon_counts.",
                                       specs$sampleInfo$sample_name[i]))

         if (is.null(exonMat)) {
            exonMat <- curExon
         } else {
            curExon$chr <- NULL
            curExon$start <- NULL
            curExon$end <- NULL
            curExon$gene_id <- NULL
            exonMat <- cbind(exonMat,curExon)
         }
      }


      intronFn <- paste0(specs$baseDir,
                "/results/RNAseq/star_align/",
                specs$sampleInfo$sample_name[i], "/",
                specs$sampleInfo$sample_name[i],".intron_coverage_pergene.txt")
      if (file.exists(intronFn)){
         if (specs$debug) print(paste("-> reading intron coverage: ", date()))

         curIntron <- read.table(intronFn,
            header=FALSE,sep="\t",
            colClasses=c("character","numeric","numeric","character","numeric")
         )
         colnames(curIntron) <- c("chr", "start", "end", "gene_id",
                                  paste0("intron_counts.",
                                         specs$sampleInfo$sample_name[i]))


         if (is.null(intronMat)) {
            intronMat <- curIntron
         } else {
            curIntron$chr          <- NULL
            curIntron$start        <- NULL
            curIntron$end          <- NULL
            curIntron$gene_id      <- NULL
            intronMat <- cbind(intronMat, curIntron)
         }
      }
      if (specs$debug) print(paste("    -> DONE: ",date()))
   }

   debugSubsetGenes <- NULL

   if (loadExons) {
   print("Reading intron/exon coverage tables")
   samplesWithEI <- gsub("intron_counts.", "",
         colnames(intronMat)[grep("intron_counts.", colnames(intronMat))])

   exonMat <- merge( geneInfo,
                     exonMat,
                     all.y      = TRUE,
                     by         = c("gene_id"))

   intronMat <- merge( geneInfo,
                       intronMat,
                       all.y = TRUE,
                       by = c("gene_id"))

   exonMat$exon.length <- c(exonMat$end - exonMat$start)
   intronMat$intron.length <- c(intronMat$end - intronMat$start)

   exonMat$chr          <- NULL
   exonMat$start        <- NULL
   exonMat$end          <- NULL
   intronMat$chr        <- NULL
   intronMat$start      <- NULL
   intronMat$end        <- NULL


# Convert from exon/intron-levelcounts to per gene intronic and exonic RPKM

   flatIntronMat <- ddply(
      intronMat, c("gene_id","gene_name"),
      function(a){
         apply(a[,
            c( colnames(intronMat)[grep("intron_counts", colnames(intronMat))],
               "intron.length")], 2, sum)
      }
   )

   for (i in 1:length(samplesWithEI)) {
      intronCol <- paste0("intron_counts.", samplesWithEI[i])
      intronFPKM <- paste0("intron_fpkm.", samplesWithEI[i])

      flatIntronMat[,intronFPKM] <- c((1E9 * flatIntronMat[, intronCol])/ 
                                      (flatIntronMat$intron.length *
                                       sum(flatIntronMat[, intronCol])))
   }


   flatExonMat <- ddply(
      exonMat, c("gene_id", "gene_name"),
      function(a){
         apply(a[,
            c( colnames(exonMat)[grep("exon_counts", colnames(exonMat))],
               "exon.length")], 2, sum)
      }
   )

   for (i in 1:length(samplesWithEI)) {
      exonCol <- paste0("exon_counts.", samplesWithEI[i])
      exonFPKM <- paste0("exon_fpkm.", samplesWithEI[i])

      flatExonMat[,exonFPKM] <- c((1E9 * flatExonMat[,exonCol])/ 
                                (flatExonMat$exon.length *
                                 sum(flatExonMat[,exonCol])))
   }


   covMat <- merge(flatIntronMat, flatExonMat,
                  all=TRUE,
                  by=c("gene_id","gene_name"))

   covMat <- na.omit(covMat)
   for (i in 1:length(samplesWithEI)) {
      intronFPKM <- paste0("intron_fpkm.", samplesWithEI[i])
      exonFPKM <- paste0("exon_fpkm.", samplesWithEI[i])
      intronCounts <- paste0("intron_counts.", samplesWithEI[i])
      relIntronFPKM <- paste0("rel_intron_fpkm.", samplesWithEI[i])

      covMat[,relIntronFPKM] <- c( covMat[,intronCounts] /
                                   covMat$intron.length) / 
                                c( sum(covMat[,intronCounts]) /
                                   sum(covMat$intron.length))

      exonCounts <- paste0("exon_counts.", samplesWithEI[i])
      relExonFPKM <- paste0("rel_exon_fpkm.", samplesWithEI[i])
      covMat[,relExonFPKM] <- c( covMat[,exonCounts] /
                                 covMat$exon.length) / 
                              c( sum(covMat[,exonCounts]) /
                                 sum(covMat$exon.length))

      ratioCol <- paste0("intron_exon_fpkm_ratio.", samplesWithEI[i])
      covMat[,ratioCol] <- c(covMat[,intronFPKM] / covMat[,exonFPKM])

      relRatioCol <- paste0("rel_intron_exon_fpkm_ratio.", samplesWithEI[i])
      covMat[,relRatioCol] <- c(covMat[,relIntronFPKM] / covMat[,relExonFPKM])
   }
   print("-> DONE")
   } else {
      covMat<-NA
   }



## Prepare transcript expression matrices for output

   allCols <- colnames(txExprMat)
   tpmCols <- allCols[grep("tpm.",allCols)]
   estCountCols <- allCols[grep("est_counts.",allCols)]

   # Reorder transcript matrix columns
   newTxExprMat <- txExprMat[, c("transcript_id", "transcript_name", "length",
                               "gene_id", "gene_name", "gene_biotype",
                               tpmCols, estCountCols)]

   # Format to max 3 decimal digits
   if (printTxExprMat) {
      x <- newTxExprMat
      x$length <- NULL
      for (i in c(tpmCols, estCountCols)) {
         x[,i] <- sprintf("%.3f", x[,i])}

      gz1 <- gzfile(paste0(specs$outDir,
                           "/transcriptExpressionMatrix.txt.gz"),"w")
      write.table(x,
                  file          = gz1,
                  col.names     = TRUE,
                  row.names     = FALSE,
                  quote         = FALSE,
                  sep           = "\t")
      close(gz1)
      x <- NULL
   }


##   biotype content
   biotypeExpr <- aggregate( txExprMat[,c(tpmCols)],
      by=list(gene_biotype = txExprMat[["gene_biotype"]]),
      FUN = sum )
   rownames(biotypeExpr)<-biotypeExpr$gene_biotype
   biotypeExpr$gene_biotype<-NULL


   # Remove ERCC, INTACT, rRNA entries; renormalize TPM to 1M total
   cleanTxExprMat <- newTxExprMat[newTxExprMat$gene_biotype %in%
                                  c("protein_coding"),]

#   cleanTxExprMat <- newTxExprMat[grep("ERCC|INTACT|rRNA",
#                                       newTxExprMat$gene_name, invert=TRUE),]

   print("Renormalizing expression table after removing ERCC, INTACT, rRNA")
   for (i in (1:length(tpmCols))) {
      cleanTxExprMat[, tpmCols[i]] <- (1E6 * cleanTxExprMat[, tpmCols[i]] / 
                                      sum(cleanTxExprMat[, tpmCols[i]]))
   }
   print("-> DONE")


   #Plot transcript length distribution
   if (plotTxLenDistr) {        # not sure if necessary, redundant with plotTxLengthDistr?
      lineColors <- topo.colors(length(tpmCols))

      lineColors[grep("T1", tpmCols)] <- "red"

      pdf(paste0(specs$outDir,
                 "/global_transcriptome_length_distribution.pdf"))
      lengthDistrs <- list()
      lengthDistrMaxY <- 0
      lengthDistrMaxX <- 0
      for (i in (1:length(tpmCols))) {
         lengthDistrs[[tpmCols[i]]] <- density(
            rep(log10(cleanTxExprMat$length),
                floor(cleanTxExprMat[,tpmCols[i]])))

         lengthDistrMaxY <- max(lengthDistrMaxY,
                                max(lengthDistrs[[tpmCols[i]]]$y))

         lengthDistrMaxX <- max(lengthDistrMaxX,
                                max(lengthDistrs[[tpmCols[i]]]$x))
      }

      plot( lengthDistrs[[tpmCols[1]]],
            main        = "",
            xlab        = "length (log10 nt)",
            ylab        = "Density",
            lwd         = 2,
            col         = lineColors[1] )

      for (i in (2:length(tpmCols))) {
         lines(lengthDistrs[[tpmCols[i]]], lwd=2, col=lineColors[i])
      }
      legend("topright",
             tpmCols,
             col        = lineColors,
             cex        = 0.4,
             text.col   = lineColors)
      dev.off()
   }

   # Format to max 3 decimal digits
   if (printTxCleanMat) {
      x <- cleanTxExprMat
      x$length <- NULL
      for (i in c(tpmCols, estCountCols)) x[, i] <- sprintf("%.3f", x[, i])

      gz1 <- gzfile(paste0(specs$outDir,
                "/transcriptExpressionMatrix.no_RRNA_ERCC_INTACT.txt.gz"), "w")
      write.table(x,
                  file          = gz1,
                  col.names     = TRUE,
                  row.names     = FALSE,
                  quote         = FALSE,
                  sep           = "\t")
      close(gz1)
      x <- NULL
   }

   print("Calculating gene expression matrix")

   geneExprMat <- aggregate( txExprMat[, c(tpmCols, estCountCols)],
                             by = list(gene_id    = txExprMat[["gene_id"]],
                                       gene_name  = txExprMat[["gene_name"]]),
                             FUN = sum )

   # Format to max 3 decimal digits
   if (printGeneExprMat) {
      x <- geneExprMat
      for (i in c(tpmCols, estCountCols)) x[,i] <- sprintf("%.3f", x[,i])

      gz1 <- gzfile(paste0(specs$outDir, "/geneExpressionMatrix.txt.gz"), "w")
      write.table(x,
                  file          = gz1,
                  col.names     = TRUE,
                  row.names     = FALSE,
                  quote         = FALSE,
                  sep           = "\t")
      close(gz1)
      x <- NULL
   }


   cleanGeneExprMat <- aggregate( cleanTxExprMat[,c(tpmCols, estCountCols)],
                           by = list(gene_id   = cleanTxExprMat[["gene_id"]],
                                     gene_name = cleanTxExprMat[["gene_name"]]),
                           FUN = sum )

   # Format to max 3 decimal digits
   if (printGeneCleanMat) {
      x <- cleanGeneExprMat
      for (i in c(tpmCols, estCountCols)) {
         x[,i] <- sprintf("%.3f", x[,i]) }

      gz1 <- gzfile(paste0(specs$outDir,
                     "/geneExpressionMatrix.no_RRNA_ERCC_INTACT.txt.gz"),"w")
      write.table(x,
                  file          = gz1,
                  col.names     = TRUE,
                  row.names     = FALSE,
                  quote         = FALSE,
                  sep           = "\t")
      close(gz1)
      x <- NULL
   }


   exprMat.forDE <- cleanGeneExprMat[, c("gene_id", estCountCols)]
   rownames(exprMat.forDE) <- exprMat.forDE$gene_id
   exprMat.forDE$gene_id <- NULL


   print("Setup transcript expression table by cell type")
   cleanTxExprMat.bycell <- cleanTxExprMat[, c("transcript_id", "gene_id",
                                               "gene_name", "length")]
   for (celltype in sort(specs$sampleInfo$celltype)){

         curSamples <- specs$sampleInfo$sample_name[
                           specs$sampleInfo$celltype == celltype]

         curTpmCols <- paste0("tpm.", curSamples)
         newTpmCol <- paste0("tpm.", celltype)
         if (length(curTpmCols) == 1) {
            curAvgTPM <- cleanTxExprMat[, curTpmCols[1]]
         } else {
            curAvgTPM <- apply(cleanTxExprMat[, curTpmCols], 1, mean)
         }
         cleanTxExprMat.bycell[, newTpmCol] <- curAvgTPM

   }


   print("Setup gene expression table by cell type")
# Set up bycell gene expression matrix
   cleanGeneExprMat.bycell <- cleanGeneExprMat[,c("gene_id","gene_name")]
   rownames(cleanGeneExprMat.bycell) <- cleanGeneExprMat.bycell$gene_name

   for (cellType in sort(specs$sampleInfo$cellType)){

      curSamples <- specs$sampleInfo$sample_name[
                        specs$sampleInfo$cellType == cellType ]

      curTpmCols <- paste0("tpm.", curSamples)
      newTpmCol <- paste0("tpm.", cellType)

      if (length(curTpmCols) == 1) {
         curAvgTPM <- cleanGeneExprMat[, curTpmCols[1]]
      } else {
         curAvgTPM <- apply(cleanGeneExprMat[, curTpmCols], 1, mean)
      }
      cleanGeneExprMat.bycell[, newTpmCol] <- curAvgTPM

   }

   print("Setup gene expression table by lobe")
# Set up bylobe gene expression matrix
   cleanGeneExprMat.bylobe <- cleanGeneExprMat[,c("gene_id","gene_name")]
   rownames(cleanGeneExprMat.bylobe) <- cleanGeneExprMat.bylobe$gene_name

   for (lobe in sort(specs$sampleInfo$lobe)){

      curSamples <- specs$sampleInfo$sample_name[
                        specs$sampleInfo$lobe == lobe ]

      curTpmCols <- paste0("tpm.", curSamples)
      newTpmCol <- paste0("tpm.", lobe)

      if (length(curTpmCols) == 1) {
         curAvgTPM <- cleanGeneExprMat[, curTpmCols[1]]
      } else {
         curAvgTPM <- apply(cleanGeneExprMat[, curTpmCols], 1, mean)
      }
      cleanGeneExprMat.bylobe[, newTpmCol] <- curAvgTPM

   }


   return(list(
      biotypeExpr       = biotypeExpr,
      geneExpr          = cleanGeneExprMat,
      geneExpr.withERCC = geneExprMat,
      geneExpr.bycell   = cleanGeneExprMat.bycell,
      geneExpr.bylobe   = cleanGeneExprMat.bylobe,
      txExpr            = cleanTxExprMat,
      txExpr.bycell     = cleanTxExprMat.bycell,

      exprCounts        = exprMat.forDE,
      transcriptInfo    = transcriptInfo,
      geneInfo          = geneInfo,
      exonIntronMat     = covMat
   ))

}




compareReplicates <- function( data,
                               plotRepScatter          = TRUE,
                               returnDatOnly           = FALSE,
                               mode                    = "bioReps",
                               selectFrom ){
# Purpose: make replicate scatterplots and quantify correlation

   library(png)

   allSamples <- data$specs$sampleInfo

   if (!missing(selectFrom)) {
         allSamples <- allSamples[allSamples$sample_name %in% selectFrom,] }

   tpmRange <- range(1 + dat$expr$geneExpr[,paste0("tpm.",allSamples$sample_name)])

   bioRepCorMat.sample1<-c()
   bioRepCorMat.sample2<-c()
   bioRepCorMat.cor<-c()
   for (driver in c(unique(allSamples$driver))){
      print(paste0("Computing biological replicate correlation for: ", driver))
      curSamples <- allSamples[allSamples$driver == driver, ]
      curRepR <- c()
      if (nrow(curSamples) == 1) next;
      for (i in (1:(nrow(curSamples) - 1))) {

         colX <- paste0("tpm.", curSamples$sample_name[i])
         curX <- data$expr$geneExpr[, colX]

         for (j in ((i + 1):nrow(curSamples))) {

            colY <- paste0("tpm.", curSamples$sample_name[j])
            curY <- data$expr$geneExpr[, colY]

            curR <- cor(log2(1 + curX), log2(1 + curY))
            curRepR <- c(curRepR, curR)

            bioRepCorMat.sample1 <- c(bioRepCorMat.sample1,
                                      curSamples$sample_name[i])
            bioRepCorMat.sample2 <- c(bioRepCorMat.sample2,
                                      curSamples$sample_name[j])
            bioRepCorMat.cor <- c(bioRepCorMat.cor, curR)

            if (plotRepScatter) {
               curX <- data$expr$geneExpr[, colX]
               curY <- data$expr$geneExpr[, colY]

   tmppngfn <- tempfile()
   png(file = tmppngfn, height = 3.1, width = 3.1, units = "in", res = 300)
   par(mar = c(0, 0,0, 0))
               plot(1 + curX,
                    1 + curY,
                    log = "xy",
                    ann=FALSE,axes=FALSE,main="",
                    las = 1,
                    xlab= paste0(curSamples$sample_name[i], " (TPM + 1)"),
                    ylab= paste0(curSamples$sample_name[j], " (TPM + 1)"),
                    xlim = tpmRange,
                    ylim = tpmRange,
                    pch = 20,
                    cex = 0.5)
   dev.off()
   pngbg <- readPNG(tmppngfn)
   pngbg <- as.raster(pngbg)

               pdf(paste0(data$specs$outDir, "/scatter_bio_replicate_",
                          driver, "_",
                          curSamples$sample_name[i], "_vs_",
                          curSamples$sample_name[j], ".pdf"),
                   height=3.5, width=3.5)

               par(mar = c(3.75, 3.75, 0.5, 0.5),
                   mgp = c(2, 0.6, 0), 
                   cex=1.3,cex.axis=0.8, cex.lab=1)

               par(cex=1.3)

               plot(1 + curX,
                    1 + curY,
                    log = "xy",
                    main= "",
                    type="n",
                    xlim = tpmRange,
                    ylim = tpmRange,
                    xlab= paste0(curSamples$sample_name[i], " (TPM + 1)"),
                    ylab= paste0(curSamples$sample_name[j], " (TPM + 1)"),
                    pch = 20,
                    cex = 0.5)
               lim <- par()
   rasterImage(pngbg, 10^lim$usr[1], 10^lim$usr[3],
                      10^lim$usr[2], 10^lim$usr[4])


               abline(a=0,b=1,lwd=2)

               legend("topleft",
                      legend= paste0("pearson r=", sprintf("%.2f", curR)),
                      bty   = "n")

               dev.off()
            }
         }
      }

   }
   range.repR <- range(bioRepCorMat.cor)
   print(paste("Biological Replicate r range: ",
               paste0(range.repR, collapse=" - ")))

   bioRepCorMat<-data.frame( rep1 = bioRepCorMat.sample1,
                             rep2 = bioRepCorMat.sample2,
                             cor  = bioRepCorMat.cor,
                             stringsAsFactors = FALSE)
   if (returnDatOnly) {return(bioRepCorMat)}

   outFn<-paste0(data$specs$outDir,
                    "/table_biological_replicate_correlation.txt")
   write.table(bioRepCorMat,
               outFn,
               col.names = TRUE,
               row.names = FALSE,
               quote     = FALSE,
               sep       = "\t")

   pdf(paste0(data$specs$outDir,
              "/hist_biological_replicate_correlations.pdf"))
   par(cex=1.5, mar=c(5, 5, 1, 1))
   hist(bioRepCorMat.cor,
        lwd  = 2,
        xlab = "Biological replicate correlation (pearson)",
        main = "")
   dev.off()

   return(TRUE)
}


makeSampleLists <- function(data) {
# purpose: make lists of samples to use for most analysis and inferState() analysis

   sampleList <- list()

   sampleList$all <- data$specs$sampleInfo$sample_name

   if ("figSamples" %in% names(data$specs)) {

   if ("skip.samples" %in% names(data$specs$figSamples)) {
      sampleList$all <- setdiff(sampleList$all,
                                data$specs$figSamples$skip.samples) }

   }

   sampleList$good <- sampleList$all[!sampleList$all %in%
                                     data$specs$deg.skipSamples]

   return(sampleList)

}



plotHeatmap <- function( data,
                         tpmMat,
                         tpmOverlay = FALSE,
                         minMaxExprThresh = 30,
                         clusterRows = FALSE,
                         clusterCols = FALSE,
                         show_rownames = TRUE,
                         colGaps = FALSE,
                         plotSamples,
                         sampleLabels,
                         minFc = 4,
                         pdfHeight = 4,
                         pdfWidth = 5,
                         geneList=NULL, #optional list of genes
                         figName ){

# If no gene_id column, make one, assume rownames are names

   library(pheatmap)
   library(RColorBrewer)

   if (missing(geneList)) { geneList <- tpmMat$gene_name }

   if (missing(plotSamples)) {
      tpmCols <- paste0("tpm.", data$specs$sampleInfo$sample_name)
   } else {
      tpmCols <- paste0("tpm.", plotSamples)
   }
   tpmMat <- tpmMat[match(geneList, tpmMat$gene_name), tpmCols]


   rownames(tpmMat) <- geneList

   tpmMat <- tpmMat[ apply(tpmMat, 1, max) >= minMaxExprThresh, ]
   print(paste0("   all genes (TPM > ", minMaxExprThresh, "): ", nrow(tpmMat)))
   colnames(tpmMat) <- gsub("tpm.", "", colnames(tpmMat))

   if (!missing(sampleLabels)) {
      colLabels <- sampleLabels
   } else {
      colLabels <- colnames(tpmMat)
   }

   curFn <- paste0("heatmap_tpm.pdf")

   curMat <- log10(1 + tpmMat)
   curMat[curMat > 3] <- 3
   curBreaks <- seq(0, 3, length.out=101)
   colorScaleFull <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)

   if (!missing(figName)) curFn <- paste0(figName, ".", curFn)

   print(paste0("NOW PHEATMAPPING n=",nrow(curMat)," rows with cluster_rows=",clusterRows))
   if (nrow(curMat) <= 2) {clusterRows <- FALSE}
   curMat <- na.omit(curMat)

   displayNumbers <- FALSE
   if (tpmOverlay) {
      displayNumbers <- round(tpmMat[rownames(curMat),], digits=0)
   }


   pdf(paste0(data$specs$outDir, "/", curFn), onefile=FALSE,
       height=pdfHeight, width=pdfWidth)
   pheatmap(curMat,
            scale        = "none",
            treeheight_row = 0,
            border_color = NA,
            cluster_cols = clusterCols,
            cluster_rows = clusterRows,
            show_rownames = show_rownames,
            labels_col = colLabels,
            gaps_col = colGaps,
            breaks = curBreaks,
            display_numbers = displayNumbers,
            number_color = "black",
            legend_breaks = c(0,1,2,3),
            legend_labels = c(0,10,100,1000),
            col             = colorScaleFull,
#            main            = "Gene expression (TPM)",
            height       = pdfHeight,
            width        = pdfWidth, 
            fontsize_col     = 6,
             fontsize_number=getHmapFontsize(nrow(curMat)) - 1,
             fontsize_row  = getHmapFontsize(nrow(curMat)))
   dev.off()


   tpmMat <- log2((1 + tpmMat) / (1 + apply(tpmMat, 1, mean)))

   tpmMat[tpmMat > 2.5] <- 2.5
   tpmMat[tpmMat < -2.5] <- -2.5
   curBreaks <- seq(-2.5, 2.5, length.out=101)

   curFn <- paste0("heatmap_relexpr.pdf")

   colorScaleFull <- colorRampPalette(c("steelBlue3", "white", "darkOrange3"))(100)

   if (!missing(figName)) curFn <- paste0(figName, "_", curFn)
   pdf(paste0(data$specs$outDir, "/",  curFn),
       height=pdfHeight, width=pdfWidth,
       onefile=FALSE)

   tpmMat <- na.omit(tpmMat)
   pheatmap( tpmMat,
             scale           = "none",
             treeheight_row = 0,
#             main            = "relative gene expression (log2 TPM/mean)",
             main            = "",
             border_color    = NA,
             cluster_cols = clusterCols,
             cluster_rows = clusterRows,
             show_rownames = show_rownames,
             gaps_col = colGaps,
             labels_col = colLabels,
             height          = pdfHeight,
             width           = pdfWidth,
             col             = colorScaleFull,
             breaks          = curBreaks,
             fontsize_col  = 3.5,
             fontsize_row  = getHmapFontsize(nrow(tpmMat)))
   dev.off()

}


getHmapFontsize <- function (n) {
# Given number of rows, returns font for pheatmap

      if (n < 10) {
         return(6)
      } else if (n < 100) {
         return(3)
      } else if (n < 200) {
         return(2.5)
      } else if (n < 400) {
         return(2)
      } else {
         return(1.5)
      }

}



writeDataTables <- function(data, figName = "msDataTable") {

# Table 1: Transcript expression in all samples dat$expr$txExpr
# Table 2: Gene expression in all samples dat$expr$geneExpr
# Table 3: Gene expression - celltype average dat$expr$geneExpr.bycell
# Table 4: Gene expression - lobe average dat$expr$geneExpr.bylobe

   
   readmeFn <- paste0(dat$specs$outDir,"/",figName,"_README.txt")
   
   outFn1 <- paste0(dat$specs$outDir,"/",figName,
                     "_table1_transcript_x_sample.txt.gz")
   outFn2 <- paste0(dat$specs$outDir,"/",figName,
                     "_table2_gene_x_sample.txt.gz")
   outFn3 <- paste0(dat$specs$outDir,"/",figName,
                     "_table3_gene_x_subtype.txt.gz")
   outFn4 <- paste0(dat$specs$outDir,"/",figName,
                     "_table4_gene_x_class.txt.gz")


   print(paste0("Writing data table 1: ", basename(outFn1)))
   outTab <- data$expr$txExpr[,paste0("tpm.",sort(data$specs$sampleInfo$sample_name))]
   outTab <- round(outTab, digits=2)
   outTab <- cbind(dat$expr$txExpr[,c("transcript_id", "transcript_name",
                                      "gene_id", "gene_name")], outTab)
   colnames(outTab) <- gsub("tpm.", "", colnames(outTab))
   gz1 <- gzfile(outFn1, "w")
   write.table(outTab, file=gz1,
               col.names=TRUE,
               row.names=FALSE,
               quote=FALSE, sep="\t")
   close(gz1)


   print(paste0("Writing data table 2: ", basename(outFn2)))
   outTab <- data$expr$geneExpr[,paste0("tpm.",sort(data$specs$sampleInfo$sample_name))]
   outTab <- round(outTab, digits=2)
   outTab <- cbind(dat$expr$geneExpr[,c("gene_id", "gene_name")], outTab)
   colnames(outTab) <- gsub("tpm.", "", colnames(outTab))
   gz1 <- gzfile(outFn2, "w")
   write.table(outTab, file=gz1,
               col.names=TRUE,
               row.names=FALSE,
               quote=FALSE, sep="\t")
   close(gz1)


   print(paste0("Writing data table 3: ", basename(outFn3)))
   outTab <- data$expr$geneExpr.bycell
   tpmCols <- colnames(outTab); tpmCols <- tpmCols[grepl("tpm.",tpmCols)]
   outTab[,tpmCols] <- round(outTab[,tpmCols], digits=2)
   gz1 <- gzfile(outFn3, "w")
   write.table(outTab, file=gz1,
               col.names=TRUE,
               row.names=FALSE,
               quote=FALSE, sep="\t")
   close(gz1)


   print(paste0("Writing data table 4: ", basename(outFn4)))
   outTab <- data$expr$geneExpr.bylobe
   tpmCols <- colnames(outTab); tpmCols <- tpmCols[grepl("tpm.",tpmCols)]
   outTab[,tpmCols] <- round(outTab[,tpmCols], digits=2)
   gz1 <- gzfile(outFn4, "w")
   write.table(outTab, file=gz1,
               col.names=TRUE,
               row.names=FALSE,
               quote=FALSE, sep="\t")
   close(gz1)

}

defineSubtypeMarkers.onemat <- function(dat){
# Define KC-subtype marker genes

   library(pander)
   library(limma)
   library(pheatmap)

   genes.tmg <- list()

   allSamples <- dat$specs$sampleInfo
   allSamples <- allSamples[!allSamples$sample_name %in%
                            dat$specs$deg.skipSamples,]
   allSamples <- allSamples[allSamples$driver != "mock",]
   allSamples <- allSamples[allSamples$driver != "mock",]

   curMat.counts <- dat$expr$geneExpr[, paste0("est_counts.", allSamples$sample_name)]
   rownames(curMat.counts) <- dat$expr$geneExpr$gene_name
   detGenes <- rownames(curMat.counts)[which(apply(curMat.counts,1,function(x) { sum(x > 1)}) > ncol(curMat.counts) / 2)]
   curMat.counts <- curMat.counts[detGenes,]
   print(paste0("Evaluating ", length(detGenes)," genes detected in all samples"))

   tpmMat <- dat$expr$geneExpr[, paste0("tpm.", allSamples$sample_name)]
   tpmMat <- round(tpmMat,digits=3)
   tpmMat$gene_name <- dat$expr$geneExpr$gene_name
   rownames(tpmMat) <- tpmMat$gene_name

# sort by lobe and cell type
   tpmCols <- c()
   tpmCol.feats <- list(cellType = c(), lobe = c())
   tpmCol.gaps <- c()
   for (lobe in sort(unique(allSamples$lobe))) {
      for (cellType in sort(unique(allSamples$cellType[allSamples$lobe == lobe]))) {
         curSamples <- allSamples$sample_name[allSamples$cellType == cellType]
         tpmCols <- c(tpmCols, paste0("tpm.",curSamples))
         tpmCol.gaps <- c(tpmCol.gaps, length(tpmCols))

         tpmCol.feats$cellType <- c(tpmCol.feats$cellType, rep(cellType,length(curSamples)))
         tpmCol.feats$lobe <- c(tpmCol.feats$lobe, rep(lobe,length(curSamples)))
      }
      tpmCol.gaps <- c(tpmCol.gaps, length(tpmCols))
   }

   tpmMat <- tpmMat[,tpmCols]
   tpmCol.annot <- data.frame(cellType = tpmCol.feats$cellType,
                              lobe = tpmCol.feats$lobe)
   rownames(tpmCol.annot) <- colnames(tpmMat)
   tpmCol.annotColors <- list(
      lobe      = c(alpha.beta          = "darkslategray4",
                    alphap.betap        = "palegreen4",
                    gamma               = "darkOrchid4"),
      cellType  = c(alpha.beta.c        = "darkslategray1",
                    alpha.beta.p        = "darkslategray2",
                    alpha.beta.s        = "darkslategray3",
                    alphap.betap.ap     = "palegreen1",
                    alphap.betap.m      = "palegreen2",
                    gamma.d             = "darkOrchid1",
                    gamma.m             = "darkOrchid2")
                              )

   diffGenes <- list()

   for (groupType in c("cellType","lobe")) {
      diffGenes[[groupType]] <- list()
      print(paste0("Finding ",groupType,"-specific genes"))

      ct <- factor(allSamples[,groupType])

      design <- model.matrix(~0 + ct)
      colnames(design) <- levels(ct)

      v <- voom(counts = curMat.counts, design = design, normalize="quantile")
      fit <- lmFit(v, design)

      allGroups <- unique(allSamples[,groupType])
      print(paste0("Comparing across ",length(allGroups)," ",groupType))
      for (group in sort(allGroups)) {
         print(paste0("Finding ",groupType," ",group,"-specific genes"))

         otherGroups <- setdiff(allGroups,group)

         ct <- relevel(ct, ref=group)
         design <- model.matrix(~ct)
         fit <- lmFit(v, design)
         fit$genes <- rownames(v)
         fit2 <- fit[,c(paste0("ct",otherGroups))]
         fit2 <- eBayes(fit2, trend=TRUE)

         results <- decideTests(fit2, lfc=1, adjust.method="BH", p.value=0.05)

         posSig <- rowSums(results<0)==length(otherGroups)
         print(paste0("-> positive signature: ", paste0(rownames(fit2[posSig,]),collapse=", ")))

         negSig <- rowSums(results>0)==length(otherGroups)
         print(paste0("-> negative signature: ", paste0(rownames(fit2[negSig,]),collapse=", ")))

         posGenes <- rownames(fit2[posSig,])
         negGenes <- rownames(fit2[negSig,])

         diffGenes[[groupType]][[group]] <- list(
            posSig = posGenes,
            negSig = negGenes)

         if (length(posGenes) + length(negGenes) > 0) {
            print("---> drawing heatmap")
            curGenes <- c(posGenes,negGenes)
            print(paste0("---> length curGenes = ",length(curGenes)))
            print(paste0("---> length posSig=",length(posGenes)))
            print(paste0("---> length negSig=",length(negGenes)))
            curMat <- tpmMat[curGenes,]

            curMat <- curMat + 1
            curMat <- log2(curMat / apply(curMat,1,mean))

            curMat[curMat < -2] <- -2
            curMat[curMat > 2] <- 2
            curBreaks <- seq(-2, 2, length.out=101)

            rowGaps <- NULL
            if (length(posGenes) > 0 & length(negGenes) > 0) {
               rowGaps <- length(posGenes) }

            print(paste0("rowGaps at: ", paste0(rowGaps, collapse=", ")))
            print(paste0("in a matrix that has ", nrow(curMat)," rows"))

            curHeight <- 7.5
            curWidth <- 5

            rowFont <- 6

            if (nrow(curMat) > 50) { rowFont <- 3.5 }
            if (nrow(curMat) > 80) { rowFont <- 3 }

            pdf(paste0(dat$specs$outDir,"/voom_onemat_sigGenes_heatmap.",
                       groupType,"_",group,".pdf"), onefile=FALSE,family="ArialMT",
                height=curHeight, width=curWidth)
            pheatmap(curMat,
                     main="relative expression log2(TPM/mean)",
                     border_color = NA,
                     show_rownames = TRUE,
                     show_colnames = TRUE,
                     gaps_col = tpmCol.gaps,
                     gaps_row = rowGaps,
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     annotation_legend = TRUE,
                     annotation_names_col = TRUE,
                     breaks = curBreaks,
                     color = colorRampPalette(c("steelBlue2","white","darkOrange2"))(100),
                     fontsize      = 8,
                     fontsize_row  = rowFont,
                     font.main = 1
            )
            dev.off()
         }

         if (0) {
            pdf(paste0(dat$specs$outDir,"/voom_onemat_venn.",
                       groupType,"_",group,".pdf"))
            vennDiagram(results)
            dev.off()
         }

         write.fit(fit2, file=paste0(dat$specs$outDir,"/voom_onemat_pvals.",groupType,"_",group,".txt"),
                   adj="BH", sep="\t",digits=3)

         if (0) {
         tx <- read.table(paste0(dat$specs$outDir,"/voom_onemat_pvals.",groupType,"_",group,".txt"),
                          header=TRUE, sep="\t", quote="")
         txCols <- colnames(tx)
         txCols <- txCols[!grepl("Genes", txCols)]
         txCols <- c("Genes", txCols)
         tx <- tx[,txCols]
      
         outMat <- tx
      
         write.table(outMat, file=paste0(dat$specs$outDir,"/allexpr_voom_pvals.",groupType,"_",group,".txt"),
                     col.names=TRUE,sep="\t", row.names=FALSE,quote=FALSE)
   
         print("---------------------------------------------------------------")
         }
      }
   }

   outTable <- NULL
   for (groupType in names(diffGenes)) {
      for (group in names(diffGenes[[groupType]])) {
         curRow <- c(groupType,
                     group,
                     length(diffGenes[[groupType]][[group]]$posSig),
                     length(diffGenes[[groupType]][[group]]$negSig))
         
         print(paste(curRow))
         outTable <- c(outTable, curRow)
      }
   }
   outTable <- as.data.frame(matrix(ncol=4,data=outTable,byrow=TRUE))
   colnames(outTable) <- c("groupType", "group", "posSig", "negSig")
   pandoc.table(outTable)
   return(list(outTable = outTable, diffGenes = diffGenes))

}

writeSubtypeMarkers <- function(data, figName1, figName2) {

   outFn <- paste0(figName1,"_marker_gene_summary.txt")
   outTab <- data$subtypeMarkers$outTable
   colnames(outTab) <- c("groupType", "group", "enriched", "depleted")
   write.table(outTab, outFn, col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

   outTab2 <- data$subtypeMarkers$outTable

   outTab <- NA
   newTab <- data.frame(groupType = character(),
                        group =  character(),
                        direction = character(),
                        genes = character(),
                        stringsAsFactors=FALSE)

   for (groupType in names(data$subtypeMarkers$diffGenes)) {
      for (group in names(dat$subtypeMarkers$diffGenes[[groupType]])) {
      for (direction in c("posSig","negSig")) {
         dirLabel <- "enriched"
         if (direction == "negSig") { dirLabel <- "depleted" }
         print("HERE")

         curGenes <- data$subtypeMarkers$diffGenes[[groupType]][[group]][[direction]]
         if (length(curGenes) == 0) {next;}

         curGeneString <- paste0(curGenes,collapse=", ")
         print("HERE1")
         newTab <- rbind(newTab,
            data.frame(groupType = groupType,
              group = group,
              direction = dirLabel,
              genes = curGeneString))
         print("HERE2")

         print(newTab)

         print("HER31")
         curTab <- data.frame(groupType = groupType,
                              group = group,
                              direction = dirLabel,
                              gene = curGenes,
                              stringsAsFactors = FALSE
                              )
         if (!is.data.frame(outTab)) {
            outTab <- curTab
         } else {
            outTab <- rbind(outTab, curTab)
         }
      }
      }
   }

   outFn <- paste0(figName2,"_marker_genes.tsv")
   write.table(newTab, outFn, col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

}


writeSleuthOutput <- function(comparisonName,
                              so = NULL,
                              sleuthResults,
                              designMat,
                              upGenes = NULL,
                              downGenes = NULL,
                              outDir){

# Purpose: Given SLEUTH return object, write DE list to file
#          ALT: save R object dump
#             name with MD5 hash of string of samples compared

   if (!file.exists(outDir)){
      dir.create(outDir, recursive = TRUE) }

   print(paste0("Saving SLEUTH results to ", outDir))

# Write comparison string to disk
   fileConn <- file(paste0(outDir, "/comparison.txt"))
   writeLines(comparisonName, fileConn)
   close(fileConn)

# Dump Sleuth objects to disk
   if (!is.null(so)){ saveRDS(so, file = paste0(outDir, "/so.rds")) }
   saveRDS(sleuthResults, file = paste0(outDir, "/sleuthResults.rds"))
   saveRDS(designMat, file = paste0(outDir, "/designMat.rds"))

# Write gene list to disk
   if (!(is.null(upGenes))){
      upFn <- paste0(outDir, "/upGenes.txt")
      write.table(upGenes, file = upFn, quote = FALSE,
                  row.names = FALSE, col.names = FALSE)
   }

   if (!(is.null(downGenes))){
      downFn <- paste0(outDir, "/downGenes.txt")
      write.table(downGenes, file = downFn, quote = FALSE,
                  row.names = FALSE, col.names = FALSE)
   }

# Put file down saying run done
   fileConn <- file(paste0(outDir, "/done.txt"))
   writeLines("1", fileConn)
   close(fileConn)

   return(1)
}

readSleuthOutput <- function(comparisonName, resultsDir){

   print(paste0("Reading SLEUTH results from ", resultsDir))

   comparisonName <- scan(file = paste0(resultsDir, "/comparison.txt"),
                        what = "character")

   so <- NA
   soFn <- paste0(resultsDir, "/so.rds")
   if (file.exists(soFn)) { so <- readRDS(file = soFn) }

   sleuthResults <- readRDS(file = paste0(resultsDir, "/sleuthResults.rds"))
   designMat <- readRDS(file = paste0(resultsDir, "/designMat.rds"))
   upGenes   <- scan(file = paste0(resultsDir, "/upGenes.txt"), what = "character")
   downGenes <- scan(file = paste0(resultsDir, "/downGenes.txt"), what = "character")

   return(list(
      comparisonName= comparisonName,
      so            = so,
      designMat     = designMat,
      sleuthResults = sleuthResults,
      upGenes       = upGenes,
      downGenes     = downGenes
   ))

}


setCellGroups <- function() {

   cellGroups <- list(c("MB607","MB131"),
                      c("MB370","MB418"),
                      c("MB371","MB185","MB594"))

   labels <- list(c(expression(paste(gamma,"d")),
                    expression(paste(gamma,"m"))),
                  c(expression(paste(alpha,"'/", beta,"'ap")),
                     expression(paste(alpha,"'/", beta,"'m"))),
                  c(expression(paste(alpha,"/", beta,"p")),
                     expression(paste(alpha,"/", beta,"s")),
                     expression(paste(alpha,"/", beta,"c"))))

   labels.simple <- list(c("g.d", "g.m"),
                         c("a'/b'ap", "a'/b'.m"),
                         c("a/b.p", "a/b.s", "a/b.c"))

   return(list(cells = cellGroups,
               labels = labels,
               labels.simple=labels.simple))

}


barplotTPM <- function(dat, geneList, outPre = "",
                       figName = "barplotTPM",
                       sumGenes=FALSE, sumLabel,
                       justPrintSamples = FALSE,
                       makePlot = TRUE
                       ) {
# Make TPM barplot for list of genes
# add dots for individual rep points.

   cellGroups <- setCellGroups()

   condSamples <- list()
   for (driver in unlist(cellGroups$cells)) {
      condSamples[[driver]] <- dat$specs$sampleInfo$sample_name[
                                  dat$specs$sampleInfo$driver == driver]
   }

   if (justPrintSamples) {
      curSampleList <- as.vector(unlist(condSamples))
      figDescrFn <- paste0(dat$specs$outDir, "/", figName, "_samples.txt")
      write.table(curSampleList, figDescrFn,
                  row.names = FALSE, col.names = FALSE, quote = FALSE)
      return(1)
   }


   tpmCols <- c()
   for (driver in names(condSamples)) {
      tpmCols <- c(tpmCols, condSamples[[driver]]) }
   tpmCols <- paste0("tpm.",tpmCols)

   eMat <- dat$expr$geneExpr[match(geneList[geneList %in%
                                            dat$expr$geneExpr$gene_name],
                                    dat$expr$geneExpr$gene_name),
                             c("gene_name", tpmCols),
                             drop=FALSE]
   if (sumGenes == TRUE) {
      new.eMat <- as.data.frame(t(apply(eMat[,tpmCols], 2, sum)))
      colnames(new.eMat) <- tpmCols
      new.eMat$gene_name <- sumLabel
      eMat <- new.eMat
      geneList <- sumLabel
   }

   newCols <- paste0("avgTPM.",names(condSamples))
   newCol2srcCol <- list()
   for (cond in names(condSamples))  {
      srcCol <- paste0("tpm.",condSamples[[cond]])

      newCol <- paste0("avgTPM.",cond)
      newCol2srcCol[[newCol]] <- srcCol

      eMat[,newCol] <- apply(eMat[,srcCol],1,mean)
   }
   print(paste("newCols is: ", newCols, collapse=", "))
   eMat<- as.data.frame(eMat, stringsAsFactors=FALSE,as.is=TRUE)
   rownames(eMat) <- eMat$gene_name

   cellLabels = c(cellGroups$labels[[1]][1],
                  cellGroups$labels[[1]][2],
                  cellGroups$labels[[2]][1],
                  cellGroups$labels[[2]][2],
                  cellGroups$labels[[3]][1],
                  cellGroups$labels[[3]][2],
                  cellGroups$labels[[3]][3])

   if (makePlot) {
   for (gene in geneList) {
      print(paste0("NOW IN on gene ", gene))

      outFn <- paste0(dat$specs$outDir,"/",figName,".barTPM_", gene,".pdf")
      outFn <- gsub(" ","_",outFn)
      outFn <- gsub("\\(","_",outFn)
      outFn <- gsub("\\)","_",outFn)

      print(paste0("-> printing to ", outFn))

      pdf(outFn, family = "ArialMT", height = 3.5, width  = 5)
      par(mar=c(4,5,2,1),xpd=TRUE)
      curRow <- t(eMat[gene,newCols])
      yRange <- range(curRow)
      yMax <- yRange[2]
      barplot(curRow,
              beside=TRUE,
              names=NA,
              space=c(0,0, 0.5, 0, 0.5, 0, 0),
              las=2,
              col=c("darkGrey"),
              border="black",
              ylim=c(0,yRange[2]),
              cex.lab=1.2,
              cex.axis=1.2,
              font.main=1,
              main=gene)
      mtext(side=2,line=3.5,"Expression (TPM)",cex=1.2)
      realX <- c(c(0:1,2.5:3.5,5:7) + 0.5)
      axis(side=1, las=1, lwd=0,
           at=realX,
           labels= cellLabels,
           cex=1.2,
           line=0)

      for (i in 1:length(newCols)) {
         points(c(rep(realX[i], length(newCol2srcCol[[newCols[i]]]))),
                eMat[gene,newCol2srcCol[[newCols[i]]]],
                pch=20)
      }

      dev.off()
   }
   }

   avgEmat <- eMat[,newCols]
   colnames(avgEmat) <- cellLabels
   return(avgEmat)

}
