#' @description
#' Vertex Wise analysis of MGH formatted Data
#' The main function to use is ...
#'
#' @author Ryan Muetzel & Sander Lamballais
#'
#' @examples
#' Put examples in here
#'
#' @references
#' Muetzel, et al, Bullying involvement and cortical morphology in children.
'_PACKAGE'

################################################################################################################
#### Wrote a function (+ methods) that standardizes imputed (and non-imputed) datasets #########################
#### This does require the mice and mi packages if imputed datasets are derived from them ######################
#### (other imputation packages do not require to be loaded)  ##################################################
################################################################################################################

imp2list <- function(x) UseMethod("imp2list", x)


imp2list.amelia <- function(x) x$imputations
imp2list.aregImpute <- function(x) stop("The aregImpute format is not supported,
                                        as it does not include the raw dataset.
                                        Please create your own design matrices.")
imp2list.data.frame <- function(x) list(x)
# imp2list.list <- function(x) {
#   s <- t(sapply(x, names))
#   if (all(is.null(s))) stop("The provided datasets (in the list) do not have column names")
#   if (any(is.null(s))) stop("Some provided datasets (in the list) do not have column names")
#   count <- sum(!duplicated(s))
#   if (count != 1) s else stop("The supplied datasets (in the list) do not share column names.")
# }
imp2list.matrix <- function(x) imp2list(as.data.frame(x))
imp2list.mi <- function(x) mi::complete(x)
imp2list.mids <- function(x) lapply(seq_len(x$m), function(y) mice::complete(x, y))
imp2list.missForest <- function(x) list(x$ximp)



################################################################################################################
#### Using this so often, just make it a separate function######################################################
################################################################################################################
message2 <- function(verbose, ...) {
  if (verbose) message(...)
}

pool_quick <- function(cm, sm){
  m <- nrow(cm)
  k <- ncol(cm)
  cp <- rep(1/m, m) %*% cm
  smn <- rep(1/m, m) %*% (sm^2)
  e <- cm - matrix(cp, nrow = m, ncol = k, byrow = TRUE)
  b <- (t(e) %*% e)/(m - 1)
  t <- smn + (1+1/m) * diag(b)
  sp <- sqrt(t)
  return(list(cp, sp))
}


#' Load an MGH file into memory
#'
#' @param input.file full path to the mgh file
#'
#' @return mgh object with data and various header elements
#'
#' @examples
#' mgh <- load.mgh('/export/local/scratch/lh.fwhm10.fsaverage.thickness.mgh')
#'
#'
load.mgh <- function(input.file) {
  # written by Heath Pardoe, heath.pardoe at nyumc.org, 09/12/2013
  to.read <- file(input.file, "rb")
  v <- readBin(to.read, integer(), endian = "big")
  ndim1 <- readBin(to.read, integer(), endian = "big")
  ndim2 <- readBin(to.read, integer(), endian = "big")
  ndim3 <- readBin(to.read, integer(), endian = "big")
  nframes <- readBin(to.read, integer(), endian = "big")
  type <- readBin(to.read, integer(), endian = "big")
  dof <- readBin(to.read, integer(), endian = "big")
  close(to.read)

  to.read <- file(input.file,"rb")
  dump <- readBin(to.read, double(), size = 4, n = 71, endian = "big")
  #THIS SEEMED TO ONLY READ IN THE FIRST CHUNK OF DATA, and ignored if additional subjects were merged in
  #added ndim1*nframes
  x <- readBin(to.read ,double(), size = 4, n = ndim1*nframes, endian = "big")
  close(to.read)
  list(x = x, v = v, ndim1 = ndim1, ndim2 = ndim2, ndim3 = ndim3, nframes =
         nframes, type = type, dof = dof)
}

#' Save out an MGH file from memory
#'
#' @param vol MGH object (as from load.mgh)
#' @param fname file name to be used to save out the data
#'
#' @examples
#' save.mgh(mgh, '/export/local/scratch/pval.mgh')
#'
save.mgh <-function(vol,fname) {

  # R translation of save_mgh.m
  # written by Heath Pardoe, heath.pardoe at nyumc.org, 09/12/2013

  MRI.UCHAR <-  0
  MRI.INT <-    1
  MRI.LONG <-   2
  MRI.FLOAT <-  3
  MRI.SHORT <-  4
  MRI.BITMAP <- 5
  MRI.TENSOR <- 6
  slices <- c(1:256)

  fid <- file(fname, open = "wb", blocking = TRUE)

  width <- vol$ndim1
  height <- vol$ndim2
  depth <- vol$ndim3
  #RLM added
  nframes <- vol$nframes

  writeBin(as.integer(1), fid, size = 4, endian = "big")
  writeBin(as.integer(width), fid, size = 4, endian = "big")
  writeBin(as.integer(height), fid, size = 4, endian = "big")
  writeBin(as.integer(depth), fid, size = 4, endian = "big")
  #here we replace the default of 1 frame
  writeBin(as.integer(nframes), fid, size = 4, endian = "big")
  #writeBin(as.integer(1),fid,size = 4, endian = "big")

  # HP note: I've ignored all the diffusion tensor stuff
  writeBin(as.integer(MRI.FLOAT), fid, size = 4, endian = "big")

  writeBin(as.integer(1), fid, size = 4, endian = "big")
  # dof = fread(fid, 1, 'int');
  ## HP note: ignored ^this line^ from save_mgh.m

  UNUSED.SPACE.SIZE <- 256
  USED.SPACE.SIZE <- (3 * 4 + 4 * 3 * 4)  # space for ras transform

  unused.space.size <- UNUSED.SPACE.SIZE - 2

  # ignored all the stuff about "M" - could probably do it if necessary so let me know
  #if (nargin > 2)
  ## fwrite(fid, 1, 'short')        # ras.good.flag <- 0
  # writeBin(1,fid,size = 2, endian = "big")
  # unused.space.size <- unused.space.size - USED.SPACE.SIZE
  ## fwrite(fid, sizes(1), 'float32')  # xsize
  ## fwrite(fid, sizes(2), 'float32')  # ysize
  ## fwrite(fid, sizes(3), 'float32')  # zsize
  #
  # fwrite(fid, M(1,1), 'float32')   # x.r
  # fwrite(fid, M(2,1), 'float32')   # x.a
  # fwrite(fid, M(3,1), 'float32')   # x.s

  # fwrite(fid, M(1,2), 'float32')   # y.r
  # fwrite(fid, M(2,2), 'float32')   # y.a
  # fwrite(fid, M(3,2), 'float32')   # y.s

  # fwrite(fid, M(1,3), 'float32')   # z.r
  # fwrite(fid, M(2,3), 'float32')   # z.a
  # fwrite(fid, M(3,3), 'float32')   # z.s

  # fwrite(fid, M(1,4), 'float32')   # c.r
  # fwrite(fid, M(2,4), 'float32')   # c.a
  # fwrite(fid, M(3,4), 'float32')   # c.s
  #else
  # fwrite(fid, 0, 'short')        # ras.good.flag <- 0
  writeBin(as.integer(0), fid, size = 2, endian = "big")

  #   } #

  writeBin(as.integer(rep.int(0, unused.space.size)), fid, size = 1)
  bpv <- 4    # bytes/voxel
  nelts <- width * height   # bytes per slice
  #writeBin(vol$x, fid, size = 4, endian = "big")
  writeBin(vol$x, fid, size = 4, endian = "big")
  close(fid)
}

#'call the FreeSurfer tool mris_preproc to merge surface data into a single mgh file
#'
#' @param DIR the output folder name, should have a trailing /
#' @param idVar The id variable in your data frame that matches the id variable in the SUBJECTS_DIR
#' @param dataFrame The data frame with the list of subjects you want to merge
#' @param target usually fsaverage, unless you made a custom template
#' @param hemi lh or rh
#' @param cacheIn This is the root output name of the qcache recon-all. example: thickness.fwhm10.fsaverage
#' @param cacheOut This is the name of the file you want to save out in mgh format
#' @param verbose defaults to FALSE to supress loads out output from the merging process.
#'
#' @examples
#' mergeMgh('/export/local/scratch/tests', 'subjid', myData, 'lh', 'thickness.fwhm10.fsaverage', 'lh.thickness.fwhm10.fsaverage.mgh', verbose = verbose)
#'
mergeMgh <- function(DIR, idVar, dataFrame, hemi, cacheIn, cacheOut, target = "fsaverage", verbose = FALSE) {
  ################################################################################################################
  #### Why is DIR not being used? ################################################################################
  ################################################################################################################
  #cacheOut <- paste(DIR, cacheOut, '.mgh', sep='')
  idList <- gsub(".mgh", "_idList.txt", cacheOut)
  message2(verbose, "Writing ID numbers to file: \n", idList, "\n", cacheOut)
  write(dataFrame[[idVar]], file = idList, ncolumns = 1)
  cmdStr <- paste("mris_preproc", "--f",  idList, "--target", target, "--hemi",
                  hemi, "--cache-in", cacheIn, "--out", cacheOut)
  message2(verbose, cmdStr)
  system(cmdStr, ignore.stdout = !verbose)
}

#' Compute smoothness of residuals for mc-z correction
#'
#' @param DIR the output folder name, should have a trailing /
#' @param hemi lh or rh
#' @param eres the residual output from a vertex wise analysis
#' @param mghBaseName to help name the output dat file with the fwhm value
#'
#' @return fwhm numeric
#'
#' @examples
#' fwhm <- calcFwhm('/export/local/scratch/test/', 'lh', '/export/local/scratch/test/test.eres.mgh', 'test.vertexWise')

calcFwhm <- function(DIR, hemi, eres, mghBaseName, mask = NULL, target = "fsaverage", verbose = FALSE) {
  message2(verbose, "Calculating smoothing for multiple testing correction.")
  mghBaseName <- paste0(hemi, '.', mghBaseName)
  oDat <- paste0(DIR, mghBaseName, ".fwhm.dat")
  outMask <- paste0(DIR, mghBaseName, ".finalMask.mgh")
  cmdStr <- paste("mris_fwhm", "--i", eres, "--hemi", hemi, "--subject", target, "--prune", "--cortex", "--dat", oDat, "--out-mask", outMask)
  if (!is.null(mask)) paste(cmdStr, "--mask", mask)
  message2(verbose, cmdStr)
  system(cmdStr, ignore.stdout = !verbose)
  fwhm <- read.table(oDat)
  fwhm <- round(fwhm)
  return(fwhm)
}

#' Apply mc-z simulation correction to a log10p map
#'
#' @param DIR the output folder name, should have a trailing /
#' @param hemi lh or rh
#' @param mghBaseName to help name the output dat file with the fwhm value
#' @param fwhm (from calcFwhm)
#' @param mask mgh mask file, default is to just use the cortex label option
#' @param cwpThr cluster wise p thresh...default 0.025...i.e., 'two hemi'
#' @param mczThr mcz threshold...default is 3
#' @param csdSign sign for significance testing...abs (default), pos or neg
#' @param verbose allow terminal output
#'

runMriSurfCluster <- function(DIR, hemi, mghBaseName, fwhm, mask = NULL, cwpThr = 0.025, mczThr = NULL, csdSign = "abs", verbose = FALSE) {
  ################################################################################################################
  #### If default mczThr = 2, then where does that 13 come from? #################################################
  ################################################################################################################
  if (is.null(mczThr)){
    mczThr = paste0("th", as.character(30))
  }
  else {
    mczThr =  paste0("th", as.character(mczThr))
  }

  mghBaseName <- paste0(hemi, '.', mghBaseName)
  oBaseName <- paste0(DIR, mghBaseName, ".cache.", mczThr, '.', csdSign, ".sig")
  log10p <- paste0(DIR, mghBaseName, ".log10p.mgh")
  if (fwhm < 10) fwhm <- paste0('0', fwhm)
  csd <- file.path(Sys.getenv("FREESURFER_HOME"), "average/mult-comp-cor/fsaverage", hemi, paste0("cortex/fwhm", fwhm), csdSign, mczThr, "mc-z.csd")
  cmdStr <- paste("mri_surfcluster",
                  "--in", log10p,
                  "--csd", csd,
                  "--cwsig", paste0(oBaseName, ".cluster.mgh"),
                  "--vwsig", paste0(oBaseName, ".voxel.mgh"),
                  "--sum", paste0(oBaseName, ".cluster.summary"),
                  "--ocn",  paste0(oBaseName, ".ocn.mgh"),
                  "--oannot", paste0(oBaseName, ".ocn.annot"),
                  "--annot", "aparc",
                  "--csdpdf",  paste0(oBaseName, ".pdf.dat"),
                  "--cwpvalthresh", cwpThr,
                  "--o", paste0(oBaseName, ".masked.mgh"),
                  "--no-fixmni",
                  "--surf", "white")
  cmdStr <- if (!is.null(mask)) paste(cmdStr, "--mask", mask) else paste(cmdStr, "--cortex")
  system(cmdStr)
}

#' Run vertex wise analysis using mgh data
#' The "_Y" denotes that the vertex data are the outcome measure in the LM
#'
#' @param DIR the output folder name, should have a trailing /
#' @param idVar The id variable in your data frame that matches the id variable in the SUBJECTS_DIR
#' @param dataFrame The data frame with the list of subjects you want to merge
#' @param target usually fsaverage, unless you made a custom template
#' @param hemi lh or rh
#' @param x variable of interest in data frame (e.g., age)
#' @param covars vector of variables you want added to the model
#' @param mask mgh file to mask analysis. default is to use the cortex label from the target
#' @param con.vec contrast vector if you want something other than just the positve contrast of x
#' @param measure thickness, pial_lgi, etc
#' @param fwhm smoothing kernel used on mgh data (defualt is 10)
#' @param nCores the number of cores you want to use (for Parallelization). Not sure how well this works with more than 10 cores)
#' @param verbose defaults to FALSE to supress loads out output from the merging process.
#'
#'

run_lm_VertexWiseParRcpp_Y <- function(idVar,
                                       dataFrame,
                                       hemi,
                                       outDir,
                                       x,
                                       covars = NULL,
                                       fitEq = NULL,
                                       mask = NULL,
                                       con.vec = NULL,
                                       target = "fsaverage",
                                       measure = "thickness",
                                       fwhm = ifelse(measure == "pial_lgi", 5, 10),
                                       mgh = NULL,
                                       mghBaseName,
                                       cleanUp = TRUE,
                                       cleanUpBM = TRUE,
                                       tmpDir = outDir,
                                       nCores = 2,
                                       clobber = FALSE,
                                       writeDesign = TRUE, ################ NEW ARGUMENT!!!
                                       verbose = FALSE) {
  #this will likely fail if someone uses another target...
  #target <- match.arg(target)
  measure <- match.arg(measure, c('area', 'area.pial', 'curv', 'jacobian_white', 'pial_lgi', 'sulc', 'thickness', 'volume', 'w-g.pct.mgh', 'white.H', 'white.K'))
  if (substr(outDir, nchar(outDir), nchar(outDir)) != "/") outDir <- paste0(outDir, "/")

  # Set up some prefixes
  mghHemiBaseName <- paste0(hemi, '.', mghBaseName)
  outDirHemiBase <- paste0(outDir, mghHemiBaseName, '/')

  # Set up some data-related things
  dataFrame2 <- imp2list(dataFrame)
  dataFrameTest <- dataFrame2[[1]]

  ################################################################################################################
  #### Requiring functions belongs in the DESCRIPTION file of the package ########################################
  ################################################################################################################
  #require("foreach")
  #require("doMC")
  #require("bigmemory")
  #require("RcppEigen")
  pkgs <- c('foreach', 'doMC', 'bigmemory', 'RcppEigen')

  if (is.null(fitEq)){
    stop('You must supply a regression equation...')
  }
  # Check if outDirHemiBase exists (and remake if allowed), otherwise create (outDir &) outDirHemiBase
  if (dir.exists(outDirHemiBase)){
    if (clobber){
      unlink(outDirHemiBase, recursive = TRUE)
      dir.create(outDirHemiBase)
    } else {
      stop("Output folder already exists, so the function will exit to avoid overwriting:\n",
           outDirHemiBase,
           "\nUse the clobber=TRUE argument to overwrite an existing folder.")
    }
  } else {
    if (!dir.exists(outDir)) dir.create(outDir)
    dir.create(outDirHemiBase)
  }

  ################################################################################################################
  #### Make new function that generates masks if they're not there  ##############################################
  ################################################################################################################
  if (is.null(mask)){
    maskDir <- file.path(Sys.getenv("SUBJECTS_DIR"), "qdec", 'R', "masks")
    mask <- paste(hemi, target, "cortex.mask.mgh", sep = '.')
    mask <- file.path(maskDir, mask)
  }
  if (!file.exists(mask)) stop("Mask file does not exist. \n", mask)

  # Load the mask
  maskData <- load.mgh(mask)

  ################################################################################################################
  #### Still need to check this part #############################################################################
  ################################################################################################################
  #if an mgh file is supplied, double check the order
  if(!is.null(mgh)) {
    if (file.exists(mgh)){
      message2(verbose, "Checking pre-defined MGH file for subject ID order.")
      idList <- gsub(".mgh", "_idList.txt", mgh)
      idVals <- read.table(idList)
      if (!identical(dataFrameTest[,idVar], idVals[,1])){
        stop("Something is wrong. Ids in data frame do not match MGH file!")
      }
    }
  }
  #otherwise, create one
  else {
    if (!is.null(measure) & !is.null(target) & !is.null(fwhm)){
      fwhmCache <- paste0("fwhm", fwhm)
      cacheIn <- paste(measure, fwhmCache, target, sep = '.')
      message2(verbose, "Creating cache file for: ", cacheIn)
    }
    else {
      stop("To create MGH file, you need to specify the measure, fwhm and target.")
    }
    message("Merging MGH file.")
    mgh <- paste0(outDirHemiBase, mghHemiBaseName, ".mgh")
    mergeMgh(DIR = outDirHemiBase, idVar = idVar, dataFrame = dataFrameTest,
             hemi = hemi, cacheIn = cacheIn, cacheOut = mgh, target = target,
             verbose = verbose)
  }

  on.exit({
    if (cleanUp){
      message2(verbose, "Cleaning up MGH files.")
      file.remove(mgh)
    }
  }, add = TRUE)

  # read the data into memory
  message("Reading mgh data into memory.")
  mghData <- load.mgh(mgh)

  #In order for the parallel processing not to puke on memory, we need to simplify what we give it
  #ie we can't give it the full mghMem object...it will fork that to each sub process.
  nVertices <- mghData$ndim1
  nFrames <- mghData$nframes
  #if an mgh file is not supplied, make one
  #i'm not sure, but the .inorder arg may be slowing things down...
  #now with big mem we don't care
  #in any case, the results seem to come back in order, strangely enough

  oldwd <- setwd(outDirHemiBase)
  on.exit(setwd(oldwd), add = TRUE)

  #convert the long vector to a matrix and then a bigmem matrix
  backFile <- paste0(mghHemiBaseName, ".bin")
  backStatsFile <- gsub(".bin", "_stats.bin", backFile)
  backEresFile <- gsub(".bin", "_eres.bin", backFile)
  descFile <- paste0(mghHemiBaseName, ".desc")
  descStatsFile <- gsub(".desc", "_stats.desc", descFile)
  descEresFile <- gsub(".desc", "_eres.desc", descFile)
  message2(verbose, "Creating big.memory matrices")
  #on some clusters, it might be helpful to point to a super fast disk (e.g., ssd drive or node-specific disk)
  #otherwise the IO slows down the parallezation (use tmpDir)
  mghData$x <- bigmemory::as.big.matrix(x = matrix(mghData$x, nrow = mghData$ndim1, ncol = mghData$nframes),
                             backingfile = backFile,
                             descriptorfile = descFile,
                             backingpath = tmpDir)
  #make an empty matrix to hold the vertex-wise pvals, t-vals and b-coefs
  vertexStatData = bigmemory::as.big.matrix(x = matrix(1, nrow = mghData$ndim1, ncol = 4),
                                 backingfile = backStatsFile,
                                 descriptorfile = descStatsFile,
                                 backingpath = tmpDir)
  vertexEresData = bigmemory::as.big.matrix(x = matrix(0, nrow = mghData$ndim1, ncol = mghData$nframes),
                                 backingfile = backEresFile,
                                 descriptorfile = descEresFile,
                                 backingpath = tmpDir)

  on.exit({
    if (cleanUpBM){
      message2(verbose, "Cleaning up big.memory backing files.")
      cleanFiles <- c(backFile, backStatsFile, backEresFile,
                      descFile, descStatsFile, descEresFile)
      file.remove(paste0(tmpDir, cleanFiles))
    }
  }, add = TRUE)

  #a bit of an attempt at parallelizing...
  message2(verbose, "Setting up parallel computing.")
  availCores <- parallel::detectCores()
  if (nCores > availCores) nCores <- availCores - 2
  if (nCores < 1) nCores <- 1
  doMC::registerDoMC(cores = nCores)

  #define the model
  message("Defining and checking the model.")
  #fit_eq <- paste(paste('~', x, " +"), paste(covars, collapse = " + "))
  X_list <- list()
  for (i in seq_along(dataFrame2)) X_list[[i]] <- model.matrix(as.formula(fitEq), data = dataFrame2[[i]])

  X <- X_list[[1]]
  N <- nrow(X)
  k <- ncol(X) - 1
  df <- N - k - 1

  ##################
  #check the data frame for completeness
  nMiss <- sum(is.na(X))
  if (nMiss > 0) stop("Found ", nMiss, " missing data point(s) in variables used in model. Please check the data frame...")
  ##################

  varNames <- attributes(X)$dimnames[[2]]
  message("Variables in X: ")
  cat(varNames, '\n')
  if(is.null(con.vec)){
    message("No contrast vector provided. Determining contrasts based on: ", x)
    #find which column is the variable of interest
    #conInd <- match(x, names(X[1,]))
    #make an empty vector
    #con.vec <- rep(0, ncol(X))
    #re-code the column of interest
    #con.vec[conInd] <- 1
    con.vec <- as.numeric(startsWith(varNames, x))
    names(con.vec) <- varNames
    #make sure only one contrast is defined
    #later we should just loop over each one...
    if (sum(as.numeric(con.vec==1)) != 1){
      cat(con.vec, '\n')
      cat(x, '\n')
      cat(varNames, '\n')
      stop('Contrast vector not properly defined')
    }
  }
  message2(verbose, "Contrast vector is: ", con.vec)
  if (writeDesign){
    message2(verbose, "Writing out design matrix and contrast vector to disk.")
    if (length(X_list) == 1) {
      write.table(X, file = paste0(outDirHemiBase, "X.dat"))
    } else {
      for (i in seq_along(X_list)) write.table(X_list[[i]], file = paste0(outDirHemiBase, "X", i, ".dat"))
    }
    write.table(as.matrix(t(con.vec)), file = paste0(outDirHemiBase, "ctx.dat"), row.names = FALSE)
  }

  message("Fitting model at each vertex.")
  pb <- txtProgressBar(min = 1, max = mghData$ndim1, style = 3)
  #tried to get this to work, but it failed
  #we must need to just require it at the beginning...
  #otherwise it cant find dopar in the namespace
  if (nCores > 1){
    `%dopar%` <- foreach::`%dopar%`
  }
  else {
    `%dopar%` <- foreach::`%do%`
  }
  #consider getting rid of packages....this was only added for the bigr cluster...
  fitMat <- foreach::foreach (j = 1:nVertices, .inorder = FALSE, .packages = pkgs) %dopar% {
  #fitMat <- foreach::foreach (j = 1:nVertices, .inorder = FALSE) %dopar% {
      if (getTxtProgressBar(pb) < j) setTxtProgressBar(pb, j)
      #point to the input surface data
      mghBm <- bigmemory::attach.big.matrix(descFile, backingpath = tmpDir)
      #point to the output matrix
      statBm <- bigmemory::attach.big.matrix(descStatsFile, backingpath = tmpDir)
      eresBm <- bigmemory::attach.big.matrix(descEresFile, backingpath = tmpDir)
      #keep track of the vertex #
      if (j == 1) write(fitEq, file = paste0(outDirHemiBase, "model.txt"))
      #pull out a given vertex, put it into your dataframe
      #fit the model
      #First test if it is within the mask, or if there are any 0s in the vertex
      if (maskData$x[j] == 0 | 0 %in% mghBm[j,]){
        #skip the vertex if not part of mask or has zeros
        bVal <- 0
        beta <- 0
        seVal <- 0
        tVal <- 0
        pVal <- 1
      }
      else{

        y <- mghBm[j,]

        #fm <- RcppBlaze::fastLmPure(X, mghBm[i,], 1)

        m <- length(X_list)
        k <- length(con.vec)
        se_m <- coef_m <- matrix(NA, m, k)
        beta_m <- matrix(NA, m, sum(con.vec))

        for (i in seq_len(m)){
          fm <- RcppEigen::fastLmPure(X_list[[i]], y, method = 2)
          coef_m[i,] <- fm$coefficients
          se_m[i,] <- fm$se

          xsd <- sd(X_list[[i]][,!!con.vec], na.rm = TRUE)
          ysd <- sd(sapply(y, as.numeric), na.rm = TRUE)

          beta_m[i,] <- as.numeric(fm$coefficients %*% con.vec) * xsd / ysd

        }
        if (m == 1){
          coef_pool <- coef_m
          se_pool <- se_m
          beta_pool <- beta_m
        } else {
          pooled <- pool_quick(coef_m, se_m)
          coef_pool <- pooled[[1]]
          se_pool <- pooled[[2]]
          beta_pool <- mean(beta_m)
        }

        bVal <- as.numeric(coef_pool %*% con.vec)
        beta <- beta_pool
        seVal <- as.numeric(se_pool %*% con.vec)
        tVal <- bVal / seVal
        pVal <- 2 * pt(abs(tVal), df = df, lower.tail = FALSE)

        ############## REWRITE TO WORK WITH MI

        if (is.nan(pVal)) {
          pVal <- 1
          tVal <- 0
          bVal <- 0
          beta <- 0
        }
      }
      #Put the stats into the big mem matrix
      statBm[j,1] <- -1 * (log10(pVal))
      statBm[j,2] <- tVal
      statBm[j,3] <- bVal
      statBm[j,4] <- beta
      if (maskData$x[j] == 0 | 0 %in% mghBm[j,]){
        eresBm[j,] <- as.vector(matrix(0, nFrames, 1))
      }
      else {
        ################################################################################################################
        #### should be pooled somehow ##################################################################################
        ################################################################################################################
        eresBm[j,] <- as.vector(fm$residuals)
      }

    }
    close(pb)


  #now save out the data
  message("Saving estimates to mgh format.")
  ##TODO
  ###we have a problem with naming things after the mgh file
  ###if you input an MGH from a different analysis with a different basename
  ###it will revert to saving these new analyses in the old folder with the old name...
  ###so we'll need to fix this in the event peoepl really want to use existing mgh '4d' files...

  eres <- gsub(".mgh", ".eres.mgh", mgh)
  nFrames <- mghData$nframes
  if (object.size(vertexEresData[,]) > 2^31){
    #save out each subject one at a time
    mghData$nframes <- 1
    eresDir <- paste(outDirHemiBase, 'eres',  sep='/')
    dir.create(eresDir)
    eresFileList <- c()
    for (s in 1:nFrames){
      eresMgh <- paste(s, 'eres.mgh', sep='_')
      eresMgh <- paste(eresDir, eresMgh, sep='/')
      mghData$x <- as.vector(vertexEresData[,s])
      save.mgh(mghData, eresMgh)
      eresFileList <- c(eresFileList, eresMgh)
    }
    write(eresFileList, 'eresFileList.txt')
    cmdStr <- paste('mri_concat', '--o', eres, '--f', 'eresFileList.txt')
    system(cmdStr)
    on.exit({
      if (cleanUp){
        message2(verbose, "Cleaning up MGH files (part 2).")
        unlink(eresDir, recursive = TRUE)
      }
    }, add = TRUE)
  }
  else{
    mghData$x <- as.vector(vertexEresData[,])
    save.mgh(mghData, eres)
  }
  on.exit({
      if (cleanUp){
      message2(verbose, "Cleaning up MGH files (part 2b).")
      file.remove(eres)
      }
  }, add = TRUE)

  # Write log10p, t, b and beta values from vertexStatData to separate files
  outStats <- c("log10p", "t", "b", "beta")
  mghData$nframes <- 1
  for (i in seq_along(outStats)){
    oFile <- gsub(".mgh", paste0(".", outStats[i], ".mgh"), mgh)
    mghData$x <- vertexStatData[,i]
    save.mgh(mghData, oFile)
  }

  # Run final mask
  fwhm <- calcFwhm(outDirHemiBase, hemi, eres, mghBaseName, mask = mask, verbose = verbose)
  finalMask <- gsub(".mgh", ".finalMask.mgh", mgh)
  runMriSurfCluster(outDirHemiBase, hemi, mghBaseName, fwhm, mask = finalMask, verbose = verbose)

  # End
  message("Done.")
  return(NULL)
}
