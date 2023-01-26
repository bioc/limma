sampleInfoFromGEO <- function(file, remove.constant.columns=TRUE)
# Parse sample information from a GEO series matrix file into character matrices
# Gordon Smyth
# Created 7 Nov 2021.
{
  con <- file(file, "r")  
  on.exit(close(con))
  i <- 0
  repeat {
    i <- i+1
    txt <- readLines(con,n=1)
    if(!length(txt)) stop("Sample information not found in file. Input should be a GEO series matrix file.")
    if(identical(substring(txt,1,8),"!Sample_")) break
  }
  Lines <- list()
  i <- 0
  repeat {
    i <- i+1
    Lines[[i]] <- txt
    txt <- readLines(con,n=1)
    if(!length(txt)) break
    if(!identical(substring(txt,1,8),"!Sample_")) break
  }
  Lines <- lapply(Lines,function(x) strsplit(x,split="\t")[[1]])
  NSamples <- lengths(Lines)
  Short <- which(NSamples < max(NSamples))
  if(length(Short)) for (i in Short) Lines[[i]] <- NULL
  NSamples <- NSamples[1] - 1L
  Names <- vapply(Lines, FUN=function(x) substring(x[1],9,nchar(x[1])), FUN.VALUE="")
  Lines <- vapply(Lines, FUN=function(x) gsub("\"","",x[-1]), FUN.VALUE=character(NSamples))
  colnames(Lines) <- Names
  rownames(Lines) <- 1:NSamples

# Parse characteristics columns
  i <- grep("characteristics_ch1",colnames(Lines))
  if(length(i)) {
    CharacteristicsCh1 <- .parseGEOSampleCharacteristics(Lines[,i,drop=FALSE])
    Lines <- Lines[,-i]
  } else {
    CharacteristicsCh1 <- NULL
  }
  i <- grep("characteristics_ch2",colnames(Lines))
  if(length(i)) {
    CharacteristicsCh2 <- .parseGEOSampleCharacteristics(Lines[,i,drop=FALSE])
    Lines <- Lines[,-i]
  } else {
    CharacteristicsCh2 <- NULL
  }

  if(remove.constant.columns && length(Lines) > 0L) {
    jl <- .columnsAllEqual(Lines)
    if(any(jl)) Lines <- Lines[,!jl,drop=FALSE]
  }
  if(remove.constant.columns && length(CharacteristicsCh1) > 0L) {
    jl <- .columnsAllEqual(CharacteristicsCh1)
    if(any(jl)) CharacteristicsCh1 <- CharacteristicsCh1[,!jl,drop=FALSE]
  }
  if(remove.constant.columns && length(CharacteristicsCh2) > 0L) {
    jl <- .columnsAllEqual(CharacteristicsCh2)
    if(any(jl)) CharacteristicsCh2 <- CharacteristicsCh2[,!jl,drop=FALSE]
  }

  list(SampleInfo=Lines,CharacteristicsCh1=CharacteristicsCh1,CharacteristicsCh2=CharacteristicsCh2)
}

.parseGEOSampleCharacteristics <- function(x)
# Parse the sample characteristics rows of a GEO series matrix file
# Gordon Smyth
# Created 7 Nov 2021.
{
  nsamples <- nrow(x)
  characteristic.value <- strsplit2(x,split=": ")
  characteristic <- characteristic.value[,1]
  value <- characteristic.value[,2]
  Characteristics <- unique(characteristic)
  Characteristics <- Characteristics[Characteristics != ""]
  ncharacteristics <- length(Characteristics)
  targets <- matrix(NA_character_,nsamples,ncharacteristics)
  rownames(targets) <- 1:nsamples
  colnames(targets) <- Characteristics
  sample <- rep(1:nsamples,ncol(x))
  for (i in 1:length(x)) {
    if(value[i] != "") targets[sample[i],characteristic[i]] <- value[i]
  }
  targets
}

.columnsAllEqual <- function(x)
# Identify which columns of a matrix have all elements the same
# Gordon Smyth
# Created 7 Nov 2021.
{
  x <- as.matrix(x)
  n <- nrow(x)
  if(n <= 1L) return(rep_len(TRUE,ncol(x)))
  Eq <- x[1:(n-1),,drop=FALSE] == x[2:n,,drop=FALSE]
  Eq[is.na(Eq)] <- FALSE
  colSums(Eq) == (n-1L)
}
