read.idat <- function(idatfiles, bgxfile, dateinfo=FALSE, annotation="Symbol", tolerance=0L, verbose=TRUE)
#	Read GenomeStudio IDAT files for Illumina gene expression BeadChips
#	Matt Ritchie and Gordon Smyth
#	Created 30 September 2013.	Last modified 14 June 2016.
{
#	Need illuminaio package
	OK <- requireNamespace("illuminaio",quietly=TRUE)
	if(!OK) stop("illuminaio package required but is not installed (or can't be loaded)")

#	Check idatfiles
	idatafiles <- as.character(idatfiles)
	nsamples <- length(idatfiles)

#	Initialize EListRaw object
	elist <- new("EListRaw")
	elist$source <- "illumina"
	elist$targets <- data.frame("IDATfile"=idatfiles,stringsAsFactors=FALSE)
	if(dateinfo) elist$targets$DecodeInfo <- elist$targets$ScanInfo <- rep_len("", nsamples)

#	Read bead manifest file
	if(verbose) cat("Reading manifest file", bgxfile, "... ")
	bgx <- illuminaio::readBGX(bgxfile)
	if(verbose) cat("Done\n")
	nregprobes <- nrow(bgx$probes)
	nctrlprobes <- nrow(bgx$control)
	nprobes <- nregprobes+nctrlprobes

#	Assemble gene annotation
	elist$genes <- rbind(bgx$probes[,c("Probe_Id","Array_Address_Id")], bgx$controls[,c("Probe_Id","Array_Address_Id")])

#	Set probe control status
	elist$genes$Status <- "regular"
	elist$genes$Status[(nregprobes+1):nprobes] <- bgx$controls[,"Reporter_Group_Name"]

#	Add optional annotation columns
	if(!is.null(annotation) && !is.na(annotation)) {
		annotation <- as.character(annotation)
		annotation <- intersect(names(bgx$probes),annotation)
	}
	if(length(annotation)) {
		ac <- annotation %in% names(bgx$controls)
		for (i in 1:length(annotation)) {
			elist$genes[[annotation[i]]] <- NA_character_
			elist$genes[[annotation[i]]][1:nregprobes] <- bgx$probes[[annotation[i]]]
			if(ac[i]) elist$genes[[annotation[i]]][(nregprobes+1L):nprobes] <- bgx$controls[[annotation[i]]]
		}
	}	

#	Initalize expression matrices
	elist$E <- matrix(0, nprobes, nsamples)
	elist$E[] <- NA
	colnames(elist$E) <- removeExt(idatfiles)
	rownames(elist$E) <- elist$genes[,"Array_Address_Id"]	
	elist$other$STDEV <- elist$other$NumBeads <- elist$E

#	Read IDAT files
	for(j in 1:nsamples) {
		if(verbose) cat("\t", idatfiles[j], "... ")
		tmp <- illuminaio::readIDAT(idatfiles[j])
		if(verbose) cat("Done\n")
		ind <- match(elist$genes$Array_Address_Id, tmp$Quants$IllumicodeBinData)

#		Check for whether values are available for all probes
		if(anyNA(ind)) {
			nna <- sum(is.na(ind))
			if(nna > tolerance)
				stop("Can't match all ids in manifest with those in idat file ", idatfiles[i], "\n", sum(is.na(ind)), 
				     " missing - please check that you have the right files, or consider setting \'tolerance\'=", sum(is.na(ind)))
			i <- which(!is.na(ind))
			ind <- ind[i]
			elist$E[i,j] <- tmp$Quants$MeanBinData[ind]
			elist$other$STDEV[i,j] <- tmp$Quants$DevBinData[ind]
			elist$other$NumBeads[i,j] <- tmp$Quants$NumGoodBeadsBinData[ind]
		} else {
			elist$E[,j] <- tmp$Quants$MeanBinData[ind]
			elist$other$STDEV[,j] <- tmp$Quants$DevBinData[ind]
			elist$other$NumBeads[,j] <- tmp$Quants$NumGoodBeadsBinData[ind]
		}

		if(dateinfo) {
			elist$targets$DecodeInfo[j] = paste(tmp$RunInfo[1,], collapse=" ")
			elist$targets$ScanInfo[j] = paste(tmp$RunInfo[2,], collapse=" ")
		}
	}

	if(verbose) cat("Finished reading data.\n")
	return(elist)
}
