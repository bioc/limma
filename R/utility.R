#  UTILITY FUNCTIONS

.matvec <- function(M,v) {
#	Multiply the columns of matrix by the elements of a vector,
#	i.e., compute M %*% diag(v)
#	Gordon Smyth
#	5 July 1999
#
	v <- as.vector(v)
	M <- as.matrix(M)
	if(length(v)!=dim(M)[2]) stop("Dimensions do not match")
	t(v * t(M))
}

.vecmat <- function(v,M) {
#	Multiply the rows of matrix by the elements of a vector,
#	i.e., compute diag(v) %*% M
#	Gordon Smyth
#	5 July 1999
#
	v <- as.vector(v)
	M <- as.matrix(M)
	if(length(v)!=dim(M)[1]) stop("Dimensions do not match")
	v * M
}

isNumeric <- function(x) {
#	Test for numeric argument or data.frame with numeric columns
#	Gordon Smyth
#	12 April 2003

	is.numeric(x) || (is.data.frame(x) && length(x)>0 && all(unlist(lapply(x,is.numeric))))
}

blockDiag <- function(...)
#	Block diagonal matrix
#	Gordon Smyth
#	8 Feb 2004
{
	e <- list(...)
	d <- matrix(unlist(lapply(e,dim)),ncol=2,byrow=TRUE)
	if(nrow(d) != length(e)) stop("all arguments must be matrices")
	dsum <- apply(d,2,sum)
	z <- array(0,dsum)
	dimnames(z) <- list(character(dsum[1]),character(dsum[2]))
	coord <- c(0,0)
	for (i in 1:length(e)) {
		coord1 <- coord[1]+1:d[i,1]
		coord2 <- coord[2]+1:d[i,2]
		z[coord1,coord2] <- e[[i]]
		rownames(z)[coord1] <- rownames(e[[i]],do.NULL=FALSE,prefix="")
		colnames(z)[coord2] <- colnames(e[[i]],do.NULL=FALSE,prefix="")
		coord <- coord+d[i,]
	}
	z
}

helpMethods <- function(genericFunction) {
#	Prompt user for help topics on methods for generic function
#	Gordon Smyth
#	21 April 2003.  Last revised 24 June 2014.

	objectclass <- class(genericFunction)
 	if(objectclass != "standardGeneric") {
		if(objectclass == "character" && isGeneric(genericFunction))
			genericFunction <- getGeneric(genericFunction)
		else {
			cat("Not a generic function\n")
			return(invisible())
		}
	}
	functionname <- genericFunction@generic
	methodnames <- names(findMethods(genericFunction))
	nmethods <- length(methodnames)
	if(nmethods == 0) {
		cat("No available methods\n")
		return(invisible())
	}
	aliasnames <- paste(functionname,",",methodnames,"-method",sep="")
	for (i in 1:nmethods) cat(i,": ",aliasnames[i],"\n",sep="")
	cat("Type number to choose help topic: ")
	n <- as.integer(readline())
	if(!is.na(n) && n > 0 && n <= nmethods)
		eval(parse(text=paste("help(\"",aliasnames[n],"\")",sep="")))
	else {
	 	cat("No topic chosen\n")
	 	return(invisible())
	}
}

limmaUsersGuide <- function(view=TRUE)
#	Find and optionally view limma User's Guide
#	Gordon Smyth
#	25 Oct 2004.
{
	f <- system.file("doc","usersguide.pdf",package="limma")
	if(view) {
		if(.Platform$OS.type == "windows") 
			shell.exec(f)
		else
			system(paste(Sys.getenv("R_PDFVIEWER"),f,"&"))
	}
	return(f)
}

changeLog <- function(n=30L, package="limma")
#	Write first n lines of package change log or NEWS file.
#	Originally written for limma package, but works for any package
#	with a change log or NEWS file.
#	Gordon Smyth
#	Created 20 Sep 2004.  Last modified 14 Apr 2020.
{
#	Allow users to supply package name as a single argument
	if(is.character(n) && nargs()==1L) {
		package <- n
		n <- 30L
	}

#	If package isn't installed there's no point in looking for the change log.
	if(!(package %in% .packages(all.available=TRUE))) {
		message("Package ",package," is not installed.")
		return(invisible())
	}

#	Look for change log file in several possible naming conventions.
	Path <- system.file("doc","changelog.txt",package=package)
	FileExists <- file.exists(Path)
	if(!FileExists) {
		Path <- system.file("doc","ChangeLog",package=package)
		FileExists <- file.exists(Path)
	}
	if(!FileExists) {
		Path <- system.file("changelog.txt",package=package)
		FileExists <- file.exists(Path)
	}
	if(!FileExists) {
		Path <- system.file("ChangeLog",package=package)
		FileExists <- file.exists(Path)
	}
	if(!FileExists) {
		Path <- system.file("CHANGELOG",package=package)
		FileExists <- file.exists(Path)
	}
	if(!FileExists) {
		Path <- system.file("NEWS",package=package)
		FileExists <- file.exists(Path)
	}
	if(!FileExists) {
		Path <- system.file("NEWS.md",package=package)
		FileExists <- file.exists(Path)
	}

#	If found, write change log to R session.
	if(FileExists)
		writeLines(readLines(Path,n=n))
	else {
#		NEWS.Rd is a structured file rather than a change log.
#		If we find that file but no change log, then suggest use of news().
		if(file.exists(system.file("NEWS.Rd",package=package))) {
			message("No change log file can be found. Use news(package='",package,"') instead.")
		} else {
			message("No change log or NEWS file can be found")
		}
	}
}
