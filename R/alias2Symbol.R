##  ALIAS2SYMBOL.R

alias2Symbol <- function(alias,species="Hs",expand.symbols=FALSE)
#  Convert a set of alias names to official gene symbols
#  via Entrez Gene identifiers
#  Di Wu, Gordon Smyth and Yifang Hu
#  4 Sep 2008. Last revised 13 April 2016.
{
	alias <- as.character(alias)
	species <- match.arg(species,c("Dm","Hs","Mm","Rn"))

#	Get access to required annotation functions
	suppressPackageStartupMessages(OK <- requireNamespace("AnnotationDbi",quietly=TRUE))
	if(!OK) stop("AnnotationDbi package required but not available")

#	Get alias to symbol mappings
	switch(species,
		Hs = {
			suppressPackageStartupMessages(OK <- requireNamespace("org.Hs.eg.db",quietly=TRUE))
			if(!OK) stop("org.Hs.eg.db package required but not available")
			ALIAS2EG <- org.Hs.eg.db::org.Hs.egALIAS2EG
			SYMBOL <- org.Hs.eg.db::org.Hs.egSYMBOL
		}, Mm = {
			suppressPackageStartupMessages(OK <- requireNamespace("org.Mm.eg.db",quietly=TRUE))
			if(!OK) stop("org.Mm.eg.db package required but not available")
			ALIAS2EG <- org.Mm.eg.db::org.Mm.egALIAS2EG
			SYMBOL <- org.Mm.eg.db::org.Mm.egSYMBOL
		}, Rn = {
			suppressPackageStartupMessages(OK <- requireNamespace("org.Rn.eg.db",quietly=TRUE))
			if(!OK) stop("org.Rn.eg.db package required but not available")
			ALIAS2EG <- org.Rn.eg.db::org.Rn.egALIAS2EG
			SYMBOL <- org.Rn.eg.db::org.Rn.egSYMBOL
		}, Dm = {
			suppressPackageStartupMessages(OK <- requireNamespace("org.Dm.eg.db",quietly=TRUE))
			if(!OK) stop("org.Dm.eg.db package required but not available")
			ALIAS2EG <- org.Dm.eg.db::org.Dm.egALIAS2EG
			SYMBOL <- org.Dm.eg.db::org.Dm.egSYMBOL
		}
	)

	if(expand.symbols) {
		alias <- intersect(alias,AnnotationDbi::Rkeys(ALIAS2EG))
		eg <- AnnotationDbi::mappedLkeys(ALIAS2EG[alias])
		AnnotationDbi::mappedRkeys(SYMBOL)[eg]
	} else {
		isSymbol <- alias %in% AnnotationDbi::Rkeys(SYMBOL) 
		alias2 <- intersect(alias[!isSymbol],AnnotationDbi::Rkeys(ALIAS2EG))
		eg <- AnnotationDbi::mappedLkeys(ALIAS2EG[alias2])
		c(alias[isSymbol],AnnotationDbi::mappedRkeys(SYMBOL[eg]))
	}
}

alias2SymbolTable <- function(alias,species="Hs")
#  Convert a vector of alias names to the vector of corresponding official gene symbols
#  via Entrez Gene identifiers
#  Di Wu, Gordon Smyth and Yifang Hu
#  Created 3 Sep 2009.  Last modified 13 April 2016.
{
	alias <- as.character(alias)
	species <- match.arg(species,c("Dm","Hs","Mm","Rn"))

#	Get access to required annotation functions
	suppressPackageStartupMessages(OK <- requireNamespace("AnnotationDbi",quietly=TRUE))
	if(!OK) stop("AnnotationDbi package required but not available")

#	Get alias to symbol mappings
	switch(species,
		Hs = {
			suppressPackageStartupMessages(OK <- requireNamespace("org.Hs.eg.db",quietly=TRUE))
			if(!OK) stop("org.Hs.eg.db package required but not available")
			ALIAS2EG <- org.Hs.eg.db::org.Hs.egALIAS2EG
			SYMBOL <- org.Hs.eg.db::org.Hs.egSYMBOL
		}, Mm = {
			suppressPackageStartupMessages(OK <- requireNamespace("org.Mm.eg.db",quietly=TRUE))
			if(!OK) stop("org.Mm.eg.db package required but not available")
			ALIAS2EG <- org.Mm.eg.db::org.Mm.egALIAS2EG
			SYMBOL <- org.Mm.eg.db::org.Mm.egSYMBOL
		}, Rn = {
			suppressPackageStartupMessages(OK <- requireNamespace("org.Rn.eg.db",quietly=TRUE))
			if(!OK) stop("org.Rn.eg.db package required but not available")
			ALIAS2EG <- org.Rn.eg.db::org.Rn.egALIAS2EG
			SYMBOL <- org.Rn.eg.db::org.Rn.egSYMBOL
		}, Dm = {
			suppressPackageStartupMessages(OK <- requireNamespace("org.Dm.eg.db",quietly=TRUE))
			if(!OK) stop("org.Dm.eg.db package required but not available")
			ALIAS2EG <- org.Dm.eg.db::org.Dm.egALIAS2EG
			SYMBOL <- org.Dm.eg.db::org.Dm.egSYMBOL
		}
	)

	isSymbol <- alias %in% AnnotationDbi::Rkeys(SYMBOL)
	Symbol <- alias
	Symbol[!isSymbol] <- NA

	OtherAliases <- alias[!isSymbol]
	isAlias <- OtherAliases %in% AnnotationDbi::Rkeys(ALIAS2EG)
	if(!any(isAlias)) return(Symbol)
	OtherAliases <- OtherAliases[isAlias]

	AliasTbl <- AnnotationDbi::toTable(ALIAS2EG[OtherAliases])
	if(anyDuplicated(AliasTbl$alias_symbol)) warning("Multiple symbols ignored for one or more aliases")
	SymbolTbl <- AnnotationDbi::toTable(SYMBOL[AliasTbl$gene_id])
	m <- match(OtherAliases,AliasTbl$alias_symbol)
	GeneID <- AliasTbl$gene_id[m]
	m <- match(GeneID,SymbolTbl$gene_id)
	Symbol[!isSymbol][isAlias] <- SymbolTbl$symbol[m]
	Symbol
}
