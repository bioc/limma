chooseLowessSpan <- function(n=1000, small.n=50, min.span=0.3, power=1/3)
#	Choose an optimal span for lowess smoothing of variance trends as a function of the number of points.
#	Large spans are used for small datasets and smaller spans for larger datasets.
#	Gordon Smyth
#	Created 12 July 2020. Last modified 3 April 2024.
{
	pmin( min.span + (1-min.span) * (small.n/n)^power, 1)
}
