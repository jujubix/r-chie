\name{Basepair/Helix Conversion}
\alias{expandHelix}
\alias{collapseHelix}
\title{
	Expand or collapse helices to and from basepairs
}
\description{
	Given a helix data frame, expands a helix of arbitrary length into helices
	of length 1 (i.e. basepairs).  Also does the reverse operation of clustering
	consecutive basepairs (or helices), and merging/collapsing them into a
	single helix.
}
\usage{
	expandHelix(helix)
	collapseHelix(helix)
}
\arguments{
	\item{helix}{
		A helix data frame.
	}
}
\details{
	During the expansion, basepairs expanded from a single helix will all be
	assigned the value of the originating helix (the same goes for all other
	columns besides i, j, and length).  During collapsing, only
	helices/basepairs of equal value will be grouped together.  The ordering
	of collapsed helices returned will be sorted by value (increasing order).
	For any other columns besides i, j, length and value, values will be
	obtained from the corresponding columns of the outer most basepair.
}
\value{
	Returns a helix data frame.
}
\author{
	Daniel Lai
}

\examples{
	# Create helix data frame
	helix <- data.frame(2, 8, 3, 0.5)
	helix[2, ] <- c(5, 15, 4, -0.5)
	helix <- as.helix(helix)
	helix$colour <- c("red", "blue")

	# Before expansion
	print(helix)
	# After expansion
	print(expanded <- expandHelix(helix))
	# Collapse back (sorted by value)
	print(collapseHelix(expanded))
}
\keyword{manip}

