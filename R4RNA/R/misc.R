## Copyright (C) 2011 Daniel Lai and Irmtraud M. Meyer (www.e-rna.org)
## Contact information: irmtraud.meyer@cantab.net

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

as.helix <- function(x, length = NULL) {
	x <- as.data.frame(x)
	if (nrow(x) == 0) {
		x <- data.frame(i = vector(), j = vector(), length = vector(),
			value = vector())
	}

	if(ncol(x) < 4) {
		stop("Expecting data.frame-like structure with at least 4 columns")
	}
	names <- colnames(x)
	names[1:4] <- c("i", "j", "length", "value")
	colnames(x) <- names
	rownames(x) <- NULL
	x$i <- as.integer(x$i)
	x$j <- as.integer(x$j)
	x$length <- as.integer(x$length)
	x$value <- as.numeric(x$value)

	if (is.null(length)) {
		length <- max(x$i, x$j)
	}

	attr(x, "length") <- length
	return(x)
}

is.helix <- function(x) {
	valid <- TRUE
	if (!is.data.frame(x)) {
		warning("Not a valid data frame")
		valid <- FALSE
	}
	if (ncol(x) < 4) {
		warning("Should have at least 4 columns")
		valid <- FALSE
	}
	if (is.null(attr(x, "length"))) {
		warning("No length attribute")
		valid <- FALSE
	}
	if (!all(lapply(x, class)[1:4] == c("integer", "integer", "integer",
			"numeric"))) {
		warning("Columns of invalid class")
		valid <- FALSE
	}
	if(!all(colnames(x)[1:4] == c("i", "j", "length", "value"))) {
		warning("Invalid column names")
		valid <- FALSE
	}
	return(valid)
}

collapseHelix <- function(helix) {
	if(!is.helix(helix)) {
		stop("Input is not a valid helix, aborting")
	} else {
		if (nrow(helix) == 0) {
			return(helix)
		}
		helix <- expandHelix(helix)
		helix <- helix[!duplicated(helix), ]
	}
	length <- attr(helix, "length")

	sums <- rowSums(helix[, 1:2])
	order <- order(helix$value, sums, helix$i)
	helix <- helix[order, ]
	sums <- sums[order]

	vend <- rep(FALSE, nrow(helix))
	vsta <- rep(FALSE, nrow(helix))

	# Change in value indicates helix end
	ends <- findInterval(unique(helix$value), helix$value)
	vsta[c(0, ends[-length(ends)]) + 1] <- TRUE
	vend[ends] <- TRUE

	# Change in basepair sum indicates helix end
	ends <- cumsum(rle(sums)$lengths)
	vend[ends] <- TRUE
	vsta[c(0, ends[-length(ends)]) + 1] <- TRUE

	# Basepair position difference != 1 indicates helix end
	ends <- which(c(diff(helix$i), 0) != 1)
	vend[ends] <- TRUE
	vsta[c(0, ends[-length(ends)]) + 1] <- TRUE

	starts <- which(vsta)
	ends <- which(vend)
	counts <- which(vend) - which(vsta) + 1
	ungap <- (helix$i[vsta] == helix$i[vend] - counts + 1) & (helix$j[vsta] == helix$j[vend] + counts - 1)

	if (length(which(ungap == FALSE)) > 0) {
		stop("This error should never occur... inform developer if it does")
	}

	output <- helix[vsta, ]
	output$length <- counts
	
	return(as.helix(output, length))
}

expandHelix <- function(helix) {
	if(!is.helix(helix)) {
		stop("Input is not a valid helix, aborting")
	} else {
		if (nrow(helix) == 0 || all(helix$length == 1)) {
			return(helix)
		}
	}
	length <- attr(helix, "length")
	
	output <- as.data.frame(matrix(NA, sum(helix$length), ncol(helix)))

	# Compute i
	diff <- rep(1, sum(helix$length))
	end_index <- cumsum(helix$length)
	start_index <- c(0, end_index[-length(end_index)]) + 1

	# Compute j
	end_pos <- helix$i + helix$length - 1
	diff[start_index] <- helix$i - c(0, end_pos[-length(end_pos)])
	output[, 1] <- cumsum(diff)
	

	diff <- rep(-1, sum(helix$length))
	end_pos <- helix$j - helix$length + 1
	diff[start_index] <- helix$j - c(0, end_pos[-length(end_pos)])
	output[, 2] <- cumsum(diff)

	output[, 3] <- 1
	
	for (i in 4:ncol(helix)) {
		output[, i] <- rep(helix[, i], helix$length)
	}

	colnames(output) <- colnames(helix)
	return(as.helix(output, length))
}

expandNumberedHelix <- function(helix) {
	number <- rep(1:nrow(helix), helix$length)
	helix <- expandHelix(helix)
	helix$number <- number
	return(helix)	
}

isConflictingHelix <- function(helix, any = TRUE) {
	if (!is.helix(helix)) {
		stop("Not a valid helix data frame, aborting...")
	}
	helices <- nrow(helix)
	if(helices == 0) {
		return (vector())
	}
	if (helices == 1) {
		return(c(FALSE))
	}
	
	helix <- expandNumberedHelix(helix)
	end <- nrow(helix)
	helix$keep <- T
	for (i in 1:(end - 1)) {
		if (helix$keep[i]) {
			kill <- which(helix$i[(i + 1):end] == helix$i[i] |
				helix$j[(i + 1):end] == helix$j[i] |
				helix$j[(i + 1):end] == helix$i[i] |
				helix$i[(i + 1):end] == helix$j[i])
			helix$keep[kill + i] <- F
		}
	}
	keep <- c()
	for (i in 1:helices) {
		if (any) {
			keep <- c(keep, all(helix$keep[which(helix$number == i)]))
		} else {
			keep <- c(keep, any(helix$keep[which(helix$number == i)]))
		}
	}	
	return(!keep)
}

isDuplicatingHelix <- function(helix, any = TRUE) {
	if (!is.helix(helix)) {
		stop("Invalid helix data frame, aborting")
	}
	helices <- nrow(helix)

	if(helices == 0) {
		return (vector());
	}

	expand <- expandNumberedHelix(helix)
	flip <- which(expand$i > expand$j)
	tmp <- expand$i[flip]
	expand$i[flip] <- expand$j[flip]
	expand$j[flip] <- tmp
	expand$keep <- !duplicated(expand[, c("i", "j")])

	if (all(helix$length == 1)) {
		return(!expand$keep)
	} else {
		indices <- lapply(lapply(1:helices, "==", expand$number), which)
		logicals <- lapply(indices, function(x, keep) { keep[x]}, expand$keep)

		if (any) {
			return(!unlist(lapply(logicals, any)))
		} else {
			return(!unlist(lapply(logicals, all)))
		}
	}
}

isOverlappingHelix <- function(helix, query, any = FALSE) {
	if(!is.helix(helix) | !is.helix(query)) {
		stop("One of the inputs was not a valid helix data frame, aborting")
	}
	helices <- nrow(helix)

	if (helices == 0) {
		return (vector());
	}
	if (nrow(query) == 0) {
		return (rep(FALSE, helices))
	}

	helix <- expandNumberedHelix(as.helix(helix[, 1:4]))
	y <- expandHelix(query)[, 1:4]
	y$length <- NULL
	y$value <- NULL
	y$hit <- T
	helix <- merge(helix, y, all.x = TRUE)
	helix$hit[is.na(helix$hit)] <- F
	overlap <- c()
	for (i in 1:helices) {
		if (any) {
			overlap <- c(overlap, any(helix$hit[which(helix$number == i)]))
		} else {
			overlap <- c(overlap, all(helix$hit[which(helix$number == i)]))
		}
	}
	return(overlap)	
}

logseq <- function(from, to, length.out) {
	if (missing(length.out)) {
		length.out <- log10(max(from, to) / min(from, to)) + 1
	}
	return(exp(seq(log(from), log(to), length.out = length.out)))
}

logfloor <- function(x) {
	return(10 ** floor(log10(x)))
}

logceiling <- function(x) {
	return(10 ** ceiling(log10(x)))
}

basepairFrequency <- function(helix) {
	if (!is.helix(helix)) {
		stop("Not a valid helix data frame, aborting")
	}
	if (nrow(helix) == 0) {
		return(helix)
	} else {
		basepairs <- expandHelix(helix)[, c("i", "j")]
	}
	counts <- table(paste(basepairs$i, basepairs$j))
	pos <- as.integer(unlist(strsplit(names(counts), " ")))
	odds <- seq(1, length(pos), by = 2)
	output <- data.frame(i = pos[odds], j = pos[odds + 1], length = 1,
		value = as.integer(counts))
	output <- output[order(-output$value, output$i, output$j), ]
	return(as.helix(output, attr(helix, "length")))
}

# Returns a logical row each row of the matrix, TRUE of pseudoknotted, else false
# Assumes these are  helices of LENGTH 1
is_pseudoknotted <- function(row, matrix) {
	return((row$i <= matrix$i & matrix$i <= row$j & row$j <= matrix$j) |
	(matrix$i <= row$i & row$i <= matrix$j & matrix$j <= row$j))
}

# Adds a group column, in which helices of each numbered group are non-pseudoknotted
unknottedGroups <- function(helix) {
	if (!is.helix(helix)) {
		stop("Invalid input")
	}
	if (any(helix$length > 1)) {
		warning("Expanding helix to basepairs")
		helix <- expandHelix(helix)
	}
	if (nrow(helix) == 0) {
		stop("No basepairs in helix")
	}
	group <- rep(0, nrow(helix));
 	group[1] <- 1
	if (nrow(helix) == 1) {
		return(group)
	}
	for (i in 2:nrow(helix)) {
		test <- 1
		while (any(is_pseudoknotted(helix[i, ], helix[which(group == test), ]))) {
			test <- test + 1
		}
		group[i] <- test
	}
	return(group)
}

