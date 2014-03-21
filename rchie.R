#!/usr/bin/env Rscript

## Copyright (C) 2012 Daniel Lai and Irmtraud M. Meyer (www.e-rna.org)
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

library(R4RNA) # Use with R4RNA v0.1.3
library(RColorBrewer)
suppressPackageStartupMessages(library(optparse))

option_list <- list(
	make_option("--format1", default = "helix", help = paste(
		'Input format, choose:',
		'"vienna", "connect", "bpseq", "helix" [default "%default"]',
		sep = "\n\t\t")),
	make_option("--decreasing1", type = "integer", help =
		"1 to sort basepairs by value in decreasing order, or 0 for increasing"),
	make_option("--minimum1", type = "double", help =
		'Ignore basepairs with values smaller than this value'),
	make_option("--maximum1", type = "double", help =
		'Ignore basepairs with values greater than this value'),
	make_option("--colour1", type = "character", help = paste(
		'A quoted comma-separated string denoting the colours to use for groups.', 
		'Expects hex, but will probably take colour names',
		'e.g. \""#1A9850,#91CF60,#D9EF8B\"', sep = "\n\t\t")),
	make_option("--palette1", type = "character", help = paste(
		'Ignored when --colour1 is set.  One of the following colour palettes:',
			'BrBG, PiYG, PRGn, PuOr, RdBu, RdGy, RdYlBu, RdYlGn, Spectral,',
			'Accent, Dark2, Paired, Pastel1, Pastel2, Set1, Set2, Set3, Blues,',
			'BuGn, BuPu, GnBu, Greens, Greys, Oranges, OrRd, PuBu, PuBuGn,',
			'PuRd, Purples, RdPu, Reds, YlGn, YlGnBu, YlOrBr, YlOrRd', sep =
			"\n\t\t\t")),
	make_option("--revpal1", "store_true", default = FALSE,
		help = "Reverses colours selected by palette1"),
	make_option("--group1", type = "integer", help =
		'Ignored when --colour1 is set.  An integer number of groups'),
	make_option("--rule1", default = 0, help = paste(
		'Specifies how to create basepair groups:',
		'0: Equal number of basepairs',
		'1: Equidistance value intervals',
		'2: Alias for rule 4 (requires --value1)',
		'3: User set number of basepairs (requires --value1)',
		'4: User set value intervals (requires --value1)',
		'5: Non-pseudoknotted groups',
		'6: Dot-bracket bracket type (requires --format vienna)',
		'7: Basepair covariation value (requires --msa)',
		'8: Basepair conservation value (requires --msa)',
		'9: Percent canonical basepairs (requires --msa)', sep = "\n\t\t\t")),
	make_option("--value1", type = "character", help = paste(
		'A quoted comma-separated string denoting the values to use for',
		'grouping. If the first symbol is a negative sign, it must be preceded',
		'with a backslash. e.g. \"\\-1.5,-1,0,1,2,3\"', sep = "\n\t\t")),
	make_option("--legend1", "store_true", default = FALSE,
		help = 'Adds legend to first structure\n\n'),

	make_option("--reference", "store_true", default = FALSE, help = 'Treat second structure as reference of an overlapping structure diagram'),
	make_option("--format2", default = "helix", help = 'Format for second file [default "%default"]'),
	make_option("--decreasing2", help = "Decreasing for second file"),
	make_option("--minimum2", type = "double", help = 'Minimum for second file'),
	make_option("--maximum2", type = "double", help = 'Maximum for second file'),
	make_option("--colour2", type = "character", help = 'Colour for second file'),
	make_option("--palette2", type = "character", help = 'Palette for second file'),
	make_option("--revpal2", "store_true", default = FALSE, help = "Reverses colours selected by palette2"),
	make_option("--group2", type = "integer", help = 'Group for second file'),
	make_option("--rule2", default = 0, help = 'Rule for second file.'),
	make_option("--value2", type = "character", help = 'Value for second file'),
	make_option("--legend2", "store_true", default = FALSE, help = 'Adds legend to second structure\n\n'),

	make_option("--msafile", default = "", help =
		paste('FASTA-format multiple sequence alignment of', 
		'structure, triggers covariance plotting')),
	make_option("--msacol", type = "character", help = paste(
		'Quoted-comma separated string for multiple sequence alignment colours', 
		'Corresponding to: conservation, covariation, one-sided covariation,',
		'\tinvalid pair, upaired, gap, and ambiguous base/pair.',
		'If --basecol is specified, then colours correspond to:',
		'\tA, U/T, G, C, Gap, and Ambigous',
		sep = "\n\t\t\t")),
	make_option("--msagrid", "store_true", default = FALSE,
		help = 'Plots multiple sequence alignment as grid instead of lines'),
	make_option("--msatext", "store_true", default = FALSE,
		help = 'Displays nucleotide base (only works with --msagrid is set)'),
	make_option("--msaspecies", "store_true", default = FALSE,
		help = 'Display FASTA headers left of each sequence'),
	make_option("--basecol", "store_true", default = FALSE,
		help = 'Colour bases by nucleotide instead of basepair status'),
	make_option("--legend3", "store_true", default = FALSE,
		help = 'Adds legend for multiple sequence alignment\n\n'),

	make_option("--noscale", "store_true", default = FALSE,
		help = 'Suppress basepair scale in output'),
	make_option("--quiet", "store_true", default = FALSE,
		help = 'Suppress processing information'),
	make_option("--pdf", "store_true", default = FALSE,
		help = 'Also outputs in pdf format'),
	make_option("--output", default = "rchie.png",
		help = 'Output file name [default "%default"]'),
	make_option("--show", "store_true", default = FALSE,
		help = 'Show output file, requires acroread and eog on command line')
)

parser <- OptionParser(usage = "%prog [options] input.txt",
	option_list = option_list)

args <- parse_args(parser, positional_arguments = TRUE);
opt <- args$options
args <- args$args
verbose <- !opt$quiet

# print(opt)

# i == 1 for first, 2 for second
parseInput <- function(i) {

	# Read input
	if (verbose) { message("Reading input file") }
	tmp <- opt[paste("format", i, sep = "")]
	if (is.na(args[i])) {
		stop("ERROR: No input file found, use --help for usage")
	}
	if (verbose) { message(paste("\tFile:", args[i], "\n\tFormat:", tmp)) }
	if (tmp %in% c("vienna", "connect", "bpseq", "helix")) {
		if (tmp == "vienna") {
			input <- readVienna(args[i])
		}
		if (tmp == "connect") {
			input <- readConnect(args[i])
		}
		if (tmp == "bpseq") {
			input <- readBpseq(args[i])
		}
		if (tmp == "helix") {
			input <- readHelix(args[i])
		}
	} else {
		stop(paste("ERROR: Invalid file format:", tmp))
	}

	# Sort input
	if (verbose) { message("Sorting input file") }
	if (paste("decreasing", i, sep = "") %in% names(opt)) {
		tmp <- opt[paste("decreasing", i, sep = "")]
		if (verbose) { message(paste("\tUser-defined decreasing:", tmp)) }
		if (tmp == 1 || tmp == 0) {
			decreasing <- as.numeric(tmp)
			input <- input[order(input$value, decreasing = decreasing), ]
		} else {
			stop("ERROR: Expecting 0 or 1 for decreasing\n")
		}
	} else {
		# Inferring directionality
		if (any(is.na(input$value)) || input$value[1] > input$value[nrow(input)]) {
			decreasing <- 1
		} else {
			decreasing <- 0
		}
		if (verbose) { message(paste("\tInferred decreasing:", decreasing)) }
	}

	# Filter values
	if (verbose) { message("Filtering input file by values") }	
	if (paste("minimum", i, sep = "") %in% names(opt)) {
		tmp <- opt[paste("minimum", i, sep = "")]
		if (verbose) { message(paste("\tRemoving all values less than:", tmp)) }
		if (is.na(tmp)) { 
			stop("ERROR: Invalid minimum value\n")
		} else {
			input <- input[which(input$value >= tmp), ]
		}
	}

	if (paste("maximum", i , sep = "") %in% names(opt)) {
		tmp <- opt[paste("maximum", i, sep = "")]
		if (verbose) { message(paste("\tRemoving all values greater than:", tmp)) }
		if (is.na(tmp)) {
			stop("ERROR: Invalid maximum value\n")
		} else {
			input <- input[which(input$value <= tmp), ]
		}
	}

	if (nrow(input) == 0) {
		stop("ERROR: No basepairs to plot, too much filtering or attempted to filter valueless basepairs?\n")
	}

	# Expand and clean input
	if (verbose) { message("Expanding helices into basepairs") }
	input <- expandHelix(input)
	if (verbose) { message("Removing duplicate basepairs") }
	input <- input[which(!isDuplicatingHelix(input)), ]

	# Parse colours
	colour <- c()
	if (verbose) { message("Parsing colour from input") }
	if (paste("colour", i , sep = "") %in% names(opt)) {
		if (verbose) { message("\tParsing --colour option") }
		tmp <- as.character(opt[paste("colour", i , sep = "")])
		colour <- unlist(strsplit(tmp, ","))
		if (verbose) { message(paste("\t", length(colour), " colours found: ", paste(colour, collapse = " "), sep = "")) }
	} else {
		if (paste("palette", i , sep = "") %in% names(opt)) {
			if (verbose) { message("\tParsing --palette option") }
			tmp <- as.character(opt[paste("palette", i, sep = "")])
			if (tmp %in% rownames(brewer.pal.info)) {
				if(paste("group", i, sep = "") %in% names(opt)) {
					group <- as.integer(opt[paste("group", i, sep = "")])
					if (is.na(group)) {
						stop("ERROR: --group expects single integer argument")
					}
					if (group > 8) {
						message("WARNING: Capping group at 8")
						group <- 8
					}
				} else {
					stop("ERROR: palette options requires --group")
				}
				if (verbose) { message(paste("\tGenerating", group, "colours from", tmp, "palette")) }
				if (group < 3) {
					colour <- brewer.pal(3, tmp)
					if (group == 1) {
						colour <- colour[3]
					} else {
						colour <- colour[c(1, 3)]
					}
				} else {
					colour <- brewer.pal(group, tmp)
				}
				if (verbose) { message(paste("\t", length(colour), " colours generated: ", paste(colour, collapse = " "), sep = "")) }
				if (paste("revpal", i, sep = "") %in% names(opt)) {
					colour <- rev(colour)
					if (verbose) { message("Reversing palette colours")}
				}
			} else {
				stop("ERROR: Invalid palette name")
			}
		} else {
			if (verbose) { message("\tNo colour detected, defaulting to black") }
		}
	}

	# Assign colours + create legend according rule
	legend <- list()
	title <- NULL
	if (length(colour) == 0) {
		if (verbose) { message("No colous specified, making everything black") }
		col <- colourByCount(input, "#000000")
	} else {
		rule <- opt[paste("rule", i, sep = "")]
		if (verbose) { message(paste("Colouring basepairs according to rule:", rule)) }
		if (!(rule %in% c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9))) {
			stop(paste("ERROR: Invalid argument to rule option:", rule))
		}
		if (rule == 0) {
			if (verbose) { message("\tColouring by groups of equal size") }
			col <- colourByCount(input, unlist(colour))
			title <- "Basepairs"
		}
		if (rule == 1) {
			if (verbose) { message("\tColouring by groups of equidistant values") }
			col <- colourByValue(input, colour)	
		}
		if (rule == 2 || rule == 3 || rule == 4) {
			value <- as.numeric(unlist(strsplit(gsub("\\\\", "", opt[paste("value", i , sep = "")]), ",")))
			if (verbose) { message("\tColouring by groups by user-defined value") }
			if (verbose) { message(paste("\t", length(value), " values: ", paste(value, collapse = " "), sep = "")) }
			if (length(value) == length(colour) && !any(is.na(value))) {
				if (rule == 3) {
					if (verbose) { message("\tUsing values as basepairs in each group") }
					if (all(value >= 0)) {
						col <- colourByCount(input, colour, value)
						title <- "Basepairs"
					} else {
						stop("ERROR: Invalid integer helix count, negative number?")				
					}
				} else { # rule == 4 || rule == 2
					if (verbose) { message("\tUsing values as range boundaries") }
					value <- c(value, ifelse(decreasing, Inf, -Inf))
					col <- colourByValue(input, colour, value)
				}
			}  else {
				stop("ERROR: Mismatch between valid value and colour lengths")
			}
		
		}
		if (rule == 5) {
			if (verbose) { message("\tColouring by unknotted groups") }
			col <- colourByUnknottedGroups(input, colour, get = FALSE)
			title <- "Unknotted Basepairs"
		}
		if (rule == 6) {
			if (verbose) { message("\tColouring by bracket type") }
			if (opt[paste("format", i, sep = "")] != "vienna") {
				stop("ERROR: Must be vienna format input for this rule")
			}
			tmp <- readVienna(args[i], palette = colour)
			tmp <- merge(input, tmp, sort = FALSE)
			col <- tmp$col
			tmp <- table(col)[colour]
			tmp[is.na(tmp)] <- 0
			attr(col, "legend") <- paste(tmp, "/", sum(tmp), sep = "")
			attr(col, "fill") <- colour
		}
		if (rule == 7 || rule == 8 || rule == 9) {
			if (verbose) { message("\tColouring by multiple sequence alignment information") }
			if (!cov) {
				stop("ERROR: Rule requires multiple sequence alignment")
			}
			if (rule == 7) {
				if (verbose) { message("\tColouring by covariation") }
				col <- colourByCovariation(input, msa, colour)
				title <- "Covariation"
			}
			if (rule == 8) {
				if (verbose) { message("\tColouring by conservation") }
				col <- colourByConservation(input, msa, colour)
				title <- "Conservation"
			}
			if (rule == 9) {
				if (verbose) { message("\tColouring by percent canonical basepairs") }
				col <- colourByCanonical(input, msa, colour)
				title <- "% Canonical Basepair"
			}
		}
	}
	if (verbose) { message("\tColouring basepairs") }
	input$col <- NA
	input$col[1:min(length(col), nrow(input))] <- col
	legend <- list(attr(col, "legend"), attr(col, "fill"), title)
	return(list(input, legend))
}

# Read MSA
cov <- FALSE
if (opt$msafile != "") {
	if (verbose) { message("\n*****PROCESSING MSA FILE*****") }
	msa <- try(readFasta(opt$msafile, filter = TRUE))
	if(class(msa) == "try-error") {
		stop("ERROR: Invalid FASTA file")
	} else {
		if (length(msa) < 1) {
			stop("ERROR: No valid nucleotide sequences detected")
		}
		cov <- TRUE
	}
	
	msacol <- NA
	if ("msacol" %in% names(opt)) {
		if (verbose) { message("\tParsing --msacol option") }
		tmp <- as.character(opt["msacol"])
		msacol <- unlist(strsplit(tmp, ","))
		if (verbose) { message(paste("\t", length(msacol), " colours found: ", paste(msacol, collapse = " "), sep = "")) }
	}
	msatext <- opt$msatext
	basecol <- opt$basecol
	msalegend <- opt$legend3
	msaspecies <- opt$msaspecies
	species <- 0
	if (msaspecies) {
		species <- max(nchar(names(msa))) * 0.8438 + 2
	}
}

input <- list()
value <- list()

# Read and parse file parameters
if (verbose) { message("\n*****PROCESSING FIRST STRUCTURE*****") }
output <- parseInput(1)
input[[1]] <- output[[1]]
input[[1]] <- input[[1]][which(!is.na(input[[1]]$col)), ]
value[[1]] <- output[[2]]

if (nrow(input[[1]]) == 0) {
	warning("WARNING: No basepairs in first structure to plot")
}

sec <- FALSE
if (length(args) > 1) {
	if (verbose) { message("\n*****PROCESSING SECOND STRUCTURE*****") }
	sec <- TRUE
	output <- parseInput(2)
	input[[2]] <- output[[1]]
	input[[2]] <- input[[2]][which(!is.na(input[[2]]$col)), ]
	value[[2]] <- output[[2]]

	if (nrow(input[[2]]) == 0) {
		warning("WARNING: No basepairs in second structure to plot")
		sec <- FALSE
	}
}

if (verbose) { message("\n*****CREATING PLOT*****") }
lwd <- 5
cex <- 1.7
scale.lwd <- 2
text.cex = 0.4
scale <- !opt$noscale

if (opt$pdf) {
	png <- NA
	pdf <- opt$output
} else {
	png <- opt$output
	pdf <- NA
}


width <- max(attr(input[[1]], "length"))

if (opt$reference) {
	opt$legend2 <- FALSE
}

pad <- c(4, 4, 4, 4)

if (opt$legend1 || opt$legend2) {
	pad[4] <- 25
}

textheight <- 2.244812 # pre-compute with strheight

if (opt$legend3) {
	pad[1] <- textheight * 2
}

if (sec) { # Two structure
	width <- max(width, attr(input[[2]], "length"))
	if (cov) {
		if (opt$reference) {
			if (verbose) { message("\tOverlapping Covariance") }
			plotOverlapCovariance(input[[1]], input[[2]], msa, png = png, pdf =
				pdf, lwd = lwd, cex = cex, scale.lwd = scale.lwd, scale = scale,
				palette = msacol, text = msatext, text.cex = text.cex,
				base.colour = basecol, legend = msalegend, species = species,
				grid = opt$msagrid, pad = pad)
		} else {
			if (verbose) { message("\tDouble Covariance") }
			plotDoubleCovariance(input[[1]], input[[2]], msa, png = png, pdf =
				pdf, lwd = lwd, cex = cex, scale.lwd = scale.lwd, scale = scale,
				palette = msacol, text = msatext, text.cex = text.cex,
				base.colour = basecol, legend = msalegend, species = species,
				grid = opt$msagrid, pad = pad)
		}
	} else {
		if (opt$reference) {
			if (verbose) { message("\tOverlapping Structures") }
			plotOverlapHelix(input[[1]], input[[2]], png = png, pdf = pdf, lwd =
				lwd, cex = cex, scale.lwd = scale.lwd, scale = scale, line =
				TRUE, arrow = TRUE, pad = pad)
		} else {
			if (verbose) { message("\tDouble Structures") }
			plotDoubleHelix(input[[1]], input[[2]], png = png, pdf = pdf, lwd =
				lwd, cex = cex, scale.lwd = scale.lwd, scale = scale, line =
				TRUE, arrow = TRUE, pad = pad)
		}
	}
	usr <- par("usr")
	if (opt$legend1) {
		legend(width, usr[4], legend = value[[1]][[1]], fill = value[[1]][[2]], border = NA, bty = "n", title = value[[1]][[3]], xjust = 0, title.adj = 0.5)
	}
	if (opt$legend2) {
		legend(width, usr[3], legend = value[[2]][[1]], fill = value[[2]][[2]], border = NA, bty = "n", yjust = 0, title = value[[2]][[3]], xjust = 0, title.adj = 0.5)
	}
} else { # One structure
	if (cov) {
		if (verbose) { message("\tSingle Covariance") }
		plotCovariance(msa, input[[1]], png = png, pdf = pdf, grid =
			opt$msagrid, lwd = lwd, cex = cex, scale.lwd = scale.lwd, scale =
			scale, palette = msacol, text =
			msatext, text.cex = text.cex, base.colour = basecol, legend =
			msalegend, species = species, pad = pad)
	} else {
		if (verbose) { message("\tSingle Structure") }
		plotHelix(input[[1]], line = TRUE, arrow = TRUE, png = png, pdf = pdf,
			lwd = lwd, cex = cex, scale.lwd = scale.lwd, scale = scale,
			pad = pad)
	}
	usr <- par("usr")
	if (opt$legend1) {
		legend(width, usr[4], legend = value[[1]][[1]], fill = value[[1]][[2]], border = NA, bty = "n", title = value[[1]][[3]], xjust = 0, title.adj = 0.5)
	}
}

dim <- par("usr")
text(dim[1] + 1, dim[3] + 1, "www.e-rna.org", adj = 0, cex = 0.75, xpd = TRUE)
garbage <- dev.off()
if (verbose) {
	message("\nPrinting warnings")
	warnings()
}

if (opt$show) {
	if (opt$pdf) {
		system(paste("acroread", opt$output, "&"))	
	} else {
		system(paste("eog", opt$output, "&"))
	}
}
