#' Draw a 2D NMR Spectrum
#' 
#' This function simulates 2D NMR spectra.  Only 1st order coupling can be
#' handled, and there is currently no capacity for doublet of doublets and
#' other such peaks.  The field strength of the "instrument" is taken into
#' account.
#' 
#' @param peaks A data frame with the following columns: delta, mult
#' (multiplicity), J, area, pw.  Multiplicity should be given by a number, so
#' use 2 for a doublet.  J is in Hz (use 0 for singlets).  pw is the peak width
#' at half-height in Hz.
#'
#' @param x.range A numeric vector of length 2 giving the ppm range desired.
#' Must be increasing.
#'
#' @param MHz Integer.  The operating frequency of the instrument, in MHz.
#'
#' @param ppHz Points per Hz: The number of data points per Hz to use in
#' calculating the spectrum.  See \code{\link{makeSpec}} for details.
#'
#' @param M An adjacency matrix indicating which peaks are coupled.
#' The order of rows and columns must be the same as in \code{peaks}.
#'
#' @param type The type of 2D spectrum desired.  One of \code{c("COSY", "TOCSY")}.
#'
#' @param levels A vector of levels for the contour plot.  Must be in [0...1].
#'
#' @param \ldots Parameters to be passed to the plotting function.
#'
#' @return Returns a matrix.
#'
#' @author Bryan A. Hanson, DePauw University. \email{hanson@@depauw.edu}
#' @seealso \code{\link{makeSpec}}
#' @keywords utilities
#' @export
#' @examples
#' # ethyl 2-ethyl-3-oxobutyrate
#' peaks1 <- data.frame(
#' 	delta = c(4.20, 3.34, 2.23, 1.88, 1.28, 0.94),
#' 	mult = c(4, 3, 1, 5, 3, 3),
#' 	J = c(14, 14, 0, 14, 14, 14),
#' 	area = c(2, 1, 3, 2, 3, 3),
#' 	pw = c(2, 2, 2, 2, 2, 2))
#' 
#' #              A, B, C, D, E, F
#' AM <- matrix(c(0, 0, 0, 0, 1, 0,  # A
#'                0, 0, 0, 1, 0, 0,  # B
#'                0, 0, 0, 0, 0, 0,  # C
#'                0, 1, 0, 0, 0, 1,  # D
#'                1, 0, 0, 0, 0, 0,  # E
#'                0, 0, 0, 1, 0, 0), # F
#' 			   ncol = 6)
#' 
#' res <- plot2DNMRspec(peaks = peaks1, x.range = c(0, 5), MHz = 500, ppHz = 1, M = AM)
#' res <- plot2DNMRspec(peaks = peaks1, x.range = c(0, 5), MHz = 500, ppHz = 1,
#' levels = c(0.85, 0.9, 0.95), type = "TOCSY")

plot2DNMRspec <- function (peaks, x.range = c(0, 12), MHz = 300, ppHz = 1,
	type = "COSY", M = NULL, levels = seq(0.5, 1.0, by = 0.1),
	...) {
	
	if ((type == "COSY") & (is.null(M))) stop("You must supply an adjacency matrix indicating the coupling")

	# Create spectra using makeSpec
	# For COSY, add together spectra which make up a multplet (collapse rows),
	# then, using the adjacency matrix, add together any peaks that should be coupled
	# For TOCSY, simply use the full spectrum
	
	# Compute the NMR spectrum
	
	pk <- plotNMRspec(peaks = peaks, x.range = x.range, MHz = MHz, ppHz = ppHz, plot = FALSE)

	# For COSY, collapse the multiplets into single rows
	
	if (type == "COSY") {
		pk2 <- pk[-c(1,2), ] # We don't need the x or y.sum rows
		npk <- matrix(NA_real_, nrow = nrow(peaks), ncol = ncol(pk2))
	
		for (i in 1:nrow(peaks)) {
			end <- sum(peaks$mult[1:i])
			st <- end - peaks$mult[i] + 1
			if (st == end)  npk[i,] <- pk2[st,]	# the rows to gather (all pieces of a multiplet)
			if (st != end)	npk[i,] <- colSums(pk2[c(st:end),])
			}
		
		M1 <- matrix(NA_integer_, nrow = nrow(M), ncol = ncol(npk))	# Prep new data matrix

		for (i in 1:nrow(M)) {
			add <- which(M[i,] == 1L)
			if (length(add) > 0) M1[i,] <- colSums(npk[add,, drop = FALSE]) + npk[i,]
			}
		M1 <- na.omit(M1) # contains just the 'spectra' needed
		}
	
	if (type == "TOCSY") M1 <- pk[2,, drop = FALSE] # Use the entire spectrum for TOCSY

	# Now, add some noisy rows to the actual data.  This gives a much better result
	
	nr <- 10*nrow(M1) # 10x rows seems to work well
	nc <- ncol(M1)
	M2 <- matrix(rnorm(nc*nr, 0, 0.001), nrow = nr, ncol = nc)
	M2[1:nrow(M1),] <- M1 # insert the real data
	
	# And plot
	
	M3 <- cor(M2)
	dd <- seq(x.range[1], x.range[2], length.out = ncol(M3))
	contour(x = dd, y = dd, z = M3, drawlabels = FALSE, levels = levels, ...)
	
	return(M3)
	} # end of plot2DNMRspec
	
	
	

