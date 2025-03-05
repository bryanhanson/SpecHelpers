#'
#' Compute & Possibly Draw an NMR Multiplet with Optional Annotations
#'
#' @param J Numeric. A vector giving the coupling constants.
#' @param pw Numeric.  The value of the peak width at half-height.  Passed to [makeSpec()],
#'        where it is the `gamma` argument.
#' @param plot Logical. Shall the multiplet be drawn?
#' @param plotJtree Logical. Shall the Jtree be drawn? `plot` must also be TRUE in this case, 
#'        and is set automatically if needed.
#' @param showJvalues Logical. Should the J values be added to the plot?  Only relevant if
#'        `plotJtree = TRUE`.
#' @param showJtreeGuides Logical. Shall dotted guides be drawn between the
#'        last leaves of the Jtree and the peak maxima?  Only relevant if
#'        `plotJtree = TRUE`.
#'
#' @return A matrix as produced by [makeSpec()].
#'
#' @author Bryan A. Hanson, DePauw University. \email{hanson@@depauw.edu}
#'
#' @reference Roman A. Valiulin *NMR Multiplet Interpretation*, 2nd Edition, de Gruyter, 2025.
#'
#' @keywords utilities
#' @export
#' @examples
#' # Example 3.1 from Valiulin, a ddt.
#' res <- multiplet(J = c(16.8, 10.1, 6.7, 6.7))
#' # Example 3.2 from Valiulin, a tt.
#' res <- multiplet(J = c(6.1, 6.1, 2.15, 2.15))
#' # Example 3.3 from Valiulin, a dddd.
#' res <- multiplet(J = c(12.7, 12.2, 10.0, 4.9))
#' # Some other nice examples
#' res <- multiplet(J = c(15, 12, 8, 7), pw = 0.25)
#' res <- multiplet(J = c(15, 8, 5, 2))

multiplet <- function(J = c(15, 12, 2), pw = 0.5,
                      plot = TRUE, plotJtree = TRUE,
                      showJvalues = TRUE, showJtreeGuides = TRUE) {
                      	
  if (plotJtree) plot <- TRUE # in case user didn't realize both must be TRUE to get the Jtree

  # helper function
  pm <- function(val) { return(sort(c(-val / 2, val / 2))) } 

  J <- sort(J, decreasing = TRUE)
  no_J <- length(J)

  # Jtree will (initially) contain the positive Hz values at each level of splitting
  # plus the 0.0 line, always present, needed for drawing Jtree if requested
  # Not sorted!  In 'natural' order.
  Jtree <- vector("list", no_J + 1)
  ljt <- length(Jtree)
  Jtree[[1]] <- 0.0
  Jtree[[2]] <- J[1] / 2 # there will always be at least one J value
  if (no_J > 1) {
    for (i in 3:ljt) {
      res <- NA_real_
      for (k in 1:length(Jtree[[i - 1]])) {
        res <- c(res, Jtree[[i - 1]][k] + pm(J[i - 1]))
      }
      Jtree[[i]] <- na.omit(res)
      attributes(Jtree[[i]]) <- NULL
    }
  }

  # the last element of Jtree contains the right half of the final peak positions
  # to be plotted; generally these are positive values but negative values are
  # possible with a lot of J values; for the purpose of drawing the spectrum,
  # order doesn't matter
  peaks <- Jtree[[ljt]]
  peaks <- sort(c(-peaks, peaks))
  pl <- data.frame(x0 = peaks, area = 1.0, gamma = pw)

  dr <- range(peaks) * 2.0
  ans <- makeSpec(pl,
    x.range = dr,
    plot = FALSE, type = "lorentz", dd = 20)

  if (plot) {
    if (!plotJtree) limy <- range(ans["y.sum", ])
    extra <- 2.0 # space above spectrum for drawing the Jtree; could be a fn of no_J
    if (plotJtree) limy <- c(0.0, max(ans["y.sum", ] * extra))

    plot(x = ans["x", ], y = ans["y.sum", ], ylim = limy,
         type = "n", yaxt = "n", bty = "n",
         ylab = "", xlab = "Hz")
    lines(x = ans["x", ], y = ans["y.sum", ])
  }

  if (plotJtree) {
    # get peak maxima locations (the last leaves of the inverted tree)
    np <- nrow(ans) - 2 # no of peaks
    DF <- data.frame(x = rep(NA_real_, np), y = rep(NA_real_, np))
    for (i in 1:np) {
      tmp <- which.max(ans[i + 2, ])
      DF$y[i] <- ans["y.sum", tmp]
      DF$x[i] <- ans["x", tmp]
    }

    # divide the upper part of the Jtree plot into layers to hold the branches
    # every other layer is potentially a different height, controlled by layer_ratio
    # units are the internal ones
    # variable "boundary" refers to the y coord where we switch between vertical and
    # diagonal segments
    layer_ratio <- 0.6 # smaller values mean less height for branches
    total_tree_height <- 0.5 * diff(limy)
    layer_height <- total_tree_height / (no_J + layer_ratio * (no_J + 1))
    boundary <- c(
      total_tree_height,
      total_tree_height +
        cumsum(rep(c(layer_height, layer_ratio * layer_height), times = no_J + 1))
    )
    boundary <- rev(boundary) # reverse as we will draw from the top down

    # reflect the tree (DO NOT sort)
    for (i in 1:ljt) {
      Jtree[[i]] <- c(-rev(Jtree[[i]]), Jtree[[i]])
    }

    # draw the tree, working from the top
    for (i in 1:ljt) {
      for (j in 1:length(Jtree[[i]])) { # j moves horizontally
        for (k in 1:(length(boundary) - 1)) {
          if ((i * 2) == k) { # draw vertical segments; x-coord does not change layer-to-layer
            segments(Jtree[[i]][j], boundary[k], Jtree[[i]][j], boundary[k + 1], col = "red")
          }
          if ((i * 2) == (k - 1)) { # draw diagonal segments
            segments(Jtree[[i]][j], boundary[k], Jtree[[i + 1]][j * 2 - 1], boundary[k + 1], col = "red")
            segments(Jtree[[i]][j], boundary[k], Jtree[[i + 1]][j * 2], boundary[k + 1], col = "red")
          }
        }
      }
    }

    if (showJtreeGuides) {
      # draw dotted lines from the last leaves the tree to the peak maxima
      # lines must end on the overal spectral envelope, not a peak that might
      # be buried within it
      lastk <- boundary[length(boundary)]
      Jtree[[ljt]] <- sort(Jtree[[ljt]]) # for this purpose we must sort
      for (i in 1:length(Jtree[[ljt]])) {
        segments(Jtree[[ljt]][i], lastk, Jtree[[ljt]][i], DF$y[i], col = "pink", lty = 2)
      }
    }

    if (showJvalues) {
      # annotate with the J values
      labs <- paste("J =", J, sep = " ")
      lab_x_pos <- 0.8 * max(ans["x", ])
      use <- 4:length(boundary)
      use <- use[use %% 2 == 0]
      bump <- 0.5 * diff(boundary)
      lab_y_pos <- boundary[use] + bump[use]
      text(x = lab_x_pos, y = lab_y_pos, labels = labs)
    }
    
  } # end of plotJtree

  invisible(ans)
}
