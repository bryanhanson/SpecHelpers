#'
#' Compute & Possibly Draw an NMR Multiplet with Optional Annotations
#'
#' Serves as a teaching and self-study tool to understand complex NMR multiplets.
#' Inspired by the Valiulin book (see the reference).  One can draw a multiplet,
#' and optionally draw the splitting tree along with annotations of the J values
#' and guides connecting the tree to the peak maxima.
#'
#' @param J Numeric. A vector giving the coupling constants.
#' @param I Numeric. Nuclear spin quantum number.  Half or whole integer.  Currently allowed
#'        values are 1/2, 1, 3/2, 5/2, 3 (note there are no stable isotopes with I = 2).
#' @param pw Numeric.  Half the peak width at half-maximum (HWHM).  Passed to [makeSpec()],
#'        and then [lorentzCurve()] where it is the `gamma` argument.
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
#' @references Roman A. Valiulin *NMR Multiplet Interpretation*, 2nd Edition, de Gruyter, 2025.
#'
#' @keywords utilities
#' @export
#' @examples
#' # Examples of I = 1/2
#' # Example 3.1 from Valiulin, a ddt.
#' res <- multiplet(J = c(16.8, 10.1, 6.7, 6.7))
#' # Example 3.2 from Valiulin, a tt.
#' res <- multiplet(J = c(6.1, 6.1, 2.15, 2.15))
#' # Example 3.3 from Valiulin, a dddd.
#' res <- multiplet(J = c(12.7, 12.2, 10.0, 4.9))
#' # Some other nice examples
#' res <- multiplet(J = c(15, 12, 8, 7), pw = 0.25)
#' res <- multiplet(J = c(15, 8, 5, 2))
#'
#' # Examples of I = 1
#' res <- multiplet(J = 32, I = 1) # CDCl3 13C -> 2H coupling
#' res <- multiplet(J = c(20, 18), I = 1)
#'
#' # Examples of I = 3/2
#' res <- multiplet(J = 1.13, I = 3 / 2, pw = 0.1) # NaBF4 19F -> 11B coupling
#' res <- multiplet(J = c(10, 7), I = 3 / 2)
multiplet <- function(J = c(15, 12, 2), I = 1 / 2, pw = 0.5,
                      plot = TRUE, plotJtree = TRUE,
                      showJvalues = TRUE, showJtreeGuides = TRUE) {
  # helper function
  compute_peaks <- function(Jtree, idx, I, J) {
    # Compute the peak positions for the current level idx by inspecting the idx - 1
    # level, except in the case of idx = 1 which is trivially initialized
    if (idx == 1L) {
      return(0.0)
    }
    res <- NA_real_
    jdx <- idx - 1 # J is one element shorter than Jtree
    for (i in 1:length(Jtree[[jdx]])) {
      if (I == 1 / 2) offset <- c(-J[jdx] / 2, J[jdx] / 2)
      if (I == 1) offset <- c(-J[jdx], 0.0, J[jdx])
      if (I == 3 / 2) offset <- c(-J[jdx] * 3 / 2, -J[jdx] / 2, J[jdx] / 2, J[jdx] * 3 / 2)
      if (I == 5 / 2) offset <- c(-J[jdx] * 5 / 2, -J[jdx] * 3 / 2, -J[jdx] / 2, J[jdx] / 2, J[jdx] * 3 / 2, J[jdx] * 5 / 2)
      if (I == 3) offset <- c(-J[jdx] * 3, -J[jdx] * 2, -J[jdx], 0.0, J[jdx], J[jdx] * 2, J[jdx] * 3)
      res <- c(res, Jtree[[jdx]][i] + offset)
    }
    res <- res[-1]
  }

  OK_I <- c(1 / 2, 1, 3 / 2, 5 / 2, 3)
  if (!I %in% OK_I) stop("Don't know how to handle that value of I")
  if (plotJtree) plot <- TRUE # in case user didn't realize both must be TRUE to get the Jtree
  J <- sort(J, decreasing = TRUE)
  no_J <- length(J)
  no_Jtree <- no_J + 1

  # Jtree will store the Hz values at each level of splitting relative to 0.0
  # as well as the 0.0 line, needed for drawing Jtree if requested

  Jtree <- vector("list", no_Jtree)
  names(Jtree) <- paste("Lvl", 1:no_Jtree, sep = "_")

  for (i in 1:length(Jtree)) {
    Jtree[[i]] <- compute_peaks(Jtree, i, I, J)
  }

  peaks <- Jtree[[no_Jtree]]
  pl <- data.frame(x0 = peaks, area = 1.0, gamma = pw)

  dr <- range(peaks) * 1.5
  ans <- makeSpec(pl,
    x.range = dr,
    plot = FALSE, type = "lorentz", dd = 40
  )

  if (plot) {
    if (!plotJtree) limy <- range(ans["y.sum", ])
    extra <- 2.0 # space above spectrum for drawing the Jtree; could be a fn of no_J
    if (plotJtree) limy <- c(0.0, max(ans["y.sum", ] * extra))

    plot(
      x = ans["x", ], y = ans["y.sum", ], ylim = limy,
      type = "n", yaxt = "n", bty = "n",
      ylab = "", xlab = "Hz"
    )
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

    # divide the upper part of the Jtree plot into layers to hold the branches and guides
    # layer_ratio specifies the space allocated for vertical segments vs diagonal segments
    # units are the internal ones
    # variable "boundary" refers to the y coord where we switch between vertical and diagonal segments
    layer_ratio <- 1.0 # smaller values mean less space allocated for diagonal dotted lines
    total_tree_height <- 0.5 * diff(limy)
    layer_height <- total_tree_height / (no_J + (no_J + 1))
    boundary <- c(
      total_tree_height,
      total_tree_height +
        cumsum(rep(c(layer_height, layer_ratio * layer_height), times = no_J + 1))
    )
    boundary <- rev(boundary) # reverse as we will draw from the top down

    # draw the tree, working from the top
    for (i in 1:no_Jtree) {
      for (j in 1:length(Jtree[[i]])) { # j moves horizontally
        for (k in 1:(length(boundary) - 1)) {
          if ((i * 2) == k) { # draw VERTICAL segments; x-coord does not change layer-to-layer
            segments(Jtree[[i]][j], boundary[k], Jtree[[i]][j], boundary[k + 1], col = "red")
          }
          if ((i * 2) == (k - 1)) { # draw DIAGONAL segments
            if (I == 1 / 2) { # in this case 2 segments have to be drawn
              segments(Jtree[[i]][j], boundary[k], Jtree[[i + 1]][j * 2 - 1], boundary[k + 1], col = "pink", lty = 3)
              segments(Jtree[[i]][j], boundary[k], Jtree[[i + 1]][j * 2], boundary[k + 1], col = "pink", lty = 3)
            }
            if (I == 1) { # in this case 3 segments have to be drawn
              segments(Jtree[[i]][j], boundary[k], Jtree[[i + 1]][j * 3 - 1], boundary[k + 1], col = "pink", lty = 3)
              segments(Jtree[[i]][j], boundary[k], Jtree[[i + 1]][j * 3 - 2], boundary[k + 1], col = "pink", lty = 3)
              segments(Jtree[[i]][j], boundary[k], Jtree[[i + 1]][j * 3], boundary[k + 1], col = "pink", lty = 3)
            }
            if (I == 3 / 2) { # in this case 4 segments have to be drawn
              segments(Jtree[[i]][j], boundary[k], Jtree[[i + 1]][j * 4 - 1], boundary[k + 1], col = "pink", lty = 3)
              segments(Jtree[[i]][j], boundary[k], Jtree[[i + 1]][j * 4 - 2], boundary[k + 1], col = "pink", lty = 3)
              segments(Jtree[[i]][j], boundary[k], Jtree[[i + 1]][j * 4 - 3], boundary[k + 1], col = "pink", lty = 3)
              segments(Jtree[[i]][j], boundary[k], Jtree[[i + 1]][j * 4], boundary[k + 1], col = "pink", lty = 3)
            }
            if (I == 5 / 2) { # in this case 6 segments have to be drawn
              segments(Jtree[[i]][j], boundary[k], Jtree[[i + 1]][j * 6 - 1], boundary[k + 1], col = "pink", lty = 3)
              segments(Jtree[[i]][j], boundary[k], Jtree[[i + 1]][j * 6 - 2], boundary[k + 1], col = "pink", lty = 3)
              segments(Jtree[[i]][j], boundary[k], Jtree[[i + 1]][j * 6 - 3], boundary[k + 1], col = "pink", lty = 3)
              segments(Jtree[[i]][j], boundary[k], Jtree[[i + 1]][j * 6 - 4], boundary[k + 1], col = "pink", lty = 3)
              segments(Jtree[[i]][j], boundary[k], Jtree[[i + 1]][j * 6 - 5], boundary[k + 1], col = "pink", lty = 3)
              segments(Jtree[[i]][j], boundary[k], Jtree[[i + 1]][j * 6], boundary[k + 1], col = "pink", lty = 3)
            }
            if (I == 3) { # in this case 7 segments have to be drawn
              segments(Jtree[[i]][j], boundary[k], Jtree[[i + 1]][j * 7 - 1], boundary[k + 1], col = "pink", lty = 3)
              segments(Jtree[[i]][j], boundary[k], Jtree[[i + 1]][j * 7 - 2], boundary[k + 1], col = "pink", lty = 3)
              segments(Jtree[[i]][j], boundary[k], Jtree[[i + 1]][j * 7 - 3], boundary[k + 1], col = "pink", lty = 3)
              segments(Jtree[[i]][j], boundary[k], Jtree[[i + 1]][j * 7 - 4], boundary[k + 1], col = "pink", lty = 3)
              segments(Jtree[[i]][j], boundary[k], Jtree[[i + 1]][j * 7 - 5], boundary[k + 1], col = "pink", lty = 3)
              segments(Jtree[[i]][j], boundary[k], Jtree[[i + 1]][j * 7 - 6], boundary[k + 1], col = "pink", lty = 3)
              segments(Jtree[[i]][j], boundary[k], Jtree[[i + 1]][j * 7], boundary[k + 1], col = "pink", lty = 3)
            }
          }
        }
      }
    }

    if (showJtreeGuides) {
      # draw dotted lines from the last leaves the tree to the peak maxima
      # lines must end on the overall spectral envelope, not a peak that might be buried within it
      lastk <- boundary[length(boundary)]
      for (i in 1:length(Jtree[[no_Jtree]])) {
        segments(Jtree[[no_Jtree]][i], lastk, Jtree[[no_Jtree]][i], DF$y[i], col = "pink", lty = 3)
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
