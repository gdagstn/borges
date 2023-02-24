#' Smooth polygons
#'
#' Uses kernel smoothing to smooth polygons
#'
#' @param poly a data frame containing ordered coordinates with polygon vertices
#' @param smoothness numeric, the extent of kernel smoothing. Higher means
#'     rounder shapes. Default is 3.
#' @param min_points numeric, the minimum number of vertices to smooth.
#'     Default is 8.
#' @param n_dense numeric, the number of points to add to the polygon for more
#'     smoothing. Default is 10.
#'
#' @return A data frame containing ordered coordinates for polygon vertices and
#'    three columns indicating whether they are a hole ("inner") or not ("outer"),
#'    to which cluster they belong, and a sub-clustering to allow \code{ggplot2}
#'    to draw them as holes.
#'
#' @details This is a refactoring of `smoothr::smooth_ksmooth()` to isolate the
#'    necessary code and avoid heavy GDAL-based dependencies. The code has been
#'    simplified assuming `wrap = TRUE` and adding some other bits to handle DF
#'    with other columns as input.
#'
#' @author Matthew Strimas-Mackey, modified by Giuseppe D'Agostino
#'
#' @importFrom stats ksmooth

smoothPolygon <- function(poly, smoothness = 3, min_points = 3, n_dense = 10) {

  if (nrow(poly) < min_points) {

    poly_sm <- poly

  } else {

    poly_coords <- as.matrix(poly[,1:2])

    d_poly = sqrt(rowSums(diff(as.matrix(poly_coords))^2))
    bandwidth = mean(d_poly) * smoothness
    dense_poly = addPoints(poly_coords, steps = n_dense)
    npt = nrow(dense_poly)
    wrapped <- rbind(dense_poly[-npt, ], dense_poly, dense_poly[-1, ])
    d_dense = sqrt(rowSums(diff(as.matrix(wrapped))^2))

    d_x = c(0, cumsum(d_dense))
    poly_sm <- NULL

    for (i in seq_len(ncol(wrapped))) {
      ks <- ksmooth(d_x, wrapped[, i], n.points = length(d_x),
                    kernel = "normal", bandwidth = bandwidth)
      poly_sm <- cbind(poly_sm, ks[["y"]])

      if (i == 1) {
        keep_rows <- (ks$x >= d_x[npt]) & (ks$x <= d_x[(2 * npt - 1)])
      }
    }

    poly_sm <- as.data.frame(poly_sm[keep_rows, ])
    poly_sm[nrow(poly_sm), ] <- poly_sm[1, ]

    #Restore/add original columns
    colnames(poly_sm) <- c("x", "y")
    poly_sm <- as.data.frame(poly_sm)

    for (i in colnames(poly)[3:ncol(poly)]) {
      poly_sm[, i] <- unique(poly[, i])
    }

  }
  return(poly_sm)
}

#' Add points
#'
#' Increases the number of points in a polygon while maintaining the shape
#'
#' @param poly a data frame containing ordered coordinates with polygon vertices
#' @param steps numeric, the number of points that should be added between
#'    each point. Default is 10.
#'
#' @return A data frame containing densified coordinates for polygon vertices and
#'    any other column (assuming original columns contained unique values).
#'
#' @details Internal use only.
#'
#' @author Giuseppe D'Agostino


addPoints <- function(poly, steps = 5) {

  colnames(poly) <- NULL
  polygon_coords = as.matrix(poly[,1:2])
  new_coords = poly[1,]

  for(i in 1:(nrow(poly)-1))
  {
    new_xy = cbind(seq(polygon_coords[i,1],
                       polygon_coords[i+1,1], length.out = steps+1),
                   seq(polygon_coords[i,2],
                       polygon_coords[i+1,2], length.out = steps+1))
    tmp_coords = polygon_coords[i,]
    new_coords = rbind(new_coords, tmp_coords, new_xy)
  }

  new_xy_last = cbind(seq(polygon_coords[nrow(polygon_coords),1],
                          polygon_coords[1,1], length.out = steps+1),
                      seq(polygon_coords[nrow(polygon_coords),2],
                          polygon_coords[1,2], length.out = steps+1))

  new_coords = new_coords[2:nrow(new_coords),]
  new_coords = new_coords[!duplicated(new_coords),]
  new_coords = as.data.frame(rbind(new_coords, new_xy_last))

  rownames(new_coords) = NULL
  colnames(new_coords) <- c("x", "y")

  for (i in colnames(poly)[3:ncol(poly)]) {
    new_coords[, i] <- unique(poly[, i])
  }

  return(new_coords)
}


#' Build a bounding box
#'
#' Creates a bounding box around a subset of square vertices on a grid
#'
#' @param bpg output of \code{boxPointsGrid()}
#' @param stepsize numeric, the size of the step for the grid
#'
#' @author Giuseppe D'Agostino
#'
#' @return A data frame containing x and y coordinates for square vertices
#'    (with z = 1) and the bounding box filled with z = 0. The box is
#'    larger than the set of vertices by \code{stepsize}.

makeBoundingBox <- function(bpg, stepsize = 1) {

  bounds = list("x" = range(bpg$x) + c(-stepsize, stepsize),
                "y" = range(bpg$y) + c(-stepsize, stepsize))

  full = expand.grid(seq(bounds$x[1], bounds$x[2], by = stepsize),
                     seq(bounds$y[1], bounds$y[2], by = stepsize))

  full = as.data.frame(full)
  colnames(full) = c("x", "y")
  full$z = 0

  bound = rbind(bpg, full)
  bound = bound[!duplicated(bound[,1:2]),]

  return(bound)
}


