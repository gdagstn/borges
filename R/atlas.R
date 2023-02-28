#' Prepare a SingleCellExperiment Atlas
#'
#' Prepares a 2D atlas for plotting
#'
#' @param x a \code{SingleCellExperiment} class object, OR a \code{matrix} with 2 columns, OR a \code{data.frame} with 2 columns. 
#' @param res numeric, the resolution of the boundary estimation via \code{oveRlay::makeOverlay()}. Default is 300
#' @param labels character, the name of the column in \code{colData(sce)} with labels
#'     OR a character vector of labels of length \code{nrow(x)} if is a \code{matrix} or a \code{data.frame}
#' @param dimred character, the name of the reduced dimension in \code{reducedDim(sce)} 
#'     to be used as basis. Only used if \code{x} is a \code{SingleCellExperiment}.
#' @param as_map logical, should the coordinates be changed to accommodate a 
#'     geographical projection? Default is FALSE
#'
#' @return a list containing the following slots:
#'    `atlas` - the atlas coordinates and labels;
#'    `dr` - a data.frame containing coordinates in 2D embedding of choice
#'    `coords` - a data.frame containing coordinates for label positioning
#'    `res` - the resolution used
#'
#' @author Giuseppe D'Agostino
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment reducedDimNames reducedDim
#' @importFrom methods is
#' @importFrom scales rescale
#' @importFrom stats median
#'
#' @export


prepAtlas <- function (x, res = 300, labels, dimred = NULL, as_map = FALSE) {
  # Sanity checks
  
  if (!is(x, "SingleCellExperiment") & !is(x, "data.frame") & 
      !(is(x, "matrix"))) 
    stop("You should provide either a `SingleCellExperiment` object, a data.frame, or a matrix.")
  if (is(x, "SingleCellExperiment") & is.null(dimred)) 
    stop("You should provide a dimensional reduction name (`dimred`) if you are using a `SingleCellExperiment` object as input.")
  if (is(x, "SingleCellExperiment") & is.null(labels)) 
    stop("You should provide a `colData` column name for labels (`labels`) if you are using a `SingleCellExperiment` object as input.")
  if (!is(res, "numeric")) 
    stop("`res` must be numeric")
  if (res < 10) 
    stop("`res` must be at least 10")
  if ((is(x, "matrix") | is(x, "data.frame")) & ncol(x) < 2) 
    stop("`x` must contain at least two columns (x and y coordinates)")
  if (is(x, "SingleCellExperiment")) {
    if (!labels %in% colnames(colData(x))) 
      stop("`labels` not found in the colData slot of the SingleCellExperiment object")
    if (!dimred %in% reducedDimNames(x)) 
      stop("`dimred` not found in the reducedDim slot of the SingleCellExperiment object")
  }
  if (!is(x, "SingleCellExperiment") & (!is(labels, "character") & 
                                        !is(labels, "factor"))) 
    stop("`labels` must be a character vector or a factor")
  if (!is(x, "SingleCellExperiment") & (is(labels, "character") | 
                                        is(labels, "factor")) & length(labels) != nrow(x)) 
    stop("`labels` must have the same length as the number of points")
  
  # End checks
  
  if (is(x, "SingleCellExperiment")) {
    dr = reducedDim(x, dimred)
  }
  else {
    dr = x
  }
  if (any(is.na(dr[, 1]))) {
    if (is(x, "SingleCellExperiment")) {
      x = x[, which(!is.na(dr[, 1]) & !is.na(dr[, 2]))]
      dr = reducedDim(x, dimred)
    }
    else {
      dr = dr[!is.na(dr[, 1]) & !is.na(dr[, 2]), ]
    }
  }
  colnames(dr) = c("x", "y")
  if (is(x, "SingleCellExperiment")) {
    labels = colData(x)[, labels]
  }
  else {
    labels = labels
  }
  if (as_map) {
    dr[, 1] = rescale(dr[, 1], to = c(-120, +120))
    dr[, 2] = rescale(dr[, 2], to = c(-70, +70))
  }
  atlas = makeBoundaries(data = dr, res = res, labels = labels)
  coords = cbind(tapply(X = dr[, 1], INDEX = labels, FUN = median), 
                 tapply(X = dr[, 2], INDEX = labels, FUN = median))
  colnames(coords) = c("x", "y")
  coords = as.data.frame(coords)
  coords$n = as.numeric(table(labels))
  ret = list(atlas = atlas, dr = dr, coords = coords, res = res)
  return(ret)
}

#' Plot a SingleCellExperiment Atlas
#'
#' Plots a 2D atlas from a SingleCellExperiment
#'
#' @param atlas_ret a list as returned by \code{prepAtlas()}
#' @param plot_cells logical, should cells/points be plotted? Default is TRUE
#' @param add_contours logical, should density contours be drawn? Default is TRUE
#' @param show_labels logical, should labels be plotted? Default is TRUE
#' @param shade_borders logical, should borders be shaded? Default is TRUE
#' @param shade_offset numeric, the length and direction of offset segments. Default is -1 (inwards).
#' @param shade_skip numeric, the number of points to be skipped between shade segments. Default is 10.
#' @param as_map logical, should the coordinates be changed to accommodate a geographical projection? Default is FALSE
#' @param map_proj character, the map projection. Accepts any one character argument to \code{ggplot2::coord_map()}
#' @param map_theme character, one of "classic" (default), "renaissance", "medieval", "modern"
#' @param pal character, a vector of colors. Default is NULL and decided by the theme
#' @param capitalize_labels logical, should labels be capitalized? Default is FALSE
#'
#' @return a \code{ggplot2} plot object with an atlas-like representation of a 2D point cloud
#'
#' @author Giuseppe D'Agostino
#'
#' @importFrom scales rescale
#' @importFrom oveRlay makeOverlay
#' @importFrom MASS kde2d
#' @importFrom colorspace darken
#' @importFrom grDevices rgb
#' @importFrom ggplot2 ggplot aes .data geom_polygon geom_point geom_density_2d element_rect element_line unit
#' @importFrom ggplot2 scale_x_continuous ylim theme_void theme coord_map coord_fixed scale_fill_manual scale_color_manual
#' @importFrom ggrepel geom_label_repel
#' @importFrom methods is
#'
#' @export

plotAtlas <- function(atlas_ret, plot_cells = TRUE, add_contours = TRUE,
                      show_labels = TRUE, shade_borders = TRUE, shade_offset = -1, 
                      shade_skip = 10, as_map = FALSE, map_proj = "lagrange",
                      map_theme = "classic", pal = NULL, capitalize_labels = FALSE) {

  atlas = atlas_ret$atlas
  dr = atlas_ret$dr
  coords = atlas_ret$coords
  res = atlas_ret$res
  
  if(shade_borders & (!is(shade_offset, "numeric") | !is(shade_skip, "numeric"))) stop("`shade_skip` and `shade_offset` must both be numeric")
  
  if(is.null(pal) & is.null(map_theme)) stop("Must provide at least one of `pal` or `map_theme`")

  if(!is.null(map_theme)) {
    maptheme = mapTheme(map_theme)
    if(!is.null(pal)) maptheme$pal = pal
  }
  if(is.null(map_theme)) {
    maptheme = mapTheme("classic")
    if(!is.null(pal)) maptheme$pal = pal
  }

  if(!as_map) {
    dr[,1] = rescale(dr[,1], to = range(dr[,2]))
    atlas[,1] = rescale(atlas[,1], to = range(atlas[,2]))
    coords[,1] = rescale(coords[,1], to = range(coords[,2]))
  }

  ov = makeOverlay(atlas[,1:2], res = res, offset_prop = 0)
  ov1 = makeOverlay(atlas[,1:2], res = res, offset_prop = 0.004)
  ov2 = makeOverlay(atlas[,1:2], res = res, offset_prop = 0.008)
  ov3 = makeOverlay(atlas[,1:2], res = res, offset_prop = 0.012)
  
  kde = kde2d(dr[,1], dr[,2], n = 50)
  b2 = max(kde$z)
  b1 = 0 + diff(range(kde$z))/3

  if(capitalize_labels) {
    coords$label = toupper(rownames(coords))
  } else {
    coords$label = rownames(coords)
  }

  p = ggplot(data = atlas, aes(x = .data[["x"]], y = .data[["y"]],
                               group = .data[["cluster_label"]], fill = .data[["label"]],
                               col = .data[["label"]])) +
    geom_polygon(linewidth = 0.5) +
    scale_fill_manual(values = maptheme$pal) +
    scale_color_manual(values = darken(maptheme$pal, amount = 0.5)) +
    geom_polygon(data = ov[ov$hole == "outer",],
                 inherit.aes = FALSE,
                 mapping = aes(x = .data[["x"]], y = .data[["y"]], group = .data[["id_hole"]]),
                 fill = NA,
                 linewidth = 0.3,
                 linetype = "solid",
                 col = rgb(0,0,0,0.5)) +
    geom_polygon(data = ov1[ov1$hole == "outer",],
                 inherit.aes = FALSE,
                 mapping = aes(x = .data[["x"]], y = .data[["y"]], group = .data[["id_hole"]]),
                 fill = NA,
                 linewidth = 0.3,
                 linetype = "solid",
                 color = rgb(0,0,0,0.4)) +
    geom_polygon(data = ov2[ov2$hole == "outer",],
                 inherit.aes = FALSE,
                 mapping = aes(x = .data[["x"]], y = .data[["y"]], group = .data[["id_hole"]]),
                 fill = NA,
                 linewidth = 0.3,
                 linetype = "22",
                 color = rgb(0,0,0,0.4)) +
    geom_polygon(data = ov3[ov3$hole == "outer",],
                 inherit.aes = FALSE,
                 mapping = aes(x = .data[["x"]], y = .data[["y"]], group = .data[["id_hole"]]),
                 fill = NA,
                 linewidth = .3,
                 linetype = "13",
                 color = rgb(0,0,0,0.4))
  
  names(p$layers) = c("waiver", "overlay1", "overlay2", "overlay3", "overlay4")
  
  if(shade_borders) {
    shade = addShading(p$data, offset = shade_offset, skip = shade_skip)
    p = p + geom_segment(data = shade, 
                         mapping = aes(x = x0, xend = x1,
                                       y = y0, yend = y1),
                         inherit.aes = FALSE,
                         color = "black",
                         #alpha = 0.2,
                         linewidth = 0.1)
    
    names(p$layers)[length(p$layers)] = "shade"
  }
  
  if(plot_cells) {
    point_alpha = (rescale(c(0, b1, b2), to = c(0, 1))/20)[2]
        p = p + geom_point(data = as.data.frame(dr),
                       inherit.aes = FALSE,
                       mapping = aes(x = .data[["x"]], y = .data[["y"]]),
                       size = 0.01,
                       color = rgb(0,0,0, 0.02))
        
    names(p$layers)[length(p$layers)] = "points"
    
  }
  if(add_contours) {
    p = p + geom_density_2d(data = as.data.frame(dr),
                            inherit.aes = FALSE,
                            mapping = aes(x = .data[["x"]], y = .data[["y"]]),
                                          col = rgb(0,0,0,0.25),
                                          linewidth = 0.2,
                                          breaks = seq(b1, b2, length.out = 20))
    
    names(p$layers)[length(p$layers)] = "kde_contours"
  }
  if(show_labels) {

    p = p + geom_label_repel(data = coords,
                             mapping = aes(x = .data[["x"]], y = .data[["y"]],
                                           label = .data[["label"]]),
                             inherit.aes = FALSE,
                             size = maptheme$labelsize,
                             label.r = 0.01,
                             label.padding = unit(0.1, "lines"),
                             fill = maptheme$labelfill,
                             alpha = 0.8,
                             label.size = maptheme$border,
                             family = maptheme$labelfont)
    names(p$layers)[length(p$layers)] = "labels"

  }
  p = p +
    theme_void() +
    theme(panel.background = element_rect(fill = maptheme$background, color = NA),
          legend.position = "none",
          panel.grid.major = element_line(color = maptheme$gridlines, linewidth = 0.6),
          panel.grid.minor = element_line(color = maptheme$gridlines, linewidth = 0.25))

  if(as_map) {
    p = p + scale_x_continuous(breaks = seq(-120, 120, by = 30),
                               minor_breaks = seq(-120, 120, by = 15)) +
      ylim(-90, 110) +
      coord_map(map_proj)
  } else {
    p = p + coord_fixed()
  }
  return(p)
}

#' @importFrom scales rescale
#' @noRd
boxPointsGrid = function(data, res, stepsize = NULL, labels = NULL, min_pts = NULL) {

  if (!is.null(stepsize) & is.null(res))
    res = 1/stepsize

  mat_r = data
  mat_r[, 1] = rescale(mat_r[, 1], to = c(0, res))
  mat_r[, 2] = rescale(mat_r[, 2], to = c(0, res))
  breaks = seq(-1, res + 2)
  xcut = factor(findInterval(mat_r[, 1], breaks), levels = breaks)
  ycut = factor(findInterval(mat_r[, 2], breaks), levels = breaks)

  if(!is.null(labels)) {
    udf = as.data.frame(data)
    colnames(udf) = c("x", "y")
    udf$label = as.character(labels)
    udf$xcut = xcut
    udf$ycut = ycut
    udf$label[is.na(udf$label)] = "NA"
  }

  boxes = table(ycut, xcut)
  boxes = as.matrix(boxes[nrow(boxes):1, ])
  bgrid = expand.grid(breaks, breaks)
  colnames(bgrid) = c("x", "y")
  bgrid$z = as.numeric(t(apply(boxes, 2, rev)))
  bgrid_ne = bgrid[bgrid$z > 0,]

  if(!is.null(labels)) {
    bgrid_xy = paste0(bgrid_ne$x, "_", bgrid_ne$y)
    names(bgrid_xy) = rownames(bgrid_ne)
    bgn = names(bgrid_xy)
    names(bgn) = bgrid_xy
    udf$xy = paste0(udf$xcut, "_", udf$ycut)
    into = udf[udf$xy %in% bgrid_xy, ]
    into$rown = bgn[into$xy]

    intol = lapply(split(into, into$rown), function(x) {
      if(length(x) == 1) return(x) else return(names(which.max(table(x$label))))
    })

    intol = unlist(intol)
    bgrid_ne$maj = intol[rownames(bgrid_ne)]
    bgrid$maj = NA
    bgrid[rownames(bgrid_ne), "maj"] = bgrid_ne$maj
  }

  if (!is.null(min_pts)) bgrid$z[bgrid$z < min_pts] = 0
  bgrid$Ax = bgrid$x
  bgrid$Ay = bgrid$y
  bgrid$Bx = bgrid$x - 1
  bgrid$By = bgrid$y
  bgrid$Cx = bgrid$x - 1
  bgrid$Cy = bgrid$y - 1
  bgrid$Dx = bgrid$x
  bgrid$Dy = bgrid$y - 1

  boxes_empty = as.data.frame(mapply(c, bgrid[bgrid$z == 0, c("Ax", "Ay")],
                                     bgrid[bgrid$z == 0, c("Bx", "By")],
                                     bgrid[bgrid$z == 0, c("Cx", "Cy")],
                                     bgrid[bgrid$z == 0, c("Dx", "Dy")]))

  colnames(boxes_empty) = c("x", "y")
  boxes_empty$z = 0
  if(!is.null(labels)) boxes_empty$maj = 0

  boxes_nempty = as.data.frame(mapply(c, bgrid[bgrid$z > 0, c("Ax", "Ay")],
                                      bgrid[bgrid$z > 0, c("Bx", "By")],
                                      bgrid[bgrid$z > 0, c("Cx", "Cy")],
                                      bgrid[bgrid$z > 0, c("Dx", "Dy")]))

  colnames(boxes_nempty) = c("x", "y")
  boxes_nempty$z = 1
  if(!is.null(labels)) boxes_nempty$maj = bgrid_ne$maj
  boxes_final = rbind(boxes_nempty, boxes_empty)
  boxes_final = boxes_final[!duplicated(boxes_final[, 1:2]), ]

  return(boxes_final)
}

#' @importFrom isoband isobands
#' @importFrom reshape2 acast
#' @importFrom scales rescale
#' @noRd
makeBoundaries = function (data, res, labels, stepsize = NULL, min_pts = NULL, minsize = 4, smooth = TRUE) {

  bpg = boxPointsGrid(data, res, labels, stepsize = NULL, min_pts)
  bpg = bpg[bpg$maj != 0,]
  bpgl = split(bpg, bpg$maj)

  bpglc = lapply(names(bpgl), function(n) {
    x = bpgl[[n]]
    bound = makeBoundingBox(x[,1:3])
    m = t(acast(bound, formula = x ~ y, value.var = "z"))
    ib = isobands(x = seq_len(ncol(m)), y = seq_len(nrow(m)),
                  z = m, levels_low = 0, levels_high = 1)
    ib_df = data.frame(x = ib[[1]]$x, y = ib[[1]]$y, cluster = ib[[1]]$id)
    ib_df <- ib_df[which(!ib_df$x %in% range(ib_df$x) & !ib_df$y %in%
                           range(ib_df$y)), ]
    ib_df$x <- rescale(ib_df$x, to = range(bound[bound$z == 1,
                                                 1]))
    ib_df$y <- rescale(ib_df$y, to = range(bound[bound$z == 1,
                                                 2]))
    ib_df <- ib_df[!duplicated(ib_df[, 1:2]), ]
    grid_range = list(x = range(x[, 1]), y = range(x[,2]))

    ib_df[, 1] = rescale(ib_df[, 1], to = grid_range$x)
    ib_df[, 2] = rescale(ib_df[, 2], to = grid_range$y)

    if (!is.null(minsize)) {
      keep = names(table(ib_df$cluster))[table(ib_df$cluster) >= minsize]
      if(length(keep) == 0) return(NULL)
      ib_df = ib_df[ib_df$cluster %in% keep,]
    }

    if(smooth) {
      ibl = split(ib_df, ib_df$cluster)
      ibs = lapply(ibl, function(x) smoothPolygon(x, min_points = 3))
      ib_df = as.data.frame(do.call(rbind, ibs))
    }

    lab = gsub(" ", "_", n)
    ib_df$label = lab
    ib_df$cluster_label = paste0(lab, "_", ib_df$cluster)

    return(ib_df)
  })

  bpglc = bpglc[!unlist(lapply(bpglc, is.null))]
  ret = do.call(rbind, bpglc)
  ret = ret[!is.na(ret[,1]),]
  ret = rescale_res(data = data, ret = ret, res = res)

  return(ret)
}

#' @noRd
within_interval <- function(x, int, closed = FALSE){
  if(closed) {
    keep = which(x < max(int) & x > min(int))
  } else {
    keep = which(x <= max(int) & x >= min(int))
  }
  return(keep)
}

#' @importFrom scales rescale
#' @noRd
rescale_res <- function(data, ret, res) {
  mat_r = data
  mat_r[, 1] = rescale(mat_r[, 1], to = c(0, res))
  mat_r[, 2] = rescale(mat_r[, 2], to = c(0, res))

  keep_x = within_interval(as.numeric(mat_r[,1]), range(ret[,1]))
  keep_y = within_interval(as.numeric(mat_r[,2]), range(ret[,2]))
  keep = intersect(keep_x, keep_y)
  ret[,1] = rescale(ret[,1], to = range(data[keep,1]))
  ret[,2] = rescale(ret[,2], to = range(data[keep,2]))

  return(ret)
}

#' @noRd
mapTheme <- function(theme) {

  if(theme == "classic"){
    pal = c("#F3F6F2", "#F2F3DE", "#E2F3DF", "#F7E18F", "#DEE6A8", "#E1DDA5",
            "#E5D49C", "#F6CE81", "#D6CB9E", "#EEBE82", "#B3CA8D", "#AABBB8",
            "#BCB482", "#ACB77F", "#ADB2AE", "#CFAB74", "#C8A882", "#D09F7C",
            "#A3AB95", "#A4A7B9", "#D09B81", "#94A79C", "#9D7A5A", "#8E7967",
            "#857A52", "#6F7579", "#767068")
    background = "#c4e2e2"
    gridlines = "#a5cccc"
    labelfont = "serif"
    labelfill = "blanchedalmond"
    labelsize = 3
    border = 0.1
    labelstroke = 0.5
  }

  if(theme == "modern"){
    pal = c("#F6EF96", "#F7E18F", "#C08968", "#C5B373", "#ADA367",
            "#F4C1AC", "#9A8754", "#C08D7F", "#EEB8AE", "#A37263", "#91CAE1",
            "#B7CE9D", "#CDE4B6", "#9BB186", "#978F67", "#AA8178",
            "#9D88AD", "#8D8BAC", "#F5EECC", "#BFA1BF", "#94B4C4", "#9BC0A8",
            "#82745A", "#CEE6CA", "#818668", "#8B7791", "#9A977D", "#B0AAC3",
            "#A396AF", "#C3B0C8", "#71808F", "#725F71", "#6F7282", "#91A3A8",
            "#E7F4F4", "#C6CFCF", "#777E79")
    background = "#81BFDF"
    gridlines = "#B5DDF4"
    labelfont = "sans"
    labelfill = "white"
    labelsize = 2.5
    border = 0.25
    labelstroke = 0

  }

  if(theme == "renaissance"){
    pal = c("#945934", "#C08968", "#D09F7C", "#B89E66", "#E1BD8E", "#BD9A6E",
            "#ADA367", "#9A8754", "#EFD5AB", "#BCB482", "#775033", "#9D7A5A",
            "#7D7A45", "#F2E5B8", "#8C724A", "#EFCBAF", "#959969", "#785D3A",
            "#F0DBBC", "#A18E70", "#B4A384", "#4E4517", "#492611", "#F5EECC",
            "#D0C3A4", "#5A5130", "#6D6F4E", "#CEB8A4", "#82745A", "#5F4C34",
            "#BDA394", "#90726A", "#5D5541", "#372A12", "#5C4A3F", "#868370",
            "#B9B2A3", "#4D3936", "#8A7E7B", "#514C46")
    background = "#edecda"
    gridlines = "#d4d3c5"
    labelfont = "serif"
    labelfill = "blanchedalmond"
    labelsize = 4
    border = 0.1
    labelstroke = 0

  }
  if(theme == "medieval"){

    pal = c("#F2E7D3", "#F2DBCF", "#D6CB9E", "#CFC8C6", "#D0C3A4", "#B4A992",
            "#BFA087", "#B4A384", "#C09C99", "#88996C", "#A18E70",
            "#9A8C92", "#A18961", "#AA8178", "#818668", "#9A7585", "#777E79",
            "#6C7D67", "#956246", "#6D6F5B", "#5C6F86", "#6D6244", "#B5301E",
            "#83543F", "#AA352D", "#535E4A", "#7B4E56", "#5D5541", "#515254",
            "#5F4C34", "#3D4958", "#66392E", "#5C3B44", "#2F4460", "#3E443C",
            "#434231", "#354335", "#362720", "#242537")
    background = "#A3AB95"
    gridlines = "#00000000"
    labelfont = "serif"
    labelsize = 4
    labelfill = "blanchedalmond"
    border = 0.25
    labelstroke = 0
  }

  return(list(pal = pal,
              background = background,
              gridlines = gridlines,
              labelfont = labelfont,
              labelfill = labelfill,
              labelsize = labelsize,
              border = border))
}


#' @details taken from https://stackoverflow.com/a/69936814/8391346
#' @noRd

polygonOffset <- function(x, y, d) {
  angle <- atan2(diff(y), diff(x)) + pi/2
  angle <- c(angle[1], angle)
  data.frame(x = d * cos(angle) + x, y = d * sin(angle) + y)
}

#' @importFrom polyclip pointinpolygon 
#' @noRd

addShading <- function(p, offset = -1, skip = 10) {
  p = p[seq(1, nrow(p), by = skip),]
  p2 = polygonOffset(p$x, p$y, offset)
  df = data.frame(x0 = p$x, x1 = p2$x, y0 = p$y, y1 = p2$y)
  ins = lapply(split(p, p$cluster_label),
               function(x) pointinpolygon(P = list(x = df$x1, y = df$y1),
                                                     A = list(x = x$x, y = x$y)))
  insdf = do.call(cbind, ins)
  inp = which(rowSums(insdf) != 0)
  df = df[inp,]
  return(df)
}