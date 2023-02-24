<<<<<<< HEAD
<img src="https://user-images.githubusercontent.com/21171362/220892837-bbabcf91-d6a5-45a9-835b-fa0a9b461f93.png" align="right" width="180"/>
=======
<img src="https://user-images.githubusercontent.com/21171362/220892837-bbabcf91-d6a5-45a9-835b-fa0a9b461f93.png" align="right" alt="" width="180" />
>>>>>>> 96764cbc1c859af827ebae0c0e4ed4f0669597fb

# borges

> Of *Exactitude in Science*
>
> ...In that Empire, the craft of Cartography attained such Perfection that the Map of a Single province covered the space of an entire City, and the Map of the Empire itself an entire Province. In the course of Time, these Extensive maps were found somehow wanting, and so the College of Cartographers evolved a Map of the Empire that was of the same Scale as the Empire and that coincided with it point for point. Less attentive to the Study of Cartography, succeeding Generations came to judge a map of such Magnitude cumbersome, and, not without Irreverence, they abandoned it to the Rigours of sun and Rain. In the western Deserts, tattered Fragments of the Map are still to be found, Sheltering an occasional Beast or beggar; in the whole Nation, no other relic is left of the Discipline of Geography.
>
> ---From Travels of Praiseworthy Men (1658) by J. A. Suarez Miranda

**`borges`** is a small data visualization package that allows you to plot your single cell RNA-seq dataset as an antique or modern cartographic atlas. It uses 2D coordinates - be it UMAP, tSNE, PCA or anything else - and depicts cell labels as continuous territories in an ocean, separated by rivers or seas.

This is all done through the use of [**`oveRlay`**](github.com/gdagstn/oveRlay), **`ggplot2`** and a few other libraries.

**`borges`** is still very much **under development** so any feedback (especially bug reports) is more than welcome.

## Install

To install **`borges`** you will need to install first **`oveRlay`**.

Use `remotes::install_github()` or `devtools::install_github()` as follows:

```{r}
remotes::install_github("gdagstn/oveRlay")
remotes::install_github("gdagstn/borges")
```

## Usage

**`borges`** has only two functions:

`prepAtlas()` to prepare the atlas coordinates from a `SingleCellExperiment` object, a `matrix` or a `data.frame`, and `plotAtlas()` to display it as a `ggplot2` plot.

For a practical demonstration, let's download a `SingleCellExperiment` object using the `scRNAseq` BioConductor package from Zeisel et al. 2018, "Molecular architecture of the mouse nervous system" ([link](https://doi.org/10.1016/j.cell.2018.06.021)).

This is quite a large file which will take a while to download.

```{r}
# BiocManager::install("scRNAseq")
zeisel = scRNAseq::ZeiselNervousData()
```

The Zeisel dataset has an "unnamed" `reducedDim` slot that contains a t-SNE embedding for cells of the nervous system. There are several labels in the `colData` slot, and we will choose one that offers a good balance between detail and redability.

```{r}
zeiselatlas = prepAtlas(zeisel, 
                        dimred = "unnamed", 
                        res = 400, 
                        labels = "TaxonomyRank3", 
                        as_map = TRUE)
```

The atlas can be plotted:

```{r}
plotAtlas(zeiselatlas)
```

<img src="https://user-images.githubusercontent.com/21171362/220776819-447f8394-0215-4a78-9f35-54c7c0968e0e.png" width="800"/>

Note that plotting can take a few seconds to a minute due to the high level of detail. To have less detailed maps, you can set the `res` argument in `prepAtlas()` to a smaller value, e.g. 250 or 300.

The arguments of `plotAtlas()` allow you to control a few graphical elements:

-   `plot_cells` (logical) to plot cells (as small, semi-transparent dots)

-   `add_contours` (logical) to add 2D kernel density contour estimates, clipped to stay within land masses (mostly)

-   `show_labels` (logical) to show labels using `geom_label_repel()` from `ggrepel`

-   `capitalize_labels` (logical) to capitalize all labels

## Plotting generic 2D point clouds

**`borges`** can also be used on any generic 2D point cloud represented as a `matrix` or `data.frame`, as long as they have two columns (the first one is taken to have coordinates in the x axis, and the second one in the y axis). Moreover, if you are supplying either `matrix` or `data.frame`, the `labels` argument must be a character vector with labels for every point.

```{r}
set.sedd(123)

mats = rbind(matrix(rnorm(1000, 0, 1), ncol = 2),
             matrix(rnorm(1000, 4, 1), ncol = 2),
             matrix(rnorm(1000, -3, 1), ncol = 2))

labels = c(rep("Cluster 1", 500),
           rep("Cluster 2", 500),
           rep("Cluster 3", 500))

atl = prepAtlas(mats, res = 100, labels = labels)

plotAtlas(atl)
```

<img src="https://user-images.githubusercontent.com/21171362/221169182-b38f56bf-f995-43b7-8cb1-909e0360ee26.png" width="400"/>

## Geographical projection

The `as_map` argument in `plotAtlas()` controls whether it will be plotted using a geographical projection, and the `map_proj` argument controls the type of projection. Any one-character argument to `ggplot2::coord_map()` is acceptable.

For instance, we can plot the atlas using a "globular" projection:

```{r}
plotAtlas(zeiselatlas, as_map = TRUE, map_proj = "globular")
```

<img src="https://user-images.githubusercontent.com/21171362/220779384-b4090017-be9e-4761-abb7-e0ba547e6086.png" width="800"/>

## Map themes

**`borges`** comes with a few different themes pre-packaged:

-   **classic**: the default theme

-   **modern**: a modern political atlas-like theme

-   **renaissance**: palette from 16th century maps

-   **medieval**: palette from 14th century maps

```{r}
plotAtlas(zeiselatlas, map_theme = "renaissance")
```

<img src="https://user-images.githubusercontent.com/21171362/220781789-95a75ff2-371a-4c8d-bd1c-6b95b86e4fea.png" width="800"/>

In the future, themes will support different fonts and additional aesthetic elements.
