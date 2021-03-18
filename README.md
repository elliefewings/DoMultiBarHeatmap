# DoMultiBarHeatmap
Create multibar annotated heatmap for Seurat object

Adapted from Seurat and work by Arjun Arkal Rao written in a [github discussion](https://github.com/satijalab/seurat/issues/2201).

## Installation

Requires R packages:  
* ggplot
* rlang
* Seurat
* magrittr

### Install with devtools
```
devtools::install_github("elliefewings/DoMultiBarHeatmap")
```
### Or download and install from parent directory
```
devtools::install("DoMultiBarHeatmap")

#Load library
library(DoMultiBarHeatmap)
```

## Usage
```
data("pbmc_small"
DoMultiBarHeatmap(object = pbmc_small, group.by="cluster", additional.group.by="cell_type")
```
### Arguments

| Argument        | Definition  |
| ------------- |:-------------|
| **object** | Seurat object |
| **features** | A vector of features to plot, defaults to VariableFeatures(object = object)   |
| **cells** | A vector of cells to plot   |
| **group.by** | A vector of variables to group cells by; defaults to cell identity classes |
| **additional.group.by** | A second vector of variables to group cells by |
| **group.bar** | Add a color bar showing group status for cells (boolean) |
| **group.colors** | Colors to use for the color bar |
| **disp.min** | Minimum display value (all values below are clipped) |
| **disp.max** | Maximum display value (all values above are clipped); defaults to 2.5 if slot is 'scale.data', 6 otherwise |
| **slot** | Data slot to use, choose from 'raw.data', 'data', or 'scale.data' |
| **assay** | Assay to pull from |
| **label** | Label the cell identies above the color bar |
| **size** | Size of text above color bar |
| **hjust** | Horizontal justification of text above color bar |
| **angle** | Angle of text above color bar |
| **raster** | If true, plot with geom_raster, else use geom_tile. geom_raster may look blurry on some viewing applications such as Preview due to how the raster is interpolated. Set this to FALSE if you are encountering that issue (note that plots may take longer to produce/render). |
| **draw.lines** | Include white lines to separate the groups |
| **lines.width** | Integer number to adjust the width of the separating white lines. Corresponds to the number of "cells" between each group. |
| **group.bar.height** | Scale the height of the color bar |
| **combine** | Combine plots into a single patchworked ggplot object. If FALSE, return a list of ggplot objects |

