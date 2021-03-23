# utility
some utility scripts


## Volcano plot

```Shell
python volcano_plot.py -d volcano_data.csv --title "Volcano plot"
```
volcano_plot.py reads volcano_data.csv and generates a figure, which looks like
<img src="https://github.com/chpngyu/utility/blob/main/volcano_plot.jpg">

Give -h will show more other options
```Shell
python volcano_plot.py -h
```


## PCA

Compute and draw Principal Component Analysis (PCA) in 2D or 3D. 


## Histogram

Compute and draw the histogram for matrix data. This script allows user to provide your own data in several formats, including of csv, tsv, excel, ..etc.
Then user can select the subset of data to plot and set the regions and bins to plot. The dimension of the histogram are avaiable for 1D and 2D.

## Heatmap and Clustering

Plot a matrix dataset as a hierarchically-clustered heatmap.

