# Set directory
setwd("/Users/claramal/Desktop/Moncla_Lab/H5N5-seabird-phylo/analyses/")
getwd()

# Load packages
library(seraphim)
library(diagram)
library(geodata)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(sf)
library(terra)
library(dplyr)

mostRecentSamplingDatum = 2026.0520547945205

# Plot the CSV and save the plot as a png to the output_plots directory (Chat GPT helped adjust this script)

mcc_file <- "./A6_concat-gen_traits_2026-08-04_MCC.csv"
if (!dir.exists("seraphim_plots")) dir.create("seraphim_plots")

for (file_path in mcc_file)
{
  k <- read.csv(file_path)
  k$isnode <- ifelse(k$tipLabel == "", "node", "tip") # This line adds coloring depending on whether there is a node or a tip
  prob = 0.95; precision = 0.025
  startDatum = min(k[,"startYear"]) 
  
  # Option 1 - color scale for hosts
  host_colors <- c("#6aa1f4", "#f69ba6", "#6a994e", "#a874d2", "#923e4d", "#e9c46a")
  host_names <- c("seab", "bdpr", "gds", "mammal", "scav", "grdbd")
  names(host_colors) <- host_names
  
  # Option 2- - color scale for tips vs. nodes
  node_colors <- c("red", "blue")
  ext_node_shapes <- c(1, 0) # open
  int_node_shapes <- c(16, 15) # filled in
  node_names <- c("node", "tip")
  names(node_colors) <- node_names
  names(ext_node_shapes) <- node_names
  names(int_node_shapes) <- node_names

  # Option 3 - color scale by ending year of a branch (time at end node of branch)
  colour_scale = colorRampPalette(brewer.pal(11,"RdYlBu"))(141)[21:121] # This line provides a gradient of 100 colors from red > yellow > green
  minYear = min(k[,"endYear"]); maxYear = max(k[,"endYear"])
  endYears_indices = (((k[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
  endYears_colours = colour_scale[endYears_indices]

  # Define a background (map) for the plot
  northern_hemi <- ne_countries(continent = c("north america", "europe"), returnclass = "sf")
  northern_hemi_vect <- vect(northern_hemi)
  custom_ext <- ext(-165, 45, 36, 83) # if you want to make custom dimensions in the raster, use this code
  r_template <- rast(custom_ext, resolution = 0.3) # This code implements any custom dimensions for your raster
  template_raster <- rasterize(northern_hemi_vect, r_template, field = 1)
  
  countries = c("USA", "Canada", "Greenland", "Iceland", "Spain", "Ireland", "France", 
                "Italy", "Sweden", "Norway", "Denmark", "Portugal", "United Kingdom",
                "Netherlands", "Belgium", "Switzerland", "Germany", "Austria",
                "Hungary", "Greece", "Russia", "Finland", "Poland", "Croatia",
                "Ukraine", "Svalbard and Jan Mayen", "Estonia", "Czechia", "Slovakia",
                "Slovenia", "Albania", "Serbia", "Montenegro", "Kosovo", "Macedonia",
                "Bulgaria", "Romania", "Moldova", "Belarus", "Lithuania", "Latvia")
  borders = gadm(country = countries, level = 0, path = "/Users/claramal/Desktop/Moncla_Lab/H5-avian-ecology/SERAPHIM/nov2021-2024_no-suliformes/")
  
  # Define output PNG filename
  output_file <- paste0("plot_", tools::file_path_sans_ext(basename(file_path)), "_jittered.pdf")
  pdf(file = file.path("seraphim_plots", output_file), width = 6, height = 6.3)
  par(mar=c(0,0,0,0), oma=c(1.2,3.5,1,0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
  
  # Now plot the map background
  plot(template_raster, col="lightgrey", box=F, axes=F, colNA="white", legend=F)
  plot(borders, add=T, lwd=0.1, border="white")
  
  # This loop plots the points and lines
  # According to CHAT GPT: color is based on endLon/endLat, not startLon/startLat, except for a one-time marker at the first iteration
  for (i in dim(k)[1]:1)
  {
    # Ssets a jitter for the x and y values for ending locations, so they're easier to distinguish
    tip_x_jit <- ifelse(k[i, "isnode"] == "node", k[i, "endLon"], jitter(k[i, "endLon"], amount = 1))
    tip_y_jit <- ifelse(k[i, "isnode"] == "node", k[i, "endLat"], jitter(k[i, "endLat"], amount = 1))
    
    # Plot curved lines that end on the jittered x, y coordinates
    curvedarrow(cbind(k[i,"startLon"],k[i,"startLat"]), cbind(tip_x_jit,tip_y_jit), arr.length=0,
                arr.width=0, lwd=0.2, lty=1, lcol="gray32", arr.col=NA, arr.pos=FALSE, curve=0.1, dr=NA, endhead=F)
    
    # Starting nodes
    #points(k[i,"startLon"], k[i,"startLat"], pch=16, col=endYears_colours[i], cex=0.3) # Louise wanted the starting points colored too, so this line and the one below do that
    points(k[i,"startLon"], k[i,"startLat"], pch=16, col = host_colors[(k[i, "startHost"])], cex=0.3) # This line colors by host group
    #points(k[i,"startLon"], k[i,"startLat"], pch=16, col = node_colors[(k[i, "isnode"])], cex=0.3) # This line colors by whether there is a tipLabel or not
    points(k[i,"startLon"], k[i,"startLat"], pch=1, col="gray10", cex=0.3) # This line of code provides a dark outline for each starting point
    
    # Ending nodes
    #points(tip_x_jit, tip_y_jit, pch=int_node_shapes[(k[i,"isnode"])], col=endYears_colours[i], cex=0.3) # These next two lines color the ending points, this one colors by time
    #points(tip_x_jit, tip_y_jit, pch=ext_node_shapes[(k[i, "isnode"])], col="gray10", cex=0.3) # This line of code provides an dark outline for each ending point
    points(tip_x_jit, tip_y_jit, pch=int_node_shapes[(k[i,"isnode"])], col = host_colors[(k[i, "endHost"])], cex = 0.3) # This line jitters the colored internal part of ending points
    points(tip_x_jit, tip_y_jit, pch=ext_node_shapes[(k[i,"isnode"])], col = "gray10", cex = 0.3) # This line jitters the external outline of ending points
    
  }
  
  # Plot legends
  rect(xmin(template_raster), ymin(template_raster), xmax(template_raster), ymax(template_raster), xpd=T, lwd=0.2)
  axis(1, c(ceiling(xmin(template_raster)), floor(xmax(template_raster))), pos=ymin(template_raster), mgp=c(0,0.2,0), cex.axis=0.5, lwd=0, lwd.tick=0.2, padj=-0.8, tck=-0.01, col.axis="gray30")
  axis(2, c(ceiling(ymin(template_raster)), floor(ymax(template_raster))), pos=xmin(template_raster), mgp=c(0,0.5,0), cex.axis=0.5, lwd=0, lwd.tick=0.2, padj=1, tck=-0.01, col.axis="gray30")
  legend("bottom", legend = names(host_colors), col = host_colors, pch = 16, bty = "n", cex = 0.7, pt.cex = 0.6) # Host legend
  
  # The following codes make a legend for the ending times of each branch
  #rast = raster(matrix(nrow=1, ncol=2)); rast[1] = min(k[,"endYear"]); rast[2] = max(k[,"endYear"]) # I want the scale to only show endYear colors
  #plot(rast, legend.only=T, add=T, col=colour_scale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.40,0.80,0.14,0.15),
  #     legend.args=list(text="", cex=0.7, line=0.3, col="gray30"), horizontal=T,
  #     axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2, tck=-0.5, col.axis="gray30", line=0, mgp=c(0,-0.02,0)))
  #legend("bottom", legend = names(node_colors), col = node_colors, pch = 16, bty = "n", cex = 0.7, pt.cex = 0.6)
  
  dev.off()
}
