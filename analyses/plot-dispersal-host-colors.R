# Set directory
setwd("/Users/claramal/Desktop/Moncla_Lab/H5N5-seabird-phylo/SERAPHIM/")
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

# Step 1 - set trees directory (this assumes that you have already extracted information for the posterior distribution)
# localTreesDirectory = "./Extracted_trees" # You don't have a distribution of trees for this, so you shouldn't need this line
# nberOfTreesToSample = 96 # See above comment
mostRecentSamplingDatum = 2026.0520547945205

# Step 2 - extract spatio-temporal information from whole MCC tree
# Use the CSV you made yourself, with host information included
mcc_tab = read.csv("./A6_concat-gen_traits_2026-08-04_MCC.csv", head=T)


# Step 3 - Plot the CSV and save the plot as a png to the output_plots directory (Chat GPT helped a lot with this)
mcc_file <- "./A6_concat-gen_traits_2026-08-04_MCC.csv"
if (!dir.exists("seraphim_plots")) dir.create("seraphim_plots")

for (file_path in mcc_file)
{
  k <- read.csv(file_path)
  k$isnode <- ifelse(k$tipLabel == "", "node", "tip") # This line adds coloring depending on whether there is a node or a tip
  #nberOfExtractionFiles = nberOfTreesToSample
  prob = 0.95; precision = 0.025
  startDatum = min(k[,"startYear"]) 
#  polygons = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob, startDatum, precision))
  
  
  # Define a color scale
  
  host_colors <- c("#6aa1f4", "#f69ba6", "#6a994e", "#a874d2", "#923e4d", "#e9c46a")
  host_names <- c("seab", "bdpr", "gds", "mammal", "scav", "grdbd")
  names(host_colors) <- host_names
  
  node_colors <- c("red", "blue")
  ext_node_shapes <- c(1, 0) # open
  int_node_shapes <- c(16, 15) # filled in
  node_names <- c("node", "tip")
  names(node_colors) <- node_names
  names(ext_node_shapes) <- node_names
  names(int_node_shapes) <- node_names

  colour_scale = colorRampPalette(brewer.pal(11,"RdYlBu"))(141)[21:121] # This line provides a gradient of 100 colors from red > yellow > green
  #minYear = min(k[,"startYear"]); maxYear = max(k[,"endYear"])
  minYear = min(k[,"endYear"]); maxYear = max(k[,"endYear"])
  endYears_indices = (((k[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
  endYears_colours = colour_scale[endYears_indices]
#  polygons_colours = rep(NA, length(polygons))
#  for (i in 1:length(polygons))
#  {
#    date = as.numeric(names(polygons[[i]]))
#    polygon_index = round((((date-minYear)/(maxYear-minYear))*100)+1)
#    polygons_colours[i] = paste0(colour_scale[polygon_index],"40")
#  }

  # Step 6 - plot HPD regions and MCC tree
  # Start by defining a background (map) for the plot
  northern_hemi <- ne_countries(continent = c("north america", "europe"), returnclass = "sf")
  northern_hemi_vect <- vect(northern_hemi)
  custom_ext <- ext(-165, 45, 36, 83) # if you want to make custom dimensions in the raster, use this code
  #custom_ext <- ext(-165, 165, 15, 83) # This is because I want to see whats going on in Russia and Japan
  r_template <- rast(custom_ext, resolution = 0.3) # This code implements any custome dimensions for your raster
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
  #output_file <- paste0("russia_test_plot.pdf")
  pdf(file = file.path("seraphim_plots", output_file), width = 6, height = 6.3)
  par(mar=c(0,0,0,0), oma=c(1.2,3.5,1,0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
  plot(template_raster, col="lightgrey", box=F, axes=F, colNA="white", legend=F)
  
#  for (i in 1:length(polygons))
#  {
#    plot(polygons[[i]], axes=F, col=polygons_colours[i], add=T, border=NA)
#  }
  plot(borders, add=T, lwd=0.1, border="white")
  
#  for (i in 1:dim(k)[1]) # I removed this loop because I moved into the points loop
#  {
    #curvedarrow(cbind(k[i,"startLon"],k[i,"startLat"]), cbind(k[i,"endLon"],k[i,"endLat"]), arr.length=0,
    #            arr.width=0, lwd=0.2, lty=1, lcol="gray32", arr.col=NA, arr.pos=FALSE, curve=0.1, dr=NA, endhead=F)
#  }
  
  # This loop plots the points
  # According to CHAT GPT: color is based on endLon/endLat, not startLon/startLat, except for a one-time marker at the first iteration
  for (i in dim(k)[1]:1)
  {
#     if (i == 1) # This "if" loop colors the first starting point
#    {
#      points(k[i,"startLon"], k[i,"startLat"], pch=16, col=colour_scale[1], cex=0.3)
#      points(k[i,"startLon"], k[i,"startLat"], pch=1, col="gray10", cex=0.3)
#    }
#    points(k[i,"startLon"], k[i,"startLat"], pch=16, col="gray10", cex=0.13) # But this line ensure no other starting points are colored
    
    tip_x_jit <- ifelse(k[i, "isnode"] == "node", k[i, "endLon"], jitter(k[i, "endLon"], amount = 1))
    tip_y_jit <- ifelse(k[i, "isnode"] == "node", k[i, "endLat"], jitter(k[i, "endLat"], amount = 1))
    
    # Lines
    # The ammendment below should end the arrow on the jittered x, y coordinates
    #curvedarrow(cbind(k[i,"startLon"],k[i,"startLat"]), cbind(k[i,"endLon_jit"],k[i,"endLat_jit"]), arr.length=0,
    #            arr.width=0, lwd=0.2, lty=1, lcol="gray32", arr.col=NA, arr.pos=FALSE, curve=0.1, dr=NA, endhead=F)
    curvedarrow(cbind(k[i,"startLon"],k[i,"startLat"]), cbind(tip_x_jit,tip_y_jit), arr.length=0,
                arr.width=0, lwd=0.2, lty=1, lcol="gray32", arr.col=NA, arr.pos=FALSE, curve=0.1, dr=NA, endhead=F)
    
    # Starting nodes
    #points(k[i,"startLon"], k[i,"startLat"], pch=16, col=endYears_colours[i], cex=0.3) # Louise wanted the starting points colored too, so this line and the one below do that
    points(k[i,"startLon"], k[i,"startLat"], pch=16, col = host_colors[(k[i, "startHost"])], cex=0.3) # This line colors by host group
    #points(k[i,"startLon"], k[i,"startLat"], pch=16, col = node_colors[(k[i, "isnode"])], cex=0.3) # This line colors by whether there is a tipLabel or not
    points(k[i,"startLon"], k[i,"startLat"], pch=1, col="gray10", cex=0.3) # This line of code provides a dark outline for each starting point
    
    # Ending nodes
    #points(tip_x_jit, tip_y_jit, pch=int_node_shapes[(k[i,"isnode"])], col=endYears_colours[i], cex=0.3) # These next two lines color the ending points, this one colors by time
    #points(k[i,"endLon"], k[i,"endLat"], pch=int_node_shapes[(k[i, "isnode"])], col = host_colors[(k[i, "endHost"])], cex=0.3) # This line colors ending points by host group
    #points(k[i,"endLon"], k[i,"endLat"], pch=16, col = node_colors[(k[i, "isnode"])], cex=0.3) # This line colors by whether there is a tipLabel or not
    #points(tip_x_jit, tip_y_jit, pch=ext_node_shapes[(k[i, "isnode"])], col="gray10", cex=0.3) # This line of code provides an dark outline for each ending point
  
    points(tip_x_jit, tip_y_jit, pch=int_node_shapes[(k[i,"isnode"])], col = host_colors[(k[i, "endHost"])], cex = 0.3) # This line jitters the colored internal part of ending points
    points(tip_x_jit, tip_y_jit, pch=ext_node_shapes[(k[i,"isnode"])], col = "gray10", cex = 0.3) # This line jitters the external outline of ending points
    
  }
  
  # The following 6 lines plot ending points so that ending nodes are plotted first, and ending tips are plotted on top of that
  # ChatGPT said to keep these outside of the loop, but I am not fully sure why
  # This doesn't work
  #nodes <- k[, "isnode"] == "node"
  #tips <- k[, "isnode"] == "tip"
  
  #points(k[nodes,"x_jit"], k[nodes,"y_jit"], pch=int_node_shapes[(k[nodes,"isnode"])], col = host_colors[(k[nodes, "endHost"])], cex = 0.3) # This line jitters the colored internal part of ending points
  #points(k[nodes,"x_jit"], k[nodes,"y_jit"], pch=ext_node_shapes[(k[nodes,"isnode"])], col = "gray10", cex = 0.3) # This line jitters the external outline of ending points
  
  #points(k[tips,"x_jit"], k[tips,"y_jit"], pch=int_node_shapes[(k[tips,"isnode"])], col = host_colors[(k[tips, "endHost"])], cex = 0.3) # This line jitters the colored internal part of ending points
  #points(k[tips,"x_jit"], k[tips,"y_jit"], pch=ext_node_shapes[(k[tips,"isnode"])], col = "gray10", cex = 0.3) # This line jitters the external outline of ending points
  
  rect(xmin(template_raster), ymin(template_raster), xmax(template_raster), ymax(template_raster), xpd=T, lwd=0.2)
  axis(1, c(ceiling(xmin(template_raster)), floor(xmax(template_raster))), pos=ymin(template_raster), mgp=c(0,0.2,0), cex.axis=0.5, lwd=0, lwd.tick=0.2, padj=-0.8, tck=-0.01, col.axis="gray30")
  axis(2, c(ceiling(ymin(template_raster)), floor(ymax(template_raster))), pos=xmin(template_raster), mgp=c(0,0.5,0), cex.axis=0.5, lwd=0, lwd.tick=0.2, padj=1, tck=-0.01, col.axis="gray30")
  legend("bottom", legend = names(host_colors), col = host_colors, pch = 16, bty = "n", cex = 0.7, pt.cex = 0.6)
  #rast = raster(matrix(nrow=1, ncol=2)); rast[1] = min(k[,"startYear"]); rast[2] = max(k[,"endYear"]) # This line ascribes the dates to the color scale
  #rast = raster(matrix(nrow=1, ncol=2)); rast[1] = min(k[,"endYear"]); rast[2] = max(k[,"endYear"]) # I want the scale to only show endYear colors
  #plot(rast, legend.only=T, add=T, col=colour_scale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.40,0.80,0.14,0.15),
  #     legend.args=list(text="", cex=0.7, line=0.3, col="gray30"), horizontal=T,
  #     axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2, tck=-0.5, col.axis="gray30", line=0, mgp=c(0,-0.02,0)))
  #legend("bottom", legend = names(node_colors), col = node_colors, pch = 16, bty = "n", cex = 0.7, pt.cex = 0.6)
  
  dev.off()
}
