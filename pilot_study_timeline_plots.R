# Load the data  
masterdata <- read.csv("/Users/jasonklemm/Desktop/other/PilotsAMREE.csv")  

# Create subsets with proper indexing  
MRP1 <- subset(masterdata, Project_code == "MRP1")  
MRP1$Sample_index <- 1:nrow(MRP1)  

SMX1 <- subset(masterdata, Project_code == "SMX1")  
SMX1$Sample_index <- 1:nrow(SMX1)  

SMX2 <- subset(masterdata, Project_code == "SMX2")  
SMX2$Sample_index <- 1:nrow(SMX2)  

SMX3 <- subset(masterdata, Project_code == "SMX3")  
SMX3$Sample_index <- 1:nrow(SMX3)  

# Assign unique colors for each mutation event
mutation_colors_mrp1 <- c("purple", "orange", "red", "forestgreen")

# Helper function to get MIC at a given Sample_index
get_mic_at_index <- function(data, idx) {
  if (idx %in% data$Sample_index) {
    return(data$MIC[data$Sample_index == idx][1])
  } else {
    return(NA)
  }
}

# Open PNG device for a high quality figure
png(filename = "/Users/jasonklemm/Desktop/combplot.png", 
    width = 3000, height = 2250, res = 300, pointsize = 12)

par(mfrow = c(2,2), mar = c(4.5, 4.5, 4, 1))

# Updated plotting function with smaller points and black outlines
plot_project <- function(data, project_name, line_color, mutation_lines = NULL, mutation_labels = NULL, mutation_colors = NULL) {
  plot(data$Sample_index, data$MIC, type = "l", col = line_color,
       lwd = 2, xlab = "Time (Days)", ylab = "MIC (ug/mL)",
       main = project_name)
  lines(data$Sample_index, data$AB_conc, col = "black", lty = 2, lwd = 2)
  
  legend_items <- c("MIC", "AB_conc")
  legend_colors <- c(line_color, "black")
  legend_types <- c(1, 2)
  legend_width <- c(2, 2)
  legend_pch <- c(NA, NA)
  legend_pt_bg <- c(NA, NA)
  
  # Add points at the MIC value for mutations, with unique colors and black outlines
  if (!is.null(mutation_lines) && !is.null(mutation_labels)) {
    for (i in seq_along(mutation_lines)) {
      mic_val <- get_mic_at_index(data, mutation_lines[i])
      if (!is.na(mic_val)) {
        # Use pch=21 for filled circles with black borders, bg for fill color
        points(mutation_lines[i], mic_val, col = "black", bg = mutation_colors[i], pch = 21, cex = 1.2)
      }
      legend_items <- c(legend_items, mutation_labels[i])
      legend_colors <- c(legend_colors, "black")
      legend_types <- c(legend_types, NA)
      legend_width <- c(legend_width, NA)
      legend_pch <- c(legend_pch, 21)  # Use same point style in legend
      legend_pt_bg <- c(legend_pt_bg, mutation_colors[i])  
    }
  }
  
  legend("bottomright", 
         legend = legend_items,
         col = legend_colors,
         pt.bg = legend_pt_bg,
         lty = legend_types,
         lwd = legend_width,
         pch = legend_pch,
         pt.cex = c(rep(NA, 2), rep(1.2, length(mutation_labels))),  # Point size in legend
         cex = 0.8,  # legend text size
         bty = "n")
}

# Plot MRP1 with all four mutation points and unique colors
plot_project(MRP1, "MRP1", "blue", mutation_lines = c(16, 8, 26, 11), mutation_labels <- c(  
  expression(italic("mrdA") ~ "p.G505R"),  
  expression(italic("mrdA") ~ "p.T547I"),  
  expression(italic("ftsI") ~ "p.R167L"),  
  expression(italic("ftsI") ~ "p.A513S")  
)  , mutation_colors = mutation_colors_mrp1)
# SMX1, SMX2: folP p.F28L label, unique color
plot_project(SMX1, "SMX1", "green", mutation_lines = 16, mutation_labels = expression(italic("folP") ~ "p.F28L"), mutation_colors = "purple")
plot_project(SMX2, "SMX2", "red", mutation_lines = c(3,13), mutation_labels = c(expression(italic("folP") ~ "p.P64S"), expression(italic("folP") ~ "p.F28L")), mutation_colors = c("orange", "purple"))
# SMX3: folP p.P64S label, unique color
plot_project(SMX3, "SMX3", "orange", mutation_lines = 5, mutation_labels = expression(italic("folP") ~ "p.P64S"), mutation_colors = "red")

dev.off()

library(png)
img <- readPNG("/Users/jasonklemm/Desktop/combplot.png")
plot.new()
rasterImage(img,0,0,1,1)