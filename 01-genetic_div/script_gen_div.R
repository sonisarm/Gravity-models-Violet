# GENETIC DIVERSITY VS ENVIRONMENTAL VAR #
    # AUTHOR: SONIA SARMIENTO #
      # PROJECT: TEIDE VIOLETS #
          # NOVEMBER 2024 #

##### Load libraries
library(LandGenCourse)
library(ggplot2)
library(corMLPE)
library(terra)
library(raster)
library(sf)
library(dplyr)
library(tidyr)
library(stringr)


##### Load coordinates and genetic data
pathtocsv <- "pathtocsv"
df <- read.csv(psate0(pathtocsv, '00-data/df.csv'), sep=';', header=TRUE)
gendiv <- read.csv(psate0(pathtocsv, '00-data/gendiv.csv'), sep=';', header=TRUE)

# get coords
coords <- df[,2:3] 

##### Load raster
pathtoraster <- "pathtoraster"
rast <- rast(paste0(pathtoraster, "/raster.tif"))

# optional: crop raster
vector <- vect(Coord, geom = c("X", "Y")) 
crs(vector) <- crs(rast)
buffer <- terra::buffer(vector, width = 5000)
rast <- terra::crop(rast, buffer)


#### Extract Node data
Node <- terra::extract(rast, Coords) # also works with vector
Node <- Node %>% dplyr::select(-ID)
Node.df <- cbind(df, Node)

# Calculate mean environmental variable per population
Node.pop <- Node.df %>%
  group_by(Pop) %>%
  summarise(across(3:23, mean, na.rm = FALSE))



#### Run linear models (all variables)
length(Node.pop)
variables <- colnames(Node.pop[3:23])# Define the variables to be tested
pval <- data.frame(Variable = character(), # Initialize an empty data frame to hold the p vals
                      P_Value = numeric(),
                      Significant = logical(),
                      stringsAsFactors = FALSE)

# Define significance level threshold
alpha <- 0.05
pval <- data.frame()

# Loop through each variable
for (variable in variables) {
  # Create the formula dynamically
  formula <- as.formula(paste("AR ~", variable))
  
  # Fit the regression model
  model <- lm(formula, data = VPop)
  
  # Get the summary
  summary_model <- summary(model)
  if (!is.null(summary_model$coefficients) && 
      nrow(summary_model$coefficients) >= 2 && 
      ncol(summary_model$coefficients) >= 4) {
    
    # Extract the p-value
    p_value <- summary_model$coefficients[2, 4]
    est <- summary_model$coefficients[2, 1]
    sterr <- summary_model$coefficients[2, 2]
    tval <- summary_model$coefficients[2, 3]
    
    # Check significance
    significant <- p_value < alpha
    
    # Append the pval to the data frame
    pval <- rbind(pval, data.frame(Variable = variable,
                                   Estimate=est,
                                   StandardError=sterr,
                                   Tvalue=tval,
                                   P_Value = p_value,
                                   Significant = significant))
  } else {
    # Handle case where coefficients table is not as expected
    cat("Error: Unable to extract coefficients for variable", variable, "\n")
  }
  
}



# Filter env. var with p < threshold
sig_pval <- pval[pval$Significant == TRUE, ]


#### Check collinearity between variables
Vsite.cor <- cor(Node.pop[, c(sig_vars)], y = NULL, 
                 use = "complete.obs", 
                 method = "pearson")
diag(Vsite.cor) <- 0			  		
p <- 0.7 # set threshold

# create dataframe with correlations
Vcor.idx <- which(Vsite.cor > p | Vsite.cor < -p, arr.ind = TRUE)
Vcor.names <- vector()
Vcor.p <- vector()
for(i in 1:nrow(Vcor.idx)) {
  Vcor.p[i] <- Vsite.cor[Vcor.idx[i,][1], Vcor.idx[i,][2]]
  Vcor.names [i] <- paste(rownames(Vsite.cor)[Vcor.idx[i,][1]],
                          colnames(Vsite.cor)[Vcor.idx[i,][2]], sep="_")
}	
# print dataframe
data.frame(parm=Vcor.names, p=Vcor.p) #pairwise correlation between length/Tmean/snow
Vnode.cor <- spatialEco::collinear(VPop[3:23], p=p) 

#### Plot
graphics.off()
# var1
model <- lm(AR ~ var1, data = Node.pop)
plot(Node.pop$var1, Node.pop$AR, xlab='var1', ylab='AR')
abline(model, col = "black", lwd = 1)
label_positions <- c(3, 1, 3, 3, 3, 3)  # 1 = below, 2 = left, 3 = above, 4 = right
text(Node.pop$var1, Node.pop$AR, labels = Node.pop$Pop, pos = label_positions, cex = 0.8, col = "black")
