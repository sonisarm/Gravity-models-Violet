# GRAVITY MODELS FOR INDIVIDUALS #
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
library(devtools)
library(AICcmodavg)
library(openxlsx)
library(tibble)
#remotes::install_github("jeffreyevans/spatialEco") # used for collinearity

##### Add required packages 
p <- c("raster", "igraph", "sp", "GeNetIt", "spatialEco", "leaflet",
       "sf", "terra", "sfnetworks", "spdep", "dplyr", "tmap", "devtools") 
if(any(!unlist(lapply(p, requireNamespace, quietly=TRUE)))) { 
  m = which(!unlist(lapply(p, requireNamespace, quietly=TRUE)))
  suppressMessages(invisible(lapply(p[-m], require,    
                                    character.only=TRUE)))
  stop("Missing library, please install ", paste(p[m], collapse = " "))
} else {
  if(packageVersion("GeNetIt") < "0.1-5") {
    remotes::install_github("jeffreyevans/GeNetIt")
  } 
  suppressMessages(invisible(lapply(p, require, character.only=TRUE)))
}

# Some needed functions
back.transform <- function(y) exp(y + 0.5 * stats::var(y))
rmse = function(p, o){ sqrt(mean((p - o)^2)) }



##### Load coordinates and genetic data
pathtocsv <- "pathtocsv"
df <- read.csv(psate0(pathtocsv, 'df.csv'), sep=';', header=TRUE)
gendif <- read.csv(psate0(pathtocsv, 'gendif.csv'), sep=';', header=TRUE)

# get coords
coords <- df[,2:3] 

# optional: get pop centroid
aggregate(cbind(X, Y) ~ Pop, data = df, FUN = mean)

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
Node.sf <- st_as_sf(Node.df, coords = c("X", "Y"), crs = prj, agr = "constant") 


#### Extract Edge data
# calculate distances
Dist <- knn.graph(Node.sf, row.names=Node.sf$ID)
Dist <- merge(Dist, st_drop_geometry(Node.sf), 
                   by.y="ID", by.x="from_ID") 

# create from-to-pop column
gendif$fromtopop <- paste(gendiv$Ind1, '-', gendiv$Ind2)
gendif$Gflow <- flow(gendif$GD) 

# Get stats from env data
Dist$from.to <- paste(Dist$from_ID, Dist$to_ID, sep=".")
Edge.df <- merge(Dist, gendif, by = "from.to")

suppressWarnings(
  Edge.stats <- graph.statistics(Edge.df, r = rast, 
                                 buffer= NULL, stats = c("median")))

# NOTES. Other stats: "min","mean","max", "var". Buffer can also be changed.

Sites <- st_sf(data.frame(as.data.frame(Edge.df), Edge.stats), 
                         geometry=Edge.df$geometry)

# write gravity model input
st_write(Sites, paste0('edge_data.gpkg'), append=FALSE)


#### Evaluate variable collinearity
vars <- c('vars') # env variables
Vs <- Sites %>% dplyr::select(all_of(vars))

# log transform
Vslog <- lapply(Vs, function(x) if (any(x < 0)) x + 10 else x)
Vslog <- as.data.frame(Vslog)
Vslog <- log(Vslog)

p = 0.7 # Set upper limit of for collinearity
Site.cor <- cor(Vslog, y = NULL, 
                 use = "complete.obs", 
                 method = "pearson")
diag(Site.cor) <- 0			  		
cor.idx <- which(Site.cor > p | Site.cor < -p, arr.ind = TRUE)
Vcor.names <- vector()
Vcor.p <- vector()
for(i in 1:nrow(cor.idx)) {
  cor.p[i] <- Vsite.cor[cor.idx[i,][1], cor.idx[i,][2]]
  cor.names [i] <- paste(rownames(Site.cor)[cor.idx[i,][1]],
                          colnames(Site.cor)[cor.idx[i,][2]], sep="_")
}	

data.frame(parm=Vcor.names, p=Vcor.p) 

# Suggestions of which variables to drop:
node.cor <- spatialEco::collinear(Vslog, p=p) # All correlations are <=0.7

#### Create final Sites df (only with non-collinear variables, choose above)
variab <- intersect(colnames(Vs), colnames(Sites))
variab <- c(variab, 'Gflow')

# log transform
Sites[variab] <- log(Sites[variab])

#### Gravity models

# Null model #
(null <- gravity(y = "Gflow", x = c("length"), d = "length", group = "from_ID", constrained=TRUE, ln=TRUE,
                 data = Sites, fit.method = "ML"))

# Node model # 

(nodevar1 <- gravity(y = "Gflow", x = c("length", "nodevar1"), 
                       d = "length", group = "from_ID", constrained=TRUE, ln=TRUE,
                       data = VgdataTeideIndv, fit.method = "ML"))

(nodevar2 <- gravity(y = "Gflow", x = c("length","nodevar2"),  
                       d = "length", group = "from_ID", constrained=TRUE, ln=TRUE,
                       data = VgdataTeideIndv, fit.method = "ML"))

# Edge model # 
grav.results <- list()

for (variable in variables) {
  colname <- gsub("median\\.", "", variable)
  colname <- paste0(colname, '_1')
  grav.results[[colname]] <- gravity(y = "Gflow", x = c("length", variable), 
                                     d = "length", group = "from_ID", constrained=TRUE, 
                                     data = Sites, fit.method = "ML", ln = FALSE)}

for (variable in variables) {
  colname <- gsub("median\\.", "", variable)
  colname <- paste0(colname, '_2')
  grav.results[[colname]] <- gravity(y = "Gflow", x = c("length", "nodevar1", variable), 
                                     d = "length", group = "from_ID", constrained=TRUE, 
                                     data = Sites, fit.method = "ML", ln = FALSE)}
for (variable in variables) {
  colname <- gsub("median\\.", "", variable)
  colname <- paste0(colname, '_3')
  grav.results[[colname]] <- gravity(y = "Gflow", x = c("length", "nodevar2", variable), 
                                     d = "length", group = "from_ID", constrained=TRUE, 
                                     data = Sites, fit.method = "ML", ln = FALSE)}
for (variable in variables) {
  colname <- gsub("median\\.", "", variable)
  colname <- paste0(colname, '_4')
  grav.results[[colname]] <- gravity(y = "Gflow", x = c("length", "nodevar1", "nodevar2", variable), 
                                     d = "length", group = "from_ID", constrained=TRUE, 
                                     data = Sites, fit.method = "ML", ln = FALSE)}

# Hypothesis and global model # 
envar <- grep("median.", colnames(Sites), value = TRUE)
hypothesis <- list(hypothesis1 = c("var1", "var2"), global=variables)

for (name in names(hypothesis)) {grav.results[[paste0(name, "_1")]] <- gravity(y = "Gflow", x = c("length", hypothesis[[name]]), constrained=TRUE, 
                                                                               d = "length", group = "from_ID", data = Sites, fit.method = "ML", ln = FALSE)}

for (name in names(hypothesis)) {grav.results[[paste0(name, "_2")]] <- gravity(y = "Gflow", x = c("length", "nodevar1", hypothesis[[name]]), constrained=TRUE, 
                                                                               d = "length", group = "from_ID", data = Sites, fit.method = "ML", ln = FALSE)}

for (name in names(hypothesis)) {grav.results[[paste0(name, "_3")]] <- gravity(y = "Gflow", x = c("length", "nodevar2", hypothesis[[name]]), constrained=TRUE, 
                                                                               d = "length", group = "from_ID", data = Sites, fit.method = "ML", ln = FALSE)}

for (name in names(hypothesis)) {grav.results[[paste0(name, "_4")]] <- gravity(y = "Gflow", x = c("length", "nodevar1", "nodevar2", hypothesis[[name]]), constrained=TRUE, 
                                                                               d = "length", group = "from_ID", data = Sites, fit.method = "ML", ln = FALSE)}



#### Compare models
models <- names(grav.results)
comparison <- compare.models(models)
# choose AIC threshold 
AIC = 4
comparison$deltaAIC <- comparison$deltaAIC + 0.001 # avoid infinite by adding constant
compAIC <- comparison %>% dplyr::filter(deltaAIC < 4)
# calculate AICw:
compAIC$likelihoods <- exp(-0.5 * comparison$deltaAIC)
compAIC$akaike_weights <- compAIC$likelihoods / sum(compAIC$likelihoods)


#### Statistics

# Stats of null
par(mfrow=c(2,3))
for (i in 1:6) { plot(null, type=i) } 

# back transform log values
Vgd <- back.transform(Sites$Gflow)


#### Fit final model
grav <- list() 

for (modname in models) {
  mod <- grav.results[[modname]] # get model
  weight <- compAIC$akaike_weights[match(modname, compAIC$model)] # get weight
  g <- trimws(unlist(strsplit(gsub("[()]", "", mod$fixed.formula[3]), " \\+ ", '\\,'))) # get formula
  global_fit <- gravity(y = "Gflow", x = g, d = "length", 
                        group = "from_ID", data = Sites, ln=FALSE, fit.method = "REML")
  #Effect size calculation: global model
  grav[[modname]] <- gravity.es(global_fit)
  par(mfrow=c(2,3))
  for (q in 1:6) { plot(global_fit, type=q) } 
  # Make individual-level (group) predictions (per slope) and show RMSE 
  Vglobal.p <- predict(global_fit, y = "Gflow", x = g,  
                       newdata=Sites, groups = Sites$from_ID,
                       back.transform = "simple")
  cat("RMSE of global", rmse(Vglobal.p, Vgd), "\n")
  
  # Aggregrate estimates and plot ###
  Vglobal.p <- data.frame(EID = Sites$from.to,
                          NID = Sites$from_ID,  
                          p=Vglobal.p) # Step 1: Creating a data frame from predictions
  Vedge.p <- tapply(Vglobal.p$p, Vglobal.p$EID, mean) # Step 2: Calculating mean predictions for each edge
  Vdist20m.graphTeideIndv_results[[paste0("global.flow_",modname)]] <- Vedge.p*weight
  
  
  Vnode.p <- tapply(Vglobal.p$p, Vglobal.p$NID, mean) # Step 3: Calculating mean and variance predictions for each node

  Vglobal.node.var <- tapply(Vglobal.p$p, Vglobal.p$NID, var)
  idx <- which(Sites$ID %in% names(Vnode.p)) # Step 4: Assigning mean and variance to "Sites"
  Sites[[paste0("global.flow_",modname)]][idx] <- Vnode.p*weight
  Sites[[paste0("global.var_",modname)]][idx] <- Vglobal.node.var*weight

    Vglobal.pTO <- predict(global_fit, y = "Gflow", x = g,  
                         newdata=VgdataTeideIndv, groups = VgdataTeideIndv$to_ID,
                         back.transform = "simple")
  Vglobal.pTO <- data.frame(EID = VgdataTeideIndv$from.to,
                            NID = VgdataTeideIndv$to_ID,  
                            p=Vglobal.pTO)
  Vnode.pTO <- tapply(Vglobal.pTO$p, Vglobal.pTO$NID, mean)
  Vglobal.node.varTO <- tapply(Vglobal.pTO$p, Vglobal.pTO$NID, var)
  Sites[[paste0("global.flow_", modname)]][5] <- Vnode.pTO[4]*weight
  Sites[[paste0("global.var_", modname)]][5] <- Vglobal.node.varTO[4]*weight
  Sites[[paste0("global.var_", modname)]][1] <- Vglobal.node.varTO[1]*weight
}


Sites <- Sites %>% filter(ID %in% df$ID) 

# Calculate final model (global flow and variance are proportions according to AICw)
Sites <- Sites %>%rowwise() %>% mutate(global.flow = sum(c_across(contains('global.flow_'))))
Sites <- Sites %>%  rowwise() %>% mutate(global.var = sum(c_across(contains('global.var')))) 
Sites <- Sites %>%  rowwise() %>% mutate(global.flow = sum(c_across(contains('global.flow'))))

# save results
st_write(Sites, paste0('Results_', AIC, '.gpkg'), append=FALSE)
st_write(Node.sf, paste0('Node_', AIC, '.gpkg'), append=FALSE)


#### Env effect size and output as excel workbook
wb <- createWorkbook()
for(sheet_name in names(grav)) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet = sheet_name, grav[[sheet_name]], rowNames = TRUE)
}

saveWorkbook(wb, paste0("wb_", AIC, ".xlsx"), overwrite = TRUE)


# Get t-values and p-values
tval <- do.call(rbind, lapply(grav, function(df) {
  df %>% rownames_to_column(var = "variable")}))

tval_mean <- tval %>% group_by(variable) %>% summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))



#### Plot
pal <- colorRampPalette(rev(c("red","orange","blue")), bias=0.15)

Map <-  tm_shape(Sites) +
  tm_lines("global.flow", palette=pal(10), title.col="Edges: global.flow") +
  tm_shape(Node.sf) +
  tm_symbols(col = "global.flow", size = "global.var", 
             shape = 20, scale = 3, palette=pal(10), 
             title.col="Nodes: global.flow", title.size="Nodes: global.var") 

# static map
tmap_mode(c("plot", "view")[1])
Map + tm_layout(legend.outside=TRUE, legend.position = c("right", "top"),
                legend.text.size=1.5, asp = 0.7,legend.outside.size = 0.25)


#### Optional: Plot high flow only

hf.sites <- Sites %>% dplyr::filter(global.flow >0.65) 
hf.node <- Node.sf %>% dplyr::filter(ID %in% hf.sites$from_ID)

table(hf.sites$Pop)

pal <- colorRampPalette(rev(c("red","orange")), bias=0.15)

Map <- 
  tm_shape(hf.sites) +
  tm_lines("global.flow", palette=pal(5), title.col="Edges: global.flow") +
  tm_shape(hf.node) +
  tm_symbols(col = "global.flow", size = "global.var", 
             shape = 20, scale = 3, palette=pal(10), 
             title.col="Nodes: global.flow", title.size="Nodes: global.var") 

tmap_mode(c("plot", "view")[1])
Map + tm_layout(legend.outside=TRUE, legend.position = c("right", "top"),
                legend.text.size=1.5, asp = 0.7,legend.outside.size = 0.25)


