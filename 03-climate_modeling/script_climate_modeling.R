# CLIMATE CHANGE PROJECTION #
  # AUTHOR: SONIA SARMIENTO #
  # PROJECT: TEIDE VIOLETS #
      # NOVEMBER 2024 #


##### Load libraries
library(terra)
library(raster)
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ade4)
library(tmap)
library(openxlsx)
library(RColorBrewer)


# Some needed functions
back.transform <- function(y) exp(y + 0.5 * stats::var(y))
rmse = function(p, o){ sqrt(mean((p - o)^2)) }

# Define variables
cmodel <- c('csiro', 'miroc')
years <- c(30, 50, 80)
rcp <- c(26, 60, 85)


##### Load climate change scenario files (separately - or within loop, step 3.)
for (c in cmodel){
  for (r in rcp){
    for (y in years){
    # load raster
    mypath <- paste0("YourPath", c, "\\", r, "\\", y)
    file <- paste0(mypath, '\\yourfile.tif')
    raster <- rast(file)
    # Assign raster to a variable in the global environment
    var_name <- paste0(c, "_RCP", r, "_Y", y)  # Create a unique variable name
    assign(var_name, raster, envir = .GlobalEnv)  # Save to global environment
    }
  }
}


#####  Load coordinates 
pathtocsv <- "pathtocsv"
df <- read.csv(psate0(pathtocsv, '00-data/df.csv'), sep=';', header=TRUE)
gendif <- read.csv(psate0(pathtocsv, '00-data/gendif.csv'), sep=';', header=TRUE)

# get coords
coords <- df[,2:3] 

# GD -> gene flow
gendif$Gflow <- flow(gendif$GD)


##### Get statistics at node & between sites (per model)
for (c in cmodel){
  for (r in rcp){
    for (y in years){
  # Load raster (climatic variables - future)
  raster <- rast("futureraster.tif") # change name for every model / rcp / year
  
  # Extract data at node:
  Node.data <- terra::extract(raster, coords)
  Node.df <- merge(Node.data, coords)
  Node.sf <- st_as_sf(Node.df, coords = c("X", "Y"), crs = prj, agr = "constant") # add X Y coords and change to sf object

  # Extract data at edge:
  # create saturated graph
  Dist <- knn.graph(Node.sf, row.names = Node.sf$ID) # ID or pops, depending on your study sample
  Dist <- merge(Dist, st_drop_geometry(Node.sf), 
                by.y="ID", by.x="from_ID")
  
  # merge with gene flow
  Dist$from.to <- paste(Dist$from_ID, Dist$to_ID, sep=".")
  Edge.df <- merge(Dist, gendif, by = "from.to")
  
  # Get statistics of edges
  suppressWarnings(
    Edge.stats <- graph.statistics(Edge.df, r = raster, 
                                   buffer= NULL, stats = c("median")))

  Sites <- st_sf(data.frame(as.data.frame(Edge.df), Edge.stats), 
                           geometry=Edge.df$geometry)
  
  # Load topographic variables (unchanged)
  Topo <- st_read('edge_data.gpkg') # see 02-script_grav_models
  topovar <- c('var1', 'var2') # see 02-script_grav_models
  Topo <- Topo %>% dplyr::select(all_of(topovar))
  
  Sitesfut <- st_join(Sites, Topo, join = st_nearest_feature, left = T) # merges by matching geometry

  # 3.5 Write final file
  st_write(Sitesfut, paste0("edge_data_future_", c, "_rcp", r, "_year", y, ".gpkg"), append=FALSE) # append=FALSE to rewrite data
    }
  }
}

##### Run gravity models
var <- c('Gflow', 'var1', 'var2') # see 02-script_grav_models

for (c in cmodel){
  for (r in rcp){
    for (y in years){
      # Load raster data
      future.stack <- rast("futureraster.tif")
      # Load St file (future)
      Sites <-  st_read("FileName.gpkg")
      # Subset data 
      Sites <- Sites %>% dplyr::select('from.to', 'from_ID',
                                       'i', 'j', 'to_ID', 'Gflow',
                                       'length', all_off(variab))
      
      # NOTE: add constant to any negative value. For example:
      Sites$env <- Sites$env + 1
      
      # Log transform values
      Sites[variab] <- log(Sites[variab])
      
      ### Gravity models #
      # see 02-script_grav_models for chosen models - in Teide violet's case is global model
      envar <- grep("median.", colnames(Sites), value = TRUE)

      
      # Define the different variable combinations to use in the gravity function
      var_comb <- list(
        c("length", envar),
        c("length", "nodevar1", envar),
        c("length", "nodevar2", envar),
        c("length", "nodevar1", "nodevar2", envar))
      
      grav.results <- list()
    
      # Run gravity function for each combination, storing results as "global_1", "global_2", etc.
      for (i in seq_along(var_comb)) {
        grav.results[[paste0("global_", i)]] <- gravity(
          y = "Gflow", x = var_comb[[i]], constrained = TRUE, 
          d = "length", group = "from_ID", data = Sites, fit.method = "ML", ln = FALSE
        )
      }
      
      
      # Compare models
      model_names <- names(grav.results)
      comparison <- compare.models(null, node1, node2, model_names)
 

      # Stats of null model
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
        Sites[[paste0("global.flow_",modname)]] <- Vedge.p*weight
        
        
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
      st_write(Sites, paste0('Results_', AIC, '_', c, "_rcp", r, "_year", y,'.gpkg'), append=FALSE)
      st_write(Node.sf, paste0('Node_', AIC, '_', c, "_rcp", r, "_year", y,'.gpkg'), append=FALSE)
      
      #### Env effect size and output as excel workbook
      wb <- createWorkbook()
      for(sheet_name in names(grav)) {
        addWorksheet(wb, sheet_name)
        writeData(wb, sheet = sheet_name, grav[[sheet_name]], rowNames = TRUE)}
      
      saveWorkbook(wb, paste0("wb_", AIC, ".xlsx"), overwrite = TRUE)
  
      # Get t-values and p-values
      tval <- do.call(rbind, lapply(grav, function(df) {
        df %>% rownames_to_column(var = "variable")}))
      
      tval_mean <- tval %>% group_by(variable) %>% summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))
      
    }
  }
}



##### Comparison (future - present)
# Load df
Edgep <- st_read('Results_4.gpkg') # see 02-script_grav_models
Nodep <- st_read('Node_4.gpkg')
# Filter
Nodep<- Nodep %>% dplyr::select(ID, global.flow, global.var) # define dfs to save the results! 
Edgep <- Edgep %>% dplyr::select(from.to, global.flow)

# Read files:
file_data <- list() # Create an empty list to store the variables
directory <- "yourdirectory"

for (c in cmodel){
  for (r in rcp){
    for (y in years){
    # Construct the file name
    files <- list.files(directory)
    file_name_edge <- paste0('Results_', AIC, '_', c, "_rcp", r, "_year", y,'.gpkg')
    file_name_node <- paste0('Node_', AIC, '_', c, "_rcp", r, "_year", y,'.gpkg')
    # Check if the file exists in the directory
    if (file_name_node %in% files) {
      # Read the file and assign it to the list
      tmp <- st_read(file.path(directory, file_name_node))
      assign(paste0("Node_", c, "_", r, "_Y",y), tmp)
      
    } else {
      warning(sprintf("File %s not found in the directory.", file_name_node))
    }
    # Check if the file exists in the directory
    if (file_name_edge %in% files) {
      # Read the file and assign it to the list
      tmp <- st_read(file.path(directory, file_name_edge))
      assign(paste0("Edge_", c, "_", r, "_Y",y), tmp)
      
    } else {
      warning(sprintf("File %s not found in the directory.", file_name_edge))
      }
    }
  }
}


# Calculate difference future - present
for (c in cmodel){
  for (r in rcp){
    for (y in years){
    w <- stringr::str_extract(c, "^.{1}") # 2 first characters of model name
    w
    #z <- stringr::str_extract(y, ".{3}$") # number of spp
    df_name <- paste0("Node_", c, "_", r, "_Y",y) #construct the df's name
    df_name
    df <- get(df_name)    # Access the dataframe using the constructed name
    df <- df %>% dplyr::filter(ID %in% Nodep$ID)
    
    Nodep[[paste0("gflow_", w, r, "_Y", y)]] <- df$global.flow - Nodep$global.flow
    Nodep[[paste0("gvar_", w, r, "_Y", y)]] <- df$global.var - Nodep$global.var
    Nodep
    df_name <- paste0("Edge_",  c, "_", r, "_Y",y) #construct the df's name
    df <- get(df_name)    # Access the dataframe using the constructed name
    
    df <- df %>% dplyr::filter(from.to %in% Edgep$from.to)
    Edgep[[paste0("gflow_", w, r, "_Y", y)]] <- df$global.flow - Edgep$global.flow
    }
  }
}

# Plot difference
pal <- colorRampPalette(rev(c("red","orange","blue")), bias=0.15)
for (c in cmodel){
  for (r in rcp){
   for (y in years){
  w <- stringr::str_extract(cmodel, "^.{1}") # 2 first characters of model name
  w
  tf <-paste0("gflow_",  w, r, "_Y", y)
  tv <-paste0("gvar_", w, r, "_Y", y)

  #Map
  Map <- 
    tm_shape(Edgep) +
    tm_lines(tf, palette=pal(10), title.col="Edges: global.flow") +
    tm_shape(Nodep) +
    tm_symbols(col = tf, size = tv, 
               shape = 20, scale = 2, palette=pal(10), 
               title.col="Nodes: global.flow", title.size="Nodes: global.var") +
    tm_layout(title = paste("RCP",  r, ", Y", y), title.snap.to.legend=TRUE,
              title.position = c("right", "top"), title.size = 1,
              legend.outside = TRUE, legend.position = c("right", "top"),
              legend.text.size = 0.7, asp = 0.7, legend.outside.size = 0.25)
  

  # Save the plot as a file
  plot_file <- paste0("plot_", c,"_AIC", AIC, "_RCP", r, "_Y", y, ".png")
  tmap_save(Map, filename = plot_file, width = 8, height = 6, units = "in", dpi = 300)
  
    }
  }
}


# plots are saved with the code above.
