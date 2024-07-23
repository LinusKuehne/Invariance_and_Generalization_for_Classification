# In this script, we group the wildfires into five and nine distinct environments, resp.
# Furthermore, we create the world maps with the visualization of the different environments.


library(ggplot2)
library(ggmap)
library(gridExtra)
library(ggpubr)



# get the path of this script
script_dir <- getwd()


# load in the dataset
load(file.path(script_dir, "../../data/exported_pyrocb.rdata"))


# insert authentification key for stadiamaps 
# you can get your own key here: https://client.stadiamaps.com/signup/
register_stadiamaps(key = "7ad859a8-6d8c-4433-a91f-b2ec1451043d")





#-------------------------------------------------------------------------------
# obtain the maps of North America (AM) and Australia (AUS)
#-------------------------------------------------------------------------------

# observations collected in AUS
ind.aus <- which(envVars[,1]>100)

# observations collected in AM
ind.am <- which(envVars[,1]<0)

# extract lat/lon for AUS wildfires
loc_aus <- envVars[ind.aus, 1:2]

# extract lat/lon for AM wildfires
loc_am <- envVars[ind.am, 1:2]

# define a bounding box for the AUS map
bbox_aus <- c(left = 112, bottom = -44.72, right = 155, top = -8.23)

# define a bounding box for the AM map
bbox_am <- c(left = -137, bottom = 15, right = -69, top = 63)

# download the AUS map
aus_map <- get_stadiamap(bbox_aus, maptype = "stamen_terrain", zoom = 5)

# download the AM map
am_map <- get_stadiamap(bbox_am, maptype = "stamen_terrain", zoom = 5)

# generate a plot of the maps if interested
#ggmap(aus_map)
#ggmap(am_map)

#-------------------------------------------------------------------------------






#-------------------------------------------------------------------------------
# group the wildfires into 9 environments respecting natural ecoregions
#-------------------------------------------------------------------------------

# number of environments
num.env <- 9

env <- rep(num.env, nrow(cube))
biome <- rep("A", nrow(cube))
env9 <- rep("A", nrow(cube))

# env1 (top NA)
i1 <- intersect(ind.am[which(loc_am$V2>50)], ind.am[which(loc_am$V1> -115)])
env[i1] <- 1
biome[i1] <- "1: Boreal forest"
env9[i1] <- "env1"


# env2 (NW NA)
i2 <- intersect(ind.am[which(loc_am$V2>46)], ind.am[which(loc_am$V1< -112)])
env[i2] <- 2
biome[i2] <- "2: Temperate conifer forests"
env9[i2] <- "env2"

# env3 (int NA)
i3 <- intersect(ind.am[which(loc_am$V2<48)], ind.am[which(loc_am$V2>44)])
i3 <- setdiff(i3,i2)
env[i3] <- 3
biome[i3] <- "3: Temperate shrublands"
env9[i3] <- "env3"

# env4 (Cali N NA)
i4 <- intersect(ind.am[which(loc_am$V2<44)], ind.am[which(loc_am$V2>39)])
env[i4] <- 4
biome[i4] <- "4: Temperate conifer forests"
env9[i4] <- "env4"

# env5 (Cali NA)
i5 <- intersect(ind.am[which(loc_am$V2<39)], ind.am[which(loc_am$V1< -116.5)])
env[i5] <- 5
biome[i5] <- "5: Mediterranean forest and scrub"
env9[i5] <- "env5"

# env6 (inland NA)
i6 <- intersect(ind.am[which(loc_am$V2<39)], ind.am[which(loc_am$V2>37)])
i6 <- setdiff(i6,i5)
env[i6] <- 6
biome[i6] <- "6: Shrublands / temperate conifer forest"
env9[i6] <- "env6"

# env6 continued (AZ NA) 
i7 <- ind.am[which(env[ind.am] == num.env)]
env[i7] <- 7
biome[i7] <- "6: Shrublands / temperate conifer forest"
env9[i7] <- "env6"


# env7 (SW AUS) 
i8 <- intersect(ind.aus[which(loc_aus$V2< -31)], ind.aus[which(loc_aus$V1<135)])
env[i8] <- 8
biome[i8] <- "7: Mediterranean forest and scrub"
env9[i8] <- "env7"


# env8 (SW AUS)
i9 <- intersect(ind.aus[which(loc_aus$V2> -31)], ind.aus[which(loc_aus$V1<135)])
i9 <- union(i9, ind.aus[which(loc_aus$V2>-20)])
env[i9] <- 9
biome[i9] <- "8: Xeric shrubland and savanna"
env9[i9] <- "env8"


# env9 (SW AUS) 
i10 <- ind.aus[which(loc_aus$V1>143)]
env[i10] <- 10
biome[i10] <- "9: Temperate broadleaf & mixed forests"
env9[i10] <- "env9"



# store the grouping in a factor
env9 <- as.factor(env9)

#-------------------------------------------------------------------------------







#-------------------------------------------------------------------------------
# plot these nine environments
#-------------------------------------------------------------------------------


# plot colors
cols <- c("brown4", "black",  "hotpink", "cyan", "yellow", "mediumorchid","blue1", "red", "darkorange")


loc_aus$Biome <- factor(biome[ind.aus])
loc_am$Biome <- factor(biome[ind.am])


# plot for AM
p_am <- ggmap(am_map) + 
  geom_point(data = loc_am, aes(x = V1, y = V2, color = Biome), size = 1.5) +
  geom_point(data = loc_aus, aes(x = V1, y = V2, color = Biome), size = 1.5) +
  scale_color_manual(values=cols[1:num.env]) +
  xlab(expression(paste("Longitude (", degree,"E)"))) +
  ylab(expression(paste("Latitude (", degree,"N)"))) +
  labs(title = "Wildfires in North America", color = "Environment and biome")+
  theme(legend.background = element_blank(), legend.key = element_blank())+
  guides(color = guide_legend(ncol = 3))

#p_am




# plot for AUS
p_aus <- ggmap(aus_map) + 
  geom_point(data = loc_aus, aes(x = V1, y = V2, color = Biome), size = 1.5) +
  scale_fill_discrete(breaks=c("8: Mediterranean Forest and Scrub", "9: Xeric Shrubland and Savanna", "10: Temperate Broadleaf & Mixed Forests")) +
  scale_color_manual(values=cols[7:9]) +
  xlab(expression(paste("Longitude (", degree,"E)"))) +
  ylab(expression(paste("Latitude (", degree,"N)"))) +
  labs(title = "Wildfires in Australia", color = "Environment and biome") +
  theme(legend.background = element_blank(), legend.key = element_blank()) +
  guides(color = guide_legend(ncol = 3))

#p_aus

# combine plots
g <- ggarrange(p_am+ theme(legend.title = element_blank()), p_aus+ theme(legend.title = element_blank()), ncol=2, nrow=1, common.legend = TRUE, legend="bottom") +
  guides(color = guide_legend(ncol = 4))

#g

ggsave(filename = file.path(script_dir, "../saved_plots/map9.pdf"), width = 8, height = 5.5)

 #-------------------------------------------------------------------------------







#-------------------------------------------------------------------------------
# group the wildfires into 5 environments respecting natural ecoregions
#-------------------------------------------------------------------------------

# number of environments
num.env <- 5

# use regional clusters already defined for pyrocast ICP paper
env <- event_df$cluster_regional+1

# switch some environments around
i.1 <- which(env == 1)
i.4 <- which(env == 4)

env[i.1] <- 4
env[i.4] <- 1

i.3 <- which(env == 3)
i.5 <- which(env == 5)

env[i.3] <- 5
env[i.5] <- 3


biome <- rep("A", nrow(cube))
env5 <- rep("A", nrow(cube))

for(cl in 1:5){
  ind <- which(env == cl)
  if(cl == 1){
    biome[ind] <- "1: Boreal and conifer forests"
    env5[ind] <- "env1"
  }
  if(cl == 2){
    biome[ind] <- "2: Temperate conifer forests"
    env5[ind] <- "env2"
  }
  if(cl == 3){
    biome[ind] <- "3: Mediterranean and conifer forests"
    env5[ind] <- "env3"
  }
  if(cl == 4){
    biome[ind] <- "4: Shrublands"
    env5[ind] <- "env4"
  }
  if(cl == 5){
    biome[ind] <- "5: Mixed forests"
    env5[ind] <- "env5"
  }
}



# store this as a factor
env5 <- as.factor(env5)


#-------------------------------------------------------------------------------







#-------------------------------------------------------------------------------
# plot these five environments
#-------------------------------------------------------------------------------


# plotting colors
cols <- c("black",  "hotpink", "darkorange", "red", "blue1")


loc_aus$Biome <- factor(biome[ind.aus])
loc_am$Biome <- factor(biome[ind.am])


# plot for AM
p_am <- ggmap(am_map) + 
  geom_point(data = loc_am, aes(x = V1, y = V2, color = Biome), size = 1.5) +
  geom_point(data = loc_aus, aes(x = V1, y = V2, color = Biome), size = 1.5) +
  scale_color_manual(values=cols[1:num.env]) +
  xlab(expression(paste("Longitude (", degree,"E)"))) +
  ylab(expression(paste("Latitude (", degree,"N)"))) +
  labs(title = "Wildfires in North America", color = "Environment and biome")+
  theme(legend.background = element_blank(), legend.key = element_blank())+
  guides(color = guide_legend(ncol = 3))

#p_am




# plot for AUS
p_aus <- ggmap(aus_map) + 
  geom_point(data = loc_aus, aes(x = V1, y = V2, color = Biome), size = 1.5) +
  scale_fill_discrete(breaks=c("8: Mediterranean Forest and Scrub", "9: Xeric Shrubland and Savanna", "10: Temperate Broadleaf & Mixed Forests")) +
  scale_color_manual(values=cols[4:5]) +
  xlab(expression(paste("Longitude (", degree,"E)"))) +
  ylab(expression(paste("Latitude (", degree,"N)"))) +
  labs(title = "Wildfires in Australia", color = "Environment and biome") +
  theme(legend.background = element_blank(), legend.key = element_blank())+
  guides(color = guide_legend(ncol = 3))
#p_aus

# combine plots
g <- ggarrange(p_am+ theme(legend.title = element_blank()), p_aus+ theme(legend.title = element_blank()), ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
#g


ggsave(filename = file.path(script_dir, "../saved_plots/map5.pdf"), width = 8, height = 5.5)

#-------------------------------------------------------------------------------



# store the environments
save(env5, env9, file = file.path(script_dir, "../saved_data/discrete_envs.rdata"))






writeLines(capture.output(sessionInfo()), file.path(script_dir, "../sessionInfo/create_environments.txt"))







