############## Network Analysis for Housing Quality Indicators and Depression #########################

## This code was adapted by Faye Sanders, from Lucy Waldren's Network Analysis Code: https://github.com/Lucy-Wal/Autism_ADHD_Adults

# Load libraries
library(bootnet)
library(dplyr)
library(ggplot2)
library(NetworkComparisonTest)
library(networktools)
library(qgraph)


# Load data

df <- load("....") # your data here

# Only keep housing and depression variables in the dataframe
df <- dplyr::select(df, Problems_house_a, Facil_house_a, Total_Num_rooms_a, depression_age30, Temp_house_a, a092, Mum_bedroom_dec_a, Kitchen_dec_a, Living_room_dec_a, Other_room_dec_a)
# Here we are creating a housing decorations variable by grouping decoration variables together
df$Decorations_house_a <- rowMeans(df[, c("Mum_bedroom_dec_a", "Kitchen_dec_a", "Living_room_dec_a", "Other_room_dec_a")])
# And removing the item-level variables we just grouped together
df <- subset(df, select = -c(Mum_bedroom_dec_a, Kitchen_dec_a,Living_room_dec_a,Other_room_dec_a))

# Shorten variable names so they are easier to display on the graphs
df <- df %>%
  rename(
    Decorations = Decorations_house_a,
    Problems = Problems_house_a,
    Facilities = Facil_house_a,
    Size = Total_Num_rooms_a,
    Depression = depression_age30,
    Temperature = Temp_house_a,
    Feelings = a092
  )


# Compute the correlation matrix
cor_matrix <- cor(df, use = "pairwise.complete.obs", method = "spearman")
print(cor_matrix)


## Running an EBICGlasso network

network <- estimateNetwork(df, default = "EBICglasso", corMethod = "spearman",tuning = 0.5, weighted = T)

edge_colors <- ifelse(cor_matrix > 0, "steelblue2", "red") # Positive edges = blue, negative edges = red
edge_widths <- abs(cor_matrix) * 5 # Scale up the edges so that they are more easily visible

# Plot the network analysis
network_plot <- qgraph(cor_matrix, 
                       layout = "spring",   
                       graph = "glasso",    
                       sampleSize = nrow(df), 
                       tuning = 0.5,        # Tuning parameter for regularization
                       labels = FALSE,
                       color = c("coral", "magenta2", "palevioletred1", "skyblue", "mistyrose1", "plum", "mediumpurple1"), # Specifying node colours
                       edge.color = edge_colors,
                       edge.width = edge_widths,
                       label.cex = 2,
                       vsize = 8)


# Bootstrapping
nonpara_boot <- bootnet(network, nBoots = 1000, nCores = 6) 
plot(nonpara_boot, order = "sample", labels = T)                                
plot(nonpara_boot, plot = "difference", onlyNonZero = TRUE, order = "sample") 

# Describe edges in the network
all_edges <- network$graph 
print("All Edges:"); print(psych::describe(all_edges)) 
print(all_edges)
edge_weights <- all_edges[all_edges != 0]
non_zero_edges <- sum(all_edges != 0)
mean_edge_weight <- mean(edge_weights)
SD_edge_weight <- sd(edge_weights)

#Central Stability 
casedrop_boot <- bootnet(network, nBoots = 1000, caseMax = 1, nCores = 8, type = "case", statistics = c("betweenness", "closeness", "expectedInfluence")) 
plot(casedrop_boot, statistics = c( "betweenness", "closeness", "expectedInfluence")) 
corStability(casedrop_boot) 

community_vector <- c(1, 1, 1, 2, 1, 1, 1) # 1's = housing variables, 2 = depression (depends on the order of your variables in the dataframe)

bridges <- bridge(network$graph, communities = community_vector)                 #Generate Bridge EI centrality 
bridges$Bridge_EI_z <- scale(bridges$`Bridge Expected Influence (1-step)`) 
plot(bridges, include=c("Bridge Expected Influence (1-step)"), zscore = TRUE, labels = F)
print(bridges)

