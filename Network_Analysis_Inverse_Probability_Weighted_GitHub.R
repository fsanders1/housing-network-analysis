# --- Inverse Probability Weighted Network analysis ---

library(dplyr)
library(weights)
library(qgraph)
library(bootnet)
library(networktools)

# 1. Load original (non-imputed) dataset
raw_df <- readRDS("df.rds")

# Variables to check selection completeness (only actual selection-relevant variables, e.g., baseline housing)
vars_for_selection <- c("Problems_house_a", "Facil_house_a", "Total_Num_rooms_a",
                        "Temp_house_a", "a092",
                        "Decorations_house_a",
                        "contextual_age30", "age_at_depression_age30")

# Create selection indicator: S = 1 if complete on these baseline variables
raw_df$S <- ifelse(complete.cases(raw_df[, vars_for_selection]), 1, 0)

# Keep only ID + S for merging
S_df <- raw_df %>% select(cidB2957, S)

# 2. Load imputed dataset
load("imputed_df.RData")
imp_df <- dataset
rm(dataset)

# 3. Rename ID in imputed dataset before merging
imp_df <- imp_df %>%
  rename(cidB2957 = cidB2957.x)

# 4. Apply inclusion criteria for analysis (e.g., <50% missing housing items)
imp_df <- imp_df %>%
  filter(housing_a_missingness < 50,
         housing_f_missingness < 50,
         housing_g_missingness < 50)

# 5. Merge selection indicator into imputed dataset
df <- merge(imp_df, S_df, by = "cidB2957", all.x = TRUE)

# 6. Fit selection model using observed covariates (predicting S)
# In this case, we use SES and age
sel_model <- glm(S ~ age_at_depression_age30 + contextual_age30, 
                 family = binomial(),
                 data = df)

# 7. Predicted probability of being observed
df$ps <- predict(sel_model, newdata = df, type = "response")

# 8. Compute stabilized IP weights
pS <- mean(df$S == 1, na.rm = TRUE)
df$ipw <- ifelse(df$S == 1, 
                 pS / df$ps, 
                 (1 - pS) / (1 - df$ps))

# --- Weighted network analysis ---

# Keep analysis variables
df_analysis <- df %>%
  select(Problems_house_a, Facil_house_a, Total_Num_rooms_a,
         depression_age30, Temp_house_a, a092,
         Decorations_house_a, contextual_age30, age_at_depression_age30)

# Rename variables for nicer plotting
df_analysis <- df_analysis %>%
  rename(
    Decorations = Decorations_house_a,
    Problems = Problems_house_a,
    Facilities = Facil_house_a,
    Size = Total_Num_rooms_a,
    Depression = depression_age30,
    Temperature = Temp_house_a,
    Feelings = a092,
    SES_Risk = contextual_age30,
    Age_at_Depressive_Symptoms = age_at_depression_age30
  )

# Rank the variables for Spearman correlations
df_ranked <- as.data.frame(lapply(df_analysis, rank, na.last = "keep"))

# Compute weighted Spearman correlation matrix
cor_matrix_w <- wtd.cor(df_ranked, weight = df$ipw)
cor_matrix_w_mat <- cor_matrix_w$cor

edge_colors <- ifelse(cor_matrix_w_mat > 0, "steelblue2", "red") # Positive edges = blue, negative edges = red
edge_widths <- abs(cor_matrix_w_mat) * 5

# Build network from weighted correlations
network_w <- qgraph(cor_matrix_w_mat, 
                    layout = "spring",   
                    graph = "glasso",    
                    sampleSize = nrow(df), 
                    tuning = 0.5,  
                    edge.color = edge_colors,
                    labels = FALSE,
                    color = "white")

network_w <- EBICglasso(cor_matrix_w_mat, n = nrow(df), gamma = 0.5)

edge_df <- as.data.frame(network_w)
diag(edge_df) <- "--"
print(edge_df)

# Mean edge weights
edge_values <- network_w[upper.tri(network_w)]
edge_values_nonzero <- edge_values[edge_values != 0]
mean_edge_weight <- mean(edge_values_nonzero)
sd_edge_weight   <- sd(edge_values_nonzero)
min_edge_weight  <- min(edge_values_nonzero)
max_edge_weight  <- max(edge_values_nonzero)
nonzero_edges <- length(edge_values_nonzero)
cat("Nonzero edges:", nonzero_edges, "\n")
cat("Mean edge weight:", mean_edge_weight, "\n")
cat("SD edge weight:", sd_edge_weight, "\n")
cat("Range:", min_edge_weight, "to", max_edge_weight, "\n")


# Bootnet requires a network object from estimateNetwork
network_for_boot <- estimateNetwork(df_ranked, 
                                    default = "EBICglasso", 
                                    corMethod = "cor_auto",  # uses Spearman if df_ranked
                                    tuning = 0.5)

# 1. Non-parametric bootstrap
nonpara_boot <- bootnet(network_for_boot, nBoots = 1000, nCores = 6, type = "nonparametric")
plot(nonpara_boot, order = "sample", labels = TRUE)
plot(nonpara_boot, plot = "difference", onlyNonZero = TRUE, order = "sample")

# 2. Case-dropping bootstrap (centrality stability)
casedrop_boot <- bootnet(network_for_boot, 
                         nBoots = 1000, 
                         caseMin = 0.1, 
                         caseMax = 0.894, 
                         nCores = 8, 
                         type = "case", 
                         statistics = c("betweenness", "closeness", "expectedInfluence")) 

plot(casedrop_boot, statistics = c("betweenness", "closeness", "expectedInfluence"))

# 3. Bridge expected influence
community_vector <- c(1, 1, 1, 2, 1, 1, 1, 3, 3) # Adjust to your actual communities
bridges <- bridge(network_w, communities = community_vector)
plot(bridges, include = c("Bridge Expected Influence (1-step)"), labels = FALSE)
print(bridges)

cs_values <- corStability(casedrop_boot)
print(cs_values)
