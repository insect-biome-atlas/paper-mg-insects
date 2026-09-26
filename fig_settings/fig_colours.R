# Common settings for figures

library(viridis)

# Create series of distinct colorblindness-friendly color sets
my_cols <- viridis_pal()(5)
mg_col <- my_cols[2]
rf_col <- my_cols[1]
df_col <- my_cols[3]
se_col <- my_cols[4]
diamond_col <- my_cols[5]
forest_cols4 <- c(rf_col, df_col, mg_col, se_col)

# Colors from viridis_pal()(6)
my_cols <- viridis_pal()(6)
rf_col <- my_cols[3]        # dark-green-ish "#2A788EFF"
df_col <- my_cols[5]        # light-green-ish "#7AD151FF"
se_col <- "deepskyblue"     # works better for box plots, also in combination with the above
diamond_col <- my_cols[6]   # yellow-ish "#FDE725FF"
forest_cols3 <- c(rf_col, df_col, se_col)

# Two-set colours (used in simple contrasts)
mg_col <- "darkgreen"
se_col <- "blue"

rf_col <- my_cols[3]        # dark-green-ish "#2A788EFF"
df_col <- my_cols[5]        # light-green-ish "#7AD151FF"

