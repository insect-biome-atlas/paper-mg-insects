# Make combined figure with different types of validation plots

# Retrieve plots from other scripts
# (generating separate figure files
# for them in the process)
source("plot_filter_res.R")
p1 <- plot_A
source("plot_placement_validation.R")
p2 <- plot_A
source("plot_validation.R")
p3 <- plot_B
p4 <- plot_C
source("plot_validation_by_order.R")
p5 <- plot_match(merge(D,E2),"right")   # To get the legend
p6 <- plot_match(merge(D,E1),"none")    # To remove the legend

# Save combined figure
ggsave("../figs/Fig_combined_validation.jpg",
       device="jpg",
       width=14.0,
       height=16.0,
       units="in",
#       dpi=300,
       plot = p1 + p2 + p3 + p4 + p5 + p6 +
            plot_layout(heights=c(0.3,0.4,0.4), ncol=2) +
            plot_annotation(tag_levels="A")
        )

