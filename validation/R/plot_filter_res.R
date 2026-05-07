# Plot filter results

library(dplyr)
library(ggplot2)
library(patchwork)

plot_filter_res <- function(df) {

    ggplot() +
        geom_line (data=df, aes(x=step_no, y=remaining_frac, color=country, linetype=country)) +
        geom_point(data=df, aes(x=step_no, y=remaining_frac, color=country, shape=country), size=3) +
        theme_minimal() +
        theme(plot.margin = unit(c(0.5, 0.5, 1.5, 0.5), "cm"),
              axis.title.x=element_text(margin=margin(t=10, r=0, b=0, l=0),size=12),
              axis.title.y=element_text(margin=margin(t=0, r=10, b=0, l=0),size=12)) +
        scale_color_manual(values = c("green","blue")) +
        guides(
            color = guide_legend(title="Country"),      # Country color legend
            shape = guide_legend(title="Country"),      # Country shape legend
            linetype = guide_legend(title="Country")    # Country linetype legend
        ) +
        scale_y_continuous(limits=c(0.0,1.0)) +
        ylab("Fraction remaining OTUs") +
        xlab("Filtering step")
}

# Read in and augment result data frame
filter_res_mg <- read.delim("../../iba_data/filter_annotation_stats_mg.tsv")
filter_res_se <- read.delim("../../iba_data/filter_annotation_stats_se.tsv")
filter_res_mg$step_no <- as.numeric(row.names(filter_res_mg))
filter_res_se$step_no <- as.numeric(row.names(filter_res_se))
filter_res_mg$country <- "Madagascar"
filter_res_se$country <- "Sweden"
filter_res <- rbind(filter_res_mg,filter_res_se)

# Generate plots
p1 <- plot_filter_res(filter_res)

ggsave("../figs/Fig_filter_res.jpg", width=6.0, height=5.0, plot=p1)
