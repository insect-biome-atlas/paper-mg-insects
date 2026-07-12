# Plot filtering results

library(dplyr)
library(ggplot2)
library(patchwork)

plot_filter_res <- function(df) {

    ggplot() +
        geom_line (data=df, aes(x=step_no, y=remaining_frac, color=country, linetype=country)) +
        geom_point(data=df, aes(x=step_no, y=remaining_frac, color=country, shape=country), size=3) +
        theme_minimal(base_size=16) +
        theme(plot.margin = unit(c(0.5, 0.5, 1.5, 0.5), "cm"),
              axis.title.x=element_text(margin=margin(t=10, r=0, b=0, l=0),size=16),
              axis.title.y=element_text(margin=margin(t=0, r=10, b=0, l=0),size=16),
              legend.position="right") +
        scale_color_manual(values = c("green","blue")) +
        guides(
            color = guide_legend(title="Country"),      # Country color legend
            shape = guide_legend(title="Country"),      # Country shape legend
            linetype = guide_legend(title="Country")    # Country linetype legend
        ) +
        scale_y_continuous(limits=c(0.0,1.0)) +
        scale_x_continuous(breaks=c(2,4,6,8)) +
        ylab("Fraction remaining OTUs") +
        xlab("Filtering step")
}

# Read in and augment result data frame
filter_res_mg <- read.delim("../../iba_data/filter_annotation_stats_mg.tsv")
filter_res_se <- read.delim("../../iba_data/filter_annotation_stats_lysate_litter_se.tsv")
filter_res_mg$step_no <- as.numeric(row.names(filter_res_mg))
filter_res_se$step_no <- as.numeric(row.names(filter_res_se))
filter_res_mg$country <- "Madagascar"
filter_res_se$country <- "Sweden"
filter_res <- rbind(filter_res_mg,filter_res_se)

# Generate plots
plot_A <- plot_filter_res(filter_res)

ggsave("../figs/Fig_filter_res.jpg", width=9.0, height=7.5, plot=plot_A)
