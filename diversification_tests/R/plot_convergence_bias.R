# Generate trace plots demonstrating convergence

library(ggplot2)
library(patchwork)

# Read in samples
D1 <- read.csv("../treeppl_bias/tdbd_bias_lez_125000-100_1.csv")
D2 <- read.csv("../treeppl_bias/tdbd_bias_lez_125000-100_2.csv")
D3 <- read.csv("../treeppl_bias/tdbd_bias_lez_125000-100_3.csv")
D4 <- read.csv("../treeppl_bias/tdbd_bias_lez_125000-100_4.csv")
D1$run <- "run1"
D2$run <- "run2"
D3$run <- "run3"
D4$run <- "run4"
D1$iteration <- seq(0,125000,by=100)
D2$iteration <- seq(0,125000,by=100)
D3$iteration <- seq(0,125000,by=100)
D4$iteration <- seq(0,125000,by=100)
D <- rbind(D1,D2)
D <- rbind(D,D3)
D <- rbind(D,D4)

# Read in likelihood values
E1 <- read.csv("../treeppl_bias/tdbd_bias_lez_125000-100_1_logLike.csv")
E2 <- read.csv("../treeppl_bias/tdbd_bias_lez_125000-100_2_logLike.csv")
E3 <- read.csv("../treeppl_bias/tdbd_bias_lez_125000-100_3_logLike.csv")
E4 <- read.csv("../treeppl_bias/tdbd_bias_lez_125000-100_4_logLike.csv")
E1$run <- "run1"
E2$run <- "run2"
E3$run <- "run3"
E4$run <- "run4"
E1$iteration <- 1:nrow(E1)
E2$iteration <- 1:nrow(E2)
E3$iteration <- 1:nrow(E3)
E4$iteration <- 1:nrow(E4)
E <- rbind(E1,E2)
E <- rbind(E,E3)
E <- rbind(E,E4)

# Define plot fxn
plot_fxn <- function(D, col, ylab) {

    D <- D[,c("run", "iteration", col)]
    colnames(D)[3] <- "val"
    D <- D[D$iteration >= 25000,]
    ggplot(D, aes(x=iteration, y=val, color=run)) +
        theme_minimal() + 
        geom_line(alpha=0.5) +
        labs(color="MCMC analysis") +
        ylab(ylab) +
        xlab("Iteration") +
        ggtitle(NULL)
}

# Render trace plots
p1 <- plot_fxn(E, "logLike", "log likelihood")
p2 <- plot_fxn(D, "lambda", "lambda")
p3 <- plot_fxn(D, "epsilon", "epsilon")
p4 <- plot_fxn(D, "z", "z")

# Save plots
ggsave(file = "../figs/Fig_convergence_bias.jpg",
       width = 7.5,
       height = 7.5,
       plot = p1 + p2 + p3 + p4 +
           plot_layout(ncol=2, guides="collect") +
           plot_annotation(tag_levels="A") &
           theme (legend.position="bottom")
       )

