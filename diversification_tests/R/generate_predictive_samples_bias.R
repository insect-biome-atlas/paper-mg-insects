# Generate samples of otus from predictive distribution of lambda, epsilon and z
# Account for sampling bias of Coleoptera
rho <- 0.6648029

# Function for transforming time
transformTime <- function(time, z) {
    return (1.0 - exp(-z*time))/z;
}

# Sample no. surviving lineages after specified time
sampleN <- function(time, lambda, mu) {

    ert <- exp(-(lambda-mu)*time)
    p0 <- mu*(1.0 - ert)/(lambda - mu*ert)
    f <- (lambda - mu) / (lambda - mu*ert);
    p <- f*f*ert
    cum <- p0 + p

    u <- runif(1, p0, 1.0)

    n <- 1
    g <- (lambda*(1.0 - ert))/(lambda - mu*ert)
    while(cum < u) {
        p <- p*g
        cum <- cum + p
        n <- n + 1
    }
    return (n)
}

# Read in data
D <- read.delim("../data/diversification_data.tsv");
D$age <- max(D$placement_age,1e-6)

# Read in samples from posterior
P1 <- read.csv("../treeppl_bias/tdbd_bias_lez_125000-100_1.csv");
P2 <- read.csv("../treeppl_bias/tdbd_bias_lez_125000-100_2.csv");
P3 <- read.csv("../treeppl_bias/tdbd_bias_lez_125000-100_3.csv");
P4 <- read.csv("../treeppl_bias/tdbd_bias_lez_125000-100_4.csv");
P <- list(P1,P2,P3,P4)

# Sample n from tdbd process
res <- data.frame()
for (run in 1:4) {
    X <- P[[run]]
    cat("Samples from run",run,"\n")
    for (iter in 250:1250 ) {

        if ((iter-250)%%100 == 0)
            cat(100*((iter-250)/1000),"% complete\n")

        lambda <- X$lambda[iter+1];
        epsilon <- X$epsilon[iter+1];
        mu <- lambda*epsilon;
        z <- X$z[iter+1];

        otus <- numeric(nrow(D))
        for (i in 1:nrow(D)) {
            time <- transformTime(D$age[i], z)
            otus[i] <- sampleN(time, lambda, mu)
            if (D$Order[i]=="Coleoptera")  {
                n <- rbinom(1,otus[i],rho)
                while (n==0) { n <- rbinom(1,otus[i],rho) } # Condition on sampling at least 1 taxon (as in inference)
                otus[i] <- n
            }
        }
        col <- paste0("sample",iter)
        if (iter==250 && run==1) {
            res <- data.frame(otus)
            colnames(res) <- col
        } else {
            res <- cbind(res,data.frame(otus))
            colnames(res)[ncol(res)] <- col
        }
    }
}
write.table(res,"../data/predictive_samples_bias.tsv",sep="\t",row.names=FALSE)

