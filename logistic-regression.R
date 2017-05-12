#!/bin/env Rscript
#    #!/usr/bin/Rscript --vanilla

args <- commandArgs(TRUE)
if(length(args)!=3) {
  cat("usage: <feature-vectors.txt> <alpha:0=ridge,1=lasso> <out.betas>\n");
  q(status=1)
}
infile <- args[1]
ALPHA <- as.numeric(args[2])
outfile <- args[3]

library(glmnet)
library(methods)

data <- read.table(infile,header=T)
y <- data$category
x <- data[2:ncol(data)]
X <- as.matrix(x)
sink(outfile) # "coefficients-cv-binomial.txt"
cvfit = cv.glmnet(X, y, family="binomial", alpha=ALPHA)
coef(cvfit, s = "lambda.min")
sink()


