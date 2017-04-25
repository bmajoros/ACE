#!/bin/env Rscript
#    #!/usr/bin/Rscript --vanilla
library(glmnet)
library(methods)

args <- commandArgs(TRUE)
if(length(args)!=2) {
  cat("usage: <feature-vectors.txt> <out-betas.txt>\n");
  q(status=1)
}
infile <- args[1]
outfile <- args[2]

data <- read.table(infile,header=T)
y <- data$category
x <- data[2:ncol(data)]
X <- as.matrix(x)
sink(outfile) # "coefficients-cv-binomial.txt"
cvfit = cv.glmnet(X, y, family="binomial")
coef(cvfit, s = "lambda.min")
sink()


