#!/bin/env Rscript
#    #!/usr/bin/Rscript --vanilla

args <- commandArgs(TRUE)
if(length(args)!=1) {
  cat("usage: <feature-vectors.txt>\n");
  q(status=1)
}
infile <- args[1];

data <- read.table(infile,header=T)
y <- data$category
x <- data[2:4097]
X <- as.matrix(x)
sink("coefficients-cv-binomial.txt")
cvfit = cv.glmnet(X, y, family="binomial")
coef(cvfit, s = "lambda.min")
sink()


