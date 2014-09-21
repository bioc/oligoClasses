## Adding toy examples
## Saving the script that creates the examples here
##   so we can improve later

NS <- 5
FS <- 9
rp <- 4
ids <- as.character(t(outer(1:FS, 1:rp, paste, sep=".")))
sns <- paste("sample", 1:NS, sep="")
fns <- paste("feature", ids, sep="")
set.seed(1)
tmpExprs <- matrix(as.integer(2^rnorm(NS*FS*rp, mean=12, sd=2)), nc=NS)
rownames(tmpExprs) <- fns
colnames(tmpExprs) <- sns

## ExpressionFeatureSet
efsExample <- new("ExpressionFeatureSet",
                  exprs=tmpExprs,
                  annotation="example")
save(efsExample, file="../../data/efsExample.rda")

## SnpFeatureSet
ids <- as.character(t(outer(1:FS, c("1-A", "2-A", "1-B", "2-B"), paste, sep=".")))
fns <- paste("SNP", ids, sep="")
rownames(tmpExprs) <- fns
sfsExample <- new("SnpFeatureSet",
                  exprs=tmpExprs,
                  annotation="example")
save(sfsExample, file="../../data/sfsExample.rda")

## SnpCnvFeatureSet

## TilingFeatureSet

## ExonFeatureSet

## GeneFeatureSet


## SnpQSet
set.seed(1)
NS <- 5
FS <- 9
ata <- matrix(rnorm(NS*FS, mean=12, sd=2), nc=NS)
atb <- matrix(rnorm(NS*FS, mean=12, sd=2), nc=NS)
sta <- matrix(rnorm(NS*FS, mean=12, sd=2), nc=NS)
stb <- matrix(rnorm(NS*FS, mean=12, sd=2), nc=NS)
rownames(ata) <- rownames(atb) <- rownames(sta) <- rownames(stb) <- paste("SNP", 1:FS, sep="")
colnames(ata) <- colnames(atb) <- colnames(sta) <- colnames(stb) <- paste("sample", 1:NS, sep="")
sqsExample <- new("SnpQSet",
                  antisenseThetaA=ata,
                  antisenseThetaB=atb,
                  senseThetaA=sta,
                  senseThetaB=stb,
                  annotation="example")
save(sqsExample, file="../../data/sqsExample.rda")

## SnpCnvQSet
scqsExample <- new("SnpCnvQSet",
                   thetaA=ata,
                   thetaB=atb,
                   annotation="example")
save(scqsExample, file="../../data/scqsExample.rda")


