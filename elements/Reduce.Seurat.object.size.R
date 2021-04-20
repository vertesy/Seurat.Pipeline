######################################################################
# Reduce.Seurat.object.size.R
######################################################################
# source('~/GitHub/Packages/Seurat.pipeline/elements/Reduce.Seurat.object.size.R')
# source("https://raw.githubusercontent.com/vertesy/Seurat.Pipeline/main/elements/Reduce.Seurat.object.size.R")
# rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)

# Functions ------------------------
memory.biggest.objects()

# Setup ------------------------
OutDir = "~/Dropbox/Abel.IMBA/AnalysisD/"
setup_MarkdownReports(OutDir = OutDir, scriptname = "")
OutDirOrig = OutDir

# Metadata ------------------------
xdata <- combined.obj@assays$RNA@scale.data
xdata2 <- round(xdata)

ls.mem <- ls( envir = .GlobalEnv)
ls.obj <- lapply(ls.mem, get)
Sizes.of.objects.in.mem <- unlapply(ls.obj, object.size)
names(Sizes.of.objects.in.mem) <- ls.mem
(topX = sort(Sizes.of.objects.in.mem,decreasing = TRUE)[1:5])

# Parameters ------------------------
xdata <- combined.obj@assays$RNA@scale.data
xdata2 <- round(xdata)

# How did it become so big? Memory size check ------------------------

inFile = "~/Dropbox/Abel.IMBA/AnalysisD/_RDS.files/SEO.rds/re_m1500c21.dSample.Organoids_2000/combined.obj__2020.12.11_15.51.Rds.gz"
# inFile = "~/Documents/RDS.files/EDA.Literature/combined.obj__2021.01.05_14.56.Rds.gz" # 20GB
combined.obj <- read_rds(inFile)
combined.obj <- RecallReduction(obj = combined.obj, dim = 2)
recall.all.genes()


isave.RDS(object = combined.obj, inOutDir = F)
# Add ScaleData
tic(); combined.obj <- ScaleData(combined.obj, assay = "RNA", verbose = T, vars.to.regress = p$'variables.2.regress'); toc()
isave.RDS(object = combined.obj, inOutDir = F)

```r
# p$'variables.2.regress'  = NULL

> print(object.size(combined.obj), units = "auto")
7.8 Gb

tic(); combined.obj <- ScaleData(combined.obj, assay = "RNA", verbose = T, vars.to.regress = p$'variables.2.regress'); toc()
> print(object.size(combined.obj), units = "auto")
24.1 Gb

tic(); combined.obj <- ScaleData(combined.obj, assay = "integrated", verbose = T, vars.to.regress = p$'variables.2.regress'); toc()
> print(object.size(combined.obj), units = "auto")
24.1 Gb
```

- It is also not clear to me why adding ScaleData/integrated does not change the size at all. (Its only 2000 genes, but still)


```r
rounded.scdat <- round(combined.obj@assays$RNA@scale.data, digits = 2)
view.head(rounded.scdat)
combined.obj@assays$RNA@scale.data <- rounded.scdat
print(object.size(combined.obj), units = "auto")
24.1 Gb
```

- I am surprised that eliminating most digits did not decrease the object size. This rounding to save size works well when saving matrices into .csv files.

```r
# BEFORE
> view.head(combined.obj@assays$RNA@scale.data)
CTGATCCCACAAGCCC-1_1 AGAATAGAGTAAGTAC-1_1
AL627309.1          -0.10839875          -0.10839875
AL669831.5          -0.35907276          -0.35907276
LINC00115           -0.09303572          -0.09303572
FAM41C              -0.02871722          -0.02871722
AL645608.7          -0.16244730          -0.16244730
AL645608.1          -0.04892921          -0.04892921
SAMD11              -0.13318995          -0.13318995
NOC2L               -0.44912416          -0.44912416
KLHL17              -0.13462563          -0.13462563
AL645608.8          -0.08671774          -0.08671774

# AFTER
> view.head(round(combined.obj@assays$RNA@scale.data, digits = 2))
CTGATCCCACAAGCCC-1_1 AGAATAGAGTAAGTAC-1_1
AL627309.1                -0.11                -0.11
AL669831.5                -0.36                -0.36
LINC00115                 -0.09                -0.09
FAM41C                    -0.03                -0.03
AL645608.7                -0.16                -0.16
AL645608.1                -0.05                -0.05
SAMD11                    -0.13                -0.13
NOC2L                     -0.45                -0.45
KLHL17                    -0.13                -0.13
AL645608.8                -0.09                -0.09
```

```r
tbx <- table(rounded.scdat[1:5000, 1:5000])
head(sort(tbx, decreasing = T))
sorted.values.scaledata <- sort(tbx, decreasing = T)
qpie(sorted.values.scaledata)
```


```r
tbx <- table(round(rounded.scdat[1:5000, 1:5000], digits = 1))
head(sort(tbx, decreasing = T))
ratio <- tbx/ sum(tbx)
head(sort(ratio, decreasing = T))
sorted.values.scaledata <- sort(tbx, decreasing = T)
qpie(sorted.values.scaledata)
```


# QC ------------------------
if (F) {
  workers = 2
  plan("multiprocess", workers = workers)
  options(future.globals.maxSize = 10 * 1024^3) # 10 * 1024 ^ 3 byte = 10 GiB
}

# ------------------------
# ------------------------
# ------------------------
# ------------------------
# ------------------------



