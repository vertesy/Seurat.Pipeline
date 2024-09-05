######################################################################
# Plot Soup.or.Cell.Stats.visualisation.R
######################################################################
# source('~/GitHub/Packages/Seurat.pipeline/misc/Soup.or.Cell.Stats.visualisation.R')
try.dev.off()

# Functions ------------------------
# source('~/Github/TheCorvinas/R/DatabaseLinke.r')
require(tidyverse);  require(cowplot)
require(MarkdownReports);
require(tictoc);
# source('~/GitHub/Packages/CodeAndRoll/CodeAndRoll.R')
require(CodeAndRoll2)
require(Stringendo)
require(ReadWriter)
# source('https://raw.githubusercontent.com/vertesy/CodeAndRoll/master/CodeAndRoll.R')
# require(ggExpress);
require(ggExpress)

# Setup ------------------------------------------------------------------------
OutDir <- OutDirOrig <- "~/Dropbox/Abel.IMBA/AnalysisD/Abel/SEO/SNP.demux.Soup.or.Cell.Stats/HCA/s2/"
setup_MarkdownReports(OutDir = OutDir, scriptname = 'Soup.or.Cell.Stats.visualisation.R')

# Read in ------------------------

# InputProjectDir = "/Volumes/SNP.deconvolution/abel.202102.vcf/SEO.R10015.124719/"
# InputProjectDir = "/Volumes/SNP.deconvolution/abel.202102.vcf/HCA_Brain_ABhashPool_Bp_TK297s1_transcriptome"
InputProjectDir = "/Volumes/SNP.deconvolution/abel.202102.vcf/HCA_Brain_ABhashPool_Bp_TK297s2_transcriptome"
# InputProjectDir = "/Volumes/SNP.deconvolution/abel.202102.vcf/HCA_Brain_ABhashPool_Bp_TK297s3_transcriptome"

path.SubDirs <- list.dirs(path = InputProjectDir)[-1]
for (i in 1:l(path.SubDirs)) {

  # Read in ------------------------
  Subdir.path <- path.SubDirs[i]
  SubDir <- basename(Subdir.path)

  create_set_SubDir(basename(Subdir.path))
  file.AR <- list.files(path = Subdir.path, pattern = "ambient_rna")
  Ambient <-  read.simple.vec(kpps(Subdir.path,file.AR))
  Ambient.RNA <- iround(as.numeric(
    str_split_fixed(
      str_split_fixed(Ambient, pattern = "ambient RNA estimated as ", n=2)[2]
      , pattern = '%', n=2)[1]
    ))
  qbarplot(Ambient.RNA, ylim = c(0,100), label.pos = "out",  label = TRUE, ylab ="% Ambient RNA estimated")

  "Manually created file containing the names of the lines is necessary: lines.txt, one name per line."
  file.lines <- list.files(path = Subdir.path, pattern = "lines.txt") ;
  cell.lines <-  read.simple.vec(kpps(Subdir.path,file.lines))
  nr.lines <- length(cell.lines)
  orig.IDs <- 0:(nr.lines-1)

  file.assignments <- list.files(path = Subdir.path, pattern = "^clusters\\..+tsv|clusters.tsv")
  stopifnot(length(file.assignments)==1)
  df.assignments <- readr::read_tsv(file = (kpps(Subdir.path,file.assignments)))

  # Parse metadata ------------------------
  df.assignments$assignment.named <- translate(vec = df.assignments$assignment
                                               , old = orig.IDs
                                               , new = ppp(orig.IDs, cell.lines))
  # Basic stats ------------------------
  status.of.assignment <- table(df.assignments$status)
  qpie(status.of.assignment)

  {
    cell.line.assignment <- table(df.assignments$assignment.named)
    qbarplot(cell.line.assignment, subtitle = SubDir
             , label = TRUE, label.pos = "out", xlab.angle = 45)
  }

  {
    idx.small.fraction <- (cell.line.assignment/sum(cell.line.assignment) ) < 0.02
    cell.line.assignment.filt <- cell.line.assignment[which(!idx.small.fraction)]
    qpie(cell.line.assignment.filt, NamedSlices = T, subtitle = paste(SubDir, "| Fractions < 2% are removed"))
  }

  {
    df.clean <- df.assignments %>%
      filter(status %in% 'singlet') %>%
      select(c(starts_with("cluster"), starts_with("assignment.named")))

    max1 <- rowMax(df.clean[,1:nr.lines])
    max2 <- apply(df.clean[,1:nr.lines], 1, MaxN, 2)

    # max1.clipped <- clip.values(max1, thr = -500, high = F)
    max1.clipped <- clip.outliers(max1, high = T, probs = 0.05, showhist = T)
    # range(max1.clipped)

    assignment.Q <- tibble(
      "Ratio.2nd.Max" = (-max2 / -max1),
      "log(probability)" = max1.clipped,
        "Line" = df.clean$assignment.named
    )
    assignment.quality.scatterplots <- qscatter(assignment.Q, subtitle= "Y axis is cut at the lowest 5% of the values"
                                                , facet.by = 'Line', ncol = nr.lines
                                                , cols = rgb(0,.75,0,.5), size = 1,) + geom_density_2d()
    qqSave(assignment.quality.scatterplots, w = nr.lines*4, h = 4, ext = 'jpg')
  }
  create_set_Original_OutDir()



}

# Make plots ------------------------------------------------------------------------

# Draw and combine plots------------------------------------------------------------------------
