# precompile static vignette since it takes too long for CRAN servers:
setwd(here::here("vignettes"))

knitr::knit("UsingSpikeSlabGAM.Rnw.orig")
file.rename("UsingSpikeSlabGAM.Rnw.txt", "UsingSpikeSlabGAM.tex")
tinytex::latexmk("UsingSpikeSlabGAM.tex")
file.rename("UsingSpikeSlabGAM.pdf", "UsingSpikeSlabGAM-precompiled.pdf")
knitr::purl("UsingSpikeSlabGAM.Rnw.orig", output = "UsingSpikeSlabGAM.R")

tools::compactPDF("UsingSpikeSlabGAM-precompiled.pdf", qpdf = "qpdf", gs_quality = "ebook")
pdfs <- list.files(pattern = glob2rx("UsingSpikeSlabGAMplot*.pdf"))
file.remove(pdfs)
file.remove("UsingSpikeSlabGAM.log")
file.remove("UsingSpikeSlabGAM.tex")

# unlink("./cache", recursive = TRUE)
