#!/usr/bin/env Rscript
library(tidyverse)
library(ggnewscale)
library(ggrepel)
library(data.table)
library(glue)
library(RColorBrewer)
library(scales)
library(cowplot)
library(argparse)
library(karyoploteR)
library(GenomicRanges)

dir <- paste0(getwd(), "/workflow/scripts")
load(glue("{dir}/chm13.karyo.RData"))


# create parser object
indir <- "~/Desktop"
parser <- ArgumentParser()
parser$add_argument("-a",
                    "--asm",
                    help = "bed file with all the asm mapping",
                    default = glue("{indir}/GM08714_1.bed"))
parser$add_argument("-b", "--asm2", help = "bed file with a second asm mapping")#, default = glue("{indir}/GM08714_2.bed"))
parser$add_argument("-k", "--karyotype", help = "karyotpye file for different genomes")
parser$add_argument("--min", help = "minimum amount of total alginemnts between a target and query for it to appear", default = 1e6)
parser$add_argument("-p", "--plot", help = "output plot, must have .pdf ext.", default = "~/Desktop/ideogram.pdf")
args <- parser$parse_args()
filename <- args$asm


asmdf <- function(filename, colors, minalnsize = args$min) {
  asmvshg <- read.table(filename, header = T, comment.char = ">")
  names(asmvshg)[1:3] <- c("chr", "start", "end")
  asmvshg <- asmvshg %>%
    group_by(query_name, chr) %>%
    summarise(
      bp_aligned = sum(end - start),
      min_start = min(start),
      max_end = max(end),
    ) %>%
    filter(bp_aligned > minalnsize) %>%
    merge(asmvshg) %>%
    mutate(group_num = group_indices(., query_name, chr)) %>%
    group_by(query_name, chr) %>%
    mutate(count_in_chr = n(),) %>%
    ungroup() %>%
    arrange(chr, min_start) %>%
    data.table()
  
  asmvshg$name <- asmvshg$query_name
  print(head(asmvshg))
  print(tail(asmvshg))
  curcolor <- 1
  lencolors <- length(colors)
  precontig <- ""
  asmcolor <- NULL
  seen <- c()
  y <- NULL
  for (i in 1:nrow(asmvshg)) {
    contig <- as.character(asmvshg$name[i])
    if (contig != precontig) {
      curcolor <- (curcolor + 1) %% lencolors
      precontig <- contig
    }
    asmcolor <- c(asmcolor, colors[curcolor + 1])
    y <- c(y, curcolor / 4)
  }
  asmvshg$color <- asmcolor
  asmvshg$y <- y
  asmvshg$y1 <- asmvshg$y + .25
  print(head(asmvshg))
  return(asmvshg)
}

asmvshg <- asmdf(filename, c("#2081f9", "#f99820"))

tables = list(asmvshg)
if (!is.null(args$asm2)) {
  asmvshg2 <- asmdf(args$asm2, c("#159934", "#99157a"))
  tables[[2]] = asmvshg2
}
tables


cex <- 0.5

print("Plotting")

pdf(file = args$plot,
    width = 9,
    height = 11)

if (is.null(args$asm2)) {
  kp <-
    plotKaryotype(genome = GENOME,
                  cytobands = CYTOFULL,
                  chromosomes = NOM)
} else {
  kp <-
    plotKaryotype(
      genome = GENOME,
      cytobands = CYTOFULL,
      chromosomes = NOM,
      plot.type = 2
    )
}

for (i in 1:length(tables)) {
  data = tables[[i]]
  data$mid = (data$y + data$y1) / 2
  
  data_u = data %>%
    group_by(chr, query_name, mid, min_start, max_end, color) %>%
    summarise()
  
  # add spanning lines
  kpRect(
    kp,
    chr = data_u$chr,
    x0 = data_u$min_start,
    x1 = data_u$max_end,
    y0 = data_u$mid,
    y1 = data_u$mid,
    col = data_u$color,
    data.panel = i
  )
  # add colored blocks
  kpRect(
    kp,
    chr = data$chr,
    x0 = data$start,
    x1 = data$end,
    y0 = data$y,
    y1 = data$y1,
    col = data$color,
    data.panel = i
  )
}

dev.off()


