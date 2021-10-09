source("workflow/scripts/plotutils.R")
f <- "/Users/mrvollger/Desktop/EichlerVolumes/assembly_breaks/nobackups/asm-to-reference-alignment/all_candidate_windows.gz"
f <- "results/CHM13_V1.1/gene-conversion/all_candidate_windows.tbl.gz"
f <- snakemake@input$tbl
print(f)
df <- fread(f, nThread = 8, sep = "\t") %>% data.table()

gc.df <- df[
    (overlap == 0 &
        perID_by_all - perID_by_all.liftover > 0.1 &
        mismatches.liftover - mismatches > 4 &
        matches + mismatches > 9000 &
        matches.liftover + mismatches.liftover > 9000),
]
print("data subset")
if (F) {
    ggplot(gc.df) +
        geom_histogram(aes(mismatches.liftover - mismatches), bins = 20) +
        scale_x_continuous(trans = "log10") +
        annotation_logticks() +
        theme_cowplot()
}

gc.df$reference_name <- gc.df$`#reference_name`
gc.df$name <- paste(
    gc.df$mismatches.liftover - gc.df$mismatches,
    ";",
    gc.df$reference_name,
    ":",
    gc.df$reference_start,
    sep = ""
)
gc.df$strand <- "."
gc.df$color <- "0,127,211"
gc.df$score <- 0
gc.df$thickStart <- gc.df$reference_start.liftover
gc.df$thickEnd <- gc.df$reference_end.liftover

odf <- gc.df[, c(
    "reference_name.liftover",
    "reference_start.liftover",
    "reference_end.liftover",
    "name",
    "score",
    "strand",
    "thickStart",
    "thickEnd",
    "color",
    "reference_name",
    "reference_start",
    "reference_end",
    "mismatches.liftover",
    "mismatches",
    "perID_by_all.liftover",
    "perID_by_all",
    "original_source"
)]

odf2 <- data.table(copy(odf))
odf2$color <- "211,144,0"

odf2$reference_name.liftover <- odf$reference_name
odf2$reference_start.liftover <- odf$reference_start
odf2$reference_end.liftover <- odf$reference_end

odf2$reference_name <- odf$reference_name.liftover
odf2$reference_start <- odf$reference_start.liftover
odf2$reference_end <- odf$reference_end.liftover

odf2$thickStart <- odf2$reference_start.liftover
odf2$thickEnd <- odf2$reference_end.liftover

fwrite(
    rbind(odf, odf2),
    file = snakemake@output$tbl,
    sep = "\t",
    row.names = F,
    quote = F
)
# make the interaction file
names <- c(
    "#chrom", "chromStart", "chromEnd",
    "name",
    "score",
    "value",
    "exp",
    "color",
    "sourceChrom", "sourceStart", "sourceEnd", "sourceName", "sourceStrand",
    "targetChrom", "targetStart", "targetEnd", "targetName", "targetStrand"
)
sdf <- copy(odf)
ndf <- data.table()

ndf$`#chrom` <- sdf$`reference_name.liftover`
ndf$chromStart <- sdf$`reference_start.liftover`
ndf[sdf$`reference_start.liftover` > sdf$reference_start]$chromStart <-
    copy(sdf[sdf$`reference_start.liftover` > sdf$reference_start]$reference_start)

ndf$chromEnd <- sdf$`reference_end.liftover`
ndf[sdf$`reference_end.liftover` < sdf$reference_end]$chromEnd <-
    copy(sdf[sdf$`reference_end.liftover` < sdf$reference_end]$reference_end)

ndf$name <- "."
ndf$score <- 500 # sdf$`mismatches.liftover` - sdf$mismatches
ndf$value <- sdf$mismatches.liftover / sdf$mismatches
ndf$exp <- "."
ndf$color <- 0

ndf$sourceChrom <- sdf$reference_name.liftover
ndf$sourceStart <- sdf$reference_start.liftover
ndf$sourceEnd <- sdf$reference_end.liftover
ndf$sourceName <- "."
ndf$sourceStrand <- "."

ndf$targetChrom <- sdf$reference_name
ndf$targetStart <- sdf$reference_start
ndf$targetEnd <- sdf$reference_end
ndf$targetName <- "."
ndf$targetStrand <- "."

# fix the columns when interchromosomal
inter <- as.character(sdf$reference_name) != as.character(sdf$reference_name.liftover)
sum(inter)
dim(odf)
dim(sdf)
dim(ndf)
print(head(ndf))
if (T) {
    print("inter: setting 1-3")
    ndf[inter]$`#chrom` <- sdf[inter]$reference_name.liftover
    ndf[inter]$chromStart <- sdf[inter]$reference_start.liftover
    ndf[inter]$chromEnd <- sdf[inter]$reference_end.liftover

    print("inter: setting source")
    ndf[inter]$sourceChrom <- sdf[inter]$reference_name
    ndf[inter]$sourceStart <- sdf[inter]$reference_start
    ndf[inter]$sourceEnd <- sdf[inter]$reference_end

    print("inter: setting target name")
    ndf[inter]$targetChrom <- sdf[inter]$reference_name.liftover
    print("inter: setting target start")
    ndf[inter]$targetStart <- sdf[inter]$reference_start.liftover
    print("inter: setting target end")
    ndf[inter]$targetEnd <- sdf[inter]$reference_end.liftover
} else {
    ndf <- ndf[!inter]
}
print(head(ndf))
print(dim(ndf))
fwrite(
    ndf,
    file = snakemake@output$interact,
    sep = "\t",
    row.names = F,
    quote = F
)
print("written")