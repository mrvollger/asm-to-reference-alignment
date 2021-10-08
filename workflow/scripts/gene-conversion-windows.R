source("workflow/scripts/plotutils.R")
f <- "/Users/mrvollger/Desktop/EichlerVolumes/assembly_breaks/nobackups/asm-to-reference-alignment/all_candidate_windows.gz"
f <- snakemake@input$tbl
print(f)
df <- fread(f, nThread = 8, sep = "\t") %>% data.table()

gc.df <- df[
    (overlap == 0 &
        perID_by_all - perID_by_all.liftover > 0.1 &
        mismatches.liftover - mismatches > 1 &
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

odf <- gc.df[, c(
    "reference_name.liftover",
    "reference_start.liftover",
    "reference_end.liftover",
    "perID_by_all.liftover",
    "mismatches.liftover",
    "#reference_name",
    "reference_start",
    "reference_end",
    "perID_by_all",
    "mismatches"
)]
write.table(odf,
    file = snakemake@output$tbl,
    sep = "\t", row.names = F, quote = F
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

odf$`#chrom` <- sdf$reference_name.liftover
odf$chromStart <- sdf$reference_start.liftover
odf[chromStart > reference_start]$chromStart <-
    sdf[reference_start.liftover > reference_start]$reference_start

odf$chromEnd <- sdf$reference_end.liftover
odf[chromEnd < reference_end]$chromEnd <-
    sdf[reference_end.liftover < reference_end]$reference_end

odf$name <- "."
odf$score <- sdf$mismatches.liftover - sdf$mismatches
odf$value <- sdf$mismatches.liftover / sdf$mismatches
odf$exp <- "."
odf$color <- 0

odf$sourceChrom <- sdf$reference_name.liftover
odf$sourceStart <- sdf$reference_start.liftover
odf$sourceEnd <- sdf$reference_end.liftover
odf$sourceName <- "."
odf$sourceStrand <- "."

odf$targetChrom <- sdf$reference_name
odf$targetStart <- sdf$reference_start
odf$targetEnd <- sdf$reference_end
odf$targetName <- "."
odf$targetStrand <- "."

# fix the columns when interchromosomal
inter <- sdf$reference_name != sdf$reference_name.liftover
# TODO fix this so inters can be shown
if (F) {
    odf[inter]$`#chrom` <- sdf[inter]$reference_name.liftover
    odf[inter]$chromStart <- sdf[inter]$reference_start.liftover
    odf[inter]$chromEnd <- sdf[inter]$reference_end.liftover

    odf[inter]$sourceChrom <- sdf[inter]$reference_name
    odf[inter]$sourceStart <- sdf[inter]$reference_start
    odf[inter]$sourceEnd <- sdf[inter]$reference_end

    odf[inter]$targetChrom <- sdf[inter]$`#chrom`
    odf[inter]$targetStart <- sdf[inter]$chromStart
    odf[inter]$targetEnd <- sdf[inter]$chromEnd
} else {
    odf[!inter]
}
write.table(odf[, ..names],
    file = snakemake@output$interact,
    sep = "\t", row.names = F, quote = F
)