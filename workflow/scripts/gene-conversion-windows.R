source("workflow/scripts/plotutils.R")
f <- "/Users/mrvollger/Desktop/EichlerVolumes/assembly_breaks/nobackups/asm-to-reference-alignment/results/CHM13_V1.1/gene-conversion/all_candidate_windows.tbl.gz"
f <- "results/CHM13_V1.1/gene-conversion/all_candidate_windows.tbl.gz"
window <- 1e4
window <- snakemake@params$window
simplify <- snakemake@params$simplify
merge <- snakemake@params$merge

f <- snakemake@input$bed
print(f)
options(scipen = 999)

df <- fread(f, nThread = 8, sep = "\t") %>%
    mutate(reference_name = `#reference_name`) %>%
    group_by(group, contig, reference_name, reference_name.liftover, sample)

if (merge) {
    df <- df %>%
        ungroup() %>%
        group_by(group, reference_name, reference_name.liftover)
}

if (simplify) {
    print("SIMPLIFY")
    df <- df %>%
        mutate(mmscore = matches - matches.liftover) %>%
        group_by(group, contig, reference_name, reference_name.liftover, sample) %>%
        top_n(1, mmscore) %>%
        dplyr::select(-mmscore)
}

df <- df %>%
    summarise(
        `#reference_name` = unique(`#reference_name`),
        reference_name = unique(reference_name),
        reference_start = min(reference_start),
        reference_end = max(reference_end),
        matches = sum(matches),
        mismatches = sum(mismatches),
        deletion_events = sum(deletion_events),
        insertion_events = sum(insertion_events),
        deletions = sum(deletions),
        insertions = sum(insertions),
        reference_name.liftover = unique(reference_name.liftover),
        reference_start.liftover = min(reference_start.liftover),
        reference_end.liftover = max(reference_end.liftover),
        matches.liftover = sum(matches.liftover),
        mismatches.liftover = sum(mismatches.liftover),
        deletion_events.liftover = sum(deletion_events.liftover),
        insertion_events.liftover = sum(insertion_events.liftover),
        deletions.liftover = sum(deletions.liftover),
        insertions.liftover = sum(insertions.liftover),
        sample = paste(unique(sample), collapse = ";"),
        # sample= unique(sample),
        contig = paste(unique(contig), collapse = ";"),
        contig_start = min(contig_start),
        contig_end = max(contig_end),
        overlap = sum(overlap)
    ) %>%
    mutate(
        perID_by_matches = 100 * matches / (matches + mismatches),
        perID_by_events =
            100 * matches / (matches + mismatches + insertion_events + deletion_events),
        perID_by_all =
            100 * matches / (matches + mismatches + insertions + deletions),
        perID_by_matches.liftover =
            100 * matches.liftover / (matches.liftover + mismatches.liftover),
        perID_by_events.liftover =
            100 * matches.liftover / (matches.liftover + mismatches.liftover + insertion_events.liftover + deletion_events.liftover),
        perID_by_all.liftover =
            100 * matches.liftover / (matches.liftover + mismatches.liftover + insertions.liftover + deletions.liftover),
        original_source =
            paste(contig, ":", contig_start, "-", contig_end, sep = "")
    ) %>%
    data.table()

remove_par <- (df$reference_name == "chrX" &
    df$reference_name.liftover == "chrY") |
    (df$reference_name == "chrY" &
        df$reference_name.liftover == "chrX")
df <- df[!remove_par, ]

df
# reference_name == reference_name.liftover &
gc.df <- df[
    (perID_by_matches >= 99.5 &
        perID_by_all > perID_by_all.liftover &
        (
            (mismatches.liftover - mismatches >= 2 * window / 1e4) |
                (mismatches.liftover / mismatches > 2)
        )
    )
]
dim(gc.df)
print("data subset")
if (F) {
    dim(gc.df)
    p <- df[perID_by_all >= 99.0 &
        perID_by_all - perID_by_all.liftover > 0.00 &
        matches - matches.liftover > 0 &
        matches + mismatches > 9e3 &
        matches.liftover + mismatches.liftover > 9e3 &
        `#reference_name` != "chrY"] %>%
        filter(perID_by_events > 99.5) %>%
        ggplot() +
        geom_histogram(aes(matches - matches.liftover),
            binwidth = 1
        ) +
        facet_zoom(xlim = c(0, 20)) +
        # scale_y_continuous(trans = "log10") +
        # scale_x_continuous(trans = "log10") +
        # annotation_logticks() +
        theme_cowplot()
    nrow(p$data)
    p <- p + annotate("text",
        x = 500, y = 500,
        label = paste("n = ", comma(nrow(p$data))),
        size = 3
    )
    ggsave("~/Desktop/gc.pdf", plot = p)
}

gc.df$name <- paste(
    gc.df$mismatches.liftover - gc.df$mismatches,
    ";",
    gc.df$reference_name,
    ":",
    gc.df$reference_start,
    sep = ""
)
gc.df$name <- gc.df$mismatches.liftover - gc.df$mismatches

gc.df$strand <- "."
gc.df$color <- "0,127,211"
gc.df$score <- pmin(gc.df$mismatches.liftover - gc.df$mismatches, 1000)
gc.df$thickStart <- gc.df$reference_start.liftover
gc.df$thickEnd <- gc.df$reference_end.liftover
gc.df$status <- "Acceptor"

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
    "status",
    "reference_name",
    "reference_start",
    "reference_end",
    "mismatches.liftover",
    "mismatches",
    "perID_by_all.liftover",
    "perID_by_all",
    "sample",
    "original_source"
)]
fwrite(
    odf %>%
        arrange(reference_name.liftover, reference_start.liftover) %>%
        dplyr::rename(`#reference_name.liftover` = reference_name.liftover) %>%
        data.table(),
    file = snakemake@output$acceptor,
    sep = "\t",
    row.names = F,
    quote = F,
    scipen = 999
)

#
# add donor sites
#
odf2 <- data.table(copy(odf))
odf2$color <- "211,144,0"
odf2$status <- "Donor"

odf2$reference_name.liftover <- odf$reference_name
odf2$reference_start.liftover <- odf$reference_start
odf2$reference_end.liftover <- odf$reference_end

odf2$reference_name <- odf$reference_name.liftover
odf2$reference_start <- odf$reference_start.liftover
odf2$reference_end <- odf$reference_end.liftover

odf2$thickStart <- odf2$reference_start.liftover
odf2$thickEnd <- odf2$reference_end.liftover

fwrite(
    rbind(odf, odf2) %>%
        arrange(reference_name.liftover, reference_start.liftover) %>%
        dplyr::rename(`#reference_name.liftover` = reference_name.liftover) %>%
        data.table(),
    file = snakemake@output$bed,
    sep = "\t",
    row.names = F,
    quote = F,
    scipen = 999
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
ndf$score <- pmin(sdf$mismatches.liftover - sdf$mismatches, 1000)
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
    quote = F,
    scipen = 999
)
print("written")