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