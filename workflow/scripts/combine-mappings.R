source("workflow/scripts/plotutils.R")
sample <- "HG002_1"
windowf <- "results/CHM13_V1.1/gene-conversion/HG002_1_windows.tbl.gz"
liftoverf <- "results/CHM13_V1.1/gene-conversion/HG002_1_liftover_windows.tbl.gz"
liftoverf <- "/Users/mrvollger/Desktop/EichlerVolumes/assembly_breaks/nobackups/asm-to-reference-alignment/results/CHM13_V1.1/gene-conversion/HG002_1_liftover_windows.tbl.gz"
windowf <- "/Users/mrvollger/Desktop/EichlerVolumes/assembly_breaks/nobackups/asm-to-reference-alignment/results/CHM13_V1.1/gene-conversion/HG002_1_windows.tbl.gz"
windowf <- snakemake@input$window
liftoverf <- snakemake@input$liftover
sample <- snakemake@wildcards$sm
window <- 10000
window <- snakemake@params$window
print(windowf)
print(liftoverf)
print(sample)

window.df <- fread(windowf) %>%
    separate(query_name,
        into = c("original_mapping", "original_source"), sep = "::"
    ) %>%
    separate(original_source,
        into = c("contig", "start", "end"),
        sep = ":|-",
        remove = FALSE
    ) %>%
    mutate(
        contig_start = as.numeric(start) + query_start,
        contig_end = as.numeric(start) + query_end
    ) %>%
    dplyr::select(-query_start, -query_end, -query_length, -start, -end) %>%
    data.table()

liftover.df <- fread(liftoverf) %>%
    mutate(
        original_mapping = paste(
            query_name, ":",
            query_start, "-", query_end,
            sep = ""
        ),
        original_source = paste(
            `#reference_name`, ":",
            reference_start, "-", reference_end,
            sep = ""
        ),
    ) %>%
    dplyr::rename(
        contig = `#reference_name`,
        contig_start = `reference_start`,
        contig_end = `reference_end`,
        reference_name.liftover = query_name,
        reference_start = query_start,
        reference_end = query_end
    ) %>%
    dplyr::select(-query_length) %>%
    data.table()

overlap_bp <- function(df) {
    # x1 <= y2 && y1 <= x2
    intersects <- df$`#reference_name` == df$reference_name.liftover &
        df$reference_start <= df$reference_end.liftover &
        df$reference_start.liftover <= df$reference_end
    overlap <- df$reference_end.liftover - df$reference_start
    intersects * overlap
}

df <- merge(
    window.df,
    liftover.df,
    by = c("original_mapping", "original_source"),
    suffixes = c("", ".liftover")
) %>%
    mutate(overlap = overlap_bp(.)) %>%
    filter(
        overlap == 0 &
            # matches + mismatches >= 0.9 * window &
            # matches.liftover + mismatches.liftover >= 0.9 * window &
            matches + mismatches >= 0.9 * (matches.liftover + mismatches.liftover) &
            matches - matches.liftover > -2 * window / 1e3
    ) %>%
    relocate(original_mapping, .after = last_col()) %>%
    relocate(original_source, .after = last_col()) %>%
    arrange(`#reference_name`, `reference_start`, `reference_end`) %>%
    data.table()
df$sample <- sample

write.table(df, file = snakemake@output$bed, sep = "\t", row.names = F, quote = F)