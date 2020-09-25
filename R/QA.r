#collect summary stats for all genes in a hylite output file
# '@ export
hylite_QA <- function(in_dir){
    summs <- list.files(in_dir, pattern = ".read.summary.txt", full.names=TRUE)
    do.call(rbind.data.frame, lapply(summs, .QA_one_file))

}

#'@ export
QA_plot <- function(QA_df){
    ggplot(QA_df, aes(p_reads, fill=assignment)) + 
        geom_histogram(colour='white') + 
        facet_grid(sample ~ assignment) +
        theme_bw()

}

.QA_one_file <- function(fn){
    read_summ <- read_hylite(fn)
    n <- rowSums(read_summ[, c(2:6,8)])

    assign_by_gene <- data.frame(GENE = read_summ[["GENE"]],
                                 assigned = rowSums(read_summ[,2:5])/n, 
                                 uninfromative= read_summ[["UNINFORMATIVE"]]/n,
                                 unknown= read_summ[["UNK"]]/n
    )
    res <- pivot_longer(data=assign_by_gene, -GENE, names_to="assignment", values_to="p_reads")
    res$sample <- str_match(fn, "\\.(.+)\\.read\\.summary\\.txt$")[1,2]
    res
}


