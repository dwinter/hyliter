#'@ export 
read_hylite <- function(file_name){
    read.table(file_name,sep="\t", stringsAsFactors=FALSE, header=1)    
}


#'@ export
parental_DE <- function(expr_data, p1 ,p2){
    spp <- c("GENE", sub("\\..+", "", colnames(expr_data)[-1]))
    parents <- c(p1,p2)
    read_counts <- as.matrix(expr_data[, spp %in% parents])   
    col_data <- data.frame(species=spp[ spp %in% parents ])
    DESeqDataSetFromMatrix( read_counts, 
                           colData = col_data, 
                           design = ~ species)
}



##ASSUMES TWO PARENTS
hybrid_DE <- function(results_dir, species, include_N=FALSE){
    fnames <- list.files(results_dir, 
                         paste0(species, ".*.", "read.summary.txt"), 
                         full.names=TRUE)
    read_counts <- do.call(cbind, lapply(fnames, process_one_h_file, include_N))
    col_data <- data.frame(parent =colnames(read_counts))
    colnames(read_counts) <- paste(colnames(read_counts), rep(1:(ncol(read_counts)/2), each=2), sep="_")
    DESeqDataSetFromMatrix( read_counts, 
                            colData = col_data, 
                            design = ~ parent)
    
    
}


unexpressed <- function(mod, level, assay_n=1, min_reads=2){
    M <- assays(mod)[[assay_n]]
    classifier <- colData(mod)[[level]]
    res <- sapply(levels(classifier), function(x) apply(M[,classifier == x], 1, function(y) all(y <= min_reads)))
    colnames(res) <- levels(classifier)
    res
}


nullified_loci <- function(parent_mod, hybrid_mod, level, assay_n=1, min_reads=2){
    p <- unexpressed(parent_mod, level, assay_n, min_reads)
    h <- unexpressed(hybrid_mod, level, assay_n, min_reads)
    which(!(p[,2]) & h[,2])
}

.classify_expression <- function(x){
    e_from <- ifelse(x[,1], "NULL", "P1")
    paste(e_from, ifelse(x[,2], "NULL", "P2"), sep="_")
}

classify_null_alleles <- function(parent_mod, H1, H2, level=1){
    data.frame( 
        H2 = .classify_expression(unexpressed(H2,level=1)),
        H1 = .classify_expression(unexpressed(H1,level=1)),
        P  = .classify_expression(unexpressed(parent_mod,level=1))
    )    
}

.parse_counts <- function(mod, idx, org_name){
    cdata <- data.frame(counts = assay(mod)[idx,], organism = org_name)
    res <- cbind(do.call(rbind.data.frame, strsplit(rownames(cdata), "_")), cdata)
    names(res)[1:2] <- c("spp", "rep")
    res
}


#classified <- classify_null_alleles(parent_mod, H1_mod, H2_mod)
# agg <- aggregate(P ~ H1, data=classified, FUN=table)



extract_read_counts <- function(H_mod, P_mod, idx){
     h <- .parse_counts(H_mod, idx, "hybrid")
     p <- .parse_counts(P_mod, idx, "parent")
     p$spp <- sapply(strsplit(as.character(p$spp), "\\."), "[[", 1)
     rbind(h,p)     
}

plot_read_counts <- function(H_mod, P_mod, idx){
    res <- extract_read_counts(H1_mod, parent_mod, idx)
    ggplot(res, aes(organism, counts, colour=spp, group=rep)) + geom_jitter(height=0, width=0.04)
}


extract_read_counts_multi_h <- function(H1_mod, H2_mod, P_mod, idx){
     h1 <- .parse_counts(H1_mod, idx, "H1")
     h2 <- .parse_counts(H2_mod, idx, "H2")
     p <- .parse_counts(P_mod, idx, "parent")
     p$spp <- sapply(strsplit(as.character(p$spp), "\\."), "[[", 1)
     rbind(p,h1,h2)
}


process_one_h_file <- function(fname, include_N){
    expr_data <- read_hylite(fname)
    if(include_N){
        return(cbind(rowSums(expr_data[,2:3]), rowSums(expr_data[,4:5])))
    }
   as.matrix(expr_data[,c(2,4)])
}
        
       
fit_and_classify <- function(mod, cutoff=2){
    fit <- results(DESeq(mod))
    parent_class <- cut(fit$log2FoldChange, breaks=c(-Inf, -cutoff,cutoff, Inf), labels=c("P1_dom", "equal", "P2_dom"))
    data.frame(fit, parent_class)
}

fit_multi_hylite <- function(...){
    fits <- lapply(list(...), function(m) results(DESeq(m)))
    fits
}

fc_from_fits <- function(fits){
}



#Assume two parents
hybrid_specific_SNVs <- function(snp_file){
    SNVs <- read_hylite(snp_file)
    child_alt <- sapply(strsplit(SNVs[,5],","), function(x) "1" %in% x)
    parents_ref <- apply(SNVs[,6:7], 1, function(x) all(x==0))
    SNVs[child_alt & parents_ref,]
}


