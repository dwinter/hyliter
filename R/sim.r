# functions
# simulate ortholog-homeolog states

# classify genes to ortholog-homeolog classes
classify <- function(vector){
	
	ortholog <- vector[1]
	homeolog <- vector[2]
	
	if     ( ortholog ==  0 && homeolog ==  0 ){ class = 'I' }
	else if( ortholog == +1 && homeolog == +1 ){ class = 'I' }
	else if( ortholog == -1 && homeolog == -1 ){ class = 'I' }
	else if( ortholog == +1 && homeolog ==  0 ){ class = 'L' }
	else if( ortholog == -1 && homeolog ==  0 ){ class = 'L' }
	else if( ortholog ==  0 && homeolog == +1 ){ class = 'S' }
	else if( ortholog ==  0 && homeolog == -1 ){ class = 'S' }
	else if( ortholog == +1 && homeolog == -1 ){ class = 'R' }
	else if( ortholog == -1 && homeolog == +1 ){ class = 'R' }
	else { class = 'U' }
	
	return(class)
}

# count numbers of genes in ortholog-homeolog classes
count <- function(vector){
	
	class.I <- length( which(vector == 'I') )
	class.L <- length( which(vector == 'L') )
	class.S <- length( which(vector == 'S') )
	class.R <- length( which(vector == 'R') )
	class.U <- length( which(vector == 'U') )
	
	return( c(class.I, class.L, class.S, class.R, class.U) )
}

# make random string for filenames
randomstring <- function(n=1, length=12){

    randomString <- c(1:n)
    for (i in 1:n){
        
        randomString[i] <- paste(sample(c(0:9, letters, LETTERS),
                                 length, replace=T),
                                 collapse="")
    }
    return(randomString)
}

# wrapper script for running simulations and analyses
runsims <- function(iterations, n.genes, scalar, bio.reps, significance, fold.change, o.n.diffexp, o.seqdepth, o.fraction.upregulated, h.n.diffexp, h.seqdepth, h.fraction.upregulated, write_files=FALSE){
	
	# make data frame
	this.res <- as.data.frame(matrix(0, ncol=26, nrow=0))
	colnames(this.res) <- c(
		"ngenes", "nreps", "o.nreads", "o.de", "o.frac.up", "h.nreads", "h.de", "h.frac.up",
		"o.up", "o.down", "h.up", "h.down", "o.up.fc", "o.down.fc", "h.up.fc", "h.down.fc",
		"I", "L", "S", "R", "U", "I.fc", "L.fc", "S.fc", "R.fc", "U.fc")
	
	# record simulation parameters
	conditions <- c(n.genes, bio.reps, o.seqdepth, o.n.diffexp, o.fraction.upregulated, h.seqdepth, h.n.diffexp, h.fraction.upregulated)
	
	# simulate orthologs
    filename1 <- if(write_files) paste0(randomstring(), ".rds") else NULL    
	o.sim <- generateSyntheticData(	dataset = "o.sim", 
									n.vars = n.genes * scalar,
									samples.per.cond = bio.reps,
									n.diffexp = o.n.diffexp * scalar,
									seqdepth = o.seqdepth * scalar,
									fraction.upregulated = o.fraction.upregulated,
									output.file = filename1)
	
	runDiffExp(data.file = filename1, result.extent = "DESeq2",
	               Rmdfunction = "DESeq2.createRmd",
	               output.directory = ".", fit.type = "parametric",
	               test = "Wald")
	
	filename2   <- paste0(fileprefix1, "_DESeq2.rds")
	o.sim.dat   <- readRDS(filename2)
	
	# simulate homeologs
	fileprefix3 <- randomstring()
	filename3   <- paste0(fileprefix3, ".rds")
	h.sim <- generateSyntheticData(	dataset = "h.sim", 
									n.vars = n.genes * scalar,
									samples.per.cond = bio.reps,
									n.diffexp = h.n.diffexp * scalar,
									seqdepth = h.seqdepth * scalar,
									fraction.upregulated = h.fraction.upregulated,
									output.file = filename3)
	
	runDiffExp(data.file = filename3, result.extent = "DESeq2",
	               Rmdfunction = "DESeq2.createRmd",
	               output.directory = ".", fit.type = "parametric",
	               test = "Wald")
	
	filename4   <- paste0(fileprefix3, "_DESeq2.rds")
	h.sim.dat <- readRDS(filename4)
	
	# remove files
	if (file.exists(filename1)) file.remove(filename1)
	if (file.exists(filename2)) file.remove(filename2)
	if (file.exists(filename3)) file.remove(filename3)
	if (file.exists(filename4)) file.remove(filename4)
	
	
	for( i in 1:iterations ){
		
		# make subsampled dataset (note: datasets can differ in size by 1)
		o.values <- sample(1:length(o.sim.dat@result.table[,1]), n.genes, replace=F)
		h.values <- sample(1:length(h.sim.dat@result.table[,1]), n.genes, replace=F)
		
		this.o.sim.dat <- o.sim.dat@result.table[o.values,]
		this.h.sim.dat <- h.sim.dat@result.table[h.values,]
		
		# determine differential expression (statistically significant)
		o.sign.up   <- which( this.o.sim.dat$adjpvalue <= significance & this.o.sim.dat$logFC >= 0 )
		o.sign.down <- which( this.o.sim.dat$adjpvalue <= significance & this.o.sim.dat$logFC < 0 )
		o.nonsign   <- setdiff(1:n.genes, c(o.sign.up, o.sign.down))
		
		h.sign.up   <- which( this.h.sim.dat$adjpvalue <= significance & this.h.sim.dat$logFC >= 0 )
		h.sign.down <- which( this.h.sim.dat$adjpvalue <= significance & this.h.sim.dat$logFC < 0 )
		h.nonsign   <- setdiff(1:n.genes, c(h.sign.up, h.sign.down))
	
		sim.nde <- c(length(o.sign.up), length(o.sign.down), length(h.sign.up), length(h.sign.down))
	
		# determine differential expression (statistically significant and >=x-fold change)	
		o.sign.up.fc   <- which( this.o.sim.dat$adjpvalue <= significance & this.o.sim.dat$logFC >= log2(fold.change) )
		o.sign.down.fc <- which( this.o.sim.dat$adjpvalue <= significance & this.o.sim.dat$logFC <= -log2(fold.change) )
		o.nonsign.fc   <- setdiff(1:n.genes, c(o.sign.up.fc, o.sign.down.fc))
		
		h.sign.up.fc   <- which( this.h.sim.dat$adjpvalue <= significance & this.h.sim.dat$logFC >= log2(fold.change) )
		h.sign.down.fc <- which( this.h.sim.dat$adjpvalue <= significance & this.h.sim.dat$logFC <= -log2(fold.change) )
		h.nonsign.fc   <- setdiff(1:n.genes, c(h.sign.up.fc, h.sign.down.fc))
		
		sim.nde.fc <- c(length(o.sign.up.fc), length(o.sign.down.fc), length(h.sign.up.fc), length(h.sign.down.fc))
		
		# build change vectors (statistically significant)
		orthologs <- vector(length=n.genes)
		orthologs[o.sign.up]   <- +1
		orthologs[o.sign.down] <- -1
		orthologs[o.nonsign]  <- 0
		
		homeologs <- vector(length=n.genes)
		homeologs[h.sign.up]   <- +1
		homeologs[h.sign.down] <- -1
		homeologs[h.nonsign]  <- 0
		
		simulation <- matrix(data=cbind(orthologs, homeologs), ncol=2)
	
		# build change vectors (statistically significant and >=x-fold change)	
		orthologs.fc <- vector(length=n.genes)
		orthologs.fc[o.sign.up.fc]   <- +1
		orthologs.fc[o.sign.down.fc] <- -1
		orthologs.fc[o.nonsign.fc]  <- 0
		
		homeologs.fc <- vector(length=n.genes)
		homeologs.fc[h.sign.up.fc]   <- +1
		homeologs.fc[h.sign.down.fc] <- -1
		homeologs.fc[h.nonsign.fc]  <- 0
		
		simulation.fc <- matrix(data=cbind(orthologs.fc, homeologs.fc), ncol=2)
		
		# calculate ortholog-homeolog expression classes
		outcome    <- apply(simulation, 1, classify)
		outcome.fc <- apply(simulation.fc, 1, classify)
		
		# summarize results
		this.res[i,] <- c( conditions, sim.nde, sim.nde.fc, count(outcome), count(outcome.fc) )	
	}
	
	return( this.res )
}
