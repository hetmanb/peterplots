#!/usr/bin/env Rscript
# List of packages for session
usePackage<-function(p){
  # load a package if installed, else load after installation.
  # Args:
  #   p: package name in quotes
  if (!is.element(p, installed.packages()[,1])){
    print(paste('Package:',p,'Not found, Installing Now...'))
    install.packages(p, dep = TRUE, repos = "http://cran.us.r-project.org")}
  print(paste('Loading Package :',p))
  require(p, character.only = TRUE)  
}
usePackage("ggplot2")
usePackage("seqinr")

# Load arguments list (fasta files) from the command line
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

#Create an output folder, and copy fastas with updated contig names into it
files <- args
system("mkdir output")
command <- paste0("awk", " '", "/>/{sub(\">\",\"&\"FILENAME\"_\");sub(/\\.fasta/,x)}1", "' ")
for(n in files){
  cmd <- paste0(command , n, " > ", "./output/renamed_", n)
  system(cmd)
}

#Create a concatenome to work from to build the contig information table
system("cat ./output/renamed_*.fasta > ./output/cat_genomes.fasta")  
d <- read.fasta("./output/cat_genomes.fasta", seqtype = "DNA")

# For each file, creates a data frame called gc_table and populates it with
# length, coverage, gc content and genome for each contig in the fasta file
gc_table <- data.frame(attr(d, "name"))
colnames(gc_table) <- "name"
for(i in 1:length(d)){
  gc_table$GC[i] <- GC(d[[i]])
  gc_table$length[i] <- length(d[[i]])
  split_name <- strsplit(as.character(gc_table[i,1]), split = "_")
  gc_table$coverage[i] <- as.numeric(split_name[[1]][(match("cov", split_name[[1]]) + 1)])
  gc_table$genome[i] <- substr(gc_table[i,1], 1, (regexpr("NODE", gc_table[i,1])[1][1])-2)
}
write.csv(gc_table, "./output/gc_table.csv", row.names = F)
system("rm ./output/cat_genomes.fasta")

#plot the results via ggplot
plot.new()
pp <- ggplot(gc_table, aes(length,coverage, colour = GC)) + geom_point()
pp <- pp + scale_y_log10()
pp <- pp + scale_x_log10()
pp <- pp + facet_wrap(~ genome, ncol = 4)
pp <- pp + scale_color_gradientn(colours = c("red", "yellow", "green"))
pp <- pp + geom_vline("500bp", xintercept = 500) + annotate(geom = "text", x = 400, y = 1000, label = "500 bp", size = 3.5, angle = 90)

#adjust the sizing of the output table to suit the number of genomes being analyzed
l <- length(unique(gc_table$genome))
roundUp <- function(x,to=4)
{
  to*(x%/%to + as.logical(x%%to))
}

h <- roundUp(l, 4)

w <- ifelse(l == 1, 4,
            ifelse(l == 2, 8, 
                   ifelse(l == 3, 12, 
                          ifelse(l >= 4, 16, 16))))

pdf("./output/peter_plot.pdf", width = w, height = h)
plot(pp)
dev.off()
system("rm Rplots.pdf")





