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
usePackage('plotly')
usePackage('htmlwidgets')

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
gc_table$Annot <- NA
gc_table$GC <- NA
gc_table$Length <- NA
gc_table$Coverage <- NA
gc_table$Genome <- NA



colnames(gc_table) <- c("Name", "Annot", "GC", "Length", "Coverage", "Genome")
for(i in 1:length(d)){
  gc_table$Annot[i] <- attr(d[[i]], "Annot")
  gc_table$GC[i] <- GC(d[[i]])
  gc_table$Length[i] <- length(d[[i]])
  split_name <- strsplit(as.character(gc_table[i,2]), split = "[_ ]+") 
  gc_table$Coverage[i] <- as.numeric(split_name[[1]][(match("cov", split_name[[1]]) + 1)])
  gc_table$Genome[i] <- strsplit(as.character(gc_table[i,1]), '_')[[1]][1]
}
write.csv(gc_table, "./output/gc_table.csv", row.names = F)
system("rm ./output/cat_genomes.fasta")

#plot the results via ggplot
plot.new()
pp <- ggplot(gc_table, aes(Length,Coverage, colour = GC, text = Name)) + 
  geom_point() +
  scale_y_log10() +
  scale_x_log10() +
  scale_color_gradientn(colours = c("red", "yellow", "green")) +
  geom_vline("500bp", xintercept = 500) + annotate(geom = "text", x = 400, y = 1000, label = "500 bp", size = 3.5, angle = 90) +
  facet_wrap(~ Genome, ncol = 4) 

#adjust the sizing of the output table to suit the number of genomes being analyzed
l <- length(unique(gc_table$Genome))
roundUp <- function(x,to=4)
{
  to*(x%/%to + as.logical(x%%to))
}

h <- roundUp(l, 4)

w <- ifelse(l == 1, 4,
            ifelse(l == 2, 8, 
                   ifelse(l == 3, 12, 
                          ifelse(l >= 4, 16, 16))))

h2<-h*100
w2<-w*100

plotly_pp <- ggplotly(pp, width = w2, height = h2, tooltip = c("Name","Length", "Coverage", "GC"))
htmlwidgets::saveWidget(plotly_pp, 'peter_plot.html')

pdf(file = "peter_plot.pdf", width = w, height = h)
plot(pp)
dev.off()
system('rm Rplots.pdf')