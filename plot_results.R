library(ggplot2)
library(grid)

clargs <- commandArgs(TRUE)

print(clargs)
feature_name <- clargs[1]
region_name <- clargs[2]
real_bps_hit <- as.numeric(clargs[3])
real_features_hit <- as.numeric(clargs[4])
real_bases_overlapped <- as.numeric(clargs[5])
directory <- clargs[6]
color="lightblue"
if (length(clargs) == 7) {
   color <- clargs[7]
}

pdf(paste(directory, "/", gsub(" ", "_", region_name), "_to_", gsub(" ", "_", feature_name), ".pdf",sep=""))
bps_with_hits <- read.table(paste(directory, "/", 'bps_with_hits.txt', sep=""))
feature_hits <- read.table(paste(directory, "/", 'feature_hits.txt', sep=""))
bases_overlapped <- read.table(paste(directory, "/", 'bases_overlapped.txt', sep=""))

filename <- paste(paste(directory, "/", "results.txt", sep=""))
unlink(filename)

plot_random_hist <- function(realnum, iterations, title, outfile, metric, region_name, feature_name, color) {
  names(iterations) <- c("Iteration", "Overlaps")
  rcdf <- ecdf(iterations$Overlaps)
  quant <- rcdf(realnum)
  write(c(metric, "REAL", region_name, feature_name, realnum), sep="\t", file=outfile, append=TRUE, ncolumns=5)
  write(c(metric, "QUANTILE", region_name, feature_name, quant), sep="\t", file=outfile, append=TRUE, ncolumns=5)
  real_quant_label <- paste("Obs. Value: ", realnum, "\n", "Quantile: ", quant, sep="")
  p <- ggplot(iterations, aes(x=Overlaps)) + labs(title=title, count="Number of Overlaps") 
  med <- median(iterations$Overlaps)
  print(paste("median=",med))
#  binwidth <- ifelse(med == 0,max(1,(max(iterations$Overlaps) - min(iterations$Overlaps))/30),max(1,med/100))
  binwidth <- max(1,(max(iterations$Overlaps) - min(iterations$Overlaps))/60)
  print(paste("binwidth: ", binwidth))
  p <- p + geom_histogram(binwidth=binwidth,colour="black", fill=color)
  h <- with(iterations, hist(Overlaps, plot=FALSE, breaks=seq(min(Overlaps)-binwidth, max(Overlaps)+binwidth, by=binwidth)))  
  arrowlen <- max(h$counts) / 5
  print(h$breaks)
  print(h$counts)
  print(paste("real: ", realnum))
  if (realnum > max(h$breaks)) {
     predAtReal <- 0
  } else {
    print(paste("break for real: ", min(which(h$breaks > realnum))))
    predAtReal <- max(h$counts[min(which(h$breaks > realnum)) - 1], h$counts[min(which(h$breaks > realnum)) - 2])
  }
  print(paste("count at real", predAtReal))
  print(paste("arrow len", arrowlen))	
  print(paste("arrow y will be ", (arrowlen + predAtReal)))
  p <- p +							geom_segment(x=realnum, y=(arrowlen + predAtReal), xend=realnum, yend=predAtReal, arrow = arrow(length = unit(0.5, "cm"))) + 
       	 							annotate("text", x=realnum, y=arrowlen*1.5 + predAtReal, label="Observed\nOverlaps") +
       	 							xlim(min(min(iterations$Overlaps) - binwidth * 2, realnum - binwidth * 2), max(max(iterations$Overlaps) + binwidth * 2, realnum + binwidth * 2)) +
       	 							theme_bw() +
								theme(plot.background = element_blank()
      							  		           ,panel.grid.major = element_blank()			  
   										   ,panel.grid.minor = element_blank()) 
  plot(p)								
}

plot_random_hist(real_bps_hit, bps_with_hits, paste("Number of ", region_name, "s that Overlap a ", feature_name, sep=""), filename, "REGION_HITS", region_name, feature_name, color)
plot_random_hist(real_features_hit, feature_hits, paste("Number of ", feature_name, "s that Overlap a ", region_name, sep=""), filename, "FEATURE_HITS", region_name, feature_name, color)
plot_random_hist(real_bases_overlapped, bases_overlapped, paste("Portion of ", region_name, "s in ", feature_name, "s (bp)", sep=""), filename, "OVERLAP", region_name, feature_name, color)

shifted_bp_hits <- read.table(paste(directory, "/", 'shifted_bp_hits.txt', sep=""))
shifted_feature_hits <- read.table(paste(directory, "/",'shifted_feature_hits.txt', sep=""))
shifted_bases_overlapped <- read.table(paste(directory, "/",'shifted_bases_overlapped.txt', sep=""))

plot_shift_hits <- function(shift_results, title) {
  names(shift_results) <- c("Shift", "Overlaps")
  p <- ggplot(shift_results, aes(x=Shift, y=Overlaps)) + opts(title=title)
  p + geom_line() +
       	 							theme_bw() +
								theme(plot.background = element_blank()
      							  		           ,panel.grid.major = element_blank()			  
   										   ,panel.grid.minor = element_blank()) 
}
      
plot_shift_hits(shifted_bp_hits, paste("Number of ", region_name, "s that Overlap a", feature_name, sep=""))
plot_shift_hits(shifted_feature_hits, paste("Number of ", feature_name, "s that Overlap a ", region_name, sep=""))
plot_shift_hits(shifted_bases_overlapped, paste("Portion of ", region_name, "s in ", feature_name, "s (bp)", sep=""))
dev.off()
