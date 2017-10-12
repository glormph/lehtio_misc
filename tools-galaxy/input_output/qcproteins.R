library(ggplot2)
library(reshape2)
library(VennDiagram)

# testing tables
#args = c('multiplot.R', 'test-data/rebuilt_peptides.txt', 'peptides', 'true')
#args = c('test-data/rebuilt_proteins.txt', 'proteins', 'true')

args = commandArgs(trailingOnly = T)
isobaric_labeled = args[4] == 'true'
featname = substring(args[3], 1, nchar(args[3])-1) # protein, peptide, gene, etc
proteins = read.table(args[2], header=T, sep='\t')
multiplot_fn = args[1]
setnames = args[5:length(args)]
source(multiplot_fn)

# VennDiagram, max Five sets
wideproteins = dcast(proteins, Accession~Set, sum)
if (length(wideproteins) < 7) {
  setlist = list()
  for (setname in names(wideproteins[2:length(wideproteins)])) {
    setprots = wideproteins$Accession[complete.cases(wideproteins[,c(setname, 'Accession')])]
    setlist[[setname]] = setprots
  }
  venn.diagram(x=setlist, filename=sprintf("venn_%ss.png", featname), imagetype='png')
}


if (isobaric_labeled) {
  pagelength = length(setnames) * 1.5
} else {
  pagelength = length(setnames) / 8
}
print(sprintf('Pagelength %s', pagelength))
pdf(sprintf('%ss.pdf', featname), height=pagelength)

if (featname == 'protein') {
    setmedian = median(subset(proteins, Set==setnames[1])$Coverage)
    ggplot(subset(proteins, Set==setnames[1]), aes(x=Coverage)) +
      geom_histogram(position='identity', alpha=0.4, bins=50) +
      ylab(sprintf('%s count', featname)) + 
      annotation_custom(grob=textGrob(label=sprintf('Median: %.3f', setmedian)))
}

ticks = element_line(size=10/length(setnames))
if (length(setnames) > 50) {
  font = element_text(size=1000/length(setnames))
} else {
  font = element_text(size=100/length(setnames))
}
ggplot(proteins, aes(Set, MS1.precursor.area)) +
  geom_boxplot() + # geom_jitter(alpha=0.2) + 
  scale_y_log10() + ylab(sprintf('Precursor area per %s', featname)) + coord_flip() +
  theme(text=font, axis.ticks=ticks, panel.grid=element_line(size=0))
  

ggplot(na.omit(proteins), aes(Set)) +
  geom_bar() + ylab(sprintf('nr of %ss', featname)) + coord_flip() +
  theme(text=font, panel.grid=element_line(size=0), axis.ticks=ticks)

max_rm = function(vec) {
  val = max(vec, na.rm=T)
  if (val == Inf || val == -Inf) {
    return(NA)
  }
  return(val)
}

if (isobaric_labeled) {
  channels = names(proteins)[grepl('plex_[0-9]+[CN]*$', names(proteins))]
  isobarics = melt(proteins, id.vars='Set', measure.vars=channels)
  plots = list()
  for (i in 1:length(setnames)) local({
    i <- i
      plots[[i]] <<- ggplot(subset(isobarics, Set==setnames[i]), aes(x=variable, y=value)) +  
                            geom_boxplot() + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + ggtitle(setnames[i]) +
                            scale_y_log10() + ylab(sprintf('isobaric ratio of %s', featname))
  })
  multiplot(plotlist=plots, cols=2)
}
if (isobaric_labeled) {
  quanteds = names(proteins)[grepl('quanted', names(proteins))]
  proteins$max_chpsms = apply(proteins[c(quanteds)], 1, max_rm)
  quantednrs = melt(proteins, id.vars='Set', measure.vars="max_chpsms")
  quantednrs = transform(quantednrs, setrank=ave(value, Set, FUN = function(x) rank(x, ties.method = "random")))
  ggplot(quantednrs, aes(y=value, x=setrank)) + 
    geom_step(aes(color=Set), size=2) + 
    scale_y_log10() + ylab(sprintf('nr of PSMs per %s', featname)) +
    xlab(sprintf('%s rank', featname))
}

dev.off()
