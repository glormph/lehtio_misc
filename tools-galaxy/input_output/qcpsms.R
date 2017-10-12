library(ggplot2)
library(reshape2)

# for testing interactive
#setwd()
#args = c('multiplot.R', test-data/preproc_psms.txt', 'test-data/amount_spectra_plates', 'WC', 'Cytsol')
#setnames = c('WC', 'Cytosol')

args = commandArgs(trailingOnly = T)
multiplot_fn = args[1]
amount_spec = read.table(args[3], header=T, sep='\t')
setnames = args[4:length(args)]

psms = read.table(args[2], header=T, sep='\t')
source(multiplot_fn)
if ('Fractions' %in% colnames(psms)) {
  is_ipg = T
} else {
  is_ipg = F
}

platenames = unique(psms$Plate_ID)
if (length(platenames) < 30) {
  pagelength = length(platenames) * 1.5
} else {
  pagelength = length(platenames) / 8
}
print(sprintf('Pagelemgth %s', pagelength))
pdf('psms.pdf', height=pagelength)

#IDed MS2 spectra
idcols = c('Plate_ID', 'amount')
if (is_ipg) {
  amount_psms = aggregate(Fractions~Plate_ID, psms, length)
} else {
  amount_psms = aggregate(Biological.set~Plate_ID, psms, length)
}

ticks = element_line(size=10/length(platenames))
if (length(platenames) > 50) {
  font = element_text(size=1000/length(platenames))
} else {
  font = element_text(size=100/length(platenames))
}

names(amount_psms) = idcols
amount_psms$counts = 'PSM-IDed'
amount_spec = aggregate(amount_ms2~Plate_ID, amount_spec, sum)
names(amount_spec) = idcols
amount_spec$counts = 'MS2-spectra'
amount_id = rbind(amount_psms, amount_spec)
ggplot(amount_id) +
  geom_bar(aes(Plate_ID, y=amount, fill=counts), stat='identity', position='dodge') +
  coord_flip() + theme(text=font, axis.ticks=ticks, panel.grid=element_line(size=0))


# empty channels of psms
if (length(grep('plex', names(psms)))) {
  channels = names(psms)[grepl('plex', names(psms))]
  psm_empty = melt(psms[c("Plate_ID", channels)], id.vars="Plate_ID")
  psm_empty = psm_empty[psm_empty$value==0,]
  psm_empty$value = 1
  psm_empty = aggregate(value~Plate_ID + variable, psm_empty, sum)
  names(psm_empty) = c('Plate_ID', 'channel', 'nr_missing_values')
  ggplot(psm_empty, aes(Plate_ID, nr_missing_values)) + 
    geom_bar(aes(fill=channel), stat='identity', position="dodge")
}

# RT per plate
if (length(platenames) < 20) {
plots = list()
for (i in 1:length(setnames)) local({
  i <- i
  plots[[i]] <<- ggplot(subset(psms, Biological.set==setnames[i]), aes(x=Retention.time.min., fill=Plate_ID)) + 
    geom_histogram(position='identity', bins=100, alpha=.2) +
    ggtitle(setnames[i])
})
multiplot(plotlist=plots, cols=2)
# mass error plot
plots = list()
for (i in 1:length(setnames)) local({
  i <- i
  plots[[i]] <<- ggplot(subset(psms, Biological.set==setnames[i]), aes(x=PrecursorError.ppm., fill=Plate_ID)) + 
    geom_histogram(position='identity', binwidth=0.4, alpha=0.4) +
    ggtitle(setnames[i])
})
multiplot(plotlist=plots, cols=2)
}

# delta pI per plate
if (is_ipg) {
  plots = list()
  for (i in 1:length(setnames)) local({
    i <- i
    plots[[i]] <<- ggplot(subset(psms[complete.cases(psms[,"Delta.pI"]),], Biological.set==setnames[i]), aes(x=Delta.pI, fill=Plate_ID)) + 
      geom_histogram(position='identity', bins=100, alpha=.2) +
      ggtitle(setnames[i])
    })
    multiplot(plotlist=plots, cols=2)

  # fraction yield
  psm_fracs = psms[order(psms$percolator.svm.score, decreasing=T),]
  psm_fracs = psm_fracs[!duplicated(psm_fracs[c("Plate_ID", "Peptide")]),]
  psm_fracs = dcast(psm_fracs[,c("Biological.set", "Plate_ID", "Fractions")], Biological.set+Plate_ID~Fractions)
  psm_fracs = melt(psm_fracs, id.vars=c('Plate_ID', 'Biological.set'))
  colnames(psm_fracs) <- c('Plate_ID', 'Biological.set', 'fraction_nr', 'amount_peptides')
  plots = list()
  for (i in 1:length(setnames)) local({
    i <- i
    plots[[i]] <<- ggplot(subset(psm_fracs, Biological.set==setnames[i]), aes(fraction_nr, amount_peptides)) + 
      geom_bar(aes(fill=Plate_ID), alpha=0.2, stat='identity', position=position_identity()) +
      ggtitle(setnames[i])
  })
  multiplot(plotlist=plots, cols=2)
}



# missed cleavages
psms_wide = dcast(psms[,c("Plate_ID", "missed_cleavage")], Plate_ID~missed_cleavage)
sum_psms = apply(as.data.frame(psms_wide[,c(2:length(psms_wide))]), 1, sum)
for (mcno in names(psms_wide)[c(2:length(psms_wide))]) {
  psms_wide[,mcno] = psms_wide[,mcno]/sum_psms * 100
}
psms_mc = melt(psms_wide, id.vars='Plate_ID')
colnames(psms_mc) <- c('Plate_ID', 'nr_missed_cleavages', 'percent_PSMs')
ggplot(psms_mc, aes(Plate_ID, percent_PSMs)) + 
  geom_bar(aes(fill=nr_missed_cleavages), stat='identity', position="dodge") +
  scale_y_continuous(breaks=seq(0,100, 10)) +
  coord_flip() + theme(text=font, axis.ticks=ticks, panel.grid=element_line(size=0))

dev.off()
