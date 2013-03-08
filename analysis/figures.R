setwd("~/research/spearman/analysis")

file.out = F
output.type = 'svg'

load.data = T
fig.abund.vs.sites = F
fig.motif.avoidance = F
fig.paralog = T

master.fname <- paste(raw.data.dirname,"scer-raw-extended.Rdata",sep='/')

if (load.data) {
	yres <- get(load(master.fname))
	x <- read.table("~/research/hsp70-binding/data/scer-hsp70-motifs-summary-vd0.txt", header=T)
	#x <- read.table("~/research/hsp70-binding/data/scer-hsp70-motifs-summary-vd5.txt", header=T)
	#x <- read.table("~/research/hsp70-binding/data/scer-hsp70-motifs-summary-r5.txt", header=T)
	z <- match(yres$bg$orf, x$orf)
	d <- data.frame(yres$bg, yres$est, x[z,])
	rc <- rcormat(d[,c('prot','mrna','cai','length','num.sites','num.motifs','max.score','prop.sites',p.0('num.motifs',1:10))])
}

if (fig.abund.vs.sites) {
	
	if (file.out) dev.out('fig.abund.propsite', width=4, height=4, output.type=output.type)
	plot(d$prot, d$prop.sites, log='x', las=1, xlab='Protein level', ylab='Proportion of sequence scored as Hsp70 binding site')
	if (file.out) dev.off()
	pcor(cortest(d$prot, d$prop.sites))
}

if (fig.motif.avoidance) {
	n.mot <- 1:10
	if (file.out) dev.out('fig.abund.propsite', width=4, height=4, output.type=output.type)
	plot(n.mot, rc$r['prot',p.0('num.motifs',n.mot)], type='l', col='red', ylim=c(-0.15,0.05))
	lines(n.mot, rc$r['mrna',p.0('num.motifs',n.mot)], col='blue')
	lines(n.mot, rc$r['cai',p.0('num.motifs',n.mot)], col='green')
	if (file.out) dev.off()
}

if (fig.paralog) {
	pairorf <- subset(d, !is.na(pair.orf))
	zp <- match(pairorf$orf, d$pair.orf)
	zo <- match(pairorf$orf, d$orf)
	
	ratio.flds <- c('prot','mrna', 'length','prop.sites','num.sites','num.motifs.3')
	all.flds <- c('orf','gene',ratio.flds)
	dp <- data.frame(d[zo,all.flds], d[zp,all.flds], d[zo,ratio.flds]/d[zp,ratio.flds])
	colnames(dp) <- c(all.flds, p.0('pair',all.flds), p.0('ratio',ratio.flds))
	plot(dp$ratio.prot, dp$ratio.num.sites, log='xy')
	abline(0,1)
	pcor(cortest(dp$ratio.prot, dp$ratio.num.motifs.3))
	pcor(cortest(dp$ratio.prot, dp$ratio.prop.sites))
	pcor(cortest(dp$ratio.prot, dp$ratio.num.sites))
	pcor(cortest(dp$ratio.mrna, dp$ratio.prop.sites))
	pcor(cortest(dp$ratio.mrna, dp$ratio.num.sites))
}
