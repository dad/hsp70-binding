setwd("~/research/spearman/analysis")

file.out = F
output.type = 'svg'

load.data = T
fig.abund.vs.sites = T
fig.motif.avoidance = T

master.fname <- paste(raw.data.dirname,"scer-raw-extended.Rdata",sep='/')

if (load.data) {
	yres <- get(load(master.fname))
	x <- read.table("~/research/hsp70-binding/data/scer-hsp70-motifs-summary-5.txt", header=T)
	z <- match(yres$bg$orf, x$orf)
	d <- data.frame(yres$bg, yres$est, x[z,])
	rc <- rcormat(d[,c('prot','mrna','cai','length','num.sites','num.motifs','min.score','prop.sites',p.0('num.motifs',1:10))])
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
