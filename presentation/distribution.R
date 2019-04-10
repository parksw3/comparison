pdf("distribution.pdf", width=6, height=4)
par(mar=c(5, 4, 2, 2) + 0.1)
curve(1/13* exp(-x/13), xlim=c(0, 25), ylim=c(0, 0.2),
	  xlab="Generation time (days)",
	  ylab="Density",
	  bty="n", lwd=2)
curve(dgamma(x, shape=25, scale=13/25), add=TRUE, col=2, lwd=2)

legend(
	x=0, y=0.2,
	legend=c("Exponential", "Realistic"),
	col=c(1, 2),
	lty=c(1, 1),
	lwd=2
)
dev.off()
