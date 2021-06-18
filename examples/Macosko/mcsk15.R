# -----------------------------------------------------------------------------
# +++ ~/ utilitities
# -----------------------------------------------------------------------------

MACOSKO_COLORS1 <- c(
	"#A5C93D",
	"#8B006B",
	"#2000D7",
	"#538CBA",
	"#8B006B",
	"#B33B19",
	"#8B006B",
	"#8B006B",
	"#8B006B",
	"#C38A1F",
	"#538CBA",
	"#FFFF00"
)

MACOSKO_COLORS2 <- c(
	"#B33B19",
	"#C38A1F",
	rep("#9AC222", 21),
	rep("#538CBA", 2),
	rep("#2000D7", 8),
	rep("#8B006B", 6)
)

MACOSKO_LABELS1 <- c(
	"Amacrine cells",
	"Astrocytes",
	"Bipolar cells",
	"Cones",
	"Fibroblasts",
	"Horizontal cells",
	"Microglia",
	"Muller glia",
	"Pericytes",
	"Retinal ganglion cells",
	"Rods",
	"Vascular endothelium"
)

mcsk15.legend <- function(pltt = MACOSKO_COLORS1, text = MACOSKO_LABELS1, ncol = 3, cex = 1.0, pt.cex = 1.0)
{
	par(mar = c(0.1, 1.5, 0.1, 1.0), oma = c(0.0, 0.1, 0.0, 0.1))
	plot(1, 1, xlab = '', ylab = '', xaxt = "n", yaxt = "n", bty = "n", type = "n", pin = c(10.0, 5.0))
	legend('left', legend = text, bty = 'n', pch = 15, cex = cex, pt.cex = pt.cex, col = pltt, ncol = ncol)
}

mcsk15.lbls <- function(l = 1)
{
	L <- as.matrix(read.csv('mcsk15_lbls.csv.gz'))
	L[, l]
}
