require('tikzDevice')
options(tikzLatexPackages=c(getOption('tikzLatexPackages'), "\\usetikzlibrary{backgrounds}"))

plot.width <- 308.43/72
plot.height <- 184.546607/72
plot.margin <- c(3.1, 3.6, 0.1, 0.1)
leg.inset <- c(0, 0, 0, 0)

myps <- function(outputFile, paper="special", width=6, height=4, family="Times", outputDir='./') {
	tikz(paste(outputDir, outputFile, ".tex", sep=""), standAlone = TRUE, width=width, height=height)
	done <- function() {
		graphics.off()
		system(paste("pdflatex -output-directory ", outputDir, " ", outputFile, ".tex", sep=""))

		system(paste('rm', paste(paste(outputDir, outputFile, c('.tex', '.aux', '.log'), sep=''), collapse=' ')))
	}
	return(done)
}

#define own function to plot arrows (i.e., confidence intervals, avoiding the problem with too small confidence intervals)
arrows <- function(x0, y0, x1=x0, y1=y0, col = 1, lwd = 2, lty = "solid", angle = 90, length = 1, code = 1) {
	segments(x0, y0, x1, y1, col=col, lwd=lwd, lty=lty)
	l <- grconvertX(length, from="inches", to="device")*2
	segments(x0-l, y0, x0+l, y0, col=col, lwd=lwd, lty=lty)
#	segments(x0-l, y1, x0+l, y1, col=col, lwd=lwd, lty=lty)
}
