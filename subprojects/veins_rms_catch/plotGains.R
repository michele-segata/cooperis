library(ggplot2)
library(viridis)
library('tikzDevice')
library(RColorBrewer)

source('plot-utils.R')


plot.width <- 308.43/72*1.5
plot.height <- 184.546607/72/1.25*1.5
plot.far.field <- function(outputFile, d, norm.db=F) {
    # p <- ggplot(d, aes(x=phi, y=theta)) + geom_raster(aes(fill=gain)) + scale_fill_viridis(option="turbo", limits = c(0, 300)) +
    lim <- c(0, 70)
    if (norm.db) {
        d$gain <- d$gain / max(d$gain)
        d$gain <- 10*log10(d$gain)
        d$gain[d$gain < -30] <- -30
        lim <- c(-30, 0)
    }
    p <- ggplot(d, aes(x=phi, y=theta)) + geom_raster(aes(fill=gain)) + scale_fill_viridis(option="turbo", limits=lim) +
    scale_x_continuous(limits = c(-180, 180), expand = c(0, 0), breaks=seq(-180, 180, by=45)) +
    scale_y_continuous(limits = c(0, 90), expand = c(0, 0), breaks=seq(0, 90, by=15)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          legend.title=element_blank(),
          legend.key.height = unit(0.5, "inch"))+
     xlab("reflection azimuth [degrees]") +
     ylab("reflection elevation [degrees]")
     ggsave(paste0(outputFile, ".pdf"), p, units="in", width=plot.width, height=plot.height, device=cairo_pdf)
}

#graphs params
plot.margin <- c(3.1, 3.3, 0.1, 0.1)
leg.inset <- c(0, 0, 0, 0)

# nicer, but with xquartz it shows horrible white lines. with x11() it is nice
#plot.far.field <- function(outputFile, d) {
#
#    done <- myps(outputFile = outputFile, width=plot.width, height=plot.height)
#    par(mar=plot.margin, xpd=F)
#    palette(turbo(100))
#    plot.new()
#
#    plot.window(xlim=c(0, (phis-1)), ylim=c(0, (thetas-1)), yaxs="i", xaxs="i")
#
#    filled.contour(x=0:(phis-1), y=0:(thetas-1), z=g, levels=pretty(c(0, 300), 100), col=palette(), grid(col=rgb(1,1,1,0)))
#
#    axis(1, lwd=0, lwd.ticks=1, las=1)
#    axis(2, lwd=0, lwd.ticks=1, las=1)
#    title(ylab='elevation $\\theta_r$ [$^{\\circ}$]', line=2.2)
#    title(xlab='azimuth $\\varphi_r$ [$^{\\circ}$]', line=2)
#    box()
#    done()
##filled.contour(z=g)
##par(fg = NA)
#x11()
#filled.contour(x=0:(phis-1), y=0:(thetas-1), z=g, zlim=c(0, 300), nlevels=100, color.palette=function(n){turbo(n)})
#filled.contour(x=0:(phis-1), y=0:(thetas-1), z=g, nlevels=100, color.palette=function(n){turbo(n)})
#}

files <- c("gains_-45_45.csv", "gains_-30_45.csv", "gains_-80_80.csv", "gains_10_10.csv", "gains_75_10.csv")

for (f in files) {
    d <- read.csv(f)
    print(max(d$gain))
    plot.far.field(gsub(f, pattern=".csv", replacement=""), d)
    plot.far.field(paste0(gsub(f, pattern=".csv", replacement=""), "_db"), d, norm.db=T)
}

#thetas <- 91
#phis <- 181
#g <- t(matrix(d$gain, nrow=thetas, ncol=phis, byrow=F))
