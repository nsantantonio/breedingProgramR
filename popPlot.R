

seel(run1)
selPlot <- function(pop1, pop2){


getPopStats <- function(pop, meanVariety = TRUE){
    VgRGSC <- sapply(run1[["RGSC"]], genicVarG)
    gvRGSC <- sapply(run1[["RGSC"]], function(x) mean(gv(x)))
    # gvVariety <- sapply(run1[["VDP"]][["variety"]], function(x) mean(gv(x)))
    gvVariety <- lapply(run1[["VDP"]][["variety"]], function(x) gv(x))
    meanVariety <- TRUE
    if(meanVariety){
    	Yvariety <- sapply(gvVariety, mean)
    	Xvariety <- RGSCyr[1:length(Yvariety)]
    } else {
    	nVariety <- sapply(gvVariety, nrow)
    	Yvariety <- unlist(gvVariety)
    	Xvariety <- rep(RGSCyr[1:length(nVariety)], times = nVariety)
    }
    SDgRGSC <- sqrt(VgRGSC)
    RGSCgen <- c(0, 1:(tail(generation * GScylcePerYr, 1)))
    RGSCyr <- c(generation) * GScylcePerYr
    RGSCyrLab <- c(generation)

    return(list(Vg = VgRGSC, gv = gvRGSC, yr = RGSCyr, yrLab = RGSCyrLab, gen, sd, vx, vy)
}

run1[["VDP"]][["variety"]][[1]]@id # note, select doesnt select individuals, it selects phenotypes. fml

length(VgRGSC)

SDgRGSC <- sqrt(VgRGSC)
RGSCgen <- c(0, 1:(tail(generation * GScylcePerYr, 1)))
RGSCyr <- c(generation) * GScylcePerYr
RGSCyrLab <- c(generation)

polyCol <- "gray"
lineCol <- "black"
ptCol <- "black"


xpoly <- c(RGSCgen, rev(RGSCgen), RGSCgen[1])
ypoly <- c(gvRGSC + SDgRGSC, rev(gvRGSC - SDgRGSC), gvRGSC[1] + SDgRGSC[1])

xlims <- range(c(0, RGSCgen))
ylims <- range(c(ypoly, Yvariety) * 1.1)
plot(NA, xlim = xlims, ylim = ylims, xaxt = "n", xlab = "generation", ylab = "standardized genetic value")
axis(1, at = c(0, RGSCyr), labels = c(0, RGSCyrLab))
polygon(x = xpoly, y = ypoly, col = polyCol)
lines(x = RGSCgen, y = gvRGSC, type = "l", col = lineCol, lwd = 2)
points(Xvariety, Yvariety, col = ptCol, pch = 16)
legend("topright", legend = c("RGSC mean", expression(paste('RGSC ', sigma[g])), "Variety mean"), 
       # fill = c(NA, polyCol, NA), 
       bty = "n",
       col = c(lineCol, "black", ptCol),
       lty = c(1, 0, 0), lwd = c(2, 0, 0),
       pch = c(NA, 22, 16),
       pt.bg = c(NA, polyCol, ptCol),
)
}

pdf("genValByGen.pdf", width = 10, height = 7)
dev.off()
