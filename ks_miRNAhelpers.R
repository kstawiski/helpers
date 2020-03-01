library(plyr)
library(dplyr)
library(edgeR)
library(epiDisplay)
library(rsq)
library(MASS)
library(Biocomb)
library(caret)
library(dplyr)
library(epiDisplay)
library(pROC)
library(ggplot2)
library(DMwR)
library(ROSE)
library(gridExtra)
library(gplots)
library(devtools)
library(stringr)
library(data.table)
library(tidyverse)

ks.wykreskorelacji = function(var1, var2, labvar1, labvar2, title, yx = T, metoda = 'pearson', gdzie_legenda = "topleft") {
  var1 = as.numeric(as.character(var1))
  var2 = as.numeric(as.character(var2))
  plot(var1, var2, pch = 19, cex=0.5, xlab=labvar1, ylab=labvar2, main=title)
  if (yx == T) { abline(0,1, col='gray') }
  abline(fit <- lm(var2 ~ var1), col='black', lty = 2)
  temp = cor.test(var1, var2, method = metoda)
  if (metoda=='pearson') {
    if (temp$p.value < 0.0001) {
      legend(gdzie_legenda, bty="n", legend=paste("r =",round(cor(var1,var2,use="complete.obs"),2),", p < 0.0001\nadj. R² =", format(summary(fit)$adj.r.squared, digits=4)))
    } else {
      legend(gdzie_legenda, bty="n", legend=paste("r =",round(cor(var1,var2,use="complete.obs"),2),", p =",round(temp$p.value, 4),"\nadj. R² =", format(summary(fit)$adj.r.squared, digits=4))) } }
  if (metoda=="spearman")
  { if (temp$p.value < 0.0001) {
    legend(gdzie_legenda, bty="n", legend=paste("rho =",round(cor(var1,var2,use="complete.obs"),2),", p < 0.0001"))
  } else {
    legend(gdzie_legenda, bty="n", legend=paste("rho =",round(cor(var1,var2,use="complete.obs"),2),", p =",round(temp$p.value, 4))) } } 
  
}

ks.heatmap.3 <- function(x,
                      Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                      distfun = dist,
                      hclustfun = hclust,
                      dendrogram = c("both","row", "column", "none"),
                      symm = FALSE,
                      scale = c("none","row", "column"),
                      na.rm = TRUE,
                      revC = identical(Colv,"Rowv"),
                      add.expr,
                      breaks,
                      symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                      col = "heat.colors",
                      colsep,
                      rowsep,
                      sepcolor = "white",
                      sepwidth = c(0.05, 0.05),
                      cellnote,
                      notecex = 1,
                      notecol = "cyan",
                      na.color = par("bg"),
                      trace = c("none", "column","row", "both"),
                      tracecol = "cyan",
                      hline = median(breaks),
                      vline = median(breaks),
                      linecol = tracecol,
                      margins = c(5,5),
                      ColSideColors,
                      RowSideColors,
                      side.height.fraction=0.3,
                      cexRow = 0.2 + 1/log10(nr),
                      cexCol = 0.2 + 1/log10(nc),
                      labRow = NULL,
                      labCol = NULL,
                      key = TRUE,
                      keysize = 1.5,
                      density.info = c("none", "histogram", "density"),
                      denscol = tracecol,
                      symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                      densadj = 0.25,
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      lmat = NULL,
                      lhei = NULL,
                      lwid = NULL,
                      ColSideColorsSize = 1,
                      RowSideColorsSize = 1,
                      KeyValueName="Value",...){
  
  invalid <- function (x) {
    if (missing(x) || is.null(x) || length(x) == 0)
      return(TRUE)
    if (is.list(x))
      return(all(sapply(x, invalid)))
    else if (is.vector(x))
      return(all(is.na(x)))
    else return(FALSE)
  }
  
  x <- as.matrix(x)
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale))
    "none"
  else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (length(col) == 1 && is.character(col))
    col <- get(col, mode = "function")
  if (!missing(breaks) && (scale != "none"))
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
  if (is.null(Rowv) || is.na(Rowv))
    Rowv <- FALSE
  if (is.null(Colv) || is.na(Colv))
    Colv <- FALSE
  else if (Colv == "Rowv" && !isTRUE(Rowv))
    Colv <- FALSE
  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2)
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote))
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                 c("both", "row"))) {
      if (is.logical(Colv) && (Colv))
        dendrogram <- "column"
      else dedrogram <- "none"
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting row dendogram.")
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                 c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv))
        dendrogram <- "row"
      else dendrogram <- "none"
      warning("Discrepancy: Colv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting column dendogram.")
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  }
  else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc)
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  }
  else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow))
    labRow <- if (is.null(rownames(x)))
      (1:nr)[rowInd]
  else rownames(x)
  else labRow <- labRow[rowInd]
  if (is.null(labCol))
    labCol <- if (is.null(colnames(x)))
      (1:nc)[colInd]
  else colnames(x)
  else labCol <- labCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
    if (missing(col) || is.function(col))
      breaks <- 16
    else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks)
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function")
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei))
    lhei <- c(keysize, 4)
  if (missing(lwid) || is.null(lwid))
    lwid <- c(keysize, 4)
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    
    if (!missing(ColSideColors)) {
      #if (!is.matrix(ColSideColors))
      #stop("'ColSideColors' must be a matrix")
      if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
        stop("'ColSideColors' must be a matrix of nrow(x) rows")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      #lhei <- c(lhei[1], 0.2, lhei[2])
      lhei=c(lhei[1], side.height.fraction*ColSideColorsSize/2, lhei[2])
    }
    
    if (!missing(RowSideColors)) {
      #if (!is.matrix(RowSideColors))
      #stop("'RowSideColors' must be a matrix")
      if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
        stop("'RowSideColors' must be a matrix of ncol(x) columns")
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
      #lwid <- c(lwid[1], 0.2, lwid[2])
      lwid <- c(lwid[1], side.height.fraction*RowSideColorsSize/2, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }
  
  if (length(lhei) != nrow(lmat))
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat))
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  
  if (!missing(RowSideColors)) {
    if (!is.matrix(RowSideColors)){
      par(mar = c(margins[1], 0, 0, 0.5))
      image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    } else {
      par(mar = c(margins[1], 0, 0, 0.5))
      rsc = t(RowSideColors[,rowInd, drop=F])
      rsc.colors = matrix()
      rsc.names = names(table(rsc))
      rsc.i = 1
      for (rsc.name in rsc.names) {
        rsc.colors[rsc.i] = rsc.name
        rsc[rsc == rsc.name] = rsc.i
        rsc.i = rsc.i + 1
      }
      rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
      image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
      if (length(rownames(RowSideColors)) > 0) {
        axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  if (!missing(ColSideColors)) {
    
    if (!is.matrix(ColSideColors)){
      par(mar = c(0.5, 0, 0, margins[2]))
      image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    } else {
      par(mar = c(0.5, 0, 0, margins[2]))
      csc = ColSideColors[colInd, , drop=F]
      csc.colors = matrix()
      csc.names = names(table(csc))
      csc.i = 1
      for (csc.name in csc.names) {
        csc.colors[csc.i] = csc.name
        csc[csc == csc.name] = csc.i
        csc.i = csc.i + 1
      }
      csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
      image(csc, col = as.vector(csc.colors), axes = FALSE)
      if (length(colnames(ColSideColors)) > 0) {
        axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr"))
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
  retval$carpet <- x
  if (exists("ddr"))
    retval$rowDendrogram <- ddr
  if (exists("ddc"))
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
          col = na.color, add = TRUE)
  }
  axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
       cex.axis = cexCol)
  if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
       cex.axis = cexRow)
  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr))
    eval(substitute(add.expr))
  if (!missing(colsep))
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  if (!missing(rowsep))
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol,
               lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        abline(h = i + hline, col = linecol, lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote))
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
         col = notecol, cex = notecex)
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }
  else plot.new()
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  }
  else plot.new()
  if (!is.null(main))
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    }
    else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    
    z <- seq(min.raw, max.raw, length = length(col))
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv)
    if (scale == "row")
      mtext(side = 1, "Row Z-Score", line = 2)
    else if (scale == "column")
      mtext(side = 1, "Column Z-Score", line = 2)
    else mtext(side = 1, KeyValueName, line = 2)
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
            lwd = 1)
      axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
      title("Color Key\nand Density Plot")
      par(cex = 0.5)
      mtext(side = 2, "Density", line = 2)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
            col = denscol)
      axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
      title("Color Key\nand Histogram")
      par(cex = 0.5)
      mtext(side = 2, "Count", line = 2)
    }
    else title("Color Key")
  }
  else plot.new()
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                  high = retval$breaks[-1], color = retval$col)
  invisible(retval)
}

# DE funkcja
ks.miRNA.de = function(ttpm_pofiltrze, klasy)
{
  wyniki = data.frame(miR = colnames(ttpm_pofiltrze))
  
  library(dplyr)
  
  for (i in 1:length(colnames(ttpm_pofiltrze)))
  {
    # Filtry
    #wyniki[i,"jakieszero"] = ifelse(sum(ttpm_pofiltrze[,i] == 0) > 0, "tak", "nie")
    
    # Średnia i SD
    wyniki[i,"mean logtpm"] = mean(ttpm_pofiltrze[,i])
    wyniki[i,"median logtpm"] = median(ttpm_pofiltrze[,i])
    wyniki[i,"SD logtpm"] = sd(ttpm_pofiltrze[,i])
    
    # Cancer
    tempx = ttpm_pofiltrze[klasy == "Cancer",]
    wyniki[i,"cancer mean"] = mean(tempx[,i])
    wyniki[i,"cancer median"] = median(tempx[,i])
    wyniki[i,"cancer SD"] = sd(tempx[,i])
    
    # Cancer
    tempx = ttpm_pofiltrze[klasy != "Cancer",]
    wyniki[i,"control mean"] = mean(tempx[,i])
    wyniki[i,"control median"] = median(tempx[,i])
    wyniki[i,"control SD"] = sd(tempx[,i])
    
    # DE
    temp = t.test(ttpm_pofiltrze[,i] ~ as.factor(klasy))
    fc = (temp$estimate[1] - temp$estimate[2])
    wyniki[i,"log10FC (subtr estim)"] = fc
    wyniki[i,"log10FC"] = wyniki[i,"cancer mean"] - wyniki[i,"control mean"]
    wyniki[i,"log2FC"] = wyniki[i,"log10FC"] / log10(2)
    
    revfc = (temp$estimate[2] - temp$estimate[1])
    wyniki[i,"reversed_log10FC"] = revfc
    wyniki[i,"reverse_log2FC"] = wyniki[i,"reversed_log10FC"] / log10(2)
    #wyniki[i,"log2FC"] = log2(fc)
    wyniki[i,"p-value"] = temp$p.value
  }
  wyniki[,"p-value Bonferroni"] = p.adjust(wyniki$`p-value`, method = "bonferroni")
  wyniki[,"p-value Holm"] = p.adjust(wyniki$`p-value`, method = "holm")
  wyniki[,"-log10(p-value Bonferroni)"] = -log10(p.adjust(wyniki$`p-value`, method = "bonferroni"))
  wyniki[,"p-value BH"] = p.adjust(wyniki$`p-value`, method = "BH")
  
  return(wyniki %>% arrange(`p-value BH`))
}

ks.wykreskorelacji = function(var1, var2, labvar1, labvar2, title) {
  plot(var1, var2, pch = 19, cex=0.5, xlab=labvar1, ylab=labvar2, main=title)
  abline(0,1, col='gray')
  abline(fit <- lm(var2 ~ var1), col='black', lty = 2)
  temp = cor.test(var1, var2, method = 'pearson')
  legend("topleft", bty="n", legend=paste("r =",round(cor(var1,var2),2),", p =",round(temp$p.value, 4),"\nR² =", format(summary(fit)$adj.r.squared, digits=4)))
}

ks.miR.formula = function(wybrane_miRy) {
  as.formula(paste0("Class ~ ",paste0(as.character(wybrane_miRy), collapse = " + ")))
}

ks.wczytajmix = function(wd = getwd(), smote_over = 10000, use_smote_not_rose = F, replace_smote = F) {
  oldwd = getwd()
  setwd(wd)
  train = dplyr::select(read.csv("mixed_train.csv", stringsAsFactors = F), starts_with("hsa"), Class)
  
  test = dplyr::select(read.csv("mixed_test.csv", stringsAsFactors = F), starts_with("hsa"), Class)
  valid = dplyr::select(read.csv("mixed_valid.csv", stringsAsFactors = F), starts_with("hsa"), Class)
  train$Class = factor(train$Class, levels = c("Control","Cancer"))
  test$Class = factor(test$Class, levels = c("Control","Cancer"))
  valid$Class = factor(valid$Class, levels = c("Control","Cancer"))
  
  # Wywalamy miRy z zerowa wariancja
  temp = train %>% filter(Class == "Cancer")
  temp2 = as.numeric(which(apply(temp, 2, var) == 0))
  temp = train %>% filter(Class == "Control")
  temp3 = as.numeric(which(apply(temp, 2, var) == 0))
  temp4 = unique(c(temp2, temp3))
  if (length(temp4) > 0) {
  train = train[,-temp4]
  test = test[,-temp4]
  valid = valid[,-temp4]
  }
  
  if(use_smote_not_rose) {
    train_smoted = DMwR::SMOTE(Class ~ ., data = train, perc.over = smote_over,perc.under=100, k=10)
  train_smoted$Class = factor(train_smoted$Class, levels = c("Control","Cancer"))
  } else {
  rosed = ROSE(Class ~ ., data = train, N = nrow(train)*10, seed = 1)
  train_smoted = rosed[["data"]]
  train_smoted$Class = factor(train_smoted$Class, levels = c("Control","Cancer"))
  }
  
  
  
  if (replace_smote == T) { train = train_smoted }
  
  trainx = dplyr::select(train, starts_with("hsa"))
  trainx_smoted = dplyr::select(train_smoted, starts_with("hsa"))
  setwd(oldwd)
  return(list(train, test, valid, train_smoted, trainx, trainx_smoted))
}

ks.selekcjamiRNA = function(wd = getwd(), m = c(1:70),
                            max_iterations = 10, code_path = "/home/konrad/snorlax/2019_PRELUDIUM/feature_selection/",
                            register_parallel = T, clx = NULL, stamp = as.numeric(Sys.time()),
                            prefer_no_features = 11, conda_path = "/home/konrad/anaconda3/bin/conda") {

  oldwd = getwd()
  setwd(wd)
  

  if(!dir.exists("temp")) { dir.create("temp") }
  
  run_id = stamp
  formulas = list()
  times = list()
  
  zz <- file(paste0("temp/",stamp,"featureselection.log"), open = "wt")
  sink(zz)
  #sink(zz, type = "message")
  pdf(paste0("temp/",stamp,"featureselection.pdf"))
  
  
  
  
  dane = ks.wczytajmix(); train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]
  
  n = 1
  if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); cat(paste0("\n\nMatched method ", n, " with those requested.. Starting.."));
  start_time <- Sys.time()
    
  formulas[["all"]] = ks.miR.formula(colnames(trainx))
  
  end_time <- Sys.time()
  saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS"))
  saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
 save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  cat("\n\nStandard DE...")
  wyniki = ks.miRNA.de(trainx, train$Class)
  istotne = filter(wyniki, `p-value BH` <= 0.05) %>% arrange(`p-value BH`)
  istotne_top = wyniki %>% arrange(`p-value BH`) %>% head(prefer_no_features)
  istotne_topBonf = wyniki %>% arrange(`p-value Bonferroni`) %>% head(prefer_no_features)
  istotne_topHolm = wyniki %>% arrange(`p-value Holm`) %>% head(prefer_no_features)
  istotne_topFC = wyniki %>% arrange(desc(abs(`log2FC`))) %>% head(prefer_no_features)

  train_sig = dplyr::select(train, as.character(istotne$miR), Class)
  trainx_sig = dplyr::select(train, as.character(istotne$miR))
  
  cat("\n\nPreparing for SMOTE...")
  #train_sig_smoted = DMwR::SMOTE(Class ~ ., data = train_sig, perc.over = 10000,perc.under=100, k=10)
  train_sig_smoted = dplyr::select(train_smoted, as.character(istotne$miR), Class)
  train_sig_smoted$Class = factor(train_sig_smoted$Class, levels = c("Control","Cancer"))
  trainx_sig_smoted = dplyr::select(train_sig_smoted, starts_with("hsa"))
  
  wyniki_smoted = ks.miRNA.de(trainx_smoted, train_smoted$Class)
  istotne_smoted = filter(wyniki_smoted, `p-value BH` <= 0.05) %>% arrange(`p-value BH`)
  istotne_top_smoted = wyniki_smoted %>% arrange(`p-value BH`) %>% head(prefer_no_features)
  istotne_topBonf_smoted = wyniki_smoted %>% arrange(`p-value Bonferroni`) %>% head(prefer_no_features)
  istotne_topHolm_smoted = wyniki_smoted %>% arrange(`p-value Holm`) %>% head(prefer_no_features)
  istotne_topFC_smoted = wyniki_smoted %>% arrange(desc(abs(`log2FC`))) %>% head(prefer_no_features)
  
  # Caret prep
  if (register_parallel) {
  cat("\n\nGetting cluster ready...")
  if(is.null(clx)) {
  library(doParallel)
  cl <- makePSOCKcluster(detectCores() - 2)
  registerDoParallel(cl) }
    else { registerDoParallel(clx) }
  }
  
  # 0. All and sig
  cat("\n\nStarting selected methods...")
  
  
  n = n + 1
  if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); 
    cat("\n\nStarting SIG")
    start_time <- Sys.time()
  
  formulas[["sig"]] = ks.miR.formula(as.character(istotne$miR))
  formulas[["sigtop"]] = ks.miR.formula(as.character(istotne_top$miR))
  formulas[["sigtopBonf"]] = ks.miR.formula(as.character(istotne_topBonf$miR))
  formulas[["sigtopHolm"]] = ks.miR.formula(as.character(istotne_topHolm$miR))
  formulas[["topFC"]] = ks.miR.formula(as.character(istotne_topFC$miR))
  formulas[["sigSMOTE"]] = ks.miR.formula(as.character(istotne_smoted$miR))
  formulas[["sigtopSMOTE"]] = ks.miR.formula(as.character(istotne_top_smoted$miR))
  formulas[["sigtopBonfSMOTE"]] = ks.miR.formula(as.character(istotne_topBonf_smoted$miR))
  formulas[["sigtopHolmSMOTE"]] = ks.miR.formula(as.character(istotne_topHolm_smoted$miR))
  formulas[["topFCSMOTE"]] = ks.miR.formula(as.character(istotne_topFC_smoted$miR))
  
  
  end_time <- Sys.time()
  saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS"))
  saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
 save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  # 1. Fold-change and sig. filter
  n = n + 1
  if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting.."));
  start_time <- Sys.time()
  
  cat("\n\nStarting FC and SIG")
  fcsig = as.character(istotne$miR[abs(istotne$log2FC)>1])
  formulas[["fcsig"]] = ks.miR.formula(fcsig)
  
  fcsig = as.character(istotne_smoted$miR[abs(istotne_smoted$log2FC)>1])
  formulas[["fcsigSMOTE"]] = ks.miR.formula(fcsig)
  
  end_time <- Sys.time()
  saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS"))
  saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
 save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  # 2. CFS
  n = n + 1
  if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting.."));
    start_time <- Sys.time()
  cat("\n\nStarting CFS")
  cfs = select.cfs(train)
  formulas[["cfs"]] = ks.miR.formula(as.character(cfs$Biomarker))
  
  cfs = select.cfs(train_smoted)
  formulas[["cfsSMOTE"]] = ks.miR.formula(as.character(cfs$Biomarker))
  
  cfs = select.cfs(train_sig)
  formulas[["cfs_sig"]] = ks.miR.formula(as.character(cfs$Biomarker))
  
  cfs = select.cfs(train_sig_smoted)
  formulas[["cfsSMOTE_sig"]] = ks.miR.formula(as.character(cfs$Biomarker))
  end_time <- Sys.time()
  saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS"))
  saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
 save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  # 3. classifier.loop
  n = n + 1
  if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting.."));
  start_time <- Sys.time()
  cat("\n\nStarting classifier loop")
  classloop = classifier.loop(train, feature.selection = "auc", method.cross="fold-crossval", classifiers=c("svm","lda","rf","nsc"), no.feat=prefer_no_features)
  f_classloop = rownames(classloop$no.selected)[classloop$no.selected[,1]>0]
  formulas[["classloop"]] = ks.miR.formula(f_classloop)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  classloop = classifier.loop(train_smoted, feature.selection = "auc", method.cross="fold-crossval", classifiers=c("svm","lda","rf","nsc"), no.feat=prefer_no_features)
  f_classloop = rownames(classloop$no.selected)[classloop$no.selected[,1]>0]
  formulas[["classloopSMOTE"]] = ks.miR.formula(f_classloop)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  classloop = classifier.loop(train_sig, feature.selection = "auc", method.cross="fold-crossval", classifiers=c("svm","lda","rf","nsc"), no.feat=prefer_no_features)
  f_classloop = rownames(classloop$no.selected)[classloop$no.selected[,1]>0]
  formulas[["classloop_sig"]] = ks.miR.formula(f_classloop)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  classloop = classifier.loop(train_sig, feature.selection = "auc", method.cross="fold-crossval", classifiers=c("svm","lda","rf","nsc"), no.feat=prefer_no_features)
  f_classloop = rownames(classloop$no.selected)[classloop$no.selected[,1]>0]
  formulas[["classloopSMOTE_sig"]] = ks.miR.formula(f_classloop)
  
  end_time <- Sys.time()
  saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS"))
  saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
 save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  # 4. select.forward.Corr
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  cat("\n\nStarting select.forward.Corr")
  fcfs = select.forward.Corr(train, disc.method="MDL")
  formulas[["fcfs"]] = ks.miR.formula(fcfs)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  fcfs = select.forward.Corr(train_smoted, disc.method="MDL")
  formulas[["fcfsSMOTE"]] = ks.miR.formula(fcfs)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  fcfs = select.forward.Corr(train_sig, disc.method="MDL")
  formulas[["fcfs_sig"]] = ks.miR.formula(fcfs)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  fcfs = select.forward.Corr(train_sig_smoted, disc.method="MDL")
  formulas[["fcfsSMOTE_sig"]] = ks.miR.formula(fcfs)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  # 5. select.forward.wrapper
  fwrap = select.forward.wrapper(train)
  formulas[["fwrap"]] = ks.miR.formula(fwrap)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  fwrap = select.forward.wrapper(train_smoted)
  formulas[["fwrapSMOTE"]] = ks.miR.formula(fwrap)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  fwrap = select.forward.wrapper(train_sig)
  formulas[["fwrap_sig"]] = ks.miR.formula(fwrap)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  fwrap = select.forward.wrapper(train_sig_smoted)
  formulas[["fwrapSMOTE_sig"]] = ks.miR.formula(fwrap)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
 save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  # 6, 7, 8.
  #select.process(dattable,method="InformationGain",disc.method="MDL",
  #               threshold=0.2,threshold.consis=0.05,attrs.nominal=numeric(),
  #               max.no.features=10)
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  cat("\n\nStarting MDL")
  formulas[["AUC_MDL"]] = ks.miR.formula(colnames(train)[select.process(train, method="auc", disc.method = "MDL", max.no.features = prefer_no_features)])
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  formulas[["SU_MDL"]] = ks.miR.formula(colnames(train)[select.process(train, method="symmetrical.uncertainty", disc.method = "MDL", max.no.features = prefer_no_features)])
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  formulas[["CorrSF_MDL"]] = ks.miR.formula(colnames(train)[select.process(train, method="CorrSF", disc.method = "MDL", max.no.features = prefer_no_features)])
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  
  formulas[["AUC_MDLSMOTE"]] = ks.miR.formula(colnames(train_smoted)[select.process(train_smoted, method="auc", disc.method = "MDL", max.no.features = prefer_no_features)])
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  
  
  formulas[["SU_MDLSMOTE"]] = ks.miR.formula(colnames(train_smoted)[select.process(train_smoted, method="symmetrical.uncertainty", disc.method = "MDL", max.no.features = prefer_no_features)])
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  formulas[["CorrSF_MDLSMOTE"]] = ks.miR.formula(colnames(train_smoted)[select.process(train_smoted, method="CorrSF", disc.method = "MDL", max.no.features = prefer_no_features)])
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  formulas[["AUC_MDL_sig"]] = ks.miR.formula(colnames(train_sig)[select.process(train_sig, method="auc", disc.method = "MDL", max.no.features = prefer_no_features)])
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  formulas[["SU_MDL_sig"]] = ks.miR.formula(colnames(train_sig)[select.process(train_sig, method="symmetrical.uncertainty", disc.method = "MDL", max.no.features = prefer_no_features)])
  
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  formulas[["CorrSF_MDL_sig"]] = ks.miR.formula(colnames(train_sig)[select.process(train_sig, method="CorrSF", disc.method = "MDL", max.no.features = prefer_no_features)])
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  formulas[["AUC_MDLSMOTE_sig"]] = ks.miR.formula(colnames(train_sig_smoted)[select.process(train_sig_smoted, method="auc", disc.method = "MDL", max.no.features = prefer_no_features)])
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  formulas[["SU_MDLSMOTE_sig"]] = ks.miR.formula(colnames(train_sig_smoted)[select.process(train_sig_smoted, method="symmetrical.uncertainty", disc.method = "MDL", max.no.features = prefer_no_features)])
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  formulas[["CorrSF_MDLSMOTE_sig"]] = ks.miR.formula(colnames(train_sig_smoted)[select.process(train_sig_smoted, method="CorrSF", disc.method = "MDL", max.no.features = prefer_no_features)])
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
 save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  # 9. bounceR - genetic
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  cat("\n\nStarting bounceR..")
  library(bounceR)
  mrmr <- bounceR::featureSelection(data = train,
                           target = "Class",
                           max_time = "15 mins",
                           selection = selectionControl(n_rounds = NULL,
                                                        n_mods = NULL,
                                                        p = prefer_no_features,
                                                        penalty = 0.5,
                                                        reward = 0.2),
                           bootstrap = "regular",
                           early_stopping = "aic",
                           boosting = boostingControl(mstop = 100, nu = 0.1),
                           cores = parallel::detectCores()-1)
  formulas[["bounceR-full"]] = mrmr@opt_formula
  formulas[["bounceR-stability"]] = ks.miR.formula(as.character(mrmr@stability[1:prefer_no_features,] %>% pull('feature')))
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
 save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  mrmr <- featureSelection(data = train_smoted,
                           target = "Class",
                           max_time = "15 mins",
                           selection = selectionControl(n_rounds = NULL,
                                                        n_mods = NULL,
                                                        p = prefer_no_features,
                                                        penalty = 0.5,
                                                        reward = 0.2),
                           bootstrap = "regular",
                           early_stopping = "aic",
                           boosting = boostingControl(mstop = 100, nu = 0.1),
                           cores = parallel::detectCores()-1)
  formulas[["bounceR-full_SMOTE"]] = mrmr@opt_formula
  formulas[["bounceR-stability_SMOTE"]] = ks.miR.formula(as.character(mrmr@stability[1:prefer_no_features,] %>% pull('feature')))
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
 save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  
  mrmr <- featureSelection(data = train_sig,
                           target = "Class",
                           max_time = "15 mins",
                           selection = selectionControl(n_rounds = NULL,
                                                        n_mods = NULL,
                                                        p = prefer_no_features,
                                                        penalty = 0.5,
                                                        reward = 0.2),
                           bootstrap = "regular",
                           early_stopping = "aic",
                           boosting = boostingControl(mstop = 100, nu = 0.1),
                           cores = parallel::detectCores()-1)
  formulas[["bounceR-full_SIG"]] = mrmr@opt_formula
  formulas[["bounceR-stability_SMOTE"]] = ks.miR.formula(as.character(mrmr@stability[1:prefer_no_features,] %>% pull('feature')))
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
 save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  
  mrmr <- featureSelection(data = train_sig_smoted,
                           target = "Class",
                           max_time = "15 mins",
                           selection = selectionControl(n_rounds = NULL,
                                                        n_mods = NULL,
                                                        p = prefer_no_features,
                                                        penalty = 0.5,
                                                        reward = 0.2),
                           bootstrap = "regular",
                           early_stopping = "aic",
                           boosting = boostingControl(mstop = 100, nu = 0.1),
                           cores = parallel::detectCores()-1)
  formulas[["bounceR-full_SIGSMOTE"]] = mrmr@opt_formula
  formulas[["bounceR-stability_SIGSMOTE"]] = ks.miR.formula(as.character(mrmr@stability[1:prefer_no_features,] %>% pull('feature')))
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
 save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  
  # 10. RFE RandomForest
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  cat("\n\nStarting RFE RF")
  ctrl <- rfeControl(functions =rfFuncs,
                     method = "cv", number = 10,
                     saveDetails = TRUE,
                     allowParallel = TRUE,
                     returnResamp = "all",
                     verbose = T)
  
  rfProfile <- rfe(trainx, train$Class,
                   sizes = 3:11,
                   rfeControl = ctrl)
  plot(rfProfile, type=c("g", "o"))
  formulas[["RandomForestRFE"]] = ks.miR.formula(predictors(rfProfile))
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
 save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  rfProfile <- rfe(trainx_smoted, train_smoted$Class,
                   sizes = 3:11,
                   rfeControl = ctrl)
  plot(rfProfile, type=c("g", "o"))
  formulas[["RandomForestRFESMOTE"]] = ks.miR.formula(predictors(rfProfile))
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
 save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  rfProfile <- rfe(trainx_sig, train_sig$Class,
                   sizes = 3:11,
                   rfeControl = ctrl)
  plot(rfProfile, type=c("g", "o"))
  formulas[["RandomForestRFE_sig"]] = ks.miR.formula(predictors(rfProfile))
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
 save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  rfProfile <- rfe(trainx_sig_smoted, train_sig_smoted$Class,
                   sizes = 3:11,
                   rfeControl = ctrl)
  plot(rfProfile, type=c("g", "o"))
  formulas[["RandomForestRFESMOTE_sig"]] = ks.miR.formula(predictors(rfProfile))
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
 save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  
  
  # 11. Genetic 
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  cat("\n\nStarting Genetic")
  ga_ctrl <- gafsControl(functions = rfGA, method = "repeatedcv", number=10, repeats=5, allowParallel=T)
  rf_ga <- gafs(x = trainx, y = train$Class,
                iters = max_iterations,
                gafsControl = ga_ctrl)
  plot(rf_ga) + theme_bw()
  print(rf_ga)
  formulas[["GeneticAlgorithmRF"]] = ks.miR.formula(rf_ga$ga$final)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
 save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  ga_ctrl <- gafsControl(functions = rfGA, method = "repeatedcv", number=10, repeats=5, allowParallel=T)
  rf_ga <- gafs(x = trainx_smoted, y = train_smoted$Class,
                iters = max_iterations,
                gafsControl = ga_ctrl)
  plot(rf_ga) + theme_bw()
  print(rf_ga)
  formulas[["GeneticAlgorithmRFSMOTE"]] = ks.miR.formula(rf_ga$ga$final)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
 save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  ga_ctrl <- gafsControl(functions = rfGA, method = "repeatedcv", number=10, repeats=5, allowParallel=T)
  rf_ga <- gafs(x = trainx_sig, y = train_sig$Class,
                iters = max_iterations,
                gafsControl = ga_ctrl)
  plot(rf_ga) + theme_bw()
  print(rf_ga)
  formulas[["GeneticAlgorithmRF_sig"]] = ks.miR.formula(rf_ga$ga$final)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
 save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  ga_ctrl <- gafsControl(functions = rfGA, method = "repeatedcv", number=10, repeats=5, allowParallel=T)
  rf_ga <- gafs(x = trainx_sig_smoted, y = train_sig_smoted$Class,
                iters = max_iterations,
                gafsControl = ga_ctrl)
  plot(rf_ga) + theme_bw()
  print(rf_ga)
  formulas[["GeneticAlgorithmRFSMOTE_sig"]] = ks.miR.formula(rf_ga$ga$final)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
 save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  
  # 12. SimulatedAnealing
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  cat("\n\nStarting SimulatedAnealing")
  sa_ctrl <- safsControl(functions = rfSA,
                         method = "repeatedcv",
                         number=5, repeats=10, allowParallel=T,
                         improve = 50)
  
  rf_sa <- safs(x = trainx, y = train$Class,
                iters = max_iterations,
                safsControl = sa_ctrl)
  print(rf_sa)
  plot(rf_sa) + theme_bw()
  formulas[["SimulatedAnnealingRF"]] = ks.miR.formula(rf_sa$sa$final)
  
  
  
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  rf_sa <- safs(x = trainx_smoted, y = train_smoted$Class,
                iters = max_iterations,
                safsControl = sa_ctrl)
  print(rf_sa)
  plot(rf_sa) + theme_bw()
  formulas[["SimulatedAnnealingRFSMOTE"]] = ks.miR.formula(rf_sa$sa$final)
  
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  rf_sa <- safs(x = trainx_sig, y = train_sig$Class,
                iters = max_iterations,
                safsControl = sa_ctrl)
  print(rf_sa)
  plot(rf_sa) + theme_bw()
  formulas[["SimulatedAnnealingRF_sig"]] = ks.miR.formula(rf_sa$sa$final)
  
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  rf_sa <- safs(x = trainx_sig_smoted, y = train_sig_smoted$Class,
                iters = max_iterations,
                safsControl = sa_ctrl)
  print(rf_sa)
  plot(rf_sa) + theme_bw()
  formulas[["SimulatedAnnealingRFSMOTE_sig"]] = ks.miR.formula(rf_sa$sa$final) 
  
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  
  
  # 14. Boruta (https://www.jstatsoft.org/article/view/v036i11)
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  library(Boruta)
  cat("\n\nStarting Boruta")
  bor = Boruta(trainx, train$Class)
  formulas[["Boruta"]] = ks.miR.formula(names(bor$finalDecision)[as.character(bor$finalDecision) == "Confirmed"])
  
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  bor = Boruta(trainx_smoted, train_smoted$Class)
  formulas[["BorutaSMOTE"]] = ks.miR.formula(names(bor$finalDecision)[as.character(bor$finalDecision) == "Confirmed"])
  
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  
  # 15. spFSR - feature selection and ranking by simultaneous perturbation stochastic approximation
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  cat("\n\nStarting spFSR")
  library(spFSR)
  knnWrapper    <- makeLearner("classif.knn", k = 5) 
  classifTask   <- makeClassifTask(data = train, target = "Class")
  perf.measure  <- acc
  spsaMod <- spFeatureSelection(
    task = classifTask,
    wrapper = knnWrapper,
    measure = perf.measure ,
    num.features.selected = prefer_no_features,
    iters.max = max_iterations,
    num.cores = detectCores() - 1)
  formulas[["spFSR"]] = ks.miR.formula(spsaMod$features)
  
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  classifTask   <- makeClassifTask(data = train_smoted, target = "Class")
  perf.measure  <- acc
  spsaMod <- spFeatureSelection(
    task = classifTask,
    wrapper = knnWrapper,
    measure = perf.measure ,
    num.features.selected = prefer_no_features,
    iters.max = max_iterations,
    num.cores = detectCores() - 2)
  formulas[["spFSRSMOTE"]] = ks.miR.formula(spsaMod$features)
  
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  # varSelRF
  cat("\n\nStarting varSelRF")
  library(varSelRF)
  var.sel <- varSelRF(trainx, train$Class, ntree = 500, ntreeIterat = max_iterations, vars.drop.frac = 0.05, whole.range = T, keep.forest = T)
  formulas[["varSelRF"]] = ks.miR.formula(var.sel$selected.vars)
  
  var.sel <- varSelRF(trainx_smoted, train_smoted$Class, ntree = 500, ntreeIterat = max_iterations, vars.drop.frac = 0.05, whole.range = T, keep.forest = T)
  formulas[["varSelRFSMOTE"]] = ks.miR.formula(var.sel$selected.vars)
  
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  # 13. WxNet (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6642261/)
  cat("\nStarting WxNet")
  library(plyr)
  library(dplyr)
  library(reticulate)
  library(tidyverse)
  library(data.table)
  library(DMwR)
  

  
  # Set the path to the Python executale file
  #use_python("/anaconda3/bin/python", required = T)
  
  
  conda_list()
  use_condaenv("tensorflow", required = T)
  py_config()
  
  dane = ks.wczytajmix(replace_smote = F); train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]
  
  #train_org = train
  #train = SMOTE(Class ~ ., data = train_org, perc.over = 10000,perc.under=100)
  
  # Przygotowanie
  trainx = dplyr::select(train, starts_with("hsa"))
  testx = dplyr::select(test, starts_with("hsa"))
  trainx = t(trainx)
  testx = t(testx)
  ids = paste0("train",1:ncol(trainx))
  ids2 = paste0("test",1:ncol(testx))
  colnames(trainx) = gsub("-", "", ids)                     
  colnames(testx) = gsub("-", "", ids2) 
  
  traindata = cbind(data.frame(fnames = rownames(trainx), trainx))
  fwrite(traindata, paste0(code_path, "wx/DearWXpub/src/train-data.csv"), row.names = F, quote = F)                  
  testdata = cbind(data.frame(fnames = rownames(testx), testx))
  fwrite(testdata, paste0(code_path, "wx/DearWXpub/src/test-data.csv"), row.names = F, quote = F)     
  
  trainanno = data.frame(id = colnames(trainx), label = train$Class)
  fwrite(trainanno, paste0(code_path, "wx/DearWXpub/src/train-anno.csv"), row.names = F, quote = F)  
  testanno = data.frame(id = colnames(testx), label = test$Class)
  fwrite(testanno, paste0(code_path, "wx/DearWXpub/src/test-anno.csv"), row.names = F, quote = F)
  
  
  # Wywolanie WX
  
  out <- tryCatch(
    {
      setwd(paste0(code_path,"wx/DearWXpub/src/"))
      system(paste0(conda_path," activate tensorflow"))
      py_run_file("wx_konsta.py", local = T)
      
      np <- import("numpy")
      w = np$load("wyniki.npy",allow_pickle=T)
      formulas[["Wx"]] = ks.miR.formula(w)
    },
    error=function(cond) {
      message("ERROR:")
      message(cond)
      # Choose a return value in case of error
    },
    warning=function(cond) {
      message("WARNING:")
      message(cond)
    },
    finally={
      setwd(wd)
    }
  )
    
 
  
  
  
  # Wx with SMOTE
  dane = ks.wczytajmix(replace_smote = T); train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]
  
  # Przygotowanie
  trainx = dplyr::select(train, starts_with("hsa"))
  testx = dplyr::select(test, starts_with("hsa"))
  trainx = t(trainx)
  testx = t(testx)
  ids = paste0("train",1:ncol(trainx))
  ids2 = paste0("test",1:ncol(testx))
  colnames(trainx) = gsub("-", "", ids)                     
  colnames(testx) = gsub("-", "", ids2) 
  
  traindata = cbind(data.frame(fnames = rownames(trainx), trainx))
  fwrite(traindata, paste0(code_path, "wx/DearWXpub/src/train-data.csv"), row.names = F, quote = F)                  
  testdata = cbind(data.frame(fnames = rownames(testx), testx))
  fwrite(testdata, paste0(code_path, "wx/DearWXpub/src/test-data.csv"), row.names = F, quote = F)     
  
  trainanno = data.frame(id = colnames(trainx), label = train$Class)
  fwrite(trainanno, paste0(code_path, "wx/DearWXpub/src/train-anno.csv"), row.names = F, quote = F)  
  testanno = data.frame(id = colnames(testx), label = test$Class)
  fwrite(testanno, paste0(code_path, "wx/DearWXpub/src/test-anno.csv"), row.names = F, quote = F)
  
  
  # Wywolanie WX
  # setwd(paste0(code_path,"wx/DearWXpub/src/"))
  # py_run_file("wx_konsta.py", local = T)
  # 
  # np <- import("numpy")
  # w = np$load("wyniki.npy",allow_pickle=T)
  # print(w)
  # setwd(wd)
  # formulas[["WxSMOTE"]] = ks.miR.formula(w)
  out <- tryCatch(
    {
      setwd(paste0(code_path,"wx/DearWXpub/src/"))
      py_run_file("wx_konsta.py", local = T)
      
      np <- import("numpy")
      w = np$load("wyniki.npy",allow_pickle=T)
      formulas[["WxSMOTE"]] = ks.miR.formula(w)
    },
    error=function(cond) {
      message("ERROR:")
      message(cond)
      # Choose a return value in case of error
    },
    warning=function(cond) {
      message("WARNING:")
      message(cond)
    },
    finally={
      setwd(wd)
    }
  )
  
  dane = ks.wczytajmix(replace_smote = F); train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]
  
  # Przygotowanie
  trainx = dplyr::select(train, starts_with("hsa"))
  testx = dplyr::select(test, starts_with("hsa"))
  trainx = t(trainx)
  testx = t(testx)
  ids = paste0("train",1:ncol(trainx))
  ids2 = paste0("test",1:ncol(testx))
  colnames(trainx) = gsub("-", "", ids)                     
  colnames(testx) = gsub("-", "", ids2) 
  
  traindata = cbind(data.frame(fnames = rownames(trainx), trainx))
  fwrite(traindata, paste0(code_path, "wx/DearWXpub/src/train-data.csv"), row.names = F, quote = F)                  
  testdata = cbind(data.frame(fnames = rownames(testx), testx))
  fwrite(testdata, paste0(code_path, "wx/DearWXpub/src/test-data.csv"), row.names = F, quote = F)     
  
  trainanno = data.frame(id = colnames(trainx), label = train$Class)
  fwrite(trainanno, paste0(code_path, "wx/DearWXpub/src/train-anno.csv"), row.names = F, quote = F)  
  testanno = data.frame(id = colnames(testx), label = test$Class)
  fwrite(testanno, paste0(code_path, "wx/DearWXpub/src/test-anno.csv"), row.names = F, quote = F)
  
  
  # Wywolanie WX
  # setwd(paste0(code_path,"wx/DearWXpub/src/"))
  # py_run_file("wx_konsta_z.py", local = T)
  # 
  # np <- import("numpy")
  # w = np$load("wyniki.npy",allow_pickle=T)
  # print(w)
  # setwd(wd)
  # formulas[["Wx_Zscore"]] = ks.miR.formula(w)
  out <- tryCatch(
    {
      setwd(paste0(code_path,"wx/DearWXpub/src/"))
      py_run_file("wx_konsta_z.py", local = T)
      
      np <- import("numpy")
      w = np$load("wyniki.npy",allow_pickle=T)
      formulas[["Wx_Zscore"]] = ks.miR.formula(w)
    },
    error=function(cond) {
      message("ERROR:")
      message(cond)
      # Choose a return value in case of error
    },
    warning=function(cond) {
      message("WARNING:")
      message(cond)
    },
    finally={
      setwd(wd)
    }
  )
  
  dane = ks.wczytajmix(replace_smote = T); train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]
  
  # Przygotowanie
  trainx = dplyr::select(train, starts_with("hsa"))
  testx = dplyr::select(test, starts_with("hsa"))
  trainx = t(trainx)
  testx = t(testx)
  ids = paste0("train",1:ncol(trainx))
  ids2 = paste0("test",1:ncol(testx))
  colnames(trainx) = gsub("-", "", ids)                     
  colnames(testx) = gsub("-", "", ids2) 
  
  traindata = cbind(data.frame(fnames = rownames(trainx), trainx))
  fwrite(traindata, paste0(code_path, "wx/DearWXpub/src/train-data.csv"), row.names = F, quote = F)                  
  testdata = cbind(data.frame(fnames = rownames(testx), testx))
  fwrite(testdata, paste0(code_path, "wx/DearWXpub/src/test-data.csv"), row.names = F, quote = F)     
  
  trainanno = data.frame(id = colnames(trainx), label = train$Class)
  fwrite(trainanno, paste0(code_path, "wx/DearWXpub/src/train-anno.csv"), row.names = F, quote = F)  
  testanno = data.frame(id = colnames(testx), label = test$Class)
  fwrite(testanno, paste0(code_path, "wx/DearWXpub/src/test-anno.csv"), row.names = F, quote = F)
  
  
  # Wywolanie WX
  # setwd(paste0(code_path,"wx/DearWXpub/src/"))
  # py_run_file("wx_konsta_z.py", local = T)
  # 
  # np <- import("numpy")
  # w = np$load("wyniki.npy",allow_pickle=T)
  # print(w)
  # setwd(wd)
  # formulas[["Wx_ZscoreSMOTE"]] = ks.miR.formula(w)
  out <- tryCatch(
    {
      setwd(paste0(code_path,"wx/DearWXpub/src/"))
      py_run_file("wx_konsta_z.py", local = T)
      
      np <- import("numpy")
      w = np$load("wyniki.npy",allow_pickle=T)
      formulas[["Wx_ZscoreSMOTE"]] = ks.miR.formula(w)
    },
    error=function(cond) {
      message("ERROR:")
      message(cond)
      # Choose a return value in case of error
    },
    warning=function(cond) {
      message("WARNING:")
      message(cond)
    },
    finally={
      setwd(wd)
    }
  )
  
  
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  # My.stepwise
  n = n + 1
  if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); 
    cat("\n\nStarting My.stepwise")
    start_time <- Sys.time()
    
    dane = ks.wczytajmix(replace_smote = F); train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]
    
    library(My.stepwise)
    temp = capture.output(My.stepwise.glm(Y = "Class", colnames(trainx), data = train, sle = 0.05, sls = 0.05, myfamily = "binomial"))
    # temp2 = temp[length(temp)-1]
    # temp3 = temp[length(temp)-3]
    # temp4 = temp[length(temp)-5]
    temp5 = paste(temp[(length(temp)-12):length(temp)], collapse = " ")
    wybrane = FALSE
    for(i in 1:length(colnames(trainx)))
    {
      temp6 = colnames(trainx)[i]
      wybrane[i] = grepl(temp6, temp5)
    }
    formulas[["Mystepwise_glm_binomial"]] = ks.miR.formula(colnames(trainx)[wybrane])
    
    wyniki = ks.miRNA.de(trainx, train$Class)
    istotne = filter(wyniki, `p-value BH` <= 0.05) %>% arrange(`p-value BH`)
    
    temp = capture.output(My.stepwise.glm(Y = "Class", as.character(istotne$miR), data = train, sle = 0.05, sls = 0.05, myfamily = "binomial"))
    # temp2 = temp[length(temp)-1]
    # temp3 = temp[length(temp)-3]
    # temp4 = temp[length(temp)-5]
    temp5 = paste(temp[(length(temp)-12):length(temp)], collapse = " ")
    wybrane = FALSE
    for(i in 1:length(colnames(trainx)))
    {
      temp6 = colnames(trainx)[i]
      wybrane[i] = grepl(temp6, temp5)
    }
    formulas[["Mystepwise_sig_glm_binomial"]] = ks.miR.formula(colnames(trainx)[wybrane])
    
    end_time <- Sys.time()
    saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS"))
    saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
    save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  n = n + 1
  if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); 
    cat("\n\nStarting My.stepwise SMOTE")
    start_time <- Sys.time()
    
    dane = ks.wczytajmix(replace_smote = F); train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]
    
    library(My.stepwise)
    temp = capture.output(My.stepwise.glm(Y = "Class", colnames(trainx_smoted), data = train_smoted, sle = 0.05, sls = 0.05, myfamily = "binomial"))
    # temp2 = temp[length(temp)-1]
    # temp3 = temp[length(temp)-3]
    # temp4 = temp[length(temp)-5]
    temp5 = paste(temp[(length(temp)-12):length(temp)], collapse = " ")
    wybrane = FALSE
    for(i in 1:length(colnames(trainx_smoted)))
    {
      temp6 = colnames(trainx_smoted)[i]
      wybrane[i] = grepl(temp6, temp5)
    }
    formulas[["Mystepwise_glm_binomialSMOTE"]] = ks.miR.formula(colnames(trainx_smoted)[wybrane])
    
    wyniki = ks.miRNA.de(trainx, train$Class)
    istotne = filter(wyniki, `p-value BH` <= 0.05) %>% arrange(`p-value BH`)
    
    temp = capture.output(My.stepwise.glm(Y = "Class", as.character(istotne$miR), data = train_smoted, sle = 0.05, sls = 0.05, myfamily = "binomial"))
    # temp2 = temp[length(temp)-1]
    # temp3 = temp[length(temp)-3]
    # temp4 = temp[length(temp)-5]
    temp5 = paste(temp[(length(temp)-12):length(temp)], collapse = " ")
    wybrane = FALSE
    for(i in 1:length(colnames(trainx_smoted)))
    {
      temp6 = colnames(trainx_smoted)[i]
      wybrane[i] = grepl(temp6, temp5)
    }
    formulas[["Mystepwise_sig_glm_binomialSMOTE"]] = ks.miR.formula(colnames(trainx_smoted)[wybrane])
    
    end_time <- Sys.time()
    saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS"))
    saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
    save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  n = n + 1
  if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); 
    cat("\n\nStarting stepAIC")
    start_time <- Sys.time()
    
    dane = ks.wczytajmix(replace_smote = F); train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]
    
    temp = glm(Class ~ ., data = train, family = "binomial")
    temp2 = stepAIC(temp)
    
    formulas[["stepAIC"]] = temp2$formula
    
    wyniki = ks.miRNA.de(trainx, train$Class)
    istotne = filter(wyniki, `p-value BH` <= 0.05) %>% arrange(`p-value BH`)
    
    train.sig = dplyr::select(train, as.character(istotne$miR), Class)
    
    temp = glm(Class ~ ., data = train.sig, family = "binomial")
    temp2 = stepAIC(temp)
    
    formulas[["stepAICsig"]] = temp2$formula
    
    end_time <- Sys.time()
    saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS"))
    saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
    save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  n = n + 1
  if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); 
    cat("\n\nStarting stepAIC SMOTE")
    start_time <- Sys.time()
    
    dane = ks.wczytajmix(replace_smote = F); train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]
    
    temp = glm(Class ~ ., data = train_smoted, family = "binomial")
    temp2 = stepAIC(temp)
    
    formulas[["stepAIC_SMOTE"]] = temp2$formula
    
    wyniki = ks.miRNA.de(trainx, train$Class)
    istotne = filter(wyniki, `p-value BH` <= 0.05) %>% arrange(`p-value BH`)
    
    train.sig = dplyr::select(train_smoted, as.character(istotne$miR), Class)
    
    temp = glm(Class ~ ., data = train.sig, family = "binomial")
    temp2 = stepAIC(temp)
    
    formulas[["stepAICsig_SMOTE"]] = temp2$formula
    
    end_time <- Sys.time()
    saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS"))
    saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
    save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  n = n + 1
  if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); 
    cat("\n\nStarting MK method (iterated RFE)")
    start_time <- Sys.time()
    
    dane = ks.wczytajmix(replace_smote = F); train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]
    
    selectedMirsCV <- mk.iteratedRFE(trainSet = train, useCV = T, classLab = 'Class', checkNFeatures = prefer_no_features)$topFeaturesPerN[[prefer_no_features]]
    selectedMirsTest <- mk.iteratedRFE(trainSet = train, testSet = test, classLab = 'Class', checkNFeatures = prefer_no_features)$topFeaturesPerN[[prefer_no_features]]
    
    formulas[["iteratedRFECV"]] = ks.miR.formula(selectedMirsCV$topFeaturesPerN[[prefer_no_features]])
    formulas[["iteratedRFETest"]] = ks.miR.formula(selectedMirsTest$topFeaturesPerN[[prefer_no_features]])
    
    end_time <- Sys.time()
    saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS"))
    saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
    save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  n = n + 1
  if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); 
    cat("\n\nStarting MK method (iterated RFE) with SMOTE")
    start_time <- Sys.time()
    
    dane = ks.wczytajmix(replace_smote = F); train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]
    
    selectedMirsCV <- mk.iteratedRFE(trainSet = train_smoted, useCV = T, classLab = 'Class', checkNFeatures = prefer_no_features)$topFeaturesPerN[[prefer_no_features]]
    selectedMirsTest <- mk.iteratedRFE(trainSet = train_smoted, testSet = test, classLab = 'Class', checkNFeatures = prefer_no_features)$topFeaturesPerN[[prefer_no_features]]
    
    formulas[["iteratedRFECV"]] = ks.miR.formula(selectedMirsCV$topFeaturesPerN[[prefer_no_features]])
    formulas[["iteratedRFETest"]] = ks.miR.formula(selectedMirsTest$topFeaturesPerN[[prefer_no_features]])
    
    end_time <- Sys.time()
    saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS"))
    saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
    save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  
  # End
  cat("\n\nEnding...")
  stopCluster(cl)
  saveRDS(formulas, paste0("temp/featureselection_formulas",stamp,".RDS"))
  saveRDS(formulas, paste0("featureselection_formulas.RDS"))
  nazwa = paste0("temp/featureselection_data",stamp,".rda")
  save(list = ls(), file = nazwa)
  dev.off() 
  sink()
  print(formulas)
  setwd(oldwd)
  return(formulas)
}

mk.iteratedRFE <- function(trainSet, testSet = NULL, initFeatures = colnames(trainSet), classLab, checkNFeatures = 25, votingIterations = 100000, useCV = F, nfolds = 10, initRandomState = 42 ) {
  
  set.seed(initRandomState)
  
  #prepare output data structures
  resAcc <- c(rep(0, checkNFeatures))
  resVotes <- data.frame(matrix(0, nrow = length(initFeatures), ncol = checkNFeatures), row.names = initFeatures)
  for(i in 1:checkNFeatures) colnames(resVotes)[i] <- toString(i)
  resTop <- list()
  
  
  for (i in 1:votingIterations) {
    
    if(useCV == F) {
      params <- rfeControl(functions = rfFuncs, saveDetails = T)
      iter <- rfeIter(x = trainSet[, initFeatures], y = as.factor(trainSet[, classLab]), testX = testSet[, initFeatures], testY = as.factor(testSet[, classLab]), sizes = 1:checkNFeatures,
                      metric = "Accuracy", rfeControl = params)
      
      for(j in 1:checkNFeatures) {
        tmp <- iter$pred[iter$pred$Variables == j, ]
        
        acc <- length(which(tmp$pred == tmp$obs)) / nrow(tmp) #calculate and add accuracy
        resAcc[j] <- resAcc[j] + acc
        
        selected <- iter$finalVariables[[j+1]]$var
        numb <- iter$finalVariables[[j+1 ]]$Variables[1]
        
        resVotes[selected, numb] <- resVotes[selected, numb] + 1
      }
      
      
    }
    else {
      
      seeds <- vector(mode = "list", length = nfolds + 1) # add random seeds for cross validation
      for(i in 1:nfolds) seeds[[i]] <- sample.int(1000000000, checkNFeatures + 1)
      seeds[nfolds + 1] <- sample.int(1000000000, 1)
      
      params <- rfeControl(functions = rfFuncs, number = nfolds, saveDetails = T)
      iter <- rfe(x = trainSet[, initFeatures], y = as.factor(trainSet[, classLab]), sizes = 1:checkNFeatures, rfeControl = params)
      
      for(j in 1:checkNFeatures) {
        tmp <- iter$variables[iter$variables$Variables == j, ]
        for(k in tmp$var) resVotes[k, j] <- resVotes[k, j] + 1 # increase a voting score for each fold
        
        resAcc[j] <- resAcc[j] + iter$results[iter$results$Variables == j, "Accuracy"]
        
      }
    }
  }
  
  resAcc <- resAcc / votingIterations #make average accuracy
  
  for(i in 1:ncol(resVotes)) resTop[[i]] <- rownames(resVotes[order(-resVotes[, i])[1:i], ])
  
  returning <- list(data.frame(resAcc), resVotes, resTop)
  names(returning) <- c("accuracyPerNFeatures", "votesPerN", "topFeaturesPerN")
  return(returning)
  
}

ks.benchmark = function(wd = getwd(), search_iters = 2000, keras_epochs = 5000, keras_threads = floor(parallel::detectCores()/2), search_iters_mxnet = 5000,
                        cores = detectCores()-1, input_formulas = readRDS("featureselection_formulas_final.RDS"), 
                        output_file = "benchmark.csv", mxnet = F, gpu = F, 
                        #algorithms = c("nnet","svmRadial", "svmLinear","rf","C5.0","mlp", "mlpML","xgbTree"), 
                        algorithms = c("mlpKerasDropout", "mlpKerasDecay", "mlpML","xgbTree","svmRadial", "svmLinear","rf","C5.0"), 
                        holdout = T, stamp = as.character(as.numeric(Sys.time()))) {
  library(plyr)
  library(dplyr)
  library(edgeR)
  library(epiDisplay)
  library(rsq)
  library(MASS)
  library(Biocomb)
  library(caret)
  library(dplyr)
  library(epiDisplay)
  library(pROC)
  library(ggplot2)
  library(DMwR)
  library(stringr)
  library(psych)
  library(C50)
  library(randomForest)
  library(nnet)
  library(reticulate)
  library(stargazer)

  if(!dir.exists("temp")) { dir.create("temp") }
  
  use_condaenv("tensorflow")
  zz <- file(paste0("temp/benchmark",stamp,".log"), open = "wt")
  pdf(paste0("temp/benchmark",stamp,".pdf"))
  sink(zz)
  sink(zz, type = "message")
  library(doParallel)
  
  
  oldwd = getwd()
  setwd(wd)
  
  formulas = input_formulas
  
  wybrane = list()
  ile_miRNA = numeric()
  for (i in 1:length(formulas)) {
    wybrane[[i]] = all.vars(as.formula(formulas[[i]]))[-1]
    ile_miRNA[i] = length(all.vars(as.formula(formulas[[i]])))-1
  }
  
  hist(ile_miRNA, breaks = 150)
  psych::describe(ile_miRNA)
  
  # Wywalamy formuły z więcej niż 15 miRNA
  #formulas_old = formulas
  #formulas = formulas[which(ile_miRNA <= 15)]
  #print(formulas)
  
  dane = ks.wczytajmix(wd = wd, replace_smote = F); train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]
  
  if(!dir.exists("models")) { dir.create("models") }
  
  wyniki = data.frame(method = names(formulas))
  wyniki$SMOTE = ifelse(grepl("SMOTE", names(formulas)),"Yes","No")

  
  if(mxnet == T) { 
    library(mxnet)
    #gpu = system("nvidia-smi", intern = T)
    if (gpu == T) {
    mxctx = mx.gpu() 
    print("MXNET is using GPU :)")
    } else {
      mxctx = mx.cpu()
      print("MXNET is using CPU :(")
    }
    
    
    algorytmy = c("glm", "mxnetAdam","mxnet", algorithms)
  } else { 
    algorytmy = c("glm",algorithms)
    }
  
  

  
  
  for (ii in 1:length(algorytmy)) {
    algorytm = algorytmy[ii]
    print(algorytm)
    
    if(grepl("Keras",algorytm)) {
      cl <- makePSOCKcluster(keras_threads)
      registerDoParallel(cl)
    } else {
      cl <- makePSOCKcluster(cores-1)
      registerDoParallel(cl)
    }
    
    for (i in 1:length(formulas)) {
      print(paste0("Testing formula: ", as.character(formulas[[i]])))
      
      wyniki$miRy[i] = as.character(formulas[[i]])
      temptrain = train 
      if (wyniki$SMOTE[i] == "Yes") { temptrain = train_smoted } 
      
      # Hold-out czy CV?
      temptrainold = temptrain
      if(holdout == T) {
        #fit_on = list(rs1 = 1:nrow(temptrain), rs2 = 1:nrow(temptrain))
        #pred_on = list(rs1 = (nrow(temptrain)+1):((nrow(temptrain)+1)+nrow(test)), rs2 = ((nrow(temptrain)+1)+nrow(test)+1):((nrow(temptrain)+1)+nrow(test)+1+nrow(valid)))
        #temptrain = rbind.fill(temptrain,test,valid)
        fit_on = list(rs1 = 1:nrow(temptrain))
        pred_on = list(rs1 = (nrow(temptrain)+1):((nrow(temptrain))+nrow(test)))
        temptrain = rbind.fill(temptrain,test)
      }
      
      
      # wyniki2 = tryCatch({
      if(algorytm == "mxnet") {
          hyperparameters = expand.grid(layer1 = seq(2,14,1), layer2 = c(0,2,4), layer3 = c(0), activation = c('relu', 'sigmoid', 'tanh', 'softrelu'),
                                        dropout = c(0,0.05),learning.rate= c(0.00001, 0.01), momentum = c(0, 0.8, 0.99))
          train_control <- trainControl(method="cv", repeats=5, number = 10, classProbs = TRUE,verboseIter = TRUE, summaryFunction = twoClassSummary, savePredictions = TRUE)
          if(holdout == T) { train_control <- trainControl(method="cv", index= fit_on, indexOut = pred_on, indexFinal = fit_on[[1]], verboseIter = TRUE, 
                                                           classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = TRUE) }
          model1 = caret::train(as.formula(formulas[[i]]), ctx = mxctx, optimizer = 'sgd', 
                                #optimizer_params=(('learning_rate',0.1),('lr_scheduler',lr_sch)),
                                #preProc = c("center", "scale"),
                                #epoch.end.callback = mx.callback.early.stop(5,10,NULL, maximize=TRUE, verbose=TRUE),
                                #eval.data = list(data=dplyr::select(test, starts_with("hsa")),label=dplyr::select(test, Class)),
                                epoch.end.callback=mx.callback.early.stop(30, 30),
                                #preProc = c("center", "scale"),
                                num.round = search_iters_mxnet, data=temptrain, trControl=train_control, method=algorytm, tuneGrid = hyperparameters)
          print(model1$finalModel)
      } else if(algorytm == "mxnetAdam") {
        hyperparameters = expand.grid(layer1 = seq(2,14,1), layer2 = c(0,2,4), layer3 = c(0), activation = c('relu', 'sigmoid', 'tanh', 'softrelu'),
                                      dropout = c(0,0.05), beta1=0.9, beta2=0.999, learningrate= c(0.00001, 0.001, 0.01))
        train_control <- trainControl(method="cv", repeats=5, number = 10, classProbs = TRUE,verboseIter = TRUE, summaryFunction = twoClassSummary, savePredictions = TRUE)
        if(holdout == T) { train_control <- trainControl(method="cv", index= fit_on, indexOut = pred_on, indexFinal = fit_on[[1]], verboseIter = TRUE, 
                                                         classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = TRUE) }
        model1 = caret::train(as.formula(formulas[[i]]), ctx = mxctx,
                              #optimizer_params=(('learning_rate',0.1),('lr_scheduler',lr_sch)),
                              #preProc = c("center", "scale"),
                              #epoch.end.callback = mx.callback.early.stop(5,10,NULL, maximize=TRUE, verbose=TRUE),
                              #eval.data = list(data=dplyr::select(test, starts_with("hsa")),label=dplyr::select(test, Class)),
                              epoch.end.callback=mx.callback.early.stop(30, 30),
                              #preProc = c("center", "scale"),
                              num.round = search_iters_mxnet, data=temptrain, trControl=train_control, method=algorytm, tuneGrid = hyperparameters)
        print(model1$finalModel)
      } else if(grepl("Keras",algorytm)) {
        train_control <- trainControl(method="cv", number = 5, search="random", classProbs = TRUE, verboseIter = TRUE,
                                      summaryFunction = twoClassSummary, savePredictions = TRUE, indexFinal = fit_on[[1]])
        #if(holdout == T) { train_control <- trainControl(method="cv", index= fit_on, indexOut = pred_on, indexFinal = fit_on[[1]], verboseIter = TRUE, 
        #                                                 classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = TRUE) }
        model1 = caret::train(as.formula(formulas[[i]]), data=temptrain, trControl=train_control, method=algorytm, tuneLength = search_iters, 
                              epochs = keras_epochs)
        print(model1$finalModel)
      } else if (algorytm == "glm") {
        train_control <- trainControl(method="repeatedcv", repeats=5, number = 10, search="random", classProbs = TRUE, verboseIter = TRUE,
                                      summaryFunction = twoClassSummary, savePredictions = TRUE)
        if(holdout == T) { train_control <- trainControl(method="cv", index= fit_on, indexOut = pred_on, indexFinal = fit_on[[1]], verboseIter = TRUE, 
                                                         classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = TRUE) }
        model1 = caret::train(as.formula(formulas[[i]]), data=temptrain, trControl=train_control, method=algorytm, family="binomial")
        print(model1$finalModel)
          } else {
        train_control <- trainControl(method="repeatedcv", repeats=5, number = 10, search="random", classProbs = TRUE, verboseIter = TRUE,
                                      summaryFunction = twoClassSummary, savePredictions = TRUE)
        if(holdout == T) { train_control <- trainControl(method="cv", index= fit_on, indexOut = pred_on, indexFinal = fit_on[[1]], verboseIter = TRUE, 
                                                         classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = TRUE) }
        model1 = caret::train(as.formula(formulas[[i]]), data=temptrain, trControl=train_control, method=algorytm, tuneLength = search_iters)
        print(model1$finalModel)
        }
      
      
      modelname = as.numeric(Sys.time())
      print(paste0("MODELID: ", modelname))
      saveRDS(model1, paste0("models/",modelname,".RDS"))
      wyniki[i,paste0(algorytm,"_modelname")] = modelname
      
      t1_roc = roc(temptrainold$Class ~ predict(model1, newdata = temptrainold , type = "prob")[,2])
      wyniki[i,paste0(algorytm,"_train_ROCAUC")] = t1_roc$auc
      wyniki[i,paste0(algorytm,"_train_ROCAUC_lower95CI")] = as.character(ci(t1_roc))[1]
      wyniki[i,paste0(algorytm,"_train_ROCAUC_upper95CI")] = as.character(ci(t1_roc))[3]
      
      t1 = caret::confusionMatrix(predict(model1, newdata = temptrainold), as.factor(temptrainold$Class), positive = "Cancer")
      wyniki[i,paste0(algorytm,"_train_Accuracy")] = t1$overall["Accuracy"]
      wyniki[i,paste0(algorytm,"_train_Sensitivity")] = t1$byClass["Sensitivity"]
      wyniki[i,paste0(algorytm,"_train_Specificity")] = t1$byClass["Specificity"]
      
      v1 = caret::confusionMatrix(predict(model1, newdata = test), as.factor(test$Class), positive = "Cancer")
      wyniki[i,paste0(algorytm,"_test_Accuracy")]  = v1$overall["Accuracy"]
      wyniki[i,paste0(algorytm,"_test_Sensitivity")]  = v1$byClass["Sensitivity"]
      wyniki[i,paste0(algorytm,"_test_Specificity")]  = v1$byClass["Specificity"]
      
      v1 = caret::confusionMatrix(predict(model1, newdata = valid), as.factor(valid$Class), positive = "Cancer")
      wyniki[i,paste0(algorytm,"_valid_Accuracy")]  = v1$overall["Accuracy"]
      wyniki[i,paste0(algorytm,"_valid_Sensitivity")]= v1$byClass["Sensitivity"]
      wyniki[i,paste0(algorytm,"_valid_Specificity")] = v1$byClass["Specificity"]
      # return(wyniki)
      # }
      # , warning = function(warning_condition) {
      #   print(paste("WARNING:  ",warning_condition))
      #   return(wyniki)
      # }, error = function(error_condition) {
      #   print(paste("WARNING:  ",error_condition))
      #   wyniki[i,paste0(algorytm,"_train_Accuracy")] = as.character(error_condition)
      #   wyniki[i,paste0(algorytm,"_train_Sensitivity")] = NA
      #   wyniki[i,paste0(algorytm,"_train_Specificity")] = NA
      # 
      #   wyniki[i,paste0(algorytm,"_test_Accuracy")]  = NA
      #   wyniki[i,paste0(algorytm,"_test_Sensitivity")]  = NA
      #   wyniki[i,paste0(algorytm,"_test_Specificity")]  = NA
      #   
      #   wyniki[i,paste0(algorytm,"_valid_Accuracy")]  = NA
      #   wyniki[i,paste0(algorytm,"_valid_Sensitivity")]= NA
      #   wyniki[i,paste0(algorytm,"_valid_Specificity")] = NA
      #   return(wyniki)
      # }, finally={
      #   stargazer(wyniki, type = 'html', out = paste0("temp/wyniki",stamp,".html"), summary = F)
      # })
      stargazer(wyniki, type = 'html', out = paste0("temp/wyniki",stamp,".html"), summary = F)
      write.csv(wyniki,paste0("temp/",output_file))
    }
    
    stopCluster(cl)
    
  }
  
  stargazer(wyniki, type = 'html', out = paste0("temp/wyniki",stamp,".html"), summary = F)
  
  write.csv(wyniki,output_file)
  setwd(oldwd)
  
  sink(type = "message")
  sink()
  dev.off()
  return(wyniki)
}

ks.korektaaliasowmiRNA = function(temp, species = "hsa") {
  library(foreach)
  library(doParallel)
  library(dplyr)
  miRbase_aliasy = fread("ftp://mirbase.org/pub/mirbase/CURRENT/aliases.txt.gz")
  colnames(miRbase_aliasy) = c("MIMAT","Aliasy")
  miRbase_aliasy_hsa = miRbase_aliasy %>% filter(str_detect(Aliasy, paste0(species,"*"))) %>% filter(str_detect(MIMAT, "MIMAT*"))
  
  #setup parallel backend to use many processors
  cores=detectCores()
  cl <- makePSOCKcluster(cores-1) #not to overload your computer
  registerDoParallel(cl)
  
  temp2 = colnames(temp)
  final <- foreach(i=1:length(temp2), .combine=c) %dopar% {
    #for(i in 1:length(temp2)) {
    naz = temp2[i]
    library(data.table)
    library(stringr)
    for (ii in 1:nrow(miRbase_aliasy_hsa)) {
      temp3 = str_split(as.character(miRbase_aliasy_hsa[ii,2]), ";")
      temp4 = temp3[[1]]
      temp4 = temp4[temp4 != ""]
      if(naz %in% temp4) { naz = temp4[length(temp4)] }
    }
    naz
  }
  
  colnames(temp) = final
  stopCluster(cl)
  return(temp)
}

ks.setup = function(install_pkgs = T, create_pyenv = T) {
  # apt install libcairo2-dev libopencv-dev build-essential git ninja-build ccache libopenblas-dev libglu1-mesa-dev freeglut3-dev mesa-common-dev libtool autossh ocl-icd-opencl-dev libfreetype6-dev libssl-dev default-jdk default-jre libcurl4-openssl-dev libxml2-dev
  if(install_pkgs) {
  install.packages("BiocManager")
  BiocManager::install(c("reticulate","devtools","plyr","dplyr","edgeR","epiDisplay","rsq","MASS","Biocomb","caret","dplyr",
                       "pROC","ggplot2","DMwR", "doParallel", "Boruta", "spFSR", "varSelRF", "stringr", "psych", "C50", "randomForest",
                       "foreach","data.table", "ROSE", "deepnet", "gridExtra", "stargazer","gplots"))
  devtools::install_github("STATWORX/bounceR")
  #system("../mxnet/docs/install/install_mxnet_ubuntu_python.sh")
  }
  library(reticulate) 
  library(devtools)
  
  # paczki nie w CRAN
  #install.packages("BiocManager")
  
  py_config()
  
  if(create_pyenv) {
  gpu = system("nvidia-smi", intern = T)
  if (length(gpu) > 1) {
    reticulate::conda_create("tensorflow",c("tensorflow-gpu","keras","mxnet-gpu"))
  } else { 
    reticulate::conda_create("tensorflow",c("tensorflow","keras","mxnet")) 
    }
  }
  
  use_condaenv("tensorflow")
  

}

ks.merge_formulas = function(wd = getwd(), max_miRNAs = 11, add = list()) {
  oldwd = getwd()
  setwd(wd)
  formulas_files = list.files("temp","formulas*", all.files = T, full.names = T)
  formulas = list()
  formulas_names = character()
  for (i in 1:length(formulas_files)) {
    temp = readRDS(formulas_files[i])
    tempn = names(temp)
    formulas = c(formulas, temp)
    formulas_names = c(formulas_names, tempn)
  }
  
  temp = data.frame(name = formulas_names, formula = unlist(as.character(formulas)), stringsAsFactors = F)
  final = temp %>% dplyr::distinct()

  final$formula
  formulas = list()
  all = as.list(final$formula)
  names(all) =  make.names(final$name, unique = T)
  saveRDS(all, "featureselection_formulas_all.RDS")
  for(i in 1:nrow(final)){
    final$ile_miRNA[i] = length(all.vars(as.formula(final$formula[i])))-1
  }
  
  for (i in 1:length(add)) {
    tempdupa = data.frame(name = names(add)[i], formula = paste0("Class ~ ", as.character(ks.miR.formula(add[[i]]))[3]), ile_miRNA = 0, stringsAsFactors = F)
    final = rbind(tempdupa, final)
  }
  
  finalold = final
  final = final %>% filter(ile_miRNA <= max_miRNAs)
  if (("fcsig" %in% final$name) == FALSE) { final = rbind(finalold[which(finalold$name=="fcsig"),], final)
  }
  if (("cfs_sig" %in% final$name) == FALSE) { final = rbind(finalold[which(finalold$name=="cfs_sig"),], final) }
  formulas_final = as.list(final$formula)
  names(formulas_final) = make.names(final$name, unique = T)
  setwd(oldwd)
  saveRDS(formulas_final, "featureselection_formulas_final.RDS")
  #sink()
  #dev.off()
  return(formulas_final)
}

ks.miRNA_signiture_overlap = function(which_formulas = c("sig","cfs"), benchmark_csv = "benchmark1578929876.21765.csv")
{
  benchmark = read.csv(benchmark_csv, stringsAsFactors = F)
  rownames(benchmark) = make.names(benchmark$method, unique = T)
  wybrane = list()
  for (i in 1:length(which_formulas)) {
    ktora_to = match(which_formulas[i], rownames(benchmark))
    temp = as.formula(benchmark$miRy[ktora_to])
    wybrane[[rownames(benchmark)[ktora_to]]] = all.vars(temp)[-1]
  }
  require("gplots")
  temp = venn(wybrane)
  temp
}

ks.get_benchmark = function(benchmark_csv = "benchmark1578929876.21765.csv"){
  benchmark = read.csv(benchmark_csv, stringsAsFactors = F)
  rownames(benchmark) = make.names(benchmark$method, unique = T)
  benchmark$method = rownames(benchmark)
  return(benchmark)
}

ks.get_benchmark_methods = function(benchmark_csv = "benchmark1578929876.21765.csv"){
  library(limma)
  benchmark = read.csv(benchmark_csv, stringsAsFactors = F)
  rownames(benchmark) = make.names(benchmark$method, unique = T)
  benchmark$method = rownames(benchmark)
  temp = dplyr::select(benchmark, ends_with("_valid_Accuracy"))
  metody = strsplit2(colnames(temp), "_")[,1]
  return(metody)
}

ks.best_signiture_proposals = function(benchmark_csv = "benchmark1578929876.21765.csv", without_train = F){
benchmark = read.csv(benchmark_csv, stringsAsFactors = F)
rownames(benchmark) = make.names(benchmark$method, unique = T)
acc = dplyr::select(benchmark, ends_with("_train_Accuracy"), ends_with("_test_Accuracy"),ends_with("_valid_Accuracy") )
if (without_train == T) { acc = dplyr::select(benchmark, ends_with("_test_Accuracy"),ends_with("_valid_Accuracy")) }
acc$metaindex = rowMeans(acc)
acc$method = rownames(benchmark)
acc$miRy = benchmark$miRy
rownames(acc) = make.names(benchmark$method, unique = T)
acc = acc %>% arrange(desc(metaindex))
return(acc)
}

ks.best_signiture_proposals_meta11 = function(benchmark_csv = "benchmark1578929876.21765.csv"){
  benchmark = read.csv(benchmark_csv, stringsAsFactors = F)
  rownames(benchmark) = make.names(benchmark$method, unique = T)
  temp1 =  dplyr::select(benchmark, ends_with("_valid_Sensitivity"))
  temp2 = dplyr::select(benchmark, ends_with("_valid_Specificity"))
  #acc = dplyr::select(benchmark, ends_with("_train_Accuracy"), ends_with("_test_Accuracy"),ends_with("_valid_Accuracy") )
  #if (without_train == T) { acc = dplyr::select(benchmark, ends_with("_test_Accuracy"),ends_with("_valid_Accuracy")) }
  temp1$temp = rowMeans(temp1)
  temp2$temp = rowMeans(temp2)
  acc = benchmark
  acc$metaindex = temp1$temp + temp2$temp - 1
  acc$method = rownames(benchmark)
  acc$miRy = benchmark$miRy
  rownames(acc) = make.names(benchmark$method, unique = T)
  acc = acc %>% arrange(desc(metaindex))
  return(acc)
}

ks.get_miRNAs_from_benchmark = function(benchmark_csv = "benchmark1578990441.6531.csv", method = "fcsig")
{
  benchmark = read.csv(benchmark_csv, stringsAsFactors = F)
  rownames(benchmark) = make.names(benchmark$method, unique = T)
  return(all.vars(as.formula(benchmark$miRy[which(rownames(benchmark) == method)]))[-1])
}

ks.best_signiture_de = function(selected_miRNAs, use_mix = F)
{
  dane = ks.wczytajmix(replace_smote = F); train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]
  
  wyniki = ks.miRNA.de(trainx, train$Class)
  #write.csv(wyniki, "DE_tylkotrain.csv")
  if (use_mix ==T) { 
   mix = rbind(train,test,valid)
   wyniki = ks.miRNA.de(dplyr::select(mix,starts_with("hsa")), mix$Class) }
  #write.csv(wynikiw, "DE_wszystko.csv")
  
  wyniki = wyniki %>% arrange(`p-value BH`)
  
  return(wyniki[match(selected_miRNAs, wyniki$miR),c("miR","log2FC","p-value BH")])
}

ks.vulcano_plot = function(selected_miRNAs, DE = ks.miRNA.de(), only_label = NULL)
{
  temp = DE[match(selected_miRNAs, DE$miR),]
  tlabels = gsub("\\.","-", selected_miRNAs)
  if (is.null(only_label))
  {
  }
  else {
    tlabels = rep(NA, length(gsub("\\.","-", selected_miRNAs)))
    
    tlabels[match(only_label, DE$miR)] = gsub("\\.","-", selected_miRNAs)[match(only_label, DE$miR)]
  }
  thres = -log10(0.05)
  library(ggplot2)
  library(ggrepel)
  x_limits <- c(-3, +3)
  pval = -log10(temp$`p-value BH`)
  ggplot(data = temp, aes(x = temp$`log2FC`, y = pval, label = tlabels)) + 
    geom_text_repel(arrow = arrow(length = unit(0.01, "npc"), type = "closed", ends = "last"),
                    force = 50) +
    geom_point(color = 'black') +
    theme_classic(base_size = 16) +
    labs(x = "Log2(FC)") +
    labs(y = "-Log10(adjusted p-value)") +
    geom_hline(yintercept=thres, linetype="dashed", color = "red") +
    geom_vline(xintercept = 0, linetype="dotted", 
               color = "blue", size=0.5) + 
    xlim(-3,3) +
    ylim(0,8)
}

ks.mydist=function(c) {dist(c,method="euclidian")}
ks.myclust=function(c) {hclust(c,method="average")}

ks.heatmap = function(x = trainx[,1:10], rlab = data.frame(Batch = dane$Batch, Class = dane$Class), zscore = F) {

  assigcolor = character()
  assigcode = character()
  kolory = rep(palette(),20)[-1]
  kolor_i = 1
  for (i in 1:ncol(rlab)){
    o_ile = as.numeric(length(unique(rlab[,i])))
    assigcode = c(assigcode, as.character(unique(rlab[,i])))
    assigcolor = c(assigcolor, kolory[kolor_i:(kolor_i+o_ile-1)])
    #levels(rlab[,i]) = topo.colors(length(unique(rlab[,i])))
    levels(rlab[,i]) = kolory[kolor_i:(kolor_i+o_ile-1)]
    kolor_i = kolor_i+o_ile
  }
  assig = data.frame(assigcode, assigcolor)
  #levels(rlab$Batch) = rainbow(length(unique(dane$Batch)))
  #levels(rlab$Class) = c("red","green") # red - cancer, green - control
  x2 = as.matrix(x)
  colnames(x2) = gsub("\\.","-", colnames(x2))
  
  if(zscore == F) {
    brks<-ks.diverge_color(x2, centeredOn = median(x2))
    
  # colors = seq(min(x2), max(x2), by = 0.01)
  # my_palette <- colorRampPalette(c("blue", "white", "red"))(n = length(colors) - 1)
  
  rlab = as.matrix(rlab)
  ks.heatmap.3(x2, hclustfun=ks.myclust, distfun=ks.mydist, 
               RowSideColors=t(rlab), 
               margins = c(10,10),
            KeyValueName="log10(TPM)",
            symm=F,symkey=F,symbreaks=T, scale="none",
            col=as.character(brks[[2]]), 
            breaks=as.numeric(brks[[1]]$brks),
            #legend = T
            #,scale="column"
  )
  legend("topright",
         assigcode, fill=assigcolor, horiz=F, bg="transparent", cex=0.5)
  } else {
    x3 = x2
    for(i in 1:ncol(x2)) {
      x3[,i] = scale(x2[,i])
    }
    
    brks<-ks.diverge_color(x3, centeredOn = median(0))
    
    # colors = seq(min(x2), max(x2), by = 0.01)
    # my_palette <- colorRampPalette(c("blue", "white", "red"))(n = length(colors) - 1)
    
    rlab = as.matrix(rlab)
    ks.heatmap.3(x3, hclustfun=ks.myclust, distfun=ks.mydist, 
                 RowSideColors=t(rlab), 
                 margins = c(10,10),
                 KeyValueName="Z-score log10(TPM)",
                 symm=F,symkey=F,symbreaks=T, scale="none",
                 col=as.character(brks[[2]]), 
                 breaks=as.numeric(brks[[1]]$brks),
                 #legend = T
                 #,scale="column"
    )
    legend("topright",
           assigcode, fill=assigcolor, horiz=F, bg="transparent", cex=0.5)
  }
}

ks.diverge_color <- function(data,centeredOn=0){
  library(classInt)
  nHalf=50
  Min <- min(data,na.rm=TRUE)
  Max <- max(data,na.rm=TRUE)
  Thresh <- centeredOn
  pal<-colorRampPalette(c("blue", "white", "red"))(n = 11)
  rc1<-colorRampPalette(colors=c(pal[1],pal[2]),space="Lab")(10)
  for(i in 2:10){
    tmp<-colorRampPalette(colors=c(pal[i],pal[i+1]),space="Lab")(10)
    rc1<-c(rc1,tmp)
  }
  rb1 <- seq(Min, Thresh, length.out=nHalf+1)
  rb2 <- seq(Thresh, Max, length.out=nHalf+1)[-1]
  rampbreaks <- c(rb1, rb2)
  cuts <- classIntervals(data, style="fixed",fixedBreaks=rampbreaks)
  return(list(cuts,rc1))
}

# filter out by default if less than 10 counts in more than 1/2 samples.
ks.counts_to_log10tpm = function(danex, metadane = metadane, ids = metadane$ID, filtr = T, 
                                 filtr_minimalcounts = 10, filtr_howmany = 1/2) {
  danex = as.matrix(danex)
  
  dane_counts = t(danex)
  colnames(dane_counts) = ids
  mode(dane_counts) = 'numeric'
  
  # Normalizacja do TPM
  dane3 = DGEList(counts=dane_counts, genes=data.frame(miR = rownames(dane_counts)), samples = metadane)
  tpm = cpm(dane3, normalized.lib.sizes=F, log = F, prior.count = 0.001)
  tpm = tpm + 0.001
  tpm = log10(tpm)
  ttpm = t(tpm)
  saveRDS(dane3,"TPM_DGEList.rds")
  cat("\nDGEList unfiltered object with TPM was saved as TPM_DGEList.rds.")
  
  if (filtr == F) {
    cat("\nReturned data are log10(TPM).")
  return(ttpm)
  } else {
    zostaw = vector()
    zostaw[1] = TRUE
    for (i in 1:nrow(dane3)) 
    {
      temp = as.numeric(dane3$counts[i,])
      if (sum(temp < filtr_minimalcounts) > (filtr_howmany)*length(temp)) { zostaw[i] = FALSE } else { zostaw[i] = TRUE }
    }
    
    dane4 = dane3[zostaw,]
    ttpm_pofiltrze = ttpm[,which(zostaw)]
    
    cat("\nDGEList filtered object with TPM was saved as TPM_DGEList_filtered.rds.")
    cat(paste0("\n(After filtering) miRNAs left: ", sum(zostaw), " | filtered out: ", sum(zostaw==FALSE), "."))
    saveRDS(dane4,"TPM_DGEList_filtered.rds")
    cat("\nReturned data are log10(TPM).")
    return(ttpm_pofiltrze)
  }
  
}

ks.PCA = function(ttpm_pofiltrze, meta) {
  dane.pca <- prcomp(ttpm_pofiltrze, scale. = TRUE)
  library(ggbiplot)
  ggbiplot(dane.pca,var.axes = FALSE,ellipse=TRUE,circle=TRUE, groups=as.factor(meta))
}

ks.PCA_3D = function(ttpm_pofiltrze, meta) {
  dane.pca <- prcomp(ttpm_pofiltrze, scale. = TRUE)
  library(plotly)
  pc = as.data.frame(dane.pca$x)
  pc = cbind(pc, meta)
  
  p <- plot_ly(data = pc, x = ~PC1, y = ~PC2, z = ~PC3, color = ~meta) %>%
    add_markers() %>%
    layout(scene = list(xaxis = list(title = 'PC1'),
                        yaxis = list(title = 'PC2'),
                        zaxis = list(title = 'PC3')))
  
  p
}

# tempp must contain Class
ks.prepare_split = function(metadane = metadane, ttpm = ttpm_pofiltrze, train_proc = 0.6)
{
  tempp = cbind(metadane, ttpm)
  # Podzial - http://topepo.github.io/caret/data-splitting.html#simple-splitting-based-on-the-outcome
  set.seed(1)
  mix = rep("unasign", nrow(tempp))
  library(caret)
  train.index <- createDataPartition(tempp$Class, p = train_proc, list = FALSE)
  mix[train.index] = "train"
  train = tempp[train.index,]
  

  tempp2 =  tempp[-train.index,]
  train.index <- createDataPartition(tempp2$Class, p = .5, list = FALSE)
  test = tempp2[train.index,]
  valid = tempp2[-train.index,]
  write.csv(train, "mixed_train.csv",row.names = F)
  write.csv(test, "mixed_test.csv",row.names = F)
  write.csv(valid, "mixed_valid.csv",row.names = F)
  
  train$mix = "train"
  test$mix = "test"
  valid$mix = "valid"
  
  mixed = rbind(train,test,valid)
  
  metadane2 = dplyr::select(mixed, -starts_with("hsa"))
  ttpm2 = dplyr::select(mixed, starts_with("hsa"))
  
  mixed2 = cbind(metadane2, ttpm2)
  write.csv(mixed2, "mixed.csv",row.names = F)
  cat("\nSaved 3 sets as csv in working directory. Retruned mixed dataset.")
  return(mixed)
}

ks.keras_create_model <- function(i, hyperparameters, how_many_features = ncol(x_train_scale)) {
  # tempmodel <- keras_model_sequential() %>% 
  #   { if(hyperparameters[i,10]==T) { layer_dense(. , units = hyperparameters[i,1], kernel_regularizer = regularizer_l2(l = 0.001),
  #                                                activation = hyperparameters[i,4], input_shape = c(ncol(x_train_scale))) } else 
  #                                  { layer_dense(. , units = hyperparameters[i,1], activation = hyperparameters[i,4], 
  #                                                                input_shape = c(ncol(x_train_scale))) } } %>%
  #   { if(hyperparameters[i,7]>0) { layer_dropout(. , rate = hyperparameters[i,7]) } else { . } } %>%
  #   { if(hyperparameters[i,2]>0) {
  #   if(hyperparameters[i,11]==T) { layer_dense(. , units = hyperparameters[i,2], activation = hyperparameters[i,5], 
  #               kernel_regularizer = regularizer_l2(l = 0.001)) } else { 
  #                 layer_dense(units = hyperparameters[i,2], activation = hyperparameters[i,5]) } } }  %>% 
  #   { if(hyperparameters[i,8]>0) { layer_dropout(rate = hyperparameters[i,8]) } else { . } } %>%
  #   { if(hyperparameters[i,3]>0) {
  #   if(hyperparameters[i,12]==T) { layer_dense(units = hyperparameters[i,3], activation = hyperparameters[i,6], 
  #                                                kernel_regularizer = regularizer_l2(l = 0.001)) } else 
  #                                   {layer_dense(units = hyperparameters[i,3], activation = hyperparameters[i,6])} } else { . } } %>% 
  #   { if(hyperparameters[i,9]>0) { layer_dropout(rate = hyperparameters[i,9]) } else { . } } %>%
  #   layer_dense(units = 1, activation = 'sigmoid')
  library(keras)
  tempmodel <- keras_model_sequential()
  if(hyperparameters[i,10]==T) { layer_dense(tempmodel , units = hyperparameters[i,1], kernel_regularizer = regularizer_l2(l = 0.001),
                                             activation = hyperparameters[i,4], input_shape = c(how_many_features)) } else 
                                             { layer_dense(tempmodel , units = hyperparameters[i,1], activation = hyperparameters[i,4], 
                                                           input_shape = c(how_many_features)) }
  if(hyperparameters[i,7]>0) { layer_dropout(tempmodel , rate = hyperparameters[i,7]) }
  if(hyperparameters[i,2]>0) {
    if(hyperparameters[i,11]==T) { layer_dense(tempmodel , units = hyperparameters[i,2], activation = hyperparameters[i,5], 
                                               kernel_regularizer = regularizer_l2(l = 0.001)) } else 
                                               {layer_dense(tempmodel, units = hyperparameters[i,2], activation = hyperparameters[i,5]) } }
  
  if(hyperparameters[i,2]>0 & hyperparameters[i,8]>0) { layer_dropout(tempmodel, rate = hyperparameters[i,8]) }
  if(hyperparameters[i,3]>0) {
    if(hyperparameters[i,12]==T) { layer_dense(tempmodel, units = hyperparameters[i,3], activation = hyperparameters[i,6], 
                                               kernel_regularizer = regularizer_l2(l = 0.001)) } else 
                                               { layer_dense(tempmodel, units = hyperparameters[i,3], activation = hyperparameters[i,6])} }
  if(hyperparameters[i,3]>0 & hyperparameters[i,9]>0) { layer_dropout(rate = hyperparameters[i,9]) }
  layer_dense(tempmodel, units = 2, activation = 'softmax')
  
  print(tempmodel)
  
  
  dnn_class_model = keras::compile(tempmodel, optimizer = hyperparameters[i,13],
                            loss = 'binary_crossentropy',
                            metrics = 'accuracy')
  
}

# autoencoder ma softmax na deep feature
# jesli zmienna autoencoder >0 budujemy autoencoder 5-wartwowy, jesli <0 budujemy autoencoder 3-warstwowy
ks.deep_learning = function(selected_miRNAs = ".", wd = getwd(),
                            SMOTE = F, keras_batch_size = 64, clean_temp_files = T,
                            save_threshold_trainacc = 0.85, save_threshold_testacc = 0.8, keras_epochae = 5000,
                            keras_epoch = 2000, keras_patience = 50, 
                            hyperparameters = expand.grid(layer1 = seq(3,11, by = 2), layer2 = c(0,seq(3,11, by = 2)), layer3 = c(0,seq(3,11, by = 2)),
                                                          activation_function_layer1 = c("relu","sigmoid"), activation_function_layer2 = c("relu","sigmoid"), activation_function_layer3 = c("relu","sigmoid"),
                                                          dropout_layer1 = c(0, 0.1), dropout_layer2 = c(0, 0.1), dropout_layer3 = c(0),
                                                          layer1_regularizer = c(T,F), layer2_regularizer = c(T,F), layer3_regularizer = c(T,F),
                                                          optimizer = c("adam","rmsprop","sgd"), autoencoder = c(0,7,-7), balanced = SMOTE, formula = as.character(ks.miR.formula(selected_miRNAs))[3], scaled = c(T,F),
                                                          stringsAsFactors = F),
                            keras_threads = ceiling(parallel::detectCores()/2), start = 1, end = nrow(hyperparameters), output_file = "deeplearning_results.csv", save_all_vars = F) 
  {
  
  codename = sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(output_file))
  #options(warn=-1)
  oldwd = getwd()
  setwd = setwd(wd)
  set.seed(1)
  temp_dir = tempdir()
  options(bitmapType = 'cairo', device = 'png')
  library(dplyr)
  library(keras)
  library(plyr)
  library(foreach)
  library(doParallel)
  #library(doSNOW)
  library(data.table)
  fwrite(hyperparameters, paste0("hyperparameters_",output_file))
  
  dane = ks.wczytajmix(wd = wd, replace_smote = F); train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]
  if (SMOTE == T) { train = train_smoted }
  
  #cores=detectCores()
  cat(paste0("\nTemp dir: ", temp_dir, "\n"))
  cat("\nStarting preparing cluster..\n")
  #cl <- makePSOCKcluster(keras_threads) #not to overload your computer
  cl = makeCluster(keras_threads, outfile=paste0("temp/", ceiling(as.numeric(Sys.time())), "deeplearning_cluster.log"))
  registerDoParallel(cl)
  on.exit(stopCluster(cl))
  cat("\nCluster prepared..\n")
  

  
  #args= names(mget(ls()))
  #export = export[!export %in% args]
  
  # tu musi isc iteracja
  cat(paste0("\nStarting parallel loop.. There are: ", end-start+1, " hyperparameter sets to be checked.\n"))
  final <- foreach(i=as.numeric(start):as.numeric(end), .combine=rbind, .verbose=T, .inorder=F
                   ,.errorhandling="remove", .export = ls(), .packages = loadedNamespaces()
                   ) %dopar% {
    start_time <- Sys.time()
    library(keras)
    library(ggplot2)
    library(dplyr)
    library(data.table)
    set.seed(1)
    
    cat("\nStarting hyperparameters..\n")
    print(hyperparameters[i,])
    source("/home/konrad/snorlax/helpers/ks_miRNAhelpers.R")
    options(bitmapType = 'cairo', device = 'png')
    model_id = paste0(format(i, scientific = FALSE), "-", ceiling(as.numeric(Sys.time())))
    if(SMOTE == T) { model_id = paste0(format(i, scientific = FALSE), "-SMOTE-", ceiling(as.numeric(Sys.time()))) }
    tempwyniki = data.frame(model_id=model_id)
    tempwyniki[1, "model_id"] = model_id

    
    if(!dir.exists(paste0(temp_dir,"/models"))) { dir.create(paste0(temp_dir,"/models"))}
    if(!dir.exists(paste0(temp_dir,"/models/keras",model_id))) { dir.create(paste0(temp_dir,"/models/keras",model_id))}
    cat(paste0("\nTraining model: ",temp_dir,"/models/keras",model_id,"\n"))
    #pdf(paste0(temp_dir,"/models/keras",model_id,"/plots.pdf"), paper="a4")
    
    con <- file(paste0(temp_dir,"/models/keras",model_id,"/training.log"))
    sink(con, append=TRUE)
    sink(con, append=TRUE, type="message")
    
    early_stop <- callback_early_stopping(monitor = "val_loss", mode="min", patience = keras_patience)
    cp_callback <- callback_model_checkpoint(
      filepath =  paste0(temp_dir,"/models/keras",model_id,"/finalmodel.hdf5"),
      save_best_only = TRUE, period = 5, monitor = "val_loss",
      verbose = 0
    )
    ae_cp_callback <- callback_model_checkpoint(
      filepath =  paste0(temp_dir,"/models/keras",model_id,"/autoencoderweights.hdf5"),
      save_best_only = TRUE, save_weights_only = T, period = 5, monitor = "val_loss",
      verbose = 0
    )
    
    x_train <- train %>% 
      { if (selected_miRNAs != ".") { dplyr::select(., selected_miRNAs) } else { dplyr::select(., starts_with("hsa")) } } %>% 
      as.matrix()
    y_train <- train %>% 
      dplyr::select("Class") %>% 
      as.matrix()
    y_train[,1] = ifelse(y_train[,1] == "Cancer",1,0)
    
    
    x_test <- test %>% 
      { if (selected_miRNAs != ".") { dplyr::select(.,selected_miRNAs) } else { dplyr::select(.,starts_with("hsa")) } } %>% 
      as.matrix()
    y_test <- test %>% 
      dplyr::select("Class") %>% 
      as.matrix()
    y_test[,1] = ifelse(y_test[,1] == "Cancer",1,0)
    
    x_valid <- valid %>% 
      { if (selected_miRNAs != ".") { dplyr::select(.,selected_miRNAs) } else { dplyr::select(.,starts_with("hsa")) } } %>% 
      as.matrix()
    y_valid <- valid %>% 
      dplyr::select("Class") %>% 
      as.matrix()
    y_valid[,1] = ifelse(y_valid[,1] == "Cancer",1,0)
    
    if(hyperparameters[i, 17] == T) {
    x_train_scale = x_train %>% scale()
    
    col_mean_train <- attr(x_train_scale, "scaled:center") 
    col_sd_train <- attr(x_train_scale, "scaled:scale")
    
    x_test_scale <- x_test %>% 
      scale(center = col_mean_train, 
            scale = col_sd_train)
    
    x_valid_scale <- x_valid %>% 
      scale(center = col_mean_train, 
            scale = col_sd_train)
    }
    else {
      x_train_scale = x_train
      x_test_scale <- x_test
      x_valid_scale <- x_valid
    }
    
    input_layer <- layer_input(shape = c(ncol(x_train_scale)))
    
    #psych::describe(x_test_scale)

    # czy autoenkoder?
    if(hyperparameters[i,14] != 0) {
      
      
      n1 <- hyperparameters[i,14]
      n3 <- ncol(x_train_scale)
      
      if (hyperparameters[i,14]>0) {
        encoder <- 
          input_layer %>%
          #layer_dense(units = ceiling(n3/2), activation = hyperparameters[i,6])    %>% 
          layer_dense(units = n1, activation = "sigmoid")  # dimensions of final encoding layer
        
        decoder <- encoder %>% 
          #layer_dense(units = ceiling(n3/2), activation = hyperparameters[i,6]) %>% 
          layer_dense(units = n3, hyperparameters[i,6])  # dimension of original variable
      }
      else {
        n1 = -n1
        encoder <- 
          input_layer %>%
          layer_dense(units = n1, activation = "sigmoid")  # dimensions of final encoding layer
        
        decoder <- encoder %>% 
          layer_dense(units = n3, hyperparameters[i,6])  # dimension of original variable
      }
      ae_model <- keras_model(inputs = input_layer, outputs = decoder)
      ae_model
      
      ae_model %>% 
        keras::compile(loss = "mean_absolute_error",
                optimizer = hyperparameters[i,13],
                metrics = c("mean_squared_error"))
      
      summary(ae_model)
      

      
      #ae_tempmodelfile = tempfile()
      ae_history <- fit(ae_model, x = x_train_scale,
                   y = x_train_scale,
                   epochs = keras_epochae, batch_size = keras_batch_size,
                   shuffle = T, verbose = 0,
                   view_metrics = FALSE,
                   validation_data = list(x_test_scale, x_test_scale),
                   callbacks = list(ae_cp_callback, early_stop))
      #file.copy(ae_tempmodelfile, paste0(temp_dir,"/models/keras",model_id,"/finalmodel.hdf5"), overwrite = T)
      
      ae_history_df <- as.data.frame(ae_history)
      fwrite(ae_history_df, paste0(temp_dir,"/models/keras",model_id,"/autoencoder_history_df.csv.gz"))
      
      compare_cx <- data.frame(
        train_loss = ae_history$metrics$loss,
        test_loss = ae_history$metrics$val_loss
      ) %>%
        tibble::rownames_to_column() %>%
        mutate(rowname = as.integer(rowname)) %>%
        tidyr::gather(key = "type", value = "value", -rowname)
      
      plot1 = ggplot(compare_cx, aes(x = rowname, y = value, color = type)) +
        geom_line() +
        xlab("epoch") +
        ylab("loss") +
        theme_bw()
      
      
      ggplot2::ggsave(file = paste0(temp_dir,"/models/keras",model_id,"/training_autoencoder.png"), grid.arrange(plot1, nrow =1, top = "Training of autoencoder"))
      
      cat("\n- Reloading autoencoder to get weights...\n")
      encoder_model <- keras_model(inputs = input_layer, outputs = encoder)
      encoder_model %>% load_model_weights_hdf5(paste0(temp_dir,"/models/keras",model_id,"/autoencoderweights.hdf5"), 
                                                skip_mismatch = T, 
                                                by_name = T)
      
      cat("\n- Autoencoder saving...\n")
      save_model_hdf5(encoder_model, paste0(temp_dir,"/models/keras",model_id,"/autoencoder.hdf5"))

      
      cat("\n- Creating deep features...\n")
      
      ae_x_train_scale <- encoder_model %>% 
        predict(x_train_scale) %>% 
        as.matrix()
      fwrite(ae_x_train_scale, paste0(temp_dir,"/models/keras",model_id,"/deepfeatures_train.csv"))
      
      ae_x_test_scale <- encoder_model %>% 
        predict(x_test_scale) %>% 
        as.matrix()
      fwrite(ae_x_test_scale, paste0(temp_dir,"/models/keras",model_id,"/deepfeatures_test.csv"))
      
      ae_x_valid_scale <- encoder_model %>% 
        predict(x_valid_scale) %>% 
        as.matrix()
      fwrite(ae_x_valid_scale, paste0(temp_dir,"/models/keras",model_id,"/deepfeatures_valid.csv"))
      
      # # podmiana żeby nie edytować kodu
      x_train_scale = as.matrix(ae_x_train_scale)
      x_test_scale = as.matrix(ae_x_test_scale)
      x_valid_scale = as.matrix(ae_x_valid_scale)
      
      cat("\n- Training model based on deep features...\n")
      dnn_class_model <- ks.keras_create_model(i, hyperparameters = hyperparameters, how_many_features = ncol(x_train_scale))
      
      #tempmodelfile = tempfile()
      history <- fit(dnn_class_model, x = x_train_scale, 
                   y = to_categorical(y_train),
                   epochs = keras_epoch, 
                   validation_data = list(x_test_scale, to_categorical(y_test)),
                   callbacks = list(cp_callback, early_stop),
                   verbose = 0,
                   view_metrics = FALSE,
                   batch_size = keras_batch_size, shuffle = T)
      print(history)
      #plot(history, col="black")
      history_df <- as.data.frame(history)
      fwrite(history_df, paste0(temp_dir,"/models/keras",model_id,"/history_df.csv.gz"))
      
      cat("\n- Saving history and plots...\n")
      compare_cx <- data.frame(
        train_loss = history$metrics$loss,
        test_loss = history$metrics$val_loss
      ) %>%
        tibble::rownames_to_column() %>%
        mutate(rowname = as.integer(rowname)) %>%
        tidyr::gather(key = "type", value = "value", -rowname)
      
      plot1 = ggplot(compare_cx, aes(x = rowname, y = value, color = type)) +
        geom_line() +
        xlab("Epoch") +
        ylab("Loss") +
        theme_bw()
      
      compare_cx <- data.frame(
        train_accuracy = history$metrics$accuracy,
        test_accuracy = history$metrics$val_accuracy
      ) %>%
        tibble::rownames_to_column() %>%
        mutate(rowname = as.integer(rowname)) %>%
        tidyr::gather(key = "type", value = "value", -rowname)
      
      plot2 = ggplot(compare_cx, aes(x = rowname, y = value, color = type)) +
        geom_line() +
        xlab("Epoch") +
        ylab("Accuracy") + 
        theme_bw()
      
      ggplot2::ggsave(file = paste0(temp_dir,"/models/keras",model_id,"/training.png"), grid.arrange(plot1, plot2, nrow =2, top = "Training of final neural network"))
      
      # ewaluacja
      cat("\n- Saving final model...\n")
      dnn_class_model = load_model_hdf5(paste0(temp_dir,"/models/keras",model_id,"/finalmodel.hdf5"))
      y_train_pred <- predict(object = dnn_class_model, x = x_train_scale)
      y_test_pred <- predict(object = dnn_class_model, x = x_test_scale)
      y_valid_pred <- predict(object = dnn_class_model, x = x_valid_scale)
      
      
      # wybranie odciecia
      pred = data.frame(`Class` = train$Class, `Pred` = y_train_pred)
      library(cutpointr)
      cutoff = cutpointr(pred, Pred.2, Class, pos_class = "Cancer", metric = youden)
      summary(cutoff)
      ggplot2::ggsave(file = paste0(temp_dir,"/models/keras",model_id,"/cutoff.png"), plot(cutoff))
      wybrany_cutoff = cutoff$optimal_cutpoint
      #wyniki[i, "training_AUC"] = cutoff$AUC
      tempwyniki[1, "training_AUC"] = cutoff$AUC
      #wyniki[i, "cutoff"] = wybrany_cutoff
      tempwyniki[1, "cutoff"] = wybrany_cutoff
      cat(paste0("\n\n---- TRAINING AUC: ",cutoff$AUC," ----\n\n"))
      cat(paste0("\n\n---- OPTIMAL CUTOFF: ",wybrany_cutoff," ----\n\n"))
      cat(paste0("\n\n---- TRAINING PERFORMANCE ----\n\n"))
      pred$PredClass = ifelse(pred$Pred.2 >= wybrany_cutoff, "Cancer", "Control")
      pred$PredClass = factor(pred$PredClass, levels = c("Control","Cancer"))
      cm_train = caret::confusionMatrix(pred$PredClass, pred$Class, positive = "Cancer")
      print(cm_train)
      
      t1_roc = pROC::roc(Class ~ as.numeric(Pred.2), data=pred)
      tempwyniki[1, "training_AUC2"] = t1_roc$auc
      tempwyniki[1, "training_AUC_lower95CI"] = as.character(ci(t1_roc))[1]
      tempwyniki[1, "training_AUC_upper95CI"] = as.character(ci(t1_roc))[3]
      saveRDS(t1_roc, paste0(temp_dir,"/models/keras",model_id,"/training_ROC.RDS"))
      
      tempwyniki[1, "training_Accuracy"] = cm_train$overall[1]
      tempwyniki[1, "training_Sensitivity"] = cm_train$byClass[1]
      tempwyniki[1, "training_Specificity"] = cm_train$byClass[2]
      tempwyniki[1, "training_PPV"] = cm_train$byClass[3]
      tempwyniki[1, "training_NPV"] = cm_train$byClass[4]
      tempwyniki[1, "training_F1"] = cm_train$byClass[7]
      saveRDS(cm_train, paste0(temp_dir,"/models/keras",model_id,"/cm_train.RDS"))
      
      cat(paste0("\n\n---- TESTING PERFORMANCE ----\n\n"))
      pred = data.frame(`Class` = test$Class, `Pred` = y_test_pred)
      pred$PredClass = ifelse(pred$Pred.2 >= wybrany_cutoff, "Cancer", "Control")
      pred$PredClass = factor(pred$PredClass, levels = c("Control","Cancer"))
      cm_test = caret::confusionMatrix(pred$PredClass, pred$Class, positive = "Cancer")
      print(cm_test)
      tempwyniki[1, "test_Accuracy"] = cm_test$overall[1]
      tempwyniki[1, "test_Sensitivity"] = cm_test$byClass[1]
      tempwyniki[1, "test_Specificity"] = cm_test$byClass[2]
      tempwyniki[1, "test_PPV"] = cm_test$byClass[3]
      tempwyniki[1, "test_NPV"] = cm_test$byClass[4]
      tempwyniki[1, "test_F1"] = cm_test$byClass[7]
      saveRDS(cm_test, paste0(temp_dir,"/models/keras",model_id,"/cm_test.RDS"))
      
      cat(paste0("\n\n---- VALIDATION PERFORMANCE ----\n\n"))
      pred = data.frame(`Class` = valid$Class, `Pred` = y_valid_pred)
      pred$PredClass = ifelse(pred$Pred.2 >= wybrany_cutoff, "Cancer", "Control")
      pred$PredClass = factor(pred$PredClass, levels = c("Control","Cancer"))
      cm_valid = caret::confusionMatrix(pred$PredClass, pred$Class, positive = "Cancer")
      print(cm_valid)
      tempwyniki[1, "valid_Accuracy"] = cm_test$overall[1]
      tempwyniki[1, "valid_Sensitivity"] = cm_test$byClass[1]
      tempwyniki[1, "valid_Specificity"] = cm_test$byClass[2]
      tempwyniki[1, "valid_PPV"] = cm_test$byClass[3]
      tempwyniki[1, "valid_NPV"] = cm_test$byClass[4]
      tempwyniki[1, "valid_F1"] = cm_test$byClass[7]
      saveRDS(cm_valid, paste0(temp_dir,"/models/keras",model_id,"/cm_valid.RDS"))
      
      mix = rbind(train,test,valid)
      mixx = rbind(x_train_scale, x_test_scale, x_valid_scale)
      y_mixx_pred <- predict(object = dnn_class_model, x = mixx)
      
      mix$Podzial = c(rep("Training",nrow(train)),rep("Test",nrow(test)),rep("Validation",nrow(valid)))
      mix$Pred = y_mixx_pred[,2]
      mix$PredClass = ifelse(mix$Pred >= wybrany_cutoff, "Cancer", "Control")
      fwrite(mix, paste0(temp_dir,"/models/keras",model_id,"/data_predictions.csv.gz"))
      fwrite(cbind(hyperparameters[i,], tempwyniki), paste0(temp_dir,"/models/keras",model_id,"/wyniki.csv"))
      
      wagi = get_weights(dnn_class_model)
      saveRDS(wagi, paste0(temp_dir,"/models/keras",model_id,"/finalmodel_weights.RDS"))
      save_model_weights_hdf5(dnn_class_model, paste0(temp_dir,"/models/keras",model_id,"/finalmodel_weights.hdf5"))
      saveRDS(dnn_class_model, paste0(temp_dir,"/models/keras",model_id,"/finalmodel.RDS"))
      
      # czy jest sens zapisywac?
      #sink() 
      #sink(type="message")
      cat(paste0("\n\n== ",model_id, ": ", tempwyniki[1, "training_Accuracy"], " / ", tempwyniki[1, "test_Accuracy"], " ==> ", tempwyniki[1, "training_Accuracy"]>save_threshold_trainacc & tempwyniki[1, "test_Accuracy"]>save_threshold_testacc))
      if(tempwyniki[1, "training_Accuracy"]>save_threshold_trainacc & tempwyniki[1, "test_Accuracy"]>save_threshold_testacc) {
      # zapisywanie modelu do właściwego katalogu
      if (save_all_vars) { save(list = ls(all=TRUE), file = paste0(temp_dir,"/models/keras",model_id,"/all.Rdata.gz"), compress = "gzip", compression_level = 9) }
      if(!dir.exists(paste0("models/",codename,"/"))) { dir.create(paste0("models/",codename,"/")) }
      zip(paste0("models/",codename,"/",codename, "_", model_id,".zip"),list.files(paste0(temp_dir,"/models/keras",model_id), full.names = T, recursive = T, include.dirs = T))
      # file.copy(list.files(paste0(temp_dir,"/models/keras",model_id), pattern = "_wyniki.csv$", full.names = T, recursive = T, include.dirs = T),paste0("temp/",codename,"_",model_id,"_deeplearningresults.csv"))
      if (!dir.exists(paste0("temp/",codename,"/"))) { dir.create(paste0("temp/",codename,"/")) }
      file.copy(list.files(paste0(temp_dir,"/models/keras",model_id), pattern = "_wyniki.csv$", full.names = T, recursive = T, include.dirs = T),paste0("temp/",codename,"/",model_id,"_deeplearningresults.csv"))
      #dev.off()
      }
    } else {
      x_train <- train %>% 
        { if (selected_miRNAs != ".") { dplyr::select(.,selected_miRNAs) } else { dplyr::select(.,starts_with("hsa")) } } %>% 
        as.matrix()
      y_train <- train %>% 
        dplyr::select("Class") %>% 
        as.matrix()
      y_train[,1] = ifelse(y_train[,1] == "Cancer",1,0)
      
      
      x_test <- test %>% 
        { if (selected_miRNAs != ".") { dplyr::select(.,selected_miRNAs) } else { dplyr::select(.,starts_with("hsa")) } } %>% 
        as.matrix()
      y_test <- test %>% 
        dplyr::select("Class") %>% 
        as.matrix()
      y_test[,1] = ifelse(y_test[,1] == "Cancer",1,0)
      
      x_valid <- valid %>% 
        { if (selected_miRNAs != ".") { dplyr::select(.,selected_miRNAs) } else { dplyr::select(.,starts_with("hsa")) } } %>% 
        as.matrix()
      y_valid <- valid %>% 
        dplyr::select("Class") %>% 
        as.matrix()
      y_valid[,1] = ifelse(y_valid[,1] == "Cancer",1,0)
      
      if(hyperparameters[i, 17] == T) {
        x_train_scale = x_train %>% scale()
        
        col_mean_train <- attr(x_train_scale, "scaled:center") 
        col_sd_train <- attr(x_train_scale, "scaled:scale")
        
        x_test_scale <- x_test %>% 
          scale(center = col_mean_train, 
                scale = col_sd_train)
        
        x_valid_scale <- x_valid %>% 
          scale(center = col_mean_train, 
                scale = col_sd_train)
      }
      else {
        x_train_scale = x_train
        x_test_scale <- x_test
        x_valid_scale <- x_valid
      }
      
      
      dnn_class_model <- ks.keras_create_model(i, hyperparameters = hyperparameters, how_many_features = ncol(x_train_scale))
      
      history <-  fit(dnn_class_model, x = x_train_scale, 
                   y = to_categorical(y_train),
                   epochs = keras_epoch, 
                   validation_data = list(x_test_scale, to_categorical(y_test)),
                   callbacks = list(
                     cp_callback, 
                   #callback_reduce_lr_on_plateau(monitor = "val_loss", factor = 0.1),
                   #callback_model_checkpoint(paste0(temp_dir,"/models/keras",model_id,"/finalmodel.hdf5")),
                     early_stop),
                   verbose = 0,
                   view_metrics = FALSE,
                   batch_size = keras_batch_size, shuffle = T)
      print(history)
      #plot(history, col="black")
      history_df <- as.data.frame(history)
      fwrite(history_df, paste0(temp_dir,"/models/keras",model_id,"/history_df.csv.gz"))
      
      # pdf(paste0(temp_dir,"/models/keras",model_id,"/plots.pdf"))
      compare_cx <- data.frame(
        train_loss = history$metrics$loss,
        test_loss = history$metrics$val_loss
      ) %>%
        tibble::rownames_to_column() %>%
        mutate(rowname = as.integer(rowname)) %>%
        tidyr::gather(key = "type", value = "value", -rowname)
      
      plot1 = ggplot(compare_cx, aes(x = rowname, y = value, color = type)) +
        geom_line() +
        xlab("epoch") +
        ylab("loss") +
        theme_bw()
      
      compare_cx <- data.frame(
        train_accuracy = history$metrics$accuracy,
        test_accuracy = history$metrics$val_accuracy
      ) %>%
        tibble::rownames_to_column() %>%
        mutate(rowname = as.integer(rowname)) %>%
        tidyr::gather(key = "type", value = "value", -rowname)
      
      plot2 = ggplot(compare_cx, aes(x = rowname, y = value, color = type)) +
        geom_line() +
        xlab("epoch") +
        ylab("loss") + 
        theme_bw()
      
      ggplot2::ggsave(file = paste0(temp_dir,"/models/keras",model_id,"/training.png"), grid.arrange(plot1, plot2, nrow =2, top = "Training of final neural network"))
      
      # ewaluacja
      dnn_class_model = load_model_hdf5(paste0(temp_dir,"/models/keras",model_id,"/finalmodel.hdf5"))
      y_train_pred <- predict(object = dnn_class_model, x = x_train_scale)
      y_test_pred <- predict(object = dnn_class_model, x = x_test_scale)
      y_valid_pred <- predict(object = dnn_class_model, x = x_valid_scale)
      
      
      # wybranie odciecia
      pred = data.frame(`Class` = train$Class, `Pred` = y_train_pred)
      library(cutpointr)
      cutoff = cutpointr(pred, Pred.2, Class, pos_class = "Cancer", metric = youden)
      print(summary(cutoff))
      ggplot2::ggsave(file = paste0(temp_dir,"/models/keras",model_id,"/cutoff.png"), plot(cutoff))
      wybrany_cutoff = cutoff$optimal_cutpoint
      #wyniki[i, "training_AUC"] = cutoff$AUC
      tempwyniki[1, "training_AUC"] = cutoff$AUC
      #wyniki[i, "cutoff"] = wybrany_cutoff
      tempwyniki[1, "cutoff"] = wybrany_cutoff
      
      cat(paste0("\n\n---- TRAINING PERFORMANCE ----\n\n"))
      pred$PredClass = ifelse(pred$Pred.2 >= wybrany_cutoff, "Cancer", "Control")
      pred$PredClass = factor(pred$PredClass, levels = c("Control","Cancer"))
      cm_train = caret::confusionMatrix(pred$PredClass, pred$Class, positive = "Cancer")
      print(cm_train)
      
      t1_roc = pROC::roc(Class ~ as.numeric(Pred.2), data=pred)
      tempwyniki[1, "training_AUC2"] = t1_roc$auc
      tempwyniki[1, "training_AUC_lower95CI"] = as.character(ci(t1_roc))[1]
      tempwyniki[1, "training_AUC_upper95CI"] = as.character(ci(t1_roc))[3]
      saveRDS(t1_roc, paste0(temp_dir,"/models/keras",model_id,"/training_ROC.RDS"))
      #ggplot2::ggsave(file = paste0(temp_dir,"/models/keras",model_id,"/training_ROC.png"), grid.arrange(plot(t1_roc), nrow =1, top = "Training ROC curve"))
      
      tempwyniki[1, "training_Accuracy"] = cm_train$overall[1]
      tempwyniki[1, "training_Sensitivity"] = cm_train$byClass[1]
      tempwyniki[1, "training_Specificity"] = cm_train$byClass[2]
      tempwyniki[1, "training_PPV"] = cm_train$byClass[3]
      tempwyniki[1, "training_NPV"] = cm_train$byClass[4]
      tempwyniki[1, "training_F1"] = cm_train$byClass[7]
      saveRDS(cm_train, paste0(temp_dir,"/models/keras",model_id,"/cm_train.RDS"))
      
      cat(paste0("\n\n---- TESTING PERFORMANCE ----\n\n"))
      pred = data.frame(`Class` = test$Class, `Pred` = y_test_pred)
      pred$PredClass = ifelse(pred$Pred.2 >= wybrany_cutoff, "Cancer", "Control")
      pred$PredClass = factor(pred$PredClass, levels = c("Control","Cancer"))
      cm_test = caret::confusionMatrix(pred$PredClass, pred$Class, positive = "Cancer")
      print(cm_test)
      tempwyniki[1, "test_Accuracy"] = cm_test$overall[1]
      tempwyniki[1, "test_Sensitivity"] = cm_test$byClass[1]
      tempwyniki[1, "test_Specificity"] = cm_test$byClass[2]
      tempwyniki[1, "test_PPV"] = cm_test$byClass[3]
      tempwyniki[1, "test_NPV"] = cm_test$byClass[4]
      tempwyniki[1, "test_F1"] = cm_test$byClass[7]
      saveRDS(cm_test, paste0(temp_dir,"/models/keras",model_id,"/cm_test.RDS"))
      
      cat(paste0("\n\n---- VALIDATION PERFORMANCE ----\n\n"))
      pred = data.frame(`Class` = valid$Class, `Pred` = y_valid_pred)
      pred$PredClass = ifelse(pred$Pred.2 >= wybrany_cutoff, "Cancer", "Control")
      pred$PredClass = factor(pred$PredClass, levels = c("Control","Cancer"))
      cm_valid = caret::confusionMatrix(pred$PredClass, pred$Class, positive = "Cancer")
      print(cm_valid)
      tempwyniki[1, "valid_Accuracy"] = cm_test$overall[1]
      tempwyniki[1, "valid_Sensitivity"] = cm_test$byClass[1]
      tempwyniki[1, "valid_Specificity"] = cm_test$byClass[2]
      tempwyniki[1, "valid_PPV"] = cm_test$byClass[3]
      tempwyniki[1, "valid_NPV"] = cm_test$byClass[4]
      tempwyniki[1, "valid_F1"] = cm_test$byClass[7]
      saveRDS(cm_valid, paste0(temp_dir,"/models/keras",model_id,"/cm_valid.RDS"))
      
      mix = rbind(train,test,valid)
      mixx = rbind(x_train_scale, x_test_scale, x_valid_scale)
      y_mixx_pred <- predict(object = dnn_class_model, x = mixx)
      
      mix$Podzial = c(rep("Training",nrow(train)),rep("Test",nrow(test)),rep("Validation",nrow(valid)))
      mix$Pred = y_mixx_pred[,2]
      mix$PredClass = ifelse(mix$Pred >= wybrany_cutoff, "Cancer", "Control")
      fwrite(mix, paste0(temp_dir,"/models/keras",model_id,"/data_predictions.csv.gz"))
      fwrite(cbind(hyperparameters[i,], tempwyniki), paste0(temp_dir,"/models/keras",model_id,"/",model_id,"_wyniki.csv"))
      
      wagi = get_weights(dnn_class_model)
      saveRDS(wagi, paste0(temp_dir,"/models/keras",model_id,"/finalmodel_weights.RDS"))
      save_model_weights_hdf5(dnn_class_model, paste0(temp_dir,"/models/keras",model_id,"/finalmodel_weights.hdf5"))
      saveRDS(dnn_class_model, paste0(temp_dir,"/models/keras",model_id,"/finalmodel.RDS"))
      
      
      # czy jest sens zapisywac?
      
      cat(paste0("\n\n== ",model_id, ": ", tempwyniki[1, "training_Accuracy"], " / ", tempwyniki[1, "test_Accuracy"], " ==> ", tempwyniki[1, "training_Accuracy"]>save_threshold_trainacc & tempwyniki[1, "test_Accuracy"]>save_threshold_testacc))
      if(tempwyniki[1, "training_Accuracy"]>save_threshold_trainacc & tempwyniki[1, "test_Accuracy"]>save_threshold_testacc) {
      # zapisywanie modelu do właściwego katalogu
      if (save_all_vars) { save(list = ls(all=TRUE), file = paste0(temp_dir,"/models/keras",model_id,"/all.Rdata.gz"), compress = "gzip", compression_level = 9) }
      if(!dir.exists(paste0("models/",codename,"/"))) { dir.create(paste0("models/",codename,"/")) }
      zip(paste0("models/",codename,"/",codename, "_", model_id,".zip"),list.files(paste0(temp_dir,"/models/keras",model_id), full.names = T, recursive = T, include.dirs = T))
      if (!dir.exists(paste0("temp/",codename,"/"))) { dir.create(paste0("temp/",codename,"/")) }
      file.copy(list.files(paste0(temp_dir,"/models/keras",model_id), pattern = "_wyniki.csv$", full.names = T, recursive = T, include.dirs = T),paste0("temp/",codename,"/",model_id,"_deeplearningresults.csv"))
    } }
    
    sink() 
    sink(type="message")
    #dev.off()
    tempwyniki2 = cbind(hyperparameters[i,],tempwyniki)
    tempwyniki2[1,"name"] = paste0(codename,"_", model_id)
    tempwyniki2[1,"worth_saving"] = as.character(tempwyniki[1, "training_Accuracy"]>save_threshold_trainacc & tempwyniki[1, "test_Accuracy"]>save_threshold_testacc)
    end_time <- Sys.time()
    tempwyniki2[1,"training_time"] = as.character(end_time - start_time)
    tempwyniki2
  }
  
  saveRDS(final, paste0(output_file,".RDS"))
  cat("\nAll done!! Ending..\n")
  if (file.exists(output_file)) { 
    tempfi = fread(output_file)
    final = rbind(tempfi, final) }
  fwrite(final, output_file)
  setwd(oldwd)
  #options(warn=0)
  # sprzątanie
  if(clean_temp_files) {
  unlink(paste0(normalizePath(temp_dir), "/", dir(temp_dir)), recursive = TRUE) }
}
