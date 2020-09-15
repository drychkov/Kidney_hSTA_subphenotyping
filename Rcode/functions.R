source("package_install_load.R")

#################
# Functions adopted from Hughey JJ, Butte AJ. 2015. Robust meta-analysis of gene expression using the elastic net. Nucleic Acids Research 43:e79–e79. doi:10.1093/nar/gkv229

fixCustomCdfGeneIds = function(geneIds) {
  return(sub('_at', '', geneIds))
}


fixGeoSampleNames = function(sampleNames) {
  sampleNames = paste0(toupper(sampleNames), '_')
  regexResult = regexpr('^GSM[0-9]+[^0-9]', sampleNames)
  sampleNamesNew = mapply(function(sampleName, matchLength) substr(sampleName, 1, matchLength-1),
                          sampleNames, attr(regexResult, 'match.length'))
  return(sampleNamesNew)}


fixCelSampleNames = function(sampleNames) {
  sampleNamesNew = gsub('\\.cel$', '', sampleNames, ignore.case=TRUE)
  return(sampleNamesNew)}

getGeneProbeMappingDirect = function(featureDf, geneColname, probeColname='ID') {
  mapping = featureDf[,c(probeColname, geneColname)]
  mapping = mapping[apply(mapping, MARGIN=1, function(x) all(!is.na(x) & x!='')),]
  mapping = data.frame(lapply(mapping, as.character), stringsAsFactors=FALSE)
  colnames(mapping) = c('probeSet', 'geneId')
  return(mapping)}


getGeneProbeMappingAnno = function(featureDf, dbName, interName) {
  mappingProbeIntermediate = featureDf[!is.na(featureDf[,interName]) & featureDf[,interName]!='', c('ID', interName)]
  colnames(mappingProbeIntermediate) = c('probeSet', 'geneInter')
  mappingProbeIntermediate$geneInter = gsub("\\.[0-9]",'', mappingProbeIntermediate$geneInter)
  mapTmp1 = eval(parse(text=dbName))
  mapTmp2 = mappedkeys(mapTmp1)
  mapTmp3 = as.list(mapTmp1[mapTmp2])
  geneId = do.call(c, mapTmp3)
  geneInter = do.call(c, mapply(function(inter, len) rep_len(inter, len), names(mapTmp3), sapply(mapTmp3, length),
                                SIMPLIFY=FALSE))
  if (dbName == 'org.Hs.egUNIGENE2EG') {
    geneInter = sub('.', '', geneInter, fixed = TRUE)}
  mappingIdInter = data.frame(geneId, geneInter, stringsAsFactors=FALSE)
  mapping = merge(mappingIdInter, mappingProbeIntermediate, by = 'geneInter', sort = FALSE)}



calcExprsByGeneEmat = function(emat, mapping) {
  cat("Mapping probes to genes...", "\n")
  geneIds = unique(mapping[,'geneId'])
  exprsByGene = matrix(nrow = length(geneIds), 
                       ncol = ncol(emat), 
                       dimnames = list(geneIds, colnames(emat)))
  for (geneId in geneIds) {
    probeSet = mapping[mapping[,'geneId'] == geneId, 'probeSet']
    probeSet = probeSet[probeSet %in% rownames(emat)]
    if (length(probeSet) == 0) {
      next
    }
    exprsTmp = emat[probeSet,, drop = FALSE]
    if (nrow(exprsTmp) == 1) {
      exprsByGene[geneId,] = exprsTmp
    } else {
      exprsByGene[geneId,] = rowMedians(t(exprsTmp), na.rm = TRUE)
    }
  }
  # exprsByGene = exprsByGene[complete.cases(exprsByGene),]
  return(exprsByGene)}

cleanStudyData = function(esetList, sampleMeta) {
  # select relevant samples, convert to matrix
  ematList = foreach(studyName = names(esetList)) %do% {
    keepIdx = (colnames(esetList[[studyName]]) %in% 
                 sampleMeta[sampleMeta[,'study'] == studyName, 'sample'])
    exprs(esetList[[studyName]])[,keepIdx]
  }
  names(ematList) = names(esetList)
  return(ematList)
}

# Convert Enterz Ids to Gene Symbols
geneId2Symbol = function(emat) {
  # Remove control probes 
  controlProbes <- grep("AFFX", rownames(emat))
  if (length(controlProbes)!=0) {emat <- emat[-controlProbes,]}
  
  # ematTrain = ematTrain[,sampleTrain$sample]
  geneIds = as.character(rownames(emat))
  geneSymbols = getSYMBOL(geneIds,'org.Hs.eg')
  
  rownames(emat) = unname(geneSymbols)
  emat = emat[!is.na(rownames(emat)),]
  
  return(emat)
}

symbol2geneID = function(genes) {
  geneIds = bitr(genes, fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = "org.Hs.eg.db")
  genes1 = geneIds$ENTREZID
  names(genes1) = geneIds$SYMBOL
  
  notFound = genes[!genes %in% names(genes1)]
  if (length(notFound) > 0) {
    message("Can't find ", paste0(notFound, collapse = ", "), " gene(s).", "\n")
  } 
  return(genes1)
}

instaScore = function(dataMat, features = NULL, model, norm.data = FALSE) {
  if (ncol(dataMat) == 0) return(numeric(0))
  if (is.null(features)) features = rownames(dataMat)
  if (norm.data) dataMat = t(scale(t(dataMat[features,, drop = F]), 
                                   center = T, scale = T))
  
  scores = colSums(model$finalModel$coefficients[-1] * dataMat) + model$finalModel$coefficients[1]
  return(scores)
}

scorePlot = function(mat, features = NULL, meta, 
                     class = "class", p_cl = c("AR", "Normal"), 
                     model = NULL, method = 1, thres = NA, norm.data = FALSE) {
  
  if (is.null(features)) features = rownames(mat)
  mat = mat[features,, drop = FALSE]
  mat1 = mat[, colnames(mat) %in% meta[meta[,class] %in% c("AR", "Normal"), "sample"]]
  mat2 = mat[, colnames(mat) %in% meta[meta[,class] == "STA", "sample"]]
  
  if (!is.null(model)) {
    scores1 = instaScore(mat1, model = model, norm.data = norm.data)
    scores2 = instaScore(mat2, model = model, norm.data = norm.data)
    scores = c(scores1, scores2)
    
  }
  
  meta = meta[colnames(mat),]
  data_xy = data.frame(scores = scores, 
                       class = meta[names(scores),class],
                       size = rep(1,length(scores)))
  pval = wilcox.test(scores[rownames(data_xy[data_xy$class == p_cl[1],])],
                     scores[rownames(data_xy[data_xy$class == p_cl[2],])])$p.value    
  cat("Wilcoxon test comparing ", p_cl[1], " to ", p_cl[2], "; p-value = ", pval, "\n", sep = "")
  
  data_xy$class = factor(data_xy$class, 
                         levels = c("AR", "Normal", "hSTA/mAR", "hSTA/mSTA", "mAR_eGFR", "mSTA_eGFR"))
  
  library("ggbeeswarm")
  col_pal = pal_nejm(alpha = 0.9)(8)[c(1, 4, 3, 2, 5, 6)]
  cols = c("AR" = col_pal[1], "Normal" = col_pal[2], "hSTA/mAR" = col_pal[3], 
           "hSTA/mSTA" = col_pal[4], "mAR_eGFR" = col_pal[5], "mSTA_eGFR" = col_pal[6])
  p = ggplot(data = data_xy, aes(x = scores, y = 0, fill = class)) +
    geom_quasirandom(aes(group = class),
                     method = 'smile',
                     # width = 0.3,
                     dodge.width = 0.01,
                     varwidth = T,
                     shape = 21,
                     color = "white", 
                     cex = 3,
                     bandwidth = .002,
                     alpha = 0.7,
                     nbins = 100,
                     show.legend = TRUE,
                     groupOnX = FALSE) +
    theme_bw() +
    scale_fill_manual(values = cols) +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          # axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          panel.background = element_rect()) +
    xlab("") +
    guides(size = "none")
  
  # if (all(c("AR", "Normal") %in% unique(meta[colnames(mat), class]))) {
  A.density <- density(subset(data_xy, class == p_cl[1])$scores, 
                       from = min(data_xy$scores), 
                       to = max(data_xy$scores), n = 2^10)
  B.density <- density(subset(data_xy, class == p_cl[2])$scores, 
                       from = min(data_xy$scores), 
                       to = max(data_xy$scores), n = 2^10)
  
  intersection.point <- A.density$x[which(diff((A.density$y - B.density$y) > 0) != 0) + 1]
  print(intersection.point)
  
  p = p #+ geom_vline(xintercept = intersection.point, color = "red", linetype = 2)
  
  ydens <- axis_canvas(p, axis = "x", coord_flip = FALSE) +
    geom_density(data = data_xy, aes(x = scores, fill = class),
                 alpha = 0.7, size = 0.2) + 
    # ggpubr::fill_palette(pal_nejm(alpha = 0.9)(8)[c(1, 4, 3, 2, 5, 6)])
    ggpubr::fill_palette(cols)
  p1 <- insert_xaxis_grob(p, ydens, grid::unit(.2, "null"), position = "top") 
  ggdraw(p1)
}


plotHeatmap = function(emat, sampleMeta, rowNames = FALSE, colCluster = TRUE, 
                       rowCluster = TRUE, classes = "class", title = "", 
                       rowFont = 10, fontsize = 10, bars = NA, colSplit = 2) {
  emat = emat[, colnames(emat) %in% rownames(sampleMeta)]
  sampleMeta = sampleMeta[colnames(emat),]
  annotation_col = sampleMeta[, classes, drop = FALSE]
  
  cl = sort(unique(unlist(sampleMeta[, classes])))
  cols = rep("", length(cl))
  names(cols) = cl
  
  # col_palette = pal_nejm(alpha = 0.8)(8)
  cols[names(cols) == "AR"] = "#BC3C29CC"
  cols[names(cols) == "Normal"] = "#20854ECC"
  cols[names(cols) == "STA"] = "#0072B5CC"
  cols[names(cols) == "mSTA"] = "#0072B5CC"
  cols[names(cols) == "mAR"] = "#E18727CC"
  cols[names(cols) == "hSTA/mSTA"] = "#0072B5CC"
  cols[names(cols) == "hSTA/mAR"] = "#E18727CC"
  
  cols[names(cols) == "BL"] = "#7876B1CC"
  cols[names(cols) == "AR+CAN"] = "#FFDC91CC"
  cols[names(cols) == "BL+ABMR"] = "#FFDC91CC"
  cols[names(cols) == "TCMR"] = "#E18727CC"
  cols[names(cols) == "Mixed"] = "#6F99ADCC"
  cols[names(cols) == "ABMR+TCMR"] = "#6F99ADCC"
  cols[names(cols) == "BL+CAN"] = "#7876B1CC"
  cols[names(cols) == "ABMR"] = "#EE4C97CC"
  
  ann_colors = lapply(classes, function(x) cols[sort(unique(sampleMeta[, x]))])
  names(ann_colors) = classes
  
  colors = colorRampPalette(rev(brewer.pal(n = 11, name = 'RdBu')))(100)
  
  if (class(bars) == "data.frame") {
    for (i in 1:length(bars)) {
      if (colnames(bars)[i] == "log2_FC") {
        cols.extra = diverge.color(start.color = "#008837", end.color = "#7b3294", 
                                   min.value = min(bars[[i]]*10), max.value = max(bars[[i]]*10), 
                                   mid.value = 0, mid.color = "white")
      } else if (colnames(bars)[i] == "correlation") {
        cols.extra = diverge.color(start.color = "#67a9cf", end.color = "#fc8d59", 
                                   min.value = min(bars[[i]]*10), max.value = max(bars[[i]]*10), 
                                   mid.value = 0, mid.color = "white")
      } else if (colnames(bars)[i] == "fold change") {
        cols.extra = diverge.color(start.color = "#008837", end.color = "#7b3294", 
                                   min.value = min(bars[[i]]*10), max.value = max(bars[[i]]*10), 
                                   mid.value = 1, mid.color = "white")
      } else {
        cols.extra = pal_material(palette = c("pink", "purple", "indigo", "cyan", "green", "lime",
                                              "amber", "deep-orange", "brown", "blue-grey")[i],
                                  n = length(bars[[i]]), alpha = 0.65)(length(bars[[i]]))
      }
      
      ann_colors = c(ann_colors,  list(cols.extra))
      names(ann_colors)[length(ann_colors)] = names(bars)[i]
    }
  }
  
  ematVars = sapply(rownames(emat), function(x) var(emat[x,]))
  if (any(ematVars == 0)) {
    message("Genes ", paste0(names(ematVars[ematVars == 0]), collapse = ","), " have zero variance")
    emat = emat[ematVars != 0, ]
  }
  
  
  breaks = seq(from = -3, to = 3, length.out = 100)
  breaks[length(breaks)] <- max(max(emat),max(breaks))
  breaks[1] <- min(min(emat),min(breaks))
  
  pheatmap(emat, 
           color = colors,
           scale = "row",
           breaks = breaks, 
           cluster_rows = rowCluster, 
           cluster_cols = colCluster, 
           show_colnames = FALSE, 
           show_rownames = rowNames, 
           border_color = NA, 
           clustering_method = "ward.D",
           cutree_rows = 2,
           cutree_cols =  colSplit,
           fontsize_row = rowFont,
           fontsize = fontsize,
           annotation_col = annotation_col,
           annotation_row = bars,
           annotation_colors = ann_colors,
           treeheight_row = 30,
           treeheight_col = 50,
           angle_col = "90",
           # cex = 0.9,
           silent = FALSE)
}

diverge.color <- function(start.color, end.color, min.value, max.value, 
                          mid.value = 0, mid.color = "ivory") {
  # based on ideas from Maureen Kennedy, Nick Povak, and Alina Cansler
  
  # creates a palette for the current session for a divergent-color
  # graphic with a non-symmetric range
  # "cuts" = the number of slices to be made in the range above and below "mid.value"
  
  ramp1 <- colorRampPalette(c(start.color, mid.color))
  ramp2 <- colorRampPalette(c(mid.color, end.color))
  
  # now specify the number of values on either side of "mid.value"
  
  max.breaks <- round(max.value - mid.value)
  # min.breaks <- round(mid.value - min.value)
  min.breaks <- ceiling(mid.value - min.value)
  
  num.breaks <- max(max.breaks,min.breaks)
  
  low.ramp <- ramp1(num.breaks)
  high.ramp <- ramp2(num.breaks)
  
  # now create a combined ramp from the higher values of "low.ramp" and 
  # the lower values of "high.ramp", with the longer one using all values 
  # high.ramp starts at 2 to avoid duplicating zero
  
  myColors <- c(low.ramp[(num.breaks - min.breaks):num.breaks], high.ramp[2:max.breaks])
  
  return(myColors)
}

# PCA Plot
plotPCA = function(emat, sampleMeta, title = "", interest = "class", 
                   alpha = 1, cex = 1, cexLab = 1.5, ellipse = TRUE) {
  sampleMeta = sampleMeta[rownames(sampleMeta) %in% colnames(emat), ]
  emat = emat[, rownames(sampleMeta)]
  classes = gsub("_[A-Z]", "", sampleMeta[, interest])
  names = factor(classes, levels = sort(unique(classes)))
  
  cols = rep("", length(classes))
  names(cols) = classes
  
  col_palette = pal_nejm(alpha = alpha)(8)
  cols = rainbow(length(unique(names)))[as.integer(names)]
  names(cols) = classes
  
  cols[names == "AR"] = "#BC3C29CC"
  cols[names == "Normal"] = "#20854ECC"
  cols[names == "STA"] = "#0072B5CC"
  cols[names == "mSTA"] = "#0072B5CC"
  cols[names == "mAR"] = "#E18727CC"
  cols[names == "hSTA/mSTA"] = "#0072B5CC"
  cols[names == "hSTA/mAR"] = "#E18727CC"
  
  cols[names == "BL"] = "#7876B1CC"
  cols[names == "AR+CAN"] = "#FFDC91CC"
  cols[names == "TCMR"] = "#E18727CC"
  cols[names == "Mixed"] = "#6F99ADCC"
  cols[names == "BL+CAN"] = "#7876B1CC"
  cols[names == "ABMR"] = "#EE4C97CC"
  cols[names == "BL+ABMR"] = "#FFDC91CC"
  cols[names == "TCMR"] = "#E18727CC"
  cols[names == "ABMR+TCMR"] = "#6F99ADCC"
  cols[names == "ABMR"] = "#EE4C97CC"
  
  
  ematVars = sapply(rownames(emat), function(x) var(emat[x,]))
  if (any(ematVars == 0)) {
    message("Genes ", paste0(names(ematVars[ematVars == 0]), collapse = ","), " have zero variance")
    emat = emat[ematVars != 0, ]
  }
  emat = emat[complete.cases(emat),]
  res.svd = svd(scale(t(emat)))
  
  ggbiplot(res.svd, ellipse = ellipse, circle = TRUE, 
           obs.scale = 1, var.scale = 1,
           var.axes = FALSE, labels = NULL, groups = classes, pointSize = cex, alpha = alpha) + 
    # geom_point(size = 10) +
    scale_size(range = c(0.8, cex)) +
    scale_colour_manual(name = interest, values = cols) +
    theme(text = element_text(size = 15)) +
    theme_bw() +
    guides(size = FALSE) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          plot.title = element_text(hjust = 0.5),
    ) + 
    # coord_fixed(1) +
    labs(title = title)
}

ggbiplot <- function(pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE, 
                     obs.scale = 1 - scale, var.scale = scale, 
                     groups = NULL, ellipse = FALSE, ellipse.prob = 0.68, 
                     labels = NULL, labels.size = 3, alpha = 1, 
                     var.axes = TRUE, 
                     circle = FALSE, circle.prob = 0.69, 
                     varname.size = 3, varname.adjust = 1.5, 
                     varname.abbrev = FALSE, pointSize = 1, ...) {
  library(ggplot2)
  library(plyr)
  library(scales)
  library(grid)
  
  stopifnot(length(choices) == 2)
  
  # Recover the SVD
  if (class(pcobj) == 'prcomp') {
    nobs.factor <- sqrt(nrow(pcobj$x) - 1)
    d <- pcobj$sdev
    u <- sweep(pcobj$x, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- pcobj$rotation
  } else if (inherits(pcobj, 'princomp')) {
    nobs.factor <- sqrt(pcobj$n.obs)
    d <- pcobj$sdev
    u <- sweep(pcobj$scores, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- pcobj$loadings
  } else if (inherits(pcobj, 'PCA')) {
    nobs.factor <- sqrt(nrow(pcobj$call$X))
    d <- unlist(sqrt(pcobj$eig)[1])
    u <- sweep(pcobj$ind$coord, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- sweep(pcobj$var$coord,2,sqrt(pcobj$eig[1:ncol(pcobj$var$coord),1]),FUN = "/")
  } else if (class(pcobj) == 'list') {
    nobs.factor <- sqrt(nrow(pcobj$u) - 1)
    d <- pcobj$d
    u <- sweep(pcobj$u, 2, 1 / nobs.factor, FUN = '*')
    # u <- sweep(u, 2, 1 / (d * nobs.factor), FUN = '*') #predict(pcobj)$x/nobs.factor
    v <- pcobj$v
    d.total <- sum(d^2)
  } else {
    stop('Expected a object of class prcomp, princomp, PCA, or lda')
  }
  
  # Scores
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[,choices], 2, d[choices]^obs.scale, FUN = '*'))
  
  # Directions
  v <- sweep(v, 2, d^var.scale, FUN = '*')
  df.v <- as.data.frame(v[, choices])
  
  names(df.u) <- c('xvar', 'yvar')
  names(df.v) <- names(df.u)
  
  if (pc.biplot) {
    df.u <- df.u * nobs.factor
  }
  
  # Scale the radius of the correlation circle so that it corresponds to 
  # a data ellipse for the standardized PC scores
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
  
  # Scale directions
  v.scale <- rowSums(v^2)
  df.v <- r * df.v / sqrt(max(v.scale))
  
  # Change the labels for the axes
  if (obs.scale == 0) {
    u.axis.labs <- paste('standardized PC', choices, sep = '')
  } else {
    u.axis.labs <- paste('PC', choices, sep = '')
  }
  
  # Append the proportion of explained variance to the axis labels
  u.axis.labs <- paste(u.axis.labs, 
                       sprintf('(%0.1f%% explained var.)', 
                               100 * d[choices]^2/sum(d^2)))
  
  # Score Labels
  if (!is.null(labels)) {
    df.u$labels <- labels
  }
  
  # Grouping variable
  if (!is.null(groups)) {
    df.u$groups <- groups
  }
  
  # Variable Names
  if (varname.abbrev) {
    df.v$varname <- abbreviate(rownames(v))
  } else {
    df.v$varname <- rownames(v)
  }
  
  # Variables for text label placement
  df.v$angle <- with(df.v, (180/pi) * atan(yvar / xvar))
  df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar)) / 2)
  
  # Base plot
  g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + 
    xlab(u.axis.labs[1]) + ylab(u.axis.labs[2]) #+ coord_fixed(1)
  
  if (var.axes) {
    # Draw circle
    if (circle) 
    {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
      circle <- data.frame(xvar = r * cos(theta), yvar = r * sin(theta))
      g <- g + geom_path(data = circle, color = muted('white'), 
                         size = 1/2, alpha = 1/3)
    }
    
    # Draw directions
    g <- g +
      geom_segment(data = df.v,
                   aes(x = 0, y = 0, xend = xvar, yend = yvar),
                   arrow = arrow(length = unit(1/2, 'picas')), 
                   color = muted('red'))
  }
  
  # Draw either labels or points
  if (!is.null(df.u$labels)) {
    if (!is.null(df.u$groups)) {
      g <- g + geom_text(aes(label = labels, color = groups), 
                         size = labels.size)
    } else {
      g <- g + geom_text(aes(label = labels), size = labels.size)      
    }
  } else {
    if (!is.null(df.u$groups)) {
      g <- g + geom_point(aes(color = groups, size = pointSize), alpha = alpha)
    } else {
      g <- g + geom_point(aes(size = pointSize), alpha = alpha)      
    }
  }
  
  # Overlay a concentration ellipse if there are groups
  if (!is.null(df.u$groups) && ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))
    
    ell <- ddply(df.u, 'groups', function(x) {
      if (nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x$xvar, x$yvar))
      mu <- c(mean(x$xvar), mean(x$yvar))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = '+'), 
                 groups = x$groups[1])
    })
    names(ell)[1:2] <- c('xvar', 'yvar')
    g <- g + geom_path(data = ell, aes(color = groups, group = groups))
  }
  
  # Label the variable axes
  if (var.axes) {
    g <- g + 
      geom_text(data = df.v, 
                aes(label = varname, x = xvar, y = yvar, 
                    angle = angle, hjust = hjust), 
                color = 'darkred', size = varname.size)
  }
  return(g)
}




# UMAP Plot
plotUMAP = function(emat, sampleMeta, interest = "class", 
                    verbose = FALSE, n = 15, alpha = 0.75, pointSize = 1, title = "") {
  library("uwot")
  sampleMeta = sampleMeta[sampleMeta$sample %in% colnames(emat), ]
  
  set.seed(309)
  umap_data = uwot::umap(t(emat), n_neighbors = n, verbose = verbose, scale = TRUE)
  
  classes = sampleMeta[colnames(emat), interest]
  names = factor(classes, levels = sort(unique(classes)))
  # col_palette = pal_nejm(alpha = alpha)(8)
  cols = rainbow(length(unique(names)), alpha = alpha)[as.integer(names)]
  # cols = rep("", length(classes))
  names(cols) = classes
  
  cols[names == "AR"] = "#BC3C29CC"
  cols[names == "Normal"] = "#20854ECC"
  cols[names == "STA"] = "#0072B5CC"
  cols[names == "mSTA"] = "#0072B5CC"
  cols[names == "mAR"] = "#E18727CC"
  cols[names == "hSTA/mSTA"] = "#0072B5CC"
  cols[names == "hSTA/mAR"] = "#E18727CC"
  
  cols[names == "BL"] = "#7876B1CC"
  cols[names == "AR+CAN"] = "#FFDC91CC"
  cols[names == "TCMR"] = "#E18727CC"
  cols[names == "Mixed"] = "#6F99ADCC"
  cols[names == "BL+CAN"] = "#7876B1CC"
  cols[names == "ABMR"] = "#EE4C97CC"
  cols[names == "BL+ABMR"] = "#FFDC91CC"
  cols[names == "TCMR"] = "#E18727CC"
  cols[names == "ABMR+TCMR"] = "#6F99ADCC"
  cols[names == "ABMR"] = "#EE4C97CC"
  
  ggplot(data = as.data.frame(umap_data), aes(x = V1, y = V2)) + 
    geom_point(aes(color = names), size = pointSize) +
    xlab("UMAP1") + ylab("UMAP2") +
    scale_colour_manual(name = interest, values = cols) +
    theme(text = element_text(size = 15)) +
    theme_bw() +
    guides(size = FALSE) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          plot.title = element_text(hjust = 0.5),
    ) + 
    labs(title = title)
  
  
}

plotROC = function(test, target, leg = "") {
  library("ROCR")

  # Negative - Positive classes are in order (e.g. 0 < 1, -1 < 1, 'a' < 'b', FALSE < TRUE)
  
  predictions = as.vector(test)
  pred = prediction(predictions = predictions, 
                    labels = target)
  
  perf_AUC = performance(pred, "auc") #Calculate the AUC value
  AUC = perf_AUC@y.values[[1]]
  perf_AUCPR = performance(pred, "aucpr") #Calculate the AUC value
  AUCPR = perf_AUCPR@y.values[[1]]
  
  perf_ROC = performance(pred,"tpr","fpr") #plot the actual ROC curve
  perf_RR = performance(pred,"prec","rec") #plot the actual PR curve
  # perf_ROC = performance(pred,"sens", "spec") #plot the actual ROC curve
  
  par(mfrow = c(1,2))
  plot(perf_ROC, #main = paste("ROC plot", leg), 
       col = "red", lwd = 3)
  segments(x0 = 0, y0 = 0, x1 = 1, y1 = 1)
  # abline(a = 0, b = 1)
  # plot(x=0:1, y=0:1, type="l")
  text(0.5, 0.5, paste("AUC = ",format(AUC, digits = 3, scientific = FALSE)))
  plot(perf_RR, #xlim = c(0,1), ylim = c(0,1),  #main = paste("ROC plot", leg), 
       col = "blue", lwd = 3)
  text((0.5+min(perf_RR@x.values[[1]][!is.nan(perf_RR@y.values[[1]])])/2), 0.77, paste("AUCPR = ",format(AUCPR, digits = 3, scientific = FALSE)))
  par(mfrow = c(1,1))
  # return(perf_RR)
}

