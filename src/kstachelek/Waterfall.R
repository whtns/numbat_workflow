library(akima)
library(ape) # MST
library(Hmisc) #extrapolate
library(dendextend)
library(FactoMineR)
#library(RHmm) 
library(ggplot2)
library(caTools)
library(signal) #library(KEGGREST)
library(pheatmap)
library(scatterplot3d)
#source("src/KEGGREST_tmp.R")
#library("RDAVIDWebService")
library(GO.db)
#library(org.Mm.eg.db)
#library(hash)
library(MASS)  # standard, no need to install
library(class)  # standard, no need to install
#library(WGCNA)
library(flashClust)
library(parallel)
library(snow)
#library(rgl)
library(RColorBrewer)

Waterfall.version = "Waterfall version 2.0"

ScaleTPM = function (TPM){
  .d = apply(TPM, 2, function(X) return(X/sum(X))) * 1000000
  return (as.data.frame(.d))
}

raw2FPKM = function (mat,start,end){
  len_vect = (end - start)
  FPKM = as.data.frame(apply(mat,2,function(cell){((cell*1e+03)/len_vect)/(sum(cell)/1e+06)}))
  return(FPKM)
  
}

raw2TPM = function (mat,start,end){
  len_vect = (end - start)
  FPKM = apply(mat,2,function(cell){((cell*1e+03)/len_vect)/(sum(cell)/1e+06)})
  TPM = ScaleTPM(FPKM)
  return(TPM)
  
}

minDist2D = function(v){
  return(apply(as.matrix(dist(v)),1,function(v){min(v[v>0])}))
}

Waterfall.Cluster = function(TPM_Data, nClusters=NULL, simplecolors = F,colorchoices = NULL,cell.cols = NULL){
  "Heirarchical clustering. plots various clusters, returns colors based on nClusters"
  dat = as.data.frame(TPM_Data)
  d_mat <- dist(t(dat));
  hcobj <- stats::hclust(d_mat, method = "ward", members=NULL)
  
  drawFit = function(h_in,k_out){
    
    edgeIndex = max(which(k_out == max(k_out)))
    A = min(k_out)
    D = (max(k_out) - A)
    seek = D*exp(-1)
    Tau = h_in[which.min(abs(k_out-(seek)))] - h_in[edgeIndex]
    
    k_est = c()
    for (h in h_in){
      ret = D*exp(-(h/Tau))+A
      k_est = append(k_est,ret)
    }
    k_est = c(rep(max(k_out),edgeIndex),k_est[1:(length(k_est)-edgeIndex)])
    lines(h_in,k_est,col="darkgreen",lty=2)
    legend('topright', c("k output","k model","best h"), lty=c(1,2,1), bty='n', cex=.75, col= c("black","green","red")) 

  }
  
  #Determine best k (assumes exponential decay)
  if (is.null(nClusters)){
    unity = max(hcobj$height)
    h_candidates = seq(0,unity,100)
    k_output = c()
    for (h in h_candidates){
      k = max(cutree(hcobj,h=h))
      k_output = append(k_output,k)
    }
   
    hIndex = (max(k_output)-min(k_output))*(exp(-3))
    best_h = h_candidates[which.min(abs(k_output-hIndex))]
    plot(h_candidates,k_output,type="l",main = "h-seek")
    abline(v = best_h,col="red")
    drawFit(h_candidates,k_output)
    
    br = cutree(hcobj,h=best_h)
  }else{
    br = cutree(hcobj,k=nClusters)
  }

  if (is.null(colorchoices)) colorchoices = c("#4882C3", "#F26A6A", "#13751B", "#FF6A00", "#E2CF00", "#980B0B", rgb(0,0,1), rgb(1,0,1), rgb(1,1,0), rgb(0,1,1), rgb(.5,0,0), rgb(0,.5,0), rgb(0,0,.5), rgb(.5,0,.5), rgb(.5,.5,0), rgb(0,.5,.5), rgb(0,0,0))
  if(simplecolors) colorchoices = colors();
  brcols = colorchoices[br]
  
  .dend = as.dendrogram(hcobj)
  labels_colors(.dend) = brcols[order.dendrogram(.dend)]
  if(!is.null(cell.cols)) labels_colors(.dend) = cell.cols[order.dendrogram(.dend)]
  plot(.dend)
  
  out = cbind( brcols, br)
  rownames(out) = names(br)
  return(list(out, dend))
}


rotate2D = function(M, angle){
  alpha = angle * pi/180
  #rotation matrix
  rotm <- matrix(c(cos(alpha),sin(alpha),-sin(alpha),cos(alpha)),ncol=2)
  #shift, rotate, shift back
  M2 <- t(rotm %*% (
    t(M)-c(M[1,1],M[1,2])
  )+c(M[1,1],M[1,2]))
}


Waterfall.Pseudotime = function( TPM_Data, angle=0, cols = "green", plot_title="PCA plot", seed=5, nclust=6, scalePCA = FALSE, label=FALSE,lines = FALSE,invert_pt = FALSE,mst=TRUE,rndY = FALSE, threeD_plot=FALSE) {
  "Returns a vector, 'pseudotime', with cell names and a PT value from 0-1"
  print(Waterfall.version)
  
  "this makes a PCA plot using all of the data supplied. \n Data should be rows=genes, cols= cells"
  dat = as.data.frame(TPM_Data)
  pca2 = PCA(t(dat), graph = F, scale.unit = scalePCA)
  scatterplot3d(pca2$ind$coord[,1],pca2$ind$coord[,3],pca2$ind$coord[,2], main="3D Scatterplot", type="h", color = cols)
  
  if (threeD_plot){
    scatterplot3d(pca2$ind$coord[,1],pca2$ind$coord[,3],pca2$ind$coord[,2], main="3D Scatterplot", type="p", col = cols)
  }
  
  pca2.coords = pca2$ind$coord[,1:2]
  pca2.coords = rotate2D(pca2.coords, angle)
  
  
  if (invert_pt){
    shift_val = min(pca2.coords[,1])
    pca2.coords[,1] = (pca2.coords[,1] - min(pca2.coords[,1]))
    pca2.coords[,1] = max((pca2.coords[,1])) - pca2.coords[,1]
    pca2.coords[,1] = pca2.coords[,1] + shift_val
  } 
  
  plot( pca2.coords, col = cols, main=plot_title, cex=2, pch=20, bty="n", asp = 1)
  if (label) text(pca2.coords, labels=colnames(dat), col = cols, pos=1, cex = .6) #label each point with the cell name?

  #points(pca2.coords[names(g1),],cex=1.75,pch = 5)
  #points(pca2.coords[names(g2),],cex=1.75,pch = 3)
  
  set.seed(seed)
  km <- kmeans(pca2.coords,nclust)
  km.centers <- km$centers
  km.centers
  km.centers <- km.centers[order(km.centers[,1]),]
  rownames(km.centers) = 1:length(rownames(km.centers))
  #points(km.centers, col="black")
  
  #14.) Make minimum spanning tree (including S0, SA, Outliers)
  m = mst(dist(km.centers))
  w <- which(m!=0)
  i <- as.vector(row(m))[w]
  j <- as.vector(col(m))[w]
  
  #segments( km.centers[i,1], km.centers[i,2], km.centers[j,1], km.centers[j,2], col="#FF0000" ,lwd=2 )
  
  #calculate pseudotime for each cell
  euc.dist = function(x1, x2) ((x1[1]-x2[1])**2 + (x1[2]-x2[2])**2) **.5
  
  #TODO: modify lower/upper bound computation to allow for expansion closer to min/max vals along 1st component
  range = abs(max(pca2.coords[,1]) - min(pca2.coords[,1])) * 1.1
  lower_bound = max(pca2.coords[,1]) - range
  upper_bound = min(pca2.coords[,1]) + range
  
  myrange = seq(lower_bound,upper_bound,range /1000 )

  myrange.len = length(myrange)
  
  ypt = as.numeric(pca2.coords[,2])
  xpt = as.numeric(pca2.coords[,1])
  y.loess <-loess(ypt~xpt, span = .75)
  ypt.predict <-predict(y.loess,data.frame(xpt = xpt,se=T))
  
  km = km.centers
  if(mst){
    ax = as.data.frame(approxExtrap(km.centers[,1], km.centers[,2], xout =  myrange))
  }else{
    ax = as.data.frame(approxExtrap(xpt, ypt.predict, xout =  myrange))
  }  

  points(ax,  col = "gold", cex=.5) #plot the axis (must have pca plot already to se3)
  ax$pseudotime = 0
  for (i in 2:nrow(ax)){
    dis = euc.dist( ax[i-1, 1:2], ax[i, 1:2]) / 100 #scaling
    ax$pseudotime[i] = dis + ax$pseudotime[i-1]
  }
  
  #Calculate PT for each cell
  #for each cell, find closest point on axis
  pseudotime = c()
  pseudotime[colnames(dat)]  = 0
  distance = pseudotime #same structure as pt
  
  
  for (cellname in colnames(dat) ){
    coords = pca2.coords[cellname,]
    i = which.min(apply(ax[,1:2], 1, function(x1) euc.dist(x1, coords) ) )
    pt = ax$pseudotime[i]
    pseudotime[cellname]  = pt
    if(lines) lines( rbind(ax[i,1:2], coords), col="red")
    
    mysign = (coords[2] > ax[i,2]) * 2 -1
    
    distance2Axis = mysign * euc.dist(coords, ax[i, 1:2])
    distance[cellname] = distance2Axis
    
  }
  
  pseudotime = pseudotime - min(pseudotime)
  pseudotime = pseudotime * 100 / max(pseudotime)

  if (rndY) distance = sample(seq(0,100,.01),length(distance),replace = TRUE)
  
  output = (cbind(pseudotime, distance))[order(pseudotime),]
  
  plot(output, col = cols[rownames(output)], cex=3, pch=20, main = paste0("Pseudotime", plot_title),xlab = "Pseudotime" )
  
  return(output)
}

Waterfall.PCA = function(TPM_Data, ind_sup_vect = NULL, cols = NULL, plot_title="PCA plot", twoDcomp = c(1,2), threeDcomp = c(1,2,3),
                         scalePCA = FALSE, label=FALSE, twoD_plot = TRUE, threeD_plot=FALSE, pdf_out = NULL, multiplex_plot=FALSE){

  if ("Pseudotime" %in% rownames(TPM_Data)){ 
    TPM_Data = TPM_Data[-which(rownames(TPM_Data) %in% "Pseudotime"),]
  }
  

  dat = as.data.frame(TPM_Data)
  
  if(is.null(cols)){ 
    cols = rep("black",ncol(dat))
  }
  
  if(missing(ind_sup_vect)){
    pca2 = PCA(t(dat), graph = F, scale.unit = scalePCA, ncp = 20)  
  } else{
    pca2 = PCA(t(dat), graph = F, scale.unit = scalePCA, ncp = 20, ind.sup = ind_sup_vect)  
  }
  
#  pca2 = PCA(t(dat), graph = F, scale.unit = scalePCA, ncp = 20, ind.sup = ind_sup_vect)
  twoD.lab1 = paste0("PC",twoDcomp[1]," (",round(pca2$eig[twoDcomp[1],2],2),"%)")
  twoD.lab2 = paste0("PC",twoDcomp[2]," (",round(pca2$eig[twoDcomp[2],2],2),"%)")
  
  threeD.lab1 = paste0("PC",threeDcomp[1]," (",round(pca2$eig[threeDcomp[1],2],2),"%)")
  threeD.lab2 = paste0("PC",threeDcomp[2]," (",round(pca2$eig[threeDcomp[2],2],2),"%)")
  threeD.lab3 = paste0("PC",threeDcomp[3]," (",round(pca2$eig[threeDcomp[3],2],2),"%)")
  
  if (twoD_plot){
  pca2.coords = pca2$ind$coord[,twoDcomp]
  pca2.supp.coords = pca2$ind.sup$coord[,twoDcomp]
  plot(pca2.coords,xlab = twoDcomp[1],ylab = twoDcomp[2], col = cols, main=plot_title, cex=0.5, pch=20, bty="n", asp = 1)
  points(pca2.supp.coords,xlab = twoDcomp[1],ylab = twoDcomp[2], col = cols, main=plot_title, cex=0.5, pch=20, bty="n", asp = 1)
  }
  
  if (threeD_plot){
   
    scatterplot3d(pca2$ind$coord[,threeDcomp[1]],pca2$ind$coord[,threeDcomp[3]],pca2$ind$coord[,threeDcomp[2]],
           xlab = threeD.lab1,ylab = threeD.lab2,zlab = threeD.lab3,main="3D Scatterplot", type="h", color = cols)
  }
  
  if (multiplex_plot){
   
    pdf(pdf_out)
    par(mfrow = c(2,2))  
    for (i in twoDcomp[1]:(20-1)){
      j=i+1
      for (j in j:20){
        c_pair = c(i, j)
        plot_title = paste("PC", i, "&", j)
        pca2.coords = pca2$ind$coord[,c_pair]
        plot(pca2.coords,xlab = i,ylab = j, col = cols, main=plot_title, cex=0.5, pch=20, bty="n", asp = 1)  
        points(pca2.supp.coords,xlab = twoDcomp[1],ylab = twoDcomp[2], col = cols, main=plot_title, cex=0.5, pch=20, bty="n", asp = 1)
      }
    }
  dev.off()  
  }
  
  pc.cor = list()
  for (i in 1:(ncol(pca2$var$cor))){
    tmp = pca2$var$cor[,i]
    p = tmp[tmp>0 & tmp<Inf & !is.na(tmp)]; p = c(p[order(p,decreasing = TRUE)])
    n = tmp[tmp<0 & tmp>-Inf & !is.na(tmp)]; n = c(n[order(n)])
    
    pc.cor[[paste0("pc",as.character(i))]] = list(p=p,n=n)
  }
  
  
  return(list(pcaN = pca2$ind$coord, pc.cor = pc.cor))
}

PTgeneplot = function(genename, TPM_data, span=.75, hmm = FALSE, s = 2, unit = NULL , gene_txt = NULL,minimum=NULL, beta=1,apply_to_graph=FALSE,...){
  
  for (i in 1:length(genename)){
    if(genename[i] %notin% rownames(TPM_data)){
      if (is.na(genename[i]) == TRUE){
        break()
      } else{
      warning(sprintf("The genename %s is not a gene. The results will exclude the gene. Go back and check the spelling",genename[i]))
      genename = genename[-i]
      }
    }
  }
  if (length(genename) == 1){
    par(mfrow = c(1,1))
  } else if (length(genename) == 2){
    par(mfrow = c(2,1))
  } else if (length(genename) <= 4 && length(genename) >= 3){
    par(mfrow = c(2,2))
  } else if (length(genename) >= 5 && length(genename) <= 9){
    par(mfrow = c(3,3))
  } else{
    par(mfrow = c(4,4))
  }
  for (i in 1:length(genename)){
  
  PT = as.numeric(TPM_data["Pseudotime",])
  exp = as.numeric(TPM_data[genename[i],])
  if(apply_to_graph) exp = exp^beta
  
  ypt = exp
  xpt = as.numeric(PT)
  y.loess <-loess(ypt~xpt, span = span)
  ypt.predict <-predict(y.loess,data.frame(x=xpt),se=T)
  
  if (is.null(gene_txt)){ 
    main = genename[i]
  }else{
    main = gene_txt
  }
  plot(PT, exp , pch=20, main=main,...)
  
  
  fit = ypt.predict$fit
  margin = ypt.predict$se.fit
  SCALAR = 1
  
  lines(xpt,fit,lwd=1)
  lines(xpt, ypt.predict$fit - SCALAR*margin, lty = 2)
  lines(xpt, ypt.predict$fit + SCALAR*margin, lty = 2)
  
  y.poly = c((fit+SCALAR*margin)[1:length(fit)], (fit-SCALAR*margin)[length(fit):1])
  
  x.poly = c(PT, PT[order(PT,decreasing = TRUE)])
  
  polygon(x.poly, y.poly,col="#00009933",border = NA)
  
  if (hmm){
    onoff = state.foo(genename[i],TPM_data,s,minimum=minimum,beta=beta)
    brewercolor = (brewer.pal(s,"YlOrRd"))
    brewercolor = brewercolor[s:1]
    brewercolor = brewercolor[-2]
    brewercolor = append(brewercolor,"#000000")
    clrs = brewercolor[s:1]
    #clrs = colfunc(s)
    #clrs = c("#000000","#FFEE00")
    
    offset = min(exp) + (min(exp)-max(exp))/30
    hmmX = seq(.5,99.5,101/length(onoff))
    hmmY = rep(offset,length(hmmX))
    
    for (pti in 2:length(hmmX)){
      x1 = hmmX[pti-1]
      y1 = hmmY[pti-1]
      x2 = hmmX[pti]
      y2 = hmmY[pti]
      lines(c(x1,x2),c(y1,y2),col=clrs[onoff[pti-1]],lwd = 9)
    } 
  }
  }
}

geneExpPTCor = function(pca_coords, cor_genes){
  print("not done")
  
}

###PT gene plots
PTgeneplot.trend = function(g, p, span=.75, est_i = NULL, PT_points = 101, ...){
  
  if (is.character(g)){
    if (!g %in% rownames(p)) return()
    time = as.numeric(p["Pseudotime",])
    exp = as.numeric(p[g,])
  }else if (is.vector(g)){
    time = p
    exp = g
  }
  
  ypt = as.numeric(exp)
  xpt = as.numeric(time)
  y.loess = loess(ypt~xpt, span = span)
  
  if (is.null(est_i)){
    ypt.predict = predict(y.loess,data.frame(xpt = seq(0,100,(101/PT_points))),se=T)
    
  }else{
    ypt.predict <-predict(y.loess,est_i,se=T)
  }
  
  out = ypt.predict$fit
  out[which(is.na(out))] = out[which(is.na(out))-1]
  
  return (out)
  
}


pluralPTgeneplots = function(gene_list, TPM_data, span=.75, colors = NULL){
  
  mmNorm = function(v,new_max = 100){
    v = (v-min(v))
    v = as.numeric(v/max(v)) * new_max
    return (v)
  }
  
  plot(c(-100,200),c(-100,200),ylim=c(0,100),xlim=c(0,100),ylab = "- Dynamic Time Warp -",xlab = "Pseudotime",yaxt='n')
  for (i in 1:length(gene_list)){
    
    genename = gene_list[i]
    ypt = as.numeric(TPM_data[genename,])
    xpt = as.numeric(TPM_data["Pseudotime",])
    
    y.loess = loess(ypt~xpt, span = span)
    ypt.predict = predict(y.loess,data.frame(x=xpt),se=T)
    
    
    fit = mmNorm(ypt.predict$fit)
    color = i
    if(!is.null(colors)) color = colors[i]
      
    lines(xpt,fit,lwd=2,col=color)
    
  }
}

multiplePTgeneplot = function(genename, TPM_data, PT, sub1, sub2, span=.75, ...){
  
  if(genename %notin% rownames(TPM_data)) return();
  
  sub1 = names(PT[names(PT)%in%sub1])
  sub2 = names(PT[names(PT)%in%sub2])
  
  exp = as.numeric(TPM_data[genename,])
  nms = names(PT)
  PT = as.numeric(PT); names(PT) = nms

  plot(PT, exp , pch=20, ...)
  
  ypt1 = as.numeric(TPM_data[genename,sub1])
  xpt1 = as.numeric(PT[sub1])
  y.loess1 = loess(ypt1~xpt1, span = span)
  ypt.predict1 = predict(y.loess1,data.frame(x=xpt1),se=T)
  
  ypt2 = as.numeric(TPM_data[genename,sub2])
  xpt2 = as.numeric(PT[sub2])
  y.loess2 = loess(ypt2~xpt2, span = span)
  ypt.predict2 = predict(y.loess2,data.frame(x=xpt2),se=T)
  
  lines(xpt1,ypt.predict1$fit,lwd=2,col="black")
  lines(xpt2,ypt.predict2$fit,lwd=2,col="red")
  
}

state.foo <-function(g, p, s=2,minimum = NULL, unit=NULL,beta=1){
  
  
  find_unit = function(PT){
    PT = as.numeric(PT)
    for (i in seq(0.1, length(PT),0.01)){
      bin_count = as.numeric(round(length(PT)/i))
      freq = hist(PT,breaks = bin_count, plot = FALSE)
      counts = freq$counts
      if(min(counts)!=0){
        return(as.numeric(round(i/2.5)))
        break
      }
    }
    warning("Could not bit to unity")
  }
  
  
  onoff_determine.foo <-function(exp,onoff,time,s){
    obs = unlist(lapply(split(exp, t(round(time/unit))),mean))
    ResFit = HMMFit(obs, nStates=s)# Baum-Welch
    state = viterbi(ResFit,obs)
    names(onoff) = time.unit
    onoff[match(as.character(names(obs)),names(onoff))] = state[[1]]
    
    #finds nearest unveiled value from HMM across vector indices 
    i = 1
    unattended_bins = TRUE 
    while (unattended_bins) {
      onoff[which(onoff==.5)] = onoff[which(onoff==.5)+i]
      onoff[which(onoff==.5)] = onoff[which(onoff==.5)-i]
      if (length(onoff[onoff==.5]) == 0) unattended_bins = FALSE
      i = i + 1
        
    }
    
    state_order <-order(ResFit$HMM$distribution$mean)
    
    onoff <-sapply(onoff,function(X){state_order[X]})
    return(onoff)
    
  }
  if (is.character(g)){
    time = as.numeric(p["Pseudotime",])
    exp = as.numeric(p[g,])^beta
  }else if (is.vector(g)){
    time = p
    exp = g^beta
  }
  

  
  if(is.null(unit)) unit = find_unit(time)
  if(is.null(minimum)) minimum = min(exp)
  time.unit = 0:max(round((time-min(time))/unit))
  
  onoff = rep(.5,length(time.unit))
  endpoint = length(exp)
  onethird = round(length(exp)*1/3)
  twothird = round(length(exp)*2/3)
  
  if (mean(exp[1:onethird]) > minimum & mean(exp[(onethird+1):twothird]) > minimum & mean(exp[(twothird+1):endpoint]) > minimum ){
    
    onoff = onoff_determine.foo(exp,onoff,time,s)
    
    s.new = 3
    state_curious = TRUE
    while ((s.new < 4) & (state_curious)){ # Loop through different s vals. If HMM raster still suspect, attempt using new s vals 
      if (all(onoff[2:length(onoff)]!=onoff[1:length(onoff)-1])) # states oscilate
      {
        onoff <-rep(.5,length(time.unit))
        onoff <-onoff_determine.foo(exp,onoff,time,s.new)
        if (!max(onoff)==1) onoff = round(onoff/max(onoff)*2)
        s.new = s.new + 1
        
      }else{
        state_curious = FALSE
      }
    }
    
  }else{
    onoff <- rep(1,length(time.unit))
    
  }


  onoff <-as.integer(ceiling(onoff))

  return(onoff)
}
  

distMap1D = function(x){
  dst = c()
  for (i in 1:length(x)){
    min_dst = min(abs(x[-i] - x[i]))
    dst = append(dst,min_dst)
  }
  return(dst)
}

##################################################################
##################################################################
##################################################################
################### GENE ONTOLOGY ANALYSIS FXNs  #################
##################################################################
##################################################################


#finds genes associated with GO term
loadGenesFromGo = function(process_name,exclusions = c(),min_genes=0){
  #process_name: process name as string
  #exclusions: vector of GO Evidence codes deemed inappropriate for the study
  #min_genes: minimum gene count 
  
  pid = as.character(term_ids[process_name]) #process ID
  egs = get(pid, org.Mm.egGO2ALLEGS) #
  genes = unlist(mget(egs,org.Mm.egSYMBOL))
  evidence = names(egs)
  rmi = which(evidence %in% exclusions) #find genes with unreliable evidence
  if (length(rmi) > 0) genes = genes[-rmi] 
  genes = unique(genes[duplicated(genes)]) #remove genes with multiple evidence codes 
  if (length(genes) < min_genes) return(NULL)
  return (genes)
}

plotGo.sin = function(process_name=NULL, TPMdata,exclusions = c(),min_genes = 0){
  if(is.null(process_name)) return(NULL)
  .loadGenes = loadGenesFromGo(process_name,exclusions=exclusions,min_genes=min_genes)
  if (is.null(.loadGenes)) return(NULL)
  
  Genes = .loadGenes[.loadGenes %in% rownames(TPMdata)]
  Desc = .loadGenes[.loadGenes %in% rownames(TPMdata)]
  
  SinScaledData = sin_scale2.foo(TPMdata[Genes,])
  
  .summed = apply(SinScaledData, 2, function(column) mean(column))
  
  xpt = as.numeric(TPMdata["Pseudotime",]) 
  ypt = as.numeric(.summed)
  #plot(xpt,ypt)
  #print(xpt)
  
  y.loess <-loess( ypt~xpt, span = .75)

  ypt.predict <-predict(y.loess,data.frame(x=xpt),se=T)
  ypt_fit = ypt.predict$fit
  ypt_fit = ypt_fit - min(ypt_fit)
  max = max(ypt_fit)
  if (max ==0){
  max = 1
  }
  ypt_fit = ypt_fit / max
  #print("10")
  return( list("x" = xpt, "y" = ypt_fit, "genes" = Genes, "ngenes" = length(Genes)))
  
  
}


MakeGoPlots.Sin = function(process_name, dataset, exclusions = c(),min_genes = 0,colors = NULL,plt = TRUE, main = NULL,give_lab = TRUE){
  pname = as.character(process_name)
  fit = list ( "x"= as.numeric(dataset["Pseudotime",]),"y"= rep(0, ncol(dataset)),"Genes"=c(),"ngenes"=0 )
  try( expr = (fit= plotGo.sin(process_name = pname, dataset, exclusions=exclusions, min_genes=min_genes)), silent = T)
  
  if (give_lab){
    if (is.null(main)){
      main = process_name
    }
    else (main = main)
  }

  
  if(plt){
    x_label = "Pseudotime"
    y_label = "Process Activity"
    if (is.null(colors)){
      plot(fit$x,fit$y,main = main,ylim=c(0,1),xlab = x_label,ylab = y_label)
    }else{
      plot(fit$x,fit$y,main= main,col = colors,ylim=c(0,1),xlab = x_label,ylab = y_label)
    }
    
  }
  return(fit) 
}

CreateGOTable = function(tpm_data,out_file="GO_activity.txt") {
  
  PT = tpm_data["Pseudotime",]
  a = matrix(PT,nrow = 1,ncol = length(PT))
  rn= c("Pseudotime")
  #for (i in 1:length(trms)){ ## for server
  for (i in 1:100){ ##for tests
    pname = as.character(trms[i])
    print(pname)
    rn = c(rn, pname)
    fit = list ( "x"= as.numeric(tpm_data["Pseudotime",]),"y"= rep(0, ncol(tpm_data)),"Genes"=c(),"ngenes"=0 )
    try( expr = (fit= plotGo.sin(process_name= pname, tpm_data, exclusions=c(), min_genes=0)), silent = T)
    enrichment=  as.numeric(fit$y)
    a = rbind(a, enrichment)

  }
  colnames(a) = colnames(tpm_data)
  rownames(a) = rn

  write.table(a, file = out_file)
  
}


getIntersectingGenesFromGo = function(geneList, pval_thresh = 0.01, write_csv = FALSE, csv_path = "~/Desktop/GOannotation.csv"){
  david = DAVIDWebService(email="pujadpat@usc.edu", url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
  ent = AnnotationDbi::select(org.Mm.eg.db, geneList, "ENTREZID", "SYMBOL")[,2]
  result = addList(david, ent, idType="ENTREZ_GENE_ID", listName="ent", listType="Gene")
  annochart = getFunctionalAnnotationChart(david)
  go_terms=c()
  for (i in 1:dim(annochart)[1])
  {
    if (annochart[i,5] <= pval_thresh & (length(grep("GOTERM", annochart[i,1])) == 1)){
      go_terms <- append(go_terms, annochart[i,2])
    }
  }
  intersect_genes = hash()
  max_length = 0
  for (i in 1:length(go_terms)){
    pval = annochart[(annochart[,2]==go_terms[i]),5]
    pid = substr(go_terms[i],0,10)
    try ({
      egs = get(pid, org.Mm.egGO2ALLEGS)
      genes = unlist(mget(egs,org.Mm.egSYMBOL))
      genes = unique(genes[duplicated(genes)])
      if (length(genes) >= 1){
        if (length(intersect(genes, geneList))>=1){
          term = substring(go_terms[i], 12)
          .set(intersect_genes, keys = term, values =c(pval, intersect(genes, geneList)))
          if (length(intersect(genes, geneList))>= max_length){
            max_length = length(intersect(genes, geneList))
          }
        }
      } })
  }
  if (write_csv)
  {
    gene_table = matrix("",nrow = length(intersect_genes),ncol = max_length+1)
    rownames(gene_table) = names(intersect_genes)
    for (i in 1:length(intersect_genes)){
      genes = intersect_genes[[names(intersect_genes)[i]]]
      for (j in 1:length(genes)){
        gene_table[i,j]=genes[j]
      }
    }
    gene_table = gene_table[order(as.double(gene_table[,1])),]
    write.csv(gene_table, file = csv_path ,row.names=TRUE)
  }
  return (intersect_genes)
  # use intersect_genes[["term"]] or values(intersect_genes, keys="term") to pull up genes for specific go term
}

corGenes = function(mat, pt,method="pearson"){
  if ("Pseudotime" %in% rownames(mat)){ mat = mat[-which(rownames(mat) %in% "Pseudotime"),]}
  
  pt = as.numeric(pt)
  c = apply(mat,1,function(v){cor(v,pt,method=method)})
  pos = c[c>0]; pos = pos[order(pos,decreasing = TRUE)]
  neg = c[c<0]; neg = neg[order(neg)]
  
  return(list(p=pos,n=neg))
  
}


##################################################################
##################################################################
##################################################################
################### KEGG ANALYSIS FXNs  ##########################
##################################################################
##################################################################
Waterfall.go.CleanData = function(dat, i = 4, cutoff=10){ #returns a TPM_data type array
  # takes data, TPM format (rows = TPM data, cols = cells)
  # applies the rule: Gene must be expressed in > 2 cells, at TPM > 1
  keep = apply(dat, 1, function(row){sort(row, decreasing = T)[i] > cutoff})
  return(dat[keep,])
}

Waterfall.go.GetCorrelatedGenes = function(TPM_data, log=F){#outputs a list of genes, assumes a pseudogene, pseudotime
  pt = TPM_data["Pseudotime",]
  
  GetGenePseudotimeCorrelation = function(Row){
    X = as.numeric(Row)
    if (log==T){ Row = log(Row + .0001)}
    Y = as.numeric(pt)
    C = cor(Y,X, method="spearman")
    return(C)
  }
  GeneCorrelation = apply(TPM_data, 1, GetGenePseudotimeCorrelation)
  GeneCorrelation = sort(GeneCorrelation, decreasing = T)
  return(GeneCorrelation)
}


loadGenesFromKegg = function(processID ){
  #20.) Get some lists of genes from KEGG database in Japan
  genes = keggGet(processID)[[1]]$GENE
  if (is.null(genes)) return(NULL)
  goi.clut <-unlist(strsplit(genes,split="; "))# mmu00190: oxidative phosphorylation

  #21.) convert Raw KEGG info into a list of genes invovled in OxPhos, store gene function in goi_group
  goi = goi.clut[ seq(from=2, to=length(goi.clut), by=3)]
  goi_group = goi.clut[ seq(from=3, to=length(goi.clut), by=3)]
  
  return (cbind(goi, goi_group) )
}


sin_scale2.foo <-function(m, defaultMedian = 50){
  scale_row.foo <-function(X){
    medX <-max(median(X), defaultMedian)
    maxX <-max(max(X),2* defaultMedian)
    Y <-rep(0,length(X))
    l.idx <-which(X>=medX)
    s.idx <-which(X<medX)
    
    Y[s.idx] <- 0.5*sin((X[s.idx]-medX)/medX*pi/2)+.5
    Y[l.idx] <- 0.5*sin((X[l.idx]-medX)/(max(X)-medX)*pi/2)+.5
    
    return(Y)
  }
  
  return(t(apply(m,1,scale_row.foo)))
}


plotKegg.sin = function(processID = "mmu00190", TPMdata, processTitle= NULL, type = 'l', color = NULL, ...){
  #"Plots a KEgg PROCESS, e.g. Metabolism"
  .loadGenes = loadGenesFromKegg(processID)
  if (is.null(.loadGenes)) return(NULL)
  Genes = .loadGenes[.loadGenes[,1] %in% rownames(TPMdata), 1]
  Desc = .loadGenes[.loadGenes[,1] %in% rownames(TPMdata), 2]
  
  scaleddata = ScaleTPM( TPMdata[Genes,] )

  SinScaledData = sin_scale2.foo( scaleddata ) #scaleddata #

  .summed = apply(SinScaledData, 2, function(column) mean(column))
  
  xpt = as.numeric(TPMdata["Pseudotime",]) #as.numeric(1:length(HMM_data.summed))
  ypt = as.numeric(.summed)

  
  y.loess <-loess( ypt~xpt, span = .75)
  ypt.predict <-predict(y.loess,data.frame(x=xpt),se=T)
  
  #plot( xpt, ypt, ylab = processTitle, ylim = c(min(ypt.predict$fit), max(ypt.predict$fit)))
  
  ypt_fit = ypt.predict$fit
  ypt_fit = ypt_fit - min(ypt_fit)
  ypt_fit = ypt_fit / max(ypt_fit)

  #plot(xpt,ypt_fit,lwd=1, cex = 1,main=processTitle, type = 'l', col = dc)
  if( !is.null(color)){
    plot(xpt,ypt_fit,lwd=1, cex = 1,main= processTitle, type = 'l', col = color)
  }else (plot(xpt,ypt_fit,lwd=1, cex = 1,main= processTitle, type = 'l'))
  
  print(Desc)
  ypt_fit = as.matrix(ypt_fit)
  colnames(ypt_fit) = processTitle

  return( list("x" = xpt, "y" = ypt_fit, "genes" = Genes, "ngenes" = length(Genes)))
}


MakeMultipleKeggPlots.Sin = function(processesIndex, dataset, graphTitle= NULL, color = "black", ...){
  #produces as many kegg plots as it can given a list of processes.
  #use the function SearchPOI to find the index of a given process by a search term, 
  #those indexes need to go into this func
  #note, by "idex" i am NOT referrign to mmu1000 type codes. I'm looking for a number from 0-467 (total # of Kegg processes)
  #note2 , this will output an array with z columns, z = total# of cells in the dataset
  
  processFits = sapply(processesIndex, function(idx){
    pid = as.character(Kegg.processes[idx,1])
    pname = as.character(Kegg.processes[idx,2])
    #assign("graphtitle", graphTitle, envir = .GlobalEnv)
    graphTitle = pname
    fit = rep(0, ncol(dataset))
    try( expr = (fit= plotKegg.sin(pid, dataset, color= color, processTitle = graphTitle)$y), silent = T)
    return(fit)
  })
  colnames(processFits) = Kegg.processes[processesIndex,2]
  return (t(processFits))
}

overlayPlots = function(TPM_data,processes,location='bottomleft', color = NULL, graphTitle = NULL){

  for(i in 1:length(processes)){
  
    if (!is.null(color)){
    color = color
  } else color = "black" 
  
    dv = try(as.matrix(MakeMultipleKeggPlots.Sin(Search_POI(processes[i]), TPM_data, color = color[i], processTitle = graphTitle)))
    print(dv)
    if (i == 1) {
      avd = dv
    }else{
      avd = rbind(avd, dv)
    }
    par(new=T)
  }
  final.proc <- c(processes, "avd")
  avd = apply(avd, 2, mean) #SCALE TO 1, MIN IS 0 MAX IS 1 (in ylim) [done]
  #dev.off()
  plot(as.numeric(jt["Pseudotime",]),as.numeric(avd), type = "l", lty = 2, ylim=c(0,1.0),main=graphTitle) #plots average line 
  
  ##creating the graph's legend
  legend.lty <- c() #empty array for line sample types
  for (i in 1:length(processes)){ #fills array with line sample type of processes lines
    legend.lty[i] = 1
  }
  legend.lty<- c(legend.lty, 2) #adds line type of average line 
  
  #if (!is.null(col)){ #coloring the legend correctly, average line will always be black
   # assign ("dlegendc", col, envir = .GlobalEnv)
  #} else (assign ("dlegendc", "black", envir = .GlobalEnv))
  
  legend(location, final.proc, lty=legend.lty, bty='n', cex=.75, col= c(color, "black")) 
  par(new=T)
  return(avd)
  
}

##OVERLAY PLOTS FOR GO###
overlayPlots2 = function(process_name = c(),dataset,location='bottomleft', color = NULL, graphTitle = NULL){
  for(i in 1:length(process_name)){

  dv = try(t(as.matrix(MakeGoPlots.Sin(process_name[i], dataset, exclusions = c(),min_genes = 0, color = color[i] ,plt = TRUE,give_lab = FALSE)$y)), silent = TRUE)
  par(new=T)
  print(dv)
  
  if (i == 1) {
    avd = dv
  }else{
    avd = rbind(avd, dv)
  }
  par(new=T)
  }
 
  final.proc <- c(process_name, "average")
  avd = apply(avd, 2, mean) #SCALE TO 1, MIN IS 0 MAX IS 1 (in ylim) [done]
  plot(as.numeric(dataset["Pseudotime",]),as.numeric(avd), type = "l", lty = 2, ylim=c(0,1.0),main=graphTitle) #plots average line 
  
  ##creating the graph's legend
  legend.lty <- c() #empty array for line sample types
  for (i in 1:length(process_name)){ #fills array with line sample type of processes lines
    legend.lty[i] = 1
  }
  legend.lty<- c(legend.lty, 2) #adds line type of average line 
  
  #if (!is.null(col)){ #coloring the legend correctly, average line will always be black
  # assign ("dlegendc", col, envir = .GlobalEnv)
  #} else (assign ("dlegendc", "black", envir = .GlobalEnv))
  
  legend(location, final.proc, lty=legend.lty, bty='n', cex=.75, col= c(color, "black")) 
  return(avd)
  
}

Search_POI = function(searchStrings) {
  #returns indexes of processes of interest
  indexes = c()
  for (p in searchStrings){
    idx = grep(p, Kegg.processes[,2])
    if (length(idx) != 1){
      print (p)
      print( Kegg.processes[idx,])
      break
    }
    indexes = c(idx, indexes)
  }
  return (indexes)
}



plotQuery = function(ProcessDataset, queryOutput, outputname = "plot1.pdf"){
  processes2plot = ProcessDataset[queryOutput,]
  pt = ProcessDataset["Pseudotime",]
  
  pdf(file=outputname, width=10, height=10)
  par(mfrow=c(3,3), mar=c(1.2,1.2,1.2,1.2) + .5)
  
  for (i in 1:nrow(processes2plot)){
    print(rownames(processes2plot))
    plot(t(pt), t(processes2plot[i,]), main=rownames(processes2plot)[i] )
  }
  dev.off()
}


writeKEGGProcessesTXT = function(tpm_data,out_file="KEGG_processes_activity.txt") {
  fits.all1 = MakeMultipleKeggPlots.Sin(1:476, tpm_data)
  colnames(fits.all1) = colnames(tpm_data)
  fits.all1 = rbind(fits.all1, tpm_data["Pseudotime",])
  
  write.table(fits.all1, file = out_file)
    
}


##########################
# Gene Networks
##########################

orderTop = function(v,top){
  if (length(v) < top){
    warning("top value given longer than length of ordered vector. all vals given ")
    top = length(v)
  } 
  v = v[order(v,decreasing = TRUE)]
  v = v[1:top]
  return(v)
}

getAssociations = function(path){
  m = read.csv(path)
  ids = as.vector(m[,1])
  genes = as.vector(m[,2])
  names(genes) = ids
  
  return(genes)
}

rmDupGenes = function(mat,blend_function = "mean"){
  mat = as.data.frame(mat)
  ID = rownames(mat)
  genes = mat[,1]
  dups = unique(genes[duplicated(genes)])
  mat = mat[,-1]
  avg_mat = as.data.frame(matrix(0,nrow = length(dups), ncol = ncol(mat)))
  rownames(avg_mat) = dups
  colnames(avg_mat) = colnames(mat)

  for (i in 1:length(dups)){
    g = as.character(dups[i])
    avg_vect = as.numeric(apply(mat[genes %in% g,],2,get(blend_function)))
    avg_mat[i,] = avg_vect
  }
  
  mat = mat[!genes %in% dups,]
  rownames(mat) = genes[!genes %in% dups]
  mat = rbind(mat,avg_mat)
  return(mat)
  
}

networkAnalysisBlock = function(datMat,minModuleSize = 100,softPower = 12,tpm_thresh = NULL, top_var = NULL,module_genes_path = "~/Desktop/default_modules.csv",block_save_dir = "~/Desktop/dat",write_modMat = FALSE,plotMatrix = FALSE,association_names=NULL){
  cndGenes = colnames(datMat)
  
  if (!is.null(tpm_thresh)){
    cndGenes = apply(datMat,2,max) #computes max expression
    cndGenes = names(cndGenes[cndGenes>tpm_thresh]) #removes genes below threshold
  }

  if (!is.null(top_var)) cndGenes = names(orderTop(apply(datMat[,cndGenes],2,var),top_var)) #keeps most variable genes
  
  datMat = datMat[,cndGenes]
  
  #enable multicore / multithreading 
  enableWGCNAThreads()
  allowWGCNAThreads()
  
  powers = c(c(1:10), seq(from = 1, to=20, by=2))
  sft = pickSoftThreshold(datMat, powerVector = powers, verbose = 5)
  
  sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  #this line corresponds to using an R^2 cut-off of h
  abline(h=0.90,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  
  adjacency = adjacency(datMat, power = softPower);
  TOM = TOMsimilarity(adjacency); rownames(TOM) = rownames(adjacency); colnames(TOM) = colnames(adjacency)
  dissTOM = 1-TOM 
  
  if(!file.exists(block_save_dir)) dir.create(block_save_dir)
  childName = paste0(block_save_dir,"/dat")
  bwnet = blockwiseModules(datMat, maxBlockSize = ncol(datMat),
                           power = softPower, minModuleSize = minModuleSize,
                           reassignThreshold = 0, mergeCutHeight = 0.25,
                           numericLabels = TRUE,
                           saveTOMs = TRUE,
                           saveTOMFileBase = childName,
                           verbose = 3)
  
  geneTree = bwnet$dendrograms[[1]];
  bwModuleColors = labels2colors(bwnet$colors)
  
  # open a graphics window
  #sizeGrWindow(6,6)
  par(new=FALSE)
  bwnet$dendrograms[[1]]$labels = as.vector(association_names[colnames(datMat)])
  # Plot the dendrogram and the module colors underneath for block 1
  plotDendroAndColors(bwnet$dendrograms[[1]], bwModuleColors[bwnet$blockGenes[[1]]],
                      "Module colors", main = "TOM dissimilarity & Module Assignment",
                      dendroLabels = NULL, hang = 0.01,
                      guideAll = TRUE, guideHang = 0.05,cex.dendroLabels = 0.1)
  
  
  cols = colorRampPalette(c("red", "orange","yellow","lightblue"),space = "rgb")
  if(plotMatrix){
    col_vals = cols(length(unique(dissTOM)))
    col_vals = col_vals[length(col_vals):1]
    adj_mat = adjacency
    adj_mat = adj_mat^(1/SOFT_PWR)
    diag(adj_mat) = min(adj_mat)
    TOMplot(adj_mat, geneTree, bwModuleColors, main = "Adjacency Heatmap (Selected Genes)",col = col_vals)
  } 
  
  if (is.logical(association_names)){
    if (association_names){
      association_names = colnames(datMat)
      names(association_names) = association_names
      
    }else{
      warning("Logical must bet TRUE to write file")
    }
    
  }
  if(is.null(association_names)){
    association_names = colnames(datMat)
    names(association_names) = association_names
  }
  member_counts = table(bwModuleColors)
  moduleMat = matrix("",nrow=max(member_counts),ncol=length(member_counts)*2)
  nms = rep("module membership",length(member_counts)*2)
  nms[seq(1,length(member_counts)*2,2)] = names(member_counts)
  colnames(moduleMat) = nms
  for (i in seq(1,length(member_counts)*2,2)){
    mod_col = names(member_counts)[seq(1,length(member_counts),.5)[i]]
    Ids = colnames(datMat)[which(bwModuleColors==mod_col)]
    
    if(length(Ids)>1){
      modMembership = apply(adjacency[Ids,Ids],1,sum)
    } else{
      modMembership = 1
    }
    modMembership = modMembership - min(modMembership)
    modMembership = modMembership/max(modMembership)
    srt = order(modMembership,decreasing = TRUE)
    
    genes = association_names[Ids[srt]]
    modMembership = modMembership[srt]
    
    moduleMat[1:member_counts[mod_col],i] = genes
    moduleMat[1:member_counts[mod_col],i+1] = modMembership
    
  }
  
  if(write_modMat) {write.csv(moduleMat,file = module_genes_path,row.names=FALSE) }
  
  
  
  return(list(datMat = datMat, adjacency = adjacency, TOM = TOM, dendo = geneTree, modColors = bwModuleColors, moduleMat = moduleMat))
}


corWindow = function(genes,tpm_data,mod_mat,pt,start=0,stop=100,qtl=1){
  genes = as.character(genes)
  pt = as.numeric(pt)
  cond_pt = which((pt>=start)&(pt<=stop))
  
  mai = as.numeric(mod_mat[1:length(genes),which(mNet$moduleMat[1,] %in% genes)+1])
  cond_mai = quantile(mai,qtl)
  genes = genes[which(mai>cond_mai)]
  
  cor_vals = apply(tpm_data[genes,cond_pt],1,function(v){cor(v,pt[cond_pt])})
  
  up = cor_vals[cor_vals>0]
  dn = cor_vals[cor_vals<0]
  
  return(list(up = up, dn = dn))
}

findBestSample = function(dat_mat, groups_lst = NULL,genes_lst=NULL,minModuleSize = 100,softPower = 12,tpm_thresh = NULL, top_var = NULL,association_names=NULL){
  cndGenes = colnames(dat_mat)
  
  if (!is.null(tpm_thresh)){
    cndGenes = apply(dat_mat,2,max) #computes max expression
    cndGenes = names(cndGenes[cndGenes>tpm_thresh]) #removes genes below threshold
  }
  
  if (!is.null(top_var)) cndGenes = names(orderTop(apply(dat_mat[,cndGenes],2,var),top_var)) #keeps most variable genes
  
  dat_mat = dat_mat[,cndGenes]; rownames(dat_mat) = as.vector(sapply(rownames(dat_mat),function(s){strsplit(s," ")[[1]][1]}))
  
  #enable multicore / multithreading 
  enableWGCNAThreads()
  allowWGCNAThreads()
  

  adjacency = adjacency(dat_mat, power = softPower);
  TOM = TOMsimilarity(adjacency); rownames(TOM) = rownames(adjacency); colnames(TOM) = colnames(adjacency)

  dif_mat = matrix(0,nrow = length(genes_lst),ncol = length(groups_lst)); rownames(dif_mat) = names(genes_lst); colnames(dif_mat) = names(groups_lst)
  for (i in 1:length(names(genes_lst))){
    module = names(genes_lst)[i]
    genes = genes_lst[[module]]
    mean_TOM = TOM[genes,genes]; diag(mean_TOM) = NA; mean_TOM = mean(mean_TOM,na.rm=TRUE)

    for (j in 1:length(names(groups_lst))){
      group = names(groups_lst)[j]
      samples = as.vector(sapply(groups_lst[[group]],function(s){strsplit(s," ")[[1]][1]}))
      rm_idx = which(rownames(dat_mat) %in% samples)
      if (length(rm_idx)==0) warning(paste("No samples were removed for",group))
      
      tmp_dat = dat_mat[-rm_idx,]
      
      adjacency_tmp = adjacency(tmp_dat, power = softPower);
      TOM_tmp = TOMsimilarity(adjacency_tmp); rownames(TOM_tmp) = rownames(adjacency_tmp); colnames(TOM_tmp) = colnames(adjacency_tmp)
      
      mean_TOM_tmp = TOM_tmp[genes,genes]; diag(mean_TOM_tmp) = NA; mean_TOM_tmp = mean(mean_TOM_tmp,na.rm=TRUE)
      
      red = mean_TOM_tmp - mean_TOM 
      print(group)
      print(paste0("Mean TOM:",mean_TOM,"\n","Percent Change:",100*red/mean_TOM))
      
      dif_mat[i,j] = 100*red/mean_TOM
      
    }
    
  }
  return(dif_mat)

}


findBestModule = function(mod_mat, tpm_dat, PT, take_top = 10, window = NULL,alpha=.05){
  PT = as.numeric(PT)
  
  if(!is.null(window)){
    windowIdx = which(PT>=window[1] & PT<window[2])
    PT = PT[windowIdx]
    tpm_dat = tpm_dat[,windowIdx]
    
  }
  
  ascIdx = c()
  for(i in seq(2,ncol(mod_mat),2)){
    geneCol = i-1
    modMemCol = i
    
    geneNames = as.character(mod_mat[,geneCol]); geneNames = geneNames[geneNames != ""]
    
    corVals = apply(tpm_dat[geneNames,],1,function(v){cor(v,PT,method = "spearman")})
    corVals[is.na(corVals)] = 0
    modMemVals = as.numeric(mod_mat[which(mod_mat[,geneCol] %in% geneNames),modMemCol]); names(modMemVals) = geneNames
    
    ascIdx = c(ascIdx,sum(abs(corVals[1:take_top]*modMemVals[1:take_top])))
  }
  
  names(ascIdx) = colnames(mod_mat)[seq(1,ncol(mod_mat),2)]
  
  prob_vals = dnorm(ascIdx,mean=mean(ascIdx),sd=sd(ascIdx))
  prob_vals = dnorm(ascIdx,mean=mean(ascIdx[prob_vals>alpha]),sd=sd(ascIdx[prob_vals>alpha]))
  
  d = density(ascIdx)
  centerY = (max(d$y) - min(d$y))/2
  plot(d, xlab = "Module Associadion Index", ylab = "Frequency",main = "")
  for (modColor in names(ascIdx)){
    abline(v=ascIdx[modColor],col=modColor)
    if (prob_vals[modColor] <= alpha) text(ascIdx[modColor]-((max(ascIdx)-min(ascIdx))/40),centerY,paste0(modColor," p ~ ",format(round(prob_vals[modColor],4),scientific=FALSE)),srt=90,col=modColor)
  }  

  return(prob_vals)
}

findBestModuleDifferential = function(mod_mat, tpm_dat, PT, take_top = 10, group_1 = NULL, group_2 = NULL){
  PT = as.numeric(PT)
  
  if(!is.null(group_1) & !is.null(group_2)){
    samples = c(group_1,group_2)
    tpm_dat = tpm_dat[,samples]
  }
  else{
      warning("Group counts != 2")
  }
  PT = as.numeric(PT)

  ascIdx = c()
  corVals = signatureAll(tpm_dat,group_1)
  for(i in seq(2,ncol(mod_mat),2)){
    geneCol = i-1
    modMemCol = i
    
    geneNames = as.character(mod_mat[,geneCol]); geneNames = geneNames[geneNames != ""]
    
    #run wilcoxian signed-rank test
    modMemVals = -10*log(as.numeric(mod_mat[which(mod_mat[,geneCol] %in% geneNames),modMemCol])); names(modMemVals) = geneNames
    modMemVals[!is.finite(modMemVals)] = max(modMemVals[is.finite(modMemVals)])
    ascIdx = c(ascIdx,sum(abs(-10*log(corVals[1:take_top,2])*modMemVals[1:take_top])))
  }
  
  names(ascIdx) = colnames(mod_mat)[seq(1,ncol(mod_mat),2)]
  
  prob_vals = dnorm(ascIdx,mean=mean(ascIdx),sd=sd(ascIdx))
  
  d = density(ascIdx)
  centerY = (max(d$y) - min(d$y))/2
  plot(d, xlab = "Module Associadion Index", ylab = "Frequency",main = "")
  for (modColor in names(ascIdx)){
    abline(v=ascIdx[modColor],col=modColor)
    if (ascIdx[modColor] > (mean(ascIdx+(2*sd(ascIdx))))) text(ascIdx[modColor]-((max(ascIdx)-min(ascIdx))/40),centerY,paste0(modColor," p ~ ",format(round(prob_vals[modColor],4),scientific=FALSE)),srt=90,col=modColor)
  }  
  
  ascIdx = -dnorm(ascIdx,mean=mean(ascIdx),sd=sd(ascIdx),log=TRUE)
  ascIdx = ascIdx[order(ascIdx,decreasing = TRUE)]
  
  return(ascIdx)
}

plotHubGenes = function(mod_mat,module_color,top=3,lgnd=FALSE,...){
  
  if (mode(module_color) == "character"){
    gene_i = which(colnames(mod_mat) == module_color)
    asso_i = gene_i + 1
  } 
  
  HI = mod_mat[,asso_i]
  HI = as.numeric(HI[!HI==""])
  genes = mod_mat[(1:length(HI)),gene_i]
  
  cols= rep("gray",length(HI))
  cols[1:top] = "blue"
  
  title = paste0(toupper(module_color)," Intramodular Hub Ranking")
  plot(HI,col=cols,pch = 19,main=title,xlab = "Gene Rank Order",ylab="Normalized Gene Hub Ranking",...)
  if(lgnd) legend(length(genes)/1.61803398875,1,legend = genes[1:top],title = "Hub Genes",bty = "n",pch=20,col="red")
}


pvolcanoPlot = function(log2FCs, pvals, genenames = NA, fcthresh = 2, alpha =.05, main = "Volcano Plot"){
  x = log2FCs
  y = -log10(pvals)
  x = x[-which(is.na(x))]; x[x==Inf] = max(x[x<Inf]); x[x==-Inf] = min(x[x>-Inf])
  y = y[-which(is.na(y))]
  if(!is.null(rownames(log2FCs))){   ###Row Names###
    rownames = rownames(log2FCs)
  }
  else(rownames = genenames)

  upIdx = intersect(which(x > fcthresh),which(y >-log10(alpha)))
  downIdx = intersect(which(x < -fcthresh),which(y >-log10(alpha)))
  vol_cols = rep("black", length(x)) ###Colors###
  vol_cols[upIdx] = "blue"
  vol_cols[downIdx] = "red"
  
  plot(x, y, col = vol_cols, xlab = "Log2FC", ylab = "-log10(Pval)",pch=19,cex=.25,main=main)
  abline(a=-log10(alpha), b = 0, lty = 5)
  
  upGenes = genes[upIdx]; upGenes = upGenes[order(fc[upGenes],decreasing = TRUE)]; write.csv(cbind(upGenes,fc[upGenes],p[upGenes]),"~/Desktop/upgenes.csv")
  dnGenes = genes[downIdx]; dnGenes = dnGenes[order(fc[dnGenes],decreasing = TRUE)]; write.csv(cbind(dnGenes,fc[dnGenes],p[dnGenes]),"~/Desktop/dngenes.csv")
  
  candidates = list("upGenes" = upGenes, "dnGenes" = dnGenes) 
  return(candidates)
}

pvolcanoPlot2 = function(log2FCs, pvals, genenames = NA, fcthresh = 2, alpha =.05, main = "Volcano Plot"){
  x = log2FCs
  y = -log10(pvals)
  rm_vals = unique(c(which(is.na(x)),which(is.na(y)),which(is.nan(x)),which(is.nan(y))))
  x = x[-rm_vals]; x[x==Inf] = max(x[x<Inf]); x[x==-Inf] = min(x[x>-Inf])
  y = y[-rm_vals]
  if(!is.null(rownames(log2FCs))){   ###Row Names###
    rownames = rownames(log2FCs)
  }
  else(rownames = genenames)
  
  upIdx = intersect(which(x > fcthresh),which(y >-log10(alpha)))
  downIdx = intersect(which(x < -fcthresh),which(y >-log10(alpha)))
  vol_cols = rep("black", length(x)) ###Colors###
  vol_cols[upIdx] = "blue"
  vol_cols[downIdx] = "red"
  
  plot(x, y, col = vol_cols, xlab = "Log2FC", ylab = "-log10(Pval)",pch=19,cex=.25,main=main)
  abline(a=-log10(alpha), b = 0, lty = 5)
  
  upGenes = genes[upIdx]; upGenes = upGenes[order(fc[upGenes],decreasing = TRUE)]; write.csv(cbind(upGenes,fc[upGenes],p[upGenes]),"~/Desktop/upgenes.csv")
  dnGenes = genes[downIdx]; dnGenes = dnGenes[order(fc[dnGenes],decreasing = TRUE)]; write.csv(cbind(dnGenes,fc[dnGenes],p[dnGenes]),"~/Desktop/dngenes.csv")
  
  candidates = list("upGenes" = upGenes, "dnGenes" = dnGenes) 
  return(candidates)
}

ptCor = function(exp_dat,pt,top=NULL,method="pearson",beta=1){
  pt=as.numeric(pt)
  coef = apply(exp_dat,1,function(v){cor(v^beta,pt,method = method)})
  if ("Pseudotime" %in% names(coef)) coef = coef[-which(names(coef) %in% "Pseudotime")]
  up = coef[coef>0]; up = up[order(up,decreasing = TRUE)]
  dn = coef[coef<0]
  
  return(list(up=up,dn=dn))
}

cascadingExpression = function(exp_data,beta=1,unit=1,decreasing=FALSE,on_prop = .1,
                               avg_mindist_prop = .4,write_img = FALSE,parent_dir = getwd(),
                               name="raster",cell_dim = 6.25,asn_table = NULL,...){
  
  orderGenes = function(mat){
    avg = c()
    keep = c()
    for (i in 1:nrow(mat)){
      v = mat[i,]
      onIdx = which(v==2)
      if ((length(onIdx) >= length(v) * on_prop) | (mean(minDist2D(onIdx)) <= length(v)*avg_mindist_prop)){
        keep = append(keep,i)
        avg = append(avg,mean(onIdx))
      }
    }
    return(keep[order(avg,decreasing = decreasing)])
  }
  
  pt = as.numeric(exp_data["Pseudotime",])
  hmm = t(apply(exp_data,1,function(v) state.foo(v,pt,beta=beta,unit=unit)))
  hmm = hmm[orderGenes(hmm),]
  if (!is.null(asn_table)) rownames(hmm) = asn_table[rownames(hmm)]
  
  if(write_img){
    full_path = file.path(parent_dir,paste0(name,".pdf"))
    pheatmap(hmm,cluster_row=F,cluster_col=F,color=c("#000000","#FFEE00"),legend=F,show_rownames=T,show_colnames=F,cellwidth=cell_dim,cellheight=cell_dim,...)
    print(paste0("Writing img to: ",full_path))
    pdf(file = full_path)#,width=cell_dim*ncol(hmm),height=cell_dim*nrow(hmm))
    pheatmap(hmm,cluster_row=F,cluster_col=F,color=c("#000000","#FFEE00"),legend=F,show_rownames=T,show_colnames=F,cellwidth=cell_dim,cellheight=cell_dim,...)
    dev.off()
  }
  
  return(hmm)
}



