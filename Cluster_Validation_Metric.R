## Metric Functions for Clustering

confusion.table <- function(mat){
  r.total <- rowSums(mat)
  c.total <- colSums(mat)
  N <- sum(r.total)
  sq.mat <- mat^2
  sq.sum <- sum(rowSums(sq.mat))
  r.sq.sum <- sum(r.total^2)
  c.sq.sum <- sum(c.total^2)
  a <- 0.5*(sq.sum-N)
  b <- 0.5*(r.sq.sum-sq.sum)
  c <- 0.5*(c.sq.sum-sq.sum)
  d <- 0.5*(sq.sum+N^2-r.sq.sum-c.sq.sum)
  c.mat <- matrix(c(a,b,c,d),nrow = 2,byrow = T)
  return(c.mat)
}

cluster.metric <- function(true,pred,method,beta=NULL){
  true <- as.factor(true)
  pred <- as.factor(pred)
  mat <- table(true,pred)
  c.mat <- confusion.table(mat)
  a <- c.mat[1,1]
  b <- c.mat[1,2]
  c <- c.mat[2,1]
  d <- c.mat[2,2]
  s <- a + b + c + d
  
  if (method == "Purity"){
    N <- length(true)
    l.pred <- length(mat[1,])
    int <- numeric()
    for (i in 1:l.pred) {
      int[i] <- max(mat[,i])
    }
    score <- sum(int)/N
  }
  else if (method == "F"){
    precision <- a/(a+c)
    recall <- a/(a+b)
    score <- (beta^2+1)*precision*recall/(beta^2*precision+recall)
  }
  else if (method == "Rand"){
    score <- (a+d)/(a+b+c+d)
  }
  else if (method == "Jaccard"){
    score <- a/(a+b+c)
  }
  else if (method == "Fowlkes-Mallows"){
    score <- a/sqrt((a+b)*(a+c))
  }
  else if (method == "Hubert-Arabie"){
    numerator <- s*(a+d)-((a+b)*(a+c)+(c+d)*(b+d))
    denominator <- s^2 - ((a+b)*(a+c)+(c+d)*(b+d))
    score <- numerator/denominator
  }
  return(score)
}
  

# external cluster validation
ExternalClusterValidation = function(true, pred, method = "ARI", beta = 1){
  N = length(true)
  
  true_lev = levels(as.factor(true))
  true_lev_num = length(true_lev)
  true = map(true, data.frame(true_lev, 1:length(true_lev)))
  pred_lev = levels(as.factor(pred))
  pred_lev_num = length(pred_lev)
  pred = map(pred, data.frame(pred_lev, 1:length(pred_lev)))
  
  ct = table(true, pred)
  r.total = rowSums(ct)
  c.total = colSums(ct)
  N = sum(r.total)
  sq.mat = ct^2
  sq.sum = sum(rowSums(sq.mat))
  r.sq.sum = sum(r.total^2)
  c.sq.sum = sum(c.total^2)
  a = 0.5*(sq.sum-N)
  b = 0.5*(r.sq.sum-sq.sum)
  c = 0.5*(c.sq.sum-sq.sum)
  d = 0.5*(sq.sum+N^2-r.sq.sum-c.sq.sum)
  s = a + b + c + d
  
  if (method == "purity" || method == "all"){
    #purity = ClusterR::external_validation(true_labels = true, clusters = pred, 
    #                                       method = "purity")
    purity = numeric()
    for (i in 1:pred_lev_num) {
      purity[i] = max(ct[,i])
    }
    purity = sum(purity)/N
  }
  else{
    purity = NA
  }
  if (method == "F" || method == "all"){
    precision = a/(a+c)
    recall = a/(a+b)
    f = (beta^2+1)*precision*recall/(beta^2*precision+recall)
  }
  else{
    f = NA
  }
  if (method == "NVD" || method == "all"){
    max_col = numeric()
    for (j in 1:ncol(ct)) {
      max_col[j] = max(ct[,j])
    }
    max_row = numeric()
    for (i in 1:nrow(ct)) {
      max_row[i] = max(ct[i,])
    }
    nvd = 1 - (sum(max_row) + sum(max_col)) / (2*N)
  }
  else{
    nvd = NA
  }
  if (method == "PSI" || method == "all"){
    ct_r = rowSums(ct)
    ct_c = colSums(ct)
    sij = matrix(0, nrow = nrow(ct), ncol = ncol(ct))
    for (i in 1:nrow(ct)) {
      for (j in 1:ncol(ct)) {
        sij[i,j] = ct[i,j]/max(ct_r[i],ct_c[j])
      }
    }
    S = sum(sij)
    
    ct_r = sort(ct_r)
    ct_c = sort(ct_c)
    E = 0
    for (k in 1:min(nrow(ct),ncol(ct))) {
      E = E + ct_r[i]*ct_c[i]/max(ct_r[i],ct_c[i])
    }
    E = E/N

    if (S >= E & max(nrow(ct),ncol(ct)) > 1){
      psi = (S - E)/(max(nrow(ct),ncol(ct)) - E)
    }
    else if (S < E){
      psi = 0
    }
    else if (nrow(ct) == 1 & ncol(ct) == 1){
      psi = 1
    }
  }
  else{
    psi = NA
  }
  
  if (method == "ARI" || method == "all"){
    #ari = ClusterR::external_validation(true_labels = true, clusters = pred, 
    #                                    method = "adjusted_rand_index")
    ari_n = s*(a+d)-((a+b)*(a+c)+(c+d)*(b+d))
    ari_d = s^2 - ((a+b)*(a+c)+(c+d)*(b+d))
    ari = ari_n/ari_d
  }
  else{
    ari = NA
  }
  if (method == "NMI" || method == "all"){
    nmi = ClusterR::external_validation(true_labels = true, clusters = pred, 
                                        method = "nmi")
  }
  else{
    nmi = NA
  }
  
  index = c(purity, f, nvd, psi, ari, nmi)
  names(index) = c("Purity", "F_Measure", "NVD", "PSI", "ARI", "NMI")
  return(index)
}

# internal cluster validation
InternalClusterValidation = function(data, cluster, method = "CDbw", 
                                     distance = "euclidean", p = 1, 
                                     neighbSize = 10, minkowski_p = 1, 
                                     threads = 1){
  map = function(x, mat){
    y = numeric()
    n = length(x)
    for (i in 1:n) {
      if (x[i] %in% mat[,1]){
        j = min(which(x[i]==mat[,1]))
        y[i] = mat[j,2]
      }
      else {
        y[i] = NA
      }
    }
    return(y)
  }
  
  cl = levels(as.factor(cluster))
  K = length(cl)
  cluster = map(cluster, data.frame(cl, 1:K))
  d = ncol(data)
  N = nrow(data)
  
  center = matrix(NA, nrow = K, ncol = d)
  for (i in 1:K) {
    center[i,] = colMeans(data[which(cluster == i),])
  }
  C = colMeans(data)
  
  dist = ClusterR::distance_matrix(data, method = distance, upper = T, diagonal = T, 
                                   minkowski_p = minkowski_p, threads = threads)
  dist_c = ClusterR::distance_matrix(center, method = distance, upper = T, diagonal = T, 
                                     minkowski_p = minkowski_p, threads = threads)
  
  data_cen = data
  for (j in 1:N) {
    data_cen[j, ] = center[cluster[j], ]
  }
  
  dist_cen = ClusterR::distance_matrix(data_cen, method = distance, upper = T, diagonal = T, 
                                       minkowski_p = minkowski_p, threads = threads)
  
  sigma = function(data){
    s = colSums((data - colMeans(data))^2)
    return(s)
  }
  centerdist = function(data, distance){
    data = rbind(data, colMeans(data))
    dist = ClusterR::distance_matrix(data, method = distance, upper = T, diagonal = T, 
                                     minkowski_p = minkowski_p, threads = threads)
    return(dist[nrow(dist), ncol(dist) - 1])
  }
  
  # icv for compactness & separatness
  
  if (method == "RMSSTD" || method == "all"){
    Wk = 0
    for (i in 1:K) {
      subdata = data[which(cluster == i), ]
      Wk = Wk + sum(sigma(subdata))
    }
    rmsstd = Wk / sqrt(N*(N-K))
    rmsstd = c(rmsstd, rmsstd, NA, NA)
  }
  else{
    rmsstd = c(NA, NA, NA, NA)
  }
  if (method == "RS" || method == "all"){
    Wk = 0
    for (i in 1:K) {
      subdata = data[which(cluster == i), ]
      Wk = Wk + sum(sigma(subdata))
    }
    W = sum(sigma(data))
    rs = 1 - Wk / W
    rs = c(rs, NA, rs, NA) 
  }
  else{
    rs = c(NA, NA, NA, NA)
  }
  if (method == "MHG" || method == "all"){
    mhg = sum(dist * dist_cen) * (2/(N*(N-1)))
    mhg = c(mhg, NA, mhg, NA) 
  }
  else{
    mhg = c(NA, NA, NA, NA)
  }
  if (method == "CH" || method == "all"){
    #icv = fpc::calinhara(data, cluster)
    Wk = 0
    for (i in 1:K) {
      subdata = rbind(data[which(cluster == i), ], center[i, ])
      dist_sub = ClusterR::distance_matrix(subdata, method = distance, upper = T, diagonal = T, 
                                           minkowski_p = minkowski_p, threads = threads)
      Wk = Wk + sum((dist_sub[nrow(dist_sub), ])^2)
    }
    Wk = Wk / (N-K)
    data_cen_c = rbind(data_cen, C)
    dist_cen_c = ClusterR::distance_matrix(data_cen_c, method = distance, upper = T, diagonal = T, 
                                           minkowski_p = minkowski_p, threads = threads)
    Bk = sum((dist_cen_c[nrow(dist_cen_c), ])^2) / (K-1)
    ch = Bk / Wk
    ch = c(ch, Wk, Bk, NA) 
  }
  else{
    ch = c(NA, NA, NA, NA)
  }
  if (method == "I" || method == "all"){
    WD_v = numeric()
    for (i in 1:K) {
      WD_vi = centerdist(data[which(cluster == i),], distance = distance)
      WD_v = c(WD_v, WD_vi)
    }
    TD_v = centerdist(data, distance = distance)
    I_c = sum(TD_v)/sum(WD_v)
    I_s = max(dist_c)
    I = (I_c * I_s / K)^p
    I = c(I, I_c, I_s, NA) 
  }
  else{
    I = c(NA, NA, NA, NA)
  }
  if (method == "D" || method == "all"){
    D = clValid::dunn(distance = dist, clusters = cluster)
    D = c(D, NA, NA, NA) 
  }
  else{
    D = c(NA, NA, NA, NA)
  }
  if (method == "S" || method == "all"){
    S = cluster::silhouette(cluster, dist = dist)
    S_sum = summary(S)
    s = c(S_sum$avg.width, NA, NA, NA) 
  }
  else{
    s = c(NA, NA, NA, NA)
  }
  if (method == "DB" || method == "all"){
    WD_s = numeric()
    for (i in 1:K) {
      WD_vi = centerdist(data[which(cluster == i),], distance = distance)
      WD_s[i] = sum(WD_vi) / length(WD_vi)
    }
    WD_s = outer(WD_s, WD_s, "+")
    WD_q = WD_s / dist_c
    WD_q_max = numeric()
    for (i in 1:K) {
      WD_q_max[i] = max(WD_q[-i,i])
    }
    db = sum(WD_q_max) / K
    db = c(db, NA, NA, NA)
  }
  else{
    db = c(NA, NA, NA, NA)
  }
  if (method == "XB" || method == "all"){
    WD_v = numeric()
    for (i in 1:K) {
      WD_vi = centerdist(data[which(cluster == i),], distance = distance)
      WD_v = c(WD_v, WD_vi)
    }
    xb_c = sum(WD_v^2)
    xb_s = numeric()
    for (i in 1:K) {
      xb_s[i] = min(dist_c[-i,i])^2
    }
    xb_s = min(xb_s)
    xb = xb_c/(N*xb_s)
    xb = c(xb, xb_c, xb_s, NA)
  }
  else{
    xb = c(NA, NA, NA, NA)
  }
  if (method == "SD" || method == "all"){
    
    dist_c_max = max(dist_c)
    dist_c_min = numeric()
    for (i in 1:K) {
      dist_c_min[i] = min(dist_c[-i,i])
    }
    dist_c_min = min(dist_c_min)
    alpha = dist_c_max / dist_c_min
    Dis = alpha * sum(1/(colSums(dist_c)))
    
    dist_max = max(dist)
    dist_min = numeric()
    for (i in 1:K) {
      dist_min[i] = min(dist[-i,i])
    }
    dist_min = min(dist_min)
    beta = dist_max / dist_min
    Dis_max = beta * sum(1/(colSums(dist)))
    
    sigma_i = 0
    for (i in 1:K) {
      subdata = data[which(cluster == i),]
      sigma_i = sigma_i + sqrt(sum((sigma(subdata))^2))/(nrow(subdata)-1)
    }
    sigma_data = sqrt(sum((sigma(data))^2))/(nrow(data)-1)
    Scat = 1/K * sigma_i / sigma_data 
    
    sd = Dis_max * Scat + Dis
    sd = c(sd, Scat, Dis, NA) 
  }
  else{
    sd = c(NA, NA, NA, NA)
  }
  if (method == "SDbw" || method == "all"){
    
    sigma_i = 0
    for (i in 1:K) {
      subdata = data[which(cluster == i),]
      sigma_i = sigma_i + sqrt(sum((sigma(subdata))^2))/(nrow(subdata)-1)
    }
    sigma_data = sqrt(sum((sigma(data))^2))/(nrow(data)-1)
    Scat = 1/K * sigma_i / sigma_data 
    
    stdev = 1/K * sqrt(sigma_i)
    sdbw_density = function(x, y){
      d = sqrt(sum((x-y)^2))
      return(ifelse(d > stdev, 0, 1))
    }
    density_ci = numeric()
    for (i in 1:K) {
      subdata = data[which(cluster == i),]
      ni = nrow(subdata)
      density_ci[i] = 0
      for (j in 1:ni) {
        density_ci[i] = density_ci[i] + sdbw_density(subdata[j,], center[i,])
      }
    }
    density_cij = matrix(NA, nrow = K, ncol = K)
    for (i in 1:(K-1)) {
      subdata_i = data[which(cluster == i),]
      for (j in (i+1):K) {
        subdata_j = data[which(cluster == j),]
        subdata_ij = rbind(subdata_i, subdata_j)
        nij = nrow(subdata_ij)
        density_cij[i,j] = 0
        cij = (center[i,] + center[j,]) / 2
        for (l in 1:nij) {
          density_cij[i,j] = density_cij[i,j] + sdbw_density(subdata_ij[l,], cij)
        }
      }
    }
    
    Dens_bw = 0
    for (i in 1:(K-1)) {
      for (j in (i+1):K) {
        Dens_bw = Dens_bw + density_cij[i,j]/max(density_ci[i], density_ci[j])
      }
    }
    Dens_bw = 2/(K*(K-1)) * Dens_bw
    
    sdbw = Scat + Dens_bw
    
    sdbw = c(sdbw, Scat, Dens_bw, NA)
  }
  else{
    sdbw = c(NA, NA, NA, NA)
  }
  if (method == "CDbw" || method == "all"){
    cdbw = fpc::cdbw(data, cluster)
    cdbw = c(cdbw$cdbw, cdbw$compactness, cdbw$sep, NA)
  }
  else{
    cdbw = c(NA, NA, NA, NA)
  }
  
  # icv for connectivity
  
  if (method == "connect" || method == "all"){
    cnnt = clValid::connectivity(distance = dist, clusters = cluster, neighbSize = neighbSize)
    cnnt = c(cnnt, NA, NA, cnnt)
  }
  else{
    cnnt = c(NA, NA, NA, NA)
  }
  
  index = data.frame(rmsstd,rs,mhg,ch,I,D,s,db,xb,sd,sdbw,cdbw,cnnt)
  colnames(index) = c("RMSSTD", "RS", "MHG", "CH", "I", "D", "S",
                      "DB", "XB", "SD", "S_dbw", "C_dbw", "Connectedness")
  rownames(index) = c("Index", "Compactness", "Separatness", "Connectedness")
  return(index)
}




