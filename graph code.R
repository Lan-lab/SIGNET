##########################################################################
### 本文档提供函数：clusterplot, featureplot, theatplot, gbarplot, highlight,
### map, levelscount, getNetElement
### clusterplot()：输入data, meta，输出按照特定方法降维后的数据矩阵
### featureplot()：输入data, feature，输出按照特定方法降维后的数据矩阵
### theatplot()：输入data, meta, metacolor，输出添加列标签的热图，颜色自定
### gbarplot()：输入data，输出可以直接提供ggplot绘图的重构数据表，注意输入
### gbarplot()的表格列必须为：regulon, group, 1k, 5k, 10k, total
### 1k,5k,10k，为0的行total栏应标注为NA，不可用斜杠代替
### highlight(): 输入字符串向量和需要高亮的字符，输出分组字符串
### map(): 输入向量与对应列表，按照映射列表替换向量内容
### levelscount(): 统计因子向量各因子水平并输出（默认输出表格格式，可选向量形式）
### getNetElement(): 输入节点映射列表，输出可igraph的网络对象，可选模式
### mode = "plane"为普通模式，mode = "crosshighlight"，将交叉元素高亮
##########################################################################

## 读取数据
data = read.csv("", header = T, row.names = 1)
meta = read.csv("", header = T, row.names = 1)
meta = as.character(meta)

## 绘图宏包
library(reshape2)
library(umap)
library(ggsci)
library(pheatmap)
library(igraph)
library(GGally)
library(ggplot2)
windowsFonts(TNM = windowsFont("Times New Roman"))
windowsFonts(HNR = windowsFont("Helvetica Neue CE 55 Roman"))

## 预处理
clusterplot = function(data, meta, method = c("umap", "tsne")){
  if(ncol(data)==length(meta)){
    data = as.data.frame(t(data))
  }
  if(method == "umap"){
    dataumap = umap::umap(data)
    dataplot = data.frame(dataumap$layout[,1:2], meta)
    colnames(dataplot) = c("UMAP_1", "UMAP_2","CellType")
    datalist = list(dataumap = dataumap, dataplot = dataplot)
  }
  if(method == "tsne"){
    datatsne = tsne::tsne(as.matrix(data), k = 2, perplexity = 50)
    dataplot = data.frame(datatsne, meta)
    colnames(dataplot) = c("tSNE_1", "tSNE_2","CellType")
    datalist = list(datatsne = datatsne, dataplot = dataplot)
  }
  return(datalist)
}

featureplot = function(data, feature, method = c("umap", "tsne")){
  r = which(rownames(data) == feature)
  f = as.vector(t(data[r,1:ncol(data)]))
  if(method == "umap"){
    dataumap = umap::umap(data)
    dataplot = data.frame(dataumap$layout[,1:2], f)
    colnames(dataplot) = c("UMAP_1", "UMAP_2","TF")
    datalist = list(dataumap = dataumap, dataplot = dataplot)
  }
  if(method == "tsne"){
    datatsne = tsne::tsne(as.matrix(data), k = 2, perplexity = 50)
    dataplot = data.frame(datatsne, f)
    colnames(dataplot) = c("tSNE_1", "tSNE_2","TF")
    datalist = list(datatsne = datatsne, dataplot = dataplot)
  }
  return(datalist)
}

# 1. 热图
theatmap = function(data, meta, metacolor, ...){
  ann_col = data.frame(CellType = meta)
  rownames(ann_col) = colnames(data)
  color_label = metacolor
  names(color_label) = levels(as.factor(meta))
  t = pheatmap(data, scale = "row", 
               clustering_distance_rows = "correlation", 
               clustering_distance_cols = "correlation",
               treeheight_row = 30, treeheight_col = 0,
               color = colorRampPalette(c("#3c5488","white","#dc0000"))(50),
               annotation_col = ann_col,
               annotation_colors = list(CellType = color_label),
               annotation_names_col = F,
               annotation_legend = T,
               show_colnames = F, border_color = NA, legend = T,
               fontsize = 9, fontsize_row = 6, ...)
  return(t)
}

## 颜色示例
# 彩虹色
color = c("#FF0000", "#FF7A00", "#FFF500", "#52FF00", "#00FFE0",
          "#00A3FF", "#1400FF", "#8F00FF", "#FF00F5", "#FF007A")
# nature配色
color = pal_npg(palette = c("nrc"))(10)
# jco配色
color = pal_jco(palette = c("default"))(10)
# 热图：颜色自定
theatmap(data, meta, metacolor = color)

# 2. 聚类
cplot = clusterplot(data, meta, method = "umap")
ggplot(cplot$dataplot, aes(UMAP_1, UMAP_2, color = CellType)) +
  geom_point(size =1, shape = 19, alpha = 0.5) + 
  scale_colour_brewer(palette="Paired") + 
  #scale_color_manual(values = color_label) +
  labs(x = "UMAP_1", y = "UMAP_2", colour ="CellType") +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_text(family = "Times New Roman", face = "bold", size = 20),
        legend.text = element_text(family = "Times New Roman",size = 16),
        #legend.position = "none",
        #panel.background = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = .5),
        panel.grid = element_blank()) 

# 3. feature map
fplot = featureplot(data, feature, method = "umap")
ggplot(fplot$dataplot, aes(UMAP_1, UMAP_2, color = TF)) +
  geom_point(size = 1, shape = 19, alpha = 0.5) + 
  scale_colour_gradient(low = "#e2e4e3", high = "#f35a79",
                        breaks = c(0.02,0.04,0.06,0.08,0.10,0.12,0.14,0.16),
                        labels = c(0.02,0.04,0.06,0.08,0.10,0.12,0.14,0.16)) + 
  labs(x = "UMAP_1", y = "UMAP_2", title = "FLOS1 regulon", colour ="AUC enrichment score") + 
  guides(colour = guide_colorbar(title.theme = element_text(family = "HNR", face = "bold", size = 16, angle = 270),
                                 title.position = "right", title.hjust = 0.5, title.vjust = 1,
                                 label.theme = element_text(family = "HNR",face = "bold", size = 12),
                                 barwidth = 0.5, barheight = 15, ticks = FALSE,
                                 frame.colour = "black", frame.linewidth = 1)) +
  theme(plot.title = element_text(family = "HNR", face = "bold",size = 20, hjust = 0.5),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.background = element_blank(),
        #panel.background = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = .5),
        panel.grid = element_blank())  

# 4. Highlight graph
data = read.csv("", header = T, row.names = 1)
dataumap = umap::umap(as.data.frame(t(data)))
# highlight: 按照term将文字向量分组
highlight = function(x, term){
  y = numeric()
  for (i in 1:length(x)) {
    if(isTRUE(x[i] %in% term)){
      y[i] = x[i]
    }
    else {
      y[i] = "Others"
    }
  }
  return(y)
}
datahl = highlight(meta, "myocyte")
dataplot = data.frame(dataumap$layout[,1:2], datahl)
colnames(dataplot) = c("UMAP_1", "UMAP_2", "CellType")
# highlight graph
ggplot(dataplot, aes(UMAP_1, UMAP_2, color = CellType)) +
  geom_point(size =1, shape = 19, alpha = 0.5) + 
  #scale_colour_brewer(palette="Paired") + 
  #scale_color_npg() +
  scale_color_manual(values = c("#3c5488", "gray80")) +
  labs(x = "UMAP_1", y = "UMAP_2", colour ="CellType") +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_text(family = "HNR", face = "bold", size = 12),
        legend.text = element_text(family = "HNR",size = 10),
        #legend.position = "none",
        #panel.background = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = .5),
        panel.grid = element_blank()) 

# 5. 柱形图
data = read.csv("", header = T)
colnames(data) = c("regulon","group","1k","5k","10k","total")

gbarplot = function(data){
  mdata = reshape2::melt(data, id = c("regulon","group","total"))
  datanum = data
  for (i in nrow(datanum)) {
    if(!is.na(datanum$total[i])){
      datanum[i,3:5] = round(datanum[i,3:5]*datanum$total[i])
    }
  }
  mdatanum = reshape2::melt(datanum, id = c("regulon","group","total"))
  return(list(data1 = mdata, data2 = mdatanum))
}

gbplot = gbarplot(data)

# 比例图
ggplot(gbplot$data1, aes(variable, value, fill = group)) + 
  geom_bar(stat = "identity", position = "dodge", width = .5,
           colour = "black", alpha = .85) +
  facet_grid(. ~ regulon) + 
  scale_y_continuous(name = "Proportion/%",
                     breaks = c(0,0.25,0.50,0.75,1.00),
                     labels = c("0","25","50","75","100")) +
  scale_x_discrete(name = "") +
  scale_fill_manual(values = c("#f39b7f","#4dbbd5"),
                    #values = c("#e64b35","#4dbbd5"),
                    #values = c("#e64b35","#3c5488"),
                    labels = c("MCDM","SCENIC"),
                    guide = guide_legend(title = "Method")) +
  theme(axis.title = element_text(family = "HNR", face = "bold", size = "10", color = "black"),
        axis.text = element_text(family = "HNR", size = "8", color = "black"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.grid = element_blank(),
        legend.title = element_text(family = "HNR", face = "bold", size = "10", color = "black"),
        legend.text = element_text(family = "HNR", size = "10", color = "black"),
        strip.text = element_text(family = "HNR", face = "bold.italic", size = "8", color = "black"), 
        strip.background = element_rect(fill = "white"))

# 绝对数图
ggplot(gbplot$data2, aes(variable, value, fill = group)) + 
  geom_bar(stat = "identity", position = "dodge", width = .5,
           colour = "black", alpha = .85) +
  facet_grid(. ~ regulon) + 
  scale_x_discrete(name = "") +
  scale_fill_manual(values = c("#f39b7f","#4dbbd5"),
                    #values = c("#e64b35","#4dbbd5"),
                    #values = c("#e64b35","#3c5488"),
                    labels = c("MCDM","SCENIC"),
                    guide = guide_legend(title = "Method")) +
  theme(axis.title = element_text(family = "HNR", face = "bold", size = "10", color = "black"),
        axis.text = element_text(family = "HNR", size = "8", color = "black"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.grid = element_blank(),
        legend.title = element_text(family = "HNR", face = "bold", size = "10", color = "black"),
        legend.text = element_text(family = "HNR", size = "10", color = "black"),
        strip.text = element_text(family = "HNR", face = "bold.italic", size = "8", color = "black"), 
        strip.background = element_rect(fill = "white"))

# 6. 网络图
## 读取数据
netdata = read.csv("data/copaired_hnscc7000_t80.csv", header = T, row.names = 1)
# 映射函数
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
# 水平统计函数
levelscount = function(x, output = "df"){
  if(!is.factor(x)){
    stop("input requires a factor")
  }
  lev = levels(x)
  l = length(lev)
  y = numeric()
  for (i in 1:l) {
    y[i] = length(which(x == lev[i]))
  }
  if (output == "v"){
    names(y) = lev
    return(y)
  }
  if (output == "df"){
    df = data.frame(levels = lev, counts = y)
    return(df)
  }
}
# 边/节点列表生成函数
getNetElement = function(data, center, color, mode = "plane", ...){
  # 普通模式
  if(mode == "plane"){
    if(length(center) + 1 > length(color)){
      stop("The length of palette is shorter than the center")
    }
    n = length(center) + 1
    colormap = data.frame(center = c(center,"NTF"), color[1:n])
    # 边列表
    row.id = data[,1] %in% center
    edgelist = data[row.id,]
    edgecolor = map(edgelist[,1], colormap)
    edgelist = data.frame(edgelist, edgecolor)
    # 节点列表
    from = as.character(levels(as.factor(edgelist[,1])))
    to = as.character(levels(as.factor(edgelist[,2])))
    dif = setdiff(to, from)
    nodeid = c(from, dif)
    nodegroup = c(from, rep("NTF", length(dif)))
    nodelist = data.frame(nodeid, nodegroup)
  }
  
  # 重节点高亮模式
  if(mode == "crosshighlight"){
    if(length(center) + 2 > length(color)){
      stop("The length of palette is shorter than the center")
    }
    n = length(center) + 2
    colormap = data.frame(center = c(center,"NTF","MNTF"), color[1:n])
    # 边列表
    row.id = data[,1] %in% center
    edgelist = data[row.id,]
    edgecolor = map(edgelist[,1], colormap)
    edgelist = data.frame(edgelist, edgecolor)
    # 节点列表
    from = as.character(levels(as.factor(edgelist[,1])))
    to = as.character(levels(as.factor(edgelist[,2])))
    dif = setdiff(to, from)
    ntf = levelscount(as.factor(edgelist[,2]))
    to2 = ntf[ntf$counts>=2,1]
    dif2 = setdiff(to2, from)
    dif = setdiff(dif, dif2)
    nodeid = c(from, dif, dif2)
    nodegroup = c(from, rep("NTF", length(dif)), rep("MNTF", length(dif2)))
    nodelist = data.frame(nodeid, nodegroup)
  }
  net = graph_from_data_frame(d=edgelist, vertices=nodelist, directed=T, ...)
  return(list(edge = edgelist, node = nodelist, net = net))
}

# 生成网络对象
start = c("FOSL1","MYC","FOS")
#netcol =  substring(pal_jco(palette = c("default"))(10), first = 1, last = 7)
netcol = c("#0073c2","#efc000","#cd534c","gray80","#003c67") # 3个tf
#netcol = c("#efc000","#cd534c","gray80","#003c67") # 2个tf
#netcol = c("#cd534c","gray80","#003c67") # 1个tf
names(netcol) = c(start, "NTF", "MNTF")
#names(netcol) = c(start, "NTF")
netsize = c(rep(20, length(start)), 8, 16)
names(netsize) = c(start, "NTF", "MNTF")
net = getNetElement(netdata, center = start, color = netcol, mode = "crosshighlight")

# 网络图
ggnet2(net$net, mode = "kamadakawai", 
       shape = 19, alpha = 0.85, 
       size = "nodegroup", size.palette = netsize, 
       color = "nodegroup", color.palette = netcol, color.legend = "Group", 
       label = T, label.size = 1.25, label.color = "white", label.alpha = 1, 
       edge.alpha = 0.75, edge.color = "edgecolor", edge.size = 0.5,
       arrow.size = 8, arrow.gap = 0.01,
       legend.size = 8, legend.position = "bottom") + 
  guides(size = F)





