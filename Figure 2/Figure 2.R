# 坚果促生项目数据分析及可视化脚本
# 本实验涉及到的代谢组和转录组数据由百迈客公司进行测试# Figure 2A-E直接在百迈客云（BMKCloud）进行绘制（https://www.biocloud.net/）

# 设置路径
setwd("E:/macadamia/promotion/Figure 2")

# 加载必要的包
library(igraph)
library(ggplot2)
library(reshape2)
library(RColorBrewer)


# 1. 数据输入
# Rhodanobacter 相关代谢物
rhodanobacter_data <- data.frame(
  metabolite = c("neg_1310", "neg_2143", "pos_2050", "pos_2053", "pos_2372", 
                 "pos_315", "pos_339", "pos_3444", "pos_351", "pos_3647", 
                 "pos_4037", "pos_436", "pos_5279", "pos_573", "pos_5990", 
                 "pos_7645", "pos_7855"),
  coefficient = c(0.840, 0.828, 0.865, 0.865, 0.821, 0.806, 0.869, -0.855, 
                  0.855, 0.879, 0.830, 0.891, -0.818, 0.860, 0.818, -0.842, 0.857),
  pvalue = c(0.00236, 0.00307, 0.00123, 0.00123, 0.00362, 0.00486, 0.00109, 
             0.00164, 0.00164, 0.00081, 0.00294, 0.00054, 0.00381, 0.00142, 
             0.00381, 0.00222, 0.00153),
  genus = "Rhodanobacter"
)

# Ferrovibrio 相关代谢物
ferrovibrio_data <- data.frame(
  metabolite = c("neg_1243", "neg_1456", "neg_1512", "neg_979", "pos_2498", "pos_351", "pos_474"),
  coefficient = c(0.914, 0.831, 0.875, 0.834, 0.821, 0.860, 0.821),
  pvalue = c(0.00021, 0.00288, 0.00090, 0.00271, 0.00359, 0.00142, 0.00359),
  genus = "Ferrovibrio"
)

# Rhodoblastus 相关代谢物
rhodoblastus_data <- data.frame(
  metabolite = c("neg_1243", "neg_1293", "neg_1564", "neg_1595", "neg_2257", 
                 "neg_2526", "neg_4182", "pos_339", "pos_351", "pos_5279", "pos_6065"),
  coefficient = c(0.877, 0.801, 0.819, 0.806, 0.863, 0.840, -0.840, 0.801, 0.833, -0.847, 0.826),
  pvalue = c(0.00087, 0.00533, 0.00372, 0.00490, 0.00130, 0.00236, 0.00236, 0.00533, 0.00277, 0.00200, 0.00322),
  genus = "Rhodoblastus"
)

# 合并所有数据
all_edges <- rbind(rhodanobacter_data, ferrovibrio_data, rhodoblastus_data)
 
# 2. 代谢物名称映射（根据之前的映射表）
metabolite_names <- data.frame(
  ID = c("neg_1069", "neg_1089", "neg_1117", "neg_1166", "neg_1243", "neg_1293",
         "neg_1301", "neg_1310", "neg_1316", "neg_1354", "neg_1367", "neg_1428",
         "neg_1439", "neg_1442", "neg_1450", "neg_1456", "neg_1512", "neg_1533",
         "neg_1564", "neg_1595", "neg_1633", "neg_177", "neg_1873", "neg_1912",
         "neg_1951", "neg_1960", "neg_2095", "neg_2143", "neg_2257", "neg_2274",
         "neg_2286", "neg_2360", "neg_2391", "neg_2467", "neg_2526", "neg_2894",
         "neg_3073", "neg_332", "neg_3357", "neg_340", "neg_362", "neg_3637",
         "neg_3639", "neg_3654", "neg_3665", "neg_3667", "neg_3669", "neg_3671",
         "neg_3713", "neg_3725", "neg_3729", "neg_3731", "neg_3732", "neg_3758",
         "neg_4018", "neg_4107", "neg_4110", "neg_4182", "neg_448", "neg_492",
         "neg_494", "neg_513", "neg_527", "neg_573", "neg_602", "neg_604",
         "neg_667", "neg_936", "neg_979",
         "pos_1854", "pos_1897", "pos_2050", "pos_2053", "pos_2095", "pos_2185",
         "pos_2296", "pos_2342", "pos_2351", "pos_2353", "pos_2354", "pos_2356",
         "pos_2372", "pos_2498", "pos_2573", "pos_2608", "pos_2831", "pos_2895",
         "pos_2912", "pos_2933", "pos_2949", "pos_315", "pos_3169", "pos_318",
         "pos_3241", "pos_328", "pos_3288", "pos_339", "pos_3444", "pos_351",
         "pos_3647", "pos_3684", "pos_3686", "pos_379", "pos_4037", "pos_436",
         "pos_441", "pos_4508", "pos_470", "pos_472", "pos_474", "pos_4967",
         "pos_5279", "pos_573", "pos_5990", "pos_6063", "pos_6065", "pos_6287",
         "pos_6533", "pos_7103", "pos_7169", "pos_7368", "pos_7501", "pos_7502",
         "pos_7519", "pos_7520", "pos_7644", "pos_7645", "pos_7654", "pos_7655",
         "pos_7661", "pos_7743", "pos_7855"),
  Name = c("Prephenate", "19-Chloroproansamitocin", "6-Oxo-1,4,5,6-tetrahydronicotinate",
           "Ent-toddalolactone", "Aspirin", "Methylimidazole acetaldehyde",
           "Olivil 4'-O-glucoside", "L-Homocysteine", "Maltotriose", "FMN",
           "Helicid", "Dihydrodehydrodiconiferyl alcohol-4-O-glucoside",
           "6-Heptenyl Glucosinolate", "4-Methoxycinnamic acid", "fumagillin",
           "Amlodipine", "5,6,7-Trimethoxycoumarin", "2-Isopropylmaleate",
           "Salicylacyl Glucuronide", "5-(3-Hydroxypropyl)-7-methoxy-2-(3',4'-methylenedioxyphenyl)benzofuran",
           "Olivil-4'-O-glucoside", "2-Oxo-3-hydroxy-4-phosphobutanoate",
           "Catechin-5-O-glucoside", "Bruceine A", "Apiopaeonoside",
           "3-alpha(S)-Strictosidine", "3alpha-dihydrocadambine", "Nafcillin",
           "2-Methyl-1-hydroxybutyl-ThPP", "1-O-p-Cumaroylglycerol", "Arctiin",
           "Ajmalan-17-one", "L-xylo-Hexulonolactone", "5'-Deoxy-5-fluorocytidine",
           "Osmanthuside H", "Coumermic acid", "Sedoheptulose",
           "Pinocembrin-7-O-neohesperidoside", "6-Keto-prostaglandin E1", "Sucrose",
           "Tenuifoliside A", "Staphyloferrin B", "Albiflorin",
           "Eicosapentaenoic Acid ethyl ester", "Diethylcarbamazine", "L-Norvaline",
           "sanshodiol", "trans-3-Oxo-alpha-ionol", "Estradiol Cypionate",
           "Roquefortine L", "Sarpagine", "2-Hexylphosphoric Acid", "Val Glu Thr Glu",
           "1-Deoxy-11beta-hydroxypentalenate", "Alisol F", "Momordicoside F2",
           "Erythromycin C", "3,28-Dihydroxylup-20(29)-ene (Betulin)",
           "5-Hydroxy-N-formylkynurenine", "Syringetin-3-O-rutinoside",
           "4-Coumaroylshikimate", "5,6,7,4'-Tetramethoxyflavone", "epsilon-Rhodomycin T",
           "N6-(L-1,3-Dicarboxypropyl)-L-lysine", "2-Hydroxy-6-oxoocta-2,4,7-trienoate",
           "4-Pyridoxolactone", "Indoxylsulfuric acid", "Rubrofusarin-6-O-beta-D-gentiobioside",
           "Tropolone", "Ethyl maltol", "3-beta-D-Galactosyl-sn-glycerol",
           "Grevilloside Q", "N-Succinyl-LL-2,6-diaminoheptanedioate", "D-Glucosamine 1-phosphate",
           "4-Hydroxycoumarin", "Matairesinol", "hydrocinnamic acid",
           "Allysine(6-Oxo DL-Norleucine)", "Pogostone", "Carboxyibuprofen",
           "Arbutin", "ISOPEONOL", "Danshensu", "Homovanillin", "Garbanzol",
           "10,16-Dihydroxyhexadecanoic acid", "Echinone", "Deoxyshikonin",
           "6-Hydroxymelatonin", "Sapidolide A", "3.alpha.-Mannobiose", "Aucubin",
           "Plaunol E", "Petroselinic acid", "Ranunculin", "Gentamicin C2",
           "Sulpiride", "Octadecanoic acid", "N(alpha)-Acetyl-L-2,4-diaminobutyrate",
           "Gentiopicrin", "Toluate", "Triptolide", "3-[(1E,4R)-4-Hydroxycyclohex-2-en-1-ylidene]-2-oxopropanoate",
           "Gibberellin A19", "C20916", "Abscisic aldehyde", "4-Hydroxybutanoic acid",
           "Fumiquinazoline A", "Urocanic acid", "(R)-5,6-Dihydrothymine",
           "3-Methyl-1-hydroxybutyl-ThPP", "Apigenin-4'-O-glucoside", "Pyochelin",
           "(S)-DNPA", "Rimantadine", "(8S,Z)-6-((S)-3-hydroxy-2-methylpropylidene)-8-methyloctahydroindolizin-8-ol",
           "7-Dehydrodesmosterol", "Nisoldipine", "Myristoyl-EA", "Phylloquinol",
           "Isosakuranetin-7-O-glucoside (Isosakuranin)", "3-[(3aS,4S,7aS)-7a-Methyl-1,5-dioxo-octahydro-1H-inden-4-yl]propanoate",
           "Cochlearine", "Methyl-6-Paradol", "2,3-Dinor-8-iso prostaglandin F1alpha",
           "Gomisin N", "Hydroxymethylbilane", "Ciprofloxacin", "Roseoside",
           "Tobramycin", "Ketosantalic acid", "13(S)-HPOT")
)

# 合并代谢物名称
all_edges$metabolite_name <- metabolite_names$Name[match(all_edges$metabolite, metabolite_names$ID)]

# 对于未匹配的代谢物，使用原始ID
all_edges$metabolite_name[is.na(all_edges$metabolite_name)] <- all_edges$metabolite[is.na(all_edges$metabolite_name)]
 
# 3. 构建 igraph 对象

# 创建节点数据框
nodes <- data.frame(
  name = unique(c(all_edges$metabolite_name, all_edges$genus)),
  type = c(rep("metabolite", length(unique(all_edges$metabolite_name))),
           rep("genus", length(unique(all_edges$genus)))),
  stringsAsFactors = FALSE
)
nodes$type[grep("Rhodanobacter|Ferrovibrio|Rhodoblastus", nodes$name)] <- "genus"

# 创建边数据框
edges <- data.frame(
  from = all_edges$genus,
  to = all_edges$metabolite_name,
  weight = abs(all_edges$coefficient),
  correlation = all_edges$coefficient,
  direction = ifelse(all_edges$coefficient > 0, "positive", "negative"),
  pvalue = all_edges$pvalue,
  signif = ifelse(all_edges$pvalue < 0.001, "***",
                  ifelse(all_edges$pvalue < 0.01, "**", "*")),
  stringsAsFactors = FALSE
)

# 创建 igraph 图
g <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
 
# 4. 设置图形参数

# 节点颜色
V(g)$color <- ifelse(V(g)$type == "genus", 
                     ifelse(V(g)$name == "Rhodanobacter", "#E41A1C",
                            ifelse(V(g)$name == "Ferrovibrio", "#377EB8",
                                   ifelse(V(g)$name == "Rhodoblastus", "#4DAF4A", "#FFFF99"))),
                     "#FFD700")

# 节点形状：菌属为方形，代谢物为圆形
V(g)$shape <- ifelse(V(g)$type == "genus", "square", "circle")

# 节点大小
V(g)$size <- ifelse(V(g)$type == "genus", 25, 15)

# 边的颜色：正相关红色，负相关蓝色
E(g)$color <- ifelse(E(g)$direction == "positive", "#E41A1C", "#377EB8")

# 边的宽度：按相关系数的绝对值
E(g)$width <- E(g)$weight * 5

# 边的透明度
E(g)$alpha <- 0.6
 
# 5. 绘图 - 方法1：基础 igraph 绘图

setwd("E:/macadamia/promotion/network")

# 保存为 PDF
pdf("three_genera_network_basic.pdf", width = 12, height = 10)

# 布局：按菌属分组
layout <- layout_with_fr(g, niter = 1000)

plot(g,
     layout = layout,
     vertex.label = V(g)$name,
     vertex.label.cex = 0.7,
     vertex.label.dist = 1,
     vertex.label.degree = 0,
     edge.arrow.size = 0,
     main = "Core Genera-Metabolite Correlation Network\n(Rhodanobacter, Ferrovibrio, Rhodoblastus)",
     sub = paste("Positive correlations (red) | Negative correlations (blue)\n",
                 "* p < 0.05, ** p < 0.01, *** p < 0.001")
)

# 添加图例
legend("topleft",
       legend = c("Rhodanobacter", "Ferrovibrio", "Rhodoblastus", "Metabolite",
                  "Positive correlation", "Negative correlation"),
       col = c("#E41A1C", "#377EB8", "#4DAF4A", "#FFD700", "#E41A1C", "#377EB8"),
       pch = c(15, 15, 15, 19, NA, NA),
       lty = c(NA, NA, NA, NA, 1, 1),
       lwd = c(NA, NA, NA, NA, 2, 2),
       bty = "n",
       cex = 0.8)

dev.off()
 
# 6. 绘图 - 方法2：使用 ggnetwork (更美观)
 
install.packages("ggnetwork")
library(ggnetwork)
library(ggrepel)

# 转换为 ggnetwork 格式
net <- ggnetwork(g, layout = layout_with_fr, niter = 1000)

# 设置边的线型
net$linetype <- ifelse(net$direction == "positive", "solid", "dashed")

# 绘图
p <- ggplot(net, aes(x = x, y = y, xend = xend, yend = yend))
p <- p + geom_edges(aes(color = direction, alpha = 0.7, size = weight),
                    arrow = arrow(length = unit(0, "cm")))
p <- p + geom_nodes(aes(shape = type, fill = type, size = type))
p <- p + geom_text_repel(aes(label = vertex.names), size = 3.5, max.overlaps = 20)
p <- p + scale_color_manual(values = c("positive" = "#E41A1C", "negative" = "#377EB8"),
                            name = "Correlation")
p <- p + scale_fill_manual(values = c("genus" = "#4DAF4A", "metabolite" = "#FFD700"),
                           name = "Node type")
p <- p + scale_shape_manual(values = c("genus" = 22, "metabolite" = 21),
                            name = "Node type")
p <- p + scale_size_manual(values = c("genus" = 8, "metabolite" = 5),
                           guide = "none")
p <- p + scale_alpha_continuous(guide = "none", range = c(0.3, 0.8))
p <- p + theme_void()
p <- p + ggtitle("Core Genera-Metabolite Correlation Network",
                 subtitle = "Rhodanobacter | Ferrovibrio | Rhodoblastus")
p <- p + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
               plot.subtitle = element_text(hjust = 0.5, size = 10),
               legend.position = "right")

# 保存
ggsave("three_genera_network_ggnetwork.pdf", p, width = 14, height = 10, dpi = 300)
ggsave("three_genera_network_ggnetwork.png", p, width = 14, height = 10, dpi = 300)
 
# 7. 绘制分组网络图

library(circlize)

# 准备 circos 数据
circos_data <- data.frame(
  from = all_edges$genus,
  to = all_edges$metabolite_name,
  value = all_edges$coefficient,
  direction = all_edges$direction
)

# 按菌属分组
circos_data$group <- circos_data$from

# 打开 PDF
pdf("three_genera_network_circos.pdf", width = 10, height = 10)

# 初始化 circos
circos.clear()
circos.par(start.degree = 90, gap.degree = 5)

# 定义颜色
genus_colors <- c("Rhodanobacter" = "#E41A1C", 
                  "Ferrovibrio" = "#377EB8", 
                  "Rhodoblastus" = "#4DAF4A")

# 创建因子
all_nodes <- unique(c(circos_data$from, circos_data$to))
all_groups <- ifelse(all_nodes %in% names(genus_colors), all_nodes, "metabolite")
names(all_groups) <- all_nodes

# 设置颜色
group_colors <- c(genus_colors, "metabolite" = "#FFD700")

# 绘制 chord diagram
chordDiagram(circos_data[, c("from", "to")],
             col = ifelse(circos_data$direction == "positive", 
                          rgb(228/255, 26/255, 28/255, 0.6),
                          rgb(55/255, 126/255, 184/255, 0.6)),
             transparency = 0.4,
             annotationTrack = "grid",
             preAllocateTracks = list(track.height = 0.1))

# 添加标签
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + 0.1, sector.name,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),
              cex = 0.6)
}, bg.border = NA)

# 添加标题
title("Core Genera-Metabolite Correlation Network (Chord Diagram)",
      cex.main = 1.2, line = -2)

dev.off()
 
# 8. 输出统计摘要

cat("\n========== 网络统计摘要 ==========\n\n")
cat("总节点数:", vcount(g), "\n")
cat("总边数:", ecount(g), "\n")
cat("代谢物节点数:", sum(V(g)$type == "metabolite"), "\n")
cat("菌属节点数:", sum(V(g)$type == "genus"), "\n\n")

cat("各菌属连接的代谢物数量：\n")
for(genus in c("Rhodanobacter", "Ferrovibrio", "Rhodoblastus")) {
  neighbor_count <- length(neighbors(g, genus))
  cat(sprintf("  %s: %d 个代谢物\n", genus, neighbor_count))
}

cat("\n正相关边数:", sum(E(g)$direction == "positive"), "\n")
cat("负相关边数:", sum(E(g)$direction == "negative"), "\n")

# 导出边列表和节点列表
write.csv(as_data_frame(g, what = "vertices"), "network_nodes.csv", row.names = FALSE)
write.csv(as_data_frame(g, what = "edges"), "network_edges.csv", row.names = FALSE)

cat("\n文件已保存至:", getwd(), "\n")
cat("  - three_genera_network_basic.pdf\n")
cat("  - three_genera_network_ggnetwork.pdf / .png\n")
cat("  - three_genera_network_circos.pdf\n")
cat("  - network_nodes.csv\n")
cat("  - network_edges.csv\n")
