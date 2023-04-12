# 論文のデータ
A1 <- tibble(
  B1 = c(31,35,31),
  B2 = c(35,40,35),
  B3 = c(35,30,35),
)
A2 <- tibble(
  B1 = c(50,40,40),
  B2 = c(60,50,55),
  B3 = c(45,45,50)
)
A3 <- tibble(
  B1 = c(40,39,34),
  B2 = c(45,46,40),
  B3 = c(42,39,36)
)

# 手動での計算
SNR_data <- tibble(
  data = list(A1,A2,A3),
  tag = c("水準A1","水準A2","水準A3")) %>% 
  rowwise() %>% 
  mutate(
    SUM = list(colSums(data)),
    r = ncol(data),
    m = nrow(data),
    ST = sum(data^2),
    Sm = sum(data)^2/(r*m),
    Snxm = sum(SUM^2)/m-Sm,
    Se = ST-Sm-Snxm,
    Ve = Se/(r*m-1-(r-1)),
    Sn = Snxm + Se,
    SNR = 10*log10(Sm/Sn),
    AVE = sum(data)/(r*m),
    LAVE = list(SUM/m),
    SUMdelta2 = sum((AVE-LAVE)^2),
    dof1 = 1,
    dof2 = r*m-1,
    ncp1 = r*m*AVE^2/Ve,
    ncp2 = m*SUMdelta2/Ve
  ) %>% 
  ungroup()
# 二重非心F分布に従うと想定したSN比のサンプリング
SNR_sample <- 
  SNR_data %>% 
  rowwise() %>% 
  mutate(
    func1 =
      list(rchisq2),
    sample = 
      list(func1(n=10000000,dof1,dof2,ncp1,ncp2)),
    sample = 
      list(10*log10(sample/(r*m)))
  )%>% 
  ungroup()


# 信頼度
sum(SNR_sample$sample[[1]] > SNR_sample$sample[[2]])/10000000 # 0.888765
sum(SNR_sample$sample[[1]] > SNR_sample$sample[[3]])/10000000 # 0.5698815
sum(SNR_sample$sample[[3]] > SNR_sample$sample[[2]])/10000000 # 0.8554172

#### 水準A2を水準A1へ変更する場合--------------------------------------
# 基準値　standard　
stn <- 2
# 比較対象　comparison
cmp <- 1

# データ作成
max_sample <- 
  max(SNR_sample$sample[[cmp]]/SNR_sample$sample[[stn]]) 
min_sample <- 
  min(SNR_sample$sample[[cmp]]/SNR_sample$sample[[stn]]) 
df1 <- tibble(
  value = SNR_sample$sample[[cmp]]/SNR_sample$sample[[stn]]
)

value = SNR_sample$sample[[cmp]]/SNR_sample$sample[[stn]]
# 確率密度分布
dens <- density(value)
probs=c(0.00001,0.025,0.05,0.5,0.975,0.99999)
quan <- quantile(value, probs = probs)
# グラフ用にデータフレームを作成
data_A <- tibble(X=dens$x,
                 Y=dens$y)
# 信頼区間描画用
data_B <- tibble(X2=dens$x[which(dens$x>quan["5%"] & dens$x<quan["99.999%"])],
                 Y2=dens$y[which(dens$x>quan["5%"] & dens$x<quan["99.999%"])])
yend <- min(dens$y[which(dens$x>quan["5%"] & dens$x<1)])
# ヒストグラム
g_A1A2 <-
  ggplot() +
  geom_histogram(df1, mapping = aes( x = value, y = after_stat(density)), 
                 position = "identity", binwidth = 0.05, alpha = 0.5 , color = "green", fill = "white") +
  geom_density(df1, mapping = aes( x = value, y = after_stat(density)), alpha = 0.5 , color = "red", linewidth = 1.0) +
  geom_ribbon(data_B, mapping = aes(x=X2, ymin=0, ymax=Y2), alpha=0.3, fill="yellow") +
  geom_vline(xintercept = 1, linetype="dotted", linewidth = 1.0) +
  geom_segment(mapping = aes(x=quan["5%"], y=0, xend=quan["5%"], yend=yend), linewidth = 1.0, color = "red") +
  theme_bw(base_family = "Japan1GothicBBB") + 
  theme(text = element_text(size = 12),
        plot.title    = element_text(color = "black", size = 12),
        plot.subtitle = element_text(color = "orange", size = 12)) +
  labs(y="density", x="向上率(A1/A2)", 
       title = "SN比向上率ヒストグラム",
       subtitle = "水準A2を水準A1へ変更した場合のSN比向上率") +
  # annotate("text", x=1, y=-0.2, label="100%", adj="center", color="darkgreen", size = 4) +
  annotate("text", x=min_sample, y=-0.1, 
           label=as.character(round(min_sample,digits = 2)), adj="center", color="darkgreen", size = 4) +
  annotate("text", x=max_sample, y=-0.1, 
           label=as.character(round(max_sample,digits = 2)), adj="center", color="darkgreen", size = 4) +
  annotate("text", x=quan["5%"], y=-0.1, label=as.character(round(quan["5%"],digits=2)),
           adj="centor", color="darkgreen", size = 4) +
  annotate("text", x=quan["50%"]*1.05, y=max(data_B$Y2)*0.2, label="95%信頼区間",
           adj="center", color="darkred", family = "Japan1GothicBBB", size = 4)
  # annotate("text", x=quan["97.5%"], y=-max(data_B$Y2)/40, label=as.character(round(quan["97.5%"],digits=2)),
  #          adj="center", color="darkgreen", size = 4)
  # scale_x_continuous(breaks=seq(min_sample,max_sample,0.5),limits = c(min_sample,max_sample))
  
plot(g_A1A2)
# グラフ保存
graph_path="./PDF/SN比向上率ヒストグラム_A1A2.pdf"  
ggsave(graph_path, plot = g_A1A2, device = cairo_pdf, dpi=300, width=8, height=6)

# 累積確率
ecdf.fun <- ecdf(SNR_sample$sample[[cmp]]/SNR_sample$sample[[stn]])
ecdf.fun(c(1.0))
ecdf.y <- sapply(seq(min_sample,max_sample,0.01),ecdf.fun)

df2 <- tibble(
  X = seq(min_sample,max_sample,0.01),
  Y = sapply(seq(min_sample,max_sample,0.01),ecdf.fun)*100
)
df3 <- tibble(
  X = seq(min_sample,1,0.01),
  Y = sapply(seq(min_sample,1,0.01),ecdf.fun)*100
)
# 累積確率グラフ
g_ecdf_A1A2 <- ggplot() +
  geom_path(df2, mapping = aes(x=X,y=Y), color="green", linewidth=1.0) +
  # geom_ribbon(df3, mapping = aes(x=X, ymin=0, ymax=Y), alpha=0.3, fill="lightgreen") +
  geom_vline(xintercept=1.0, linetype="dotted") +
  geom_hline(yintercept=ecdf.fun(c(1.0))*100, linetype="dotted") +
  theme_bw(base_family = "Japan1GothicBBB") + 
  theme(text = element_text(size = 12),
        plot.title    = element_text(color = "black", size = 12),
        plot.subtitle = element_text(color = "orange", size = 12)) +
  labs(y="累積確率(%)", x="向上率(A1/A2)", 
       title = "SN比向上率の累積確率",
       subtitle = "水準A2を水準A1へ変更した場合のSN比向上率の累積確率") +
  # scale_x_continuous(breaks=seq(0,max_sample,0.5),limits = c(0,max_sample)) +
  annotate("text", x=min_sample, y=(ecdf.fun(c(1.0))+0.03)*100, 
           label=paste0(as.character(round(ecdf.fun(c(1.0))*100,digits=2)),"%"),
           adj="left", color="darkgreen", size = 4) +
  annotate("text", x=1, y=-5, label="100%", adj="center", color="darkgreen", size = 4)
plot(g_ecdf_A1A2)
# グラフ保存
graph_path="./PDF/SN比の累積確率分布_A1A2.pdf"  
ggsave(graph_path, plot = g_ecdf_A1A2, device = cairo_pdf, dpi=300, width=8, height=6)

g_A1A2 + g_ecdf_A1A2
# 合成グラフ保存
graph_path="./PDF/SN比のヒストグラム+累積確率分布_A1A2-2.pdf"  
ggsave(graph_path, plot = g_A1A2 + g_ecdf_A1A2, device = cairo_pdf, dpi=300, width=12, height=6)

#### 水準A3を水準A1へ変更する場合--------------------------------------
# 基準値　standard　
stn <- 3
# 比較対象　comparison
cmp <- 1

# データ作成
max_sample <- 
  max(SNR_sample$sample[[cmp]]/SNR_sample$sample[[stn]]) 
min_sample <- 
  min(SNR_sample$sample[[cmp]]/SNR_sample$sample[[stn]]) 
df1 <- tibble(
  value = SNR_sample$sample[[cmp]]/SNR_sample$sample[[stn]]
)

value = SNR_sample$sample[[cmp]]/SNR_sample$sample[[stn]]
# 確率密度分布
dens <- density(value)
probs=c(0.00001,0.025,0.05,0.5,0.975,0.99999)
quan <- quantile(value, probs = probs)
# グラフ用にデータフレームを作成
data_A <- tibble(X=dens$x,
                 Y=dens$y)
# 信頼区間描画用
data_B <- tibble(X2=dens$x[which(dens$x>quan["5%"] & dens$x<quan["99.999%"])],
                 Y2=dens$y[which(dens$x>quan["5%"] & dens$x<quan["99.999%"])])
yend <- min(dens$y[which(dens$x>quan["5%"] & dens$x<1)])
# ヒストグラム
g_A1A3 <-
  ggplot() +
  geom_histogram(df1, mapping = aes( x = value, y = after_stat(density)), 
                 position = "identity", binwidth = 0.05, alpha = 0.5 , color = "green", fill = "white") +
  geom_density(df1, mapping = aes( x = value, y = after_stat(density)), alpha = 0.5 , color = "red", linewidth = 1.0) +
  geom_ribbon(data_B, mapping = aes(x=X2, ymin=0, ymax=Y2), alpha=0.3, fill="yellow") +
  geom_vline(xintercept = 1, linetype="dotted", linewidth = 1.0) +
  geom_segment(mapping = aes(x=quan["5%"], y=0, xend=quan["5%"], yend=yend), linewidth = 1.0, color = "red") +
  theme_bw(base_family = "Japan1GothicBBB") + 
  theme(text = element_text(size = 12),
        plot.title    = element_text(color = "black", size = 12),
        plot.subtitle = element_text(color = "orange", size = 12)) +
  labs(y="density", x="向上率(A1/A3)", 
       title = "SN比向上率ヒストグラム",
       subtitle = "水準A3を水準A1へ変更した場合のSN比向上率") +
  # annotate("text", x=1, y=-0.2, label="100%", adj="center", color="darkgreen", size = 4) +
  annotate("text", x=min_sample, y=-0.1, 
           label=as.character(round(min_sample,digits = 2)), adj="center", color="darkgreen", size = 4) +
  annotate("text", x=max_sample, y=-0.1, 
           label=as.character(round(max_sample,digits = 2)), adj="center", color="darkgreen", size = 4) +
  annotate("text", x=quan["5%"], y=-0.1, label=as.character(round(quan["5%"],digits=2)),
           adj="centor", color="darkgreen", size = 4) +
  annotate("text", x=quan["50%"]*1.05, y=max(data_B$Y2)*0.2, label="95%信頼区間",
           adj="center", color="darkred", family = "Japan1GothicBBB", size = 4)
# annotate("text", x=quan["97.5%"], y=-max(data_B$Y2)/40, label=as.character(round(quan["97.5%"],digits=2)),
#          adj="center", color="darkgreen", size = 4)
# scale_x_continuous(breaks=seq(min_sample,max_sample,0.5),limits = c(min_sample,max_sample))

plot(g_A1A3)
# グラフ保存
graph_path="./PDF/SN比向上率ヒストグラム_A1A3.pdf"  
ggsave(graph_path, plot = g_A1A3, device = cairo_pdf, dpi=300, width=8, height=6)

# 累積確率
ecdf.fun <- ecdf(SNR_sample$sample[[cmp]]/SNR_sample$sample[[stn]])
ecdf.fun(c(1.0))
ecdf.y <- sapply(seq(min_sample,max_sample,0.01),ecdf.fun)

df2 <- tibble(
  X = seq(min_sample,max_sample,0.01),
  Y = sapply(seq(min_sample,max_sample,0.01),ecdf.fun)*100
)
df3 <- tibble(
  X = seq(min_sample,1,0.01),
  Y = sapply(seq(min_sample,1,0.01),ecdf.fun)*100
)
# 累積確率グラフ
g_ecdf_A1A3 <- ggplot() +
  geom_path(df2, mapping = aes(x=X,y=Y), color="green", linewidth=1.0) +
  # geom_ribbon(df3, mapping = aes(x=X, ymin=0, ymax=Y), alpha=0.3, fill="lightgreen") +
  geom_vline(xintercept=1.0, linetype="dotted") +
  geom_hline(yintercept=ecdf.fun(c(1.0))*100, linetype="dotted") +
  theme_bw(base_family = "Japan1GothicBBB") + 
  theme(text = element_text(size = 12),
        plot.title    = element_text(color = "black", size = 12),
        plot.subtitle = element_text(color = "orange", size = 12)) +
  labs(y="累積確率(%)", x="向上率(A1/A3)", 
       title = "SN比向上率の累積確率",
       subtitle = "水準A3を水準A1へ変更した場合のSN比向上率の累積確率") +
  # scale_x_continuous(breaks=seq(0,max_sample,0.5),limits = c(0,max_sample)) +
  annotate("text", x=min_sample, y=(ecdf.fun(c(1.0))+0.03)*100, 
           label=paste0(as.character(round(ecdf.fun(c(1.0))*100,digits=2)),"%"),
           adj="left", color="darkgreen", size = 4) +
  annotate("text", x=1, y=-5, label="100%", adj="center", color="darkgreen", size = 4)
plot(g_ecdf_A1A3)
# グラフ保存
graph_path="./PDF/SN比の累積確率分布_A1A3.pdf"  
ggsave(graph_path, plot = g_ecdf_A1A3, device = cairo_pdf, dpi=300, width=8, height=6)

g_A1A3 + g_ecdf_A1A3
# 合成グラフ保存
graph_path="./PDF/SN比のヒストグラム+累積確率分布_A1A3-2.pdf"  
ggsave(graph_path, plot = g_A1A3 + g_ecdf_A1A3, device = cairo_pdf, dpi=300, width=12, height=6)
