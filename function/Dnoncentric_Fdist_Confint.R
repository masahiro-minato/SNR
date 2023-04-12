Dnoncentric_Fdist_Confint <- function(dof1,dof2,ncp1,ncp2,r,m,
                                      tag = "",
                                      SNR = NA,
                                      mode = "both", # "SNR","F2"
                                      x.scale_F2 = "exclusive", # "free"
                                      x.scale_SNR = "free", # "exclusive"
                                      text_F2 = "on", # "off"
                                      graph_path = "./PDF/二重非心F分布95％信頼区間.pdf",
                                      graph_save = FALSE,
                                      graph_width = 8, 
                                      graph_height = 6,
                                      n = 100000,
                                      rep_num = 1000){
  # 二重非心F分布とSN比換算の信頼区間の推定とグラフ化を行う
  
  # dof1,dof2:自由度
  # ncp1,ncp2:非心度
  # SNR:SN比の計算値 表示不要の場合はデフォルトのNAとする
  # mode:二重非心F分布とSN比のグラフ選択 "both","SNR","F2"
  # x.scale_F2:二重非心F分布のx軸範囲を0.001％～99.999％へ限定するかどうか 限定する⇒"exclusive","しない⇒free"
  # x.scale_SNR:SN比確率分布のx軸範囲を0.001％～99.999％へ限定するかどうか 限定する⇒"exclusive","しない⇒free"
  # text_F2:二重非心F分布グラフに信頼区間の数値表示をするかどうかの選択　表示する⇒"on"
  # tag:グラフ左上部にタグを表示する場合はテキストを記載する　ex. tag = "A"
  # graph_path:グラフの保存先
  # graph_save:グラフの保存有無
  # graph_width:保存グラフの幅
  # graph_height:保存グラフの高さ
  # n:モンテカルロ法における乱数の数
  # rep_num:モンテカルロ法の繰り返し数 但し、Rのforループは遅い
  # 戻り値 p_Confdensity:グラフ出力
  
  # フォント設定
  par(family="Noto Sans")
  # windowsFonts(Japan1GothicBBB = windowsFont("Japan1GothicBBB"))
  seed <- 123
  # 分位数
  probs = c(0.00001,0.025,0.5,0.975,0.99999)
  # 信頼度
  Conf <- probs[[4]]-probs[[2]]
  # 2重非心F分布の信頼区間算出
  q1_F2 <- chisq2_quantile(dof1,dof2,ncp1,ncp2,
                        probs=probs,
                        n=n,
                        rep_num=rep_num)
  print(q1_F2[c(1:6)])
  # 確率密度分布
  dens_F2 <- density(q1_F2[c(-1,-2,-3,-4,-5,-6)])
  # グラフ用にデータフレームを作成
  data_A <- tibble(X=dens_F2$x,
                   Y=dens_F2$y)
  # 信頼区間描画用
  data_B <- tibble(X2=dens_F2$x[which(dens_F2$x>q1_F2["2.5%"] & dens_F2$x<q1_F2["97.5%"])],
                   Y2=dens_F2$y[which(dens_F2$x>q1_F2["2.5%"] & dens_F2$x<q1_F2["97.5%"])])
  # 目盛り設定
  Step_F2 <- 10^(GetDigit(q1_F2["99.999%"]-q1_F2["0.001%"])-1)
  if (q1_F2["0.001%"] < Step_F2){
    Start_F2 <- 0
  }else{
    Start_F2 <- (q1_F2["0.001%"]%/%10^(GetDigit(q1_F2["0.001%"])-1))*10^(GetDigit(q1_F2["0.001%"])-1)
  }
  End_F2 <- q1_F2["99.999%"]
  if ((End_F2 / Step_F2) <= 2 ){
    Step_F2 <- Step_F2/2
  }else if ((End_F2-Start_F2)/Step_F2 < 4){
    Step_F2 <- Step_F2/2
  }
  
  # テキスト位置
  height = 0.8
  # 二重非心F分布 確率密度分布グラフ
  p_F2_Confdensity <- ggplot() +
    geom_path(data_A, mapping = aes(x=X,y=Y), color="green", linewidth=1.0) +
    # geom_path(data_A, mapping = aes(x=X,y=Y), color="green", size=1.0) +
    geom_ribbon(data_B,mapping = aes(x=X2, ymin=0, ymax=Y2), alpha=0.3, fill="lightgreen") +
    geom_vline(xintercept=c(q1_F2["mean"]), linetype="dotted") +
    theme_bw() + 
    theme(text = element_text(size = 12),
          plot.title    = element_text(color = "black", size = 12),
          plot.subtitle = element_text(color = "orange", size = 12),
          plot.tag      = element_text(color = "blue", size = 18)) +
    labs(y="density", x="F''", 
         title = paste0("二重非心Ｆ分布",as.character(Conf*100),"％信頼区間"),
         subtitle = paste0("自由度ν1=",dof1,", ν2=",dof2, " 非心度λ1=",round(ncp1,digits=2),", λ2=",round(ncp2,digits=2)),
         tag = tag)
  if (x.scale_F2 == "exclusive"){
    p_F2_Confdensity <- 
      p_F2_Confdensity + scale_x_continuous(breaks=seq(Start_F2,End_F2,Step_F2),limits = c(Start_F2,End_F2))
  }
  if (text_F2 == "on"){
    p_F2_Confdensity <- 
      p_F2_Confdensity + 
      annotate("text", x=q1_F2["2.5%"], y=-max(dens_F2$y)/40, label=as.character(round(q1_F2["2.5%"],digits=2)),
               adj="center", color="darkgreen", size = 4) +
      annotate("text", x=q1_F2["mean"], y=max(dens_F2$y)/40, label=paste0("平均 = ",as.character(round(q1_F2["mean"],digits=2))),
               adj="center", color="darkgreen", size = 4) +
      annotate("text", x=q1_F2["97.5%"], y=-max(dens_F2$y)/40, label=as.character(round(q1_F2["97.5%"],digits=2)),
               adj="center", color="darkgreen", size = 4)
  }
  
  # SN比 確率密度分布グラフ
  if (mode == "both" | mode == "SNR"){
    # SN比へ換算
    q1_SN <- 10*log10(q1_F2/(r*m))
    # 確率密度分布
    dens_SN <- density(q1_SN[c(-1,-2,-3,-4,-5,-6)])
    # グラフ用にデータフレームを作成
    data_C <- tibble(X=dens_SN$x,
                     Y=dens_SN$y)
    # 信頼区間描画用
    data_D <- tibble(X2=dens_SN$x[which(dens_SN$x>q1_SN["2.5%"] & dens_SN$x<q1_SN["97.5%"])],
                     Y2=dens_SN$y[which(dens_SN$x>q1_SN["2.5%"] & dens_SN$x<q1_SN["97.5%"])])
    # 目盛り設定
    Step_SN <- 10^(GetDigit(q1_SN["99.999%"]-q1_SN["0.001%"])-1)
    if (q1_SN["0.001%"] < Step_SN){
      Start_SN <- 0
    }else{
      Start_SN <- (q1_SN["0.001%"]%/%10^(GetDigit(q1_SN["0.001%"])-1))*10^(GetDigit(q1_SN["0.001%"])-1)
    }
    End_SN <- q1_SN["99.999%"]
    if ((End_SN / Step_SN) <= 2 ){
      Step_SN <- Step_SN/2
    }else if ((End_SN-Start_SN)/Step_SN < 4){
      Step_SN <- Step_SN/2
    }
    # SN比 確率密度分布グラフ
    p_SNR_Confdensity <- ggplot() +
      geom_path(data_C, mapping = aes(x=X,y=Y), color="red", linewidth=1.0) +
      geom_ribbon(data_D,mapping = aes(x=X2, ymin=0, ymax=Y2), alpha=0.3, fill="mistyrose") +
      geom_vline(xintercept=c(q1_SN["mean"]), linetype="dotted") +
      theme_bw() + 
      theme(text = element_text(size = 12),
            plot.title    = element_text(color = "black", size = 12),
            plot.subtitle = element_text(color = "orange", size = 12)) +
      labs(y="density", x="10log10(F''/rm)", 
           title = paste0("SN比",as.character(Conf*100),"％信頼区間"),
           subtitle = paste0("自由度ν1=",dof1,", ν2=",dof2, " 非心度λ1=",round(ncp1,digits=2),", λ2=",round(ncp2,digits=2))) +
      annotate("text", x=q1_SN["50%"], y=max(dens_SN$y)/2, label=paste0(as.character(Conf*100),"%信頼区間"),
               color="darkgreen", family = "Japan1GothicBBB", size = 4) +
      annotate("text", x=q1_SN["2.5%"], y=-max(dens_SN$y)/40, label=as.character(round(q1_SN["2.5%"],digits=2)),
               adj="center", color="darkgreen", size = 4) +
      annotate("text", x=q1_SN["mean"], y=max(dens_SN$y)/40, label=paste0("平均 = ",as.character(round(q1_SN["mean"],digits=2))),
               adj="center", color="darkgreen", size = 4) +
      annotate("text", x=q1_SN["97.5%"], y=-max(dens_SN$y)/40, label=as.character(round(q1_SN["97.5%"],digits=2)),
               adj="center", color="darkgreen", size = 4)
    if (!is.na(SNR)){
      p_SNR_Confdensity <-
        p_SNR_Confdensity + annotate("text", x=SNR, y=max(dens_SN$y)*1.03, label=paste0("SN比 = ",as.character(round(SNR,digits=2)),"db"),
                                     adj="center", color="darkgreen", size = 4)
    }
    if (x.scale_SNR == "exclusive"){
      p_SNR_Confdensity <- 
        p_SNR_Confdensity + scale_x_continuous(breaks=seq(Start_SN,End_SN,Step_SN), limits = c(Start_SN,End_SN))
    }
    if (mode == "SNR"){
      p_SNR_Confdensity <- 
        p_SNR_Confdensity + theme(plot.tag = element_text(color = "blue", size = 18)) +
        labs(tag = tag)
    }
  }
  if (mode == "both"){
    p_Confdensity <- p_F2_Confdensity + p_SNR_Confdensity
  }else if (mode == "SNR"){
    p_Confdensity <- p_SNR_Confdensity
  }else{
    p_Confdensity <- p_F2_Confdensity
  }
  # グラフ保存
  if (graph_save == TRUE){
    ggsave(graph_path, plot = p_Confdensity, device = cairo_pdf, dpi=300, 
           width=graph_width, height=graph_height)
  }
  return(p_Confdensity)
}
