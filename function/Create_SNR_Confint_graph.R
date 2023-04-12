Create_SNR_Confint_graph <- function(
    data, 
    tag = c(NA),
    graph_path = "./PDF/望目特性SN比95％信頼区間.pdf",
    graph_save = FALSE,
    graph_width = 10, 
    graph_height = 15,
    text_F2 = "off",
    n = 100000,
    rep_num = 100
    ){
  # 二重非心F分布に従うと想定して望目特性SN比の確率密度分布グラフを作成する
  # data:tibble形式の対象データをリストで渡す　ex. list(A1,A2,A3)
  # tag:グラフ左上部にタグを表示する場合はテキストを記載する　
  #     データリストの要素数タグ数は一致させる必要あり　ex. tag = c("水準A1","水準A2","水準A3")
  # graph_path:グラフの保存先
  # graph_save:グラフの保存有無
  # graph_width:保存グラフの幅
  # graph_height:保存グラフの高さ
  # text_F2:二重非心F分布グラフに信頼区間の数値表示をするかどうかの選択　表示する⇒"on"
  # n:モンテカルロ法における乱数の数
  # rep_num:モンテカルロ法の繰り返し数 但し、Rのforループは遅い
  # 戻り値 SNR_fig:グラフ出力
  
  # データ作成
  if (anyNA(tag)){
    tag <- rep("", length(data))
  }
  
  SNR_data <- tibble(
    data = data,
    tag = tag
    ) %>% 
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
  
  # 二重非心F分布の計算
  SNR_graph <- 
    SNR_data %>% 
    rowwise() %>% 
    mutate(
      func =
        list(Dnoncentric_Fdist_Confint),
      F2_fig =
        list(func(dof1,dof2,ncp1,ncp2,r,m,tag,SNR,
                  text_F2 = text_F2,
                  n=n,
                  rep_num=rep_num))
    )
  # グラフの抽出
  SNR_fig <- 
    SNR_graph %>% 
    pull(F2_fig) %>% 
    patchwork::wrap_plots(ncol = 1)
  
  # グラフ保存
  if (graph_save == TRUE){
    ggsave(graph_path, plot = SNR_fig, device = cairo_pdf, dpi=300, 
           width=graph_width, height=graph_height)
  }
  return(SNR_fig)
}
