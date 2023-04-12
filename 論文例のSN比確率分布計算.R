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

# 関数での計算
system.time(SNR_graph <- 
              Create_SNR_Confint_graph(data = list(A1,A2,A3), tag = c("水準A1","水準A2","水準A3")))


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

SNR_graph <- 
  SNR_data %>% 
  rowwise() %>% 
  mutate(
    func =
      list(Dnoncentric_Fdist_Confint),
    F2_fig =
      list(func(dof1,dof2,ncp1,ncp2,r,m,tag,SNR,
                text_F2 = "off",
                n=100000,
                rep_num=100))
  )

SNR_fig <- 
  SNR_graph %>% 
  pull(F2_fig) %>% 
  patchwork::wrap_plots(ncol = 1)

# グラフ描画
plot(SNR_fig)

# グラフ保存
graph_path="./PDF/SN比の確率分布.pdf"  
ggsave(graph_path, plot = SNR_fig, device = cairo_pdf, dpi=300, width=10, height=15)

