rchisq2 <- function(n,df1,df2,ncp1,ncp2){
  # 2重非心F分布の乱数
  # n:生成する乱数の数
  # df1,df2:自由度
  # ncp1,ncp2:非心度
  seed <- 123
  return((rchisq(n,df1,ncp1)/df1)/(rchisq(n,df2,ncp2)/df2))
}

chisq2_quantile <- function(df1,df2,ncp1,ncp2,
                            probs=c(0.00001,0.025,0.5,0.975,0.99999),
                            n=100000,
                            rep_num=100){
  # 2重非心F分布の信頼区間をモンテカルロ法にて算出する
  # df1,df2:自由度
  # ncp1,ncp2:非心度
  # probs:求める分位数
  # n:モンテカルロ法における乱数の数
  # rep_num:モンテカルロ法の繰り返し数
  # 戻り値 quan:probsで指定した分位数と乱数
  seed <- 123
  quan <- c(0,0,0,0,0)
  mean <- 0
  rcn <- c()
  pb <- progress_bar$new(total = rep_num,
                         format = "[:bar] :percent 終了まで: :eta",
                         clear = TRUE, width= 60)
  for (i in 1:rep_num){
    pb$tick()
    rc <- rchisq2(n,df1,df2,ncp1,ncp2)
    qu <- quantile(rc,probs = probs)
    quan <- quan + qu
    me <- mean(rc)
    mean <- mean + me
    rcn <- append(rcn, rc)
    Sys.sleep(1 / rep_num)
  }
  quan <- round(quan/rep_num,digits = 2)
  mean <- round(mean/rep_num,digits = 2)
  names(mean) = c("mean")

  return(c(quan,mean,rcn))
}

GetDigit <- function(num){
  # 数字の桁数を返す関数
  return(floor(log10(num)+1))
}
