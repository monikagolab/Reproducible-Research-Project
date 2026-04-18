library(quantmod)
library(vars)     
library(urca)
library(xts)

set.seed(123)

# 1. Dane ceny -> log-zwroty (BTC/ETH/SOL)
tickers   <- c("BTC-USD","ETH-USD","SOL-USD")
startDate <- as.Date("2017-01-01")
endDate <- as.Date("2025-09-06")  # data użyta w raporcie


get_px <- function(sym) {
  # 1 kolumna: 
  suppressWarnings(Ad(getSymbols(sym, src = "yahoo",
                                 from = startDate, to = endDate,
                                 auto.assign = FALSE)))
}

px <- Reduce(merge, lapply(tickers, get_px))
colnames(px) <- c("BTC","ETH","SOL")

ret3 <- na.omit(diff(log(px)))   # wspólny przekrój po merge + na.omit
cat("Zakres 3-coin:", as.character(start(ret3)), "→", as.character(end(ret3)), "\n")

# Rysunek 1. Log-zwroty BTC/ETH/SOL (wspólny przekrój)
rysunek1 <- function(ret3){
  op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
  par(mar = c(5,5,2.2,2), cex.axis = 0.9, cex.lab = 1.0,
      mgp = c(2.2, 0.6, 0), lend = 1, xaxs = "i", yaxs = "i")
  
  rng  <- range(ret3, na.rm = TRUE)
  cols <- c(BTC = adjustcolor("black",      alpha.f = 0.60),
            ETH = adjustcolor("red3",       alpha.f = 0.55),
            SOL = adjustcolor("forestgreen",alpha.f = 0.55))
  
  plot(index(ret3), coredata(ret3$BTC), type = "h", lwd = 0.6, col = cols["BTC"],
       xlab = "", ylab = "log-zwrot",
       main = "Log-zwroty BTC/ETH/SOL (wspólny przekrój)", ylim = rng)
  lines(index(ret3), coredata(ret3$ETH), type = "h", lwd = 0.6, col = cols["ETH"])
  lines(index(ret3), coredata(ret3$SOL), type = "h", lwd = 0.6, col = cols["SOL"])
  
  abline(h = 0, col = "grey60", lwd = 1.2)
  legend("topright", bty = "n", inset = 0.01, lwd = 1.5,
         col = cols, legend = c("BTC","ETH","SOL"))
}

rysunek1(ret3)

# 2. Stacjonarność (ADF, PP i KPSS na zwrotach)
print(summary(ur.df(ret3$BTC, type="drift", selectlags="AIC")))
print(summary(ur.pp(ret3$BTC, type="Z-tau", model="constant", lags="short")))
print(summary(ur.kpss(ret3$BTC, type="mu", lags="short")))

print(summary(ur.df(ret3$ETH, type="drift", selectlags="AIC")))
print(summary(ur.pp(ret3$ETH, type="Z-tau", model="constant", lags="short")))
print(summary(ur.kpss(ret3$ETH, type="mu", lags="short")))

print(summary(ur.df(ret3$SOL, type="drift", selectlags="AIC")))
print(summary(ur.pp(ret3$SOL, type="Z-tau", model="constant", lags="short")))
print(summary(ur.kpss(ret3$SOL, type="mu", lags="short")))


# 3. Wybór p (BIC->AIC) i estymacja VAR na pełnej próbie

# 3.1 Wybór rzędu p (BIC -> AIC)
lag.max <- min(10L, floor(12 * (nrow(ret3)/100)^(1/4)))  # opcjonalna heurystyka
sel <- VARselect(ret3, lag.max = lag.max, type = "const")
p <- sel$selection[["SC(n)"]]
if (is.na(p) || p < 1) p <- sel$selection[["AIC(n)"]]
p <- ifelse(is.na(p) || p < 1, 1L, as.integer(p))
cat("Wybrany rząd p (BIC->AIC):", p, "(lag.max =", lag.max, ")\n")
print(round(sel$criteria, 3))  # AIC, HQ, SC, FPE dla przejrzystości

# 3.2 Estymacja VAR(p) ze stałą
var_model <- VAR(ret3, p = p, type = "const")
sm <- summary(var_model)
print(sm)

N_total <- nrow(ret3); N_eff <- N_total - p
cat(sprintf("Obserwacje: total=%d, użyte=%d (utracone %d przez opóźnienia), p=%d, deterministycznie: const\n",
            N_total, N_eff, p, p))

# 3.3 Stabilność
rts <- roots(var_model)
cat(sprintf("Max |root| = %.5f  %s\n",
            max(abs(rts)), ifelse(max(abs(rts)) < 1, "→ stabilny", "-> niestabilny")))

## 4. Diagnostyka

# 4.1 Autokorelacja reszt (Portmanteau)
print(serial.test(var_model, lags.pt = 10, type = "PT.asymptotic"))

# 4.2 Heteroskedastyczność (ARCH)
print(arch.test(var_model, lags.multi = 5))

# 4.3 Normalność (JB: univariate + multivariate)
print(normality.test(var_model, multivariate.only = FALSE))

# 4.4 Stabilność (CUSUM)
stab <- stability(var_model)
plot(stab)

# 4.5 Macierz korelacji reszt
R <- round(cor(residuals(var_model)), 4)
print(R)

# 5. IRF (Cholesky & uogólnione) + FEVD + Granger

H     <- 3
vars3 <- colnames(ret3)

# 5.1 IRF (Cholesky)

# BTC – większy lewy margines
op <- par(no.readonly = TRUE)
par(ask = FALSE, mar = c(5, 8, 3, 2), cex.axis = 0.8, cex.lab = 1.1, mgp = c(3, 1.6, 0))
plot(irf(var_model, impulse = "BTC", response = vars3,
         n.ahead = 3, ortho = TRUE, boot = TRUE, runs = 500, ci = 0.95))
par(op)  # przywróć poprzednie ustawienia graficzne

# ETH – trochę inne marginesy/cexy
op <- par(no.readonly = TRUE)
par(ask = FALSE, mar = c(5, 9, 3, 2), cex.axis = 0.75, cex.lab = 1.1, mgp = c(3, 2, 0))
plot(irf(var_model, impulse = "ETH", response = vars3,
         n.ahead = H, ortho = TRUE, boot = TRUE, runs = 500, ci = 0.95))

par(op)  # przywróć poprzednie ustawienia graficzne

# SOL – domyślne
plot(irf(var_model, impulse = "SOL", response = vars3,
         n.ahead = H, ortho = TRUE, boot = TRUE, runs = 500, ci = 0.95))

# 5.2 IRF uogólnione (Pesaran–Shin) – domyślne
plot(irf(var_model, impulse = "BTC", response = vars3,
         n.ahead = H, ortho = FALSE, boot = TRUE, runs = 500, ci = 0.95))

plot(irf(var_model, impulse = "ETH", response = vars3,
         n.ahead = H, ortho = FALSE, boot = TRUE, runs = 500, ci = 0.95))

plot(irf(var_model, impulse = "SOL", response = vars3,
         n.ahead = H, ortho = FALSE, boot = TRUE, runs = 500, ci = 0.95))

# 5.3 FEVD
fevd_res <- fevd(var_model, n.ahead = H)
plot(fevd_res)

# 5.4 Granger 
g_btc <- causality(var_model, cause = "BTC")$Granger
g_eth <- causality(var_model, cause = "ETH")$Granger
g_sol <- causality(var_model, cause = "SOL")$Granger
gr_tb <- data.frame(
  cause   = c("BTC","ETH","SOL"),
  F       = c(g_btc$statistic,   g_eth$statistic,   g_sol$statistic),
  df1     = c(g_btc$parameter[1],g_eth$parameter[1],g_sol$parameter[1]),
  df2     = c(g_btc$parameter[2],g_eth$parameter[2],g_sol$parameter[2]),
  p_value = c(g_btc$p.value,     g_eth$p.value,     g_sol$p.value)
)
print(gr_tb)

# 6. Walidacja OOS (expanding 1-step)
stopifnot(exists("ret3"), exists("p"))
set.seed(123)

# 6.1 Parametry
n <- nrow(ret3)
split  <- max(p + 50, floor(0.8 * n))      # start testu: 80% próby lub min. p+50
series <- c("BTC","ETH","SOL")
idx_test <- (split+1):n

# 6.2 Pomocnicze funkcje
fit_var_and_forecast_1 <- function(train_mat, p) {
  fit <- VAR(train_mat, p = p, type = "const")
  fc  <- predict(fit, n.ahead = 1)$fcst
  c(BTC = fc$BTC[1,1], ETH = fc$ETH[1,1], SOL = fc$SOL[1,1])
}
fit_ar1_and_forecast_1 <- function(y) {
  if (length(y) < 5) return(NA_real_)
  y1 <- y[-length(y)]; y2 <- y[-1]
  fit <- lm(y2 ~ y1)
  as.numeric(coef(fit)[1] + coef(fit)[2] * tail(y, 1))
}
MAE  <- function(e) mean(abs(e), na.rm = TRUE)
RMSE <- function(e) sqrt(mean(e^2, na.rm = TRUE))

# 6.3 Rolling OOS: prognozy VAR / AR(1) / Zero
oos <- data.frame(
  date    = as.Date(index(ret3)[idx_test]),
  BTC     = as.numeric(ret3[idx_test, "BTC"]),
  ETH     = as.numeric(ret3[idx_test, "ETH"]),
  SOL     = as.numeric(ret3[idx_test, "SOL"]),
  VAR_BTC = NA_real_, VAR_ETH = NA_real_, VAR_SOL = NA_real_,
  AR1_BTC = NA_real_, AR1_ETH = NA_real_, AR1_SOL = NA_real_,
  ZERO_BTC = 0, ZERO_ETH = 0, ZERO_SOL = 0
)

row <- 1L
for (t in split:(n-1)) {
  train <- ret3[1:t, ]
  # VAR
  oos[row, c("VAR_BTC","VAR_ETH","VAR_SOL")] <- fit_var_and_forecast_1(train, p)
  # AR(1)
  for (s in series) oos[row, paste0("AR1_", s)] <- fit_ar1_and_forecast_1(as.numeric(train[, s]))
  if (row %% 500 == 0) cat("OOS step:", row, "z", length(idx_test), "\n")
  row <- row + 1L
}

# 6.4 Metryki: per seria + średnio
models <- c("VAR","AR1","ZERO")
met_list <- list()
for (m in models) for (s in series) {
  y <- oos[[s]]; yhat <- oos[[paste0(m,"_",s)]]
  met_list[[length(met_list)+1]] <- data.frame(model=m, series=s,
                                               mae=MAE(y-yhat), rmse=RMSE(y-yhat))
}
metrics <- do.call(rbind, met_list)
avg <- aggregate(cbind(mae, rmse) ~ model, data = metrics, FUN = mean)

# 6.5 Wydruk
pretty <- function(df) within(df, { mae = round(as.numeric(mae),6); rmse = round(as.numeric(rmse),6) })

cat("Rolling OOS (expanding), 1-step ahead\n")
cat(sprintf("Train/Test split: %d / %d obserwacji\n", split, n - split))
cat(sprintf("VAR lag p = %d\n\n", p))

cat("Metryki per seria (MAE / RMSE):\n")
print(pretty(metrics)[order(metrics$series, metrics$model), ])

cat("\nŚrednio po seriach (MAE / RMSE):\n")
print(pretty(avg))

best_rmse <- avg$model[which.min(avg$rmse)]
best_mae  <- avg$model[which.min(avg$mae)]
cat(sprintf("\nWinner (RMSE, avg): %s = %.6f\n", best_rmse, min(avg$rmse)))
cat(sprintf("Winner (MAE,  avg): %s = %.6f\n", best_mae,  min(avg$mae)))
