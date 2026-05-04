rm(list = ls())
cat("\014")

library(readxl)
library(ggplot2)
library(ggfortify)
library(dynlm)
library(patchwork)
library(tstools)
library(stargazer)
library(strucchange)
library(MASS)
library(car)


# PR1
# Nacteni dat a vizualizace
my_path <- file.path(getwd(), "EKON_USA_data.xlsx")
my_data <- read_xlsx(my_path)

# Prevedeni dat na casove rady
ts_data <- ts(my_data, start=c(1955,1), frequency=4)

p1 <- autoplot(ts_data[,'Log_Real_GDP'], size=1, colour = 'red', facets = FALSE)
p2 <- autoplot(ts_data[,'Unemployment'], size=1, colour = 'blue', facets = FALSE)
p1/p2



log_GDP <- ts_data[,'Log_Real_GDP']
unemp <- ts_data[, 'Unemployment']


bp_gdp <- breakpoints(log_GDP ~ time(log_GDP))
summary(bp_gdp)
plot(log_GDP)
lines(bp_gdp)


bp_unemp <- breakpoints(unemp ~ time(unemp))
summary(bp_unemp)
plot(unemp)
lines(bp_unemp)


# Tvorba trendu a dummy promenne
 D08 <- create_dummy_ts(
   end_basic = c(2025,4),
   dummy_start = c(2008,3),
   dummy_end = NULL,
   sp = FALSE,
   start_basic = c(1955, 1),
   basic_value = 0,
   dummy_value = 1,
   frequency = 4
 )

# tvorba casovy trend
Time <- ts(seq(1:length(D08)), start=c(1955,1), frequency=4)

break_point <- window(Time,start=c(2008,3),end=c(2008,3))
Time_rel <- Time - break_point[1]




# Odhad mezer vystupu
pom_data <- ts.union(log_GDP, D08, Time)

pom <- lm(log_GDP ~ Time + I(D08 * Time_rel), data = pom_data)

summary(pom)
yc <- ts(residuals(pom)*100,start=c(1955,1), frequency=4)

autoplot(yc, size=1, facets = FALSE) +
  ggtitle("Cyclical output") +
   xlab("Time") + ylab("yc")




# Odhad cyklicke nezamestnanosti
pom <- lm(unemp ~ D08)
summary(pom)
uc <- ts(residuals(pom),start=c(1955,1), frequency=4)

autoplot(uc, size=1, colour = 'black', facets = FALSE)+
   ggtitle("Cyclical unemployment")+
   xlab("Time") + ylab("uc")


# PR2
# Odhad staticke verze Okunova vztahu (rovnice 3)
Okun_1 <- lm(uc ~ yc -1)
summary(Okun_1)



# Chow test strukturalniho zlomu
cat("\n--- Chow test strukturalniho zlomu ---\n")
cat("Zlom 2008 Q3:\n")
break_point <- window(Time,start=c(2008,3),end=c(2008, 3))
sctest(uc ~ yc-1 , type = "Chow", point = break_point)



# CUSUM test
cusum_test <- efp(uc ~ yc - 1, type = "Rec-CUSUM")
plot(cusum_test, main = "Recursive CUSUM")

# Chowuv predpovedni test - pomocna funkce
chow_predictive_test <- function(y, x, break_time, freq = 4) {
  t_all <- time(y)

  # Najdi index zlomu
  idx <- which(abs(t_all - break_time) < (0.5 / freq))
  if (length(idx) == 0) stop("Zlomovy bod nebyl nalezen v casove rade.")

  n1 <- idx
  n  <- length(y)
  n2 <- n - n1
  k  <- 1  # pocet parametru (model bez interceptu)

  model_full <- lm(y ~ x - 1)
  model_sub  <- lm(y[1:n1] ~ x[1:n1] - 1)

  rss_full <- sum(residuals(model_full)^2)
  rss_sub  <- sum(residuals(model_sub)^2)

  f_stat  <- ((rss_full - rss_sub) / n2) / (rss_sub / (n1 - k))
  p_value <- pf(f_stat, df1 = n2, df2 = n1 - k, lower.tail = FALSE)

  cat("Chow Predictive Test | Zlom:", break_time, "\n")
  cat("  n1:", n1, " n2:", n2, " k:", k, "\n")
  cat("  F-stat:", round(f_stat, 4), "\n")
  cat("  p-value:", round(p_value, 4), "\n\n")

  return(invisible(list(f = f_stat, p = p_value, n1 = n1, n2 = n2)))
}

cat("\n--- Chow Predictive Test ---\n")
chow_predictive_test(uc, yc, break_time = 2008.5)



Okun_3 <- lm(uc ~ yc + I(D08*yc) -1)
summary(Okun_3)

# Okunuv koeficient po zlomu
a_D <- coefficients(Okun_3)[1]+coefficients(Okun_3)[2]
# kovariancni matice odhadu koeficientu
Sigma <- vcov(Okun_3)
std_a_D <- sqrt(Sigma[1,1]+Sigma[2,2]+2*Sigma[1,2])

linearHypothesis(Okun_3, "I(D08 * yc) = 0")
writeLines(sprintf("Okunuv koeficient po zlomu (sm. odchylka): %.4f (%.4f)\n", a_D, std_a_D))





yc_pos <- ifelse(yc > 0, yc, 0)
yc_neg <- ifelse(yc < 0, yc, 0)
model_asym <- lm(uc ~ yc_pos + yc_neg - 1)
summary(model_asym)
linearHypothesis(model_asym, "yc_pos - yc_neg = 0")



# PR 3
# vyber AIC/BIC
aic_vals <- c()
bic_vals <- c()
for (i in 1:6) {
  model <- dynlm(uc ~ L(uc, 1:i) + L(yc, 1:i))
  aic_vals[i] <- AIC(model)
  bic_vals[i] <- BIC(model)
}

best_lag_aic <- which.min(aic_vals)
best_lag_bic <- which.min(bic_vals)
cat("AIC best lag:", best_lag_aic, "\n")
cat("BIC best lag:", best_lag_bic, "\n")



best_lag_model <- dynlm(uc ~ L(uc, 1:best_lag_bic) + L(yc, 1:best_lag_bic) -1)
summary(best_lag_model)

delta <- coefficients(best_lag_model)

Sigma <- vcov(best_lag_model)

# vypocet dlouhodobeho Okunova koeficientu
a_LR <- sum(delta[3:4])/(1-sum(delta[1:2]))
writeLines(sprintf("Dlouhodoby Okunuv koeficient je %.4f \n", a_LR))




model_structural_dyn <- dynlm(uc ~ L(uc, 1:best_lag_bic) + L(yc, 1:best_lag_bic) + I(D08 * L(yc, 1:best_lag_bic)) -1)
summary(model_structural_dyn)

delta <- coefficients(model_structural_dyn)
delta2 <- delta[1:4]

Sigma <- vcov(model_structural_dyn)
Sigma2 <-Sigma[1:4, 1:4]

# vypocet dlouhodobeho Okunova koeficientu před zlomem
a_LR_3 <- sum(delta[3:4])/(1-sum(delta[1:2]))
a_LR_4 <- sum(delta[3:6])/(1-sum(delta[1:2]))
writeLines(sprintf("Dlouhodoby Okunuv koeficient před zlomem je %.4f \n", a_LR_3))
writeLines(sprintf("Dlouhodoby Okunuv koeficient po zlomu je %.4f \n", a_LR_4))




model_structural_dyn2 <- dynlm(
    uc ~ L(uc, 1:best_lag_bic) + L(yc, 1:best_lag_bic) +
        I(D08 * L(uc, 1:best_lag_bic)) + I(D08 * L(yc, 1:best_lag_bic)) -1)
summary(model_structural_dyn2)

delta <- coefficients(model_structural_dyn2)
delta2 <- delta[1:4]

Sigma <- vcov(model_structural_dyn2)
Sigma2 <- Sigma[1:4, 1:4]

# vypocet dlouhodobeho Okunova koeficientu před zlomem
a_LR_5 <- sum(delta[3:4])/(1-sum(delta[1:2]))
a_LR_6 <- (sum(delta[3:4])+sum(delta[7:8]))/(1-sum(delta[1:2])-sum(delta[5:6]))
writeLines(sprintf("Dlouhodoby Okunuv koeficient před zlomem je %.4f \n", a_LR_5))
writeLines(sprintf("Dlouhodoby Okunuv koeficient po zlomu je %.4f \n", a_LR_6))




calc_delta_se <- function(model, expr) {
  res <- deltaMethod(model, expr)
  return(res$SE)
}

# pomocna funkce pro simulacni metodu
calc_sim_se <- function(model, func_LR, n_sim = 10000) {
  set.seed(123)
  beta_hat <- coef(model)
  Sigma_hat <- vcov(model)

  # Simulace koeficientů z multivariačního normálního rozdělení
  sim_betas <- mvrnorm(n_sim, mu = beta_hat, Sigma = Sigma_hat)

  # Výpočet a_LR pro každou simulaci
  sim_a_LR <- apply(sim_betas, 1, func_LR)

  return(sd(sim_a_LR, na.rm = TRUE))
}


# zakladni dyn model (best_lag_model)
names_m1 <- names(coef(best_lag_model))
expr1 <- sprintf("(`%s` + `%s`) / (1 - `%s` - `%s`)",
                 names_m1[3], names_m1[4], names_m1[1], names_m1[2])

func_1 <- function(b) (b[3] + b[4]) / (1 - b[1] - b[2])

se_delta_1 <- deltaMethod(best_lag_model, expr1)$SE
se_sim_1   <- calc_sim_se(best_lag_model, func_1)

cat(sprintf("Model 1 (Základní): a_LR = %.4f | SE Delta = %.4f | SE Simul. = %.4f\n", a_LR, se_delta_1, se_sim_1))


# structural model 1 (jen yc - model_structural_dyn)
names_m2 <- names(coef(model_structural_dyn))

# pred zlomem
func_3 <- function(b) (b[3] + b[4]) / (1 - b[1] - b[2])
expr3 <- sprintf("(`%s` + `%s`) / (1 - `%s` - `%s`)",
                 names_m2[3], names_m2[4], names_m2[1], names_m2[2])
se_delta_3 <- deltaMethod(model_structural_dyn, expr3)$SE
se_sim_3   <- calc_sim_se(model_structural_dyn, func_3)

# po zlomu
func_4 <- function(b) (b[3] + b[4] + b[5] + b[6]) / (1 - b[1] - b[2])
expr4 <- sprintf("(`%s` + `%s` + `%s` + `%s`) / (1 - `%s` - `%s`)",
                 names_m2[3], names_m2[4], names_m2[5], names_m2[6], names_m2[1], names_m2[2])
se_delta_4 <- deltaMethod(model_structural_dyn, expr4)$SE
se_sim_4   <- calc_sim_se(model_structural_dyn, func_4)

cat(sprintf("Model 2 (Před): a_LR = %.4f | SE Delta = %.4f | SE Simul. = %.4f\n", a_LR_3, se_delta_3, se_sim_3))
cat(sprintf("Model 2 (Po):    a_LR = %.4f | SE Delta = %.4f | SE Simul. = %.4f\n", a_LR_4, se_delta_4, se_sim_4))


# structural model 2
names_m3 <- names(coef(model_structural_dyn2))

# pred zlomem
func_5 <- function(b) (b[3] + b[4]) / (1 - b[1] - b[2])
expr5 <- sprintf("(`%s` + `%s`) / (1 - `%s` - `%s`)",
                 names_m3[3], names_m3[4], names_m3[1], names_m3[2])
se_delta_5 <- deltaMethod(model_structural_dyn2, expr5)$SE
se_sim_5   <- calc_sim_se(model_structural_dyn2, func_5)

# po zlomu
func_6 <- function(b) (b[3] + b[4] + b[7] + b[8]) / (1 - (b[1] + b[2] + b[5] + b[6]))
expr6 <- sprintf("(`%s` + `%s` + `%s` + `%s`) / (1 - (`%s` + `%s` + `%s` + `%s`))",
                 names_m3[3], names_m3[4], names_m3[7], names_m3[8],
                 names_m3[1], names_m3[2], names_m3[5], names_m3[6])
se_delta_6 <- deltaMethod(model_structural_dyn2, expr6)$SE
se_sim_6   <- calc_sim_se(model_structural_dyn2, func_6)

cat(sprintf("Model 3 (Před): a_LR = %.4f | SE Delta = %.4f | SE Simul. = %.4f\n", a_LR_5, se_delta_5, se_sim_5))
cat(sprintf("Model 3 (Po):    a_LR = %.4f | SE Delta = %.4f | SE Simul. = %.4f\n", a_LR_6, se_delta_6, se_sim_6))