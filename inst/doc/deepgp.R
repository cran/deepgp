## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center",
  fig.width = 7,
  fig.height = 5,
  echo = TRUE,
  eval = TRUE
)
set.seed(0)

## -----------------------------------------------------------------------------
library(deepgp)

## -----------------------------------------------------------------------------
higdon <- function(x) {
  i <- which(x <= 0.6)
  x[i] <- 2 * sin(pi * 0.8 * x[i] * 4) + 0.4 * cos(pi * 0.8 * x[i] * 16)
  x[-i] <- 2 * x[-i] - 1
  return(x)
}

## ---- fig.width = 6, fig.height = 4-------------------------------------------
# Training data
n <- 24
x <- seq(0, 1, length = n)
y <- higdon(x)

# Testing data
np <- 100
xp <- seq(0, 1, length = np)
yp <- higdon(xp)

plot(xp, yp, type = "l", col = 4, xlab = "X", ylab = "Y", main = "Higdon function")
points(x, y)

## ---- fig.height = 3.5, fig.width = 6-----------------------------------------
fit1 <- fit_one_layer(x, y, nmcmc = 10000, verb = FALSE)
plot(fit1)

## -----------------------------------------------------------------------------
fit1 <- trim(fit1, 5000, 2) # remove 5000 as burn-in, thin by half

## ---- fig.width = 6, fig.height = 4.5-----------------------------------------
fit1 <- predict(fit1, xp, lite = FALSE)
plot(fit1)

## ---- fig.height = 3----------------------------------------------------------
fit2 <- fit_two_layer(x, y, nmcmc = 8000, verb = FALSE)
plot(fit2)

## ---- fig.height = 4----------------------------------------------------------
fit2 <- continue(fit2, 2000, verb = FALSE)

## -----------------------------------------------------------------------------
fit2$nmcmc

## ---- fig.width = 6, fig.height = 4-------------------------------------------
plot(fit2, trace = FALSE, hidden = TRUE) 

## ---- fig.width = 6, fig.height = 4.5-----------------------------------------
fit2 <- trim(fit2, 5000, 2)
fit2 <- predict(fit2, xp, lite = FALSE)
plot(fit2)

## ---- fig.width = 6, fig.height = 4.5-----------------------------------------
fit3 <- fit_three_layer(x, y, nmcmc = 10000, verb = FALSE)
fit3 <- trim(fit3, 5000, 2)
fit3 <- predict(fit3, xp, lite = FALSE)
plot(fit3)

## -----------------------------------------------------------------------------
metrics <- data.frame("RMSE" = c(rmse(yp, fit1$mean), 
                                 rmse(yp, fit2$mean),
                                 rmse(yp, fit3$mean)),
                      "CRPS" = c(crps(yp, fit1$mean, diag(fit1$Sigma)),
                                 crps(yp, fit2$mean, diag(fit2$Sigma)),
                                 crps(yp, fit3$mean, diag(fit3$Sigma))),
                      "SCORE" = c(score(yp, fit1$mean, fit1$Sigma),
                                  score(yp, fit2$mean, fit2$Sigma),
                                  score(yp, fit3$mean, fit3$Sigma)))
rownames(metrics) <- c("One-layer GP", "Two-layer DGP", "Three-layer DGP")
metrics

## ---- eval = FALSE------------------------------------------------------------
#  fit2_vec <- fit_two_layer(x, y, nmcmc = 10000, vecchia = TRUE, m = 20)

## ---- eval = FALSE------------------------------------------------------------
#  fit2_no_g <- fit_two_layer(x, y, nmcmc = 10000, true_g = 1e-4)

## ---- eval = FALSE------------------------------------------------------------
#  fit1_sep <- fit_one_layer(x, y, nmcmc = 10000, sep = TRUE)

## -----------------------------------------------------------------------------
fit2 <- predict(fit2, xp, EI = TRUE)

## ---- fig.width = 5, fig.height = 4-------------------------------------------
plot(fit2)
par(new = TRUE)
plot(xp, fit2$EI, type = "l", lwd = 2, col = 3, axes = FALSE, xlab = "", ylab = "")
points(xp[which.max(fit2$EI)], max(fit2$EI), pch = 17, cex = 1.5, col = 3)

## ---- fig.width = 5, fig.height = 4-------------------------------------------
fit2 <- predict(fit2, xp, entropy_limit = 0)
plot(fit2)
par(new = TRUE)
plot(xp, fit2$entropy, type = "l", lwd = 2, col = 3, axes = FALSE, xlab = "", ylab = "")

## ---- echo = FALSE------------------------------------------------------------
set.seed(0)

## ---- fig.width = 5, fig.height = 4-------------------------------------------
f <- function(x) as.numeric(x > 0.5)

# Training data
x <- seq(0, 1, length = 8) 
y <- f(x)

# Testing data
xp <- seq(0, 1, length = 100)
yp <- f(xp)

plot(xp, yp, type = "l", col = 4, xlab = "X", ylab = "Y", main = "Step function")
points(x, y)

## ---- fig.width = 5, fig.height = 4-------------------------------------------
fit1 <- fit_one_layer(x, y, nmcmc = 10000, cov = "exp2", true_g = 1e-4, verb = FALSE)
fit1 <- trim(fit1, 5000, 5)
fit1 <- predict(fit1, xp)
plot(fit1)

## ---- fig.width = 5, fig.height = 4-------------------------------------------
fit2 <- fit_two_layer(x, y, nmcmc = 10000, cov = "exp2", true_g = 1e-4, verb = FALSE)
fit2 <- trim(fit2, 5000, 5)
fit2 <- predict(fit2, xp)
plot(fit2)

## ---- fig.width = 7, fig.height = 4-------------------------------------------
imse1 <- IMSE(fit1, xp)
imse2 <- IMSE(fit2, xp)
par(mfrow = c(1, 2))
plot(xp, imse1$value, type = "l", ylab = "IMSE", main = "One-layer")
points(xp[which.min(imse1$value)], min(imse1$value), pch = 17, cex = 1.5, col = 4)
plot(xp, imse2$value, type = "l", ylab = "IMSE", main = "Two-layer")
points(xp[which.min(imse2$value)], min(imse2$value), pch = 17, cex = 1.5, col = 4)

## ---- eval = FALSE------------------------------------------------------------
#  alc1 <- ALC(fit1, xp)
#  alc2 <- ALC(fit2, xp)

