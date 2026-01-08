## ----include = FALSE----------------------------------------------------------
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
booth <- function(x) {
  y <- rep(NA, length(x))
  for (i in 1:length(x)) {
    if (x[i] <= 0.58) {
      y[i] <- sin(pi * x[i] * 6) + cos(pi * x[i] * 12)
    } else y[i] <- 5 * x[i] - 4.9
  }
  return(y)
}

## ----fig.width = 6, fig.height = 4.5------------------------------------------
# Training data
n <- 20
x_booth <- seq(0, 1, length = n)
y_booth <- booth(x_booth)

# Testing data
np <- 100
xp_booth <- seq(0, 1, length = np)
yp_booth <- booth(xp_booth)

plot(xp_booth, yp_booth, type = "l", col = 4, xlab = "X", ylab = "Y", 
     main = "Booth function")
points(x_booth, y_booth)

## ----fig.height = 3.5, fig.width = 6------------------------------------------
gp_booth <- fit_one_layer(x_booth, y_booth, nmcmc = 10000, true_g = 1e-6, 
                          verb = FALSE)
plot(gp_booth)

## -----------------------------------------------------------------------------
gp_booth <- trim(gp_booth, 5000, 2) # remove 5000 as burn-in, thin by half

## ----fig.width = 6, fig.height = 4.5------------------------------------------
gp_booth <- predict(gp_booth, xp_booth, lite = FALSE)
plot(gp_booth)

## ----fig.height = 3-----------------------------------------------------------
dgp_booth <- fit_two_layer(x_booth, y_booth, nmcmc = 8000, true_g = 1e-6, 
                           verb = FALSE)
plot(dgp_booth)

## ----fig.height = 4-----------------------------------------------------------
dgp_booth <- continue(dgp_booth, 2000, verb = FALSE)

## -----------------------------------------------------------------------------
dgp_booth$nmcmc

## ----fig.width = 6, fig.height = 4--------------------------------------------
plot(dgp_booth, trace = FALSE, hidden = TRUE) 

## ----fig.width = 6, fig.height = 4.5------------------------------------------
dgp_booth <- trim(dgp_booth, 5000, 2)
dgp_booth <- predict(dgp_booth, xp_booth, lite = FALSE)
plot(dgp_booth)

## ----fig.width = 6, fig.height = 4.5------------------------------------------
dgp3_booth <- fit_three_layer(x_booth, y_booth, nmcmc = 10000, true_g = 1e-6, 
                              verb = FALSE)
dgp3_booth <- trim(dgp3_booth, 5000, 2)
dgp3_booth <- predict(dgp3_booth, xp_booth, lite = FALSE)
plot(dgp3_booth)

## -----------------------------------------------------------------------------
metrics <- data.frame("RMSE" = c(rmse(yp_booth, gp_booth$mean), 
                                 rmse(yp_booth, dgp_booth$mean),
                                 rmse(yp_booth, dgp3_booth$mean)),
                      "CRPS" = c(crps(yp_booth, gp_booth$mean, diag(gp_booth$Sigma)),
                                 crps(yp_booth, dgp_booth$mean, diag(dgp_booth$Sigma)),
                                 crps(yp_booth, dgp3_booth$mean, diag(dgp3_booth$Sigma))),
                      "SCORE" = c(score(yp_booth, gp_booth$mean, gp_booth$Sigma),
                                  score(yp_booth, dgp_booth$mean, dgp_booth$Sigma),
                                  score(yp_booth, dgp3_booth$mean, dgp3_booth$Sigma)))
rownames(metrics) <- c("One-layer GP", "Two-layer DGP", "Three-layer DGP")
metrics

## ----fig.width = 6, fig.height = 4.5------------------------------------------
dgp_booth_noisy <- fit_two_layer(x_booth, y_booth, nmcmc = 10000, verb = FALSE)
dgp_booth_noisy <- trim(dgp_booth_noisy, 5000, 2)
dgp_booth_noisy <- predict(dgp_booth_noisy, xp_booth, lite = FALSE)
plot(dgp_booth_noisy)

## -----------------------------------------------------------------------------
tray <- function(x) {
  x <- x * 4 - 2
  p1 <- abs(100 - sqrt(apply(x^2, 1, sum)) / pi)
  p2 <- abs(apply(sin(x), 1, prod) * exp(p1)) + 1
  y <- -0.0001 * (p2)^(0.1)
  return((y + 1.9) / 0.2)
}

## -----------------------------------------------------------------------------
grid <- seq(0, 1, length = 30)
xp_tray <- as.matrix(expand.grid(grid, grid))
yp_tray <- tray(xp_tray)
par(mar = c(1, 1, 1, 1))
persp(grid, grid, matrix(yp_tray, nrow = length(grid)), xlab = "x1", ylab = "x2",
      zlab = "y", theta = 30, phi = 30, r = 30)

## -----------------------------------------------------------------------------
x_tray <- matrix(runif(50 * 2), ncol = 2)
y_tray <- tray(x_tray)
dgp_tray <- fit_two_layer(x_tray, y_tray, nmcmc = 10000, monowarp = TRUE, 
                          true_g = 1e-6, verb = FALSE)
dgp_tray <- trim(dgp_tray, 5000, 2)

## ----fig.height = 3.5---------------------------------------------------------
plot(dgp_tray, trace = FALSE, hidden = TRUE)

## ----fig.height = 4-----------------------------------------------------------
dgp_tray <- predict(dgp_tray, xp_tray, lite = TRUE)
plot(dgp_tray)

## -----------------------------------------------------------------------------
step <- function(x) {
  y <- pnorm((x - 0.5)/0.04)
  dy <- dnorm((x - 0.5)/0.04)/0.04
  return(list(y = y, dy = dy))
}

# Training data
x_step <- seq(0, 1, length = 6) 
y_step <- step(x_step)

# Testing data
xp_step <- seq(0, 1, length = 100)
yp_step <- step(xp_step)

par(mfrow = c(1, 2))
plot(xp_step, yp_step$y, type = "l", col = 4, xlab = "X", ylab = "Y",
     main = "Step function")
points(x_step, y_step$y)
plot(xp_step, yp_step$dy, type = "l", col = 4, xlab = "X", ylab = "Y",
     main = "Step function Gradient")
points(x_step, y_step$dy)

## -----------------------------------------------------------------------------
gp_step <- fit_one_layer(x_step, y_step$y, y_step$dy, nmcmc = 10000, 
                          true_g = 1e-6, cov = "exp2", verb = FALSE)
gp_step <- trim(gp_step, 5000, 50) # leave 100 iterations
gp_samples <- post_sample(gp_step, xp_step, grad = TRUE)

par(mfrow = c(1, 2))
matplot(xp_step, t(gp_samples$y), type = "l", xlab = "x", ylab = "y", main = "GP")
points(x_step, y_step$y, pch = 20)
matplot(xp_step, t(gp_samples$dy), type = "l", xlab = "x", ylab = "dy", main = "GP")
points(x_step, y_step$dy, pch = 20)

## -----------------------------------------------------------------------------
dgp_step <- fit_two_layer(x_step, y_step$y, y_step$dy, nmcmc = 10000, 
                          true_g = 1e-6, cov = "exp2", verb = FALSE)
dgp_step <- trim(dgp_step, 5000, 50) # leave 100 iterations
dgp_samples <- post_sample(dgp_step, xp_step, grad = TRUE)

par(mfrow = c(1, 2))
matplot(xp_step, t(dgp_samples$y), type = "l", xlab = "x", ylab = "y", main = "DGP")
points(x_step, y_step$y, pch = 20)
matplot(xp_step, t(dgp_samples$dy), type = "l", xlab = "x", ylab = "dy", main = "DGP")
points(x_step, y_step$dy, pch = 20)

## -----------------------------------------------------------------------------
dgp_booth <- predict(dgp_booth, xp_booth, EI = TRUE)

## ----fig.width = 5, fig.height = 4--------------------------------------------
plot(dgp_booth)
par(new = TRUE)
plot(xp_booth, dgp_booth$EI, type = "l", lwd = 2, col = 3, axes = FALSE, 
     xlab = "", ylab = "")
points(xp_booth[which.max(dgp_booth$EI)], max(dgp_booth$EI), pch = 17, 
       cex = 1.5, col = 3)

## ----fig.width = 5, fig.height = 4--------------------------------------------
dgp_booth <- predict(dgp_booth, xp_booth, entropy_limit = 0)
plot(dgp_booth)
par(new = TRUE)
plot(xp_booth, dgp_booth$entropy, type = "l", lwd = 2, col = 3, axes = FALSE, 
     xlab = "", ylab = "")
points(xp_booth[which.max(dgp_booth$entropy)], max(dgp_booth$entropy), pch = 17, 
       cex = 1.5, col = 3)

## ----fig.width = 5, fig.height = 4--------------------------------------------
gp_step <- fit_one_layer(x_step, y_step$y, nmcmc = 10000, true_g = 1e-6, 
                         cov = "exp2", verb = FALSE)
gp_step <- trim(gp_step, 5000, 2)
gp_step <- predict(gp_step, xp_step, lite = TRUE)
plot(gp_step)

## ----fig.width = 5, fig.height = 4--------------------------------------------
dgp_step <- fit_two_layer(x_step, y_step$y, nmcmc = 10000, true_g = 1e-6, 
                         cov = "exp2", verb = FALSE)
dgp_step <- trim(dgp_step, 5000, 2)
dgp_step <- predict(dgp_step, xp_step, lite = TRUE)
plot(dgp_step)

## ----fig.width = 7, fig.height = 4--------------------------------------------
gp_imse <- IMSE(gp_step, xp_step)
dgp_imse <- IMSE(dgp_step, xp_step)
par(mfrow = c(1, 2))
plot(xp_step, gp_imse$value, type = "l", ylab = "IMSE", main = "One-layer")
points(xp_step[which.min(gp_imse$value)], min(gp_imse$value), pch = 17, 
       cex = 1.5, col = 4)
plot(xp_step, dgp_imse$value, type = "l", ylab = "IMSE", main = "Two-layer")
points(xp_step[which.min(dgp_imse$value)], min(dgp_imse$value), pch = 17, cex = 1.5, col = 4)

## ----eval = FALSE-------------------------------------------------------------
# gp_alc <- ALC(gp_step, xp_step)
# dgp_alc <- ALC(dgp_step, xp_step)

