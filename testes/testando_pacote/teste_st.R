##-- Pacotes ----
library(TGMRF)
library(tmvtnorm)
library(truncnorm)
library(OpenMPController)

##-- Configurações ----
rowid <- 8
colid <- 5
n_vars <- 2
N <- rowid*colid*n_vars

betas <- 0.5
intercept <- -.5
nu <- 0.3

P <- length(betas)
X <- scale(matrix(c(rnorm(n = rowid*colid*n_vars, mean = 0)), ncol = P))

type_data <- 'log-normal'
family <- 'poisson'
seed <- sample(x = 1:10000, 1)#123456

rho_s <- 0.95*2.25
rho_t <- 0.00001
rho_st <- 0.00001

##-- Dados ----
data_poisson <- rtgmrf(rowid = rowid, colid = colid, X = X, n_var = n_vars,
                       rho_s = rho_s, rho_t = rho_t, rho_st = rho_st,
                       betas = betas, nu = nu, intercept = intercept,
                       type_data = type_data, family = family,
                       seed = seed)

hist(data_poisson$y, breaks = 20)

bd <- data.frame(y = data_poisson$y,
                 regiao = data_poisson$reg, grupo = data_poisson$var,
                 X1 = data_poisson$X,
                 eps = data_poisson$eps,
                 check.names = F)
neigh <- data_poisson$neigh

##-- Configurações ----
nsim <- 20000
burnin <- 0
thin <- 1
formula <- paste0('y ~ X1')

##-- Metropolis rcpp ----
vars_beta <- diag(c(0.03, 0.02), nrow = P+1, ncol = P+1)
vars_eps <- diag(0.0001, nrow = N, ncol = N)
var_log_nu <- 0.85
var_rho <- diag(c(0.1, 0.3, 0.05), ncol = 3, nrow = 3)

omp_set_num_threads(n = 1)
type_model <- type_data

c_beta <- 1#(2.38^2)/P
c_eps <- 1#(2.38^2)/N
c_nu <- 1#(2.38^2)
c_rho <- 1#(2.38^2)/3

system.time(
  out_metcpp <- tgmrf(data = bd, formula = formula,
                      spatial_var = "regiao", group_var = "grupo",
                      beta = c(intercept, betas), eps = data_poisson$eps,
                      nu = nu,
                      rho_s = 0.95, rho_t = 0, rho_st = 0,
                      family = family, type = type_model, mat_type = "car", method = "metropolis",
                      nsim = nsim, burnin = burnin, thin = thin,
                      prior_param = list("nu" = list("shape" = 0.2, "rate" = 0.2)),
                      E = NULL, neigh = neigh,
                      fix_rho = list("rho_s" = FALSE, "rho_t" = FALSE, "rho_st" = FALSE),
                      MCMC_config = list('metropolis' = list('var_beta' = vars_beta,
                                                             'var_eps' = vars_eps,
                                                             'var_log_nu' = var_log_nu,
                                                             'var_rho' = var_rho)),
                      range = list("rho_s" = c(-1, 1), "rho_t" = c(-1, 1), "rho_st" = c(-1, 1)),
                      scale = FALSE, verbose = TRUE,
                      c_beta = c_beta, c_eps = c_eps, c_nu = c_nu, c_rho = c_rho)
)

dplyr::n_distinct(out_metcpp$nu)*100/nsim
dplyr::n_distinct(out_metcpp$eps[,1])*100/nsim

acf(out_metcpp$beta, lag.max = 100)
acf(out_metcpp$rho, lag.max = 100)
acf(out_metcpp$nu, lag.max = 100)

##-- Comparando resultados ----
coef(glm(formula = formula, family = "poisson", data = bd))
coef(out_metcpp)

par(mfrow = c(1, 2))
plot.ts(out_metcpp$beta[, 1])
abline(h = intercept, col = "red")
plot.ts(out_metcpp$beta[, 2])
abline(h = betas, col = "red")

plot.ts(out_metcpp$eps[, 1:10])
par(mfrow = c(1, 2))
plot(x = colMeans(out_metcpp$eps), y = bd$eps, asp = 0)
abline(a = 0, b = 1, col = "red")
plot(x = colMeans(out_metcpp$mu), y = data_poisson$mu, asp = 0)
abline(a = 0, b = 1, col = "red")

par(mfrow = c(2, 3))
plot.ts(out_metcpp$rho[,1])
abline(h = rho_s, col = "red")
plot.ts(out_metcpp$rho[,2])
abline(h = rho_t, col = "red")
plot.ts(out_metcpp$rho[,3])
abline(h = rho_st, col = "red")

plot(density(out_metcpp$rho[,1]))
abline(v = rho_s, col = "red")
plot(density(out_metcpp$rho[,2]))
abline(v = rho_t, col = "red")
plot(density(out_metcpp$rho[,3]))
abline(v = rho_st, col = "red")

par(mfrow = c(1, 2))
plot.ts(out_metcpp$nu)
abline(h = nu, col = "red")
plot(density(out_metcpp$nu))
abline(v = nu, col = "red")

summary(out_metcpp)

