##-- Packages ----
library(tidyr)
library(spdep)
library(plot3D)
library(plotly)
library(TGMRF)

max_range_cover <- function(Ws = Ws, Wt = Wt, rho_s = NULL, rho_t = NULL, rho_st = NULL){
  res <- try(TGMRF:::max_range(Ws = Ws, Wt = Wt, rho_s = rho_s, rho_t = rho_t, rho_st = rho_st), silent = TRUE)

  if(class(res) == "try-error"){
    res <- NA
  } else{
    res <- res$max_value[1]
  }

  return(res)
}
is.pd <- function(Ws, Wt, rho_s, rho_t, rho_st){
  Q <- TGMRF:::buildQ(Ws = Ws, Wt = Wt, tau = 1, rho_s = rho_s, rho_t = rho_t, rho_st = rho_st)

  sim <- isSymmetric(Q)
  eig <- all(eigen(Q)$values > 0)

  if(sim & eig){
    return(TRUE)
  } else{
    return(FALSE)
  }
}

##-- Configurações ----
rowid <- 6
colid <- 5
n_var <- 10

neigh <- cell2nb(rowid, colid)
Ws <- nb2mat(neigh, style = 'B')
Wt <- abs(outer(1:n_var, 1:n_var, "-")) == 1

##-- Max values ----
max_rho <- TGMRF:::max_range(Ws = Ws, Wt = Wt)

seq_rho_s <- seq(-max_rho$max_value[1], max_rho$max_value[1], length.out = 30)
seq_rho_t <- seq(-max_rho$max_value[2], max_rho$max_value[2], length.out = 30)
seq_rho_st <- seq(-max_rho$max_value[3], max_rho$max_value[3], length.out = 30)

df_rho_max <- expand.grid(rho_s = seq_rho_s, rho_t = seq_rho_t)
df_rho_pd <- expand.grid(rho_s = seq_rho_s, rho_t = seq_rho_t, rho_st = seq_rho_st)

df_rho_pd$pd <- apply(X = df_rho_pd, MARGIN = 1, FUN = function(x) is.pd(Ws = Ws, Wt = Wt, rho_s = x[1], rho_t = x[2], rho_st = x[3]))

df_rho_max$rho_st <- apply(X = df_rho_max, MARGIN = 1, FUN = function(x) max_range_cover(Ws = Ws, Wt = Wt, rho_s = x[1], rho_t = x[2]))
df_rho_max$rho_st_n <- -df_rho_max$rho_st
df_rho_max <- gather(df_rho_max, key = signal, value = rho_st, -rho_s, -rho_t) %>%
  select(-signal)

plot_ly(x = df_rho_max$rho_s, y = df_rho_max$rho_t, z = df_rho_max$rho_st, mode = "scatter3d") %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'rho_s'),
                      yaxis = list(title = 'rho_t'),
                      zaxis = list(title = 'rho_st')))
