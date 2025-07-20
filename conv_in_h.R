## ----xaringanExtra-clipboard, echo = FALSE--------------------------------------------------------
htmltools::tagList(
  xaringanExtra::use_clipboard(
    button_text = "<i class=\"fa-solid fa-clipboard\" style=\"color: #00008B\"></i>",
    success_text = "<i class=\"fa fa-check\" style=\"color: #90BE6D\"></i>",
    error_text = "<i class=\"fa fa-times-circle\" style=\"color: #F94144\"></i>"
  ),
  rmarkdown::html_dependency_font_awesome()
)


## body .main-container {
##   max-width: 100% !important;
##   width: 100% !important;
## }
## body {
##   max-width: 100% !important;
## }
## 
## body, td {
##    font-size: 16px;
## }
## code.r{
##   font-size: 14px;
## }
## pre {
##   font-size: 14px
## }
## .custom-box {
##   background-color: #f5f7fa; /* Light grey-blue background */
##   border-color: #e1e8ed; /* Light border color */
##   color: #2c3e50; /* Dark text color */
##   padding: 15px; /* Padding inside the box */
##   border-radius: 5px; /* Rounded corners */
##   margin-bottom: 20px; /* Spacing below the box */
## }
## .caption {
##   margin: auto;
##   text-align: center;
##   margin-bottom: 20px; /* Spacing below the box */
## }

## -------------------------------------------------------------------------------------------------
# Set seed for reproducibility
set.seed(1982) 
# Set global options for all code chunks
knitr::opts_chunk$set(
  # Disable messages printed by R code chunks
  message = TRUE,    
  # Disable warnings printed by R code chunks
  warning = TRUE,    
  # Show R code within code chunks in output
  echo = TRUE,        
  # Include both R code and its results in output
  include = TRUE,     
  # Evaluate R code chunks
  eval = TRUE,       
  # Enable caching of R code chunks for faster rendering
  cache = FALSE,      
  # Align figures in the center of the output
  fig.align = "center",
  # Enable retina display for high-resolution figures
  retina = 2,
  # Show errors in the output instead of stopping rendering
  error = TRUE,
  # Do not collapse code and output into a single block
  collapse = FALSE
)
# Start the figure counter
fig_count <- 0
# Define the captioner function
captioner <- function(caption) {
  fig_count <<- fig_count + 1
  paste0("Figure ", fig_count, ": ", caption)
}
# Define the function to truncate a number to two decimal places
truncate_to_two <- function(x) {
  floor(x * 10000) / 10000
}

m1table <- rSPDE:::m1table
m2table <- rSPDE:::m2table
m3table <- rSPDE:::m3table
m4table <- rSPDE:::m4table


## -------------------------------------------------------------------------------------------------
# install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)
# inla.upgrade(testing = TRUE)
# remotes::install_github("inlabru-org/inlabru", ref = "devel")
# remotes::install_github("davidbolin/rspde", ref = "devel")
# remotes::install_github("davidbolin/metricgraph", ref = "devel")
# library(INLA)
# library(inlabru)
library(rSPDE)
library(MetricGraph)
library(grateful)

library(ggplot2)
library(reshape2)
library(plotly)


## -------------------------------------------------------------------------------------------------
# For each m and beta, this function returns c_m/b_{m+1} and the roots of rb and rc
my.get.roots <- function(order, beta) {
  mt <- get(paste0("m", order, "table"))
  rb <- rep(0, order + 1)
  rc <- rep(0, order)
  if(order == 1) {
      rc = approx(mt$beta, mt[[paste0("rc")]], beta)$y
    } else {
      rc = sapply(1:order, function(i) {
        approx(mt$beta, mt[[paste0("rc.", i)]], beta)$y
      })
    }
    rb = sapply(1:(order+1), function(i) {
      approx(mt$beta, mt[[paste0("rb.", i)]], xout = beta)$y
    })
    factor = approx(mt$beta, mt$factor, xout = beta)$y
  return(list(rb = rb, rc = rc, factor = factor))
}

# Function the polynomial coefficients in increasing order like a+bx+cx^2+...
poly.from.roots <- function(roots) {
  
  coef <- 1
  for (r in roots) {coef <- convolve(coef, c(1, -r), type = "open")}
  return(coef)
}
# Function to compute the partial fraction parameters
compute.partial.fraction.param <- function(factor, pr_roots, pl_roots, time_step, scaling) {
  
  pr_coef <- c(0, poly.from.roots(pr_roots)) 
  pl_coef <- poly.from.roots(pl_roots) 
  factor_pr_coef <- pr_coef
  pr_plus_pl_coef <- factor_pr_coef + ((scaling * time_step)/factor) * pl_coef
  res <- gsignal::residue(factor_pr_coef, pr_plus_pl_coef)
  return(list(r = res$r, p = res$p, k = res$k)) 
}
# Function to compute the fractional operator
my.fractional.operators.frac <- function(L, beta, C, scale.factor, m = 1, time_step) {
  
  C <- Matrix::Diagonal(dim(C)[1], rowSums(C)) 
  Ci <- Matrix::Diagonal(dim(C)[1], 1 / rowSums(C)) 
  I <- Matrix::Diagonal(dim(C)[1])
  L <- L / scale.factor 
  LCi <- L %*% Ci
  if(beta == 1){
    L <- L * scale.factor^beta
    return(list(Ci = Ci, C = C, LCi = LCi, L = L, m = m, beta = beta, LHS = C + time_step * L))
  } else {
    scaling <- scale.factor^beta
    roots <- my.get.roots(m, beta)
    poles_rs_k <- compute.partial.fraction.param(roots$factor, roots$rc, roots$rb, time_step, scaling)

    partial_fraction_terms <- list()
    for (i in 1:(m+1)) {partial_fraction_terms[[i]] <- (LCi - poles_rs_k$p[i] * I)/poles_rs_k$r[i]}
    partial_fraction_terms[[m+2]] <- ifelse(is.null(poles_rs_k$k), 0, poles_rs_k$k) * I
    return(list(Ci = Ci, C = C, LCi = LCi, L = L, m = m, beta = beta, partial_fraction_terms = partial_fraction_terms))
  }
}
# Function to solve the iteration
my.solver.frac <- function(obj, v){
  beta <- obj$beta
  m <- obj$m
  C <- obj$C
  Ci <- obj$Ci
  if (beta == 1){
    return(solve(obj$LHS, v))
  } else {
    partial_fraction_terms <- obj$partial_fraction_terms
    output <- partial_fraction_terms[[m+2]] %*% v
    for (i in 1:(m+1)) {output <- output + solve(partial_fraction_terms[[i]], v)}
    return(Ci %*% output)
  }
}
# Function to build a tadpole graph and create a mesh
gets.graph.tadpole <- function(h){
  edge1 <- rbind(c(0,0),c(1,0))
  theta <- seq(from=-pi,to=pi,length.out = 10000)
  edge2 <- cbind(1+1/pi+cos(theta)/pi,sin(theta)/pi)
  edges <- list(edge1, edge2)
  graph <- metric_graph$new(edges = edges, verbose = 0)
  graph$set_manual_edge_lengths(edge_lengths = c(1,2))
  graph$build_mesh(h = h)
  return(graph)
}
# Function to compute the eigenfunctions 
tadpole.eig <- function(k,graph){
x1 <- c(0,graph$get_edge_lengths()[1]*graph$mesh$PtE[graph$mesh$PtE[,1]==1,2]) 
x2 <- c(0,graph$get_edge_lengths()[2]*graph$mesh$PtE[graph$mesh$PtE[,1]==2,2]) 

if(k==0){ 
  f.e1 <- rep(1,length(x1)) 
  f.e2 <- rep(1,length(x2)) 
  f1 = c(f.e1[1],f.e2[1],f.e1[-1], f.e2[-1]) 
  f = list(phi=f1/sqrt(3)) 
  
} else {
  f.e1 <- -2*sin(pi*k*1/2)*cos(pi*k*x1/2) 
  f.e2 <- sin(pi*k*x2/2)                  
  
  f1 = c(f.e1[1],f.e2[1],f.e1[-1], f.e2[-1]) 
  
  if((k %% 2)==1){ 
    f = list(phi=f1/sqrt(3)) 
  } else { 
    f.e1 <- (-1)^{k/2}*cos(pi*k*x1/2)
    f.e2 <- cos(pi*k*x2/2)
    f2 = c(f.e1[1],f.e2[1],f.e1[-1],f.e2[-1]) 
    f <- list(phi=f1,psi=f2/sqrt(3/2))
  }
}
return(f)
}
# Function to order the vertices for plotting
plotting.order <- function(v, graph){
  edge_number <- graph$mesh$VtE[, 1]
  pos <- sum(edge_number == 1)+1
  return(c(v[1], v[3:pos], v[2], v[(pos+1):length(v)], v[2]))
}
# Function to plot in 3D
graph.plotter.3d <- function(graph, 
                             U_true, 
                             U_approx1, 
                             #U_approx2, 
                             time_seq, 
                             time_step){
  x <- graph$mesh$V[, 1]
  y <- graph$mesh$V[, 2]
  x <- plotting.order(x, graph)
  y <- plotting.order(y, graph)
  weights <- graph$mesh$weights
  
  cumsum1 <- sqrt(time_step * cumsum(t(weights) %*% (U_true - U_approx1)^2))
  #cumsum2 <- sqrt(time_step * cumsum(t(weights) %*%(U_true - U_approx2)^2))
  #cumsum3 <- sqrt(time_step * cumsum(t(weights) %*% (U_approx1 - U_approx2)^2))
  
  U_true <- apply(U_true, 2, plotting.order, graph = graph)
  U_approx1 <- apply(U_approx1, 2, plotting.order, graph = graph)
  #U_approx2 <- apply(U_approx2, 2, plotting.order, graph = graph)
  
  plot_data <- data.frame(
  x = rep(x, times = ncol(U_true)),
  y = rep(y, times = ncol(U_true)),
  z_true = as.vector(U_true),
  z_approx1 = as.vector(U_approx1),
  #z_approx2 = as.vector(U_approx2),
  frame = rep(time_seq, each = length(x)))
  
  x_range <- range(x)
  y_range <- range(y)
  z_range <- range(c(U_true, 
                     #U_approx2, 
                     U_approx1))
  
  p <- plot_ly(plot_data, frame = ~frame) %>%
  add_trace(x = ~x, y = ~y, z = ~z_true, type = "scatter3d", mode = "lines", name = "T",
            line = list(color = "green", width = 2)) %>%
  add_trace(x = ~x, y = ~y, z = ~z_approx1, type = "scatter3d", mode = "lines", name = "A1",
            line = list(color = "red", width = 2)) %>%
  # add_trace(x = ~x, y = ~y, z = ~z_approx2, type = "scatter3d", mode = "lines", name = "A2", 
  #           line = list(color = "blue", width = 2)) %>%
  layout(
    scene = list(
      xaxis = list(title = "x", range = x_range),
      yaxis = list(title = "y", range = y_range),
      zaxis = list(title = "z", range = z_range),
      aspectratio = list(x = 2.4, y = 1.2, z = 1.2),
      camera = list(eye = list(x = 1.5, y = 1.5, z = 1), center = list(x = 0, y = 0, z = 0))),
    updatemenus = list(list(type = "buttons", showactive = FALSE,
                            buttons = list(
          list(label = "Play", method = "animate",
               args = list(NULL, list(frame = list(duration = 100, redraw = TRUE), fromcurrent = TRUE))),
          list(label = "Pause", method = "animate",
               args = list(NULL, list(mode = "immediate", frame = list(duration = 0), redraw = FALSE)))
        )
      )
    ),
    title = "Time: 0"
  ) %>% plotly_build()
  for (i in seq_along(p$x$frames)) {
    t <- time_seq[i]
    err1 <- signif(cumsum1[i], 6)
    #err2 <- signif(cumsum2[i], 4)
    #err3 <- signif(cumsum3[i], 4)
    p$x$frames[[i]]$layout <- list(title = paste0("Time: ", t, 
                                                  # " | cum E=T-A2: ", err2, 
                                                  # " | cum E=A1-A2: ", err3,
                                                  " | cum E=T-A1: ", err1))
  }
  return(p)
}
# Function to plot the error at each time step
error.at.each.time.plotter <- function(graph, U_true, U_approx1, U_approx2, time_seq, time_step) {
  weights <- graph$mesh$weights
  error_at_each_time1 <- t(weights) %*% (U_true - U_approx1)^2
  error_at_each_time2 <- t(weights) %*% (U_true - U_approx2)^2
  error_between_both_approx <- t(weights) %*% (U_approx1 - U_approx2)^2
  
  error1 <- sqrt(as.double(t(weights) %*% (U_true - U_approx1)^2 %*% rep(time_step, ncol(U_true))))
  error2 <- sqrt(as.double(t(weights) %*% (U_true - U_approx2)^2 %*% rep(time_step, ncol(U_true))))
  errorb <- sqrt(as.double(t(weights) %*% (U_approx1 - U_approx2)^2 %*% rep(time_step, ncol(U_true))))
  
  fig <- plot_ly() %>% 
  add_trace(
  x = ~time_seq, y = ~error_at_each_time1, type = 'scatter', mode = 'lines+markers',
  line = list(color = 'red', width = 2),
  marker = list(color = 'red', size = 4),
  name = paste0("E=T-A1: ", sprintf("%.3e", error1))
) %>% 
  add_trace(
  x = ~time_seq, y = ~error_at_each_time2, type = 'scatter', mode = 'lines+markers',
  line = list(color = 'blue', width = 2, dash = "dot"),
  marker = list(color = 'blue', size = 4),
  name = paste0("E=T-A2: ", sprintf("%.3e", error2))
) %>% 
  add_trace(
  x = ~time_seq, y = ~error_between_both_approx, type = 'scatter', mode = 'lines+markers',
  line = list(color = 'orange', width = 2),
  marker = list(color = 'orange', size = 4),
  name = paste0("E=A1-A2: ", sprintf("%.3e", errorb))
) %>% 
  layout(
  title = "Error at Each Time Step",
  xaxis = list(title = "Error at Each Time Step"),
  yaxis = list(title = "Error"),
  legend = list(x = 0.1, y = 0.9)
)
  return(fig)
}
# Function to compute the eigenvalues and eigenfunctions
gets.eigen.params <- function(N_finite = 4, kappa = 1, alpha = 0.5, graph){
  EIGENVAL_ALPHA <- NULL
  EIGENFUN <- NULL
  for (j in 0:N_finite) {
      lambda_j_alpha <- (kappa^2 + (j*pi/2)^2)^(alpha/2)
      e_j <- tadpole.eig(j,graph)$phi
      EIGENVAL_ALPHA <- c(EIGENVAL_ALPHA, lambda_j_alpha)         
      EIGENFUN <- cbind(EIGENFUN, e_j)
      if (j>0 && (j %% 2 == 0)) {
        lambda_j_alpha <- (kappa^2 + (j*pi/2)^2)^(alpha/2)
        e_j <- tadpole.eig(j,graph)$psi
        EIGENVAL_ALPHA <- c(EIGENVAL_ALPHA, lambda_j_alpha)         
        EIGENFUN <- cbind(EIGENFUN, e_j)
      }
  }
  return(list(EIGENVAL_ALPHA = EIGENVAL_ALPHA,
              EIGENFUN = EIGENFUN))
}

g_linear <- function(r, A, lambda_j_alpha_half) {
  return(A * exp(-lambda_j_alpha_half * r))
  }
G_linear <- function(t, A) {
  return(A * t)
  }
g_exp <- function(r, A, mu) {
  return(A * exp(mu * r))
  }
G_exp <- function(t, A, lambda_j_alpha_half, mu) {
  exponent <- lambda_j_alpha_half + mu
  return(A * (exp(exponent * t) - 1) / exponent)
  }
g_poly <- function(r, A, n) {
  return(A * r^n)
}
G_poly <- function(t, A, lambda_j_alpha_half, n) {
  t <- as.vector(t)
  k_vals <- 0:n
  sum_term <- sapply(t, function(tt) {
    sum(((-lambda_j_alpha_half * tt)^k_vals) / factorial(k_vals))
  })
  coeff <- ((-1)^(n + 1)) * factorial(n) / (lambda_j_alpha_half^(n + 1))
  return(A * coeff * (1 - exp(lambda_j_alpha_half * t) * sum_term))
}
g_sin <- function(r, A, omega) {
  return(A * sin(omega * r))
}
G_sin <- function(t, A, lambda_j_alpha_half, omega) {
  denom <- lambda_j_alpha_half^2 + omega^2
  numerator <- exp(lambda_j_alpha_half * t) * (lambda_j_alpha_half * sin(omega * t) - omega * cos(omega * t)) + omega
  return(A * numerator / denom)
}
g_cos <- function(r, A, theta) {
  return(A * cos(theta * r)) 
}
G_cos <- function(t, A, lambda_j_alpha_half, theta) {
  denom <- lambda_j_alpha_half^2 + theta^2
  numerator <- exp(lambda_j_alpha_half * t) * (lambda_j_alpha_half * cos(theta * t) + theta * sin(theta * t)) - lambda_j_alpha_half
  return(A * numerator / denom)
}
construct_piecewise_projection <- function(projected_U_approx, time_seq, overkill_time_seq) {
  projected_U_piecewise <- matrix(NA, nrow = nrow(projected_U_approx), ncol = length(overkill_time_seq))
  
  # Assign value at t = 0
  projected_U_piecewise[, which(overkill_time_seq == 0)] <- projected_U_approx[, 1]
  
  # Assign values for intervals (t_{k-1}, t_k]
  for (k in 2:length(time_seq)) {
    idxs <- which(overkill_time_seq > time_seq[k - 1] & overkill_time_seq <= time_seq[k])
    projected_U_piecewise[, idxs] <- projected_U_approx[, k]
  }
  
  return(projected_U_piecewise)
}
loglog_line_equation <- function(x1, y1, slope) {
  b <- log10(y1 / (x1 ^ slope))
  
  function(x) {
    (x ^ slope) * (10 ^ b)
  }
}
compute_guiding_lines <- function(h_vector, alpha_vector, errors_projected, theoretical_rate) {
  guiding_lines <- matrix(NA, nrow = length(h_vector), ncol = length(alpha_vector))
  
  for (j in seq_along(alpha_vector)) {
    guiding_lines_aux <- matrix(NA, nrow = length(h_vector), ncol = length(h_vector))
    
    for (k in seq_along(h_vector)) {
      point_x1 <- h_vector[k]
      point_y1 <- errors_projected[k, j]
      slope <- theoretical_rate[j]
      
      line <- loglog_line_equation(x1 = point_x1, y1 = point_y1, slope = slope)
      guiding_lines_aux[, k] <- line(h_vector)
    }
    
    guiding_lines[, j] <- rowMeans(guiding_lines_aux)
  }
  
  return(guiding_lines)
}



## -------------------------------------------------------------------------------------------------
# Parameters
T_final <- 2
kappa <- 15
N_finite = 4 # choose even
adjusted_N_finite <- N_finite + N_finite/2 + 1
# Coefficients for u_0 and f
coeff_U_0 <- 50*(1:adjusted_N_finite)^-1
coeff_U_0[-5] <- 0
coeff_FF <- rep(0, adjusted_N_finite)
coeff_FF[7] <- 10

AAA = 1
MU = 0.3
NNN = 2
OMEGA = pi
THETA = pi

# Time step and mesh size
POWERS <- seq(14, 9, by = -1)

time_steps <- 0.1 * 2^-POWERS
h_vector <- time_steps^(1/2)

# Overkill parameters
overkill_time_step <- 0.1 * 2^-14
overkill_h <- (0.1 * 2^-14)^(1/2)

# Finest time and space mesh
overkill_time_seq <- seq(0, T_final, length.out = ((T_final - 0) / overkill_time_step + 1))
overkill_graph <- gets.graph.tadpole(h = overkill_h)

# Compute the weights on the finest mesh
overkill_graph$compute_fem() # This is needed to compute the weights
overkill_weights <- overkill_graph$mesh$weights

m_values <- c()
alpha_vector <- seq(1, 1.8, by = 0.2)
# Create a matrix to store the errors


## -------------------------------------------------------------------------------------------------
# Load the errors_projected data
load(here::here("data_files/errors_projected_h.RData"))


## ----fig.align='center', fig.dim= c(4,5), fig.cap = captioner("Caption")--------------------------
observed_rates <- numeric(length(alpha_vector))
for (u in 1:length(alpha_vector)) {
  observed_rates[u] <- coef(lm(log10(errors_projected[, u]) ~ log10(h_vector)))[2]
}

theoretical_rates <- alpha_vector

guiding_lines <- compute_guiding_lines(h_vector, alpha_vector, errors_projected, theoretical_rate = alpha_vector)
default_colors <- scales::hue_pal()(length(alpha_vector))
plot_lines <- lapply(1:ncol(guiding_lines), function(i) {
  geom_line(
    data = data.frame(x = h_vector, y = guiding_lines[, i]),
    aes(x = x, y = y),
    color = default_colors[i],
    linetype = "dashed",
    show.legend = FALSE
  )
})

df <- as.data.frame(cbind(h_vector, errors_projected))
colnames(df) <- c("h_vector", alpha_vector)
df_melted <- melt(df, id.vars = "h_vector", variable.name = "column", value.name = "value")
alpha_formatted <- formatC(alpha_vector, format = "f", digits = 2)

custom_labels <- paste0(alpha_formatted, 
                        " | ", 
                        formatC(theoretical_rates, format = "f", digits = 4), 
                        " | ", 
                        formatC(observed_rates, format = "f", digits = 4))

df_melted$column <- factor(df_melted$column, levels = alpha_vector, labels = custom_labels)

p <- ggplot() +
  geom_line(data = df_melted, aes(x = h_vector, y = value, color = column)) +
  geom_point(data = df_melted, aes(x = h_vector, y = value, color = column)) +
  plot_lines +
  labs(
    title = expression("Convergence in " * italic(h)),
    x = expression(italic(h)~(log[10]~scale)),
    y = expression(Error~(log[10]~scale)),
    color = "          α  | theo  | obs"
  ) +
  scale_x_log10(breaks = h_vector, labels = round(h_vector, 4)) +
  scale_y_log10(labels = function(x) formatC(x, format = "f", digits = 4)) +
  theme_minimal() +
  theme(text = element_text(family = "Palatino"),
        legend.position = "bottom",
        legend.direction = "vertical",
        plot.margin = margin(0, 0, 0, 0),
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))
p


## -------------------------------------------------------------------------------------------------
ggsave(here::here("data_files/conv_rates_h.png"), width = 4, height = 5, plot = p, dpi = 300)


## -------------------------------------------------------------------------------------------------
cite_packages(output = "paragraph", out.dir = ".")

