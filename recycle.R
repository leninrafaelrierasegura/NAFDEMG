tau_from_h <- function(h, alpha_star) {
  
  stopifnot(h > 0, h < 1)
  stopifnot(alpha_star > 0)
  
  a <- h^alpha_star
  
  f <- function(tau) {
    tau - a * abs(log(tau))
  }
  
  uniroot(f, interval = c(.Machine$double.eps, 1 - 1e-12))$root
}

tau_from_h(0.1, 1)
# 0.1776866

tau_from_h(0.01, 1)
# 0.0567167

tau_from_h(1e-4, 1.5)
# 0.000924...

h <- 0.01
alpha_star <- 1.5

tau <- tau_from_h(h, alpha_star)

c(
  tau = tau,
  rhs = h^alpha_star * abs(log(tau)),
  error = tau - h^alpha_star * abs(log(tau))
)

library(lamW)

tau_from_h_exact <- function(h, alpha_star) {
  return(h^alpha_star * lamW::lambertW0(h^(-alpha_star)))
}




tau_from_h(0.1, 1)
tau_from_h_exact(0.1,1)


h <- 1e-8
alpha_star <- 1.5

tau <- tau_from_h(h, alpha_star)
TAU <- tau_from_h_exact(h, alpha_star)

c(
  tau = tau,
  TAU = TAU,
  rhs = h^alpha_star * abs(log(tau)),
  RHS = h^alpha_star * abs(log(TAU)),
  error = tau - h^alpha_star * abs(log(tau)),
  ERROR = TAU - h^alpha_star * abs(log(TAU))
)


solve_htau <- function(alpha, m, tol = 1e-15) {
  
  stopifnot(alpha > 0, alpha < 2)
  stopifnot(m > 0)
  
  C <- exp(-2*pi*sqrt((1 - alpha/2) * m))
  
  f <- function(h) {
    -log(C * h^(alpha - 2)) - C / h^2
  }
  
  ## Find a valid bracket automatically
  h.grid <- 10^seq(-16, -0.01, length.out = 500)
  
  vals <- f(h.grid)
  ind <- which(diff(sign(vals)) != 0)
  
  if(length(ind) == 0)
    stop("No root found.")
  
  h <- uniroot(
    f,
    lower = h.grid[ind[1]],
    upper = h.grid[ind[1] + 1],
    tol = tol
  )$root
  
  tau <- C * h^(alpha - 2)
  
  ## residuals of both equations
  err1 <- abs(tau - h^alpha * abs(log(tau)))
  err2 <- abs(tau - C * h^(alpha - 2))
  
  list(
    h = h,
    tau = tau,
    error1 = err1,
    error2 = err2
  )
}
sol <- solve_htau(alpha = 1.5, m = 30)

sol$h
sol$tau
sol$error1
sol$error2



library(Rmpfr)

solve_htau_mpfr <- function(alpha, m,
                            precBits = 512,
                            tol = 1e-80,
                            maxit = 100) {
  
  ## MPFR constants
  alpha <- mpfr(alpha, precBits)
  m     <- mpfr(m, precBits)
  
  one  <- mpfr(1, precBits)
  two  <- mpfr(2, precBits)
  piMP <- Const("pi", precBits)
  
  ## C = exp(-2*pi*sqrt((1-alpha/2)m))
  C <- exp(-two*piMP * sqrt((one - alpha/two) * m))
  
  ## Initial guess (double precision formula)
  Cd <- asNumeric(C)
  ad <- asNumeric(alpha)
  
  f0 <- function(h)
    -log(Cd * h^(ad - 2)) - Cd/h^2
  
  h0 <- uniroot(f0, c(1e-16, 0.999))$root
  
  ## Convert to MPFR
  h <- mpfr(h0, precBits)
  
  ## Newton iterations
  for (k in seq_len(maxit)) {
    
    f  <- -log(C * h^(alpha - two)) - C/h^2
    df <- -(alpha - two)/h + two*C/h^3
    
    h.new <- h - f/df
    
    if (abs(h.new - h) < mpfr(tol, precBits)) {
      h <- h.new
      break
    }
    
    h <- h.new
  }
  
  tau <- C * h^(alpha - two)
  
  err1 <- abs(tau - h^alpha * abs(log(tau)))
  err2 <- abs(tau - C * h^(alpha - two))
  
  list(
    h = h,
    tau = tau,
    error1 = err1,
    error2 = err2,
    iterations = k
  )
}

sol <- solve_htau_mpfr(
  alpha = 1.5,
  m = 40,
  precBits = 512,
  tol = 1e-100
)

formatMpfr(sol$h, digits = 100)
formatMpfr(sol$tau, digits = 100)

formatMpfr(sol$error1, digits = 30)
formatMpfr(sol$error2, digits = 30)

sol$iterations