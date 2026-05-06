library(MetricGraph)
library(ggplot2)

get_letter_dict <- function() {
  
  arc <- function(cx, cy, r, from, to, n = 80) {
    th <- seq(from, to, length.out = n)
    cbind(cx + r * cos(th), cy + r * sin(th))
  }
  
  list(
    
    A = function() list(
      rbind(c(0,0), c(1,4)),
      rbind(c(1,4), c(2,0)),
      rbind(c(0.5,2), c(1.5,2))
    ),
    
    B = function() list(
      rbind(c(0,0), c(0,2)),
      rbind(c(0,2), c(0,4)),
      rbind(c(0,4), c(1,4)),
      arc(1,3,1, pi/2, -pi/2),
      rbind(c(1,2), c(0,2)),
      arc(1,1,1, pi/2, -pi/2),
      rbind(c(0,0), c(1,0))
    ),
    
    C = function() list(
      arc(2,2,2, pi/2, 3*pi/2)
    ),
    
    D = function() list(
      rbind(c(0,0), c(0,4)),
      arc(0,2,2, -pi/2, pi/2)
    ),
    
    E = function() list(
      rbind(c(0,0), c(0,4)),
      rbind(c(0,4), c(2,4)),
      rbind(c(0,2), c(1,2)),
      rbind(c(0,0), c(2,0))
    ),
    
    F = function() list(
      rbind(c(0,0), c(0,4)),
      rbind(c(0,4), c(2,4)),
      rbind(c(0,2), c(1,2))
    ),
    
    G = function() list(
      arc(2,2,2, pi/2, 3*pi/2),
      rbind(c(1,2), c(2,2)),
      rbind(c(2,2), c(2,0))
    ),
    
    H = function() list(
      rbind(c(0,0), c(0,4)),
      rbind(c(0,2), c(2,2)),
      rbind(c(2,0), c(2,4))
    ),
    
    I = function() list(
      rbind(c(0,4), c(2,4)),
      rbind(c(1,4), c(1,0)),
      rbind(c(0,0), c(2,0))
    ),
    
    J = function() list(
      rbind(c(0,4), c(2,4)),
      rbind(c(1,4), c(1,0.5)),
      arc(1.5,0.5,0.5, pi, 2*pi)
    ),
    
    K = function() list(
      rbind(c(0,0), c(0,4)),
      rbind(c(0,2), c(2,4)),
      rbind(c(0,2), c(2,0))
    ),
    
    L = function() list(
      rbind(c(0,0), c(0,4)),
      rbind(c(0,0), c(2,0))
    ),
    
    M = function() list(
      rbind(c(0,0), c(0,4)),
      rbind(c(0,4), c(1,2)),
      rbind(c(1,2), c(2,4)),
      rbind(c(2,4), c(2,0))
    ),
    
    N = function() list(
      rbind(c(0,0), c(0,4)),
      rbind(c(0,4), c(2,0)),
      rbind(c(2,0), c(2,4))
    ),
    
    O = function() list(
      rbind(c(0,1), c(0,3)),
      arc(1,1,1, pi, 2*pi),
      arc(1,3,1, 0, pi),
      rbind(c(2,1), c(2,3))
    ),
    
    P = function() list(
      rbind(c(0,0), c(0,4)),
      rbind(c(0,4), c(1,4)),
      arc(1,3,1, pi/2, -pi/2),
      rbind(c(1,2), c(0,2))
    ),
    
    Q = function() list(
      rbind(c(0,1), c(0,3)),
      arc(1,1,1, pi, 2*pi),
      arc(1,3,1, 0, pi),
      rbind(c(2,1), c(2,3)),
      rbind(c(1,1), c(1+1/sqrt(2),1-1/sqrt(2))),
      rbind(c(1+1/sqrt(2),1-1/sqrt(2)), c(2,0))
    ),
    
    R = function() list(
      rbind(c(0,0), c(0,4)),
      rbind(c(0,4), c(1,4)),
      arc(1,3,1, pi/2, -pi/2),
      rbind(c(1,2), c(0,2)),
      rbind(c(1,2), c(2,0))
    ),
    
    S = function() list(
      arc(1,3,1, pi/2, 3*pi/2),
      arc(1,1,1, -pi/2, pi/2),
      rbind(c(0,0), c(1,0)),
      rbind(c(1,4), c(2,4))
    ),
    
    T = function() list(
      rbind(c(0,4), c(2,4)),
      rbind(c(1,4), c(1,0))
    ),
    
    U = function() list(
      rbind(c(0,4), c(0,1)),
      arc(1,1,1, pi, 2*pi),
      rbind(c(2,1), c(2,4))
    ),
    
    V = function() list(
      rbind(c(0,4), c(1,0)),
      rbind(c(1,0), c(2,4))
    ),
    
    W = function() list(
      rbind(c(0,4), c(0,0)),
      rbind(c(0,0), c(1,2)),
      rbind(c(1,2), c(2,0)),
      rbind(c(2,0), c(2,4))
    ),
    
    X = function() list(
      rbind(c(0,0), c(1,2)),
      rbind(c(1,2), c(2,4)),
      rbind(c(0,4), c(1,2)),
      rbind(c(1,2), c(2,0))
    ),
    
    Y = function() list(
      rbind(c(0,4), c(1,2)),
      rbind(c(2,4), c(1,2)),
      rbind(c(1,2), c(1,0))
    ),
    
    Z = function() list(
      rbind(c(0,4), c(2,4)),
      rbind(c(2,4), c(0,0)),
      rbind(c(0,0), c(2,0))
    ),
    "?" = function() list(
      arc(1,3,1, -pi/2, pi),
      rbind(c(1,2), c(1,0))
    )
    
  )
}

shift_letter <- function(letter_edges, xt, yt) {
  lapply(letter_edges, function(e) {
    e + matrix(c(xt, yt), nrow(e), 2, byrow = TRUE)
  })
}

dict <- get_letter_dict()

make_word <- function(word, start_x = 0, start_y = 0, spacing = 2) {
  edges <- list()
  
  for (i in seq_along(strsplit(word, "")[[1]])) {
    letter <- strsplit(word, "")[[1]][i]
    
    letter_edges <- dict[[letter]]()
    xt <- start_x + (i - 1) * spacing
    yt <- start_y
    
    edges <- c(edges, shift_letter(letter_edges, xt, yt))
  }
  
  edges
}

edges <- c(
  make_word("THANKYOU", 0, 5),
  make_word("FOR", 0, 1),
  make_word("YOUR", 6, 1),
  make_word("ATTENTION", 0, -3)
)

edges <- c(make_word("CAN", 0, 0) ,make_word("I", 8, 0), make_word("WRITE", 12, 0),
           make_word("EVERYTHING", 1, -4), 
           make_word("I", 4, -8), make_word("WANT?", 8, -8))



edges <- c(make_word("WANNA", 0, 0),
           make_word("GO", 0, -4),
           make_word("FOR", 6, -4), 
           make_word("A", 0, -8), make_word("WALK?", 4, -8))

make_layout <- function(lines, offset = NULL,
                        x_start = 0, y_start = 0, 
                        x_gap = 2, y_gap = 4, word_spacing = 2) {
  
  if (is.null(offset)) {
    offset <- rep(0, length(lines))
  }
  
  if (length(offset) != length(lines)) {
    stop("offset must have same length as lines")
  }
  
  edges <- list()
  y <- y_start
  
  for (i in seq_along(lines)) {
    words <- strsplit(lines[i], " ")[[1]]
    
    x <- x_start + offset[i]   # <-- key change
    
    for (w in words) {
      edges <- c(edges, list(make_word(w, x, y)))
      x <- x + nchar(w) * word_spacing + x_gap
    }
    
    y <- y - y_gap
  }
  
  return(unlist(edges, recursive = FALSE))
}


lines <- c("CAN I WRITE", "EVERYTHING", "I WANT?")
offset <- c(0,1,4)

edges <- make_layout(
  lines = lines,
  offset = offset,
  x_start = 0,
  y_start = 0,
  x_gap = 2,
  y_gap = 4,
  word_spacing = 2
)

make_square_layout_words <- function(phrase,
                                     x_start = 0,
                                     y_start = 0,
                                     x_gap = 2,
                                     y_gap = 4,
                                     word_spacing = 2) {
  
  words <- strsplit(phrase, " ")[[1]]
  n <- length(words)
  
  # number of rows = sqrt(largest square <= n)
  k <- floor(sqrt(n))
  k <- floor(k*1.5) # 1.5 makes it more rectangular, which looks better for short phrases
  
  # distribute words evenly
  base <- floor(n / k)
  extra <- n %% k
  
  sizes <- rep(base, k)
  if (extra > 0) {
    sizes[1:extra] <- sizes[1:extra] + 1
  }
  
  # build lines
  idx <- 1
  lines <- character(k)
  
  for (i in seq_len(k)) {
    lines[i] <- paste(words[idx:(idx + sizes[i] - 1)], collapse = " ")
    idx <- idx + sizes[i]
  }
  
  # ✅ CENTER USING CHARACTER LENGTH
  line_lengths <- nchar(lines)
  max_len <- max(line_lengths)
  
  offset <- (max_len - line_lengths) / 2
  offset <- offset * word_spacing   # scale to your layout
  
  edges <- make_layout(
    lines = lines,
    offset = offset,
    x_start = x_start,
    y_start = y_start,
    x_gap = x_gap,
    y_gap = y_gap,
    word_spacing = word_spacing
  )
  
  list(
    edges = edges,
    lines = lines,
    offset = offset
  )
}

clean_phrase <- function(phrase) {
  
  # 1. remove accents (é -> e, etc.)
  phrase <- iconv(phrase, from = "UTF-8", to = "ASCII//TRANSLIT")
  
  # 2. replace anything that is NOT a letter with a space
  phrase <- gsub("[^A-Za-z]", " ", phrase)
  
  # 3. uppercase
  phrase <- toupper(phrase)
  
  # 4. collapse multiple spaces
  phrase <- gsub("\\s+", " ", phrase)
  
  # 5. trim leading/trailing spaces
  phrase <- trimws(phrase)
  
  return(phrase)
}



phrase <- "ALL I WANT IS TO WRITE THINGS IN A WAY THAT LOOKS LIKE A SQUARE THE ONLY PROBLEM IS THAT I CANNOT ADD COMMAS"


phrase <- "The Matérn covariance (or Matérn kernel) is a widely used function in spatial statistics, machine learning, and geostatistics that determines the covariance—or similarity—between two points based on the distance separating them"

phrase <- "This thesis develops a unified framework for the modeling, analysis, and computation of fractional differential operators on metric graphs, with applications to Gaussian random fields, diffusion processes, and optimal control. Metric graphs provide a natural setting for representing network-structured data, such as transportation or river systems, where classical Euclidean models are inadequate.

We first introduce a new class of non-stationary Gaussian fields with arbitrary smoothness on compact metric graphs, extending the classical Whittle-Matérn construction. The proposed model allows for spatially varying parameters and greater flexibility, while preserving key analytical properties. Regularity results are established, and an efficient computational framework based on finite element discretization and rational approximation of fractional operators is developed, enabling scalable inference.

Next, we study fractional diffusion equations on metric graphs governed by fractional powers of the shifted Kirchhoff-Laplacian. A fully discrete numerical scheme is proposed, combining backward Euler time-stepping, finite element discretization, and rational approximations to reduce the fractional operator to a sequence of sparse elliptic solves. Rigorous error estimates are derived for all discretization components and validated through experiments.

Building on this framework, we consider optimal control problems governed by fractional diffusion equations on metric graphs, subject to control constraints. A forward-backward sweep method is employed to compute the optimal control, and contractivity of the iteration is proved, ensuring convergence. Error estimates for the control variable are established across temporal, spatial, and rational discretizations.

Finally, we develop an exponentially convergent Markov approximation for Gaussian Whittle-Matérn fields with non-integer smoothness parameters. The approach expresses the field as a sum of independent Gaussian Markov random fields, enabling efficient simulation and likelihood-based inference without finite element discretization.

Overall, this work provides new theoretical insights and computational tools for fractional models on metric graphs, bridging stochastic modeling, numerical analysis, and optimal control, and enabling scalable methods for network-based applications."
phrase <- clean_phrase(phrase)


res <- make_square_layout_words(phrase)
res$lines
edges <- res$edges


graph <- metric_graph$new(edges = edges,
                          longlat = FALSE,
                          perform_merges = TRUE)

graph$prune_vertices()


p_PhDProposalAbstract <- graph$plot(vertex_size = 0.1, vertex_color = "blue") + 
  theme_bw() +
  theme(text = element_text(family = "Palatino"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),  # x-axis numbers
        axis.text.y = element_text(size = 12),
        plot.margin = margin(0, 0, 0, 0)) 

scaler <- 1.5
ggsave(here::here("data_files/PhDProposalAbstract.png"), 
       plot = p_PhDProposalAbstract, width = 9*scaler, height = 4.25*scaler, dpi = 300)

ggsave("~/Desktop/leninPresentations/data_files/PhDProposalAbstract.png", 
       plot = p_PhDProposalAbstract, width = 9*scaler, height = 4.25*scaler, dpi = 300)





















