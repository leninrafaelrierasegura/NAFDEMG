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
      rbind(c(0,0), c(0,4)),
      arc(0.8,3,0.8, pi/2, -pi/2),
      arc(0.8,1,0.8, pi/2, -pi/2)
    ),
    
    C = function() list(
      arc(1,2,1.2, pi/4, 7*pi/4)
    ),
    
    D = function() list(
      rbind(c(0,0), c(0,4)),
      arc(0.8,2,1.2, pi/2, -pi/2)
    ),
    
    E = function() list(
      rbind(c(0,0), c(0,4)),
      rbind(c(0,4), c(2,4)),
      rbind(c(0,2), c(1.5,2)),
      rbind(c(0,0), c(2,0))
    ),
    
    F = function() list(
      rbind(c(0,0), c(0,4)),
      rbind(c(0,4), c(2,4)),
      rbind(c(0,2), c(1.5,2))
    ),
    
    G = function() list(
      arc(1,2,1.2, pi/4, 7*pi/4),
      rbind(c(1,2), c(2,2))
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
      rbind(c(1,4), c(1,1)),
      arc(0.5,1,0.5, 0, pi)
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
      arc(1,3,1, pi/2, -pi/2)
    ),
    
    Q = function() list(
      rbind(c(0,1), c(0,3)),
      arc(1,1,1, pi, 2*pi),
      arc(1,3,1, 0, pi),
      rbind(c(2,1), c(2,3)),
      rbind(c(1.2,1.2), c(2,0))
    ),
    
    R = function() list(
      rbind(c(0,0), c(0,4)),
      rbind(c(0,4), c(1,4)),
      arc(1,3,1, pi/2, -pi/2),
      rbind(c(1,2), c(0,2)),
      rbind(c(1,2), c(2,0))
    ),
    
    S = function() list(
      arc(1,3,1, pi, 2*pi),
      arc(1,1,1, 0, pi)
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
      rbind(c(0,4), c(0.5,0)),
      rbind(c(0.5,0), c(1,3)),
      rbind(c(1,3), c(1.5,0)),
      rbind(c(1.5,0), c(2,4))
    ),
    
    X = function() list(
      rbind(c(0,0), c(2,4)),
      rbind(c(0,4), c(2,0))
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

edges <- c(make_word("ABCDEFGHI", 0, 8),
           make_word("JKLMNOPQ", 0, 4),
           make_word("RSTUVWXYZ", 0, 0))

graph <- metric_graph$new(edges = edges,
                          longlat = FALSE,
                          perform_merges = TRUE)

graph$plot()

