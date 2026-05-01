library(MetricGraph)
library(ggplot2)
gets.graph.thankyou <- function() {
  
  arc <- function(cx, cy, r, from, to, n = 80) {
    th <- seq(from, to, length.out = n)
    cbind(cx + r * cos(th), cy + r * sin(th))
  }
  
  # ============================================================
  # ROW 1 (top): y = 5..9
  # Letters are contiguous: width = 2 each, no gaps
  # ============================================================
  
  # T (0..2)
  eT1 <- rbind(c(0,9), c(2,9))
  eT2 <- rbind(c(1,9), c(1,5))
  
  # H (2..4)
  eH1 <- rbind(c(2,5), c(2,9))
  eH2 <- rbind(c(2,7), c(4,7))
  eH3 <- rbind(c(4,5), c(4,9))
  
  # A (4..6)
  eA1 <- rbind(c(4,5), c(5,9))
  eA2 <- rbind(c(5,9), c(6,5))
  eA3 <- rbind(c(4.5,7), c(5.5,7))
  
  # N (6..8)
  eN1 <- rbind(c(6,5), c(6,9))
  eN2 <- rbind(c(6,9), c(8,5))
  eN3 <- rbind(c(8,5), c(8,9))
  
  # K (8..10)
  eK1 <- rbind(c(8,5), c(8,9))
  eK2 <- rbind(c(8,7), c(10,9))
  eK3 <- rbind(c(8,7), c(10,5))
  
  # Y (10..12)
  eY1 <- rbind(c(10,9), c(11,7)) + rbind(c(2,0), c(2,0))
  eY2 <- rbind(c(12,9), c(11,7)) + rbind(c(2,0), c(2,0))
  eY3 <- rbind(c(11,7), c(11,5)) + rbind(c(2,0), c(2,0))
  
  # O (12..14)
  eO1 <- rbind(c(12,6), c(12,8)) + rbind(c(2,0), c(2,0))
  eO2 <- arc(13,6,1, pi, 2*pi) 
  eO2 <- eO2 + cbind(rep(2, nrow(eO2)), rep(0, nrow(eO2)))
  eO3 <- arc(13,8,1, 0, pi)
  eO3 <- eO3 + cbind(rep(2, nrow(eO3)), rep(0, nrow(eO3)))
  eO4 <- rbind(c(14,6), c(14,8)) + rbind(c(2,0), c(2,0))
  
  # U (14..16)
  eU1 <- rbind(c(14,9), c(14,6)) + rbind(c(2,0), c(2,0))
  eU2 <- arc(15,6,1, pi, 2*pi)
  eU2 <- eU2 + cbind(rep(2, nrow(eU2)), rep(0, nrow(eU2)))
  eU3 <- rbind(c(16,6), c(16,9)) + rbind(c(2,0), c(2,0))
  
  # FOR
  forxt <- -19
  forty <- -4
  # F (18..20)
  xt <- 14 + forxt
  yt <- 4 + forty
  
  eFF1 <- rbind(c(6,1), c(6,5)) + rbind(c(xt,yt), c(xt,yt))
  eFF2 <- rbind(c(6,5), c(8,5)) + rbind(c(xt,yt), c(xt,yt))
  eFF3 <- rbind(c(6,3), c(7.5,3)) + rbind(c(xt,yt), c(xt,yt))
  
  # O
  xt <- 10 + forxt
  yt <- 0 + forty
  eOO1 <- rbind(c(12,6), c(12,8)) + rbind(c(xt,yt), c(xt,yt))
  eOO2 <- arc(13,6,1, pi, 2*pi) 
  eOO2 <- eOO2 + cbind(rep(xt, nrow(eOO2)), rep(yt, nrow(eOO2)))
  eOO3 <- arc(13,8,1, 0, pi)
  eOO3 <- eOO3 + cbind(rep(xt, nrow(eOO3)), rep(yt, nrow(eOO3)))
  eOO4 <- rbind(c(14,6), c(14,8)) + rbind(c(xt,yt), c(xt,yt))
  
  # R
  xt <- 22 + forxt
  yt <- 4 + forty
  eRR1 <- rbind(c(2,1), c(2,5)) + rbind(c(xt,yt), c(xt,yt))
  eRR2 <- rbind(c(2,5), c(3,5)) + rbind(c(xt,yt), c(xt,yt))
  eRR3 <- arc(3,4,1, pi/2, -pi/2)
  eRR3 <- eRR3 + cbind(rep(xt, nrow(eRR3)), rep(yt, nrow(eRR3)))
  eRR4 <- rbind(c(3,3), c(2,3)) + rbind(c(xt,yt), c(xt,yt))
  eRR5 <- rbind(c(3,3), c(4,1)) + rbind(c(xt,yt), c(xt,yt))
  # ============================================================
  # ROW 2 (bottom): y = 1..5  (touches row 1 at y=5)
  # ============================================================
  
  # YOUR
  yourxt <- 9
  youryt <- 0
  # Y (10..12)
  xt <- -10 + yourxt
  yt <- -4 + youryt
  eYY1 <- rbind(c(10,9), c(11,7)) + rbind(c(xt,yt), c(xt,yt))
  eYY2 <- rbind(c(12,9), c(11,7)) + rbind(c(xt,yt), c(xt,yt))
  eYY3 <- rbind(c(11,7), c(11,5)) + rbind(c(xt,yt), c(xt,yt))
  
  # O (4..6)
  xt <- -2 + yourxt
  yt <- 0 + youryt
  eO5 <- rbind(c(4,2), c(4,4)) + rbind(c(xt,yt), c(xt,yt))
  eO6 <- arc(5,2,1, pi, 2*pi)
  eO6 <- eO6 + cbind(rep(xt, nrow(eO6)), rep(yt, nrow(eO6)))
  eO7 <- arc(5,4,1, 0, pi)
  eO7 <- eO7 + cbind(rep(xt, nrow(eO7)), rep(yt, nrow(eO7)))
  eO8 <- rbind(c(6,2), c(6,4)) + rbind(c(xt,yt), c(xt,yt))
  
  # U (14..16)
  xt <- -10 + yourxt
  yt <- -4 + youryt
  eUU1 <- rbind(c(14,9), c(14,6)) + rbind(c(xt,yt), c(xt,yt))
  eUU2 <- arc(15,6,1, pi, 2*pi)
  eUU2 <- eUU2 + cbind(rep(xt, nrow(eUU2)), rep(yt, nrow(eUU2)))
  eUU3 <- rbind(c(16,6), c(16,9)) + rbind(c(xt,yt), c(xt,yt))
  
  # R
  xt <- 4 + yourxt
  yt <- 0 + youryt
  eR1 <- rbind(c(2,1), c(2,5)) + rbind(c(xt,yt), c(xt,yt))
  eR2 <- rbind(c(2,5), c(3,5)) + rbind(c(xt,yt), c(xt,yt))
  eR3 <- arc(3,4,1, pi/2, -pi/2)
  eR3 <- eR3 + cbind(rep(xt, nrow(eR3)), rep(yt, nrow(eR3)))
  eR4 <- rbind(c(3,3), c(2,3)) + rbind(c(xt,yt), c(xt,yt))
  eR5 <- rbind(c(3,3), c(4,1)) + rbind(c(xt,yt), c(xt,yt))
  
  # ATTENTION
  attxt <- -10
  attyt <- -4
  # A
  xt <- 6 + attxt
  yt <- -4 + attyt
  eAA1 <- rbind(c(4,5), c(5,9)) + rbind(c(xt,yt), c(xt,yt))
  eAA2 <- rbind(c(5,9), c(6,5)) + rbind(c(xt,yt), c(xt,yt))
  eAA3 <- rbind(c(4.5,7), c(5.5,7)) + rbind(c(xt,yt), c(xt,yt))
  
  # T
  xt <- 12 + attxt
  yt <- -4 + attyt
  eTT1 <- rbind(c(0,9), c(2,9)) + rbind(c(xt,yt), c(xt,yt))
  eTT2 <- rbind(c(1,9), c(1,5)) + rbind(c(xt,yt), c(xt,yt))
  
  # T
  xt <- 14 + attxt
  yt <- -4 + attyt
  eTTT1 <- rbind(c(0,9), c(2,9)) + rbind(c(xt,yt), c(xt,yt))
  eTTT2 <- rbind(c(1,9), c(1,5)) + rbind(c(xt,yt), c(xt,yt))
  
  # E (6..8) 
  xt <- 10 + attxt
  yt <- 0 + attyt
  eF1 <- rbind(c(6,1), c(6,5)) + rbind(c(xt,yt), c(xt,yt))
  eF2 <- rbind(c(6,5), c(8,5)) + rbind(c(xt,yt), c(xt,yt))
  eF3 <- rbind(c(6,3), c(7.5,3)) + rbind(c(xt,yt), c(xt,yt))
  eF4 <- rbind(c(6,1), c(8,1)) + rbind(c(xt,yt), c(xt,yt))
  
  # N (6..8)
  xt <- 12 + attxt
  yt <- -4 + attyt
  eNN1 <- rbind(c(6,5), c(6,9)) + rbind(c(xt,yt), c(xt,yt))
  eNN2 <- rbind(c(6,9), c(8,5)) + rbind(c(xt,yt), c(xt,yt))
  eNN3 <- rbind(c(8,5), c(8,9)) + rbind(c(xt,yt), c(xt,yt))

  # T
  xt <- 20 + attxt
  yt <- -4 + attyt
  eTTTT1 <- rbind(c(0,9), c(2,9)) + rbind(c(xt,yt), c(xt,yt))
  eTTTT2 <- rbind(c(1,9), c(1,5)) + rbind(c(xt,yt), c(xt,yt))
  
  # I (14..16)
  xt <- 8 + attxt
  yt <- 0 + attyt
  eI1 <- rbind(c(14,5), c(16,5)) + rbind(c(xt,yt), c(xt,yt))
  eI2 <- rbind(c(15,5), c(15,1)) + rbind(c(xt,yt), c(xt,yt))
  eI3 <- rbind(c(14,1), c(16,1)) + rbind(c(xt,yt), c(xt,yt))
  
  # O (4..6)
  xt <- 20 + attxt
  yt <- 0 + attyt
  eOO5 <- rbind(c(4,2), c(4,4)) + rbind(c(xt,yt), c(xt,yt))
  eOO6 <- arc(5,2,1, pi, 2*pi)
  eOO6 <- eOO6 + cbind(rep(xt, nrow(eOO6)), rep(yt, nrow(eOO6)))
  eOO7 <- arc(5,4,1, 0, pi)
  eOO7 <- eOO7 + cbind(rep(xt, nrow(eOO7)), rep(yt, nrow(eOO7)))
  eOO8 <- rbind(c(6,2), c(6,4)) + rbind(c(xt,yt), c(xt,yt))
  
  # N (6..8)
  xt <- 20 + attxt
  yt <- -4 + attyt
  eNNN1 <- rbind(c(6,5), c(6,9)) + rbind(c(xt,yt), c(xt,yt))
  eNNN2 <- rbind(c(6,9), c(8,5)) + rbind(c(xt,yt), c(xt,yt))
  eNNN3 <- rbind(c(8,5), c(8,9)) + rbind(c(xt,yt), c(xt,yt))
  
  
  #ky_bridge <- rbind(c(10,9), c(12,9)) 
  # ============================================================
  # EXTRA: vertical bridge (guarantees connectivity)
  # ============================================================
  
  edges <- list(
    # top
    eT1,eT2,eH1,eH2,eH3,eA1,eA2,eA3,eN1,eN2,eN3,
    eK1,eK2,eK3,eY1,eY2,eY3,eYY1,eYY2,eYY3,
    eO1,eO2,eO3,eO4,eOO1,eOO2,eOO3,eOO4,
    eU1,eU2,eU3,eFF1,eFF2,eFF3,
    
    # bottom
    eAA1,eAA2,eAA3,
    eRR1,eRR2,eRR3,eRR4,eRR5,
    eO5,eO6,eO7,eO8,eF1,eF2,eF3,eF4,
    eI1,eI2,eI3,
    
    eUU1,eUU2,eUU3,eTT1,eTT2,eTTT1,eTTT2,
    eR1,eR2,eR3,eR4,eR5,eNN1,eNN2,eNN3,
    eTTTT1,eTTTT2,eOO5,eOO6,eOO7,eOO8,
    eNNN1,eNNN2,eNNN3
    
    # bridge
    # ky_bridge
  )
  
  
  return(edges)
}

# ── Usage ─────────────────────────────────────────────
edges <- gets.graph.thankyou()

graph <- metric_graph$new(edges = edges,
                          longlat = FALSE,
                          perform_merges = TRUE)
#graph$prune_vertices()

#graph <- graph_components$new(edges = edges)


# thanksall <- graph$plot(vertex_size = 2, edge_color = "darkblue", edge_width = 1) +
#   theme_void() +
#   theme(plot.margin = margin(0, 0, 0, 0))
# thanksall
# ggsave(thanksall, filename = "data_files/thanksall.png", width = 9, height = 6, dpi = 300)

graph$build_mesh(h = 0.05)

lon <- graph$mesh$V[, 1]
lat <- graph$mesh$V[, 2]
f <- lat + 3

edge_number <- graph$mesh$VtE[,1]
distance_on_edge <- graph$mesh$VtE[,2]

graph$add_observations(
  data = data.frame(
    edge_number = edge_number,
    distance_on_edge = distance_on_edge,
    f = f),
  clear_obs = TRUE,
  normalized = TRUE)

library(plotly)
thanksall3d <- graph$plot(data = "f",
  vertex_size = 2, 
  type = "plotly",
  edge_color = "darkblue",
  edge_width = 1) |>
  layout(scene = list(
    camera = list(eye = list(x=1.5,y=0,z=1)),
    xaxis = list(
      title = "",
      showticklabels = FALSE,
      ticks = "",
      showgrid = FALSE,
      zeroline = FALSE,
      showbackground = FALSE
    ),
    yaxis = list(
      title = "",
      showticklabels = FALSE,
      ticks = "",
      showgrid = FALSE,
      zeroline = FALSE,
      showbackground = FALSE
    ),
    zaxis = list(
      title = "",
      showticklabels = FALSE,
      ticks = "",
      showgrid = FALSE,
      zeroline = FALSE,
      showbackground = FALSE
    )
    ))

thanksall3d$x$attrs <- lapply(thanksall3d$x$attrs, function(tr) {
  if (!is.null(tr$marker)) {
    tr$marker$showscale <- FALSE
    tr$marker$colorbar <- NULL
  }
  tr$showscale <- FALSE
  tr
})

#save(thanksall3d, file = here::here("old/thanksall3d.RData"))
save_dual_for_presentation(thanksall3d)
save_plotly_figure_fixed(thanksall3d, dpi = 600, scale = 2, viewer_change = 1)




