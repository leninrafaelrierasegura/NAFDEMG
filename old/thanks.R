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
  
  # ============================================================
  # ROW 2 (bottom): y = 1..5  (touches row 1 at y=5)
  # ============================================================
  
  # P (0..2)
  eP1 <- rbind(c(0,1), c(0,5))
  eP2 <- rbind(c(0,5), c(1,5))
  eP3 <- arc(1,4,1, pi/2, -pi/2)
  eP4 <- rbind(c(1,3), c(0,3))
  
  # R (2..4)
  eR1 <- rbind(c(2,1), c(2,5))
  eR2 <- rbind(c(2,5), c(3,5))
  eR3 <- arc(3,4,1, pi/2, -pi/2)
  eR4 <- rbind(c(3,3), c(2,3))
  eR5 <- rbind(c(3,3), c(4,1))
  
  # O (4..6)
  eO5 <- rbind(c(4,2), c(4,4))
  eO6 <- arc(5,2,1, pi, 2*pi)
  eO7 <- arc(5,4,1, 0, pi)
  eO8 <- rbind(c(6,2), c(6,4))
  
  # F (6..8)
  eF1 <- rbind(c(6,1), c(6,5))
  eF2 <- rbind(c(6,5), c(8,5))
  eF3 <- rbind(c(6,3), c(7.5,3))
  
  # D (8..10)
  tt <- (7 - sqrt(17)) / 2   # ≈ 1.438
  eD1 <- rbind(c(8,1), c(8,5)) + rbind(c(tt,0), c(tt,0))
  eD3 <- arc(8,3,2, pi/2, -pi/2)
  eD3 <- eD3 + cbind(rep(tt, nrow(eD3)), rep(0, nrow(eD3)))
  
  # A (10..12) 
  eA4 <- rbind(c(10,1), c(11,5)) + rbind(c(1,0), c(1,0))
  eA5 <- rbind(c(11,5), c(12,1)) + rbind(c(1,0), c(1,0))
  eA6 <- rbind(c(10.5,3), c(11.5,3)) + rbind(c(1,0), c(1,0))
  
  # V (12..14)
  eV1 <- rbind(c(12,5), c(13,1))
  eV2 <- rbind(c(13,1), c(14,5))
  
  # I (14..16)
  eI1 <- rbind(c(14,5), c(16,5))
  eI2 <- rbind(c(15,5), c(15,1))
  eI3 <- rbind(c(14,1), c(16,1))
  
  # D (16..18)
  eD5 <- rbind(c(16,1), c(16,5))
  eD6 <- rbind(c(16,5), c(16,5))
  eD7 <- arc(16,3,2, pi/2, -pi/2)
  eD8 <- rbind(c(16,1), c(16,1))
  
  ky_bridge <- rbind(c(10,9), c(12,9)) 
  yv_bridge <- rbind(c(13,5), c(14,5))
  # ============================================================
  # EXTRA: vertical bridge (guarantees connectivity)
  # ============================================================
  e_bridge <- rbind(c(8,5), c(8,5))  # shared vertex (redundant but safe)
  
  edges <- list(
    # top
    eT1,eT2,eH1,eH2,eH3,eA1,eA2,eA3,eN1,eN2,eN3,
    eK1,eK2,eK3,eY1,eY2,eY3,eO1,eO2,eO3,eO4,eU1,eU2,eU3,
    
    # bottom
    eP1,eP2,eP3,eP4,eR1,eR2,eR3,eR4,eR5,
    eO5,eO6,eO7,eO8,eF1,eF2,eF3,
    eD1,eD3,eA4,eA5,eA6,
    eV1,eV2,eI1,eI2,eI3,eD5,eD6,eD7,eD8,
    
    # bridge
    #ky_bridge,
    yv_bridge
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


# thanks <- graph$plot(vertex_size = 2, edge_color = "darkblue", edge_width = 1)
# 
# ggsave(thanks, filename = "data_files/thanks.png", width = 8, height = 4, dpi = 300)

# graph <- graph$get_largest()
# graph$plot()
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
thanks3d <- graph$plot(data = "f",
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

thanks3d$x$attrs <- lapply(thanks3d$x$attrs, function(tr) {
  if (!is.null(tr$marker)) {
    tr$marker$showscale <- FALSE
    tr$marker$colorbar <- NULL
  }
  tr$showscale <- FALSE
  tr
})


save(thanks3d, file = here::here("old/thanks3d.RData"))
