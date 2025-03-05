## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(rollupTree)
library(massProps)
suppressPackageStartupMessages({library(igraph)})

## -----------------------------------------------------------------------------
test_table

## -----------------------------------------------------------------------------
test_tree

## ----echo = FALSE-------------------------------------------------------------
plot(test_tree,layout=layout_as_tree(test_tree, 2, mode="in"), vertex.shape = 'none', edge.arrow.mode = 0)

## -----------------------------------------------------------------------------
rollup_mass_props(test_tree, test_table)

## ----echo = FALSE-------------------------------------------------------------
na_mass_props_and_unc <- function(d, t, v) {
  xyz <- c("x", "y", "z")
  list(
    mass = NA,
    center_mass = c(x = NA, y = NA, z = NA),
    inertia = matrix(nrow = 3, ncol = 3, dimnames = list(xyz, xyz)),
    poi_conv = "+",
    point = FALSE,
    sigma_mass = NA,
    sigma_center_mass = c(x = NA, y = NA, z = NA),
    sigma_inertia = matrix(nrow = 3, ncol = 3, dimnames = list(xyz, xyz))
  )
}
na_mass_props_and_unc_update <- function(d, t, s) {
  update_mass_props_and_unc(d, t, s, override = na_mass_props_and_unc)
}
sawe_input <- rollup(sawe_tree, sawe_table, update = na_mass_props_and_unc_update, validate_ds = validate_mass_props_and_unc_table)

## -----------------------------------------------------------------------------
sawe_input

## -----------------------------------------------------------------------------
rollup_mass_props_and_unc(sawe_tree, sawe_input)

## ----echo = FALSE-------------------------------------------------------------
sawe_input <- sawe_table[which(sawe_table$id != "Combined"), ]
sawe_input

## -----------------------------------------------------------------------------
t <- rollup_radii_of_gyration_unc(sawe_tree,
  add_radii_of_gyration(
    rollup_mass_props_and_unc(sawe_tree, sawe_table)
  )
)
sawe_result <- t[t$id == "Combined", ]
sawe_result

## ----echo = FALSE-------------------------------------------------------------
agree <- function(l) if (l) "agrees" else stop("disagreement")

## -----------------------------------------------------------------------------
mass <- sum(sawe_input$mass)

## ----echo = FALSE-------------------------------------------------------------
mass

## -----------------------------------------------------------------------------
C <- apply(sawe_input$mass / mass * sawe_input[, c("Cx", "Cy", "Cz")], 2, sum)

## ----echo = FALSE-------------------------------------------------------------
C

## -----------------------------------------------------------------------------
moi <- function(I, v1, v2, m, c1, c2) {
  sum(I + m * ((v1^2 + v2^2) - (c1^2 + c2^2)))
}
MOI <- c(
  Ixx = moi(sawe_input$Ixx, sawe_input$Cy, sawe_input$Cz, sawe_input$mass, C["Cy"], C["Cz"]),
  Iyy = moi(sawe_input$Iyy, sawe_input$Cx, sawe_input$Cz, sawe_input$mass, C["Cx"], C["Cz"]),
  Izz = moi(sawe_input$Izz, sawe_input$Cx, sawe_input$Cy, sawe_input$mass, C["Cx"], C["Cy"])
)

## ----echo = FALSE-------------------------------------------------------------
MOI

## -----------------------------------------------------------------------------
poi <- function(I, v1, v2, m, c1, c2) {
  sum(I + m * (v1 * v2 - c1 * c2))
}
POI <- c(
  Ixy = poi(sawe_input$Ixy, sawe_input$Cx, sawe_input$Cy, sawe_input$mass, C["Cx"], C["Cy"]),
  Ixz = poi(sawe_input$Ixz, sawe_input$Cx, sawe_input$Cz, sawe_input$mass, C["Cx"], C["Cz"]),
  Iyz = poi(sawe_input$Iyz, sawe_input$Cy, sawe_input$Cz, sawe_input$mass, C["Cy"], C["Cz"])
)

## ----echo = FALSE-------------------------------------------------------------
POI

## -----------------------------------------------------------------------------
rog <- function(I, m) sqrt(I / m)
ROG <- c(
  kx = rog(sawe_result$Ixx, sawe_result$mass),
  ky = rog(sawe_result$Iyy, sawe_result$mass),
  kz = rog(sawe_result$Izz, sawe_result$mass)
)

## ----echo = FALSE-------------------------------------------------------------
ROG

## -----------------------------------------------------------------------------
sigma_mass <- sqrt(sum(sawe_input$sigma_mass^2))

## ----echo = FALSE-------------------------------------------------------------
sigma_mass

## -----------------------------------------------------------------------------
sigma_cm <- function(m, sigma_v, sigma_m, v, c, mass) {
  sqrt(sum((m * sigma_v)^2 + (sigma_m * (v - c))^2)) / mass
}
sigma_C <- c(
  sigma_Cx = sigma_cm(sawe_input$mass, sawe_input$sigma_Cx, sawe_input$sigma_mass, sawe_input$Cx, C["Cx"], mass),
  sigma_Cy = sigma_cm(sawe_input$mass, sawe_input$sigma_Cy, sawe_input$sigma_mass, sawe_input$Cy, C["Cy"], mass),
  sigma_Cz = sigma_cm(sawe_input$mass, sawe_input$sigma_Cz, sawe_input$sigma_mass, sawe_input$Cz, C["Cz"], mass)
)

## ----echo = FALSE-------------------------------------------------------------
sigma_C

## -----------------------------------------------------------------------------
sigma_moi <- function(sigma_I, mass, sigma_mass, v1, v2, c1, c2, sigma_v1, sigma_v2) {
  sqrt(sum(
    sigma_I^2 +
    (2 * mass * (v1 - c1) * sigma_v1)^2 +
    (2 * mass * (v2 - c2) * sigma_v2)^2 +
    (((v1 - c1)^2 + (v2 - c2)^2) * sigma_mass)^2
  ))
}
sigma_MOI <- c(
  sigma_Ixx = sigma_moi(sawe_input$sigma_Ixx, sawe_input$mass, sawe_input$sigma_mass, sawe_input$Cy,
                                  sawe_input$Cz, C["Cy"], C["Cz"], sawe_input$sigma_Cy, sawe_input$sigma_Cz),
  sigma_Iyy = sigma_moi(sawe_input$sigma_Iyy, sawe_input$mass, sawe_input$sigma_mass, sawe_input$Cx,
                                  sawe_input$Cz, C["Cx"], C["Cz"], sawe_input$sigma_Cx, sawe_input$sigma_Cz),
  sigma_Izz = sigma_moi(sawe_input$sigma_Izz, sawe_input$mass, sawe_input$sigma_mass, sawe_input$Cx,
                                  sawe_input$Cy, C["Cx"], C["Cy"], sawe_input$sigma_Cx, sawe_input$sigma_Cy)
)

## ----echo = FALSE-------------------------------------------------------------
sigma_MOI

## -----------------------------------------------------------------------------
sigma_poi <- function(sigma_I, mass, sigma_mass, v1, v2, c1, c2, sigma_v1, sigma_v2) {
  sqrt(sum(
    sigma_I^2 +
    ((v1 - c1) * mass * sigma_v2)^2 +
    ((v1 - c1) * (v2 - c2) * sigma_mass)^2 +
    ((v2 - c2) * mass * sigma_v1)^2
  ))
}
sigma_POI <- c(
  sigma_Ixy = sigma_poi(sawe_input$sigma_Ixy, sawe_input$mass, sawe_input$sigma_mass, sawe_input$Cx,
                                  sawe_input$Cy, C["Cx"], C["Cy"], sawe_input$sigma_Cx, sawe_input$sigma_Cy),
  sigma_Ixz = sigma_poi(sawe_input$sigma_Ixz, sawe_input$mass, sawe_input$sigma_mass, sawe_input$Cx,
                                  sawe_input$Cz, C["Cx"], C["Cz"], sawe_input$sigma_Cx, sawe_input$sigma_Cz),
  sigma_Iyz = sigma_poi(sawe_input$sigma_Iyz, sawe_input$mass, sawe_input$sigma_mass, sawe_input$Cy,
                                  sawe_input$Cz, C["Cy"], C["Cz"], sawe_input$sigma_Cy, sawe_input$sigma_Cz)
)

## ----echo = FALSE-------------------------------------------------------------
sigma_POI

## -----------------------------------------------------------------------------
sigma_rog <- function(mt, It, sigma_mt, sigma_It, v1, c1, v2, c2, sigma_m) {
  sqrt(sigma_It^2 / (4 * mt * It) + (It * sigma_mt^2) / (4 * mt^3) - sum(sigma_m^2 * ((v1 - c1)^2 + (v2 - c2)^2)) / (2 * mt^2))
}
sigma_ROG <- c(
  sigma_kx = sigma_rog(sawe_result$mass, sawe_result$Ixx, sawe_result$sigma_mass, sawe_result$sigma_Ixx,
                       sawe_input$Cy, sawe_result$Cy, sawe_input$Cz, sawe_result$Cz, sawe_input$sigma_mass),
  sigma_ky = sigma_rog(sawe_result$mass, sawe_result$Iyy, sawe_result$sigma_mass, sawe_result$sigma_Iyy,
                       sawe_input$Cx, sawe_result$Cx, sawe_input$Cz, sawe_result$Cz, sawe_input$sigma_mass),
  sigma_kz = sigma_rog(sawe_result$mass, sawe_result$Izz, sawe_result$sigma_mass, sawe_result$sigma_Izz,
                       sawe_input$Cx, sawe_result$Cx, sawe_input$Cy, sawe_result$Cy, sawe_input$sigma_mass)
)

## ----echo = FALSE-------------------------------------------------------------
sigma_ROG

## ----echo = FALSE-------------------------------------------------------------
result <- rollup_radii_of_gyration_unc(sawe_tree, add_radii_of_gyration(rollup_mass_props_and_unc(sawe_tree, sawe_table)))
parts <- result[1:2, ]
aggr <- result[3, ]
# does not agree
sigma_rog_sawe<- function(k, kt, v1, c1, v2, c2, sigma_m, m, mt, sigma_v1, sigma_v2, sigma_k) {
  sqrt(sum(
    ((k^2 - kt^2 + (v1 - c1)^2 + (v2 - c2)^2) * sigma_m)^2 +
      (2 * m * (v1 - c1) * sigma_v1)^2 +
      (2 * m * (v2 - c2) * sigma_v2)^2 +
      (2 * m * k * sigma_k)^2
  )) /  (2 * kt * mt)
}
sigma_ROG_sawe <- c(
  sigma_kx = sigma_rog_sawe(parts$kx, aggr$kx, parts$Cy, aggr$Cy, parts$Cz, aggr$Cz,
                            parts$sigma_mass, parts$mass, aggr$mass, parts$sigma_Cy, parts$sigma_Cz, parts$sigma_kx),
  sigma_ky = sigma_rog_sawe(parts$ky, aggr$ky, parts$Cx, aggr$Cx, parts$Cz, aggr$Cz,
                            parts$sigma_mass, parts$mass, aggr$mass, parts$sigma_Cx, parts$sigma_Cz, parts$sigma_ky),
  sigma_kz = sigma_rog_sawe(parts$kz, aggr$kz, parts$Cx, aggr$Cx, parts$Cy, aggr$Cy,
                            parts$sigma_mass, parts$mass, aggr$mass, parts$sigma_Cx, parts$sigma_Cy, parts$sigma_kz)
)
# agrees
alt_moi_sawe <- function(m, mt, k, v1, c1, v2, c2) {
  sum(m * (k^2 + v1^2 + v2^2)) - mt * (c1^2 + c2^2)
}
alt_MOI_sawe <- c(
  Ixx = alt_moi_sawe(parts$mass, aggr$mass, parts$kx, parts$Cy, aggr$Cy, parts$Cz, aggr$Cz),
  Iyy = alt_moi_sawe(parts$mass, aggr$mass, parts$ky, parts$Cx, aggr$Cx, parts$Cz, aggr$Cz),
  Izz = alt_moi_sawe(parts$mass, aggr$mass, parts$kz, parts$Cx, aggr$Cx, parts$Cy, aggr$Cy)
)
# does not agree
alt_sigma_moi_sawe <- function(k, v1, c1, v2, c2, sigma_m, m, mt, sigma_v1, sigma_v2, sigma_k) {
  sqrt(sum(
    ((k^2 + (v1 - c1)^2 + (v2 - c2)^2) * sigma_m)^2 +
      (2 * m * (v1 - c1) * sigma_v1)^2 +
      (2 * m * (v2 - c2) * sigma_v2)^2 +
      (2 * m * k * sigma_k)^2
  ))
}
alt_sigma_MOI_sawe <- c(
  sigma_Ixx = alt_sigma_moi_sawe(parts$kx, parts$Cy, aggr$Cy, parts$Cz, aggr$Cz,
                            parts$sigma_mass, parts$mass, aggr$mass, parts$sigma_Cy, parts$sigma_Cz, parts$sigma_kx),
  sigma_Iyy = alt_sigma_moi_sawe(parts$ky, parts$Cx, aggr$Cx, parts$Cz, aggr$Cz,
                            parts$sigma_mass, parts$mass, aggr$mass, parts$sigma_Cx, parts$sigma_Cz, parts$sigma_ky),
  sigma_Izz = alt_sigma_moi_sawe(parts$kz, parts$Cx, aggr$Cx, parts$Cy, aggr$Cy,
                            parts$sigma_mass, parts$mass, aggr$mass, parts$sigma_Cx, parts$sigma_Cy, parts$sigma_kz)
)

## ----echo = FALSE-------------------------------------------------------------
mp_tree_depth <- max(dfs(mp_tree, 2, mode = "in", order=FALSE, dist=TRUE)$dist) + 1
nv <- length(igraph::V(mp_tree))
ne <- length(igraph::E(mp_tree))
nl <- length(which(igraph::degree(mp_tree, mode="in") == 0))
ni <- (nv - 1) * 20
no <- (nv - nl) * 20

## ----eval=FALSE---------------------------------------------------------------
# benchmark('mp + unc no validation' = rollup_mass_props_and_unc_fast(mp_tree, mp_table),
#           'mp + unc    validation' = rollup_mass_props_and_unc(mp_tree, mp_table),
#           'mp       no validation' = rollup_mass_props_fast(mp_tree, mp_table),
#           'mp          validation' = rollup_mass_props(mp_tree, mp_table),
#           replications = 1,
#           columns = c("test", "replications", "elapsed", "user.self", "sys.self")
# )

