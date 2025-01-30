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
sawe_result <- rollup_mass_props_and_unc(sawe_tree, sawe_table)[3, ]
sawe_result

## ----echo = FALSE-------------------------------------------------------------
agree <- function(l) if (l) "agrees" else "does not agree"

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
  I + m * ((v1^2 + v2^2) - (c1^2 + c2^2))
}
MOI <- c(
  Ixx = sum(moi(sawe_input$Ixx, sawe_input$Cy, sawe_input$Cz, sawe_input$mass, C["Cy"], C["Cz"])),
  Iyy = sum(moi(sawe_input$Iyy, sawe_input$Cx, sawe_input$Cz, sawe_input$mass, C["Cx"], C["Cz"])),
  Izz = sum(moi(sawe_input$Izz, sawe_input$Cx, sawe_input$Cy, sawe_input$mass, C["Cx"], C["Cy"]))
)

## ----echo = FALSE-------------------------------------------------------------
MOI

## -----------------------------------------------------------------------------
poi <- function(I, v1, v2, m, c1, c2) {
  I + m * (v1 * v2 - c1 * c2)
}
POI <- c(
  Ixy = sum(poi(sawe_input$Ixy, sawe_input$Cx, sawe_input$Cy, sawe_input$mass, C["Cx"], C["Cy"])),
  Ixz = sum(poi(sawe_input$Ixz, sawe_input$Cx, sawe_input$Cz, sawe_input$mass, C["Cx"], C["Cz"])),
  Iyz = sum(poi(sawe_input$Iyz, sawe_input$Cy, sawe_input$Cz, sawe_input$mass, C["Cy"], C["Cz"]))
)

## ----echo = FALSE-------------------------------------------------------------
POI

## -----------------------------------------------------------------------------
sigma_mass <- sqrt(sum(sawe_input$sigma_mass^2))

## ----echo = FALSE-------------------------------------------------------------
sigma_mass

## -----------------------------------------------------------------------------
sigma_cm <- function(m, sigma_v, sigma_m, v, c) {
  (m * sigma_v)^2 + (sigma_m * (v - c))^2
}
sigma_C <- c(
  sigma_Cx = sqrt(sum(sigma_cm(sawe_input$mass, sawe_input$sigma_Cx, sawe_input$sigma_mass, sawe_input$Cx, C["Cx"]))) / mass,
  sigma_Cy = sqrt(sum(sigma_cm(sawe_input$mass, sawe_input$sigma_Cy, sawe_input$sigma_mass, sawe_input$Cy, C["Cy"]))) / mass,
  sigma_Cz = sqrt(sum(sigma_cm(sawe_input$mass, sawe_input$sigma_Cz, sawe_input$sigma_mass, sawe_input$Cz, C["Cz"]))) / mass
)

## ----echo = FALSE-------------------------------------------------------------
sigma_C

## -----------------------------------------------------------------------------
sigma_moi <- function(sigma_I, mass, sigma_mass, v1, v2, c1, c2, sigma_v1, sigma_v2) {
  sigma_I^2 +
    (2 * mass * (v1 - c1) * sigma_v1)^2 +
    (2 * mass * (v2 - c2) * sigma_v2)^2 +
    (((v1 - c1)^2 + (v2 - c2)^2) * sigma_mass)^2
}
sigma_MOI <- c(
  sigma_Ixx = sqrt(sum(sigma_moi(sawe_input$sigma_Ixx, sawe_input$mass, sawe_input$sigma_mass, sawe_input$Cy,
                                  sawe_input$Cz, C["Cy"], C["Cz"], sawe_input$sigma_Cy, sawe_input$sigma_Cz))),
  sigma_Iyy = sqrt(sum(sigma_moi(sawe_input$sigma_Iyy, sawe_input$mass, sawe_input$sigma_mass, sawe_input$Cx,
                                  sawe_input$Cz, C["Cx"], C["Cz"], sawe_input$sigma_Cx, sawe_input$sigma_Cz))),
  sigma_Izz = sqrt(sum(sigma_moi(sawe_input$sigma_Izz, sawe_input$mass, sawe_input$sigma_mass, sawe_input$Cx,
                                  sawe_input$Cy, C["Cx"], C["Cy"], sawe_input$sigma_Cx, sawe_input$sigma_Cy)))
)

## ----echo = FALSE-------------------------------------------------------------
sigma_MOI

## -----------------------------------------------------------------------------
sigma_poi <- function(sigma_I, mass, sigma_mass, v1, v2, c1, c2, sigma_v1, sigma_v2) {
  sigma_I^2 +
    ((v1 - c1) * mass * sigma_v2)^2 +
    ((v1 - c1) * (v2 - c2) * sigma_mass)^2 +
    ((v2 - c2) * mass * sigma_v1)^2
}
sigma_POI <- c(
  sigma_Ixy = sqrt(sum(sigma_poi(sawe_input$sigma_Ixy, sawe_input$mass, sawe_input$sigma_mass, sawe_input$Cx,
                                  sawe_input$Cy, C["Cx"], C["Cy"], sawe_input$sigma_Cx, sawe_input$sigma_Cy))),
  sigma_Ixz = sqrt(sum(sigma_poi(sawe_input$sigma_Ixz, sawe_input$mass, sawe_input$sigma_mass, sawe_input$Cx,
                                  sawe_input$Cz, C["Cx"], C["Cz"], sawe_input$sigma_Cx, sawe_input$sigma_Cz))),
  sigma_Iyz = sqrt(sum(sigma_poi(sawe_input$sigma_Iyz, sawe_input$mass, sawe_input$sigma_mass, sawe_input$Cy,
                                  sawe_input$Cz, C["Cy"], C["Cz"], sawe_input$sigma_Cy, sawe_input$sigma_Cz)))
)

## ----echo = FALSE-------------------------------------------------------------
sigma_POI

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

