test_that("plot_peak_network() labels a node by the requested family, not the globally strongest one", {
  skip_if_not_installed("visNetwork")

  # Feature 1 belongs to two adduct families: family 1, where its role has the
  # higher mean_edge_score (the one a global, family_id-agnostic pick would
  # choose), and family 2, where its role is different and its score lower.
  # Requesting family_id = 2 should label feature 1 with family 2's role.
  family_members <- data.frame(
    family_id = c(1L, 1L, 2L, 2L),
    idx = c(1L, 2L, 1L, 3L),
    mz = c(100.000, 138.963, 100.000, 144.961),
    adduct = c("[M+H]+", "[M+K]+", "[M+Na]+", "[M+2Na-H]+"),
    neutral_mass_from_feature = c(99.000, 99.000, 99.000, 99.000),
    mean_edge_score = c(0.9, 0.9, 0.3, 0.3),
    stringsAsFactors = FALSE
  )

  family_edges <- data.frame(
    family_id = c(1L, 2L),
    idx_i = c(1L, 1L),
    idx_j = c(2L, 3L),
    mz_i = c(100.000, 100.000),
    mz_j = c(138.963, 144.961),
    adduct_i = c("[M+H]+", "[M+Na]+"),
    adduct_j = c("[M+K]+", "[M+2Na-H]+"),
    score_adduct = c(0.9, 0.3),
    is_valid_adduct = c(TRUE, TRUE),
    neutral_mass_mean = c(99.000, 99.000),
    neutral_err_ppm = c(1, 1),
    stringsAsFactors = FALSE
  )

  wf <- list(adduct_families = list(
    family_edges = family_edges,
    family_members = family_members
  ))
  class(wf) <- "peakguider_workflow"

  node_role <- function(net, id) {
    nodes <- net$x$nodes
    list(
      label = nodes$label[nodes$id == id],
      group = nodes$group[nodes$id == id]
    )
  }

  # Requesting family 2 explicitly: feature 1 must be labelled with its
  # family-2 role ([M+Na]+), matching what plot_adduct_family(family_id = 2)
  # would show for the same feature - not its family-1 role, even though
  # family 1 has the higher mean_edge_score.
  net_family2 <- plot_peak_network(
    relation_table = wf, family_id = 2, min_score = 0.1, show_edge_labels = FALSE
  )
  role2 <- node_role(net_family2, "1")
  expect_match(role2$label, "[M+Na]+", fixed = TRUE)
  expect_false(grepl("[M+H]+", role2$label, fixed = TRUE))
  expect_identical(role2$group, "family_2")

  # Requesting family 1 explicitly: feature 1 must be labelled with its
  # family-1 role.
  net_family1 <- plot_peak_network(
    relation_table = wf, family_id = 1, min_score = 0.1, show_edge_labels = FALSE
  )
  role1 <- node_role(net_family1, "1")
  expect_match(role1$label, "[M+H]+", fixed = TRUE)
  expect_identical(role1$group, "family_1")

  # No family_id (full network): unchanged pre-existing behaviour - the
  # globally strongest family (by mean_edge_score) is used.
  net_full <- plot_peak_network(
    relation_table = wf, min_score = 0.1, show_edge_labels = FALSE
  )
  role_full <- node_role(net_full, "1")
  expect_match(role_full$label, "[M+H]+", fixed = TRUE)
  expect_identical(role_full$group, "family_1")
})
