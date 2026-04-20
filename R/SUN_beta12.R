#' SUN beta 1.2 — plateau-cache intervals as official intervals
#'
#' This wrapper sources beta 1.1 and overrides `SUN()` so that:
#' - plateau cache-aligned intervals replace `@tools[[tools_key]]$intervals`
#' - recommendation / top-k / silhouette interval set follows replaced intervals
#' - cache reuse and fallback behavior are preserved

# Beta 1.1 helpers are loaded from package file `R/SUN_beta11.R`.

sun_beta12_apply_plateau_intervals <- function(obj, tools_key = "all_interval_analysis", graph.name = NULL) {
  tk <- tools_key
  obj <- run_plateau_analysis_beta11(
    obj,
    tools_key = tk,
    slot_key = "plateau_intervals_beta11",
    r_min = NULL,
    r_max = NULL,
    drop_n_le_1 = TRUE,
    graph.name = graph.name,
    validate = TRUE
  )

  pl2 <- sun_plateau_intervals_cache_aligned(obj, tools_key = tk)
  iv_classic <- obj@tools[[tk]]$intervals

  obj@tools[[tk]]$intervals_classic_legacy <- iv_classic
  if (is.data.frame(pl2) && nrow(pl2) > 0L) {
    obj@tools[[tk]]$intervals <- pl2
    obj@tools[[tk]]$intervals_definition <- "plateau_cache_aligned_r_left"
  } else {
    obj@tools[[tk]]$intervals <- iv_classic
    obj@tools[[tk]]$intervals_definition <- "classic_stable_intervals_fallback"
  }
  obj
}

sun_beta12_write_label_from_cache <- function(obj, recommended_resolution, tools_key = "all_interval_analysis") {
  cache <- obj@tools[[tools_key]]$analysis_obj$cache
  ck <- as.character(as.integer(round(as.numeric(recommended_resolution) * 100)))
  if (!ck %in% names(cache)) stop("SUN beta 1.2: cache key missing for rec resolution: ", ck, call. = FALSE)
  ids <- cache[[ck]]$ids
  if (is.null(ids)) stop("SUN beta 1.2: cache ids missing for key ", ck, call. = FALSE)
  labs <- align_ids_to_cells(ids, colnames(obj))
  obj@meta.data[["SUN_Label"]] <- as.character(labs)
  obj@meta.data[["SUN_Label_beta11"]] <- as.character(labs)
  obj
}

#' One-click SUN workflow (beta 1.2 official).
#'
#' Replaces `@tools[[tools_key]]$intervals` with plateau cache-aligned intervals.
#' All downstream recommendation objects are based on the replaced interval set.
SUN <- function(seurat_obj, r_min = 0, r_max = 1, ...) {
  if (!inherits(seurat_obj, "Seurat")) stop("seurat_obj must be a Seurat object.", call. = FALSE)

  dots <- list(...)
  tk <- if ("tools_key" %in% names(dots)) dots$tools_key else "all_interval_analysis"
  gn <- if ("graph.name" %in% names(dots)) dots$graph.name else NULL

  if (!"initial_cache" %in% names(dots)) {
    slot_obj <- seurat_obj@tools[[tk]]
    if (!is.null(slot_obj) && is.list(slot_obj$analysis_obj) && length(slot_obj$analysis_obj$cache) > 0L) {
      cand <- slot_obj$analysis_obj$cache
      ok <- tryCatch(validate_sun_initial_cache(seurat_obj, cand, graph.name = gn), error = function(e) e)
      if (!inherits(ok, "error")) dots$initial_cache <- cand
    }
  }

  obj <- tryCatch(
    do.call(run_all_interval_analysis, c(list(seurat_obj = seurat_obj, r_min = r_min, r_max = r_max), dots)),
    error = function(e) {
      if (is.null(dots$initial_cache)) stop(e)
      dots_retry <- dots
      dots_retry$initial_cache <- NULL
      do.call(run_all_interval_analysis, c(list(seurat_obj = seurat_obj, r_min = r_min, r_max = r_max), dots_retry))
    }
  )

  obj <- tryCatch(
    sun_beta12_apply_plateau_intervals(obj, tools_key = tk, graph.name = gn),
    error = function(e) {
      warning(
        "SUN beta 1.2: plateau replacement failed, using classic intervals. ",
        conditionMessage(e),
        call. = FALSE
      )
      obj
    }
  )

  iv <- obj@tools[[tk]]$intervals
  rec <- recommend_resolution_topk_silhouette(
    seurat_obj = obj,
    tools_key = tk,
    intervals_df = iv
  )

  rr <- rec$recommended_resolution
  if (length(rr) != 1L || !is.finite(rr)) {
    if (is.data.frame(iv) && nrow(iv) > 0L) {
      j <- which.max(iv$interval_width %||% (iv$r_right - iv$r_left))[1L]
      rec$recommended_resolution <- as.numeric(iv$r_left[j])
      rec$recommendation_basis <- "fallback_widest_interval_r_left_beta12"
    } else {
      stop("SUN beta 1.2: no finite recommended resolution.", call. = FALSE)
    }
  } else if (identical(obj@tools[[tk]]$intervals_definition, "plateau_cache_aligned_r_left")) {
    rec$recommendation_basis <- "plateau_cache_intervals_topk_silhouette_beta12"
  }

  obj <- sun_beta12_write_label_from_cache(obj, rec$recommended_resolution, tools_key = tk)

  obj@tools[["SUN"]] <- list(
    recommended_resolution = rec$recommended_resolution,
    recommendation_basis = rec$recommendation_basis,
    top_k_intervals = rec$top_k_intervals,
    candidates = rec$candidates,
    intervals_definition = obj@tools[[tk]]$intervals_definition,
    intervals_replaced = TRUE
  )

  list(
    seurat_obj = obj,
    recommended_resolution = rec$recommended_resolution,
    recommendation = rec
  )
}
