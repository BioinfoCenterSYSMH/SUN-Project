#!/usr/bin/env Rscript
#' SUN beta 1.3 — **standalone** (`SUN_beta13.R`)
#'
#' Same behavior as sourcing `SUN_beta1.1.R` then `SUN_beta_1.3.R`: plateau intervals,
#' `SUN()` uses **bluster::approxSilhouette** on an embedding (not graph-distance silhouette).
#' This file is self-contained: **no `source()` of other R scripts**.
#'
#' @param silhouette_reduction,silhouette_dims Passed through `SUN(...)`; see `SUN_beta_1.3.R` docs.
#'
#' Requires: Seurat, Matrix, igraph, **Bioconductor bluster** (`BiocManager::install("bluster")`).

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(igraph)
})
if (!requireNamespace("bluster", quietly = TRUE)) {
  stop(
    "SUN beta 1.3 needs package 'bluster'. Install with:\n",
    "  BiocManager::install(\"bluster\")",
    call. = FALSE
  )
}
suppressPackageStartupMessages({
  library(bluster)
})

#' @keywords internal
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' @keywords internal
round_res2 <- function(x) round(as.numeric(x), 2L)

#' @keywords internal
to_centi <- function(x) as.integer(round(as.numeric(x) * 100))

#' @keywords internal
from_centi <- function(i) as.numeric(as.integer(i)) / 100

#' @keywords internal
centi_grid <- function(cmin, cmax, step_centi) {
  g <- seq.int(cmin, cmax, by = step_centi)
  if (tail(g, 1L) != cmax) g <- c(g, cmax)
  unique(sort(as.integer(g)))
}

#' Get Seurat SNN graph
#'
#' @param object A Seurat object.
#' @param graph.name Optional graph name. Defaults to `<DefaultAssay>_snn`.
#'
#' @return A Seurat graph object.
#' @keywords internal
get_snn_graph <- function(object, graph.name = NULL) {
  graph.name <- graph.name %||% paste0(Seurat::DefaultAssay(object), "_snn")
  if (!graph.name %in% names(object@graphs)) {
    stop("Graph not found: ", graph.name, call. = FALSE)
  }
  g <- object@graphs[[graph.name]]
  if (is.null(g)) stop("Graph is NULL: ", graph.name, call. = FALSE)
  g
}

#' @keywords internal
compute_clusters_on_snn <- function(
  snn,
  resolution,
  modularity.fxn = 1L,
  algorithm = 1L,
  n.start = 10L,
  n.iter = 10L,
  random.seed = 42L,
  group.singletons = TRUE,
  verbose = FALSE
) {
  ids <- Seurat:::RunModularityClustering(
    SNN = snn,
    modularity = modularity.fxn,
    resolution = round_res2(resolution),
    algorithm = algorithm,
    n.start = n.start,
    n.iter = n.iter,
    random.seed = random.seed,
    print.output = verbose,
    temp.file.location = NULL,
    edge.file.name = NULL
  )
  cn <- colnames(snn)
  if (!is.null(cn)) names(ids) <- cn
  ids <- Seurat:::GroupSingletons(
    ids = ids,
    SNN = snn,
    group.singletons = group.singletons,
    verbose = verbose
  )
  ids <- as.factor(ids)
  list(ids = ids, n_clusters = length(unique(ids)))
}

#' @keywords internal
.filter_intervals_left_tail_n_match <- function(intervals_df, tol, get_cached) {
  if (nrow(intervals_df) == 0L) return(intervals_df)
  tol_r <- round_res2(tol)

  n_tail_at_right_c <- function(right_c) {
    r_t <- round_res2(from_centi(right_c) - tol_r)
    as.integer(get_cached(to_centi(r_t))$n_clusters)
  }

  min_right_c_for_width <- function(rL, hi_c) {
    c0 <- to_centi(round_res2(as.numeric(rL) + tol_r + 0.01))
    c0 <- max(c0, to_centi(rL) + 1L)
    c <- c0
    while (c <= hi_c) {
      rRv <- from_centi(c)
      if (rRv - as.numeric(rL) > tol_r) {
        rt <- round_res2(rRv - tol_r)
        if (rt > as.numeric(rL)) return(c)
      }
      c <- c + 1L
    }
    NA_integer_
  }

  keep <- rep(FALSE, nrow(intervals_df))
  out_r_right <- rep(NA_real_, nrow(intervals_df))

  for (i in seq_len(nrow(intervals_df))) {
    rL <- as.numeric(intervals_df$r_left[i])
    rR <- as.numeric(intervals_df$r_right[i])
    nL <- as.integer(intervals_df$n_clusters[i])
    hi_c <- to_centi(rR)
    lo_c <- min_right_c_for_width(rL, hi_c)
    if (is.na(lo_c) || lo_c > hi_c) next

    ans <- NA_integer_
    for (m in seq.int(lo_c, hi_c)) {
      if (as.integer(n_tail_at_right_c(m)) == nL) ans <- m
    }
    if (!is.na(ans)) {
      keep[i] <- TRUE
      out_r_right[i] <- round_res2(from_centi(ans))
    }
  }

  out <- intervals_df[keep, , drop = FALSE]
  out$r_right <- out_r_right[keep]
  rownames(out) <- NULL
  out
}

#' @keywords internal
#' Map cache list names to centi keys used by \\code{get_cached}.
#' Avoids false cache misses when RDS stores names like \\code{"0.35"} vs \\code{"35"}.
normalize_sun_cache_key <- function(k) {
  k <- as.character(k)[1L]
  if (!nzchar(k)) return(k)
  ic <- suppressWarnings(as.integer(k))
  if (!is.na(ic)) return(as.character(ic))
  x <- suppressWarnings(as.numeric(k))
  if (is.finite(x) && x == floor(x) && x >= 0 && x <= 100000) {
    return(as.character(as.integer(x)))
  }
  if (is.finite(x) && x >= 0 && x <= 1) {
    return(as.character(as.integer(round(x * 100))))
  }
  k
}

#' @keywords internal
same_sun_cache_entry <- function(a, b) {
  if (is.null(a) || is.null(b)) return(FALSE)
  if (!identical(as.integer(a$n_clusters), as.integer(b$n_clusters))) return(FALSE)
  ia <- a$ids
  ib <- b$ids
  if (is.null(names(ia)) && is.null(names(ib))) {
    return(identical(as.character(ia), as.character(ib)))
  }
  if (is.null(names(ia)) || is.null(names(ib))) return(FALSE)
  oa <- order(names(ia))
  ob <- order(names(ib))
  identical(as.character(ia[oa]), as.character(ib[ob]))
}

#' @keywords internal
.unwrap_first_seurat <- function(obj, max_depth = 20L) {
  if (inherits(obj, "Seurat")) return(obj)
  if (!is.list(obj) || max_depth <= 0L) return(NULL)
  for (i in seq_along(obj)) {
    u <- .unwrap_first_seurat(obj[[i]], max_depth - 1L)
    if (!is.null(u)) return(u)
  }
  NULL
}

#' Extract SUN analysis cache from a Seurat or nested SUN container
#'
#' Walks nested lists (e.g. `fb_normal$seurat_obj`) to find the first `Seurat`,
#' then reads `seurat_obj@tools[[tools_key]]$analysis_obj$cache`.
#'
#' @param obj A `Seurat` or a list wrapping one (SUN result objects).
#' @param tools_key Tools slot name, default `"all_interval_analysis"`.
#' @return List: `seurat_obj`, `cache` (named list), `tools_key`.
get_sun_cache_from_object <- function(obj, tools_key = "all_interval_analysis") {
  seu <- .unwrap_first_seurat(obj)
  if (is.null(seu) || !inherits(seu, "Seurat")) {
    stop("Could not find a Seurat inside obj (max nesting 20).", call. = FALSE)
  }
  ta <- seu@tools[[tools_key]]
  if (is.null(ta) || !is.list(ta)) {
    stop("tools[[\"", tools_key, "\"]] missing or not a list.", call. = FALSE)
  }
  ca <- ta[["analysis_obj"]][["cache"]]
  if (is.null(ca) || length(ca) == 0L) {
    stop("No analysis_obj$cache under tools[[\"", tools_key, "\"]].", call. = FALSE)
  }
  list(seurat_obj = seu, cache = ca, tools_key = tools_key)
}

#' Plateau stable intervals from a SUN cache (RLE on sorted keys)
#'
#' **Definition (beta 1.1):** maximal runs where `n_clusters` is constant as
#' resolution increases along **cached** centi keys only. `r_left` / `r_right`
#' are the smallest / largest resolution in each run among those keys.
#'
#' @param cache Named list as in `analysis_obj$cache` (centi-like names).
#' @param r_min,r_max Optional bounds (inclusive) on resolution after rounding.
#' @param drop_n_le_1 Drop runs with `n_clusters <= 1`.
#' @return `data.frame` with `r_left`, `r_right`, `n_clusters`, `interval_width`.
plateau_intervals_from_cache <- function(
    cache,
    r_min = NULL,
    r_max = NULL,
    drop_n_le_1 = TRUE) {
  if (is.null(cache) || length(cache) == 0L) {
    return(data.frame(
      r_left = numeric(0), r_right = numeric(0), n_clusters = integer(0),
      interval_width = numeric(0)
    ))
  }
  keys <- names(cache)
  if (is.null(keys) || !length(keys)) {
    stop("cache must be a named list.", call. = FALSE)
  }
  centi <- suppressWarnings(as.integer(vapply(keys, normalize_sun_cache_key, "")))
  if (anyNA(centi) || any(!is.finite(centi))) {
    stop("Could not parse all cache names as centi integers.", call. = FALSE)
  }
  nvec <- vapply(seq_along(keys), function(i) {
    as.integer(cache[[keys[i]]][["n_clusters"]])
  }, integer(1L))
  o <- order(centi, nvec)
  centi <- centi[o]
  nvec <- nvec[o]
  dupk <- unique(centi[duplicated(centi)])
  if (length(dupk)) {
    for (k in dupk) {
      nv <- unique(nvec[centi == k])
      if (length(nv) > 1L) {
        warning(
          "Duplicate cache key ", k, " with conflicting n_clusters; keeping first occurrence.",
          call. = FALSE
        )
      }
    }
  }
  keep_first <- !duplicated(centi)
  centi <- centi[keep_first]
  nvec <- nvec[keep_first]
  if (!is.null(r_min)) {
    c0 <- to_centi(round_res2(r_min))
    keep <- centi >= c0
    centi <- centi[keep]
    nvec <- nvec[keep]
  }
  if (!is.null(r_max)) {
    c1 <- to_centi(round_res2(r_max))
    keep <- centi <= c1
    centi <- centi[keep]
    nvec <- nvec[keep]
  }
  if (length(centi) == 0L) {
    return(data.frame(
      r_left = numeric(0), r_right = numeric(0), n_clusters = integer(0),
      interval_width = numeric(0)
    ))
  }
  ch <- c(TRUE, nvec[-1L] != nvec[-length(nvec)])
  br <- c(which(ch), length(nvec) + 1L)
  out <- data.frame(
    r_left = numeric(0), r_right = numeric(0), n_clusters = integer(0)
  )
  for (j in seq_len(length(br) - 1L)) {
    s <- br[j]
    e <- br[j + 1L] - 1L
    n0 <- as.integer(nvec[s])
    if (isTRUE(drop_n_le_1) && n0 <= 1L) next
    rl <- round_res2(from_centi(centi[s]))
    rr <- round_res2(from_centi(centi[e]))
    out <- rbind(out, data.frame(r_left = rl, r_right = rr, n_clusters = n0))
  }
  rownames(out) <- NULL
  if (nrow(out) > 0L) out$interval_width <- out$r_right - out$r_left
  else out$interval_width <- numeric(0)
  out
}

#' Store plateau intervals on a Seurat (beta 1.1)
#'
#' Reads cache from an existing SUN tools slot, computes plateau intervals,
#' writes `seurat_obj@tools[[slot_key]]`. **No Louvain** unless you pre-fill
#' cache elsewhere.
#'
#' @inheritParams get_sun_cache_from_object
#' @param slot_key Where to store results, default `"plateau_intervals_beta11"`.
#' @inheritParams plateau_intervals_from_cache
#' @param graph.name Passed to [validate_sun_initial_cache()] when `validate=TRUE`.
#' @param validate If `TRUE`, run [validate_sun_initial_cache()] on the cache.
#' @return Updated `Seurat`.
run_plateau_analysis_beta11 <- function(
    seurat_obj,
    tools_key = "all_interval_analysis",
    slot_key = "plateau_intervals_beta11",
    r_min = NULL,
    r_max = NULL,
    drop_n_le_1 = TRUE,
    graph.name = NULL,
    validate = TRUE) {
  ext <- get_sun_cache_from_object(seurat_obj, tools_key = tools_key)
  if (isTRUE(validate)) {
    validate_sun_initial_cache(ext$seurat_obj, ext$cache, graph.name = graph.name)
  }
  pl <- plateau_intervals_from_cache(
    ext$cache,
    r_min = r_min,
    r_max = r_max,
    drop_n_le_1 = drop_n_le_1
  )
  ext$seurat_obj@tools[[slot_key]] <- list(
    intervals = pl,
    meta = list(
      definition = "plateau: n_clusters constant along increasing cached resolutions",
      source_tools_key = tools_key,
      r_min = r_min,
      r_max = r_max,
      n_rows = nrow(pl)
    )
  )
  ext$seurat_obj
}

#' Discover stable resolution intervals on SNN graph
#'
#' @param snn A Seurat SNN graph.
#' @param r_min Minimum resolution.
#' @param r_max Maximum resolution.
#' @param modularity.fxn Seurat modularity function id.
#' @param algorithm Seurat algorithm id.
#' @param n.start Louvain starts.
#' @param n.iter Louvain iterations.
#' @param random.seed Random seed.
#' @param group.singletons Group singleton option.
#' @param coarse_by Coarse grid step.
#' @param refine_with_bisection Whether to refine jump boundaries.
#' @param tol Left-tail tolerance.
#' @param show_progress Verbose progress.
#' @param auto_halve_coarse_by Whether to halve coarse_by when empty.
#' @param coarse_by_min Minimal coarse_by when halving.
#' @param max_halving_iters Max halving iterations.
#' @param require_left_tail_n_match Require `n(r_left)==n(r_right-tol)`.
#' @param initial_cache Optional named cache list from previous run.
#'
#' @return A list with `intervals`, `cache`, and `meta`.
stable_resolution_intervals_louvain_snn <- function(
  snn,
  r_min = 0,
  r_max = 1,
  modularity.fxn = 1L,
  algorithm = 1L,
  n.start = 10L,
  n.iter = 10L,
  random.seed = 42L,
  group.singletons = TRUE,
  coarse_by = 0.05,
  refine_with_bisection = TRUE,
  tol = 0.01,
  show_progress = TRUE,
  auto_halve_coarse_by = TRUE,
  coarse_by_min = NULL,
  max_halving_iters = 30L,
  require_left_tail_n_match = TRUE,
  initial_cache = NULL
) {
  empty_intervals <- data.frame(r_left = numeric(0), r_right = numeric(0), n_clusters = integer(0))
  r_min <- round_res2(r_min); r_max <- round_res2(r_max); tol <- round_res2(tol); coarse_by <- round_res2(coarse_by)
  if (is.null(coarse_by_min)) coarse_by_min <- max(tol, 0.01)
  coarse_by_min <- round_res2(coarse_by_min)

  cache <- new.env(parent = emptyenv())
  n_cache_seeded <- 0L
  n_seed_skipped_redundant <- 0L
  n_seed_conflict_kept_first <- 0L
  if (!is.null(initial_cache)) {
    for (k in names(initial_cache)) {
      kk <- normalize_sun_cache_key(k)
      ent <- initial_cache[[k]]
      if (exists(kk, envir = cache, inherits = FALSE)) {
        old <- get(kk, envir = cache, inherits = FALSE)
        if (same_sun_cache_entry(old, ent)) {
          n_seed_skipped_redundant <- n_seed_skipped_redundant + 1L
        } else {
          n_seed_conflict_kept_first <- n_seed_conflict_kept_first + 1L
          warning(
            "SUN cache seed: key ", kk, " already occupied; keeping first entry, not overwriting with slot name ",
            k,
            call. = FALSE
          )
        }
        next
      }
      assign(kk, ent, envir = cache)
      n_cache_seeded <- n_cache_seeded + 1L
    }
  }

  n_louvain_new_calls <- 0L
  new_call_resolutions <- numeric(0)
  new_call_keys <- character(0)
  get_cached <- function(centi) {
    key <- as.character(as.integer(centi))
    if (exists(key, envir = cache, inherits = FALSE)) return(get(key, envir = cache, inherits = FALSE))
    n_louvain_new_calls <<- n_louvain_new_calls + 1L
    r_miss <- from_centi(centi)
    new_call_resolutions <<- c(new_call_resolutions, r_miss)
    new_call_keys <<- c(new_call_keys, key)
    obj <- compute_clusters_on_snn(
      snn = snn,
      resolution = r_miss,
      modularity.fxn = modularity.fxn,
      algorithm = algorithm,
      n.start = n.start,
      n.iter = n.iter,
      random.seed = random.seed,
      group.singletons = group.singletons,
      verbose = FALSE
    )
    assign(key, obj, envir = cache)
    obj
  }

  first_centi_where_count_exceeds <- function(lo_excl, hi_incl, n_floor, tol_centi) {
    lo <- as.integer(lo_excl); hi <- as.integer(hi_incl)
    while ((hi - lo) > tol_centi) {
      mid <- (lo + hi + 1L) %/% 2L
      n_mid <- get_cached(mid)$n_clusters
      if (n_mid > n_floor) hi <- mid else lo <- mid
    }
    hi
  }

  attempt <- 0L
  cur_by <- coarse_by
  final_intervals <- empty_intervals

  repeat {
    attempt <- attempt + 1L
    cmin <- to_centi(r_min); cmax <- to_centi(r_max)
    step_centi <- to_centi(cur_by)
    tol_centi <- max(1L, to_centi(tol))
    gr <- centi_grid(cmin, cmax, step_centi)
    n_vec <- vapply(gr, function(g) as.integer(get_cached(g)$n_clusters), integer(1))

    if (!refine_with_bisection) {
      breaks <- c(1L, which(n_vec[-1L] != n_vec[-length(n_vec)]) + 1L, length(n_vec) + 1L)
      out <- data.frame(r_left = numeric(0), r_right = numeric(0), n_clusters = integer(0))
      for (k in seq_len(length(breaks) - 1L)) {
        s <- breaks[k]; e <- breaks[k + 1L] - 1L
        out <- rbind(out, data.frame(r_left = from_centi(gr[s]), r_right = from_centi(gr[e]), n_clusters = n_vec[s]))
      }
    } else {
      jumps <- integer(0)
      for (i in seq_len(length(gr) - 1L)) {
        if (n_vec[i] < n_vec[i + 1L]) {
          jumps <- c(jumps, first_centi_where_count_exceeds(gr[i], gr[i + 1L], n_vec[i], tol_centi))
        }
      }
      boundaries <- sort(unique(c(cmin, jumps, cmax)))
      if (length(boundaries) <= 1L) out <- empty_intervals else {
        out <- data.frame(
          r_left = from_centi(boundaries[-length(boundaries)]),
          r_right = from_centi(boundaries[-1L]),
          n_clusters = integer(length(boundaries) - 1L)
        )
        for (k in seq_len(nrow(out))) out$n_clusters[k] <- as.integer(get_cached(boundaries[k])$n_clusters)
      }
    }

    out <- out[out$n_clusters > 1L, , drop = FALSE]
    out$r_left <- round_res2(out$r_left); out$r_right <- round_res2(out$r_right); rownames(out) <- NULL
    final_intervals <- out

    if (isTRUE(require_left_tail_n_match) && nrow(final_intervals) > 0L) {
      final_intervals <- .filter_intervals_left_tail_n_match(final_intervals, tol = tol, get_cached = get_cached)
    }
    if (nrow(final_intervals) > 0L) break
    if (!auto_halve_coarse_by || attempt >= max_halving_iters) break
    next_by <- round_res2(cur_by / 2)
    if (next_by < coarse_by_min || next_by < 0.01) break
    cur_by <- next_by
  }

  cache_keys <- ls(envir = cache, all.names = TRUE)
  cache_list <- setNames(vector("list", length(cache_keys)), cache_keys)
  for (k in seq_along(cache_keys)) cache_list[[k]] <- get(cache_keys[k], envir = cache, inherits = FALSE)
  dup_final_names <- any(duplicated(names(cache_list)))

  list(
    intervals = final_intervals,
    cache = cache_list,
    meta = list(
      r_min = r_min, r_max = r_max, coarse_by_final = cur_by, tol = tol,
      modularity.fxn = modularity.fxn, algorithm = algorithm,
      n.start = n.start, n.iter = n.iter, random.seed = random.seed,
      group.singletons = group.singletons, refine_with_bisection = refine_with_bisection,
      require_left_tail_n_match = require_left_tail_n_match,
      n_cache_keys_seeded = n_cache_seeded,
      n_seed_skipped_redundant = n_seed_skipped_redundant,
      n_seed_conflict_kept_first = n_seed_conflict_kept_first,
      n_louvain_new_calls = n_louvain_new_calls,
      louvain_new_call_resolutions = new_call_resolutions,
      louvain_new_call_keys = new_call_keys,
      n_cache_keys_after = length(cache_keys),
      cache_list_has_duplicate_names = dup_final_names
    )
  )
}

#' @keywords internal
#' Validate that cached Louvain results belong to the same cells / SNN as `seurat_obj`.
#' @return Invisibly `TRUE`, or `stop()` with an English message if incompatible.
validate_sun_initial_cache <- function(seurat_obj, cache, graph.name = NULL) {
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object.", call. = FALSE)
  }
  if (is.null(cache) || length(cache) == 0L) {
    return(invisible(TRUE))
  }

  cells <- colnames(seurat_obj)
  if (length(cells) == 0L) {
    stop("Seurat object has no cells.", call. = FALSE)
  }

  snn <- get_snn_graph(seurat_obj, graph.name = graph.name)
  cn_snn <- colnames(snn)
  if (!is.null(cn_snn) && length(cn_snn) > 0L) {
    if (!identical(length(cn_snn), length(cells))) {
      stop("SNN graph size does not match the number of cells in the Seurat object.", call. = FALSE)
    }
    if (!identical(sort(cn_snn), sort(cells))) {
      stop("SNN graph cell names do not match Seurat::colnames().", call. = FALSE)
    }
  }

  for (k in names(cache)) {
    ent <- cache[[k]]
    if (is.null(ent) || is.null(ent[["ids"]])) next
    ids <- ent$ids
    nm <- names(ids)
    if (!is.null(nm)) {
      if (anyNA(nm) || any(!nzchar(nm))) {
        stop("Cached clustering has missing or empty cell names.", call. = FALSE)
      }
      if (!all(nm %in% cells)) {
        stop("Cached clustering cell names are not all present in current colnames.", call. = FALSE)
      }
      if (length(unique(nm)) != length(nm)) {
        stop("Cached clustering contains duplicate cell names.", call. = FALSE)
      }
    } else {
      if (length(ids) != length(cells)) {
        stop("Unnamed cached clustering length does not match the number of cells.", call. = FALSE)
      }
    }
  }

  invisible(TRUE)
}

#' Run stable interval analysis and store results in Seurat tools slot
#'
#' @param seurat_obj A Seurat object with SNN graph already computed.
#' @param graph.name Optional graph name. Defaults to `<DefaultAssay>_snn`.
#' @param r_min Minimum resolution.
#' @param r_max Maximum resolution.
#' @param modularity.fxn Seurat modularity function id.
#' @param algorithm Seurat algorithm id.
#' @param n.start Louvain starts.
#' @param n.iter Louvain iterations.
#' @param random.seed Random seed.
#' @param group.singletons Group singleton option.
#' @param coarse_by Coarse grid step.
#' @param refine_with_bisection Whether to refine jump boundaries.
#' @param tol Left-tail tolerance.
#' @param show_progress Verbose progress.
#' @param auto_halve_coarse_by Whether to halve coarse_by when empty.
#' @param coarse_by_min Minimal coarse_by when halving.
#' @param max_halving_iters Max halving iterations.
#' @param require_left_tail_n_match Require `n(r_left)==n(r_right-tol)`.
#' @param initial_cache Optional named cache list from previous run.
#' @param tools_key Seurat tools key to write.
#'
#' @return Updated Seurat object.
run_all_interval_analysis <- function(
  seurat_obj,
  graph.name = NULL,
  r_min = 0,
  r_max = 1,
  modularity.fxn = 1L,
  algorithm = 1L,
  n.start = 10L,
  n.iter = 10L,
  random.seed = 42L,
  group.singletons = TRUE,
  coarse_by = 0.05,
  refine_with_bisection = TRUE,
  tol = 0.01,
  show_progress = TRUE,
  auto_halve_coarse_by = TRUE,
  coarse_by_min = NULL,
  max_halving_iters = 30L,
  require_left_tail_n_match = TRUE,
  initial_cache = NULL,
  tools_key = "all_interval_analysis"
) {
  if (!inherits(seurat_obj, "Seurat")) stop("seurat_obj must be a Seurat object.", call. = FALSE)
  snn <- get_snn_graph(seurat_obj, graph.name = graph.name)
  analysis <- stable_resolution_intervals_louvain_snn(
    snn = snn,
    r_min = r_min, r_max = r_max, modularity.fxn = modularity.fxn, algorithm = algorithm,
    n.start = n.start, n.iter = n.iter, random.seed = random.seed, group.singletons = group.singletons,
    coarse_by = coarse_by, refine_with_bisection = refine_with_bisection, tol = tol,
    show_progress = show_progress, auto_halve_coarse_by = auto_halve_coarse_by,
    coarse_by_min = coarse_by_min, max_halving_iters = max_halving_iters,
    require_left_tail_n_match = require_left_tail_n_match, initial_cache = initial_cache
  )
  intervals <- analysis$intervals
  if (nrow(intervals) > 0L) intervals$interval_width <- intervals$r_right - intervals$r_left
  analysis$intervals <- intervals
  seurat_obj@tools[[tools_key]] <- list(
    intervals = intervals,
    analysis_obj = analysis,
    params = list(
      graph.name = graph.name, r_min = r_min, r_max = r_max,
      modularity.fxn = modularity.fxn, algorithm = algorithm,
      n.start = n.start, n.iter = n.iter, random.seed = random.seed,
      group.singletons = group.singletons, coarse_by = coarse_by,
      refine_with_bisection = refine_with_bisection, tol = tol,
      require_left_tail_n_match = require_left_tail_n_match,
      initial_cache_provided = !is.null(initial_cache)
    )
  )
  seurat_obj
}

#' @keywords internal
align_ids_to_cells <- function(ids, cells_obj) {
  if (!is.null(names(ids))) {
    mi <- match(cells_obj, names(ids))
    as.character(ids[mi])
  } else {
    if (length(ids) != length(cells_obj)) {
      stop("Clustering vector has no names and length != length(cells_obj).", call. = FALSE)
    }
    as.character(ids)
  }
}

#' @keywords internal
.to_centi_r_for_key <- function(r) {
  as.character(as.integer(round(as.numeric(r) * 100)))
}

#' Plateau intervals aligned to analysis cache (internal).
#'
#' Subset of `plateau_intervals_beta11$intervals` whose `r_left` maps to a
#' cache key with non-`NULL` `ids` — same rows used for plateau-first
#' recommendation and Shiny silhouette plots.
#'
#' @param obj A `Seurat` object after [run_plateau_analysis_beta11()] (or [SUN()]).
#' @param tools_key Tools slot with `analysis_obj$cache` (default `all_interval_analysis`).
#' @return A `data.frame` with plateau columns, possibly zero rows.
#' @keywords internal
sun_plateau_intervals_cache_aligned <- function(obj, tools_key = "all_interval_analysis") {
  empty <- data.frame(
    r_left = numeric(0),
    r_right = numeric(0),
    n_clusters = integer(0),
    interval_width = numeric(0),
    stringsAsFactors = FALSE
  )
  pl <- obj@tools[["plateau_intervals_beta11"]][["intervals"]]
  ts <- obj@tools[[tools_key]]
  if (is.null(ts) || is.null(ts$analysis_obj) || is.null(ts$analysis_obj$cache)) {
    return(empty)
  }
  cache <- ts$analysis_obj$cache
  if (length(cache) == 0L || is.null(pl) || !is.data.frame(pl) || nrow(pl) == 0L) {
    return(empty)
  }
  keys <- vapply(seq_len(nrow(pl)), function(i) .to_centi_r_for_key(pl$r_left[i]), character(1L))
  ok_row <- keys %in% names(cache) & vapply(keys, function(k) !is.null(cache[[k]]$ids), logical(1L))
  pl[ok_row, , drop = FALSE]
}
#' @keywords internal
mean_silhouette_embedding_approx <- function(emb, pred_chr) {
  ok <- !is.na(pred_chr) & nzchar(as.character(pred_chr))
  if (sum(ok) < 3L) return(NA_real_)
  emb2 <- emb[ok, , drop = FALSE]
  lab2 <- as.character(pred_chr[ok])
  if (length(unique(lab2)) < 2L) return(NA_real_)
  out <- tryCatch(
    bluster::approxSilhouette(x = emb2, clusters = lab2),
    error = function(e) NA_real_
  )
  if (identical(out, NA_real_) || is.null(out$width)) return(NA_real_)
  mean(as.numeric(out$width), na.rm = TRUE)
}

#' `recommend_resolution_topk_silhouette` (beta 1.3): embedding-based silhouette.
#' Keeps output column name `silhouette_mean` for compatibility, but values are
#' **mean(bluster::approxSilhouette)** on the chosen reduction (not graph distance).
recommend_resolution_topk_silhouette <- function(
    seurat_obj,
    tools_key = "all_interval_analysis",
    k = 3L,
    max_cells_silhouette = 2000L,
    subsample_seed = 1L,
    intervals_df = NULL,
    silhouette_reduction = NULL,
    silhouette_dims = NULL) {
  if (is.null(silhouette_reduction) || length(silhouette_reduction) == 0L ||
    !nzchar(as.character(silhouette_reduction[[1L]]))) {
    silhouette_reduction <- Sys.getenv("SUN_SILHOUETTE_REDUCTION", "pca")
  }
  silhouette_reduction <- as.character(silhouette_reduction[[1L]])

  if (!inherits(seurat_obj, "Seurat")) stop("seurat_obj must inherit from Seurat.", call. = FALSE)
  k <- as.integer(k)
  if (k < 1L) stop("k must be >= 1.", call. = FALSE)

  tools_slot <- seurat_obj@tools[[tools_key]]
  if (is.null(tools_slot) || is.null(tools_slot$analysis_obj) || is.null(tools_slot$analysis_obj$cache)) {
    stop("Need cache at seurat_obj@tools[[tools_key]]$analysis_obj$cache.", call. = FALSE)
  }
  cache <- tools_slot$analysis_obj$cache
  iv <- intervals_df %||% tools_slot$intervals
  if (is.null(iv) || !is.data.frame(iv) || nrow(iv) == 0L) {
    stop("No stable intervals (intervals_df / tools$intervals empty).", call. = FALSE)
  }

  em <- tryCatch(
    Seurat::Embeddings(seurat_obj, reduction = silhouette_reduction),
    error = function(e) {
      stop(
        "Cannot read Embeddings(seurat_obj, reduction='", silhouette_reduction,
        "'). Build this reduction with the same settings used for FindNeighbors.",
        call. = FALSE
      )
    }
  )
  cells_obj <- colnames(seurat_obj)
  miss_em <- setdiff(cells_obj, rownames(em))
  if (length(miss_em) > 0L) {
    stop(
      length(miss_em), " cells in object missing from Embeddings rownames for reduction '",
      silhouette_reduction, "'.",
      call. = FALSE
    )
  }
  em <- em[cells_obj, , drop = FALSE]
  if (is.null(silhouette_dims)) {
    silhouette_dims <- seq_len(min(25L, ncol(em)))
  } else {
    silhouette_dims <- as.integer(silhouette_dims)
    silhouette_dims <- silhouette_dims[is.finite(silhouette_dims) &
      silhouette_dims >= 1L & silhouette_dims <= ncol(em)]
    if (!length(silhouette_dims)) {
      stop("silhouette_dims empty or out of range for ncol(embedding)=", ncol(em), call. = FALSE)
    }
  }
  em <- em[, silhouette_dims, drop = FALSE]

  n_sub <- min(length(cells_obj), as.integer(max_cells_silhouette))
  set.seed(as.integer(subsample_seed))
  cells_sub_names <- cells_obj[sample.int(length(cells_obj), n_sub, replace = FALSE)]
  cells_sub_names2 <- cells_sub_names[cells_sub_names %in% rownames(em)]
  if (length(cells_sub_names2) < 3L) stop("Too few sampled cells for silhouette.", call. = FALSE)
  em_sub <- em[cells_sub_names2, , drop = FALSE]

  to_centi_key <- function(r) as.character(as.integer(round(as.numeric(r) * 100)))
  top <- iv
  top$width <- top$r_right - top$r_left
  missing_keys <- character(0)
  sils <- rep(NA_real_, nrow(top))
  for (i in seq_len(nrow(top))) {
    key <- to_centi_key(top$r_left[i])
    if (!key %in% names(cache) || is.null(cache[[key]]$ids)) {
      missing_keys <- c(missing_keys, key)
      next
    }
    pred_chr <- align_ids_to_cells(cache[[key]]$ids, cells_sub_names2)
    sils[i] <- mean_silhouette_embedding_approx(em_sub, pred_chr)
  }
  top$silhouette_mean <- sils
  missing_keys <- unique(missing_keys)

  finite_mask <- is.finite(top$silhouette_mean) & is.finite(top$width)
  top$score <- NA_real_
  if (any(finite_mask)) {
    sil_vals <- top$silhouette_mean[finite_mask]
    width_vals <- top$width[finite_mask]
    sil_rng <- max(sil_vals) - min(sil_vals)
    w_rng <- max(width_vals) - min(width_vals)
    sil_norm <- if (sil_rng == 0) rep(0.5, length(sil_vals)) else (sil_vals - min(sil_vals)) / sil_rng
    w_norm <- if (w_rng == 0) rep(0.5, length(width_vals)) else (width_vals - min(width_vals)) / w_rng
    top$score[finite_mask] <- 0.5 * sil_norm + 0.5 * w_norm
    max_score <- max(top$score[finite_mask], na.rm = TRUE)
    cand <- top[is.finite(top$score) & top$score == max_score, , drop = FALSE]
    recommended_resolution <- min(cand$r_left)
  } else {
    cand <- top[FALSE, , drop = FALSE]
    recommended_resolution <- NA_real_
  }

  idx2 <- which(is.finite(top$score))
  ord <- order(-top$score[idx2], -top$silhouette_mean[idx2], -top$width[idx2], top$r_left[idx2])
  top_idx <- idx2[ord]
  top_idx <- top_idx[seq_len(min(k, length(top_idx)))]
  top_k_intervals <- top[top_idx, , drop = FALSE]

  list(
    recommended_resolution = recommended_resolution,
    recommendation_basis = "embedding_approx_bluster_topk_by_score_equal_sil_width_beta13",
    top_k_intervals = top_k_intervals,
    candidates = cand,
    missing_cache_keys = missing_keys
  )
}

sun_beta13_apply_plateau_intervals <- function(obj, tools_key = "all_interval_analysis", graph.name = NULL) {
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

sun_beta13_write_label_from_cache <- function(obj, recommended_resolution, tools_key = "all_interval_analysis") {
  cache <- obj@tools[[tools_key]]$analysis_obj$cache
  ck <- as.character(as.integer(round(as.numeric(recommended_resolution) * 100)))
  if (!ck %in% names(cache)) stop("SUN beta 1.3: cache key missing for rec resolution: ", ck, call. = FALSE)
  ids <- cache[[ck]]$ids
  if (is.null(ids)) stop("SUN beta 1.3: cache ids missing for key ", ck, call. = FALSE)
  labs <- align_ids_to_cells(ids, colnames(obj))
  obj@meta.data[["SUN_Label"]] <- as.character(labs)
  obj@meta.data[["SUN_Label_beta11"]] <- as.character(labs)
  obj
}

#' One-click SUN workflow (beta 1.3).
#'
#' Same plateau interval replacement as beta 1.2, but recommendation uses
#' **bluster::approxSilhouette** on `silhouette_reduction` / `silhouette_dims`.
SUN <- function(seurat_obj, r_min = 0, r_max = 1, ...) {
  if (!inherits(seurat_obj, "Seurat")) stop("seurat_obj must be a Seurat object.", call. = FALSE)

  dots <- list(...)
  tk <- if ("tools_key" %in% names(dots)) dots$tools_key else "all_interval_analysis"
  gn <- if ("graph.name" %in% names(dots)) dots$graph.name else NULL

  sil_red <- if ("silhouette_reduction" %in% names(dots)) dots$silhouette_reduction else NULL
  sil_dims <- if ("silhouette_dims" %in% names(dots)) dots$silhouette_dims else NULL
  dots_run <- dots[!names(dots) %in% c("silhouette_reduction", "silhouette_dims")]

  if (!"initial_cache" %in% names(dots_run)) {
    slot_obj <- seurat_obj@tools[[tk]]
    if (!is.null(slot_obj) && is.list(slot_obj$analysis_obj) && length(slot_obj$analysis_obj$cache) > 0L) {
      cand <- slot_obj$analysis_obj$cache
      ok <- tryCatch(validate_sun_initial_cache(seurat_obj, cand, graph.name = gn), error = function(e) e)
      if (!inherits(ok, "error")) dots_run$initial_cache <- cand
    }
  }

  obj <- tryCatch(
    do.call(run_all_interval_analysis, c(list(seurat_obj = seurat_obj, r_min = r_min, r_max = r_max), dots_run)),
    error = function(e) {
      if (is.null(dots_run$initial_cache)) stop(e)
      dots_retry <- dots_run
      dots_retry$initial_cache <- NULL
      do.call(run_all_interval_analysis, c(list(seurat_obj = seurat_obj, r_min = r_min, r_max = r_max), dots_retry))
    }
  )

  obj <- tryCatch(
    sun_beta13_apply_plateau_intervals(obj, tools_key = tk, graph.name = gn),
    error = function(e) {
      warning(
        "SUN beta 1.3: plateau replacement failed, using classic intervals. ",
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
    intervals_df = iv,
    silhouette_reduction = sil_red,
    silhouette_dims = sil_dims
  )

  rr <- rec$recommended_resolution
  if (length(rr) != 1L || !is.finite(rr)) {
    if (is.data.frame(iv) && nrow(iv) > 0L) {
      j <- which.max(iv$interval_width %||% (iv$r_right - iv$r_left))[1L]
      rec$recommended_resolution <- as.numeric(iv$r_left[j])
      rec$recommendation_basis <- "fallback_widest_interval_r_left_beta13"
    } else {
      stop("SUN beta 1.3: no finite recommended resolution.", call. = FALSE)
    }
  } else if (identical(obj@tools[[tk]]$intervals_definition, "plateau_cache_aligned_r_left")) {
    rec$recommendation_basis <- "plateau_cache_intervals_topk_silhouette_embedding_approx_beta13"
  }

  obj <- sun_beta13_write_label_from_cache(obj, rec$recommended_resolution, tools_key = tk)

  obj@tools[["SUN"]] <- list(
    recommended_resolution = rec$recommended_resolution,
    recommendation_basis = rec$recommendation_basis,
    top_k_intervals = rec$top_k_intervals,
    candidates = rec$candidates,
    intervals_definition = obj@tools[[tk]]$intervals_definition,
    intervals_replaced = TRUE,
    silhouette_mode = "bluster_approxSilhouette_on_reduction",
    silhouette_reduction = if (is.null(sil_red) || !nzchar(as.character(sil_red)[1])) {
      Sys.getenv("SUN_SILHOUETTE_REDUCTION", "pca")
    } else {
      as.character(sil_red)[1]
    },
    silhouette_dims = sil_dims
  )

  list(
    seurat_obj = obj,
    recommended_resolution = rec$recommended_resolution,
    recommendation = rec
  )
}

#' Alias for TOP3 recommendation
#'
#' @inheritParams recommend_resolution_topk_silhouette
#' @return Same as [recommend_resolution_topk_silhouette()].
recommend_resolution_top3_silhouette <- function(...) {
  recommend_resolution_topk_silhouette(..., k = 3L)
}

#' Optional alias: forwards to [recommend_resolution_topk_silhouette()] (embedding silhouette in beta 1.3) with default `max_cells_silhouette = 5000`.
recommend_resolution_topk_silhouette_graph_distance <- function(
    seurat_obj,
    tools_key = "all_interval_analysis",
    k = 3L,
    max_cells_silhouette = 5000L,
    subsample_seed = 1L,
    intervals_df = NULL,
    silhouette_reduction = NULL,
    silhouette_dims = NULL
) {
  recommend_resolution_topk_silhouette(
    seurat_obj = seurat_obj,
    tools_key = tools_key,
    k = k,
    max_cells_silhouette = max_cells_silhouette,
    subsample_seed = subsample_seed,
    intervals_df = intervals_df,
    silhouette_reduction = silhouette_reduction,
    silhouette_dims = silhouette_dims
  )
}

message(
  "SUN beta 1.3 loaded: interval score uses bluster::approxSilhouette on ",
  "Embeddings(reduction = SUN_SILHOUETTE_REDUCTION or silhouette_reduction arg)."
)
