#' SUN beta 1.1 — plateau-first recommendation + classic SUN workflow
#'
#' Self-contained copy of `近似图距离.R` with **SUN beta 1.1** additions:
#' - **Plateau intervals** (definition: resolution increases while cluster count
#'   is constant on **cached** Louvain results): [r_left, r_right] from sorted
#'   cache keys via run-length encoding of `n_clusters`. No new Louvain when
#'   only using existing cache.
#' - Helpers: [get_sun_cache_from_object()], [plateau_intervals_from_cache()],
#'   [run_plateau_analysis_beta11()].
#' - **Official recommendation (`SUN()`):** after [run_all_interval_analysis()],
#'   runs [run_plateau_analysis_beta11()], then [recommend_resolution_topk_silhouette()]
#'   with **`intervals_df` = plateau intervals** (same strategy as
#'   `run_eval_SUN_beta11_vs_label.R`). Falls back to classic stable intervals, then
#'   widest plateau, if needed. Writes **`meta.data$SUN_Label`** and **`SUN_Label_beta11`**
#'   from the cache at the chosen resolution.
#' - Unchanged: classic jump-based `stable_resolution_intervals_louvain_snn` for
#'   building `all_interval_analysis` cache and intervals.
#'
#' @import Seurat
#' @import Matrix
#' @import igraph

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

#' Plateau-first resolution recommendation + label columns (internal).
#'
#' Mirrors `run_eval_SUN_beta11_vs_label.R`: plateau intervals →
#' `recommend_resolution_topk_silhouette(..., intervals_df = pl2)` →
#' `SUN_Label` / `SUN_Label_beta11`.
#'
#' @keywords internal
sun_recommend_resolution_plateau_then_classic <- function(obj, tools_key, graph.name) {
  tk <- tools_key
  cache <- obj@tools[[tk]]$analysis_obj$cache
  if (is.null(cache) || length(cache) == 0L) {
    stop("SUN: empty analysis cache after run_all_interval_analysis.", call. = FALSE)
  }

  obj <- tryCatch(
    run_plateau_analysis_beta11(
      obj,
      tools_key = tk,
      slot_key = "plateau_intervals_beta11",
      r_min = NULL,
      r_max = NULL,
      drop_n_le_1 = TRUE,
      graph.name = graph.name,
      validate = TRUE
    ),
    error = function(e) {
      warning(
        "SUN: run_plateau_analysis_beta11 failed; recommendation will use classic stable intervals. ",
        conditionMessage(e),
        call. = FALSE
      )
      obj
    }
  )

  pl2 <- sun_plateau_intervals_cache_aligned(obj, tools_key = tk)

  rec <- NULL
  if (nrow(pl2) > 0L) {
    rec_pl <- recommend_resolution_topk_silhouette(
      obj,
      tools_key = tk,
      k = 3L,
      max_cells_silhouette = 2000L,
      subsample_seed = 1L,
      intervals_df = pl2
    )
    rr <- rec_pl$recommended_resolution
    if (length(rr) == 1L && is.finite(rr)) {
      rec <- rec_pl
      rec$recommendation_basis <- "plateau_intervals_topk_silhouette"
    }
  }

  if (is.null(rec)) {
    rec <- recommend_resolution_topk_silhouette(seurat_obj = obj, tools_key = tk)
  }

  rr <- rec$recommended_resolution
  if (length(rr) != 1L || !is.finite(rr)) {
    if (nrow(pl2) > 0L) {
      best_i <- which.max(pl2$interval_width)[1L]
      rec_r <- as.numeric(pl2$r_left[best_i])
      ordw <- order(-pl2$interval_width, seq_len(nrow(pl2)))
      top_idx <- ordw[seq_len(min(3L, length(ordw)))]
      rec <- list(
        recommended_resolution = rec_r,
        recommendation_basis = "fallback_widest_plateau_r_left",
        top_k_intervals = pl2[top_idx, , drop = FALSE],
        candidates = pl2[best_i, , drop = FALSE],
        missing_cache_keys = character(0L)
      )
    } else {
      stop("SUN: could not determine a finite recommended_resolution.", call. = FALSE)
    }
  }

  ck <- .to_centi_r_for_key(rec$recommended_resolution)
  if (!ck %in% names(cache)) {
    stop("SUN: cache missing key for recommended resolution: ", ck, call. = FALSE)
  }
  ent <- cache[[ck]]
  if (is.null(ent$ids)) {
    stop("SUN: cache entry has no ids for key ", ck, call. = FALSE)
  }
  labs <- align_ids_to_cells(ent$ids, colnames(obj))
  obj@meta.data[["SUN_Label"]] <- as.character(labs)
  obj@meta.data[["SUN_Label_beta11"]] <- as.character(labs)

  list(obj = obj, rec = rec)
}

#' @keywords internal
get_graph_adjacency <- function(snn_graph) {
  if (inherits(snn_graph, "Matrix")) return(snn_graph)
  if (inherits(snn_graph, "Graph")) {
    sn <- try(methods::slotNames(snn_graph), silent = TRUE)
    if (!inherits(sn, "try-error") && "graph" %in% sn) {
      adj <- methods::slot(snn_graph, "graph")
      if (inherits(adj, "Matrix")) return(adj)
    }
  }
  adj <- try(Matrix::as(snn_graph, "dgCMatrix"), silent = TRUE)
  if (!inherits(adj, "try-error") && inherits(adj, "Matrix")) return(adj)
  stop("Unable to extract sparse adjacency matrix from graph object.", call. = FALSE)
}

#' @keywords internal
#' Non-finite graph distances (unreachable pairs) are imputed as \\code{max(finite D) + 1}
#' before averaging, so silhouette is well-defined on disconnected induced subgraphs.
mean_silhouette_graph_distance_from_D <- function(D, cl) {
  n <- length(cl)
  if (n < 3L) return(NA_real_)
  cl <- as.factor(cl)
  k <- length(levels(cl))
  if (k < 2L) return(NA_real_)
  D <- as.matrix(D)
  finite_vals <- D[is.finite(D)]
  if (length(finite_vals) == 0L) return(NA_real_)
  fill_val <- max(finite_vals, na.rm = TRUE) + 1
  D[!is.finite(D)] <- fill_val
  cl_idx <- as.integer(cl)
  cluster_sizes <- tabulate(cl_idx, nbins = k)
  mean_to_c <- matrix(NA_real_, nrow = n, ncol = k)
  for (c in seq_len(k)) {
    idx_c <- which(cl_idx == c)
    if (length(idx_c) == 0L) next
    v <- rowMeans(D[, idx_c, drop = FALSE], na.rm = TRUE)
    v[is.nan(v)] <- NA_real_
    mean_to_c[, c] <- v
  }
  a <- rep(NA_real_, n); b <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    cs <- cl_idx[i]; m_same <- cluster_sizes[cs]
    if (m_same > 1L) {
      ma <- mean_to_c[i, cs]
      if (is.finite(ma)) a[i] <- (ma * m_same) / (m_same - 1L)
    }
    vals <- mean_to_c[i, setdiff(seq_len(k), cs), drop = TRUE]
    vals <- vals[is.finite(vals)]
    if (length(vals) > 0L) b[i] <- min(vals)
  }
  denom <- pmax(a, b)
  s <- (b - a) / denom
  s[!is.finite(s)] <- NA_real_
  mean(s, na.rm = TRUE)
}

#' Recommend resolution by TOP-k graph-distance silhouette score
#'
#' `silhouette_mean` uses approximate graph distances on a sampled induced SNN subgraph.
#' Final score is parameter-free:
#' `score = 0.5 * sil_norm + 0.5 * width_norm`.
#'
#' @param seurat_obj A Seurat object that already has
#'   `@tools[[tools_key]]$analysis_obj$cache`.
#' @param tools_key Tools slot key. Default `"all_interval_analysis"`.
#' @param k Number of top intervals to return.
#' @param max_cells_silhouette Max sampled cells for graph-distance silhouette.
#' @param subsample_seed Sampling seed.
#' @param intervals_df Optional custom intervals data.frame.
#'
#' @return A list with `recommended_resolution`, `recommendation_basis`,
#'   `top_k_intervals`, `candidates`, and `missing_cache_keys`.
recommend_resolution_topk_silhouette <- function(
  seurat_obj,
  tools_key = "all_interval_analysis",
  k = 3L,
  max_cells_silhouette = 2000L,
  subsample_seed = 1L,
  intervals_df = NULL
) {
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

  snn <- get_snn_graph(seurat_obj, graph.name = NULL)
  A <- get_graph_adjacency(snn)
  cells_obj <- colnames(seurat_obj)
  n_sub <- min(length(cells_obj), as.integer(max_cells_silhouette))
  set.seed(as.integer(subsample_seed))
  cells_sub_names <- cells_obj[sample.int(length(cells_obj), n_sub, replace = FALSE)]

  rn <- rownames(A); cn <- colnames(A)
  if (is.null(rn) || is.null(cn)) {
    stop("SNN adjacency matrix has no row/col names; cannot align to cell names.", call. = FALSE)
  }
  keep <- cells_sub_names %in% rn & cells_sub_names %in% cn
  cells_sub_names2 <- cells_sub_names[keep]
  if (length(cells_sub_names2) < 3L) stop("Too few sampled cells for silhouette.", call. = FALSE)
  idx <- match(cells_sub_names2, rn)
  A_sub <- A[idx, idx, drop = FALSE]
  g_sub <- igraph::graph_from_adjacency_matrix(A_sub, mode = "undirected", weighted = NULL, diag = FALSE)
  D <- igraph::distances(g_sub, v = igraph::V(g_sub), to = igraph::V(g_sub), weights = NA)

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
    ok <- !is.na(pred_chr) & nzchar(pred_chr)
    if (sum(ok) < 3L) next
    sils[i] <- mean_silhouette_graph_distance_from_D(D[ok, ok, drop = FALSE], pred_chr[ok])
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
    recommendation_basis = "graph_distance_topk_by_score_equal_sil_width_r_left",
    top_k_intervals = top_k_intervals,
    candidates = cand,
    missing_cache_keys = missing_keys
  )
}

#' Alias for TOP3 recommendation
#'
#' @inheritParams recommend_resolution_topk_silhouette
#' @return Same as [recommend_resolution_topk_silhouette()].
recommend_resolution_top3_silhouette <- function(...) {
  recommend_resolution_topk_silhouette(..., k = 3L)
}

#' Optional alias: same as [recommend_resolution_topk_silhouette()] with default `max_cells_silhouette = 5000`.
recommend_resolution_topk_silhouette_graph_distance <- function(
    seurat_obj,
    tools_key = "all_interval_analysis",
    k = 3L,
    max_cells_silhouette = 5000L,
    subsample_seed = 1L,
    intervals_df = NULL
) {
  recommend_resolution_topk_silhouette(
    seurat_obj = seurat_obj,
    tools_key = tools_key,
    k = k,
    max_cells_silhouette = max_cells_silhouette,
    subsample_seed = subsample_seed,
    intervals_df = intervals_df
  )
}

#' One-click SUN workflow
#'
#' Requires a Seurat object with a precomputed SNN graph. Runs
#' [run_all_interval_analysis()], then **plateau-first** recommendation (same
#' strategy as `run_eval_SUN_beta11_vs_label.R`: [run_plateau_analysis_beta11()]
#' → [recommend_resolution_topk_silhouette()] with plateau `intervals_df`,
#' classic intervals as fallback, then widest plateau if needed). Writes
#' `meta.data$SUN_Label` and `SUN_Label_beta11` from the chosen cache resolution,
#' and stores a short summary in `obj@tools$SUN`.
#'
#' If `initial_cache` is not passed in `...` and a non-empty
#' `seurat_obj@tools[[tools_key]]$analysis_obj$cache` exists (e.g. from a prior
#' [run_all_interval_analysis()]), that cache is validated and passed as
#' `initial_cache` so [FindClusters]-equivalent work is skipped for hits inside
#' the resolution search. When validation fails, issues an English
#' \code{warning()} that the SUN cache does not match and runs without reusing it.
#' If analysis still errors while reusing cache, retries once without it and warns.
#'
#' @param seurat_obj A `Seurat` object.
#' @param r_min,r_max Resolution scan bounds passed to [run_all_interval_analysis()]
#'   (defaults `0` and `1`). Use e.g. `r_max = 5` to sweep a wider range (Shiny demo).
#' @param ... Optional arguments forwarded to [run_all_interval_analysis()].
#'   Pass `initial_cache = NULL` explicitly to disable auto-reuse from `@tools`.
#'
#' @return List with `seurat_obj`, `recommended_resolution`, and `recommendation`
#'   (side effect: `SUN_Label` / `SUN_Label_beta11` on `seurat_obj`).
SUN <- function(seurat_obj, r_min = 0, r_max = 1, ...) {
  if (!inherits(seurat_obj, "Seurat")) stop("seurat_obj must be a Seurat object.", call. = FALSE)

  dots <- list(...)
  tk <- if ("tools_key" %in% names(dots)) dots$tools_key else "all_interval_analysis"
  gn <- if ("graph.name" %in% names(dots)) dots$graph.name else NULL

  if (!"initial_cache" %in% names(dots)) {
    slot_obj <- seurat_obj@tools[[tk]]
    if (
      !is.null(slot_obj) && is.list(slot_obj$analysis_obj) &&
        length(slot_obj$analysis_obj$cache) > 0L
    ) {
      cand <- slot_obj$analysis_obj$cache
      ok <- tryCatch(
        validate_sun_initial_cache(seurat_obj, cand, graph.name = gn),
        error = function(e) e
      )
      if (inherits(ok, "error")) {
        warning(
          "SUN cache does not match this Seurat object; clusterings will be recomputed as needed. ",
          conditionMessage(ok),
          call. = FALSE
        )
      } else {
        dots$initial_cache <- cand
      }
    }
  }

  obj <- tryCatch(
    do.call(
      run_all_interval_analysis,
      c(list(seurat_obj = seurat_obj, r_min = r_min, r_max = r_max), dots)
    ),
    error = function(e) {
      if (is.null(dots$initial_cache)) {
        stop(e)
      }
      warning(
        "SUN cache does not match this Seurat object; retrying analysis without the stored cache.",
        call. = FALSE
      )
      dots_retry <- dots
      dots_retry$initial_cache <- NULL
      do.call(
        run_all_interval_analysis,
        c(list(seurat_obj = seurat_obj, r_min = r_min, r_max = r_max), dots_retry)
      )
    }
  )

  out_lab <- sun_recommend_resolution_plateau_then_classic(obj, tools_key = tk, graph.name = gn)
  obj <- out_lab$obj
  rec <- out_lab$rec

  obj@tools[["SUN"]] <- list(
    recommended_resolution = rec$recommended_resolution,
    recommendation_basis = rec$recommendation_basis,
    top_k_intervals = rec$top_k_intervals,
    candidates = rec$candidates
  )

  list(
    seurat_obj = obj,
    recommended_resolution = rec$recommended_resolution,
    recommendation = rec
  )
}

