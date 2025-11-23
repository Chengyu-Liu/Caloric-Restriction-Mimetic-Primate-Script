#Monkey scRNAseq age-model LOOCV -----by LCY
library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
library(future)
library(RColorBrewer)
library(patchwork)
library(ggalluvial)
library(ggridges)
library(cowplot)
library(reshape)
library(reshape2)
library(ggpubr)
library(ggrepel)
library(glmnet)
set.seed(5)
#----------------------------------------------------------------------------------------------#
plan("multicore", workers = 25)
options(future.globals.maxSize = 400000 * 1024^2)
#-------------------------------------#
ts <- "liver"

analysis_dir <- 'liuchengyu/'
sc_clk_dir <- paste0(analysis_dir,'snClock/')


out_root <- paste0(sc_clk_dir ,'/training/')
dir.create(out_root)

srt <- readRDS(paste0(analysis_dir,'liver_celltype.rds'))
# srt refers to the seurat object obtained after integrating samples and performing cell-type annotation.
# The metadata file must include the following fields: sample name (sample), cell type (cell_abbr), age (age_number), age group (age), and sex(gender).
srt@meta.data$bcd <- rownames(srt@meta.data)
srt@meta.data$bcd <- rownames(srt@meta.data)
o_donors <- unique(srt@meta.data[srt@meta.data$age == "O",]$bcd)
train_id <- sample(o_donors, length(o_donors) / 2)
test_id  <- setdiff(o_donors, train_id)

srt$O_split <- srt$age
srt$O_split <- as.vector(srt$O_split)
srt@meta.data[srt$bcd %in% train_id,]$O_split <- "O_train"
srt@meta.data[srt$bcd %in% test_id, ]$O_split <- "O_test"

combination <- srt

## ===================== 1) Creat meta cells =====================
vc1 <- SetupForWGCNA(
  combination,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.1, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "tutorial" # the name of the hdWGCNA experiment
)
vc1=MetacellsByGroups(
  vc1,
  group.by = c('sample','gender','age_number','age',"cell_abbr",'O_split'),
  ident.group = "cell_abbr",
  k = 25,
  reduction = "umap",
  assay = NULL,
  slot = "counts",
  mode = "average",
  wgcna_name = "tutorial",
  cells.use = NULL,
  min_cells = 50,
  max_shared = 15,
  target_metacells = 100,
  max_iter = 5000,
  verbose = FALSE
)

metadat <- GetMetacellObject(vc1)
metadat$age_number <- as.numeric(metadat$age_number)


combination <- metadat
combination <- NormalizeData(combination)
slt <- 6
## ===================== 2) Hyperparameters =====================
ALPHAS             <- c(0.1,0.3,0.5,0.7,0.9)
N_FOLDS_INNER      <- 7
TYPE_MEASURE       <- "mae"
USE_CALIBRATION    <- TRUE
TARGET_MIN_GENES   <- 20
MAX_ROUNDS         <- 50
FALLBACK_LMIN      <- TRUE

B_REPEATS          <- 10           
SUBSAMPLE_FRAC     <- 0.8          
RESAMPLE_REPLACE   <- FALSE        
TOP_NONZERO_FRAC   <- 0.90         
FREQ_THR_INIT      <- 0.65         
FREQ_THR_STEP      <- 0.05

USE_SAMPLE_WEIGHTS <- TRUE

GROUPED_INNER_CV   <- TRUE

MAX_POINTS_PLOT    <- 150000

ASSAY              <- "RNA"        
DATA_SLOT          <- "data"       
SAMPLE_COL         <- "sample"
GROUP_COL          <- "age"     
AGE_COL            <- "age_number"
CELLTYPE_COL       <- "cell_abbr"  

MIN_PCT_CELLS      <- 0.05         
MIN_VAR_GENE       <- 1e-10        


## ===================== 3) Utility functions =====================
to_x <- function(M) if (inherits(M, "dgCMatrix")) M else as.matrix(M)

# Split by sample to generate a foldid vector (with length equal to the number of cells).
make_foldid_by_sample <- function(sample_id, nfolds, seed=NULL) {
  if (!is.null(seed)) {
    old <- if (exists(".Random.seed", envir=.GlobalEnv)) get(".Random.seed", envir=.GlobalEnv) else NULL
    on.exit({ if (!is.null(old)) assign(".Random.seed", old, envir=.GlobalEnv) }, add=TRUE)
    set.seed(seed)
  }
  s_u <- unique(sample_id)
  ns  <- length(s_u)
  nfolds_eff <- max(2, min(nfolds, ns - 1))
  s_u <- sample(s_u)  # shuffle
  # Distribute samples across folds as evenly as possible.
  bins <- rep(1:nfolds_eff, length.out = ns)
  fold_map <- setNames(bins, s_u)
  as.integer(fold_map[sample_id])
}

inner_cv_fit <- function(X, y, sample_id=NULL, weights=NULL,
                         alphas=ALPHAS, nfolds=N_FOLDS_INNER,
                         type.measure=TYPE_MEASURE, seed=NULL,
                         grouped_by_sample=GROUPED_INNER_CV) {
  if (!is.null(seed)) {
    old <- if (exists(".Random.seed", envir=.GlobalEnv)) get(".Random.seed", envir=.GlobalEnv) else NULL
    on.exit({ if (!is.null(old)) assign(".Random.seed", old, envir=.GlobalEnv) }, add=TRUE)
    set.seed(seed)
  }
  n <- length(y)
  if (grouped_by_sample && !is.null(sample_id)) {
    foldid <- make_foldid_by_sample(sample_id, nfolds=nfolds, seed=seed)
    nfolds_eff <- length(unique(foldid))
  } else {
    nfolds_eff <- max(2, min(nfolds, n - 1))
    foldid <- sample(rep(1:nfolds_eff, length.out = n))
  }
  best <- list(score=Inf, alpha=NA, cvfit=NULL, foldid=foldid)
  for (a in alphas) {
    cvfit <- cv.glmnet(
      x=to_x(X), y=y, family="gaussian",
      alpha=a, nfolds=nfolds_eff, foldid=foldid,
      type.measure=type.measure, standardize=TRUE,
      nlambda=200, lambda.min.ratio=0.01,
      weights=weights
    )
    sc <- cvfit$cvm[match(cvfit$lambda.1se, cvfit$lambda)]
    if (sc < best$score) best <- list(score=sc, alpha=a, cvfit=cvfit, foldid=foldid)
  }
  best
}

coef_1se_or_min <- function(cvfit) {
  co_1se <- coef(cvfit, s="lambda.1se")
  beta_1se <- as.numeric(co_1se[-1, , drop=TRUE]); names(beta_1se) <- rownames(co_1se)[-1]
  if (sum(beta_1se != 0) == 0 && FALLBACK_LMIN) {
    co_min <- coef(cvfit, s="lambda.min")
    beta_min <- as.numeric(co_min[-1, , drop=TRUE]); names(beta_min) <- rownames(co_min)[-1]
    list(beta=beta_min, s="lambda.min", coefmat=co_min)
  } else list(beta=beta_1se, s="lambda.1se", coefmat=co_1se)
}

one_sd_importance <- function(beta_vec, X_train) {
  nz <- which(beta_vec != 0)
  if (length(nz) == 0) return(data.frame())
  feats <- names(beta_vec)[nz]
  sdX <- apply(as.matrix(X_train[, feats, drop=FALSE]), 2, sd)
  eff <- beta_vec[feats] * sdX
  out <- data.frame(feature=feats, beta=beta_vec[feats],
                    effect_per_1SD=eff, abs_effect_1SD=abs(eff),
                    stringsAsFactors=FALSE)
  out[order(-out$abs_effect_1SD), , drop=FALSE]
}

# Sample-level calibration: first average the predicted values across training cells within each sample, then fit the model y_true_sample ~ y_pred_sample.
calibrate_model_sample <- function(y_true_cell, y_pred_cell, sample_id) {
  df <- data.frame(s=sample_id, y=y_true_cell, p=y_pred_cell)
  agg <- df %>% group_by(s) %>% summarise(y=first(y), p=mean(p, na.rm=TRUE), .groups="drop")
  y_true <- agg$y; y_pred <- agg$p
  ok <- is.finite(y_true) & is.finite(y_pred)
  y_true <- y_true[ok]; y_pred <- y_pred[ok]
  if (length(y_true) < 3 || sd(y_pred) < 1e-8) {
    a <- mean(y_true) - mean(y_pred); b <- 1; return(c(a=a,b=b))
  }
  fit <- try(lm(y_true ~ y_pred), silent=TRUE)
  if (inherits(fit,"try-error")) {
    b <- cov(y_true, y_pred) / var(y_pred); if (!is.finite(b)) b <- 1
    a <- mean(y_true) - b*mean(y_pred); return(c(a=a,b=b))
  }
  co <- coef(fit); if (any(!is.finite(co))) {
    b <- cov(y_true, y_pred) / var(y_pred); if (!is.finite(b)) b <- 1
    a <- mean(y_true) - b*mean(y_pred); return(c(a=a,b=b))
  }
  c(a=unname(co[1]), b=unname(co[2]))
}
apply_calibration <- function(yhat, calib) {
  a <- suppressWarnings(as.numeric(calib["a"])); if (!is.finite(a)) a <- 0
  b <- suppressWarnings(as.numeric(calib["b"])); if (!is.finite(b)) b <- 1
  a + b * yhat
}

# Training subsample: draw cells proportionally by sample (to prevent any single large sample from dominating).
cell_subsample_by_sample <- function(all_idx, sample_id, frac=SUBSAMPLE_FRAC, replace=RESAMPLE_REPLACE, seed=NULL) {
  if (!is.null(seed)) {
    old <- if (exists(".Random.seed", envir=.GlobalEnv)) get(".Random.seed", envir=.GlobalEnv) else NULL
    on.exit({ if (!is.null(old)) assign(".Random.seed", old, envir=.GlobalEnv) }, add=TRUE)
    set.seed(seed)
  }
  sel <- integer(0)
  for (s in unique(sample_id)) {
    idx_s <- all_idx[sample_id == s]
    k <- if (replace) max(1, round(length(idx_s) * frac)) else max(1, floor(length(idx_s) * frac))
    sel <- c(sel, sample(idx_s, k, replace=replace))
  }
  sort(unique(sel))
}

# Sample equal weighting: the weight of each cell within a sample is set to 1 / (number of cells in that sample).
make_sample_weights <- function(sample_id_vec) {
  tab <- table(sample_id_vec)
  1 / as.numeric(tab[as.character(sample_id_vec)])
}


predict_cells_with_model <- function(model, X_block) {
  feats <- intersect(model$feats, colnames(X_block))
  if (length(feats) == 0) return(rep(NA_real_, nrow(X_block)))
  pr_raw <- as.numeric(predict(model$cvfit, newx = to_x(X_block[, feats, drop=FALSE]), s=model$s_choice))
  apply_calibration(pr_raw, model$calib)
}

## ===================== 4) Data preparation (looping over cell types) =====================
md <- combination@meta.data

mat <- GetAssayData(combination, assay = ASSAY, slot = DATA_SLOT)  # genes x cells
stopifnot(inherits(mat, "dgCMatrix"))

ct_levels <- unique(as.character(md[[CELLTYPE_COL]]))[slt]

message("number of cell types: ", length(ct_levels))

for (ct in ct_levels) {
  #ct <- 'BC'
  message("\n===== Running cell type: ", ct, " =====")
  out_dir <- file.path(out_root, paste0("celltype_", make.names(ct)))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  cells_ct <- rownames(md)[md[[CELLTYPE_COL]] == ct]
  if (length(cells_ct) < 10) {
    message("[Skip] too few cells：", ct)
    next
  }
  
  X_ct <- Matrix::t(mat[, cells_ct, drop=FALSE])  # cells x genes
  meta_ct <- md[cells_ct, , drop=FALSE]
  
  # Global gene QC (within this cell type)
  det_rate <- Matrix::colMeans(X_ct > 0) 
  keep_gene <- which(det_rate >= MIN_PCT_CELLS)
  X_ct <- X_ct[, keep_gene, drop=FALSE]
  var_ok <- apply(X_ct, 2, var) > MIN_VAR_GENE
  X_ct <- X_ct[, var_ok, drop=FALSE]
  message("GeneNumber (after QC): ", ncol(X_ct))
  
  samp_vec <- as.character(meta_ct[[SAMPLE_COL]])
  grp_vec  <- as.character(meta_ct[[GROUP_COL]])
  age_vec  <- as.numeric(meta_ct[[AGE_COL]])
  stopifnot(all(grp_vec %in% c("Y","M","O","CR")))
  
  idx_ctrl <- which(grp_vec %in% c("Y","M","O"))
  idx_cr   <- which(grp_vec == "CR")
  
  X_ctrl    <- X_ct[idx_ctrl, , drop=FALSE]
  y_ctrl    <- age_vec[idx_ctrl]
  grp_ctrl  <- grp_vec[idx_ctrl]
  samp_ctrl <- samp_vec[idx_ctrl]
  
  X_cr      <- if (length(idx_cr)>0) X_ct[idx_cr, , drop=FALSE] else NULL
  y_cr      <- if (length(idx_cr)>0) age_vec[idx_cr] else NULL
  samp_cr   <- if (length(idx_cr)>0) samp_vec[idx_cr] else character(0)
  
  o_samples <- unique(samp_ctrl[grp_ctrl == "O"])
  if (length(o_samples) < 2) {
    message("[Skip] O-group samples are too few to perform LOOCV: ", ct)
    next
  }
  message("Number of O-group samples:", length(o_samples))
  
  table_train <- list()
  table_test  <- list()
  model_store <- list()
  
  # Outer layer: leave-one-O-sample-out (the test set consists of all cells from that O sample)
  for (fold_idx in seq_along(o_samples)) {
    test_samp <- o_samples[fold_idx]
    
    te_idx <- which(samp_ctrl == test_samp)
    tr_idx <- setdiff(seq_len(nrow(X_ctrl)), te_idx)
    
    X_tr_full <- X_ctrl[tr_idx, , drop=FALSE]
    y_tr_full <- y_ctrl[tr_idx]
    g_tr_full <- grp_ctrl[tr_idx]
    s_tr_full <- samp_ctrl[tr_idx]
    
    X_te_full <- X_ctrl[te_idx, , drop=FALSE]
    y_te_full <- y_ctrl[te_idx]
    s_te_full <- samp_ctrl[te_idx]
    
    feats_current <- colnames(X_tr_full)
    
    for (round_id in seq_len(MAX_ROUNDS)) {
      bag_features <- list()
      bag_effects  <- list()
      
      for (b in seq_len(B_REPEATS)) {
        seed_sub <- 10100000L + fold_idx*10000L + round_id*100L + b
        seed_cv  <- 20200000L + fold_idx*10000L + round_id*100L + b
        
        idx_all <- seq_len(nrow(X_tr_full))
        sel     <- cell_subsample_by_sample(idx_all, s_tr_full, frac=SUBSAMPLE_FRAC,
                                            replace=RESAMPLE_REPLACE, seed=seed_sub)
        
        Xb <- X_tr_full[sel, feats_current, drop=FALSE]
        yb <- y_tr_full[sel]
        sb <- s_tr_full[sel]
        
        wb <- if (USE_SAMPLE_WEIGHTS) make_sample_weights(sb) else NULL
        
        best <- inner_cv_fit(Xb, yb, sample_id=sb, weights=wb, seed=seed_cv)
        co   <- coef_1se_or_min(best$cvfit)
        imp  <- one_sd_importance(co$beta, Xb)
        
        # Training side: cell-level prediction → sample mean → sample-level calibration
        pred_tr_raw_cell <- as.numeric(predict(best$cvfit, newx=to_x(Xb), s=co$s))
        calib <- if (USE_CALIBRATION) {
          calibrate_model_sample(y_true_cell=yb, y_pred_cell=pred_tr_raw_cell, sample_id=sb)
        } else c(a=0,b=1)
        
        # Testing side: predict each cell of the O sample individually
        pred_te_raw_cell <- as.numeric(predict(best$cvfit, newx=to_x(X_te_full[, feats_current, drop=FALSE]), s=co$s))
        pred_te_cal_cell <- apply_calibration(pred_te_raw_cell, calib)
        # Compute the model-selection error for the O sample using its sample-mean prediction
        p_te_mean_raw <- mean(pred_te_raw_cell, na.rm=TRUE)
        p_te_mean_cal <- mean(pred_te_cal_cell, na.rm=TRUE)
        y_te_true     <- unique(y_te_full)[1]
        abs_err_raw   <- abs(p_te_mean_raw - y_te_true)
        abs_err_cal   <- abs(p_te_mean_cal - y_te_true)
        
        # Record
        table_train[[length(table_train)+1]] <- data.frame(
          fold=fold_idx, round=round_id, b=b,
          n_train_cells=nrow(Xb), n_train_samples=n_distinct(sb),
          inner_cv_mae=best$score,
          n_features=length(feats_current), alpha=best$alpha, lambda_choice=co$s
        )
        table_test[[length(table_test)+1]] <- data.frame(
          fold=fold_idx, round=round_id, b=b, test_sample=test_samp,
          age_true=y_te_true,
          age_pred= p_te_mean_raw, age_pred_cal=p_te_mean_cal,
          abs_err_raw=abs_err_raw, abs_err_cal=abs_err_cal,
          n_features=length(feats_current), alpha=best$alpha, lambda_choice=co$s
        )
        
        # save model
        key <- sprintf("fold%02d_round%02d_b%02d", fold_idx, round_id, b)
        model_store[[key]] <- list(
          fold=fold_idx, round=round_id, b=b, o_test_sample=test_samp,
          feats=feats_current, alpha=best$alpha, s_choice=co$s,
          cvfit=best$cvfit, calib=calib, inner_cv_mae=best$score,
          seed_sub=seed_sub, seed_cv=seed_cv
        )
        
        
        if (nrow(imp) > 0) {
          k_keep <- max(1, floor(nrow(imp) * TOP_NONZERO_FRAC))
          bag_features[[length(bag_features)+1]] <- imp$feature[seq_len(k_keep)]
          bag_effects[[length(bag_effects)+1]]  <- imp[seq_len(k_keep), c("feature","abs_effect_1SD")]
        }
      } # end B
      
      if (length(bag_features) == 0) {
        feats_next <- if (length(feats_current) > TARGET_MIN_GENES)
          feats_current[seq_len(max(TARGET_MIN_GENES, floor(length(feats_current)/2)))]
        else feats_current
      } else {
        all_feats <- unlist(bag_features, use.names=FALSE)
        freq_tab  <- sort(table(all_feats), decreasing=TRUE)
        freq      <- as.numeric(freq_tab); names(freq) <- names(freq_tab)
        freq      <- freq / length(bag_features)
        eff_df    <- do.call(rbind, bag_effects)
        med_eff   <- tapply(eff_df$abs_effect_1SD, eff_df$feature, median, na.rm=TRUE)
        med_eff   <- med_eff[names(freq)]
        stats_r   <- data.frame(feature=names(freq), freq=freq, med_abs_eff=as.numeric(med_eff))
        stats_r   <- stats_r[order(-stats_r$freq, -stats_r$med_abs_eff), ]
        
        pi_thr <- FREQ_THR_INIT
        sel <- stats_r$feature[stats_r$freq >= pi_thr]
        while (length(sel) < TARGET_MIN_GENES && pi_thr > 0) {
          pi_thr <- max(0, pi_thr - FREQ_THR_STEP)
          sel <- stats_r$feature[stats_r$freq >= pi_thr]
        }
        if (length(sel) < TARGET_MIN_GENES) sel <- head(stats_r$feature, TARGET_MIN_GENES)
        if (length(sel) > length(feats_current)) sel <- head(sel, length(feats_current))
        feats_next <- sel
      }
      
      if (length(feats_next) <= TARGET_MIN_GENES) break
      feats_current <- feats_next
    } # end ROUNDS
  } # end FOLDS
  
  
  table1_train_df <- if (length(table_train)) do.call(rbind, table_train) else data.frame()
  table2_test_df  <- if (length(table_test))  do.call(rbind, table_test)  else data.frame()
  
  # Select the optimal (round, b) for each O sample”
  tbl_err <- table2_test_df %>% mutate(abs_err_sel = ifelse(is.finite(abs_err_cal), abs_err_cal, abs_err_raw))
  selected_pairs_df <- tbl_err %>%
    group_split(test_sample) %>%
    lapply(function(df) {
      df2 <- df %>% filter(round >= 2)
      if (nrow(df2) > 0) df2 %>% arrange(abs_err_sel) %>% slice(1)
      else df %>% arrange(abs_err_sel) %>% slice(1)
    }) %>%
    bind_rows() %>%
    select(fold, round, b, test_sample, age_true, age_pred_cal, abs_err_cal, abs_err_sel)
  
  
  selected_keys <- sprintf("fold%02d_round%02d_b%02d", selected_pairs_df$fold, selected_pairs_df$round, selected_pairs_df$b)
  selected_models <- lapply(selected_keys, function(k) model_store[[k]])
  names(selected_models) <- selected_keys
  keep <- !vapply(selected_models, is.null, logical(1))
  selected_models   <- selected_models[keep]
  selected_pairs_df <- selected_pairs_df[keep, , drop=FALSE]
  stopifnot(length(selected_models) > 0)
  
  ### cell-level prediction -----------------
  # A) O-group cells: predict each cell using the selected model for that sample
  o_cells <- rownames(meta_ct)[meta_ct[[GROUP_COL]]=="O"]
  pred_o_rows <- list()
  if (length(o_cells) > 0) {
    o_samples_all <- unique(meta_ct[o_cells, SAMPLE_COL])
    for (s in o_samples_all) {
      hit <- which(selected_pairs_df$test_sample == s)
      if (length(hit) == 0) next
      mdl_key <- sprintf("fold%02d_round%02d_b%02d",
                         selected_pairs_df$fold[hit], selected_pairs_df$round[hit], selected_pairs_df$b[hit])
      mdl <- model_store[[mdl_key]]
      if (is.null(mdl)) next
      
      idx_cells_s <- which(meta_ct[[SAMPLE_COL]] == s & meta_ct[[GROUP_COL]]=="O")
      X_s <- X_ct[idx_cells_s, , drop=FALSE]
      yhat <- predict_cells_with_model(mdl, X_s)
      
      pred_o_rows[[length(pred_o_rows)+1]] <- data.frame(
        cell   = rownames(X_ct)[idx_cells_s],
        sample = as.character(s),
        group  = "O",
        age_true = as.numeric(meta_ct[[AGE_COL]][idx_cells_s]),
        age_pred = yhat,
        model_key = mdl_key,
        stringsAsFactors = FALSE
      )
    }
  }
  pred_o_cells_df <- if (length(pred_o_rows)) bind_rows(pred_o_rows) else data.frame()
  
  # B) CR-group cells: predict each cell using the mean of all selected models
  cr_cells <- rownames(meta_ct)[meta_ct[[GROUP_COL]]=="CR"]
  pred_cr_cells_df <- data.frame()
  if (length(cr_cells) > 0) {
    X_cr_all <- X_ct[cr_cells, , drop=FALSE]
    pred_mat <- sapply(selected_models, function(m) predict_cells_with_model(m, X_cr_all))
    if (is.null(dim(pred_mat))) pred_mat <- matrix(pred_mat, ncol=1)
    yhat_cr <- rowMeans(pred_mat, na.rm=TRUE)
    pred_cr_cells_df <- data.frame(
      cell   = cr_cells,
      sample = as.character(meta_ct[cr_cells, SAMPLE_COL]),
      group  = "CR",
      age_true = as.numeric(meta_ct[cr_cells, AGE_COL]),
      age_pred = yhat_cr,
      model_key = NA_character_,
      stringsAsFactors = FALSE
    )
  }
  
  # C) Y/M-group cells: predict each cell using the mean of the selected model set (for visualization only, not for performance claims)
  ym_cells <- rownames(meta_ct)[meta_ct[[GROUP_COL]] %in% c("Y","M")]
  pred_ym_cells_df <- data.frame()
  if (length(ym_cells) > 0) {
    X_ym_all <- X_ct[ym_cells, , drop=FALSE]
    pred_mat <- sapply(selected_models, function(m) predict_cells_with_model(m, X_ym_all))
    if (is.null(dim(pred_mat))) pred_mat <- matrix(pred_mat, ncol=1)
    yhat_ym <- rowMeans(pred_mat, na.rm=TRUE)
    pred_ym_cells_df <- data.frame(
      cell   = ym_cells,
      sample = as.character(meta_ct[ym_cells, SAMPLE_COL]),
      group  = as.character(meta_ct[ym_cells, GROUP_COL]),
      age_true = as.numeric(meta_ct[ym_cells, AGE_COL]),
      age_pred = yhat_ym,
      model_key = NA_character_,
      stringsAsFactors = FALSE
    )
  }
  
  # Merge cell-level prediction tables
  per_cell_pred <- bind_rows(pred_ym_cells_df, pred_o_cells_df, pred_cr_cells_df) %>%
    mutate(group = factor(group, levels=c("Y","M","O","CR")),
           abs_err = abs(age_pred - age_true),
           celltype = ct)
  
  mae_cells_O <- per_cell_pred %>%
    filter(group == "O") %>%
    summarise(MAE_cells_O = mean(abs_err, na.rm=TRUE), n_cells = n()) %>%
    mutate(celltype = ct)
  
  mae_cells_O_by_sample <- per_cell_pred %>%
    filter(group == "O") %>%
    group_by(sample) %>%
    summarise(MAE_cells_O = mean(abs_err, na.rm=TRUE), n_cells = n(), .groups="drop") %>%
    mutate(celltype = ct)
  
  ## ----------------- save results -----------------
  write.csv(table1_train_df, file.path(out_dir, "table1_train_log_cellwise.csv"), row.names=FALSE)
  write.csv(table2_test_df,  file.path(out_dir, "table2_test_O_meanPred_all_rb_cellwise.csv"), row.names=FALSE)
  write.csv(selected_pairs_df, file.path(out_dir, "selected_pairs_perO_cellwise.csv"), row.names=FALSE)
  write.csv(per_cell_pred, file.path(out_dir, "per_cell_predictions_YMOCR.csv"), row.names=FALSE)
  write.csv(mae_cells_O, file.path(out_dir, "MAE_cells_O_overall.csv"), row.names=FALSE)
  write.csv(mae_cells_O_by_sample, file.path(out_dir, "MAE_cells_O_by_sample.csv"), row.names=FALSE)
  
  ## save models----------
  model_index_rows <- list()
  for (i in seq_along(selected_models)) {
    mdl <- selected_models[[i]]
    key <- names(selected_models)[i]
    file_i <- file.path(out_dir, paste0("model_", key, ".rds"))
    mdl$note <- "SC per-cell training; per-O best (round,b); sample-level calibration; features fixed within round"
    mdl$n_features <- length(mdl$feats)
    saveRDS(mdl, file_i)
    a <- tryCatch(as.numeric(mdl$calib["a"]), error=function(e) NA_real_)
    b <- tryCatch(as.numeric(mdl$calib["b"]), error=function(e) NA_real_)
    model_index_rows[[length(model_index_rows)+1]] <- data.frame(
      fold=mdl$fold, round=mdl$round, b=mdl$b,
      test_sample = if (!is.null(mdl$o_test_sample)) mdl$o_test_sample else NA_character_,
      n_features  = length(mdl$feats),
      alpha       = if (!is.null(mdl$alpha)) mdl$alpha else NA_real_,
      lambda      = if (!is.null(mdl$s_choice)) as.character(mdl$s_choice) else NA_character_,
      calib_a     = a, calib_b = b,
      inner_cv_mae= if (!is.null(mdl$inner_cv_mae)) mdl$inner_cv_mae else NA_real_,
      rds_path    = file_i,
      stringsAsFactors = FALSE
    )
  }
  if (length(model_index_rows)) {
    write.csv(bind_rows(model_index_rows), file.path(out_dir, "index_models_perO_best.csv"), row.names=FALSE)
  }
  
  ## ----------------- Display -----------------
  plot_df <- per_cell_pred %>% filter(group %in% c("O","CR")) %>%
    mutate(.plot_id = row_number())
  if (nrow(plot_df) > MAX_POINTS_PLOT) {
    set.seed(1)
    plot_df <- plot_df %>% slice_sample(n = MAX_POINTS_PLOT)
  }
  
  R_spearman_O <- per_cell_pred %>%
    filter(group == "O") %>%
    summarise(R = suppressWarnings(cor(age_true, age_pred, method="spearman", use="complete.obs"))) %>%
    pull(R)
  lab_mae  <- sprintf("MAE_cells_O = %.2f", mae_cells_O$MAE_cells_O[1])
  lab_R    <- sprintf("R_O = %.2f", R_spearman_O)
  
  p_scatter <- ggplot(plot_df, aes(x=age_true, y=age_pred, color=group)) +
    geom_abline(slope=1, intercept=0, linetype="dashed") +
    geom_point(alpha=0.35, size=0.6) +
    annotate("label", x = min(plot_df$age_true, na.rm=TRUE),
             y = max(plot_df$age_pred, na.rm=TRUE),
             hjust = 0, vjust = 1,
             label = paste(lab_R, lab_mae, sep=" | "),
             size=3, fill="white", color="black", label.size = 0) +
    labs(title = paste0("Cell-wise age prediction — ", ct),
         x = "Chronological age (year)",
         y = "Predicted age (year, calibrated)") +
    theme_bw()
  ggsave(file.path(out_dir, paste0("scatter_cellwise_", make.names(ct), ".png")),
         plot=p_scatter, width=4.5, height=4.5, dpi=200)
  
  plot_df2 <- per_cell_pred %>% filter(group %in% c("O","CR")) %>%
    mutate(abs_err = (age_pred - age_true))
  p_violin <- ggplot(plot_df2, aes(x=group, y=abs_err, fill=group)) +
    geom_violin(trim=FALSE, alpha=0.7) +
    geom_point(alpha=0.35, size=0.6) +
    geom_boxplot(width=0.12, outlier.shape=NA) +
    labs(title = paste0("Cell-wise |error| — ", ct),
         x = NULL, y = "|Predicted - True| (year)") +
    theme_bw()
  ggsave(file.path(out_dir, paste0("abs_error_violin_", make.names(ct), ".png")),
         plot=p_violin, width=4.5, height=4.2, dpi=200)
  
  message("Saved results under: ", out_dir)
}

message("\nAll done.")




