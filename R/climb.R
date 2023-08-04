#' CLIMB deconvolution 
#'
#' By using a scRNA-seq reference dataset (sc) and a bulk target dataset (bulk), CLIMB will generate cell-type abundance
#' and cell-type gene expression for each mixture given in the bulk data. 
#'
#' @param sc ExpressionSet object containing the scRNA-seq reference dataset
#' @param bulk ExpressionSet object containing the mixtures to be deconvoluted
#' @param ratio_cancer_cells numeric vector, same size as number of mixtures. if we have data about the fraction of cancer cells present in mixture, this can improve deconvolution accuracy. If ratio_cancer_cells is provided, then CLIMB-BA method will be used. It is important that cancer cell-types in the scRNA-seq dataset contains the pattern "-like" to be recognized as cancer cells.
#' @param norm_coefs boolean indicating whether coefficients should be normalized by total RNA content
#' @param dwls_weights boolean indicating whether a 2-pass procedure should be applied, with DWLS-like gene-specific weights applied on the 2nd pass.
#' @param verbose boolean indicating whether to print message during running
#' @param up.lim numeric scalar, impose a l-infinity norm to the linear model (upper bound for coefficient)
#' @param lambda Regularization factor
#' @param norm_factor if ratio_cancer_cells is provided, indicate strength of smooth normalization
#' @param cancer_pattern a string pattern present in all cancer cell-type. Only for these cell-types CLIMB will assume the presence of differentially expressed genes
#' @param mode indicate which mode to use between: ['all'] run bulk deconvolution of cell-type proportions - cell-type expression - and differential expression analysis (without per sample DE analysis). ['all+'] same as 'all' but adds per-sample DE analysis (takes long time to run). ['abundance'] only performs cell-type proportions inference. ['expression'] performs cell-type proportions together with cell-type expression prediction.
#' @param predict_abundance boolean indicating whether cell-type proportions should be assessed. This needs to be true for cell-type expression prediction to run.
#' @param predict_expression boolean indicating whether cell-type expression should be predicted (takes longer to run)
#' @param DE_analysis boolean indicating whether to perform DE analysis between two conditions. Requires the vector 'conditions' to be given as input.
#' @param conditions vector associated with each bulk samples. Should be the same size as the number of bulks, and should contains two values only: 'condition' and 'control'. For instance, disease samples (condition) versus healthy (control), or KO samples (condition) versus wild type (control).
#' @param patient_specific_DE boolean indicating whether to perform DE for each sample individually. In that case, each sample will be compared to the 'control' group in the vector 'conditions'. If no vector conditions is provided, every other samples will be used as control samples. 
#' @param final_res previous CLIMB result to use for DE analysis.
#' @param min_common_genes minimum number of genes to be in common between bulk and scRNA-seq datsets 
#' @export
climb <- function (sc, bulk, cancer_pattern = "none", mode = "abundance", norm_coefs = TRUE, dwls_weights=TRUE,
    ratio_cancer_cells = NA, up.lim = Inf, lambda = 0, norm_factor = 0.1, verbose=TRUE,
    conditions = NA, predict_abundance = TRUE, predict_expression = TRUE, 
    DE_analysis = FALSE, patient_specific_DE = FALSE, final_res = list(), 
    min_common_genes = 100) 
{
    if (mode == "all") {
        if(verbose){message("ALL mode: prediction of cell-type abundance / high-resolution cell-type expression / DE analysis between conditions")}
        predict_abundance = TRUE
        predict_expression = TRUE
        DE_analysis = TRUE
        patient_specific_DE = FALSE
        stopifnot(!all(is.na(conditions)))
    }
    if (mode == "all+") {
        if(verbose){message("ALL+ mode: prediction of cell-type abundance / high-resolution cell-type expression / DE analysis between conditions if provided AND at sample level")}
        if(verbose){message("WARNING: sample-level DE can take long!")}
        predict_abundance = TRUE
        predict_expression = TRUE
        DE_analysis = TRUE
        patient_specific_DE = TRUE
    }
    if (mode == "abundance") {
        if(verbose){message("ABUNDANCE mode: predicting cell-type proportions in bulks")}
        predict_abundance = TRUE
        predict_expression = FALSE
        DE_analysis = FALSE
        patient_specific_DE = FALSE
    }
    if (mode == "expression") {
        if(verbose){message("EXPRESSION mode: predicting cell-type expression in bulks - requires single-cell coefficients fitted by CLIMB")}
        predict_abundance = TRUE
        predict_expression = TRUE
        DE_analysis = FALSE
        patient_specific_DE = FALSE
    }
    if (mode == "DE.only") {
        if(verbose){message("Running DE analysis based on existing CLIMB object")}
        predict_abundance = FALSE
        predict_expression = FALSE
        DE_analysis = TRUE
        patient_specific_DE = FALSE
        stopifnot(length(final_res) > 0)
    }
    num <- function(x) {
        return(as.numeric(as.character(x)))
    }
    ct.props = list()
    ct.exprs = list()
    sc.mat = exprs(sc)
    cell_expr = colSums(exprs(sc))
    save_coefs = list()
    sc$cellType = factor(sc$cellType)
    cellTypes = levels(sc$cellType)
    common_genes = intersect(rownames(bulk), rownames(sc))
    if(verbose){message(paste0(length(common_genes), " common genes found between scRNA-seq refererence and bulk datasets"))}
    if (length(common_genes) < min_common_genes) {
        stop("too few genes found between scRNA-seq refererence and bulk dataset")
    }
    sc = sc[common_genes, ]
    scmat = exprs(sc)
    bulk = bulk[common_genes, ]
    N = dim(bulk)[2]
    G = dim(bulk)[1]
    K = length(cellTypes)
    S_pred_mapping_n = array(rep(0, N * G * K), c(N, G, K))
    if (predict_abundance) {
        if(verbose){message("Bulk to single-cell mapping for prediction of cell-type abundance / expression")}
        for (i in 1:N) {
            y = num(exprs(bulk)[, i])
            if (dwls_weights){
                # CLIMB first pass
                fit = glmnet(scmat, y, lower.limits = 0, upper.limits = up.lim, standardize=T)
                coefs = coef(fit)[-1, dim(coef(fit))[2]]
                # Implement DWLS-like weights using single-cell expression matrix
                alpha_tilde_tm1 = coefs/(sum(coefs))
                weights_tm1 = 1 / (scmat %*% num(alpha_tilde_tm1))^2
                weights_tm1[weights_tm1 == Inf] <- max(weights_tm1[weights_tm1 != Inf])
                q_weights_tm1 = quantile(weights_tm1, probs = seq(0,1,0.01))
                weights_tm1[weights_tm1 > q_weights_tm1[95]] <- q_weights_tm1[95]
                # Run second pass of CLIMB
                fit = glmnet(scmat, y, standardize = T, lower.limits = 0.0, lambda=0.0, weights = weights_tm1)
                coefs = coef(fit)[-1,dim(coef(fit))[2]]
            } else {
                # CLIMB one-pass original (no weights)
                fit = glmnet(scmat, y, lower.limits = 0.0, lambda=0.0, upper.limits = up.lim, standardize=T)
                coefs = coef(fit)[-1, dim(coef(fit))[2]]
            }
            if (norm_coefs) { 
                agg = aggregate(coefs / cell_expr, list(sc$cellType), sum, drop = F)
            } else {
                agg = aggregate(coefs, list(sc$cellType), sum, drop = F)
            }
            agg$x[is.na(agg$x)] <- 0
            ppred = (agg$x)/sum(agg$x)
            names(ppred) = agg$Group.1
            ct.props[[i]] = ppred
            if (predict_expression) {
                pcor_expr_pred = list()
                all_celltypes = levels(sc$cellType)
                pred_exprs = list()
                for (k in 1:length(all_celltypes)) {
                  this_ct = all_celltypes[k]
                  sel_ct = sc$cellType == this_ct
                  pred_expr = (t(coefs)[sel_ct] %*% t(scmat[, 
                    sel_ct]))
                  pred_expr[is.na(pred_expr)] = 0
                  pred_exprs[[k]] = pred_expr
                  S_pred_mapping_n[, , k] = pred_expr
                }
                ct_exprs_pred = do.call(rbind, pred_exprs)
                ct.exprs[[i]] = ct_exprs_pred
            }
            save_coefs[[i]] = coefs
        }
        final_res$props = do.call(rbind, ct.props)
        if(verbose){message("Cell-type abundance prediction done. ")}
    }
    if (predict_expression) {
        if(verbose){message("Starting high-resolution expression deconvolution")}
        if ( cancer_pattern == 'none' ){
            for (g in 1:G) {
                for (n in 1:N) {
                    S_pred_mapping_n[n, g, ] = ct.exprs[[n]][,g]
                }
            }
            dimnames(S_pred_mapping_n)[[1]] = colnames(bulk)
            dimnames(S_pred_mapping_n)[[2]] = rownames(bulk)
            dimnames(S_pred_mapping_n)[[3]] = cellTypes
            final_res$expr.highres = S_pred_mapping_n
            final_res$expr.mapping = S_pred_mapping_n
            final_res$expr.overall = colSums(S_pred_mapping_n, dims = 1) 
            final_res$coefs = save_coefs
        } else {
            normal_sel = !grepl(cancer_pattern, sc$cellType)
            cancer_sel = grepl(cancer_pattern, sc$cellType)
            cancer_ct_sel = grepl(cancer_pattern, cellTypes)
            alpha_overal = do.call(rbind, save_coefs)
            alpha_cancer = do.call(rbind, save_coefs)[, cancer_sel]
            sc.cancer = sc[, cancer_sel]
            C_overal = exprs(sc)
            Y_hat_overal = alpha_overal %*% t(C_overal)
            Y_true_bulk = t(exprs(bulk))
            Epsilon_ng = Y_true_bulk - Y_hat_overal
            S_pred_n = array(rep(0, N * G * K), c(N, G, K))
            for (g in 1:G) {
                Epsilon_g = num(Epsilon_ng[, g])
                if (sd(Epsilon_g) != 0) {
                    fit = glmnet(alpha_cancer, Epsilon_g, intercept = TRUE)
                    C_diff_cancer = num(coef(fit)[-1, dim(coef(fit))[2]])
                    for (n in 1:N) {
                      Epsilon_cancer_n = aggregate(C_diff_cancer * 
                        alpha_cancer[n, ], list(sc.cancer$cellType), 
                        sum)$x
                      Epsilon_n = rep(0, K)
                      Epsilon_n[cancer_ct_sel] = Epsilon_cancer_n
                      S_pred_mapping_n[n, g, ] = ct.exprs[[n]][, 
                        g]
                      S_pred_n[n, g, ] = S_pred_mapping_n[n, g, ] + 
                        Epsilon_n
                      S_pred_n[n, g, ][S_pred_n[n, g, ] < 0] = 0
                    }
                }
                else {
                    S_pred_n[n, g, ] = rep(0, K)
                }
                if (g%%1000 == 0) {
                    if(verbose){message(paste0("High-Resolution expression prediction: ", 
                      g, " genes processed..."))}
                }
            }
            dimnames(S_pred_n)[[1]] = dimnames(S_pred_mapping_n)[[1]] = colnames(bulk)
            dimnames(S_pred_n)[[2]] = dimnames(S_pred_mapping_n)[[2]] = rownames(bulk)
            dimnames(S_pred_n)[[3]] = dimnames(S_pred_mapping_n)[[3]] = cellTypes
            final_res$expr.highres = S_pred_n
            final_res$expr.mapping = S_pred_mapping_n
            final_res$expr.overall = colSums(S_pred_mapping_n, dims = 1)
            final_res$coefs = save_coefs
        }
    }
    else {
        final_res$expr.highres = S_pred_mapping_n
        final_res$expr.mapping = S_pred_mapping_n
        final_res$expr.overall = colSums(S_pred_mapping_n, dims = 1) 
        final_res$coefs = save_coefs
    }
    if (!all(is.na(conditions))) {
        DE_analysis = TRUE
    }
    if (DE_analysis == TRUE) {
        if(verbose){message("Starting DE analysis")}
        define_signif <- function(p_) {
            signif_p_ = p_
            signif_p_[p_ < 0.05] = "*"
            signif_p_[p_ < 0.01] = "**"
            signif_p_[p_ < 0.001] = "***"
            signif_p_[p_ < 1e-04] = "****"
            signif_p_[p_ >= 0.05] = "n.s"
            return(signif_p_)
        }
        S_mat = round(final_res$expr.pred)
        ct_prop = final_res$props
        tot_expr = num(colSums(exprs(bulk)))/1e+06
        N = dim(S_mat)[1]
        G = dim(S_mat)[2]
        K = dim(S_mat)[3]
        pvals = array(rep(1, G * K), c(G, K))
        padjs = array(rep(1, G * K), c(G, K))
        fcs = array(rep(0, G * K), c(G, K))
        pvals.N = array(rep(1, N * G * K), c(N, G, K))
        padjs.N = array(rep(1, N * G * K), c(N, G, K))
        fcs.N = array(rep(1, N * G * K), c(N, G, K))
        ct.pvals = array(rep(1, K), c(K))
        ct.padjs = array(rep(1, K), c(G, K))
        ct.fcs = array(rep(0, K), c(K))
        if (!all(is.na(conditions))) {
            if(verbose){message(paste0("DE analysis of cell-type proportions"))}
            colData = data.frame(condition = conditions, tot_expr = tot_expr)
            w_k_per1000cells = round(t(1000 * ct_prop)) + 1
            dds.ct <- DESeqDataSetFromMatrix(countData = w_k_per1000cells, 
                colData = colData, design = (~tot_expr + condition))
            dds.ct <- tryCatch(expr = {
                DESeq(dds.ct)
            }, error = function(cond) {
                if(verbose){message("Error with DE analysis, using fit with mean instead")}
                return(DESeq(dds.ct, fitType = "mean"))
            })
            df_res.ct = results(dds.ct)[order(results(dds.ct)$pvalue, 
                decreasing = F), ]
        }
        for (k in 1:K) {
            if(verbose){message(paste0("DE analysis of cell-type ", dimnames(S_mat)[[3]][k]))}
            if (!all(is.na(conditions))) {
                if(verbose){message(paste0("DE analysis between two conditions: ",
                  unique(conditions)[1], " vs ", unique(conditions)[2]))}
                S_k = S_mat[, , k]
                w_k = ct_prop[, k]
                if (rankMatrix(S_k)[1] == 0) {
                  if(verbose){message("matrix not full ranked, skipping cell-type")}
                  next
                }
                colData = data.frame(condition = conditions, 
                  celltype_prop = w_k, tot_expr = tot_expr)
                dds <- DESeqDataSetFromMatrix(countData = t(S_k) + 
                  1, colData = colData, design = (~tot_expr + 
                  celltype_prop + condition))
                dds <- tryCatch(expr = {
                  DESeq(dds)
                }, error = function(cond) {
                  if(verbose){message("Error with DE analysis, using fit with mean instead")}
                  return(DESeq(dds, fitType = "mean"))
                })
                fcs[, k] = num(-1 * results(dds, tidy = TRUE)$log2FoldChange)
                pvals[, k] = num(results(dds, tidy = TRUE)$pvalue)
                padjs[, k] = num(results(dds, tidy = TRUE)$padj)
            }
            if (patient_specific_DE) {
                if(verbose){message("DE analysis per sample")}
                for (n in 1:N) {
                  cond_temp = conditions
                  cond_temp[n] <- "condition_n"
                  sel.n = cond_temp != "condition"
                  conditions.n = cond_temp[sel.n]
                  S_k.n = S_k[sel.n, ]
                  w_k.n = w_k[sel.n]
                  if (rankMatrix(S_k.n)[1] == 0) {
                    if(verbose){message("matrix not full ranked, skipping cell-type")}
                    next
                  }
                  tot_expr.n = num(total_expr/1e+06)[sel.n]
                  colData = data.frame(condition = conditions.n, 
                    celltype_prop = w_k.n, tot_expr = tot_expr.n)
                  dds <- DESeqDataSetFromMatrix(countData = t(S_k.n) + 
                    1, colData = colData, design = (~tot_expr + 
                    celltype_prop + condition))
                  dds <- DESeq(dds)
                  dds <- tryCatch(expr = {
                    DESeq(dds)
                  }, error = function(cond) {
                    if(verbose){message("Error with DE analysis, using fit with mean instead")}
                    return(DESeq(dds, fitType = "mean"))
                  })
                  fcs.N[n, , k] = num(-1 * results(dds, tidy = TRUE)$log2FoldChange)
                  pvals.N[n, , k] = num(results(dds, tidy = TRUE)$pvalue)
                  padjs.N[n, , k] = num(results(dds, tidy = TRUE)$padj)
                }
            }
        }
        log_pvals = array(-log10(pvals), c(G, K))
        log_padj = array(-log10(padjs), c(G, K))
        signif_pred = define_signif(padjs)
        dimnames(pvals) = dimnames(log_pvals) = dimnames(padjs) = dimnames(log_padj) = dimnames(signif_pred) = dimnames(fcs) = dimnames(S_mat[1, 
            , ])
        df_res = cbind(melt(pvals), melt(log_pvals)[, 3], melt(padjs)[, 
            3], melt(log_padj)[, 3], melt(fcs)[, 3], melt(signif_pred)[, 
            3])
        colnames(df_res) = c("gene", "celltype", "pval", "log10_pval", 
            "padj", "log10_padj", "log2_FC", "signif")
        df_res = df_res[order(df_res$log10_pval, decreasing = T), 
            ]
        log_pvals.N = array(-log10(pvals.N), c(N, G, K))
        log_padj.N = array(-log10(padjs.N), c(N, G, K))
        signif_pred.N = define_signif(padjs.N)
        dimnames(pvals.N) = dimnames(log_pvals.N) = dimnames(padjs.N) = dimnames(log_padj.N) = dimnames(signif_pred.N) = dimnames(fcs.N) = dimnames(S_mat)
        df_res.N = cbind(melt(pvals.N), melt(log_pvals.N)[, 4], 
            melt(padjs.N)[, 4], melt(log_padj.N)[, 4], melt(fcs.N)[, 
                4], melt(signif_pred.N)[, 4])
        colnames(df_res.N) = c("sample", "gene", "celltype", 
            "pval", "log10_pval", "padj", "log10_padj", "log2_FC", 
            "signif")
        df_res.N = df_res.N[order(df_res.N$log10_pval, decreasing = T), 
            ]
        final_res$DE.expr.conditions = df_res
        final_res$DE.props.conditions = df_res.ct
        final_res$DE.expr.persample
    }
    return(final_res)
} 
