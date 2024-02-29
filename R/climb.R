#' CLIMB deconvolution 
#'
#' By using a scRNA-seq reference dataset (sc) and a bulk target dataset (bulk), CLIMB will generate cell-type abundance
#' and cell-type gene expression for each mixture given in the bulk data. 
#'
#' @param sc ExpressionSet object containing the scRNA-seq reference dataset
#' @param bulk ExpressionSet object containing the mixtures to be deconvoluted
#' @param verbose boolean indicating whether to print message during running
#' @param up.lim numeric scalar, impose a l-infinity norm to the linear model (upper bound for coefficient)
#' @param lambda Regularization factor
#' @param cancer_pattern a string pattern present in all cancer cell-type. Only for these cell-types CLIMB will assume the presence of differentially expressed genes
#' @param mode indicate which mode to use between: ['all'] run bulk deconvolution of cell-type proportions - cell-type expression - and differential expression analysis (without per sample DE analysis). ['all+'] same as 'all' but adds per-sample DE analysis (takes long time to run). ['abundance'] only performs cell-type proportions inference. ['expression'] performs cell-type proportions together with cell-type expression prediction.
#' @param min_common_genes minimum number of genes to be in common between bulk and scRNA-seq datsets 
#' @param ratio_cell_increase percentage increase at each step of the embirical bayes procedure
#' @param n.iter.subsampling number of subsampling that will be performed on the single-cell reference (results from each subsample are then averaged) 
#' @param min.n.cells minimum number of cells per cell type to subsample. If a cell type has less cells in reference, then sampling is done with replacement.
#' @param n.top_mean.genes number of genes used for bulk-specific gene selection.
#' @export

climb <- function (sc, bulk, mode = "abundance", up.lim = Inf, lambda = 0, 
    verbose = TRUE, cancer_pattern = "*", conditions = NA, final_res = list(), 
    min_common_genes = 100, ratio_cell_increase = 0.02, n.iter.subsampling = 5, 
    min.n.cells = 50, n.top_mean.genes = 500) 
{
    if (mode == "all") {
        if (verbose) {
            message("ALL mode: prediction of cell-type abundance / high-resolution cell-type expression / DE analysis between conditions")
        }
        predict_abundance = TRUE
        predict_expression = TRUE
        DE_analysis = TRUE
        patient_specific_DE = FALSE
        stopifnot(!all(is.na(conditions)))
    }
    if (mode == "all+") {
        if (verbose) {
            message("ALL+ mode: prediction of cell-type abundance / high-resolution cell-type expression / DE analysis between conditions if provided AND at sample level")
        }
        if (verbose) {
            message("WARNING: sample-level DE can take long!")
        }
        predict_abundance = TRUE
        predict_expression = TRUE
        DE_analysis = TRUE
        patient_specific_DE = TRUE
    }
    if (mode == "abundance") {
        if (verbose) {
            message("ABUNDANCE mode: predicting cell-type proportions in bulks")
        }
        predict_abundance = TRUE
        predict_expression = FALSE
        DE_analysis = FALSE
        patient_specific_DE = FALSE
    }
    if (mode == "expression") {
        if (verbose) {
            message("EXPRESSION mode: predicting cell-type expression in bulks - requires single-cell coefficients fitted by CLIMB")
        }
        predict_abundance = TRUE
        predict_expression = TRUE
        DE_analysis = FALSE
        patient_specific_DE = FALSE
    }
    if (mode == "DE.only") {
        if (verbose) {
            message("Running DE analysis based on existing CLIMB object")
        }
        predict_abundance = FALSE
        predict_expression = FALSE
        DE_analysis = TRUE
        patient_specific_DE = FALSE
        stopifnot(length(final_res) > 0)
    }
    num <- function(x) {
        return(as.numeric(as.character(x)))
    }
    reformat_strings2 <- function(vector_string) {
        vector_string <- gsub("\\-$", "", vector_string)
        vector_string <- gsub("\\+", "", vector_string)
        vector_string <- gsub("minus", "", vector_string)
        vector_string <- gsub("plus", "", vector_string)
        vector_string <- gsub("\\ ", "\\.", vector_string)
        vector_string <- gsub("[^[:alnum:] ]", "", vector_string)
        return(vector_string)
    }
    reformat_celltypes <- function(celltype_labels) {
        celltype_labels <- reformat_strings2(as.vector(celltype_labels))
        celltype_labels <- factor(celltype_labels)
        return(celltype_labels)
    }
    ct.props = list()
    ct.exprs = list()
    sc.mat = exprs(sc)
    cell_expr = colSums(exprs(sc))
    save_coefs = list()
    sc$cellType = factor(sc$cellType)
    cellTypes = levels(sc$cellType)
    common_genes = intersect(rownames(bulk), rownames(sc))
    if (verbose) {
        message(paste0(length(common_genes), " common genes found between scRNA-seq refererence and bulk datasets"))
    }
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
        if (verbose) {
            message("Bulk to single-cell mapping for prediction of cell-type abundance / expression")
        }
        for (i in 1:N) {
            y = num(exprs(bulk)[, i])
            fit = glmnet(scmat, y, lower.limits = 0, lambda = 0, 
                upper.limits = up.lim, standardize = T)
            coefs = coef(fit)[-1, dim(coef(fit))[2]]
            agg = aggregate(coefs, list(sc$cellType), sum, drop = F)
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
        climb_props.init = data.frame(do.call(rbind, ct.props))
        rownames(climb_props.init) = colnames(bulk)
        colnames(climb_props.init) = levels(sc$cellType)
        final_res$props.init = climb_props.init
        if (verbose) {
            message("First pass of cell-type abundance prediction done. Start second pass...")
        }
        max_size = round(dim(sc)[2]/2)
        sc_mat = exprs(sc)
        median_counts = median(colSums(sc_mat))
        sc_mat_norm = log(t(t(sc_mat)/colSums(sc_mat)) * median_counts + 
            1)
        K = length(unique(sc$cellType))
        ratio_cell_increase = min(0.01 * K, 0.1)
        bulk_mat = exprs(bulk)
        bulk_mat_norm = log(t(t(bulk_mat)/colSums(bulk_mat)) * 
            median_counts + 1)
        N = dim(bulk_mat_norm)[2]
        inter.genes = intersect(rownames(sc_mat_norm), rownames(bulk_mat_norm))
        sc_mat_norm = sc_mat_norm[inter.genes, ]
        sc = sc[inter.genes, ]
        bulk_mat_norm = bulk_mat_norm[inter.genes, ]
        bulk = bulk[inter.genes, ]
        sc_mat_avg = aggregate(t(sc_mat_norm), list(sc$cellType), 
            mean)
        celltypes = levels(sc$cellType)
        celltypes = reformat_celltypes(celltypes)
        rownames(sc_mat_avg) = sc_mat_avg$Group.1
        rownames(sc_mat_avg) = reformat_celltypes(rownames(sc_mat_avg))
        sc_mat_avg = sc_mat_avg[, -1]
        sc_mat_avg = sc_mat_avg[as.character(celltypes), ]
        sel.genes = apply(as.matrix(sc_mat_avg), 2, sd) != 0
        sc_mat_avg = sc_mat_avg[, sel.genes]
        bulk_mat_norm = bulk_mat_norm[sel.genes, ]
        sc_mat_norm = sc_mat_norm[sel.genes, ]
        bulk = bulk[sel.genes, ]
        sc = sc[sel.genes, ]
        all(rownames(sc_mat_avg) == celltypes)
        fc_to_average = t(sc_mat_avg)/colMeans(sc_mat_avg)
        fc_to_second_l = list()
        for (c in 1:ncol(sc_mat_avg)) {
            col_c = sc_mat_avg[, c]
            if (col_c[order(-col_c)][2] > 0) {
                fc_to_second_l[[c]] = col_c/col_c[order(-col_c)][2]
            }
            else {
                fc_to_second_l[[c]] = col_c * 1e-07
            }
        }
        fc_to_second = do.call(rbind, fc_to_second_l)
        double_condition_matrix = fc_to_average > num(quantile(flatten(fc_to_average), 
            p = 0.95)) & fc_to_second > num(quantile(flatten(fc_to_average), 
            p = 0.75))
        double_condition_matrix = double_condition_matrix[rowSums(double_condition_matrix) > 
            0, ]
        celltype_counts = matrix(0, ncol = K, nrow = 1)
        list_genes = list()
        r_ = 1
        for (r in 1:nrow(double_condition_matrix)) {
            row_ = double_condition_matrix[r, ]
            if (celltype_counts[1, which.max(row_)] <= 100) {
                celltype_counts[1, which.max(row_)] = celltype_counts[1, 
                  which.max(row_)] + 1
                list_genes[[r_]] = rownames(double_condition_matrix)[r]
                r_ = r_ + 1
            }
        }
        top_var = unique(unlist(list_genes))
        n_top_var = length(top_var)
        if ((n_top_var/2) < n.top_mean.genes) {
            n.top_mean.genes = round(n_top_var/2)
        }
        sc_mat_norm_v = sc_mat_norm[top_var, ]
        bulk_mat_norm_v = bulk_mat_norm[top_var, ]
        sc_mat_avg_sub = sc_mat_avg[, top_var]
        max_celltype_size = max(num(table(sc$cellType)))
        if (max_celltype_size > max_size) {
            max_celltype_size <- max_size
        }
        min.n.cells = min(min.n.cells, round(0.02 * max_size))
        ref.counts = table(sc$cellType)
        prior_proportions = final_res$props.init
        for (q in 1:2) {
            ct.props = list()
            for (n in 1:N) {
                top_means = -1 * sort(-bulk_mat_norm_v[, n] * 
                  apply(sc_mat_avg_sub, 2, max))
                top_means = top_means[1:n.top_mean.genes]
                gene_topExpr = rev(names(top_means))
                sc_mat_norm_ = sc_mat_norm_v[gene_topExpr, ]
                bulk_mat_norm_ = bulk_mat_norm_v[gene_topExpr, 
                  ]
                y = num(bulk_mat_norm_[, n])
                celltype_counts = matrix(prior_proportions[n, 
                  ] * (max_size), nrow = 1, ncol = K)
                colnames(celltype_counts) = colnames(prior_proportions)
                celltype_counts[celltype_counts < min.n.cells] <- min.n.cells
                sum_props = matrix(0, nrow = 1, ncol = K)
                for (x in 1:n.iter.subsampling) {
                  all_samples = list()
                  for (k in 1:length(celltypes)) {
                    ct = as.character(celltypes[k])
                    set.seed(x * k)
                    all_samples[[k]] = sample(grep(paste0("^", 
                      ct, "$"), sc$cellType), round(num(celltype_counts[, 
                      ct])), replace = T)
                  }
                  sample_ = unlist(all_samples)
                  sc.sub = sc[gene_topExpr, sample_]
                  scmat = sc_mat_norm_[, sample_]
                  fit = glmnet(scmat, y, lower.limits = 0, upper.limits = 0.001)
                  coefs = coef(fit)[-1, dim(coef(fit))[2]]
                  agg = aggregate(coefs, list(sc.sub$cellType), 
                    sum, drop = F)
                  agg$x[is.na(agg$x)] <- 0
                  ppred = (agg$x)/sum(agg$x)
                  names(ppred) = agg$Group.1
                  sum_props = sum_props + ppred
                }
                pred_prop = sum_props/rowSums(sum_props)
                ct.props[[n]] = pred_prop
            }
            climb_props.corrected = data.frame(do.call(rbind, 
                ct.props))
            rownames(climb_props.corrected) = colnames(bulk)
            colnames(climb_props.corrected) = levels(sc$cellType)
            prior_proportions = climb_props.corrected
        }
        final_res$props.corrected = data.frame(climb_props.corrected)
        weights_ = num(cor(table(sc$cellType), t(final_res$props.init)))
        weights_[weights_ < 0] <- 0
        final_res$props.corrected = as.matrix(final_res$props.init) * 
            num(weights_) + as.matrix(final_res$props.corrected) * 
            (1 - num(weights_))
        rownames(final_res$props.corrected) = colnames(bulk)
        colnames(final_res$props.corrected) = levels(sc$cellType)
        if (verbose) {
            message("Second pass is done. Cell-type abundance prediction is over.")
        }
    }
    if (predict_expression) {
        if (verbose) {
            message("Starting high-resolution expression deconvolution")
        }
        if (cancer_pattern == "none") {
            for (g in 1:G) {
                for (n in 1:N) {
                  S_pred_mapping_n[n, g, ] = ct.exprs[[n]][, 
                    g]
                }
            }
            dimnames(S_pred_mapping_n)[[1]] = colnames(bulk)
            dimnames(S_pred_mapping_n)[[2]] = rownames(bulk)
            dimnames(S_pred_mapping_n)[[3]] = cellTypes
            final_res$expr.highres = S_pred_mapping_norm
            final_res$expr.mapping = S_pred_mapping_norm
            final_res$expr.overall = colSums(S_pred_mapping_norm, 
                dims = 1)
            final_res$coefs = save_coefs
        }
        else {
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
                    S_pred_n[n, g, ] = S_pred_mapping_n[n, g, 
                      ] + Epsilon_n
                    S_pred_n[n, g, ][S_pred_n[n, g, ] < 0] = 0
                  }
                }
                else {
                  S_pred_n[n, g, ] = rep(0, K)
                }
                if (g%%1000 == 0) {
                  if (verbose) {
                    message(paste0("High-Resolution expression prediction: ", 
                      g, " genes processed..."))
                  }
                }
            }
            dimnames(S_pred_n)[[1]] = dimnames(S_pred_mapping_n)[[1]] = colnames(bulk)
            dimnames(S_pred_n)[[2]] = dimnames(S_pred_mapping_n)[[2]] = rownames(bulk)
            dimnames(S_pred_n)[[3]] = dimnames(S_pred_mapping_n)[[3]] = cellTypes
            final_res$expr.highres = S_pred_n
            final_res$expr.mapping = S_pred_mapping_norm
            final_res$expr.overall = colSums(S_pred_mapping_norm, 
                dims = 1)
            final_res$coefs = save_coefs
        }
    }
    else {
        final_res$expr.highres = S_pred_mapping_norm
        final_res$expr.mapping = S_pred_mapping_norm
        final_res$expr.overall = colSums(S_pred_mapping_norm, dims = 1)
        final_res$coefs = save_coefs
    }
    return(final_res)
}

### CLIFF functions ### 
sigmoid_ <- function(x, a){
    return( 1 / (1 + exp(-x + a)))
}
mse <- function(tr, pr){ return( sum((tr - pr)^2) / length(tr) ) }
rmse <- function(tr, pr){ return( sqrt( sum((tr - pr)^2) / length(tr) ) ) }
cliff_gene_selection <- function(bulk_mat, drug_d, min.genes=10){
    stopifnot(dim(bulk_mat)[1] == length(drug_d))
    for(i in 1:100){
        la = 10^(-i)
        fit = glmnet(bulk_mat, y=drug_d, lambda=la)
        coefs_ = coef(fit)[-1]
        n_genes = sum(coefs_>0)
        if ( n_genes > min.genes ){
            return(colnames(bulk_mat)[coefs_ > 0])
        }
    }
}
cliff <- function (climb_output, drug_data, mutation_data = NULL, min.mutation = 0,
                    max.em.steps = 100, mode = "highres", regularization = "none", cancer_pattern='like')
{
    message("Prepare CLIFF input from CLIMB output, mutation data, and drug sensitivity data")
	climb_expr = climb_output$expr.highres
	climb_prop = climb_output$props
	rownames(climb_prop) = dimnames(climb_expr)[[1]]
	climb_expr_overall = climb_output$expr.overall
    if (is.null(mutation_data) ) {
		message("No mutation data provided")
		mutation_data = matrix(0, ncol = 2, nrow = dim(climb_prop)[1])
		colnames(mutation_data) = c("a", "b") ; rownames(mutation_data) = rownames(climb_prop)
	}
	mean_auc = mean(drug_data$auc)
	drug_data$auc = drug_data$auc + (0.5 - mean(mean(drug_data$auc)))
	sel.sample = Reduce(intersect, list(dimnames(climb_expr)[[1]], rownames(mutation_data), drug_data$sample))
	order.sample.climb = match(sel.sample, dimnames(climb_expr)[[1]])
	drug_data = drug_data[match(sel.sample, drug_data$sample),]
	climb_expr = climb_expr[order.sample.climb, , ]
	climb_prop = climb_prop[order.sample.climb, ]
    mutation_data = mutation_data[sel.sample, ]
	N = dim(climb_prop)[1]
    if (sum(colSums(mutation_data) >= min.mutation) == 0){
        message("Not enough mutation in selected samples")
        mutation_data = matrix(0, ncol = 2, nrow = N)
        colnames(mutation_data) = c("a", "b") ; rownames(mutation_data) = rownames(climb_prop)
    } else {
        sel.mutation = colSums(mutation_data) >= min.mutation
        mutation_data = mutation_data[, sel.mutation]
        mutation_data = as.matrix(mutation_data)
    }
    # Check that we have the same sample names in good order
    stopifnot(all(dimnames(climb_expr)[[1]] == rownames(mutation_data)))
    stopifnot(all(dimnames(climb_expr)[[1]] == drug_data$sample))
    avg_auc = mean(drug_data$auc)
    message("Prepare CLIFF input from CLIMB output, mutation data, and drug sensitivity data")
    K = num(dim(climb_expr)[3])
    N = num(dim(climb_expr)[1])
    sample_names = dimnames(climb_expr)[[1]]
    tabs_ = list()
    for (n in 1:N) {
        if (mode == "overall") {
            climb_expr[n, , ] = climb_expr_overall
        } else {
            q99 = quantile(climb_expr[n, , ], p=0.95)
            climb_expr[n, , ][climb_expr[n, , ] > q99] <- q99
        }
        df_expr = data.frame(t(log2(climb_expr[n, , ] + 1))/rowSums(t(log2(climb_expr[n,
        , ] + 1)))) * climb_prop[n, ]
        df_expr[is.na(df_expr)] <- 0
        rownames(df_expr) = paste0(rownames(df_expr), "_", sample_names[n])
        df_expr$props = climb_prop[n, ]
        mat.mut = matrix(0, ncol = length(mutation_data[n, ]),
            nrow = length(colnames(climb_prop)))
        mat.mut[grepl(cancer_pattern, colnames(climb_prop)),
            ] = mutation_data[n, ]
        colnames(mat.mut) = paste0("mut.", colnames(mutation_data))
        df_expr_0 = cbind(df_expr, mat.mut)
        df_expr_1 = cbind(df_expr, mat.mut)
        rownames(df_expr_0) = paste0(rownames(df_expr_0), ".0")
        rownames(df_expr_1) = paste0(rownames(df_expr_1), ".1")
        tabs_[[n]] = rbind(df_expr_0, df_expr_1)
    }
    cliff_input = do.call(rbind, tabs_)
    prop_col = num(cliff_input[, grepl("props", colnames(cliff_input))])
    cliff_input = cliff_input[, !grepl("props", colnames(cliff_input))]
    drug_data = drug_data[rep(seq_len(nrow(drug_data)), each = 2 *
        K), ]
    drug_data$y = drug_data$auc
    drug_data$y_bin = drug_data$auc
    for (i in 0:(2 * N - 1)) {
        start_ = i * K + 1
        end_ = (i + 1) * K
        if (i%%2 == 1) {
            drug_data$y_bin[start_:end_] <- 0
        }
        else {
            drug_data$y_bin[start_:end_] <- 1
        }
    }
    message("Launch CLIFF algorithm and iterate the EM algorithm")
    drug_data$pi_hat = mean_auc
    drug_data$y_hat = 0
    drug_data$q_0 = 0
    drug_data$q_1 = 0
    drug_data$w_0 = 0
    drug_data$w_1 = 0
    prev_rmse = 1e+05
    lambda = 0
    drug_data$y = drug_data$auc
    it.increasing.rmse = -1
    min_rmse = 1000
    for (e in 1:max.em.steps) {
        for (i in 0:(2 * N - 1)) {
            start_ = i * K + 1
            end_ = (i + 1) * K
            div_factor_0 = sum(prop_col[start_:end_] * (1 - drug_data$pi_hat[start_:end_]))
            div_factor_1 = sum(prop_col[start_:end_] * drug_data$pi_hat[start_:end_])
            drug_data$q_0[start_:end_] <- ((prop_col[start_:end_] *
                (1 - drug_data$pi_hat[start_:end_]))/div_factor_0)
            drug_data$q_1[start_:end_] <- ((prop_col[start_:end_] *
                drug_data$pi_hat[start_:end_])/div_factor_1)
            drug_data$w_0[start_:end_] <- (1 - drug_data$y[start_:end_]) *
                drug_data$q_0[start_:end_]
            drug_data$w_1[start_:end_] <- drug_data$y[start_:end_] *
                drug_data$q_1[start_:end_]
        }
        y = drug_data$y_bin
        drug_data$weights = drug_data$w_0
        drug_data$weights[drug_data$y_bin == 1] <- drug_data$w_1[drug_data$y_bin ==
            1]
        w_ = drug_data$weights
        set.seed(1)
        if (regularization == "L2") {
            fit = glmnet(as.matrix(cliff_input), y, weights = w_,
                family = "binomial", intercept = F, scale = F,
                alpha = 0)
        }
        else {
            fit = glmnet(as.matrix(cliff_input), y, weights = w_,
                family = "binomial", intercept = F, scale = F,
                lambda = 0)
        }
        cliff_coefs = coef(fit)[-1, dim(coef(fit))[2]]
        drug_data$pi_hat = sigmoid_(as.matrix(cliff_input) %*%
            cliff_coefs, log((1 - mean_auc)/mean_auc))
        for (i in 0:(2 * N - 1)) {
            start_ = i * K + 1
            end_ = (i + 1) * K
            drug_data$y_hat[start_:end_] <- sum(prop_col[start_:end_] *
                drug_data$pi_hat[start_:end_])
        }
        drug_data.sub = drug_data[drug_data$y_bin == 1, ]
        drug_data.sub$cellType = rep(colnames(climb_prop), N)
        pred_mean_pi = aggregate(drug_data.sub$pi_hat, list(drug_data.sub$cellType),
            mean)
        pi_hat_nk = data.frame(drug_data.sub$pi_hat)
        pi_hat_nk$sample = gsub(".*_", "", rownames(pi_hat_nk))
        pi_hat_nk$cellType = gsub("_.*", "", rownames(pi_hat_nk))
        colnames(pi_hat_nk) = c("pi_hat_nk", "sample", "celltype")
        PI_hat_nk = dcast(pi_hat_nk, sample ~ celltype, value.var = "pi_hat_nk")
        rownames(PI_hat_nk) = PI_hat_nk[, 1]
        PI_hat_nk = PI_hat_nk[, -1]
        rownames(PI_hat_nk) = gsub("\\.0", "", rownames(PI_hat_nk))
        PI_hat_nk = PI_hat_nk[rownames(mutation_data), ]
        auc_values = drug_data.sub$auc[seq(1,length(drug_data.sub$auc),K)]
        rmse_ = rmse(rowSums(as.matrix(PI_hat_nk) * as.matrix(climb_prop)),
            auc_values)
        if ( e > 1 & rmse_ < min_rmse & abs(rmse_ - min_rmse) > 1e-3 ){
            min_rmse = rmse_
            min_coefs = cliff_coefs
            min_PI = PI_hat_nk
            it.increasing.rmse = 0
        }
        else {
            it.increasing.rmse = it.increasing.rmse + 1
            if (it.increasing.rmse >= 3) {
                message(paste0("early stopping of EM algorithm at step ",
                  e))
                break
            }
        }
    }
    return(list(PI_hat_nk, mutation_data, climb_prop))
}

