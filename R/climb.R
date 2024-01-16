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
#' @export

climb <- function (sc, bulk, mode = "abundance", 
    up.lim = Inf, lambda = 0, verbose = TRUE, cancer_pattern = "*",
    conditions = NA, final_res = list(), min_common_genes = 100, ratio_cell_increase=0.02, 
    n.iter.subsampling=5, min.n.cells=50, n.top_mean.genes=500 ) 
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
    reformat_strings2 <- function(vector_string){
        # replace plus and minus (e.g. useful for CD34+, CD34- populations)
        vector_string <- gsub('\\-$', '', vector_string) ; vector_string <- gsub('\\+', '', vector_string)
        vector_string <- gsub('minus', '', vector_string) ;        vector_string <- gsub('plus', '', vector_string)
        vector_string <- gsub('\\ ', '\\.', vector_string) ; vector_string <- gsub('[^[:alnum:] ]','',vector_string)
        return(vector_string)
    }
    reformat_celltypes <- function(celltype_labels){
        celltype_labels <- reformat_strings2(as.vector(celltype_labels)) ;   celltype_labels <- factor(celltype_labels)
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
            agg = aggregate(coefs, list(sc$cellType), sum, 
              drop = F)
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
        # Make maximum sizes :
        max_size = round(dim(sc)[2] / 2)
        # sc normalization
        sc_mat = exprs(sc)
        median_counts = median(colSums(sc_mat))
        sc_mat_norm = log(t(t(sc_mat) / colSums(sc_mat))*median_counts + 1)
        K = length(unique(sc$cellType))
        ratio_cell_increase=min(0.01*K, 0.1)
        # bulk normalization
        bulk_mat = exprs(bulk)
        bulk_mat_norm = log(t(t(bulk_mat) / colSums(bulk_mat))*median_counts + 1)
        N = dim(bulk_mat_norm)[2]
        # gene intersection 
        inter.genes = intersect(rownames(sc_mat_norm), rownames(bulk_mat_norm))
        sc_mat_norm = sc_mat_norm[inter.genes,] ; sc = sc[inter.genes,]
        bulk_mat_norm = bulk_mat_norm[inter.genes,] ; bulk = bulk[inter.genes,]
        # compute average expression
        sc_mat_avg = aggregate(t(sc_mat_norm), list(sc$cellType), mean)
        celltypes = levels(sc$cellType)
        celltypes = reformat_celltypes(celltypes)
        rownames(sc_mat_avg) = sc_mat_avg$`Group.1`
        rownames(sc_mat_avg) = reformat_celltypes(rownames(sc_mat_avg))
        sc_mat_avg = sc_mat_avg[,-1]
        sc_mat_avg = sc_mat_avg[as.character(celltypes),]
        sel.genes = apply(as.matrix(sc_mat_avg), 2, sd) != 0
        sc_mat_avg = sc_mat_avg[,sel.genes]
        bulk_mat_norm = bulk_mat_norm[sel.genes,]
        sc_mat_norm = sc_mat_norm[sel.genes,]
        bulk = bulk[sel.genes,]
        sc = sc[sel.genes,]
        all(rownames(sc_mat_avg) == celltypes)
        fc_to_average = t(sc_mat_avg) / colMeans(sc_mat_avg)
        fc_to_second_l = list()
        for(c in 1:ncol(sc_mat_avg)){
            col_c = sc_mat_avg[,c]
            if (col_c[order(-col_c)][2] > 0){
                fc_to_second_l[[c]] = col_c / col_c[order(-col_c)][2]
            } else {
                fc_to_second_l[[c]] = col_c*0.0000001
            }
        }
        fc_to_second = do.call(rbind, fc_to_second_l)
        double_condition_matrix = fc_to_average > num(quantile(flatten(fc_to_average), p=0.95)) & fc_to_second > num(quantile(flatten(fc_to_average), p=0.75))
        double_condition_matrix = double_condition_matrix[rowSums(double_condition_matrix) > 0,]
        celltype_counts = matrix(0, ncol=K, nrow=1)
        list_genes = list()
        r_=1
        for(r in 1:nrow(double_condition_matrix)){
            row_ = double_condition_matrix[r,]
            if (celltype_counts[1,which.max(row_)] <= 100){
                celltype_counts[1,which.max(row_)] = celltype_counts[1,which.max(row_)]+1
                list_genes[[r_]] = rownames(double_condition_matrix)[r]
                r_ = r_+1
            }
        }
        top_var = unique(unlist(list_genes))
        n_top_var = length(top_var)
        # if not even variable genes, limit the number of top genes in bulks
        if( (n_top_var/2) < n.top_mean.genes){
            n.top_mean.genes = round(n_top_var/2)
        }
        sc_mat_norm_v = sc_mat_norm[top_var,]
        bulk_mat_norm_v = bulk_mat_norm[top_var,]
        sc_mat_avg_sub = sc_mat_avg[,top_var]
        max_celltype_size = max(num(table(sc$cellType)))
        if(max_celltype_size > max_size){ max_celltype_size <- max_size }
        min.n.cells = min(min.n.cells, round(0.02*max_size))
        ref.counts = table(sc$cellType)
        prior_proportions = final_res$props.init
        for(q in 1:2){
            ct.props = list()
            for(n in 1:N){
                top_means = -1*sort(-bulk_mat_norm_v[, n]*apply(sc_mat_avg_sub, 2, max))
                top_means = top_means[1:n.top_mean.genes]
                gene_topExpr = rev(names(top_means))

                sc_mat_norm_ = sc_mat_norm_v[gene_topExpr,]
                bulk_mat_norm_ = bulk_mat_norm_v[gene_topExpr,]
                y = num(bulk_mat_norm_[, n])
                
                celltype_counts = matrix(prior_proportions[n,]*(max_size), nrow=1, ncol=K)
                colnames(celltype_counts) = colnames(prior_proportions)
                celltype_counts[celltype_counts < min.n.cells] <- min.n.cells 
                # Run deconvolution on 5 sub-sampled single-cell datasets
                sum_props = matrix(0, nrow=1, ncol=K)
                for(x in 1:n.iter.subsampling){
                    # Sub-sampling cells from single-cell dataset
                    all_samples = list()
                    for(k in 1:length(celltypes)){
                        ct = as.character(celltypes[k])
                        set.seed(x*k)
                        all_samples[[k]] = sample(grep(paste0('^',ct,'$'),sc$cellType), round(num(celltype_counts[,ct])), replace = T)
                    }
                    sample_ = unlist(all_samples)
                    sc.sub = sc[gene_topExpr,sample_]
                    scmat = sc_mat_norm_[,sample_]
                    # Launch deconvolution
                    fit = glmnet(scmat, y, lower.limits = 0, upper.limits = 0.001) 
                    coefs = coef(fit)[-1, dim(coef(fit))[2]]

                    agg = aggregate(coefs, list(sc.sub$cellType), sum, 
                                    drop = F)
                    agg$x[is.na(agg$x)] <- 0
                    ppred = (agg$x)/sum(agg$x)
                    names(ppred) = agg$Group.1
                    sum_props = sum_props+ppred
                }
                pred_prop = sum_props / rowSums(sum_props)
                ct.props[[n]] = pred_prop
            }
            climb_props.corrected = data.frame(do.call(rbind, ct.props))
            rownames(climb_props.corrected) = colnames(bulk)
            colnames(climb_props.corrected) = levels(sc$cellType)
            prior_proportions = climb_props.corrected
        }
        final_res$props.corrected = data.frame(climb_props.corrected)
        minmax <- function(x, na.rm = TRUE) { return((x- min(x)) /(max(x)-min(x))) }
        weights_ = minmax(num(cor(table(sc$cellType), t(final_res$props.init))))
        final_res$props.final = as.matrix(final_res$props.init) * num(weights_) + as.matrix(final_res$props.corrected) * (1 - num(weights_))
        rownames(final_res$props.corrected) = colnames(bulk)
        colnames(final_res$props.corrected) = levels(sc$cellType)
        rownames(final_res$props.final) = colnames(bulk)
        colnames(final_res$props.final) = levels(sc$cellType)
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
            final_res$expr.highres = S_pred_mapping_n
            final_res$expr.mapping = S_pred_mapping_n
            final_res$expr.overall = colSums(S_pred_mapping_n, 
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
            final_res$expr.mapping = S_pred_mapping_n
            final_res$expr.overall = colSums(S_pred_mapping_n, 
                dims = 1)
            final_res$coefs = save_coefs
        }
    }
    else {
        final_res$expr.highres = S_pred_mapping_n
        final_res$expr.mapping = S_pred_mapping_n
        final_res$expr.overall = colSums(S_pred_mapping_n, dims = 1)
        final_res$coefs = save_coefs
    }
    return(final_res)
} 
