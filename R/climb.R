#' CLIMB deconvolution 
#'
#' By using a scRNA-seq reference dataset (sc) and a bulk target dataset (bulk), CLIMB will generate cell-type abundance
#' and cell-type gene expression for each mixture given in the bulk data. 
#'
#' @param sc ExpressionSet object containing the scRNA-seq reference dataset
#' @param bulk ExpressionSet object containing the mixtures to be deconvoluted
#' @param ratio_cancer_cells numeric vector, same size as number of mixtures. if we have data about the fraction of cancer cells present in mixture, this can improve deconvolution accuracy. If ratio_cancer_cells is provided, then CLIMB-BA method will be used. It is important that cancer cell-types in the scRNA-seq dataset contains the pattern "-like" to be recognized as cancer cells.
#' @param up.lim numeric scalar, impose a l-infinity norm to the linear model (upper bound for coefficient)
#' @param lambda Regularization factor
#' @param norm_factor if ratio_cancer_cells is provided, indicate strength of smooth normalization
#' @param cancer_pattern a string pattern present in all cancer cell-type. Only for these cell-types CLIMB will assume the presence of differentially expressed genes
#' @export
climb <- function(sc, bulk, cancer_pattern = "-like", predict_expression=TRUE, ratio_cancer_cells=NA, up.lim=Inf, lambda=0, norm_factor=0.1){
    num <- function(x){ return(as.numeric(as.character(x)))}
    ct.props = list() ; ct.exprs = list()
    sc.mat = exprs(sc)
    cell_expr = colSums(exprs(sc))
    save_coefs = list() ; save_ncoefs = list()
    cellTypes = levels(sc$cellType)
    N = dim(bulk)[2] ; G = dim(bulk)[1] ; K = length(cellTypes)
    S_pred_mapping_n = array(rep(0,N*G*K), c(N,G,K))
    display('Bulk to single-cell mapping for prediction of cell-type abundance')
    for(i in 1:N){
        y = num(exprs(bulk)[,i])
        fit = glmnet(sc.mat, y, lower.limits = 0.0, upper.limits = up.lim, lambda = lambda)
        coefs = coef(fit)[-1,dim(coef(fit))[2]]
        norm_coefs = coefs / cell_expr # Normalize coefs with total expression per cell
        agg = aggregate(norm_coefs, list(sc$cellType), sum, drop=F)
        # If information about cancer cell ratio is given, then a smooth normalization is done hereafter
        if(!is.na(sum(num(ratio_cancer_cells)))){
            colnames(agg) = c('celltype', 'sum_coefs')
            # Correct proportions with ratio cancer cells
            b_n = num(ratio_cancer_cells[i]) # True ratio cancer cells
            b_hat_n = sum(agg[grepl(cancer_pattern,agg$celltype),]$sum_coefs) / sum(agg$sum_coefs) # predicted blast count
            # normalization factors
            f_n = (b_n + norm_factor) / (b_hat_n + norm_factor)
            f_n_n = (1 - b_n + norm_factor) / (1 - b_hat_n + norm_factor)
            agg_norm = agg
            agg_norm[grepl(cancer_pattern,agg$celltype),]$sum_coefs <- agg_norm[grepl(cancer_pattern,agg$celltype),]$sum_coefs * f_n
            agg_norm[!grepl(cancer_pattern,agg$celltype),]$sum_coefs <- agg_norm[!grepl(cancer_pattern,agg$celltype),]$sum_coefs * f_n_n
            agg_norm$sum_coefs = num(agg_norm$sum_coefs)
            agg_norm = agg_norm[match(levels(sc$cellType), agg_norm$celltype),]
            # Compute proportions and return results
            ppred = (agg_norm$sum_coefs) / sum(agg_norm$sum_coefs)
            ct.props[[i]] = ppred
        } else {
            ppred = (agg$x) / sum(agg$x)      
            names(ppred) = agg$`Group.1`
            ct.props[[i]] = ppred
        }
        if(predict_expression){
            pcor_expr_pred = list()
            all_celltypes = levels(sc$cellType)
            pred_exprs = list()
            for(k in 1:length(all_celltypes)){
                this_ct = all_celltypes[k]
                sel_ct = sc$cellType == this_ct
                pred_expr = (t(coefs)[sel_ct] %*% t(sc.mat[,sel_ct])) #/ sum(t(coefs)[sel_ct])
                pred_expr[is.na(pred_expr)] = 0
                pred_exprs[[k]] = pred_expr
                S_pred_mapping[,,k] = pred_expr
            }
            ct_exprs_pred = do.call(rbind, pred_exprs)
            ct.exprs[[i]] = ct_exprs_pred
        }
        save_coefs[[i]] = coefs
        save_ncoefs[[i]] = norm_coefs
    }
    climb.prop = do.call(rbind,ct.props)
    display('Cell-type abundance done. ')
    if(predict_expression && any(grepl(cancer_pattern, cellTypes))){
        display('Starting high-resolution expression deconvolution')
        normal_sel = !grepl(cancer_pattern,sc$cellType) ; cancer_sel = grepl(cancer_pattern,sc$cellType)
        cancer_ct_sel = grepl(cancer_pattern,cellTypes)
        alpha_overal = do.call(rbind,save_coefs)
        alpha_cancer = do.call(rbind,save_coefs)[,cancer_sel]
        sc.cancer = sc[,cancer_sel]
        C_overal = exprs(sc)
        Y_hat_overal = alpha_overal %*% t(C_overal)
        Y_true_bulk  = t(exprs(bulk))
        Epsilon_ng = Y_true_bulk - Y_hat_overal
        ### Compute Prediction of Cell-type Expression 
        S_pred_n = array(rep(0,N*G*K), c(N,G,K))
        S_pred_mapping_n = array(rep(0,N*G*K), c(N,G,K))
        # iterate over genes
        for(g in 1:G){
            Epsilon_g = num(Epsilon_ng[,g])
            if(sd(Epsilon_g) != 0){
                fit = glmnet(alpha_cancer, Epsilon_g, intercept = TRUE)
                C_diff_cancer = num(coef(fit)[-1,dim(coef(fit))[2]])
                # iterate over bulks
                for(n in 1:N){
                    Epsilon_cancer_n = aggregate(C_diff_cancer*alpha_cancer[n,], list(sc.cancer$cellType), sum)$x
                    Epsilon_n = rep(0, K)
                    Epsilon_n[cancer_ct_sel] = Epsilon_cancer_n
                    S_pred_mapping_n[n,g,] = ct.exprs[[n]][,g]   # cell-type expression from mapping
                    S_pred_n[n,g,] = S_pred_mapping_n[n,g,] + Epsilon_n  # cell-type expression from prediction
                    S_pred_n[n,g,][S_pred_n[n,g,]<0] = 0
                } 
            } else {
                S_pred_n[n,g,] = rep(0,K)
            }
            if( g %% 1000 == 0){
                display(paste0('High-Resolution expression prediction: ', g, ' genes processed...'))
            }
        }
        dimnames(S_pred_n)[[1]] = dimnames(S_pred_mapping_n)[[1]] = colnames(bulk)
        dimnames(S_pred_n)[[2]] = dimnames(S_pred_mapping_n)[[2]] = rownames(bulk)
        dimnames(S_pred_n)[[3]] = dimnames(S_pred_mapping_n)[[3]] = cellTypes
        final_res = list(as.matrix(climb.prop), S_pred_n, S_pred_mapping_n, save_coefs, save_ncoefs)
        names(final_res) = c('props', 'expr.pred', 'expr.mapping', 'coefs', 'coefs.norm')
        return(final_res)
    } else if (predict_expression) {
        dimnames(S_pred_mapping_n)[[1]] = colnames(bulk)
        dimnames(S_pred_mapping_n)[[2]] = rownames(bulk)
        dimnames(S_pred_mapping_n)[[3]] = cellTypes
        final_res = list(as.matrix(climb.prop), S_pred_mapping_n, save_coefs, save_ncoefs)
        names(final_res) = c('props', 'expr.mapping', 'coefs', 'coefs.norm')
        return(final_res)
    } else { 
        final_res = list(as.matrix(climb.prop), save_coefs, save_ncoefs)
        names(final_res) = c('props', 'coefs', 'coefs.norm')
        return(final_res)
    }
}
