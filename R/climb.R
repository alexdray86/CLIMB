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
#' @export
climb <- function(sc, bulk, predict_expression=TRUE, ratio_cancer_cells=NA, up.lim=Inf, lambda=0, norm_factor=0.1){
    ct.props = list() ; ct.exprs = list()
    sc.mat = exprs(sc)
    N = dim(bulk)[2]
    cell_expr = colSums(exprs(sc))
    for(i in 1:N){
        y = num(exprs(bulk)[,i])
        fit = glmnet::glmnet(sc.mat, y, lower.limits = 0.0, upper.limits = up.lim, lambda = lambda)
        coefs = coef(fit)[-1,dim(coef(fit))[2]]
        norm_coefs = coefs / cell_expr # Normalize coefs with total expression per cell
        agg = aggregate(norm_coefs, list(sc$cellType), sum)
        # If information about cancer cell ratio is given, then a smooth normalization is done hereafter
        if(!is.na(sum(num(ratio_cancer_cells)))){
            colnames(agg) = c('celltype', 'sum_coefs')
            # Correct proportions with ratio cancer cells
            b_n = num(ratio_cancer_cells[i]) # True ratio cancer cells
            b_hat_n = sum(agg[grepl("-like",agg$celltype),]$sum_coefs) / sum(agg$sum_coefs) # predicted blast count
            # normalization factors
            f_n = (b_n + norm_factor) / (b_hat_n + norm_factor)
            f_n_n = (1 - b_n + norm_factor) / (1 - b_hat_n + norm_factor)
            agg_norm = agg
            agg_norm[grepl("-like",agg$celltype),]$sum_coefs <- agg_norm[grepl("-like",agg$celltype),]$sum_coefs * f_n
            agg_norm[!grepl("-like",agg$celltype),]$sum_coefs <- agg_norm[!grepl("-like",agg$celltype),]$sum_coefs * f_n_n
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
                pred_expr = (t(coefs)[sel_ct] %*% t(sc.mat[,sel_ct])) / sum(t(coefs)[sel_ct])
                pred_expr[is.na(pred_expr)] = 0
                pred_exprs[[k]] = pred_expr
            }
            ct_exprs_pred = do.call(rbind, pred_exprs)
            ct.exprs[[i]] = ct_exprs_pred
        }
        
    }
    climb.prop = do.call(rbind,ct.props)
    return( list(as.matrix(climb.prop), ct.exprs) )
}
