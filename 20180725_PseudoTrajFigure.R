#20180725_PseudoTrajFigure.R
#July 25, 2018

library(monocle)
#must load these for ddply() function to work
library(plyr)
library(reshape2)

#Load CellDataSet for PC = 3
setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180619_Monocle/")
gbm_current <- readRDS(file = "gbm_cds_subset_pseudotraj_PC3.RData")

#Write modified plot_genes_in_pseudotime() function to remove glia from pseudotime.
plot_genes_in_pseudotime_modified<- function (cds_subset, min_expr = NULL, cell_size = 0.75, nrow = NULL, 
          ncol = 1, panel_order = NULL, color_by = "State", trend_formula = "~ sm.ns(Pseudotime, df=3)", 
          label_by_short_name = TRUE, relative_expr = TRUE, vertical_jitter = NULL, 
          horizontal_jitter = NULL) 
{
  f_id <- NA
  Cell <- NA
  if (cds_subset@expressionFamily@vfamily %in% c("negbinomial", 
                                                 "negbinomial.size")) {
    integer_expression <- TRUE
  }
  else {
    integer_expression <- FALSE
    relative_expr <- TRUE
  }
  if (integer_expression) {
    cds_exprs <- exprs(cds_subset)
    if (relative_expr) {
      if (is.null(sizeFactors(cds_subset))) {
        stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
      }
      cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/sizeFactors(cds_subset))
    }
    cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
  }
  else {
    cds_exprs <- reshape2::melt(as.matrix(exprs(cds_subset)))
  }
  if (is.null(min_expr)) {
    min_expr <- cds_subset@lowerDetectionLimit
  }
  colnames(cds_exprs) <- c("f_id", "Cell", "expression")
  cds_pData <- pData(cds_subset)
  cds_fData <- fData(cds_subset)
  cds_exprs <- merge(cds_exprs, cds_fData, by.x = "f_id", by.y = "row.names")
  cds_exprs <- merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")
  if (integer_expression) {
    cds_exprs$adjusted_expression <- cds_exprs$expression
  }
  else {
    cds_exprs$adjusted_expression <- log10(cds_exprs$expression)
  }
  if (label_by_short_name == TRUE) {
    if (is.null(cds_exprs$gene_short_name) == FALSE) {
      cds_exprs$feature_label <- as.character(cds_exprs$gene_short_name)
      cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
    }
    else {
      cds_exprs$feature_label <- cds_exprs$f_id
    }
  }
  else {
    cds_exprs$feature_label <- cds_exprs$f_id
  }
  cds_exprs$f_id <- as.character(cds_exprs$f_id)
  cds_exprs$feature_label <- factor(cds_exprs$feature_label)
  new_data <- data.frame(Pseudotime = pData(cds_subset)$Pseudotime)
  model_expectation <- genSmoothCurves(cds_subset, cores = 1, 
                                       trend_formula = trend_formula, relative_expr = T, new_data = new_data)
  colnames(model_expectation) <- colnames(cds_subset)
  expectation <- ddply(cds_exprs, .(f_id, Cell), function(x) data.frame(expectation = model_expectation[x$f_id, 
                                                                                                        x$Cell]))
  cds_exprs <- merge(cds_exprs, expectation)
  cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
  cds_exprs$expectation[cds_exprs$expectation < min_expr] <- min_expr
  if (is.null(panel_order) == FALSE) {
    cds_exprs$feature_label <- factor(cds_exprs$feature_label, 
                                      levels = panel_order)
  }
  q <- ggplot(aes(Pseudotime, expression), data = cds_exprs)
  if (is.null(color_by) == FALSE) {
    q <- q + geom_point(aes_string(color = color_by), size = I(cell_size), 
                        position = position_jitter(horizontal_jitter, vertical_jitter))
  }
  else {
    q <- q + geom_point(size = I(cell_size), position = position_jitter(horizontal_jitter, 
                                                                        vertical_jitter))
  }
  q <- q + geom_line(aes(x = Pseudotime, y = expectation), 
                     data = cds_exprs)
  q <- q + scale_y_log10() + facet_wrap(~feature_label, nrow = nrow, 
                                        ncol = ncol, scales = "fixed")  ###Change to fixed scale
  if (min_expr < 1) {
    q <- q + expand_limits(y = c(min_expr, 1))
  }
  if (relative_expr) {
    q <- q + ylab("Relative Expression")
  }
  else {
    q <- q + ylab("Absolute Expression")
  }
  q <- q + xlab("Pseudo-time")
  q <- q + scale_x_continuous(limits=c(13, 35)) ##This is the line I added to modify x axis
  q <- q + monocle_theme_opts()
  q
}

###Modified monocle_theme_opts() written for Monocle to manipulate pplot theme
monocle_theme_opts <- function()
{
    theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_line(size=0.25, color="black")) +
    theme(axis.line.y = element_line(size=0.25, color="black")) +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
    theme(panel.background = element_rect(fill='white')) +
    theme(legend.key=element_blank()) +
    theme(strip.text = element_text(size=20)) + ##Added to change label size
    theme(text = element_text(size=20))  ##Added to change text size
}

#plot genes over pseudotime for neurons and progenitors ONLY (pseudotime = 13-35)
my_genes <- row.names(subset(fData(gbm_current),
                             gene_short_name %in% c("ase","nSyb", "Nlg2", "Nlg3",  
                                                  "Syt1", "Syt4", "cpx")))

gbm_current <- gbm_current[my_genes,]
setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180725_PseudoTrajFigure/")
png(file="TrajectoryValidatingGenes_fixed.png", width = 2000, height = 300)
print(plot_genes_in_pseudotime_modified(gbm_current, color_by = "Pseudotime", nrow = 2, ncol = 8, 
                               panel_order = c("ase","nSyb", "Nlg2", "Nlg3",  "Syt1", "Syt4", "cpx")), 
                               cell_size = 1.5)
dev.off()





