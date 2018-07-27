#20180727_modifiedvgmfunction.R
#

visualize_gene_markers_modified <- function (gbm, gene_probes, projection, limits = c(0, 10), marker_size = 0.1, 
          title = NULL) 
{
  gbm_trunc <- trunc_gbm_by_genes(gbm, gene_probes)
  gene_values <- t(as.matrix(exprs(gbm_trunc)))
  gene_values[gene_values < limits[1]] <- limits[1]
  gene_values[gene_values > limits[2]] <- limits[2]
  colnames(gene_values) <- gene_probes
  projection_names <- colnames(projection)
  colnames(projection) <- c("Component.1", "Component.2")
  proj_gene <- data.frame(cbind(projection, gene_values))
  proj_gene_melt <- melt(proj_gene, id.vars = c("Component.1", 
                                                "Component.2"))
  p <- ggplot(proj_gene_melt, aes(Component.1, Component.2)) + 
    geom_point(aes(colour = value), size = marker_size) + 
    facet_wrap(~variable, nrow = 2) + scale_colour_gradient(low = "grey", 
                                                  high = "red", name = "val") + labs(x = projection_names[1], 
                                                                                     y = projection_names[2])
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  p <- p + theme_bw() + theme(plot.title = element_text(hjust = 0.5), 
                              panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  return(p)
}

#Neuron vs. glia
genes <-c("ase", "repo", "eya", "hth", "hbn", "Lim1")
quartz("Neuron vs. glia", 8,2)
visualize_gene_markers_modified(gbm_log, genes, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0, 4))