library(pathfindR)
setwd('~/Documents/desktop-tutorial/Sonja/SAN_GSE30936179/pathway')

### necessary files downloaded in KEGG_signaling_pathways.R
mmu_kegg_genes <- readRDS("mmu_kegg_genes.RDS")
mmu_kegg_descriptions <- readRDS("mmu_kegg_descriptions.RDS")
path2SIF <- file.path('~/Documents/desktop-tutorial/Sonja/SAN_GSE30936179/pathway', "mmusculusPIN.sif")
path2SIF = normalizePath(path2SIF)

getPlot = function(de,mmu_kegg_genes,mmu_kegg_descriptions,path2SIF,color = 'gray'){
  output_df <- run_pathfindR(de,convert2alias = FALSE,
                             gene_sets = "Custom",custom_genes = mmu_kegg_genes,
                             custom_descriptions = mmu_kegg_descriptions,pin_name_path = path2SIF)
  sig = output_df[grep('signaling',output_df$Term_Description),c('Term_Description','highest_p')]
  #row.names(sig) = sig$Term_Description
  sig$log10p = -log10(sig$highest_p)
  sig = sig[order(sig$log10p),]
  par(mar=c(3,15,2,2))
  barplot(sig$log10p,horiz=T,names.arg=sig$Term_Description,las=2,cex.names =0.5,col = color)
  # clustered_df <- cluster_enriched_terms(output_df)
  #kegg <- get_kegg(species = "mmu", path = "~/Documents/desktop-tutorial/Sonja/data/")
  #combine_pathfindR_results
  return(output_df)
}

justPlot = function(output_df,color = 'gray'){
  sig = output_df[grep('signaling',output_df$Term_Description),c('Term_Description','highest_p')]
  #row.names(sig) = sig$Term_Description
  sig$log10p = -log10(sig$highest_p)
  sig = sig[order(sig$log10p),]
  par(mar=c(3,15,2,2))
  barplot(sig$log10p,horiz=T,names.arg=sig$Term_Description,las=2,cex.names =0.5,col = color)
}
load("/Users/zhezhenwang/Documents/desktop-tutorial/Sonja/frCarlos/HTseq_DESeq_LA_padj0.05LFC0.5.rData")
up = subset(resHTdata,log2FoldChange>0.5 & padj<0.05)
uppawy = getPlot(up[,c('gene_id','log2FoldChange','padj')],mmu_kegg_genes,mmu_kegg_descriptions,path2SIF)
save(uppawy,file='uppawy_tbx5_mouse.rData')
justPlot(uppawy,color = 'darkviolet',mmu_kegg_genes,mmu_kegg_descriptions,path2SIF)
dev.copy2pdf(file = 'Tbx5_up_pathway_mouse.pdf')

dn = subset(resHTdata,log2FoldChange< -0.5 & padj<0.05)
dnpawy = getPlot(dn[,c('gene_id','log2FoldChange','padj')],mmu_kegg_genes,mmu_kegg_descriptions,path2SIF)
save(dnpawy,file='dnpawy_tbx5_mouse.rData')
justPlot(dnpawy,color = 'mediumpurple')
dev.copy2pdf(file = 'Tbx5_down_pathway_mouse.pdf')

nodalref = read.table('~/Documents/desktop-tutorial/Sonja/data/nodal_Srivastava/305913r1_online_table_ii.txt',sep = '\t',head = T)
nodalref = nodalref[-1,]
colnames(nodalref)[2:3] = c('log2FC','log2CPM')

up = subset(nodalref,log2FC>0.5 & FDR<0.05)
uppawy = getPlot(up[,c('Gene','log2FC','FDR')],color = 'goldenrod',mmu_kegg_genes,mmu_kegg_descriptions,path2SIF)
save(uppawy,file = 'SAN_up_pathway_mouse.rData')
dev.copy2pdf(file = 'SAN_up_pathway_mouse.pdf')

down = subset(nodalref,log2FC< -0.5 & FDR<0.05)
dnpawy = getPlot(down[,c('Gene','log2FC','FDR')],mmu_kegg_genes,mmu_kegg_descriptions,path2SIF)
justPlot(dnpawy,color = 'khaki')
save(dnpawy,file = 'SAN_down_pathway_mouse.rData')
dev.copy2pdf(file = 'SAN_down_pathway_mouse.pdf')

load("/Users/zhezhenwang/Documents/desktop-tutorial/Sonja/combined/cod/final/DE_tac_lfc0.5padj0.05.rData")
de = cbind(data.frame(gene_id = row.names(DE_tac)),DE_tac[,c('log2FoldChange','padj')]) 
depawy = getPlot(de,color = 'blue',mmu_kegg_genes,mmu_kegg_descriptions,path2SIF)
save(depawy,file='depawy_tac_mouse.rData')
#justPlot(uppawy,color = 'blue')
dev.copy2pdf(file = 'TAC_de_pathway_mouse.pdf')

up = subset(de, log2FoldChange>0)
uppawy = getPlot(up,color = 'darkblue',mmu_kegg_genes,mmu_kegg_descriptions,path2SIF)
save(uppawy,file='uppawy_tac_mouse.rData')
dev.copy2pdf(file = 'TAC_up_pathway_mouse.pdf')

dn = subset(de, log2FoldChange<0)
dnpawy = getPlot(dn,color = 'dodgerblue',mmu_kegg_genes,mmu_kegg_descriptions,path2SIF)
save(dnpawy,file='dnpawy_tac_mouse.rData')
dev.copy2pdf(file = 'TAC_dn_pathway_mouse.pdf')

load("/Users/zhezhenwang/Documents/desktop-tutorial/Sonja/AngII/results/rt.deseq2result.rData")
DE_angii = subset(rt,abs(log2FoldChange)>0.5 & padj<0.05)
de = cbind(data.frame(gene_id = row.names(DE_angii)),DE_angii[,c('log2FoldChange','padj')]) 
depawy = getPlot(de,color = 'red',mmu_kegg_genes,mmu_kegg_descriptions,path2SIF)
save(depawy,file='depawy_angII_mouse.rData')
dev.copy2pdf(file = 'angII_de_pathway_mouse.pdf')

up = subset(de, log2FoldChange>0)
uppawy = getPlot(up,color = 'darkred',mmu_kegg_genes,mmu_kegg_descriptions,path2SIF)
save(uppawy,file='uppawy_angII_mouse.rData')
dev.copy2pdf(file = 'angII_up_pathway_mouse.pdf')

dn = subset(de, log2FoldChange<0)
dnpawy = getPlot(dn,color = 'salmon',mmu_kegg_genes,mmu_kegg_descriptions,path2SIF)
save(dnpawy,file='dnpawy_angII_mouse.rData')
dev.copy2pdf(file = 'angII_dn_pathway_mouse.pdf')



RA_clustered_fuzzy <- cluster_enriched_terms(dnpawy, method = "fuzzy")
term_gene_heatmap(result_df = dnpawy, genes_df = dn)
# Upset
UpSet_plot(result_df = dnpawy, genes_df = dn)
dev.copy2pdf(file = 'AngII_down_pathway.pdf')
## score
score_matrix <- score_terms(enrichment_table = RA_clustered[RA_clustered$Status == "Representative", ],
                            exp_mat = RA_exp_mat,
                            cases = cases,
                            use_description = TRUE, # default FALSE
                            label_samples = FALSE, # default = TRUE
                            case_title = "RA",  # default = "Case"
                            control_title = "Healthy", # default = "Control"
                            low = "#f7797d", # default = "green"
                            mid = "#fffde4", # default = "black"
                            high = "#1f4037")
## combine result
# =combine_pathfindR_results(result_A = , 
#                           result_B = , 
#                           plot_common = T)
# combined_results_graph

load("/Users/zhezhenwang/Documents/desktop-tutorial/Sonja/frCarlos/HTseq_DESeq_LA_padj0.05LFC0.5.rData")
de = subset(resHTdata,abs(log2FoldChange)>0.5 & padj<0.05)
depawy = getPlot(de[,c('gene_id','log2FoldChange','padj')])
