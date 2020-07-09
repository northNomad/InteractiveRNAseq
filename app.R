wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)


# load packages -----------------------------------------------------------
suppressPackageStartupMessages({
  library(shiny)
  library(tidyverse)
  library(RColorBrewer)
  library(viridis)
  library(fgsea)
  library(AnnotationDbi)
  library(DESeq2)
  library(patchwork)
  library(gt)
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
  library(stats)
})

# Set up summarizedExperiment ---------------------------------------------
raw <- read_csv("www/readcount_genename.csv")
tpm <- read_csv("www/TPM_genename.csv")
raw <- inner_join(raw[, -21], dplyr::select(tpm, gene_id, gene_chr), by = "gene_id")
feature_data <- raw[, c(1, 20:28)]
raw_s1s2 <- raw[, 1:13]
raw_all <- raw[, 1:19]

feature_gr <- GRanges(
  seqnames = feature_data$gene_chr,
  ranges = IRanges(start = feature_data$gene_start, end = feature_data$gene_end),
  strand = feature_data$gene_strand,
  symbol = feature_data$gene_name,
  biotype = feature_data$gene_biotype,
  description = feature_data$gene_description,
  tf_family = feature_data$tf_family
)
names(feature_gr) <- feature_data$gene_id

# col_data_s1s2 <- data.frame(treatment = rep(c("DMSO", "IVO", "ENTO", "ENTO.IVO", "R406", "R406.IVO"), 2))
# col_data_s1s2$treatment <- factor(col_data_s1s2$treatment, levels = c("DMSO", "IVO", "ENTO", "ENTO.IVO", "R406", "R406.IVO"))
# col_data_s1s2$sample <- factor(rep(c("sample1", "sample2"), each=6))
# rownames(col_data_s1s2) <- colnames(raw_s1s2)[-1]

# col_data_all <- data.frame(treatment = rep(c("DMSO", "IVO", "ENTO", "ENTO.IVO", "R406", "R406.IVO"), 3))
# col_data_all$treatment <- factor(col_data_all$treatment, levels = c("DMSO", "IVO", "ENTO", "ENTO.IVO", "R406", "R406.IVO"))
# col_data_all$sample <- factor(rep(c("sample1", "sample2", "sample3"), each=6))
# rownames(col_data_all) <- colnames(raw_all)[-1]


assay_s1s2 <- raw_s1s2[, -1] %>% as.matrix() 
assay_all <- raw_all[, -1] %>% as.matrix



# dds_s1s2 <- DESeqDataSet(se_s1s2, design = ~treatment)
# dds_s1s2 <- dds_s1s2[rowSums(counts(dds_s1s2)) >= 10, ] #remove low counts
# dds_s1s2 <- DESeq(dds_s1s2)
# 
# dds_all <- DESeqDataSet(se_all, design = ~treatment)
# dds_all <- dds_all[rowSums(counts(dds_all)) >= 10, ]
# dds_all <- DESeq(dds_all)

pathways_all <- list.files("www/MSigDB_mmu_version", pattern = "rds")
treatment <- c("DMSO", "IVO", "ENTO", "R406", "ENTO.IVO", "R406.IVO")

# ui ----------------------------------------------------------------------

ui <- navbarPage(
  tabPanel("GSEA",
    selectInput(inputId = "input_pathway", 
                label = "Select Mouse Version of MSigDB",
                choices = pathways_all
    ),
    sliderInput(inputId = "input_fill_fdr_cutoff",
                label = "Select FDR Cutoff in fGSEA for Colouring Barplots",
                min=0, max=1, step = 0.01, value = .1
                ),
    #setting reference treatment
    selectInput(inputId = "input_reference_treatment",
                label = "Select Reference Samples",
                choices = treatment,
                selected = "DMSO"
    ),
    selectInput(inputId = "input_interest_treatment",
                label = "Select Samples of Interest",
                choices = treatment,
                selected = "IVO"
    ),
    actionButton(inputId = "input_samples",
                 label = "UPDATE SAMPLES"
                 ),
    plotOutput("temp_plot")
  ),
  
  #tab panel for venn
  tabPanel("Venn",
    
  )
  
)

server <- function(input, output){
  interest_treatment <- eventReactive(input$input_samples,{input$input_interest_treatment})
  reference_treatment <- eventReactive(input$input_samples, {input$input_reference_treatment})
  

  #colData for summarizedExperiment
  col_data_s1s2 <- reactive({
    x <- data.frame(treatment = factor(rep(c("DMSO", "IVO", "ENTO", "ENTO.IVO", "R406", "R406.IVO"), 2)))
    x$treatment <- relevel(x$treatment, ref = reference_treatment())
    x$sample <- factor(rep(c("sample1", "sample2"), each=6))
    rownames(x) <- colnames(raw_s1s2)[-1]
    x
  })
  
  col_data_all <- reactive({
    x <- data.frame(treatment = factor(rep(c("DMSO", "IVO", "ENTO", "ENTO.IVO", "R406", "R406.IVO"), 3)))
    x$treatment <- relevel(x$treatment, ref = reference_treatment())
    x$sample <- factor(rep(c("sample1", "sample2", "sample3"), each=6))
    rownames(x) <- colnames(raw_all)[-1]
    x
  })
  #summarizedExperiment
  se_s1s2 <- reactive({
    x <- SummarizedExperiment(
          assays = assay_s1s2,
          rowRanges = feature_gr,
          colData = col_data_s1s2()
         )
    x
  })
    
  se_all <- reactive({
    x <- SummarizedExperiment(
      assays = assay_all,
      rowRanges = feature_gr,
      colData = col_data_all()
    )
    x
  })
  
  #call dds
  dds_s1s2 <- reactive({
    x <- DESeqDataSet(se_s1s2(), design = ~treatment)
    x <- x[rowSums(counts(x)) >= 10, ] #remove low counts
    x <- DESeq(x)
    x
  })

  dds_all <- reactive({
    x <- DESeqDataSet(se_all(), design = ~treatment)
    x <- x[rowSums(counts(x)) >= 10, ] #remove low counts
    x <- DESeq(x)
    x
  })
  

  #shrink lfc
  reslfc_s1s2 <- reactive(
    {
      coef <- paste0("treatment_", interest_treatment(), "_vs_", reference_treatment())
      x <- lfcShrink(dds_s1s2(), coef = coef)
      x <- x[order(x$pvalue), ]
      x$symbol <- mcols(se_s1s2())$symbol[match(rownames(x), rownames(mcols(se_s1s2())))]
      x
    }
  )
  reslfc_all <- reactive(
    {
      coef <- paste0("treatment_", interest_treatment(), "_vs_", reference_treatment())
      x <- lfcShrink(dds_all(), coef = coef)
      x <- x[order(x$pvalue), ]
      x$symbol <- mcols(se_all())$symbol[match(rownames(x), rownames(mcols(se_all())))]
      x
    }
  )
  
  #Obtain vector of log2fold change with ENTREZID of gene as names
  named_vector_ENTREZID_s1s2 <- reactive(
    {
      x <- reslfc_s1s2() %>% as.data.frame() %>% mutate(ENSEMBL = rownames(.))
      y <- AnnotationDbi::select(org.Mm.eg.db,
                                 keys = x$ENSEMBL,
                                 columns = "ENTREZID",
                                 keytype = "ENSEMBL"
                                 )
      x <- right_join(y, x, by = "ENSEMBL") %>% drop_na(ENTREZID)
      z <- x$log2FoldChange
      names(z) <- x$ENTREZID
      z
    }
  )
  named_vector_ENTREZID_all <- reactive(
    {
      x <- reslfc_all() %>% as.data.frame() %>% mutate(ENSEMBL = rownames(.))
      y <- AnnotationDbi::select(org.Mm.eg.db,
                                 keys = x$ENSEMBL,
                                 columns = "ENTREZID",
                                 keytype = "ENSEMBL"
      )
      x <- right_join(y, x, by = "ENSEMBL") %>% drop_na(ENTREZID)
      z <- x$log2FoldChange
      names(z) <- x$ENTREZID
      z
    }
  )
  
  list_pathway <- reactive({readRDS(file.path("./www/MSigDB_mmu_version", input$input_pathway))})
  
  
  temp_res <- reactive({fgsea(pathways = list_pathway(),
                    stats=named_vector_ENTREZID_s1s2(),
                    nperm = 1000
                    )})
  
  #temp-genesets in c7 that contain mecom
  ls <- readRDS("www/MSigDB_mmu_version/Mm.c7.all.v7.1.entrez.rds") 
  
  output$temp_plot <- renderPlot({
    ggplot(temp_res()
           , aes(reorder(pathway, NES), NES)) + 
      geom_col(aes(fill = padj < input$input_fill_fdr_cutoff)) + 
      coord_flip()
  }
  )
  
}


shinyApp(ui=ui, server=server)
