library(tidyverse)
library(shiny)
library(DT)

# Load data
Cellline = c("All_filtered","ES","Nas2","iPS")

dataset = "All_filtered"
# Gene names
GeneID_option <- read_csv(paste0('data/',dataset,'/geneID.txt'),col_names = FALSE,col_types=c(col_character()))
colnames(GeneID_option) <- 'Ensembl gene ID|Gene name'

# load all viz found in viz subdirectory
myviz_option <- tibble(path=list.files(paste0('data/',dataset,'/viz'), full.names = TRUE)) %>% 
  mutate(name=basename(path)) %>%
  tidyr::separate(name, 'name', '\\.', extra='drop') %>%   # reorder viz
  select(-path) %>%
  arrange(name)

shinyUI(fluidPage(
  titlePanel('Zahra single-cell RNA-seq - Human embryonic stem cells'),
  sidebarLayout(
    sidebarPanel(
      radioButtons('dataset', 'Cellline', Cellline),
	  radioButtons('viz', 'Visualization', myviz_option$name)
    ),
    mainPanel()
  ),
  
  sidebarLayout(
    sidebarPanel(
	  	selectInput(inputId = "gene", label = "Gene name", choices = GeneID_option)
	),
    mainPanel(
    	plotOutput('plotGene')
    )
  ),
  
  sidebarLayout(
      sidebarPanel(
      radioButtons('cond', 'Cell Condition', c("Sample","Cellline","Timepoint"))
    ),
    mainPanel(
      plotOutput('plotCondition'),
	  DT::dataTableOutput('auctable_cond')
    )
  ),
  
  sidebarLayout(
    sidebarPanel(
      radioButtons('clustering', 'Clustering Algorithm', c('Ward','Louvain')),
   	  conditionalPanel(condition = "input.clustering == 'Ward'",
      	sliderInput(inputId = "n_cluster",
                  label = "Number of clusters", 
                  value = 12, min = 1, max = 30)
	  ),
   	  conditionalPanel(condition = "input.clustering == 'Louvain'",
		radioButtons('knn', 'Nb. of nearest neighbours', c(5,10,20,30,50,100))
	  )
    ),
    mainPanel(
      plotOutput('plotCluster'),
	  DT::dataTableOutput('auctable_clust')
    )
  )
))
