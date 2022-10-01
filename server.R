library(tidyverse)
library(shiny)
library(DT)

clust_colors <- readLines('data/cluster_colors.txt')
time_colors <- readLines('data/time_colors.txt')
cond_colors <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

Cellline <- c("ES","Nas2","iPS")

# declare data lists
myltq <- list()
mycells <- list()
Louvain <- list()
Ward <- list()
myviz <- list()
myauc <- list()

# Load data for All_filtered
dataset = "All_filtered"

# Gene expression
load(paste0("data/",dataset,"/log_transcription_quotients.RData"))
myltq[[dataset]] <- mydata
rm(mydata)

# cell subgroups
mycells[[dataset]] <- myltq[[dataset]] %>% colnames() %>% .[-1] %>%
	tibble(id=.) %>% 
	tidyr::separate(id, c("Batch","Cellline","Timepoint"), extra='drop', remove=FALSE) %>% 
	group_by(id)
mycells[[dataset]]$Sample = paste(paste(mycells[[dataset]]$Batch,mycells[[dataset]]$Cellline),mycells[[dataset]]$Timepoint)

# Get Clusterings
Louvain[[dataset]] <- read_csv(paste0('data/',dataset,'/Clustering/Louvain.txt'),col_types=cols(.default=col_integer())) %>%
	mutate(Cellline=mycells[[dataset]]$Cellline) %>%
	mutate(id=mycells[[dataset]]$id) %>%
	group_by(id)
Ward[[dataset]] <- read_csv(paste0('data/',dataset,'/Clustering/Ward.txt'), col_names=paste0('k', 1:30), col_types=cols(.default=col_integer())) %>%
	mutate(Cellline=mycells[[dataset]]$Cellline) %>%
	mutate(id=mycells[[dataset]]$id) %>%
	group_by(id)
	
# load all viz found in viz subdirectory
myviz[[dataset]] <- tibble(path=list.files(paste0('data/',dataset,'/viz'), full.names = TRUE)) %>%   # extract viz name
	mutate(name=basename(path)) %>% 
	tidyr::separate(name, 'name', '\\.', extra='drop') %>%   # reorder viz
	mutate(space=purrr::map(path, ~ read_csv(..1, col_names=c('x', 'y'), col_types=cols(.default=col_double())) %>% 
	mutate(cell=mycells[[dataset]]$id) %>% mutate(Cellline=mycells[[dataset]]$Cellline) )) %>% 
	select(-path) %>% 
	arrange(name)

# Get data in cellline subgroups
for(cellline in Cellline){
	# filter ltq matrix
	myltq[[cellline]] <- myltq[[dataset]][,c(TRUE,mycells[[dataset]]$Cellline==cellline)]
	# filter mycells
	mycells[[cellline]] <- mycells[[dataset]] %>% filter(Cellline==cellline)
	# filter clusterings
	Louvain[[cellline]] <- Louvain[[dataset]] %>% filter(Cellline==cellline)
	Ward[[cellline]] <- Ward[[dataset]] %>% filter(Cellline==cellline)
}

for(cellline in Cellline){
	# load all viz found in viz subdirectory and filter by cellline
	myviz[[cellline]] <- tibble(path=list.files(paste0('data/',dataset,'/viz'), full.names = TRUE)) %>%   # extract viz name
		mutate(name=basename(path)) %>% 
		tidyr::separate(name, 'name', '\\.', extra='drop') %>%   # reorder viz
		mutate(space=purrr::map(path, ~ read_csv(..1, col_names=c('x', 'y'), col_types=cols(.default=col_double())) %>% 
		mutate(cell=mycells[[dataset]]$id) %>% mutate(Cellline=mycells[[dataset]]$Cellline) %>% filter(Cellline==cellline) )) %>% 
		select(-path) %>% 
		arrange(name)
}


# load all AUC found in AUC subdirectory
for(cellline in c(dataset,Cellline)){
	myauc[[cellline]] <- tibble(path=list.files(paste0('data/',cellline,'/AUC'), full.names = TRUE)) %>% #extract names
		mutate(name=basename(path)) %>%
		tidyr::separate(name, 'name', '\\.', extra='drop') %>%
		mutate(mytable=purrr::map(path, ~ read_csv(..1, col_types=cols("Row"=col_character(),.default=col_double())) )) %>%
		select(-path)
}

shinyServer(function(input, output, session) {
 
	shiny_dataset <- reactive({
		input$dataset
	})

	shiny_viz <- reactive({
		input$viz
	})

  shiny_gene <- reactive({
		input$gene
  })

  shiny_cond <- reactive({
    input$cond
  })

  shiny_clustering <- reactive({
    input$clustering
  })

  shiny_n_cluster <- reactive({
	  input$n_cluster
  })

  shiny_knn <- reactive({
	  input$knn
  })

  output$plotGene <- renderPlot({
    myviz[[shiny_dataset()]] %>% filter(name==shiny_viz()) %>% 
      unnest(space) %>% 
      left_join(
        myltq[[shiny_dataset()]] %>% 
          filter(GeneID==shiny_gene()) %>%
          tidyr::pivot_longer(-GeneID, names_to = 'cell')
      ) %>% 
      ggplot() +
      geom_point(aes(x, y, col=value), alpha=.5, stroke=0, size=1) +
      #scale_color_viridis_c(option = "plasma") +
	  scale_colour_gradient(low="grey90",high="red") +
      labs(col=str_split(shiny_gene(), '\\|')[[1]][2]) +
      theme_void() +
      theme(legend.position = 'bottom') +
      annotate('point', 0, 0, pch=3, size=4) +
	  coord_fixed(ratio = 1/2) +
      NULL
  })
  
  output$plotCondition <- renderPlot({
	if (shiny_cond() == 'Timepoint'){
	  	my_colors = time_colors
	}else{
		my_colors = cond_colors 
	}
    myviz[[shiny_dataset()]] %>% filter(name==shiny_viz()) %>% 
      unnest(space) %>% 
      left_join(
        mycells[[shiny_dataset()]] %>% select(cell=id, cond=shiny_cond())
      ) %>% 
      ggplot() +
      geom_point(aes(x, y, col=factor(cond)), alpha=.5, stroke=0, size=1) +
      scale_color_manual(values=my_colors, guide=guide_legend(override.aes=list(alpha=1, size=4))) +
      labs(col=shiny_cond()) +
      theme_void() +
      theme(legend.position = 'bottom') +
      annotate('point', 0, 0, pch=3, size=4) +
	  coord_fixed(ratio = 1/2) +
      NULL
  })
  
  output$plotCluster <- renderPlot({
	  if (shiny_clustering()=='Louvain'){
		my_colors = clust_colors[ sort(unique( Louvain[[shiny_dataset()]][[paste0('knn',shiny_knn())]] )) ]
		myviz[[shiny_dataset()]] %>% filter(name==shiny_viz()) %>% 
		  unnest(space) %>% 
		  left_join(
			Louvain[[shiny_dataset()]] %>% select(cell=id, cl=paste0('knn',shiny_knn()) )
		  ) %>% 
		  ggplot() +
		  geom_point(aes(x, y, col=factor(cl)), alpha=.5, stroke=0, size=1) +
		  scale_color_manual(values=my_colors, guide=guide_legend(override.aes=list(alpha=1, size=4))) +
		  labs(col='clusters') +
		  theme_void() +
		  theme(legend.position = 'bottom') +
		  annotate('point', 0, 0, pch=3, size=4) +
		  coord_fixed(ratio = 1/2) +
		  NULL
	  }else if (shiny_clustering()=='Ward'){
		my_colors = clust_colors[ sort(unique( Ward[[shiny_dataset()]][[paste0('k',shiny_n_cluster())]] )) ]
		myviz[[shiny_dataset()]] %>% filter(name==shiny_viz()) %>% 
		  unnest(space) %>% 
		  left_join(
			 Ward[[shiny_dataset()]] %>% select(cell=id, cl=paste0('k',shiny_n_cluster()) )
		  ) %>% 
		  ggplot() +
		  geom_point(aes(x, y, col=factor(cl)), alpha=.5, stroke=0, size=1) +
		  scale_color_manual(values=my_colors, guide=guide_legend(override.aes=list(alpha=1, size=4))) +
		  labs(col='clusters') +
		  theme_void() +
		  theme(legend.position = 'bottom') +
		  annotate('point', 0, 0, pch=3, size=4) +
		  coord_fixed(ratio = 1/2) +
		  NULL
	  }
  })

    output$auctable_cond <- DT::renderDataTable({
		myauc[[shiny_dataset()]] %>% filter(name==shiny_cond()) %>% unnest(mytable) %>% select(-name) %>% column_to_rownames("Row")
	})

    output$auctable_clust <- DT::renderDataTable({
		if (shiny_clustering()=='Louvain'){	
			myauc[[shiny_dataset()]] %>% filter(name==paste0('Louvain_knn',shiny_knn())) %>% unnest(mytable) %>% select(-name) %>% column_to_rownames("Row")
		}else if (shiny_clustering()=='Ward'){
			myauc[[shiny_dataset()]] %>% filter(name=='Ward_4') %>% unnest(mytable) %>% select(-name) %>% column_to_rownames("Row")
		}
	}) 
})

