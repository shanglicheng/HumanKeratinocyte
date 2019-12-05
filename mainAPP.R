#BiocManager::install("shiny")
library(shiny)
plotTSNEByType = function(project.list) {
  library(ggplot2)
  g = ggplot(data = project.list$tSNE, mapping = aes(x = tSNE_1, y = tSNE_2, color = project.list$Type)) + 
    scale_colour_manual(name = "Type", values = c("mediumspringgreen", "lightskyblue", "indianred1")) + 
    geom_point(size = 1) + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"))
  return(g)
}
plotTSNEByCluster = function(project.list) {
  library(ggplot2)
  g = ggplot(data = project.list$tSNE, mapping = aes(x = tSNE_1, y = tSNE_2, color = project.list$Cluster)) + 
    scale_colour_manual(name = "Cluster", values = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")) + 
    geom_point(size = 1) + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"))
  return(g)
}
plotMarkerWithShape_skin = function(project.list, marker){
  df_tmp = data.frame(project.list$tSNE[grep('Normal skin', project.list$Type),], 
                      Log2RPKM = log2(as.numeric(project.list$RPKM[marker, grep('Normal skin', project.list$Type)])+1),
                      Cluster = project.list$Cluster[grep('Normal skin', project.list$Type)])
  library(ggplot2)
  g = ggplot(mapping = aes(x = tSNE_1, y = tSNE_2, shape = Cluster), data = df_tmp) + 
    geom_point(mapping = aes (colour = Log2RPKM), size = 1, data = df_tmp) + 
    scale_colour_gradient2(low = "gray", mid = "pink" ,high = "red", midpoint = ((max(df_tmp$Log2RPKM)+min(df_tmp$Log2RPKM))/2),
                           space = "Lab", na.value = "midnightblue", guide = "colourbar",
                           limits=c(min(df_tmp$Log2RPKM), max(df_tmp$Log2RPKM))) + 
    scale_shape_manual(values = 1:(nlevels((df_tmp$Cluster)))) + theme_bw()+ 
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
    xlim(c(min(project.list$tSNE$tSNE_1), max(project.list$tSNE$tSNE_1))) + 
    ylim(c(min(project.list$tSNE$tSNE_2), max(project.list$tSNE$tSNE_2))) + 
    ggtitle(paste0(marker," in normal skin"))
  return(g)
}
plotMarkerWithShape_acute = function(project.list, marker){
  df_tmp = data.frame(project.list$tSNE[grep('Acute wound', project.list$Type),], 
                      Log2RPKM = log2(as.numeric(project.list$RPKM[marker, grep('Acute wound', project.list$Type)])+1),
                      Cluster = project.list$Cluster[grep('Acute wound', project.list$Type)])
  library(ggplot2)
  g = ggplot(mapping = aes(x = tSNE_1, y = tSNE_2, shape = Cluster), data = df_tmp) + 
    geom_point(mapping = aes (colour = Log2RPKM), size = 1, data = df_tmp) + 
    scale_colour_gradient2(low = "gray", mid = "pink" ,high = "red", midpoint = ((max(df_tmp$Log2RPKM)+min(df_tmp$Log2RPKM))/2),
                           space = "Lab", na.value = "midnightblue", guide = "colourbar",
                           limits=c(min(df_tmp$Log2RPKM), max(df_tmp$Log2RPKM))) + 
    scale_shape_manual(values = 1:(nlevels((df_tmp$Cluster)))) + 
    theme_bw()+ 
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black")) + 
    xlim(c(min(project.list$tSNE$tSNE_1), max(project.list$tSNE$tSNE_1))) + 
    ylim(c(min(project.list$tSNE$tSNE_2), max(project.list$tSNE$tSNE_2))) + 
    ggtitle(paste0(marker," in acute wound"))
  return(g)
}
plotMarkerWithShape_pu = function(project.list, marker){
  df_tmp = data.frame(project.list$tSNE[grep('PU', project.list$Type),], 
                      Log2RPKM = log2(as.numeric(project.list$RPKM[marker, grep('PU', project.list$Type)])+1),
                      Cluster = project.list$Cluster[grep('PU', project.list$Type)])
  library(ggplot2)
  g = ggplot(mapping = aes(x = tSNE_1, y = tSNE_2, shape = Cluster), data = df_tmp) + 
    geom_point(mapping = aes (colour = Log2RPKM), size = 1, data = df_tmp) + 
    scale_colour_gradient2(low = "gray", mid = "pink" ,high = "red", midpoint = ((max(df_tmp$Log2RPKM)+min(df_tmp$Log2RPKM))/2),
                           space = "Lab", na.value = "midnightblue", guide = "colourbar",
                           limits=c(min(df_tmp$Log2RPKM), max(df_tmp$Log2RPKM))) + 
    scale_shape_manual(values = 1:(nlevels((df_tmp$Cluster)))) + theme_bw()+ 
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
    xlim(c(min(project.list$tSNE$tSNE_1), max(project.list$tSNE$tSNE_1))) + 
    ylim(c(min(project.list$tSNE$tSNE_2), max(project.list$tSNE$tSNE_2))) + 
    ggtitle(paste0(marker," in PU"))
  return(g)
}
plotMarkerBar_skin = function(project.list, marker){
  od = c(grep('KC_1', project.list$Cluster), 
         grep('KC_2', project.list$Cluster), 
         grep('KC_3', project.list$Cluster),
         grep('KC_4', project.list$Cluster),
         grep('Melanocyte', project.list$Cluster),
         grep('Immune cell', project.list$Cluster))
  rpkm.od = project.list$RPKM[,od]
  type.od = project.list$Type[od]
  cluster.od = project.list$Cluster[od]
  df_plot = data.frame(x = c(1:length(as.numeric(rpkm.od[marker, grep('Normal skin', project.list$Type)]))),
                       RPKM = as.numeric(rpkm.od[marker, grep('Normal skin', project.list$Type)]),
                       Cluster = cluster.od[grep('Normal skin', project.list$Type)])
  
  library(ggplot2)
  g = ggplot(df_plot, mapping = aes(x = x, y = RPKM, fill = Cluster, color = Cluster)) + 
    geom_bar(stat="identity", size=0.01) + 
    scale_color_manual(values = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")) +
    scale_fill_manual(values = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")) + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"), panel.border = element_blank(), 
          panel.background = element_blank()) + 
    ggtitle(paste0(marker," in normal skin"))
  return(g)
}
plotMarkerBar_acute = function(project.list, marker){
  od = c(grep('KC_1', project.list$Cluster), 
         grep('KC_2', project.list$Cluster), 
         grep('KC_3', project.list$Cluster),
         grep('KC_4', project.list$Cluster),
         grep('Melanocyte', project.list$Cluster),
         grep('Immune cell', project.list$Cluster))
  rpkm.od = project.list$RPKM[,od]
  type.od = project.list$Type[od]
  cluster.od = project.list$Cluster[od]
  df_plot = data.frame(x = c(1:length(as.numeric(rpkm.od[marker, grep('Acute wound', project.list$Type)]))),
                       RPKM = as.numeric(rpkm.od[marker, grep('Acute wound', project.list$Type)]),
                       Cluster = cluster.od[grep('Acute wound', project.list$Type)])
  
  library(ggplot2)
  g = ggplot(df_plot, mapping = aes(x = x, y = RPKM, fill = Cluster, color = Cluster)) + 
    geom_bar(stat="identity", size=0.01) + 
    scale_color_manual(values = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")) +
    scale_fill_manual(values = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")) + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"), 
          panel.border = element_blank(), 
          panel.background = element_blank()) + 
    ggtitle(paste0(marker," in normal skin"))
  return(g)
}
plotMarkerBar_pu = function(project.list, marker){
  od = c(grep('KC_1', project.list$Cluster), 
         grep('KC_2', project.list$Cluster), 
         grep('KC_3', project.list$Cluster),
         grep('KC_4', project.list$Cluster),
         grep('Melanocyte', project.list$Cluster),
         grep('Immune cell', project.list$Cluster))
  rpkm.od = project.list$RPKM[,od]
  type.od = project.list$Type[od]
  cluster.od = project.list$Cluster[od]
  df_plot = data.frame(x = c(1:length(as.numeric(rpkm.od[marker, grep('PU', project.list$Type)]))),
                       RPKM = as.numeric(rpkm.od[marker, grep('PU', project.list$Type)]),
                       Cluster = cluster.od[grep('PU', project.list$Type)])
  
  library(ggplot2)
  g = ggplot(df_plot, mapping = aes(x = x, y = RPKM, fill = Cluster, color = Cluster)) + 
    geom_bar(stat="identity", size=0.01) + 
    scale_color_manual(values = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")) +
    scale_fill_manual(values = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")) + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"), 
          panel.border = element_blank(), 
          panel.background = element_blank()) + 
    ggtitle(paste0(marker," in normal skin"))
  return(g)
}

ui = shinyUI(
  fluidPage(
    ###
    #Title
    titlePanel("Expressions of genes in human skin cells"),
    
    #Part 1 
    sidebarLayout(position = "left",
                  sidebarPanel("",
                               textInput(inputId = "GeneName", label = "Gene symbol", value = "KRT10", width = "100px")
                  ),
                  mainPanel()
    ),
    mainPanel(
      tabsetPanel(
        type = 'tabs',
        tabPanel('Cells in tSNE',
                 fluidRow('',
                          column(width = 4, plotOutput('tsne.cluster', width = 390, height = 300)),
                          column(width = 4, plotOutput('tsne.type', width = 390, height = 300))
                 ),
                 fluidRow('',
                          column(width = 4, plotOutput('skin', width = 390, height = 300)),
                          column(width = 4, plotOutput('acute', width = 390, height = 300)),
                          column(width = 4, plotOutput('pu', width = 390, height = 300))
                 ),
                 fluidRow('',
                          column(width = 4, plotOutput('exp.skin', width = 300, height = 300)),
                          column(width = 4, plotOutput('exp.acute', width = 300, height = 300)),
                          column(width = 4, plotOutput('exp.pu', width = 300, height = 300))
                 )
        )
      )
    )
  )
)

server <- shinyServer(
  function(input, output) {
    ###
    #load data
    load("./data/data.Rd")
    
    output$tsne.type = renderPlot({
      plotTSNEByType(project.list = project.list)
    })
    
    output$tsne.cluster <- renderPlot({
      plotTSNEByCluster(project.list = project.list)
    })
    
    output$skin <- renderPlot({
      plotMarkerWithShape_skin(project.list = project.list, marker = input$GeneName)
    })
    
    output$acute <- renderPlot({
      plotMarkerWithShape_acute(project.list = project.list, marker = input$GeneName)
    })
    
    output$pu <- renderPlot({
      plotMarkerWithShape_pu(project.list = project.list, marker = input$GeneName)
    })
    
    output$exp.skin <- renderPlot({
      plotMarkerBar_skin(project.list = project.list, marker = input$GeneName)
    })
    
    output$exp.acute <- renderPlot({
      plotMarkerBar_acute(project.list = project.list, marker = input$GeneName)
    })
    
    output$exp.pu <- renderPlot({
      plotMarkerBar_pu(project.list = project.list, marker = input$GeneName)
    })
  }
)

shinyApp(ui = ui, server = server)
