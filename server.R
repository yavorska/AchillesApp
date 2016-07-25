library(shiny)

shinyServer(function(input, output) {

  # Cancer Pie Chart
  output$cancersPieChart <- renderPlot({
    df_can <- data.frame(Cancer = cancers, Number = no_cancers)

    ggplot(df_can, aes(x = "", y = no_cancers, fill = factor(Cancer))) +
      geom_bar(stat = "identity") +
      coord_polar(theta = "y") +
      scale_fill_discrete(name = "Entity") +
      ggtitle("Cell Lines in Achilles") +
      xlab("") + ylab("") +
      theme_bw()
  })

  # List of cell lines in dataset
  output$cancersList <- DT::renderDataTable(data.frame(RNAiObject@cancers), server = F, selection = 'single')

  output$variation.within.cell.line.plot <- renderPlot({
    validate(
      need(input$cancersList_rows_selected, "Choose a cell line.")
    )
    cancerIndex <- input$cancersList_rows_selected
    cancerID <- RNAiObject@cancers[cancerIndex]
    all <- variation.within.cell.line(RNAiObject, cancerID, input$geneName_variation.within.cell.line, "all")
    median <- variation.within.cell.line(RNAiObject, cancerID, input$geneName_variation.within.cell.line, "median")

    # View plots side by side
    multiplot(all, median, cols = 2)
  })

  output$variation.within.entity.plot <- renderPlot({
    geneName <- input$geneName_variation.within.entity
    entityName <- input$entityName_variation.within.entity
    variation.within.entity(RNAiObject, entityName = entityName, geneName = geneName)[[2]]
  })

  output$variation.within.entities.plot <- renderPlot({
    geneName <- input$geneName_variation.within.entity
    variation.within.entities(RNAiObject, geneName, type = "value", statistic = "median")[[2]]
  })

  output$consistency.of.phenotypes.plot <- renderPlot({
    geneName <- input$geneName_consistency.of.phenotypes
    entityName <- input$entityName_consistency.of.phenotypes
    box <- consistency.of.phenotypes(RNAiObject, entityName = entityName, geneName = geneName, type = "box")[[2]]
    line <- consistency.of.phenotypes(RNAiObject, entityName = entityName, geneName = geneName, type = "line")[[2]]
    multiplot(box, line, cols = 1)
  })

  output$consistency.of.phenotypes.all.plot <- renderPlot({
    geneName <- input$geneName_consistency.of.phenotypes
    entityName <- "all"
    box <- consistency.of.phenotypes(RNAiObject, entityName = entityName, geneName = geneName, type = "box")[[2]] + 
            theme(legend.position = "none")
    line <- consistency.of.phenotypes(RNAiObject, entityName = entityName, geneName = geneName, type = "line")[[2]] + 
            theme(legend.position = "none")
    multiplot(box, line, cols = 1)
  })

  output$p53.dependency.by.entity.plot <- renderPlot({
    geneName <- input$geneName_p53.dependency.by.entity
    entityName <- input$entityName_p53.dependency.by.entity
    status <- input$status_p53.dependency.by.entity

    line <- p53.dependency.by.entity(RNAiObject, entityName = entityName, geneName = geneName, p53 = status, type = "line")[[2]]
    box <- p53.dependency.by.entity(RNAiObject, entityName = entityName, geneName = geneName, p53 = status, type = "box")[[2]]
    multiplot(line, box, cols = 1)
  })

  output$p53.dependency.by.entity.all.plot <- renderPlot({
    geneName <- input$geneName_p53.dependency.by.entity
    entityName <- "all"
    status <- input$status_p53.dependency.by.entity

    line <- p53.dependency.by.entity(RNAiObject, entityName = entityName, geneName = geneName, p53 = status, type = "line")[[2]] + 
              theme(legend.position = "none")
    box <- p53.dependency.by.entity(RNAiObject, entityName = entityName, geneName = geneName, p53 = status, type = "box")[[2]]
    multiplot(line, box, cols = 1)
  })

  output$riger.gene.comparison.all.plot <- renderPlot({
    geneName <- input$geneName_GENEE

    if(input$data_GENEE) {
      riger.gene.comparison(Achilles2.0, geneName, "all")
    } else {
      riger.gene.comparison(Achilles2.0, geneName, "all")
    }
  })

  output$riger.gene.comparison.entities.plot <- renderPlot({
    geneName <- input$geneName_GENEE

    if(input$data_GENEE) {

      multiplot(riger.gene.comparison(Achilles2.0, geneName, "colon"),
                riger.gene.comparison(Achilles2.0, geneName, "lung"),
                riger.gene.comparison(Achilles2.0, geneName, "ovary"),
                riger.gene.comparison(Achilles2.0, geneName, "pancreas"),
                cols = 2
                )
    } else {

      multiplot(riger.gene.comparison(Achilles2.4, geneName, "breast"),
                riger.gene.comparison(Achilles2.4, geneName, "cns"),
                riger.gene.comparison(Achilles2.4, geneName, "colon"),
                riger.gene.comparison(Achilles2.4, geneName, "haem"),
                riger.gene.comparison(Achilles2.4, geneName, "lung"),
                riger.gene.comparison(Achilles2.4, geneName, "ovary"),
                riger.gene.comparison(Achilles2.4, geneName, "pancreas"),
                riger.gene.comparison(Achilles2.4, geneName, "skin"),
                cols = 3
                )
    }
  })

  output$riger.gene.comparison.datasets.all.plot <- renderPlot({
    geneName <- input$geneName_GENEE

    riger.gene.comparison.datasets(Achilles2.0, Achilles2.4, "Achilles v2.0", "Achilles v 2.4", geneName, "all", 1, "Comparing RIGER results for all cell lines")

  })

  output$riger.gene.comparison.datasets.entities.plot <- renderPlot({
    geneName <- input$geneName_GENEE

    multiplot(riger.gene.comparison.datasets(Achilles2.0, Achilles2.4, "Achilles v2.0", "Achilles v 2.4", geneName, "colon", 1, "Comparing RIGER results for colon cell lines"),
              riger.gene.comparison.datasets(Achilles2.0, Achilles2.4, "Achilles v2.0", "Achilles v 2.4", geneName, "lung", 1, "Comparing RIGER results for lung cell lines"),
              riger.gene.comparison.datasets(Achilles2.0, Achilles2.4, "Achilles v2.0", "Achilles v 2.4", geneName, "ovary", 1, "Comparing RIGER results for ovary cell lines"),
              riger.gene.comparison.datasets(Achilles2.0, Achilles2.4, "Achilles v2.0", "Achilles v 2.4", geneName, "pancreas", 1, "Comparing RIGER results for pancreas cell lines"),
              cols = 2)
  })

})
