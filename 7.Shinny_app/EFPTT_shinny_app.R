library(shinydashboard)
library(shiny)
library(data.table)
library(fastDummies)
library(ggplot2)
library(multiROC)
library(caret)
library(pROC)
library(GGally)
library(glmnet)
library(randomForest)
library(reshape2)
library(DT)
library(plotly)
library(reactable)
library(ggpubr)
library(gghalves)

set.seed(123)

# create path
path <- "/home/minhhoang/Documents/Thyroid cancer/TOTAL/Shinny_app"

# input data 
traindata <- as.data.frame(fread(paste0(path,"/Data/traindata_model.tsv",header=T))
data <- rbind(Plot = NA, traindata)


ui <- dashboardPage(
                    dashboardHeader(title = "Basic dashboard"),
                    dashboardSidebar(
                    sidebarMenu(
                    menuItem("Training data", tabName = "training_data", icon = icon("table")),
                    menuItem("Classification", icon = icon("th"), tabName = "classification"),
      
                    fileInput("file1", "Choose CSV File",
                             multiple = TRUE,
                             accept = c("text/csv","text/comma-separated-values,text/plain",".csv")))
                  ),

      dashboardBody(
                tabItems(
                tabItem(tabName = "training_data",
                    reactable(
                        data,
                        defaultColDef = colDef(cell = function(value, index, name) {
                            # Render a "show plot" button in the first row - this button won't
                            # do anything by itself, but it will trigger the custom cell click action
                            # when clicked.
                            if (index == 1 && name %in% names(traindata)) {
                              htmltools::tags$button(shiny::icon("line-chart"), `aria-label` = "Show plot")
                            } else {
                              value
                            }
                          }),
                          # Custom cell click action on the "show plot" row to send an event to Shiny.
                          # input$plot_column will be used to access the clicked column name/ID.
                          onClick = JS("
                function(rowInfo, colInfo) {
                  if (rowInfo.index === 0) {
                    if (window.Shiny) {
                      Shiny.setInputValue('plot_column', colInfo.id, { priority: 'event' })
                    }
                  }
                }
              ")
                        ),
             
                plotlyOutput("graph1")
              ),
      
      tabItem(tabName = "classification",
          column(width = 12, 
            DTOutput('trace_table', height = 400)),
              fluidRow(
                box(plotlyOutput("graph")),
                box(title = "Case Analyses Details", height = 400,
                    DTOutput("plot2")
                ))
      )
    )
    # Boxes need to be put in a row (or column)
  )
)

server <- function(input, output) {
  
  output$graph1 <- renderPlotly({
    req(input$plot_column %in% names(traindata))
    traindata2<- traindata[,c(input$plot_column,"subtype")]
    mtraindata2<- melt(traindata2)
    print(mtraindata2)
    mtraindata2 %>%
      plot_ly() %>% 
      add_trace(x = ~subtype,y = ~value, color = ~subtype, type = "box", 
                hoverinfo = 'name+y') %>%
      add_markers(x = ~subtype, y = ~value, color = ~subtype,
                  marker = list(size = 10),
                  hoverinfo = "text",
                  text = ~paste0("Group: ",subtype,
                                 "<br>xval: ",value),
                  showlegend = FALSE) %>% 
      
      layout(legend = list(orientation = "h",
                           x =0.5, xanchor = "center",
                           y = 1, yanchor = "bottom"
      ),
      # xaxis = list(title = "Group",
      #              showticklabels = FALSE))
      xaxis = list(title =  paste("Plot for", input$plot_column),
                   showticklabels = FALSE))
    
  })
  set.seed(122)
  histdata <- rnorm(500)
  
  testdata <- reactive({
    inFile <- input$file1
    if (is.null(inFile)) return(NULL)     
    else read.csv(input$file1$datapath,
                  header = TRUE,
                  sep = "\t")
  })
  
  dotplot1 <- reactive({
  testdata <- testdata()
  traindata <- as.data.frame(fread(paste0(path,"/Data/traindata_model.tsv",header=T))
    
    # creat model
    set.seed(123)
    traindata_rf <- randomForest(as.factor(subtype)~., data=traindata, ntree=1000, proximity=TRUE)
    methylclass <- predict(traindata_rf, newdata=traindata)
    rf_score <- as.data.frame(predict(traindata_rf, newdata=traindata, type="prob"))
    x <- as.matrix(rf_score) 
    y <- as.factor(methylclass)  
    cvfit <- cv.glmnet(x, y, family = "multinomial", type.multinomial = "grouped")
    testPred <- as.matrix(predict(traindata_rf, newdata=testdata, type="prob"))
    predic_new <- as.factor(predict(cvfit, newx = testPred, s = "lambda.min", type = "class"))
    predic_new_score  <- predict(cvfit, newx = testPred, type='response')
    
    
    # creat data frame 
    predict <- as.data.frame(predic_new_score)
    colnames(predict) <- c("FA","FC","fvPTC", "NIFTP")
    predict$Predict <- names(predict)[1:4][apply(predict[,1:4], 1, which.max)]
    new_data <-predict %>% mutate(highest_degree = pmax(FA, FC, fvPTC,NIFTP,na.rm = TRUE))
    # draw 
    library(ggplot2)
    new_data$number <- seq(1,nrow(testdata),1)
    h <- rep(new_data$predict,each=4)
    h <- rep(paste0(seq(1,nrow(testdata),1),"_", new_data$Predict),each=4)
    c <- t(new_data[,c(1:4)])
    library(reshape2)
    mc <- melt(c)
    dotplot1 <- as.data.frame(cbind(mc,h))
    dotplot1 <- dotplot1[,-2]
    colnames(dotplot1) <- c("Score_of_Classes", "Scores", "Predict")
    samples <- paste0(seq(1,nrow(testdata),1),"_", predict$Predict)
    dotplot1$Predict<- factor(dotplot1$Predict,levels=samples)
    list(new_data=new_data,dotplot1=dotplot1)
  })
  
  
  output$graph <- renderPlotly({
    req(input$file1)
    dotplot1 <- dotplot1()$dotplot1
    print(dotplot1)
    plot_ly(
      data = dotplot1,
      x = ~Predict,
      y = ~Scores,
      marker = list(size = 10),
      type = "scatter",
      mode = "markers",
      color = ~Score_of_Classes
    )%>%
      layout(barmode="overlay",
             title = "Scores of Calibration model",
             xaxis = list(tickangle = 90),
             legend = list(orientation = "h",  
                           xanchor = "center",
                           x = 0.5, y=- 0.4))
    
  })
  
  
  output$plot2 <- renderDataTable({
    req(input$file1)
    new_data <-dotplot1()$new_data
    datatable(
      new_data [,1:5],
      options = list(
        scrollX = TRUE,
        scrollY = "250px"
      )        
    )
  })
  
  
  output$trace_table <- renderDataTable({
    req(input$file1)
    
    datatable(
      testdata(),
      options = list(
        scrollX = TRUE,
        scrollY = "250px"
      )        
    )
  })
}



shinyApp(ui, server)