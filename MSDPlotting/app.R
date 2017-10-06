#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library('shiny')
library('dplyr')
library('ggplot2')
library('RColorBrewer')
library('DT')
library('viridis')
library('shinythemes')

filename <- "Simulated1"
deploy <- TRUE

if (deploy) {
  load(paste0(filename,"_shiny.Rdata"))
} else {
  setwd('..')
  load(paste0("Output/",filename,"_shiny.Rdata"))
}
w <- 1 #Deff distance in um
data <- data %>% mutate(Deff = w^2/4*(4*D/(alpha*w^2))^(1/alpha))
data.numlength <- data %>% select(ID,SizeClass,lengthtime,D,alpha) %>% unique()
data.display <- data %>%  mutate(type = dplyr::if_else(SizeClass == "S","granule","scrum")) %>% 
  mutate(alpha = round(alpha, digits = 4), D = round(D, digits = 6), lengthtime = round(lengthtime, digits = 2)) %>%
  select(ID, type, alpha, gamma = D, length = lengthtime, uniqueID) %>% unique()

#data.realalpha <- data %>% filter(alpha >= 0) %>% filter(alpha <= 2)
deltarange <- c(0.1,150)
msdrange <- c(-6,6)
trajrange <- c(0.2,300)
plotbreaks <- c(10^-6,10^-4,10^-2,10^0,10^2,10^4,10^6)
plotlabels <- c(expression("10"^{"-6"}), expression("10"^{"-4"}), expression("10"^{"-2"}), expression("10"^{"0"}),
                expression("10"^{"2"}), expression("10"^{"4"}), expression("10"^{"6"}))
mypalette <- brewer.pal(3,"Set1")
deffrange <- c(floor(log10(min(data.numlength$D, na.rm = TRUE))), ceiling(log10(max(data.numlength$D, na.rm = TRUE))))

theme.dist <- theme_bw() + theme(panel.grid.major.x = element_line(color = 'grey75', size = 0.3),
                                 panel.grid.minor.x = element_line(color = 'grey85', size = 0.25), 
                                 panel.grid.major.y = element_line(color = 'grey75', size = 0.3), 
                                 panel.grid.minor.y = element_line(color = 'grey85', size = 0.25), 
                                 panel.border = element_rect(fill = NA, color = 'black', size = 0.75),
                                 legend.key = element_blank())
theme.traj <- theme(panel.background = element_rect(fill = 'black'), line = element_blank(), axis.title = element_blank(), axis.text = element_blank(), 
                    strip.background = element_blank(), strip.text = element_blank(), legend.position = "bottom")
scale <- 500 # in nm
scalelength <- scale / pixelsize
overlaycolor <- 'forestgreen' #'chartreuse'

ui <- navbarPage(
  title = "Trajectory Analysis", selected = "MSD", theme = shinytheme("cerulean"), id ="navbar",
  tabPanel("MSD",
           withMathJax(),
           fluidRow(
             column(7,
                    h4(paste0("Filename: ",filename)),
                    "Equations:   ",
                    tags$ul(
                      tags$li("\\( MSD = 4 \\Gamma \\Delta^{\\alpha} \\)  ")
                      #tags$li("\\( D_{eff} \\equiv D(t = \\tau) = \\frac{\\omega^{2}}{4}\\Big(\\frac{4\\Gamma}{\\omega^{2}}\\Big)^{1/\\alpha}  \\),",
                      #        " where \\( \\omega^{2} = \\) ",w,"\\( \\mu m^{2} \\)")
                    ),
                    h4(textOutput("numTraj1")),
                    h4(textOutput("numTraj2")),
                    tabsetPanel(
                      tabPanel("Data Display",
                               column(6,
                                      sliderInput("traj_length", label = "Trajectory Length [s]",
                                                  min = trajrange[1], max = trajrange[2], value = c(100,200), step = 0.1, round = -1),
                                      selectInput("type", label = "Type:",
                                                  choices = c("All","Granules","Scrums"))
                               ),
                               column(6,
                                      sliderInput("delta_adjust", label = "\\( \\Delta \\) Range [s]",
                                                  min = deltarange[1], max = deltarange[2], value = 10, step = 0.1, round = -1),
                                      sliderInput("msd_adjust", label = "MSD Range [\\(\\mu m^2\\)] (Powers of 10)",
                                                  min = msdrange[1], max = msdrange[2], value = c(-4,1), step = 0.1,
                                                  round = -1)
                               )
                      ),
                      tabPanel("Style",
                               column(6,
                                      sliderInput("fontsize_adjust", label = "Font Size",
                                                  min = 8, max = 72, value = 20, step = 1, round = TRUE),
                                      sliderInput("alpha_adjust", label = "Transparancy",
                                                  min = 0, max = 1, value = 0.8, step = 0.01, round = -2)
                               ),
                               column(6,
                                      sliderInput("line_adjust", label = "Line Width",
                                                  min = 0.25, max = 5, value = 0.5, step = 0.25),
                                      sliderInput("point_adjust", label = "Point Size",
                                                  min = 0.25, max = 5, value = 0.5, step = 0.25)
                               )
                      )
                    )
                    
             ),
             column(5,
                    tabsetPanel(
                      tabPanel("CDF", plotOutput("cdfPlot") ),
                      tabPanel("Image", img(src = paste0(filename,".png"), height = "100%", width = "100%", align = "middle") )
                    )
             )
           ),
           fluidRow(
             column(8,
                    tabsetPanel(
                      tabPanel("MSD: Type", plotOutput("msdPlot.gs", brush = "msd.brush", click = "msd.click")),
                      tabPanel("MSD: \\( \\alpha \\)", plotOutput("msdPlot.alpha", brush = "msd.brush", click = "msd.click")),
                      tabPanel("MSD: \\( \\Gamma \\)", plotOutput("msdPlot.Deff", brush = "msd.brush", click = "msd.click")),
                      tabPanel("Distribution: \\( \\alpha \\)",
                               plotOutput("Dist.alpha"),
                               column(width = 6, offset = 3, sliderInput("alpha_bin", label = "\\( \\alpha \\) Bin Size",
                                                                         min = 0.01, max = 0.5, value = 0.1, step = 0.01, round = -2) )
                      ),
                      tabPanel("Distribution: \\( \\Gamma \\)",
                               plotOutput("Dist.deff"),
                               column(width = 6, offset = 3, sliderInput("deff_bin", label = "\\( \\Gamma \\) Bin Size",
                                                                         min = 0.01, max = 1, value = 0.1, step = 0.01, round = -2) )
                      )
                    )
             ),
             column(4,
                    tabsetPanel(
                      tabPanel("Select Plot", plotOutput("select.traj"),
                               "Scale bar: ",scale," nm"
                      ),
                      tabPanel("Select Table", tableOutput("select.info") )
                    )
             )
           ),
           fluidRow(
             column(8,
                    DT::dataTableOutput('display.table')
             ),
             column(4,
                    plotOutput("table.traj"),
                    "Scale bar: ",scale," nm"
             )
           )
  )
  #tabPanel(title = "Compare"),
  #tabPanel(title = "Data Sets"),
  #tabPanel(title = "Close", value = "close")
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  if(v < 0.2) {
    stopApp()
  }
  
  observe({
    if (input$navbar == "close") {
      stopApp() 
    }
  })   
  
  
  data.filter <- reactive(
    if(input$type == "Granules") {
      filter(data,lengthtime >= as.numeric(input$traj_length[1])) %>% filter(lengthtime <= as.numeric(input$traj_length[2])) %>% 
        filter(SizeClass == "S")
    } else if(input$type == "Scrums") {
      filter(data,lengthtime >= as.numeric(input$traj_length[1])) %>% filter(lengthtime <= as.numeric(input$traj_length[2])) %>% 
        filter(SizeClass == "L")
    } else {
      filter(data,lengthtime >= as.numeric(input$traj_length[1])) %>% filter(lengthtime <= as.numeric(input$traj_length[2]))
    }
  )
  
  data.numlength.filter <- reactive(
    filter(data.numlength, lengthtime >= as.numeric(input$traj_length[1])) %>% filter(lengthtime <= as.numeric(input$traj_length[2]))
  )
  
  data.ID.filter <- reactive(
    if(input$type == "Granules") {
      filter(data.ID,lengthtime >= as.numeric(input$traj_length[1])) %>% filter(lengthtime <= as.numeric(input$traj_length[2])) %>% 
        filter(SizeClass == "S")
    } else if(input$type == "Scrums") {
      filter(data.ID,lengthtime >= as.numeric(input$traj_length[1])) %>% filter(lengthtime <= as.numeric(input$traj_length[2])) %>% 
        filter(SizeClass == "L")
    } else {
      filter(data.ID,lengthtime >= as.numeric(input$traj_length[1])) %>% filter(lengthtime <= as.numeric(input$traj_length[2]))
    }
  )
  
  output$numTraj1 <- renderText(
    paste0(
      "Displaying ", prettyNum(nrow(data.numlength.filter()), big.mark = ","), " out of ",
      prettyNum(nrow(data.numlength), big.mark = ","), " trajectories. [",
      round(nrow(data.numlength.filter())/nrow(data.numlength) * 100, digits = 2), "%]"
    )
  )
  output$numTraj2 <- renderText(
    paste0("Granules: ", 
           prettyNum(nrow(data.numlength.filter() %>% filter(SizeClass == "S")), big.mark = ","), " of ",
           prettyNum(nrow(filter(data.numlength, SizeClass == "S")), big.mark = ","), " [",
           round(nrow(data.numlength.filter() %>% filter(SizeClass == "S")) / nrow(filter(data.numlength, SizeClass == "S")) * 100, digits = 2),
           "%]; Scrums: ",
           prettyNum(nrow(data.numlength.filter() %>% filter(SizeClass == "L")), big.mark = ","), " of ",
           prettyNum(nrow(filter(data.numlength, SizeClass == "L")), big.mark = ","), " [", 
           round(nrow(data.numlength.filter() %>% filter(SizeClass == "L")) / nrow(filter(data.numlength, SizeClass == "L")) * 100, digits = 2),
           "%]."
    )
  )
  
  output$cdfPlot <- renderPlot(
    ggplot(data.numlength,aes(x=lengthtime, color = SizeClass)) + 
      stat_ecdf(geom = "step", size = 1.5) + 
      scale_color_brewer(palette = "Set1", labels = c("scrum", "granule"), name = NULL) + 
      #geom_vline(xintercept = 1+as.numeric(input$traj_length), size = 1) + 
      annotate("rect", xmin = as.numeric(input$traj_length[1]), xmax = as.numeric(input$traj_length[2]), 
               ymin = -Inf, ymax = Inf, fill = 'deepskyblue1', alpha = 0.4) +
      #scale_x_log10() + 
      labs(x = "trajectory length [s]", y = "CDF") + theme.dist + 
      theme(text = element_text(size = as.numeric(input$fontsize_adjust)), legend.position = c(0.85,0.12), legend.background = element_blank())
  )
  
  color.type <- reactive(
    if(input$type == "Granules") {
      scale_color_manual(values = mypalette[2], labels = c("granule"), name = NULL)
    } else if(input$type == "Scrums") {
      scale_color_manual(values = mypalette[1], labels = c("scrum"), name = NULL)
    } else {
      scale_color_manual(values = mypalette, labels = c("scrum", "granule"), name = NULL)
    }
  )
  
  output$msdPlot.gs <- renderPlot({
    p <- ggplot(data.filter(),aes(x=lagtime, y=MSDum)) + 
      geom_line(aes(color = SizeClass, group = uniqueID), alpha = as.numeric(input$alpha_adjust), size = as.numeric(input$line_adjust)) + 
      geom_point(aes(color = SizeClass), size = as.numeric(input$point_adjust), alpha = as.numeric(input$alpha_adjust)) + 
      scale_x_log10() + scale_y_log10(breaks = plotbreaks, labels = plotlabels) + color.type() +
      coord_cartesian(xlim = c(0.1,as.numeric(input$delta_adjust)), ylim = 10^as.numeric(input$msd_adjust)) +
      theme.dist + theme(text = element_text(size = as.numeric(input$fontsize_adjust)), legend.position = "bottom") + 
      labs(x = expression(Delta~"[s]"), y = "MSD ["~mu~"m"^{"2"}~"]")
    if(!is.null(input$display.table_row_last_clicked)) {
      seldata <- data %>% filter(uniqueID == data.display$uniqueID[input$display.table_row_last_clicked])
      p <- p + geom_line(data = seldata, aes(x=lagtime, y=MSDum), color = overlaycolor, alpha = as.numeric(input$alpha_adjust), size = as.numeric(input$line_adjust)) + 
        geom_point(data = seldata, aes(x=lagtime, y=MSDum), color = overlaycolor, size = as.numeric(input$point_adjust), alpha = as.numeric(input$alpha_adjust))
    }
    p}
  )
  
  output$msdPlot.alpha <- renderPlot({
    p <- ggplot(data.filter(),aes(x=lagtime, y=MSDum)) + 
      geom_line(aes(color = alpha, group = uniqueID), alpha = as.numeric(input$alpha_adjust), size = as.numeric(input$line_adjust)) + 
      geom_point(aes(color = alpha), size = as.numeric(input$point_adjust), alpha = as.numeric(input$alpha_adjust)) + 
      scale_x_log10() + scale_y_log10(breaks = plotbreaks, labels = plotlabels) + 
      coord_cartesian(xlim = c(0.1,as.numeric(input$delta_adjust)), ylim = 10^as.numeric(input$msd_adjust)) +
      scale_color_gradient2(low = 'red', high = 'blue', mid = 'white', midpoint = 1, na.value = 'purple', limits = c(0,2), 
                            breaks = c(0,0.5,1,1.5,2), name = expression(alpha), guide = guide_colorbar(barwidth = 15, barheight = 1.5)) + 
      theme_dark() + theme(text = element_text(size = as.numeric(input$fontsize_adjust)), legend.position = "bottom") + 
      labs(x = expression(Delta~"[s]"), y = "MSD ["~mu~"m"^{"2"}~"]")
    if(!is.null(input$display.table_row_last_clicked)) {
      seldata <- data %>% filter(uniqueID == data.display$uniqueID[input$display.table_row_last_clicked])
      p <- p + geom_line(data = seldata, aes(x=lagtime, y=MSDum), color = overlaycolor, alpha = as.numeric(input$alpha_adjust), size = as.numeric(input$line_adjust)) + 
        geom_point(data = seldata, aes(x=lagtime, y=MSDum), color = overlaycolor, size = as.numeric(input$point_adjust), alpha = as.numeric(input$alpha_adjust))
    }
    p}
  )
  
  output$msdPlot.Deff <- renderPlot({
    p <- ggplot(data.filter(),aes(x=lagtime, y=MSDum)) + 
      geom_line(aes(color = D, group = uniqueID), alpha = as.numeric(input$alpha_adjust), size = as.numeric(input$line_adjust)) + 
      geom_point(aes(color = D), size = as.numeric(input$point_adjust), alpha = as.numeric(input$alpha_adjust)) + 
      scale_x_log10() + scale_y_log10(breaks = plotbreaks, labels = plotlabels) + 
      coord_cartesian(xlim = c(0.1,as.numeric(input$delta_adjust)), ylim = 10^as.numeric(input$msd_adjust)) +
      scale_color_gradient(trans = "log10", na.value = 'purple', #limits = 10^deffrange, 
                           name = expression(Gamma~"["~mu~"m"^{"2"}~"/s]"), guide = guide_colorbar(barwidth = 15, barheight = 1.5)) + 
      theme.dist + theme(text = element_text(size = as.numeric(input$fontsize_adjust)), legend.position = "bottom") + 
      labs(x = expression(Delta~"[s]"), y = "MSD ["~mu~"m"^{"2"}~"]")
    if(!is.null(input$display.table_row_last_clicked)) {
      seldata <- data %>% filter(uniqueID == data.display$uniqueID[input$display.table_row_last_clicked])
      p <- p + geom_line(data = seldata, aes(x=lagtime, y=MSDum), color = overlaycolor, alpha = as.numeric(input$alpha_adjust), size = as.numeric(input$line_adjust)) + 
        geom_point(data = seldata, aes(x=lagtime, y=MSDum), color = overlaycolor, size = as.numeric(input$point_adjust), alpha = as.numeric(input$alpha_adjust))
    }
    p}
  )
  
  output$Dist.alpha <- renderPlot(
    ggplot(data.numlength.filter(),aes(x=alpha)) + 
      geom_histogram(alpha = 0.5, color = 'black', position = 'identity', binwidth = as.numeric(input$alpha_bin), aes(y = ..density.., fill = SizeClass)) + 
      geom_density(aes(color = SizeClass), size = 1+as.numeric(input$line_adjust)) +
      scale_fill_brewer(palette = "Set1", labels = c("scrum", "granule"), name = NULL) + 
      scale_color_brewer(palette = "Set1", guide = FALSE) + scale_x_continuous(limits = c(0,2)) +
      labs(x = expression(alpha), y = "pdf") + theme.dist + 
      theme(text = element_text(size = as.numeric(input$fontsize_adjust)), 
            legend.position = "right", legend.background = element_blank())
  )
  
  output$Dist.deff <- renderPlot(
    ggplot(data.numlength.filter(),aes(x=D)) + 
      geom_histogram(alpha = 0.5, color = 'black', position = 'identity', binwidth = as.numeric(input$deff_bin), aes(y = ..density.., fill = SizeClass)) + 
      geom_density(aes(color = SizeClass), size = 1+as.numeric(input$line_adjust)) +
      scale_fill_brewer(palette = "Set1", labels = c("scrum", "granule"), name = NULL) + 
      scale_color_brewer(palette = "Set1", guide = FALSE) + scale_x_log10() +
      labs(x = expression(Gamma), y = "pdf") + theme.dist + 
      theme(text = element_text(size = as.numeric(input$fontsize_adjust)), 
            legend.position = "right", legend.background = element_blank())
  )
  
  
  
  output$select.info <- renderTable({
    if(!is.null(input$msd.brush)) {
      sel <- brushedPoints(data.filter(), input$msd.brush) 
    } else {
      sel <- nearPoints(data.filter(), input$msd.click, maxpoints = 5)
    }
    sel %>%  mutate(type = dplyr::if_else(SizeClass == "S","granule","scrum")) %>% 
      mutate(ID = as.integer(ID)) %>% select(ID, type, alpha, gamma = D) %>% unique()
  }, include.rownames=FALSE, digits = 3, display = c("s","d","s","f","E"))
  
  output$select.traj <- renderPlot({
    sel <- data.frame()
    if(!is.null(input$msd.brush)) {
      sel <- brushedPoints(data.filter(), input$msd.brush) %>% select(uniqueID,ID,SizeClass) %>% unique()
      if(nrow(sel) > 4) {
        sel <- sel[1:4,]
      }
    } else {
      sel <- nearPoints(data.filter(), input$msd.click, maxpoints = 4) %>% select(uniqueID,ID,SizeClass) %>% unique()
    }
    if(nrow(sel) > 0) {
      selpos <- positions %>% filter(uniqueID == sel$uniqueID[1]) %>% 
        mutate(realtime = realtime-min(realtime), x = x-mean(x), y = y-mean(y)) 
      for(i in 2:nrow(sel)) {
        temp <- positions %>% filter(uniqueID == sel$uniqueID[i]) %>% 
          mutate(realtime = realtime-min(realtime), x = x-mean(x), y = y-mean(y))
        selpos <- bind_rows(selpos,temp)
      }
      limit <- max( abs(c(selpos$x, selpos$y)) )
      selpos$uniqueID <- factor(selpos$uniqueID, ordered = TRUE, levels = sel$uniqueID)
      scale <- data.frame(x = c(limit-scalelength,limit), y = c(-limit,-limit), uniqueID = selpos$uniqueID[1])
      sel2 <- sel %>%  mutate(type = dplyr::if_else(SizeClass == "S","granule","scrum"))
      label <- data.frame(x = -limit, y = limit, uniqueID = unique(sel2$uniqueID), lab = paste(sel2$ID,sel2$type))
      ps <- ggplot(selpos, aes(x = x, y = y, color = realtime)) + geom_path(size = 1.5) + 
        scale_color_viridis(name = "time (s)", option = "D", guide = guide_colorbar(barheight = 2, barwidth = 14)) +
        theme.traj + theme(legend.text = element_text(size = as.numeric(input$fontsize_adjust)-2), 
                           legend.title = element_text(size = as.numeric(input$fontsize_adjust)-4)) + 
        geom_path(data = scale, color = 'white', size = 3) + 
        geom_text(data = label, aes(label = lab), color = 'white', size = as.numeric(input$fontsize_adjust)/4, hjust = 0) +
        coord_fixed(xlim = c(-limit,limit), ylim = c(-limit,limit)) + facet_wrap(~uniqueID,ncol = 2)
    } else {
      ps <- ggplot() + geom_blank()
    }
    ps
  })
  
  output$display.table <- DT::renderDataTable(
    datatable(data.display, rownames = FALSE, filter = 'top', selection = 'single',
              options=list(dom = 'ltipr', 
                           pageLength = 5,
                           lengthMenu = c(5, 10, 25, 50),
                           columnDefs = list(list(visible=FALSE, targets=5),
                                             list(className = 'dt-center', targets = 1))
              )) %>% 
      formatRound('alpha', digits = 4) %>% formatRound('gamma', digits = 5)
  )
  
  output$table.traj <- renderPlot({
    if(!is.null(input$display.table_row_last_clicked)) {
      selpos <- positions %>% filter(uniqueID == data.display$uniqueID[input$display.table_row_last_clicked]) %>% 
        mutate(realtime = realtime-min(realtime), x = x-mean(x), y = y-mean(y)) 
      limit <- max( abs(c(selpos$x, selpos$y)) )
      scale <- data.frame(x = c(limit-scalelength,limit), y = c(-limit,-limit))
      pt <- ggplot(selpos, aes(x = x, y = y, color = realtime)) + geom_path(size = 1.5) + 
        scale_color_viridis(name = "time (s)", option = "D", guide = guide_colorbar(barheight = 2, barwidth = 14)) +
        theme.traj + theme(legend.text = element_text(size = as.numeric(input$fontsize_adjust)-2), 
                           legend.title = element_text(size = as.numeric(input$fontsize_adjust)-4)) +
        geom_path(data = scale, color = 'white', size = 3) +
        coord_fixed(xlim = c(-limit,limit), ylim = c(-limit,limit)) 
    } else{
      pt <- ggplot() + geom_blank()
    }
    pt
  })
  
}


# Run the application 
shinyApp(ui = ui, server = server)

