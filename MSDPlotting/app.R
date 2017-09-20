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

setwd('..')
filename <- "170427_3B11M_P13_Plate2a_Top"
load(paste0("Output/",filename,"_shiny.Rdata"))
data.numlength <- data %>% select(ID,SizeClass,lengthtime,D,alpha) %>% unique()
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
ctrl.width <- 3

ui <- fluidPage(
  withMathJax(),
  fluidRow(
    column(11,
           h3(paste0("Filename: ",filename))
    ),
    column(1,
           actionButton("close", "Close")
    )
  ),
  fluidRow(h6("Style Control")),
  fluidRow(
    column(ctrl.width,
           sliderInput("fontsize_adjust", label = "Font Size",
                       min = 8, max = 72, value = 24, step = 1, round = TRUE)
    ),
    column(ctrl.width,
           sliderInput("alpha_adjust", label = "Transparancy",
                       min = 0, max = 1, value = 0.8, step = 0.01, round = -2)
    ),
    column(ctrl.width,
           sliderInput("line_adjust", label = "Line Width",
                       min = 0.25, max = 5, value = 0.5, step = 0.25)
    ),
    column(ctrl.width,
           sliderInput("point_adjust", label = "Point Size",
                       min = 0.25, max = 5, value = 0.5, step = 0.25)
    )
  ),
  fluidRow(h6("Data Display Control")),
  fluidRow(
    column(ctrl.width,
           sliderInput("traj_length", label = "Trajectory Length [s]",
                       min = trajrange[1], max = trajrange[2], value = 25, step = 0.1, round = -1)
    ),
    column(ctrl.width,
           selectInput("type", label = "Type:",
                       choices = c("All","Granules","Scrums"))
    ),
    column(ctrl.width,
           sliderInput("delta_adjust", label = "\\( \\Delta \\) Range [s]",
                       min = deltarange[1], max = deltarange[2], value = 10, step = 0.1, round = -1)
    ),
    column(ctrl.width,
           sliderInput("msd_adjust", label = "MSD Range [\\(\\mu m^2\\)] (Powers of 10)",
                       min = msdrange[1], max = msdrange[2], value = c(-4,1), step = 0.1,
                       round = -1)
    )
  ),
  fluidRow(
    column(7,
           h4(textOutput("numTraj1")),
           h4(textOutput("numTraj2"))
    )
  ),
  fluidRow(
    column(7,
           plotOutput("msdPlot.gs")
    ),
    column(5,
           plotOutput("cdfPlot")
    )
  ),
  fluidRow(
    column(7,
           plotOutput("msdPlot.alpha")
    ),
    column(5,
           plotOutput("Dist.alpha")
    )
  ),
  fluidRow(
    column(width = 3, offset = 9,
           sliderInput("alpha_bin", label = "\\( \\alpha \\) Bin Size",
                       min = 0.01, max = 0.5, value = 0.1, step = 0.01, round = -2)
    )
  ),
  fluidRow(
    column(7,
           plotOutput("msdPlot.Deff")
    ),
    column(5,
           plotOutput("Dist.deff")
    )
  ),
  fluidRow(
    column(width = 3, offset = 9,
           sliderInput("deff_bin", label = "\\( D_{eff} \\) Bin Size",
                       min = 0.01, max = 1, value = 0.1, step = 0.01, round = -2)
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  observe({
    if (input$close > 0) {
      stopApp() 
    }
  })
  
  
  data.filter <- reactive(
    if(input$type == "Granules") {
      filter(data,lengthtime >= as.numeric(input$traj_length)) %>% filter(SizeClass == "S")
    } else if(input$type == "Scrums") {
      filter(data,lengthtime >= as.numeric(input$traj_length)) %>% filter(SizeClass == "L")
    } else {
      filter(data,lengthtime >= as.numeric(input$traj_length))
    }
  )
  
  data.numlength.filter <- reactive(
    filter(data.numlength, lengthtime >= as.numeric(input$traj_length))
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
  
  color.type <- reactive(
    if(input$type == "Granules") {
      scale_color_manual(values = mypalette[2], labels = c("granule"), name = NULL)
    } else if(input$type == "Scrums") {
      scale_color_manual(values = mypalette[1], labels = c("scrum"), name = NULL)
    } else {
      scale_color_manual(values = mypalette, labels = c("scrum", "granule"), name = NULL)
    }
  )
  
  output$msdPlot.gs <- renderPlot(
    ggplot(data.filter(),aes(x=lagtime, y=MSDum)) + 
      geom_line(aes(color = SizeClass, group = uniqueID), alpha = as.numeric(input$alpha_adjust), size = as.numeric(input$line_adjust)) + 
      geom_point(aes(color = SizeClass), size = as.numeric(input$point_adjust), alpha = as.numeric(input$alpha_adjust)) + 
      scale_x_log10() + scale_y_log10(breaks = plotbreaks, labels = plotlabels) + color.type() +
      coord_cartesian(xlim = c(0.1,as.numeric(input$delta_adjust)), ylim = 10^as.numeric(input$msd_adjust)) +
      theme.dist + theme(text = element_text(size = as.numeric(input$fontsize_adjust)), legend.position = "bottom") + 
      labs(x = expression(Delta~"[s]"), y = "MSD ["~mu~"m"^{"2"}~"]")
  )
  
  output$cdfPlot <- renderPlot(
    ggplot(data.numlength,aes(x=lengthtime, color = SizeClass)) + 
      stat_ecdf(geom = "step", size = 1+as.numeric(input$line_adjust)) + 
      scale_color_brewer(palette = "Set1", labels = c("scrum", "granule"), name = NULL) + 
      geom_vline(xintercept = 1+as.numeric(input$traj_length), size = as.numeric(input$line_adjust)) + 
      scale_x_log10() + labs(x = "trajectory length [s]", y = "CDF") + theme.dist + 
      theme(text = element_text(size = as.numeric(input$fontsize_adjust)), legend.position = "bottom", legend.background = element_blank())
  )
  
  output$msdPlot.alpha <- renderPlot(
    ggplot(data.filter(),aes(x=lagtime, y=MSDum)) + 
      geom_line(aes(color = alpha, group = uniqueID), alpha = as.numeric(input$alpha_adjust), size = as.numeric(input$line_adjust)) + 
      geom_point(aes(color = alpha), size = as.numeric(input$point_adjust), alpha = as.numeric(input$alpha_adjust)) + 
      scale_x_log10() + scale_y_log10(breaks = plotbreaks, labels = plotlabels) + 
      coord_cartesian(xlim = c(0.1,as.numeric(input$delta_adjust)), ylim = 10^as.numeric(input$msd_adjust)) +
      scale_color_gradient2(low = 'red', high = 'blue', mid = 'white', midpoint = 1, na.value = 'grey50', limits = c(0,2), 
                            breaks = c(0,0.5,1,1.5,2), name = expression(alpha), guide = guide_colorbar(barwidth = 15, barheight = 1.5)) + 
      theme_dark() + theme(text = element_text(size = as.numeric(input$fontsize_adjust)), legend.position = "bottom") + 
      labs(x = expression(Delta~"[s]"), y = "MSD ["~mu~"m"^{"2"}~"]")
  )
  
  output$msdPlot.Deff <- renderPlot(
    ggplot(data.filter(),aes(x=lagtime, y=MSDum)) + 
      geom_line(aes(color = D, group = uniqueID), alpha = as.numeric(input$alpha_adjust), size = as.numeric(input$line_adjust)) + 
      geom_point(aes(color = D), size = as.numeric(input$point_adjust), alpha = as.numeric(input$alpha_adjust)) + 
      scale_x_log10() + scale_y_log10(breaks = plotbreaks, labels = plotlabels) + 
      coord_cartesian(xlim = c(0.1,as.numeric(input$delta_adjust)), ylim = 10^as.numeric(input$msd_adjust)) +
      scale_color_gradient(trans = "log10", na.value = 'grey50', #limits = 10^deffrange, 
                           name = expression("D"["eff"]~"["~mu~"m"^{"2"}~"/s]"), guide = guide_colorbar(barwidth = 15, barheight = 1.5)) + 
      theme.dist + theme(text = element_text(size = as.numeric(input$fontsize_adjust)), legend.position = "bottom") + 
      labs(x = expression(Delta~"[s]"), y = "MSD ["~mu~"m"^{"2"}~"]")
  )
  
  output$Dist.alpha <- renderPlot(
    ggplot(data.numlength.filter(),aes(x=alpha)) + 
      geom_histogram(alpha = 0.5, color = 'black', position = 'identity', binwidth = as.numeric(input$alpha_bin), aes(y = ..density.., fill = SizeClass)) + 
      geom_density(aes(color = SizeClass), size = 1+as.numeric(input$line_adjust)) +
      scale_fill_brewer(palette = "Set1", labels = c("scrum", "granule"), name = NULL) + 
      scale_color_brewer(palette = "Set1", guide = FALSE) + scale_x_continuous(limits = c(0,2)) +
      labs(x = expression(alpha), y = "pdf") + theme.dist + 
      theme(text = element_text(size = as.numeric(input$fontsize_adjust)), 
            legend.position = "bottom", legend.background = element_blank())
  )
  
  output$Dist.deff <- renderPlot(
    ggplot(data.numlength.filter(),aes(x=D)) + 
      geom_histogram(alpha = 0.5, color = 'black', position = 'identity', binwidth = as.numeric(input$deff_bin), aes(y = ..density.., fill = SizeClass)) + 
      geom_density(aes(color = SizeClass), size = 1+as.numeric(input$line_adjust)) +
      scale_fill_brewer(palette = "Set1", labels = c("scrum", "granule"), name = NULL) + 
      scale_color_brewer(palette = "Set1", guide = FALSE) + scale_x_log10() +
      labs(x = expression("D"["eff"]), y = "pdf") + theme.dist + 
      theme(text = element_text(size = as.numeric(input$fontsize_adjust)), 
            legend.position = "bottom", legend.background = element_blank())
  )
}


# Run the application 
shinyApp(ui = ui, server = server)

