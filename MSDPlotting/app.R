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

theme.dist <- theme_bw() + theme(panel.grid.major.x = element_line(color = 'grey75', size = 0.3),
                               panel.grid.minor.x = element_line(color = 'grey85', size = 0.25), 
                               panel.grid.major.y = element_line(color = 'grey75', size = 0.3), 
                               panel.grid.minor.y = element_line(color = 'grey85', size = 0.25), 
                               panel.border = element_rect(fill = NA, color = 'black', size = 0.75),
                               legend.title = element_blank(), legend.key = element_blank())

ui <- fixedPage(
  withMathJax(),
  fixedRow(width = 12,
           column(12,
                  h3(paste0("Filename: ",filename))
           )
  ),
  fixedRow(width = 12,
           column(3,
                  sliderInput("fontsize_adjust", label = "Font Size",
                              min = 8, max = 72, value = 24, step = 1, round = TRUE)
           ),
           column(3,
                  sliderInput("alpha_adjust", label = "Transparancy",
                              min = 0, max = 1, value = 0.8, step = 0.01, round = -2)
           ),
           column(3,
                  sliderInput("line_adjust", label = "Line Width",
                              min = 0.25, max = 5, value = 0.5, step = 0.25)
           ),
           column(3,
                  sliderInput("point_adjust", label = "Point Size",
                              min = 0.25, max = 5, value = 0.5, step = 0.25)
           )
  ),
  fixedRow(width = 12,
           column(6,
                  sliderInput("delta_adjust", label = "Delta Range [s]",
                              min = deltarange[1], max = deltarange[2], value = 10, step = 0.1, round = -1)
           ),
           column(6,
                  sliderInput("msd_adjust", label = "MSD Range [\\(\\mu m^2\\)] (Powers of 10)",
                              min = msdrange[1], max = msdrange[2], value = c(-4,1), step = 0.1,
                              round = -1)
           )
  ),
  fixedRow(width = 12,
           column(6,
                  sliderInput("traj_length", label = "Trajectory Length [s]",
                              min = trajrange[1], max = trajrange[2], value = 25, step = 0.1, round = -1)
           ),
           column(3,
                  sliderInput("alpha_bin", label = "\\( \\alpha \\) Bin Size",
                              min = 0.01, max = 0.5, value = 0.1, step = 0.01, round = -2)
           ),
           column(3,
                  sliderInput("deff_bin", label = "\\( D_{eff} \\) Bin Size",
                              min = 0.001, max = 1, value = 0.01, step = 0.001, round = -2)
           )
  ),
  fixedRow(width = 12,
           column(12,
                  h4(textOutput("numTraj1")),
                  h4(textOutput("numTraj2"))
           )
  ),
  fixedRow(
    column(7,
           plotOutput("msdPlot.gs")
    ),
    column(5,
           plotOutput("cdfPlot")
    )
  ),
  fixedRow(
    column(7,
           plotOutput("msdPlot.alpha")
    ),
    column(5,
           plotOutput("Dist.alpha")
    )
  ),
  fixedRow(
    column(7,
           plotOutput("msdPlot.Deff")
    ),
    column(5,
           plotOutput("Dist.deff")
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

  data.filter <- reactive(
    filter(data,lengthtime >= as.numeric(input$traj_length))
  )
  
  data.numlength.filter <- reactive(
    filter(data.numlength, lengthtime >= as.numeric(input$traj_length))
  )

  output$numTraj1 <- renderText(
    paste0(
      "Displaying ", nrow(data.numlength.filter()), " out of ",
      nrow(data.numlength), " trajectories. [",
      round(nrow(data.numlength.filter())/nrow(data.numlength) * 100, digits = 2), "%]"
      )
  )
  output$numTraj2 <- renderText(
    paste0("Granules: ", 
      nrow(data.numlength.filter() %>% filter(SizeClass == "S")), " of ",
      nrow(filter(data.numlength, SizeClass == "S")), " [",
      round(nrow(data.numlength.filter() %>% filter(SizeClass == "S")) / nrow(filter(data.numlength, SizeClass == "S")) * 100, digits = 2),
      "%]; Scrums: ",
      nrow(data.numlength.filter() %>% filter(SizeClass == "L")), " of ",
      nrow(filter(data.numlength, SizeClass == "L")), " [", 
      round(nrow(data.numlength.filter() %>% filter(SizeClass == "L")) / nrow(filter(data.numlength, SizeClass == "L")) * 100, digits = 2),
      "%]."
    )
  )
  
  output$msdPlot.gs <- renderPlot(
    ggplot(data.filter(),aes(x=lagtime, y=MSDum)) + 
      geom_line(aes(color = SizeClass, group = uniqueID), alpha = as.numeric(input$alpha_adjust), size = as.numeric(input$line_adjust)) + 
      geom_point(aes(color = SizeClass), size = as.numeric(input$point_adjust), alpha = as.numeric(input$alpha_adjust)) + 
      scale_x_log10(limits = c(0.1,as.numeric(input$delta_adjust))) + 
      scale_y_log10(limits = 10^as.numeric(input$msd_adjust), breaks = plotbreaks, labels = plotlabels) + 
      scale_color_brewer(palette = "Set1", labels = c("scrum", "granule"), name = NULL) +
      theme.dist + theme(text = element_text(size = as.numeric(input$fontsize_adjust)), legend.position = "bottom") + 
      labs(x = expression(Delta~"[s]"), y = "MSD ["~mu~"m"^{"2"}~"]")
  )
  
  output$cdfPlot <- renderPlot(
    ggplot(data.numlength,aes(x=lengthtime, color = SizeClass)) + 
      stat_ecdf(geom = "step", size = 1+as.numeric(input$line_adjust)) + scale_color_brewer(palette = "Set1", labels = c("scrum", "granule")) + 
      geom_vline(xintercept = as.numeric(input$traj_length)) + 
      scale_x_log10() + labs(x = "trajectory length [s]", y = "CDF") + theme.dist + 
      theme(text = element_text(size = as.numeric(input$fontsize_adjust)), legend.position = "bottom", legend.background = element_blank())
  )
  
   output$msdPlot.alpha <- renderPlot(
     ggplot(data.filter(),aes(x=lagtime, y=MSDum)) + 
       geom_line(aes(color = alpha, group = uniqueID), alpha = as.numeric(input$alpha_adjust), size = as.numeric(input$line_adjust)) + 
       geom_point(aes(color = alpha), size = as.numeric(input$point_adjust), alpha = as.numeric(input$alpha_adjust)) + 
       scale_x_log10(limits = c(0.1,as.numeric(input$delta_adjust))) + 
       scale_y_log10(limits = 10^as.numeric(input$msd_adjust), breaks = plotbreaks, labels = plotlabels) + 
       scale_color_gradient2(low = 'red', high = 'blue', mid = 'white', midpoint = 1, na.value = 'grey50', limits = c(0,2), 
                             breaks = c(0,0.5,1,1.5,2), name = expression(alpha), guide = guide_colorbar(barwidth = 15, barheight = 1.5)) + 
       theme_dark() + theme(text = element_text(size = as.numeric(input$fontsize_adjust)), legend.position = "bottom") + 
       labs(x = expression(Delta~"[s]"), y = "MSD ["~mu~"m"^{"2"}~"]")
   )
   
   output$msdPlot.Deff <- renderPlot(
     ggplot(data.filter(),aes(x=lagtime, y=MSDum)) + 
       geom_line(aes(color = D, group = uniqueID), alpha = as.numeric(input$alpha_adjust), size = as.numeric(input$line_adjust)) + 
       geom_point(aes(color = D), size = as.numeric(input$point_adjust), alpha = as.numeric(input$alpha_adjust)) + 
       scale_x_log10(limits = c(0.1,as.numeric(input$delta_adjust))) + 
       scale_y_log10(limits = 10^as.numeric(input$msd_adjust), breaks = plotbreaks, labels = plotlabels) + 
       scale_color_gradient(na.value = 'grey50',  name = expression("D"["eff"]), guide = guide_colorbar(barwidth = 15, barheight = 1.5)) + 
       theme_dark() + theme(text = element_text(size = as.numeric(input$fontsize_adjust)), legend.position = "bottom") + 
       labs(x = expression(Delta~"[s]"), y = "MSD ["~mu~"m"^{"2"}~"]")
   )
   
   output$Dist.alpha <- renderPlot(
     ggplot(data.numlength.filter(),aes(x=alpha)) + 
       geom_histogram(alpha = 0.5, color = 'black', position = 'identity', binwidth = as.numeric(input$alpha_bin), aes(y = ..density.., fill = SizeClass)) + 
       geom_density(aes(color = SizeClass), size = 1+as.numeric(input$line_adjust)) +
       scale_fill_brewer(palette = "Set1", labels = c("scrum", "granule")) + 
       scale_color_brewer(palette = "Set1", guide = FALSE) + scale_x_continuous(limits = c(0,2)) +
       labs(x = expression(alpha), y = "pdf") + theme.dist + 
       theme(text = element_text(size = as.numeric(input$fontsize_adjust)), 
             legend.position = "bottom", legend.background = element_blank())
   )
   
   output$Dist.deff <- renderPlot(
     ggplot(data.numlength.filter(),aes(x=D)) + 
       geom_histogram(alpha = 0.5, color = 'black', position = 'identity', binwidth = as.numeric(input$deff_bin), aes(y = ..density.., fill = SizeClass)) + 
       geom_density(aes(color = SizeClass), size = 1+as.numeric(input$line_adjust)) +
       scale_fill_brewer(palette = "Set1", labels = c("scrum", "granule")) + 
       scale_color_brewer(palette = "Set1", guide = FALSE) + scale_x_continuous() +
       labs(x = expression("D"["eff"]), y = "pdf") + theme.dist + 
       theme(text = element_text(size = as.numeric(input$fontsize_adjust)), 
             legend.position = "bottom", legend.background = element_blank())
   )
}


# Run the application 
shinyApp(ui = ui, server = server)

