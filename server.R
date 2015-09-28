library(shiny)
library(data.table)
library(ggplot2)

##### 
# To do:
# -- paper reference on user interface
# -- switch to view by orf DONE
# -- specify orf by url

ps_dt <- data.table(read.table("scer_aggregation_psup_long.txt",stringsAsFactors=FALSE,header=TRUE))

theme_set(theme_minimal(base_size=14) %+replace% theme(legend.position="none"))

scale_y_pSup <- function(name=NULL) {
    list(scale_y_continuous(name,expand=c(0.01,0),limits=c(0,1),
                       #breaks=seq(0,1,.25),labels=c("0","","0.5","","1"))
                       breaks=c(0,1)),
         theme(axis.title.y=element_text(angle=0)))
}

scale_time <- function(name=expression("time at "*46*degree~C*" (mins)"),
                       text=TRUE) {
    if (text) {
        mylim=c(0,9.9)
    } else {
        mylim=c(0,8)
    }
    scale_x_continuous(name,expand=c(0.01,0.1),limits=mylim,
                       breaks=c(0,2,4,8))
}

plotmygenes <- function(mygenes,data=ps_dt,
                        tempexps=c("30C.rep2","37C.8min","42C.8min","46C.8min"),
                        temps=c(30,37,42,46), tempbreaks=temps,
                        timeexps=c("30C.rep1","46C.2min","46C.4min","46C.8min"),
                        times=c(0,2,4,8), 
                        errorbars=FALSE,linesize=0.8,
                        idType=c("gene","orf")) {
    names(temps) <- tempexps
    names(times) <- timeexps
    
    if(idType=="gene") {
        ps_dt_temp <- subset(ps_dt, gene %in% mygenes & experiment %in% tempexps)
        ps_dt_time <- subset(ps_dt, gene %in% mygenes & experiment %in% timeexps)
    } else if(idType=="orf") {
        ps_dt_temp <- subset(ps_dt, orf %in% mygenes & experiment %in% tempexps)
        ps_dt_time <- subset(ps_dt, orf %in% mygenes & experiment %in% timeexps)
    } else {
        stop("idType must be gene or orf")
    }
    ps_dt_temp$temp <- temps[ps_dt_temp$experiment ]
    plot_temp <- ggplot(data=ps_dt_temp,
                        aes_string(x="temp",y="psup",ymin="psup.lo",ymax="psup.hi",
                                   colour=idType,label=idType)) +
        geom_line(size=linesize) + 
        geom_text(size=4,data=subset(ps_dt_temp,temp==max(temp)),
                  aes(x=max(temp)+2,y=psup)) +
        scale_x_continuous(expression("temperature "*(degree*C)*" of 8 min shock"),breaks=tempbreaks,labels=tempbreaks) +
        scale_y_pSup("pSup")
    
    ps_dt_time$time <- times[ps_dt_time$experiment ]
    
    plot_time <- ggplot(data=ps_dt_time,
                        aes_string(x="time",y="psup",ymin="psup.lo",ymax="psup.hi",
                                   colour=idType,label=idType)) +
        geom_line(size=linesize) +
        geom_text(size=4,data=subset(ps_dt_time,time==max(time)),
                  aes(x=max(time)+1.1,y=psup)) +
        scale_y_pSup("pSup") + scale_time()
    
    if (errorbars) {
        plot_temp <- plot_temp + geom_pointrange()
        plot_time <- plot_time + geom_pointrange()
    }
    
    return(list(plot_time=plot_time,plot_temp=plot_temp))

}

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
    
    # Expression that generates a histogram. The expression is
    # wrapped in a call to renderPlot to indicate that:
    #
    #  1) It is "reactive" and therefore should re-execute automatically
    #     when inputs change
    #  2) Its output type is a plot
    
    output$plot <- renderPlot({
        # time plot of 46C psup
        ids <- toupper(strsplit(gsub(" ", "", input$ids, fixed = TRUE),",")[[1]])
        plots <- plotmygenes(ids,errorbars=input$interval,idType=input$idType)
        if (input$plotType == "time") {
            return(plots$plot_time)
        } else if (input$plotType == "temperature") {
            return(plots$plot_temp)
        }
    })
#     output$tempPlot <- renderPlot({
#         # temperature plot of 46C psup
#         genes <- strsplit(input$genes,",")[[1]]
#         plotmygenes(genes,errorbars=input$interval)$plot_temp
#     })
})