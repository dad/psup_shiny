# install.packages(c("shiny","data.table","ggplot2","devtools"))
# devtools::install_github("aoles/shinyURL")
# query strings:
# psup_shiny/?ids=PGK1,OLA1,PMA1&idType=gene&interval=F&plotType=time
# psup_shiny/?ids=YER165W&idType=orf&interval=T&plotType=time

# To publish:
# runGitHub( "psup_shiny", "dad")
# and Publish

library(shiny)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(shinyURL)

ps_dt <- read_tsv("scer_aggregation_psup_long.txt", comment='#') #,stringsAsFactors=FALSE,header=TRUE)

theme_set(theme_minimal(base_size=14) %+replace% theme(legend.position="none"))

splitstring = function(s) {
  # Take string formatted as "ribosomal=(RPS1,RPS2); glycolytic=(GLK1,TDH3)" and split into a
  # list(ribosomal=('RPS1','RPS2'), glycolytic=c('GLK1','TDH3'))
  #if (T) {  s = "ribosomal=(rps1,RPS2); glycolytic=(GLK1,TDH3); agg=(ded1,pab1); pma1" }
  classes = str_split(s,';')[[1]]
  
  x = sapply(classes, function(ss) {
    res = str_split(ss,'=')[[1]]
    if (length(res)==2){
      label=gsub("[[:space:]]",'',res[1])
      genenames = res[2]
    } else {
      # No label; use gene as label
      trimname = gsub("[[:space:]]",'',res[1])
      label=toupper(trimname)
      genenames = trimname
    }
    genes=toupper(str_split(gsub("\\((.*)\\)", "\\1",genenames),',')[[1]])
    g = list(g=genes)
    names(g)[1] = label
    #print(g)
    g
  }, USE.NAMES = F)
  x
}

scale_y_pSup <- function(name=NULL) {
    list(scale_y_continuous(name,expand=c(0.01,0),limits=c(0,1),
                       #breaks=seq(0,1,.25),labels=c("0","","0.5","","1"))
                       breaks=c(0,1)),
         theme(axis.title.y=element_text(angle=0, vjust=0.5)))
}

scale_time <- function(name=expression("Minutes at "*46*degree~C*""),
                       text=TRUE) {
    if (text) {
        mylim=c(0,9.9)
    } else {
        mylim=c(0,8)
    }
    scale_x_continuous(name,expand=c(0.01,0.1),limits=mylim,
                       breaks=c(0,2,4,8))
}

plotmygenes = function(mygenes,data=ps_dt,
                       tempexps=c("30C.rep2","37C.8min","42C.8min","46C.8min"),
                       temps=c(30,37,42,46), tempbreaks=temps,
                       timeexps=c("30C.rep1","46C.2min","46C.4min","46C.8min"),
                       times=c(0,2,4,8), 
                       errorbars=FALSE,linewidth=0.8,
                       idType=c("gene","orf")) {
  names(temps) = tempexps
  names(times) = timeexps
  
  if(idType=="gene") {
    ps_dt_temp = ps_dt |> filter(gene %in% mygenes & experiment %in% tempexps)
    ps_dt_time = ps_dt |> filter(gene %in% mygenes & experiment %in% timeexps)
  } else if(idType=="orf") {
    ps_dt_temp = ps_dt |> filter(orf %in% mygenes & experiment %in% tempexps)
    ps_dt_time = ps_dt |> filter(orf %in% mygenes & experiment %in% timeexps)
  } else {
    stop("idType must be gene or orf")
  }
  ps_dt_temp$temp = temps[ps_dt_temp$experiment ]
  plot_temp = ggplot(data=ps_dt_temp,
                     aes(x=.data[["temp"]],y=.data[["psup"]],ymin=.data[["psup.lo"]],ymax=.data[["psup.hi"]],
                         colour=idType,label=idType)) +
    geom_line(linewidth=linewidth) + 
    geom_text_repel(size=4,data=ps_dt_temp |> filter(temp==max(temp)),
                    aes(x=max(temp)+0.5,y=psup), xlim=c(46,52)) +
    coord_cartesian(xlim=c(30,52)) +
    #        scale_x_continuous(expression("Temperature "*(degree*C)*" of 8 min. treatment"),
    scale_x_continuous("Temperature (\u00b0C) of 8 min. treatment",
                       breaks=tempbreaks,labels=tempbreaks, expand=c(0,0)) +
    scale_y_pSup("Proportion\nin\nsupernatant")
  
  ps_dt_time$time = times[ps_dt_time$experiment ]
  
  plot_time = ggplot(data=ps_dt_time,
                     aes_string(x="time",y="psup",ymin="psup.lo",ymax="psup.hi",
                                colour=idType,label=idType)) +
    geom_line(size=linewidth) +
    geom_text_repel(size=4,data=subset(ps_dt_time,time==max(time)),
                    aes(x=max(time)+0.1,y=psup), xlim=c(8,12)) +
    scale_y_pSup("Proportion\nin\nsupernatant") + scale_time()
  
  if (errorbars) {
    plot_temp = plot_temp + geom_pointrange()
    plot_time = plot_time + geom_pointrange()
  }
  
  return(list(plot_time=plot_time,plot_temp=plot_temp))
  
}

plotgenes_categories = function(gene_cat_list,data=ps_dt,
                                tempexps=c("30C.rep2","37C.8min","42C.8min","46C.8min"),
                                temps=c(30,37,42,46), tempbreaks=temps,
                                timeexps=c("30C.rep1","46C.2min","46C.4min","46C.8min"),
                                times=c(0,2,4,8), 
                                errorbars=FALSE,linewidth=0.8,
                                idType=c("gene","orf")) {
  names(temps) = tempexps
  names(times) = timeexps
  
  if (!is_list(gene_cat_list)) {
    gene_cat_list = list("<unnamed>"=gene_cat_list)
  }
  
  ps_dt_temp_all = bind_rows(lapply(names(gene_cat_list), function(cat) {
    mygenes = gene_cat_list[[cat]]
    ps_dt |> filter(.data[[idType]] %in% mygenes & experiment %in% tempexps) |> mutate(category=cat)
  }))
  
  ps_dt_time_all = bind_rows(lapply(names(gene_cat_list), function(cat) {
    mygenes = gene_cat_list[[cat]]
    ps_dt |> filter(.data[[idType]] %in% mygenes & experiment %in% timeexps) |> mutate(category=cat)
  }))
  
  # Summarize by category
  ps_dt_temp = ps_dt_temp_all |> group_by(experiment, category) |> summarise(sd=sd(psup,na.rm=T), se=sd/sqrt(n()), psup=mean(psup), psup.lo=psup-sd, psup.hi=psup+sd, category=first(category))
  ps_dt_time = ps_dt_time_all |> group_by(experiment, category) |> summarise(sd=sd(psup,na.rm=T), se=sd/sqrt(n()), psup=mean(psup), psup.lo=psup-sd, psup.hi=psup+sd, category=first(category))
  
  ps_dt_temp$temp = temps[ps_dt_temp$experiment ]
  plot_temp = ggplot(data=ps_dt_temp,
                     aes(x=.data[["temp"]],y=.data[["psup"]],ymin=.data[["psup.lo"]],ymax=.data[["psup.hi"]],
                         colour=category,label=category)) +
    geom_line(linewidth=linewidth) + 
    geom_text_repel(size=4,data=ps_dt_temp |> filter(temp==max(temps)),
                    aes(x=max(temps)+0.5,y=psup), xlim=c(46,52)) +
    coord_cartesian(xlim=c(30,52)) +
    #        scale_x_continuous(expression("Temperature "*(degree*C)*" of 8 min. treatment"),
    scale_x_continuous("Temperature (\u00b0C) of 8 min. treatment",
                       breaks=tempbreaks,labels=tempbreaks, expand=c(0,0)) +
    scale_y_pSup("Proportion\nin\nsupernatant")
  
  ps_dt_time$time = times[ps_dt_time$experiment ]
  
  plot_time = ggplot(data=ps_dt_time,
                     aes(x=.data[["time"]],y=.data[["psup"]],ymin=.data[["psup.lo"]],ymax=.data[["psup.hi"]],
                         colour=category,label=category)) +
    geom_line(size=linewidth) +
    geom_text_repel(size=4,data=subset(ps_dt_time,time==max(times)),
                    aes(x=max(times)+0.1,y=psup), xlim=c(8,12)) +
    scale_y_pSup("Proportion\nin\nsupernatant") + scale_time()
  
  if (errorbars) {
    plot_temp = plot_temp + geom_pointrange()
    plot_time = plot_time + geom_pointrange()
  }
  
  return(list(plot_time=plot_time,plot_temp=plot_temp))
  
}

plot_from_input = function(s, errorbars, idType) {
  if (grepl('=',s,fixed=TRUE)) { 
    # This is a named list
    ids = splitstring(s)
    print(paste("named categories:",ids))
    fun = plotgenes_categories
  } else {
    ids = toupper(strsplit(gsub(" ", "", s, fixed = TRUE),",")[[1]])
    print(paste("flat list:",ids))
    fun = plotmygenes
  }
  fun(ids, errorbars=errorbars, idType=idType)
}

# Define server logic required to draw a plot
shinyServer(function(input, output,session) {
    
    # Expression that generates a timecourse. The expression is
    # wrapped in a call to renderPlot to indicate that:
    #
    #  1) It is "reactive" and therefore should re-execute automatically
    #     when inputs change
    #  2) Its output type is a plot
    
    output$plot <- renderPlot({
        # time plot of 46C psup
        
        plots = plot_from_input(input$ids, errorbars=input$interval,idType=input$idType)
        if (input$plotType == "time") {
            return(plots$plot_time)
        } else if (input$plotType == "temperature") {
            return(plots$plot_temp)
        }
    })
    shinyURL.server(session)
})