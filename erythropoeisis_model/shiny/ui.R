# downloadButton

library(shiny)
#print('entered UI.r')
shinyUI(pageWithSidebar(
  headerPanel=headerPanel('Erythropoeisis model'),
  
  sidebarPanel=sidebarPanel(
    
    submitButton('Apply changes ...'), br(),
    
    h6("Set parameters"),
    
    sliderInput(inputId='ss',label='Stem cell to stem cell cmpt',min=0,max=0.1,value=0.02),
    
    sliderInput(inputId='sc',label='Stem cell to C cmpt',min=0,max=0.1,value=0.01),
    
    sliderInput(inputId='sd',label='Stem cell death rate',min=0,max=0.1,value=0.01), br(),
    
    sliderInput(inputId='cp',label='C cmpt to  P cmpt',min=0,max=1,value=0.1),
    
    sliderInput(inputId='cc',label='C cmpt recycling rate',min=0,max=1,value=0.05),
    
    sliderInput(inputId='h',label='C cpmt death controller',min=0,max=0.2,value=0.1), br(),
    
    sliderInput(inputId='pp',label='P cmpt recycling rate',min=0,max=1,value=0.1),
    
    sliderInput(inputId='pd',label='P cmpt death rate',min=0,max=0.2,value=0.05),
    
    sliderInput(inputId='tau.p',label='P cmpt to R cmpt delay time',min=0,max=10,value=5), br(),
    
    sliderInput(inputId='tau.r',label='R cmpt death delay time',min=0,max=200,value=100),
    
    numericInput(inputId='rh',label='R to Hb converter',value=1e-8,min=0), br(),
    
    h6('Set initial values'),
    
    numericInput(inputId='S',label='Stem cell numbers',value=1e10,min=0), 
    
    numericInput(inputId='C',label='C cmpt cell numbers',value=0,min=0),
    
    numericInput(inputId='cd',label='Initial C cmpt death rate',value=0,min=0),
    
    numericInput(inputId='P',label='P cmpt cell numbers',value=0,min=0), 
    
    numericInput(inputId='R',label='R cmpt cell numbers',value=0,min=0), 
    
    numericInput(inputId='newP',label='New P cmpt cell numbers at start',value=0,min=0),
    
    numericInput(inputId='newR',label='New R cmpt cell numbers at start',value=0,min=0),
    
    numericInput(inputId='tf',label='Final time',value=500,min=0),
    
    numericInput(inputId='dt',label='Time units',value=1,min=0), br(),
    
    h6('Enter events details'),
    
    textInput(inputId='vars',label='Variables to be changed (comma sep)'),
    
    textInput(inputId='time',label='Times at which events occur'),
    
    textInput(inputId='val',label='Values of the variables'),
    
    textInput(inputId='method',label='Methods (add, mult, rep)'), br(),
    
    textInput(inputId='dPrint',label='Times at which state variable values need to be saved'), br(),
    
    downloadButton(outputId='downloadParam',label='Download'),
    helpText('Download initial parameters and state variables at end and at specified times.')
        
  ),
  mainPanel=mainPanel(
    tabsetPanel(
      tabPanel("Plot",plotOutput(outputId='mPlot',height=600)),
      tabPanel("Text",verbatimTextOutput(outputId='initParam_stateVar'))
    )
    
  )
))

#print('exited ui.R')