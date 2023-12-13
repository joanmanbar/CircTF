# Load required libraries
library(shiny)
library(ggplot2)
library(data.table)
library(tidyr)
library(dplyr)

# Load data
all_TPM <- fread("/Users/joanbarreto/Downloads/All_Genotypes_TPM.csv")
CircadianGenes <- read.csv("/Users/joanbarreto/Downloads/circadianANDnoncirc.csv")
# sort the data by JTK_BH.Q (corrected p-value)
CircadianGenes <- CircadianGenes[order(CircadianGenes$JTK_BH.Q),]
CircadianGenes <- CircadianGenes[,c("Genotype", "CycID", "JTK_BH.Q", "Label")]
setnames(CircadianGenes, "CycID", "Geneid") # rename CycID to Geneid

# # Define parameters
# From <- 10
# To <- 15
# Rank <- 1
# SortByGenotype <- TRUE
# Genotype <- "A03"
# Means <- TRUE
# SD <- TRUE


# Single gene function
SingleGene <- function(Rank=1, InvertRank=FALSE, Genotype=NULL){
  Candidate <- CircadianGenes
  if(InvertRank == TRUE){
      Candidate <- Candidate[order(Candidate$JTK_BH.Q, decreasing = TRUE),]
    }
  # check if Genotype is provided
  if(!is.null(Genotype)){
    Candidate <- Candidate[which(Candidate$Genotype == Genotype),]
    Candidate <- Candidate$Geneid
    df_subset <- all_TPM[which(all_TPM$Geneid %in% Candidate),]
    df_subset <- df_subset[!duplicated(df_subset$Geneid), ] # Remove duplicated rows
    df_subset <- df_subset[Rank,]
  } else {
    Candidate <- CircadianGenes$Geneid[Rank]
    df_subset <- all_TPM[which(all_TPM$Geneid %in% Candidate),]
    df_subset <- df_subset[!duplicated(df_subset$Geneid), ] # Remove duplicated rows
  }
  # Melt the data
  df <- df_subset %>%
    gather(key, value, starts_with("WTP")) %>%
    separate(key, into = c("Sample", "Rep"), sep = "_R")
  # Line plot
  p <- ggplot(df, aes(x = Sample, y = value)) +
    geom_line(aes(color = Rep, group = Rep)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Timepoint", y = "TPM")
  if(InvertRank == TRUE){
    p + ggtitle(sprintf("Bottom %s gene: %s", Rank, df_subset$Geneid))
  } else (
    p + ggtitle(sprintf("Top %s gene: %s", Rank, df_subset$Geneid))
  )

  # # Example
  # SingleGene(1, "A03", InvertRank = FALSE)
}


# Multi-gene function
MultiGene <- function(From=1, To=10, InvertRank=FALSE, Genotype=NULL, SD=TRUE){
  Genes <- CircadianGenes
  # Check Genotype is provided
  if(!is.null(Genotype)){
    Genes <- Genes[Genes$Genotype == Genotype,]
  }
  # Check InvertRank is TRUE
  if(InvertRank == TRUE){
    Genes <- Genes[order(Genes$JTK_BH.Q, decreasing = TRUE),]
  }
  Genes <- Genes[From:To,] # Subset the data
  geneIDs <- Genes$Geneid # Get the gene IDs
  df_subset <- all_TPM[which(all_TPM$Geneid %in% geneIDs),] # Subset the data
  df_subset <- df_subset[!duplicated(df_subset$Geneid), ] # Remove duplicated rows
  # Convert to long format
  df_subset_long <- df_subset %>%
    gather(key, value, starts_with("WTP")) %>%
    separate(key, into = c("Sample", "Rep"), sep = "_R")
  # Prepare to plot
  if(SD == TRUE){
    # Get means and sd
    df_subset_summary <- df_subset_long %>%
      group_by(Geneid, Sample) %>%
      summarise(mean_value = mean(value),
                sd_value = sd(value))
    # Line plot with standard deviation bars
    ggplot(df_subset_summary, aes(x = Sample, y = mean_value, group = Geneid, color = Geneid)) +
      geom_line() +
      geom_point() +
      geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), width = 0.2) +
      labs(title = "Mean values with standar deviation bars across replicates for each WTP",
          x = "Timepoint",
          y = "Mean TPM") +
      theme_minimal()
    } else {
      # Get means only
      df_subset_summary <- df_subset_long %>%
        group_by(Geneid, Sample) %>%
        summarise(mean_value = mean(value),
                  sd_value = sd(value))
      # Line plot without standard deviation bars
      ggplot(df_subset_summary, aes(x = Sample, y = mean_value, group = Geneid, color = Geneid)) +
        geom_line() +
        geom_point() +
        labs(title = "Mean values across replicates for each WTP",
            x = "Timepoint",
            y = "Mean TPM") +
        theme_minimal()
  }

  # # Example
  # MultiGene(From=1, To=10, InvertRank=TRUE, Genotype=NULL, SD=TRUE)
}




# UI
ui <- fluidPage(
  titlePanel("Gene Plotting App"),
  sidebarLayout(
    sidebarPanel(
      selectInput("Function", "Select Function", choices = c("SingleGene", "MultiGene")),
      numericInput("rank", "Rank:", 1),
      checkboxInput("invertRank", "Invert Rank", FALSE),
      textInput("genotype", "Genotype (optional):", ""),
      numericInput("from", "From:", 1),
      numericInput("to", "To:", 10),
      checkboxInput("sd", "Include Standard Deviation", TRUE),
      actionButton("plotButton", "Generate Plot")
    ),
    mainPanel(
      plotOutput("genePlot")
    )
  )
)

# Server
server <- function(input, output) {
  geneData <- eventReactive(input$plotButton, {
    if (input$Function == "SingleGene") {
      SingleGene(input$rank, input$invertRank, if (input$genotype == "") NULL else input$genotype)
    } else {
      MultiGene(input$from, input$to, input$invertRank, if (input$genotype == "") NULL else input$genotype, input$sd)
    }

  })

  output$genePlot <- renderPlot({
    req(geneData())  # Ensure geneData is available
    geneData()
  })
}

# Run the app
shinyApp(ui, server)

