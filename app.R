#
# Jenna Everard
#
# Last Updated: June 8, 2026
#

library(shiny)
library(zip)
library(ggplot2)
library(dplyr)
library(tidyr)

# Function to read XRD output file with BTX style header
read_xrd_file <- function(filepath) {
  
  lines <- readLines(filepath, warn = FALSE)
  header_lns <- lines[startsWith(lines, "#")]
  data_lns <- lines[!startsWith(lines, "#")]
  data <- tryCatch(
    {
      read.table(text = paste(data_lns, collapse="\n"),
                 header = FALSE,
                 col.names = c("two_theta", "intensity"))
    },
    error = function(e) {NULL}
  )
  
  list(header = header_lns, data = data)
}

# Function to read XRF output file with BTX style header
read_xrf_file <- function(filepath) {
  
  lines <- readLines(filepath, warn = FALSE)
  header_lns <- lines[startsWith(lines, "#")]
  data_lns <- lines[!startsWith(lines, "#")]
  data <- tryCatch(
    {
      read.table(text = paste(data_lns, collapse="\n"),
                 header = FALSE,
                 col.names = c("index", "Energy", "Intensity"))
    },
    error = function(e) {NULL}
  )
  
  list(header = header_lns, data = data)
}

# Function to save background corrected output
#save

# Define UI
ui <- fluidPage(

    # Application title
    titlePanel("BTX II XRD/XRF: Background Subtraction"),

    # Sidebar with inputs
    sidebarLayout(
        sidebarPanel(
          
          # Input data run as a ZIP (direct output from BTX)
          fileInput(
            inputId="datazipfile",
            label="Upload your data as a ZIP file",
            accept=".zip"
          ),
          
          # Input blank run as a ZIP (direct output from BTX)
          fileInput(
            inputId="blankzipfile",
            label="Upload your blank/background as a ZIP file",
            accept=".zip"
          ),
          
          # Prompt for a file name to name all files to save
          textInput(
            inputId="filename",
            label="Sample Name (for naming output files): "
          ),
          
          radioButtons(
            "units",
            "X-axis units for XRD display plot (* Feature coming soon!)",
            choices = c("Angstrom"="Angstrom", "2Theta"="2Theta"),
            selected="2Theta",
            inline=TRUE
          )
        ),

        mainPanel(
          plotOutput("XRD_both_patterns"),
          plotOutput("XRD_background_corrected"),
          plotOutput("XRF_both_patterns"),
          plotOutput("XRF_background_corrected")
        )
    )
)

# Define server logic
server <- function(input, output) {

    plots <- reactive({
      
      req(input$datazipfile, input$blankzipfile, input$filename)
      
      data_fp <- input$datazipfile$datapath
      blank_fp <- input$blankzipfile$datapath
      
      # Verify that both uploaded files are actually zip files
      if(tolower(tools::file_ext(input$datazipfile$name)) != "zip") {
        showNotification("The uploaded data file is not a ZIP file", type="error")
        return()
      }
      if(tolower(tools::file_ext(input$blankzipfile$name)) != "zip") {
        showNotification("The uploaded blank file is not a ZIP file", type="error")
        return()
      }
      
      # Create temporary directories for unzipped inputs
      data_dir <- tempfile("data_")
      blank_dir <- tempfile("blank_")
      dir.create(data_dir)
      dir.create(blank_dir)
      
      # Unzip inputs into temporary directories
      unzip(data_fp, exdir = data_dir)
      unzip(blank_fp, exdir = blank_dir)
      
      # Quick lil function to find files :)
      find_file <- function(base, pattern) {
        all_files <- list.files(base, pattern=pattern, recursive=TRUE, full.names=TRUE)
        if(length(all_files) == 0) return(NULL)
        all_files[1]
      }
      
      # find the data files
      data_txt <- find_file(data_dir, "-film.txt$")
      blank_txt <- find_file(blank_dir, "-film.txt$")
      data_xrf_txt <- find_file(data_dir, "-xrf.txt$")
      blank_xrf_txt <- find_file(blank_dir, "-xrf.txt$")
      
      # If any of the data files are missing, CRASH
      if(any(sapply(list(data_txt, blank_txt, data_xrf_txt, blank_xrf_txt), is.null))) {
        showNotification("ZIP folder is missing expected files.", type="error")
        return()
      }
      
      # Try to read the TXT files
      data <- read_xrd_file(data_txt)
      blank <- read_xrd_file(blank_txt)
      data_xrf <- read_xrf_file(data_xrf_txt)
      blank_xrf <- read_xrf_file(blank_xrf_txt)
      
      if(is.null(data$data) || is.null(blank$data) || is.null(data_xrf$data) || is.null(blank_xrf$data)) {
        showNotification("Unable to read data files", type="error")
      }
      
      #######
      # XRD #
      #######
      
      # Interpolate both onto the same two-theta grid
      min_twotheta <- max(min(data$data$two_theta), min(blank$data$two_theta))
      max_twotheta <- min(max(data$data$two_theta), max(blank$data$two_theta))
      common_twotheta <- seq(min_twotheta, max_twotheta, by=0.05)
      
      data_interp <- approx(data$data$two_theta, data$data$intensity, xout=common_twotheta)$y
      blank_interp <- approx(blank$data$two_theta, blank$data$intensity, xout=common_twotheta)$y
      
      # Subtract blank from data
      corrected <- data.frame(two_theta = common_twotheta, intensity = data_interp - blank_interp)
      
      # Plot 1: XRD data and XRD blank
      df1 <- tibble(
        x = common_twotheta,
        data = data_interp,
        blank = blank_interp
      ) %>%
        pivot_longer(cols = c(data, blank),
                     names_to = "Series",
                     values_to = "Value")
      
      p1 <- ggplot(df1, aes(x = x, y = Value, color = Series)) +
        geom_line(size = 1) +
        labs(title = "Raw XRD Data",
             x = "Two-Theta",
             y = "Intensity") +
        theme_minimal()
      
      # Plot 2: Background-corrected XRD
      df2 <- tibble(
        x = common_twotheta,
        data = data_interp - blank_interp
      )
      
      p2 <- ggplot(df2, aes(x = x, y = data)) +
        geom_line(size = 1) +
        labs(title = "Background Subtracted XRD Data",
             x = "Two-Theta",
             y = "Intensity") +
        theme_minimal()
      
      #######
      # XRF #
      #######
      
      # Interpolate both onto the same two-theta grid
      min_energy <- max(min(data_xrf$data$Energy), min(blank_xrf$data$Energy))
      max_energy <- min(max(data_xrf$data$Energy), max(blank_xrf$data$Energy))
      common_energy <- seq(min_energy, max_energy, by=0.0036437)
      
      data_xrf_interp <- approx(data_xrf$data$Energy, data_xrf$data$Intensity, xout=common_energy)$y
      blank_xrf_interp <- approx(blank_xrf$data$Energy, blank_xrf$data$Intensity, xout=common_energy)$y
      
      # Crop energy range
      mask <- common_energy >= 3 & common_energy <= 8
      common_energy <- common_energy[mask]
      data_xrf_interp <- data_xrf_interp[mask]
      blank_xrf_interp <- blank_xrf_interp[mask]
      
      # Subtract blank from data
      corrected_xrf <- data.frame(energy = common_energy, intensity = data_xrf_interp - blank_xrf_interp)
      
      # Plot 3: XRF data and XRF blank
      df3 <- tibble(
        x = common_energy,
        data = data_xrf_interp,
        blank = blank_xrf_interp
      ) %>%
        pivot_longer(cols = c(data, blank),
                     names_to = "Series",
                     values_to = "Value")
      
      p3 <- ggplot(df3, aes(x = x, y = Value, color = Series)) +
        geom_line(size = 1) +
        labs(title = "Raw XRF Data",
             x = "Energy (keV)",
             y = "Intensity (au)") +
        theme_minimal()
      
      # Plot 4: Background-corrected XRF
      df4 <- tibble(
        x = common_energy,
        data = data_xrf_interp - blank_xrf_interp
      )
      peak_lines <- tibble(
        element = c("K", "Ca", "Ti", "Mn", "Fe"),
        energy = c(3.31, 3.69, 4.51, 5.90, 6.40)
      )
      
      p4 <- ggplot(df4, aes(x = x, y = data)) +
        geom_line(size = 1) +
        geom_vline(data = peak_lines, aes(xintercept = energy),
                   linetype = "dashed", color = "seagreen", alpha = 0.7) +
        geom_text(data = peak_lines, 
                  aes(x = energy, y = max(df4$data, na.rm=TRUE) * 0.95,
                      label = element),
                  angle = 90, vjust = -0.4, hjust = 0, size = 5, color = "seagreen") +
        labs(title = "Background Subtracted XRF Data",
             x = "Energy (keV)",
             y = "Intensity (au)") +
        theme_minimal()
      
      list(p1 = p1, p2 = p2, p3 = p3, p4 = p4)
      
    })
    
    # Push all four plots to the output
    
    output$XRD_both_patterns <- renderPlot({plots()$p1})
    
    output$XRD_background_corrected <- renderPlot({plots()$p2})
    
    output$XRF_both_patterns <- renderPlot({plots()$p3})
    
    output$XRF_background_corrected <- renderPlot({plots()$p4})
}

# Run the application 
shinyApp(ui = ui, server = server)
