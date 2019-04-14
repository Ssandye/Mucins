library(tidyverse)
library(readxl)
library(ggrepel)
library(data.table)
library(RColorBrewer)

MassCal <- function(sequence){
  for (sequence_number in 1:length(sequence)) {
    peptide_vector <- strsplit(sequence[sequence_number], 
                               split = "")[[1]]
    peptide_length <- length(peptide_vector)
    C <- 12.010736
    H <- 1.007941
    O <- 15.999405
    S <- 32.064787
    N <- 14.006703
    proton <- 1.007276466
    electron <- 0.00054857990943
    residueMass <- function(residue) {
      if (residue == "A") 
        mass = C * 3 + H * 5 + N + O
      if (residue == "R") 
        mass = C * 6 + H * 12 + N * 4 + O
      if (residue == "N") 
        mass = C * 4 + H * 6 + N * 2 + O * 2
      if (residue == "D") 
        mass = C * 4 + H * 5 + N + O * 3
      if (residue == "E") 
        mass = C * 5 + H * 7 + N + O * 3
      if (residue == "Q") 
        mass = C * 5 + H * 8 + N * 2 + O * 2
      if (residue == "G") 
        mass = C * 2 + H * 3 + N + O
      if (residue == "H") 
        mass = C * 6 + H * 7 + N * 3 + O
      if (residue == "I") 
        mass = C * 6 + H * 11 + N + O
      if (residue == "L") 
        mass = C * 6 + H * 11 + N + O
      if (residue == "K") 
        mass = C * 6 + H * 12 + N * 2 + O
      if (residue == "M") 
        mass = C * 5 + H * 9 + N + O + S
      if (residue == "m") 
        mass = C * 5 + H * 9 + N + O * 2 + S
      if (residue == "F") 
        mass = C * 9 + H * 9 + N + O
      if (residue == "P") 
        mass = C * 5 + H * 7 + N + O
      if (residue == "S") 
        mass = C * 3 + H * 5 + N + O * 2
      if (residue == "T") 
        mass = C * 4 + H * 7 + N + O * 2
      if (residue == "W") 
        mass = C * 11 + H * 10 + N * 2 + O
      if (residue == "Y") 
        mass = C * 9 + H * 9 + N + O * 2
      if (residue == "V") 
        mass = C * 5 + H * 9 + N + O
      if (residue == "C") 
        mass = C * 5 + H * 8 + N * 2 + O * 2 + S
      return(mass)
    }
    masses <- sapply(peptide_vector, residueMass)
    pm <- sum(masses, 2 * H, O)
    return(pm)
  }
}

FullPlot <- function(File){
  p <-   ggplot(File, aes(x = Mass, y = Intensity)) + 
    geom_bar(stat = "identity",width = diff(range(File$Mass)) / 500, fill = "black") + 
    theme_bw() + ggtitle("Deconvoluted MS1 spectra")
  return(p)
}

AnnoPlot <- function(File){
  p <-   ggplot(File, aes(x = Mass, y = Intensity)) + 
    geom_bar(stat = "identity",width = diff(range(File$Mass)) / 500, fill = "black") + theme_bw() +
    geom_label_repel(aes(label = Assignment, color = Assigned), box.padding = unit(0.3, "lines"), 
                     point.padding = unit(0.4, "lines"), show.legend = F, size = 3) + 
    scale_color_manual(values = c("red", "blue"))
  return(p)
}