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

MUC1 <- "DIGLNGIPAPDNKPAPGSTAPPAHGVTSAPDTRPAPGSTAPPAHGVTSAPDTRPAPGSTAPPAHGVTSAPDTRPAPGSTAPPAHGVTSAPDTRPAPGSTAPPAHGVTSAPDTRPAPGSTAPPAHGVTSAPDTRPAPGSTAPPAHGVTSAAAHHHHHH"

MUC1 <- "PAPGSTAPPAHGVTSAPDTRPAPGSTAPPAHGVTSAPDTRPAPGSTAPPAHGVTSAPDTRPAPGSTAPPAHGVTSAPDTRPAPGSTAPPAHGVTSAPDTRPAPGSTAPPAHGVTSAPDTRPAPGSTAPPAHGVTSAAAHHHHHH"

AMass <- MassCal(MUC1)


#Hex <- 162.140848

# dHex <- 146.141443
# 
# HexNAc <- 203.192845
# 
# HexHexNAc <- 203.192845 + 162.140848
# 
# NeuAc <- 291.25503
# 
# SO3 <- 80.063002

HexNum <- as.data.frame(c(0:20))
GlyNum <- HexNum
names(GlyNum) <- "Hex"

AddSugar <- function(MinGly, MaxGly, GlyName){
  GN1 <- as.data.frame(rep(c(MinGly:MaxGly), each = length(GlyNum[,1])))
  names(GN1) <- GlyName
  GlyNum <- cbind(GlyNum,GN1)
  return(GlyNum)
}

GlyNum <- AddSugar(0,34,"HexNAc")
GlyNum <- AddSugar(0,34,"NeuAc")
GlyNum <- AddSugar(0,20,"dHex")
GlyNum <- AddSugar(0,1,"Phos")

GlyMass <- data.table(GlyNum)%>%
  .[, Average := Hex * 162.140848 + dHex * 146.141443 + HexNAc * 203.192845 + NeuAc * 291.25503 + Phos * 79.979917]%>%
  .[, Structure := str_c("Hex(", Hex, ")HexNAc(", HexNAc, ")NeuAc(", NeuAc, ")dHex(", dHex,")Phos(", Phos, ")")]%>%
  .[, Structure := str_replace(Structure, "NeuAc\\(0\\)", "")]%>%
  .[, Structure := str_replace(Structure, "dHex\\(0\\)", "")]%>%
  .[, Structure := str_replace(Structure, "Phos\\(0\\)", "")]%>%
  .[, Structure := str_replace(Structure, "HexNAc\\(0\\)", "")]%>%
  .[, Structure := str_replace(Structure, "Hex\\(0\\)", "")]

GlyMass <- GlyMass[dHex == 0 & Phos == 0 & Hex == 0]


## MUC1_T

MUC1_T <- read_excel("data/300319_Lingbo_Muc1_T(v2).xlsx", 
                     sheet = "Lingbo-TdeSA(14-17)")%>%data.table()

MUC1_T <- MUC1_T[,.(Mass,Intensity,Name)]%>%
  .[order(Intensity, decreasing = TRUE)]

MUC1_T[,MassDif := Mass - MUC1_T[[1,1]]]


FullPlot(MUC1_T[1:25])


Temp_MUC <- GlyMass[abs(Average - abs(MUC1_T[[1,1]] - AMass)) <= 1]

MUC1_T[, Assignment := round(MassDif)%>%as.character()]

MUC1_T$A2 <- MUC1_T$Assignment

MUC1_T[[1,5]] <- Temp_MUC[2,Structure]

MUC1_T[[2,5]] <- "-1 dHex"

MUC1_T[[3,5]] <- "-1 HexHexNAc"

MUC1_T[[4,5]] <- "-1 HexHexNAc -1 dHex"

MUC1_T[[5,5]] <- "+1 dHex"

MUC1_T[,Assigned := !(Assignment == A2)]

AnnoPlot(MUC1_T)


## MUC1 ST

MUC1_ST <- read_excel("data/300319_Lingbo_Muc1_T(v2).xlsx", 
                     sheet = "Lingbo-T(14-18)")%>%data.table()

MUC1_ST <- MUC1_ST[,.(Mass,Intensity,Name)]%>%
  .[order(Intensity, decreasing = TRUE)]

MUC1_ST[,MassDif := Mass - MUC1_ST[[1,1]]]

FullPlot(MUC1_ST[1:25])


MUC1_ST[, Assignment := round(MassDif)%>%as.character()]

MUC1_ST$A2 <- MUC1_ST$Assignment

GlyMass[abs(Average - abs(MUC1_ST[[1,1]] - AMass)) <= 0.5]

Temp_MUC <- GlyMass[abs(Average - abs(MUC1_ST[[1,1]] - AMass)) <= 0.5]

MUC1_ST[[1,5]] <- Temp_MUC[1,Structure]

MUC1_ST[[2,5]] <- "HexHexNAc(34)NeuAc(1)"

MUC1_ST[[3,5]] <- "HexHexNAc(32)NeuAc(3)dHex(4)"

MUC1_ST[[4,5]] <- "HexHexNAc(32)NeuAc(2)dHex(4)"

MUC1_ST[[5,5]] <- "HexHexNAc(34)NeuAc(3)"

MUC1_ST[[6,5]] <- "HexHexNAc(33)NeuAc(2)dHex(1)"

MUC1_ST[[7,5]] <- "HexHexNAc(34)NeuAc(3)dHex(2)"

MUC1_ST[[8,5]] <- "HexHexNAc(33)NeuAc(1)dHex(2)"

MUC1_ST[[9,5]] <- "HexHexNAc(32)NeuAc(2)dHex(1)"

GlyMass[abs(Average + AMass - MUC1_ST[[1,1]])/ (Average + AMass) * 1000000 <= 3]

MUC1_ST[,Assigned := !(Assignment == A2)]

AnnoPlot(MUC1_ST)

GlyMass[abs(Average + AMass - 27277.9)/ (Average + AMass) * 1000000 <= 10]
