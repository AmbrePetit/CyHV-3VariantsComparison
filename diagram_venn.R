# load Venn diagram package
library("grid")
library("futile.logger")
library("VennDiagram")
library(dplyr)

data <- read.csv("/Users/ambre/Downloads/vcf-compare-VN.txt", header = FALSE, sep = "\t")


####### Variants unique à chaque échantillon : 
# On récupère seulement les lignes de l'échantillon tout seul c'est-à-dire seulement les lignes où il y a les deux dernières colonnes sont vide
data_one_set <- data[data[, ncol(data) - 1] == "" & data[, ncol(data)] == "",  ]
# Extraction des nombres pour l'échantillon seul 
number_one_set <- as.numeric(gsub("\\D", "", data_one_set[, 1]))
# Extraction des noms de la deuxieme colonne 
name_one_set <- data_one_set[,2]
# Associer les nombres aux noms 
df_one_set <- data.frame(Number_one_set = number_one_set, Name_one_set = name_one_set)

#Variables comprenant le nombre pour chaque échantillon unique
P15_only <- grep("^P15", df_one_set$Name_one_set) 
Number_P15_only <- df_one_set$Number_one_set[P15_only]
P30_only <- grep("^P30", df_one_set$Name_one_set) 
Number_P30_only <- df_one_set$Number_one_set[P30_only]
P50_only <- grep("^P50", df_one_set$Name_one_set) 
Number_P50_only <- df_one_set$Number_one_set[P50_only]


####### Variants communs à deux échantillons : 
# On récupère seulement les lignes où il y a la dernières colonnes vide
last_column_empty <- data[, ncol(data)] == ""
data_two_set <- data %>%
  filter(last_column_empty) %>%
  anti_join(data %>%
              filter(last_column_empty & data[, ncol(data) - 1] == ""), by = NULL)


# Extraction des nombres pour les variants communs à deux échantillons :
number_two_set <- as.numeric(gsub("\\D", "", data_two_set[, 1]))
# Extraction des noms de la deuxieme colonne : 
name_two_set <- data_two_set[,2]
# Extraction des noms de la troisème colonne : 
name_other_two_set <- data_two_set[,3]

# Associer les nombres aux noms 
df_two_set <- data.frame(Number_two_set = number_two_set, Name_two_set = name_two_set, Name_two_set2 = name_other_two_set)
#Variables comprenant le nombre pour chaque échantillon communs deux à deux : 
P15_P30 <- df_two_set[grepl("^P15", df_two_set$Name_two_set) & grepl("^P30", df_two_set$Name_two_set2), ]
Number_P15_P30 <- P15_P30$Number_two_set
P15_P50 <- df_two_set[grepl("^P15", df_two_set$Name_two_set) & grepl("^P50", df_two_set$Name_two_set2), ]
Number_P15_P50 <- P15_P50$Number_two_set
P30_P50 <- df_two_set[grepl("^P30", df_two_set$Name_two_set) & grepl("^P50", df_two_set$Name_two_set2), ]
Number_P30_P50 <- P30_P50$Number_two_set


####### Variants communs aux trois échantillons : 
# On récupère seulement les lignes où il n'y a aucune colonne vide
data_three_set <- function(row) {
  all(sapply(row, function(x) !is.na(x) && nchar(x) > 0))
}
# Appliquer la fonction à chaque ligne
data_three_set <- data[apply(data, 1, data_three_set), ]

# Extraction des nombres pour les variants communs tous les échantillons :
number_three_set <- as.numeric(gsub("\\D", "", data_three_set[, 1]))
# Extraction des noms de la deuxieme colonne : 
name_three_set <- data_three_set[,2]
# Extraction des noms de la troisème colonne : 
name_other_three_set <- data_three_set[,3]
# Extraction des noms de la dernière colonne : 
name_last_three_set <- data_three_set[,4]

# Associer les nombres aux noms 
df_three_set <- data.frame(Number_three_set = number_three_set, Name_three_set = name_three_set, Name_three_set2 = name_other_three_set, Name_three_set3 = name_last_three_set)
#Variables comprenant le nombre pour chaque échantillon communs à tous : 
P15_P30_P50 <- df_three_set[grepl("^P15", df_three_set$Name_three_set) & grepl("^P30", df_three_set$Name_three_set2) & grepl("^P50", df_three_set$Name_three_set3), ]
Number_P15_P30_P50 <- P15_P30_P50$Number_three_set


Number_P15 <- sum(Number_P15_only, Number_P15_P30, Number_P15_P50, Number_P15_P30_P50, na.rm = TRUE)
Number_P30 <- sum(Number_P30_only, Number_P30_P50, Number_P15_P30, Number_P15_P30_P50, na.rm = TRUE)
Number_P50 <- sum(Number_P50_only, Number_P30_P50, Number_P15_P50, Number_P15_P30_P50, na.rm = TRUE)


# Création du diagramme de Venn
grid.newpage()
draw.triple.venn(area1=Number_P15, area2=Number_P30, area3=Number_P50,
                 n12=sum(Number_P15_P30, Number_P15_P30_P50, na.rm = TRUE), n23=sum(Number_P30_P50, Number_P15_P30_P50, na.rm = TRUE), n13=sum(Number_P15_P50, Number_P15_P30_P50, na.rm = TRUE), n123=Number_P15_P30_P50,
                 category=c("P15","P30","P50"),
                 lty = "blank",
                 fill=c("Green","Red","Blue"))







##### Diagramme de Venn pour comparaison des variants caller : 

# load Venn diagram package
library("grid")
library("futile.logger")
library("VennDiagram")
library(dplyr)

# Récupération du nombre total de variants de chaque fichier
nanosv <- "/Users/ambre/Desktop/MASTER/M2/Bioanalyse, transcriptomique/combiSV/P15-2.trimed1000.sv_nanosv.vcf"
lignes_nanosv <- readLines(nanosv)
area2 <- sum(grepl("^DQ657948.1", lignes_nanosv))

cutesv <- "/Users/ambre/Desktop/MASTER/M2/Bioanalyse, transcriptomique/combiSV/P15-2.trimed1000.sv_cutesv.vcf"
lignes_cutesv <- readLines(cutesv)
area1 <- sum(grepl("^DQ657948.1", lignes_cutesv))

svim <- "/Users/ambre/Desktop/MASTER/M2/Bioanalyse, transcriptomique/combiSV/P15-2.trimed1000.sv_svim.vcf"
lignes_svim <- readLines(svim)
area4 <- sum(grepl("^DQ657948.1", lignes_svim))

nanosv <- "/Users/ambre/Desktop/MASTER/M2/Bioanalyse, transcriptomique/combiSV/P15-2.trimed1000.sv_sniffles.vcf"
lignes_sniffles <- readLines(nanosv)
area3 <- sum(grepl("^DQ657948.1", lignes))


# Création du diagramme de Venn
grid.newpage()
draw.quad.venn(area1=area1, area2=area2, area3=area3, area4=area4, n12=3, n13=2, n14=0, n23=0, n24=0,
               n34=0, n123=1, n124=0, n134=0, n234=0, n1234=0, 
               category=c("CuteSV","NanoSV","Sniffles", "SVIM"),
               lty = "blank",
               fill=c("Green","Red","Blue","Yellow"))












