# Path analysis code to test more complex versions than simple mediations

# Set working directory
setwd("D:/R_code")

# if you need to install the library: install.packages("readxl")
# activate excel reader
library("readxl")

# load data
df1 <- read_excel("df_path_analysis_plasma.xlsx",sheet = "Sheet1",col_types = "numeric")

# if you need to install the library: #install.packages("lavaan")
# activate lavaan for path analysis
library(lavaan)


## multiple models to be tested, activate only one at a time, these are all different path analysis

# Specify model, can also write it in a single line
specmod <- "
dmpfc_lac ~ plasma_lac
bold_dmpfc ~ dmpfc_lac
PHE ~ dmpfc_lac
PHE ~ bold_dmpfc
"

specmod <- "
dmpfc_lac ~ plasma_lac
bold_dmpfc ~ dmpfc_lac
kEp ~ dmpfc_lac
kEp ~ bold_dmpfc
"

specmod <- "
dmpfc_lac ~ plasma_lac
bold_dmpfc ~ dmpfc_lac
MHE ~ dmpfc_lac
MHE ~ bold_dmpfc
"

specmod <- "
dmpfc_lac ~ plasma_lac
bold_dmpfc ~ dmpfc_lac
kEm ~ dmpfc_lac
kEm ~ bold_dmpfc
"

specmod <- "
dmpfc_lac ~ plasma_lac
bold_dmpfc ~ dmpfc_lac
THE ~ dmpfc_lac
THE ~ bold_dmpfc
"

specmod <- "
dmpfc_lac ~ plasma_lac
bold_dmpfc ~ dmpfc_lac
PHE ~ dmpfc_lac
PHE ~ bold_dmpfc
PHE ~ plasma_lac
"

# Estimate our model
fitmod <- sem(specmod,data=df)

# Summarize the results since it is degree of freedom = 0, it is just identified, so fit.measures won't give more
summary(fitmod,fit.measures=TRUE, standardized = TRUE, rsquare=TRUE)
