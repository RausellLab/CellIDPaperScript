library(tidyverse)
Bench1 <- readRDS("data/IntestinalBenchCell.rds")%>%  spread(methods, value) %>%  rename(`Reference Data` = data)
Bench2 <- readRDS("data/IntestinalBenchOverall.rds")%>%  spread(methods, value) %>%  rename(`Reference Data` = data)
xlsx::write.xlsx(Bench1, file = "../FinalTable/SupTable5.xlsx", sheetName = "Cell_Population", append = T)
xlsx::write.xlsx(Bench2, file = "../FinalTable/SupTable5.xlsx", sheetName = "Overall", append = T)

