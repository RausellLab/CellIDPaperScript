library(tidyverse)
PancreasBenchCell <- readRDS("data/PancreasBenchCell.rds") %>%  mutate(value = round(value, digits =3)) %>%  spread(methods, value) %>%  as_tibble() %>%  rename(query_data = data)
PancreasBenchOverall <- readRDS("data/PancreasBenchOverall.rds") %>%  mutate(value = round(value, digits =3)) %>%  spread(methods, value) %>%  as_tibble() %>%  rename(query_data = data)
xlsx::write.xlsx(x = PancreasBenchCell, file = "../FinalTable/SupTable9.xlsx", sheetName = "AssessmentPerCellType_Pancreas", append = F)
xlsx::write.xlsx(x = PancreasBenchOverall, file = "../FinalTable/SupTable9.xlsx", sheetName = "OverallAssessment_Pancreas", append = T)


EpithelialBenchCell <- readRDS("data/EpithelialBenchCell.rds") %>% mutate(value = round(value, digits =3)) %>%  spread(methods, value) %>%  as_tibble() %>%  rename(query_data = data)
EpithelialBenchOverall <- readRDS("data/EpithelialBenchOverall.rds") %>% mutate(value = round(value, digits =3)) %>%   spread(methods, value) %>%  as_tibble() %>%  rename(query_data = data)
xlsx::write.xlsx(x = EpithelialBenchCell, file = "../FinalTable/SupTable9.xlsx", sheetName = "AssessmentPerCellType_Epithelial", append = T)
xlsx::write.xlsx(x = EpithelialBenchOverall, file = "../FinalTable/SupTable9.xlsx", sheetName = "OverallAssessment_Epithelial", append = T)
