Bench1 <- readRDS("data/ATACBenchCell.rds")%>% mutate(value = round(value, digits =3)) %>%   spread(methods, value) %>%  rename(`Reference Data` = data)
Bench2 <- readRDS("data/ATACBenchOverall.rds")%>% mutate(value = round(value, digits =3)) %>%   spread(methods, value) %>%  rename(`Reference Data` = data)
xlsx::write.xlsx(Bench1, file = "../FinalTable/SupTable7.xlsx", sheetName = "Cell_Population", append = T)
xlsx::write.xlsx(Bench2, file = "../FinalTable/SupTable7.xlsx", sheetName = "Overall", append = T)
