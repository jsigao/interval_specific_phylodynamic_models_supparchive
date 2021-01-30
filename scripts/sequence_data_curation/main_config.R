
library(RColorBrewer)

reference_genome_name_string <- "EPI_ISL_402125"
nuc_unambig <- c("a", "c", "g", "t")

# coarse discretization
regions_coarse <- list(ChinaMainland = c("China", "China_Hubei", 
                                         "China_Guangdong", "China_Hunan",
                                         "China_Anhui", "China_Fujian", "China_Jiangsu", "China_Jiangxi", "China_Shanghai", "China_Zhejiang",
                                         "China_Beijing", "China_Heilongjiang", "China_Shandong", 
                                         "China_Chongqing", "China_Henan", "China_Sichuan", "China_Yunnan"), 
                       AsiaEast = c("Hong Kong", "Japan", "South Korea", "Taiwan"),
                       AsiaSouthSoutheast = c("India", "Nepal", "Pakistan", "Sri Lanka",
                                              "Brunei", "Cambodia", "Indonesia", "Malaysia", 
                                              "Philippines", "Singapore", "Thailand", "Vietnam"), 
                       MiddleEastAfrica = c("Bahrain", "Georgia", "Egypt", "Iran", "Israel", "Jordan", "Kuwait", "Lebanon", "Oman", 
                                            "Saudi Arabia", "Turkey", "United Arab Emirates", "Palestine",
                                            "Algeria", "Democratic Republic of the Congo", "Gambia", "Morocco", "Nigeria", "Senegal", "South Africa", "Tunisia"),
                       
                       EuropeNorth = c("Denmark", "Finland", "Iceland", "Norway", "Sweden", 
                                       "Ireland", "United Kingdom"),
                       EuropeWestCentralEast = c("France", "Belgium", "Netherlands",
                                                 "Austria", "Croatia", "Czech Republic", "Germany", "Hungary", "Luxembourg", 
                                                 "Poland", "Slovakia", "Slovenia", "Switzerland",
                                                 "Greece", "Lithuania", "Latvia", "Russia"),
                       
                       Italy = c("Italy", "Italy"),
                       PortugalSpain = c("Portugal", "Spain"),
                       
                       CanadaUSAWest = c("Canada_British Columbia", "Canada_Manitoba", "Canada_Saskatchewan", "USA_Alaska",
                                         "USA_Washington", "USA_Oregon", "USA_California", "USA_Idaho", "USA_Nevada", "USA_Utah",
                                         "USA_Colorado", "USA_Wyoming", "USA_Arizona", "USA_Texas", "USA_Diamond Princess", "USA_Grand Princess",
                                         "USA_Grand Princess", "USA_Hawaii"), 
                       CanadaUSAEast = c("Canada_New Brunswick", "Canada_Newfoundland and Labrador", "Canada_Nova Scotia",
                                         "Canada_Ontario", "Canada_Prince Edward Island", "Canada_Quebec",
                                         "USA_New York", "USA_Connecticut", "USA_Massachusetts", "USA_Vermont",
                                         "USA_New Hampshire", "USA_New Jersey", "USA_Pennsylvania", "USA_Rhode Island", 
                                         "USA_Illinois", "USA_Indiana", "USA_Iowa", "USA_Kansas", "USA_Minnesota", "USA_Missouri",
                                         "USA_Nebraska", "USA_Ohio", "USA_Wisconsin", "USA_Washington DC",
                                         "USA_District of Columbia", "USA_Florida", "USA_Georgia", "USA_Louisiana", "USA_Mississippi",
                                         "USA_Maryland", "USA_North Carolina", "USA_South Carolina", "USA_Virginia"),
                       LatinAmerica = c("Mexico", "Costa Rica", "Panama", "USA_Puerto Rico", "Belize",
                                        "Argentina", "Brazil", "Chile", "Colombia", "Ecuador", "Peru", "Saint Barthélemy"),
                       
                       Oceania = c("Australia", "New Zealand"))

# fine discretization
regions_fine <- list(ChinaHubei = c("China_Hubei", "China_Hubei"),
                     ChinaSouth = c("China_Guangdong", "China_Hunan"),
                     ChinaEast = c("China_Anhui", "China_Fujian", "China_Jiangsu", "China_Jiangxi", "China_Shanghai", "China_Zhejiang"),
                     ChinaNorthWest = c("China_Beijing", "China_Heilongjiang", "China_Henan", "China_Shandong", 
                                        "China_Chongqing", "China_Sichuan", "China_Yunnan"),
                      
                      HongKongTaiwan = c("Hong Kong", "Taiwan"),
                      JapanKorea = c("Japan", "South Korea"),
                      
                      AsiaSoutheast = c("Brunei", "Cambodia", "Indonesia", "Malaysia", 
                                        "Philippines", "Singapore", "Thailand", "Vietnam"),
                      AsiaSouth = c("India", "Nepal", "Pakistan", "Sri Lanka"), 
                      
                      MiddleEast = c("Bahrain", "Georgia", "Egypt", "Iran", "Israel", "Jordan", "Kuwait", 
                                     "Lebanon", "Oman", "Saudi Arabia", "Turkey", "United Arab Emirates", "Palestine"),
                      Africa = c("Algeria", "Democratic Republic of the Congo", "Gambia", "Morocco", "Nigeria", "Senegal", "South Africa", "Tunisia"),
                      
                      NorthEuropeContinental = c("Denmark", "Finland", "Norway", "Sweden"),
                      Iceland = c("Iceland", "Iceland"),
                      
                      BritishIsles = c("Ireland", "United Kingdom"),
                      
                      EuropeWest = c("France", "Belgium", "Netherlands"),
                      
                      EuropeCentralEast = c("Austria", "Croatia", "Czech Republic", "Germany", "Hungary", "Luxembourg", 
                                        "Poland", "Slovakia", "Slovenia", "Switzerland", "Greece", "Lithuania", "Latvia", "Russia"),
                      
                      Italy = c("Italy", "Italy"),
                      PortugalSpain = c("Portugal", "Spain"),
                     
                     Washington = c("USA_Washington", "USA_Washington"),
                     
                     CanadaUSAWest = c("Canada_British Columbia", "Canada_Manitoba", "Canada_Saskatchewan", "USA_Alaska",
                                       "USA_Oregon", "USA_California", 
                                       "USA_Colorado", "USA_Idaho", "USA_Nevada", "USA_Utah", "USA_Wyoming", 
                                       "USA_Arizona", "USA_Texas", 
                                       "USA_Diamond Princess", "USA_Grand Princess", "USA_Hawaii"), 
                     
                     NewYork = c("USA_New York", "USA_New York"),
                     
                     CanadaUSAEast = c("Canada_New Brunswick", "Canada_Newfoundland and Labrador", "Canada_Nova Scotia",
                                       "Canada_Ontario", "Canada_Prince Edward Island", "Canada_Quebec",
                                       "USA_Connecticut", "USA_Massachusetts", "USA_Vermont",
                                       "USA_New Hampshire", "USA_New Jersey", "USA_Pennsylvania", "USA_Rhode Island", 
                                       "USA_Illinois", "USA_Indiana", "USA_Iowa", "USA_Kansas", "USA_Minnesota", "USA_Missouri",
                                       "USA_Nebraska", "USA_Ohio", "USA_Wisconsin", 
                                       "USA_Washington DC", "USA_District of Columbia", "USA_Florida", "USA_Georgia", "USA_Louisiana", "USA_Mississippi",
                                       "USA_Maryland", "USA_North Carolina", "USA_South Carolina", "USA_Virginia"),
                     
                     LatinAmerica = c("Mexico", "Costa Rica", "Panama", "USA_Puerto Rico", "Belize",
                                      "Argentina", "Brazil", "Chile", "Colombia", "Ecuador", "Peru", "Saint Barthélemy"),
                      
                      Oceania = c("Australia", "New Zealand"))

# area names to be printed
state_name_formals <- list(coarse = c("Mainland China", "East Asia", "Southeast and South Asia", "Middle East and Africa", "North Europe", 
                                      "West, Central, and East Europe", "Italy", "Spain and Portugal", "West USA and Canada", "East USA and Canada", 
                                      "Latin America", "Oceania"),
                           fine = c("Hubei China", "South China", "East China", "North and West China", "Hong Kong and Taiwan", "Japan and Korea", 
                                    "Southeast Asia", "South Asia", "Middle East", "Africa", "North Europe", "Iceland", 
                                    "British Isles", "West Europe", "Central Europe", "Italy", "Spain and Portugal", "Washington", 
                                    "West USA and Canada", "New York", "East USA and Canada", "Latin America", "Oceania"))

# grouping areas into regions
region_divisions <- list(coarse = list(AsiaPacific = c(1:3, 12), EuropeAfrica = c(4:8), America = c(9:11)),
                         fine = list(AsiaPacific = c(1:8, 23), EuropeAfrica = c(9:17), America = c(18:22)))

# color setting for plotting purposes
state_color_asiapacific_coarse <- rev(brewer.pal("YlOrRd", n = length(region_divisions$coarse$AsiaPacific)))
state_color_asiapacific_fine <- c("#67000D", "#990000", "#CC0000", "#FF0000", "#FF3D00", "#FF6D00", "#FFA500", "#311B92", "#FFCA28")

state_color_europeafrica_coarse <- brewer.pal("BuPu", n = length(region_divisions$coarse$EuropeAfrica))
state_color_europeafrica_fine <- c("#673AB7", "#9C27B0", "#66BB6A", "#CCFF00", "#388E3C", "#C0CA33", "#9E9D24", "#827717", "#33691E")

state_color_america_coarse <- brewer.pal("PiYG", n = 5)[c(1, 5, 4)]
state_color_america_fine <- c("#1A237E", "#42A5F5", "#0D47A1", "#90CAF9", "#3D5AFE")

state_colors <- list(coarse = adjustcolor(c(state_color_asiapacific_coarse[1:3], 
                                            state_color_europeafrica_coarse, state_color_america_coarse, state_color_asiapacific_coarse[4]), alpha.f = 0.9),
                     fine = adjustcolor(c(state_color_asiapacific_fine[1:8], 
                                          state_color_europeafrica_fine, state_color_america_fine, state_color_asiapacific_fine[9]), alpha.f = 0.9))

state_orders <- list(coarse = c(1, 2, 3, 12, 4, 6, 7, 8, 5, 11, 10, 9), 
                     fine = c(1, 2, 3, 4, 5, 6, 7, 23, 8, 9, 10, 15, 16, 14, 17, 13, 11, 12, 22, 21, 20, 19, 18))

state_focal_indices <- list(coarse = list(chinamain = c(1), asiapacific = c(2:3, 12), europe = 5:8, us = 9:10), 
                            fine = list(chinamain = 1:4, asiapacific = c(5:8, 23), europe = 11:17, us = 18:21))


