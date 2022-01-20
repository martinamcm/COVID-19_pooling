############################################################
# Epi curve for Hong Kong COVID-19 - updated for 5th wave 
#
# Date: 2022-01-11
############################################################


# Helper functions --------------------------------------------------------

# Inset epi curve

gg_inset <- function (grob, 
                      xmin = -Inf, 
                      xmax = Inf, 
                      ymin = -Inf, 
                      ymax = Inf, 
                      data) {
  layer(data = data, 
        stat = StatIdentity, 
        position = PositionIdentity, 
        geom = ggplot2::GeomCustomAnn,
        inherit.aes = TRUE, 
        params = list(grob = grob, 
                      xmin = xmin, 
                      xmax = xmax, 
                      ymin = ymin, 
                      ymax = ymax)
        )
}




# Panel A - Epi curve with RCHE inset -------------------------------------


## Main epi curve


# read most recent publicly available CHP data 

case <- read.csv(
  "~/SPH Dropbox/HK_COVID/data_CHP_public/enhanced_surveillance/hkcase_20220111.csv"
  )


# case numbers related to RCHE - staff and residents

Tot_case <- c(1298, 1302, 1307, 1310, 1311, 1312, 1316, 1321, 1323, 1334, 1335, 
              1336, 1337, 1338, 1339, 1340, 1341, 1342, 1343, 1344, 1345, 1346,
              1347, 1348, 1349, 1350, 1351, 1352, 1353, 1354, 1355, 1360, 1395, 
              1407, 1408, 1410, 1473, 1485, 1560, 1563, 1606, 1659, 1688, 1691,
              1700, 1708, 1908, 2052, 2224, 2246, 2254, 2305, 2343, 2351, 2404, 
              2405, 2477, 2506, 2543, 2586, 2587, 2588, 2589, 2590, 2591, 2592, 
              2604,2695, 2698, 2701, 2702, 2703, 2704, 2705, 2708, 2745, 2746, 
              2747, 2748, 2749, 2750, 2756, 2793, 2794, 2795, 2820, 2842, 2882, 
              2893, 2938, 2940, 2941, 2942, 2955, 2970, 2975, 2991, 3068, 3090,
              3153, 3215, 3260, 3501, 3502, 3553, 3569, 3584, 3585, 3637, 3683, 
              3707, 3733, 3734, 3735, 3739, 3771, 3781, 3791, 3802, 3860, 3903, 
              3962, 3963, 3990, 4065, 4072, 4085, 4179, 4367, 4565, 4630, 4654, 
              5083, 5928, 6224, 6275, 6342, 6370, 6381, 6524, 6606, 6691, 6724, 
              6745, 6774, 6795, 6798, 6800,6812, 6831, 6839, 6846, 6849, 6850, 
              6923, 7053, 7127, 7128, 7129, 7130, 7131, 7132, 7133, 7320, 7584, 
              7673, 8603, 8720, 8741, 8927, 9911, 9951, 10100, 10258, 10302, 
              10503, 10525, 10535) 


# map RCHE case numbers to case data
case_Eld <- case[which(case$case.no %in% Tot_case),]

# clean dates
case_Eld$onset.date <- as.Date(case_Eld$onset.date,"%d/%m/%Y")
case_Eld$confirm.date <- as.Date(case_Eld$confirm.date,"%d/%m/%Y")


## Define by location 

##Kong Tai Elderly Home outbreak - by resident, staff
KongTaiCase <- c(1298,1302,1307,1310,1311,1312,1316,1321,1323,1334,1335,1336,
                 1337,1338,1339,1340,1341,1342,1343,1344,1345,1346,1347,1348,
                 1349,1350,1351,1352,1353,1354,1355,1360,1395,1407,1408,1410,
                 1473,1485,1560,1606,1659,1688,2343,3860)

KongTaiLink <- c(rep("Resident",3),
                 rep("Staff",3),
                 "Resident",
                 "Staff",
                 rep("Resident",25),
                 rep("Staff",3),
                 rep("Resident",4),
                 "Staff",
                 "Resident",
                 rep("Staff",2)
                 )


##Harmony Villa Elderly Home 
HVCase <- c(1691,1700,1708,1908,3781,3990)
HVLink <- c("Staff",
            "Resident",
            rep("Staff",2),
            rep("Resident",2)
            )


##Cornwall Elderly Home 
CornwallCase <- c(2254,2477,2543,2586,2587, 2588, 2589, 2590, 2591, 2592, 2701, 
                  2702, 2703, 2704, 2705, 2708, 2745, 2746, 2747, 2748, 2749, 
                  2750, 2756, 2793, 2794, 2820, 2842, 2938, 2940, 2941,
                  2942, 2955, 2975, 3090, 3153, 3215, 3260, 3501,3683,3802)
CornwallLink <- c("Staff",
                  "Resident",
                  rep("Staff",2),
                  rep("Resident",23),
                  "Staff",
                  rep("Resident",4),
                  "Staff",
                  rep("Resident",7)
                  )


##Lung Hang Elderly Home 
LungHangCase <- c(2052, 2224, 2305, 2404, 2405, 2506, 2795, 2882,3569,4367)
LungHangLink <- c(rep("Staff",2),
                  "Resident",
                  rep("Resident",3),
                  rep("Staff",2),
                  rep("Resident",2)
                  )

## King Fok Limited 
KingFokCase <- c(3068,3502,3584,3585,3637,3707,3733,3734,3735,3739,3771,3791,
                 3962,3963)
KingFokLink <- c(rep("Resident",9),
                 "Staff",
                 rep("Resident",2),
                 rep("Staff",2)
                 )

## SAGE Kai Yip
KaiYipCase <- c(3903,4065,4072,4085,4179,4565)
KaiYipLink <- c(rep("Resident",6))


## Tai Po - Hong Wo Home for the Elderly

HongWoCase <- c(2351,2695)
HongWoLink <- c(rep("Resident",2))


## Po Chung Cheun Ying Home for the Elderly

PoChungCase <- c(2698,2991)
PoChungLink <- c("Resident",
                 "Staff"
                 )


## HoYukChing 

HoYukChingCase <- c(5928, 6224, 6275,  6370, 6381, 6691, 6724, 6745, 6774, 6795,
                    6798, 6800,6812, 6831, 6839, 6846, 6850, 6923, 7053, 7127, 
                    7128, 7129,7131, 7132, 7133, 7320) 
HoYukChingLink <- c(rep("Staff",4),
                    rep("Resident",4),
                    rep("Staff",2),
                    rep("Resident",2),
                    "Staff",
                    rep("Resident",12), 
                    "Staff"
                    )

##  Pak Lok Nursing Home

PakLokCase <- c(6342, 6524, 6849, 7130)
PakLokLink <- c(rep("Resident",4))

## San Po Kong

SanPoKongCase <- c(7584, 7673)
SanPoKongLink <- c(rep("Staff",2))


## Kowloon Bay Health Centre

KowloonBayCase <- c(8603, 8720, 8741, 8927)
KowloonBayLink <- c("Staff", 
                    rep("Resident",2), 
                    "Staff"
                    )


## Sheung On Elderly Home YMT

SheungOnCase <- c(9911, 9951)
SheungOnLink <- c(rep("Resident",2))

## Kin Ling Elderly Home 

KinLingCase <- c(10100, 10258, 10302)
KinLingLink <- c(rep("Resident",3))


## Salvation Army Nam Shan

NamShanCase <- c(10503, 10525, 10535)
NamShanLink <- c("Resident", 
                 "Staff", 
                 "Resident"
                 )



# create data frame

Eld_Link <- data.frame(case.no = c(KongTaiCase, HVCase, CornwallCase, LungHangCase, 
                                   KingFokCase, KaiYipCase, HongWoCase, PoChungCase, 
                                   HoYukChingCase, PakLokCase, SanPoKongCase, 
                                   KowloonBayCase, SheungOnCase, KinLingCase, 
                                   NamShanCase
                                   ), 
                       Homes = c(rep("KongTai",44), rep("HV",6), rep("Cornwall",40),
                                 rep("LH",10), rep("KingFok",14), rep("KaiYip",6),
                                 rep("HongWo",2), rep("PoChung",2), rep("HoYukChing",26),
                                 rep("PakLok",4), rep("SanPoKong",2), rep("KowloonBay",4),
                                 rep("SheungOn",2), rep("KinLing",3), rep("NamShan",3)
                                 ),
                       Links = c(KongTaiLink, HVLink, CornwallLink, LungHangLink,
                                 KingFokLink, KaiYipLink, HongWoLink, PoChungLink, 
                                 HoYukChingLink, PakLokLink, SanPoKongLink, 
                                 KowloonBayLink, SheungOnLink, KinLingLink,NamShanLink
                                 )
                       )

# merge 
case_merge <- merge(case_Eld, Eld_Link, by = "case.no")


# Plot epi curve 

## Inset RCHE plots 

### RCHE epi in third wave
inc_conf <- incidence(case_merge$confirm.date, 
                      interval = 1, 
                      groups = case_merge$Homes
                      )

case_merge_1 <- case_merge[case_merge$confirm.date<"2020/08/23",]
inc_conf_1 <- incidence(case_merge_1$confirm.date, 
                        interval = 1, 
                        groups = case_merge_1$Homes
                        )

colorpanelsEP1 <- c("Cornwall" = "#440154FF",
                    "HongWo" = "#443A83FF", 
                    "HV" = "#3B528BFF", 
                    "KaiYip" = "#31688EFF", 
                    "KingFok" = "#21908CFF", 
                    "KongTai" = "#27AD81FF", 
                    "LH" = "#5DC863FF", 
                    "PoChung"="#C7E020FF"
                    )

ggRCHEepi_1 <- plot(inc_conf_1,
                    alpha = 0.9,
                    stack = TRUE,
                    border = FALSE) + 
  theme_bw() +
  scale_fill_manual(values = colorpanelsEP1
                    ) +
  ggtitle("") + 
  ylim(c(0,25)) + 
  xlab("") + 
  ylab("") + 
  removeGridX() + 
  removeGridY()+
  scale_x_date(date_breaks = "months", 
               date_labels = "%d-%m-%y"
               ) + 
  theme(legend.position = "none",
        axis.text = element_text(size = 7)
        )



### RCHE epi in fourth wave

case_merge_2 <- case_merge[case_merge$confirm.date>"2020/10/02" & 
                             case_merge$confirm.date<"2021/01/07" ,]
inc_conf_2 <- incidence(case_merge_2$confirm.date, 
                        interval = 1, 
                        groups = case_merge_2$Homes
                        )


colorpanelsEP2 <-  c("HoYukChing" = "#481F70FF", 
                     "KinLing" = "#287C8EFF",  
                     "KowloonBay" = "35B779FF", 
                     "NamShan" = "#8FD744FF", 
                     "PakLok" = "#AADC32FF",
                     "SanPoKong" = "#E3E418FF", 
                     "SheungOn" = "#FDE725FF"
                     )

ggRCHEepi_2 <- plot(inc_conf_2, 
                    alpha = 0.9,
                    stack = TRUE,
                    border = FALSE
                    ) + 
  theme_bw() +
  ggtitle("") + 
  ylim(c(0,25)) + 
  xlab("") + 
  ylab("") + 
  removeGridX() + 
  removeGridY() + 
  scale_fill_manual(values = colorpanelsEP2 
                    ) +
  scale_x_date(date_breaks = "months", 
               date_labels = "%d-%m-%y"
               ) + 
  theme(legend.position = "none",
        axis.text = element_text(size=7)
        )


## RCHE > 07/01/2021  

case_merge_3 <- case_merge[case_merge$confirm.date>"2021/01/07" ,]
inc_conf_3 <- incidence(case_merge_3$confirm.date, 
                        interval = 1, 
                        groups = case_merge_3$Homes
                        )


colorpanelsEP2 <-  c("HoYukChing" = "#481F70FF", 
                     "KinLing" = "#287C8EFF",  
                     "KowloonBay" = "35B779FF", 
                     "NamShan" = "#8FD744FF", 
                     "PakLok" = "#AADC32FF",
                     "SanPoKong" = "#E3E418FF", 
                     "SheungOn" = "#FDE725FF"
                     )

ggRCHEepi_3 <- plot(inc_conf_3,
                    alpha = 0.9,
                    stack = TRUE,
                    border = FALSE
                    ) + 
  theme_bw() +
  ggtitle("") + 
  ylim(c(0,25)) + 
  xlab("") + 
  ylab("") + 
  removeGridX() + 
  removeGridY() + 
  scale_fill_manual(values = colorpanelsEP2 
                    ) +
  scale_x_date(date_breaks = "months", 
               date_labels = "%d-%m-%y"
               ) +
  theme(legend.position = "none",
        axis.text = element_text(size=7)
        )


## All HK Epi Curve 

case$confirm.date <- as.Date(case$confirm.date,"%d/%m/%Y")
case$RCHEcase <- ifelse(case$case.no %in% Tot_case, 1, 0)
case$RCHEcase <- factor(case$RCHEcase,
                        levels = c(0,1),
                        labels = c("Other","RCHE")
                        )

# make in to an incidence object
inc_all <- incidence(case$confirm.date, 
                     interval = "1 week",
                     groups = case$RCHEcase
                     )

# plot
ggepiHK <- plot(inc_all,
                stack = TRUE,
                border = "white"
                ) + 
  theme_bw() +
  ggtitle("") + 
  ylab("Number of confirmed cases") + 
  scale_fill_manual(values = c("#D55E00","#999999"),
                    breaks = c("RCHE","Other")
                    ) +
  xlab("Date of confirmation") + 
  scale_x_date(date_breaks = "months", date_labels = "%m-%y") +
  theme(legend.position = c(0.94,0.85),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(4,"mm"),
        axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   size = 8),
        axis.title = element_text(size=9)
        ) + 
  removeGridX() + 
  removeGridY() + 
  labs("RCHE") +
  annotation_custom(ggplotGrob(ggRCHEepi_1), 
                    ymin = 400, 
                    ymax = 1010, 
                    xmin=18230, 
                    xmax=18455
                    ) + 
  labs(fill="Cases") +
  annotation_custom(ggplotGrob(ggRCHEepi_2), 
                    ymin = 400, 
                    ymax = 1010,  
                    xmin=18652, 
                    xmax=18840
                    ) + 
  labs(fill="Cases") +
  annotation_custom(ggplotGrob(ggRCHEepi_3), 
                    ymin = 400, 
                    ymax = 1010, 
                    xmin=18840, 
                    xmax=18955
                    ) + 
  labs(fill="Cases")






# Panel B - onset of symptoms in symptomatic cases ------------------------


# define difference in onset and confirmation date
case_merge$date_diff <- as.Date(case_merge$confirm.date, format="%Y/%m/%d")-
  as.Date(case_merge$onset.date, format="%Y/%m/%d")


#### Time series from onset to confirmation

# second wave
case_merge_seg <- subset(case_merge,is.na(onset.date) == FALSE)
case_merge_seg <- case_merge_seg[case_merge_seg$onset.date<"2020-08-21",]
case_merge_seg <- case_merge_seg[order(case_merge_seg$Homes),]
case_merge_seg$ID <- seq(1,87,1)


ggons <- ggplot(case_merge_seg, 
                aes(x = ID,
                    ymin = onset.date,
                    ymax = confirm.date,
                    color = Homes,
                    alpha = Links)
                ) + 
  geom_linerange(size=2, 
                 position = position_dodge(0.3)
                 ) +
  theme_bw() + 
  scale_color_viridis(discrete = T, 
                      option = "D",
                      alpha = 1,
                      begin = 0,  
                      labels=c("Cornwall",
                               "Harmony Villa",
                               "Kai Yip",
                               "King Fok",
                               "Kong Tai",
                               "Lung Hang") 
                      ) +
  ggtitle("") + 
  ylab("Onset-to-Confirmation") + 
  xlab("") + 
  scale_y_date(date_breaks = "days", 
               date_labels = "%d-%b") + 
  scale_alpha_discrete(range=c(1,0.5)
                       ) +
  theme(legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        axis.text = element_text(size=8),
        axis.title = element_text(size=9),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   size = 8)
        ) + 
  labs(fill="Link") +
  coord_flip() + 
  labs(color = "RCHE") + 
  removeGridY()

case_merge_seg$Facetvals <- factor(case_merge_seg$Homes,
                                   levels = c('KaiYip',
                                              'KingFok',
                                              'PoChung',
                                              'HongWo',
                                              'Cornwall',
                                              'LH',
                                              'HV',
                                              'KongTai')
                                   )

testaxis <- seq(as.Date("2020-07-04"),
                as.Date("2020-08-23"),
                1)

labeldates <- c("04-07", "", "06-07","", "08-07","", "10-07","", "12-07","",
                "14-07","", "16-07","","18-07", "", "20-07", "", "22-07", "", 
                "24-07","", "26-07","", "28-07","", "30-07","", "01-08", "", 
                "03-08", "", "05-08","", "07-08", "", "09-08","", "11-08","",
                "13-08", "","15-08","", "17-08", "", "19-08", "", "21-08", "",
                "23-08")


colorpanels <- c("Cornwall" = "#440154FF", 
                 "HoYukChing" = "#481F70FF",
                 "HongWo" = "#443A83FF",
                 "HV" = "#3B528BFF", 
                 "KaiYip" = "#31688EFF", 
                 "KinLing" = "#287C8EFF",  
                 "KingFok" = "#21908CFF", 
                 "KongTai" = "#27AD81FF", 
                 "KowloonBay" = "35B779FF",
                 "LH" = "#5DC863FF", 
                 "NamShan" = "#8FD744FF", 
                 "PakLok" = "#AADC32FF", 
                 "PoChung" = "#C7E020FF", 
                 "SanPoKong" = "#E3E418FF", 
                 "SheungOn" = "#FDE725FF")

ggons.facet <- ggplot(case_merge_seg,
                      aes(x = ID,
                          ymin = onset.date,
                          ymax = confirm.date,
                          color = Homes)
                      ) + 
  geom_linerange(size = 2, 
                 position = position_dodge(1)
                 ) +
  facet_grid(Facetvals~.,
             drop = TRUE,
             scale = "free_y",
             space = "free") + 
  theme_bw() +
  ggtitle("") + 
  ylab("Onset-to-Confirmation") + 
  xlab("") + 
  scale_y_date(breaks = testaxis, 
               labels=labeldates
               ) +
  theme(legend.position = "none",
        panel.spacing = unit(0.1, "lines"),  
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 11),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y = element_blank(),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   size = 8)
        ) + 
  labs(fill = "Link") +
  coord_flip() + 
  labs(color = "RCHE") + 
  removeGridY() + 
  scale_color_manual(values = colorpanels
                     ) + 
  scale_alpha_discrete(range = c(1,0.5))


#  third wave 

case_merge_seg <- subset(case_merge,is.na(onset.date)==FALSE)
case_merge_seg1 <- case_merge_seg[case_merge_seg$onset.date>"2020-08-21",]
case_merge_seg1 <- case_merge_seg1[order(case_merge_seg1$Homes),]
case_merge_seg1$ID <- seq(1,24,1)

case_merge_seg1$Facetvals <- factor(case_merge_seg1$Homes,levels = 
                                      c('HoYukChing',
                                        'PakLok', 
                                        'KinLing', 
                                        'KowloonBay',
                                        'NamShan',
                                        'SanPoKong',
                                        'SheungOn')
                                    )

colorpanels_1 <- c("HoYukChing" = "#481F70FF", 
                   "PakLok" = "#AADC32FF", 
                   "KinLing" = "#287C8EFF", 
                   "KowloonBay" = "35B779FF", 
                   "NamShan" = "#8FD744FF", 
                   "SanPoKong" = "#E3E418FF", 
                   "SheungOn" = "#FDE725FF")

testaxis_1 <- seq(as.Date("2020-11-12"),
                  as.Date("2021-02-02"),
                  1)

labeldates_1 <- c("12-11","", "14-11", "","16-11","", "18-11", "", "20-11", "", 
                  "22-11", "", "24-11","", "26-11","", "28-11", "","30-11","", 
                  "02-12", "", "04-12", "", "06-12","","08-12", "", "10-12", 
                  "","12-12", "", "14-12", "", "16-12", "", "18-12", "", "20-12",
                  "", "22-12", "", "24-12", "", "26-12", "", "28-12", "", "30-12", 
                  "", "01-01", "", "03-01", "", "05-01", "", "07-01", "", "09-01",
                  "11-01", "", "13-01", "", "15-01", "", "17-01", "", "19-01", "", 
                  "21-01", "", "23-01", "", "25-01", "", "27-01", "", "29-01", 
                  "", "31-01", "", "02-02", "")

ggons.facet2 <- ggplot(case_merge_seg1, 
                       aes(x = ID,
                           ymin = onset.date,
                           ymax = confirm.date,
                           color = Homes)
                       ) + 
  geom_linerange(size = 2, 
                 position = position_dodge(1)
                 ) +
  facet_grid(Facetvals~.,
             drop = TRUE,
             scale = "free_y",
             space = "free") + 
  theme_bw() + 
  scale_color_manual(values = colorpanels_1
                     ) +
  ggtitle("") + 
  ylab("Onset-to-Confirmation") + 
  xlab("") + 
  scale_y_date(breaks = testaxis_1,
               labels = labeldates_1
               ) + 
  scale_alpha_discrete(range = c(1,0.5)
                       ) +
  theme(legend.position = "none",
        panel.spacing = unit(0.1,"lines"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 11),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y = element_blank(),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   size = 8)
        ) + 
  labs(fill = "Link") +
  coord_flip() + 
  labs(color = "RCHE") + 
  removeGridY()




# Panel C - Vaccination in HK ---------------------------------------------


# read most recent vax data
vax_data <- read_csv(
  "~/Documents/GitHub/COVID-19_pooling/data/vax_time_HK.csv"
  )

# make in to long format

vax_clean <- 
  gather(vax_data, dose, vax, c(`Sinovac 1st dose`, 
                                `Sinovac 2nd dose`,
                                `Sinovac 3rd dose`,
                                `BioNTech 1st dose`,
                                `BioNTech 2nd dose`,
                                `BioNTech 3rd dose`),
         factor_key = TRUE
         )  %>% 
  mutate(date_vax = as.Date(Date, "%d/%m/%Y"))


# combine M and F to give totals each day
  
n_pop = 962300 # population of 70s and over in HK

ggvax <-
  vax_clean %>%
  filter(`Age Group` == "80 and above" | `Age Group` == "70-79",
         dose == "BioNTech 2nd dose" | dose == "Sinovac 2nd dose"
         ) %>%
  group_by(dose, date_vax) %>%
  summarise_if(is.numeric, funs(sum)) %>%
  padr::pad(group = "dose") %>%
  mutate_at(
    vars(vax),
    ~replace(., 
             is.na(.),
             0)
    ) %>%
  mutate(id = seq(1, length(dose))
         ) %>%
  arrange(dose, date_vax) %>%
  group_by(dose, date_vax) %>%
  mutate(
    av_vax = slider::slide_index_dbl(.x = vax, 
                             .i = id, 
                             .f = mean,
                             .before = 6)
         ) %>% 
  ggplot() +
  geom_col(aes(x = date_vax, 
               y = (vax/n_pop)*100,
               group = dose,
               fill = dose),
           position = position_dodge(3),
           alpha = 0.4) +
  geom_line(aes(x = date_vax, 
                y = (av_vax/n_pop)*100,
                group = dose, 
                colour = dose)
  ) +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  theme_minimal() +
  ylab("Second doses (%)") +
  labs(fill = "Vaccine type",
       colour = "Vaccine type") + 
  theme(
    legend.position = "bottom",
    text = element_text(size = 8),
    legend.key.size = unit(4, 'mm'), 
    #legend.key.height = unit(4, 'mm'), 
    #legend.key.width = unit(4, 'mm'), 
    legend.text = element_text(size = 9)
    )




# Number of beds  ---------------------------------------------------------


RCHE <- c("Cornwall", "Ho Yuk Ching", "Hong Wo","Harmony Villa","Kai Yip", 
          "Kin Ling", "King Fok", "Kong Tai", "Kowloon Bay", "Lung Hang", 
          "Nam Shan", "Pak Lok", "Po Chung", "San Po Kong", "Sheung On"
          )
RCHE <- factor(RCHE, levels = c("Cornwall", "Ho Yuk Ching", "Hong Wo",
                                "Harmony Villa","Kai Yip", "Kin Ling",
                                "King Fok", "Kong Tai", "Kowloon Bay", 
                                "Lung Hang", "Nam Shan", "Pak Lok", 
                                "Po Chung", "San Po Kong", "Sheung On")
               )

beds <- c(rep(20,length(RCHE)))

data.beds <- data.frame(RCHE,beds)

colorpanels_2 <- c("Cornwall" = "#440154FF", 
                   "Ho Yuk Ching" = "#481F70FF",
                   "Hong Wo" = "#443A83FF", 
                   "Harmony Villa" = "#3B528BFF", 
                   "Kai Yip" = "#31688EFF", 
                   "Kin Ling" = "#287C8EFF",  
                   "King Fok" = "#21908CFF", 
                   "Kong Tai" = "#27AD81FF", 
                   "Kowloon Bay" = "35B779FF",
                   "Lung Hang" = "#5DC863FF", 
                   "Nam Shan" = "#8FD744FF", 
                   "Pak Lok" = "#AADC32FF", 
                   "Po Chung" = "#C7E020FF", 
                   "San Po Kong" = "#E3E418FF", 
                   "Sheung On" = "#FDE725FF"
                   )

ggbeds <- 
  ggplot(
    data.beds
    ) + 
  geom_bar(
    aes(
      x = RCHE,
      y = beds,
      fill = RCHE
      ),
    width = 0.4,
    stat = "identity"
    ) + 
  scale_fill_manual(
    values = colorpanels_2
    ) + 
  theme_classic() + 
  theme(legend.position = "top",
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5,"line")
        ) + 
  xlab("") + 
  ylim(c(0,200)) + 
  ylab("Number of beds") + 
  guides(fill = 
           guide_legend(nrow = 4,
                        byrow = TRUE)
         )

leg <- get_legend(ggbeds)
legplot <- as_ggplot(leg)


# Combine panels on plot --------------------------------------------------



ggarrange(ggepiHK,
          ggarrange(ggons.facet,
                    ggarrange(ggons.facet2,
                              legplot,
                              NULL,
                              ggvax,
                              heights = c(1.2, 0.5, 0.1, 1 ),
                              nrow = 4,
                              labels=c("", "", "C")
                              ), 
                    widths = c(1, 1.4),
                    labels = c("B", "")
                    ),
          nrow=2,
          heights = c(1,2),
          labels = "A"
          ) 



  ggsave(
    "Figure1.pdf",
    plot = last_plot(),
    width = 20,
    height = 26,
    units = "cm"
  )


