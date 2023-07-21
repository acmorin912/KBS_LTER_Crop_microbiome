setwd("~/Desktop/Microbiome/MSA23/fungaltraits")

data<- read.delim("Overall_fungaltraits.txt", header= TRUE)

data_leaf<- filter(data, Origin == "leaf")
  data_leaf$Primary_guild_unique <- reorder(data_leaf$Primary_guild_unique, -data_leaf$Primary_guild_percent)
  data_leaf <- filter(data_leaf, Primary_guild_percent >= 1)
data_root<- filter(data, Origin == "root")
    data_root$Primary_guild_unique <- reorder(data_root$Primary_guild_unique, -data_root$Primary_guild_percent)
    data_root <- filter(data_root, Primary_guild_percent >= 1)
data_soil<- filter(data, Origin == "soil")
  data_soil$Primary_guild_unique <- reorder(data_soil$Primary_guild_unique, -data_soil$Primary_guild_percent)
  data_soil <- filter(data_soil, Primary_guild_percent >= 1)
data_stem<- filter(data, Origin == "stem")
    data_stem$Primary_guild_unique <- reorder(data_stem$Primary_guild_unique, -data_stem$Primary_guild_percent)
    data_stem <- filter(data_stem, Primary_guild_percent >= 1)
library(ggplot2)

# Create the bar plot
leaf_bar<- ggplot(data_leaf, aes(x = Primary_guild_unique, y = Primary_guild_percent, fill = Management)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(values = c("Organic" = "#06a651", "Conventional" = "#feffbe", "No-Till" = "#feff40")) +
  theme_minimal() +
  xlab("Primary Guild") +
  ylab("Percentage") +
  ggtitle("Leaves")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_continuous(breaks = seq(0, 40, by = 2))+
  guides(fill = FALSE)
leaf_bar

root_bar<- ggplot(data_root, aes(x = Primary_guild_unique, y = Primary_guild_percent, fill = Management)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(values = c("Organic" = "#06a651", "Conventional" = "#feffbe", "No-Till" = "#feff40")) +
  theme_minimal() +
  xlab("Primary Guild") +
  ylab("Percentage") +
  ggtitle("Roots")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_continuous(breaks = seq(0, 50, by = 2))+
  guides(fill = FALSE)
root_bar

soil_bar<- ggplot(data_soil, aes(x = Primary_guild_unique, y = Primary_guild_percent, fill = Management)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(values = c("Organic" = "#06a651", "Conventional" = "#feffbe", "No-Till" = "#feff40")) +
  theme_minimal() +
  xlab("Primary Guild") +
  ylab("Percentage") +
  ggtitle("Soil")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_continuous(breaks = seq(0, 55, by = 2))+
  guides(fill = FALSE)
soil_bar

stem_bar<- ggplot(data_stem, aes(x = Primary_guild_unique, y = Primary_guild_percent, fill = Management)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(values = c("Organic" = "#06a651", "Conventional" = "#feffbe", "No-Till" = "#feff40")) +
  theme_minimal() +
  xlab("Primary Guild") +
  ylab("Percentage") +
  ggtitle("Stems")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_continuous(breaks = seq(0, 40, by = 2))
stem_bar

library(ggpubr)
ggarrange(leaf_bar,root_bar,soil_bar,stem_bar,
          labels = c("A", "B","C","D"),
          align = c("h"),
          ncol = 4, nrow = 1)
