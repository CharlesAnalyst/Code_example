setwd("/home/galaxy/project/m6AQTL/data/data_for_plot/")
library(ggplot2)
library(grid)
theme_set(
  theme_bw() + 
    theme(legend.position = "right")
)

# Load data
df <- read.table("data_for_bubble.txt", sep = "\t", header = TRUE)
head(df)

ggplot(df, aes(x = log2_Number, y = QTL_Type_pos)) + 
  geom_point(aes(colour = QTL_Type, size=log2_Number, shape=Feature_Type), alpha = 0.5) + 
  theme(legend.title = element_blank()) + theme(legend.text=element_text(size=18)) + guides(shape=guide_legend(override.aes = list(size=4))) + guides(colour=guide_legend(override.aes = list(size=4))) + 
  scale_colour_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#FF69B4", "#7B68EE")) + 
  scale_size(range = c(2, 15), guide = FALSE) +  # Adjust the range of points size
  theme(axis.title = element_blank(), axis.text.x = element_text(face = "bold"), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.border = element_blank()) + # , axis.line = element_line() 
  annotate("segment", x=0.0, xend = 0.0, y=0, yend = 6.5, linetype="dotted") +
  annotate("segment", x=5.0, xend = 5.0, y=0, yend = 6.5, linetype="dotted") +
  annotate("segment", x=10.0, xend = 10.0, y=0, yend = 6.5, linetype="dotted") + 
  annotate("segment", x=15.0, xend = 15.0, y=0, yend = 6.5, linetype="dotted") + 
  
  # annotate("segment", x=0.01, xend = 5.8, y=6, yend = 6, alpha=0.3) +
  annotate("text", label = expression(paste("Pseudouridine (", Psi, ")")),  x = 8, y = 5.5, size=6, face="bold") +
  # geom_text(x=-0.5, y=6, label = "Ψ", size=8) # +scale_x_continuous(expand = c(.1, .1)) # , vjust="inward",hjust="inward"
  annotate("text", label = expression("N"^6*"-Methyladenosine"~"(m"^6*"A)"), x = 7, y = 3.6, size=6, face="bold") +
  annotate("text", label = expression("5-methylcytosine"~"(m"^5*"C)"), x = 8, y = 1.5, size=6, face="bold") + 
  # ggtitle("QTL-associated RNA modifications") + theme(plot.title = element_text(size=20, face="bold", hjust = 0.6))
  annotate("segment", x=0.0, xend = 15.5, y=0.25, yend =0.25, arrow=arrow(), linetype="dotted") + 
  annotate("text", label = expression("log"[2]~"count"), x = 16.5, y = 0.25, size=6, face="bold")

           