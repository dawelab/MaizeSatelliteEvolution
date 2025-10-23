library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)

setwd("/Users/mingyuwang/Desktop/Figure_edit/output_QC_assembly")
chr7 <- read.table("CG108.CCS.Chr7.F2308.bed")
chr8 <-  read.table("CG108.CCS.Chr8.F2308.bed")
agp <- read.table("CG108_out_HALF_ragtag.scaffold.agp")
agp <- agp %>% mutate(across(c( V2, V3), as.numeric)) 
cols <-  c("chr", "start", "end", "first", "second")
colnames(chr7) <- cols
colnames(chr8) <- cols

chr7$pos <- (chr7$end + chr7$start)/2
chr8$pos <- (chr8$end + chr8$start)/2

# reshape into long format
chr7_long <- chr7 %>%pivot_longer(cols = c(first, second),names_to = "type",values_to = "value")
chr8_long <- chr8 %>%pivot_longer(cols = c(first, second),names_to = "type",values_to = "value")

# plot chr7, knob 153601741-179841948, CentC 53696541-58454692
agp7 <- agp%>% filter(V1=="chr7_RagTag")%>% filter(V5=="W")
cov7 <- ggplot(chr7_long, aes(x = pos, y = value, color = type)) +
  geom_rect(aes(xmin=153601741,xmax=179841948,ymin=0,ymax=100),fill="lightgrey",color="lightgrey", alpha=0.9)+
  geom_rect(aes(xmin=53696541,xmax=58454692,ymin=0,ymax=100),fill="lightyellow",color="lightyellow", alpha=0.9)+
  geom_rect(data = agp7,
            aes(xmin = as.numeric(V2), xmax = as.numeric(V3), ymin = 95, ymax = 100),
            fill = "black", color = "white", inherit.aes = FALSE) +
  geom_point(size = 0.5, alpha = 0.7) +
  ylim(0,100) +
  labs(x = "Chromosome7 coordinates (bp)",
       y = "Read Depth",
       color = "Alignment") +
  scale_color_manual(
    values = c("first" = "blue", "second" = "red"),
    labels = c("first" = "Primary", "second" = "Secondary")
  ) +
  theme_bw(base_size = 20)

cov7
ggsave("QC_chr7.pdf",cov7,height = 6,width = 15,dpi=300)

# plot chr8, CentC 49722591-52115502, knob 162843190-191205897
agp8 <- agp%>% filter(V1=="chr8_RagTag") %>% filter(V5=="W")
cov8 <- ggplot(chr8_long, aes(x = pos, y = value, color = type)) +
  geom_rect(aes(xmin=162843190,xmax=191205897,ymin=0,ymax=100),fill="lightgrey",color="lightgrey", alpha=0.9)+
  geom_rect(aes(xmin=49722591,xmax=52115502,ymin=0,ymax=100),fill="lightyellow",color="lightyellow", alpha=0.9)+
  geom_rect(data = agp8,
            aes(xmin = as.numeric(V2), xmax = as.numeric(V3), ymin = 95, ymax = 100),
            fill = "black", color = "white", inherit.aes = FALSE) +
  geom_point(size = 0.5, alpha = 0.7) +
  ylim(0,100) +
  labs(x = "Chromosome8 coordinates (bp)",
       y = "Read Depth",
       color = "Alignment") +
  scale_color_manual(
    values = c("first" = "blue", "second" = "red"),
    labels = c("first" = "Primary", "second" = "Secondary")
  ) +
  theme_bw(base_size = 20)

cov8
ggsave("QC_chr8.pdf",cov8,height = 6,width = 15,dpi=300)
