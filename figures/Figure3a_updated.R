library(stringr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(fuzzyjoin)
library(readr)
library(tidyr)
library(rstatix)       
library(multcompView) 
library(purrr)
'%!in%'<- Negate('%in%')
setwd("/Users/mingyuwang/Desktop/Rebecca_revision_summary/Figure_edit/out/all_HOR_bed_csv")

all_array<- read.table("All_Arrays_Total")
colnames(all_array) <- c("line","rep","class","mon_count")
all_array$data<- "Chromosomes"

#Ab10
Ab10_bed <- read.csv("AB10_HOR_bed.csv", sep = ",", row.names = 1)
colnames(Ab10_bed ) <- c("bin","Pattern","Start","End")
Ab10_string <- read.csv("AB10_string_out.csv", sep = ",", row.names = 1)
colnames(Ab10_string) <- c("chr","bin","clust_val","start","rep","string","purity")
Ab10_df <- left_join(Ab10_bed, Ab10_string,by="bin")

Ab10_add_bed <- read.csv("AB10_ab10_HOR_bed.csv", sep = ",", row.names = 1)
colnames(Ab10_add_bed) <- c("bin","Pattern","Start","End")
Ab10_add_string <- read.csv("AB10_ab10_string_out.csv",sep = ",", row.names = 1)
colnames(Ab10_add_string) <- c("chr","bin","clust_val","start","rep","string","purity")
Ab10_add <- left_join(Ab10_add_bed, Ab10_add_string,by="bin")

Ab10_all <- bind_rows(Ab10_df, Ab10_add)
Ab10_all$num_chara <- nchar(Ab10_all$string)
Ab10_all$num_hormon <- Ab10_all$num_chara * Ab10_all$purity
Ab10_select_df <- Ab10_all %>% select(bin,rep, num_chara,num_hormon) %>% distinct()

Ab10_sum_df <- Ab10_select_df %>% group_by(rep) %>% summarise(total= sum(num_chara, na.rm = TRUE))
Ab10_sum_df_hor <- Ab10_select_df %>% group_by(rep) %>% summarise(total_hor= sum(num_hormon, na.rm = TRUE)) %>% mutate(line="Ab10")
Ab10_array <- all_array %>% filter(line=="AB10")
sum_Ab10_array <- Ab10_array %>% group_by(rep) %>% summarise(total_mon= sum(mon_count, na.rm = TRUE)) %>%  mutate(line="Ab10")
prop_hor_Ab10 <- merge(sum_Ab10_array,Ab10_sum_df_hor, by=c("line","rep")) %>% mutate(prop = total_hor / total_mon)


# Fix TIL25
TIL25_bed <- read.csv("TIL25_HOR_bed.csv",sep = ",", row.names = 1)
TIL25_string <- read.csv("TIL25_string_out.csv",sep = ",", row.names = 1)
colnames(TIL25_string) <- c("chr","bin","clust_val","start","rep","string","purity")

TIL25_string$num_chara <- nchar(TIL25_string$string)
TIL25_string$num_hormon <- TIL25_string$num_chara * TIL25_string$purity
TIL25_string$poss_end <- TIL25_string$start + TIL25_string$num_chara*156
TIL25_check <- TIL25_string %>% select(chr,start,poss_end)
#write.table(file="TIL25_check.bed",TIL25_check,sep = "\t", quote = FALSE, row.names = FALSE)

# 1) Read
TIL25_bins    <- read_tsv("check_bin.sorted.txt", col_names = "bin", show_col_types = FALSE)
TIL25_TR1     <- read_tsv("TR1_summar.txt",
                          col_names = c("bin","count_tr1","repeat_tr1"),
                          col_types = cols(bin=col_character(), count_tr1=col_double(), repeat_tr1=col_character()),
                          show_col_types = FALSE)
TIL25_CentC   <- read_tsv("CentC_summar.txt",
                          col_names = c("bin","count_centc","repeat_centc"),
                          col_types = cols(bin=col_character(), count_centc=col_double(), repeat_centc=col_character()),
                          show_col_types = FALSE)
TIL25_knob180 <- read_tsv("knob_summar.txt",
                          col_names = c("bin","count_knob","repeat_knob"),
                          col_types = cols(bin=col_character(), count_knob=col_double(), repeat_knob=col_character()),
                          show_col_types = FALSE)

# 2) Merge
TIL25_merged <- TIL25_bins %>%
  left_join(TIL25_TR1,     by = "bin") %>%
  left_join(TIL25_CentC,   by = "bin") %>%
  left_join(TIL25_knob180, by = "bin") %>%
  mutate(across(starts_with("count_"), ~ coalesce(., 0)))

# (optional but recommended) ensure genomic order before up-fill
TIL25_merged <- TIL25_merged %>%
  separate(bin, into = c("chr",".range"), sep = ":", remove = FALSE) %>%
  mutate(start = as.numeric(sub("-.*", "", .range))) %>%
  arrange(chr, start) %>%
  select(-.range)

# 3) Up-fill rows where ALL repeat_* are NA
rep_cols   <- c("repeat_tr1","repeat_centc","repeat_knob")
count_cols <- c("count_tr1","count_centc","count_knob")

filled <- TIL25_merged %>%
  mutate(all_rep_na = is.na(repeat_tr1) & is.na(repeat_centc) & is.na(repeat_knob)) %>%
  mutate(across(all_of(count_cols), ~ if_else(all_rep_na, NA_real_, .))) %>%
  fill(all_of(c(rep_cols, count_cols)), .direction = "up") %>%
  mutate(across(all_of(count_cols), ~ coalesce(., 0))) %>%
  select(-all_rep_na)

# 4) Winner per bin
top_per_bin <- filled %>%
  pivot_longer(
    cols = matches("^(count|repeat)_(tr1|centc|knob)$"),
    names_to   = c(".value","type"),
    names_pattern = "(count|repeat)_(tr1|centc|knob)"
  ) %>%
  group_by(bin) %>%
  slice_max(order_by = count, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(bin, top_repeat = type, top_count = count)

# 5) chr + chr_start bin; standardize names
cleaned <- top_per_bin %>%
  separate(bin, into = c("chr","range"), sep = ":", remove = FALSE) %>%
  mutate(
    start = as.numeric(sub("-.*","", range)),
    bin   = paste0(chr, "_", start),
    top_repeat = case_when(
      top_repeat == "centc" ~ "CentC",
      top_repeat == "tr1"   ~ "TR1",
      top_repeat == "knob"  ~ "knob180",
      TRUE ~ top_repeat
    )
  ) %>%
  select(chr, bin, top_repeat, top_count)

# 6) Override window on chr4 → Cent4
cleaned <- cleaned %>%
  mutate(
    start_num = as.numeric(sub(".*_", "", bin)),
    top_repeat = if_else(chr == "chr4" & between(start_num, 111314884, 111398518),
                         "Cent4", top_repeat)
  ) %>%
  select(-start_num) %>% mutate(rep=top_repeat) %>% select(chr, bin, rep)

correct_TIL25_string <- merge(TIL25_string, cleaned, by=c("chr","bin"))  %>% rename(rep = rep.y)

TIL25_sum_df <- correct_TIL25_string %>% select(bin, rep, num_chara, num_hormon) %>% group_by(rep) %>% summarise(total= sum(num_chara, na.rm = TRUE))
TIL25_sum_df_hor <- correct_TIL25_string %>% select(bin, rep, num_chara, num_hormon) %>% group_by(rep) %>% summarise(total_hor= sum(num_hormon, na.rm = TRUE)) %>% mutate(line="TIL25")
TIL25_array <- all_array %>% filter(line=="TIL25")
sum_TIl25_array <- TIL25_array %>% group_by(rep) %>% summarise(total_mon= sum(mon_count, na.rm = TRUE)) %>%  mutate(line="TIL25")
prop_hor_TIL25 <- merge(sum_TIl25_array,TIL25_sum_df_hor, by=c("line","rep")) %>% mutate(prop = total_hor / total_mon)

# non_Ab10, non_TIL25
# discover lines from *_HOR_bed.csv (e.g., "AB10", "B73", "CG108", ...)
lines <- list.files(pattern = "_HOR_bed\\.csv$") |>
  sub("_HOR_bed\\.csv$", "", x = _)
lines <- lines[3:14]
lines <- gsub("^CG119_fix$", "CG119", lines)
lines <- lines[lines != "TIL25"]

# --- helper: read bed/string pair ---
.read_pair <- function(base) {
  bed_file <- paste0(base, "_HOR_bed.csv")
  string_file <- paste0(base, "_string_out.csv")
  
  bed <- read.csv(bed_file, sep = ",", row.names = 1)
  colnames(bed) <- c("bin","Pattern","Start","End")
  
  string <- read.csv(string_file, sep = ",", row.names = 1)
  colnames(string) <- c("chr","bin","clust_val","start","rep","string","purity")
  
  left_join(bed, string, by = "bin")
}

# --- processor for each line ---
process_line <- function(line_name) {
  message("Processing ", line_name, " ...")
  
  if (line_name == "AB10") {
    df1 <- .read_pair("AB10")
    df2 <- .read_pair("AB10_ab10")
    df_all <- bind_rows(df1, df2)
  } else if (line_name == "CG119") {
    # explicitly read CG119_fix files
    df_all <- .read_pair("CG119_fix")
  } else {
    df_all <- .read_pair(line_name)
  }
  
  df_all <- df_all %>%
    mutate(num_chara = nchar(string),
           num_hormon = num_chara * purity)
  
  select_df <- df_all %>%
    select(bin, rep, num_chara, num_hormon) %>%
    distinct()
  
  sum_df_hor <- select_df %>%
    group_by(rep) %>%
    summarise(total_hor = sum(num_hormon, na.rm = TRUE), .groups = "drop") %>%
    mutate(line = line_name)
  
  sum_array <- all_array %>%
    filter(line == line_name) %>%
    group_by(rep) %>%
    summarise(total_mon = sum(mon_count, na.rm = TRUE), .groups = "drop") %>%
    mutate(line = line_name)
  
  full_join(sum_array, sum_df_hor, by = c("line","rep")) %>%
    mutate(prop = total_hor / total_mon)
}

# --- run for all lines ---
prop_hor_all <- map_dfr(unique(lines), process_line) %>% select(line, rep, total_mon, total_hor,prop)
prop_hor_all <- bind_rows(prop_hor_all, prop_hor_Ab10,prop_hor_TIL25)

prop_hor_all$line <- factor(prop_hor_all$line,
                            levels = c(
                              "B73", "Ab10", "Mo17", "CG108", "CG119", "CG44",
                              "Tx777", "Tx779", "CML442", "K64", "TIL01", "TIL11", "TIL25"
                            )
)

plot <- ggplot(prop_hor_all, aes(x = line, y = prop, fill = rep)) +
  geom_bar(stat = "identity", color = "black") +
  facet_wrap(~ rep, ncol = 1, strip.position = "right", scales = "free_y") +
  scale_fill_manual(values = c(
    "Cent4"   = "#88CCEE",
    "CentC"   = "#6699CC",
    "knob180" = "#882255",
    "TR1"     = "#CC6677"
  )) +
  ylim(0,1) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    legend.position = "none",
    strip.background = element_rect(fill = "gray90", color = "gray50"),
    strip.text = element_text(size = 12),
    panel.spacing.y = unit(0.6, "lines"),
    panel.grid = element_blank()
  ) +
  labs(
    x = "Line",
    y = "Proportion",
    title = "Proportion of HORs per repeat type across maize lines"
  )

plot

ggsave("HOR_13_lines_prop.pdf", plot, height = 10, width = 12, dpi = 300)
