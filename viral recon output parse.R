#Check files and Load packages---- 
list.files()
pacman::p_load(
  rio,        # importing data  
  here,       # relative file pathways  
  janitor,    # data cleaning and tables
  lubridate,  # working with dates
  matchmaker, # dictionary-based cleaning
  epikit,     # age_categories() function
  tidyverse   # data management and visualization
)

#loading data for analysis and data exploration ----

NPHL_data <- read.csv("Data_for_R.csv")
dim(NPHL_data)
str(NPHL_data)
summary(NPHL_data)
names(NPHL_data)

## Data cleaning ----
summary(NPHL_data$Age.years.)
linelist <- NPHL_data %>% 
  janitor::clean_names()
names(linelist)
# linelist data set is piped through select() command, and names() prints just the column names
linelist_vars <- linelist %>% 
  select(seq_name, age, age_years, sex, collection_week, year, region, nextclade_pango, clade, clade_display)

linelist_vars <- linelist_vars %>% 
  mutate(
    # Create categories
    age_group = dplyr::case_when(
      age_years <= 5            ~ "0-5",
      age_years > 6 & age_years <= 17 ~ "6-17",
      age_years > 17 & age_years <= 45 ~ "17-45",
      age_years > 45             ~ "> 45"
    ),
    # Convert to factor
    age_group = factor(
      age_group,
      level = c("0-5", "6-17","17-45", "> 45")
    )
  )
## changing variables to factor
library(aweek)
library(forcats)
library(janitor)
library(lubridate)
library(tidyr)
linelist_vars <- linelist_vars %>%
  mutate(age_group = fct_relevel(age_group))
linelist_vars <- linelist_vars %>%
  mutate(sex = fct_relevel(sex))
linelist_vars <- linelist_vars %>%
  mutate(region = fct_relevel(region))
linelist_vars <- linelist_vars %>%
  mutate(nextclade_pango = fct_relevel(nextclade_pango))
linelist_vars <- linelist_vars %>%
  mutate(clade = fct_relevel(clade))
linelist_vars <- linelist_vars %>%
  mutate(clade_display = fct_relevel(clade_display))


#plotting
ggplot(data = linelist_vars)+
  geom_bar(mapping = aes(x = region, fill = age_group)) +
  scale_fill_discrete(drop = FALSE)+                        # show all age groups in the legend, even those not present
  labs(
    title = "Age groups by Region")

ggplot(data = linelist_vars)+
  geom_bar(mapping = aes(x = region, fill = clade)) +
  scale_fill_discrete(drop = FALSE)+                        # show all age groups in the legend, even those not present
  labs(
    title = "clades by Region")

ggplot(data = linelist_vars)+
  geom_bar(mapping = aes(x = region, fill = clade_display)) +
  scale_fill_discrete(drop = FALSE)+                        # show all age groups in the legend, even those not present
  labs(
    title = "clades by Region")


ggplot(data = linelist_vars)+
  geom_bar(mapping = aes(x = age_group, fill = clade)) +
  scale_fill_discrete(drop = FALSE)+                        # show all age groups in the legend, even those not present
  labs(
    title = "Age groups by clades")

ggplot(data = linelist_vars)+
  geom_bar(mapping = aes(x = region, fill = nextclade_pango)) +
  scale_fill_discrete(drop = FALSE)+                        # show all age groups in the legend, even those not present
  labs(
    title = "clades by Region")

ggplot(data = linelist_vars)+
  geom_bar(mapping = aes(x = age_group, fill = nextclade_pango)) +
  scale_fill_discrete(drop = FALSE)+                        # show all age groups in the legend, even those not present
  labs(
    title = "lineages by age groups")
# analyzing the variants----

Variants <- read.csv("NPHL_variants_long_table.csv")


# filter for missense variants in gene S 
missense_variants <- Variants %>% 
  filter(EFFECT == "missense_variant" & GENE == "S" & FILTER == "PASS")


# count the number of times each variant occurs and sort the top 10 variants

missense_variants_count <- missense_variants %>% 
  count(POS, HGVS_P_1LETTER) %>% arrange(desc(n)) %>% top_n(30)


# select SAMPLE in missense_variants df where HGVS_P_1LETTER is equal to p.K417N
K417N <- missense_variants %>% filter(HGVS_P_1LETTER == "p.K417N" & AF > 0.5)



