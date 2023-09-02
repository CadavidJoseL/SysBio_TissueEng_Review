# Load stop words
library(tidytext)
library(ggrepel)
library(cowplot)
library(scales)
data(stop_words)

source("extra_functions.R")

# Get IDs of papers that mention the word "systems" in their title
IDs <- NULL
for (year in seq(1999,2023)){
term <- paste0('systems[title] AND review[pt] AND eng[lang] AND hasabstract')

ID_new <- get_IDs_batched(query = term, year = year, data_base = "pubmed")
IDs <- c(IDs, ID_new)
}




# Get titles and abstracts and look for words next to "systems"
titles_list <- NULL
for (ii in seq(1, length(IDs), 100)){
  
  search <- entrez_summary(db = "pubmed", id = IDs[ii: min(ii+99, length(IDs))])
  
  titles <- lapply(search, FUN = function(x) x$title %>% tolower()) %>% unlist() 
  

  titles_list <- c(titles_list, titles)
  
}


words <- titles_list %>% str_replace_all("[[:punct:]]"," ") %>% str_replace_all("  "," ") %>%
        str_extract_all("(?<=systems ).*?(?= )") %>% unlist()

words_df <- data.frame(word = words) %>% 
            filter(!(word %in% stop_words$word)) %>%
            filter(!(word %in% c(" ", ""))) %>%
            group_by(word) %>% summarize(n = n())

words_df$word <- words_df$word %>% str_replace_all(" ","")

# Check whether pairs (systems + word) are indeed in the titles or if they were only
# included after removing a punctuation sign

words_df$is.title <- sapply(words_df$word, FUN = function(x) str_count(titles_list, paste("systems",x)) %>% sum())

# By manual inspection the interesting terms are:

systems_terms <- data.frame(term = c("biology", "pharmacology", "medicine", "neuroscience", "genetics",
                 "metabolic engineering", "immunology", "toxicology", "vaccinology",
                 "chemistry", "physiology", "pathology", "biotechnology", "microbiology",
                 "biomedicine", "biochemistry", "glycobiology", "epidemiology", "genomics",
                 "oncology", "serology", "neurobiology", "biocatalysis", "virology"))

table <- NULL

# Get article counts for each term
for (ii in seq(1, nrow(systems_terms))){
  query <- paste0('"systems ', systems_terms$term[ii], '"[tiab] AND eng[la]')
  query_PMC <- paste0('"systems ', systems_terms$term[ii], '"[article]')
  
  for (year in seq(1999, 2022)){
    
  table <- rbind(table,
                          data.frame(term = systems_terms$term[ii],
                                      year = year,
                                      count = entrez_search(db = "pubmed", term = query, mindate = year, maxdate = year)$count,
                                      count_PMC = entrez_search(db = "pmc", term = query_PMC, mindate = year, maxdate = year)$count))


  }
}

table$term <- as.factor(table$term)
table <- table %>% group_by(term) %>% mutate(count_sum = cumsum(count), count_sum_PMC = cumsum(count_PMC)) %>% ungroup()

# Write table
#write.table(table, file = "Results2023/systems_words_table.txt", quote = F, row.names = F, sep = "\t")

table <- read.table("Results2023/systems_words_table.txt", header = T, sep = "\t")

# PLot
plt1 <- table %>% group_by(term) %>% filter(count_sum > 1) %>% ungroup() %>%
        ggplot(aes(x = year, y = count_sum, color = term)) + 
        geom_line(linewidth = 1.2) + 
        facet_wrap(facets = vars(term), ncol = 4) + 
        theme(legend.position = "none") + 
        scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x))) +
        scale_x_continuous(breaks = c(2000, 2010, 2020)) +
        labs(y = "Cumulative number of articles in Pubmed", x = "Year") +
        theme(text = element_text(size = size_text),
          axis.text.y = element_text(color = "black"),
          axis.text.x = element_text(color = "black"),
          legend.key.size = unit(size_legend, "cm"))



# SUmmary of table (year in which 5 total papers were crossed, and total number of papers)
table_sum <- table %>% group_by(term) %>% summarize(time_5 = year[count_sum > 5][1], 
                                                    total = max(count_sum), 
                                                    total_PMC = max(count_sum_PMC))

plt2 <- table_sum %>% mutate(term = fct_reorder(term, total)) %>%
        ggplot(aes(x = total, y = term)) + #, fill = log10(total)/(2022 - time_5))) +
        geom_col(fill = "steelblue") +
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x))) +
        scale_fill_viridis_c() + 
        labs(y = NULL, x = "Total count") +
        theme(text = element_text(size = size_text),
              axis.text.y = element_text(color = "black"),
              axis.text.x = element_text(color = "black"),
              legend.key.size = unit(size_legend, "cm"))


plt_grid <- plot_grid(plt1, plt2, 
                      labels = c('a', 'b'), label_size = size_text, rel_widths = c(2, 1))

ggsave("systems_words_plot.pdf", plot = plt_grid, units = "cm", width = 17, height = 12)
ggsave("systems_words_plot.png", plot = plt_grid, units = "cm", width = 17, height = 12)

# Do dot plot to keep simple
plt3 <- table_sum %>% 
        ggplot(aes(x = time_5, y = total_PMC, label = term, color = log10(total)/(2022 - time_5))) + 
        geom_point(size = 2*size_dot) + 
        geom_text_repel(color = "black", size = size_annotation) + 
        scale_y_log10() + 
        scale_color_viridis_c()

table_recent <- table %>% filter(year >= 2017) %>% group_by(term) %>% 
                mutate(count_sum = cumsum(count), count_sum_PMC = cumsum(count_PMC)) %>% 
                summarize(total = max(count_sum), total_PMC = max(count_sum_PMC))
