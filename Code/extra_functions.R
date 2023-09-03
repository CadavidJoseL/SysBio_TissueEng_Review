# Additional functions and code

### Load packages
library(rentrez, quietly = T)
library(stringr, quietly = T)
library(dplyr, quietly = T)
library(tidyr, quietly = T)
library(ggplot2, quietly = T)
library(RColorBrewer, quietly = T)
library(ggvenn, quietly = T)

# Handy infix function to concatenate long lines (for queries) into a single string
# From https://stackoverflow.com/questions/6329962/split-code-over-multiple-lines-in-an-r-script
`%+%` = function(x,y) return(paste0(x,y))


# Read api key if it exists
api_key <- "../Data/api_key.txt"
if (file.exists(api_key)){
  api_key <- readChar(api_key, file.info(api_key)$size)
  set_entrez_key(api_key)
} else {
  api_key <- NULL
}


### Functions and parameters for plots

# Palettes and styles for plot elements
size_dot = 1 # dot size for plots
size_error = 0.4 # size (thickness of line) for error bars
size_text = 9
size_col = 0.1 # Thickness of line for barplots
size_legend = 0.3 # size of legend boxes in cm
size_line = 0.6 # thickness of lines
size_annotation = 3 # size for asterisks and other annotations
color_palette = "Set2" # Palette to be used. Choose from the Color Brewer palettes

# Function for standardizing plot aesthetics
add_theme <- function(plt){
  plt <- plt +
    theme_classic() +
    theme(text = element_text(size = size_text),
          axis.text.y = element_text(color = "black"),
          axis.text.x = element_text(color = "black"),
          legend.key.size = unit(size_legend, "cm"))
}

### Modify order of axes in plots with discrete levels. Add additional theme formatting
# without changing the theme itself (as done by add_theme). Top legend position

mod_axes <- function(plt, lims_x, lims_y){
  plt <- plt +
          scale_x_discrete(limits = lims_x) +
          scale_y_discrete(limits = lims_y) +
          theme(text = element_text(size = size_text),
                axis.text.y = element_text(color = "black"),
                axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
                legend.key.size = unit(size_legend, "cm"))

  return(plt)
}

### Wrapper functions for rentrez to make looping and processing a bit cleaner

# Function for querying and processing pubmed data so that we can run different queries
get_pubmed_data <- function(my_query, n_fetch = NULL, save_name = NULL, read_abstract = T, api_key = NULL, get_keywords = F){
  # Query
  my_pmid <- get_pubmed_ids(my_query, api_key)
  my_pmid$Count <- as.numeric(my_pmid$Count)
  print(my_pmid$Count)
  # Fetch data in chunks
  if (is.null(n_fetch)){n_fetch = my_pmid$Count}
  # Make sure we didn't request too many refs
  n_fetch <- min(n_fetch, my_pmid$Count)
  # Loop in blocks of 5000 references since that is what fetch_pubmed_data handles
  # nIDs could be set as an argument if we make this a function. Then we can either retrieve
  # everything or only a couple
  retstart <- 1
  retmax <- min(5000, n_fetch) # maximum number of PMID to fetch at a time
  text_df_full <- NULL
  while (retstart <= n_fetch){
    # Fetch data - starting index is default 0 so we subtract 1
    my_abstracts_xml <- fetch_pubmed_data(my_pmid, retstart =  retstart - 1, retmax = retmax)
    # Store Pubmed Records as elements of a list
    all_xml <- articles_to_list(my_abstracts_xml)
    # Perform operation (use lapply here, no further parameters) to convert each article
    # into a row of a dataframe. Skip authors to speed up
    text_df_dummy <- do.call(rbind, lapply(all_xml, article_to_df,
                                           max_chars = -1, getAuthors = F, getKeywords = get_keywords))

    # Concatenate dummy data frame to complete data frame
    text_df_full <- rbind(text_df_full, text_df_dummy)
    # Update indices for fetching
    retstart <- retstart + retmax
    retmax <- min(5000, n_fetch - retstart + 1)
  }

  # Save data frame if needed
  if (!(is.null(save_name))){
    write.table(text_df_full, file = paste0(save_name, '.txt'), sep = "\t")
  }
  # If we don't want abstracts (to reduce size of data), remove
  if (!read_abstract & !is.null(text_df_full)){
    text_df_full <- subset(text_df_full, select = -c(abstract))
  }
  # Return
  return(text_df_full)
}

# Function for batching date ranges - Pubmed only returns up to 10000 records at a time
# so we change the month ranges to get chunks of 10000 hits. We could just query month by month
# but the fewer queries, the faster this runs. This can be made faster by batching dates in a smarter way
get_IDs_batched <- function(query, year, data_base = "pubmed", pdat_char = "[pdat])"){
  # Output
  ids <- NULL
  # Month number
  months <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
  # Perform initial query to see if we can get away with passing the full year at once
  mix_query <- paste0('(',query,') AND (', year, ':', year, pdat_char) 
  # Search
  nIDs <- entrez_search(db = data_base, term = mix_query)$count
  #print(nIDs)
  if (nIDs < 10000){
    ids <- entrez_search(db = data_base, term = mix_query, retstart = 0, retmax = 10000)$ids
  }
  else {
    # Get month-by-month
    for (ii in seq(1,12)){
      #print(paste0("month: ", months[ii]))
      mix_query <- paste0('(',query,') AND (', year, '/', months[ii], ':', year, '/', months[ii], pdat_char) 
      search_query <- entrez_search(db = data_base, term = mix_query, retstart = 0, retmax = 10000)
      nIDs <- search_query$count
      # Split month by days if needed. if this split is not fine enough, we just leave it be before 
      # this becomes spaghetti code
      if (nIDs < 10000){
        if(nIDs >0){
          ids <- c(ids, search_query$ids)
          }
      }
      else {
        #print("split month")
        # Number of days in that month-year
        n_days <- lubridate::days_in_month(as.Date(paste0(year,"-",months[ii],"-01"), 
                                                   "%Y-%m-%d"))
        # Sequence of days - convert to date
        days_seq <- seq(1,n_days) 
        days_seq <- lapply(days_seq, function(x) ifelse(nchar(x) == 1, paste0("0",x), x)) %>% unlist()
        
        for (dd in seq(1, n_days-1)){
          #print(paste0("day: ", dd))
          mix_query <- paste0('(',query,') AND (', year, '/', months[ii], '/', days_seq[dd], ':', year, '/', months[ii], '/', days_seq[dd+1], pdat_char)
          search_query <- entrez_search(db = data_base, term = mix_query, retstart = 0, retmax = 10000)
          if(nIDs >0){
            ids <- c(ids, search_query$ids)
            }
        }
      } # end day split
    } #end months
  }
  # Return
  ids <- unique(ids) 
  if (length(ids) == 0){
    ids <- 0}
  return(ids)
}

###### Loop through time and categories

# Database should be pubmed (default) or pmc. If pmc, we get the PMID as well
# if we found nothing, we store the NA just to get a count of 0 for hits for a
# given year. This function is just a nice loop to keep things tidy

# Update: According to the E-utilities website, the entrez search can only retrieve
# 10000 records for pubmed. So, to avoid this issue, we iterate over months as well. 

get_queries_range <- function(query_table, add_query, years, data_base = "pubmed", 
                               name_save = NULL, n_subsample = NULL, dir_out){
  
  table_out <- NULL
  months <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
  
  for (year in years){
    #print(year)
    for (cat in seq(1, nrow(query_table))){
      # Define query
      query <- query_table$term[cat]
      category <- query_table$category[cat]
      #print(category)
      
      # Publication date might yield duplicates but we can handle those later one
      mix_query <- paste0('(',query,') AND (', add_query, ')') 
      # If data_base is PMC, we release constraints on title and abstract so that we
      # can find more hits in the text body. We also need to change other field IDs
      
      # Because there is no proximity search, some queries only work well because they are
      # constrained to the abstract. To maintain those, the queries have a "%" wildcard.
      # If the database is pubmed, we just remove it and keep the [tiab]. If pmc, we 
      # change the field to abstract
      if (data_base == "pmc"){
        mix_query <- gsub('\\[tiab%\\]','\\[ab\\]', mix_query)
        mix_query <- gsub('\\[tiab\\]','\\[article\\]', mix_query)
        # PMC handles truncation differently, so skip for now
        mix_query <- gsub('\\*','', mix_query)
        # Character for dates
        pdat_char <- "[pubdate])"
      }
      else {
        mix_query <- gsub('%','', mix_query)
        # Character for dates
        pdat_char <- "[pdat])"
      }
      # Query IDs: This function will split into months or days if necessary
      ids <- get_IDs_batched(mix_query, year, data_base = data_base, pdat_char = pdat_char)
      #print("Found " %+% length(ids))
      
      # Write dataframe - save as PMID or PMC depending on what we queried
      if (data_base == "pmc"){
        table_out <- rbind(table_out,
                           data.frame(year = year, category = category, PMC = ids))
      }
      else {
        table_out <- rbind(table_out,
                           data.frame(year = year, category = category, PMID = ids))
      }
      
    } # end category
  } # end year
  
  # Add category as a factor
  table_out$category <- as.factor(table_out$category)
  
  # Remove duplicates - will keep by default the lowest year
  table_out <- table_out[!duplicated(table_out),]
  
  # Save if name was given
  if (!is.null(name_save)){
    write.table(table_out, paste0(dir_out,name_save,".txt"), sep = "\t", quote = F, row.names = F)
  }
  
  # If number to subsample is given, subsample table and store. slice_sample is causing
  # issues since it wants the same number for all categories. An alternative would be to get a fixed
  # proportion of the whole thing, but that leads to overrepresentation of some categories
  if (!is.null(n_subsample)){
    table_sub <- NULL
    table_out <- table_out %>% na.omit()
    for (cat in unique(table_out$category)){ #using unique and not levels since some categories might be empty
      dummy <- table_out %>% filter(category == cat)
      table_sub <- rbind(table_sub,
                         dummy %>% slice_sample(n = min(n_subsample, nrow(dummy)))
      )
    }
    write.table(table_sub, paste0(dir_out,name_save,"_subsampled.txt"), sep = "\t", quote = F, row.names = F)
  }
  
  # Return
  return(table_out)
  
}


#### Make overlap plots between different tables: This function is a wrapper and calls
# get_queries_overlap for pairs of tables and returns a list with overlap tables 
# and plot objects. 

get_queries_overlap_lists <- function(table_list, save_name = NULL){
  list_out <- NULL
  n_tables <- length(table_list)
  table_names <- names(table_list)
  # Loop for all pairs of tables, including self-intersection
  for (ii in seq(1, n_tables)){
    for (jj in seq(ii, n_tables)){
      name_pair <- paste0(table_names[ii],"_", table_names[jj])
      # Handle self-intersection a bit differently
      if (ii == jj){
        list_out[[name_pair]] <- get_queries_overlap(table_list[[ii]])
        # Mod plot 
        list_out[[name_pair]]$plot <- list_out[[name_pair]]$plot +
                                      theme(text = element_text(size = size_text),
                                            axis.text.y = element_text(color = "black"),
                                            axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
                                            legend.key.size = unit(size_legend, "cm"))
      }
      else {
        list_out[[name_pair]] <- get_queries_overlap(table_list[[ii]], table_list[[jj]], use_jaccard = F)
      }
      
    }
  }
  # Make Venn diagram: Pass unique IDs (third column) after removing zeros and na
  list_out[["Venn"]] <- ggvenn(lapply(table_list, 
                                      FUN = function(x) x[x[,3] != 0,3] %>% na.omit() %>% unique()),
                               show_percentage = F,
                               fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
                               stroke_size = 0.5, set_name_size = 4)
  
  # If a save name was given, save list
  if (!is.null(save_name)){
    saveRDS(list_out, file = save_name, compress = TRUE, refhook = NULL)
  }
  # Return
  return(list_out)
}

#### Get overlap in PMID between query categories

get_queries_overlap <- function(table1, table2 = NULL, use_jaccard = T){
  # Single input is to look for overlaps within table
  if (is.null(table2)){
    table2 <- table1
    }
  # Clean inputs and change the third column to ID
  colnames(table1)[3] <- "ID"
  colnames(table2)[3] <- "ID"
  table1 <- table1 %>% na.omit() %>% filter(ID != "0")
  table2 <- table2 %>% na.omit() %>% filter(ID != "0")
  # Sanity check: Make sure categories are factors
  table1$category <- as.factor(table1$category)
  table2$category <- as.factor(table2$category)

  # Calculate overlaps
  table_out <- NULL
  for (l1 in levels(table1$category)){
    # Get subset 1
    sub1 <- subset(table1, category == l1)

    for (l2 in levels(table2$category)){
      # Get subset 2
      sub2 <- subset(table2, category == l2)

      # Get counts in each set - get unique PMIDs in case something is repeated
      # because of the pubdate
      n1 <- sub1$ID %>% unique() %>% length()
      n2 <- sub2$ID %>% unique() %>% length()
      n_both <- intersect(sub1$ID, sub2$ID) %>% length()
      n_either <- n1 + n2 - n_both

      table_out <- rbind(table_out,
                         data.frame(category1 = l1, category2 = l2,
                                    count1 = n1,
                                    count2 = n2,
                                    count_both = n_both,
                                    count_either = n_either,
                                    jaccard = n_both/n_either))

    }
  }

  if (use_jaccard){
    # Use Jaccard index
  # Make overlap plot - remove high jaccard values that will arise when checking self-intersections
  plt <- ggplot(data = table_out %>% filter(count_both > 0 & jaccard < 0.999),
                aes(x = category1, y = category2, size = 100*jaccard)) +
             geom_point() +
             labs(x = NULL, y = NULL)
  }
  else {
    # Use counts
    nMax <- max(table_out$count_both) %>% log10() %>% ceiling()

    plt <- ggplot(data = table_out %>% filter(count_both > 0),
                aes(x = category1, y = category2, color = log2(count_both))) +
             geom_point(size = 4) +
            # scale_color_continuous(type = "viridis",
             #                       breaks = trans_breaks("log2", function(x) 2^x),
              #                      labels = trans_format("log2", math_format(2^.x)),
               #                     limits = c(1, nMax)) +
              #scale_color_continuous(type = "viridis",
               #                      labels = trans_format("identity", math_format(2^.x))) + # trick to get the 2^x in label
             scale_color_viridis_c() +
             labs(x = NULL, y = NULL) +
             labs(color = "log10(Count)")
  }

  # Return output
  return(list(table = table_out, plot = plt))

}

#### Plot trends - a wrapper function. Assumes x axis is "year"

plot_trends <- function(data, var_count = "PMID", var_group1 = "PMID", var_group2 = NULL){
  
  # We can group based on a category, or not. Using an if statement for now. Can be made more general
  # I suppose
  if (!is.null(var_group2)){
    plot <- data %>% group_by(!!sym(var_group1), !!sym(var_group2)) %>% 
              summarize(year = min(year)) %>%
              group_by(!!sym(var_group2), year) %>% summarize(n = n()) %>%
              ggplot(aes(x = year, y = n, color = !!sym(var_group2))) +
              geom_line(size = size_line) +
              scale_y_log10() +
              facet_wrap(facets = vars(!!sym(var_group2)), ncol = 5) +
              labs(x = "Year", y = "Count", color = NULL) +
              theme(legend.position = "none")
  }
  else {
    plot <- data %>% group_by(!!sym(var_group1)) %>% 
      summarize(year = min(year)) %>%
      group_by(year) %>% summarize(n = n()) %>%
      ggplot(aes(x = year, y = n)) +
      geom_line(size = size_line) +
      scale_y_log10() +
      labs(x = "Year", y = "Count") +
      theme(legend.position = "none")
  }
  
  return(plot)
  
}


# Faster version - fetching multiple papers at once and processing as batch. Batches of
# 100 seem fine. Parsing abstract myself to deal with abstracts with subsections (background, methods, etc)
get_title_abstract <- function(PMID, save_name = NULL){
  
  numel <- length(PMID)
  pmid_title_abstract <- NULL

  for (retstart in seq(1,numel,100)){
    # Fetch info (abstracts and titles)
    PMID_batch <- PMID[retstart:min(numel, retstart + 99)]
    xml_text <- entrez_fetch(db = "pubmed", id = PMID_batch, retmode = "abstract", rettype = "xml") 
    
    # Some papers won't have an abstract but might still be relevant based on title
    # which is why they aren't exluded. But here that will cause problems.
    # We split the xml into papers and then process each. Sometimes there is a book in there, so we split
    articles <- str_extract_all(xml_text, "(?<=<PubmedArticle>).+?(?=</PubmedArticle>)|(?<=<PubmedBookArticle>).+?(?=</PubmedBookArticle>)") %>% unlist()
    # For some papers the above line doesn't work and I have no idea why (the pattern matches well
    # in other regex debuggers). If that's the case, we can brute force it and manually find the 
    # PubmedArticle delimiters
    if (length(articles) != length(PMID_batch)){
      articles <- NULL
      pos_start <- str_locate_all(xml_text,"<PubmedArticle>|<PubmedBookArticle>")[[1]]
      pos_end <- str_locate_all(xml_text,"</PubmedArticle>|</PubmedBookArticle>")[[1]]
      for (ii in seq(1, length(PMID_batch))){
        articles[ii] <- substr(xml_text, pos_start[ii,2]+1, pos_end[ii,1]-1)
      }
    }
    
    # For loop not so elegant but works
    abstract_list <- NULL
    title_list <- NULL
    for (ii in seq(1, length(articles))){
      # Extract abstracts (. matches any character, + looks for any number (1+), ? sets it to lazy - non greedy, to find shorter matches)
      # Gotta learn the REGEX tricks. The lookarounds omit the <Abstract> part
      abstract <- str_extract_all(articles[ii], "(?<=<Abstract>).+?(?=</Abstract>)") %>% unlist()
      # Remove<AbstractText> (* looks for any number (0+))
      abstract <- gsub("<AbstractText.*?>|</AbstractText.*?>", "", abstract)
      # Copyright info sometimes sneaks into abstract
      abstract <- gsub("<CopyrightInformation>.*?</CopyrightInformation>","", abstract)
      if(length(abstract) == 0){abstract = NA}
      # Get title
      title <- str_extract_all(articles[ii], "(?<=<ArticleTitle>).+?(?=</ArticleTitle>)") %>% unlist()
      if(length(title) == 0){title = NA}
      # Append to list
      abstract_list[ii] <- abstract
      title_list[ii] <- title
    }

    # Put together - concatenate
    pmid_title_abstract <- rbind(pmid_title_abstract,
                                 data.frame(PMID = PMID_batch,
                                      abstract = abstract_list,
                                      title = title_list))
    
    # Save at end of batch if name was given
    if (!is.null(save_name)){
      write.table(pmid_title_abstract, save_name, row.names = F, sep = "\t")
    }
  }
  
  return(pmid_title_abstract)
}

##### Functions for converting from PMC to PMID.
# Function entrez_link works, but it's slow, so it's simpler to get
# the document summaries, and we get the publication type from the retrieved PMIDs

# Get PMIDs by getting article summary (faster then entrez_links)
get_PMID_PMC <- function(PMCs, do.clean = T){

  summary <- entrez_summary(db="pmc", id = PMCs, always_return_list = T)
  # Extract PMIDs - some PMC won't map and we likely want to get rid of them
  # PMID seems to be first and then doi and PMC but these last two are not always
  # in that position
  PMIDs <- sapply(summary, FUN = function(x) x[["articleids"]][["value"]][1])
  # Build dataframe
  dictionary <- data.frame(PMC = PMCs, PMID = PMIDs)

  # If do.clean, we will remove reviews and other articles, and things not in english
  if (do.clean){
    # Drop unmatched IDs
    dictionary <- dictionary %>% filter(PMIDs != 0)
    # Get summary of PMIDs to get pub type
    # This will ignore things that don't map to PMID, so we'll have to wrangle a bit
    summary <- entrez_summary(db="pubmed", id=PMIDs, always_return_list = T)
    pubtype <- sapply(summary, FUN = function(x) x[["pubtype"]])
    lang <- sapply(summary, FUN = function(x) x[["lang"]])
    # Find if the article is just a journal article or a review. By default everything seems to
    # be indexed as Journal Article
    is.article <- sapply(pubtype, FUN = function(x) all(x == 'Journal Article'))
    dictionary <- dictionary[is.article & lang == "eng",]
  }
  return(dictionary)
}

# We can only do batches of 300 at a time or so, so this function loops batches.
# Could be a single function, but I think two functions is easier to read
get_PMID_PMC_batched <- function(PMCs, do.clean = T, save_name = "../Data/dictionary.txt"){
  # Remove IDs with 0 (that come from empty queries)
  PMCs <- PMCs[PMCs != "0"]
  # Check if a dictionary is already in file. If yes, load and only process
  # new IDs
  if (file.exists(save_name)){
    dictionary <- read.table(save_name, header = T)
    PMCs <- PMCs[!(PMCs %in% dictionary$PMC)]
    print("Previous file found")
  } else {
    dictionary = NULL
  }
  # Pass chunks of 300 at a time
  batch_size <- 300
  numel <- length(PMCs)
  ii_min <- 1
  ii_max <- batch_size
  while (ii_min <= numel){
    # Cap ii_max
    ii_max <- min(ii_max, numel)
    batch <- PMCs[ii_min:ii_max]
    dummy <- get_PMID_PMC(batch, do.clean = do.clean)
    # Concatenate
    dictionary <- rbind(dictionary, dummy)
    # Update indices
    ii_min <- ii_min + batch_size
    ii_max <- ii_max + batch_size
    # If save_name, save!
    if (!is.null(save_name)){
      write.table(dictionary, file = save_name, sep = "\t", row.names = F, quote = F)
    }
  }
  rownames(dictionary) <- NULL
  return(dictionary)
}


##### Function for looping through article PMIDs and ask for manual inspection
# We display the title of the paper one by one and ask the user to input whether 
# the paper should be reviewed in more detail later on (by hand). This could be done
# manually but it's nice to loop titles in R

verify_titles <- function(PMID, file_prechecked = NULL, file_store = NULL, batch_size = NULL){
 
  # Make data frame with unique PMID (the full table has repeated PMID that pop up)
  data <- data.frame(PMID = unique(PMID)) %>% filter(PMID != 0)
  # Add column for verification
  data$is.good <- NA
  data$is.review <- NA
  # Number of elements
  numel <- nrow(data)
  
  # If given, read pre-checked IDs so that we don't have to repeat
  if (!is.null(file_prechecked)){
    # Check if the name given is valid
    if (file.exists(file_prechecked)){
      print("Loading pre-existent file...")
      data_checked <- read.table(file_prechecked, sep = "\t", header = T)
      # Nice prompt
      print(paste0(sum(PMID %in% data_checked$PMID),
                   ' out of ',
                   length(PMID),
                   ' pre-checked IDs!'))
      # Exclude PMIDs that we've checked
      data <- data %>% filter(!(PMID %in% data_checked$PMID))
      numel <- nrow(data)
      # Put new dataframe on top of old and we only alter the new part
      data <- rbind(data, data_checked)
    }
    else {print("Invalid pre-checked file name, starting from scratch!!")}
  }
  else {print("No pre-checked file name, starting from scratch!!")}
  
  # Loop for all PMID (up to batch_size) and display title and prompt user for Y/N
  if (!is.null(batch_size)){numel <- min(numel,batch_size)}
  # Loop is there is something to loop for
  if (numel > 0){
    for (ii in seq(1, numel)){
      print(paste0(ii, ' out of ',numel))
      # Check if the abstract of the paper says "review"
      data$is.review[ii] <- entrez_search(db = "pubmed", 
                                                        term = paste0(data$PMID[ii],
                                                                      "[uid] AND review[tiab]"))$count > 0
      # Print title and ask for confirmation
      print(entrez_summary(db = "pubmed", id = data$PMID[ii])$title)
      # Query user
      data$is.good[ii] <- readline(prompt = "Check? (y/n):")
      
      # Save table every 10 entries - eliminate whatever hasn't been filled
      if ((ii %% 10) & !is.null(file_store)){
        write.table(data %>% na.omit(), file = file_store,
                    quote = F, sep = "\t", row.names = F)
      }
    }
  }
  
  # Store: Remove anything that is NA (not checked)
  numel <- nrow(data)
  data <- data %>% na.omit()
  if(!is.null(file_store)){
    write.table(data, file = file_store,
                quote = F, sep = "\t", row.names = F)
  }
  # Prompt
  print(paste0(numel - nrow(data), ' out of ', numel, ' left to check!'))
  # Return as invisible to not clutter window
  invisible(data)
}
