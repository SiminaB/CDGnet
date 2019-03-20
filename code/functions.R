##function that maps to various attributes
source("http://michael.hahsler.net/SMU/ScientificCompR/code/map.R")
##see also http://michael.hahsler.net/SMU/LearnROnYourOwn/code/igraph.html for examples

##function that inputs a list of genes and outputs a vector matching the names to the standardized names
##for the genes that were matched up using DGI (and just returns the original names for the ones that were not)
gene_names_standard <- function(genes)
{
  stg <- searchTermSummary(queryDGIdb(genes))
  ##make list with standardized gene names as names, initial gene names as entries
  gene_names_list <- as.list(as.character(stg$SearchTerm))
  names(gene_names_list) <- as.character(stg$Matches)
  ##now split on comma
  gene_names_sep_list <- lapply(gene_names_list, 
                                function(g){strsplit(g,", ")[[1]]})
  ##now create vector for id matching
  ##first repeat the standard names the correct number of times
  stand_names <- rep(names(gene_names_sep_list), sapply(gene_names_sep_list, length))
  name2stand_name <- stand_names
  names(name2stand_name) <- unlist(gene_names_sep_list)
  
  ##now print out all the genes that do not get mapped to standard names
  print("Genes that do not get mapped to standard names:")
  not_mapped <- setdiff(genes, names(name2stand_name))
  print(not_mapped)
  ##remove NA
  not_mapped <- not_mapped[!is.na(not_mapped)]

  ##for all these genes, add in the IDs at the end (except for NAs)    
  name2stand_name <- c(name2stand_name, not_mapped)
  names(name2stand_name) <- c(unlist(gene_names_sep_list), not_mapped)
  
  name2stand_name
}

##function that adds other genes to the molecular profiling (MP) data frame
MP_add_genes <- function(MP_df, add_genes)
{
  keep_genes_not_MP <- setdiff(add_genes, MP_df$Gene_protein)
  MP_plus <- MP_df
  if(length(keep_genes_not_MP) > 0)
  {
    MP_plus <- rbind(MP_df, data.frame(Gene_protein = keep_genes_not_MP,
                                       Data_type = NA,
                                       Alteration = NA))
  }
  MP_plus
}

##function that remove wild type alterations from an MP data frame
##MP_df = data frame with molecular alterations
##wild_type = ways in which the wild type can be coded in MP_df
remove_WT <- function(MP_df, 
                      wild_type = c("wild type","wildtype","wild","wt"))
{
  MP_df[!(tolower(MP_df$Alteration) %in% wild_type),]
}

##function that checks whether any of the drugs have biomarkers that are wild-type
##if any of the biomarkers are wild-type, then check the MP report
##if they are not wild-type there, then remove them from the list
##drugs_bio = data frame with drugs and biomarkers
##MP_df = data frame with molecular alterations
check_WT_biomarkers <- function(drugs_bio,MP_df){
  
  ##if any of these are biomarkers for "wild type", check against that in MP_df data
  biomarkers_WT <- unique(drugs_bio$Gene.Protein[drugs_bio$Alteration == "wild type"])
  
  if(length(biomarkers_WT)>=1)
  {
    ##what is the alteration for these biomarkers in the MP_df report? If it is not "wild type" or "WT," remove it from the list
    biomarkers_WT_MP_alt <- MP_df %>% dplyr::filter(Gene_protein %in% biomarkers_WT)
    biomarkers_WT_MP_alt <-
      biomarkers_WT_MP_alt[!(tolower(biomarkers_WT_MP_alt$Alteration) %in%
                               c("wild type","wildtype","wild","wt")),]
    keep_MP_bio <- setdiff(MP_df$Gene_protein, biomarkers_WT_MP_alt$Gene_protein)
    
    drugs_bio <- dplyr::filter(drugs_bio, Gene.Protein %in% keep_MP_bio)
  }

  drugs_bio
}

##function that calculates the shortest paths from inputs to everything else that
##is in a pathway but not in MP
##calls on calculate_shortest_path
calculate_shortest_path_from_inputs <- function(cancer_path,
                                                MP,
                                                inputs_in_path)
{
  ##convert the connections to an adjacency matrix
  g <- graph.data.frame(cancer_path[,c("from","to")])
  
  ##get the length of the shortest path from the inputs to all the other genes
  all_genes <- unique(c(as.character(cancer_path$from),
                        as.character(cancer_path$to)))
  inputs <- unique(as.character(MP$Gene_protein))
  all_genes_no_inputs <- setdiff(all_genes, inputs)
  
  shortest_paths_from_inputs <- calculate_shortest_path(inputs_in_path,
                                                        all_genes_no_inputs,
                                                        g)
}
  
##function that calculates the shortest paths from some inputs to some outputs
##returns a matrix with the inputs as rows, the outputs as columns, and the shortest path length as entries
##inputs and outputs are both character vectors
##igraph_obj is created by converting a data frame of edges using graph.data.frame
calculate_shortest_path <- function(inputs, outputs, igraph_obj)
{
  shortest_paths_from_inputs <- 
    matrix(nrow=length(inputs), ncol=length(outputs))
  rownames(shortest_paths_from_inputs) <- inputs
  colnames(shortest_paths_from_inputs) <- outputs
  
  for(input in inputs)
  {
    shortest_paths_from_inputs[input,] <-
      sapply(outputs,
             function(to,from,igraph_obj){length(attributes(shortest_paths(igraph_obj,from,to)$vpath[[1]])$names)},
             input,igraph_obj)
  }
  shortest_paths_from_inputs
}

##function that pastes the shortest paths between some inputs and some outputs
##returns a character vector with the paths marked with "->", separated by ";" if there are multiple paths
##inputs, inputs_in_path and outputs are all character vectors
##(inputs_in_path are downstream of the oncogenes)
##igraph_obj is created by converting a data frame of edges using graph.data.frame
paste_shortest_path <- function(inputs, inputs_in_path, outputs, igraph_obj)
{
  shorts <- rep(NA, length(outputs))
  names(shorts) <- outputs
  for(gene in outputs)
  {
    if(gene %in% inputs){
      shorts[gene] <- gene
    } else {
      gene_path_from_input <- list()
      for(input in inputs_in_path)
      {
        gene_path_from_input[[input]] <-
          paste(attributes(shortest_paths(igraph_obj, input, gene)$vpath[[1]])$names,
                sep="", collapse=" -> ")
      }
      gene_path_from_input_vect <- unlist(gene_path_from_input)
      gene_path_from_input_vect <-
        gene_path_from_input_vect[gene_path_from_input_vect != ""]
      shorts[gene] <- paste(gene_path_from_input_vect,
                            sep="",collapse=";")
      
    }
  }
  shorts
}

##function that prettifies output of drugs_mat
pretty_drugs_mat <- function(drugs_mat)
{
  colnames(drugs_mat)[4:6] <- c("Predicted effect",
                                "Data type",
                                "Alteration")
  
  ##remove/simplify some things further
  drugs_mat <- drugs_mat[,c("Drug","Gene.Protein","Data type",
                            "Alteration","Path","Disease","Predicted effect")]
  drugs_mat
}

##function to get nodes downstream from given nodes
##starts from the entire pathway - calls calculate_shortest_path_from_inputs and
##get_downstream_nodes
get_downstream_from_inputs <- function(cancer_path,
                                       MP,
                                       inputs_in_path)
{
  ##get the length of the shortest path between the genes/proteins in
  ##inputs_in_path and all the other genes (except for the ones that are WT in the
  ##actual MP)
  shortest_paths_from_inputs <- 
    calculate_shortest_path_from_inputs(cancer_path,
                                        MP,
                                        inputs_in_path)
  
  ##keep genes that have a shortest path > 0 from any of the oncogenes that are also positive in the MP report
  ##get just these genes
  ##want genes downstream of _all_ the inputs
  ##first get list of all downstream things
  get_downstream_nodes(inputs_in_path,
                       shortest_paths_from_inputs)
}

##function to get nodes downstream from given nodes (which are the inputs)
##inputs is a character vector, shortest_paths_from_inputs is matrix generated by calculate_shortest_path
##output of function is a list of downstream nodes for each input node
get_downstream_nodes <- function(inputs, shortest_paths_from_inputs)
{
  downstream_from_inputs <- list()
  for(input in inputs){
    downstream_from_inputs[[as.character(input)]] <-
      names(which(shortest_paths_from_inputs[input,]>0))
  }
  ##since want the shortest path from _any_ of the inputs, take the column sums
  sum_shortest_path <- colSums(shortest_paths_from_inputs)
  ##the ones that are downstream from any of the inputs have the sum > 0
  downstream_from_inputs <- names(which(sum_shortest_path>0))
  downstream_from_inputs
}
  
##filter data frames based on the gene/protein being in keep_genes
##keep only the unique drug names
filter_drugs_select_genes <- function(drugs_genes_df,
                                      keep_genes)
{
  dplyr::filter(drugs_genes_df,
                Gene.Protein %in% keep_genes)
}

##function that creates the data frame specifying the edges from a data frame like drugs_PO_FDA_keep
make_df_edges_targeted <- function(df_targeted, subtype = "sensitivity")
{
  if(nrow(df_targeted)>0)
  {
    df_edges_targeted <- data.frame(from=as.character(df_targeted$Drug),
                                    to=as.character(df_targeted$Gene.Protein),
                                    subtype=subtype,
                                    direction=-1)
  } else {
    df_edges_targeted <- data.frame(from=character(),to=character(),
                                    subtype=character(),direction=character())
  }
  df_edges_targeted
}

##function that combines all drugs from all categories
##Cat1, Cat2, Cat3, Cat4 correspond to the data frames from the 4 categories (Cat 4 may be null)
combine_drugs <- function(Cat1, Cat2, Cat3, Cat4 = NULL)
{
  drugs <- c()
  if(ncol(Cat1)>1)
  {
    drugs <- c(drugs, as.character(Cat1$Drug))
  }
  if(ncol(Cat2)>1)
  {
    drugs <- c(drugs, as.character(Cat2$Drug))
  }
  if(ncol(Cat3)>1)
  {
    drugs <- c(drugs, as.character(Cat3$Drug))
  }
  if(!is.null(Cat4))
  {
    if(ncol(Cat4)>1)
    {
      drugs <- c(drugs, as.character(Cat4$Drug))
    }
  }
  unique(drugs)
}

##function that generates pathway plots for categories 3 and 4
##edge is the matrix of edges, drugs and inputs are character vectors
graph_pathways_cats_3_4 <- function(edges, drugs, inputs)
{
  graph <- graph_from_data_frame(edges)
  V(graph)$size <- 3
  V(graph)$label.cex <- 0.4
  V(graph)$color <- "#1f77b4"
  
  V(graph)$color[V(graph)$name %in% drugs] <- "#ff7f0e"
  V(graph)$color[V(graph)$name %in% inputs] <- "#9467bd"
  
  # ##get smallest distance between each vertex and each drug
  # ##only need to consider drugs that are vertex names
  # drugs <- intersect(drugs, V(graph)$name)
  # small_dist_drug <- calculate_shortest_path(drugs, V(graph)$name, graph)
  # ##replace 0 with number of vertices (probably means that there's no direct path)
  # small_dist_drug[small_dist_drug == 0] <- length(V(graph)$name)
  # ##get smallest distance to any drug
  # small_dist_any_drug <- apply(small_dist_drug, 2, min)
  # ##put in a very small number for the actual drugs
  # small_dist_any_drug[drugs] <- 1
  # write.csv(small_dist_drug, file="small_dist_drug.csv")
  # write.csv(small_dist_any_drug, file="small_dist_any_drug.csv")
  
  plot(graph,
       edge.arrow.size=0,
       layout=layout_with_fr,
       edge.arrow.size=0.5, 
       vertex.label.cex=0.7,##0.5*map(1/small_dist_any_drug[V(graph)$name],c(1,1.5)), 
       vertex.label.family="Helvetica",
       vertex.label.font=2,
       vertex.shape="circle", 
       vertex.label.dist=1,
       vertex.size=4,##map(1/small_dist_any_drug[V(graph)$name],c(1,10)), 
       vertex.label.color="black", 
       edge.width=0.5)
}

##function that makes data frame of nodes
##assoc_df is the data frame of drug-gene/protein assocations
make_node_df <- function(assoc_df)
{
  node_names <- unique(c(as.character(assoc_df$from),
                         as.character(assoc_df$to)))
  nodes_df <- data.frame(Name = node_names,
                         Group = "Gene/Protein")
  nodes_df$Group <- as.character(nodes_df$Group)
  nodes_df$Group[grep("\\(",as.character(nodes_df$Name))] <- "Drug"
  
  nodes_df
}

##function that gets connections from one or more pathways which are within a list like the one of KEGG pathways
##list_of_paths is the list of pathways
##IDs is a vector of IDs for the names of the pathways that should be pulled out
##subtype may specify "activation" or just be NULL, in which case everything is included
##output is a data frame of connections
get_connections_paths <- function(list_of_paths, IDs, subtype=NULL)
{
  cancer_path <- c()
  
  for(ID in IDs)
  {
    cancer_path <- rbind(cancer_path,
                         list_of_paths[[ID]])
  }
  ##filter based on subtype (if applicable)
  if(!is.null(subtype))
  {
    cancer_path <- dplyr::filter(cancer_path,
                                 subtype == subtype)
  }
  
  ##remove any duplicates
  cancer_path <- cancer_path[!duplicated(cancer_path),]
  
  cancer_path
}

##get edges for categories 3 and 4
##Cat3, Cat4 are the reactive objects corresponding to categories 3 and 4
##Cat4 may be null
get_edges_cats_3_4 <- function(Cat3, Cat4 = NULL)
{
  all_edges <- as.data.frame(Cat3$edges)
  
  if(!is.null(Cat4))
  {
    all_edges <- rbind(all_edges,
                       as.data.frame(Cat4$edges))
  }
  
  drugs_mat <- Cat3$drugs_mat
  
  if(!is.null(Cat4))
  {
    drugs_mat <- rbind(drugs_mat,
                       as.data.frame(Cat4$drugs_mat))
  }
  
  ##drugs that are not used should be the set difference between everything in an edge and
  ##everything in drugs_mat that is either a drug or a gene/protein
  # drugs_excluded <- unique(setdiff(c(as.character(all_edges[,1]), 
  #                                    as.character(all_edges[,2])),
  #                                  c(as.character(drugs_mat[,"Drug"]),
  #                                    as.character(drugs_mat[,"Gene or Protein"]))))
  # 
  # all_edges <- all_edges[!(all_edges$Source %in% drugs_excluded) &
  #                          !(all_edges$Target %in% drugs_excluded),]
  # 
  all_edges <- all_edges[!duplicated(all_edges),]
  
  all_edges  
}

##for category 3 and 4 recommendations, take out the drugs already recommended above
##(i.e. recommended in categories 1 and 2 for category 3 and in categories 1, 2, or 3 
##for category 4)
##if Type3_df == NULL, then this means just take out category 1 and 2 drugs + anything that has 
##(otherwise, also take out category 3 drugs)
##drugs_biomarkers_targets is the initial input data frame, from which known recommendations
##will be removed
remove_known_recs <- function(drugs_biomarkers_targets,Type1_df, Type2_df, Type3_df = NULL)
{
  known_recs <- NULL
  
  if(ncol(Type1_df) > 1)
  {
    known_recs <- Type1_df[,c("Gene or Protein","Type","Alteration")]
  }
  if(ncol(Type2_df) > 1)
  {
    known_recs <- Type2_df[,c("Gene or Protein","Type","Alteration")]
  }
  if(ncol(Type1_df) > 1 & ncol(Type2_df) > 1)
  {
    known_recs <- rbind(Type1_df[,c("Gene or Protein","Type","Alteration")],
                        Type2_df[,c("Gene or Protein","Type","Alteration")])
  }
  if(!is.null(Type3_df))
  {
    if(ncol(Type3_df)>1)
    {
      known_recs <- rbind(known_recs[,c("Gene or Protein","Type","Alteration")],
                          Type3_df[,c("Gene or Protein","Type","Alteration")])
    }
  }
  if(!is.null(known_recs))
  {
    drugs_biomarkers_targets <- 
      drugs_biomarkers_targets[!(drugs_biomarkers_targets$Gene.Protein %in% known_recs[,1]),]
  }
  drugs_biomarkers_targets
}
  
##function to add in the shortest path to the drugs_biomarkers_targets data frame
add_shortest_path <- function(drugs_biomarkers_targets, 
                              inputs, inputs_in_path,
                              cancer_path)
{
  ##convert the connections to an adjacency matrix
  g <- graph.data.frame(cancer_path[,c("from","to")])
  
  ##get shortest path from input genes to keep_genes
  shorts <- paste_shortest_path(inputs, inputs_in_path,
                                unique(drugs_biomarkers_targets$Gene.Protein),
                                g)
  
  dplyr::mutate(drugs_biomarkers_targets, Path=shorts[Gene.Protein])
}

dataParser <- function(input, fda, fda_other, rec) {
  
  input <- distinct(input)
  fda <- distinct(fda)
  fda_other <- distinct(fda_other)
  rec <- distinct(rec)
  
  ###
  ## groups denotes the id and title for each section in the sankey diagram.
  ###
  json <- list(
    "links"=list(), 
    "nodes"=list(), 
    "groups"=list(
      "mp"=list(
        "type"="process",
        "title"="Molecular Profile",
        "bundle"=NULL,
        "id"="mp",
        "nodes"=c(),
        "def_pos"=NULL
      ),
      "fda"=list(
        "type"="process",
        "title"="FDA Approved Drugs",
        "bundle"=NULL,
        "id"="fda",
        "nodes"=c(),
        "def_pos"=NULL
      ),
      "it"=list(
        "type"="process",
        "title"="Inferred Targets",
        "bundle"=NULL,
        "id"="it",
        "nodes"=c(),
        "def_pos"=NULL
      ),
      "rec"=list(
        "type"="process",
        "title"="Recommended Therapies",
        "bundle"=NULL,
        "id"="rec",
        "nodes"=c(),
        "def_pos"=NULL
      )
    )
  )
  
  link <- list(      
    "color"="rgb(204, 235, 197)",
    "source"=NULL,
    "value"=10,
    "type"=NULL,
    "target"=NULL,
    "opacity"=1,
    "time"="*",
    "title"=NULL,
    "id"=NULL
  )
  
  node <- list(      
    "bundle"=NULL,
    "title"=NULL,
    "visibility"="visible",
    "def_pos"=NULL,
    "id"=NULL,
    "style"="process",
    "direction"="r"
  )
  
  node_count <- 1
  link_count <- 1
  
  for (i in rownames(input)) {
    row <- input[i, ]
    mp_node <- node
    mp_node$title <- row$Gene_protein
    mp_node$id <- paste0("mp^", row$Gene_protein)
    
    json$nodes[[node_count]] <- mp_node
    node_count <- node_count + 1
    
    json$groups$mp$nodes <- c(json$groups$mp$nodes, paste0("mp^", row$Gene_protein))
  }
  
  if(nrow(fda) >= 1 && !length(fda$Note) >= 1 ) {
    
    unique_fda <- aggregate(fda, by=list(fda$Drug, fda$`Gene or Protein`), FUN=paste)
    
    
    for (i in rownames(unique_fda)) {
      row <- unique_fda[i, ]
      
      fda_node <- node
      fda_node$title <- row$Drug
      fda_node$id <- paste0("fda^", row$Drug)
      
      json$nodes[[node_count]] <- fda_node
      node_count <- node_count + 1
      
      json$groups$fda$nodes <- c(json$groups$fda$nodes, paste0("fda^", row$Drug))
      
      fda_link <- link
      fda_link$source <- paste0("mp^", row[["Gene or Protein"]]) 
      fda_link$target <- paste0("fda^", row$Drug)
      fda_link$type <- "FDA approved Drugs"
      fda_link$title <- paste0("alteration: ", row$Alteration)
      fda_link$id <- paste0(fda_link$source, " -> ", fda_link$target)
      
      json$links[[link_count]] <- fda_link
      link_count <- link_count + 1
    }
  }
  
  if(nrow(fda_other) >= 1) {
    
    unique_fda_other <- aggregate(fda_other, by=list(fda_other$Drug, fda_other$`Gene or Protein`), FUN=paste)
    
    
    for (i in rownames(unique_fda_other)) {
      row <- unique_fda_other[i, ]
      
      fda_other_node <- node
      fda_other_node$title <- row[["Group.1"]]
      fda_other_node$id <- paste0("fda^", row[["Group.1"]], " (other tumor type)")
      
      json$nodes[[node_count]] <- fda_other_node
      node_count <- node_count + 1
      
      json$groups$fda$nodes <- c(json$groups$fda$nodes, 
                                 paste0("fda^", row[["Group.1"]], " (other tumor type)"))
      
      fda_link <- link
      fda_link$source <- paste0("mp^", row[["Group.2"]]) 
      fda_link$target <- paste0("fda^", row[["Group.1"]], " (other tumor type)")
      fda_link$type <- "FDA approved Drugs"
      fda_link$title <- paste0("alteration: ", row$Alteration, "\n <br> tumor drug is approved for: ", row[["Tumor in which it is approved"]])
      fda_link$id <- paste0(fda_link$source, " -> ", fda_link$target)
      
      json$links[[link_count]] <- fda_link
      link_count <- link_count + 1
    }
  }
  
  if(nrow(rec) >= 1) {
    unique_rec <- aggregate(rec, by=list(rec$Drug, rec$`Gene or Protein`), FUN=paste)
    
    for (i in rownames(unique_rec)) {
      row <- unique_rec[i, ]
      
      # rec_node <- node
      # rec_node$title <- row$Drug
      # rec_node$id <- paste0("rec^", row$Drug)
      
      # json$nodes[[node_count]] <- rec_node
      # node_count <- node_count + 1
      
      # json$groups$rec$nodes <- c(json$groups$rec$nodes, paste0("rec^", row$Drug))
      
      fda_link <- link
      fda_link$source <- paste0("it^", row[["Group.2"]]) 
      fda_link$target <- paste0("rec^", row[["Group.1"]])
      fda_link$type <- "Recommended therapies"
      fda_link$title <- paste0("alteration: ", row$Alteration, "\n <br> tumor drug is approved for: ", row[["Tumor in which it is approved"]],
                               "\n <br> pathway: ", row$Path)
      fda_link$id <- paste0(fda_link$source, " -> ", fda_link$target)
      
      json$links[[link_count]] <- fda_link
      link_count <- link_count + 1
    }
    
    for (i in unique(rec[["Drug"]])) {
      rec_node <- node
      rec_node$title <- i
      rec_node$id <- paste0("rec^", i)
      
      json$nodes[[node_count]] <- rec_node
      node_count <- node_count + 1
      
      json$groups$rec$nodes <- c(json$groups$rec$nodes, paste0("rec^", i))
    }
    
    for (i in unique(rec[["Gene or Protein"]])) {
      rec_node <- node
      rec_node$title <- i
      rec_node$id <- paste0("it^", i)
      
      json$nodes[[node_count]] <- rec_node
      node_count <- node_count + 1
      
      json$groups$it$nodes <- c(json$groups$it$nodes, paste0("it^", i))
    }
    
    for (i in unique(rec[["Path"]])) {
      
      paths <- strsplit(i, "->")
      
      rec_link <- link
      rec_link$source <- paste0("mp^", trimws(paths[[1]][1]))
      rec_link$target <- paste0("it^", trimws(paths[[1]][length(paths[[1]])]))
      rec_link$type <- "Pathway"
      rec_link$title <- i
      rec_link$id <- i
      
      json$links[[link_count]] <- rec_link
      link_count <- link_count + 1
    }
  }
  
  gmp <- json$groups$mp
  gfda <- json$groups$fda
  git <- json$groups$it
  grec <- json$groups$rec
  
  # order <- list(list(gmp$nodes, gfda$nodes), list(list(git$nodes)), list(list(grec$nodes)))
  # json$order <- order
  
  groups <- list(gmp, gfda, git, grec)
  json$groups <- groups
  json
}

##get Category 1,2 recommendations
##if cat2 == "yes", then it's category 2, otherwise it's category 1
##(difference is in using given cancer type, vs. all other cancer types)
get_cat_1_2 <- function(MP,
                        cancer_type,
                        drugs_PO_FDA_biomarkers,
                        cat2 = "no")
{
  ##change cat2 to lower case, in case someone typed "Yes," "YES" etc
  cat2 <- tolower(cat2)
  if(cat2 == "y" | cat2 == "ye")
  {
    cat2 <- "yes"
  }
  
  Type1_2_df <- NULL
  if(cat2 != "yes")
  {
    Type1_2_df <- dplyr::filter(drugs_PO_FDA_biomarkers,
                                Gene.Protein %in% MP$Gene_protein,
                                Disease == cancer_type)
  } else {
    Type1_2_df <- dplyr::filter(drugs_PO_FDA_biomarkers,
                                Gene.Protein %in% MP$Gene_protein,
                                Disease != cancer_type)
  }
  
  Type1_2_df <- check_WT_biomarkers(Type1_2_df, MP)
  
  Type1_2_df <- Type1_2_df[,c("Drug","Gene.Protein","Type","Alteration","Disease")]
  colnames(Type1_2_df)[2] <- "Gene or Protein"
  ##change "Disease" to "Tumor in which it is approved"
  colnames(Type1_2_df)[colnames(Type1_2_df) == "Disease"] <- 
    "Tumor in which it is approved"
  
  if(nrow(Type1_2_df) == 0)
  {
    Type1_2_df <- 
      data.frame(Note = "There are no recommended therapies in this category.")  
  }
  
  Type1_2_df
}
  
##get Category 3,4 recommendations
##if cat4 == "yes", then it's category 4, otherwise it's category 3
##(difference is in using given cancer type, vs. all other cancer types)
get_cat_3_4 <- function(MP, ##input data frame
                        cancer_type, ##character string
                        cat_drugs, ##character string - whether only targeted cancer therapies, all FDA approved therapies, or all drugs in DrugBank
                        list_paths, ##list of pathways
                        Onc_df, ##data frame with oncogenes and cancer type
                        drugs_PO_FDA_biomarkers, ##FDA-approved precision oncology drugs with listed biomarkers
                        drugs_PO_FDA_targets, ##FDA-approved precision oncology drugs with targets
                        drug_targets, ##data frame of targets for specific genes/proteins - right now, have DrugBank_targets
                        FDA_approved_drugs, ##list of FDA-approved drugs)
                        Type1, ##results of get_cat_1_2
                        Type2, ##results of get_cat_1_2
                        Type3 = NULL, ##results of get_cat_3_4 (only use in cat4 == "yes)
                        cat4 = "no",
                        subtype = "activation") ##choose subtype of connections, if have multiple ones in list_paths
{
  ##change cat4 to lower case, in case someone typed "Yes," "YES" etc
  cat4 <- tolower(cat4)
  if(cat4 == "y" | cat4 == "ye")
  {
    cat4 <- "yes"
  }
  
  cancer_path <- NULL
  if(cat4 != "yes")
  {
    ##get the pathway for this cancer
    ##get just the connections in that pathway - note that only keeping activations!
    cancer_path <- get_connections_paths(list_paths, cancer_type, subtype=subtype)
  } else {
    ##get the connections from all the other cancer pathways
    other_paths <- setdiff(unique(Onc_df$Name),
                           cancer_type)
    
    cancer_path <- get_connections_paths(list_paths, other_paths, subtype=subtype)
  }
  
  if(cat4 != "yes")
  {
    ##get the oncogenes from that pathway
    Onc_df <- dplyr::filter(Onc_df,
                            Name == cancer_type)
    ##make sure these oncogenes are not wild type in the MP 
    ##(which means they would not in fact be oncogenic)
    ##get molecular alterations that are not wild-type in the MP
    MP_not_WT <- remove_WT(MP)
  } else
  {
    ##for this category, consider all oncogenes from all cancer types
    ##make sure these oncogenes are not wild type in the MP 
    ##(which means they would not in fact be oncogenic)
    ##get molecular alterations that are not wild-type in the MP
    MP_not_WT <- remove_WT(MP)
  }
  
  inputs <- unique(as.character(MP_not_WT$Gene_protein))
  
  ##get genes among the MP that are both 1) oncogenes and 2) not wild type
  ##these will be the only inputs considered further
  ##(will look downstream of them)
  inputs_in_path <- dplyr::filter(Onc_df,
                                  Oncogene %in% MP_not_WT$Gene_protein)
  inputs_in_path <- unique(as.character(inputs_in_path$Oncogene)) 
  
  ##keep just the genes that are downstream from inputs + the inputs
  downstream_genes <- get_downstream_from_inputs(cancer_path,
                                                 MP,
                                                 inputs_in_path)
  keep_genes <- unique(c(downstream_genes, inputs))
  
  ##just keep connections to and from these genes + inputs
  cancer_path <- dplyr::filter(cancer_path,
                               (to %in% keep_genes) | 
                                 (from %in% keep_genes))

  ##get the drugs which either target genes/proteins in keep_genes or for which
  ##the genes/proteins in keep_genes are biomarkers
  drugs_PO_FDA_biomarkers_keep <- filter_drugs_select_genes(drugs_PO_FDA_biomarkers,
                                                            keep_genes)
  drugs_PO_FDA_targets_keep <- filter_drugs_select_genes(drugs_PO_FDA_targets,
                                                         keep_genes)
  drugs_targets_keep <- filter_drugs_select_genes(drug_targets,
                                                  keep_genes)
  
  ##now get the drugs based on the chosen category
  ##by default, include everything
  drugs_biomarkers_targets <- 
    rbind(drugs_PO_FDA_biomarkers_keep[,1:6],
          drugs_PO_FDA_targets_keep[,1:6],
          drugs_targets_keep)
  if(cat_drugs == "Only FDA-approved targeted therapies for cancer")
  {
    drugs_biomarkers_targets <- 
      rbind(drugs_PO_FDA_biomarkers_keep[,1:6],
            drugs_PO_FDA_targets_keep[,1:6])
  }
  if(cat_drugs == "Only FDA-approved therapies")
  {
    drugs_biomarkers_targets <- 
      drugs_biomarkers_targets %>% filter(Drug %in% FDA_approved_drugs)
  }
  
  ##remove rows with wild-type biomarkers if those biomarkers are among the (non-wild-type) inputs
  drugs_biomarkers_targets <-
    drugs_biomarkers_targets %>% 
    filter(!(tolower(Alteration) %in% c("wild type","wildtype","wild","wt") &
               Gene.Protein %in% inputs))
  
  ##add shortest path from inputs to targets/biomarkers in drugs_biomarkers_targets
  drugs_biomarkers_targets <- add_shortest_path(drugs_biomarkers_targets, 
                                                inputs, inputs_in_path,
                                                cancer_path)
  
  ##simplify/remove some column names for display purposes
  drugs_biomarkers_targets <- pretty_drugs_mat(drugs_biomarkers_targets)
  
  ##create data frame of edges by combining pathway and gene-drug info
  edges_df <- data.frame(Source = c(as.character(cancer_path$from), as.character(drugs_biomarkers_targets$Drug)),
                         Target = c(as.character(cancer_path$to), as.character(drugs_biomarkers_targets$Gene.Protein)))
  
  ##change "Disease" to "Tumor in which it is approved"
  colnames(drugs_biomarkers_targets)[colnames(drugs_biomarkers_targets) == "Disease"] <- "Tumor in which it is approved"
  
  ##take out all genes/proteins already present above
  if(cat4 != "yes")
  {
    drugs_biomarkers_targets <- remove_known_recs(drugs_biomarkers_targets,
                                                  Type1_df=Type1, 
                                                  Type2_df=Type2)
  } else {
    drugs_biomarkers_targets <- remove_known_recs(drugs_biomarkers_targets,
                                                  Type1_df=Type1,
                                                  Type2_df=Type2,
                                                  Type3_df=Type3)
                                                  
  }
  
  colnames(drugs_biomarkers_targets)[c(2,3)] <- c("Gene or Protein","Type")
  
  list(drugs_mat = drugs_biomarkers_targets,
       edges = edges_df)
}
