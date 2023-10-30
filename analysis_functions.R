# Making plots from graphviz files
library(readr)
library(data.table)
library(dplyr)
library(limma)
library(tidyr)
library(tools)
library(igraph)
library(EMT)
library(cowplot)
library(tidyverse)
library(caret)
library(combinat, exclude = c("combn")) 
library(randomForest)
library(ggpubr)
library(networkD3)
library(readxl)
library(reshape)
library(forcats)
library(dplyr)
library(ggplot2)
library(stringr)
library(survival)

mergevars = function(vec) {
  if(1 %in% vec) return(1)
  else if(2 %in% vec & 0 %in% vec) return(1)
  else if(2 %in% vec) return(2)
  else if(0 %in% vec) return(0)
  else return(3)
}

# Function to fix the Stanford sample names
fix_stanford_samples = function(sampvec) {
  sampvec = gsub("SU313_C-R","SU313C-R",sampvec)
  sampvec = gsub("SU320A","SU320",sampvec)
  sampvec = as.character(sapply(sampvec, function(x) 
    if(!startsWith(x,"SU")) paste0("SU",x) else x))
  sampvec = as.character(sapply(sampvec, function(x) 
    strsplit(x, "_|-")[[1]][1]))
  return(sampvec)
}

# Function to fix the MSK sample names
fix_msk_samples = function(sampvec) {
  sampvec = gsub("*","",gsub("^","",sampvec,fixed=TRUE),fixed=TRUE)
  return(sampvec)
}

# Function to get the reported doublet rate
get_doublet_rate = function(output) {
  outputfile = read_lines(output)
  doubletchars = gsub("best relevant doublet rate:\t","",
                      outputfile[grepl("best relevant doublet rate",outputfile)],fixed=TRUE)
  return(as.numeric(doubletchars))
}

# Function to get variant from a line in a gv file
getgene = function(string) {
  var = strsplit(string,split='\"',fixed=TRUE)[[1]][2]
  gene = strsplit(var,":p.",fixed=TRUE)[[1]][1]
  # Account of times when we have a splicing variant
  if(grepl("splicing",gene)) gene = strsplit(gene,"_")[[1]][1]
  # And account for times when we have FLT3-ITD
  if(grepl("FLT3-ITD",gene)) gene = strsplit(gene, "-")[[1]][1]
  # Can return gene or variant alone if preferred
  return(c(gene,var))
}

# Function to get link between two IDs
getlink = function(string) {
  link = paste0("v",strsplit(string, " -> ")[[1]])
  return(gsub(";","",link))
}

# Function to get a matrix with the gene and its links from gv file
parsegv = function(gvfile) {
  treefile = read_lines(gvfile)
  labs = treefile[grepl("label",treefile)]
  edges = treefile[grepl("->",treefile)]
  genes = data.frame(t(sapply(labs, getgene)),1)
  colnames(genes) = c("Gene","Variant","Count")
  genes = suppressMessages(genes %>% left_join(mutclass))
  # IDs for genes
  ID = gsub("[","",paste0("v",sapply(labs, substr, start = 1, stop = 2)),fixed=TRUE)
  genes = data.frame(ID, genes)
  links = data.frame(t(sapply(edges, getlink)), 1)
  colnames(links) = c("from","to","weight")
  return(list(genes = genes, links = links))
}

# Function to get the list of paths from a tree that involve the root
get_paths = function(graphtmp) {
  l <- unlist(lapply(V(graphtmp), function(x) all_simple_paths(graphtmp, from=x)), recursive = F)
  allpathroot = lapply(1:length(l), function(x) as_ids(l[[x]]))
  rootname = vertex_attr(graphtmp)$name[vertex_attr(graphtmp)$Gene == "Root"]
  rootpres = sapply(allpathroot, function(x) rootname %in% x)
  return(allpathroot[rootpres])
}

# Function to identify duplicate paths from a vector of paths
get_dup_path = function(pathtmp, pathvec) {
  return(any(grepl(pathtmp, pathvec[pathvec != pathtmp])))
}

# Function to get the weight of a specific link
get_link_weights = function(vec, g, totpaths) {
  # Check if specific vector is present in all paths
  link_present = any(sapply(lapply(totpaths, all.equal, 
                                   current = vec), paste, collapse=",") == TRUE)
  # If a link is present, return its weight from the graph
  if(link_present) return(E(g)$weight[E(g, P = vec)])
  # If the link isn't present, get the path that have this link
  grepfrom = sapply(lapply(totpaths, function(x) x == vec[1]),any)
  grepto = sapply(lapply(totpaths, function(x) x == vec[2]),any)
  # Select the first path (this is arbitrary since doing this with the v# names)
  efull = totpaths[grepfrom & grepto][[1]]
  edgepath = efull[grep(TRUE, efull == vec[1]):grep(TRUE, efull == vec[2])]
  # Now get the weights within that path
  allweights = NULL
  for(i in 1:(length(edgepath)-1)) {
    allweights[i] = E(g)$weight[E(g, P = edgepath[c(i,i+1)])]
  }
  # return the minimum weight in that path
  return(min(allweights))
}

# Function to get all links within each path and return data frame of links
get_all_directed_links = function(g) {
  allpaths = get_paths(g)
  
  # Summarize each path as a character string and then limit to those not duplicated
  allpathstring = sapply(allpaths, function(x) paste(x, collapse = "->"))
  allpaths = allpaths[!sapply(allpathstring, get_dup_path, pathvec = allpathstring)]
  
  # Need to collapse the same paths within a single graph
  allcomb = do.call(rbind, lapply(allpaths, function(x) t(combn(x, m=2))))
  
  # Get each path weights
  weights = apply(allcomb, 1, get_link_weights, g = g, totpaths = allpaths)
  
  return(data.frame(from = allcomb[,1], to = allcomb[,2], 
                    weight = weights) %>% distinct())
}

# Function to get all links within each path and return data frame of links, without weights
get_all_directed_links_no_weight = function(g) {
  allpaths = get_paths(g)
  
  # Summarize each path as a character string and then limit to those not duplicated
  allpathstring = sapply(allpaths, function(x) paste(x, collapse = "->"))
  allpaths = allpaths[!sapply(allpathstring, get_dup_path, pathvec = allpathstring)]
  
  # Need to collapse the same paths within a single graph
  allcomb = do.call(rbind, lapply(allpaths, function(x) t(combn(x, m=2))))
  
  return(data.frame(from = allcomb[,1], to = allcomb[,2]) %>% distinct())
}

# Function to count number of FLT3-ITD variants in a file
count_itd = function(gvfile) {
  parsed = parsegv(gvfile)
  itd_count = sum(grepl("FLT3-ITD",parsed$genes$Variant))
  return(itd_count)
}

# function to get sample name from a file name
get_sample_name = function(gvfile) {
  if(!grepl("MD_A",gvfile)) sampname = paste(strsplit(gvfile, "_")[[1]][1:2], collapse = "_")
  else sampname = paste(strsplit(gvfile, "_")[[1]][1:3], collapse = "_")
  return(sampname)
}

# Make a data frame out of a vector for joining
make_df = function(namedvec) {
  return(data.frame(t(as.matrix(data.frame(namedvec)))))
}

# Function to remove a node and reconnect the graph
remove_node = function(g, node) {
  gdf = as_long_data_frame(g)
  tonodes = unique(gdf$to_name[gdf$from_name == node])
  fromnode = unique(gdf$from_name[gdf$to_name == node])
  # If this node goes to another node, then address that
  if(length(tonodes) > 0) {
    edgevec = as.vector(rbind(fromnode,tonodes))
    weights = E(g)[node %--% tonodes]$weight
    #print(weights)
    g_new = delete_edges(g, E(g)[node %--% tonodes]) %>%
      delete_vertices(node) %>%
      add_edges(edgevec)
    newedges = E(g_new)[fromnode %--% tonodes]
    for(i in 1:length(newedges)) {
      g_new = set_edge_attr(g_new, "weight",index = newedges[i], weights[i])
    }
  }
  else g_new = g %>% delete_vertices(node)
  return(g_new)
}

# Function to remove non-pathogenic variants
# Using the "driver_vars_tab" object below for a given tree
remove_nonpath = function(g, gvfile) {
  # get the non-pathogenic variants
  drivervars = driver_vars_tab$AAchange[driver_vars_tab$Tree == gvfile]
  pathlims = V(g)$Variant %in% drivervars | V(g)$Variant == "Root"
  nonpath = V(g)$name[!pathlims]
  nonpathvars = V(g)$Variant[!pathlims]
  if(length(nonpath)>0) {
    for(i in 1:length(nonpath)) g = remove_node(g, nonpath[i])
  }
  return(g)
}

# Function to get the parents of the FLT3-ITD variants
get_founder_itd = function(gvfile) {
  parsedgv = parsegv(gvfile)
  g = graph_from_data_frame(parsedgv$links, directed=TRUE, 
                            vertices=parsedgv$genes)
  itdvars = V(g)[grepl("FLT3-ITD",V(g)$Variant)]
  itdparents = NULL
  itdchildren = vector(mode = "list")
  for(i in 1:length(itdvars)) {
    itdparents[i] = V(g)$Variant[neighbors(g, itdvars[i], mode = "in")]
    itdchildren[[i]] = V(g)$Variant[neighbors(g, itdvars[i], mode = "out")]
  }
  # Get boolean for whether there are empty children (other than other ITD variants)
  itdchildren = sapply(itdchildren, function(x) x[!grepl("FLT3-ITD",x)])
  emptychildren = sapply(itdchildren, function(x) length(x) == 0)
  return(any(itdparents == "Root" & !emptychildren))
}

# Using this for the FLT3-ITD stuff below
make_one_graph_nonpath = function(gfile, analysisvar, varsize) {
  parsedgv = parsegv(gfile)
  g = graph_from_data_frame(parsedgv$links, 
                            directed=TRUE, vertices=parsedgv$genes)
  if(analysisvar == "Gene") {
    png(paste0(gfile,".png"), width = 1000, height = 1000)
    plot(g, vertex.color = genecolor[V(g)$Gene],
         vertex.label = V(g)$Variant, layout=layout_(g, as_tree()),
         vertex.size = varsize, vertex.label.cex = 2.5,
         vertex.label.family = "Helvetica",
         vertex.label.font = 2)
    dev.off()
  }
  else if(analysisvar == "Pathway") {
    png(paste0(gfile,".png"), width = 1000, height = 1000)
    plot(g, vertex.color = pathwaycolor[V(g)$Pathway],
         vertex.label = V(g)$Gene, layout=layout_(g, as_tree()),
         vertex.size = varsize, vertex.label.cex = 2.5,
         vertex.label.family = "Helvetica",
         vertex.label.font = 2)
    dev.off()
  }
  else stop("analysisvar should be Gene or Pathway")
}

# Function to make basic graph from file
make_one_graph = function(gfile, analysisvar, varsize) {
  parsedgv = parsegv(gfile)
  g = remove_nonpath(graph_from_data_frame(parsedgv$links, 
                                           directed=TRUE, vertices=parsedgv$genes), gvfile = gfile)
  if(analysisvar == "Gene") {
    plot(g, vertex.color = genecolor[V(g)$Gene],
         vertex.label = V(g)$Variant, layout=layout_(g, as_tree()),
         vertex.size = varsize)
  }
  else if(analysisvar == "Pathway") {
    plot(g, vertex.color = pathwaycolor[V(g)$Pathway],
         vertex.label = V(g)$Gene, layout=layout_(g, as_tree()),
         vertex.size = varsize)
  }
  else stop("analysisvar should be Gene or Pathway")
}


# function to get a data.frame with the number of ITD variants for all ITD merging patterns
# Using the gv files in a directory
get_itd_count_df = function(dirtmp) {
  allfiles = list.files(dirtmp)[endsWith(list.files(dirtmp),"map0.gv")]
  libtrees = allfiles[grepl("liberal", allfiles)]
  nonetrees = allfiles[grepl("none", allfiles)]
  contrees = allfiles[grepl("conservative", allfiles)]
  alltrees = allfiles[grepl("all", allfiles)]
  libcounts = sapply(libtrees, count_itd)
  names(libcounts) = sapply(libtrees, get_sample_name)
  nonecounts = sapply(nonetrees, count_itd)
  names(nonecounts) = sapply(nonetrees, get_sample_name)
  concounts = sapply(contrees, count_itd)
  names(concounts) = sapply(contrees, get_sample_name)
  allcounts = sapply(alltrees, count_itd)
  names(allcounts) = sapply(alltrees, get_sample_name)
  df_all = full_join(make_df(libcounts), make_df(concounts)) %>% 
    full_join(make_df(allcounts)) %>% full_join(make_df(nonecounts)) %>%
    replace(is.na(.), 1)
  df_all = data.frame(merge_type = 
                        c("liberal","conservative","all","none"),df_all)
  return(df_all)
}

# Function to get the parents of the FLT3-ITD variants
get_itd_parents = function(gvfile) {
  parsedgv = parsegv(gvfile)
  g = graph_from_data_frame(parsedgv$links, directed=TRUE, 
                            vertices=parsedgv$genes)
  itdvars = V(g)[grepl("FLT3-ITD",V(g)$Variant)]
  itdparents = NULL
  for(i in 1:length(itdvars)) {
    itdparents[i] = V(g)$Variant[neighbors(g, itdvars[i], mode = "in")]
  }
  return(unique(itdparents))
}

# Function to get the parents of the FLT3-ITD variants
get_itd_children = function(gvfile) {
  parsedgv = parsegv(gvfile)
  g = graph_from_data_frame(parsedgv$links, directed=TRUE, 
                            vertices=parsedgv$genes)
  itdvars = V(g)[grepl("FLT3-ITD",V(g)$Variant)]
  itdchildren = NULL
  for(i in 1:length(itdvars)) {
    kids = V(g)$Variant[neighbors(g, itdvars[i], mode = "out")]
    itdchildren[i] = paste(kids[!grepl("FLT3-ITD",kids)], collapse = ", ")
  }
  return(na.omit(itdchildren[itdchildren != ""]))
}

# Function to get support for a specific pair of variants
get_order_support = function(gfile, mutorder) {
  mutmat = mutmatlist[[gfile]]
  mut1 = mutmat[mutmat$AAchange == mutorder$from,-c(1:2)]
  mut2 = mutmat[mutmat$AAchange == mutorder$to,-c(1:2)]
  if(nrow(mut1) > 1) mut1 = apply(mut1, 2, mergevars)
  if(nrow(mut2) > 1) mut2 = apply(mut2, 2, mergevars)
  mut1 = as.numeric(mut1)
  mut2 = as.numeric(mut2)
  lims = mut1 != 3 & mut2 != 3
  numerator = sum(mut1[lims] != 0 & mut2[lims] != 0)
  denominator = sum(mut2[lims] != 0)
  return(c(numerator/denominator, numerator))
}

# Function to get support vector for an entire tree from a gv file
get_tree_support = function(gfile) {
  parsed = parsegv(gfile)
  g = graph_from_data_frame(parsed$links, directed=TRUE, 
                            vertices=parsed$genes)
  if(length(V(g)) > 1) {
    allpaths = get_all_directed_links_no_weight(g)
    # Drop paths with the root node
    rootnode = V(g)$name[V(g)$Gene == "Root"]
    allpaths = allpaths %>% filter(!apply(allpaths, 1, function(x) any(x == rootnode)))
    if(nrow(allpaths) > 0) {
      allpathvec = rep(NA, nrow(allpaths))
      cells = rep(NA, nrow(allpaths))
      fraction = rep(NA, nrow(allpaths))
      for(i in 1:nrow(allpaths)) {
        mutorder = allpaths[i,1:2]
        mutorder[1] = V(g)$Variant[V(g)$name == allpaths[i,1]]
        mutorder[2] = V(g)$Variant[V(g)$name == allpaths[i,2]]
        allpathvec[i] = paste(mutorder,collapse="->")
        pathsupport = get_order_support(gfile, mutorder)
        fraction[i] = pathsupport[1]
        cells[i] = pathsupport[2]
      }
      return(data.frame(Path = allpathvec, Fraction = fraction, Cells = cells))
    }
    else {
      return("No paths")
    }
  }
  else return("No paths")
}

# Function to get support vector for an entire tree from a gv file
# Main difference is that we haven't removed VUS's and including root node
get_itd_graph_support = function(gfile) {
  parsed = parsegv(gfile)
  g = graph_from_data_frame(parsed$links, directed=TRUE, 
                            vertices=parsed$genes)
  if(length(V(g)) > 1) {
    allpaths = get_all_directed_links(g)
    if(nrow(allpaths) > 0) {
      allpathvec = rep(NA, nrow(allpaths))
      cells = rep(NA, nrow(allpaths))
      fraction = rep(NA, nrow(allpaths))
      for(i in 1:nrow(allpaths)) {
        mutorder = allpaths[i,1:2]
        mutorder[1] = V(g)$Variant[V(g)$name == allpaths[i,1]]
        mutorder[2] = V(g)$Variant[V(g)$name == allpaths[i,2]]
        allpathvec[i] = paste(mutorder,collapse="->")
        pathsupport = get_order_support(gfile, mutorder)
        fraction[i] = pathsupport[1]
        cells[i] = pathsupport[2]
      }
      return(data.frame(Path = allpathvec, Fraction = fraction, Cells = cells))
    }
    else {
      return("No paths")
    }
  }
  else return("No paths")
}

# Function to identify the tree of choice for itds
get_optimal_itd_tree = function(gvfiles) {
  # First, drop the tree where we don't do any merging since we won't be using that
  gvfiles = gvfiles[!grepl("none",gvfiles)]
  parsedgvs = lapply(gvfiles, parsegv)
  gs = lapply(parsedgvs, function(x) 
    graph_from_data_frame(x$links, directed=TRUE, 
                          vertices=x$genes))
  n_itds = sapply(gs, function(x) 
    sum(grepl("FLT3-ITD",V(x)$Variant)))
  # If there's only one ITD, then we'll just "merge all"
  if(all(n_itds == 1)) return(c(gvfiles[grepl("all",gvfiles)], "one ITD"))
  # Check if there is only one parent across all graphs, not including the FLT3-ITD parents
  parents = unique(unlist(lapply(gvfiles, get_itd_parents)))
  parents_same = length(parents[!grepl("FLT3-ITD",parents)]) <= 1
  # Check if there are no children of the FLT3-ITD variants across all graphs
  no_children = length(unlist(lapply(gvfiles, get_itd_children))) == 0
  # Only do all these calculations if we already know that the parents 
  # are the same and there are no child notes
  if(parents_same & no_children) {
    all_links = lapply(gs, get_all_directed_links)
    # Change links to variant names
    for(i in 1:length(all_links)) {
      node_tmp = names(V(gs[[i]]))
      for(j in 1:length(node_tmp)) {
        all_links[[i]][all_links[[i]] == node_tmp[j]] = 
          V(gs[[i]])$Variant[node_tmp == node_tmp[j]]
      }
    }
    # Get the links without ITDs and paste them into individual strings
    no_itd_links = lapply(lapply(all_links, function(x)
      x[!grepl("FLT3-ITD",x[,1]) & !grepl("FLT3-ITD",x[,2]),]),
      function(x) paste(x[,1], x[,2]))
    # Get whether any of the links are different
    other_links_same = length(unlist(lapply(1:length(no_itd_links), 
                                            function(n) setdiff(no_itd_links[[n]], 
                                                                unlist(no_itd_links[-n]))))) == 0
  }
  else other_links_same = FALSE
  # If all the conditions are met, return a merged_all graph
  if(parents_same & no_children & other_links_same) {
    return(c(gvfiles[grepl("all",gvfiles)],"same parents, no kids, others same"))
  }
  # Otherwise, get the number of low support variants and then return
  # the graph with the most conservative merging strategy of the graphs
  # that have the minimum number of low support variants
  tree_support = sapply(lapply(gvfiles, get_itd_graph_support),
                        function(x) sum(is.nan(x$Fraction)) + sum(na.omit(x$Fraction) < 0.5))
  good_trees = gvfiles[tree_support == min(tree_support)]
  if(any(grepl("conserv",good_trees))) 
    return(c(good_trees[grepl("conserv",good_trees)],"best support"))
  else if(any(grepl("lib",good_trees))) 
    return(c(good_trees[grepl("lib",good_trees)],"best support"))
  else return(c(good_trees[grepl("all",good_trees)],"best support"))
}

# Function to remove low support variants from the driver_vars_tab
remove_low_support = function(tree_tmp, driver_vars_tmp) {
  # Check if there are low-support paths
  badvars = NULL
  tree_support = get_tree_support(tree_tmp)
  low_support_paths = NULL
  if(is.data.frame(tree_support)) {
    low_support_paths = tree_support[tree_support$Fraction < 0.5 | 
                                       is.na(tree_support$Fraction),]
    # Only interested in low support paths that also involve drivers
    allvars = unique(unlist(sapply(low_support_paths$Path, strsplit, split = "->")))
    nonpathvars = setdiff(allvars, driver_vars_tmp$AAchange[driver_vars_tmp$Tree == tree_tmp])
    if(length(nonpathvars) > 0) {
      low_support_paths = low_support_paths[!grepl(paste(c(paste0(nonpathvars, "->"), 
                                                           paste0("->", nonpathvars)), collapse = "|"), low_support_paths$Path),]
    }
    low_support_paths_orig = low_support_paths
    if(nrow(low_support_paths) > 0) {
      # Make graph
      parsedgv = parsegv(tree_tmp)
      g = graph_from_data_frame(parsedgv$links, directed=TRUE, 
                                vertices=parsedgv$genes)
      # Get distance from root for each node
      rootnode = V(g)$name[V(g)$Variant == "Root"]
      rootdistance = distances(g, rootnode)[1,]
      # rename the elements of this vector to their variants
      for(i in 1:length(rootdistance)) {
        names(rootdistance)[i] = V(g)$Variant[V(g)$name == names(rootdistance)[i]]
      }
      # Now keep removing variants until the low_support_paths df is empty
      while(nrow(low_support_paths) > 0) {
        # Get all variants involved in low-support connections
        varlist = sapply(low_support_paths$Path,strsplit,split="->")
        # Get the variants involved in the greatest number of bad support connections
        vartab = table(unlist(varlist))
        badvars_tmp = names(vartab)[vartab == max(vartab)]
        # If multiple tied, limit to those that have the greatest distance from root
        if(length(badvars_tmp) == 1) badvar = badvars_tmp
        else {
          distant_lims = rootdistance[badvars_tmp] == max(rootdistance[badvars_tmp])
          badvar = names(rootdistance[badvars_tmp][distant_lims])
        }
        # If still multiple, arbitrarily choose one for this iteration
        if(length(badvar) > 1) badvar = sample(badvar, 1)
        # Remove all occasions when bad variant occurs start or end of a path
        low_support_paths = low_support_paths %>% filter(!(grepl(paste0(badvar, "->"),Path,fixed=TRUE) |
                                                             grepl(paste0("->",badvar),Path,fixed=TRUE)))
        badvars = c(badvars, badvar)
      }
      # If there are multiple variants with the same amino acid change, figure out which should be dropped
      if(any(table(driver_vars_tmp %>% filter(Tree == tree_tmp & AAchange %in% 
                                              badvars) %>% select(Tree, varvec, AAchange) %>% 
                   distinct() %>% pull(AAchange)) > 1)) {
        print(paste0(tree_tmp, " has multiple of same variant"))
        # Hard code a fix for the one time this is true
        badvars = badvars[badvars != "WT1:p.V371fs"]
        driver_vars_tmp = driver_vars_tmp %>% filter(!(Tree == tree_tmp & AAchange %in% badvars))
        driver_vars_tmp = driver_vars_tmp %>% filter(!(Tree == tree_tmp &
                                                         varvec == "WT1_32417941_C_CA"))
        print("Fixed this one case")
      }
      else {
        driver_vars_tmp = driver_vars_tmp %>% filter(!(Tree == tree_tmp & AAchange %in% badvars))
      }
    }
    low_support_paths = low_support_paths_orig
  }
  return(list(driver_vars_tmp, badvars, low_support_paths))
}


# Function to replace vertex names with name of analysis variant
replace_vertex_names = function(graphtmp, strings, analysisvar) {
  namevec = vertex_attr(graphtmp)$name
  if(analysisvar == "Gene") varvec = vertex_attr(graphtmp)$Gene
  else if(analysisvar == "Pathway") varvec = vertex_attr(graphtmp)$Pathway
  else stop("analysisvar should be either 'Gene' or 'Pathway'")
  names(varvec) = namevec
  
  # Replace nodes
  stringsnew=strings
  for(i in 1:length(stringsnew)) {
    splitvec = strsplit(stringsnew[i],"->")[[1]]
    stringsnew[i] = paste(varvec[splitvec],collapse="->")
  }
  
  names(stringsnew) = strings
  # This returns a dictionary of paths with the names being the v0->v1 etc
  # and the values being the equivalent path but 
  # using the biological unit of interest (e.g. TET2->TP53)
  return(stringsnew)
}


# Function to get last vertex of a path
get_last_vertex = function(pathvec) {
  return(tail(strsplit(pathvec,"->")[[1]], n=1))
}

# Function to get alternative vertex name pased on a vector of paths
# Again, pathvec1 is the path vectors for the full graph
get_alt_vertex = function(vert, pathvec1, pathvec2) {
  pathvec2split = lapply(names(pathvec2), function(x) strsplit(x, "->")[[1]])
  pathvec2lims = sapply(pathvec2split, function(x) any(x == vert))
  # Select the overlapping path
  pathvec2tmp = pathvec2[pathvec2lims & pathvec2 %in% pathvec1][1]
  pathvec1tmp = pathvec1[pathvec1 == pathvec2tmp]
  pathvec1split = strsplit(names(pathvec1tmp),"->")[[1]]
  pathveclim2split = strsplit(names(pathvec2tmp), "->")[[1]]
  return(pathvec1split[which(pathveclim2split == vert)])
}

# Function to rename the vertices of graph2 
# so that they're compatible with graph 1 and return graph2.
# Note that pathvec1 is the path vectors for the full graph
rename_paths = function(pathvec1, pathvec2, graph1, graph2) {
  # Get the maximum index of the vertices currently in graph1
  vertnums = c(gsub("v","",V(graph1)$name), gsub("v","",V(graph2)$name))
  maxvertnum = max(as.numeric(vertnums))
  
  # Get the new paths and then the last vertex of each new path
  newverts = unique(sapply(names(pathvec2)[pathvec2 %in% 
                                             setdiff(pathvec2, pathvec1)], get_last_vertex))
  
  # Then get the vertices that aren't new, excluding the root
  rootvert2 = strsplit(names(pathvec2)[1],"->")[[1]][1]
  oldverts = setdiff(V(graph2)$name,newverts)
  oldverts = oldverts[oldverts != rootvert2]
  
  # Now replace the new vertex names with new names for the new graph
  if(length(newverts) > 0) {
    for(i in 1:length(newverts)) {
      V(graph2)$name[V(graph2)$name == newverts[i]] = paste0("v",i+maxvertnum)
    }
  }
  
  # Reassign the already existing vertices via a temp vertex name
  
  # First, get and reassign the root vertices
  rootvert1 = strsplit(names(pathvec1)[1],"->")[[1]][1]
  rootvert1tmp = gsub("v","t",rootvert1)
  V(graph2)$name[V(graph2)$name == rootvert2] = rootvert1tmp
  
  # Then if there are any old vertices that aren't the root, change their names to the temp name
  if(length(oldverts) > 0) {
    for(i in 1:length(oldverts)) {
      alttmp = gsub("v","t",get_alt_vertex(oldverts[i], pathvec1, pathvec2))
      V(graph2)$name[V(graph2)$name == oldverts[i]] = alttmp
    }
  }
  
  # Now rename the vertices so that we're not using the tmp names
  V(graph2)$name = gsub("t","v",V(graph2)$name)
  
  return(graph2)
}


# Function to convert a boolean of vertices into a numeric vector for the contract function
get_contract_num = function(boolean) {
  numvec = rep(NA, length(boolean))
  truepos = NA
  falsepos = 1
  for(j in 1:length(boolean)) {
    if(boolean[j]) {
      if(is.na(truepos)) {
        numvec[j] = j
        truepos = j
      }
      else numvec[j] = truepos
    }
    else {
      if(!is.na(truepos) & falsepos == truepos) {
        numvec[j] = falsepos + 1
        falsepos = falsepos + 2
      }
      else {
        numvec[j] = falsepos
        falsepos = falsepos + 1
      }
    }
  }
  return(numvec)
}

# Function to get duplicated vertices and collapse that part of the graph
collapse_duplicates = function(pathvec, g, analysisvar) {
  pathdup = table(pathvec)[table(pathvec) > 1]
  if(length(pathdup) == 0) return(g)
  for(i in 1:length(pathdup)) {
    # for a given duplicated path, get the duplicated vertex variants
    pathstringdups = names(pathvec)[pathvec == names(pathdup)[i]]
    dupverts = sapply(pathstringdups, get_last_vertex)
    #parentvert = tail(strsplit(pathstringdups[1],"->")[[1]], n=2)[1]
    
    # Change the attributes for the vertices and edges
    V(g)$Variant[V(g)$name %in% dupverts] = 
      paste(unique(V(g)$Variant[V(g)$name %in% dupverts]),collapse=",")
    if(analysisvar == "Pathway") {
      V(g)$Gene[V(g)$name %in% dupverts] = 
        paste(unique(V(g)$Gene[V(g)$name %in% dupverts]),collapse=",")	
    }
    # Subtracting 1 so that we are only adding the number of times duplicated
    V(g)$Count[V(g)$name %in% dupverts] = as.numeric(pathdup[i]) -1 + 
      V(g)$Count[V(g)$name %in% dupverts]
    #for(j in 1:length(dupverts)) {
    #	E(g)[parentvert %--% dupverts[j]]$weight = as.numeric(pathdup[i]) -1 + 
    #		E(g)[parentvert %--% dupverts[j]]$weight#
    #}
    
    # Make a vector we can use to contract the graph
    contract_num = get_contract_num(V(g)$name %in% dupverts)
    
    # Contract the graph, getting the first of the attributes
    # attributes were all set to the same above)
    g = contract(g, contract_num, vertex.attr.comb = "first")
  }
  return(g)
}

# Function to merge rows for string attributes
merge_rows = function(x, cnt) {
  if(!cnt) {
    x = gsub(" ","",x)
    x1 = NULL
    for(i in 1:length(x)) x1 = c(x1, strsplit(x[i], ",")[[1]])
    x1 = na.omit(unique(x))
    if(length(x1) == 1) return(x1)
    return(paste(unique(x1), collapse=","))
  }
  return(sum(as.numeric(na.omit(x))))
}

# Function to delete all vertex attributes in a vector
delete_all_attributes = function(g, vec) {
  for(i in 1:length(vec)) g = g %>% delete_vertex_attr(vec[i])
  return(g)
}

# Merge two directed keeping context for graph
# Drop the Variants attribute before doing this
merge_graphs = function(graph1, graph2, analysisvar) {
  # Get all the paths of big graph and graph to be added
  bigpaths = get_paths(graph1)
  newpaths = get_paths(graph2)
  
  # Summarize each path as a character string
  bigpathstring = sapply(bigpaths, function(x) paste(x, collapse = "->"))
  newpathstring = sapply(newpaths, function(x) paste(x, collapse = "->"))
  
  # Also replace the vertices by the analysis variable, i.e. gene or pathway
  pathvec1 = replace_vertex_names(graph1, bigpathstring, analysisvar)
  pathvec2 = replace_vertex_names(graph2, newpathstring, analysisvar)
  
  # Need to collapse the same paths within a single graph
  graph1 = collapse_duplicates(pathvec1, graph1, analysisvar)
  graph2 = collapse_duplicates(pathvec2, graph2, analysisvar)
  
  # function renames the vertices of the new graph to match (or add to) the bigger graph
  graph2renamed = rename_paths(pathvec1, pathvec2, graph1, graph2)
  
  # Now merge the two graphs
  graphnew = graph1 %u% graph2renamed
  
  # Merge the edge attributes
  weight = apply(data.frame(E(graphnew)$weight_1,E(graphnew)$weight_2),1,merge_rows,cnt=TRUE)
  graphnew = graphnew %>% set_edge_attr("weight",value=weight) %>%
    delete_edge_attr("weight_1") %>% delete_edge_attr("weight_2")
  
  # Merge the vertex attributes
  Gene = apply(data.frame(V(graphnew)$Gene_1, V(graphnew)$Gene_2),1,merge_rows,cnt=FALSE)
  Pathway = apply(data.frame(V(graphnew)$Pathway_1,V(graphnew)$Pathway_2),1,merge_rows,cnt=FALSE)
  Variant = apply(data.frame(V(graphnew)$Variant_1,V(graphnew)$Variant_2),1,merge_rows,cnt=FALSE)
  Count = apply(data.frame(V(graphnew)$Count_1,V(graphnew)$Count_2),1,merge_rows,cnt=TRUE)
  
  graphnew = graphnew %>% set_vertex_attr("Pathway",value=Pathway) %>%
    set_vertex_attr(name="Gene",value=Gene) %>% set_vertex_attr(name="Count",value=Count) %>%
    set_vertex_attr(name="Variant",value=Variant)
  # Delete the prior vertex attributes
  delattributes = c("Pathway","Count","Variant","Gene")
  delattributes = c(paste(delattributes,"1",sep="_"),paste(delattributes,"2",sep="_"))
  graphnew = delete_all_attributes(graphnew, delattributes)
  #print(V(graphnew)$Count)
  # Finally, return the merged graph
  return(graphnew)
}


# Function to get all graphs that begin with a certain feature
gene_graphs_feature = function(graphfiles, feature, analysisvar) {
  graphlist = vector(mode = "list", length = length(graphfiles))
  for(i in 1:length(graphlist)) {
    parsedgv = parsegv(graphfiles[i])
    graphlist[[i]] = remove_nonpath(graph_from_data_frame(parsedgv$links, directed=TRUE, 
                                                          vertices=parsedgv$genes), gvfile = graphfiles[i])
  }
  graphsublist = vector(mode = "list")
  j = 1
  for(i in 1:length(graphlist)) {
    g = graphlist[[i]]
    rootnode = V(g)$name[V(g)$Gene == "Root"]
    df = as_long_data_frame(g) %>% filter(from_name == rootnode)
    if(analysisvar == "Gene" & feature %in% df$to_Gene) {	
      graphsublist[[j]] = g
      names(graphsublist)[[j]] = graphfiles[i]
      j = j+1
    }
    if(analysisvar == "Pathway" & feature %in% df$to_Pathway) {
      graphsublist[[j]] = g
      names(graphsublist)[[j]] = graphfiles[i]
      j = j+1
    }
    if(!analysisvar %in% c("Gene","Pathway")) 
      stop("analysisvar should be Gene or Pathway")
  }
  return(graphsublist)
}

# Function to get a subgraph starting with a certain feature (pathway or gene)
get_subgraph_root = function(g, gdf, feature, analysisvar) {
  rootnode = V(g)$name[V(g)$Gene == "Root"]
  if(analysisvar == "Gene") featnode = unique(gdf$to_name[gdf$to_Gene == feature])
  else if(analysisvar == "Pathway") featnode = gdf$to_name[gdf$to_Pathway == feature]
  else stop("analysisvar should be Gene or Pathway")
  all_nodes = c()
  for(i in 1:length(featnode)) {
    g_sub = subcomponent(g, featnode[i], "out")
    all_nodes = c(all_nodes, as.numeric(g_sub))
  }
  subg = induced_subgraph(g, vids = unique(all_nodes))
  return(subg)
}

# Function to get a subgraph starting with a certain feature (pathway or gene)
get_subgraph_any = function(g, gdf, feature, analysisvar) {
  if(analysisvar == "Gene") featnode = unique(gdf$from[gdf$from_Gene == feature])
  else if(analysisvar == "Pathway") featnode = unique(gdf$from[gdf$from_Pathway == feature])
  else stop("analysisvar should be Gene or Pathway")
  all_nodes = c()
  for(i in 1:length(featnode)) {
    g_sub = subcomponent(g, featnode[i], "out")
    all_nodes = c(all_nodes, as.numeric(g_sub))
  }
  subg = induced_subgraph(g, vids = unique(all_nodes))
  if(length(E(subg)) > 1) {
    return(subg)
  }
}

# Make weighted and color plot
plot_weighted_color = function(g,vertmult,edgemult,arrowsize,arrowwidth,
                               legendloc,analysisvar,interact,trim,legendsize) {
  # make root a low weight
  V(g)$Count[V(g)$Gene=="Root"] = sort(V(g)$Count,decreasing=TRUE)[2]
  # Trim the terminal nodes if Count less than trim
  trimverts = V(g)$name[degree(g) == 1 & V(g)$Count <= trim]
  g = delete_vertices(g, trimverts)
  if(!interact) {
    if(analysisvar == "Pathway") {
      plot(g, vertex.color = pathwaycolor[V(g)$Pathway],
           vertex.label = "", vertex.size = vertmult*V(g)$Count,
           edge.width=edgemult*E(g)$weight,layout=layout_(g, as_tree()),
           edge.arrow.width=arrowwidth,edge.arrow.size=arrowsize)
      legend(legendloc, legend = intersect(names(pathwaycolor), unique(V(g)$Pathway)), 
             pch = 16, col = pathwaycolor[intersect(names(pathwaycolor), unique(V(g)$Pathway))], 
             bty = "n",cex=legendsize)
    }
    else if(analysisvar == "Gene") {
      plot(g, vertex.color = genecolor[V(g)$Gene],
           vertex.label = "", vertex.size = vertmult*V(g)$Count,
           edge.width=edgemult*E(g)$weight,layout=layout_(g, as_tree()),
           edge.arrow.width=arrowwidth,edge.arrow.size=arrowsize)
      legend(legendloc, legend = intersect(names(genecolor), unique(V(g)$Gene)), 
             pch = 16, 
             col = genecolor[intersect(names(genecolor), unique(V(g)$Gene))], 
             bty = "n",cex=legendsize)
    }
    else print("analysisvar should be 'Pathway' or 'Gene'")
  }
  if(interact) {
    if(analysisvar == "Pathway") {
      tkplot(g, vertex.color = pathwaycolor[V(g)$Pathway],
             vertex.label = "", vertex.size = vertmult*V(g)$Count,
             edge.width=edgemult*E(g)$weight,layout=layout_(g, as_tree()),
             edge.arrow.width=arrowwidth,edge.arrow.size=arrowsize)
      #legend(legendloc, legend = pathways, pch = 16, col = pathwaycolor, bty = "n",cex=legendsize)
    }
    else if(analysisvar == "Gene") {
      tkplot(g, vertex.color = genecolor[V(g)$Gene],
             vertex.label = "", vertex.size = vertmult*V(g)$Count,
             edge.width=edgemult*E(g)$weight,layout=layout_(g, as_tree()),
             edge.arrow.width=arrowwidth,edge.arrow.size=arrowsize)
      #legend(legendloc, legend = genes, pch = 16, col = genecolor, bty = "n",cex=legendsize)
    }
    else print("analysisvar should be 'Pathway' or 'Gene'")
  }
}

get_max_length = function(pathlist, node) {
  pathsub = pathlist[sapply(pathlist, function(x)
    any(x == node))]
  return(max(sapply(pathsub, length))-1)
}

get_node_pos = function(pathlist, node) {
  pathsub = pathlist[sapply(pathlist, function(x)
    any(x == node))]
  return(which(pathsub[[1]] == node)-1)
}

# Function to check for opposite path in a data frame
check_opposite_path = function(pathtmp, allpathtmp) {
  pathrev = rev(pathtmp)
  pathbool = apply(allpathtmp, 1, function(x) all(x == pathrev))
  return(any(pathbool))
}

# Function to get the weight of a given path, used in make_undirected_from_directed
get_path_weight = function(pathtmp, allpathtmp) {
  # limit to path vector without the weight column
  pathsub = pathtmp[1:2]
  pathrev = rev(pathsub)
  pathrevbool = apply(allpathtmp[,1:2], 1, function(x) all(x == pathrev))
  pathbool = apply(allpathtmp, 1, function(x) all(x == pathtmp))
  # return the sum of the weights for both paths
  return(sum(allpathtmp$weight[pathrevbool | pathbool]))
}

make_undirected_from_directed = function(g, analysisvar, only_bidirect) {
  # First get all the paths and convert to the relevant analysis variable (gene or pathway)
  allpaths = get_all_directed_links(g)
  if(analysisvar == "Gene") vardict = V(g)$Gene
  else if(analysisvar == "Pathway") vardict = V(g)$Pathway
  else stop("analysisvar should be either 'Gene' or 'Pathway'")
  names(vardict) = V(g)$name
  for(i in 1:length(vardict)) allpaths[allpaths == names(vardict)[i]] = vardict[i]
  
  # Eliminate the paths involving root
  # Adding the path weights together via an apply function
  if(!only_bidirect) {
    newdf = allpaths %>% filter(from != "Root") 
    newdf[,1:2] = t(apply(newdf[,1:2], 1, sort))
    newdf = newdf %>% group_by(from, to) %>% summarize(weight = sum(weight)) %>% ungroup()
    return(graph_from_data_frame(newdf, directed=FALSE))
  }
  
  # Get the subset of paths where the opposite path also exists in the dataset
  pathsrev = allpaths[apply(allpaths[,1:2], 1, check_opposite_path, 
                            allpathtmp = allpaths[,1:2]),]
  pathsrev$weight = apply(pathsrev, 1, get_path_weight, allpathtmp = pathsrev)
  pathsrev[,1:2] = t(apply(pathsrev[,1:2], 1, sort))
  
  return(graph_from_data_frame(pathsrev %>% distinct(),directed=FALSE))
}


# Function to create mutation matrix from tree
get_mutation_matrix = function(g) {
  allpaths = get_paths(g)
  gdf = as_long_data_frame(g)
  leaves = V(g)$name[!V(g)$name %in% gdf$from_name]
  clones = paste0("clone",1:length(leaves))
  clonemat = matrix(NA, ncol = length(clones), nrow = length(V(g)), 
                    dimnames = list(V(g)$name, clones))
  for(i in 1:length(leaves)) {
    # Added a [[1]] after this because this function returns
    # a list with one element
    clonepath = allpaths[sapply(allpaths, function(x) 
      if(tail(x, 1) == leaves[i]) TRUE else FALSE)][[1]]
    clonemat[,i] = as.numeric(rownames(clonemat) %in% clonepath)
  }
  return(clonemat)
}

# Get before/after matrix for pairs of genes
get_before_after_matrix = function(g) {
  allpaths = get_paths(g)
  beforeafter = matrix(NA, ncol = length(V(g)), nrow = length(V(g)), 
                       dimnames = list(V(g)$name, V(g)$name))
  for(i in 1:nrow(beforeafter)) {
    for(j in 1:ncol(beforeafter)) {
      if(j != i) {
        # Limiting to cases where two genes in same path
        pathbool = sapply(allpaths, function(x) 
          V(g)$name[i] %in% x & V(g)$name[j] %in% x)
        if(any(pathbool)) {
          pathtmp = allpaths[pathbool][[1]]
          beforeafter[i,j] = as.numeric(which(pathtmp == V(g)$name[i]) < 
                                          which(pathtmp == V(g)$name[j]))
        }
      }
    }
  }
  return(beforeafter)
}

# Function to count when two mutations occur in same or different clone for a graph
count_same_different = function(g,gene1,gene2) {
  mutmat = get_mutation_matrix(g)
  gene1name = V(g)$name[V(g)$Gene == gene1]
  gene2name = V(g)$name[V(g)$Gene == gene2]
  if(gene1 != gene2) {
    matsub = mutmat[sort(c(gene1name,gene2name)),]
    # Added this in case matsub is a vector because only one clone was involved
    if(is.null(dim(matsub))) matsub = as.matrix(matsub)
    # Limiting to clones with these genes in them
    matsub = matsub[,apply(matsub, 2, sum) != 0]
    # Added this in case matsub is a vector because only one clone was involved
    if(is.null(dim(matsub))) matsub = as.matrix(matsub)
    gene1sub = matsub[gene1name,]
    gene2sub = matsub[gene2name,]
    # Again, converting back to matrix as needed
    if(is.null(dim(gene1sub))) {
      gene1sub = as.matrix(gene1sub)
      if(any(grepl("clone",rownames(gene1sub)))) {
        gene1sub = t(gene1sub)
        rownames(gene1sub) = gene1name
      }
    }
    if(is.null(dim(gene2sub))) {
      gene2sub = as.matrix(gene2sub)
      if(any(grepl("clone",rownames(gene2sub)))) {
        gene2sub = t(gene2sub)
        rownames(gene2sub) = gene2name
      }
    }
    gene1sum = apply(gene1sub, 2, sum)
    gene2sum = apply(gene2sub, 2, sum)
    # First, checking if two genes are completely disjoint, and if so,
    # return their vector result, where total for each = separate for each
    disjoint = sum(gene1sum[gene2sum >= 1]) == 0
    if(disjoint) {
      summ = c(length(gene1name), length(gene2name), length(gene1name), length(gene2name))
      names(summ) = c(paste0(gene1, "_total"), paste0(gene2,"_total"),
                      paste0(gene1, "_separate"), paste0(gene2,"_separate"))
      return(summ)
    }
    # Now keep track of the number of times each gene is disjoint from another
    disjoint1 = 0
    disjoint2 = 0
    # Also have commented code where I keep track of whether in the same path
    #inpath1 = 0
    #inpath2 = 0
    for(i in 1:nrow(gene1sub)) {
      # Check if a specific instance of gene1 is disjoint from all of gene2
      if(sum(gene1sub[i,gene2sum >= 1]) == 0) disjoint1 = disjoint1 + 1
    }
    for(i in 1:nrow(gene2sub)) {
      if(sum(gene2sub[i,gene1sum >= 1]) == 0) disjoint2 = disjoint2 + 1
    }
    summ = c(length(gene1name), length(gene2name), disjoint1, disjoint2)
    names(summ) = c(paste0(gene1, "_total"), paste0(gene2,"_total"),
                    paste0(gene1, "_separate"), paste0(gene2,"_separate"))
    return(summ)
  }
  #print(paste("Checking duplicate",gene1))
  else if(gene1 == gene2 & length(gene1name) > 1) {
    matsub = mutmat[gene1name,]
    # If there's only one clone with both variants, then we need to make this a matrix
    if(is.null(dim(matsub))) matsub = t(as.matrix(matsub))
    ndisjoint = 0
    for(i in 1:nrow(matsub)) {
      # Only checking the variants that we haven't already checked
      othergene = matsub[-i,]
      if(is.null(dim(othergene))) othergene = t(as.matrix(othergene))
      matsubsum = apply(othergene,2, sum)
      if(sum(matsub[i,][matsubsum >= 1]) == 0) {
        ndisjoint = ndisjoint + 1
      }
    }
    summ = c(nrow(matsub), nrow(matsub), ndisjoint, ndisjoint)
    names(summ) = c(paste0(gene1, "_total"), paste0(gene2,"_total"),
                    paste0(gene1, "_separate"), paste0(gene2,"_separate"))
    return(summ)
  }
}

# Function to count when two mutations occur in same or different clone for a graph
count_before_after = function(g,gene1,gene2) {
  mutmat = get_before_after_matrix(g)
  gene1name = V(g)$name[V(g)$Gene == gene1]
  gene2name = V(g)$name[V(g)$Gene == gene2]
  mat1sub = mutmat[gene1name,gene2name]
  mat2sub = mutmat[gene2name,gene1name]
  # Convert matrix to vector as needed
  if(!is.null(dim(mat1sub))) mat1sub = as.vector(mat1sub)
  if(!is.null(dim(mat2sub))) mat2sub = as.vector(mat2sub)
  # Filtering for if the two genes aren't completely disjoint
  if(any(!is.na(mat1sub))) {
    summ = c(sum(na.omit(mat1sub)), sum(na.omit(mat2sub)))
    names(summ) = c(paste0(gene1, "_first"), paste0(gene2,"_first"))
  }
  else summ = rep(NA,2)
  names(summ) = c(paste0(gene1, "_first"), paste0(gene2,"_first"))
  return(summ)
}


get_pairwise_mats = function(graphfiles) {
  graphlist = vector(mode = "list", length = length(graphfiles))
  for(i in 1:length(graphlist)) {
    parsedgv = parsegv(graphfiles[i])
    g = remove_nonpath(graph_from_data_frame(parsedgv$links, directed=TRUE, 
                                             vertices=parsedgv$genes), gvfile = graphfiles[i])
    graphlist[[i]] = g
  }
  treegenes = lapply(graphlist, function(x) V(x)$Gene)
  allgenes = unique(unlist(treegenes))
  allgenes = allgenes[allgenes != "Root"]
  genecombos = t(combn(allgenes, m=2))
  # Add pairwise for each type of gene
  genecombos = rbind(genecombos, matrix(c(allgenes, allgenes),ncol = 2))
  
  in_same_graph = matrix(0, nrow = length(allgenes), ncol = length(allgenes))
  rownames(in_same_graph) = allgenes
  colnames(in_same_graph) = allgenes
  
  in_same_path = matrix(0, nrow = length(allgenes), ncol = length(allgenes))
  rownames(in_same_path) = allgenes
  colnames(in_same_path) = allgenes
  
  earlier_in_path = matrix(NA, nrow = length(allgenes), ncol = length(allgenes))
  rownames(earlier_in_path) = allgenes
  colnames(earlier_in_path) = allgenes
  
  for(i in 1:nrow(genecombos)) {
    gene1 = genecombos[i,1]
    gene2 = genecombos[i,2]
    if(gene1 != gene2) {
      graphsub = graphlist[sapply(treegenes, function(x) 
        gene1 %in% x & gene2 %in% x)]
    }
    else graphsub = graphlist[sapply(treegenes, function(x) 
      length(which(x == gene1)) > 1)]
    
    if(length(graphsub) > 0) {
      for(j in 1:length(graphsub)) {
        samediff = count_same_different(graphsub[[j]],gene1,gene2)
        in_same_graph[gene1,gene2] = in_same_graph[gene1,gene2] + 
          samediff[paste0(gene1,"_total")]
        if(gene2 != gene1) in_same_graph[gene2,gene1] = in_same_graph[gene2,gene1] + 
          samediff[paste0(gene2,"_total")]
        in_same_path[gene1,gene2] = in_same_path[gene1,gene2] + 
          samediff[paste0(gene1,"_total")] - 
          samediff[paste0(gene1,"_separate")]
        if(gene2 != gene1) in_same_path[gene2,gene1] = in_same_path[gene2,gene1] + 
          samediff[paste0(gene2,"_total")] - 
          samediff[paste0(gene2,"_separate")]
        pathcount = count_before_after(graphsub[[j]],gene1,gene2)
        if(any(!is.na(pathcount))) {
          # Restricting only to cases where the two genes are different
          if(gene1 != gene2) {
            if(is.na(earlier_in_path[gene1,gene2])) {
              earlier_in_path[gene1,gene2]=0
              earlier_in_path[gene2,gene1]=0
            }
            earlier_in_path[gene1,gene2] = earlier_in_path[gene1,gene2] + 
              pathcount[paste0(gene1, "_first")]
            earlier_in_path[gene2,gene1] = earlier_in_path[gene2,gene1] + 
              pathcount[paste0(gene2, "_first")]
          }
        }
      }
    }
  }
  # also get the fraction of earlier in path
  earlier_in_path_perc = earlier_in_path
  for(i in 1:nrow(earlier_in_path)) {
    for(j in i:ncol(earlier_in_path)) {
      if(j != i & !is.na(earlier_in_path[i,j])) {
        allevents = earlier_in_path[i,j] + earlier_in_path[j,i]
        earlier_in_path_perc[i,j] = earlier_in_path[i,j]/allevents
        earlier_in_path_perc[j,i] = 1 - earlier_in_path_perc[i,j]
      }
    }
  }
  # And the fraction of in the same path
  in_same_path_perc = matrix(-1, nrow = nrow(in_same_path), 
                             ncol = ncol(in_same_path), dimnames = list(allgenes, allgenes))
  for(i in 1:nrow(in_same_path)) {
    for(j in i:ncol(in_same_path)) {
      if(j != i & in_same_graph[i,j] != 0) {
        allevents = in_same_graph[i,j]
        in_same_path_perc[i,j] = in_same_path[i,j]/allevents
        in_same_path_perc[j,i] = in_same_path_perc[i,j]
      }
    }
  }
  in_same_path_perc[in_same_path_perc == -1] = NA
  return(list(in_same_graph=in_same_graph,
              in_same_path=in_same_path,
              in_same_path_perc=in_same_path_perc, 
              earlier_in_path=earlier_in_path,
              earlier_in_path_perc=earlier_in_path_perc))
}

# Function to exclude low weight nodes from a graph
exclude_low_weight = function(g, minweight) {
  lownodes = V(g)$name[V(g)$Count <= minweight]
  for(i in 1:length(lownodes)) {
    g = remove_node(g, lownodes[i])
  }
  return(g)
}

import_merge_graphs = function(graphfiles, analysisvar) {
  graphlist = vector(mode = "list", length = length(graphfiles))
  for(i in 1:length(graphlist)) {
    parsedgv = parsegv(graphfiles[i])
    graphlist[[i]] = remove_nonpath(graph_from_data_frame(parsedgv$links, 
                                                          directed=TRUE, vertices=parsedgv$genes), gvfile = graphfiles[i])
  }
  graphmerged = merge_graphs(graphlist[[1]], graphlist[[2]],analysisvar)
  #nvars = data.frame(graphfiles, number = NA)
  #nvars$number[1] = length(V(graphlist[[1]]))
  #nvars$number[2] = length(V(graphlist[[2]]))
  for(i in 3:length(graphlist)) {
    #	nvars$number[i] = length(V(graphlist[[i]]))
    if(length(V(graphlist[[i]])) > 1) {
      graphmerged = merge_graphs(graphmerged, graphlist[[i]],analysisvar)
    }
  }
  #return(list(graphmerged = graphmerged, nvars = nvars))
  return(graphmerged)
}


# And matrix of first features
get_first_var = function(gfile, analysisvar) {
  parsed = parsegv(gfile)
  g = remove_nonpath(graph_from_data_frame(parsed$links, directed=TRUE, 
                                           vertices=parsed$genes), gvfile = gfile)
  rootnode = V(g)$name[V(g)$Gene == "Root"]
  firstnodes = adjacent_vertices(g, rootnode, mode = "out")[[1]]$name
  if(analysisvar == "Gene") firstvars = V(g)$Gene[V(g)$name %in% firstnodes]
  else if(analysisvar == "Pathway") firstvars = V(g)$Pathway[V(g)$name %in% firstnodes]
  else stop("analysisvar should be Gene or Pathway")
  return(firstvars)
}

get_first_vars = function(gfiles, analysisvar) {
  firstvars = lapply(gfiles, get_first_var, analysisvar=analysisvar)
  allvars = unique(unlist(firstvars))
  firstvarmat = matrix(0, nrow = length(gfiles), ncol = length(allvars),
                       dimnames = list(gfiles, allvars))
  for(i in 1:nrow(firstvarmat)) {
    if(length(firstvars[[i]]) == 1) firstvarmat[i,firstvars[[i]]] = 1
    else if(length(firstvars[[i]]) > 1) { 
      for(j in 1:length(firstvars[[i]])) {
        firstvarmat[i,firstvars[[i]][j]] = firstvarmat[i,firstvars[[i]][j]]+1
      }
    }
  }
  return(firstvarmat)
}


# function to get % cells with mutation among non-missing
# and also to calculate overlap between two genes
get_overlap_genes = function(mutmattmp, gene1, gene2, vartmp) {
  if(!vartmp %in% c("ncells","percs")) stop("vartmp should be ncells or percs")
  gene1vals = as.numeric(mutmattmp[mutmattmp$AAchange == gene1,-c(1:2)])
  gene2vals = as.numeric(mutmattmp[mutmattmp$AAchange == gene2,-c(1:2)])
  gene1mut = sum(gene1vals %in% c(1,2))/sum(gene1vals %in% c(0,1,2))
  gene2mut = sum(gene2vals %in% c(1,2))/sum(gene2vals %in% c(0,1,2))
  lims = gene2vals != 3 & gene1vals != 3
  gene1edit = as.numeric(gsub(2, 1, gene1vals[lims]))
  gene2edit = as.numeric(gsub(2, 1, gene2vals[lims]))
  genecomp = sum((gene2vals[lims] == 1) & 
                   (gene1vals[lims] == 1))/sum((gene2vals[lims] == 1))
  returnvec = c(gene1mut, gene2mut, genecomp)
  names(returnvec) = c(gene1, gene2, paste0(gene2, " in ", gene1))
  if(vartmp == "percs") return(returnvec)
  if(vartmp == "ncells") return(c(returnvec[1]*sum(gene1vals != 3),
                                  returnvec[2]*sum(gene2vals != 3), returnvec[3]*sum(gene2vals[lims] == 1)))
}

# Compare patterns vs. VAF
get_vaf_mut_perc = function(gvfile) {
  driver_vars_tmp = driver_vars_tab %>% filter(Tree == gvfile & !is.na(Bulk_VAF))
  if(nrow(driver_vars_tmp) > 0) {
    parsedgv = parsegv(gvfile)
    g = graph_from_data_frame(parsedgv$links, 
                              directed=TRUE, vertices=parsedgv$genes)
    allpaths = get_all_directed_links(g)
    vardict = V(g)$Variant
    names(vardict) = V(g)$name
    for(i in 1:length(vardict)) allpaths[allpaths == names(vardict)[i]] = vardict[i]
    allpaths = allpaths %>% filter(to %in% driver_vars_tmp$AAchange &
                                     from %in% driver_vars_tmp$AAchange)
    if(nrow(allpaths) > 0) {
      allpaths = allpaths %>% 
        dplyr::rename(AAchange = from) %>%
        left_join(driver_vars_tmp %>% 
                    select(AAchange, ncells, Mutated_Percentage, Missing_perc, Bulk_VAF) %>%
                    filter(!is.na(Bulk_VAF)) %>% rename_with(~ paste0(.x, "_from"), 
                                                             !(starts_with("nce") | starts_with("AAc") |
                                                                 starts_with("Sam")))) %>%
        dplyr::rename(from = AAchange) %>% dplyr::rename(AAchange = to) %>%
        left_join(driver_vars_tmp %>% 
                    select(Sample, AAchange, ncells, Mutated_Percentage, Missing_perc, Bulk_VAF) %>%
                    rename_with(~ paste0(.x, "_to"), 
                                !(starts_with("nce") | starts_with("AAc") |
                                    starts_with("Sam")))) %>%
        dplyr::rename(to = AAchange)
      return(data.frame(gvfile, allpaths))
    }
    return(NULL)
  }
  return(NULL)
}

# get descendents of a branch point
get_branch_descendants = function(g) {
  gdf = as_long_data_frame(g)
  branchpoints = unique(gdf$from_name[duplicated(gdf$from_name)])
  if(length(branchpoints > 0)) {
    alldesc = NULL
    for(i in 1:length(branchpoints)) {
      newdesc = subcomponent(g, 
                             branchpoints[i], mode = "out")$name
      newdesc = newdesc[newdesc != branchpoints[i]]
      alldesc = c(alldesc, newdesc)
    }
    alldesc = unique(alldesc)
    newgenes = V(g)$Gene[V(g)$name %in% alldesc]
    return(newgenes)
  }
  else return(NULL)
}

# Rename starting node to root
rename_to_root = function(g) {
  V(g)[degree(g, mode="in") == 0]$Pathway = "Root"
  V(g)[degree(g, mode="in") == 0]$Gene = "Root"
  return(g)
}

# Function to get trees where the root node's neighbors meet a certain condition
bool_root_neighbor_condition = function(g, variant, analysisvar) {
  root_node = V(g)[degree(g, mode = "in") == 0]
  # Identify neighbors of the root node
  neighbors_of_root = neighbors(g, v = root_node, mode = "out")
  # Filter neighbors by the condition Gene == "X"
  if(analysisvar == "Gene") return(any(neighbors_of_root$Gene == variant))
  return(any(neighbors_of_root$Pathway == variant))
}

# get descendants of a branch point
get_merge_all_subtrees = function(gfiles, feature, analysisvar) {
  graphlist = vector(mode = "list")
  j = 1
  for(i in 1:length(gfiles)) {
    parsedgv = parsegv(gfiles[i])
    g = remove_nonpath(graph_from_data_frame(parsedgv$links, 
                                             directed=TRUE, vertices=parsedgv$genes), gvfile = gfiles[i])
    gdf = as_long_data_frame(g)
    if(feature %in% c(gdf$from_Gene, gdf$from_Pathway)) {
      graphlist[[j]] = get_subgraph_any(g, gdf, feature, analysisvar)
      j = j+1
    }
  }
  graphlist = graphlist[!sapply(graphlist, is.null)]
  # Rename starting node to root to take advantage of prior code
  graphlist = lapply(graphlist, rename_to_root)
  # Check times when the feature of interest is also a neighbor
  sapply(graphlist, bool_root_neighbor_condition, 
         variant = feature, analysisvar = analysisvar)
  # Separate multiple disjoint graphs into distinct graph objects
  graphlist = unlist(lapply(graphlist, decompose), recursive=FALSE)
  graphmerged = merge_graphs(graphlist[[1]], graphlist[[2]],analysisvar)
  for(i in 3:length(graphlist)) {
    if(length(V(graphlist[[i]])) > 1) {
      graphmerged = merge_graphs(graphmerged, graphlist[[i]],analysisvar)
    }
  }
  if(analysisvar == "Pathway") {
    V(graphmerged)[degree(graphmerged, mode="in") == 0]$Pathway = feature
  }
  if(analysisvar == "Gene") {
    V(graphmerged)[degree(graphmerged, mode="in") == 0]$Gene = feature
  }
  return(graphmerged)
}

get_merge_root_subtrees = function(gfiles, feature, analysisvar) {
  graphlist = vector(mode = "list")
  j = 1
  for(i in 1:length(gfiles)) {
    parsedgv = parsegv(gfiles[i])
    g = remove_nonpath(graph_from_data_frame(parsedgv$links, 
                                             directed=TRUE, vertices=parsedgv$genes), gvfile = gfiles[i])
    gdf = as_long_data_frame(g)
    gdf_root = gdf %>% filter(from_Gene == "Root")
    if(feature %in% c(gdf_root$to_Gene, gdf_root$to_Pathway)) {
      new_g = get_subgraph_root(g, gdf, feature, analysisvar)
      if(!is.null(new_g)) {
        if(bool_root_neighbor_condition(new_g,
                                        variant = feature, analysisvar = analysisvar)) {
          print(gfiles[i])
        }
      }
      graphlist[[j]] = new_g
      j = j+1
    }
  }
  graphlist = graphlist[!sapply(graphlist, is.null)]
  # Rename starting node to root to take advantage of prior code
  graphlist = lapply(graphlist, rename_to_root)
  # Separate multiple disjoint graphs into distinct graph objects
  graphlist = unlist(lapply(graphlist, decompose), recursive=FALSE)
  graphmerged = merge_graphs(graphlist[[1]], graphlist[[2]],analysisvar)
  for(i in 3:length(graphlist)) {
    
    if(length(V(graphlist[[i]])) > 1) {
      graphmerged = merge_graphs(graphmerged, graphlist[[i]],analysisvar)
    }
  }
  if(analysisvar == "Pathway") {
    V(graphmerged)[degree(graphmerged, mode="in") == 0]$Pathway = gsub(" ","", feature)
  }
  if(analysisvar == "Gene") {
    V(graphmerged)[degree(graphmerged, mode="in") == 0]$Gene = gsub(" ","", feature)
  }
  return(graphmerged)
}


# Function to get all the branch points
get_all_branch_points = function(gfiles) {
  branchpoint = NULL
  for(i in 1:length(gfiles)) {
    parsedgv = parsegv(gfiles[i])
    g = remove_nonpath(graph_from_data_frame(parsedgv$links, 
                                             directed=TRUE, vertices=parsedgv$genes), gvfile = gfiles[i])
    gdf = as_long_data_frame(g)
    branchpoints = unique(gdf$from_name[duplicated(gdf$from_name)])
    branchgenes = V(g)$Gene[V(g)$name %in% branchpoints]
    branchpoint = c(branchpoint, branchgenes)
  }
  return(branchpoint)
}

# Function to get all the branch points
get_all_npm1 = function(gfiles) {
  npm1_branch = c()
  npm1_children = vector(mode = "list")
  for(i in 1:length(gfiles)) {
    parsedgv = parsegv(gfiles[i])
    g = remove_nonpath(graph_from_data_frame(parsedgv$links, 
                                             directed=TRUE, vertices=parsedgv$genes), gvfile = gfiles[i])
    if("NPM1" %in% V(g)$Gene) {
      gdf = as_long_data_frame(g)
      npm1_kids = gdf %>% filter(from_Gene == "NPM1") %>% pull(to_Gene)
      if(length(npm1_kids) >= 1) {
        if(length(npm1_kids) > 1) {
          npm1_branch[i] = 1
        }
        else {
          npm1_branch[i] = 0
        }
        npm1_children[[i]] = npm1_kids
      }
    }
  }
  return(list(npm1_branch, npm1_children))
}

convert_pathway_names = function(pathway) {
  if(pathway == "RTK/RAS/MAP kinase pathway") return("Signaling")
  else if(!pathway %in% c("DNA methylation", "NPM1")) return("Other")
  else return(pathway)
}

get_diff_mut_perc = function(treefiles, genes1, genes2, diagnosis = FALSE) {
  graphlist = vector(mode = "list", length = length(treefiles))
  linkslist = vector(mode = "list", length = length(treefiles))
  # Be able to filter the driver variants
  if(diagnosis) driver_vars_tab_tmp = driver_vars_tab %>% filter(diagnosis == 1)
  else driver_vars_tab_tmp = driver_vars_tab
  
  for(i in 1:length(graphlist)) {
    # Focus on a specific tree
    driver_vars_tmp = driver_vars_tab_tmp %>% filter(Tree == treefiles[i])
    # Skip if the tree doesn't have the genes of interest
    if(!any(grepl(paste(genes1, collapse="|"), driver_vars_tmp$varvec)) |
       !any(grepl(paste(genes2, collapse="|"), driver_vars_tmp$varvec))) {
      graphlist[[i]] = NULL
      linkslist[[i]] = NULL
    }
    else {
      parsedgv = parsegv(treefiles[i])
      graphlist[[i]] = remove_nonpath(graph_from_data_frame(parsedgv$links, 
                                                            directed=TRUE, vertices=parsedgv$genes), gvfile = treefiles[i])
      g = graphlist[[i]]
      all_links_tmp = get_all_directed_links(g) %>% select(-weight)
      # Create the directed links by AA change
      for(k in 1:length(g)) {
        all_links_tmp$from[all_links_tmp$from == V(g)$name[k]] = V(g)$Variant[k]
        all_links_tmp$to[all_links_tmp$to == V(g)$name[k]] = V(g)$Variant[k]
      }
      # Only focus on those links that are in our driver variants
      all_links_tmp = all_links_tmp %>% 
        filter(from %in% driver_vars_tmp$AAchange & 
                 to %in% driver_vars_tmp$AAchange)
      all_links_tmp1 = all_links_tmp %>% 
        filter(grepl(paste(c(genes1), collapse="|"), from) &
                 grepl(paste(c(genes2), collapse="|"), to)) %>% 
        mutate(gene1first = 1)
      all_links_tmp = all_links_tmp %>% 
        filter(grepl(paste(c(genes2), collapse="|"), from) &
                 grepl(paste(c(genes1), collapse="|"), to)) %>%
        mutate(gene1first = 0) %>% full_join(all_links_tmp1)
      # Only proceed if we have remaining links
      if(nrow(all_links_tmp) > 0) {
        all_samps = unique(driver_vars_tmp$Sample)
        all_links_df = as.data.frame(matrix(nrow = 0, ncol = 7))
        colnames(all_links_df) = c("from", "to", "gene1first", "gene1_perc", 
                                   "gene1_vaf", "gene2_perc", "gene2_vaf")
        # Now for each sample, creating data frame with links and 
        # variables above, and rbind all the samples together
        for(samp_tmp in all_samps) {
          all_links_df = rbind(all_links_df, all_links_tmp %>% dplyr::rename(AAchange = from) %>%
                                 inner_join(driver_vars_tmp %>% filter(Sample == samp_tmp) %>%
                                              select(AAchange, Mutated_Percentage, Bulk_VAF, Missing_perc)) %>%
                                 dplyr::rename(gene1_perc = Mutated_Percentage, gene1_vaf = Bulk_VAF, 
                                               from = AAchange, AAchange = to) %>%
                                 # Denominator is 1 minus the missing percentage of cells
                                 mutate(gene1_perc = gene1_perc/(1-Missing_perc)) %>% select(-Missing_perc) %>%
                                 inner_join(driver_vars_tmp %>% filter(Sample == samp_tmp) %>%
                                              select(AAchange, Mutated_Percentage, Bulk_VAF, Missing_perc)) %>%
                                 dplyr::rename(gene2_perc = Mutated_Percentage, gene2_vaf = Bulk_VAF, 
                                               to = AAchange) %>%
                                 mutate(gene2_perc = gene2_perc/(1-Missing_perc)) %>% select(-Missing_perc))
        }
        linkslist[[i]] = data.frame(Tree = treefiles[i], all_links_df)
      }
    }
  }
  # At the end, rbind all of these data frames together
  return(do.call(rbind, linkslist))
}


get_early_late_none = function(tree_file, gene_set_list) {
  parsed = parsegv(tree_file)
  g = remove_nonpath(graph_from_data_frame(parsed$links, 
                                           directed=TRUE, vertices=parsed$genes), gvfile = tree_file)
  root_neighbors = neighbors(g, which(V(g)$Gene == "Root"))$Gene
  other_genes = setdiff(V(g)$Gene, root_neighbors)
  gene_set_vec = c()
  for(i in 1:length(gene_set_list)) {
    if(length(intersect(gene_set_list[[i]], root_neighbors))>0) {
      gene_set_vec = c(gene_set_vec, "Early")
    }
    else if(length(intersect(gene_set_list[[i]], other_genes))>0) {
      gene_set_vec = c(gene_set_vec, "Late")
    }
    else {
      gene_set_vec = c(gene_set_vec, "None")
    }
  }
  names(gene_set_vec) = names(gene_set_list)
  return(gene_set_vec)
}

write_regression_output = function(linear_mod, file_name) {
  write.table(data.frame(summary(linear_mod)$coefficient) %>%
                mutate(p = signif(Pr...t.., 3), 
                       Estimate = signif(Estimate, 3)), 
              file_name, sep = "\t")
}

# Finally, compare to clinical variables
# First create matrix of all the mutation orders
get_links_matrix = function(gfiles, analysisvar) {
  allords = NULL
  graphfiles = vector(mode = "list")
  graphlinks = vector(mode = "list")
  for(i in 1:length(gfiles)) {
    parsed = parsegv(gfiles[i])
    graphfiles[[i]] = remove_nonpath(graph_from_data_frame(parsed$links, 
                                                           directed=TRUE, vertices=parsed$genes), gvfile = gfiles[i])
    g = graphfiles[[i]]
    rootnode = V(g)$name[V(g)$Gene == "Root"]
    if(length(V(g)) > 1) {
      graphlinks[[i]] = get_all_directed_links(g) %>% filter(from != rootnode)
      if(analysisvar == "Gene") vardict = V(g)$Gene
      else if(analysisvar == "Pathway") vardict = V(g)$Pathway
      else stop("analysisvar should be either 'Gene' or 'Pathway'")
      names(vardict) = V(g)$name
      graphlinks[[i]]$from = vardict[graphlinks[[i]]$from]
      graphlinks[[i]]$to = vardict[graphlinks[[i]]$to]
    }
    else {	
      graphlinks[[i]] = as.data.frame(matrix(nrow=0,ncol=3))
      names(graphlinks[[i]]) = c("to","from","weight")
    }
  }
  
  linksdf = do.call(rbind,graphlinks) %>% group_by(from,to) %>% 
    summarize(weight=sum(weight)) %>% ungroup()
  allords = apply(linksdf %>% select(from, to),1,paste,collapse="->")
  
  ordmat = matrix(NA, nrow = length(gfiles), ncol = length(allords))
  colnames(ordmat) = allords
  rownames(ordmat) = gfiles
  for(i in 1:nrow(ordmat)) {
    if(nrow(graphlinks[[i]]) > 0) {
      glinks = graphlinks[[i]]
      g = graphfiles[[i]]
      linkvecs = apply(glinks %>% select(from, to), 1, paste, collapse="->")
      linkstab = table(linkvecs)
      # Now that we have the links of interest, get all other links that could have occurred
      if(analysisvar == "Gene") linksall = rbind(t(combn(V(g)$Gene, 2)), 
                                                 t(apply(combn(V(g)$Gene, 2), 2, rev)))
      else if(analysisvar == "Pathway") linksall = rbind(t(combn(V(g)$Pathway, 2)), 
                                                         t(apply(combn(V(g)$Pathway, 2), 2, rev)))
      linksall = apply(linksall, 1, paste, collapse="->")
      linksnot = linksall[linksall %in% colnames(ordmat) & 
                            !linksall %in% linkvecs]
      if(length(linkvecs) > 0) ordmat[i,names(linkstab)] = as.numeric(linkstab)
      if(length(linksnot) > 0) ordmat[i,linksnot] = 0
    }
  }
  return(ordmat)
}
