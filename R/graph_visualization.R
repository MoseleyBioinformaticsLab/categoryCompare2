#' unique combinations
#'
#' determine the unique combinations of annotations that exist in the
#' significant matrix of the \code{\link{cc_graph}} and assign each node in the graph
#' to a group.
#'
#' @param object the \code{\link{cc_graph}} to work on
#'
#' @return node_assignment
#' @exportMethod annotation_combinations
#' @rdname annotation_combinations
setMethod("annotation_combinations",
          signature = list(object = "cc_graph"),
          function(object) .annotation_combinations(object@significant))

check_package_installed = function(packagename){
  is_install = try(find.package(packagename))
  
  if (inherits(is_install, "try-error")) {
    stop_message = paste0(packagename, " is not installed, but is required for this operation!")
    stop(stop_message)
  } else {
    return(NULL)
  }
}


#' unique combinations
#'
#' determine the unique combinations of annotations that exist in the
#' significant matrix of the \code{\link{combined_statistics}} and assign each
#' annotation to a group.
#'
#' @param object the \code{\link{combined_statistics}} to work on
#'
#' @return node_assignment
#' @exportMethod annotation_combinations
#' @rdname annotation_combinations
setMethod("annotation_combinations",
          signature = list(object = "significant_annotations"),
          function(object) .annotation_combinations(object@significant))

.annotation_combinations <- function(sig_matrix){
  if (nrow(sig_matrix) == 0) {
    return(NULL)
  }

  unique_combinations <- unique(sig_matrix)
  # after generating the unique_combinations, we then want them ordered such
  # that a FALSE, FALSE, ... would be last, because this is often not that
  # interesting. This does that.
  # We can only do this if there is more than one row.
  if ((nrow(unique_combinations) > 1) && (ncol(unique_combinations) > 1)) {
    
    n_col <- ncol(unique_combinations)
    uniq_order <- do.call(order, c(lapply(1:n_col, function(i) unique_combinations[, i]), decreasing=TRUE))
    unique_combinations <- unique_combinations[uniq_order, ]
  } else {
    uniq_order <- order(unique_combinations, decreasing = TRUE)
    unique_combinations <- unique_combinations[uniq_order, , drop = FALSE]
  }
  
  name_combinations <- paste("G", seq(1, nrow(unique_combinations)), sep = "")

  rownames(unique_combinations) <- name_combinations

  # also generate a textual description of each group that can be added 
  # to the graph visualization
  tmp_names <- colnames(unique_combinations)
  group_description <- vapply(rownames(unique_combinations), function(in_comb){
    use_comb <- unique_combinations[in_comb, ]
    paste(tmp_names[use_comb], collapse = ",")
  }, character(1))
  
  # initialize the things that store our assignments to GO terms
  combination_assign <- rep("G", nrow(sig_matrix))
  names(combination_assign) <- rownames(sig_matrix)
  
  description_assign <- rep("G", nrow(sig_matrix))
  names(description_assign) <- rownames(sig_matrix)

  for (in_comb in name_combinations){
    has_match <- apply(sig_matrix, 1, function(in_sig){
      in_sig == unique_combinations[in_comb, ]
    })
    combination_assign[has_match] <- in_comb
    description_assign[has_match] <- group_description[in_comb]
  }
  
  
  new("node_assign", groups = unique_combinations, assignments = combination_assign,
      description = description_assign)
}

#' generate colors
#'
#' given a bunch of items, generate a set of colors for either single node colorings
#' or pie-chart annotations. Colors are generated using the \emph{hcl} colorspace,
#' and for \code{n_color >= 5}, the colors are re-ordered in an attempt to create
#' the largest contrasts between colors, as they result from being picked on a
#' circle in \emph{hcl} space.
#'
#' @param n_color how many colors to generate
#'
#' @export
#' @importFrom colorspace rainbow_hcl
generate_colors <- function(n_color){
  out_color <- rainbow_hcl(n_color, c = 100)

  if (n_color <= 4){
    return(out_color)
  }

  out_index <- seq(1, n_color)
  tmp_index <- out_index

  if ((n_color %% 2) == 1){
    swap_index_1 <- seq(2, n_color - 1, 2)
    swap_index_2 <- swap_index_1 + 1

    for (i_swap in seq(1, length(swap_index_1))){
      out_index[swap_index_1[i_swap]] <- tmp_index[swap_index_2[i_swap]]
      out_index[swap_index_2[i_swap]] <- tmp_index[swap_index_1[i_swap]]
      tmp_index <- out_index
    }

    swap_index_1 <- seq(3, n_color - 1, 2)
    swap_index_2 <- swap_index_1 + 1

    for (i_swap in seq(1, length(swap_index_1))){
      out_index[swap_index_1[i_swap]] <- tmp_index[swap_index_2[i_swap]]
      out_index[swap_index_2[i_swap]] <- tmp_index[swap_index_1[i_swap]]
      tmp_index <- out_index
    }

  } else {
    swap_index_1 <- seq(1, n_color, 2)
    swap_index_2 <- swap_index_1 + 1

    for (i_swap in seq(1, length(swap_index_1))){
      out_index[swap_index_1[i_swap]] <- tmp_index[swap_index_2[i_swap]]
      out_index[swap_index_2[i_swap]] <- tmp_index[swap_index_1[i_swap]]
      tmp_index <- out_index
    }

    swap_index_1 <- seq(2, n_color - 2, 3)
    swap_index_2 <- seq(n_color - 1, 3, -3)

    for (i_swap in seq(1, length(swap_index_1))){
      out_index[swap_index_1[i_swap]] <- tmp_index[swap_index_2[i_swap]]
      out_index[swap_index_2[i_swap]] <- tmp_index[swap_index_1[i_swap]]
      tmp_index <- out_index
    }

  }

  out_color <- out_color[out_index]
  return(out_color)
}

#' assign colors
#'
#' given a \code{\link{node_assign}}, assign colors to either the independent groups
#' of unique annotations, or to each of the experiments independently.
#'
#' @param in_assign the \code{\link{node_assign}} object generated from a \code{\link{cc_graph}}
#' @param type either "group" or "experiment"
#'
#' @export
#' @return node_assign with colors
assign_colors <- function(in_assign, type = "experiment"){
  grp_matrix <- in_assign@groups

  if (type == "experiment"){
    n_color <- ncol(grp_matrix)
    use_color <- generate_colors(n_color)
    names(use_color) <- colnames(grp_matrix)

    in_assign@colors <- use_color
    in_assign@color_type <- "pie"
    in_assign@pie_locs <- generate_piecharts(grp_matrix, use_color)
  } else {
    n_color <- nrow(grp_matrix)
    use_color <- generate_colors(n_color)
    names(use_color) <- rownames(grp_matrix)
    in_assign@colors <- use_color
    in_assign@color_type <- "solid"
  }

  return(in_assign)
}

#' create piecharts for visualization
#'
#' given a group matrix and the colors for each experiment, generate the pie graphs
#' that will be used as glyphs in Cytoscape
#'
#' this should \emph{not be exported in the final version}
#'
#' @param grp_matrix the group matrix
#' @param use_color the colors for each experiment
#'
#' @export
#' @return list of png files that are pie graphs
#' @importFrom colorspace desaturate
generate_piecharts <- function(grp_matrix, use_color){
  check_package_installed("Cairo")
  n_grp <- nrow(grp_matrix)
  n_color <- length(use_color)
  
  # defines how many pie segments are needed, common to all the pie-charts
  pie_area <- rep(1 / n_color, n_color)
  names(pie_area) <- rep("", n_color) # add blank names so nothing gets printed
  
  # use desaturated version of colors when there is non-significance
  desat_color <- desaturate(use_color)
  names(desat_color) <- names(use_color)
  piechart_strings <- purrr::map_dfr(rownames(grp_matrix), function(i_grp){
    tmp_logical <- grp_matrix[i_grp, ]
    tmp_color <- use_color
    
    # add the proper desaturated versions of the colors
    tmp_color[!tmp_logical] <- desat_color[!tmp_logical]
    
    out_str = paste0('piechart: attributelist="',
                     paste(colnames(grp_matrix), collapse = ','),
                     '" ',
                     'colorlist="',
                     paste(tmp_color, collapse = ','),
                     '" ',
                     'arcstart=-90 showlabels=false')
    
    out_file <- tempfile(i_grp, fileext = ".png")
    Cairo::Cairo(width = 640, height = 640, file = out_file, type = "png", bg = "transparent")
    par(mai = c(0, 0, 0, 0))
    pie(pie_area, col = tmp_color, clockwise = TRUE)
    dev.off()
    
    if (Sys.info()['sysname'] == "Windows") {
      out_file <- gsub("\\", "/", out_file, fixed = TRUE)
      out_file <- paste0("file:///", out_file)
    } else {
      out_file <- paste0("file://localhost", out_file)
    }
    data.frame(colorlist = out_str, piechart = out_file)
    
  })
  
  tmp_matrix = as.data.frame(matrix(1, nrow = nrow(piechart_strings),
                      ncol = ncol(grp_matrix)))
  names(tmp_matrix) = colnames(grp_matrix) 
  piechart_df = cbind(tmp_matrix, piechart_strings)
  piechart_df$visattr = rownames(grp_matrix)
  return(piechart_df)
}


#' add tooltip
#' 
#' before passing to Cytoscape, add a tooltip attribute to the graph
#' 
#' @param in_graph the graph to work with
#' @param node_data which pieces of node data to use
#' @param description other descriptive text to use
#' @param separator what separator to use for the tooltip
#' 
#' @return the graph with a new nodeData member "tooltip"
#' 
add_tooltip <- function(in_graph, node_data = c("name", "description"), description, separator = "\n"){
  use_nodes <- graph::nodes(in_graph)
  n_nodes <- length(use_nodes)
  tooltips <- vapply(use_nodes, function(in_node){
    out_tooltip <- ""
    for (i_dat in node_data) {
      out_tooltip <- paste0(out_tooltip, graph::nodeData(in_graph, in_node, i_dat), "\n")
    }
    out_tooltip <- paste0(out_tooltip, description[in_node])
    out_tooltip
  }, character(1))
  
  graph::nodeDataDefaults(in_graph, "tooltip") <- "NA"
  attr(graph::nodeDataDefaults(in_graph, "tooltip"), "class") <- "STRING"
  graph::nodeData(in_graph, use_nodes, "tooltip") <- tooltips
  
  in_graph
}

#' visualize in cytoscape
#'
#' given a graph, and the node assignments, visualize the graph in cytoscape
#' for manipulation
#'
#' @param in_graph the cc_graph to visualize
#' @param in_assign the node_assign generated
#' @param description something descriptive about the vis (useful when lots of different visualizations)
#'
#' @export
#' @return something
vis_in_cytoscape <- function(in_graph, in_assign, description = "cc2 enrichment"){

  check_package_installed("RCy3")
  in_graph <- add_tooltip(in_graph, description = in_assign@description)
  # initialize and add the visual attribute so we can color according to the
  # data that lives in in_assign
  graph::nodeDataDefaults(in_graph, "visattr") <- ""
  attr(graph::nodeDataDefaults(in_graph, "visattr"), "class") <- 'STRING'
  graph::nodeData(in_graph, names(in_assign@assignments), "visattr") <- in_assign@assignments

  cyt_window <- RCy3::createNetworkFromGraph(graph = in_graph, title = description)
  
  RCy3::layoutNetwork('force-directed', network = cyt_window)
  
  RCy3::setNodeTooltipMapping("tooltip", network = cyt_window)

  if (in_assign@color_type == "solid"){
    RCy3::setNodeColorMapping("visattr", names(in_assign@colors), colors = in_assign@colors, mapping.type = "discrete", network = cyt_window)
  } else if (in_assign@color_type == "pie"){
    RCy3::setNodeShapeDefault("ELLIPSE")
    RCy3::setNodeSizeDefault(35)
    all_nodes = graph::nodes(in_graph)
    pie_data = in_assign@pie_locs
    pie_names = names(pie_data)
    pie_names = pie_names[!(pie_names %in% "piechart")]
    node_vis = graph::nodeData(in_graph, all_nodes, "visattr")
    node_vis_df = data.frame(name = names(node_vis), visattr = unlist(node_vis))
    node_vis_df = dplyr::left_join(node_vis_df, pie_data[, pie_names], by = "visattr")
    RCy3::loadTableData(node_vis_df, data.key.column = "name", table = "node", table.key.column = "name")
    RCy3::updateStyleMapping("default", RCy3::mapVisualProperty("NODE_CUSTOMGRAPHICS_1", "colorlist", "p"))

  }
  
  return(cyt_window)
}

#' remove edges
#'
#' given a RCy3 network connection, remove edges according to provided
#' values.
#'
#' @param edge_obj a network ID to an RCy3 Cytoscape connection
#' @param cutoff what cutoff to use to remove edges
#' @param edge_attr what attribute has the values
#' @param value_direction remove those edges "under" or "over" the value
#'
#' @export
#' @return nothing
#' @rdname remove_edges
#' @aliases remove_edges
setMethod("remove_edges", signature=list(edge_obj="character", cutoff="numeric"), function(edge_obj, cutoff, edge_attr, value_direction)
  .remove_edges_cw(edge_obj, cutoff, edge_attr, value_direction))

.remove_edges_cw <-	function(edge_obj, cutoff, edge_attr = "weight", value_direction = "under"){
  edge_data <- RCy3::getTableColumns('edge', network = edge_obj)

  switch(value_direction,
         under = to_delete <- edge_data[(as.numeric(edge_data[, edge_attr]) < cutoff),],
         over = to_delete <- edge_data[(as.numeric(edge_data[, edge_attr]) > cutoff),]
  )

  RCy3::selectEdges(to_delete$SUID, by.col = "SUID", network = edge_obj)
  RCy3::deleteSelectedEdges(network = edge_obj)

  RCy3::layoutNetwork('force-directed', network = edge_obj)

  message("Removed ", nrow(to_delete), " edges from graph\n")
}

#' remove graph edges
#'
#' @param edge_obj cc_graph
#' @param cutoff the cutoff to use
#' @param edge_attr which attribute to use
#' @param value_direction remove edges with value under or over
#'
#' @export
#' @return cc_graph
#' @rdname remove_edges
setMethod("remove_edges", signature=list(edge_obj="cc_graph", cutoff="numeric"), function(edge_obj, cutoff, edge_attr, value_direction)
  .remove_edges_ccgraph(edge_obj, cutoff, edge_attr, value_direction))

.remove_edges_ccgraph <-	function(in_graph, cutoff, edge_attr = "weight", value_direction = "under"){
  edge_data <- unlist(graph::edgeData(in_graph, , , edge_attr))

  switch(value_direction,
         under = del_edges <- names(edge_data)[edge_data < cutoff],
         over = del_edges <- names(edge_data)[edge_data > cutoff]
  )

  if (length(del_edges) > 0){
    del_edges <- strsplit(del_edges, "|", fixed = TRUE)
    from_node <- sapply(del_edges, function(x){x[1]})
    to_node <- sapply(del_edges, function(x){x[2]})
    in_graph <- graph::removeEdge(from_node, to_node, in_graph)
  }

  message("Removed ", length(del_edges), " edges from graph\n")
  return(in_graph)
}

#' generate a legend
#'
#' it often helps to have a legend displayed for reference.
#'
#' @param in_assign the assign object from \code{annotation_combinations}
#' @param upper_names whether to make names uppercase for easier viewing
#' @param img should a base64 encoded data uri be returned for embedding?
#' @param width how wide should the image be if saving to an image
#' @param height how high should it be
#' @param pointsize the pointsize parameter for Cairo, determines textsize in the image
#' @param ... any other parameter to \code{pie}
#'
#' @return NULL
#' @export
generate_legend <- function(in_assign, upper_names = TRUE, img = FALSE,
                            width = 800, height = 400, pointsize = 70, ...){
  
  if (in_assign@color_type == "pie") {
    use_color <- in_assign@colors
    n_color <- length(use_color)

    # defines how many pie segments are needed, common to all the pie-charts
    pie_area <- rep(1 / n_color, n_color)

    use_labels <- names(use_color)
    if (upper_names) {
      use_labels <- toupper(use_labels)
    }

    if (!img) {
      par(mai = c(0, 0, 0, 0), ps = pointsize)
      pie(pie_area, labels = use_labels, col = use_color, clockwise = TRUE, ...)
    } else {
      check_package_installed("Cairo")
      out_file <- tempfile(pattern = "legendfile", fileext = ".png")
      Cairo::CairoPNG(file = out_file, bg = "white", width = width, height = height, pointsize = pointsize)
      par(mai = c(0, 0, 0, 0))
      pie(pie_area, labels = use_labels, col = use_color, clockwise = TRUE, ...)
      dev.off()
      base64_encode <- base64enc::dataURI(file = out_file)
      base64_encode <- sub("data:", "data:image/png", base64_encode, fixed = TRUE)
      cat(paste0('<img src="', base64_encode, '", width="200px">'))
    }
  }
}

base64_encode_images <- function(in_assign){
  if (Sys.info()['sysname'] == "Windows") {
    sub_pattern = "file:///"
  } else {
    sub_pattern = "file://localhost"
  }
  image_locs <- in_assign@pie_locs$piechart
  names(image_locs) = in_assign@pie_locs$visattr
  image_locs <- gsub(sub_pattern, "", image_locs)
  base_64_imgs <- vapply(image_locs, function(x){
    base64enc::dataURI(file = x)
  }, character(1))
  base_64_imgs
}


#' cc_graph to visnetwork
#' 
#' takes a \code{cc_graph} object and transforms it into something that can
#' be visualized using \code{visNetwork}
#' 
#' @param in_graph the cc_graph object
#' @param in_assign the colors generated by \code{assign_colors}
#' @param node_communities the communities generated by \code{label_communities}
#' @param use_nodes the list of nodes to actually use
#' 
#' @importFrom base64enc dataURI
#' @importFrom rlang .data
#' 
#' @export
#' @return list
graph_to_visnetwork <- function(in_graph, in_assign, node_communities = NULL, use_nodes = NULL){
  check_package_installed("DiagrammeR")
  in_graph <- add_tooltip(in_graph, description = in_assign@description, separator = "<br>")
  graph_nodes <- graph::nodes(in_graph)
  edge_list <- graph::edgeMatrix(in_graph)
  edge_weight <- unlist(graph::edgeData(in_graph, , , "weight"))
  
  if (is.null(use_nodes)) {
    use_nodes <- graph_nodes
  } else {
    use_nodes <- intersect(use_nodes, graph_nodes)
  }
  
  if (!is.null(node_communities)) {
    comm_nodes <- unique(unlist(lapply(node_communities, function(x){x$members})))
    use_nodes <- intersect(use_nodes, comm_nodes)
  }
  
  from_list <- graph_nodes[edge_list["from", ]]
  to_list <- graph_nodes[edge_list["to", ]]
  
  #in_nodes <- intersect(in_nodes, use_nodes)
  from_to <- data.frame(from = from_list, to = to_list)
  from_to <- dplyr::filter(from_to, (.data$from %in% use_nodes) & (.data$to %in% use_nodes))
  
  from_to$edgeid <- paste0(from_to$from, "|", from_to$to)
  from_to$weight <- edge_weight[from_to$edgeid]
  
  g_nodes <- DiagrammeR::create_node_df(n = length(use_nodes),
                            label = use_nodes,
                            group = in_assign@assignments[use_nodes]
  )
  
  if (in_assign@color_type == "pie") {
    g_nodes$shape = "image"
                           
    web_locs <- base64_encode_images(in_assign)
    
    for (igroup in g_nodes$group) {
      g_nodes[g_nodes$group == igroup, "image"] <- web_locs[igroup]
    }
  } else {
    g_nodes$shape = "circle"
    
    for (igroup in g_nodes$group) {
      g_nodes[g_nodes$group == igroup, "color"] <- in_assign@colors[igroup]
    }
  }
  
  for (inode in g_nodes$label) {
    which_node <- which(g_nodes$label %in% inode)
    g_nodes[which_node, "tooltip"] <- g_nodes[which_node, "title"] <- gsub("\\n", "<br>", graph::nodeData(in_graph, inode, "tooltip")[[1]])
  }
  
  
  if (!is.null(node_communities)) {
    g_nodes$community <- ""
    
    for (icomm in seq(1, length(node_communities))) {
      comm_data <- node_communities[[icomm]]
      g_nodes[(g_nodes$label %in% comm_data$members), "community"] <- comm_data$label
    }
  }
  
  
  from_to$from_id <- g_nodes$id[match(from_to$from, g_nodes$label)]
  from_to$to_id <- g_nodes$id[match(from_to$to, g_nodes$label)]
  
  g_edges <- DiagrammeR::create_edge_df(from = from_to$from_id,
                            to = from_to$to_id,
                            weight = from_to$weight)
  out_graph <- DiagrammeR::create_graph(g_nodes, g_edges, directed = FALSE)
  return(list(graph = out_graph, nodes = g_nodes, edges = g_edges))
}

#' assign communities
#' 
#' given a \code{cc_graph}, find communities of nodes based on their connectivity
#' and weights. 
#' 
#' @param in_graph the \code{cc_graph} object to use
#' @export
#' 
#' @return list
assign_communities <- function(in_graph){
  igraph_graph <- igraph::graph_from_graphnel(in_graph)
  walk_membership <- igraph::cluster_walktrap(igraph_graph)
  walk_communities <- igraph::membership(walk_membership)
  split_comms <- split(names(walk_communities), walk_communities)
  names(split_comms) <- NULL
  split_comms
}


#' label communities
#' 
#' Determine the label of a community based on the most generic member of
#' each community, which is defined as being the one with the most
#' annotations.
#' 
#' @param community_defs the communities from \code{assign_communities}
#' @param annotation the annotation object used for enrichment
#' 
#' @export
#' @return list
label_communities <- function(community_defs, annotation){
  n_members <- vapply(community_defs, length, numeric(1))
  
  community_defs <- community_defs[n_members > 1]
  
  all_members <- unique(unlist(community_defs))
  member_annotation_counts <- annotation@counts[all_members]
  
  get_rep_member <- lapply(community_defs, function(in_def){
    def_counts <- member_annotation_counts[names(member_annotation_counts) %in% in_def]
    max_member <- names(def_counts)[which.max(def_counts)]
    max_member[1]
  })
  
  community_info <- lapply(seq(1, length(community_defs)), function(i_def){
    if (length(annotation@description) != 0) {
      label <- annotation@description[[get_rep_member[[i_def]]]]
    } else {
      label <- get_rep_member[[i_def]]
    }
    list(rep = get_rep_member[[i_def]],
         label = label,
         members = community_defs[[i_def]])
  })
  
  community_info
}


#' vis in visNetwork
#' 
#' Visualize a \code{cc_graph} in \code{visNetwork}, with selection for communities
#' if that exists.
#' 
#' @param in_graph_info the graph structure from \code{graph_to_visnetwork}
#' 
#' @export
#' @return NULL
vis_visnetwork <- function(in_graph_info){
  check_package_installed("visNetwork")
  if (!is.null(in_graph_info$nodes$community)) {
    visNetwork::visOptions(visNetwork::visNetwork(edges = in_graph_info$edges,
                                                  nodes = in_graph_info$nodes),
                           selectedBy = "community")
  } else {
    visNetwork::visNetwork(edges = in_graph_info$edges,
                           nodes = in_graph_info$nodes)
  }
}

#' table from graph
#' 
#' Creates a table from the annotation graph, and if provided, adds the
#' community information to the table.
#' 
#' @param in_graph the \code{cc_graph} object
#' @param in_assign the \code{node_assign} object
#' @param community_info the \code{community_info} object
#' 
#' @export
#' @return data.frame
table_from_graph <- function(in_graph, in_assign = NULL, community_info = NULL){
  node_data <- graph::nodeData(in_graph, , )
  node_data <- lapply(node_data, as.data.frame, stringsAsFactors = FALSE)
  node_table <- do.call(rbind, node_data)
  
  if (!is.null(in_assign)) {
    sig_loc <- grep("sig$", names(node_table), value = TRUE)
    meas_loc <- grep("meas$", names(node_table), value = TRUE)
    
    node_table[, c(sig_loc, meas_loc)] <- NULL
    
    sig_group <- in_assign@description[node_table$name]
    
    if ("description" %in% names(node_table)) {
      cut_table_loc <- grep("description", names(node_table))
    } else {
      cut_table_loc <- grep("name", names(node_table))
    }
    
    tmp_table_1 <- cbind(node_table[, seq(1, cut_table_loc)], sig_group)
    node_table_2 <- cbind(tmp_table_1, node_table[, seq(cut_table_loc + 1, ncol(node_table))])
    
  } else {
    node_table_2 <- node_table
  }
  
  if (!is.null(community_info)) {
    in_members <- unique(unlist(lapply(community_info, function(x){x$members})))
    out_members <- setdiff(node_table_2$name, in_members)
    n_comm <- length(community_info)
    community_info[[n_comm + 1]] <- list(label = "other", members = out_members)
    
    community_info = purrr::map(seq(1, length(community_info)), function(comm_id){
      community_info[[comm_id]]$group = comm_id
      community_info[[comm_id]]
    })
    
    null_table <- node_table_2[1, ]
    rownames(null_table) <- NULL
    null_table <- lapply(null_table, function(x){
      if (is.character(x)) {
        out_value <- ""
      } else {
        out_value <- as.numeric(NA)
      }
      out_value
    })
    null_table <- as.data.frame(null_table, stringsAsFactors = FALSE)
    
    node_table_split <- lapply(community_info, function(in_info){
      header_table <- null_table
      if ("description" %in% names(node_table_2)) {
        header_table$description <- paste0("**", in_info$label, "**")
      } else {
        header_table$name <- in_info$label
      }
      tmp_table <- node_table_2[(node_table_2$name %in% in_info$members), ]
      sort_cols <- character(0)
      if ("sig_group" %in% names(tmp_table)) {
        sort_cols <- c(sort_cols, "sig_group")
      } else if (length(grep("padjust$", names(tmp_table))) != 0) {
        padjust_cols <- grep("padjust$", names(tmp_table), value = TRUE)
        sort_cols <- c(sort_cols, padjust_cols)
      }
      
      if (length(sort_cols) > 0) {
        tmp_table <- dplyr::arrange(tmp_table, .data[[sort_cols]])
      }
      
      out_table <- rbind(header_table, tmp_table)
      rownames(out_table) <- NULL
      out_table$group <- in_info$group
      out_table
    })
    out_node_table <- do.call(rbind, node_table_split)
  } else {
    out_node_table <- node_table_2
  }
  out_node_table
}

