log_to_file <- function(message) {
  write(message, file = "debug_log.txt", append = TRUE)
}

#'Read the lrobject and generate the comparative tables
#'
#'@param data LRObj with single condition
#'@param out_path output path
#'@return LRObject
#'@import dplyr
#'@importFrom tidyr %>%
create_diff_table1 <- function(data, out_path, comparison = NULL) {
  if (!is.null(comparison)) {
    for (pair in comparison) {
      log_to_file(paste("Processing comparison: ", pair[1], " vs ", pair[2]))
      ctr_name <- pair[2]
      ctr_table <- data@tables[[ctr_name]]
      exp_name <- pair[1]
      cmp_name <- paste0(exp_name, "_x_", ctr_name)
      exp_table <- data@tables[[exp_name]]
      tmp_data <- merge(exp_table, ctr_table, by = 'allpair', all = TRUE)
      
      log_to_file(paste("Columns after merge for", cmp_name, ":", paste(colnames(tmp_data), collapse = ", ")))
      
      tmp_data <- tmp_data %>%
        tidyr::separate(.data$allpair, c('ligpair', 'recpair'), sep = '@', remove = F)
      tmp_data$LRScore.x[is.na(tmp_data$LRScore.x)] <- 0
      tmp_data$LRScore.y[is.na(tmp_data$LRScore.y)] <- 0
      
      final_data <- tmp_data %>%
        dplyr::mutate(LRScore = .data$LRScore.x - .data$LRScore.y) %>%
        dplyr::select(c(.data$source.x, .data$source.y, .data$target.x, .data$target.y, .data$LRScore,
                        .data$ligpair, .data$recpair, .data$allpair, .data$type_gene_A.x, .data$type_gene_B.x,
                        .data$type_gene_A.y, .data$type_gene_B.y, .data$gene_A.x, .data$gene_B.x,
                        .data$gene_A.y, .data$gene_B.y)) %>%
        dplyr::arrange(LRScore)
      
      final_data <- final_data[final_data$LRScore != 0,]
      final_data <- final_data %>%
        dplyr::mutate(type_gene_A = dplyr::coalesce(.data$type_gene_A.x, .data$type_gene_A.y)) %>%
        dplyr::mutate(type_gene_B = dplyr::coalesce(.data$type_gene_B.x, .data$type_gene_B.y)) %>%
        dplyr::mutate(gene_A = dplyr::coalesce(.data$gene_A.x, .data$gene_A.y)) %>%
        dplyr::mutate(gene_B = dplyr::coalesce(.data$gene_B.x, .data$gene_B.y)) %>%
        dplyr::mutate(source = dplyr::coalesce(.data$source.x, .data$source.y)) %>%
        dplyr::mutate(target = dplyr::coalesce(.data$target.x, .data$target.y)) %>%
        dplyr::select(-.data$type_gene_A.x, -.data$type_gene_A.y, -.data$type_gene_B.x, -.data$type_gene_B.y,
                      -.data$gene_A.x, -.data$gene_A.y, -.data$gene_B.x, -.data$gene_B.y,
                      -.data$source.x, -.data$source.y, -.data$target.x, -.data$target.y) %>%
        dplyr::mutate(cellpair = paste0(.data$source, "@", .data$target))
      
      log_to_file(paste("Columns in final_data for", cmp_name, ":", paste(colnames(final_data), collapse = ", ")))
      
      final_data$interaction_type <- paste(final_data$type_gene_A, final_data$type_gene_B, sep = "")
      final_data <- final_data %>%
        dplyr::mutate(interaction_type = ifelse(interaction_type == "LigandReceptor", "LR", interaction_type)) %>%
        dplyr::mutate(interaction_type = ifelse(interaction_type == "ReceptorTranscription Factor", "RTF", interaction_type)) %>%
        dplyr::mutate(interaction_type = ifelse(interaction_type == "Transcription FactorLigand", "TFL", interaction_type))
      data@tables[[cmp_name]] <- final_data
      
      final <- final_data %>%
        dplyr::mutate(ccitype = paste(.data$type_gene_A, .data$type_gene_B)) %>%
        dplyr::filter(!(stringr::str_detect(.data$ccitype, "Transcription Factor"))) %>%
        dplyr::group_by(.data$cellpair) %>%
        dplyr::summarise(LRScore = sum(.data$LRScore))
      
      final <- final %>%
        tidyr::separate(.data$cellpair, c("u", "v"), sep = "@", remove = F)
      final <- na.omit(final)
      
      filtervar <- grepl('Transcription', final_data[['type_gene_A']]) | grepl('Transcription', final_data[['type_gene_B']])
      raw_inter <- table(final_data$cellpair[!filtervar])
      freq <- as.array(raw_inter)[final$cellpair] - min(as.array(raw_inter)[final$cellpair])
      freq <- freq / (max(as.array(raw_inter)[final$cellpair]) - min(as.array(raw_inter)[final$cellpair])) + 0.1
      final$freq <- as.array(freq)[final$cellpair]
      final$pair <- final$cellpair
      final <- dplyr::arrange(final, LRScore)
      
      graph1 <- igraph::graph_from_data_frame(final[, c("u", "v", "LRScore")])
      igraph::E(graph1)$inter <- final$freq
      igraph::E(graph1)$inter.raw <- as.array(raw_inter)[final$cellpair]
      igraph::E(graph1)$weight <- igraph::E(graph1)$LRScore
      igraph::E(graph1)$mean <- igraph::E(graph1)$LRScore
      data@graphs[[cmp_name]] <- graph1
      
      graph1 <- igraph::graph_from_data_frame(final_data[, c("ligpair", "recpair", "LRScore")])
      igraph::E(graph1)$weight <- igraph::E(graph1)$LRScore
      igraph::E(graph1)$mean <- igraph::E(graph1)$LRScore
      data@graphs_ggi[[cmp_name]] <- graph1
    }
  } else {
    ctr_name <- names(data@tables)[1]
    ctr_table <- data@tables[[ctr_name]]
    for (i in 2:length(data@tables)) {
      log_to_file(paste("Processing comparison: ", names(data@tables)[i]))
      exp_name <- names(data@tables)[i]
      cmp_name <- paste0(exp_name, "_x_", ctr_name)
      exp_table <- data@tables[[exp_name]]
      tmp_data <- merge(exp_table, ctr_table, by = 'allpair', all = TRUE)
      
      log_to_file(paste("Columns after merge for", cmp_name, ":", paste(colnames(tmp_data), collapse = ", ")))
      
      tmp_data <- tmp_data %>%
        tidyr::separate(.data$allpair, c('ligpair', 'recpair'), sep = '@', remove = F)
      tmp_data$LRScore.x[is.na(tmp_data$LRScore.x)] <- 0
      tmp_data$LRScore.y[is.na(tmp_data$LRScore.y)] <- 0
      
      final_data <- tmp_data %>%
        dplyr::mutate(LRScore = .data$LRScore.x - .data$LRScore.y) %>%
        dplyr::select(c(.data$source.x, .data$source.y, .data$target.x, .data$target.y, .data$LRScore,
                        .data$ligpair, .data$recpair, .data$allpair, .data$type_gene_A.x, .data$type_gene_B.x,
                        .data$type_gene_A.y, .data$type_gene_B.y, .data$gene_A.x, .data$gene_B.x,
                        .data$gene_A.y, .data$gene_B.y)) %>%
        dplyr::arrange(LRScore)
      
      final_data <- final_data[final_data$LRScore != 0,]
      final_data <- final_data %>%
        dplyr::mutate(type_gene_A = dplyr::coalesce(.data$type_gene_A.x, .data$type_gene_A.y)) %>%
        dplyr::mutate(type_gene_B = dplyr::coalesce(.data$type_gene_B.x, .data$type_gene_B.y)) %>%
        dplyr::mutate(gene_A = dplyr::coalesce(.data$gene_A.x, .data$gene_A.y)) %>%
        dplyr::mutate(gene_B = dplyr::coalesce(.data$gene_B.x, .data$gene_B.y)) %>%
        dplyr::mutate(source = dplyr::coalesce(.data$source.x, .data$source.y)) %>%
        dplyr::mutate(target = dplyr::coalesce(.data$target.x, .data$target.y)) %>%
        dplyr::select(-.data$type_gene_A.x, -.data$type_gene_A.y, -.data$type_gene_B.x, -.data$type_gene_B.y,
                      -.data$gene_A.x, -.data$gene_A.y, -.data$gene_B.x, -.data$gene_B.y,
                      -.data$source.x, -.data$source.y, -.data$target.x, -.data$target.y) %>%
        dplyr::mutate(cellpair = paste0(.data$source, "@", .data$target)) %>%
        dplyr::arrange(LRScore)
      
      log_to_file(paste("Columns in final_data for", cmp_name, ":", paste(colnames(final_data), collapse = ", ")))
      
      final_data$interaction_type <- paste(final_data$type_gene_A, final_data$type_gene_B, sep = "")
      final_data <- final_data %>%
        dplyr::mutate(interaction_type = ifelse(interaction_type == "LigandReceptor", "LR", interaction_type)) %>%
        dplyr::mutate(interaction_type = ifelse(interaction_type == "ReceptorTranscription Factor", "RTF", interaction_type)) %>%
        dplyr::mutate(interaction_type = ifelse(interaction_type == "Transcription FactorLigand", "TFL", interaction_type))
      data@tables[[cmp_name]] <- final_data
      
      final <- final_data %>%
        dplyr::mutate(ccitype = paste(.data$type_gene_A, .data$type_gene_B)) %>%
        dplyr::filter(!(stringr::str_detect(.data$ccitype, "Transcription Factor"))) %>%
        dplyr::group_by(.data$cellpair) %>%
        dplyr::summarise(LRScore = sum(.data$LRScore))
      
      final <- final %>%
        tidyr::separate(.data$cellpair, c("u", "v"), sep = "@", remove = F)
      final <- na.omit(final)
      
      filtervar <- grepl('Transcription', final_data[['type_gene_A']]) | grepl('Transcription', final_data[['type_gene_B']])
      raw_inter <- table(final_data$cellpair[!filtervar])
      freq <- as.array(raw_inter)[final$cellpair] - min(as.array(raw_inter)[final$cellpair])
      freq <- freq / (max(as.array(raw_inter)[final$cellpair]) - min(as.array(raw_inter)[final$cellpair])) + 0.1
      final$freq <- as.array(freq)[final$cellpair]
      final$pair <- final$cellpair
      final <- dplyr::arrange(final, LRScore)
      
      graph1 <- igraph::graph_from_data_frame(final[, c("u", "v", "LRScore")])
      igraph::E(graph1)$inter <- final$freq
      igraph::E(graph1)$inter.raw <- as.array(raw_inter)[final$cellpair]
      igraph::E(graph1)$weight <- igraph::E(graph1)$LRScore
      igraph::E(graph1)$mean <- igraph::E(graph1)$LRScore
      data@graphs[[cmp_name]] <- graph1
      
      graph1 <- igraph::graph_from_data_frame(final_data[, c("ligpair", "recpair", "LRScore")])
      igraph::E(graph1)$weight <- igraph::E(graph1)$LRScore
      igraph::E(graph1)$mean <- igraph::E(graph1)$LRScore
      data@graphs_ggi[[cmp_name]] <- graph1
    }
  }
  saveRDS(data, file.path(out_path, "LR_data_step2.Rds"))
  return(data)
}
