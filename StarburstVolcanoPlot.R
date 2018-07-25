## This program plot TCHAvisualize_volcano plot
### Run command
#volocavoplot(x, y, names = rownames(genes))
#volocavoplot(x, y, names = rownames(genes), x.cut = 2, y.cut = 0.01, y1.cut=1, y1.cut=0.01)
## x1.cut and y1.cut is cutoff for plot
#x1.cut and y1.cut is value for dotted line on X and Y axis
genes <- read.table("results.txt", header=TRUE)
rownames(genes) <- genes$Gene
x <- genes$log2FoldChange; y <- genes$padj
############################
volocavoplot <- function (x, y, filename = "volcano1.pdf", ylab = expression(paste(-Log[10], " (FDR corrected -P values)")), 
          xlab = NULL, title = NULL, legend = NULL, label = NULL, xlim = NULL, 
          ylim = NULL, color = c("black", "green4", "red"), names = NULL, names.fill = TRUE, 
          show.names = "significant", x.cut = 1, y.cut = 0.01, x1.cut=1, y1.cut=0.05, height = 10, width = 10, 
          highlight = NULL, highlight.color = "orange", names.size = 3, dpi = 500) 
{
  if (!is.null(names)) {
    if (all(grepl("\\|", names))) {
      names <- strsplit(names, "\\|")
      names <- unlist(lapply(names, function(x) x[1]))
    }
  }
  .e <- environment()
  threshold <- rep("1", length(x))
  names(color) <- as.character(1:3)
  if (is.null(label)) {
    label = c(`1` = "Not Significant", `2` = "Up regulared", 
              `3` = "Down regulated")
  }
  else {
    names(label) <- as.character(1:3)
  }
  sig <- y < y.cut
  sig[is.na(sig)] <- FALSE
  up <- x > x.cut
  up[is.na(up)] <- FALSE
  if (any(up & sig)) 
    threshold[up & sig] <- "2"
  down <- x < (-x.cut)
  down[is.na(down)] <- FALSE
  if (any(down & sig)) 
    threshold[down & sig] <- "3"
  if (!is.null(highlight)) {
    idx <- which(names %in% highlight)
    if (length(idx) > 0) {
      print(idx)
      threshold[which(names %in% highlight)] <- "4"
      color <- c(color, highlight.color)
      names(color) <- as.character(1:4)
    }
  }
  df <- data.frame(x = x, y = y, threshold = threshold)
  if (!is.null(highlight)) {
    order.idx <- order(df$threshold)
    down <- down[order.idx]
    sig <- sig[order.idx]
    up <- up[order.idx]
    df <- df[order.idx, ]
    names <- names[order.idx]
  }
  p <- ggplot(data = df, aes(x = x, y = -1 * log10(y), colour = threshold), 
              environment = .e) + geom_point() + ggtitle(title) + ylab(ylab) + 
    xlab(xlab) + geom_vline(aes(xintercept = -x1.cut), colour = "black", 
                            linetype = "dashed") + geom_vline(aes(xintercept = x1.cut), 
                                                              colour = "black", linetype = "dashed") + geom_hline(aes(yintercept = -1 * 
                                                                                                                        log10(y1.cut)), colour = "black", linetype = "dashed") + 
    scale_color_manual(breaks = as.numeric(names(label)), 
                       values = color, labels = label, name = legend) + 
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), legend.text = element_text(size = 10), 
                       axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), 
                       legend.position = "top", legend.key = element_rect(colour = "white"))
  if (!is.null(names)) {
    if (show.names == "significant") {
      idx <- (up & sig) | (down & sig)
      important <- c("2", "3")
    }
    else if (show.names == "highlighted") {
      if (!is.null(highlight)) {
        idx <- (names %in% highlight)
        important <- c("4")
      }
      else {
        message("Missing highlight argument")
        return(NULL)
      }
    }
    else if (show.names == "both") {
      if (!is.null(highlight)) {
        idx <- (up & sig) | (down & sig) | (names %in% 
                                              highlight)
        important <- c("2", "3", "4")
      }
      else {
        message("Missing highlight argument")
        return(NULL)
      }
    }
    else {
      message("Wrong highlight argument")
      return(NULL)
    }
    if (any(threshold %in% important)) {
      if (names.fill) {
        p <- p + geom_label_repel(data = subset(df, threshold %in% 
                                                  important), aes(label = names[idx], fill = threshold), 
                                  size = names.size, show.legend = FALSE, fontface = "bold", 
                                  color = "black", box.padding = unit(0.35, "lines"), 
                                  point.padding = unit(0.3, "lines")) + scale_fill_manual(values = color[as.numeric(important)])
      }
      else {
        p <- p + geom_text_repel(data = subset(df, threshold %in% 
                                                 important), aes(label = names[idx]), size = names.size, 
                                 show.legend = FALSE, fontface = "bold", color = "black", 
                                 point.padding = unit(0.3, "lines"), box.padding = unit(0.5, 
                                                                                        "lines"))
      }
    }
  }
  if (!is.null(filename)) {
    ggsave(p, filename = filename, width = width, height = height, 
           dpi = dpi)
  }
  else {
    return(p)
  }
}