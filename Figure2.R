# libraries
library(ggplot2)
library(patchwork)
library(ggpubr)
library(scales)
library(ggdist)
library(dplyr)
library(cowplot)
library(pdftools)
library(png)
library(ggplot2)
library(pdftools)
library(png)
library(grid)
library(grafify)
library(ggsignif)
library(scales)
## load data for Figure 2 B. 

dataset_melt <- readr::read_tsv("F2B.tsv") %>% as.data.frame()
dataset_melt$value <- as.numeric(dataset_melt$value)
dataset_melt <- dataset_melt[!is.na(dataset_melt$value),]
data_p1 <- dataset_melt %>% dplyr::filter(grepl("mut[.]",variable)) %>%
  dplyr::mutate(germline = gsub("full","Full", germline),
                germline = gsub("subset","Subset", germline))


p1 <- ggplot(data_p1, 
             aes(germline, value, fill = germline)) + 
  ## add half-violin from {ggdist} package
  ggdist::stat_halfeye( 
    ## custom bandwidth
    adjust = .5, 
    ## adjust height
    width = .6, 
    ## move geom to the right
    justification = -.2, 
    ## remove slab interval
    .width = 0, 
    point_colour = NA
  ) + 
  geom_boxplot(
    width = .12, 
    ## remove outliers
    outlier.color = NA ## `outlier.shape = NA` works as well
  ) + ## remove white space on the left
  coord_cartesian(xlim = c(1.2, NA)) +
  ggpubr::theme_pubclean(base_size = 30)  + guides(fill="none") + 
  labs(y="# of Mutations",x="IGHV Germline reference set") +
  theme(axis.title = element_text(size = 32)) +
  scale_fill_manual(values = c("#E69F00","#56B4E9")) +
  geom_signif(comparisons = list(c("Full", "Subset")), 
              map_signif_level=TRUE)


# load data for Figure 2 A
data_maskprimers <- readr::read_csv("F2A.csv")

number2 <- function (x, accuracy = NULL, scale = 1, prefix = "", suffix = "", 
                     big.mark = " ", decimal.mark = ".", style_positive = c("none", 
                                                                            "plus"), style_negative = c("hyphen", "minus", "parens"), 
                     scale_cut = NULL, trim = TRUE, ...) 
{
  if (length(x) == 0) {
    return(character())
  }
  style_positive <- rlang::arg_match(style_positive)
  style_negative <- rlang::arg_match(style_negative)
  if (!is.null(scale_cut)) {
    cut <- scales:::scale_cut(x, breaks = scale_cut, scale = scale, 
                     accuracy = accuracy, suffix = suffix)
    scale <- cut$scale
    suffix <- cut$suffix
    accuracy <- cut$accuracy
  }
  accuracy <- scales:::`%||%`(accuracy,scales:::precision(x * scale))
  
  x <- scales:::round_any(x, accuracy/scale)
  nsmalls <- -floor(log10(accuracy))
  nsmalls <- pmin(pmax(nsmalls, 0), 20)
  sign <- sign(x)
  sign[is.na(sign)] <- 0
  x <- abs(x)
  x_scaled <- scale * x
  ret <- character(length(x))
  suffix_ <- sapply(x_scaled, function(x) if(x==0 & !is.na(x)) "" else suffix)
  for (nsmall in unique(nsmalls)) {
    idx <- nsmall == nsmalls
    ret[idx] <- format(x_scaled[idx], big.mark = big.mark, 
                       decimal.mark = decimal.mark, trim = trim, nsmall = nsmall, 
                       scientific = FALSE, ...)
  }
  
  ret <- paste0(prefix, ret, suffix_)
  ret[is.infinite(x)] <- as.character(x[is.infinite(x)])
  if (style_negative == "hyphen") {
    ret[sign < 0] <- paste0("-", ret[sign < 0])
  }
  else if (style_negative == "minus") {
    ret[sign < 0] <- paste0("−", ret[sign < 0])
  }
  else if (style_negative == "parens") {
    ret[sign < 0] <- paste0("(", ret[sign < 0], ")")
  }
  if (style_positive == "plus") {
    ret[sign > 0] <- paste0("+", ret[sign > 0])
  }
  ret[is.na(x)] <- NA
  names(ret) <- names(x)
  ret
}


label_number2 <- function (accuracy = NULL, scale = 1, prefix = "", suffix = "", 
          big.mark = " ", decimal.mark = ".", style_positive = c("none", 
                                                                 "plus"), style_negative = c("hyphen", "minus", "parens"), 
          scale_cut = NULL, trim = TRUE, ...) 
{
  scales:::force_all(accuracy, scale, prefix, suffix, big.mark, decimal.mark, 
            style_positive, style_negative, scale_cut, trim, ...)
  function(x) {
    number2(x, accuracy = accuracy, scale = scale, prefix = prefix, 
           suffix = suffix, big.mark = big.mark, decimal.mark = decimal.mark, 
           style_positive = style_positive, style_negative = style_negative, 
           scale_cut = scale_cut, trim = trim, ...)
  }
}

# vanhidden -  pipeline3 - SRR4026043
# stern - pipeline1 - SRR1383456
# yaari2 - pipeline6 - SRR19445534	


p3<-ggplot(data_maskprimers %>%
             mutate(pipeline = recode(pipeline, "pipeline1" = "SRR4026043", "pipeline3" = "SRR1383456", "pipeline6" = "SRR19445534")), 
           aes(x=as.factor(error), 
               y=num_pass_read, 
               fill=as.factor(error))) + 
  geom_bar(stat = "identity") +
  facet_wrap(.~pipeline, scales = "free") +
  guides(fill = FALSE) +
  theme(strip.text.x = element_blank()) +
  scale_y_continuous(labels = label_number2(suffix = "K", scale = 1e-3)) +
  labs(y= "# of passed reads", x = "Error rate threshold") +
  ggpubr::theme_pubclean(base_size = 30) +
  theme(axis.title = element_text(size = 32), 
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  scale_fill_manual(values = c("#56B4E9","#E69F00","#56B4E9"))


v <- read.csv("F2CA.csv",header =TRUE);v$X=NULL
v_my <- read.csv("F2CB.csv",header =TRUE);v_my$X=NULL

chron_v_results <- reshape2::melt(v[,seq(1, ncol(v), 2) ])
chron_v_results$variable <- gsub("[.]","-", chron_v_results$variable)
chron_v_results$method <- "control"

chron_v_results_dfn <- reshape2::melt(v_my[,seq(1, ncol(v_my), 2) ])
chron_v_results_dfn$variable <- gsub("[.]","-", chron_v_results_dfn$variable)
chron_v_results_dfn$method <- "control-DolphinNext"

plot_data <- rbind(chron_v_results, chron_v_results_dfn)


p4<- ggplot(plot_data, aes(x=as.factor(variable), y=value, fill=method)) + 
  geom_boxplot() +
  ggpubr::theme_pubclean(base_size = 30) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
        legend.position=c(.9,0.95), axis.title = element_text(size = 32)) +
  labs(x = "V gene", y = "Frequency", fill = "") +
  scale_fill_manual(values = c("#61D04F","#DF536B"))
  





library(gridExtra)

p1 <- p1 + 
  labs(tag = "B")

p3 <- p3 + 
  labs(tag = "A")

pp <- p4 + 
  labs(tag = "C")

row1 <- cowplot::plot_grid(p3, p1, ncol = 2, rel_widths = c(0.6,0.4))#grid.arrange(p3, p1, ncol=2)
row2 <- pp

pdf("Figure2.pdf",width = 20,height = 20)
grid.arrange(row1, row2, nrow=2)
dev.off()

