ggarrange <- function (..., plotlist = NULL, ncol = NULL, nrow = NULL, labels = NULL, 
					label.x = 0, label.y = 1, hjust = -0.5, vjust = 1.5, font.label = list(size = 14, 
					color = "black", face = "bold", family = NULL), align = c("none", "h", "v", "hv"), 
					widths = 1, heights = 1, legend = NULL, common.legend = FALSE, plot_legend = NULL) 
{
	plots <- c(list(...), plotlist)
	align <- match.arg(align)
	nb.plots <- length(plots)
	nb.plots.per.page <- .nbplots_per_page(ncol, nrow)
	if (is.null(legend) & common.legend) 
		legend <- "top"
	legend <- .check_legend(legend)
	if (!is.null(legend)) 
		plots <- purrr::map(plots, function(x) {
			if (!is.null(x)) 
				x + theme(legend.position = legend)
			else x
		})
	leg <- NULL
	if (common.legend) {
		if (!is.null(plot_legend) & plot_legend <= nb.plots)
		leg <- get_legend(plots[plot_legend])
		plots <- purrr::map(plots, function(x) {
			if (!is.null(x)) 
				x + theme(legend.position = "none")
			else x
		})
	}
	if (nb.plots > nb.plots.per.page) {
		plots <- split(plots, ceiling(seq_along(plots)/nb.plots.per.page))
	}
	else plots <- list(plots)
	.lab <- .update_label_pms(font.label, label.x = label.x, 
														label.y = label.y, hjust = hjust, vjust = vjust)
	res <- purrr::map(plots, .plot_grid, ncol = ncol, nrow = nrow, 
										labels = labels, label_size = .lab$size, label_fontfamily = .lab$family, 
										label_fontface = .lab$face, label_colour = .lab$color, 
										label_x = .lab$label.x, label_y = .lab$label.y, hjust = .lab$hjust, 
										vjust = .lab$vjust, align = align, rel_widths = widths, 
										rel_heights = heights, legend = legend, common.legend.grob = leg)
	if (length(res) == 1) 
		res <- res[[1]]
	class(res) <- c(class(res), "ggarrange")
	res
}

