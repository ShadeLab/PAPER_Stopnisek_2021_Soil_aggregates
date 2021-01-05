submission/manuscript.pdf: submission/manuscript.Rmd
	R -e 'library(markdown); render("submission/manuscript.Rmd")'
