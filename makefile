all: ../DNAPlotR_0.1.tar.gz

.PHONY: all install

install:
	R -e 'devtools::install_github("sherrillmix/DNAPlotR")'

man: R/*.R
	R -e 'devtools::document()'
	touch man


inst/doc: vignettes/*.Rnw
	R -e 'devtools::build_vignettes()'
	touch inst/doc

	
../DNAPlotR_0.1.tar.gz: man R/*.R DESCRIPTION inst/doc
	R -e 'devtools::check();devtools::build()'
