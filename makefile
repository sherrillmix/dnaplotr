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

data: data-raw/makeAmpliciationLookup.R
	Rscript data-raw/makeAminoColors.R
	touch data
	
../DNAPlotR_0.1.tar.gz: man R/*.R DESCRIPTION inst/doc
	sed -i "s/^Date:.*$$/Date: `date +%Y-%m-%d`/" DESCRIPTION
	R -e 'devtools::check();devtools::build()'
