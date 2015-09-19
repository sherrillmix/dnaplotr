all: docs package

install:
	R -e 'devtools::install_github("sherrillmix/DNAPlotR")'

docs: R/*.R
	R -e 'devtools::document()'
	R -e 'devtools::build_vignettes()'


package: docs R/*.R DESCRIPTION
	R -e 'devtools::check();devtools::build()'
