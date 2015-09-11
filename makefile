all: docs package

install:
	R -e 'devtools::install_github("sherrillmix/plotDNA")'

docs: R/*.R
	R -e 'devtools::document()'

package: docs R/*.R DESCRIPTION
	R -e 'devtools::check();devtools::build()'
