#!/bin/bash
#
# Builds the manual
#

rm -f `dirname $0`/../man/*.Rd

LANG="en.ASCII" echo 'library(roxygen2); roxygenize("..")' | R --vanilla

rm -Rf mgsa.pdf mgsa-internal.pdf
R CMD Rd2pdf --no-preview --pdf -o mgsa.pdf `dirname $0`/..
R CMD Rd2pdf --internals --pdf -o mgsa-internal.pdf `dirname $0`/..
