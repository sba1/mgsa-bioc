#!/bin/bash
#
# Builds the manual
#

rm -f `dirname $0`/../man/*.Rd

LANG="en.ASCII" R CMD roxygen -d -s `dirname $0`/..

rm -Rf mgsa.pdf mgsa-internal.pdf
R CMD Rd2dvi --no-preview --pdf -o mgsa.pdf `dirname $0`/..
R CMD Rd2dvi --internals --pdf -o mgsa-internal.pdf `dirname $0`/..
