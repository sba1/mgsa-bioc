#!/bin/bash
#
# Builds the manual
#

R CMD roxygen -d `dirname $0`/..

R CMD Rd2dvi --pdf -o mgsa.pdf `dirname $0`/..
