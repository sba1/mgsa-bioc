#!/bin/bash
#
# Builds the manual
#

R CMD roxygen -d `dirname $0`/..
