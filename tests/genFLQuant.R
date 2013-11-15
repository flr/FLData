#=====================================================================
#
# Date: 30/08/2012
# Version: 0.1-0
# Authors: Ernesto Jardim & Colin Millar
#
# Short description: tests for genFLQuant
#
# ToDo:
#
# References (bibtex):
#
#!Notes:
#
#=====================================================================

library(FLCore)
data(ple4)

# start test
setCon()
zz <- startTest("genFLQuantTests.txt")
tagTest("genFLQuant testing ...")

flq <- harvest(ple4)
checkRun(flqSim <- genFLQuant(harvest(ple4)))
checkTrue(is(flqSim, "FLQuant"))
checkIdentical(units(flqSim), units(flq))

finishTest()

