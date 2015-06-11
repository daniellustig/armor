#!/usr/bin/env python

################################################################################
# one.py                                                                       #
#                                                                              #
# Script for generating a single ArOR DBT shim FSMs                            #
################################################################################
# The ArMOR project                                                            #
#                                                                              #
# Copyright (C) 2015 Daniel Lustig, Princeton University                       #
#                                                                              #
# This software may be modified and distributed under the terms                #
# of the MIT license.  See the LICENSE file for details.                       #
################################################################################

import sys
import argparse
import logging

# Argument parsing

parser = argparse.ArgumentParser(description="Generate ArMOR Shim FSM Summary")
parser.add_argument('-o', dest="graphfilepath", type=str, required=True,
        help=("Directory in which to store output files."))
parser.add_argument('-l', dest="logfile", type=str, default=None,
        help="Print the log to a file")
parser.add_argument('-v', dest="verbose", required=False,
        action="store_true", help="Verbose")
args = parser.parse_args()

# create logger
logger = logging.getLogger()
logger.setLevel(logging.NOTSET)

formatter = logging.Formatter('%(levelname)08s - %(filename)s:%(lineno)s - %(message)s')

# create file handler
if args.logfile:
    fh = logging.FileHandler(args.logfile, 'w')
    if args.verbose:
        fh.setLevel(logging.DEBUG)
    else:
        fh.setLevel(logging.INFO)
    fh.setFormatter(formatter)
    logger.addHandler(fh)

# create console handler
ch = logging.StreamHandler()
ch.setLevel(logging.WARNING)
ch.setFormatter(formatter)
logger.addHandler(ch)

# This needs to be imported *after* the loggers are set up
import armor

################################################################################
# Usage examples                                                               #
################################################################################

# Simplified version of SC upstream and TSO downstream

# Define the MOSTs relevant to the architectures in question
sc_ppo = armor.MOST({
    "ld": {"ld": 'S', "st": 'S'},
    "st": {"ld": 'S', "st": 'S'}})
tso_ppo = armor.MOST({
    "ld": {"ld": 'S', "st": 'S'},
    "st": {"ld": '-', "st": 'M'}})

# Define the source and target architectures in terms of the above MOSTs
sc  = armor.Architecture("sc2x2" , sc_ppo , {},                 ["ld", "st"])
tso = armor.Architecture("tso2x2", tso_ppo, {"mfence": sc_ppo}, ["ld", "st", "mfence"])

# Example 1: SC upstream, TSO downstream (fences need to be inserted)
armor.GenerateFSMAndPrintDOTGraph(args.graphfilepath, sc, tso)

# Example 2: TSO upstream, SC downstream (fences may be removed)
armor.GenerateFSMAndPrintDOTGraph(args.graphfilepath, tso, sc)

################################################################################

# Versions of SC upstream and TSO downstream addressing same vs. different
# addresses

sc_ppo = armor.MOST({
     "ld": {"ld&sa": 'S', "ld&da": 'S', "st": 'S'},
     "st": {"ld&sa": 'S', "ld&da": 'S', "st": 'S'}})
tso_ppo = armor.MOST({
     "ld": {"ld&sa": 'S', "ld&da": 'S', "st": 'S'},
     "st": {"ld&sa": 'L', "ld&da": '-', "st": 'M'}})

# Define the source and target architectures in terms of the above MOSTs
sc  = armor.Architecture("sc2x3" , sc_ppo , {},                 ["ld", "st"])
tso = armor.Architecture("tso2x3", tso_ppo, {"mfence": sc_ppo}, ["ld", "st", "mfence"])

# Example 1: SC upstream, TSO downstream (fences need to be inserted)
armor.GenerateFSMAndPrintDOTGraph(args.graphfilepath, sc, tso)

# Example 2: TSO upstream, SC downstream (fences may be removed)
armor.GenerateFSMAndPrintDOTGraph(args.graphfilepath, tso, sc)

################################################################################

# Versions of SC upstream and TSO downstream addressing cumulativity

sc_ppo = armor.MOST({
     "ld": {"ld": 'S', "cld": 'S', "st": 'S', "cst": 'S'},
    "cld": {"ld": 'S', "cld": 'S', "st": 'S', "cst": 'S'},
     "st": {"ld": 'S', "cld": 'S', "st": 'S', "cst": 'S'},
    "cst": {"ld": 'S', "cld": 'S', "st": 'S', "cst": 'S'}})
tso_ppo = armor.MOST({
     "ld": {"ld": 'S', "cld": 'S', "st": 'S', "cst": 'S'},
    "cld": {"ld": 'S', "cld": 'S', "st": 'S', "cst": 'S'},
     "st": {"ld": '-', "cld": 'S', "st": 'M', "cst": 'S'},
    "cst": {"ld": 'S', "cld": 'S', "st": 'S', "cst": 'S'}})

# Define the source and target architectures in terms of the above MOSTs
sc  = armor.Architecture("sc4x4" , sc_ppo , {},                 ["ld", "st"])
tso = armor.Architecture("tso4x4", tso_ppo, {"mfence": sc_ppo}, ["ld", "st", "mfence"])

# Example 1: SC upstream, TSO downstream (fences need to be inserted)
armor.GenerateFSMAndPrintDOTGraph(args.graphfilepath, sc, tso)

# Example 2: TSO upstream, SC downstream (fences may be removed)
armor.GenerateFSMAndPrintDOTGraph(args.graphfilepath, tso, sc)

################################################################################

# For more examples, see:
# - architectures.py: defines a wide range of architectures, including those
#   analyzed in the ISCA paper
# - report.py: generates a report summarizing the shim FSMs for a variety of
#   use cases, including all of those analyzed in the ISCA paper
