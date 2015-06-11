#!/usr/bin/env python

################################################################################
#                                                                              #
# architectures.py                                                             #
#                                                                              #
# A set of pre-defined MOSTs and architecture definitions                      #
#                                                                              #
################################################################################
#                                                                              #
# The ArMOR project                                                            #
#                                                                              #
# Copyright (C) 2015 Daniel Lustig, Princeton University                       #
#                                                                              #
# This software may be modified and distributed under the terms                #
# of the MIT license.  See the LICENSE file for details.                       #
#                                                                              #
################################################################################

import sys
import copy
import unittest
import logging
from armor import MOST, Architecture

################################################################################
# Baseline Model Definitions                                                   #
#                                                                              #
# These are the simplest forms.  More precise forms are derived further down.  #
################################################################################

# Sequential consistency
sc2ppo = MOST({
    "ld": {"ld": 'S', "st": 'S'},
    "st": {"ld": 'S', "st": 'S'}})

# Total Store Ordering (TSO)
tso2ppo = MOST({
    "ld": {"ld": 'S', "st": 'S'},
    "st": {"ld": '-', "st": 'M'}})

# SPARC Partial Store Ordering
pso2ppo = MOST({
    "ld": {"ld": 'S', "st": 'S'},
    "st": {"ld": '-', "st": '-'}})

# "Partial Load Ordering": by analogy to PSO
plo2ppo = MOST({
    "ld": {"ld": '-', "st": 'S'},
    "st": {"ld": '-', "st": 'M'}})

# "Load-Store Ordering only"
lso2ppo = MOST({
    "ld": {"ld": '-', "st": 'S'},
    "st": {"ld": '-', "st": '-'}})

# SPARC Relaxed Memory Ordering
rmo2ppo = MOST({
    "ld": {"ld": '-', "st": '-'},
    "st": {"ld": '-', "st": '-'}})

# The baseline versions of these architectures only have "mfence"
sc2 = Architecture("sc#2", sc2ppo, {}, ["ld", "st"], set())
tso2 = Architecture("tso#2", tso2ppo, {"mfence": sc2ppo}, ["ld", "st", "mfence"], set())
plo2 = Architecture("plo#2", plo2ppo, {"mfence": sc2ppo}, ["ld", "st", "mfence"], set())
pso2 = Architecture("pso#2", pso2ppo, {"mfence": sc2ppo}, ["ld", "st", "mfence"], set())
lso2 = Architecture("lso#2", lso2ppo, {"mfence": sc2ppo}, ["ld", "st", "mfence"], set())
rmo2 = Architecture("rmo#2", rmo2ppo, {"mfence": sc2ppo}, ["ld", "st", "mfence"], set())

################################################################################

powerppo = MOST({
        "ld" : {"ld": '-', "cld": '-', "st": '-', "cst": '-'},
        "cld": {"ld": '-', "cld": '-', "st": '-', "cst": '-'},
        "st" : {"ld": '-', "cld": '-', "st": '-', "cst": '-'},
        "cst": {"ld": '-', "cld": '-', "st": '-', "cst": '-'}})

powerlwsync = MOST({
        "ld" : {"ld": 'S', "cld": 'S', "st": 'S', "cst": 'S'},
        "cld": {"ld": 'S', "cld": 'S', "st": 'S', "cst": 'S'},
        "st" : {"ld": '-', "cld": '-', "st": 'N', "cst": 'N'},
        "cst": {"ld": '-', "cld": '-', "st": 'N', "cst": 'N'}})

powersync = MOST({
        "ld" : {"ld": 'S', "cld": 'S', "st": 'S', "cst": 'S'},
        "cld": {"ld": 'S', "cld": 'S', "st": 'S', "cst": 'S'},
        "st" : {"ld": 'S', "cld": 'S', "st": 'S', "cst": 'S'},
        "cst": {"ld": 'S', "cld": 'S', "st": 'S', "cst": 'S'}})

power = Architecture("power#", powerppo,
        {"lwsync": powerlwsync, "sync": powersync},
        ["ld", "st", "lwsync", "sync"], {"cld", "cst"})

armppo = MOST({
        "ld" : {"ld": '-', "cld": '-', "st": '-', "cst": '-'},
        "cld": {"ld": '-', "cld": '-', "st": '-', "cst": '-'},
        "st" : {"ld": '-', "cld": '-', "st": '-', "cst": '-'},
        "cst": {"ld": '-', "cld": '-', "st": '-', "cst": '-'}})

armdmb = MOST({
        "ld" : {"ld": 'S', "cld": 'S', "st": 'S', "cst": 'S'},
        "cld": {"ld": 'S', "cld": 'S', "st": 'S', "cst": 'S'},
        "st" : {"ld": 'S', "cld": 'S', "st": 'S', "cst": 'S'},
        "cst": {"ld": 'S', "cld": 'S', "st": 'S', "cst": 'S'}})

arm = Architecture("arm#", armppo,
        {"dmb": armdmb},
        ["ld", "st", "dmb"], {"cld", "cst"})

################################################################################
# ISCA paper-specific Model Definitions                                        #
#                                                                              #
# These define variants of the above models which use the architectures        #
# described in the ISCA paper (i.e., for the shims implemented within the gem5 #
# OoO pipeline                                                                 #
################################################################################

# The following two fences are defined based on the behavior of the gem5 OoO
# pipeline's microarchitectural implementation of fences.
mlfence = MOST({
    "ld": {"ld": 'S', "st": 'S'},
    "st": {"ld": 'S', "st": '-'}})

msfence = MOST({
    "ld": {"ld": '-', "st": 'S'},
    "st": {"ld": '-', "st": 'S'}})

# "allfences" defines a fine-grained set of fences that aim to give the shim
# more fine-grained control over which fences to insert
allfences = {}
for sl in ['-', 'S']:
    for ll in ['-', 'S']:
        for ss in ['-', 'S']:
            for ls in ['-', 'S']:
                label = "fence"
                if ll == 'S':
                    label += '_LL'
                if ls == 'S':
                    label += '_LS'
                if sl == 'S':
                    label += '_SL'
                if ss == 'S':
                    label += '_SS'
                if label == "fence": # nothing enforced
                    continue
                allfences[label] = MOST({
                    "ld": {"ld": ll, "st": ls},
                    "st": {"ld": sl, "st": ss}})

################################################################################

# Architecture definitions
tso2d = Architecture("tso#2d", tso2ppo, {"mfence": sc2ppo, "mlfence": mlfence, "msfence": msfence},
        ["ld", "st", "mfence", "msfence", "mlfence"], set())
plo2d = Architecture("plo#2d", plo2ppo, {"mfence": sc2ppo, "mlfence": mlfence, "msfence": msfence},
        ["ld", "st", "mfence", "msfence", "mlfence"], set())
pso2d = Architecture("pso#2d", pso2ppo, {"mfence": sc2ppo, "mlfence": mlfence, "msfence": msfence},
        ["ld", "st", "mfence", "msfence", "mlfence"], set())
lso2d = Architecture("lso#2d", lso2ppo, {"mfence": sc2ppo, "mlfence": mlfence, "msfence": msfence},
        ["ld", "st", "mfence", "msfence", "mlfence"], set())
rmo2d = Architecture("rmo#2d", rmo2ppo, {"mfence": sc2ppo, "mlfence": mlfence, "msfence": msfence},
        ["ld", "st", "mfence", "msfence", "mlfence"], set())

# A variant of RMO with lots of different fine-grained fences, for comparison
rmo16_2d = Architecture("rmo16#2d", rmo2ppo, allfences,
        ["ld", "st"] + list(allfences.keys()), set())

################################################################################

# Define a simplified (non-cumulative) form of "lwsync" to go along with the
# simplified form of Power
lwsync2 = MOST({
    "ld": {"ld": 'S', "st": 'S'},
    "st": {"ld": '-', "st": 'N'}})

# A simplified, multiple copy atomic variant of Power
powera = Architecture("powera#2", rmo2ppo, {"sync": sc2ppo, "lwsync": lwsync2},
        ["ld", "st", "lwsync", "sync"], set())

################################################################################
# Functions for automatically extending the above simpler definitions into     #
# more precise forms                                                           #
################################################################################

def SplitAddrVersion(arch, st_stsa, ld_ldsa):
    """Given an architecture with "B" instructions "ld" and "st", split each
    into "same address" vs. "different address" categories.  Also, add the
    baseline "coherence" ordering to the PPO.  For architectures such as RMO
    which do not enforce all of these orderings, the PPO will have to be
    manually weakened.
    
    Take the ordering for "st -> st&sa" as an argument.  It may be 'M' on
    multi-copy atomic architectures, and it may be 'N' on non-MCA architectures,
    so don't assume either one.
    
    Take the ordering for "ld -> ld&sa" as an argument.  It is usually 'S' on
    most architectures, but some (e.g., RMO) may occasionally weaken it, so
    allow that to be passed in.
    """
    def SplitMOST(most):
        baseline = {
                "ld" : {"ld&sa": ld_ldsa, "st&sa": 'S'    },
                "st" : {"ld&sa": 'L'    , "st&sa": st_stsa}}
        result = {}
        for a, bs in most.items():
            result.setdefault(a, {})
            for b, s in bs.items():
                if b in ["ld", "st"]:
                    result[a][b + "&da"] = s
                    if s == '-':
                        result[a][b + "&sa"] = baseline.get(a, {}).get(b + "&sa", '-')
                    else:
                        result[a][b + "&sa"] = s
                else:
                    result[a][b] = s
        return MOST(result)

    newppo = SplitMOST(arch.ppo.most)
    visible_ops = arch.visible_ops[:] + ["ld&sa", "ld&da", "st&sa", "st&da"]
    visible_ops.remove("ld")
    visible_ops.remove("st")
    result = Architecture(arch.code_name + "s", newppo, {}, visible_ops,
            arch.invisible_ops, verify=False)
    result.print_name = arch.print_name

    # Copy the values from each old MOST into the equivalent new MOST.
    for k, v in arch.mosts.items():
        result.mosts[k] = SplitMOST(v.most)

    logging.debug("Split address version of...")
    for ln in arch.Dump().split("\n"):
        logging.debug(ln)
    logging.debug("...is:")
    for ln in result.Dump().split("\n"):
        logging.debug(ln)

    result.Verify()
    return result

def CumulativeVersion(arch):
    """Given a split address version of an architecture, generate a version in
    which cumulativity is both explicitly included and enforced by default."""
    newppo = {}
    for a, bs in arch.ppo.most.items():
        newppo[      a] = {"cld": 'S', "cst": 'S'}
        newppo["c" + a] = {"cld": 'S', "cst": 'S'}
        for b, s in bs.items():
            newppo[      a][b] = s
            newppo["c" + a][b] = 'S'
    result = Architecture(arch.code_name + "C", MOST(newppo), {},
            arch.visible_ops, arch.invisible_ops | {"cld", "cst"}, verify=False)
    result.print_name = arch.print_name

    # Copy the non-cumulativity values from the old PPO into the new PPO.
    for a, bs in arch.ppo.most.items():
        for b, s in bs.items():
            result.ppo.most[a][b] = s

    # Copy the non-cumulativity values from each old MOST into the
    # corresponding new MOST.
    for k, v in arch.mosts.items():
        new_most = copy.deepcopy(newppo)
        for a, bs in v.most.items():
            new_most.setdefault(a, {})
            for b, s in bs.items():
                new_most[a][b] = s
        result.mosts[k] = MOST(new_most)

    logging.debug("Cumulative version of...")
    for ln in arch.Dump().split("\n"):
        logging.debug(ln)
    logging.debug("...is:")
    for ln in result.Dump().split("\n"):
        logging.debug(ln)

    result.Verify()
    return result

def CumulativeSplitAddrVersion(arch, st_stsa, ld_ldsa):
    """Return a version of arch which addresses cumulativity and same vs.
    different addresses."""
    return CumulativeVersion(SplitAddrVersion(arch, st_stsa, ld_ldsa))
