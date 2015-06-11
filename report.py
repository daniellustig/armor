#!/usr/bin/env python

################################################################################
# report.py                                                                    #
#                                                                              #
# Script for generating a report summarizing the ArMOR DBT Shim FSMs generated #
# for a large number of upstream and downstream MCMs                           #
################################################################################
# The ArMOR project                                                            #
#                                                                              #
# Copyright (C) 2015 Daniel Lustig, Princeton University                       #
#                                                                              #
# This software may be modified and distributed under the terms                #
# of the MIT license.  See the LICENSE file for details.                       #
################################################################################

import sys
import os
import argparse
import logging

################################################################################

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

################################################################################

# These need to be imported *after* the loggers are set up
import armor
from architectures import *

################################################################################

def GenerateLaTeXSummary(upstreams, downstreams, path, header=None):
    """Generate a summary of the ArMOR DBT shim FSM generated for each
    upstream->downstream combination"""
    try:
        os.makedirs(path)
    except os.error:
        pass

    # The script used to compile each of the graphviz files
    f = open("%s/build_graphs.sh" % path, 'w')
    f.write("#!/bin/bash\n\n")

    num_fsm_nodes = {}
    unmin_fsms = {}
    min_fsms = {}
    for up in upstreams:
        for down in downstreams:
            logging.debug("Considering table entry:")
            logging.debug("Upstream:")
            for ln in up.Dump().split('\n'):
                logging.debug(ln)
            logging.debug("Downstream:")
            for ln in down.Dump().split('\n'):
                logging.debug(ln)

            if (up.print_name == down.print_name):
                # No need to translate an MCM into itself
                continue

            if (up.ppo < down.ppo):
                # If upstream is weaker than downstream, this would correspond
                # to a fence removal scenario.  This works, and some fences are
                # removed but it isn't yet efficient...
                continue

            sys.stderr.write("(%d x %d) %s --> %s\n" %
                    (len(up.ppo.most[list(up.ppo.most.keys())[0]].keys()),
                        len(up.ppo.most.keys()),
                        up.print_name, down.print_name))
            logging.debug("Upstream %s - Downstream %s =" %
                    (up.print_name, down.print_name))
            for ln in armor.TableString((up.ppo - down.ppo).most).split('\n'):
                logging.debug(ln)

            # Generate the FSM for this (up, down) pair
            min_fsm, unmin_fsm = armor.GenerateFSMAndPrintDOTGraph(path, up, down)

            # Save the results
            min_fsms.setdefault(up.print_name, {})
            min_fsms[up.print_name][down.print_name] = min_fsm
            unmin_fsms.setdefault(up.print_name, {})
            unmin_fsms[up.print_name][down.print_name] = unmin_fsm
            num_fsm_nodes.setdefault(up.print_name, {})
            num_fsm_nodes[up.print_name][down.print_name] = len(min_fsm.nodes.keys())

            # Tell the script to build these 
            f.write("echo Drawing %s/%s_%s_minimized.pdf\n" %
                    (path, up.print_name, down.print_name))
            f.write("dot -Tpdf %s/%s_%s_minimized.gv -o %s/%s_%s_minimized.pdf\n" %
                    (path, up.print_name, down.print_name,
                        path, up.print_name, down.print_name))
            if unmin_fsm:
                f.write("echo Drawing %s/%s_%s_nonminimized.pdf\n" %
                        (path, up.print_name, down.print_name))
                f.write("dot -Tpdf %s/%s_%s_nonminimized.gv -o %s/%s_%s_nonminimized.pdf\n" %
                        (path, up.print_name, down.print_name,
                            path, up.print_name, down.print_name))

    f.write("\n\n")
    f.close()

    # Log the table with the number of nodes per FSM
    for ln in (armor.TableString(num_fsm_nodes,
        map(lambda x: x.print_name, upstreams),
        map(lambda x: x.print_name, downstreams),
        newline='\n')).split('\n'):
        logging.info(ln)

    # Emit the LaTeX summary of the FSMs

    f = open("%s/summary.tex" % path, 'w')

    f.write("""
    \\documentclass{article}
    \\usepackage[T1]{fontenc}
    \\usepackage{array}
    \\usepackage{graphicx}
    \\usepackage{fullpage}
    \\usepackage{rotating}
    \\usepackage{subcaption}
    \\usepackage{chngcntr}
    \\usepackage[table]{xcolor}
    \\usepackage[export]{adjustbox}
    \\usepackage{hyperref}
    \\renewcommand{\\thesection}{A.\\arabic{section}}
    \\counterwithin{figure}{section}
    \\renewcommand{\\thefigure}{A.\\arabic{section}.\\alph{figure}}
    \\begin{document}
    """)

    if header:
        f.write(header)

    f.write("""
    \\section{Summary}

    """)

    f.write("\\begin{figure*}[h!]\n")
    f.write("\\centering\n")
    f.write("\\begin{tabular}{|c|")
    f.write(''.join(["c|" for _ in downstreams]))
    f.write("}\n  ")

    # top left cell
    f.write("\\hline\n  ")

    # column headings
    for down in downstreams:
        f.write("& %s" % down.print_name)
    f.write(" \\\\\\hline\n")

    # rows
    for up in upstreams:
        f.write("  %s \n" % up.print_name)
        for down in downstreams:
            f.write(" & ")
            if down.print_name in num_fsm_nodes.get(up.print_name, {}):
                f.write("%d" % num_fsm_nodes[up.print_name][down.print_name])
            else:
                f.write("- \n")
        f.write("  \\\\\\hline\n")

    # end
    f.write("\\end{tabular}")
    f.write("""
        \caption{FSM node count summary.}
    \end{figure*}

    \clearpage

    """)

    f.write("\\begin{sidewaysfigure}\n")
    f.write("\\centering\n")
    f.write("\\begin{tabular}{|>{\\centering\\arraybackslash} m{1.8cm}|")
    f.write(''.join([">{\\centering\\arraybackslash} m{1.6cm}|" for _ in downstreams]))
    f.write("}\n  ")

    # top left cell
    f.write("\\hline\n  ")

    # column headings
    for down in downstreams:
        f.write("& %s" % down.print_name)
    f.write(" \\\\\\hline\n")

    # rows
    for up in upstreams:
        f.write("  %s \n" % up.print_name)
        for down in downstreams:
            if down.print_name in num_fsm_nodes.get(up.print_name, {}):
                f.write("  & ")
                f.write("\\includegraphics[width=1.6cm]")
                f.write("{%s_%s_minimized}\n" %
                        (up.print_name, down.print_name))
            else:
                f.write("  & - \n")
        f.write("  \\\\\\hline\n")

    # end
    f.write("\\end{tabular}")
    f.write("""
        \caption{FSM summary.}
    \end{sidewaysfigure}

    \clearpage

    """)

    # Upstream MOSTs
    f.write("\\section{Upstream MCMs}\n")
    for up in upstreams:
        f.write("  \\begin{figure}[h!]\n")
        f.write("    \\centering\n")
        f.write("    %s\n" % up.ppo.LaTeXRepr("PPO").replace("\n", "\n    "))
        for name, m in up.mosts.items():
            f.write("    ~\n")
            f.write("    %s\n" % m.LaTeXRepr(name).replace("\n", "\n    "))
        f.write("    \\caption{MOSTs for Upstream %s}\n" % up.print_name)
        f.write("  \\end{figure}\n")
    f.write("\\clearpage\n")

    # Downstream MOSTs
    f.write("\\section{Downstream MCMs}\n")
    for down in downstreams:
        f.write("  \\begin{figure}[h!]\n")
        f.write("    \\centering\n")
        f.write("    %s\n" % down.ppo.LaTeXRepr("PPO").replace("\n", "\n    "))
        for name, m in down.mosts.items():
            f.write("    ~\n")
            f.write("    %s\n" % m.LaTeXRepr(name).replace("\n", "\n    "))
        f.write("    \\caption{MOSTs for Downstream %s}\n" % down.print_name)
        f.write("  \\end{figure}\n")
    f.write("\\clearpage\n")

    for up in upstreams:
        for down in downstreams:
            if down.print_name in num_fsm_nodes.get(up.print_name, {}):
                f.write("\\section{%s Upstream, %s Downstream}\n" %
                        (up.print_name, down.print_name))

                # PPO Difference
                f.write("  \\begin{figure}[h!]\n")
                f.write("  \\centering \\footnotesize\n")
                f.write("  \\begin{minipage}{0.3\\textwidth}\n")
                f.write("    \\centering\n")
                diff = up.ppo - down.ppo
                f.write("    %s\n" % diff.LaTeXRepr("PPO Diff.").replace("\n", "\n    "))
                f.write("    \\caption{PPO of (Upstream %s - Downstream %s)}\n" %
                        (up.print_name, down.print_name))
                f.write("  \\end{minipage}\n")

                f.write("  \\qquad\n")

                # FSM Transition Table
                fsm = min_fsms[up.print_name][down.print_name]
                f.write("  \\begin{minipage}{0.6\\textwidth}\n")
                f.write("    \\centering")
                f.write("    %s\n" % fsm.LaTeXRepr().replace("\n", "\n    "))
                f.write("    \\caption{FSM Transition Table}\n")
                f.write("  \\end{minipage}\n")
                f.write("  \\end{figure}\n")

                if unmin_fsms[up.print_name][down.print_name]:
                    # FSM Picture (pre-minimization)
                    f.write("  \\begin{figure}[h!]\n")
                    f.write("    \\centering")
                    f.write("    \\begin{minipage}{0.45\\textwidth}\n")
                    f.write("    \\includegraphics[max width=0.95\\textwidth,max height=5in]")
                    f.write("{%s_%s_nonminimized}\n" % (up.print_name, down.print_name))
                    f.write("    \\caption{FSM (Pre-minimization)}\n")
                    f.write("    \\end{minipage}\n")

                    f.write("    \\qquad\n")

                    # FSM Picture
                    f.write("    \\begin{minipage}{0.45\\textwidth}\n")
                    f.write("    \\includegraphics[max width=0.95\\textwidth,max height=5in]")
                    f.write("{%s_%s_minimized}\n" % (up.print_name, down.print_name))
                    f.write("    \\caption{FSM}\n")
                    f.write("    \\end{minipage}\n")
                    f.write("  \\end{figure}\n")
                else:
                    # FSM Picture
                    f.write("  \\begin{figure}[h!]\n")
                    f.write("    \\centering")
                    f.write("    \\includegraphics[max width=5in,max height=5in]")
                    f.write("{%s_%s_minimized}\n" % (up.print_name, down.print_name))
                    f.write("    \\caption{FSM}\n")
                    f.write("  \\end{figure}\n")

                f.write("\\clearpage\n")

    f.write("""

    \end{document}
    """)

################################################################################

upstreams   = [sc2, tso2 , plo2 , pso2 , lso2 , rmo2 ,           powera]
downstreams = [     tso2d, plo2d, pso2d, lso2d, rmo2d, rmo16_2d, powera]

GenerateLaTeXSummary(upstreams, downstreams, args.graphfilepath + "/isca-2x2")

################################################################################

def f(ld_ldsa):
    return lambda arch: SplitAddrVersion(arch, ld_ldsa=ld_ldsa, st_stsa="M")
upstreams   = map(f('S'), [sc2, tso2 , plo2 , pso2 , lso2 ]) + [f('-')(rmo2 ),                   f('S')(powera)]
downstreams = map(f('S'), [     tso2d, plo2d, pso2d, lso2d]) + [f('-')(rmo2d), f('-')(rmo16_2d), f('S')(powera)]

GenerateLaTeXSummary(upstreams, downstreams, args.graphfilepath + "/isca-2x4")

################################################################################

# The version used for the ISCA paper!
header = """
    \\section*{Appendix}

    The following appendix presents a gallery of MOSTs and shim FSMs used in the evaluation of Sections~6 and~7.  Because addresses may not all be resolved by the time instructions reach the issue queue, we ignore same-address orderings for this particular implementation.

    \\noindent The code which automatically generated these shims (and this appendix) is available online:

    \url{https://github.com/daniellustig/armor}.
    """
upstreams   = map(CumulativeVersion, [sc2, tso2 , plo2 , pso2 , lso2 , rmo2 ,           powera]) + [power, arm]
downstreams = map(CumulativeVersion, [     tso2d, plo2d, pso2d, lso2d, rmo2d, rmo16_2d, powera]) + [power, arm]

GenerateLaTeXSummary(upstreams, downstreams, args.graphfilepath + "/isca-4x4", header)

################################################################################

def f(ld_ldsa):
    return lambda arch: CumulativeSplitAddrVersion(arch, ld_ldsa=ld_ldsa, st_stsa="M")
def g(arch):
    return SplitAddrVersion(arch, ld_ldsa='S', st_stsa='N')

upstreams   = map(f('S'), [sc2, tso2 , plo2 , pso2 , lso2 ]) + [f('-')(rmo2 ),                   f('S')(powera)] + map(g, [power, arm])
downstreams = map(f('S'), [     tso2d, plo2d, pso2d, lso2d]) + [f('-')(rmo2d), f('-')(rmo16_2d), f('S')(powera)] + map(g, [power, arm])

GenerateLaTeXSummary(upstreams, downstreams, args.graphfilepath + "/isca-4x6")

################################################################################

def f(ld_ldsa):
    return lambda arch: CumulativeSplitAddrVersion(arch, ld_ldsa=ld_ldsa, st_stsa="M")
def g(arch):
    return SplitAddrVersion(arch, ld_ldsa='S', st_stsa='N')

arches = map(f('S'), [sc2, tso2, plo2, pso2, lso2]) + [f('-')(rmo2), f('S')(powera)] + map(g, [power, arm])

GenerateLaTeXSummary(arches, arches, args.graphfilepath + "/normal")
