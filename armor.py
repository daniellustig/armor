#!/usr/bin/env python

################################################################################
#                                                                              #
# armor.py                                                                     #
#                                                                              #
# Core ArMOR analysis and maniupulation routines.                              #
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
import re
import unittest
import logging

class ArMORError(Exception):
    def __init__(self, msg):
        self.msg = msg
    def __str__(self):
        return self.msg

class InsufficientMOR(ArMORError):
    pass

class StrengthLevel():
    """For now, assume a total order on strength levels.  This could be changed
    if necessary, but for now it makes the coding easier"""
    strengths = '-LNMS'

    def __init__(self, level):
        self.level = level

    # Primitive comparison operators:
    def __eq__(self, other):
        return self.level == other.level
    def __le__(self, other):
        return (StrengthLevel.strengths.index(self.level) <=
                StrengthLevel.strengths.index(other.level))
    # Non-primitive comparison operators (defined in terms of the primitive
    # operators)
    def __ne__(self, other):
        return not self == other
    def __ge__(self, other):
        return other <= self
    def __lt__(self, other):
        return self <= other and self != other
    def __gt__(self, other):
        return self >= other and self != other

    def __add__(self, other):
        """The join (/\) operator, represented using the '+' in Python.  Since
        we are assuming a total order on strength levels for now, the join is
        implemented as the max of the input strength levels."""
        return StrengthLevel.strengths[
                max(StrengthLevel.strengths.index(self.level),
                    StrengthLevel.strengths.index(other.level))]

    def __sub__(self, other):
        """The subtraction operator, represented using the '-' in Python.
        Subtraction here is a conservative approximation: for 'a - b', if b
        is weaker than a, return all of a (i.e., ignore the fact that some
        parts of a have been enforced by b, and consider them still pending
        anyway)."""
        if self <= other:
            return StrengthLevel.strengths[0]
        else:
            return self.level

def APart(op):
    """A and B parts of operations:

    We provide the ability for multiple B operation types to map onto the same
    A operation type.  This is most useful (so far) when describing "same
    address" vs. "different address" mappings.  For example, the MOST for TSO
    looks like this:

            +-----------------------------------------|
      TSO:  | Load/Same Addr | Load/Diff Addr | Store |
    +-------+----------------+----------------+-------|
    | Load  |       S        |        S       |   S   |
    |-------+----------------+----------------+-------|
    | Store |       L        |        -       |   M   |
    +-------+-----------------------------------------|

    ...and we want to map both types of B load onto the same type of A load.
    The syntax for this is as follows:

    If a B operation has an ampersand (&) in its name, then it is mapped onto
    the A operation which has the name matching the portion of the name coming
    before the ampersand.  For example, the above table would look like this:

          +--------------------+
     TSO: | ld&sa | ld&da | st |
    +-----+-------+-------+----|
    | ld  |   S   |   S   | S  |
    |-----+-------+-------+----|
    | st  |   L   |   -   | M  |
    +-----+--------------------+
    """
    return re.sub("&.*", "", op)

class MOST():
    """A memory ordering specification table.

    A MOST is built using nested dictionaries.  Each entry in the outer
    dictionary represents a particular row within the MOST (the key is the row
    heading, and the value is the set of inner dictionaries).  Each inner
    dictionary represents cell in the MOST (the key is the column heading, and
    the value is the MOST cell strength value).
    """
    def __init__(self, most):
        self.most = most

    def __hash__(self):
        """Allow MOSTs to serve as dictionary keys"""
        try:
            return self.memoized_hash
        except AttributeError:
            pass
        # In Python, can only take the hash of an immutable object.  Therefore,
        # create an immutable equivalent to the MOST before hashing.  This
        # requires recursion into the second level of dicts.
        subhashes = {}
        for k, v in self.most.items():
            subhashes[k] = hash(frozenset(v.items()))
        self.memoized_hash = hash(frozenset(subhashes.items()))
        return self.memoized_hash

    def __repr__(self):
        """Return a string concatenating the orderings in the the table, row
        by row."""
        try:
            return self.memoized_str
        except:
            pass
        self.memoized_str = ''
        for _, bs in self.most.items():
            for _, v in bs.items():
                self.memoized_str += v
            self.memoized_str += ' '

        # drop the trailing space
        self.memoized_str = self.memoized_str[:-1]
        return self.memoized_str

    def LaTeXRepr(self, name=None):
        """Return a LaTeX tabular environment depicting the MOST."""
        s = "\\begin{tabular}{|c|"
        columns = set()
        for row in self.most.values():
            columns |= set(row.keys())
        for _ in columns:
            s += "c|"
        s += "}\n"

        if name is not None:
            s += ("\\multicolumn{%d}{c}{%s} \\\\\n" %
                    (len(columns) + 1,
                        name.replace("_", "\\_")))

        s += "\cline{2-%d}\n" % (len(columns) + 1)

        s += "\\multicolumn{1}{c|}{} "
        for b in columns:
            s += "& %s " % b.replace("&", "\\&")
        s += "\\\\\\hline\n"

        for a, bs in self.most.items():
            s += "%s " % a
            for b in columns:
                s += "& %s " % bs.get(b, "?")
            s += "\\\\\\hline\n"

        s += "\\end{tabular}"

        return s

    # Primitive comparison operators:
    def __eq__(self, other):
        for (a, bs) in self.most.items():
            for (b, s) in bs.items():
                if s != other.most[a][b]:
                    return False
        return True
    def __le__(self, other):
        for (a, bs) in self.most.items():
            for (b, s) in bs.items():
                if not StrengthLevel(s) <= StrengthLevel(other.most[a][b]):
                    return False
        return True
    # Non-primitive comparison operators (defined in terms of the primitive
    # operators)
    def __ne__(self, other):
        return not self == other
    def __ge__(self, other):
        return other <= self
    def __lt__(self, other):
        return self <= other and self != other
    def __gt__(self, other):
        return self >= other and self != other

    def KeepColumn(self, op):
        """Return a MOST of the same shape (i.e., with the same rows and
        columns) but with no ordering requirements in columns other than the
        one specified."""
        result = {}
        for a, bs in self.most.items():
            result[a] = {}
            for b, s in bs.items():
                if b == op:
                    result[a][b] = s
                else:
                    result[a][b] = '-'
        return MOST(result)

    def KeepRow(self, op):
        """Return a MOST of the same shape (i.e., with the same rows and
        columns) but with no ordering requirements in rows other than the
        one specified."""
        result = {}
        for a, bs in self.most.items():
            result[a] = {}
            for b, s in bs.items():
                if a == APart(op):
                    result[a][b] = s
                else:
                    result[a][b] = '-'
        return MOST(result)

    def __add__(self, other):
        """The join (/\) operator, represented using the '+' in Python."""
        result = {}
        for a, bs in self.most.items():
            result[a] = {}
            for b, s in bs.items():
                result[a][b] = StrengthLevel(s) + StrengthLevel(other.most[a][b])
        return MOST(result)

    def __sub__(self, other):
        """The subtraction operator, represented using the '-' in Python."""
        result = {}
        for a, bs in self.most.items():
            result[a] = {}
            for b, s in bs.items():
                result[a][b] = StrengthLevel(s) - StrengthLevel(other.most[a][b])
        return MOST(result)

    def EmptyMOSTOfSameShape(self):
        """Return a MOST of the same shape (i.e., with the same rows and
        columns) but with no ordering requirements."""
        result = {}
        for a, bs in self.most.items():
            result[a] = {}
            for b, _ in bs.items():
                result[a][b] = '-'
        return MOST(result)

    def Superset(self, other):
        """Return a MOST that describes the merging of two states.  If a pair
        of corresponding cells is the same, keep the value.  If a pair of
        corresponding cells is different, make the resulting cell '*' to make
        it clear that the original input differ in that cell."""
        result = {}
        for a, bs in self.most.items():
            result.setdefault(a, {})
            for b, s in bs.items():
                if s != other.most.get(a, {}).get(b, '*'):
                    result[a][b] = '*'
                else:
                    result[a][b] = s
        for a, bs in other.most.items():
            result.setdefault(a, {})
            for b, s in bs.items():
                if s != self.most.get(a, {}).get(b, '*'):
                    result[a][b] = '*'
                else:
                    result[a][b] = s
        return MOST(result)

def TableString(t, a_keys=None, b_keys=None, newline='\n', sep=' ', start='', end='', print_keys=True):
    """Return the contents of "t" nicely formatted as a text table"""
    def spaces(n):
        return ''.join([' ' for _ in range(n)])

    result = start

    if not a_keys:
        a_keys = t.keys()
    max_a_length = max(map(len, a_keys))
    a_keys = list(a_keys)
    a_keys.sort()

    if not b_keys:
        b_keys = set()
        for a_values in t.values():
            b_keys |= set(a_values.keys())
    b_keys = list(b_keys)
    b_keys.sort()

    # Calculate width of each column
    max_b_length = {}
    for a in a_keys:
        for b in b_keys:
            max_b_length.setdefault(b, len(b))
            max_b_length[b] = max(
                    max_b_length[b], len(str(t.get(a, {}).get(b, ' '))))

    # Top left corner
    # to make the hash (which may be negative) printable
    if print_keys:
        result += spaces(max_a_length)

    # Column headings
    if print_keys:
        for b in b_keys:
            result += sep + b

    # Rows
    for a in a_keys:
        # Put the newline here so that there's no newline at the very end
        result += newline;

        # Row heading
        if print_keys:
            result += a + spaces(max_a_length - len(a))

        # Cells
        for b in b_keys:
            cell = str(t.get(a, {}).get(b, ' '))
            result += sep + cell + spaces(max_b_length[b] - len(cell))

    result += end

    return result

class Architecture():
    def __init__(self, name, ppo, mosts, visible_ops, invisible_ops=set(), verify=True):
        """The name of an architecture has two parts: a printable part and a
        non-printable part.  There does not need to be a non-printable part,
        but if there is, it is separated from the printable part by a '#".
        The non-printable part is useful for distinguishing two different
        representations for the same architecture: e.g., "tso#2x2" vs. "tso#4x4"
        """
        self.code_name = name
        self.print_name = re.sub("#.*", "", name).upper()
        self.ppo = ppo
        self.mosts = mosts
        self.visible_ops = visible_ops
        self.invisible_ops = invisible_ops
        if verify:
            self.Verify()
    def __repr__(self):
        result  = "Architecture %s (%s)\n" % (self.print_name, self.code_name)
        result += "\t%s\tPPO" % self.ppo
        for k, v in self.mosts.items():
            result += "\n\t%s\t%s" % (v, k)
        return result
    def Dump(self):
        result  = "Architecture %s (%s)\n" % (self.print_name, self.code_name)
        result += "PPO:\n%s\n" % TableString(self.ppo.most)
        for k, v in self.mosts.items():
            result += "%s:\n%s\n" % (k, TableString(v.most))
        return result
    def Verify(self):
        for n, m in self.mosts.items():
            if not m >= self.ppo:
                logging.info("Architecture %s (%s) MOST %s is not as strong as PPO" %
                        (self.print_name, self.code_name, n))
                logging.info("PPO:")
                for ln in TableString(m.most).split('\n'):
                    logging.info(ln)
                logging.info("MOST %s:" % n)
                for ln in TableString(self.ppo.most).split('\n'):
                    logging.info(n)
                diff = self.ppo - m
                logging.info("PPO - MOST %s" % n)
                for ln in TableString(diff.most).split('\n'):
                    logging.info(ln)
                newmost = self.ppo + m
                logging.info("Updating to:")
                for ln in TableString(newmost.most).split('\n'):
                    logging.info(ln)
                self.mosts[n] = newmost

class Transition():
    """A transition from one shim FSM state to another"""
    def __init__(self, src, dst, mor_in, mors_out):
        self.src = src
        self.dst = dst
        self.mor_in = mor_in
        self.mors_out = mors_out
    def __eq__(self, other):
        return self.src == other.src and self.dst == other.dst and \
                self.mor_in == other.mor_in and \
                self.mors_out == other.mors_out
    def __hash__(self):
        return hash((self.src, self.dst, self.mor_in, frozenset(self.mors_out)))
    def __repr__(self):
        return "%s --(%s/%s)--> %s" % (
                self.src, self.mor_in, str(self.mors_out), self.dst)

def MinimalMOST(l, orderingsToEnforce):
    """Find the minimal MOST among the MOSTs in l, where l is a dictionary of
    (name, MOST) pairs"""
    minimal = None
    # TODO: Very naive inefficient algorithm!
    for i_name, i_most in l.items():
        not_min = False
        for j_name, j_most in l.items():
            if i_most > j_most:
                logging.debug("%s not minimal: greater than %s" %
                        (i_most, j_most))
                not_min = True
                break
            else:
                logging.info("%s maybe minimal: not greater than %s" %
                        (i_most, j_most))
        if not not_min:
            if minimal:
                logging.info("More than one minimal MOST within set!")
                logging.info("\t%s (requirement)" % orderingsToEnforce)
                logging.info("\t%s %s" % (minimal[1], minimal[0]))
                logging.info("\t%s %s" % (i_most, i_name))
                logging.info("\t(and maybe more...)")
            else:
                minimal = i_name, i_most
    return minimal

def WeakestSufficientMOR(orderingsToEnforce, mosts):
    """Find the minimal MOST among the MOSTs in mosts, where mosts is a
    dictionary of (name, MOST) pairs"""
    valid_mosts = {}
    for name, most in mosts.items():
        if orderingsToEnforce <= most:
            logging.info("\t  sufficient: %s %s" % (most, name))
            valid_mosts[name] = most
        else:
            logging.info("\tinsufficient: %s %s" % (most, name))
            logging.debug("\t\t%s" % (orderingsToEnforce - most))
    if not valid_mosts:
        logging.warning("No sufficiently strong MOST!")
        raise InsufficientMOR("No sufficiently strong MOST!")
    return MinimalMOST(valid_mosts, orderingsToEnforce)

def NextState(upstream, downstream, assumedReqs, currentState, op):
    """Given the current state of a shim, the upstream and downstream PPO, the
    assumedReqs, and the input operation type, calculate the next state to go
    to and the set of downstream MORs (if any) to insert"""

    logging.info("")
    logging.info("Next State Iteration")
    logging.debug("Upstream:")
    for ln in str(upstream).split('\n'):
        logging.debug(ln)
    logging.debug("Downstream:")
    for ln in str(downstream).split('\n'):
        logging.debug(ln)
    logging.info("Current state: %s" % str(currentState))
    logging.info("MOR: %s" % str(op))

    # The orderings that need to be enforced are:
    # 1) those in the PPO MOST column corresponding to the input operation
    #    (for input operation types without an explicit MOST, e.g., loads
    #     and stores)
    # 2) those in the MOST associated with the input operation
    #    (for input operation types with an explicit MOST, e.g., fences)
    # 3) the assumedReqs (i.e., any accesses which can't be observed directly
    #    and which therefore must be conservatively assumed to be required
    orderingsToEnforce = currentState.KeepColumn(op)
    if op in upstream.mosts.keys():
        orderingsToEnforce += upstream.mosts[op]
    orderingsToEnforce += assumedReqs

    # The orderings that will be marked pending once this input operation has
    # been processed are those which:
    # 1) are in the row of the upstream MOST corresponding the the input
    #    operation
    # 2) are not enforced by the corresponding entries in the downstream MOST
    #    row corresponding to the input operation
    newOrderings = (upstream.ppo - downstream.ppo).KeepRow(op)

    # If PPO is sufficiently strong by itself to enforce all of the necessary
    # orderings, then don't insert anything.  Otherwise, insert the weakest
    # sufficiently strong MOST
    if orderingsToEnforce <= downstream.ppo:
        logging.info("\tPPO sufficient")
        insertedMORs = []
        insertedMORMOST = downstream.ppo
        nextState = currentState + newOrderings + assumedReqs
    else:
        logging.info("\tPPO insufficient")
        logging.info("\t\t  %s" % orderingsToEnforce)
        logging.info("\t\t- %s" % downstream.ppo)
        logging.info("\t\t= %s" % (orderingsToEnforce - downstream.ppo))
        insertedMORName, insertedMORMOST = WeakestSufficientMOR(
                orderingsToEnforce, downstream.mosts)
        insertedMORs = [insertedMORName]
        nextState = currentState - insertedMORMOST + newOrderings + assumedReqs

    # One hard-coded optimization:
    # (st->st, check_S) - (st->st, check_M) = (st->ld, check)
    # In other words, the difference between single- and multi-copy atomicity
    # is only observable by a store following a load to the same address.
    # Single-copy atomicity prevents store buffer forwarding, while multi-copy
    # atomicity allows it.  For cases in which store->load ordering already
    # needs to be enforced, then the (check_S - check_M) store->store ordering
    # becomes redundant.
    for a in nextState.most.keys():
        if a != "st":
            continue
        for b in nextState.most[a].keys():
            if APart(b) == "st" and \
                    nextState.most[a][b] == 'S' and \
                    insertedMORMOST.most[a][b] == 'M':
                found_cell = False
                for k in nextState.most[a].keys():
                    if APart(k) == 'ld':
                        if nextState.most[a][k] != 'S':
                            logging.info(("No optimization: %s->%s is %s, not S") %
                                (a, k, nextState.most[a][k]))
                            found_cell = False
                            break
                        nextState.most[a][k] = 'S'
                        found_cell = True
                        logging.info(("Optimization! (%s->%s, check_S) - " +
                            "(%s->%s, check_M) = (%s->%s, check)") %
                            (a, b, a, b, a, k))

                if found_cell:
                    nextState.most[a][b] = '-'
                    logging.info(("Optimization! (%s->%s, check_S) - " +
                        "(%s->%s, check_M) = (%s->%s, -)") %
                        (a, b, a, b, a, b))

    # If the input operation has a corresponding downstream operation (e.g.,
    # loads and stores), then insert it (after the other necessary MOST(s) (if
    # any) have been inserted).
    if APart(op) in downstream.ppo.most.keys():
        insertedMORs.append(APart(op))

    logging.info("Inserted MORs: %s" % str(insertedMORs))
    logging.info("Next state: %s" % str(nextState))

    return Transition(currentState, nextState, op, insertedMORs)

def AssumedReqs(upstream, downstream):
    """Return the assumedReqs: the MOST corresponding to operations which cannot
    be directly observed and must therefore conservatively be assumed to always
    be pending enforcement"""
    result = upstream.ppo.EmptyMOSTOfSameShape()
    for op in downstream.invisible_ops:
        result += upstream.ppo.KeepRow(op)
        result += upstream.ppo.KeepColumn(op)
    logging.info("AssumedReqs for %s -> %s: %s" % (
        upstream.code_name, downstream.code_name, result))
    return result

class FSM():
    def __init__(self, upstream, downstream, edges):
        self.upstream = upstream
        self.downstream = downstream
        self.nodes = {}
        for e in edges:
            self.nodes.setdefault(e.src, {})
            self.nodes.setdefault(e.dst, {})
            self.nodes[e.src][e.mor_in] = (e.dst, e.mors_out)
        self.edges = edges

    def StateID(self, node):
        return self.nodes.keys().index(node)

    def __eq__(self, other):
        # FIXME: check if the upstreams and downstreams are actually the same
        if self.upstream.code_name != other.upstream.code_name:
            return False
        if self.downstream.code_name != other.downstream.code_name:
            return False

        for e in self.edges:
            if e not in other.edges:
                return False

        for e in other.edges:
            if e not in self.edges:
                return False

        return True

    def DOTGraph(self, label=""):
        """Return a DOT graph representing the FSM"""
        s = ""

        s += "digraph armor {\n"
        s += "labelloc=t;\n"
        if label:
            label = " (" + label + ")"
        s += ('label="%s -> %s%s";\n' %
                (self.upstream.print_name, self.downstream.print_name, label))

        # Use courier so that the FSMs are properly aligned
        s += '\tnode [fontname="courier"];\n'
        s += '\tedge [fontname="courier"];\n'

        # Print the edges
        for e in self.edges:
            # Print a user-readable form of the edge as a comment
            s += '\t//%s\n' % str(e)

            # Edge color scheme:
            # Red = edge for which at least one extra MOR had to be inserted
            # Green = edge for which no extra MOR had to be inserted
            # Blue = an edge with no MORs at all (uncommon?)
            if len(e.mors_out) == 0:
                logging.info("No MORs inserted...redundant fence? %s" % e)
                edge_color = "blue"
                mors_out = '-'
            elif len(e.mors_out) == 1:
                edge_color = "black"
                mors_out = e.mors_out
            else:
                edge_color = "red"
                mors_out = e.mors_out

            # Print the edge in DOT notation
            s += ('\tn%s -> n%s [label="%s/%s";color=%s;penwidth=5];\n' %
                    (str(hash(e.src)).replace('-', '_'),
                     str(hash(e.dst)).replace('-', '_'),
                     e.mor_in,
                     ';'.join(mors_out),
                     edge_color))

        # Print the nodes
        for n in self.nodes.keys():
            # The node shows the pending orderings table
            label = TableString(n.most, sep=' ', newline='\\n')
            s += ('\tn%s [label="%s";penwidth=5]\n' %
                    (str(hash(n)).replace("-", "_"), label))

        # End the graph
        s += "}\n"

        return s

    def LaTeXRepr(self):
        """Return a LaTeX tablular environment showing the FSM transition
        table"""
        s = ""
        s += "\\begin{tabular}{|c|c|c|c|c|}\n"
        s += "\\hline\n"
        s += "\\multicolumn{3}{|c|}{Input} & \\multicolumn{2}{c|}{Output} \\\\\\hline\n"
        s += "State & MOST & Op. & Op(s). & Next State \\\\\\hline\n"
        for src, dsts in self.nodes.items():
            for op, edge in dsts.items():
                s += "%d & \\texttt{%s} & %s & " % (
                        self.StateID(src), str(src).replace("-", "{-}"),
                        op.replace("_", "\\_").replace("&", "\\&"))
                s += '; '.join(
                        map(lambda x: x.replace("_", "\\_").replace("&", "\\&"),
                            edge[1])) if edge[1] else '-'
                s += " & %s \\\\\\hline\n" % self.StateID(edge[0])
        s += "\\end{tabular}\n"
        return s

    def TextRepr(self):
        """Return a text form of the FSM transition table"""
        d = {"0": {"1": "State", "2": "MOST", "3": "Input", "4": "|",
            "5": "Output", "6": "Next State"}}
        for src, dsts in self.nodes.items():
            for op, edge in dsts.items():
                d[str(len(d.keys()))] = {
                        "1": str(self.StateID(src)),
                        "2": src,
                        "3": op.replace("_", "\\_").replace("&", "\\&"),
                        "4": "|",
                        "5": '; '.join(edge[1]) if edge[1] else '-',
                        "6": self.StateID(edge[0])}
        return TableString(d, sep=" | ", print_keys=False)

    def MinimizedFSM(self, merge):
        """Given an FSM, merge any states which are identical in behavior"""

        logging.debug("MIN\tStart")

        # Calculate a pairwise set of conditions under which each pair of two
        # states may be identical.  comparison[i][j] is either 1) an assertion
        # that for i to be equal to j, a must be equal to b, or 2) an assertion
        # that i cannot be equivalent to j
        comparison = {}
        for ik, iv in self.nodes.items():
            for jk, jv in self.nodes.items():
                if hash(ik) >= hash(jk):
                    continue
                comparison.setdefault(ik, {})
                comparison[ik][jk] = set()
                for op, dst in iv.items():
                    if dst[1] != jv[op][1]:
                        # The transitions for input "op" is not the same;
                        # hence, the two source states cannot be the same
                        logging.info("MIN\t%s != %s because on op %s:" %
                                (ik, jk, op))
                        logging.info("MIN\t\t%s: %s %s" %
                                (ik, dst[1], iv[op][1]))
                        logging.info("MIN\t\t%s: %s %s" %
                                (jk, dst[1], jv[op][1]))
                        comparison[ik][jk] = None
                        break
                    else:
                        # The transitions for input "op" are the same; hence,
                        # the two source states may be the same, but only if
                        # the two destination states (a, b) are also the same
                        a = dst[0]
                        b = jv[op][0]
                        if a != b:
                            assert not a == b
                            if hash(a) < hash(b):
                                comparison[ik][jk].add((a, b))
                            else:
                                comparison[ik][jk].add((b, a))

        # Sanity check: there should have been (n*(n-1))/2 comparisons
        # performed
        len_comparison = 0
        for k, v in comparison.items():
            len_comparison += len(v.keys())
        assert len_comparison == (len(self.nodes) * (len(self.nodes) - 1)) / 2

        # Loop through the state equivalence conditions calculated above and
        # check whether any of them fail.  If they do, then mark the potential
        # state equalities depending on those conditions as invalid.  Loop
        # until convergence.
        iterate = True
        while iterate:
            iterate = False
            for i, js in comparison.items():
                for j, conditions in js.items():
                    if conditions is None:
                        continue
                    for (a, b) in conditions:
                        if comparison[a][b] is None:
                            logging.info("MIN\t%s != %s because %s != %s" %
                                    (i, j, a, b))
                            comparison[i][j] = None
                            iterate = True

        # At this point, any remaining conditional equivalences hold true.
        # Group them.
        equal_states = set()
        for i, js in comparison.items():
            for j, conditions in js.items():
                if conditions is None:
                    continue
                # i and j are equal.  Find any other equivalences involving i
                # or j and merge them together with (i=j).
                new_eq = {i, j}
                # Track the old equivalences that are superseded by the new
                # merged equivalence
                old_eqs = set()
                for eq in equal_states:
                    if i in eq or j in eq:
                        # i and j are equivalent to all of the states in eq
                        new_eq |= eq
                        # eq is now superseded by new_eq; remove it
                        old_eqs |= {eq}
                equal_states -= old_eqs
                equal_states |= {frozenset(new_eq)}
        logging.info("MIN\tEqual states list: %s" % equal_states)

        # Map each node to a canonical representative
        mapping = {}
        replacement = {}
        for n in self.nodes.keys():
            mapping[n] = n
        for eq in equal_states:
            choice = min(eq) # it doesn't really matter which one is chosen
            replacement[choice] = choice
            for x in eq:
                replacement[choice] = merge(replacement[choice], x)
                mapping[x] = choice
        logging.debug("MIN\tReplacement mapping:")
        for k, v in replacement.items():
            logging.debug("MIN\t\t%s -> %s" % (k, v))
        logging.info("MIN\tEqual states mapping:")
        for k, v in mapping.items():
            logging.info("MIN\t\t%s -> %s -> %s" %
                    (k, v, replacement.get(v, v)))

        # Merge the nodes
        new_edges = set()
        for a, bs in self.nodes.items():
            if mapping[a] != a:
                # this edge is now redundant: skip it
                continue
            for b, v in bs.items():
                dst, mors = self.nodes[a][b]
                new_edges.add(Transition(replacement.get(a, a),
                    replacement.get(mapping[dst], mapping[dst]), b, mors))

        return FSM(self.upstream, self.downstream, new_edges)

    def RemoveTransientNodes(self):
        """Remove any transient nodes: nodes which are not reachable from all
        other nodes in the graph"""

        # Calculate the adjacency list of the FSM
        adj = {}
        for e in self.edges:
            adj.setdefault(e.src, {})
            adj[e.src][e.dst] = True

        # Take the transitive closure of the adjacency list
        for k in self.nodes:
            for i in self.nodes:
                for j in self.nodes:
                    if adj.get(i, {}).get(j, False) and \
                            adj.get(j, {}).get(k, False):
                        adj[i][k] = True

        # If j is not reachable from i, then mark j as transient
        transient_nodes = set()
        for i in self.nodes:
            for j in self.nodes:
                if not adj.get(i, {}).get(j, False):
                    logging.info("MIN\tTransient node %s" % j)
                    transient_nodes.add(j)

        # Calculate the FSM with transient nodes removed
        new_edges = set()
        for e in self.edges:
            if e.src not in transient_nodes and e.dst not in transient_nodes:
                new_edges.add(e)

        return FSM(self.upstream, self.downstream, new_edges)

def GenerateFSM(upstream, downstream):
    """Given upstream and downstream architectures, generate the (not yet
    minimal) ArMOR shim FSM for translating from upstream to downstream."""

    assumedReqs = AssumedReqs(upstream, downstream)

    edges = set()
    visitedStates = set()

    # The starting state would be the empty state, except that the assumedReqs
    # is always considered pending, so the starting state is actually equal to
    # assumedReqs
    pendingStates = [assumedReqs]

    logging.info("")
    logging.info("")
    logging.info("Generating FSM for %s --> %s" %
            (upstream.code_name, downstream.code_name))
    for ln in upstream.Dump().split('\n'):
        logging.debug(ln)
    for ln in downstream.Dump().split('\n'):
        logging.debug(ln)

    while pendingStates:
        logging.debug("%d states remaining" % len(pendingStates))

        # Pick a pending state; it doesn't matter which one
        currentState = pendingStates.pop()
        # Add it to the list of visited states (so we don't enumerate it twice)
        visitedStates.add(hash(currentState))

        # Calculate the outgoing edge corresponding to each input operation
        for op in upstream.visible_ops:
            # Calculate the next state
            transition = NextState(upstream, downstream, assumedReqs,
                    currentState, op)
            edges.add(transition)

            # If this is a not-yet-visited state, add it to the list of
            # states pending exploration
            if hash(transition.dst) not in visitedStates and \
                    hash(transition.dst) not in map(hash, pendingStates):
                logging.info("New state %s" % str(transition.dst))
                pendingStates.append(transition.dst)
            else:
                logging.info("Repeated next state %s" % str(transition.dst))

    return FSM(upstream, downstream, edges)

################################################################################

def GenerateFSMAndPrintDOTGraph(graphfilepath, upstream, downstream):
    fsm = GenerateFSM(upstream, downstream)

    for e in fsm.edges:
        logging.info("\t%s" % str(e))

    # Check whether the minimized version of the graph is actually smaller
    minimized_fsm_is_smaller = False
    def f_merge(x, y):
        return x.Superset(y)
    fsm_min = fsm.MinimizedFSM(f_merge).RemoveTransientNodes()
    for e in fsm.edges:
        if e not in fsm_min.edges:
            minimized_fsm_is_smaller = True
            logging.info("Minimized graph is smaller!")
            logging.info("Minimized FSM:")
            for e in fsm_min.edges:
                logging.info("\t%s" % str(e))
            break

    if graphfilepath is not None:
        if minimized_fsm_is_smaller:
            # If the minimized FSM is actually smaller, print the before and
            # the after as two separate graphs
            if graphfilepath == '-':
                f1 = sys.stdout
                f2 = sys.stdout
            else:
                f1 = open(graphfilepath + "/%s_%s_nonminimized.gv" %
                        (upstream.print_name, downstream.print_name), 'w')
                f2 = open(graphfilepath + "/%s_%s_minimized.gv" %
                        (upstream.print_name, downstream.print_name), 'w')
            f1.write(fsm.DOTGraph("Non-minimized"))
            f2.write(fsm_min.DOTGraph("Minimized"))
        else:
            # Otherwise, if the before and the after are the same, print just
            # the one graph
            if graphfilepath == '-':
                f1 = sys.stdout
                f2 = sys.stdout
            else:
                f1 = open(graphfilepath + "/%s_%s_nonminimized.gv" %
                        (upstream.print_name, downstream.print_name), 'w')
                f2 = open(graphfilepath + "/%s_%s_minimized.gv" %
                        (upstream.print_name, downstream.print_name), 'w')
            f1.write(fsm.DOTGraph())
            f2.write(fsm.DOTGraph())

        f = open(graphfilepath + "/%s_%s_minimized.txt" %
                (upstream.print_name, downstream.print_name), 'w')
        f.write(fsm_min.TextRepr())

    # Return the generated graphs
    return fsm_min, fsm if minimized_fsm_is_smaller else None

################################################################################

class ArMORUnitTests(unittest.TestCase):
    def test_minimize(self):
        """from http://en.wikibooks.org/wiki/Digital_Circuits/Optimization
            {frozenset({3, 7}), frozenset({2, 5, 6}), frozenset({0, 1, 4})}"""
        self.assertEqual(
                FSM(Architecture("upstream", MOST({}), {}, set()),
                    Architecture("downstream", MOST({}), {}, set()),
                    {   Transition(0, 7, 0, [0]), Transition(0, 2, 1, [0]),
                        Transition(1, 7, 0, [0]), Transition(1, 5, 1, [0]),
                        Transition(2, 7, 0, [1]), Transition(2, 0, 1, [0]),
                        Transition(3, 0, 0, [1]), Transition(3, 7, 1, [0]),
                        Transition(4, 3, 0, [0]), Transition(4, 6, 1, [0]),
                        Transition(5, 3, 0, [1]), Transition(5, 1, 1, [0]),
                        Transition(6, 3, 0, [1]), Transition(6, 4, 1, [0]),
                        Transition(7, 4, 0, [1]), Transition(7, 3, 1, [0])})
                .MinimizedFSM(min),
                FSM(Architecture("upstream", MOST({}), {}, set()),
                    Architecture("downstream", MOST({}), {}, set()),
                    {
                        Transition(0, 3, 0, [0]), Transition(0, 2, 1, [0]),
                        Transition(2, 3, 0, [1]), Transition(2, 0, 1, [0]),
                        Transition(3, 0, 0, [1]), Transition(3, 3, 1, [0])}))

if __name__ == "__main__":
    unittest.main()

