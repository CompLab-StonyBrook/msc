#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from itertools import chain, combinations
import pickle


########################
#  Simple hierarchies  #
########################

number_hrcs = {'sp': {('s', 'p')},
               'ps': {('p', 's')}}

person_hrcs = {'123': {('1', '2'), ('2', '3'), ('1', '3')},
               '12|3': {('1', '2'), ('1', '3')},
               '1|23': {('1', '3'), ('2', '3')}}


########################
#  Build Crossproduct  #
########################

def extract_elements(hierarchy: set) -> set:
    return set([element for ranking in hierarchy for element in ranking])


def add_ranking(node1: str, node2: str, person: set, number: set) -> tuple:
    person1, number1 = node1
    person2, number2 = node2
    if node1 != node2 and\
       (person1 == person2 or (person1, person2) in person) and\
       (number1 == number2 or (number1, number2) in number):
        return (node1, node2)


def crossalgebra(person: set, number: set) -> set:
    persons = extract_elements(person)
    numbers = extract_elements(number)
    nodes = set([person+number
                 for person in persons
                 for number in numbers])

    output = set(add_ranking(node, other_node, person, number)
                 for node in nodes
                 for other_node in nodes)
    output.remove(None)
    return output


######################################
#  Compute Algebras with Added Arcs  #
######################################

def reachability_closure(algebra: set) -> set:
    new_algebra = algebra.copy()
    for ranking in algebra:
        for other_ranking in algebra:
            if ranking[1] == other_ranking[0] and\
               ranking[0] != other_ranking[1]:
                new_algebra.add((ranking[0], other_ranking[1]))
    if new_algebra == algebra:
        return new_algebra
    else:
        return reachability_closure(new_algebra)


def algebra_expansions(algebra: set, possible_additions: iter) -> iter:
    for additions in possible_additions:
        # print("adding", additions)
        yield algebra.copy().union(additions)


def arc_addition(algebra: set, person: set, number: set) -> set:
    # compute elements
    persons = extract_elements(person)
    numbers = extract_elements(number)

    # set of all arcs
    nodes = set(p+n for p in persons for n in numbers)
    all_arcs = set((n1, n2) for n1 in nodes for n2 in nodes
                   if n1 != n2)

    # remove arcs that are already in the algebra
    all_arcs -= algebra

    # the powerset of all remaining arcs,
    # i.e. the set of all possible algebra extensions;
    # this is an iterator, so memory usage should be manageable
    possible_additions = chain.from_iterable(combinations(all_arcs, r)
                                             for r in range(len(all_arcs)+1))

    # compute iterator of all algebras with 0 or more added arcs
    algebras = algebra_expansions(algebra, possible_additions)

    # close each algebra under reachability
    closed_algebras = set(frozenset(reachability_closure(algebra))
                          for algebra in algebras)

    # return set of all remaining algebras
    return closed_algebras


######################################
#  Compute Syncretisms from Algebra  #
######################################

def algebra_interpretation(algebra: set) -> tuple:
    presentation = ['1s', '2s', '3s', '1p', '2p', '3p']
    syncretism = [0]
    for cell in presentation[1:]:
        threshold = len(syncretism)
        for pos in range(threshold):
            other_cell = presentation[pos]
            if (other_cell, cell) in algebra and\
               (cell, other_cell) in algebra:
                value = syncretism[pos]
                syncretism.append(value)
                break
        if len(syncretism) == threshold:
            syncretism.append(max(syncretism)+1)
    return tuple(syncretism)


def number_to_letter(syncretism: tuple) -> tuple:
    return tuple([chr(65 + number) for number in syncretism])


def syn_patterns(algebras: set) -> set:
    return set(number_to_letter(algebra_interpretation(algebra))
               for algebra in algebras)


def all_syncretisms() -> dict:
    syncretisms = {}
    for num_label, num_hier in number_hrcs.items():
        for pers_label, pers_hier in person_hrcs.items():
            algebra = crossalgebra(pers_hier, num_hier)
            new_algebras = arc_addition(algebra, pers_hier, num_hier)
            syncretisms[(pers_label, num_label)] = \
                syn_patterns(new_algebras)
    syncretisms['total'] = set().union(*syncretisms.values())
    return syncretisms


try:
    syncretisms = pickle.load(open("syncretisms.p", "rb"))
except:
    syncretisms = all_syncretisms()
    pickle.dump(syncretisms, open("syncretisms.p", "wb"))
