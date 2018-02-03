#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from itertools import chain, combinations, product
from pprint import pprint
import pickle
import csv


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


##################
#  Analyze Data  #
##################

def read_data(filename: str='data.csv') -> list:
    patterns = []
    with open('data.csv') as datafile:
        data = csv.reader(datafile)
        for row in data:
            patterns.append(tuple(row[2:8]))
    return patterns


def overgeneration(algebra: set, data: list) -> list:
    return [pattern for pattern in algebra if pattern not in data]


def undergeneration(algebra: set, data: list) -> list:
    return [pattern for pattern in data if pattern not in algebra]


def all_generation(function=overgeneration) -> dict:
    data = list(set(read_data()))
    return {key: set(function(syncretisms[key], data))
            for key in syncretisms.keys()}


def show_generation(function=overgeneration):
    generated = all_generation(function=function)

    for key, val in generated.items():
        print(key, ":", len(val))
        pprint(sorted(list(val)))
        print("--------------")


def all_overgeneration() -> dict:
    return all_generation(function=overgeneration)


def all_undergeneration() -> dict:
    return all_generation(function=undergeneration)


def show_overgeneration():
    return show_generation(function=overgeneration)


def show_undergeneration():
    return show_generation(function=undergeneration)


def show_generators(algebras: dict, patterns: list) -> dict:
    overview = {}
    for pattern in patterns:
        overview[pattern] = set()
        for algebra, generated in algebras.items():
            if pattern in generated and algebra != 'total':
                overview[pattern].add(algebra)
    return overview


#############################
#  Monotonicity Calculator  #
#############################

base_algebras = {(plab, nlab): reachability_closure(crossalgebra(p, n))
                 for plab, p in person_hrcs.items()
                 for nlab, n in number_hrcs.items()}


def relabel_algebra(algebra, syncretism):
    presentation = ['1s', '2s', '3s', '1p', '2p', '3p']
    labeling = dict(zip(presentation, syncretism))
    return {(labeling[a], labeling[b])
            for (a, b) in algebra}


def test_monotonicity(algebra, syncretism):
    algebra = relabel_algebra(algebra, syncretism)
    for a, b in algebra:
        if a != b:
            for c, d in algebra:
                if (a, b) == (d, c):
                    return False
    return True


def is_monotonic(algebras, syncretism):
    return [key for key, val in algebras.items()
            if test_monotonicity(val, syncretism)]


def show_monotonicity(algebras=base_algebras, data=read_data()):
    return {syncretism: is_monotonic(algebras, syncretism)
            for syncretism in data}


def all_patterns(cells=6):
    letters = [[chr(65 + number) for number in range(cells)]
               for _ in range(cells)]
    return [i for i in product(*letters)]


try:
    patterns = pickle.load(open("patterns.p", "rb"))
except:
    patterns = all_patterns()
    pickle.dump(patterns, open("patterns.p", "wb"))


def all_monotonicity():
    return show_monotonicity(base_algebras, patterns)


def letter_to_number(xs):
    if not xs:
        return xs

    ys = [1]
    for pos in range(1, len(xs)):
        x = xs[pos]
        for left in range(pos):
            if x == xs[left]:
                ys.append(ys[left])
                break
        else:
            ys.append(max(ys) + 1)
    return ys


def isomorphic(xs, ys):
    return letter_to_number(xs) == letter_to_number(ys)
        

def monotonicity_overgeneration():
    all_monotonic = {key: val
                     for key, val in all_monotonicity().items()
                     if val}

    unique_monotonic = {}
    for key, val in all_monotonic.items():
        for key2, val2 in unique_monotonic.items():
            if isomorphic(key, key2):
                break
        else:
            unique_monotonic[key] = val

    attested_monotonic = show_monotonicity()

    return {key: val
            for key, val in unique_monotonic.items()
            for key2 in attested_monotonic
            if not isomorphic(key, key2)}

try:
    monotonic_overgen = pickle.load(open("monotonic_overgen.p", "rb"))
except:
    monotonic_overgen = monotonicity_overgeneration()
    pickle.dump(monotonic_overgen, open("monotonic_overgen.p", "wb"))
