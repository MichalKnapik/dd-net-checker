#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import math

try:
    import dd
except ImportError:
    print("Cannot import DD module! Exiting.")
    sys.exit()

try:
    import dd.cudd as _bdd
    print("Loaded CUDD from DD.")
except ImportError as e:
    print("Failed to load CUDD, but I really need it (autoref won't do). Exiting.")
    sys.exit()

# Helper functions.

def read_actions(sync_fname):
    with open(sync_fname) as f:
        actions = [f.strip() for f in f.readlines()]
        return actions

def read_model(model_fname):
    with open(model_fname) as f:
        lines = [f.strip() for f in f.readlines()]

        if lines[0] != 'states':
            print("Model file error: no leading 'states'.")
            sys.exit()
        try:
            state_sec_end = lines.index('transitions')
        except ValueError:
            print("Model file error: no 'transitions'.")
            sys.exit()
        states = lines[1:state_sec_end]

        transitions = [l[1:-1] for l in lines[state_sec_end+1:]]
        transitions = [t.split(',') for t in transitions]
        transitions = [tuple([l.strip() for l in t]) for t in transitions]

        return states,transitions

def encode_list_of_labels(label_list, bdd_mgr, variable_prefix=''):
    """Generate new bdd variables for label_list and return a dict from
    labels in label_list to their bdd encodings using these new bdd vars."""

    needed_vars = math.ceil(math.log2(len(label_list)))
    for i in range(needed_vars):
        mgr.add_var(variable_prefix + str(i))

    encodings_dict = {}
    i = 0
    for label in label_list:
        raw_encoding = ('0' * needed_vars + bin(i)[2:])[-needed_vars:]
        bool_encoding = [bool(int(bit)) for bit in raw_encoding]

        bdd_encoding = bdd_mgr.true        
        for j in range(len(bool_encoding)):
            if bool_encoding[j]:
                bdd_encoding = bdd_encoding & bdd_mgr.var(variable_prefix + str(j))
            else:
                bdd_encoding = bdd_encoding & (~bdd_mgr.var(variable_prefix + str(j)))
        
        encodings_dict[label] = bdd_encoding
        i += 1

    return encodings_dict

# The Automaton collects all the raw data and bdd structures.

class Automaton:

    def __init__(self, name=''):
        self.name = name
        self.mgr = None
        self.known_actions = None
        self.states = None
        self.transitions = None

    def read_automaton(self, model_fname, actions_list):
        self.known_actions = []
        self.states, self.transitions = read_model(model_fname)

    def __str__(self):
        ret = f'automaton {self.name}\n'
        ret += '' if self.known_actions is None or len(self.known_actions) == 0 else 'with known actions: ' + ' '.join(self.known_actions)

        return ret

if __name__ == '__main__':
    mgr = _bdd.BDD()
    actions = read_actions('tests/case_w4d3c2/sync.modgraph')
    action_label_to_bdd_encoding = encode_list_of_labels(actions, mgr, 'act')

    wuch = Automaton('wuch')
    wuch.read_automaton('tests/case_w4d3c2/AND1.modgraph', actions)
    print(wuch)
