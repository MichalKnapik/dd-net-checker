#!/usr/bin/env python3
# -*- coding: utf-8 -*-

################################################################################
# This is a script for counting the number of states and (approximately) the   #
# number of transitions in a network of synchronizing automata given in a mod- #
# graph file (see tests dir). It uses dd-tulip CUDD BDD bindings.              #
################################################################################

import sys
import math
import functools

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

def encode_list_of_labels(label_list, bdd_mgr, variable_prefix='', generate_primed=False):
    """Generate new bdd variables for label_list and return:
    - if generate_primed = False: a list of var names and dict from labels in 
      label_list to their bdd encodings using these new bdd vars.
    - if generate_primed = True: as above, plus a dict from nonprimed to primed
      var names (for let-substitutions in bdd operations)."""

    needed_vars = math.ceil(math.log2(len(label_list)))
    bdd_var_names = []

    if generate_primed:
        primed_prefix = 'primed'
        nonprimed_to_primed_dict = {}

    for i in range(needed_vars):
        
        new_var_name = variable_prefix + str(i)
        mgr.add_var(new_var_name)
        bdd_var_names.append(new_var_name)

        if generate_primed:
            new_primed_var_name = primed_prefix + new_var_name
            nonprimed_to_primed_dict[new_var_name] = new_primed_var_name
            mgr.add_var(new_primed_var_name)

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

    if not generate_primed:
        return bdd_var_names,encodings_dict
    else:
        return bdd_var_names,encodings_dict,nonprimed_to_primed_dict

# The Automaton collects all the raw data and bdd structures.

class Automaton:

    def __init__(self, name=''):
        self.name = name
        self.known_actions = None
        self.states = None
        self.transitions = None
        self.tau_label = 'tau'

        self.mgr = None
        self.state_bdd_vars = None
        self.state_bdd_vars_nonprimed_to_primed_dict = None
        self.state_to_nonprimed_bdd_encoding = None
        self.action_bdd_vars = None
        self.action_to_bdd_encoding = None
        self.transition_relation = None

    def read_automaton(self, model_fname, actions_list):
        self.known_actions = []
        self.states, self.transitions = read_model(model_fname)

        # note: labels of all the transitions that do not belong to actions_list
        # are replaced with tau_label (i.e., non-synchronized transition).
        for i in range(len(self.transitions)):
            origin,trans_label,target = self.transitions[i]
            if trans_label not in actions_list:
                self.transitions[i] = (origin, self.tau_label, target)
            else:
                self.known_actions.append(trans_label)

    def __str__(self):
        ret = f'automaton {self.name} with: \n'
        ret += 'states: ' + str(self.states)
        ret += '\n' if self.known_actions is None or len(self.known_actions) == 0 else '\nknown actions: ' + str(self.known_actions)
        ret += '\ntransitions: ' + str(list(self.transitions))
        
        return ret

    def print_bdd_debug_structs(self):
        print('-'*80)        
        print('self.state_bdd_vars = ' + str(self.state_bdd_vars))
        print('self.state_bdd_vars_nonprimed_to_primed_dict = ' + str(self.state_bdd_vars_nonprimed_to_primed_dict))
        print('self.state_to_nonprimed_bdd_encoding = ' + str(self.state_to_nonprimed_bdd_encoding))
        print('self.action_bdd_vars = ' + str(self.action_bdd_vars))
        print('self.action_to_bdd_encoding = ' + str(self.action_to_bdd_encoding))
        print('-'*80)
    
    def encode_model(self, mgr, action_bdd_vars, action_to_bdd_encoding):
        self.mgr = mgr
        
        # these encodings are shared for all the automata, so need to be provided externally
        self.action_bdd_vars = action_bdd_vars
        self.action_to_bdd_encoding = action_to_bdd_encoding

        # encode states
        self.state_bdd_vars,self.state_to_nonprimed_bdd_encoding,self.state_bdd_vars_nonprimed_to_primed_dict \
            = encode_list_of_labels(self.states, self.mgr, 'state', True)

        # tau_bdd: negated conjunction of encodings of all the actions, needed to prevent tau from syncing with
        # known actions of other agents
        all_action_bdds = [self.action_to_bdd_encoding[act] for act in self.known_actions]
        tau_bdd = ~functools.reduce(lambda x,y: x|y, all_action_bdds)

        # encode transitions
        self.transition_relation = self.mgr.false
        for source,label,target in self.transitions:

            source_bdd = self.state_to_nonprimed_bdd_encoding[source]
            label_bdd = tau_bdd if label == self.tau_label else self.action_to_bdd_encoding[label]

            target_bdd = self.state_to_nonprimed_bdd_encoding[target]
            target_bdd_primed = self.mgr.let(self.state_bdd_vars_nonprimed_to_primed_dict, target_bdd)

            self.transition_relation = self.transition_relation | (source_bdd & label_bdd & target_bdd_primed)

if __name__ == '__main__':
    mgr = _bdd.BDD()

    # read actions (common for all automata)
    actions = read_actions('tests/case_w4d3c2/sync.modgraph')
    action_bdd_var_names, action_bdd_encodings = encode_list_of_labels(actions, mgr, 'act')
    
    # make, read and encode automaton
    wuch = Automaton('wuch')
    wuch.read_automaton('tests/case_w4d3c2/AND1.modgraph', actions)
    wuch.encode_model(mgr, action_bdd_var_names,action_bdd_encodings)
    
    print(wuch)
#    wuch.print_bdd_debug_structs()
    
