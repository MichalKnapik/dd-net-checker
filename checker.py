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
import itertools

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
        lines = [f.strip() for f in f.readlines() if len(f.strip()) > 0]

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
            nonprimed_to_primed_dict[new_primed_var_name] = new_var_name
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
        self.init_state = None
        self.transitions = None
        self.tau_label = 'tau'

        self.mgr = None
        self.state_bdd_vars = None
        self.state_bdd_vars_nonprimed_to_primed_dict = None
        self.state_to_nonprimed_bdd_encoding = None
        self.action_bdd_vars = None
        self.action_to_bdd_encoding = None
        self.init_state_bdd = None

        # a map from action name (incl. tau_label) to encoding:
        # \/ nonprime_encoded_source /\ action /\ nonprime_encoded_target
        self.transition_relation = {}
        # encoding of \/_x (nonprime_encoded_state_x/\prime_encoded_state_x)
        self.identity = None

    def read_automaton(self, model_fname, actions_list):
        self.known_actions = [self.tau_label]
        self.states, self.transitions = read_model(model_fname)
        self.init_state = self.states[0]

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
            = encode_list_of_labels(self.states, self.mgr, f'{self.name}state', True)
        self.init_state_bdd = self.state_to_nonprimed_bdd_encoding[self.init_state]

        # the special silent action does not sync
        tau_bdd = self.mgr.true

        # prepare the state identity function
        self.identity = self.mgr.false
        for state,bdd_encoding in self.state_to_nonprimed_bdd_encoding.items():
            bdd_primed_encoding = self.mgr.let(self.state_bdd_vars_nonprimed_to_primed_dict, bdd_encoding)
            self.identity = self.identity | (bdd_encoding & bdd_primed_encoding)

        # encode transitions
        for action in self.known_actions:
            self.transition_relation[action] = self.mgr.false

        for source,label,target in self.transitions:
            source_bdd = self.state_to_nonprimed_bdd_encoding[source]
            label_bdd = tau_bdd if label == self.tau_label else self.action_to_bdd_encoding[label]

            target_bdd = self.state_to_nonprimed_bdd_encoding[target]
            target_bdd_primed = self.mgr.let(self.state_bdd_vars_nonprimed_to_primed_dict, target_bdd)

            # debug
            print('@'*112)
            for state,encoding in self.state_to_nonprimed_bdd_encoding.items():
                primed_encoding = self.mgr.let(self.state_bdd_vars_nonprimed_to_primed_dict, encoding)
                print(target, state, primed_encoding == target_bdd_primed)
            # --debug

            self.transition_relation[label] = self.transition_relation[label] | (source_bdd & label_bdd & target_bdd_primed)

# The Network collects automata, builds the global statespace, and has methods for analyzing it.

class Network:

    def __init__(self, automata, actions, name=''):
        self.name = name
        self.automata = automata
        self.actions = actions
        
        self.mgr = None
        self.state_bdd_vars = []
        self.state_bdd_vars_nonprimed_to_primed_dict = {}
        self.action_bdd_vars = None
        self.action_to_bdd_encoding = None

        self.init_state_bdd = None
        self.transition_relation = None

    def get_automata_that_know_action(self, action):
        return [auto for auto in self.automata if action in auto.known_actions]

    def encode_model(self, mgr, action_bdd_vars, action_to_bdd_encoding):
        self.mgr = mgr
        self.action_bdd_vars = action_bdd_vars
        self.action_to_bdd_encoding = action_to_bdd_encoding

        # call encode_model of the underlying automata, build the global initial state 
        self.init_state_bdd = mgr.true
        for automaton in self.automata:
            automaton.encode_model(self.mgr, self.action_bdd_vars, self.action_to_bdd_encoding)
            self.init_state_bdd = self.init_state_bdd & automaton.init_state_bdd
            self.state_bdd_vars.extend(automaton.state_bdd_vars)
            self.state_bdd_vars_nonprimed_to_primed_dict.update(automaton.state_bdd_vars_nonprimed_to_primed_dict)            

        # build the global transition relation

        # big debug
        sync_automata = self.get_automata_that_know_action('c')
        sync_transitions = self.mgr.true
        for auto in sync_automata:
            sync_transitions = sync_transitions & auto.transition_relation['c']

        nonprimed_state_and_action_bdd_var_names = self.state_bdd_vars + self.action_bdd_vars
        primed_state_and_action_bdd_var_names = list(self.state_bdd_vars_nonprimed_to_primed_dict.values()) + self.action_bdd_vars

        sources = self.mgr.quantify(sync_transitions, primed_state_and_action_bdd_var_names, forall=False)

        targetsprimed = self.mgr.quantify(sync_transitions, nonprimed_state_and_action_bdd_var_names, forall=False)
        reversefucker = {}
        for x,y in self.state_bdd_vars_nonprimed_to_primed_dict.items():
            reversefucker[y] = x
        print(reversefucker)
        targets = self.mgr.let(reversefucker, targetsprimed)

        print('sources:')
        self.print_bdd_states_debug(sources)
        print('targets:')
        self.print_bdd_states_debug(targets)
        print('-'*100)

        sys.exit()

        # end big debug

        self.transition_relation = self.mgr.false
        # build transitions w.r.t. self.actions (possibly synchronising)
        for action in self.actions:
            print(f'encoding {action}') # debug
            sync_automata = self.get_automata_that_know_action(action)

            if len(sync_automata) > 0: # if someone reacts to the action...
                # ...take the conjunction of all the transitions for the automata that know the action...
                sync_transitions = functools.reduce(lambda x,y: x&y, \
                                               [auto.transition_relation[action] for auto in sync_automata])
                
                other_automata = [auto for auto in self.automata if auto not in sync_automata]
                # ...with the state identity for the automata that do not know the action
                if len(other_automata) > 0:
                    identity = functools.reduce(lambda x,y: x&y, [auto.identity for auto in other_automata])
                else:
                    identity = self.mgr.true

                sync_transitions = sync_transitions & identity

                # -- debug
                nonprimed_state_and_action_bdd_var_names = self.state_bdd_vars + self.action_bdd_vars
                primed_state_and_action_bdd_var_names = list(self.state_bdd_vars_nonprimed_to_primed_dict.values()) + self.action_bdd_vars
                sources = self.mgr.quantify(sync_transitions, primed_state_and_action_bdd_var_names, forall=False)
                targetsprimed = self.mgr.quantify(sync_transitions, nonprimed_state_and_action_bdd_var_names, forall=False)
                targets = self.mgr.let(self.state_bdd_vars_nonprimed_to_primed_dict, targetsprimed)
                print('sources:')
                self.print_bdd_states_debug(sources)
                print('targets:')
                self.print_bdd_states_debug(targets)
                print('-'*100)

                # wniosek - cos nie tak z targetem!!111, ale sources są ok..
                

                # -- debug end

                
                self.transition_relation = self.transition_relation | sync_transitions

        # add local (tau) transitions
        for auto in self.automata:
            if auto.transition_relation[auto.tau_label] == self.mgr.false:
                continue
            other_automata = [oauto for oauto in self.automata if oauto != auto]
            identity = functools.reduce(lambda x,y: x&y, [auto.identity for auto in other_automata])

            local_transitions = auto.transition_relation[auto.tau_label] & identity
            self.transition_relation = self.transition_relation | local_transitions

    def compute_reachable_space(self, verbose=False):
        # call after encoding model only
        if verbose:
            print('computing reachable statespace')
        
        reachable_states_bdd = self.init_state_bdd
        prev_bdd = self.mgr.false

        nonprimed_state_and_action_bdd_var_names = self.state_bdd_vars + self.action_bdd_vars

        i = 1
        while reachable_states_bdd != prev_bdd:

            prev_bdd = reachable_states_bdd
            next_states_bdd_primed = self.mgr.quantify((reachable_states_bdd & self.transition_relation), \
                                                nonprimed_state_and_action_bdd_var_names, forall=False)

            next_states_bdd_nonprimed = self.mgr.let(self.state_bdd_vars_nonprimed_to_primed_dict, next_states_bdd_primed)

            reachable_states_bdd = reachable_states_bdd | next_states_bdd_nonprimed

            if verbose:
#                print(f'iteration {i}: reached {self.mgr.count(reachable_states_bdd)} state(s)')
                self.print_bdd_states_debug(reachable_states_bdd)

            i += 1

        return reachable_states_bdd

    def print_bdd_states_debug(self, bdd_states):
        # don't try to use it for anything but small debugging
        for gstate in itertools.product(*[auto.states for auto in self.automata]):
            gstate_encoding = functools.reduce(lambda x,y: x&y,[autom.state_to_nonprimed_bdd_encoding[loc] \
                                                                for loc,autom in zip(gstate, self.automata)])
            if gstate_encoding & bdd_states == gstate_encoding:
                print(gstate)
            
    def print_bdd_debug_structs(self):
        print('-'*80)        
        print('self.state_bdd_vars = ' + str(self.state_bdd_vars))
        print('self.state_bdd_vars_nonprimed_to_primed_dict = ' + str(self.state_bdd_vars_nonprimed_to_primed_dict))
        print('self.action_bdd_vars = ' + str(self.action_bdd_vars))
        print('self.action_to_bdd_encoding = ' + str(self.action_to_bdd_encoding))
        print('-'*80)

    def __str__(self):
        return f'Network {self.name} with automata:\n' + '\n'.join(['>> '+str(a) for a in self.automata])

if __name__ == '__main__':
    mgr = _bdd.BDD()

    # read actions (common for all automata)
    actions = read_actions('sync.modgraph')
    action_bdd_var_names, action_bdd_encodings = encode_list_of_labels(actions, mgr, 'act')
    
    # make, read and encode automatons
    wuch = Automaton('wuch')
    wuch.read_automaton('X.modgraph', actions)

    wub = Automaton('wub')
    wub.read_automaton('Y.modgraph', actions)

    wuk = Automaton('wuk')
    wuk.read_automaton('Z.modgraph', actions)    
    
    # make a network
    net = Network([wuch, wub, wuk], actions, 'pulwa')
    net.encode_model(mgr, action_bdd_var_names, action_bdd_encodings)

    print(net)
    net.print_bdd_debug_structs()
    net.compute_reachable_space(verbose=True)
