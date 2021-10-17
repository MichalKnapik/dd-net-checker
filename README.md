# Model checker for networks of synchronizing automata

This is an extremely rudimentary model checker for networks of automata in modgraph format. Currently only supports building statespace and estimating its size.
It uses CUDD BDD library via [Tulip-control Python bindings](https://github.com/tulip-control/dd).

## Example

Run:
```sh
./checker.py examples/small_system/sync examples/small_system/*.modgraph --verbose 1
```

to get:
```
Loaded CUDD from DD.
Read and encoded network of 3 automata.
Network net with automata:
>> automaton wuch0 with: 
states: ['1', '2', '3']
known actions: ['tau', 'a', 'b', 'c']
transitions: [('1', 'a', '2'), ('2', 'b', '3'), ('3', 'c', '1')]
>> automaton wuch1 with: 
states: ['1', '2', '3']
known actions: ['tau', 'a', 'c']
transitions: [('1', 'a', '2'), ('2', 'c', '3'), ('3', 'tau', '1')]
>> automaton wuch2 with: 
states: ['1', '2', '3']
known actions: ['tau', 'b', 'c']
transitions: [('1', 'b', '2'), ('2', 'c', '3'), ('3', 'tau', '1')]
** Computing reachable statespace. This might take a while... **
iteration 1: reached 1.0 state(s)
iteration 2: reached 2.0 state(s)
iteration 3: reached 3.0 state(s)
iteration 4: reached 4.0 state(s)
iteration 5: reached 6.0 state(s)
iteration 6: reached 7.0 state(s)
Reached 7.0 state(s) and 9.0 transition(s).
Done.
```

Now we know the size of the statespace, (approximately) the number of transitions, and how the statespace grows as the fixpoint is computed. Play with *verbose* to get more information about the system.

## The small system:

**see examples/small_system/sync:**

```
a
b
c
```

**see examples/small_system/X.modgraph:**

```
states
1
2
3
transitions
(1, a ,2)
(2, b ,3)
(3, c ,1)
```

**see examples/small_system/Y.modgraph:**

```
states
1
2
3
transitions
(1, a ,2)
(2, c ,3)
(3, u ,1)
```

**examples/small_system/Z.modgraph:**

```
states
1
2
3
transitions
(1, b ,2)
(2, c ,3)
(3, l ,1)
```
