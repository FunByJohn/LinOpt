from graph import *

'''
    Ford-Fulkerson algorithm (p. 305)

    Step 1. Start with a feasible flow f.

        Here we simply pick the zero flow, since this is a feasible flow.

    Step 2. Search for an augmenting path

        How to do this?

    Step 3. If no augmenting path can be found, the algorithm terminates.

    Step 4. If an augmenting path P is found, then

        (a) If delta(P) < inf, push delta(P) units of flow along P, and go to Sep 2.

        (b) If delta(P) = inf, the algorithm terminates.
'''

def ford_fulkerson(G):
    