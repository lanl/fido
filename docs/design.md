Design Tradeoffs
====

There are numerous tradeoffs that could be made in how one drives a
multi-objective optimizer and couples it to a pde solver for evaluation of some objective function.  One could:
- launch subprocess and communicate with the solver via stdin/stdout
- require the solver to specify an interface such that it can be used as a library for repeated simulations.

The first option could add a significant amount of overhead due to the repeated launching of many subprocess.  It would also require some way to modify the solver's input file and write it to disk in various `tmp` directories so the input files don't get clobbered when doing parallel runs.

The second option seems more restrictive but is actually more general, since one could implement the first option in terms of the second but not vice/versa.

What should the interface look like?  Ideally, we would like to initialize things once and not be required to touch the filesystem each time we run a simulation.  This suggests some kind of C-api where we get back some kind of context pointer that we pass around.  

Say our optimization requires $N$ objective functions $\{O_0, O_1, \ldots,O_{N-1}\}$.  Each objective function may require a different number of simulation runs and contexts.  For each objective function, $O_i$, the $q^{th}$ evaluation, $O^q_i$, requires $\{(C_i, \Delta^q_i)_k\}$ simulation contexts. $C_{i,k}$ is the *constant* portion of the context, such as grid size, cfl, initial conditions, etc. and is independent of the iteration number, $q$. $\Delta^q_{i,k}$ is the portion of the context that changes with each iteration.

This suggests that each objective function should own or be associated with it's various contexts.  If we assume that they are not too large (or too numerous) then they can all be pre-allocated in the context.  That could also be a per-wrapper decision but makes sense for `shoccs`.

Either way, there will need to be an `update_context` function.  In the pre-allocated case, it will only need $\Delta^q_{i,k}$.  In the other, it will need $(C_i, \Delta^q_{i,k})$
