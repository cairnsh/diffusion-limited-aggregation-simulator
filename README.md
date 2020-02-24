A faster DLA simulator that can deposit about 3 particles per second and does not get stuck walking around indefinitely.
Walks are accelerated by taking many steps at once, using Chernoff bounds to make it very unlikely that any collisions are missed.

It makes some approximations but it's pretty close to an exact walk on the lattice.
