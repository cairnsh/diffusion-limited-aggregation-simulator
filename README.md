A faster DLA simulator that can deposit about 3 particles per second and does not get stuck walking around indefinitely.
Walks are accelerated by taking many steps at once, using Chernoff bounds to make it very unlikely that any collisions are missed.

It makes some approximations but it's pretty close to an exact walk on the lattice.

Usage:
```
usage: simulator.py [-h] (--start | --continue CONTINUE | --plot PLOT) [--fps FPS] [--quiet] [--output_after_every OUTPUT_AFTER_EVERY]
                    [--plot-size PLOT_SIZE]

Simulate DLA in an accelerated and slightly approximate way

optional arguments:
  -h, --help            show this help message and exit
  --start               Start the simulation
  --continue CONTINUE   Continue the simulation from a CSV list of sites
  --plot PLOT           Just plot a list of sites
  --fps FPS             Frames per second for the lifelike ASCII graphics
  --quiet               Turn off the lifelike ASCII graphics
  --output_after_every OUTPUT_AFTER_EVERY
                        The number of simulation steps between plots.
  --plot-size PLOT_SIZE
                        The size of the plots in matplotlib inches, which are 100 pixels.
```

Typical usage: `python simulator.py --start`

Some screenshots. The simulator displays some Unicode graphics and statistics when it's running:
```
                                  ▀█
                                 ▄▄█▀▀
                                   ██▄▄
                                    ▄█▄▄▄▄
                            ▄▄   ▄ ▄▄█  ▀
                             █▀▀█▀▀▀██
                                ▀  ▄ █ ▄
                                  ▀█▀▀▀▀▀
                                  ▀█
                                █▄ █ ▄   ▄▄  ▀██▄▄█
                            ▄    ▀▀█▀█  ▄█▄▄▄█▀ ▀
                           ▀█▄▄▄  ▀█▄██▄██   ▄   ▄ ▄
                          ▄▄▄▄████▀▀▀  █▄ ▄█▄█▄▄▄█▀▀
                             █   ▀     ▀▀██▀
                            ▀█          ▀█▀
                           ▀▀▀       █▀▀▀█▄▄    ▄
                                    ▄█   ███▀█▀██▀
                                     ▀   ▄▄█▄▄  ██▀█▀█▄
                                          ▄██▄  ▀    ▀
                                        ▀▀▀▄█▄▄
                                            █▄▄


       70          238 particles    0.15 sec/particle
```

It also outputs plots every few steps. (By default, every 1,000,000 steps.)

![An image output by the DLA simulator.](https://github.com/cairnsh/diffusion-limited-aggregation-simulator/example_plot.png)
