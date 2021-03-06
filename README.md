A faster DLA simulator that can deposit 1-3 particles per second and does not get stuck walking around indefinitely.
Walks are accelerated by taking many steps at once, using Chernoff bounds to make it very unlikely that any collisions are missed. (The total probability of a mistake in any individual run is at most .000001.)

It creates particles by rounding a Gaussian and keeping only points with radius > 10000, and it kills particles at radius > 40000. Aside from those approximations it should be roughly an exact simulation on the lattice with probability 1 - 1e-6, up to the accuracy of numpy's `numpy.random.binomial`.

It uses ideas from "Internal DLA: Efficient Simulation of a Physical Growth Model" (https://people.mpi-inf.mpg.de/~kbringma/paper/2014ICALP.pdf).

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

Typical usage:
```
python simulator.py --start
python simulator.py --start --fps=6
python simulator.py --start --fps=29.97
python simulator.py --start --quiet &
python simulator.py --continue output/dla-final-2020-02-25-14:47:02-001.csv --output_after_every 1000000
python simulator.py --plot output/dla-final-2020-02-25-14:47:02-001.csv --plot-size 20
```

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

It also outputs plots every few steps. (By default, every 10,000,000 steps.)

The names of the images are `output/dla-(timestamp at start of run)-(number of image).png`. It writes out a CSV file containing the occupied sites together with each plot.

The program will make one final plot if you interrupt it from the keyboard, which will be called `dla-final-(timestamp)-(number).png`.

The plots become larger as the occupied set grows.![An image output by the DLA simulator.](https://github.com/cairnsh/diffusion-limited-aggregation-simulator/blob/master/example-plot.png)

Note: Y. L. Loh apparently wrote a much faster simulator for "Bias-free simulation of diffusion-limited aggregation on a square lattice." (https://arxiv.org/pdf/1407.2586.pdf)
