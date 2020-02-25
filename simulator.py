import numpy as np
import numpy.random
import time
import os
import sys
import argparse
import re
from shutil import get_terminal_size
from matplotlib import lines, pyplot as p

OUTPUTDIR = "output"

PARTICLES_START_OUTSIDE = 1e4

PARTICLES_RESET_OUTSIDE = 4e4

LONG_JUMP_MIN_DISTANCE = 20

def plot_scale_decision(plotradius):
    "decide what size a plot should be given the radius of the occupied set"
    smallradius = 10
    largeradius = 1000
    smallscale = 4
    largescale = 24
    if plotradius < smallradius:
        return smallscale
    if plotradius > largeradius:
        return largescale
    # interpolate in log
    return int(smallscale * (plotradius / smallradius) ** (np.log(largescale / smallscale) / np.log(largeradius / smallradius)))

"""
    NOTE: our random walk has diagonal steps.
    We always start at a position (x, y) with x + y even.
    This makes the jump probability a bit easier to think about,
    because the two coordinates are independent.

    When we draw the random walk in .plot() or .ascii(), we correct for this
    with the map (x, y) -> ((x+y)/2, (x-y)/2).
"""

invzeta = 6 / np.pi**2

def _slow():
    "a sequence that sums to 1, but decreases slowly"
    i, l, ll = 1, 2, 1
    while True:
        yield 2 * invzeta / l / ll**2
        i += 1
        if i == l:
            l <<= 1
            ll += 1

def _drop_first_and_rescale(iterator, drop):
    scale = 1
    for i in range(drop):
        scale -= next(iterator)
    while True:
        yield next(iterator) / scale

def centred_binomial_rv(m, size):
    return 2 * np.random.binomial(m, 0.5, size) - m

class regulatedjump:
    """
    We can walk for time pretty close to R^2 and still probably not leave
    the circle of radius R.  More precisely, the probability that a simple
    random walk of length T will go outside R is bounded above by
        8 * exp(-R^2/2T).

    The probability that our diagonal random walk will leave the
    square [-R, R]^2 is <= 16 * exp(-R^2 / 2T), so the probability
    that it will leave the circle of radius R is <=
        16 * exp(-R^2 / 4T).
    
    We choose the jump lengths so that the total probability we
    ever miss a departure is less than p, using the above sequence
    of numbers that add to 1.
    """
    
    def __init__(self, p=1e-6):
        self.p = p
        self.seq = _slow()
        self.seq = _drop_first_and_rescale(self.seq, 128)
    
    def get_allowed_p(self):
        return self.p * next(self.seq)
    
    @staticmethod
    def get_scale(p):
        # p <= 16 * exp(-(R/sqrt(2))^2 / 2T)
        # T = radius**2 / 4 / ln(8/p)
        # extra factor of sqrt(2)
        return 0.25 / np.log(16 / p)
    
    def jump(self, distance):
        p = self.get_allowed_p()
        self.scale = regulatedjump.get_scale(p)
        length = distance**2 * self.scale
        length = int(np.floor(length))
        if length < 1:
            length = 1
        try:
            return centred_binomial_rv(length, 2)
        except OverflowError:
            print("Overflow error: couldn't do a binomial of length", length)
            print("We were trying to jump inside distance", distance)
            raise
        
def large_random_pos():
    while True:
        x, y = np.random.normal(0, PARTICLES_START_OUTSIDE, 2)
        if np.sqrt(x**2 + y**2) > PARTICLES_START_OUTSIDE:
            x, y = int(x), int(y)
            if (x + y) & 1 == 0:
                return np.array([x, y], dtype=np.int64)

def stepgenerator():
    COUNT = 100
    pass
    while True:
        steps = 2 * np.random.randint(0, 2, (COUNT, 2)) - 1
        for i in range(COUNT):
            yield steps[i, :]

def DIAGONAL_TO_SQUARE_LATTICE(x, y):
    "Rotate pi/4 radians clockwise and scale by 1/sqrt(2)."
    return (x+y) >> 1, (y-x) >> 1

def SQUARE_TO_DIAGONAL_LATTICE(x, y):
    "The opposite of DIAGONAL_TO_SQUARE_LATTICE."
    return x-y, y+x

DELIMITER = re.compile(r"[ ,]")

def read_from_csv(f):
    o = open(f)
    data = []
    for line in o:
        data += [map(int, DELIMITER.split(line.strip()))]
    return data

class dla_walk:
    def __init__(self, occupied_set_csv=None):
        self.regj = regulatedjump()
        self.stepgenerator = stepgenerator()
        self.radius = 1
        self.reset()
        self.stopat = set()
        self.occ = set()
        if occupied_set_csv:
            for site in read_from_csv(occupied_set_csv):
                x, y = SQUARE_TO_DIAGONAL_LATTICE(*site)
                self.fillin(x, y)
        else:
            self.default_starting_configuration()
        self.saved_image_counter = 0
        self.timestamp = time.strftime("%Y-%m-%d-%H:%M:%S")
    
    def reset(self):
        self.pos = large_random_pos()

    def default_starting_configuration(self):
        # just one occupied site at the origin
        self.fillin(0, 0)
    
    def fillin(self, x, y):
        self.stopat.update({(x+1, y+1), (x+1, y-1), (x-1, y+1), (x-1, y-1)})
        self.occ.add((x, y))
        self.radius = max(self.radius, np.sqrt(x**2 + y**2) + 1)
        
    def walk(self):
        self.pos += 2 * np.random.randint(0, 2, 2) - 1
        if tuple(self.pos) in self.stopat:
            self.fillin(self.pos[0], self.pos[1])
            self.reset()
    
    def log(self, txt):
        if False:
            if getattr(self, "logfile", None) is None:
                self.logfile = open("log", "w")
            self.logfile.write(txt)
            self.logfile.write("\n")
            self.logfile.flush()
    
    def step(self):
        x, y = self.pos
        r = np.sqrt(x**2 + y**2)
        if r > PARTICLES_RESET_OUTSIDE: # just reset
            self.log("reset from %d %d" % (x, y))
            self.pos = large_random_pos()
            pass
        elif r > self.radius + LONG_JUMP_MIN_DISTANCE: # long jump
            self.log("long jump %d %d" % (x, y))
            self.pos += self.regj.jump(int(r - self.radius - 0.5))
            pass
        else:
            self.log("walk from %d %d" % (x, y))
            self.walk()
        return self.pos

    def currently_occupied_sites(self):
        """Return the list of currently occupied sites.

        The result is rotated into the square lattice
        by the map x, y -> (x + y) >> 1, (x - y) >> 1."""
        occ = np.array(list(self.occ))
        x, y = occ[:, 0], occ[:, 1]
        return DIAGONAL_TO_SQUARE_LATTICE(x, y)

    def _save_fig(self, filename, plotsize = 7, scatterdotsize = 2):
        x, y = self.currently_occupied_sites()

        fig, ax = p.subplots(figsize=(plotsize, plotsize))
        r = self.radius / np.sqrt(2) * 1.1
        ax.set_xlim(-r, r)
        ax.set_ylim(-r, r)

        p.scatter(x, y, scatterdotsize)

        """
        def addline(x1, y1, x2, y2):
            fig.lines.extend([ lines.Line2D([x1, x2], [y1, y2], transform=ax.transData, figure=fig) ])

        for x, y in self.occ:
            X, Y = DIAGONAL_TO_SQUARE_LATTICE(x, y)
            if (x+1, y+1) in self.occ:
                addline(X, Y, X+1, Y)
                print("line", X, Y, X+1, Y)
                pass
            if (x+1, y-1) in self.occ:
                addline(X, Y, X, Y-1)
                print("line", X, Y, X, Y-1)
        """
                
        p.savefig(filename)
        p.close(fig)

    def _save_csv(self, filename):
        x, y = self.currently_occupied_sites()

        csv = open(filename, "w")
        for j in range(x.shape[0]):
            csv.write("%d,%d\n" % (x[j], y[j]))
        csv.close()

    def generate_next_filenames(self, filename_identifier=""):
        if not os.path.isdir(OUTPUTDIR):
            raise Exception("output directory doesn't exist, but we should have made it in make_output_directory")
        self.saved_image_counter += 1
        fn = "%s/dla%s-%s-%03d" % (OUTPUTDIR, filename_identifier, self.timestamp, self.saved_image_counter)
        return {
                "png": fn + ".png",
                "csv": fn + ".csv"
        }

    def make_output_directory(self):
        try:
            os.mkdir(OUTPUTDIR)
        except FileExistsError:
            pass

    #def plot(self, plotsize = 7, scatterdotsize = 2):
    def plot(self, plotsize = None, scatterdotsize = 2, filename_identifier=""):
        self.make_output_directory()

        filenames = self.generate_next_filenames(filename_identifier)

        if plotsize is None:
            plotsize = plot_scale_decision(self.radius)

        self._save_fig(filenames['png'], plotsize, scatterdotsize)
        self._save_csv(filenames['csv'])
        return filenames['png']

    def ascii(self, insert):
        CHARS = " ▄▀█"
        size = get_terminal_size()
        COLS = size[0] - 1
        ROWS = 2*size[1]
        occ = np.zeros((COLS, ROWS))
        for j in range(COLS):
            x = j - COLS//2
            for i in range(ROWS):
                y = i - ROWS//2
                occ[j, i] = (x+y, x-y) in self.occ
        offset = ROWS & 1
        end = (ROWS + 1) // 2
        put = []
        for i in range(end):
            l1 = occ[:, 2*i - offset] if 2*i - offset >= 0 else np.zeros(COLS)
            l2 = occ[:, 2*i - offset + 1]
            string = [None] * COLS
            for j in range(COLS):
                ix = (2 if l1[j] else 0) + (1 if l2[j] else 0)
                string[j] = CHARS[ix]
            out = "".join(string)
            if i == end - 1:
                out = insert + out[len(insert):]
            put += [out]
        none = ""
        print("\n" + "\n".join(put), end=none)

    def sitecount(self):
        return len(self.occ)

#FORMATSTRING = "%9d    %9d particles    %.2f sec/particle%s"
FORMATSTRING = "%9d    %9d particles    %.2f sec/particle    %7d steps/particle    scale: %.0e%s"

def walk(csv, fps, quiet, output_after_every, plot_size, justplot):
    dwal = dla_walk(occupied_set_csv = csv)
    try:
        if not justplot:
            starting = dwal.sitecount()
            t = time.time()
            i = 0
            steps = 0
            imageat = steps + output_after_every
            fn = ""
            while True:
                timer = time.time() + 1.0/fps
                while time.time() < timer:
                    steps += 1
                    dwal.step()

                i += 1
                if not quiet:
                    sitecount = dwal.sitecount() - starting or 1
                    elapsed = time.time() - t
                    scale = np.sqrt(dwal.regj.scale)
                    extra = FORMATSTRING % (i, sitecount, elapsed / sitecount, steps / sitecount, scale, fn)
                    #extra = FORMATSTRING % (i, dwal.sitecount(), (time.time() - t) / max(1, (dwal.sitecount() - starting)), fn)
                    dwal.ascii(extra)

                if imageat <= steps:
                    imageat += output_after_every
                    fn = " " * 4 + dwal.plot(plot_size, 2)
    finally:
        dwal.plot(plot_size, 2, "-final")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Simulate DLA in an accelerated and slightly approximate way')

    mode = parser.add_mutually_exclusive_group(required=1)
    mode.add_argument('--start', action='store_true', help="Start the simulation")
    mode.add_argument('--continue', type=str, help="Continue the simulation from a CSV list of sites")
    mode.add_argument('--plot', type=str, help="Just plot a list of sites")
    parser.add_argument('--fps', type=float, default=2, help="Frames per second for the lifelike ASCII graphics")
    parser.add_argument('--quiet', action='store_true', help="Turn off the lifelike ASCII graphics")
    parser.add_argument('--output_after_every', type=int, default=int(1e7), 
            help="The number of simulation steps between plots.")
    parser.add_argument('--plot-size', type=int, help="The size of the plots in matplotlib inches, which are 100 pixels.")

    arg = vars( parser.parse_args(sys.argv[1:]) )

    if arg['fps'] <= 0:
        raise Exception("fps should be > 0, but it is %f" % arg['fps'])

    #print(arg)

    fn = arg["continue"] or arg["plot"]

    walk(
        csv = fn,
        fps = arg['fps'],
        quiet = arg['quiet'],
        output_after_every = arg['output_after_every'],
        plot_size = arg["plot_size"],
        justplot = arg["plot"] is not None
    )

    # yay
