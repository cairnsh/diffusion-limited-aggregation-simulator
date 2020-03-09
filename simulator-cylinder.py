import numpy as np
import numpy.random
import time
import os
import sys
import argparse
import re
from shutil import get_terminal_size
from matplotlib import lines, pyplot as p

import util_sequence

import walk_cylinder as walk

OUTPUTDIR = "output"

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

DELIMITER = re.compile(r"[ ,]")

CYLINDERWIDTH = 1000

DEFAULT_CONFIG = {(i, 0) for i in range(CYLINDERWIDTH)}

def MAP(x, y):
    x %= CYLINDERWIDTH
    return x, y

from util_csv import load_from_csv

class dla_walk:
    def __init__(self, occupied_set_csv=None):
        self.pos = walk.walk_handler()
        self.height = 1
        self.stop_at = set()
        self.occupy = set()

        self.pos.reset()

        if occupied_set_csv:
            self._initialize_from_csv(occupied_set_csv)
        else:
            self._initialize_from_configuration(DEFAULT_CONFIG)

        self.saved_image_counter = 0
        self.timestamp = time.strftime("%Y-%m-%d-%H:%M:%S")

    def _initialize_from_csv(self, csv_file):
        the_configuration = set()
        for site in read_from_csv(csv_file):
            x, y = site
            the_configuration.add((x, y))
        self._initialize_from_configuration(the_configuration)

    def _initialize_from_configuration(self, configuration):
        for x, y in configuration:
            self._fill(x, y)
    
    def _fill(self, x, y):
        self.stop_at.update(walk.NEAR(x, y))
        self.occupy.add((x, y))
        assert y >= -1
        self.height = max(self.height, y + 1)

    def step(self):
        x, y = self.pos.step(self.height)
        if (x,y) in self.stop_at:
            self._fill(x, y)
            self.pos.reset()
        pass
        
    def currently_occupied_sites(self):
        occ = np.array(list(self.occupy))
        return occ[:, 0], occ[:, 1]

    def _save_fig(self, filename, plotsize = 7, scatterdotsize = 2):
        x, y = self.currently_occupied_sites()

        fig, ax = p.subplots(figsize=(plotsize, plotsize))
        r = self.radius
        ax.set_xlim(-r, r)
        ax.set_ylim(-r, r)

        p.scatter(x, y, scatterdotsize)
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

    def plot(self, plotsize = None, scatterdotsize = 2, filename_identifier=""):
        self.make_output_directory()

        filenames = self.generate_next_filenames(filename_identifier)

        #if plotsize is None:
        #    plotsize = plot_scale_decision(self.radius)

        #self._save_fig(filenames['png'], plotsize, scatterdotsize)
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
                occ[j, i] = MAP(x, y) in self.occupy
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
        return len(self.occupy)

    def scale(self):
        return self.pos.scale()

FORMATSTRING = "%9d    %9d particles    %.2f sec/particle    %7d steps/particle    scale: %.0e%s"

def runs(csv, fps, quiet, output_after_every, plot_size, justplot):
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
                    scale = np.sqrt(dwal.scale())
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

    runs(
        csv = fn,
        fps = arg['fps'],
        quiet = arg['quiet'],
        output_after_every = arg['output_after_every'],
        plot_size = arg["plot_size"],
        justplot = arg["plot"] is not None
    )

    # yay
