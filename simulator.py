import numpy as np

import numpy.random

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
    
    def get_allowed_p(self):
        return self.p * next(self.seq)
    
    @staticmethod
    def get_allowed_length(p, radius):
        # p <= 16 * exp(-(R/sqrt(2))^2 / 2T)
        # T = radius**2 / 4 / ln(8/p)
        # extra factor of sqrt(2)
        return radius**2 / 4 / np.log(16 / p)
    
    def jump(self, distance):
        length = regulatedjump.get_allowed_length(
            self.get_allowed_p(),
            distance
        )
        try:
            return centred_binomial_rv(int(length - 0.5), 2)
        except OverflowError:
            print(distance, length)
            raise
        
def large_random_pos():
    while True:
        x, y = np.random.normal(0, 10000, 2)
        if max(x, y, -x, -y) > 10000:
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

class dla_walk:
    def __init__(self, starting_occupied=None):
        self.regj = regulatedjump()
        self.stepgenerator = stepgenerator()
        self.radius = 1
        self.reset()
        self.stopat = set()
        self.occ = set()
        if starting_occupied is None:
            starting_occupied = {(0, 0)}
        for site in starting_occupied:
            self.fillin(site[0], site[1])
    
    def reset(self):
        self.pos = large_random_pos()
    
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
        if r > 1e5: # just reset
            self.log("reset from %d %d" % (x, y))
            self.pos = large_random_pos()
            pass
        elif r > self.radius + 20: # long jump
            self.log("long jump %d %d" % (x, y))
            self.pos += self.regj.jump(int(r - self.radius - 0.5))
            pass
        else:
            self.log("walk from %d %d" % (x, y))
            self.walk()
        return self.pos

    def plot(self, nam, ff):
        from matplotlib import pyplot as p
        occ = np.array(list(self.occ))
        x, y = occ[:, 0], occ[:, 1]
        x, y = (x + y) >> 1, (x - y) >> 1
        p.figure(figsize=(ff, ff))
        p.scatter(x, y, s = 2)
        import os
        fn = "dla-%s.png" % nam
        i = 0
        while os.path.exists(fn):
            i += 1
            fn = "dla-%s-%d.png" % (nam, i)
        print("saving to %s" % fn)
        p.savefig(fn)
        csv = open(fn + ".csv", "w")
        for j in range(x.shape[0]):
            csv.write("%d %d\n" % (x[j], y[j]))
        csv.close()

    def ascii(self, insert):
        CHARS = " ▄▀█"
        import shutil
        size = shutil.get_terminal_size()
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

if __name__ == "__main__":
    import time
    dwal = dla_walk()
    t = time.time()
    try:
        i = 0
        while True:
            for j in range(100000):
                dwal.step()
            rat = (time.time() - t) / max(1, len(dwal.occ))
            dwal.ascii("%9d    %9d particles    %.2f sec/particle" % (i, max(0, len(dwal.occ)), rat))
            if i % 100 == 0 and i > 0:
                dwal.plot("%06d" % i, 7)
            i += 1
    finally:
        dwal.plot("final", 7)
    print()
    print("%d occupied sites" % len(dwal.occ))
