import numpy as np
import numpy.random
import util_sequence

PARTICLES_START_OUTSIDE = 1e4

PARTICLES_RESET_OUTSIDE = 4e4

LONG_JUMP_MIN_DISTANCE = 20

ALLOWED_PROBABILITY_OF_ERROR = 1e-6

def centred_binomial_rv(m, size):
    return 2 * np.random.binomial(m, 0.5, size) - m

"""
We can walk for time pretty close to R^2 and still probably not leave
the circle of radius R.  More precisely, the probability that a simple
random walk of length T will go outside R is bounded above by
    8 * exp(-R^2/2T).

A simple random walk in two dimensions can be thought of as a sum of
two independent simple random walks in one dimension X_1 and X_2:
x_n, y_n = DIAGONAL_TO_SQUARE_LATTICE(X_1, X_2)

The probability that a simple random walk starting at (0, 0) will leave
the diamond |x| + |y| <= R in T steps is the same as the probability
that one of the walks X_1, X_2 will leave the interval [-R, R].
This is bounded above by 16 * exp(-R^2 / 2T).

That diamond is contained in the circle of radius R, so we have the
same bound for the probability of leaving the circle.

So, we want to have T <= 0.5 * radius**2 / np.log(16/p).
"""

def _diagonal_walk_safe_length(radius, p):
    return 0.5 * radius**2 / np.log(16 / p)

def DIAGONAL_TO_SQUARE_LATTICE(x, y):
    "Rotate pi/4 radians clockwise and scale by 1/sqrt(2)."
    return (x+y) >> 1, (y-x) >> 1

def SQUARE_TO_DIAGONAL_LATTICE(x, y):
    "The opposite of DIAGONAL_TO_SQUARE_LATTICE."
    return x-y, y+x

def walkfor(length):
    "Do a simple random walk for time <length>."
    x, y = centred_binomial_rv(length, 2)
    return DIAGONAL_TO_SQUARE_LATTICE(x, y)

class jump_regulator:
    def __init__(self, p=ALLOWED_PROBABILITY_OF_ERROR):
        self.probabilities = util_sequence.slowly_increasing_sequence(add_up_to = p)

    def get_allowed_length(self, radius):
        length = _diagonal_walk_safe_length(radius, next(self.probabilities))
        if length < 1:
            return 1
        else:
            return int(np.floor(length))

    def jump(self, radius):
        length = self.get_allowed_length()
        try:
            return walkfor(length)
        except OverflowError:
            print("Overflow error: couldn't do a binomial of length", length)
            print("We were trying to jump inside distance", distance)
            raise

def large_random_position():
    while True:
        x, y = np.random.normal(0, PARTICLES_START_OUTSIDE, 2)

        if x**2 + y**2 < PARTICLES_START_OUTSIDE**2:
            continue

        x, y = int(x), int(y)

        if (x + y) & 1:
            continue

        return np.array([x, y])


def NEAR(x, y):
    return {(x+1, y+1), (x+1, y-1), (x-1, y+1), (x-1, y-1)}

class walk_handler:
    def __init__(self):
        self.jump = jump_regulator()
        self.reset()

    def reset(self):
        self.pos = large_random_position()

    def step_slowly(self):
        self.pos += 2 * np.random.randint(0, 2, 2) - 1

    def step_quickly(self, radius):
        self.pos += self.jump.jump(radius)

    #def step(self, radius_to_walk_slowly):
    def step(self, radius_to_walk_slowly, margin):
        r = np.sqrt(self.pos[0]**2 + self.pos[1]**2)
        if r > PARTICLES_RESET_OUTSIDE:
            self.reset()
        elif r > radius_to_walk_slowly + margin:
            self.step_quickly(r - radius_to_walk_slowly - 0.5)
        else:
            self.step_slowly()
        return self.pos
