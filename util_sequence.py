import numpy as np

# this is 1/ζ(2) = 6/pi^2 = 1/(1 + 1/4 + 1/9 + ...).
INVERSE_ZETA_2 = 6 / np.pi**2

def _slowly_increasing_sequence():
    "A sequence that sums to 1, but decreases slowly."

    index = 1
    next_highest_power_of_2 = 2
    lg_nhp2 = 1

    while True:
        yield 2 * INVERSE_ZETA_2 / next_highest_power_of_2 / lg_nhp2**2
        index += 1
        if index == next_highest_power_of_2:
            next_highest_power_of_2 <<= 1
            lg_nhp2 += 1

def slowly_increasing_sequence(add_up_to = 1, smoothness = 128):
    iterator = _slowly_increasing_sequence()

    total = 1
    for i in range(smoothness):
        total -= next(iterator)

    scale = add_up_to / total

    while True:
        yield scale * next(iterator)
