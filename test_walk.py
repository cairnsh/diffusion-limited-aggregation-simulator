import unittest  
import numpy as np
import walk
import random

class TestCentredBinomialRv(unittest.TestCase):
    def test_centred_binomial_rv(self):
        def f(size):
            j = random.randint(1000, 2000)
            rv = walk.centred_binomial_rv(j, size)
            self.assertTrue(rv.shape == (size,))
            self.assertTrue(np.max(rv) <= j)
            self.assertTrue(np.min(rv) >= -j)
            print("%d %d-binomials: mean %f" % (size, j, np.mean(rv)))
        for size in [1, 2, 1000, 2000]:
            f(size)

    def test_safe_length(self):
        DISTANCE = 200
        length = walk.diagonal_walk_safe_length(DISTANCE, 0.01)
        def f():
            shape = (2, int(np.floor(length)))
            wa = 2*np.random.randint(0, 2, shape) - 1
            wa = np.cumsum(wa, axis=1) # :(
            te = np.linalg.norm(wa, axis=0)
            return np.max(te) <= DISTANCE
        count = 0
        TOTAL = 10000
        for i in range(TOTAL):
            if not f():
                count += 1
        print("fraction:", count / TOTAL)
        self.assertTrue(count < 0.01 * TOTAL)

    def test_rotate(self):
        dtsl = walk.DIAGONAL_TO_SQUARE_LATTICE
        self.assertEqual(dtsl(1, 1), (1, 0))
        self.assertEqual(dtsl(1, -1), (0, -1))
        for i in range(40000):
            x, y = np.random.randint(-20, 21, 2)
            a, b = walk.SQUARE_TO_DIAGONAL_LATTICE(x, y)
            X, Y = walk.DIAGONAL_TO_SQUARE_LATTICE(a, b)
            self.assertEqual((x, y), (X, Y))

    def test_jump_regulator(self):
        regulator = walk.jump_regulator()
        pass

    def test_large_random_position(self):
        pass

    def test_near(self):
        pass

    def test_walk(self):
        handler = walk.walk_handler()
        for i in range(10000):
            handler.step(100)
        for i in range(40):
            handler.moveto(0, 0)
            pos = handler.step(100)
            self.assertTrue(abs(pos[0]) + abs(pos[1]) == 1)
        for i in range(1000):
            pos = handler.pos.copy()
            handler.step_quickly(100)
            pos -= handler.pos
            self.assertTrue(np.linalg.norm(pos) < 100)
