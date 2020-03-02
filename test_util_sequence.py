import unittest
from util_sequence import slowly_increasing_sequence

class TestSlowlyIncreasingSequence(unittest.TestCase):
    def test_nonnegative(self):
        z = slowly_increasing_sequence()
        for i in range(1000000):
            self.assertTrue(next(z) >= 0)

        z = slowly_increasing_sequence(128)
        for i in range(1000000):
            self.assertTrue(next(z) >= 0)

    def test_sum_is_less_than_one(self):
        z = slowly_increasing_sequence(1, 0)
        total = 0
        for i in range(1000000):
            total += next(z)
            self.assertTrue(total < 1)
        print("Sum of the first 1000000 elements:", total)

        z = slowly_increasing_sequence(1, 128)
        total = 0
        for i in range(1000000):
            total += next(z)
            if total >= 1:
                print("failure at", i, total)
            self.assertTrue(total < 1)
        print("Sum of the first 1000000 elements, dropping the first 128 and rescaling:", total)

if __name__ == "__main__":
    unittest.main()
