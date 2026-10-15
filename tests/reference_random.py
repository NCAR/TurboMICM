"""Portable reproduction of the reference host's seeded random test data.

The source tests leave default_random_engine/distribution algorithms to the C++
implementation. Retain the supplied double-precision build's minstd_rand,
MT19937, canonical uniform, and cached polar-normal sequences in Python so the
same matrices are tested with either C compiler and on other operating systems.
"""
import math


class Minstd:
    minimum, span = 1, 2147483646

    def __init__(self, seed=1):
        self.state = seed % 2147483647 or 1

    def __call__(self):
        self.state = self.state * 48271 % 2147483647
        return self.state


class MT19937:
    minimum, span = 0, 4294967296

    def __init__(self, seed=12345):
        self.state = [seed & 0xffffffff]
        for i in range(1, 624):
            x = self.state[-1]
            self.state.append((1812433253 * (x ^ (x >> 30)) + i) & 0xffffffff)
        self.index = 624

    def __call__(self):
        if self.index == 624:
            for i in range(624):
                x = (self.state[i] & 0x80000000) | (self.state[(i+1) % 624] & 0x7fffffff)
                self.state[i] = self.state[(i+397) % 624] ^ (x >> 1) ^ (0x9908b0df if x & 1 else 0)
            self.index = 0
        x = self.state[self.index]
        self.index += 1
        x ^= x >> 11
        x ^= (x << 7) & 0x9d2c5680
        x ^= (x << 15) & 0xefc60000
        x ^= x >> 18
        return x


def uniform(engine):
    first = float(engine() - engine.minimum)
    second = float(engine() - engine.minimum)
    base = float(engine.span)
    return (first + second * base) / (base * base)


def booleans(engine=None):
    engine = engine or Minstd()
    while True:
        yield (engine() - engine.minimum) & 1


def bernoulli(engine=None):
    engine = engine or MT19937()
    while True:
        yield uniform(engine) < .5


def lognormals(sigma=2., engine=None, mean=-2.):
    engine = engine or Minstd()
    while True:
        while True:
            u = 2 * uniform(engine) - 1
            v = 2 * uniform(engine) - 1
            s = u*u + v*v
            if 0 < s <= 1:
                break
        factor = math.sqrt(-2 * math.log(s) / s)
        yield math.exp(u * factor * sigma + mean)
        yield math.exp(v * factor * sigma + mean)


def uniform_integers(maximum, engine=None):
    engine = engine or Minstd()
    span = maximum+1
    width = (span-1).bit_length()
    mask = (1 << width)-1
    limit = (engine.span >> width) << width
    while True:
        value = engine()-engine.minimum
        if value < limit and (value & mask) < span:
            yield value & mask
