"""
Implements a wrapper to help count the number fo evaluations per run
"""


class CountingEvaluator:
    def __init__(self, udp):
        self.udp = udp
        self.count = 0

    def fitness(self, x):
        self.count += 1
        return self.udp.fitness(x)
