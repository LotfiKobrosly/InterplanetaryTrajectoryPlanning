"""
Implements a wrapper to help count the number fo evaluations per run (by Claude)
"""

class CountingEvaluator:
    def __init__(self, udp):
        # Use object.__setattr__ to avoid triggering __setattr__ overrides
        object.__setattr__(self, '_udp', udp)
        object.__setattr__(self, 'count', 0)

    def fitness(self, x):
        object.__setattr__(self, 'count', self.count + 1)
        return self._udp.fitness(x)

    def __getattr__(self, name):
        # Only called when normal attribute lookup fails
        # At this point _udp is guaranteed to exist
        return getattr(object.__getattribute__(self, '_udp'), name)
