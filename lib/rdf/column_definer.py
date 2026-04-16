import ROOT


class ColumnDefiner:
    """Batch-define new columns on an RDataFrame node.

    Usage:
        cols = ColumnDefiner(rdf)
        cols.add("combined_Pi0",
                 "sqrt(pow(pi01Fit[5] - PhysicsConstants::mPi0, 2)"
                 " + pow(pi02Fit[5] - PhysicsConstants::mPi0, 2))")
        cols.add("totalE", "pi01Fit[3] + pi02Fit[3]")
        rdf = cols.apply()          # returns the new node with all columns
        print(cols.names())         # ["combined_Pi0", "totalE"]
    """

    def __init__(self, rdf):
        self._rdf_base = rdf
        self._definitions = []

    def add(self, name, expression):
        """Register a new column (name, C++ expression string)."""
        self._definitions.append((name, expression))
        return self

    def add_many(self, defs):
        """Register multiple columns at once.

        defs - dict {name: expression} or list of (name, expression) tuples.
        """
        items = defs.items() if isinstance(defs, dict) else defs
        for name, expr in items:
            self._definitions.append((name, expr))
        return self

    def apply(self):
        """Apply all registered Define() calls and return the new RDF node."""
        node = self._rdf_base
        for name, expr in self._definitions:
            node = node.Define(name, expr)
        return node

    def names(self):
        """Return the list of defined column names."""
        return [name for name, _ in self._definitions]

    def __len__(self):
        return len(self._definitions)

    def __repr__(self):
        return f"ColumnDefiner({len(self._definitions)} columns: {self.names()})"
