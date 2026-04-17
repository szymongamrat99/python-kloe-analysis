import ROOT
import json


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

    Or load from JSON:
        cols = ColumnDefiner(rdf).from_json("config/columns.json")
        rdf = cols.apply()
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

    def from_json(self, path):
        """Load column definitions from a JSON file.

        JSON format:
            {
                "combined_Pi0": "sqrt(pow(pi01Fit[5] - PhysicsConstants::mPi0, 2) + ...)",
                "totalE": "pi01Fit[3] + pi02Fit[3]"
            }
        """
        with open(path, 'r') as f:
            data = json.load(f)
        self.add_many(data)
        return self

    @staticmethod
    def _load_json(path):
        """Load column definitions dict from JSON (without needing an RDF)."""
        with open(path, 'r') as f:
            return json.load(f)

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
