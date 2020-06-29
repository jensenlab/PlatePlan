"""Ingredients for building PlatePlan Experiments and Assays."""
import datetime

from . import makeids
from . import units
from . import utils


class Ingredient(object):
    def __init__(self, id_):
        self.id = id_

    def __hash__(self):
        return hash(self.id)

    def __eq__(self, other):
        return hash(self) == hash(other)

    def __repr__(self):
        return str(self)


class Stock(Ingredient):
    def __init__(
        self,
        ingredient,
        id_=None,
        quantity=None,
        labware=None,
        location=None,
        date_made=None,
        date_expires=None,
    ):
        if id_ is None:
            id_ = makeids.unique_id(prefix=ingredient.id)
        super().__init__(id_)

        self.ingredient = ingredient
        if date_made is None:
            date_made = datetime.date.today()
        if quantity is not None:
            self.quantity = units.parse(quantity)
        else:
            self.quantity = units.parse("0 l")
        self.labware = labware
        self.location = location
        self.date_made = date_made
        self.date_expires = date_expires

    def __str__(self):
        return "Stock({}): {} of {}".format(
            self.id, str(self.quantity), repr(self.ingredient)
        )

    def __hash__(self):
        # stock_id is a hex hash key but hash() must return an integer
        return int(self.id, 16)

    @property
    def ingredient_ids(self):
        return self.ingredient.ingredient_ids


class Reagent(Ingredient):
    def __init__(self, id_, name=None, molecular_weight=None, **kwargs):
        super().__init__(id_)
        self.name = name
        self.molecular_weight = molecular_weight

    def duplicate(self, other):
        """Copy attributes when making ReagentStock objects."""
        self.id = other.id
        self.name = other.name
        self.molecular_weight = other.molecular_weight

    def __str__(self):
        return "Reagent(" + self.id + ")"

    @property
    def ingredient_ids(self):
        return [self.id]


class Solution(Ingredient):
    def __init__(self, reagents, solvent=None, id_=None):
        # reagents is a dict of {reagent_id: concentration}
        # eventually we should check that the reagent ID's are valid
        # for now let solvent be a string
        if id_ is None:
            id_ = makeids.unique_id()
        super().__init__(id_)

        self.reagents = {r: units.parse(u) for r, u in reagents.items()}
        self.solvent = solvent

    def __str__(self):
        s = "Solution {} in {} with:".format(self.id, self.solvent)
        for reagent in self.reagents:
            s += "\n   {} {}".format(str(self.reagents[reagent]), reagent)
        return s

    def __repr__(self):
        return "Solution({})".format(self.id)

    def concentrated(self, factor):
        reagents = {r: factor * c for r, c in self.reagents.items()}
        id_ = self.id + " x" + str(factor)
        return Solution(reagents, solvent=self.solvent, id_=id_)

    def without(self, reagent_ids):
        # support when called as solution.without('ingredient') by converting
        # to solution.without(['ingredient'])
        reagent_ids = utils.assert_list(reagent_ids)
        reagents = {r: c for r, c in self.reagents.items() if r not in reagent_ids}
        if reagent_ids:
            id_ = self.id + " -- " + " -- ".join(reagent_ids)
        else:
            id_ = self.id
        return Solution(reagents, solvent=self.solvent, id_=id_)

    def with_only(self, reagent_ids):
        to_drop = [r for r in self.reagents if r not in reagent_ids]
        return self.without(to_drop)

    def updated(self, reagents):
        new = self.without([])
        new.reagents.update({r: units.parse(reagents[r]) for r in reagents})
        return new

    def remove_reagents(self, reagent_ids):
        new = self.without([])
        reagent_ids = utils.assert_list(reagent_ids)
        for r in reagent_ids:
            del new.reagents[r]
        return new

    @property
    def ingredient_ids(self):
        return list(self.reagents.keys())


if __name__ == "__main__":
    r1 = Reagent("reagent1")
    r2 = Reagent("reagent2")
    sol = Solution(dict(reagent1="10 g/l", reagent2="23 ug/ul"), id_="r12")
    print(r1)
    print(r2)
    print(sol)
    print(Stock(sol))
