"""Unit parsing and arithmetic.

Almost all quantities for PlatePlan require units. Rather than store
units separately and manage conversions case-by-case, we use a set of *Unit
classes (MaterialUnit, VolumeUnit, ConcentrationUnit, etc.). These classes
define arithmetic operations to abstract unit conversion details.

The only way to create *Unit objects is the `parse` function. After *Unit
object are created, they can be used in calculations. The following operations
are defined for *Unit objects:
    U = Unit (either MaterialUnit, VolumeUnit, ConcentrationUnit, etc.)
    M = MaterialUnit
    V = VolumeUnit
    C = ConcentrationUnit
    MM = MolarMassUnit
    <num> is a number

   <num> * U
   U * <num>
   U / <num>
   U{T} / U{T} (= <num>)
   U{T} + U{T} (= U{T})
   M / V (= C)
   C * V (= M)
   V * C (= M)
   M / C (= V)
   MM * M (= M)
   
   
Units should be treated as immutable objects (just like numbers). All 
arithmetic operations return new *Unit objects. You should not access
the attributes, and you should have no need to. Doing so is a sign that 
you're not treating the units of the quantity properly.
"""

import copy
import math

import pyparsing as pp

PREFIX_FACTORS = {
    "k": 1.0e3,
    None: 1.0,
    "c": 1.0e-2,
    "m": 1.0e-3,
    "u": 1.0e-6,
    "n": 1.0e-9,
    "p": 1.0e-12,
}


def find_prefix(value, prefixes=PREFIX_FACTORS):
    """Find the closest prefix for a number.
    
    Selects the largest prefix such that the scaled number is > 1.0
    Example: find_prefix(0.1) -> 'c' since 0.1 m = 10 cm
    """
    prefixes = sorted(prefixes.items(), key=lambda x: -x[1])
    for prefix, factor in prefixes:
        if value >= factor:
            break
    return prefix


def _make_unit_str(prefix, unit, volprefix=None, volunit=None):
    """Format a unit string using a (possibly None) prefix."""
    if volunit is None:
        return "{}{}".format(
            "" if prefix is None else prefix, "" if unit is None else unit
        )
    else:
        return "{}{}/{}{}".format(
            "" if prefix is None else prefix,
            "" if unit is None else unit,
            "" if volprefix is None else volprefix,
            volunit,
        )


class ConversionError(Exception):
    def __init__(self, from_unit, to_unit, message="Error converting."):
        self.from_unit = from_unit
        self.to_unit = to_unit
        self.message = message


class Unit(object):
    def __init__(self, value, unit, prefix=None):
        self.value = value
        self.unit = unit
        self.prefix = prefix

    def __mul__(self, other):
        other = parse(other)
        if isinstance(other, Unit):
            return NotImplemented
        new = copy.copy(self)
        new.value = other * self.value
        return new

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        other = parse(other)
        if isinstance(other, Unit):
            # U{T} / U{T} (= <num>)
            try:
                other = other.convert(self)
            except ConversionError:
                return NotImplemented
            return self.value / other.value
        else:
            # U / <num>
            return self * (1 / other)

    def __add__(self, other):
        try:
            new = other.convert(self)
        except ConversionError:
            return NotImplemented
        new.value += self.value
        return new.auto_prefix()

    def __repr__(self):
        return "Unit({} {})".format(self.value, self.unit_str)

    def convert(self, other):
        other = parse(other)
        if self.unit != other.unit:
            raise ConversionError(self, other)
        factor = math.log10(PREFIX_FACTORS[self.prefix]) - math.log10(
            PREFIX_FACTORS[other.prefix]
        )
        new = copy.copy(self)
        new.value *= 10 ** factor
        new.prefix = other.prefix
        new.unit = other.unit
        return new

    @property
    def base_value(self):
        return self.convert(self.unit).value

    def auto_prefix(self):
        prefix = find_prefix(self.base_value)
        return self.convert(prefix + self.unit)

    @property
    def unit_str(self):
        return _make_unit_str(self.prefix, self.unit)

    def __str__(self):
        if self.value is None:
            return self.unit_str
        else:
            return "{} {}".format(self.value, self.unit_str)


class MaterialUnit(Unit):
    def __truediv__(self, other):
        other = parse(other)
        if isinstance(other, VolumeUnit):
            # M / V (= C)
            return ConcentrationUnit(
                self.value / other.value,
                self.unit,
                self.prefix,
                den_unit=other.unit,
                den_prefix=other.prefix,
            )
        elif isinstance(other, ConcentrationUnit):
            # M / C (= V)
            new = other.volume_unit
            new.value = self.value / other.material_unit.convert(self).value
            return new
        else:
            return NotImplemented

    def __mul__(self, other):
        if isinstance(other, MolarMassUnit):
            # M * MM (= M)
            new = other.mass_unit
            new.value *= self.convert(other.material_unit).value
            return new

        other = parse(other)
        new = copy.copy(self)
        new.value = other * self.value
        return new

    def encode_json(self):
        return {
            "type": "MaterialUnit",
            "value": self.value,
            "unit": self.unit,
            "prefix": self.prefix,
        }

    @classmethod
    def decode_json(cls, d):
        return cls(d["value"], d["unit"], d["prefix"])


class VolumeUnit(Unit):
    def __mul__(self, other):
        other = parse(other)
        if isinstance(other, ConcentrationUnit):
            # V * C (= M)
            new = other.material_unit
            new.value *= self.convert(other.volume_unit).value
            return new
        else:
            # U * <num>
            return super().__mul__(other)

    def __rmul__(self, other):
        return self * other

    def encode_json(self):
        return {
            "type": "VolumeUnit",
            "value": self.value,
            "unit": self.unit,
            "prefix": self.prefix,
        }

    @classmethod
    def decode_json(cls, d):
        return cls(d["value"], d["unit"], d["prefix"])


class ConcentrationUnit(Unit):
    def __init__(self, value, unit, prefix, den_unit, den_prefix=None):
        self.value = value
        self.unit = unit
        self.prefix = prefix
        self.den_unit = den_unit
        self.den_prefix = den_prefix

    @property
    def base_value(self):
        den_unit = Unit(self.value, self.unit, self.prefix)
        return den_unit.convert(self.unit).value

    def auto_prefix(self):
        prefix = find_prefix(self.base_value)
        return self.convert(prefix + self.unit + "/" + self.den_prefix + self.den_unit)

    @property
    def material_unit(self):
        return MaterialUnit(self.value, self.unit, self.prefix)

    @property
    def volume_unit(self):
        return VolumeUnit(1, self.den_unit, self.den_prefix)

    @property
    def unit_str(self):
        return _make_unit_str(self.prefix, self.unit, self.den_prefix, self.den_unit)

    def convert(self, other, molar_mass=None):
        other = parse(other)

        if not isinstance(other, ConcentrationUnit):
            ConversionError(self, other)
        factor = (
            math.log10(PREFIX_FACTORS[self.prefix])
            - math.log10(PREFIX_FACTORS[other.prefix])
            - math.log10(PREFIX_FACTORS[self.den_prefix])
            + math.log10(PREFIX_FACTORS[other.den_prefix])
        )
        new = copy.copy(self)
        new.value *= 10 ** factor
        new.prefix = other.prefix
        new.den_prefix = other.den_prefix
        new.unit = other.unit
        new.den_unit = other.den_unit

        if self.unit == "g" and other.unit == "mol":
            if not molar_mass:
                raise ConversionError(self, other, "Missing molar mass.")
            new.value /= molar_mass
        elif self.unit == "mol" and other.unit == "g":
            if not molar_mass:
                raise ConversionError(self, other, "Missing molar mass.")
            new.value *= molar_mass
        elif self.unit != other.unit or self.den_unit != other.den_unit:
            raise ConversionError(self, other, "Units don't match.")
        return new

    def __str__(self):
        return "{} {}/{}".format(
            str(self.value),
            _make_unit_str(self.prefix, self.unit),
            _make_unit_str(self.den_prefix, self.den_unit),
        )

    def encode_json(self):
        return {
            "type": "ConcentrationUnit",
            "value": self.value,
            "unit": self.unit,
            "prefix": self.prefix,
            "den_unit": self.den_unit,
            "den_prefix": self.den_prefix,
        }

    @classmethod
    def decode_json(cls, d):
        return cls(d["value"], d["unit"], d["prefix"], d["den_unit"], d["den_prefix"])


class TemperatureUnit(Unit):
    def encode_json(self):
        return {
            "type": "TemperatureUnit",
            "value": self.value,
            "unit": self.unit,
            "prefix": self.prefix,
        }

    @classmethod
    def decode_json(cls, d):
        return cls(d["value"], d["unit"], d["prefix"])


class MolarMassUnit(Unit):
    def __init__(self, value, unit, prefix, den_unit, den_prefix=None):
        self.value = value
        self.unit = unit
        self.prefix = prefix
        self.den_unit = den_unit
        self.den_prefix = den_prefix

    @property
    def base_value(self):
        den_unit = Unit(self.value, self.unit, self.prefix)
        return den_unit.convert(self.unit).value

    def auto_prefix(self):
        prefix = find_prefix(self.base_value)
        return self.convert(prefix + self.unit + "/" + self.den_prefix + self.den_unit)

    @property
    def mass_unit(self):
        return MaterialUnit(self.value, self.unit, self.prefix)

    @property
    def material_unit(self):
        return MaterialUnit(1, self.den_unit, self.den_prefix)

    @property
    def unit_str(self):
        return _make_unit_str(self.prefix, self.unit, self.den_prefix, self.den_unit)

    def convert(self, other):
        other = parse(other)
        if not isinstance(other, MolarMassUnit):
            raise ConversionError(self, other)
        if self.unit != other.unit or self.den_unit != other.den_unit:
            raise ConversionError(self, other)
        factor = (
            math.log10(PREFIX_FACTORS[self.prefix])
            - math.log10(PREFIX_FACTORS[other.prefix])
            - math.log10(PREFIX_FACTORS[self.den_prefix])
            + math.log10(PREFIX_FACTORS[other.den_prefix])
        )
        new = copy.copy(self)
        new.value *= 10 ** factor
        new.prefix = other.prefix
        new.den_prefix = other.den_prefix
        new.unit = other.unit
        new.den_unit = other.den_unit
        return new

    def __str__(self):
        return "{} {}/{}".format(
            str(self.value),
            _make_unit_str(self.prefix, self.unit),
            _make_unit_str(self.den_prefix, self.den_unit),
        )

    def encode_json(self):
        return {
            "type": "MolarMassUnit",
            "value": self.value,
            "unit": self.unit,
            "prefix": self.prefix,
            "den_unit": self.den_unit,
            "den_prefix": self.den_prefix,
        }

    @classmethod
    def decode_json(cls, d):
        return cls(d["value"], d["unit"], d["prefix"], d["den_unit"], d["den_prefix"])


# =============== Unit Parsing ===============

# these need to be entered separately into the parser definition below;
# we can't automate it since 'l' is caseless but 'g' and 'mol' are not
MATERIAL_UNITS = ["g", "mol"]
VOLUME_UNITS = ["l"]
CONCENTRATION_UNITS = ["M"]
TEMPERATURE_UNITS = ["C"]


def _is_material_unit(unit):
    return unit in MATERIAL_UNITS


def _is_volume_unit(unit):
    return unit.lower() in VOLUME_UNITS


def _is_concentration_unit(unit):
    return unit in CONCENTRATION_UNITS


def _is_temperature_unit(unit):
    return unit.upper() in TEMPERATURE_UNITS


_plusminus = pp.oneOf("+ -")
_digitpart = pp.Word(pp.nums)
_signed_integer = pp.Optional(_plusminus) + _digitpart
_point = pp.Literal(".")
_floating = (_signed_integer + _point + pp.Optional(_digitpart)) | (
    pp.Optional(_plusminus) + _point + _digitpart
)
_exponent = pp.CaselessLiteral("e") + _signed_integer
_number = pp.Combine((_signed_integer ^ _floating) + pp.Optional(_exponent))
_base_unit = pp.oneOf("mol g M C") ^ pp.CaselessLiteral("l")
_prefixes = [k for k in PREFIX_FACTORS.keys() if k is not None]
_prefix = ~_base_unit + pp.oneOf(" ".join(_prefixes))
_unit = pp.Group(pp.Optional(_prefix, default=None) + _base_unit)
_units = pp.Group(_unit + pp.Optional(pp.Literal("/").suppress() + _unit, default=None))

_unit_expr = (
    pp.Optional(_number, default=None)
    + pp.Optional(_units, default=None)
    + pp.StringEnd()
)


class UnknownUnitError(Exception):
    def __init__(self, unit):
        self.unit = unit


def parse(s, unit_str=None):
    """Add units to a number or parse a string with units.
    
    Parses numbers with units and returns a *Unit object. Can be called
    in six ways:
        1. With a number and a unit string:
            parse(10, "mg/l") or parse("10", "mg/l")
        2. With a single string containing both a number and a unit:
            parse("10 mg/l")
        There is no need for a space between the number and unit.
        3. With only a unit (no value):
            parse("mg") # returns Unit(None mg)
        4. With a *Unit object:
            unit = parse(10, "mg")
            parse(unit)  # returns the unit unaltered
        This form is convenient if you're unsure if a value has already
        been parsed. Parsing it again won't hurt, so go ahead.
        5. With a number, in which case a float is returned:
            parse(10) # returns 10.0
        6. With None, in which case None is returned:
            parse(None) # returns None
    """

    if s is None:
        # Case 6
        if unit_str is not None:
            return parse(unit_str)
        else:
            return None
    if not isinstance(s, str) and unit_str is None:
        # Case 5
        # called as parse(<number>); return just the number
        return s
    if isinstance(s, Unit):
        # Case 4
        # s is already a parsed unit; just return it
        return s
    if not isinstance(s, str) and unit_str is not None:
        # Case 1
        # called as parse(<number>, "units"); combine to single str for parsing
        s = str(s) + unit_str
    elif isinstance(s, str) and unit_str is not None:
        s = s + unit_str

    # Case 2 or 3; just parse s

    try:
        num, units = _unit_expr.parseString(s.strip())
        # print(s, "--", num, units)
    except pp.ParseException:
        raise UnknownUnitError(s.strip())
    if num is not None:
        num = float(num)
    if units is None:
        return num

    if units[1] is None:
        if _is_material_unit(units[0][1]):
            return MaterialUnit(num, unit=units[0][1], prefix=units[0][0])
        elif _is_volume_unit(units[0][1]):
            return VolumeUnit(num, unit=units[0][1], prefix=units[0][0])
        elif _is_concentration_unit(units[0][1]):
            # molar concentration are stored as xmol/l
            return ConcentrationUnit(
                num, unit="mol", prefix=units[0][0], den_unit="l", den_prefix=None
            )
        elif _is_temperature_unit(units[0][1]):
            return TemperatureUnit(num, unit=units[0][1], prefix=units[0][0])
        else:
            raise UnknownUnitError(units[0][1])
    else:
        if units[0][1] == "g" and units[1][1] == "mol":
            return MolarMassUnit(
                num,
                unit=units[0][1],
                prefix=units[0][0],
                den_unit=units[1][1],
                den_prefix=units[1][0],
            )
        else:
            return ConcentrationUnit(
                num,
                unit=units[0][1],
                prefix=units[0][0],
                den_unit=units[1][1],
                den_prefix=units[1][0],
            )


if __name__ == "__main__":
    print(parse(10, "mg/l"))
    print(parse("10", "mg/l"))
    print(parse("10 mg/l"))
    print(parse("mg"))
    print(parse(parse(10, "mg")))
    print(parse(10))
    print(parse(None))
    print(parse("37 C"))

    m = parse("1 mg")
    v = parse("2 l")
    c = parse("3 mg/l")

    print(4 * m)
    print(4 * v)
    print(4 * c)

    print(m * 4)
    print(v * 4)
    print(c * 4)

    # print(m / 4)
    print(v / 4)
    print(c / 4)

    # print(m / m)
    print(v / v)
    print(c / c)

    print(m + m)
    # print(v + v)
    # print(c + c)

    print(m / v)
    print(c * v)
    print(v * c)
    print(m / c)

    print(parse("1.314 mg/l").convert("ug/ul").auto_prefix())
    print(parse("1.314 mM") * 10)
    print(parse("5 mol/l") * parse("1 l"))
    print(parse("1 mg").convert("ug").auto_prefix())
    print(parse("1 mg") + parse("538 ug"))
    # print(parse("1 L") / parse("50 ml"))
    print(parse("1 mg") / parse("50 ml"))
    print(parse("1 mg") / parse("50 mg/ml"))
    print(parse("1 ml") * parse("50 mg/ul"))
    print(parse("50 mg/ul") * parse("1 ml"))

    print(parse("5 g/mol") * parse("1000 mmol"))
    print(parse("5 g/l").convert("M", molar_mass=5))
    print(parse("5 g/l").convert("mol/l", molar_mass=5))
    print(parse("5 g/l").convert("mmol/l", molar_mass=5))
    print(parse("5 mol/l").convert("g/l", molar_mass=5))
    print(parse("5 mol/l").convert("mg/ul", molar_mass=5))
    print(parse("5 M").convert("mg/mol", molar_mass=5))
    print(parse("5 mol/l").convert("mol/l"))

