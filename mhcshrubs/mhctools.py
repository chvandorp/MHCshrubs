"""@package mhctools
A class MhcObject to represent the MHC nomenclature.

TODO: use private members (with underscores) and use property
for setters and getters.
"""

from __future__ import print_function
from builtins import (str, object)
from warnings import warn
import re ## for parsing MHC names

## auxiliary functions

def head(xs): return xs[0]
def tail(xs): return xs[1:]

## TODO: let MhcObject deal with these: possibly give parsing hints at init
def fromNetMhcFormat(hlastr): ## add a star before the numbers...
    exceptions = ["H-2-Db"] ## FIXME: handle this type
    if hlastr in exceptions: return hlastr
    digits = re.findall("\d+", hlastr) ## FIXME: class II uses _ as separator between locus and fields
    preamble = re.sub("[\d+:*]", '', hlastr)
    return preamble + '*' + ':'.join(digits)

def fromTreeFormat(hlastr):
    warn("use of fromTreeFormat is deprecated, use MhcObject(fmt='Newick') instead")
    hlastr = re.sub("[i]", ":", hlastr)
    hlastr = re.sub("[x]", "*", hlastr)
    return hlastr

class MhcObject(object):
    """
    A MHC object that handles MHC parsing,
    and can give various representations of the allele.
    See also http://hla.alleles.org/nomenclature/naming.html for naming conventions
    The MhcObject is also intended for translations between different notations
    """
    _species_sep = '-'
    _locus_sep = '*'
    _field_sep = ':'
    _locus_sep_newick = "x" ## looks a bit like *
    _field_sep_newick = 'i' ## looks a bit like :
    _unknown_field = '??'
    _unknown_locus = '?'
    _max_fields = 4
    _legal_fmts = ["default", "NetMHC", "Newick"]
    def __init__(self, init_name=None, fmt="default", species="HLA",
                 locus=None, field1=None, field2=None, field3=None,
                 field4=None, suffix=None):
        if fmt not in self._legal_fmts:
            raise ValueError("fmt must be one of " + ", ".join(self._legal_fmts))
        self._fields = [None for _ in range(self._max_fields)]
        if init_name is not None:
            ## remove any trailing whitespace
            self._init_name = init_name.strip()
            ## try to parse the init name
            name = self._init_name ## get a copy to 'tear down'
            ## chop off the species part
            xs = name.split(self._species_sep)
            if len(xs) > 1:
                ## allow species_sep (-) is the species identifier, to include e.g. "H-2"
                self._species = self._species_sep.join(xs[:-1])
                name = xs[-1]
            else:
                self._species = None
                name = xs[0]
            ## chop off the locus part
            if fmt == "NetMHC":
                raise Exception("parsing NetMHC format not implemented in MhcObject")
            elif fmt == "Newick":
                locus_sep = self._locus_sep_newick
            else: ## default
                locus_sep = self._locus_sep
            xs = name.split(locus_sep)
            if len(xs) > 1:
                self._locus = xs[0]
                name = self._locus_sep.join(xs[1:])
            else:
                self._locus = None
                name = xs[0]
            ## split the digits
            if fmt == "Newick":
                field_sep = self._field_sep_newick
            else: ## default or NetMHC
                field_sep = self._field_sep
            xs = name.split(field_sep)
            digits = re.compile("[0-9]+")
            for i, x in enumerate(xs[:self._max_fields]): ## if len(xs) < 4, then xs[:4] is just the shorter part
                dx = head(digits.findall(x))
                if i+1 != len(xs) and dx != x:
                    warn("ignoring invalid characters in field {0} of '{1}'".format(i+1, self._init_name))
                self._fields[i] = dx
            if len(xs) > self._max_fields:
                wstr = "the number of MHC fields exceeds {0}. Ignoring fields {1}."
                warn(wstr.format(self._max_fields, field_sep.join(xs[self._max_fields:])))
            ## find the suffix
            letters = re.compile("[a-zA-Z]+") ## TODO restrict to known suffixes
            suf = letters.findall(xs[-1])
            if len(suf) == 1:
                self._suffix = suf[0]
            else:
                self._suffix = None
                if len(suf) > 1:
                    warn("suffix is not well-formed. Ignoring suffix {0}.".format(xs[-1]))
            ## override some fields with explicitly specified values
            if species is not None: self._species = species
            if locus is not None: self._locus = locus
            if field1 is not None: self._fields[0] = field1
            if field2 is not None: self._fields[1] = field2
            if field3 is not None: self._fields[2] = field3
            if field4 is not None: self._fields[3] = field4
            if suffix is not None: self._suffix = suffix
        else:
            self._species = species
            self._locus = locus
            self._fields[0] = field1
            self._fields[1] = field2
            self._fields[2] = field3
            self._fields[3] = field4
            self._suffix = suffix
    @property
    def species(self):
        return self._species
    @property
    def locus(self):
        return self._locus
    @property
    def field1(self):
        return self._fields[0]
    @property
    def field2(self):
        return self._fields[1]
    @property
    def field3(self):
        return self._fields[2]
    @property
    def field4(self):
        return self._fields[3]
    @property
    def suffix(self):
        return self._suffix
    def __len__(self):
        """
        The len function returns the resolution
        i.e. the number of fields.
        """
        for i, field in enumerate(self._fields[::-1]):
            if field is not None:
                break
        return self._max_fields - i

    def __repr__(self):
        return "<MhcObject {}>".format(self.__str__())
    def __str__(self):
        return self.full_str()
    def __format__(self, fmt):
        """
        Used for string formatting.

        Todo:
            - what is the correct behaviour for the __format__ method?
            - add a 1f, 2f etc fmt character

        Example:
            allele = MhcObject("HLA-A*01:01:01:01")
            msg = "The most common HLA A subtype is {:st}".format(allele)
            print(msg) ## The most common HLA A subtype is HLA-A*01:01
        """
        if fmt == "tp":
            return self.type_str()
        elif fmt == "st":
            return self.subtype_str()
        elif fmt == "nm":
            return self.NetMHC_str()
        elif fmt == "nw":
            return self.Newick_str()
        else:
            return self.full_str()
    def __key__(self):
        """
        used for hashing and comparing.
        We need e.g. 01:100 to be greater than 01:33
        """
        return ("" if self.species is None else self.species,
                "" if self.locus is None else self.locus,
                0 if self.field1 is None else int(self.field1),
                0 if self.field2 is None else int(self.field2),
                0 if self.field3 is None else int(self.field3),
                0 if self.field4 is None else int(self.field4),
                "" if self.suffix is None else self.suffix)
    def __hash__(self):
        return hash(self.__key__())
    def __eq__(self, other):
        return (self.species == other.species and
                self.locus == other.locus and
                self.field1 == other.field1 and
                self.field2 == other.field2 and
                self.field3 == other.field3 and
                self.field4 == other.field4 and
                self.suffix == other.suffix)
    def __ne__(self, other):
        return not (self == other)
    def __lt__(self, other):
        return self.__key__() < other.__key__()
    def __contains__(self, other):
        """
        The members of other must be equal to selfs members or None.
        """
        if not isinstance(other, type(self)):
            raise TypeError
        if self.species is not None and other.species != self.species:
            return False
        elif self.locus is not None and other.locus != self.locus:
            return False
        elif self.field1 is not None and other.field1 != self.field1:
            return False
        elif self.field2 is not None and other.field2 != self.field2:
            return False
        elif self.field3 is not None and other.field3 != self.field3:
            return False
        elif self.field4 is not None and other.field4 != self.field4:
            return False
        elif self.suffix is not None and other.suffix != self.suffix:
            return False
        else: return True
    @property
    def type(self):
        return MhcObject(species=self.species, locus=self.locus, field1=self.field1)
    @property
    def subtype(self):
        return MhcObject(species=self.species, locus=self.locus,
                         field1=self.field1, field2=self.field2)
    @property
    def protein(self):
        """
        protein is identical to subtype, but the expression suffix is added.
        """
        return MhcObject(species=self.species, locus=self.locus,
                         field1=self.field1, field2=self.field2, suffix=self.suffix)
    def make_str(self, print_species=True, fmt="default"):
        """
        Format the MhcObject as a string.

        Kwargs:
            print_species (bool): if False, do not print the species attribute
                for brevity
            fmt (str): choose between
                - "default" -- use the default nomenclature
                - "NetMHC" -- do not print the locus separator
                - "Newick" -- a format that can be used in Newick trees and filenames

        Returns:
            a string representation of the MHC object
        """
        if fmt not in self._legal_fmts:
            raise ValueError("fmt must be one of " + ", ".join(self._legal_fmts))
        ## ignore most significant None fields
        field_strs = []
        for field in self._fields[::-1]:
            if field is None:
                if len(field_strs) > 0:
                    field_strs.append(self._unknown_field)
            else:
                field_strs.append(field)
        field_strs.reverse()
        field_sep = self._field_sep_newick if fmt == "Newick" else self._field_sep
        mhc_str = field_sep.join(field_strs)
        ## add remaining attributes
        if self._suffix is not None:
            mhc_str += str(self._suffix)
        locus_str = self._unknown_locus if self._locus is None else str(self._locus)
        if fmt == "default":
            locus_sep = self._locus_sep
        elif fmt == "Newick":
            locus_sep = self._locus_sep_newick
        else: ## fmt == "NetMHC"
            locus_sep = ""
        mhc_str = locus_str + locus_sep + mhc_str
        if self._species is not None and print_species:
            mhc_str = str(self._species) + self._species_sep + mhc_str
        return mhc_str
    def full_str(self):
        return self.make_str(print_species=True)
    def short_str(self):
        """Do not print the species attribute (e.g. 'HLA')"""
        return self.make_str(print_species=False)
    def type_str(self):
        """Print only the first field"""
        return self.type.full_str()
    def subtype_str(self):
        """Print only the first two fields"""
        return self.subtype.full_str()
    def protein_str(self):
        """Print only the first two fields and the suffix"""
        return self.protein.full_str()
    def NetMHC_str(self):
        """NetMHC doesn't use the locus separator (*)"""
        return self.make_str(fmt="NetMHC")
    def Newick_str(self):
        """
        Putting MHC allele names in a Newick tree can be problematic
        since the allele separator (:) is used in the newick format.
        Also, the locus separator (*) is not used.
        Can also be useful for filenames.
        """
        return self.make_str(fmt="Newick")
