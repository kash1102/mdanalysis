.. -*- coding: utf-8 -*-

.. working with bonds angles and torsions

Working with bonds, angles and torsions
=======================================

Chemical bonds between atoms can be manipulated using MDAnalysis in an
object oriented fashion similar to AtomGroups.  MDAnalysis supports bonds,
angles, dihedrals (also known as torsions) and improper dihedrals.
These are made available as attributes of an AtomGroup
::

     >>> import MDAnalysis as mda
     >>> from MDAnalysisTests.datafiles import PSF, DCD
     >>> u = mda.Universe(PSF, DCD)
     >>> u.atoms.bonds
     <TopologyGroup containing 3365 bonds>
     >>> u.atoms.angles
     <TopologyGroup containing 6123 angles>
     >>> u.atoms.dihedrals
     <TopologyGroup containing 8921 dihedrals>
     >>> u.atoms.impropers
     <TopologyGroup containing 541 impropers>

If bonding information was not included in the topology file, it is also
possible to guess the bonds present in the system based upon the distances
between atoms.
::

     >>> import MDAnalysis as mda
     >>> u = mda.Universe('myfile.gro')
     >>> u.atoms.guess_bonds(vdwradii={'C': 1.61, 'N':1.64}

For large systems this will become very slow as it searches all pairwise
combinations.  To circumvent this you can instead only search within a
smaller group of atoms, for example looping through the segments in a
system.
::

     for segment in u.segments:
	 segment.atoms.guess_bonds()


TopologyGroups
--------------

Regardless of the type of bond that is being dealt with, the main working
object for handling them is the TopologyGroup.

.values method (uses Cython so should be fast)


Selecting bonds
---------------

For examples working with a box of ethanol::

    >>> import MDAnalysis as mda
    >>> u = mda.Universe('ethanol.gro', guess_bonds=True)
    >>> u.bonds
    <TopologyGroup containing 11784 Bonds>
    >>> u.bonds.types()  # view available types
    [('O', 'H'), ('C', 'O'), ('C', 'H'), ('C', 'C')]
    >>> u.bonds.select_bonds(('C', 'O'))  # return all C-O bonds from the group
    <TopologyGroup containing 1473 Bonds>

Bonds are categorised based on the types of atoms.  This is done in a way
so that type (a, b, c) is equivalent to (c, b, a) ie. bonds are reversible.
For example::

    >>> u.angles.types()
    [('C', 'C', 'H'),
     ('H', 'C', 'H'),
     ('C', 'O', 'H'),
     ('C', 'C', 'O'),
     ('H', 'C', 'O')]

There are only C-C-H bonds and no H-C-C bonds.  Selection however is
aware that sometimes types are reversed::

    >>> u.angles.select_bonds(('H', 'C', 'C'))  # note reversal of type
    <TopologyGroup containing 7365 Angles>

TopologyGroups can be combined and indexed::

    >>> tg = u.angles.select_bonds(('C', 'C', 'O')) + u.angles.select_bonds(('C', 'O', 'H'))
    >>> tg.types()
    [('C', 'O', 'H'), ('C', 'C', 'O')]
    >>> tg[:100]
    <TopologyGroup containing 100 Angles>


Intersections
-------------

TopologyGroup.atomgroup_intersection method explanation
