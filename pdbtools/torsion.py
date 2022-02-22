#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.

__description__ = \
"""
Determines the dihedral angles (phi,psi) for each residue in a protein.
"""

__author__ = "Michael J. Harms"
__date__ = "070601"

import os, sys
from .helper import cmdline, geometry

def pdbTorsion(pdb):
    """
    Calculate the backbone torsion angles for a pdb file.
    """

    residue_list = []
    N = []
    CO = []
    CA = []

    resid_contents = {}
    current_residue = None
    to_take = ["N  ","CA ","C  "]
    for line in pdb:
        if line[0:4] == "ATOM" or (line[0:6] == "HETATM" and line[17:20] == "MSE"):

            if line[13:16] in to_take:

                # First residue
                if current_residue == None:
                    current_residue = line[17:26]

                # If we're switching to a new residue, record the previously
                # recorded one.
                if current_residue != line[17:26]:

                    try:
                        N.append([float(resid_contents["N  "][30+8*i:38+8*i])
                                  for i in range(3)])
                        CO.append([float(resid_contents["C  "][30+8*i:38+8*i])
                                   for i in range(3)])
                        CA.append([float(resid_contents["CA "][30+8*i:38+8*i])
                                   for i in range(3)])
                        residue_list.append(current_residue)

                    except KeyError:
                        err = "Residue %s has missing atoms: skipping.\n" % current_residue
                        sys.stderr.write(err)

                    # Reset resid contents dictionary
                    current_residue = line[17:26]
                    resid_contents = {}

                # Now record N, C, and CA entries.  Take only a unique one from
                # each residue to deal with multiple conformations etc.
                if line[13:16] not in resid_contents:
                    resid_contents[line[13:16]] = line
                else:
                    err = "Warning: %s has repeated atoms!\n" % current_residue
                    sys.stderr.write(err)

    # Record the last residue
    try:
        N.append([float(resid_contents["N  "][30+8*i:38+8*i])
                  for i in range(3)])
        CO.append([float(resid_contents["C  "][30+8*i:38+8*i])
                   for i in range(3)])
        CA.append([float(resid_contents["CA "][30+8*i:38+8*i])
                   for i in range(3)])
        residue_list.append(current_residue)

    except KeyError:
        err = "Residue %s has missing atoms: skipping.\n" % current_residue
        sys.stderr.write(err)


    # Calculate phi, psi and omega for each residue.  If the calculation fails, write
    # that to standard error and move on.
    labels = []
    dihedrals = []

    # the first residue, which does not have phi
    try:
        dihedrals.append(geometry.calcDihedrals(None,N[0],CA[0],CO[0],
                                                    N[1],CA[1]))
        labels.append(residue_list[0])
    except ValueError:
        err = "Dihedral calculation failed for %s\n" % residue_list[i]
        sys.stderr.write(err)

    # all intermediate residues                                            
    for i in range(1,len(residue_list)-1):
        try:
            dihedrals.append(geometry.calcDihedrals(CO[i-1],N[i],CA[i],CO[i],
                                                    N[i+1],CA[i+1]))
            labels.append(residue_list[i])
        except ValueError:
            err = "Dihedral calculation failed for %s\n" % residue_list[i]
            sys.stderr.write(err)

    # the last residue, which does not have psi and omega
    try:
        dihedrals.append(geometry.calcDihedrals(CO[-2],N[-1],CA[-1],CO[-1],
                                                    None,None))
        labels.append(residue_list[-1])
    except ValueError:
        err = "Dihedral calculation failed for %s\n" % residue_list[i]
        sys.stderr.write(err)
    return dihedrals, labels
