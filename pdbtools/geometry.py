#!/usr/bin/env python

# Copyright 2022, Jie Li
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.

__author__ = "Jie Li"
__date__ = "220314"

import os, sys
from .helper import cmdline, geometry
from math import nan

def pdbGeom(pdb):
    """
    Calculate the backbone bond lengths and bond angles for the given pdb file
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


    # Calculate N-CA, CA-C, C-N bond lengths and N-CA-C, CA-C-N, C-N-CA bond angles
    #  for each residue.  If the calculation fails, write
    # that to standard error and move on.
    labels = []
    bond_lengths = []
    bond_angles = []

    # all residues other than the last one
    for i in range(0,len(residue_list)-1):
        try:
            n_ca_bl = geometry.dist(N[i],CA[i])
            ca_c_bl = geometry.dist(CA[i],CO[i])
            c_n_bl = geometry.dist(CO[i],N[i+1])
            n_ca_c_angle = geometry.calcAngle(N[i],CA[i],CO[i])
            ca_c_n_angle = geometry.calcAngle(CA[i],CO[i],N[i+1])
            c_n_ca_angle = geometry.calcAngle(CO[i],N[i+1],CA[i+1])
            bond_lengths.append([n_ca_bl, ca_c_bl, c_n_bl])
            bond_angles.append([n_ca_c_angle, ca_c_n_angle, c_n_ca_angle])
            labels.append(residue_list[0])
        except ValueError:
            err = "Dihedral calculation failed for %s\n" % residue_list[i]
            sys.stderr.write(err)


    # the last residue, which do not have C-N bond length and CA-C-N, C-N-CA bond angles
    try:
        n_ca_bl = geometry.dist(N[i],CA[i])
        ca_c_bl = geometry.dist(CA[i],CO[i])
        n_ca_c_angle = geometry.calcAngle(N[i],CA[i],CO[i])
        bond_lengths.append([n_ca_bl, ca_c_bl, nan])
        bond_angles.append([n_ca_c_angle, nan, nan])
        labels.append(residue_list[0])

    except ValueError:
        err = "Dihedral calculation failed for %s\n" % residue_list[i]
        sys.stderr.write(err)
    return bond_lengths, bond_angles, labels