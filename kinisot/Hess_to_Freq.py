#!/usr/bin/python

# Comments and/or additions are welcome (send e-mail to:
# robert.paton@colostate.edu

import sys, math
import numpy as np

def read_hess(file, iso):
   # The force constant matrix is read from g09 ouptut
   # The matrix values are mass-weighted according to the isotopic masses
   # Vibrational scaling factors are not applied to matrix elements at this stage:
   # the resulting frequencies can be scaled after diagonalization

   # Read gaussian output
   g09_output = open(file, 'r')
   inlines = g09_output.readlines()

   mass_list = []
   start_force = len(inlines) < 100
   for i in range(0,len(inlines)):
      if inlines[i].strip().find('and mass') > -1:
         mass_list.append(float(inlines[i].strip().split()[8]))

      if inlines[i].strip().startswith('NAtoms='):
         d_o_f = int(inlines[i].strip().split()[1])*3

      if inlines[i].strip().find('NImag') > -1: start_force = i

   # Hessian and mass-weighted Hessian (3N x 3N matrices, N = no. atoms)
   hess_mat = np.ndarray(shape=(d_o_f,d_o_f))
   mw_hess_mat = np.ndarray(shape=(d_o_f,d_o_f))

   # Isotopic substitution will consider 1H/2H, 12C/13C and 16O/17O. More can be
   # added, but it hasn't been necessary so far...
   isolist = iso.split()
   for i in range(0,len(mass_list)):
      for atom in isolist:
         if i == int(atom)-1:
            if mass_list[i] == 1.00783: mass_list[i] = 2.0141 # protium - deuterium
            if mass_list[i] == 12.00000: mass_list[i] = 13.00335 # 12C - 13C
            if mass_list[i] == 15.99491: mass_list[i] = 16.9991 # 16O - 17O
   longline = ""
   for i in range(start_force,len(inlines)):
      longline = longline + inlines[i].rstrip().lstrip()

   forces = longline.split("NImag")[1].split('\\')[2].split(',')
   n = 0; l = -1
   for m in range (0, d_o_f):
      for n in range(0,m+1):
         l = l + 1
         Hmn = forces[l]
         sqrt_Mmn = (mass_list[m//3] * mass_list[n//3]) ** 0.5
         hess_mat[n,m] = float(Hmn)
         mw_hess_mat[m,n] = ((float(Hmn) / sqrt_Mmn))
         mw_hess_mat[n,m] = ((float(Hmn) / sqrt_Mmn))

   return mw_hess_mat
