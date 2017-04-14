#!/usr/bin/python
from __future__ import print_function

# Comments and/or additions are welcome (send e-mail to:
# robert.paton@chem.ox.ac.uk

#####################################
#           Kinisot.py              #
#####################################
###  Written by:  Rob Paton #########
###  Last modified:  Mar 27, 2016 ###
#####################################

import sys, math, time
import numpy as np
from glob import glob
from Hess_to_Freq import *

# PHYSICAL CONSTANTS (SI apart from speed of light)
# These are identical to those used by G09/G16 (http://gaussian.com/constants/)
PLANCK_CONSTANT = 6.62606957e-34 #m2 kg / s
BOLTZMANN_CONSTANT = 1.3806488e-23 #m2 kg s-2 K-1
SPEED_OF_LIGHT = 2.99792458e10 #cm s-1
ENERGY_AU = 4.35974434e-18 #J
BOHR_RADIUS = 5.2917721092e-11 #m
ATOMIC_MASS_UNIT = 1.660538921e-27 #kg

# Enables output to terminal and to text file
class Logger:
   # Designated initializer
   def __init__(self,filein,suffix,append):
      # Create the log file at the input path
      self.log = open(filein+"_"+append+"."+suffix, 'w' )

   # Write a message to the log
   def Write(self, message):
      # print the message
      print(message, end=' ')
      # Write to log
      self.log.write(message)

   # Write a message only to the log and not to the terminal
   def Writeonlyfile(self, message):
      # Write to log
      self.log.write("\n"+message+"\n")

   # Write a fatal error, finalize and terminate the program
   def Fatal(self, message):
      # print the message
      print(message+"\n")
      # Write to log
      self.log.write(message + "\n")
      # Finalize the log
      self.Finalize()
      # End the program
      sys.exit(1)

   # Finalize the log file
   def Finalize(self):
      self.log.close()

def calc_product_factor(frequency_wn, freq_scale_factor):
   """
   Calculates the product of (scaled) vibrational frequencies in order to
   obtain the Teller-Redlich product factor. There is no temperture dependence
   to this term in the B-M equation
   """
   product = 1.0
   for freq in frequency_wn:
      product = product * freq * freq_scale_factor
   return product

def calc_zpe_factor(frequency_wn, temperature, freq_scale_factor):
   """
   Calculates the vibrational ZPE of (scaled) vibrational frequencies. ZPEs
   themselves are not temperature dependent although the exponential form of
   this term in the B-M equation is.
   """
   product = 1.0
   frequency = [entry * SPEED_OF_LIGHT for entry in frequency_wn]
   for entry in frequency:
      hv_over_kt = ((PLANCK_CONSTANT*entry*freq_scale_factor)/(BOLTZMANN_CONSTANT * temperature))
      product = product * math.exp(0.5 * hv_over_kt)
   return product

def calc_excitation_factor(frequency_wn, temperature, freq_scale_factor):
   """
   Calculates the excitation factor term of the RPFR from the (scaled) vibrational
   frequencies. This term is temperature dependent.
   """
   product = 1.0
   frequency = [entry * SPEED_OF_LIGHT for entry in frequency_wn]
   for entry in frequency:
      factor = ((PLANCK_CONSTANT*entry*freq_scale_factor)/(BOLTZMANN_CONSTANT*temperature))
      product = product * (1-math.exp(-factor))
   return product

def is_linear(file):
   # Check if a molecule is linear or not
   # This affects the number of rotational d.o.f.
   # Read gaussian output
   symm = 'none'
   g09_output = open(file, 'r')
   inlines = g09_output.readlines()

   for i in range(0,len(inlines)):
      if inlines[i].strip().find('prolate symmetric top') > -1: symm = 'linear'
      if inlines[i].strip().find('asymmetric top') > -1: symm = 'none'
   return symm

class calc_rpfr:
   #Computes the Reduced Isotopic Partition Function Ratio from a structure and a given isotopic substitution
   def __init__(self, file, isomer, *args):
      # Frequencies in waveunmbers
      frequency_wn = []

      # Extract the Force constants from a g09 logfile and generate the
      # mass-weighted Hessian matrix in Hartree/(amu Bohr^2)
      mw_hessmat = read_hess(file, isomer)

      # Convert from atomic units - a bit ugly
      unit_conversion = ENERGY_AU / (BOHR_RADIUS**2 * ATOMIC_MASS_UNIT) / ((SPEED_OF_LIGHT * 2 * np.pi)**2)
      eigs = np.linalg.eigvalsh(mw_hessmat * unit_conversion)
      freqs = [ np.copysign(np.sqrt(np.abs(freq)),freq) for freq in eigs ]

      # Keep frequencies larger in magnitude than the predefined cut-off
      # This should remove 5 or 6 normal modes (linear or non-linear molecule)
      trans_rot_modes = 0
      for freq in freqs:
         if np.abs(freq) > freq_cutoff:
             if freq > 0.0: frequency_wn.append(freq)
             else:  self.im_frequency_wn = -1.0 * freq
         else: trans_rot_modes = trans_rot_modes + 1

      # Check the number of low modes = 5/6 if molecule is non-linear/linear
      if is_linear(file) == 'linear':
          if trans_rot_modes != 5: log.Writeonlyfile("\nWARNING: removed", trans_rot_modes, "low modes. Expected 5 for linear molecule")
      if is_linear(file) == 'none':
          if trans_rot_modes != 6: log.Writeonlyfile("\nWARNING: removed", trans_rot_modes, "low modes. Expected 6 for non-linear molecule")

      # Calculate the excitation factor (EXC), the ZPE (ZPE) and Teller-Redlich product factor (PF)
      self.PF = calc_product_factor(frequency_wn, freq_scale_factor)
      self.ZPE = calc_zpe_factor(frequency_wn, temperature, freq_scale_factor)
      self.EXC = calc_excitation_factor(frequency_wn, temperature, freq_scale_factor)

if __name__ == "__main__":
   # Takes arguments: g09_output_files and optional temperature and vibrational scaling factor
   files = []
   # Default values
   temperature = 298.15; freq_scale_factor = 1.0; freq_cutoff = 50.0
   # Write an output file
   log = Logger("Kinisot","dat", "output")

   if len(sys.argv) > 1:
      for i in range(1,len(sys.argv)):
         if sys.argv[i] == "-t": temperature = float(sys.argv[i+1])
         elif sys.argv[i] == "-s": freq_scale_factor = float(sys.argv[i+1])
         elif sys.argv[i] == "-iso": label = (sys.argv[i+1])
         elif sys.argv[i] == "-cutoff": freq_cutoff = (sys.argv[i+1])
         else:
            if len(sys.argv[i].split(".")) > 1:
               if sys.argv[i].split(".")[1] == "out" or sys.argv[i].split(".")[1] == "log":
                   for file in glob(sys.argv[i]): files.append(file)

      space = "   "; dash = "--"; dash_line = space * 18 + dash * 40
      log.Write("\n" + (space * 18) + "KINISOT.py v 1.0.1" + time.strftime("%Y-%m-%d %H:%M"))
      log.Write("\n" + (space * 18) + "Temp. =" + str(temperature) + "K / Vib. scale factor =" + str(freq_scale_factor) + " / Low mode cut-off = " + str(freq_cutoff) + " cm-1")
      log.Write('\n' + dash_line)
      log.Write("\n    ".ljust(51) + "  V-ratio".rjust(12) + "  ZPE".rjust(12) + "   EXC".rjust(12) + "   TRPF".rjust(12) + "   KIE".rjust(12) + "   1D-tunn".rjust(12) + "   corr-KIE".rjust(12))

   else:
      log.Fatal("\nWrong number of arguments used. Correct format: Kinisot.py 2 x g09_output_files -iso 1 2 3 (-t <temp>) (-s <scalefactor>) (-cutoff <value>)\n")

   log.Write("\n   Structure".ljust(55) + dash *40)
   # Conversion from wavenumbers to SI energy units; then divide by kT
   tofreq = SPEED_OF_LIGHT * PLANCK_CONSTANT / BOLTZMANN_CONSTANT / temperature
   KIE = []

   # Calculates the RPFR for each for the reactant and its isotopomer and then the transition structure and the same isotopomer
   for iso in ["0",label]:
      for file in files:
         rpfr = calc_rpfr(file, iso, temperature, freq_scale_factor, freq_cutoff)
         log.Write("\no ",)
         if iso != "0": log.Write((file.split(".")[0]+": iso @ "+iso).ljust(47),)
         else: log.Write(file.split(".")[0].ljust(47),)
         if hasattr(rpfr, "im_frequency_wn"):
             log.Write('{:11.1f} {:11.2e} {:11.2e} {:11.2e}'.format(rpfr.im_frequency_wn, rpfr.ZPE, rpfr.EXC, rpfr.PF))
         else:
             log.Write('            {:11.2e} {:11.2e} {:11.2e}'.format(rpfr.ZPE, rpfr.EXC, rpfr.PF))
         KIE.append(rpfr)

   # Fancy log.Writeing
   log.Write('\n' + dash_line)
   log.Write(("\n   TST-KIE @ "+str(temperature)+" K").ljust(51),)

   if hasattr(KIE[1], "im_frequency_wn") and hasattr(KIE[3], "im_frequency_wn"):
      imfreq_fac = KIE[1].im_frequency_wn/KIE[3].im_frequency_wn
      log.Write('{:11.6f}'.format(imfreq_fac))

      # A correction factor for QM-tunneling (Bell infinite parabola)
      parabolic_tunn_corr = imfreq_fac * math.sin(0.5 * tofreq * KIE[3].im_frequency_wn) / math.sin(0.5 * tofreq * KIE[1].im_frequency_wn)
   else:
      log.Write("\no  Kinisot requires a transition structure with an imaginary frequencies!")
      sys.exit()

   # Application of the Bigeleisen-Mayer equation
   # Interesting: https://en.wikipedia.org/wiki/Maria_Goeppert-Mayer
   # Compute separately the rate factors associated with differences in ZPE, vibrational entropy and then enthalpy...
   ZPE = (KIE[0].ZPE/KIE[2].ZPE) / (KIE[1].ZPE/KIE[3].ZPE)
   EXC = (KIE[0].EXC/KIE[2].EXC) / (KIE[1].EXC/KIE[3].EXC)
   TRPF = (KIE[1].PF/KIE[3].PF) / (KIE[0].PF/KIE[2].PF)
   #log.Write("   %.6f" % (ZPE) + "   %.6f" % (EXC) + "   %.6f" % (TRPF))

   # log.Write out (a) the Bigeleisen-Mayer KIE with classical nuclei and then (b) a value corrected to include quantum tunneling effects...
   KIE = imfreq_fac * ZPE * EXC * TRPF
   KIE_Tunnel = KIE * parabolic_tunn_corr
   #log.Write("   %.6f" % (KIE) + "   %.6f" % (parabolic_tunn_corr) + "   %.6f" % (KIE_Tunnel))
   log.Write('{:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f}'.format(ZPE, EXC, TRPF, KIE, parabolic_tunn_corr, KIE_Tunnel))
   log.Write('\n' + dash_line + '\n')
