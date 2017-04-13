#!/usr/bin/python

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

# PHYSICAL CONSTANTS
PLANCK_CONSTANT = 6.62606957e-34
BOLTZMANN_CONSTANT = 1.3806488e-23
SPEED_OF_LIGHT = 2.99792458e10
Eh = 4.359745e-18
a0 = 5.291772e-11
amu = 1.660468e-27

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

class calc_rpfr:
   #Computes the Reduced Isotopic Partition Function Ratio from a structure and a given isotopic substitution
   def __init__(self, file, isomer, *args):
      # Frequencies in waveunmbers
      frequency_wn = []

      # Extract the Force constants from a g09 logfile and generate the
      # mass-weighted Hessian matrix
      mw_hessmat = read_hess(file, isomer)

      # Convert from atomic units - a bit ugly
      conv_factor = Eh/a0**2/amu
      eigs = np.linalg.eigvalsh(mw_hessmat*conv_factor)
      fac = 5.30883671512e-12
      freqs = [ np.copysign(np.sqrt(np.abs(freq)),freq) * fac for freq in eigs ]

      # Keep frequencies larger in magnitude than the predefined cut-off
      # This should remove 5 or 6 normal modes (linear or non-linear molecule)
      # Add sanity check for removal?
      for freq in freqs:
         if np.abs(freq) > freq_cutoff:
             if freq > 0.0: frequency_wn.append(freq)
             else:  self.im_frequency_wn = -1.0 * freq

      # Calculate the excitation factor (EXC), the ZPE (ZPE) and Teller-Redlich product factor (PF)
      self.PF = calc_product_factor(frequency_wn, freq_scale_factor)
      self.ZPE = calc_zpe_factor(frequency_wn, temperature, freq_scale_factor)
      self.EXC = calc_excitation_factor(frequency_wn, temperature, freq_scale_factor)

if __name__ == "__main__":
   # Takes arguments: g09_output_files and optional temperature and vibrational scaling factor
   files = []
   # Default values
   temperature = 298.15; freq_scale_factor = 1.0; freq_cutoff = 50.0

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

      print "   KINISOT.py v 1.0.1", time.strftime("%Y-%m-%d %H:%M")
      print "   Temperature =", temperature, "K;   Frequency scale factor =", freq_scale_factor
      #print "   Frequency cut-off =", freq_cutoff, "wavenumbers"

   else:
      print "\nWrong number of arguments used. Correct format: Kinisot.py 2 x g09_output_files -iso 1 2 3 (-t <temp>) (-s <scalefactor>) (-cutoff <value>)\n"
      sys.exit()

   print "\n  ",
   print "Structure".ljust(32), "Imag. Freq (cm-1)".rjust(10)

   # Conversion from wavenumbers to SI energy units; then divide by kT
   tofreq = SPEED_OF_LIGHT * PLANCK_CONSTANT / BOLTZMANN_CONSTANT / temperature
   KIE = []

   # Calculates the RPFR for each for the reactant and its isotopomer and then the transition structure and the same isotopomer
   for iso in ["0",label]:
      for file in files:
         rpfr = calc_rpfr(file, iso, temperature, freq_scale_factor, freq_cutoff)
         print "o ",
         print file.split(".")[0],
         if iso != "0": print "- isotopic @", iso,
         if hasattr(rpfr, "im_frequency_wn"): print '{:27.2f}'.format(rpfr.im_frequency_wn)
         else: print "------".rjust(26)
         KIE.append(rpfr)

   # Fancy printing
   print "\n  ","".ljust(30), "V-ratio".rjust(10), "  ZPE".rjust(11), "   EXC".rjust(11), "   TRPF".rjust(11), "   KIE".rjust(11), "   1D-tunn".rjust(11), "   corr-KIE".rjust(11)
   print ("   TST-KIE @ "+str(temperature)+" K").ljust(30),

   if hasattr(KIE[1], "im_frequency_wn") and hasattr(KIE[3], "im_frequency_wn"):
      imfreq_fac = KIE[1].im_frequency_wn/KIE[3].im_frequency_wn
      print "     %.6f" % (imfreq_fac),

      # A correction factor for QM-tunneling (Bell infinite parabola)
      parabolic_tunn_corr = imfreq_fac * math.sin(0.5 * tofreq * KIE[3].im_frequency_wn) / math.sin(0.5 * tofreq * KIE[1].im_frequency_wn)
   else:
      print "o  Kinisot requires a transition structure with an imaginary frequencies!"
      sys.exit()

   # Application of the Bigeleisen-Mayer equation
   # https://en.wikipedia.org/wiki/Maria_Goeppert-Mayer
   # Compute separately the rate factors associated with differences in ZPE, vibrational entropy and then enthalpy...
   ZPE = (KIE[0].ZPE/KIE[2].ZPE) / (KIE[1].ZPE/KIE[3].ZPE)
   EXC = (KIE[0].EXC/KIE[2].EXC) / (KIE[1].EXC/KIE[3].EXC)
   TRPF = (KIE[1].PF/KIE[3].PF) / (KIE[0].PF/KIE[2].PF)
   print "   %.6f" % (ZPE), "   %.6f" % (EXC), "   %.6f" % (TRPF),

   # print out (a) the Bigeleisen-Mayer KIE with classical nuclei and then (b) a value corrected to include quantum tunneling effects...
   KIE = imfreq_fac * ZPE * EXC * TRPF
   KIE_Tunnel = KIE * parabolic_tunn_corr
   print "   %.6f" % (KIE), "   %.6f" % (parabolic_tunn_corr), "   %.6f" % (KIE_Tunnel)
