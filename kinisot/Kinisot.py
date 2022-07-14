#!/usr/bin/python
from __future__ import print_function

# Comments and/or additions are welcome (send e-mail to:
# robert.paton@colostate.edu

import pathlib, sys, math, time
import numpy as np
from glob import glob
from argparse import ArgumentParser
from kinisot.Hess_to_Freq import *
from kinisot.vib_scale_factors import scaling_data, scaling_refs

# version
__version__ = "2.0.1"

# PHYSICAL CONSTANTS (SI apart from speed of light)
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
   product = 0.0
   frequency = [entry * SPEED_OF_LIGHT for entry in frequency_wn]
   for entry in frequency:
      # changing to vibrational temperatures avoids big numbers. Everything is done logarithmically for the same reason throughout
      hv_over_k = PLANCK_CONSTANT * entry * freq_scale_factor / BOLTZMANN_CONSTANT
      product = product + np.log(hv_over_k)
   return product

def calc_zpe_factor(frequency_wn, temperature, freq_scale_factor):
   """
   Calculates the vibrational ZPE of (scaled) vibrational frequencies. ZPEs
   themselves are not temperature dependent although the exponential form of
   this term in the B-M equation is.
   """
   product = 0.0
   frequency = [entry * SPEED_OF_LIGHT for entry in frequency_wn]
   for entry in frequency:
      hv_over_kt = ((PLANCK_CONSTANT*entry*freq_scale_factor)/(BOLTZMANN_CONSTANT * temperature))
      product = product + np.log(math.exp(0.5 * hv_over_kt))
   return product

def calc_excitation_factor(frequency_wn, temperature, freq_scale_factor):
   """
   Calculates the excitation factor term of the RPFR from the (scaled) vibrational
   frequencies. This term is temperature dependent.
   """
   product = 0.0
   frequency = [entry * SPEED_OF_LIGHT for entry in frequency_wn]
   for entry in frequency:
      factor = ((PLANCK_CONSTANT*entry*freq_scale_factor)/(BOLTZMANN_CONSTANT*temperature))
      product = product + np.log(1-math.exp(-factor))
   return product

def level_of_theory(file):
   # Read gaussian output for the level of theory and basis set used
   g09_output = open(file, 'r')
   inlines = g09_output.readlines()
   level = "none"

   for i in range(0,len(inlines)):
      if inlines[i].strip().find('\\Freq\\') > -1:
          if len(inlines[i].strip().split("\\")) > 5:
              level = (inlines[i].strip().split("\\")[4])
              bs = (inlines[i].strip().split("\\")[5])
   return level+"/"+bs

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
   def __init__(self, file, isomer, temperature=298.15, freq_scale_factor=1.0, freq_cutoff=50.0):
      # Frequencies in waveunmbers
      self.frequency_wn = []

      # Extract the Force constants from a g09 logfile and generate the
      # mass-weighted Hessian matrix in Hartree/(amu Bohr^2)
      mw_hessmat = read_hess(file, isomer)

      # Convert from atomic units - a bit ugly
      unit_conversion = ENERGY_AU / (BOHR_RADIUS**2 * ATOMIC_MASS_UNIT) / ((SPEED_OF_LIGHT * 2 * np.pi)**2)
      eigs = np.linalg.eigvalsh(mw_hessmat * unit_conversion)
      freqs = [ np.copysign(np.sqrt(np.abs(freq)),freq) for freq in eigs ]

      # 5 or 6 small normal modes will be removed (depending on whether the molecule is linear or non-linear)
      if is_linear(file) == 'linear': trans_rot_modes = 5
      else: trans_rot_modes = 6

      # Keep a single imaginary frequency. It should be larger than the predefined cut-off
      if np.abs(freqs[0]) > freq_cutoff:
         self.im_frequency_wn = -1.0 * freqs[0]
         trans_rot_modes = trans_rot_modes + 1
      for freq in freqs[trans_rot_modes:]: self.frequency_wn.append(freq)

      # Calculate the excitation factor (EXC), the ZPE (ZPE) and Teller-Redlich product factor (PF)
      # returns a 1D-array of all terms
      self.PF = calc_product_factor(self.frequency_wn, freq_scale_factor)
      self.ZPE = calc_zpe_factor(self.frequency_wn, temperature, freq_scale_factor)
      self.EXC = calc_excitation_factor(self.frequency_wn, temperature, freq_scale_factor)


def main():
    # Parse Arguments
   parser = ArgumentParser()
   parser.add_argument("-t", dest="temperature", action="store", type=float, default=298.15, help="temperature in Kelvin (default 298.15K)")
   parser.add_argument("-s", dest="freq_scale_factor", action="store", type=float, default=False, help="scale factor for vibrations (default 1)")
   parser.add_argument("--iso", dest="label", action="store", default=False, help="atom number(s) of interest")
   parser.add_argument("--cutoff", dest="freq_cutoff", action="store", type=float, default=50.0, help="Frequency cutoff (default = 50 cm-1)")
   options, args = parser.parse_known_args()

   # Write an output file
   log = Logger("Kinisot","dat", "output")
   space = "   "; dash = "--"; dash_line = space * 17 + " " + dash * 37
   log.Write("\n  " + "KINISOT.py v " + __version__ + ": " + time.strftime("%Y-%m-%d %H:%M"))

   # Takes arguments: output_files and optional temperature and vibrational scaling factor
   files = []
   for arg in sys.argv:
       if pathlib.Path(arg).suffix in ['.out', '.log']:
           for file in glob(arg): files.append(file)

   if len(files) != 2:
      log.Fatal('   Kinisot requires two output files found! Exiting ...')
      sys.exit()

   if not options.label:
      log.Fatal('   Kinisot requires at least one atom to be labelled! Exiting ...')
      sys.exit()

   # Check the level of theory matches for the two files and then try to find
   # the relevant vibrational scaling factor
   l_o_t = []
   for file in files:
       l_o_t.append(level_of_theory(file))
   if l_o_t[0] != l_o_t[1]:
       log.Writeonlyfile("\nWARNING: found different levels of theory for reactant " + l_o_t[0] + " and TS " + l_o_t[1])
   else:
       level = l_o_t[0]
       for scal in scaling_data:
           if not options.freq_scale_factor:
               if level.upper().find(scal['level'].decode("utf-8") .upper()) > -1 or level.upper().find(scal['level'].decode("utf-8") .replace("-","").upper()) > -1:
                   log.Write("\n  " + "Found vibrational scaling factor " + str(scal['zpe_fac']) + " for " + level + " level of theory")
                   options.freq_scale_factor = scal['zpe_fac']
                   ref = scaling_refs[scal['zpe_ref']]
                   log.Write("\n  REF: " + ref)

   # If no match could be found then use 1.0 as default
   if not options.freq_scale_factor:
       options.freq_scale_factor = 1.00
       log.Write("\n" + (space * 18) + "Unable to find vibrational scaling factor for " + level + "; using value of 1.0")

   log.Write("\n\n" + (space * 17) + "  Temp = " + str(options.temperature) + "K / Vib. scale factor = " + str(options.freq_scale_factor))
   log.Write(("\n ").ljust(50))
   log.Write('{:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}'.format("V-ratio", "ZPE", "EXC", "TRPF", "KIE", "1D-tunn", "corr-KIE"))

# Calculates the RPFR for each for the reactant and its isotopomer and then the transition structure and the same isotopomer
   #options.freq_scale_factor = float(options.freq_scale_factor)
   KIE = []
   for iso in ["0",options.label]:
      for file in files:
         rpfr = calc_rpfr(file, iso, options.temperature, options.freq_scale_factor, options.freq_cutoff)
         KIE.append(rpfr)

   # Check for the presence of an imaginary frequency in second structure; exit gracefully if not.
   if hasattr(KIE[1], "im_frequency_wn") and hasattr(KIE[3], "im_frequency_wn"):
      imfreq_fac = KIE[1].im_frequency_wn/KIE[3].im_frequency_wn
   else: log.Fatal("\no  Kinisot requires a transition structure with an imaginary frequencies!")

   log.Write("\no " + files[0].split(".")[0].ljust(47) + "   " + dash * 37)
   log.Write("\no " + files[1].split(".")[0].ljust(47))
   log.Write('{:10.1f}'.format(KIE[1].im_frequency_wn))
   log.Write("\no " + (files[0].split(".")[0]+": iso @ "+iso).ljust(47))
   log.Write('           {:10.3e} {:10.3e} {:10.3e}'.format(np.e ** (KIE[0].ZPE - KIE[2].ZPE), np.e **(KIE[0].EXC - KIE[2].EXC), np.e ** (KIE[1].PF - KIE[3].PF)))
   log.Write("\no " + (files[1].split(".")[0]+": iso @ "+iso).ljust(47))
   log.Write('{:10.1f} {:10.3e} {:10.3e} {:10.3e}'.format(KIE[3].im_frequency_wn, np.e ** (KIE[1].ZPE - KIE[3].ZPE), np.e ** (KIE[1].EXC - KIE[3].EXC), np.e ** (KIE[0].PF - KIE[2].PF)))

   # Application of the Bigeleisen-Mayer equation
   # Interesting: https://en.wikipedia.org/wiki/Maria_Goeppert-Mayer
   # Compute separately the rate factors associated with differences in ZPE, vibrational entropy and then enthalpy...
   ZPE = np.e ** (KIE[0].ZPE - KIE[2].ZPE - KIE[1].ZPE + KIE[3].ZPE)
   EXC = np.e **(KIE[0].EXC - KIE[2].EXC - KIE[1].EXC + KIE[3].EXC)
   TRPF = np.e ** (KIE[1].PF - KIE[3].PF - KIE[0].PF + KIE[2].PF)

   # A correction factor for QM-tunneling (Bell infinite parabola)
   # Conversion from wavenumbers to SI energy units; then divide by kT
   tofreq = SPEED_OF_LIGHT * PLANCK_CONSTANT / BOLTZMANN_CONSTANT / options.temperature
   parabolic_tunn_corr = imfreq_fac * math.sin(0.5 * tofreq * KIE[3].im_frequency_wn) / math.sin(0.5 * tofreq * KIE[1].im_frequency_wn)

   # (a) the Bigeleisen-Mayer KIE with classical nuclei and (b) a value corrected to include quantum tunneling effects...
   KIE = imfreq_fac * ZPE * EXC * TRPF
   KIE_Tunnel = KIE * parabolic_tunn_corr

  # Fancy log.Writing
   log.Write('\n' + dash_line)
   log.Write(("\n  KIE @ "+str(options.temperature)+" K").ljust(50))
   log.Write('{:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f}'.format(imfreq_fac, ZPE, EXC, TRPF, KIE, parabolic_tunn_corr, KIE_Tunnel))
   log.Write('\n' + dash_line + '\n')
   log.Finalize()

if __name__ == "__main__":
   main()
