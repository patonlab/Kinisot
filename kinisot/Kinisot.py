#!/usr/bin/python
from __future__ import print_function

# Comments and/or additions are welcome (send e-mail to:
# robert.paton@colostate.edu

import pathlib, sys, math, time
import numpy as np
from glob import glob
from argparse import ArgumentParser

# Importing regardless of relative import
try:
    from .vib_scale_factors import scaling_data, scaling_refs
    from .Hess_to_Freq import *
except:
    from vib_scale_factors import scaling_data, scaling_refs
    from Hess_to_Freq import *

# version
__version__ = "2.0.2"

# PHYSICAL CONSTANTS (SI apart from speed of light)
PLANCK_CONSTANT = 6.62606957e-34 #m2 kg / s
BOLTZMANN_CONSTANT = 1.3806488e-23 #m2 kg s-2 K-1
SPEED_OF_LIGHT = 2.99792458e10 #cm s-1
ENERGY_AU = 4.35974434e-18 #J
BOHR_RADIUS = 5.2917721092e-11 #m
ATOMIC_MASS_UNIT = 1.660538921e-27 #kg

# PHYSICAL CONSTANTS (SI apart from speed of light)
PLANCK_CONSTANT = 6.626070e-34 #m2 kg / s
BOLTZMANN_CONSTANT = 1.380649e-23 #m2 kg s-2 K-1
SPEED_OF_LIGHT = 2.997925e10 #cm s-1
ENERGY_AU = 4.359745e-18 #J
BOHR_RADIUS = 5.291772e-11 #m
ATOMIC_MASS_UNIT = 1.660468e-27 #kg

# print formatting
space = "   "; dash = "--"; dash_line = space * 17 + " " + dash * 37

# Enables output to terminal and to text file
class Logger:
   # Designated initializer
   def __init__(self,filein,suffix,append):
      # Create the log file at the input path
      self.log = open(filein+"_"+append+"."+suffix, 'w' )

   # Write a message to the log
   def Write(self, message):
      # print the message
      print(message, end='')
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

def get_frequency_scaling(files, log):
   # Check the level of theory matches for all files and then try to find
   # the relevant vibrational scaling factor
   freq_scale_factor, level = 1.00, "unknown"     
   l_o_t = []
   for file in files:
      l_o_t.append(level_of_theory(file))
   if l_o_t[0] != l_o_t[1]:
      log.Writeonlyfile("\nWARNING: found different levels of theory for reactant " + l_o_t[0] + " and TS " + l_o_t[1])
   else:
      level = l_o_t[0]
      for scal in scaling_data:
         
         if level.upper().find(scal['level'].decode("utf-8") .upper()) > -1 or level.upper().find(scal['level'].decode("utf-8") .replace("-","").upper()) > -1:
            log.Write("\n  " + "Found vibrational scaling factor " + str(scal['zpe_fac']) + " for " + level + " level of theory")
            freq_scale_factor = scal['zpe_fac']
            ref = scaling_refs[scal['zpe_ref']]
            log.Write("\n  REF: " + ref)
   
   if freq_scale_factor == 1.00:
      log.Write("\n  Unable to find vibrational scaling factor for " + level + "; using value of 1.0")
  
   return freq_scale_factor

def calc_product_factor(frequency_wn):
   """
   Calculates the product of (scaled) vibrational frequencies in order to
   obtain the Teller-Redlich product factor. There is no temperature dependence
   to this term in the BM equation
   """
   product = 0.0
   frequency = [entry * SPEED_OF_LIGHT for entry in frequency_wn]

   for entry in frequency:
      # changing to vibrational temperatures avoids big numbers. Everything is done logarithmically for the same reason throughout
      hv_over_k = PLANCK_CONSTANT * entry / BOLTZMANN_CONSTANT
      product += np.log(hv_over_k)
   return product

def calc_zpe_factor(frequency_wn, temperature):
   """
   Calculates the vibrational ZPE of (scaled) vibrational frequencies. ZPEs
   themselves are not temperature dependent although the exponential form of
   this term in the BM equation is.
   """
   product = 0.0
   frequency = [entry * SPEED_OF_LIGHT for entry in frequency_wn]
   for entry in frequency:
      hv_over_kt = ((PLANCK_CONSTANT*entry)/(BOLTZMANN_CONSTANT * temperature))
      product += np.log(math.exp(0.5 * hv_over_kt))
   return product

def calc_excitation_factor(frequency_wn, temperature):
   """
   Calculates the excitation factor term of the RPFR from the (scaled) vibrational
   frequencies. This term is temperature dependent.
   """
   product = 0.0
   frequency = [entry * SPEED_OF_LIGHT for entry in frequency_wn]
   for entry in frequency:
      factor = ((PLANCK_CONSTANT*entry)/(BOLTZMANN_CONSTANT*temperature))
      product += np.log(1-math.exp(-factor))
   return product

class calc_rpfr:
   #Computes the Reduced Isotopic Partition Function Ratio from a structure and a given isotopic substitution
   def __init__(self, files, isomer, temperature=298.15, freq_scale_factor=1.0, freq_cutoff=50.0):

      self.PF, self.ZPE, self.EXC = 0.0, 0.0, 0.0

      for i, file in enumerate(files):
          # Frequencies in waveunmbers
          self.frequency_wn = []

          # Extract the Force constants from a g09 logfile and generate the
          # mass-weighted Hessian matrix in Hartree/(amu Bohr^2)
          mw_hessmat = read_hess(file, isomer[i])
          
          # Convert from atomic units - a bit ugly
          unit_conversion = ENERGY_AU / (BOHR_RADIUS**2 * ATOMIC_MASS_UNIT) / ((SPEED_OF_LIGHT * 2 * np.pi)**2)
          eigs = np.linalg.eigvalsh(mw_hessmat * unit_conversion)
          freqs = [ np.copysign(np.sqrt(np.abs(freq)),freq) * freq_scale_factor for freq in eigs ]

          # 5 or 6 small normal modes will be removed (depending on whether the molecule is linear or non-linear)
          if is_linear(file) == 'linear': trans_rot_modes = 5
          else: trans_rot_modes = 6

          # Keep a single imaginary frequency. It should be larger than the predefined cut-off
          if np.abs(freqs[0]) > freq_cutoff:
             self.im_frequency_wn = -1.0 * freqs[0]
             trans_rot_modes = trans_rot_modes + 1
          for freq in freqs[trans_rot_modes:]: self.frequency_wn.append(freq)

          wns = [("%0.2f") % wn for wn in self.frequency_wn]
          wns = [float(wn) for wn in wns]
          #print(file, isomer, wns)

          # Calculate the excitation factor (EXC), the ZPE (ZPE) and Teller-Redlich product factor (PF)
          # returns a 1D-array of all terms
          self.PF += calc_product_factor(self.frequency_wn)
          self.ZPE += calc_zpe_factor(self.frequency_wn, temperature)
          self.EXC += calc_excitation_factor(self.frequency_wn, temperature)

def compute_isotope_effect(rct, ts, prd, label, temperature=298.15, freq_scale_factor=1.0, freq_cutoff=50.0):
   # Calculates the RPFR for each species and its isotopomer
   KIE = []

   for iso in [['0'] * len(rct), label[0:len(rct)]]:
       rpfr = calc_rpfr(rct, iso, temperature, freq_scale_factor, freq_cutoff)
       KIE.append(rpfr)

   if ts != None:
       for iso in [['0'] * len(ts), label[len(rct):]]:
           rpfr = calc_rpfr(ts, iso, temperature, freq_scale_factor, freq_cutoff)
           KIE.append(rpfr)

       # Check for the presence of an imaginary frequency in second structure; exit gracefully if not.
       if hasattr(KIE[2], "im_frequency_wn") and hasattr(KIE[3], "im_frequency_wn"):
           freq_fac = KIE[2].im_frequency_wn/KIE[3].im_frequency_wn
       else: log.Fatal("\no  Kinisot requires a transition structure with an imaginary frequencies!")

   elif prd != None:
     for iso in [['0'] * len(prd), label[len(rct):]]:
         rpfr = calc_rpfr(prd, iso, temperature, freq_scale_factor, freq_cutoff)
         KIE.append(rpfr)
     freq_fac = 1.0

   # Application of the Bigeleisen-Mayer equation
   ZPE = np.e ** (KIE[0].ZPE - KIE[1].ZPE - KIE[2].ZPE + KIE[3].ZPE)
   EXC = np.e **(KIE[0].EXC - KIE[1].EXC - KIE[2].EXC + KIE[3].EXC)
   TRPF = np.e ** (KIE[2].PF - KIE[3].PF - KIE[0].PF + KIE[1].PF)

   # A correction factor for QM-tunneling (Bell infinite parabola)
   # Conversion from wavenumbers to SI energy units; then divide by kT
   tofreq = SPEED_OF_LIGHT * PLANCK_CONSTANT / BOLTZMANN_CONSTANT / temperature
   if ts != None: parabolic_tunn_corr = freq_fac * math.sin(0.5 * tofreq * KIE[3].im_frequency_wn) / math.sin(0.5 * tofreq * KIE[2].im_frequency_wn)
   else: parabolic_tunn_corr = 1.0

   # (a) the Bigeleisen-Mayer KIE with classical nuclei and (b) a value corrected to include quantum tunneling effects...
   KIE_no_tunnel = freq_fac * ZPE * EXC * TRPF
   KIE_tunnel = KIE_no_tunnel * parabolic_tunn_corr

   return KIE, ZPE, EXC, TRPF, KIE_no_tunnel, KIE_tunnel, parabolic_tunn_corr, freq_fac


def main():
    # Parse Arguments
   parser = ArgumentParser()
   parser.add_argument("-t", dest="temperature", action="store", type=float, default=298.15, help="temperature in Kelvin (default 298.15K)")
   parser.add_argument("-s", dest="freq_scale_factor", action="store", type=float, default=False, help="scale factor for vibrations (default 1)")
   parser.add_argument("--iso", dest="label", action='append', required=True, help="atom number(s) of interest")
   parser.add_argument("--cutoff", dest="freq_cutoff", action="store", type=float, default=50.0, help="Frequency cutoff (default = 50 cm-1)")
   parser.add_argument("--rct", dest="rct", action='append', required=True, help="Reactant logfile")
   parser.add_argument("--prd", dest="prd", action='append', help="Product logfile (for EQE calculation)")
   parser.add_argument("--ts", dest="ts", action='append', help="TS logfile (for KIE calculation)")
   
   options, args = parser.parse_known_args()
   log = Logger("Kinisot","dat", "output")

   if options.ts == None and options.prd == None:
       log.Fatal('   Kinisot requires either a TS for KIE or a product for EQE! Exiting ...')
       sys.exit()

   # Write an output file
   log.Write("\n  " + "KINISOT.py v " + __version__ + ": " + time.strftime("%Y-%m-%d %H:%M") + "\n")

   # Takes arguments: output_files and optional temperature and vibrational scaling factor
   if options.ts != None: files = options.rct + options.ts
   elif options.prd != None: files = options.rct + options.prd

   # if only one set of labels is provided, assume that the atom numbering is the same for rct and ts or rct and prd
   if len(options.label) == 1: options.label = options.label * 2
   
   for i, species in enumerate(files):
      print("  Species: {} isotopologue: {}".format(species, options.label[i]))
    
   if len(files) != len(options.label):
       log.Fatal("\no  For multiple reactants you need to specify the labels in each!")
   
   # if not specified try to automatically determine the vibrational scaling factor
   if not options.freq_scale_factor: 
      options.freq_scale_factor = get_frequency_scaling(files, log)

   log.Write("\n\n" + (space * 17) + "  Temp = " + str(options.temperature) + "K / Vib. scale factor = " + str(options.freq_scale_factor))
   log.Write(("\n  ").ljust(50))
   log.Write(' {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} \n'.format("V-ratio", "ZPE", "EXC", "TRPF", "KIE", "1D-tunn", "corr-KIE"))

   # Here are the ingredients and final predictions of the isotope effect
   KIE, ZPE, EXC, TRPF, KIE_no_tunnel, KIE_tunnel, parabolic_tunn_corr, freq_fac = compute_isotope_effect(options.rct, options.ts, options.prd, options.label, options.temperature, options.freq_scale_factor, options.freq_cutoff)

   # Fancy log.Writing
   log.Write("\no " + options.rct[0].split(".")[0].ljust(47) + "   " + dash * 37)
   if options.ts != None: log.Write("\no " + options.ts[0].split(".")[0].ljust(47))
   elif options.prd != None: log.Write("\no " + options.prd[0].split(".")[0].ljust(47))
   if options.ts != None: log.Write('{:10.1f}'.format(KIE[2].im_frequency_wn))
   log.Write("\no " + (options.rct[0].split(".")[0]+": iso @ "+' / '.join(options.label[0:len(options.rct)])).ljust(47))
   log.Write('           {:10.3e} {:10.3e} {:10.3e}'.format(np.e ** (KIE[0].ZPE - KIE[1].ZPE), np.e **(KIE[0].EXC - KIE[1].EXC), np.e ** (KIE[2].PF - KIE[3].PF)))
   if options.ts != None: log.Write("\no " + (options.ts[0].split(".")[0]+": iso @ "+' / '.join(options.label[len(options.rct):])).ljust(47))
   elif options.prd != None: log.Write("\no " + (options.prd[0].split(".")[0]+": iso @ "+' / '.join(options.label[len(options.rct):])).ljust(47))
   if options.ts != None: log.Write('{:10.1f} {:10.3e} {:10.3e} {:10.3e}'.format(KIE[3].im_frequency_wn, np.e ** (KIE[2].ZPE - KIE[3].ZPE), np.e ** (KIE[2].EXC - KIE[3].EXC), np.e ** (KIE[0].PF - KIE[1].PF)))
   else: log.Write('{:21.3e} {:10.3e} {:10.3e}'.format(np.e ** (KIE[2].ZPE - KIE[3].ZPE), np.e ** (KIE[2].EXC - KIE[3].EXC), np.e ** (KIE[0].PF - KIE[1].PF)))

   log.Write('\n' + dash_line)
   log.Write(("\n  KIE @ "+str(options.temperature)+" K").ljust(50))
   if options.ts != None: log.Write('{:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f}'.format(freq_fac, ZPE, EXC, TRPF, KIE_no_tunnel, parabolic_tunn_corr, KIE_tunnel))
   else: log.Write('{:21.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f}'.format(ZPE, EXC, TRPF, KIE_no_tunnel, parabolic_tunn_corr, KIE_tunnel))
   log.Write('\n' + dash_line + '\n')

if __name__ == "__main__":
   main()
