#!/usr/bin/python

# THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# Comments and/or additions are welcome (send e-mail to:
# robert.paton@chem.ox.ac.uk

#####################################
#           Kinisot.py              #
#####################################
###  Written by:  Rob Paton #########
###  Last modified:  Mar 27, 2016 ###
#####################################

import sys, math
import numpy as np
from Hess_to_Freq import *

# PHYSICAL CONSTANTS
GAS_CONSTANT = 8.3144621
PLANCK_CONSTANT = 6.62606957e-34
BOLTZMANN_CONSTANT = 1.3806488e-23
SPEED_OF_LIGHT = 2.99792458e10
AVOGADRO_CONSTANT = 6.0221415e23
AMU_to_KG = 1.66053886E-27
autokcal = 627.509541
kjtokcal = 4.184
atmos = 101.325

# vibrational energy evaluation (depends on frequencies, temp and scaling factor)
def calc_vibrational_energy(frequency_wn, temperature,freq_scale_factor):
   """
   Calculates the vibrational energy contribution (kcal/mol)
   Includes ZPE (0K) and thermal contributions
   Evib = R * Sum(0.5 hv/k + 1/(e^(hv/KT)-1))
   """
   energy = 1.0
   frequency = [entry * SPEED_OF_LIGHT for entry in frequency_wn]
   for entry in frequency:
      factor = ((PLANCK_CONSTANT*entry*freq_scale_factor)/(BOLTZMANN_CONSTANT*temperature))
      energy = energy * factor
   return energy


# vibrational zero point energy evaluation (depends on frequencies and scaling factor)
def calc_zeropoint_energy(frequency_wn,freq_scale_factor):
   """
   Calculates the vibrational ZPE (kcal/mol)
   EZPE = Sum(0.5 hv/k)
   """
   energy = 0.0
   frequency = [entry * SPEED_OF_LIGHT for entry in frequency_wn]
   for entry in frequency:
      factor = ((PLANCK_CONSTANT*entry*freq_scale_factor)/(BOLTZMANN_CONSTANT))
      temp = 0.5*factor
      temp = temp*GAS_CONSTANT
      energy = energy + temp
   energy = energy/kjtokcal/1000.0
   return energy


# rigid-rotor harmonic oscillator entropy evaluation
def calc_harmonic_entropy(frequency_wn, temperature,freq_scale_factor):
   """
   Calculates the entropic contribution (cal/(mol*K)) of a harmonic oscillator for
   a list of frequencies of vibrational modes
   Sv = R(hv/(k(e^(hv/KT)-1) - ln(1-e^(-hv/kT)))   
   """
   entropy = 1.0
   frequency = [entry * SPEED_OF_LIGHT for entry in frequency_wn]
   for entry in frequency:
      factor = ((PLANCK_CONSTANT*entry*freq_scale_factor)/(BOLTZMANN_CONSTANT*temperature))
      temp = 1-math.exp(-factor)
      entropy = entropy * temp
   return entropy


class calc_rpfr:   
   def __init__(self, file, temperature, freq_scale_factor,isomer):

      # Frequencies in waveunmbers
      frequency_wn = []
      
      # Read commandline arguments
      g09_output = open(file, 'r')
   
      mw_hessmat = read_hess(file, freq_scale_factor,isomer)
      eigs = np.linalg.eigvalsh(mw_hessmat)
      for eig in eigs:
         if eig > 0.00147: frequency_wn.append((eig / 5.8890E-07) ** 0.5)
         if eig < -0.00147: self.im_frequency_wn = (-1.0 * eig / 5.8890E-07) ** 0.5
      
      #print frequency_wn
   
      # Calculate vibrational contributions to the energy
      Uvib = calc_vibrational_energy(frequency_wn, temperature,freq_scale_factor)
      ZPE = calc_zeropoint_energy(frequency_wn, freq_scale_factor)
      # Calculate harmonic entropy for each frequency - functions defined above
      Svib = calc_harmonic_entropy(frequency_wn, temperature,freq_scale_factor)
      
      # Collect all terms 
      self.enthalpy = Uvib/autokcal
      self.zpe = ZPE/autokcal
      self.entropy = Svib/autokcal/1000.0
      self.gibbs_free_energy = self.enthalpy - temperature * self.entropy
   
if __name__ == "__main__":
   
   # Takes arguments: g09_output_files and optional temperature and vibrational scaling factor
   files = []
   temperature = "none"; freq_scale_factor = "none"

   if len(sys.argv) > 1:
      for i in range(1,len(sys.argv)):
         if sys.argv[i] == "-t": temperature = float(sys.argv[i+1])
         elif sys.argv[i] == "-s": freq_scale_factor = float(sys.argv[i+1])
         elif sys.argv[i] == "-iso": label = (sys.argv[i+1])
         else:
            if len(sys.argv[i].split(".")) > 1:
               if sys.argv[i].split(".")[1] == "out" or sys.argv[i].split(".")[1] == "log": files.append(sys.argv[i])
         
      if temperature != "none": print "   Temperature =", temperature, "K;",
      else: print "   Temperature (default) = 298.15K;",; temperature = 298.15
      if freq_scale_factor != "none": print "   Frequency scale factor =", freq_scale_factor
      else: print "   Frequency scale factor (default) = 1.0"; freq_scale_factor = 1.0
      
   else:
      print "\nWrong number of arguments used. Correct format: Kinisot.py (-t temp) (-s scalefactor) 2 x g09_output_files\n"
      sys.exit()
   
   print "\n  ",
   print "Structure".ljust(60), "Im Freq".rjust(10), "   ZPE".rjust(11)
   Standard = []
   beta = 1000.0 * autokcal * kjtokcal / GAS_CONSTANT / temperature

   isotopomers = ["0",label]
   KIE = []

   for iso in isotopomers:
      for file in files:
         rpfr = calc_rpfr(file, temperature, freq_scale_factor,iso)
         print "o ",
         print (file.split(".")[0]+": label @ "+str(iso)).ljust(60),
         if hasattr(rpfr, "im_frequency_wn"): print "%.6f" % (rpfr.im_frequency_wn),
         else: print "N/A".rjust(10),
         if hasattr(rpfr, "zpe"): print "   %.6f" % (rpfr.zpe)
         else: print "N/A".rjust(10)
         KIE.append(rpfr)
      
   
   print "\n  ","".ljust(30), "V-ratio".rjust(10), "  ZPE-ratio".rjust(11), "   Exc".rjust(11), "   TR-Prod".rjust(11), "   KIE".rjust(11), "   Tunn".rjust(11), "   corr-KIE".rjust(11)
   print ("   TST-KIE @ "+str(temperature)+" K").ljust(30),

   if hasattr(KIE[1], "im_frequency_wn") and hasattr(KIE[3], "im_frequency_wn"):
      imfreq_fac = KIE[1].im_frequency_wn/KIE[3].im_frequency_wn
      print "     %.6f" % (imfreq_fac),
   
      tofreq = SPEED_OF_LIGHT * PLANCK_CONSTANT / BOLTZMANN_CONSTANT / temperature
   
      # A correction factor for QM-tunneling (infinite parabola)
      parabolic_tunn_corr = imfreq_fac * math.sin(0.5 * tofreq * KIE[3].im_frequency_wn) / math.sin(0.5 * tofreq * KIE[1].im_frequency_wn)
   else:
      print "o  Require two transition structures with imaginary frequencies!"
      sys.exit()      

   # Compute separately the rate factors associated with differences in ZPE, vibrational entropy and then enthalpy...
   zpe_fac = math.exp(beta*((KIE[0].zpe-KIE[2].zpe)-(KIE[1].zpe-KIE[3].zpe)))
   print "   %.6f" % (zpe_fac),
   exp_fac = (KIE[0].entropy/KIE[2].entropy) / (KIE[1].entropy/KIE[3].entropy)     
   print "   %.6f" % (exp_fac),
   freq_fac = (KIE[1].enthalpy/KIE[3].enthalpy) / (KIE[0].enthalpy/KIE[2].enthalpy)
   print "   %.6f" % (freq_fac),
   
   # print out the KIE from TST with classical nuclei and then a value corrected to include quantum tunneling effects...
   kie_tot = imfreq_fac * zpe_fac * exp_fac * freq_fac
   corr_kie_tot = kie_tot * parabolic_tunn_corr
   print "   %.6f" % (kie_tot), "   %.6f" % (parabolic_tunn_corr), "   %.6f" % (corr_kie_tot)

