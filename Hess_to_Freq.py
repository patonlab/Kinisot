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

#######################################################################
#                             Hess_to_Freq                            #
#######################################################################
#######  Written by:  Rob Paton #######################################
#######  Last modified:  Mar 20, 2013 #################################
#######################################################################

import sys, math
import numpy as np

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


def read_hess(file, freq_scale_factor,iso):

		# Force constants are scaled by freq_scale_factor squared
		force_const_scale_factor = freq_scale_factor ** 2

		# Read gaussian output
		g09_output = open(file, 'r')
		
		inlines = g09_output.readlines()
		
                mass_list = []
		start_force = len(inlines) < 100
		for i in range(0,len(inlines)):
                       	if inlines[i].strip().find('and mass') > -1:
                                mass_list.append(float(inlines[i].strip().split()[8]))

			if inlines[i].strip().startswith('Deg. of freedom'):
				d_o_f = int(inlines[i].strip().split()[3]) + 6

			if inlines[i].strip().startswith('The second derivative matrix:'):  start_hess = i + 2
			if inlines[i].strip().startswith('ITU='): end_hess = i
	
			if inlines[i].strip().find('NImag') > -1: start_force = i
		
		hess_mat = np.ndarray(shape=(d_o_f,d_o_f))
	
		isolist = iso.split()
		for i in range(0,len(mass_list)):
			for atom in isolist:
				if i == int(atom)-1:
					#print "increasing the mass of atom", iso 
					if mass_list[i] == 1.00783: mass_list[i] = 2.0141
					if mass_list[i] == 12.00000: mass_list[i] = 13.00335
					if mass_list[i] == 15.99491: mass_list[i] = 16.9991
		#print start_force, len(inlines)
		longline = ""
		for i in range(start_force,len(inlines)):
			longline = longline + inlines[i].rstrip().lstrip()
		
		forces = longline.split("NImag")[1].split('\\')[2].split(',')
		
		n = 0; l = -1
		for m in range (0, d_o_f):
			for n in range(0,m+1):
				l = l + 1
				#print m, n, forces[l],
				Hmn = forces[l]
	
				sqrt_Mmn = mass_list[m/3] * mass_list[n/3] 
				sqrt_Mmn = sqrt_Mmn ** 0.5
				GDF = 435.942 / ( 5.29167 * 5.29167 )
				hess_mat[m,n] = GDF * ((float(Hmn) / sqrt_Mmn))
				hess_mat[n,m] = GDF * ((float(Hmn) / sqrt_Mmn))
		
			
		#print hess_mat
		return hess_mat

if __name__ == "__main__":
	
	# Takes arguments: g09_output_files and optional temperature and vibrational scaling factor
	files = []
	freq_scale_factor = 1.0

	if len(sys.argv) > 1:
		for i in range(1,len(sys.argv)):
			if sys.argv[i] == "-s": freq_scale_factor = float(sys.argv[i+1])
			else:
				if len(sys.argv[i].split(".")) > 1:
					if sys.argv[i].split(".")[1] == "out" or sys.argv[i].split(".")[1] == "log": files.append(sys.argv[i])
			
		print "   Frequency scale factor =", freq_scale_factor

	else:
		print "\nWrong number of arguments used. Correct format: Hess_to_Freq.py (-s scalefactor) g09_output_file(s)\n"
		sys.exit()
	

	for file in files:
		SPEED_OF_LIGHT = 2.99792458e10
		AVOGADRO_CONSTANT = 6.0221415e23
		AMU_to_KG = 1.66053886E-27
		CNV = 5.8890E-07

		mw_hessmat = read_hess(file, freq_scale_factor, 0)
		#print (mw_hessmat)
		eigs = np.linalg.eigvalsh(mw_hessmat)
		for eig in eigs: 
			if eig > 0: print (eig / CNV) ** 0.5 
			else: print -1 * (-eig/ CNV) ** 0.5 

