# -*- coding: utf-8 -*-
"""
Created on Sat May  9 15:12:11 2020

@author: willi
"""

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from numba import jitclass
import numba
from scipy.optimize import curve_fit

def calculate_reducde_chi(e,o,err):
	"""
	

	Parameters
	----------
	e : 1d np array
		expected value
	o : 1d np array
		obatined value
	err : 1d np array
		standard deviation

	Returns
	-------
	TYPE
		reduced chi squared

	"""
	total = 0
	for i in range(e.size):
		total += ((e[i]-o[i])/err[i])**2
	return total/(o.size-1)

#sets up numba variable conversion
spec = [('absorbtion_area_1', numba.float64),
        ('scattering_area_1', numba.float64),
		('density_1', numba.float64),
		('mass_1', numba.float64),
		('absorbtion_area', numba.float64[:]),
        ('scattering_area', numba.float64[:]),
		('density', numba.float64[:]),
		('mass', numba.float64[:]),
		('absorbtion_area_prime', numba.float64),
        ('scattering_area_prime', numba.float64),
		('density_prime', numba.float64),
		('mass_prime', numba.float64),
		('material_count',numba.int_)]

@jitclass(spec)
class experiment_neutrinoes:
	"""
	Sets up materials to conduct experients on
	"""
	
	def __init__(self, absorbtion_area_1, scattering_area_1, density_1, mass_1):
		self.material_count = 0
		self.absorbtion_area = np.zeros(3)
		self.scattering_area = np.zeros(3)
		self.density = np.zeros(3)
		self.mass = np.zeros(3)
		
		self.absorbtion_area[self.material_count] = absorbtion_area_1
		self.scattering_area[self.material_count] = scattering_area_1
		self.density[self.material_count]= density_1
		self.mass[self.material_count] = mass_1
		self.material_count = 1
	
	def define_another_material(self, absorbtion_area_prime,
							 scattering_area_prime, density_prime, mass_prime):
		"""
		Allows defining up to 3 materials

		Parameters
		----------
		absorbtion_area_prime : float
		scattering_area_prime : float
		density_prime : float
		mass_prime : float

		Returns
		-------
		None.

		"""
		if self.material_count == 3:
			print("Only supports 3 materials")
			return
		self.absorbtion_area[self.material_count] = absorbtion_area_prime
		self.scattering_area[self.material_count] = scattering_area_prime
		self.density[self.material_count] = density_prime
		self.mass[self.material_count] = mass_prime
		self.material_count += 1
	
	def get_cartesians(self,vec):
		"""
		Calculates x,y,z from polar cordiantes

		Parameters
		----------
		vec : 2d np array
			Array of polar vectors

		Returns
		-------
		vec : 2d np array
			Array of cartesian vectors

		"""
		for i in range(vec[:,0].size):
			x = vec[i][0]*np.sin(vec[i][1])*np.cos(vec[i][2])
			y = vec[i][0]*np.sin(vec[i][1])*np.sin(vec[i][2])
			z = vec[i][0]*np.cos(vec[i][1])
			vec[i] = np.array([x,y,z])
			
		return vec
	
	def generate_distances(self,num, mean_free_path):
		"""
		Generates a random distances

		Parameters
		----------
		num : int
			number of varibles to be generated
		mean_free_path : float

		Returns
		-------
		1d np array
			Random distance

		"""
		rand_samples = np.array(
				[np.random.sample() for x in range(num)])
		return -mean_free_path*np.log(rand_samples)
	
	def generate_polar_vectors(self, num):
		"""
		Generates a random unit vectors

		Parameters
		----------
		num : int
			number of varibles to be generated

		Returns
		-------
		2d np array
			Random unit vectors

		"""
		vectors = np.zeros((num,3))
		
		for i in range(num):
			vectors[i][0] = 1
			vectors[i][1] = np.arccos(2*np.random.sample()-1)
			vectors[i][2] = 2*np.pi*np.random.sample()
		
		return vectors
	
	def generate_mean_free_path_polar_vectors(self, num,mean_free_path):
		"""
		Generates a random polar vectors

		Parameters
		----------
		num : int
			number of varibles to be generated
		mean_free_path : float

		Returns
		-------
		2d np array
			Random polar vectors

		"""
		vectors = np.zeros((num,3))
		for i in range(num):
			vectors[i][0] = -mean_free_path*np.log(np.random.sample())
			vectors[i][1] = np.arccos(2*np.random.sample()-1)
			vectors[i][2] = 2*np.pi*np.random.sample()
			
		return vectors
	
	def fictious_step(self,mean_free_path,thicknesses,s,d,v):
		"""
		Completes a fictious step by repeating a the last step directions with
		the smallest material length

		Parameters
		----------
		mean_free_path : float
		thicknesses : 1d np array
		s : int
			A variable showing which material the particle is in
		d : float
			Distance jumped
		v : 1d np array
			A vector

		Returns
		-------
		x : float
			Distance stepeed in

		"""
		x = 0.0
		f = 0
		path_temp = 0.0
		if x > thicknesses[0] and s == 0:
			f = 1
			path_temp = min(mean_free_path[0:2])
		elif x < thicknesses[0] and s == 1:
			f = 1
			path_temp = min(mean_free_path[0:2])
		elif x > np.sum(thicknesses[0:2]) and s == 0:
			f = 1
			path_temp = min(mean_free_path[0:3])
		elif x > np.sum(thicknesses[0:2]) and s == 1:
			f = 1
			path_temp = min(mean_free_path[1:3])
		elif x < np.sum(thicknesses[0:2]) and s == 2:
			f = 1
			path_temp = min(mean_free_path[1:3])
		elif x < thicknesses[0] and s == 2:
			f = 1
			path_temp = min(mean_free_path[0:3])
		if f == 1:
			x -= v[0] * d
			while f == 1:
				d = self.generate_distances(1,path_temp)[0]
				x += d * v[0]
				if s == 0 and x >= thicknesses[0]:
					f =  0
				elif s == 1 and x < thicknesses[0]:
					f =  0
				elif s == 1 and x > np.sum(thicknesses[0:2]):
					f =  0
				elif s == 2 and x < np.sum(thicknesses[0:2]):
					f =  0
				elif x < 0:
					f = 0
				elif x > np.sum(thicknesses):
					f = 0
		return x
	
	def simulate_three_materials(self, thicknesses, num):
		"""
		Tracks multiple path histories through materials.
		If using < 3 materials pass 0 into the extra thicknesses

		Parameters
		----------
		thicknesses : 1d np array size 3
		num : int
			number of neutrons simulated

		Returns
		-------
		list
			DESCRIPTION.

		"""
		#thinkness in cm
		N_a = 6.02214076e23
		
		num_density = np.zeros(3)
		prob_absorbtion = np.zeros(3)
		mean_free_path = np.zeros(3)
		
		for i in range(self.material_count):
			num_density[i] = self.density[i]*N_a/self.mass[i]
			prob_absorbtion[i] = self.absorbtion_area[i]/(
				self.scattering_area[i]+self.absorbtion_area[i])
			mean_free_path[i] = 1/(
				num_density[i]*(self.scattering_area[i]
					  +self.absorbtion_area[i])*10**-24)
		
		relf = 0
		trans = 0
		absor = 0
		
		for i in range(num):
			importance = 1
			s = 0
			d = self.generate_distances(1,mean_free_path[0])[0]
			v = np.array([1,0,0])
			x = d
			x += self.fictious_step(mean_free_path, thicknesses, s, d, v)
			
			prob_temp = prob_absorbtion[0]
			path_temp = mean_free_path[0]
			
			while importance == 1:
				if x < 0:
					relf += 1
					importance = 0
					break
				elif x > np.sum(thicknesses):
					trans += 1
					importance = 0
					break
				if x < thicknesses[0]:
					s = 0
					prob_temp = prob_absorbtion[0]
					path_temp = mean_free_path[0]
				elif x >= thicknesses[0] and x <= np.sum(thicknesses[0:2]):
					s = 1
					prob_temp = prob_absorbtion[1]
					path_temp = mean_free_path[1]
				elif x > np.sum(thicknesses[0:2]):
					s = 2
					prob_temp = prob_absorbtion[2]
					path_temp = mean_free_path[2]
				if np.random.sample() > prob_temp:
					v_ = self.generate_polar_vectors(1)
					d = self.generate_distances(1,path_temp)[0]
					v__ = self.get_cartesians(v_)[0]
					x += v__[0] * d
					x += self.fictious_step(mean_free_path, thicknesses, s, d, v__)
				else:
					absor += 1
					importance = 0
					break
					
		return [trans/num, relf/num, absor/num]
	
	def full_error_simulation(self,num_of_exp,thickness,num_particules):
		"""
		Conduct mutiple simulations to calculte error

		Parameters
		----------
		num_of_exp : int
			Number of simulations conducted
		thickness : 1d np array size 3
		num_particules : int
			Number of particules simulated

		Returns
		-------
		list
			Results of experients

		"""
		
		t = np.zeros(num_of_exp)
		r = np.zeros(num_of_exp)
		a = np.zeros(num_of_exp)
		
		for i in range(num_of_exp):
			arr = self.simulate_three_materials(thickness,num_particules)
			t[i] = arr[0]
			r[i] = arr[1]
			a[i] = arr[2]
		
		mean_t = np.sum(t)/num_of_exp
		mean_r = np.sum(r)/num_of_exp
		mean_a = np.sum(a)/num_of_exp
		
		std_t = (np.sum((mean_t-t)**2)/(num_of_exp-1))**0.5
		std_r = (np.sum((mean_r-r)**2)/(num_of_exp-1))**0.5
		std_a = (np.sum((mean_a-a)**2)/(num_of_exp-1))**0.5
		
		return [[mean_t,std_t],[mean_r,std_r],[mean_a,std_a]]
	
	def visual_sample_three_mateirals(self,thicknesses):
		#thinkness in cm
		N_a = 6.02214076e23
		
		collisions = np.zeros((2**16,3))
		
		num_density = np.zeros(3)
		prob_absorbtion = np.zeros(3)
		mean_free_path = np.zeros(3)
		
		for i in range(self.material_count):
			num_density[i] = self.density[i]*N_a/self.mass[i]
			prob_absorbtion[i] = self.absorbtion_area[i]/(
				self.scattering_area[i]+self.absorbtion_area[i])
			mean_free_path[i] = 1/(
				num_density[i]*(self.scattering_area[i]
					  +self.absorbtion_area[i])*10**-24)
		
		relf = 0
		trans = 0
		absor = 0
		count = 0
		collisions[count] = np.array([0,0,0])
		count += 1
		
		importance = 1
		s = 0
		d = self.generate_distances(1,mean_free_path[0])[0]
		v = np.array([1,0,0])
		x = d
		x += self.fictious_step(mean_free_path, thicknesses, s, d, v)
		collisions[count] = np.array([x,0,0])
		
		prob_temp = prob_absorbtion[0]
		path_temp = mean_free_path[0]
		
		while importance == 1:
			if x < 0:
				relf += 1
				importance = 0
				break
			elif x > np.sum(thicknesses):
				trans += 1
				importance = 0
				break
			if x < thicknesses[0]:
				s = 0
				prob_temp = prob_absorbtion[0]
				path_temp = mean_free_path[0]
			elif x >= thicknesses[0] and x <= np.sum(thicknesses[0:2]):
				s = 1
				prob_temp = prob_absorbtion[1]
				path_temp = mean_free_path[1]
			elif x > np.sum(thicknesses[0:2]):
				s = 2
				prob_temp = prob_absorbtion[2]
				path_temp = mean_free_path[2]
			if np.random.sample() > prob_temp:
				v_ = self.generate_polar_vectors(1)
				d = self.generate_distances(1,path_temp)[0]
				v__ = self.get_cartesians(v_)[0]
				x += v__[0] * d
				x += self.fictious_step(mean_free_path, thicknesses, s, d, v__)
				collisions[count] = [x,d*v__[1],d*v__[2]]
				count += 1
			else:
				absor += 1
				importance = 0
				break
					
		return collisions[0:count]

################### Main Scripting area #####################################

exp_water = experiment_neutrinoes(0.6652,103.0,1.00,18.01528)
exp_lead = experiment_neutrinoes(0.158,11.221,11.35,207.2)
exp_graphite = experiment_neutrinoes(0.0045,4.74,1.67,12.011)

exp_water_graphite = experiment_neutrinoes(0.6652,103.0,1.00,18.01528)
exp_water_graphite.define_another_material(0.0045,4.74,1.67,12.011)
exp_water_graphite.define_another_material(0.6652,103.0,1.00,18.01528)

results_more = exp_water_graphite.full_error_simulation(10,np.array([1,2,1]),10000)
print(results_more)

## Draws 3d grapth of particule simulations
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

y = np.linspace(-4,4,10)
z = np.linspace(-4,4,10)

Y,Z = np.meshgrid(y,z)
X= 0*Y + 0*Z + 0
ax.plot_surface(X,Y,Z, alpha=0.2)
ax.plot_surface(X+1,Y,Z, alpha=0.2)
ax.plot_surface(X+3,Y,Z, alpha=0.2)
ax.plot_surface(X+4,Y,Z, alpha=0.2)

for i in range(30):
	v = exp_water_graphite.visual_sample_three_mateirals(np.array([1,2,1]))
	v = v.transpose()
	ax.plot(v[0],v[1],v[2])

plt.show()

def remove_zeroes(tranmissions,tranmissions_errors,thicknesses):
	"""
	Removes data that can't be used, zeroes and infinits

	Parameters
	----------
	tranmissions : 1d np array
	tranmissions_errors :  1d np array
	thicknesses :  1d np array

	Returns
	-------
	tranmissions :  1d np array
	tranmissions_errors :  1d np array
	thicknesses :  1d np array
	log_errors :  1d np array
	log_weights :  1d np array

	"""
	#removes for linear
	index_non_zero = np.transpose(np.argwhere(tranmissions!=0))
	tranmissions = tranmissions[index_non_zero][0]
	thicknesses = thicknesses[index_non_zero][0]
	tranmissions_errors = tranmissions_errors[index_non_zero][0]
	
	index_non_zero = np.transpose(np.argwhere(tranmissions_errors!=0))
	tranmissions = tranmissions[index_non_zero][0]
	thicknesses = thicknesses[index_non_zero][0]
	tranmissions_errors = tranmissions_errors[index_non_zero][0]
	
	log_errors = tranmissions_errors/tranmissions
	log_weights = 1/log_errors
	index_non_inf = np.transpose(np.argwhere(log_weights!=np.inf))
	log_weights = log_weights[index_non_inf][0]
	thicknesses = thicknesses[index_non_inf][0]
	tranmissions = tranmissions[index_non_inf][0]
	log_errors = log_errors[index_non_inf][0]
	
	return tranmissions, tranmissions_errors,thicknesses, log_errors, log_weights

## Calculates attenuation lengths
#sets up arrays
thicknesses = np.linspace(1,40,num=100)
thicknesses_init = thicknesses
tranmissions_water = np.array([])
tranmissions_water_errors = np.array([])
tranmissions_lead = np.array([])
tranmissions_lead_errors = np.array([])
tranmissions_graphite = np.array([])
tranmissions_graphite_errors = np.array([])

#simulation data
no_neutrons = 10000
no_sims = 10

#produces all simulated data
for x in thicknesses:
	temp = exp_water.full_error_simulation(no_sims,[x,0,0],no_neutrons)[0]
	tranmissions_water = np.append(tranmissions_water,temp[0])
	tranmissions_water_errors = np.append(tranmissions_water_errors,temp[1])
	
	temp = exp_lead.full_error_simulation(no_sims,[x,0,0],no_neutrons)[0]
	tranmissions_lead = np.append(tranmissions_lead,temp[0])
	tranmissions_lead_errors = np.append(tranmissions_lead_errors,temp[1])
	
	temp = exp_graphite.full_error_simulation(no_sims,[x,0,0],no_neutrons)[0]
	tranmissions_graphite = np.append(tranmissions_graphite,temp[0])
	tranmissions_graphite_errors = np.append(tranmissions_graphite_errors,temp[1])

plt.figure()

#sets up loopable object
transmissions = [tranmissions_water,tranmissions_lead,tranmissions_graphite]
transmission_errors = [tranmissions_water_errors,
					   tranmissions_lead_errors,tranmissions_graphite_errors]
colours = ['blue','red','green']
names = ['water','lead','graphite']
#fits and plots
for i in range(3):
	thicknesses = thicknesses_init
	transmissions[i], transmission_errors[i],thicknesses, log_errors,log_weights = remove_zeroes(
				transmissions[i], transmission_errors[i],thicknesses)
	[m,c],cov = np.polyfit(thicknesses,-np.log(transmissions[i]),1,w=log_weights,cov=True)
	print("Charatcheristic legnth")
	print(1/m)
	print("Std")
	print(np.sqrt(cov[0][0]/transmissions[i].size)/m**2)
	print(calculate_reducde_chi(-np.log(transmissions[i]),m*thicknesses+c,log_errors))
	plt.scatter(thicknesses,-np.log(transmissions[i]),s=2,c=colours[i])
	plt.plot(thicknesses, m*thicknesses+c,c=colours[i],label=names[i])
	plt.errorbar(thicknesses,-np.log(transmissions[i]),yerr=log_errors,fmt='none',c=colours[i])

plt.title("Log of Survivability of neutrons through three materials")
plt.xlabel("Thickness (cm)")
plt.ylabel("-ln(Transmission)")
plt.legend()
plt.show()
