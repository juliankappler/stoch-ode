# distutils: language = c++
# distutils: extra_compile_args = -std=c++11

import numpy as np
cimport numpy as np
import os
import time
cimport cpython
import types
from libcpp cimport bool
DTYPE   = np.float

cdef extern from "math.h":
	double sqrt(double x) nogil

cdef extern from "math.h":
	double exp(double x) nogil

from libcpp.vector cimport vector

# wrapper for C++11 pseudo-random number generator
cdef extern from "<random>" namespace "std":
	cdef cppclass mt19937:
		mt19937()
		mt19937(unsigned long seed)

	cdef cppclass normal_distribution[T]:
		normal_distribution()
		normal_distribution(T a, T b)
		T operator()(mt19937 gen)




cdef class integrator_base:

	cdef:
		np.ndarray trajectory
		int N_dim, N_steps
		double dt, sqrt_dt
		double t, x, x0
		double a, b, b_bp
		mt19937 gen
		long seed
		bool b_constant, verbose
		bool milstein

	def __init__(self, parameters):
		cdef:
			int i

		self.N_dim      = 1
		self.dt         = parameters['dt']
		self.N_steps    = parameters['N_steps']
		self.sqrt_dt    = sqrt(self.dt)

		self.verbose = False

		self.x    = 0
		self.x0   = 0
		self.a    = 0
		self.b    = 0
		self.b_bp = 0
		self.trajectory = np.zeros( [self.N_steps+1], dtype=DTYPE)

		# set seed for pseudo-random number generator (if provided)
		try:
			self.initialize_random_number_generator(
			                      supplied_seed=parameters['seed'])
		except KeyError:
			self.initialize_random_number_generator()

	cdef random_standard_normal(self):
		cdef:
			normal_distribution[double] dist = normal_distribution[double](0.0,1.0)
		return dist(self.gen)

	cpdef initialize_random_number_generator(self,long supplied_seed=-1):
		cdef:
			long max_long = 9223372036854775807
		if supplied_seed < 0:
			self.seed = (abs(os.getpid()) + long(time.time()*1000)) % max_long
		else:
			self.seed = supplied_seed % max_long
		self.gen = mt19937(self.seed)

	'''
	cdef set_x0(self,x0):
		cdef:
			double [:] trajectory  = self.trajectory
			double x = self.x
			int N_dim = self.N_dim

		# write initial condition to trajectory
		trajectory[0] = x0
		x = x0
		return 0
	''';

	'''
	cdef update_a(self,a_function):
		cdef:
			double x = self.x
			double a = self.a
			int N_dim = self.N_dim
			np.ndarray a_eval

		a = a_function(x)
		return 0

	cdef update_b(self,b_function):
		cdef:
			double [:] x = self.x
			double t = self.t
			double [:,:] b = self.b
			int N_dim = self.N_dim
			np.ndarray b_eval

		b_eval = b_function(np.array(x),t)
		b[i,j] = b_eval[i]
		return 0
	''';


	cdef return_output_dictionary(self):
		output_dictionary = {'N_dim':self.N_dim,
                            'N_steps':self.N_steps,
                            'dt':self.dt,
                            'x':self.trajectory,
                            't':np.arange(self.N_steps+1)*self.dt,
                            }
        #
		return output_dictionary



cdef class integrator(integrator_base):

	'''
	cdef integration_step(self):
		cdef:
			double b = self.b
			double b_bp = self.b_bp
			double a = self.a
			double x = self.x
			double dt = self.dt
			double sqrt_dt = self.sqrt_dt
			double cur_normal
			bool milstein = self.milstein
			int i, j
		#
		x += a*dt
		#
		cur_normal = self.random_standard_normal()
		x += b * cur_normal * sqrt_dt
		if milstein:
			x += b_bp * ( cur_normal**2 - 1. ) * dt/2.
		#
		return 0
	''';

	cpdef simulate(self, x0, a_func, b_func,
						b_bp_func = None, # b*db/dx .. to be used for milstein method
	            ):
		cdef:
			int N_dim = self.N_dim
			int N_steps = self.N_steps
			double [:] trajectory = self.trajectory
			double x = self.x
			double a = self.a
			double b = self.b
			double t
			bool milstein = self.milstein
			int i, j

		#self.set_x0(x0)
		x = x0
		trajectory[0] = x0

		if b_bp_func is not None:
			milstein = True

		for i in range(N_steps):
			t = i*self.dt
			# update force and diffusivity matrix
			#self.update_a(a_function=a)
			a = a_func(x,t)
			b = b_func(x,t)
			if milstein:
				b_bp = b_bp_func(x,t)
			#
			# perform integration step
			#self.integration_step()
			x += a*self.dt
			#
			cur_normal = self.random_standard_normal()
			x += b * cur_normal * self.sqrt_dt
			if milstein:
				x += b_bp * ( cur_normal**2 - 1. ) * self.dt/2.

			# update output trajectory array
			trajectory[i+1] = x

		return self.return_output_dictionary()
