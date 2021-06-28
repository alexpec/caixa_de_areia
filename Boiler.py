import warnings
import numpy


from scipy.integrate import odeint
from scipy.integrate import solve_ivp


start_paster = False
start_fermentation = False

paster_interval = 0.5*3600.0
fermentation_interval = 8*3600.0



def SetPoint(t):
	'''
		Function that returns the setpoint based in time
		
		@param: t {float} - The time in seconds, since the start of the process
	'''

	Sp0 = 173.0 + 15;

	if t <= 0.0: return Sp0
	elif t>0  and t<= 2000+1800: return 241.0
	elif t>3800 and t <= 25400+21200: return 218.0
	else: return Sp0
	#return 241.0

def Error(t, T_measured):
	'''
		Function that returns the deviation from the setpoin on the
		measured temperature

		@param: t {float} - The time in seconds, since the start of process
		@param: T_measured {float} - The measured temperature in the reactor
	'''
	#import pdb; pdb.set_trace()
	val = SetPoint(t) - T_measured

	return val



def Control(t, T_measured, *args):
	'''
		Control element that returns a float value from 0 - 1 where 0 is the inactivated actuator
		and 1 is the actuator full activated

		@param: t {float} - The time in seconds, since the start of process
		@param: T_measured - The neasured temperature in the reactor
		@param: *args - array{float} - The Kp, Ki and Kd when using PID control; Null to bang-bang
	'''
	error = Error(t, T_measured)

	Sp = SetPoint(t)

	normalized_error = error / Sp

	#if normalized_error > 1.0: 
		#warnings.warn('The actuator reached its limit at %.4f s with an normalized error of %.4f' %(t, normalized_error))
	#	normalized_error = 1.0
	#elif normalized_error < 0:
		#warnings.warn('The actuator reached its limit at %.4f s with an normalized error of %.4f' %(t, normalized_error))
	#	normalized_error = 0.0

	return normalized_error


def Generator(error, GeneratorMax):
	if error > 0.001: return GeneratorMax
	else: return 0.0
	#return GeneratorMax


def TA(t):
	return 173+15 #[K]


def Ei_j(k,L,rho,Cp,V,A,Ti,Tj):
	#Heat from media i to media j
	val = 0.0
	val = A*(k/L)*(1/(rho*V*Cp))*(Ti-Tj)
	
	return val
	


def Process(t, y, 
	    Eg0, 
	    Jk, Ja, Jt, Cd, Jv, Chc,
	    Rk, Ra, Rt, Rd, Rv, Rhc):

	#import pdb; pdb.set_trace()

	Tj, Tr = y
	e = Error(t, Tr)


	Eg = Generator(e, Eg0)
	Ta = TA(t)

	A_tampa = 3.14 * 0.30**2/4

	dydt = [
		Ei_j(Rk,Rt,Cd,Chc,Jv,Ra,Tr,Tj) + (1/(Cd*Jv*Chc))*Eg + Ei_j(Jk,Jt,Cd,Chc,Jv,Ja,Ta,Tj),
		Ei_j(Rk,Rt,Rd,Rhc,Rv,Ra,Tj,Tr)
	]


	
	#print(t, dydt, y)

	return dydt



if __name__ == '__main__':
	jacket_volume 			= (3.14 * 0.34**2 / 4)*0.30 - (3.14 * 0.28**2 / 4)*0.25 		#[m3] - 0.01184
	jacket_thikness		= 0.005							  	#[m]
	jacket_thermal_conductivity	= 0.04						  			#[W/m.K] la de rocha!!
	jacket_transfer_area		= 3.14*0.30*0.34 + 3.14*0.34**2/4			  		#m2 - 1.1765

	reactor_volume 			= (3.14 * 0.28**2 / 4)*0.25				  		#[m3] - 0.01539
	reactor_thikness		= 0.050						  		#[m]
	reactor_thermal_conductivity	= 15.0						  			#[W/m.K]
	reactor_transfer_area		= 3.14*0.25*0.28 + 3.14*0.28**2/4 		  			#m2 - 0.9407

	energy_generator		= 3000					  			#W

	coolant_density		= 1000.0			 			  		#[kg/m3]
	coolant_heat_capacity		= 4200.0		 			  			#[J/kg.K]

	reactant_density		= 1032.0		 			  			#[kg/m3]
	reactant_heat_capacity		= 4000.0		 			  			#[J/kg.K]



	Jv  = jacket_volume
	Jt  = jacket_thikness
	Jk  = jacket_thermal_conductivity
	Ja  = jacket_transfer_area

	Rv  = reactor_volume
	Rt  = reactor_thikness
	Rk  = reactor_thermal_conductivity
	Ra  = reactor_transfer_area

	Cd  = coolant_density
	Chc = coolant_heat_capacity

	Rd  = reactant_density
	Rhc = reactant_heat_capacity 


	t0 = 0.0
	tf = 100000.0
	t = numpy.linspace(t0, tf, 1000)

	y0 = numpy.array([173+4, TA(0)])


	sol = solve_ivp(Process, [t0,tf], y0, t_eval=t, args=( 
		  	energy_generator,
		  	Jk, Ja, Jt, Cd, Jv, Chc,
			Rk, Ra, Rt, Rd, Rv, Rhc
		 		        )
		)


	e_t = numpy.zeros(len(t))
	u_t = numpy.zeros(len(t))
	Eg_t = numpy.zeros(len(t))
	Sp = numpy.zeros(len(t))
	
	y = numpy.transpose(sol.y)
	for idx, i in enumerate(t):
		Tr, Tj = y[idx]
		u_t[idx] = Control(sol.t[idx], Tr)
		Eg_t[idx] = Generator(u_t[idx], energy_generator)
		Sp[idx] = SetPoint(sol.t[idx]) 
		e_t[idx] = Error(sol.t[idx], Tr)

	import matplotlib.pyplot as plt
	total_plots = 3
	fig, ax = plt.subplots(total_plots, sharex=True, sharey=False)
	ax[0].plot(t, y[:,0], t, Sp)
	ax[1].plot(t, y[:,1])
	ax[2].plot(t, Eg_t, label='Eg')
	#ax[3].plot(t, Sp, label='SP')
	#ax[4].plot(t, u_t, label='u(t)')
	#ax[5].plot(t, e_t, label='e(t)')

	plt.show()

	#import pdb; pdb.set_trace()
