import warnings
import numpy
import math


from scipy.integrate import odeint
from scipy.integrate import solve_ivp




set_point = 0.0

def Error(t, T_measured):
	'''
		Function that returns the deviation from the setpoin on the
		measured temperature

		@param: t {float} - The time in seconds, since the start of process
		@param: T_measured {float} - The measured temperature in the reactor
	'''
	#import pdb; pdb.set_trace()
	val = set_point - T_measured

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

	Sp = set_point

	normalized_error = error / Sp

	#if normalized_error > 1.0: 
		#warnings.warn('The actuator reached its limit at %.4f s with an normalized error of %.4f' %(t, normalized_error))
	#	normalized_error = 1.0
	#elif normalized_error < 0:
		#warnings.warn('The actuator reached its limit at %.4f s with an normalized error of %.4f' %(t, normalized_error))
	#	normalized_error = 0.0

	return normalized_error


def Generator(error, GeneratorMax):
	if error > 0.1: return GeneratorMax
	else: return 0.0
	#return GeneratorMax

def MassInflux(error, h, MaxReactorHeight):
	if error < -1.0 and h < MaxReactorHeight: return True
	elif h >= MaxReactorHeight: 
		warnings.warn('Max Jacket Height reached! Waiting until coolant flow')
		return False 
	elif h < 0.0:
		warnings.warn('How does h reached lesser than zero!?')
	else: return False


def TA(t):
	return 173+15 #[K]




def Ei_j(k,L,rho,Cp,V,A,Ti,Tj):
	#Heat from media i to media j
	val = 0.0
	val = A*(k/L)*(1/(rho*V*Cp))*(Ti-Tj)

	return val


def m_flux(A, h):
	if h < 0.0: return 0.0
	m = A * math.sqrt(2.0 * 9.8 * h ) 
	return m

def Temperature_in():
	return 173.15 + 15


def Process(t, y, 
	    Eg0, 
	    Jk, Ja, Jt, Cd, Jv, Chc,
	    Rk, Ra, Rt, Rd, Rv, Rhc, Mjh):

	#import pdb; pdb.set_trace()


	Tj, Tr, h = y
	e = Error(t, Tr)
	Eg = Generator(e, Eg0)
	Ta = TA(t)
	Tin = Temperature_in()
	H = 0.03 #[m]


	A_tampa = 3.14 * 0.30**2/4 #[m2]
	Ain = 3.14 * 0.0127**2/4 #[m2]
	Aout = 3.14 * 0.0127**2/4 #[m2] 
	
	if MassInflux(e, h, Mjh): m_in = m_flux(Ain,H)
	else: m_in = 0.0
	

	
	dydt = [
		Ei_j(Rk,Rt,Cd,Chc,Jv,Ra,Tr,Tj) + (1/(Cd*Jv*Chc))*Eg + Ei_j(Jk,Jt,Cd,Chc,Jv,Ja,Ta,Tj) + (1/(Cd*Jv))*(Tin*m_in-Tj*m_flux(Aout, h)),
		Ei_j(Rk,Rt,Rd,Rhc,Rv,Ra,Tj,Tr),
		(1.0 / (0.34**2 - 0.28**2))*(m_in - m_flux(Aout,h))
	]


	
	#print(t, dydt, y)

	return dydt



if __name__ == '__main__':
	jacket_volume 			= (3.14 * 0.34**2 / 4)*0.30 - (3.14 * 0.28**2 / 4)*0.25 		#[m3] - 0.01184
	jacket_thikness			= 0.005								  	#[m]
	jacket_thermal_conductivity	= 0.04						  			#[W/m.K] la de rocha!!
	jacket_transfer_area		= 3.14*0.30*0.34 + 3.14*0.34**2/4			  		#m2 - 1.1765

	reactor_volume 		= (3.14 * 0.28**2 / 4)*0.25				  		#[m3] - 0.01539
	reactor_thikness		= 0.003						  		#[m]
	reactor_thermal_conductivity	= 15.0						  			#[W/m.K]
	reactor_transfer_area		= 3.14*0.25*0.28 + 3.14*0.28**2/4 		  			#m2 - 0.9407

	energy_generator		= 3000					  				#W

	coolant_density		= 1000.0			 			  		#[kg/m3]
	coolant_heat_capacity		= 4200.0		 			  			#[J/kg.K]

	reactant_density		= 1032.0		 			  			#[kg/m3]
	reactant_heat_capacity		= 4000.0		 			  			#[J/kg.K]

	max_jacket_height		= 5.0 *10**-2								#[m]

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

	Mjh = max_jacket_height
	
	
	time_start_process = 0.0
	time_start_pasteur = 0.0
	time_end_paster = 0.0
	time_start_cooling = 0.0
	time_start_fermentation = 0.0
	time_end_fermentation = 0.0
	time_end_final_cooling = 0.0


	pasteur_time = 30*60 #[s]
	fermentation_time = 8*3600 #[s]

	heating_setpoint = 173.15 + 68 #K
	fermentation_setpoint = 173.15 + 45 #K
	end_cooling_setpoint = TA(0.0) #K
	

	def evt_end_heating(t, y, *args):
		Tj, Tr, h = y		
		return (heating_setpoint - Tr) / heating_setpoint
	evt_end_heating.terminal = True

	def evt_end_st_cooling(t, y, *args):
		Tj, Tr, h = y
		return (fermentation_setpoint - Tr) / fermentation_setpoint
	evt_end_st_cooling.terminal = True

	def evt_end_cooling(t, y, *args):
		Tj, Tr, h = y		
		return (end_cooling_setpoint - Tr) / end_cooling_setpoint
	evt_end_cooling.terminal = True


			  #Step               #duration time [s], #SP, #evt
	process_steps = {'Initital Heating': 	[None, heating_setpoint, evt_end_heating],
			  'Pasteur': 		[30*60, heating_setpoint, None],
			  'First Cooling': 	[None, fermentation_setpoint, evt_end_st_cooling],
			  'Fermenting': 	[8.0*3600.0, fermentation_setpoint, None],
			  'Second Cooling': 	[None, end_cooling_setpoint, evt_end_cooling]
			 }
			 
	processes_solvers = {'Initital Heating': None,
			  'Pasteur': None,
			  'First Cooling': None,
			  'Fermenting': None,
			  'Second Cooling': None
			 }

	processes_times = {'Initital Heating': 0.0,
			  'Pasteur': 0.0,
			  'First Cooling': 0.0,
			  'Fermenting': 0.0,
			  'Second Cooling': 0.0
			 }

	
	last_step = None
	
	for key, value in process_steps.items():		
					
		if key == 'Initital Heating':			
			y0 = numpy.array([173+4, TA(0), 0.0])
			
		else:
			y0 = processes_solvers[last_step].y.T[-1]
		evt = process_steps[key][2]
		duration = process_steps[key][0]

		if evt == None: evt = []
		elif duration == None: duration = 10000

		t0 = 0.0
		tf = duration
		t = numpy.linspace(t0, tf, 100)

		set_point = value[1]	


		processes_solvers[key] = solve_ivp(Process, [t0,tf], y0, t_eval=t, args=( 
			  	energy_generator,
			  	Jk, Ja, Jt, Cd, Jv, Chc,
				Rk, Ra, Rt, Rd, Rv, Rhc,
				Mjh
			 		        ), 
				events=process_steps[key][2],
			)
		
		last_step = key
		processes_times[key] = processes_solvers[key].t[-1]

	time = numpy.array([0.0])
	Tr   = numpy.array([0.0])
	Tj   = numpy.array([0.0])
	h    = numpy.array([0.0])
	Sp   = numpy.array([0.0])
	points_of_interest = []

	for key, item in processes_solvers.items():
		
		points = len(item.y.T[:,0])
		last_t = time[-1]
		points_of_interest.append(last_t)
		time = numpy.concatenate([time, numpy.linspace(last_t, last_t + processes_times[key], points)])
		Tj   = numpy.concatenate([Tj, item.y.T[:,0]])
		Tr   = numpy.concatenate([Tr, item.y.T[:,1]])
		h    = numpy.concatenate([h, item.y.T[:,2]])
		Sp   = numpy.concatenate([Sp, process_steps[key][1]*numpy.ones(points)])

	#import pdb; pdb.set_trace()
	
	time = (1/3600)*time 

	import matplotlib.pyplot as plt	
	fig, ax = plt.subplots(3, sharex=False, sharey=False)
	ax[0].plot(time, Tr, time, Sp)
	for i in points_of_interest:
		ax[0].axvline(x=i/3600.0, color='r', linestyle='--')
	ax[1].plot(time, Tj)
	ax[2].plot(time, h)
	
	plt.show()

		

	'''
	e_t = numpy.zeros(len(t))
	u_t = numpy.zeros(len(t))
	Eg_t = numpy.zeros(len(t))
	Sp = numpy.zeros(len(t))
	
	y = numpy.transpose(sol.y)
	for idx, i in enumerate(t):
		Tr, Tj, h = y[idx]
		u_t[idx] = Control(sol.t[idx], Tr)
		Eg_t[idx] = Generator(u_t[idx], energy_generator)
		Sp[idx] = SetPoint(sol.t[idx]) 
		e_t[idx] = Error(sol.t[idx], Tr)

	import matplotlib.pyplot as plt
	total_plots = 4
	fig, ax = plt.subplots(total_plots, sharex=True, sharey=False)
	ax[0].plot(t, y[:,0], t, Sp)
	ax[1].plot(t, y[:,1])
	ax[2].plot(t, y[:,2])
	ax[3].plot(t, Eg_t, label='Eg')
	#ax[3].plot(t, Sp, label='SP')
	#ax[4].plot(t, u_t, label='u(t)')
	#ax[5].plot(t, e_t, label='e(t)')

	plt.show()
	'''
	#import pdb; pdb.set_trace()
