"""
Created on Fri May 12 12:39:38 2023

@author: l.anderson
"""
#scripting the FlexPDE file in python
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LogNorm
from scipy.ndimage.filters import gaussian_filter
from scipy.optimize import curve_fit
import math
from mpl_toolkits.mplot3d import axes3d

plt.rcParams['font.family'] = "times"

#Flex file goes in quotations below
flexcode = """
TITLE '1D autochemotaxis for a heat source particle in a harmonic potential, Dirichlet b.c.s'     { the problem identification }

SELECT 
penwidth = 7

COORDINATES cartesian1  { coordinate system, 1D,2D,3D, etc }
VARIABLES        	    { system variables }
  u(threshold=100)      {GRADIENT OF the temperature - threshold clause is needed to quantify error estimates}
  	
GLOBAL VARIABLES	 {note : if global variables are not used for rp,vp the system becomes unstable}
  r_p(threshold=100)    {position of the heat source, i.e. the particle, which moves}

  v_p(threshold=100) {velocity of particle}

DEFINITIONS     { parameter definitions }
kappa = %s		{spring constant of the harmonic potential}
D = 0.1	{Diffusion coefficient (thermal conductivity )}
Ta = 25			{ambient temperature}
w = 0.3			{gaussian source radius}
A = 1/sqrt(4*pi*w^2)			{amplitude of gaussian heat source}
source = A*exp(-(x-r_p)^2/(4*w^2))	{defining a gaussian heat source}
grad_source = (A*(r_p - x)*exp(-(x-r_p)^2/(4*w^2)))/(2*w^2)

l = %s		{soret parameter (strength of the chemotactic coupling), l>0 for negative chemotaxis/"chemorepulsion"}
					 {note: if we set l<0 for chemoattraction the solutions are stable but uninteresting, the particle just sticks to the wall}
grad_p = eval(XCOMP(grad(u)),r_p)		{XCOMP Returns a scalar whose value is the first component of the vector argument, we must do this here to avoid issues employing the grad operator in one dimension}
gamma = 1 	{Stokes' drag coefficient - 6*pi*viscosity*radius - typical value for air/water in later simlns?}                    
M = 0			{inertia term neglected in the overdamped regime}
!v_b = l/(2*D)		{steady state velocity predicted for a point-like heat source (true delta function) }
h = 10      	{half width of the box}                    
mu = 0.1		{heat absorption rate}
Ctime = %s

INITIAL VALUES
u = -1/(2*D) * TANH(r_p)*EXP(-SQRT(mu/D)*ABS(r_p))
r_p = 0.01	{release source from initial position}
v_p = 0 {..with initial velocity} 

EQUATIONS        		{SINCE U=GRAD(TEMP) MUST USE DIRICHLET B.CS }						
 u:  D*div(grad(u)) + grad_source - dt(u) - mu*u =0 			{the PDE with source term} 
 v_p:  M*dt(v_p) + gamma*v_p = -l*eval(u,r_p) - kappa*r_p				{describes the effective chemotactic coupling force, assuming an overdamped particle with stokes drag coeff gamma}
 r_p: dt(r_p) = v_p

BOUNDARIES       { The domain definition }

  REGION 1       { For each material region }
    START(-h)   { Walk the domain boundary }
    POINT VALUE(u) = 0	 LINE TO (+h)	POINT VALUE(u) = 0	{applied Dirichlet boundary condition at either end of the box}
 !Note: no close statement for 1D since the boundary path defines the domain
    
TIME 0 TO Ctime by 0.2  { if time dependent }

MONITORS         { shows progress of the simulations while running }
PLOTS            
	FOR   T = 0 BY 0.2 TO Ctime 
        !History( r_p ) 	
        History(v_p , r_p) PrintOnly Export Format "#t#b#1#b#2" file="harmonicdata.txt"
        !ELEVATION(u) FROM (-h) TO (+h)		{ plotting an elevation of temperature along the line }
     
END

"""
"""
FlexFileName = "flexharmonicscript.pde"
import subprocess

k_values = np.arange(0,2,0.05)
l_values = np.arange(0,1,0.05)
phase_data = np.zeros(shape=(np.size(k_values),np.size(l_values)))
frequency_data = np.zeros(shape=(np.size(k_values),np.size(l_values)))
print(phase_data)

def msd_1d(x):
    #Calculates the root mean square amplitude for a 1 dimensional trajectory
    dx_sq = np.square(x)
    msd = np.mean(dx_sq)
    return np.sqrt(msd)


SAoutputdat = np.loadtxt('stabilityphasediagram.txt')
rpSA = SAoutputdat[:,0]
vpSA = SAoutputdat[:,1]

#for k in (k_values):
 #   for i in l_values:    
#***comment out below***       
k = 0.2
i = 0.6
with open(FlexFileName, "w") as f:
    print(flexcode%(k,i,"400"),file=f) 
subprocess.run(["C:\Program Files\FlexPDE7\FlexPDE7.exe", "-S", FlexFileName])
        
flexoutputrawdata = np.loadtxt('flexharmonicscript_output\\harmonicdata.txt', skiprows=8)        
t = flexoutputrawdata[1000:1500,0] #[50:,0]
r_p = flexoutputrawdata[1000:1500,2]
integer_k = int(k*20)
integer_l = int(i*20)
r_data = flexoutputrawdata[:800,2] 
t_data = flexoutputrawdata[:800,0]
v_p = flexoutputrawdata[1000:1500,1]

***comment out above***
#array[i,j] picks out the j+1 th column of the i+1 th row 
phase_data[integer_k,integer_l] = msd_1d(r_p)
print(k)
num_samples = len(t)
samplingrate = 10
# apply fast fourier transform and take absolute values
freq_spectrum=np.abs(np.fft.fft(r_p))/num_samples  
# get the list of frequencies
frequencies = np.fft.fftfreq(num_samples, 1/samplingrate)
       
#find the element of the spectrum array with max value        
a = np.unravel_index(freq_spectrum.argmax(), freq_spectrum.shape)
#extract corresponding frequency - this is the fundamental(?)/dominant frequency
frequency_data[integer_k,integer_l] = frequencies[a]
        
#print(frequency_data)
       
#INDENT EVERYTHING ABOVE THIS
    
np.savetxt('frequencymaprawdata.txt', frequency_data)
np.savetxt('heatmaprawdata.txt', phase_data)

UNCO??ENT BELOZ UNTIL 192
phase_data = np.loadtxt('heatmaprawdata.txt', skiprows=0, usecols = (0,1,2,3,4,5,6,7,8,9,10))
test = phase_data[:20,]
frequency_data = np.loadtxt('frequencymaprawdata.txt', skiprows= 0, usecols =(0,1,2,3,4,5,6,7,8,9,10))
test2 = 2*np.pi*frequency_data[:20,]+0.000001
log_norm = LogNorm(vmin=test.min().min(), vmax=test.max().max())
cbar_ticks = [math.pow(10, i) for i in range(math.floor(math.log10(test.min().min())), 1+math.ceil(math.log10(test.max().max())))]
amp_smooth = gaussian_filter(test, sigma = 0.5)
freq_smooth = gaussian_filter(test2, sigma= 0.5)
#plt.title('Heatmap in {$\kappa$, $\lambda$ } parameter space: \n Harmonic system, D = 0.1, $\mu$ = 0.1, $\sigma$ = 0.3')
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(11,6))
plt.subplot(121)
sns.heatmap(amp_smooth,norm=log_norm, vmin = np.min(amp_smooth[1:20,1:10]), vmax = np.max(amp_smooth[1:20,1:10]),
   cbar_kws={"ticks": cbar_ticks})

#sns.heatmap(amp_smooth,vmin = np.min(amp_smooth[1:20,1:10]), vmax = np.max(amp_smooth[1:20,1:10]))

#plt.yticks([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39],[0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95])
plt.yticks([0,2,4,6,8,10,12,14,16,18,20],[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
plt.xticks([0,1,2,3,4,5,6,7,8,9,10],[0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5] )
plt.title('Particle rms amplitude', loc='center', style='italic')
plt.xlabel('Soret coefficient, $\lambda$', fontsize=12)
plt.ylabel('Spring constant, $\kappa$', fontsize=12)
plt.gca().invert_yaxis()


#uploading the mathematica parametric plot data for critical parameters at phase transition
mathematicaoutputdata = np.loadtxt('mathematicaexport.txt')        
lambda_= mathematicaoutputdata[:,0]
kappa = mathematicaoutputdata[:,1]
#scaling up the data to match the heatmap plot
scaled_lamb =(lambda_*20 )
scaled_kap = (kappa*20)
plt.plot(scaled_lamb,scaled_kap, color='aqua', linewidth=3)
plt.xlim([1,10])
plt.ylim([1,20])

plt.subplot(122)

cbar_ticks2 = [math.pow(10, i) for i in range(math.floor(math.log10(test2.min().min())), 1+math.ceil(math.log10(test2.max().max())))]

sns.heatmap(freq_smooth, vmin = np.min(freq_smooth[1:20,1:10]), vmax = np.max(freq_smooth[1:20,1:10]))
plt.yticks([0,2,4,6,8,10,12,14,16,18,20],[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
plt.xticks([0,1,2,3,4,5,6,7,8,9,10],[0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5] )
plt.title('Particle oscillation frequency',style='italic')
plt.xlabel('Soret coefficient, $\lambda$', fontsize=12)
plt.ylabel('Spring constant, $\kappa$', fontsize=12)
plt.gca().invert_yaxis()
plt.plot(scaled_lamb,scaled_kap, color='aqua', linewidth=3)

plt.xlim([1,10])
plt.ylim([1,20])
plt.show()


fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(15.2,4.8))
plt.subplot(131)
fig2 = plt.plot(r_p,v_p, label="Simulation")
#plt.plot(rpSA,vpSA, label="Stability Analysis")
plt.title('(a)')
plt.xlabel('Position, X(t)', fontsize=12, fontweight='bold', style='italic')
plt.ylabel(r"Velocity, $\mathbf{ \dot X(t)}$", fontsize=12, fontweight='bold', style='italic')
plt.legend()

plt.subplot(132)
fig3 = plt.plot(t_data, r_data)
plt.title('(b)')
#plt.title("Particle trajectory in 1d: Chemotaxis + Harmonic potential \n  D = 0.1, $\mu$ = 0.1, $\sigma$=0.3,  $\kappa$= 0.6, $\lambda$ = 0.25")
plt.ylabel('Position, X(t)', fontsize=12, fontweight='bold', style='italic')
plt.xlabel('Time, t', fontsize=12, fontweight='bold', style='italic')
#plt.xlim([0,300])



uncomment below until
mathematicaoutput_amplitudes = np.loadtxt('mathematicaexport_3d.txt')  
mathematicaoutput_frequencies = np.loadtxt('mathematicaexport_3d_freqs.txt')      
amp_smooth_SA =  gaussian_filter(mathematicaoutput_amplitudes[:,:10], sigma= 0)
freq_smooth_SA = gaussian_filter(mathematicaoutput_frequencies[:,:10], sigma= 0)

fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(11,6))
plt.subplot(121)

sns.heatmap(mathematicaoutput_amplitudes[:,:10],vmin = np.min(amp_smooth_SA[1:20,1:10]), vmax = np.max(amp_smooth_SA[1:20,1:10]))
# predicted by stability analysis \n  D = 0.1, $\mu$ = 0.1, $\sigma$=0.3
plt.title('Amplitude of first harmonic mode oscillations')
plt.xlabel('Soret coefficient, $\lambda$')
plt.ylabel('Spring constant, $\kappa$')
plt.gca().invert_yaxis()
plt.yticks([0,2,4,6,8,10,12,14,16,18,20],[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
plt.xticks([0,1,2,3,4,5,6,7,8,9,10],[0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5] )
#plt.xticks([0,1,2,3,4,5,6,7,8,9,10],[0.00,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5] )
#plt.yticks([0,2,4,6,8,10,12,14,16,18],[0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95])
plt.plot(scaled_lamb,scaled_kap, color='lime', linewidth=3) 
plt.xlim([1,10])
plt.ylim([1,20])

plt.subplot(122)
sns.heatmap(mathematicaoutput_frequencies[:,:10],vmin = np.min(freq_smooth_SA[1:20,1:10]), vmax = np.max(freq_smooth_SA[1:20,1:10]))
plt.title('Frequency of first harmonic mode oscillations')
plt.xlabel('Soret coefficient, $\lambda$')
plt.ylabel('Spring constant, $\kappa$')
plt.gca().invert_yaxis()
plt.yticks([0,2,4,6,8,10,12,14,16,18,20],[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
plt.xticks([0,1,2,3,4,5,6,7,8,9,10],[0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5])
plt.plot(scaled_lamb,scaled_kap, color='lime', linewidth=3) 
plt.xlim([1,10])
plt.ylim([1,20])
plt.show()


signal_amplitudes = np.loadtxt('mathematicaexport_3d_signal.txt')  
#signalamp_smooth_SA =  gaussian_filter(signal_amplitudes[:,:10], sigma= 0)

sns.heatmap(signal_amplitudes,vmin = np.min(signal_amplitudes), vmax = np.max(signal_amplitudes))
plt.title('Amplitude of oscillations predicted by Cubic order stability analysis')
plt.xlabel('Soret coefficient, $\lambda$')
plt.ylabel('Spring constant, $\kappa$')
plt.gca().invert_yaxis()
plt.yticks([0,2,4,6,8,10,12,14,16,18,20],[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
plt.xticks([0,1,2,3,4,5,6,7,8,9,10],[0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5] )
plt.plot(scaled_lamb,scaled_kap, color='lime', linewidth=3) 
plt.xlim([1,10])
plt.ylim([1,20])


duration = 100
samplingrate= 5
num_samples = int(duration*samplingrate)
#num_samples = len(r_data)

# apply fast fourier transform and take absolute values\\
freq_spectrum=np.abs(np.fft.fft(r_p))/num_samples 
 
#freq_spectrum=np.abs(np.fft.fft(r_data))
# get the list of frequencies
frequencies = np.fft.fftfreq(num_samples, 1/samplingrate)
plt.subplot(133)
# plot frequency spectrum, note may need to use a semilog scale
plt.plot(frequencies*2*np.pi,freq_spectrum)
plt.title('(c)')
plt.xlim([-3.5,3.5])
plt.xlabel(r"Angular frequency, $ \mathbf {\omega_0 [s^{-1}]}$", fontsize=12, fontweight='bold', style='italic')
plt.ylabel("Spectral amplitude [m s]", fontsize=12, fontweight='bold', style='italic')
plt.suptitle( r"Stable dynamics ($ \mathbf {D= 0.1, \mu= 0.1, \sigma=0.3, \kappa= 0.2, \lambda= 0.6}$)", fontsize=12, fontweight='bold')
plt.show()

def freq_above_threshold(freq_spectrum, threshold,frequencies):
    x_coordinates = []
    for i in range(len(freq_spectrum)):
        if freq_spectrum[i] > threshold :
            x_coordinates.append(frequencies[i])
    return x_coordinates

#find the index of the maximum amplitude (excludiing the DC component - this represents the mean position value and not the oscillatory behaviour)

max_ampl = np.max(freq_spectrum)
max_ampl_index = np.argmax(freq_spectrum)

#print("The root mean square amplitude of oscillation is", msd_1d(r_data))
print("The dominant (linear) frequency of oscillation is", frequencies[max_ampl_index], "Hz")

print("The dominant (angular) frequency of oscillation is", 2*np.pi*frequencies[max_ampl_index], "Hz")
threshold = 1

#***until here***
#frequency_modes =freq_above_threshold(freq_spectrum, threshold, frequencies)
#print("The frequencies observed in the spectrum are:",frequency_modes)

#numerical integration of the fourier frequency spectrum using trapezoidal rule
#defining integral bounds from w_0 - w_0/2 to w_0 + w_0/2


#fa = frequencies[max_ampl_index] - frequencies[max_ampl_index]/2
#fb = frequencies[max_ampl_index] + frequencies[max_ampl_index]/2
fa = 0
fb = 0.5
#defining bounds from 3w_0 - 3w_0/2 to 3w_0 + 3w_0/2
#fc = frequencies[3*max_ampl_index] - frequencies[3*max_ampl_index]/2
#fd = frequencies[3*max_ampl_index] + frequencies[3*max_ampl_index]/2

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

a = find_nearest(frequencies, fa)
b = find_nearest(frequencies, fb)
#c = find_nearest(frequencies, fc)
#d = find_nearest(frequencies, fd)

integral_frequencies = frequencies[a:b]
spectrum = freq_spectrum[a:b]
#cubintegral_frequencies = frequencies[c:d]
#cubspectrum = freq_spectrum[c:d]

integralvalue = np.trapz(spectrum, integral_frequencies)
#cubicintegralvalue = np.trapz(cubspectrum, cubintegral_frequencies)

print("The numerical integration result across the w0 peak is:", integralvalue)
#print("The numerical integration result across the 3w0 peak is:", cubicintegralvalue)
#***include these in below**
print("rp amplitude is ",np.max(r_data)  )
print("freq spectrum peak at w0 is of height:", np.max(freq_spectrum))


def signal(t,A1,w0,phi1, A3, phi2):
    return A1*np.cos(w0*t +phi1) + A3*np.cos(3*w0*t +phi2)
#ydata = signal(t,0.2,0.0001,0.44,phi1,phi2)

plt.plot(t,r_p,'b-', label='simulation')
popt, pcov = curve_fit(signal, t, r_p, p0 = [0.2,0.44,0, 0.001, 1])
popt
#plt.plot(t,signal(t,*popt), 'r-', label = 'fit: A1=%5.3f, w0=%5.3f,phi1=%5.3f' % tuple(popt))

plt.plot(t,signal(t,*popt), 'r-', label = 'fit: A1=%5.3f, w0=%5.3f, phi1=%5.3f, A3=%5.3f, phi2=%5.3f' % tuple(popt))
#include label for the fitted parameters!!
plt.xlabel('t')
plt.ylabel('rp')
plt.legend(loc='upper center')
plt.show()
print("rp amplitude is ",np.max(r_p)  )

----
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(L,K,Z, rstride=1, cstride=1,cmap='viridis', edgecolour='none')
plt.xlabel('Soret coefficient, $\lambda$')
plt.ylabel('Spring constant, $\kappa$')
ax.set_title('surface title')
"""

from matplotlib import cm
from matplotlib.ticker import LinearLocator

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Make data.
K = np.flip(np.arange(0.01, 1, 0.05))
#lambda below
L = np.arange(0.01, 0.5, 0.05)
L, K = np.meshgrid(L, K)
Z =  np.loadtxt('mathematicaexport_3d_signal.txt') 

print(np.meshgrid(L, K))
print(Z)

# Plot the surface.
surf = ax.plot_surface(L, K, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

# Customize the z axis.
ax.set_zlim(0, 1.01)
#ax.zaxis.set_major_locator(LinearLocator(10))
# A StrMethodFormatter is used automatically
#ax.zaxis.set_major_formatter('{x:.02f}')

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()