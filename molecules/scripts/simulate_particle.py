"""
Symplectic Euler and Velocity Verlet time integration of
a particle moving in a double well potential.

Produces plots of the position of the particle
and its energy, both as function of time. Also
saves both to file.

The potential is V(x) = a*x^4 - b*x^2, where
a and b are hard-coded in the main() method
and passed to the functions that
calculate force and potential energy.

Author: Nick Mansell
Number: s2017896

"""

import sys
import math
import numpy as np
import matplotlib.pyplot as pyplot
from particle3D import Particle3D
from scipy.signal import find_peaks


def force_morse(p1, p2, alpha, D_e, r_e):
    """
    Return the morse force of two particles.

    The force is given by
        F(x) = -dV/dx = 

    Parameters
    ----------
    p1: Particle3D
    alpha: float
    D_e: float
    r_e: float

    Returns
    -------
    force_1: float
    force_2: float
    """
    r1 = p1.position
    r2 = p2.position
    r_12 = np.linalg.norm(r2-r1)
    rhat_12 = (r2-r1)/r_12
    
    force = (2*alpha*D_e)*(1-math.exp(-alpha*(r_12-r_e)))*math.exp(-alpha*(r_12-r_e))*rhat_12
    return force


def potential_morse(p1, p2, alpha, D_e, r_e):
    """
    Method to return potential energy 
    of particle in double-well potential
    V(x) = a*x^4 - b*x^2

    Parameters
    -----------
    p1: Particle3D
    a: float
    b: float

    Returns
    -------
    potential: float
    """
    r1 = p1.position
    r2 = p2.position
    r_12 = np.linalg.norm(r2-r1)
    potential = D_e*(((1-math.exp(-alpha*(r_12-r_e)))**2)-1)
    return potential

   
def find_frequency(y_values, time_values):
    """
    Function to return the wavenumber and frequency of an oscillating particle.
    wavenumber = f/c = 1/c*period
    frequency = c*wavenumber = 1/period
    
    Parameters
    ---------
    y_values: 1D numpy array
    time_values: 1D numpy array
    
    Returns
    ------
    wavenumber: float
    frequency: float
    """
    peaks = find_peaks(y_values)                 
    peaks = list(peaks[0])                          
    peaks_range = time_values[peaks[len(peaks)-1]] - time_values[peaks[0]]
    period = (peaks_range/(len(peaks)-1))*(1.018050571e-14)               
    
    wavenumber = 1/(period*(3e10))
    frequency = 1/period
    return wavenumber, frequency


def initial_frequency(y_values, time_values):
    peaks = find_peaks(y_values)
    peaks = list(peaks[0])
    f_initial = 1/((time_values[peaks[1]] - time_values[peaks[0]])*1.018050571e-14)

    return f_initial



def main():
    # Read inputs from command line
    # The variable sys.argv contains whatever you typed on the command line
    # or after %run on the ipython console to launch the code.  We can use
    # it to get user inputs.
    # Here we expect three things:
    #    the name of this file
    #    euler or verlet
    #    the name of the output file the user wants to write to
    # So we start by checking that all three are specified and quit if not,
    # after giving the user a helpful error message.
    if len(sys.argv) != 4 :
        print("You left out the name of either the input or output file (or both) when running.")
        print("In spyder, run like this instead:")
        print(f"    %run {sys.argv[0]} <euler or verlet> <desired input file> <desired output file>")
        sys.exit(1)
    else:
        mode = sys.argv[1]
        outfile_name = sys.argv[3]
        infile_name = sys.argv[2]

    # Open the output file for writing ("w")
    outfile = open(outfile_name, "w")
    
    # Open the input file for reading ("r")
    infile = open(infile_name, "r")

    # Choose our simulation parameters
    lines = infile.readlines()
    dt, numstep = lines[1].split() 
    D_e, r_e, alpha = lines[3].split()
    time = 0.0
    
    dt = float(dt)
    numstep = int(numstep)
    D_e = float(D_e)
    r_e = float(r_e)
    alpha = float(alpha)
    
    # Set up particle initial conditions:
    p1 = Particle3D.read_line(lines[4])
     
    p2 = Particle3D.read_line(lines[5])
    
    
    # Get initial force
    force = force_morse(p1, p2, alpha, D_e, r_e)

    # Write out starting time, position, and energy values
    # to the output file.
    energy = p1.kinetic_energy() + potential_morse(p1, p2, alpha, D_e, r_e) + p2.kinetic_energy()
    outfile.write(f"{time}    {p1.position}   {p2.position}   {energy}\n")
    

    # Initialise numpy arrays that we will plot later, to record
    # the trajectories of the particles.
    times = np.zeros(numstep)
    energies = np.zeros(numstep)
    separation = np.zeros(numstep)
    # Start the time integration loop
    for i in range(numstep):

        # Update the positions and velocities.
        # This will depend on whether we are doing an Euler or verlet integration
        if mode == "euler":
            # Update particle position
            p1.update_position_1st(dt)
            
            # Update particle position
            p2.update_position_1st(dt)
            
            # Calculate force
            force = force_morse(p1, p2, alpha, D_e, r_e)

            # Update particle velocity 
            p1.update_velocity(dt, force)
            
            # Update particle velocity 
            p2.update_velocity(dt, -force)
            
            
        elif mode == "verlet":
            # Update particle position using previous force
            p1.update_position_2nd(dt, force)
            
            # Update particle position using previous force
            p2.update_position_2nd(dt, -force)
            
            # Get the force value for the new positions
            force_new = force_morse(p1, p2, alpha, D_e, r_e)

            # Update particle velocity by averaging
            # current and new forces
            p1.update_velocity(dt, 0.5*(force+force_new))
            
            # Update particle velocity by averaging
            # current and new forces
            p2.update_velocity(dt, 0.5*(-(force+force_new)))
            
            # Re-define force value for the next iteration
            force = force_new
        else:
            raise ValueError(f"Unknown mode {mode} - should be euler or verlet")

        # Increase time
        time += dt
        
        # Output particle information
        energy = p1.kinetic_energy() + potential_morse(p1, p2, alpha, D_e, r_e) + p2.kinetic_energy()
        outfile.write(f"{time} {p1.position} {p2.position} {energy}\n")
        

        # Store the things we want to plot later in our arrays
        times[i] = time
        separation[i] = np.linalg.norm(p2.position-p1.position)
        energies[i] = energy
        

    # Now the simulation has finished we can close our output file
    outfile.close()
    
    #Sets the label based on the particle being observed
    if infile_name == "oxygen.dat":
        linelabel = "Oxygen"
        
    elif infile_name == "nitrogen.dat":
        linelabel = "Nitrogen"
    
    elif infile_name == "oxygenspin.dat":
        linelabel = "Oxygen with spin"
    
    else:
        linelabel = "Nitrogen with spin"

    #Set the title of the plot depending on whch integration method is in use    
    if mode == "euler":
        title = 'Symplectic Euler: '
    
    else:
        title = 'Verlet: '
        

    # Plot particle trajectory to screen. There are no units
    # here because it is an abstract simulation, but you should
    # include units in your plot labels!
    pyplot.title(title + 'Separation vs Time')
    pyplot.xlabel('Time (ang*amu/eV)')
    pyplot.ylabel('Separation (ang)')
    pyplot.plot(times, separation, 'k', label = linelabel, linewidth = 0.7)
    pyplot.legend(bbox_to_anchor =(0.3,-0.27), loc='lower right')
    pyplot.show()

    # Plot particle energy to screen
    pyplot.title(title + 'Total Energy vs Time')
    pyplot.xlabel('Time (ang*amu/eV)')
    pyplot.ylabel('Energy (eV)')
    pyplot.plot(times, energies, 'k', label = linelabel, linewidth = 0.7)
    pyplot.legend(bbox_to_anchor =(0.3,-0.27), loc='lower right')
    pyplot.show()
  
    # Calculate and print the frequency 
    wavenumber, frequency = find_frequency(separation,times)    
    print("The frequency of oscillations of this molecule is: " + str(frequency) + " Hz")


    # Calculate and print the energy innacuracy of the simulation
    emax = np.max(energies)
    emin = np.min(energies)
    e0 = energies[0]
    e_innacuracy = (emax - emin)/e0
    print("The energy innacuracy of this simulation is: " + str(abs(e_innacuracy)))

    print("The wavenumber of this molecule is: " + str(wavenumber) + " cm^-1")
    
    print(str(initial_frequency(separation,times)))
# This strange but standard python idiom means that the main function will
# only be run if we run this file, not if we just import it from another
# python file. It is good practice to include it whenever your code can be
# run as a program.
if __name__ == "__main__":
    main()

