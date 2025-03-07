"""
Nick Mansell
s2017896
"""

from particle3D import Particle3D 
import numpy as np
from basic_functions import compute_separations, compute_forces_potential
import sys


def read_initial_conditions(infile):
    """
    Reads in a list of particles from a text file and creates an initial list of Particle3D particles
    
    Parameters
    ----------
    infile: text file containing a set of
            particles
        
    Returns
    -------
    particles: List of Particle3D instances
    """
    #open file
    bodies = open(infile, "r")
    #create P3D instance for each line
    lines = bodies.readlines()
    particles = []
    #add each instance to list of particles
    for i in range(len(lines)):
        particles.append(Particle3D.read_line(lines[i]))

    return particles


def writeout_xyz(particles, i, outfile):    
    """
    Writes the trajectory of all particles to a file in the XYZ format
    
    Parameters
    ----------
    particles: List of Particle3D instances
    i : Integer value, represents the current step of the simulation
    outfile: Name of file trajectories will be saved to

    Returns
    -------
    None.
    
    """
    n = len(particles)
    #write no. particles and point to file
    outfile.write(f"{n}\npoint = {i+1}\n")
    #loop over all particles
    for j in range(n):
        #write to file
        outfile.write(f"{particles[j].__str__()}\n")
    

def extract_xvalues(separations, numstep, z):
    """
    Creates an array of all the x separation values for a given particle relative to the sun

    Parameters
    ----------
    separations: Array of shape (numstep, #particles, #particles, 3) containing all the separations of the particles for a given simulation
                 relatice to eachother                
    numstep: Integer, number of timesteps taken in simulation
    z: Integer, index of particle in the particles list

    Returns
    -------
    x_values: 1D array 

    """
    #initialise empty array
    x_values = np.zeros(numstep)
    #loop over all time iterations
    for i in range(numstep):
        #add x separation of current timestep to array
        x_values[i] = separations[i, z, 0, 0]
        
    return x_values
        
    
def compute_observables(separations, numstep, particle, particles):
    """
    Computes the Aphelion, Perihelion and Period of a particle for a given simulation. First checks to see if a full orbit has 
    been completed, if not returns N/A

    Parameters
    ----------
    separations: Array of shape (numstep, #particles, #particles, 3) containing all the separations of the particles for a given 
                 simulation relative to eachother
    numstep: Integer, number of timesteps taken in simulation
    particle: Integer, index of particle in the particles list
    particles: List of Particle3D instances

    Returns
    -------
    period: Float
    aphelion: Float
    perihelion: Float

    """
    result = 0
    #get x values from separations array
    x = extract_xvalues(separations, numstep, particle)
    step = 1
    
    #loop over all x values except the first and last
    for i in range(1, numstep-1):
        #stop if step reaches numstep
        if step >= numstep:
            break
        
        #check to see if current point is passing over first point > must use a range rather than an exact value as due to the
        #nature of the simulation it will never be exactly the same
        if (x[i-1] < x[0] and x[i+1] > x[0]):
            #if yes, add 1 to ressult, 5 to step and break loop. add 5 to step to jump ahead 5 spaces in the list of x values
            #this is done as it is possible for 2 values in a row can satisfy the if loop which would count as two successes insted
            #of just 1
            result += 1
            step += 5
            break
        
        #same as above but inversed
        elif (x[i-1] > x[0] and x[i+1] < x[0]):
            result += 1
            step += 5
            break
        
        #if no, add 1 to step and continue
        else:
            step += 1
            pass
        
    #same loops as above but from step to numstep-1,     
    for i in range(step, numstep-1):
        
        if step >= numstep:
            break 
        #if it crosses again, break > once it crosses twice a full orbit has been completed
        if (x[i-1] < x[0] and x[i+1] > x[0]):
            result += 1
            break
        
        elif (x[i-1] > x[0] and x[i+1] < x[0]):
            result += 1
            break
        
        else:
            step += 1
            pass

    #set orbit to true or false based on if a full orbit has been completed
    if result < 2:
        orbit = False
    else:
        orbit = True
        
    
    if orbit == True:    
        #if true, period = whatever value step was when loop ended
        period = step
        
        #initiate empty array
        distances = np.zeros(numstep)
        #loop over numstep
        for i in range(numstep):
            #calculate distance from sun and add to array for given particle
            distances[i] = np.linalg.norm(separations[i, particle, 0])
        #calculate aphelion and perihelion
        aphelion = np.max(distances)
        perihelion = np.min(distances)

    #if false, all observables = N/A        
    else:
        period = str("N/A")
        aphelion = str("N/A")
        perihelion = str("N/A")
    
    return period, aphelion, perihelion    


def energy_innac(energies):
    e_file = open(energies, "r")
    e_list = e_file.readlines()
    e_file.close()
    return (float(max(e_list)) - float(min(e_list)))/float(e_list[0])

def main():
    #if command line not formatted correctly, print message showing correct format
    if len(sys.argv) != 7:
        print("Incorrect command, try: %run full_source_code.py <Time step> <Total steps> <Input file> <Output file> <Energy file> <Calculate Observables? (yes/no)>")
        sys.exit(1)
    
    #read inputs from command line and open all output files
    infile = sys.argv[3]
    outfile_name = sys.argv[4]
    energies_file_name = sys.argv[5]
    outfile = open(outfile_name, "w")
    energies_file = open(energies_file_name, "w")
    #energies_file.write("#TotaL energy of simulation (Earth masses * AU^2 * Days^-2)\n")
    
    #read in particles
    particles = read_initial_conditions(infile)
    #get centre of mass velocity of particles
    com_vel = Particle3D.com_velocity(particles)
    #subtract centre of mass velocity from all particles 
    for i in range(len(particles)):
        particles[i].velocity -= com_vel
    
    #generate variables and initial arrays
    numstep = int(sys.argv[2])
    dt = float(sys.argv[1])
    n = len(particles)
    separations = np.zeros([numstep,n,n,3])
    times = np.zeros(numstep)
    time = 0.0
    force, potential = compute_forces_potential(particles, compute_separations(particles))
    
    #loop over given numstep
    for i in range(numstep):
        #add current trajectory of particles to file
        writeout_xyz(particles, i, outfile)
        
        #add current separations and time to arrays
        separations[i] = compute_separations(particles)
        times[i] = time
        
        #write current total energy of system to energy file
        energies_file.write(f"{potential + (Particle3D.total_kinetic_energy(particles))}\n")       

        #update position of all particles        
        for j in range(n):
            
            particles[j].update_position_2nd(dt, force[j])
        
        #calculate new force
        force_new, potential = compute_forces_potential(particles, compute_separations(particles))

        #update velocity of all particles
        for k in range(n):
            
            particles[k].update_velocity(dt, 0.5*(force[k] + force_new[k]))
            
        #make force = new_force for next iteration
        force = force_new
        #increase time by given dt
        time += dt
    
    #close files
    outfile.close()
    energies_file.close()
    
    #determines if use wants to calculate and print apsides and periods
    if sys.argv[6] == "yes":
        #initiate observable list
        observables = []
        #loop over all particles except sun as measurements are taken relative to the sun. technically the sun does orbit the centre
        #of mass of the solar system but this simulation in its current state cannot calculate this
        for i in range(1, n):
            period, aphelion, perihelion = compute_observables(separations, numstep, i, particles)
            
            #if the observable are N/A, add one N/A instead of doing it 3 times
            if period == str("N/A"):
                observables.append(f"{particles[i].label}: N/A")
            
            #add observables to list
            else:    
                observables.append(f"{particles[i].label}: Aphelion = {aphelion} AU  Perihelion = {perihelion} AU  Period = {period*dt} Days")
        
        #print observables to terminal
        print(observables)
        
    #print(energy_innac(energies_file_name))
   

main()