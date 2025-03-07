"""
Nick Mansell
s2017896
"""


from particle3D import Particle3D
import numpy as np



def compute_separations(particles):
    """
    Calculates the vector separations between a list of particle3D class 
    particles.

    Parameters
    ----------
    particles: List of particle3D particles.
        

    Returns
    -------
    separations: 3D numpy array containing vector separations
    of particles in the input list.
        

    """
    n = len(particles)
    #initialise output array
    separations = np.zeros([n,n,3])
    #start loops to cycle through different particle combinations
    for i in range(n):
        for j in range(n):
            #if function to avoid calculations involving the same particle twice
            #as well as making the same calculation twice
            #==> If j is ever less than i then the inverse of that calculation 
            #will have been done already
            if i < j:
                #calculate separation vector
                separation = particles[i].position - particles[j].position
                #add calculated separations to array
                separations[i,j] = separation
                separations[j,i] = -separation

    return separations




def compute_forces_potential(particles,separations):
    """
    

    Parameters
    ----------
    particles: List of particle 3D particles.
   
    separations: 3D numpy array containing vector separations
    of particles in the input list.

    Returns
    -------
    forces: 3D array containing the gravitational forces between the particles.
    
    potential: Total gravitational potential energy of the system of particles.

    """
    G = 8.887692593e-10
    n = len(particles)
    #initialise forces array and potential energy
    forces = np.zeros([n,3])
    potential = 0
    #start loops to cyle through the different combinations of particles
    for i in range(n):
        for j in range(n):
            #if function same as find_separations
            if i < j:
                #calculate gravitational force
                force = -G*particles[i].mass*particles[j].mass*separations[i,j]/(np.linalg.norm(separations[i,j])**3)
                #add calculated force to forces array
                forces[i] += force
                #add inverse of calculated force to array
                forces[j] -= force
    
                #calculate potential energy and add to potential
                potential += -G*particles[i].mass*particles[j].mass/np.linalg.norm(separations[i,j])
                
    return forces, potential

#a,b = compute_forces_potential(particles, compute_separations(particles))
#print(a)
#print(b)
