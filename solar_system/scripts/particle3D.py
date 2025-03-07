"""
CompMod Ex2: Particle3D, a class to describe point particles in 3D space

An instance describes a particle in Euclidean 3D space: 
velocity and position are [3] arrays

Author: Nick Mansell
Number: s2017896
        

"""
import numpy as np


class Particle3D(object):
    """
    Class to describe point-particles in 3D space

    Attributes
    ----------
    label: name of the particle
    mass: mass of the particle
    position: position of the particle
    velocity: velocity of the particle

    Methods
    -------
    __init__
    __str__
    kinetic_energy: computes the kinetic energy
    momentum: computes the linear momentum
    update_position_1st: updates the position to 1st order
    update_position_2nd: updates the position to 2nd order
    update_velocity: updates the velocity

    Static Methods
    --------------
    read_file: initializes a P3D instance from a file handle
    total_kinetic_energy: computes total K.E. of a list of particles
    com_velocity: computes centre-of-mass velocity of a list of particles
    """

    def __init__(self, label, mass, position, velocity):
        self.label = str(label)
        self.mass = float(mass)
        self.position = np.array(position)
        self.velocity = np.array(velocity)
        
        
        """
        Initialises a particle in 3D space.

        Parameters
        ----------
        label: str
            name of the particle
        mass: float
            mass of the particle
        position: [3] float array
            position vector
        velocity: [3] float array
            velocity vector
        """

    def __str__(self):
        """
        Return an XYZ-format string. The format is
        label    x  y  z

        Returns
        -------
        str
        """
        return str(self.label)+" "+str(self.position[0])+" "+str(self.position[1])+" "+str(self.position[2])

    def kinetic_energy(self):
        """
        Returns the kinetic energy of a Particle3D instance

        Returns
        -------
        ke: float
            1/2 m v**2
        """
        ...
        ke = 0.5*(self.mass)*((np.linalg.norm(self.velocity))**2)
        
        return ke   # ...

    def momentum(self):
        """
        Returns the momentum of a Particle3D instance
        
        Returns
        -------
        momentum: numpy float array
                  m*v
        """
        momentum = self.mass*self.velocity
        return momentum

    def update_position_1st(self, dt):
        """
        Uodates the position of a Particle3D instance
        
        Parameters
        -------
        dt: float
            time interval between original and updated position
        """
        self.position = self.position + dt*self.velocity
        
    
    def update_position_2nd(self, dt, f):         # f = force
        """
        Uodates the position of a Particle3D instance being acted on by a force
        
        Parameters
        -------
        dt: float
            time interval between original and updated position
        f: numpy float array
           time dependant force
        """
        self.position = self.position + dt*self.velocity + (dt**2)*(f/(2*self.mass))
        
    def update_velocity(self, dt, f):            # f = force
        """
        Uodates the velocity of a Particle3D instance being acted on by a force
        
        Parameters
        ----------
        dt: float
            time interval between original and updated velocity
        f: numpy float array
           time dependant force
        """
        self.velocity = self.velocity + dt*(f/self.mass)

    @staticmethod
    def read_line(line):
        """
        Creates a Particle3D instance given a line of text.

        The input line should be in the format:
        label   <mass>  <x> <y> <z>    <vx> <vy> <vz>

        Parameters
        ----------
        filename: str
            Readable file handle in the above format

        Returns
        -------
        p: Particle3D
        """
        
        variables = line.split()
        
        label = str(variables[0])
        mass = float(variables[1])
        position = np.array([float(variables[2]),float(variables[3]),float(variables[4])])
        velocity = np.array([float(variables[5]),float(variables[6]),float(variables[7])])               
        
        p = Particle3D(label,mass,position,velocity)
        
        return p

    @staticmethod
    def total_kinetic_energy(particles):
        """
        Computes the total kinetic energy of a list of Particle3D's
        
        Parameters
        ----------
        particles: list
                   a list of Particle3D instances
        
        Returns
        -------
        ke_total: float
                  total kinetic energy
        """
        ke_total = 0
        for particles in particles:
            ke = particles.kinetic_energy()
            ke_total += ke
            
        return ke_total

    @staticmethod
    def com_velocity(particles):           #v_com =(1/total mass) * sum of momentums
        """
        Computes the CoM velocity of a list of P3D's

        Parameters
        ----------
        particles: list
            A list of Particle3D instances

        Returns
        -------
        com_vel: array
            Centre-of-mass velocity
        """
        mv_tot = 0
        for particles in particles:
            mv = particles.momentum()
            mv_tot += mv
            return mv_tot
        
        mass_tot = 0
        for particles in particles:
            mass = particles.mass
            mass_tot += mass
            return mass_tot
        
        com_vel = mv_tot/mass_tot
        return com_vel
        
        
        
