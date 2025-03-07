The simulation takes a .txt file as an input of the form: name mass x y z vx vy vz
(for example: solar_system.txt)

to run: python full_source_code.py <timestep> <total steps> <input filename> <trajectory output filename> <energy output filename> <yes/no>
yes/no indicates if you would like to calcukate apside and periods and print to terminal

trajectory_output.XYZ can be animated using: plot_xyz.py trajectory_output.XYZ

simulation units:
time - days
mass - earth masses
distance - astronomical units