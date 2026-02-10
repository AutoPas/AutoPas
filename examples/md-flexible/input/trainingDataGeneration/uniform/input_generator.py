from string import Template
import numpy as np
import os

# Uses the template_input.yaml to creates various input files that emulate a
# simulation domain distributed uniformly with particles, with different thread
# counts and skin sizes. For each scenario, 5 repeat runs are made.

inputTemplateFile = open("template_input.yaml", "r")

inputTemplate = Template(inputTemplateFile.read())

verlet_skins = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
thread_counts = np.array([6, 12, 18, 24, 30, 36])
num_particles = np.array([100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 51200, 102400, 204800])

# Create directory structure:
# <num_particles>/<thread_count>/<skin>/<run>
for num_particle in num_particles:
    os.mkdir(f"./noParticles_{num_particle:0>6}")
    os.chdir(f'./noParticles_{num_particle:0>6}')
    
    for thread_count in thread_counts:
        os.mkdir(f"./{thread_count:0>2}_threads")
        os.chdir(f'./{thread_count:0>2}_threads')
        for verlet_skin in verlet_skins:
            os.mkdir(f"./skin_{verlet_skin}")
            os.chdir(f'./skin_{verlet_skin}')
            
            for run in range(5):
                os.mkdir(f"./run_{run:0>1}")
                os.chdir(f"./run_{run:0>1}")
                
                dictionary = {
                    'skin' : verlet_skin,
                    'numParticles' : num_particle
                }

                
                        
                f = open("./input.yaml", "w")
                    
                f.write(inputTemplate.substitute(dictionary))
                    
                f.close()
                    
                os.chdir('..')
            os.chdir('..')
        os.chdir("..")
    os.chdir("..")
    
    
