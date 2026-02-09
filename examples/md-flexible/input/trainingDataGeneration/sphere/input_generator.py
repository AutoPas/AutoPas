from string import Template
import numpy as np
import os

# Uses the template_input.yaml to creates various input files that emulate a
# simulation domain empty except for a sphere of particles, with different
# thread counts and skin sizes. For each scenario, 5 repeat runs are made.

inputTemplateFile = open("template_input.yaml", "r")

inputTemplate = Template(inputTemplateFile.read())

verlet_skins = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
thread_counts = np.array([6, 12, 18, 24, 30, 36])
radii = np.array([5, 10, 15, 20])
spacings = np.array([0.5, 1, 1.5, 2])

# Create directory structure:
# <radius>/<spacing>/<thread_count>/<verlet_skin>/<run>
for radius in radii:
    os.mkdir(f"./radius_{radius:0>2}")
    os.chdir(f'./radius_{radius:0>2}')

    for spacing in spacings:
        os.mkdir(f"./spacing_{spacing:.1f}")
        os.chdir(f'./spacing_{spacing:.1f}')
    
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
                        'skin'            : verlet_skin,
                        'radius'          : radius,
                        'spacing'         : spacing
                    }

                
                        
                    f = open("./input.yaml", "w")
                    
                    f.write(inputTemplate.substitute(dictionary))
                    
                    f.close()

                    os.chdir('..')
             

                os.chdir('..')
            os.chdir('..')
        os.chdir("..")
    os.chdir("..")
    
    
