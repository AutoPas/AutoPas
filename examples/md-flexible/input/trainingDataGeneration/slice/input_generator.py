from string import Template
import numpy as np
import os

# Uses the template_input.yaml to creates various input files that emulate a
# simulation domain with a big slice in the middle, for different thread counts
# and skin sizes. For each scenario, 5 repeat runs are made.

inputTemplateFile = open("template_input.yaml", "r")

inputTemplate = Template(inputTemplateFile.read())

verlet_skins = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
thread_counts = np.array([6, 12, 18, 24, 30, 36])
slice_sizes = np.array([5, 10, 15, 20])

# Create folder structure:
# <slice_size>/<thread_count>/<verlet_skin>/<run>
for slice_size in slice_sizes:
    os.mkdir(f"./slice_size_{slice_size:0>2}")
    os.chdir(f'./slice_size_{slice_size:0>2}')
    slice_start = 200 - slice_size / 2
    
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
                    'slice_size'      : slice_size,
                    'slice_start'     : slice_start
                }

                
                        
                f = open("./input.yaml", "w")
                    
                f.write(inputTemplate.substitute(dictionary))
                    
                f.close()
                    
                
                os.chdir('..')
            os.chdir('..')
        os.chdir("..")
    os.chdir("..")
    
    
