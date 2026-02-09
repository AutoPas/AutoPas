from string import Template
import numpy as np
import os

# Uses the template_input.yaml to creates various input files that emulate 
# heterogeneous scenarios using uniformly distributed clusters of gaussian-ly
# distributed particles, with different thread counts and skin sizes. For each
# scenario, 5 repeat runs are made.

inputTemplateFile = open("template_input.yaml", "r")

inputTemplate = Template(inputTemplateFile.read())

number_of_clusters = np.array([10, 20, 40, 80, 160])
verlet_skins = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
thread_counts = np.array([6, 12, 18, 24, 30, 36])

# Create directory structure:
# <num_clusters>/<thread_count>/<verlet_skin>/<run>
for num_cluster in number_of_clusters:
    os.mkdir(f"./noClusters_{num_cluster:0>3}")
    os.chdir(f'./noClusters_{num_cluster:0>3}')
    
    for thread_count in thread_counts:
        os.mkdir(f"./{thread_count:0>2}_threads")
        os.chdir(f'./{thread_count:0>2}_threads')
        for verlet_skin in verlet_skins:
            os.mkdir(f"./skin_{verlet_skin}")
            os.chdir(f'./skin_{verlet_skin}')
            
            for run in range(5):
                os.mkdir(f"./run_{run:0>1}")
                os.chdir(f"./run_{run:0>1}")


                # Generate Object String
                objectStr = ''
                # For each cluster, add a gaussian object
                for cluster in range(num_cluster):
                    # randomly generate a centre not too close to any wall
                    centre = np.random.uniform(2, 18, 3)

                    objectStr += f'''    {cluster}:
      distribution-mean              :  [{centre[0]}, {centre[1]}, {centre[2]}]
      distribution-stddeviation      :  [2, 2, 2]
      numberOfParticles              :  100
      box-length                     :  [20, 20, 20]
      bottomLeftCorner               :  [0, 0, 0]
      velocity                       :  [0, 0, 0]
      particle-type-id               :  0
'''
                
                dictionary = {
                    'skin' : verlet_skin,
                    'objects' : objectStr
                }

                
                        
                f = open("./input.yaml", "w")
                    
                f.write(inputTemplate.substitute(dictionary))
                    
                f.close()
                    
                os.chdir('..')
            os.chdir('..')
        os.chdir("..")
    os.chdir("..")
    
    
