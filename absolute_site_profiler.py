import subprocess


container_prefixes = ["LC"] #["VL","LC"]
data_layouts = ["SoA"] #,"AoS"]
site_counts = ["1", "2", "5"]
densities = ["1", "0.75", "0.5"]
functors = ["NewFunc"] #["NewFunc","OldFunc"]

#files_to_sample = ["smallscaleExplodingLiquidSites1", "smallscaleExplodingLiquidSites2", "smallscaleExplodingLiquidSites5", "smallscaleExplodingLiquidSites10"]
#files_to_sample = ["smallscaleExplodingLiquidSites5"]

samples = 3

for container_prefix in container_prefixes:
    for data_layout in data_layouts:
        for site_count in site_counts:
            for density in densities:
                for functor in functors:
                    for sample in range(samples):
                        path_to_input_file = "./build/examples/md-flexible/"
                        filename = container_prefix + data_layout + "Sites" + site_count + "Density" + density + functor
                        print("handling " + filename)
                        result = subprocess.run(["./build/examples/md-flexible/md-flexible", "--yaml-filename", path_to_input_file + filename+".yaml"], capture_output=True, text = True)
                        path_to_output_folder = "./profile_results_new_functor/" + container_prefix + data_layout + "/"
                        file = open(path_to_output_folder + filename + ".txt", 'a+')
                        file.write(result.stdout)

#for (container_prefix, data_layout, site_count, density, functor, sample) in zip(container_prefixes, data_layouts, site_counts, densities, functors, samples):
#    path_to_input_file = "./build/examples/md-flexible/"
#    filename = container_prefix + data_layout + "Sites" + site_count + "Density" + density + functor
#    print("handling " + filename)
#    result = subprocess.run(["./build/examples/md-flexible/md-flexible", "--yaml-filename", path_to_input_file + filename+".yaml"], capture_output=True, text = True)
#    path_to_output_folder = "./profile_results_new_functor/" + container_prefix + data_layout + "/"
#    file = open(path_to_output_folder + filename + ".txt", 'a+')
#    #file.write(result.stdout)

#for file_to_sample in files_to_sample:
#    for sample in range(number_of_samples):
#        result = subprocess.run(["./build/examples/md-flexible/md-flexible", "--yaml-filename", "./build/examples/md-flexible/" + file_to_sample + ".yaml"], capture_output=True, text = True)
#        file = open("./profile_results_new_functor/" + file_to_sample + "_sample" + str(sample), 'a+')
#        file.write(result.stdout)

#result = subprocess.run(["./build/examples/md-flexible/md-flexible", "--yaml-filename", "./build/examples/md-flexible/" + files_to_sample[0] + ".yaml"], capture_output=True, text = True)
#print(result)