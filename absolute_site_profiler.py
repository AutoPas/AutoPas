import subprocess

#files_to_sample = ["smallscaleExplodingLiquidSites1", "smallscaleExplodingLiquidSites2", "smallscaleExplodingLiquidSites5", "smallscaleExplodingLiquidSites10"]
files_to_sample = ["smallscaleExplodingLiquidSites5"]
number_of_samples = 2

for file_to_sample in files_to_sample:
    for sample in range(number_of_samples):
        result = subprocess.run(["./build/examples/md-flexible/md-flexible", "--yaml-filename", "./build/examples/md-flexible/" + file_to_sample + ".yaml"], capture_output=True, text = True)
        file = open("./profile_results_new_functor/" + file_to_sample + "_sample" + str(sample), 'a+')
        file.write(result.stdout)

#result = subprocess.run(["./build/examples/md-flexible/md-flexible", "--yaml-filename", "./build/examples/md-flexible/" + files_to_sample[0] + ".yaml"], capture_output=True, text = True)
#print(result)