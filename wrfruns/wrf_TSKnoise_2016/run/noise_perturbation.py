import numpy as np
from netCDF4 import Dataset

# Open wrfinput file in read and write mode
ds = Dataset('wrfinput_d03', 'r+')
#print(ds.variables)
# Get TSK (skin temperature) - shape is (Time, south_north, west_east)
tsk = ds.variables['TSK']
# Generate random noise with magnitude 1e-5
noise = np.random.uniform(-1e-5, 1e-5, size = tsk.shape)
# Apply noise
tsk[:] = tsk[:] + noise
ds.close()
print('TSK perturbed with random noise of magnitude 1e-5 K')