import numpy as np
import matplotlib.pyplot as plt
import sys

#try:
#    folder = input('Please input folder: ')
#except ValueError:
#    print("input folder number")
folder = '.'

# Read the SMD.cv file
data_file = folder+'/SMD.cv'  # Replace with the actual filename or path to your SMD.cv file
data = np.loadtxt(data_file, skiprows=1)  # Skip the header line

# Extract the relevant columns
startIndex = 0
#endIndex = 326000
endIndex = -1
timestep = [ (i -data[startIndex,0])/1000 for i in data[startIndex:endIndex, 0]]  # in ps
actual_distance = data[startIndex:endIndex, 3]  # Actual CV
leading_distance = data[startIndex:endIndex, 4] # Objective leading CV
work = data[startIndex:endIndex, -1]  # work done

# Compute the free energy profile
free_energy = -np.log(np.exp(-work))

# Plot the following situation
plt.subplot(2, 1, 1)
plt.plot(timestep, actual_distance, label='actual SN2 CV value ')
plt.plot(timestep, leading_distance, label='leading SN2 CV value')
plt.xlabel('Time elapse (ps)')
plt.ylabel('Reaction coordinate (Angstroms)')
plt.legend()

# Plot the free energy profile with x-axis flipped
plt.subplot(2, 1, 2)
plt.plot(leading_distance, free_energy) 
plt.xlabel('Reaction coordinate (Angstroms)')
plt.ylabel(r'$\Delta F$ (kJ/mol)')
plt.title('Free Energy Profile')
# Flip the x-axis tick labels
# plt.gca().invert_xaxis()

# Adjust the spacing between the subplots
plt.subplots_adjust(hspace=0.5)
plt.show()
