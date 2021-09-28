import csv
import numpy as np
import matplotlib.pyplot as plt

gpudata = np.zeros((5,6))

with open('exhongpugen9.csv', 'r') as csvfile:
    csvreader = csv.reader(csvfile, delimiter=',')
    line_count = 0
    for row in csvreader:
        if line_count == 0:
            gpunames = row
            line_count += 1
        else:
            gpudata[line_count-1,:] = row
            line_count += 1

#<codecell>
plt.plot(gpudata[:,0], gpudata[:,1], 'bo-', label=gpunames[1])
plt.plot(gpudata[:,0], gpudata[:,2], 'go-', label=gpunames[2])
plt.plot(gpudata[:,0], gpudata[:,3], 'ro-', label=gpunames[3])
plt.plot(gpudata[:,0], gpudata[:,4], 'mo-', label=gpunames[4])
plt.plot(gpudata[:,0], gpudata[:,5], 'co-', label=gpunames[5])
plt.legend()
plt.grid()
plt.title('Intel(R) UHD Graphics P630')
plt.xlabel('Number of nodes')
plt.ylabel('Calculation time [s]')
plt.semilogy()
plt.savefig('ExhaustiveOnGPU.svg', format='svg', dpi=1200)
#<codecell>

cpudata = np.zeros((5,6))
with open('exhonxeon.csv', 'r') as csvfile:
    csvreader = csv.reader(csvfile, delimiter=',')
    line_count = 0
    for row in csvreader:
        if line_count == 0:
            cpunames = row
            line_count += 1
        else:
            cpudata[line_count-1,:] = row
            line_count += 1


#<codecell>
plt.plot(cpudata[:,0], cpudata[:,1], 'bo-', label=cpunames[1])
plt.plot(cpudata[:,0], cpudata[:,2], 'go-', label=cpunames[2])
plt.plot(cpudata[:,0], cpudata[:,3], 'ro-', label=cpunames[3])
plt.plot(cpudata[:,0], cpudata[:,4], 'mo-', label=cpunames[4])
plt.plot(cpudata[:,0], cpudata[:,5], 'co-', label=cpunames[5])
plt.legend()
plt.grid()
plt.title('Intel(R) Xeon(R) E-2176G CPU @ 3.50GHz')
plt.xlabel('Number of nodes')
plt.ylabel('Calculation time [s]')
plt.semilogy()
plt.savefig('ExhaustiveOnCPU.svg', format='svg', dpi=1200)
#<codecell>



#<codecell>
plt.plot(gpudata[:,0], gpudata[:,1], 'bo-', label=gpunames[1]+'(gpu)')
plt.plot(gpudata[:,0], gpudata[:,2], 'go-', label=gpunames[2]+'(gpu)')
plt.plot(gpudata[:,0], gpudata[:,3], 'ro-', label=gpunames[3]+'(gpu)')
plt.plot(gpudata[:,0], gpudata[:,4], 'mo-', label=gpunames[4]+'(gpu)')
plt.plot(gpudata[:,0], gpudata[:,5], 'co-', label=gpunames[5]+'(gpu)')

plt.plot(cpudata[:,0], cpudata[:,1], 'b*--', label=cpunames[1]+'(cpu)')
plt.plot(cpudata[:,0], cpudata[:,2], 'g*--', label=cpunames[2]+'(cpu)')
plt.plot(cpudata[:,0], cpudata[:,3], 'r*--', label=cpunames[3]+'(cpu)')
plt.plot(cpudata[:,0], cpudata[:,4], 'm*--', label=cpunames[4]+'(cpu)')
plt.plot(cpudata[:,0], cpudata[:,5], 'c*--', label=cpunames[5]+'(cpu)')
plt.legend(ncol=2)
plt.grid()
plt.title('Intel(R) UHD Graphics P630 vs Intel(R) Xeon(R) E-2176G CPU @ 3.50GHz')
plt.xlabel('Number of nodes')
plt.ylabel('Calculation time [s]')
plt.semilogy()
plt.savefig('ExhaustiveGPUvsCPU.svg', format='svg', dpi=1200)
#<codecell>
