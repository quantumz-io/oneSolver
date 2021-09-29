import csv
import numpy as np
import matplotlib.pyplot as plt

lingpudata = np.zeros((5,3))

with open('annlinearongpugen9.csv', 'r') as csvfile:
    csvreader = csv.reader(csvfile, delimiter=',')
    line_count = 0
    for row in csvreader:
        if line_count == 0:
            linschedulertype = row
            line_count += 1
        elif line_count == 1:
            linnames = row
            line_count += 1
        else:
            lingpudata[line_count-2,:] = row
            line_count += 1

geomgpudata = np.zeros((5,3))

with open('anngeometricongpugen9.csv', 'r') as csvfile:
    csvreader = csv.reader(csvfile, delimiter=',')
    line_count = 0
    for row in csvreader:
        if line_count == 0:
            geomschedulertype = row
            line_count += 1
        elif line_count == 1:
            geomnames = row
            line_count += 1
        else:
            geomgpudata[line_count-2,:] = row
            line_count += 1

#<codecell>
plt.plot(lingpudata[:,0], lingpudata[:,1], 'bo-', label=linschedulertype[1])
plt.plot(geomgpudata[:,0], geomgpudata[:,1], 'go-', label=geomschedulertype[1])
plt.legend(title="Scheduler type")
plt.grid()
plt.title('Intel(R) UHD Graphics P630')
plt.xlabel('Number of nodes')
plt.ylabel('Energy')
plt.savefig('AnnealingEnergyOnGPU.svg', format='svg', dpi=1200)
#<codecell>

#<codecell>
plt.plot(lingpudata[:,0], lingpudata[:,2], 'bo-', label=linschedulertype[1])
plt.plot(geomgpudata[:,0], geomgpudata[:,2], 'go-', label=geomschedulertype[1])
plt.legend(title="Scheduler type")
plt.grid()
plt.title('Intel(R) UHD Graphics P630')
plt.xlabel('Number of nodes')
plt.ylabel('Calculation time [s]')
plt.savefig('AnnealingTimeOnGPU.svg', format='svg', dpi=1200)
#<codecell>


lincpudata = np.zeros((5,3))

with open('annlinearonxeon.csv', 'r') as csvfile:
    csvreader = csv.reader(csvfile, delimiter=',')
    line_count = 0
    for row in csvreader:
        if line_count == 0:
            cpulinschedulertype = row
            line_count += 1
        elif line_count == 1:
            cpulinnames = row
            line_count += 1
        else:
            lincpudata[line_count-2,:] = row
            line_count += 1

geomcpudata = np.zeros((5,3))

with open('anngeomonxeon.csv', 'r') as csvfile:
    csvreader = csv.reader(csvfile, delimiter=',')
    line_count = 0
    for row in csvreader:
        if line_count == 0:
            cpugeomschedulertype = row
            line_count += 1
        elif line_count == 1:
            cpugeomnames = row
            line_count += 1
        else:
            geomcpudata[line_count-2,:] = row
            line_count += 1
#<codecell>
plt.plot(lincpudata[:,0], lincpudata[:,1], 'bo-', label=cpulinschedulertype[1])
plt.plot(geomcpudata[:,0], geomcpudata[:,1], 'go-', label=cpugeomschedulertype[1])
plt.legend(title="Scheduler type")
plt.grid()
plt.title('Intel(R) Xeon(R) E-2146G CPU @ 3.50GHz')
plt.xlabel('Number of nodes')
plt.ylabel('Energy')
plt.savefig('AnnealingEnergyOnCPU.svg', format='svg', dpi=1200)
#<codecell>

#<codecell>
plt.plot(lincpudata[:,0], lincpudata[:,2], 'bo-', label=linschedulertype[1])
plt.plot(geomcpudata[:,0], geomcpudata[:,2], 'go-', label=geomschedulertype[1])
plt.legend(title="Scheduler type")
plt.grid()
plt.title('Intel(R) Xeon(R) E-2146G CPU @ 3.50GHz')
plt.xlabel('Number of nodes')
plt.ylabel('Calculation time [s]')
plt.savefig('AnnealingTimeOnCPU.svg', format='svg', dpi=1200)
#<codecell>

#<codecell>
plt.plot(lingpudata[:,0], lingpudata[:,1], 'bo-', label=linschedulertype[1]+'(gpu)')
plt.plot(geomgpudata[:,0], geomgpudata[:,1], 'go-', label=geomschedulertype[1]+'(gpu)')
plt.plot(lincpudata[:,0], lincpudata[:,1], 'b*--', label=cpulinschedulertype[1]+'(cpu)')
plt.plot(geomcpudata[:,0], geomcpudata[:,1], 'g*--', label=cpugeomschedulertype[1]+'(cpu)')
plt.legend(title="Scheduler type", ncol=2)
plt.grid()
plt.title('Intel(R) UHD Graphics P630 vs Intel(R) Xeon(R) E-2146G CPU @ 3.50GHz')
plt.xlabel('Number of nodes')
plt.ylabel('Energy')
plt.savefig('AnnealingEnergyOnGPUvsCPU.svg', format='svg', dpi=1200)
#<codecell>

#<codecell>
plt.plot(lingpudata[:,0], lingpudata[:,2], 'bo-', label=linschedulertype[1]+'(gpu)')
plt.plot(geomgpudata[:,0], geomgpudata[:,2], 'go-', label=geomschedulertype[1]+'(gpu)')
plt.plot(lincpudata[:,0], lincpudata[:,2], 'b*--', label=cpulinschedulertype[1]+'(cpu)')
plt.plot(geomcpudata[:,0], geomcpudata[:,2], 'g*--', label=cpugeomschedulertype[1]+'(cpu)')
plt.legend(title="Scheduler type", ncol=2)
plt.grid()
plt.title('Intel(R) UHD Graphics P630 vs Intel(R) Xeon(R) E-2146G CPU @ 3.50GHz')
plt.xlabel('Number of nodes')
plt.ylabel('Calculation time [s]')
plt.savefig('AnnealingTimeOnGPUvsCPU.svg', format='svg', dpi=1200)
#<codecell>
