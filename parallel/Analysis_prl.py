##python script for paralelization analysis
import matplotlib.pyplot as plt
#SECTION1: getdata
def read_SERIAL(file_SERIAL):
    with open(file_SERIAL, 'r') as file:
        first_line = file.readline().strip()
        # Split the line into parts
        parts = first_line.split()
        try:
            SERIE = float(parts[-1])
        except ValueError:
            print("Error: Unable to convert the last part to a float.")
            return None
        return SERIE

def extract_data(file_PARALEL):
    nprocs=[]
    times=[]
    with open(file_PARALEL, 'r') as file:
        for line in file:
            parts2 = line.strip().split()
            try:
                num_procs = int(parts2[1])
                run_time = float(parts2[-1])
            except ValueError:
                print("Error: Unable to convert values on line:", line)
                quit
            nprocs.append(num_procs)
            times.append(run_time)
    return nprocs, times

##      ->extract SERIAL TIME
file_SERIAL = "timeSERIE.dat"  
SERIE = read_SERIAL(file_SERIAL)
print(f'serial time ts (s)={SERIE}')
##      -> Now for PARALEL sets of NPROCS - TIME
file_PARALEL = "PARALELIZATION.dat"  
nprocs, times = extract_data(file_PARALEL)
for num_procs, run_time in zip(nprocs, times):
    print(f"{num_procs}procs - Time: {run_time}")

#print(nprocs[1]+nprocs[-1])
#print(times[1]+times[-1])
##SECTION2: funtion definition of speedup, efficiency

def speedup(tp,ts):
    return float(ts/tp)

def eff(speed, np):
    return float(speed/np)


#SECTION3: get datas2plot

SPEEDUP=[]
EFFICIENCY=[]

#speedup
if (len(nprocs)==len(times)):
    for x in times:
        SPEEDUP.append(speedup(float(x),SERIE))
        print('speedup=', SPEEDUP[-1])

#efficiency
for sp,n in zip(SPEEDUP,nprocs):
    EFFICIENCY.append(eff(sp, n))
    print('eff=',EFFICIENCY[-1])


##SECTION3: plots

plt.figure(1)
plt.plot(nprocs, EFFICIENCY, label='Eff(P)')
plt.xlabel('P')
plt.ylabel('Efficiency=f(P)')
plt.savefig('EffvsP.jpg')
plt.show()

plt.figure(2)
plt.plot(nprocs, SPEEDUP, label='Speedup(P)')
plt.xlabel('P')
plt.ylabel('Speedup=f(P)')
plt.savefig('SpeedupvsP.jpg')
plt.show()

plt.figure(3)
plt.plot(nprocs, times)
plt.ylabel('Time (s)')
plt.xlabel('P')
plt.savefig('timevsP.jpg')
plt.show()
