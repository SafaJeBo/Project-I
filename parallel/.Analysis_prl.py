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

tp2=times[0] #time 2 proc paralelized vs SERIE time serie2proc
#print(nprocs[1]+nprocs[-1])
#print(times[1]+times[-1])
##SECTION2: funtion definition of speedup, efficiency

def speedup(tp,ts):
    return float(ts/tp)

def eff(speed, np):
    return float(speed/np)

def serialf(SERIE,tp2):
    f=float((SERIE-tp2)/SERIE)
    return abs(f)
print(times[0])

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
plt.plot(nprocs, EFFICIENCY, color='crimson', marker='o', markersize='2', label='Eff(P)')
plt.xlabel('P')
plt.ylabel('Efficiency=f(P)')
plt.savefig('EffvsP.png')
plt.show()

plt.figure(2)
plt.plot(nprocs, SPEEDUP, color='darkolivegreen', marker='d', markersize='2', label='Speedup(P)')
plt.xlabel('P')
plt.ylabel('Speedup=f(P)')
plt.savefig('SpeedupvsP.png')
plt.show()

plt.figure(3)
plt.plot(nprocs, times, color='g', marker='^', markersize='2', label='Time/proc')
plt.ylabel('Time (s)')
plt.xlabel('P')
plt.savefig('timevsP.png')
plt.show()


## SERIAL FRACTION:

f=serialf(SERIE,tp2)  #f in seconds
ff= f/SERIE * 100 #f en %
t_par= ff*SERIE/100 #s
t_nopar = SERIE-t_par
print(f'Code has {ff:.1f}% of parallelizable code')
print('Hence, from the total execution time in serie t_s(s)={:.4f}'.format(SERIE), '\n'
      'which means: {:.4f} s serial execution (parallelizable)\n'
      '\t\t {:.4f} s serial code (NO parallelizable)'.format(t_par, t_nopar))

import numpy as np
y_id= np.array([SERIE, 1, 1/SERIE])
x_id= np.flip(y_id)
plt.figure(4)
plt.loglog(nprocs, [t_par/(n+t_nopar) for n in nprocs], color='firebrick', linestyle='-.',\
    marker='x', markersize='2', label=f't(s)={t_par:.4f}r/nproc+{t_nopar:.4f}')
plt.loglog(x_id, y_id, color='navy', label='ideal')
plt.xlabel('P')
plt.legend()
plt.title('Serial fraction: t(s) vs P')
plt.ylabel('Time t(s)')
plt.savefig('TimeFvsP.png')
plt.show()

