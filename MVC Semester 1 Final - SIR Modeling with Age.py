#Code based on following paper: https://www.nature.com/articles/s41598-021-94609-3

import matplotlib.pyplot as plt
import numpy as np


gamma=0.1
beta=0.4


#data for Pleasanton! (https://censusreporter.org/profiles/16000US0657792-pleasanton-ca/)
susceptible=np.array([0.11, 0.14, 0.8, 0.11, 0.18, 0.15, 0.12, 0.07, 0.04])*74660 #total number susceptible
infected=np.array([10, 10, 10, 10, 10, 10, 10, 10, 10]) #total number infected
recovered=np.array([0, 0, 0, 0, 0, 0, 0, 0, 0]) #total number recovered
N_i=susceptible+infected+recovered

age_matrix=np.array([[19.2, 4.8, 3.0, 7.1, 3.7, 3.1, 2.3, 1.4, 1.4], 
            [4.8, 42.4, 6.4, 5.4, 7.5, 5.0, 1.8, 1.7, 1.7], 
            [3.0, 6.4, 20.7, 9.2, 7.1, 6.3, 2.0, 0.9, 0.9], 
            [7.1, 5.4, 9.2, 16.9, 10.1, 6.8, 3.4, 1.5, 1.5], 
            [3.7, 7.5, 7.1, 10.1, 13.1, 7.4, 2.6, 2.1, 2.1], 
            [3.1, 5.0, 6.3, 6.8, 7.4, 10.4, 3.5, 1.8, 1.8], 
            [2.3, 1.8, 2.0, 3.4, 2.6, 3.5, 7.5, 3.2, 3.2], 
            [1.4, 1.7, 0.9, 1.5, 2.1, 1.8, 3.2, 7.2, 7.2], 
            [1.4, 1.7, 0.9, 1.5, 2.1, 1.8, 3.2, 7.2, 7.2]])/100.0


S_over_time=[]
I_over_time=[]
R_over_time=[]
times=[]

for day in range(350):

    #record all data at this current time stamp
    S_over_time.append(susceptible.copy())
    I_over_time.append(infected.copy())
    R_over_time.append(recovered.copy())
    times.append(day)

    #modify susceptible, infected, recovered variables
    for i in range(9): #looping through i values (see Nature paper)
        total=0
        for j in range(9):
            total+=age_matrix[i][j]*infected[j]
        total/=N_i[i]

        
        dS_dt=-beta*susceptible[i]*total
        dI_dt=beta*susceptible[i]*total-gamma*infected[i]
        dR_dt=gamma*infected[i]

        susceptible[i]+=dS_dt
        infected[i]+=dI_dt
        recovered[i]+=dR_dt 


plt.plot(times, [coln[0]/N_i[0]*100 for coln in I_over_time], label="0-9")
plt.plot(times, [coln[1]/N_i[1]*100 for coln in I_over_time], label="10-19")
plt.plot(times, [coln[2]/N_i[2]*100 for coln in I_over_time], label="20-29")
plt.plot(times, [coln[3]/N_i[3]*100 for coln in I_over_time], label="30-39")
plt.plot(times, [coln[4]/N_i[4]*100 for coln in I_over_time], label="40-49")
plt.plot(times, [coln[5]/N_i[5]*100 for coln in I_over_time], label="50-59")
plt.plot(times, [coln[6]/N_i[6]*100 for coln in I_over_time], label="60-69")
plt.plot(times, [coln[7]/N_i[2]*100 for coln in I_over_time], label="70-79")
plt.plot(times, [coln[8]/N_i[8]*100 for coln in I_over_time], label="80+")
plt.title("SIR Modeling with Age Consideration: Infected Proportions over Time")
plt.xlabel("Time (days)")
plt.ylabel("% of Subpopulation Infected")
plt.legend()
plt.show()

plt.plot(times, [coln[0]/N_i[0]*100 for coln in R_over_time], label="0-9")
plt.plot(times, [coln[1]/N_i[1]*100 for coln in R_over_time], label="10-19")
plt.plot(times, [coln[2]/N_i[2]*100 for coln in R_over_time], label="20-29")
plt.plot(times, [coln[3]/N_i[3]*100 for coln in R_over_time], label="30-39")
plt.plot(times, [coln[4]/N_i[4]*100 for coln in R_over_time], label="40-49")
plt.plot(times, [coln[5]/N_i[5]*100 for coln in R_over_time], label="50-59")
plt.plot(times, [coln[6]/N_i[6]*100 for coln in R_over_time], label="60-69")
plt.plot(times, [coln[7]/N_i[2]*100 for coln in R_over_time], label="70-79")
plt.plot(times, [coln[8]/N_i[8]*100 for coln in R_over_time], label="80+")
plt.title("SIR Modeling with Age Consideration: Recovered Proportions over Time")
plt.xlabel("Time (days)")
plt.ylabel("% of Subpopulation Recovered")
plt.legend()
plt.show()

