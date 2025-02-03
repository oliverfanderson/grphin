import numpy as np
import sys
import ast

# Reads in Altaf's files stress_orbit_significance.csv and re-calculates significance values. 
NAMES = ['bsub','elegans','fly','drerio','cerevisiae']
FILES = {name:'final_output/%s/stress_orbit_significance.csv' % name for name in NAMES}
NONZEROS = {name:[] for name in NAMES}
SIG = {name:[] for name in NAMES}
p = 0.01  ## p-value is set to 0.01.
def main():
    sig_orbits = set()
    for name in NAMES:
        f = FILES[name]
        print(name)
        print(f)
        header = False
        with open(f) as fin:
            for line in fin:
                if not header: # skip first line (header)
                    header = True
                    continue
                row = line.strip().split('\t')
                orbit = int(row[0])
                sig = int(row[1])
                obs_median = float(row[2])
                if obs_median == 0:
                    continue
                random_sampling_str = row[3].replace("np.float64(", "").replace(")", "")
                random_sampling = ast.literal_eval(random_sampling_str)
                rand_biggerthan_obs = sum([x>=obs_median for x in random_sampling])
                NONZEROS[name].append([orbit,obs_median,rand_biggerthan_obs/1000])
                if rand_biggerthan_obs/1000 < p:
                    SIG[name].append([orbit,obs_median,rand_biggerthan_obs/1000,random_sampling])
                    print('--> SIGNIFICANT! (<%.2f)'% (p),[orbit,obs_median,rand_biggerthan_obs/1000])
                    sig_orbits.add(orbit)
    

    for name in NAMES:
        print(name,len(SIG[name]))

    print('ORBITS:')
    orbits_to_sort = []
    for orbit in sorted(sig_orbits):
        species = []
        for name in NAMES:
            if any(x[0] == orbit for x in SIG[name]):
                species.append(name)

        #print(orbit,len(species),species)
        orbits_to_sort.append([orbit,len(species),species])

    print('COMMON ORBITS:')
    for o,num,species in sorted(orbits_to_sort, key=lambda x: x[1], reverse=True):
        if num < 2:
            break
        print(o,num,species)
        for name in NAMES:
            for val in NONZEROS[name]:
                if val[0] == o:
                    print('  ',name,val)



if __name__ == "__main__":
    main()
