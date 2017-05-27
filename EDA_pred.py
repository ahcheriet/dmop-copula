from PyGMO import *
from pandas import DataFrame
from math import *
from PyGMO.problem import base
import matplotlib.pylab as plt
from CEDA import *
from copula_estimator import *
import sys
import math as mat
import csv
import subprocess
import dt

from dproblems import *
prob_str = sys.argv[1]
n_individus = eval(sys.argv[2])


def plot_front(pop, fig):
    p_list = pop.compute_pareto_fronts()
    f = p_list[0]
    plt.figure(fig)
    plt.plot([pop[ind].cur_f[0] for ind in f], [pop[ind].cur_f[1]
                                                for ind in f], 'o')
def plot_file(xt,fig):
    filename = '../true_pf/'+prob_str+'/POF-nt'+str(
        dt.nt) + '-taut' + str(dt.taut) + '-'+prob_str+'-' + str(xt) + '.txt'
    print filename
    fn = open(filename,'rb')
    reader = csv.reader(fn,delimiter='\t')
    ltmp = []
    for i in reader:
        ltmp.append([eval(i[0]),eval(i[1])])
    plt.figure(fig)        
    plt.plot([item[0] for item in ltmp],[item[1] for item in ltmp],'o')                                              


def to_file(pop, filename):
    prob = pop.problem
    f_dim = prob.f_dimension
    best_idx1 = pop.compute_pareto_fronts()[0]
    functionValue = []
    filename = filename
    for i in best_idx1:
        functionValue.append(pop[i].cur_f)
    with open(filename + ".txt", "wb") as f:
        writer = csv.writer(f, delimiter=' ')
        List = [list(functionValue[i]) for i in range(len(functionValue))]
        List.sort()
        writer.writerows(List)
    f.close()


def Hypervolume(pop, dim, xt):
    to_file(pop, 'pop_tmp')
    p = subprocess.Popen(['./Hypervolume_main', 'pop_tmp.txt', '../true_pf/DMOP2-nt10-taut10/POF-nt' + str(
        dt.nt) + '-taut' + str(dt.taut) + '-DMOP2-' + str(xt) + '.txt', str(dim)], stdout=subprocess.PIPE)
    out, err = p.communicate()
    return (eval(out))


def hv(pop, filename, xt):
    to_file(pop, filename)
    p = subprocess.Popen(['./hv', filename + '.txt'], stdout=subprocess.PIPE)
    out, err = p.communicate()
    return (eval(out))



def IGD(pop, dim, filename, xt):
    to_file(pop, filename)
    p = subprocess.Popen(['./InvertedGenerationalDistance_main', filename + '.txt', '../true_pf/'+prob_str+'/POF-nt'+str(
        dt.nt) + '-taut' + str(dt.taut) + '-'+prob_str+'-' + str(xt) + '.txt', str(dim)], stdout=subprocess.PIPE)
    out, err = p.communicate()
    return (eval(out))

def Backup(pop,pop_all):
    idx = pop.compute_pareto_fronts()[0]
    for i in idx:
        pop_all.push_back(pop[i].cur_x)
    return pop_all


def Archive(pop, pop_new, prob):
    algo = algorithm.moead(gen=1)
    pop_next = population(prob, 1)

    idx = pop.compute_pareto_fronts()[0]
    for i in idx:
        pop_new.push_back(pop[i].cur_x)
    pop_new = algo.evolve(pop_new)
    idx = pop_new.get_best_idx(n_individus)
    pop_next.set_x(0, pop_new[idx[0]].cur_x)
    for i in idx[1:]:
        pop_next.push_back(pop_new[i].cur_x)
    return pop_next


prob = eval(prob_str)()
pt = 0.0
sum_igd = 0

DataListX_eda = []
DataListY_eda = []
cop_est = copula_estimator(n=100)


alg2 = ceda_moea(gen=2, n=300)
algo = algorithm.moead(gen=2)
pop = population(prob, n_individus)
pop1 = population(pop)
pop_all = population(prob, 1)
rate = 0

"""
for xt in range(1, 201, 1):
    print xt
    dt.tau = xt
    t = float(1) / float(dt.nt)
    t = t * floor(float(dt.tau) / float(dt.taut))
    if (xt + 1) % dt.taut == 0:
        plot_file(xt,1)
        plt.savefig('../fig/PF_at_' + str(xt) + '.png')
      #  plot_front(pop, 1)
      #  plot_file(xt)
      #  plt.show()

"""

for xt in range(1, dt.Taut, 1):
    print xt
    dt.tau = xt
    t = float(1) / float(dt.nt)
    t = t * floor(float(dt.tau) / float(dt.taut))
    
    if (xt + 1) % dt.taut == 0:
        quality_eda = IGD(pop, pop.problem.f_dimension, 'eda', xt)
        sum_igd = sum_igd + quality_eda
        DataListX_eda.append(xt)
        DataListY_eda.append(quality_eda)
        
    if t != pt:
        pt = t
        rate += 1 # Count the first change
        if rate == 1:
            print "I saved it"
            pop_all = Backup(pop,pop_all)
        if rate == 2:
            rate = 0
            print "Am using it"
            pop1 = cop_est.evolve(pop_all,pop)
            pop = population(Archive(pop, pop1, prob))
            pop_all = population(prob, 1)
        
    pop = algo.evolve(pop)

#print sum_igd/len(DataListX)
plt.plot(DataListX_eda[1:], DataListY_eda[1:],'-o')

plt.show()

