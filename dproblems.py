from PyGMO import population, problem
from PyGMO.problem import base


#from  CEDA import *
import csv
import subprocess


import math as mat
import dt
import dynamic_benchmark



def plot_front(pop, nfig, a=40, comp=[0, 1, 2]):
    """
    Generic plot-method for multi-objective optimization problems with more then 2 objectives
    """
    from mpl_toolkits.mplot3d import axes3d
    import matplotlib.pyplot as plt
    import numpy as np
    
    if pop.problem.f_dimension == 2:
        p_list = pop.compute_pareto_fronts()
        f = p_list[0]
        fig = plt.figure(nfig)
        ax = fig.add_subplot(111)
        ax.plot([pop[ind].cur_f[0] for ind in f],[pop[ind].cur_f[1] for ind in f],'o')	
    else:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        fit = np.transpose([ind.cur_f for ind in pop])
        ax.plot(fit[comp[0]], fit[comp[1]], fit[comp[2]], 'ro')
        ax.view_init(azim=a)
    return ax


def to_file(pop,filename):
   prob = pop.problem
   f_dim = prob.f_dimension
   best_idx1 = pop.compute_pareto_fronts()[0]
   functionValue = []
   filename = filename
   for i in best_idx1:
       functionValue.append(pop[i].cur_f)
   with open(filename+".txt", "wb") as f:
      writer = csv.writer(f,delimiter=' ')
      List = [list(functionValue[i]) for i in  range(len(functionValue))]
      List.sort()
      writer.writerows(List)    
   f.close()
        


def Archive(pop,pop_new,prob):
   algo = algorithm.moead(gen=1)
   pop_next = population(prob,1)
   
   idx = pop.compute_pareto_fronts()[0]
   for i in idx:
      pop_new.push_back(pop[i].cur_x)
   pop_new = algo.evolve(pop_new)
   
   idx = pop_new.get_best_idx(n_individus)
   pop_next.set_x(0,pop_new[idx[0]].cur_x)
   for i in idx[1:]:
      pop_next.push_back(pop_new[i].cur_x)
   return pop_next




class FDA4(base):
    """FDA4 dynamic benchmark class for pygmo problem
    """
    def __init__(self, dim = 12):
        super(FDA4,self).__init__(dim,0,3)
        lb = []
        ub = []
        for i in range(dim-1):
            lb.append(0.0)
            ub.append(1.0)
        self.set_bounds(0.0,1.0)                
        self.__dim = dim

    def _objfun_impl(self,x):
        print dt.nt
        f1,f2,f3 = dynamic_benchmark.FDA4(x,dt.tau,dt.nt,dt.taut)
        return (f1,f2,f3,)
      
    def human_readable_extra(self):
        return "\n\t Problem dimension: " + str(self.__dim)


class FDA5(base):
    """FDA5 dynamic benchmark class for pygmo problem
    """
    def __init__(self, dim = 12):
        super(FDA5,self).__init__(dim,0,3)
        lb = []
        ub = []
        for i in range(dim-1):
            lb.append(0.0)
            ub.append(1.0)
        self.set_bounds(0.0,1.0)                
        self.__dim = dim

    def _objfun_impl(self,x):
        f1,f2,f3 = dynamic_benchmark.FDA5(x,dt.tau,dt.nt,dt.taut)
        return (f1,f2,f3,)
      
    def human_readable_extra(self):
        return "\n\t Problem dimension: " + str(self.__dim)



class DIMP2(base):
    """DIMP2 dynamic benchmark class for pygmo problem
    """
    def __init__(self, dim = 10):
        super(DIMP2,self).__init__(dim,0,2)
        lb = []
        ub = []
        lb.append(0.0)
        ub.append(1.0)
        for i in range(1, dim):
            lb.append(-2.0)
            ub.append(2.0)
        self.set_bounds(lb,ub)
        self.__dim = dim

    def _objfun_impl(self,x):
        f1,f2 = dynamic_benchmark.DIMP2(x,dt.tau,dt.nt,dt.taut)
        return (f1,f2,)
      
    def human_readable_extra(self):
        return "\n\t Problem dimension: " + str(self.__dim)

class DMOP2(base):
    """dMOP2 dynamic benchmark class for pygmo problem
    """
    def __init__(self, dim = 10):
        super(DMOP2,self).__init__(dim,0,2)
        lb = []
        ub = []
        for i in range(dim-1):
            lb.append(0.0)
            ub.append(1.0)
        self.set_bounds(0.0,1.0)                
        self.__dim = dim

    def _objfun_impl(self,x):
        f1,f2 = dynamic_benchmark.dMOP2(x,dt.tau,dt.nt,dt.taut)
        return (f1,f2,)
      
    def human_readable_extra(self):
        return "\n\t Problem dimension: " + str(self.__dim)

class DMOP3(base):
    """dMOP2 dynamic benchmark class for pygmo problem
    """
    def __init__(self, dim = 10):
        super(DMOP3,self).__init__(dim,0,2)
        lb = []
        ub = []
        for i in range(dim-1):
            lb.append(0.0)
            ub.append(1.0)
        self.set_bounds(0.0,1.0)                
        self.__dim = dim
        self.r = -1
        self.rIteration = -1

    def _objfun_impl(self,x):
        f1,f2,self.r,self.rIteration = dynamic_benchmark.dMOP3(x,dt.tau,dt.nt,dt.taut,self.r,self.rIteration)
        return (f1,f2,)
      
    def human_readable_extra(self):
        return "\n\t Problem dimension: " + str(self.__dim)


class HE2(base):
    """HE2 dynamic benchmark class for pygmo problem
    """
    def __init__(self, dim = 30):
        super(HE2,self).__init__(dim,0,2)
        lb = []
        ub = []
        for i in range(dim-1):
            lb.append(0.0)
            ub.append(1.0)
        self.set_bounds(0.0,1.0)                
        self.__dim = dim

    def _objfun_impl(self,x):
        f1,f2 = dynamic_benchmark.HE2(x,dt.tau,dt.nt,dt.taut)
        return (f1,f2,)
      
    def human_readable_extra(self):
        return "\n\t Problem dimension: " + str(self.__dim)

class HE7(base):
    """HE7 dynamic benchmark class for pygmo problem
    """
    def __init__(self, dim = 10):
        super(HE7,self).__init__(dim,0,2)
        lb = []
        ub = []
        lb.append(0.0)
        ub.append(1.0)
        for i in range(1, dim):
            lb.append(-1.0)
            ub.append(1.0)
        self.set_bounds(lb,ub)
        self.__dim = dim

    def _objfun_impl(self,x):
        f1,f2 = dynamic_benchmark.HE7(x,dt.tau,dt.nt,dt.taut)
        return (f1,f2,)
      
    def human_readable_extra(self):
        return "\n\t Problem dimension: " + str(self.__dim)


class HE9(base):
    """HE9 dynamic benchmark class for pygmo problem
    """
    def __init__(self, dim = 10):
        super(HE9,self).__init__(dim,0,2)
        self.set_bounds(0,1)
        self.__dim = dim

    def _objfun_impl(self,x):
        f1,f2 = dynamic_benchmark.HE9(x,dt.tau,dt.nt,dt.taut)
        return (f1,f2,)
      
    def human_readable_extra(self):
        return "\n\t Problem dimension: " + str(self.__dim)



class DB1a(base):
    """DB1a dynamic benchmark class for pygmo problem
    """
    def __init__(self, dim=21):
        super(DB1a, self).__init__(dim, 0, 2)
        lb = []
        ub = []
        lb.append(0.0)
        ub.append(1.0)
        for i in range(1, dim):
            lb.append(-1.0)
            ub.append(1.0)
        self.set_bounds(lb,ub)            
        self.__dim = dim

    def _objfun_impl(self, x):
        f1, f2 = dynamic_benchmark.DB1a(x, dt.t)
        return (f1, f2,)

    def human_readable_extra(self):
        return "\n\t Problem dimension: " + str(self.__dim)


class DB2a(base):
    """DB2a dynamic benchmark class for pygmo problem
    """
    def __init__(self, dim=21):
        super(DB2a, self).__init__(dim, 0, 2)
        lb = []
        ub = []
        lb.append(0.0)
        ub.append(1.0)
        for i in range(1, dim):
            lb.append(-1.0)
            ub.append(1.0)
        self.set_bounds(lb,ub)            
        self.__dim = dim

    def _objfun_impl(self, x):
        f1, f2 = dynamic_benchmark.DB2a(x, dt.t)
        return (f1, f2,)

    def human_readable_extra(self):
        return "\n\t Problem dimension: " + str(self.__dim)
        
        
class DB3a(base):
    """DB3a dynamic benchmark class for pygmo problem
    """
    def __init__(self, dim=21):
        super(DB3a, self).__init__(dim, 0, 2)
        lb = []
        ub = []
        lb.append(0.0)
        ub.append(1.0)
        for i in range(1, dim):
            lb.append(-1.0)
            ub.append(1.0)
        self.set_bounds(lb,ub)
        self.__dim = dim

    def _objfun_impl(self, x):
        f1, f2 = dynamic_benchmark.DB3a(x, dt.t)
        return (f1, f2,)

    def human_readable_extra(self):
        return "\n\t Problem dimension: " + str(self.__dim)        
        
class DB4a(base):
    """DB4a dynamic benchmark class for pygmo problem
    """
    def __init__(self, dim=21):
        super(DB4a, self).__init__(dim, 0, 2)
        lb = []
        ub = []
        lb.append(0.0)
        ub.append(1.0)
        for i in range(1, dim):
            lb.append(-1.0)
            ub.append(1.0)
        self.set_bounds(lb,ub)
        self.__dim = dim

    def _objfun_impl(self, x):
        f1, f2 = dynamic_benchmark.DB4a(x, dt.t)
        return (f1, f2,)

    def human_readable_extra(self):
        return "\n\t Problem dimension: " + str(self.__dim)

class DB5a(base):
    """DB5a dynamic benchmark class for pygmo problem
    """
    def __init__(self, dim=21):
        super(DB5a, self).__init__(dim, 0, 2)
        lb = []
        ub = []
        lb.append(0.0)
        ub.append(1.0)
        for i in range(1, dim):
            lb.append(-1.0)
            ub.append(1.0)
        self.set_bounds(lb,ub)
        self.__dim = dim

    def _objfun_impl(self, x):
        f1, f2 = dynamic_benchmark.DB5a(x, dt.t)
        return (f1, f2,)

    def human_readable_extra(self):
        return "\n\t Problem dimension: " + str(self.__dim)

class DB6a(base):
    """DB6a dynamic benchmark class for pygmo problem
    """
    def __init__(self, dim=21):
        super(DB6a, self).__init__(dim, 0, 2)
        lb = []
        ub = []
        lb.append(0.0)
        ub.append(1.0)
        for i in range(1, dim):
            lb.append(-1.0)
            ub.append(1.0)
        self.set_bounds(lb,ub)
        self.__dim = dim

    def _objfun_impl(self, x):
        f1, f2 = dynamic_benchmark.DB6a(x, dt.t)
        return (f1, f2,)

    def human_readable_extra(self):
        return "\n\t Problem dimension: " + str(self.__dim)

class DB7a(base):
    """DB7a dynamic benchmark class for pygmo problem
    """
    def __init__(self, dim=21):
        super(DB7a, self).__init__(dim, 0, 2)
        lb = []
        ub = []
        lb.append(0.0)
        ub.append(1.0)
        for i in range(1, dim):
            lb.append(-1.0)
            ub.append(1.0)
        self.set_bounds(lb,ub)
        self.__dim = dim

    def _objfun_impl(self, x):
        f1, f2 = dynamic_benchmark.DB7a(x, dt.t)
        return (f1, f2,)

    def human_readable_extra(self):
        return "\n\t Problem dimension: " + str(self.__dim)

class DB8a(base):
    """DB8a dynamic benchmark class for pygmo problem
    """
    def __init__(self, dim=21):
        super(DB8a, self).__init__(dim, 0, 2)
        lb = []
        ub = []
        lb.append(0.0)
        ub.append(1.0)
        for i in range(1, dim):
            lb.append(-1.0)
            ub.append(1.0)
        self.set_bounds(lb,ub)
        self.__dim = dim

    def _objfun_impl(self, x):
        f1, f2 = dynamic_benchmark.DB8a(x, dt.t)
        return (f1, f2,)

    def human_readable_extra(self):
        return "\n\t Problem dimension: " + str(self.__dim)


