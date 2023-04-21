# In[ ]:


import numpy as np
np.set_printoptions(precision=4)

import time
import copy

lam_s = 1/(3600*175) #1.587 * 10^(-6)
mu_s = 1/(60*30) #5.556 * 10^(-4)

lam_v = 1/(2654*3600) #1.047 * 10^(-7)
mu_v = 1/(60*100) #1.667 * 10^(-4)
lam_h = 1/(60000*3600) #4.630 * 10^(-9)
mu_h = 1/(8*3600) #3.472 * 10^(-5)


# In[ ]:


def calculate_availability(servers):
    M = len(servers)
    
    ul = [0]*M
    
    var("z")

    for m in range(M):
        N = k_list[m]+3

        Q = [[0] * N for i in range(N)]
        g = [0] * N
        
        for i in range(-2,N-2):
            if i == -2:
                Q[i][i] += -mu_h
                Q[i][-3] += mu_h
            elif i == -1:
                Q[i][i] += -mu_v
                Q[i][-3] += mu_v
                Q[i][-2] += lam_h
                Q[i][i] += -lam_h
            else:
                Q[i][-2] += lam_h
                Q[i][i] += -lam_h
                Q[i][-1] += lam_v
                Q[i][i] += -lam_v

                if i < k_list[m]:
                    Q[i][i+1] += (k_list[m]-i)*mu_s
                    Q[i][i] += -(k_list[m]-i)*mu_s
                if i > 0:
                    Q[i][i-1] += i*lam_s
                    Q[i][i] += - i*lam_s

                    g[i] = g_list[m]* i
        
        l_1 = [[1]]*N
        
        Q1 = matrix(np.append(np.array(Q),l_1,axis=1))

        y = vector([0]*(N)+[1])
        
        p = list(Q1.solve_left(y))
        
        um = 0
        
        for gj, pj in zip(g,p):
            um += pj * z ^ gj
        
        ul[m] = pi([um]*servers[m])
    
    u = sigma(ul)
    
    u_m = u.coefficients(z)
    
    availability = 0
    for i in u_m:
        if i[1] >= w:
            availability += i[0]
            
    return availability


# In[ ]:


def calculate_cost(servers):
    M = len(servers)
    return sum([C_list[m]*servers[m] for m in range(M)])


# In[ ]:


def pi(l_u):
    u = z^0
    for f in l_u:
        u = u * f
    return u

def sigma(list_u):
    first = True
    for f in list_u:
        if first:
            u = f
            first = False
        else:
            l_u = u.coefficients(z)
            l_f = f.coefficients(z)
            u = 0
            for i in l_u:
                for j in l_f:
                    u += i[0]*j[0]*z^min([i[1],j[1]])
    return u


# In[ ]:


def all_search(n_max,n_parameter,servers_last, best_cost, server_best):
    servers = copy.copy(servers_last)
    
    for i in range(n_max):
        servers[n_parameter-1] = i+1
        if n_parameter:
            best_cost, server_best = all_search(n_max,n_parameter-1, servers, best_cost, server_best)
        else:
            availability = calculate_availability(servers)
            if availability >= a_goal:
                
                cost = calculate_cost(servers)
                if cost < best_cost:
                    server_best = copy.copy(servers)
                    best_cost = cost
        
    return best_cost, copy.copy(server_best)


# In[ ]:


l_max_n = [4,5]


l_M = [4,6,8]

l_k_list = [[6,4,4,2],
            [6,4,4,3,3,2],
            [6,6,4,4,3,3,2,2]]
l_g_list = [[1800,2500,3000,5000],
            [1800,2500,3000,3500,4000,6000],
            [1700,1800,2500,3000,3500,4000,5000,6000]]
l_C_I_list = [[0.2,0.4,0.5,0.8],
              [0.2,0.3,0.5,0.5,0.6,0.8],
              [0.2,0.3,0.4,0.5,0.5,0.6,0.6,0.8]]

l_w = [10000,10000,10000]

search_list = []


for nines in [2,3,4,5,6]:
    search_list.append([1,nines,5])

for C_S in [0,0.1,1]:
    search_list.append([1,5,C_S])

search_list.append([2,5,5])
search_list.append([3,5,5])

import csv
for max_n in l_max_n:
    for case, nines, C_S  in search_list:

        a_goal = 1-10^(-nines)

        M = l_M[case-1]
        k_list = l_k_list[case-1]
        g_list = l_g_list[case-1]
        C_I_list = l_C_I_list[case-1]
        w = l_w[case-1]

        C_list = [C_S + C_I_list[m]*k for m in range(M)]


        start = time.process_time()

        servers = [0]*M
        server_best = [0]*M
        best_cost = math.inf

        best_cost, server_best = all_search(max_n,M,servers,best_cost,server_best)

        print(f"max_n={max_n}, case={case}, nines={nines}, C_S={C_S}")
        print("placement="+str(server_best))
        print("availability="+str(calculate_availability(server_best)))
        print("cost="+str(best_cost))
        print(time.process_time()-start)