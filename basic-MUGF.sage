# In[ ]:


import numpy as np
import math
np.set_printoptions(precision=4)
import time

lam_i = [1/(3600*175)] * 8
mu_i = [1/(60*30)]* 8

lam_v = 1/(2654*3600) #1.047 * 10^(-7)
mu_v = 1/(60*100) #1.667 * 10^(-4)
lam_h = 1/(60000*3600) #4.630 * 10^(-9)
mu_h = 1/(8*3600) #3.472 * 10^(-5)


# In[ ]:


def calculate_availability(server):
    L = len(server)
    M = len(server[0])
    
    ul = [0]*L
    
    var("z1 z2 z3 z4 z5 z6 z7 z8 z0")

    zl = [z1,z2,z3,z4,z5,z6,z7,z8, z0]
    
    t=0
    
    for serveri in server:
        N = int(np.prod(np.array(serveri)+1)+2)

        Q = [[0] * N for i in range(N)]
        g = [[0] * M for i in range(N)]

        n_list = [int(np.prod(np.array(serveri[j:])+1)) for j in range(M)]
        n_list.append(1)


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
                power = i

                for j in range(M):
                    if i%n_list[j] >= n_list[j+1]:
                        Q[i][i-n_list[j+1]] += (i%n_list[j])//n_list[j+1]*lam_i[j]
                        Q[i][i] += -(((i%n_list[j])//n_list[j+1])*lam_i[j])

                    if i%n_list[j] < n_list[j+1]*serveri[j]:
                        Q[i][i+n_list[j+1]] += (serveri[j]-(i%n_list[j])//n_list[j+1])*mu_i[j]
                        Q[i][i] += -(serveri[j]-(i%n_list[j])//n_list[j+1])*mu_i[j]

                    g[i][j] = g_list[j]* (power // n_list[j+1])
                    power %= n_list[j+1]  

        Q = [column[:N-1]+[column[N-1]+1] for column in Q]
        
        Q1 = matrix(np.array(Q))

        y = vector([0]*(N-1)+[1])
        
        p = list(Q1.solve_left(y))

        u=0

        for i in range(N):
            z = 1
            for j in range(M):
                z *= zl[j]^g[i][j]
            u += p[i]*z
        ul[t]=u  
        t+=1

    u_m = 1
    
    for u in ul:
        u_m *= u
        
    pg_last = []
    
    expansion(u_m, M, [0]*M, zl, M, pg_last)
    
    availability = 0
    for i in pg_last:
        if min(i[1:]) >= w:
            availability += i[0]
    
            
    return availability


# In[ ]:


def expansion(u, m, l_g, zl, M, pg_last):
    u_ = u.coefficients(zl[M-m])
    for a in u_:
        if m > 0:
            l_g[M-m] = a[1]
            expansion(a[0], m-1, l_g, zl, M, pg_last)
        else:
            pg_last.append([a[0]]+l_g)


# In[ ]:


w = 10000

g_list =[2500,2500,2500,2500,2500,2500,2500,2500]

l_servers = [
    [[1,1,1,1], [1,1,1,1], [1,1,1,1], [1,1,1,1]],
    [[2,2,0,0], [0,2,2,0], [0,0,2,2], [2,0,0,2]],
    [[1,1,1,1,1], [1,1,1,1,1], [1,1,1,1,1], [1,1,1,1,1]],
    [[1,1,1,1,0], [1,1,1,0,1], [1,1,0,1,1],[1,0,1,1,1], [0,1,1,1,1]],
    [[1,1,1,1,0], [0,1,1,1,1], [1,0,1,1,1],[1,1,0,1,1], 
     [1,1,1,0,1], [2,1,1,0,0], [0,2,1,1,0], [0,0,2,1,1], 
     [1,0,0,2,1], [1,1,0,0,2]]
]


for servers in l_servers:
    l = []
    for _ in range(10):
        start = time.process_time()

        calculate_availability(servers)

        t = time.process_time()-start
        l.append(t)

    print(sum(l)/len(l))
