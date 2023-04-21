import numpy as np
import math
np.set_printoptions(precision=4)
import time
import copy
import random

lam_s = 1/(3600*175) #1.587 * 10^(-6)
mu_s = 1/(60*30) #5.556 * 10^(-4)

lam_i = [1/(3600*175)] * 6
mu_i = [1/(60*30)]* 6

lam_v = 1/(2654*3600) #1.047 * 10^(-7)
mu_v = 1/(60*100) #1.667 * 10^(-4)
lam_h = 1/(60000*3600) #4.630 * 10^(-9)
mu_h = 1/(8*3600) #3.472 * 10^(-5)

def calculate_availability(server):
    L = len(server)
    M = len(server[0])
    
    ul = [0]*L
    
    var("z1 z2 z3 z4 z5 z6 z0")

    zl = [z1,z2,z3,z4,z5,z6, z0]
    
    t=0
    
    for serveri in server:
        N = int(np.prod(np.array(serveri)+1)+2)
        #N = math.prod(map(lambda val: val+1,serveri))+2

        Q = [[0] * N for i in range(N)]
        g = [[0] * M for i in range(N)]
        
        #何回で循環するか
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

    
    not_already = set(range(L))
    u_m = zl[-1]^math.inf
    m_yet = set(range(M))
    
    
    for _ in range(M):
        mini = math.inf
        for m_ in m_yet:
            count = 0
            for l in not_already:
                if server[l][m_]:
                    count += 1        
            if count < mini:
                m = m_
                mini = count
        m_yet.remove(m)
        ll = copy.copy(not_already)
        for l in ll:
            if server[l][m]:
                u_m *= ul[l]
                not_already.remove(l)
        u_m = sigma(u_m,zl[m])
    
    if u_m != 0:
        u_m = sigma(u_m,zl[-1])
            
    return u_m

def calculate_cost(server):
    cost = 0
    M = len(server[0])
    for i in server:
        if i == [0] * M:
            continue
        cost += C_server + sum([i[em]*C_list[em] for em in range(M)])

    return cost

def sigma(u_l, z):
    u_return = 0
    if u_l != 0:
        coefs = u_l.coefficients(z)
        for coef in coefs:
            if coef[1] >= w:
                u_return += coef[0]

    return u_return

def random_addition(servers, n):
    L = len(servers)
    M = len(servers[0])
    
    choices = [[i,j] for i in range(L) for j in range(M)]
    
    random.shuffle(choices)
    
    count = 0
    
    for choice in choices :
        if sum([servers[choice[0]][em]*c_list[em] for em in range(M)])+c_list[choice[1]] <= c_server:
            servers[choice[0]][choice[1]] += 1
            count += 1
        if count >= n:
            break
    
    return copy.deepcopy(servers)

def remove_random(server_last, availability_last, cost_last):

    check = (availability_last >= a_goal)

    L = len(server_last)
    M = len(server_last[0])
    
    rand = 0#random.randint(0,0)
    if rand == 0 and check:
        cost_best = math.inf
        check_best = False
        for serl in range(L):
            for vnfm in range(M):
                serve = copy.deepcopy(server_last)
                if serve[serl][vnfm]:
                    serve[serl][vnfm] -= 1
                    availability = calculate_availability(serve)
                    cost = calculate_cost(serve)
                    check = availability >= a_goal
                    if check_best and (not check):
                        continue
                    if (not check_best) and check:
                        server_best = copy.deepcopy(serve)
                        cost_best = cost
                        availability_best = availability
                        check_best = True
                        continue
                    if cost_best > cost or (cost == cost_best and availability > availability_best):
                        server_best = copy.deepcopy(serve)
                        availability_best = availability
                        cost_best = cost
        return copy.deepcopy(server_best), availability_best, cost_best
    
    score_best = math.inf
    check_best = False
    for serl in range(L):
        for vnfm in range(M):
            serve = copy.deepcopy(server_last)
            if serve[serl][vnfm]:
                serve[serl][vnfm] -= 1
                availability = calculate_availability(serve)
                cost = calculate_cost(serve) 
                check = availability >= a_goal
                if check_best and not check:
                    continue
                if (not check_best) and check:
                    server_best = copy.deepcopy(serve)
                    score_best = (availability_last - availability) /(cost_last-cost)
                    availability_best = availability
                    cost_best = cost
                    check_best = True
                    continue
                #if score_best > (availability_last - availability) /(cost_last-cost):
                if availability == 0:
                    continue
                if score_best > (availability_last - availability) /(cost_last-cost):
                    server_best = copy.deepcopy(serve)
                    availability_best = availability
                    cost_best = cost
                    score_best = (availability_last - availability) /(cost_last-cost)
    #print("    Remove")
    #print("   ",server_best)
    #print("   ",availability_best)
    
    #print(server_best)
    
    return copy.deepcopy(server_best), availability_best, cost_best

def remove_best(server_last, availability_last, cost_last):
    
    L = len(server_last)
    M = len(server_last[0])
    
    score_best = math.inf
    check_best = False
    for serl in range(L):
        for vnfm in range(M):
            serve = copy.deepcopy(server_last)
            if serve[serl][vnfm]:
                serve[serl][vnfm] -= 1
                availability = calculate_availability(serve)
                cost = calculate_cost(serve) 
                check = availability >= a_goal
                if check_best and not check:
                    continue
                if (not check_best) and check:
                    server_best = copy.deepcopy(serve)
                    score_best = (availability_last - availability) /(cost_last-cost)
                    availability_best = availability
                    cost_best = cost
                    check_best = True
                    continue
                #if score_best > (availability_last - availability) /(cost_last-cost):
                if availability == 0:
                    continue
                if score_best > (availability_last - availability) /(cost_last-cost):
                    server_best = copy.deepcopy(serve)
                    availability_best = availability
                    cost_best = cost
                    score_best = (availability_last - availability) /(cost_last-cost)
    #print("    Remove")
    #print("   ",server_best)
    #print("   ",availability_best)
    
    #print(server_best)
    
    return copy.deepcopy(server_best), availability_best, cost_best

def first_place():
    
    check = True
    while check:
        servers = [[0] * M for _ in range(L)]
        for m in range(M):
            num = math.ceil(w / g_list[m])
            for _ in range(num):
                l = random.randrange(0,L)
                servers[l][m] += 1
        check = False
        for l in servers:
            if c_server < sum([l[m]*c_list[m] for m in range(M)]):
                check = True
                break
    
    return servers

def random_greedy(servers,T,N, remove_method):
    L = len(servers)
    M = len(servers[0])
    
    best_cost = math.inf
    best_availability = 0
    best_servers = []
    
    count = 0
    while count < T:
        servers = random_addition(servers, N)
        
        availability = calculate_availability(servers)
        
        #print(count,"times")
        #print(" ",servers)
        #print(" ",availability)
        
        count += 1

        if availability == 0:
            continue
        if availability >= a_goal:
            cost = calculate_cost(servers)
            while availability >= a_goal:
                if cost < best_cost or (cost == best_cost and availability > best_availability):
                    best_cost = cost
                    best_availability = availability
                    best_servers = copy.deepcopy(servers)
                servers, availability, cost = remove_method(servers, availability, cost)
        else:
            num = random.randint(max(0,N-2),N)
            for _ in range(num):
                servers, _, __ = remove_method(servers, calculate_availability(servers), calculate_cost(servers))
        
        #%memit
        
        if (count) % 5 == 0:         
            l_result[count / 5 + 1] = [count,best_cost,best_availability,time.time()-start,best_servers]
            print([count,best_cost,best_availability,time.time()-start,best_servers])
    
    return best_servers, best_availability, best_cost

l_N = [1,3,5,10]
l_nines = [2,3,4,5,6]

l_remove_method = [remove_random,remove_best]

l_C_server = [0,0.1,0.5,1,5,10]


T=100



case = [1,2,3]

l_M = [4,6,4]
l_L = [11,12,15]

l_c_server = [12,12,24]

l_g_list = [[1800,2500,3000,5000],[1800,2500,3000,3500,4000,6000],[1800,2500,3000,5000]]

l_C_list = [[0.2,0.4,0.5,0.8],[0.2,0.3,0.5,0.5,0.6,0.8],[0.2,0.4,0.5,0.8]]

l_c_list = [[2,3,3,6],[2,3,3,4,4,6],[2,3,3,6]]


l_w = [10000,10000,20000]


search_list = []

for N in l_N:
    search_list.append([1,N,5,5,remove_random,"Senario-1"])

for nines in l_nines:
    search_list.append([1,5,nines,5,remove_random,"Senario-2"])

for C_S in l_C_server:
    search_list.append([1,5,5,C_S,remove_random,"Senario-3"])

for N ,nines in [(N,nines) for N in [5,10] for nines in [5]]:
    search_list.append([2,N,nines,5,remove_random,"Senario-4"])

for N ,nines in [(N,nines) for N in [5,10] for nines in [5]]:
    search_list.append([3,N,nines,5,remove_random,"Senario-5"])


import itertools
import csv
import os

os.makedirs("./result/", exist_ok=True)

search_list.pop(0)
search_list.pop(0)

for case,N,nines,C_server, remove_method, senario in search_list:

    a_goal = 1 - 10^(- nines)
    export_file = f"case{case}_N{N}_nines{nines}_server{C_server}_{remove_method.__name__}.tsv"
    
    l_result = [0] * int(T/5 +2)
    l_result[0] = ["T","cost","availability","time","arrangement"]
    
    w = l_w[case-1]
    M = l_M[case-1]
    L = l_L[case-1]
    c_list = l_c_list[case-1]
    g_list = l_g_list[case-1]
    C_list = l_C_list[case-1]
    c_server = l_c_server[case-1]

    random.seed(813)

    start = time.time()

    servers = first_place()

    l_result[1] = ["Initial",calculate_cost(servers),calculate_availability(servers),time.time()-start,servers]

    servers, availability, sum_cost = random_greedy(servers,T,N, remove_method)
    
    os.makedirs(f"./result/{senario}", exist_ok=True)

    with open(f"./result/{senario}/"+export_file,"w") as f:
        tsv = csv.writer(f,delimiter="\t")
        tsv.writerows(l_result)
