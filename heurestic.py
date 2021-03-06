import math
import random
import cplex
from cplex.exceptions import CplexSolverError
import itertools
import pandas as pd


def heurestic(param):
    
        
    # Parameters

    N = param.N 
    P = param.P
    L = sef.L
    K = param.K 
    H = param.H 
    ql = param.ql 
    qk = param.qk 
    rk = param.rk
    sk = param.sk 
    rh = param.rh 
    sh = param.sh 
    V = param.V 
    S = param.S 
    W = param.W 
    A = param.A 
    Com = param.Com 
    B = param.B
    h = param.distance_matrix_PN 
    a = param.distance_matrix_SN 
    demand = param.demand 
    m = param.m 
    g = param.g 
    f = param.f
    Btc = param.Btc 
    Ch = param.Ch
    Pa= param.Pa
    subtrees = param.subtrees
    uprime = param.uprime 
    lprime = param.lprime
    
    
    # Initiating variables
    
    xijkt= [[[[0 for x in range(P)] for x in range(K*H)] for x in V] for x in V]
    zijktc= [[[[[0 for x in Com] for x in range(P)] for x in range(K*H)] for x in V] for x in V]
    ditc= [[[0 for x in Com] for x in range(P)] for x in V]
    diktc=[[[[0 for x in Com] for x in range(P)] for x in range(K*H)] for x in V]
    rikt=[[[0  for x in range(P)] for x in range(K*H)] for x in V]
    
    # Approx Algorithm
    
    # The following code solves the delivery model i.e. delivery problems for secondary vehicles and finds the approximate values of the 
    # delivery variable at the optimal and finds the optimal primary objective.
    
    for t in range(P):
    
        if t == 0 :
            total_delivery = min(K*H*qk,L*ql,sum([m[i][t][c] for i in range(V[1],len(V)) for c in Com]))
        else:
            total_delivery = min(K*H*qk,L*ql,sum([m[i][t][c] for i in range(V[1],len(V)) for c in Com],(sum([Btc[t][c] for c in Com])-sum([ditc[i][t-1][c] for i in range(V[1],len(V)) for c in Com]))/N))
        mean_unsatisfied = (sum([m[i][t][c] for i in range(V[1],len(V)) for c in Com]) - total_delivery)/float(N)
        Least_int = math.ceil(mean_unsatisfied)
        Greatest_int=math.floor(mean_unsatisfied)
        Ubound = uprime*mean_unsatisfied
        Lbound=lprime*mean_unsatisfied 
        
        if Least_int > Ubound or Greatest_int < Lbound:
            print "change uprime and lprime, original model infeasible"
            for i in range(V[1],len(V)):
                for c in Com:
                    ditc[i][t][c] = m[i][t][c]-Least_int
        
        else:
            for i in range(V[1],len(V)):
                for c in Com:
                    ditc[i][t][c] = m[i][t][c]-Least_int
        
        leftover = total_delivery - sum([ditc[i][t][c] for i in range(V[1],len(V)) for c in Com]) 
        
        
        vehicles_upto_number =[0 for x in V]
        vehicles_for_number=[0 for x in V]
        demand_upto=[0 for x in V]
        demand_for= [0 for x in V]
        dedicated_vehicles= [0 for x in V]
        sharable_vehicles= [0 for x in V]
        sharing_potential= [0 for x in V]
        demand_sharing_potential_up= [0 for x in V]
        demand_sharing_potential_down= [0 for x in V]
        unshared_demand_up=[0 for x in V]
        
        vehicles_upto_number.append(0)
        vehicles_for_number.append(0)
        for i in range(len(subtrees)):
            if subtrees[i] != None :
                dem1=0
                dem2=0
                for c in Com:
                    dem2 = dem2 + ditc[i][t][c]
                for j in range(len(subtrees[i])):
                    for c in Com:
                        dem1= dem1 + ditc[subtrees[i][j]][t][c]
                demand_upto[i]= dem1 +dem2
                demand_for[i]= dem2
                vehicles_upto_number[i] = math.ceil((dem1+dem2)/qk)
                vehicles_for_number[i]= math.ceil(dem2/qk)
                dedicated_vehicles[i]=math.floor(dem2+dem1/qk)
                sharable_vehicles[i]= vehicles_upto_number[i] - dedicated_vehicles[i]
                sharing_potential[i]= vehicles_for_number[i]- dedicated_vehicles[i]
                demand_sharing_potential_down[i]= dem2 - (math.floor(dem2/qk))*qk  
                demand_sharing_potential_up[i]= dem2 + dem1 - (math.floor(dem2+ dem1/qk))*qk
                unshared_demand_up[i]= (math.floor((dem2 + dem1)/qk))*qk
            else:
                if i!=0:
                    dem2=0
                    for c in Com:
                        dem2 =dem2 + ditc[i][t][c]
                    demand_upto[i]= dem1 +dem2
                    demand_for[i]= dem2
                    vehicles_upto_number[i]= math.ceil(dem2/qk)
                    vehicles_for_number[i]= math.ceil(dem2/qk)
                    dedicated_vehicles[i]=math.floor(dem2/qk)
                    sharable_vehicles[i]= vehicles_upto_number[i] - dedicated_vehicles[i]
                    sharing_potential[i]= vehicles_for_number[i]- dedicated_vehicles[i]
                    demand_sharing_potential_up[i]= dem2  - (math.floor(dem2/qk))*qk 
                    demand_sharing_potential_down[i]= dem2  - (math.floor(dem2/qk))*qk 
                    unshared_demand_up[i]= (math.floor(dem2/qk))*qk
        
        node_taken= [0 for j in V] 
        i=1
        while leftover > 0 and sum(node_taken) < N:
            i=1
            while i < len(V) and leftover > 0:
                c=random.choice(Com)
                if math.floor((demand_upto[i]+1)/float(qk)) == math.floor((demand_upto[i])/float(qk))  and ditc[i][t][c] < m[i][t][c] :
                    ditc[i][t][c] = ditc[i][t][c] + 1
                    leftover = leftover-1
                    i=i+1
                else:
                    node_taken[i]==1
                    i=i+1
            
        
        
        distance_from_depot = [0 for x in V] 
        for i in V:
            if i==0:
                distance_from_depot[i] = 0 
            elif Pa[i]==None:
                distance_from_depot[i] = 0 
            else:
                distance_from_depot[i]= distance_matrix_tree[Pa[i]][i] + distance_from_depot[Pa[i]]
        sorted_distance= sorted(distance_from_depot)
        sorted_node=[]
     
        for i in sorted_distance:
            if i !=0:
                j=0
                while j < len(distance_from_depot) :
                    if distance_from_depot[j]==i and sorted_node.count(j)==0 :
                           sorted_node.append(j)  
                           
                           j=len(distance_from_depot)
                    else:
                        j=j+1
        i=1               
        while leftover!=0:
            c= random.choice(Com)
            if ditc[i][t][c] < m[i][t][c] :
                ditc[sorted_node[i]][t][c]= ditc[sorted_node[i]][t][c]+1
                leftover= leftover - 1
                if i == len(sorted_node):
                    i=1
                else:
                    i=i+1
        
        for i in range(len(subtrees)):
            if subtrees[i] != None :
                dem1=0
                dem2=0
                for c in Com:
                    dem2 = dem2 + ditc[i][t][c]
                for j in range(len(subtrees[i])):
                    for c in Com:
                        dem1= dem1 + ditc[subtrees[i][j]][t][c]
                demand_upto[i]= dem1 +dem2
                demand_for[i]= dem2
                vehicles_upto_number[i] = math.ceil((dem1+dem2)/qk)
                vehicles_for_number[i]= math.ceil(dem2/qk)
                dedicated_vehicles[i]=math.floor(dem2+dem1/qk)
                sharable_vehicles[i]= vehicles_upto_number[i] - dedicated_vehicles[i]
                sharing_potential[i]= vehicles_for_number[i]- dedicated_vehicles[i]
                demand_sharing_potential_down[i]= dem2 - (math.floor(dem2/qk))*qk  
                demand_sharing_potential_up[i]= dem2 + dem1 - (math.floor(dem2+ dem1/qk))*qk
                unshared_demand_up[i]= (math.floor((dem2 + dem1)/qk))*qk
            else:
                if i!=0:
                    dem2=0
                    for c in Com:
                        dem2 =dem2 + ditc[i][t][c]
                    vehicles_upto_number[i]= math.ceil(dem2/qk)
                    vehicles_for_number[i]= math.ceil(dem2/qk)
                    dedicated_vehicles[i]=math.floor(dem2/qk)
                    sharable_vehicles[i]= vehicles_upto_number[i] - dedicated_vehicles[i]
                    sharing_potential[i]= vehicles_for_number[i]- dedicated_vehicles[i]
                    demand_sharing_potential_up[i]= dem2  - (math.floor(dem2/qk))*qk 
                    demand_sharing_potential_down[i]= dem2  - (math.floor(dem2/qk))*qk 
                    unshared_demand_up[i]= (math.floor(dem2/qk))*qk
     
    
    # The folllowing code solves the secondary vehicle routing model
    # It gives  p=aprroximate solution for xijkt, zijktc, diiktc  (these values are optimal for the given ditc. However since ditc is 
    # appriximate so these values are not gloabally optimal for the model)              
    
        pos_of_vehicle=[0 for x in range(K*H)] 
        vehicle_capacity_left=[qk for x in range(K*H)]
        weight_carried = [0 for x in range(K*H)]
        vehicles_at= [[] for x in V]
        vehicles_for=[[] for x in V]
        vehicle_path=[[]for x in range(K*H)]
        demand_left=[[0 for x in Com] for x in V]
        
    
        k=0
               
        for i in xrange(len(Ch)-1,0,-1):
            if i!=0:
                if Ch[i]!= None:
                    for x in Ch[i]:
                        for c in Com:
                            demand_left[x][c] = ditc[x][t][c] 
                        a=0
                        while demand_left[x][c] >0 and a < len(vehicles_at[x]):
                            for c in Com: 
                                diktc[x][vehicles_at[x][a]][t][c] = min(demand_left[x][c],vehicle_capacity_left[vehicles_at[x][a]]) 
                                weight_carried[vehicles_at[x][a]]= weight_carried[vehicles_at[x][a]] + diktc[x][vehicles_at[x][a]][t][c]
                                vehicle_capacity_left[vehicles_at[x][a]]= vehicle_capacity_left[vehicles_at[x][a]]-diktc[x][vehicles_at[x][a]][t][c]  # no need for this
                                demand_left[x][c]= demand_left[x][c] - diktc[x][vehicles_at[x][a]][t][c]
                
                            a=a+1
                        while sum([demand_left[x][c] for c in Com]) >0:
                            (vehicles_at[x]).append(k)
                            for c in Com: 
                                diktc[x][k][t][c] = min(demand_left[x][c],vehicle_capacity_left[k]) 
                                weight_carried[k]= weight_carried[k] + diktc[x][k][t][c]
                                vehicle_capacity_left[k]= vehicle_capacity_left[k]-diktc[x][k][t][c]  # no need for this
                                demand_left[x][c]= demand_left[x][c] - diktc[x][k][t][c]
                            k=k+1
                        b=0
                        d=b+1
                        while len(vehicles_at[x]) != vehicles_upto_number[x]  and b < len(vehicles_at[x]):
                            if d < len(vehicles_at[x]) :
                                if weight_carried[vehicles_at[x][b]] + weight_carried[vehicles_at[x][c]] < qk :
                                    redundant= vehicles_at[x][d]
                                    shared= vehicles_at[x][b]
                                    for c in Com:
                                        diktc[x][shared][t][c]= diktc[x][shared][t][c] + diktc[x][redundant][t][c]
                                        weight_carried[shared]= weight_carried[shared] + weight_carried[redundant]
                                        vehicle_capacity_left[shared]= qk- weight_carried[shared]
                                        diktc[x][redundant][t][c]=0
                                        weight_carried[redundant]=0
                                        vehicle_capacity_left[redundant]=qk
                                    for j in range(len(vehicle_path[redundant])-1):
                                        xijkt[vehicle_path[redundant][j]][vehicle_path[redundant][j+1]][shared][t]=1
                                        xijkt[vehicle_path[redundant][j+1]][vehicle_path[redundant][j]][shared][t]=1
                                        xijkt[vehicle_path[redundant][j]][vehicle_path[redundant][j+1]][redundant][t]=0
                                        xijkt[vehicle_path[redundant][j+1]][vehicle_path[redundant][j]][redundant][t]=0
                                        zijktc[vehicle_path[redundant][j]][vehicle_path[redundant][j+1]][shared][t][c]= zijktc[vehicle_path[redundant][j]][vehicle_path[redundant][j+1]][redundant][t][c]
                                        zijktc[vehicle_path[redundant][j]][vehicle_path[redundant][j+1]][redundant][t][c]=0
                                        if j !=0 :
                                            diktc[j][shared][t][c]= diktc[j][redundant][t][c]
                                            diktc[j][redundant][t][c]=0
                                    b=0
                                    d=b+1
                                else:
                                    d=d+1
                            else:
                                b=b+1 
                                d=b+1
                        
                        for y in range(len(vehicles_at[x])):    
                            xijkt[x][i][vehicles_at[x][y]][t]=1
                            xijkt[i][x][vehicles_at[x][y]][t]=1
                            vehicles_at[i].append(vehicles_at[x][y])
                            for c in Com:
                                if len(vehicle_path[vehicles_at[x][y]])==0:
                                    zijktc[i][x][vehicles_at[x][y]][t][c]=  diktc[x][vehicles_at[x][y]][t][c]
                                else:
                                    zijktc[i][x][vehicles_at[x][y]][t][c]=  diktc[x][vehicles_at[x][y]][t][c] + zijktc[x][vehicle_path[vehicles_at[x][y]][0]][vehicles_at[x][y]][t][c]
                            vehicle_path[vehicles_at[x][y]].insert(0,x)
                        pop=len(vehicles_at[x])
                        for z in range(pop):
                            vehicles_at[x].pop()
                            
        for i in vehicle_path:
            if i!=None and len(i)!=0: 
                i.insert(0,Pa[i[0]])
        
        for i in range(len(Pa)):
            if i!=0:
                if Pa[i]==None:
                    x=i
                    for c in Com:
                        demand_left[x][c] = ditc[x][t][c] 
                    a=0
                    while demand_left[x][c] >0 and a < len(vehicles_at[x]):
                        for c in Com: 
                            diktc[x][vehicles_at[x][a]][t][c] = min(demand_left[x][c],vehicle_capacity_left[vehicles_at[x][a]]) 
                            weight_carried[vehicles_at[x][a]]= weight_carried[vehicles_at[x][a]] + diktc[x][vehicles_at[x][a]][t][c]
                            vehicle_capacity_left[vehicles_at[x][a]]= vehicle_capacity_left[vehicles_at[x][a]]-diktc[x][vehicles_at[x][a]][t][c]  # no need for this
                            demand_left[x][c]= demand_left[x][c] - diktc[x][vehicles_at[x][a]][t][c]
                
                        a=a+1
                        
                        
                        
                    while sum([demand_left[x][c] for c in Com]) >0:
                        if k < K*H:
                            (vehicles_at[x]).append(k)
                            for c in Com: 
                                diktc[x][k][t][c] = min(demand_left[x][c],vehicle_capacity_left[k]) 
                                weight_carried[k]= weight_carried[k] + diktc[x][k][t][c]
                                vehicle_capacity_left[k]= vehicle_capacity_left[k]-diktc[x][k][t][c]  # no need for this
                                demand_left[x][c]= demand_left[x][c] - diktc[x][k][t][c]
                            k=k+1 
                        else:
                            for j in range(len(vehicle_capacity_left)):
                                for c in Com: 
                                    diktc[x][j][t][c] = min(demand_left[x][c],vehicle_capacity_left[j]) 
                                    weight_carried[j]= weight_carried[j] + diktc[x][j][t][c]
                                    vehicle_capacity_left[j]= vehicle_capacity_left[j]-diktc[x][j][t][c]  # no need for this
                                    demand_left[x][c]= demand_left[x][c] - diktc[x][j][t][c]
                                
                    
    # In the following code we formulate the primary vehicle routing model in cplex and solve it                       
    
    model = cplex.Cplex()
    
    for i in V:
        yijlt.append([])
        for j in V:
            if S.count(i)!=0 and S.count(j)!=0 : 
                if i != j:
                    yijlt[i].append([])
                    for l in range(L):
                        if yijlt[i][j] != None :
                            yijlt[i][j].append([])
                            for t in range(P):
                                varName = "y."+str(i)+"."+str(j)+"."+str(l)+"."+str(t)
                                if len(yijlt[i][j])!= 0:                      
                                    yijlt[i][j][l].append(varName)
                                    model.variables.add(obj = [h1[i][j]], 
                                                            lb = [0.0], 
                                                            ub = [1.0],
                                                            types= [model.variables.type.integer], 
                                                            names = [varName])
                else:
                    yijlt[i].append(None)
                    
            else:
                 yijlt[i].append(None)
    
    
    for i in V:
        vijltc.append([])
        for j in V:
            if i != j:
                vijltc[i].append([])
                for l in range(L):
                    if vijltc[i][j] != None :
                        vijltc[i][j].append([])
                        for t in range(P):
                            vijltc[i][j][l].append([])
                            for c in Com:
                                varName = "v."+str(i)+"."+str(j)+"."+str(l)+"."+str(t)+"."+str(c)
                                if len(vijltc[i][j][l])!= 0: 
                                    vijltc[i][j][l][t].append(varName)
                                    model.variables.add(obj = [0.0], 
                                                    lb = [0.0], 
                                                    ub = [cplex.infinity],
                                                    types= [model.variables.type.integer], 
                                                    names = [varName])
            else:
                vijltc[i].append(None)   
                
            
    for i in V:
        Ditc.append([])
        for t in range(P):
            Ditc[i].append([])
            for c in Com:
                varName = "D."+str(i)+"."+ str(t)+"."+str(c)
                Ditc[i][t].append(varName)
                model.variables.add(obj = [0.0], 
                                        lb = [0.0], 
                                        ub = [cplex.infinity],
                                        types= [model.variables.type.integer], 
                                        names = [varName])                                    
    
                                    
    for i in V:
        Sitc.append([])
        for t in range(P):
            Sitc[i].append([])
            for c in Com:
                varName = "S."+str(i)+"."+str(t)+"."+str(c)
                Sitc[i][t].append(varName)
                model.variables.add(obj = [0.0], 
                                    lb = [0.0], 
                                    ub = [cplex.infinity], 
                                    types= [model.variables.type.integer],
                                        names = [varName]) 
                        
                        
    for j in S:
        for t in range(P):
            for l in range(L):
                X=[]
                Y=[]
                for i in S:
                    if i != j:
                        X.append(yijlt[i][j][l][t])
                        X.append(yijlt[j][i][l][t])
                        Y.append(1)
                        Y.append(-1)
                print X
                print Y
                model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=X, val=Y)],senses = ["E"], rhs = [0]);
     
    i=0
    for t in range(P):
        for l in range(L):
            X=[]
            for j in S:
                if i !=j:
                    X.append(yijlt[i][j][l][t])
            model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=X, val=[1]*len(X))],senses = ["E"], rhs = [1]);
                   
    for t in range(P):
        for l in range(L):
            for r in range(2,len(A)+1):
                for p in itertools.combinations(A,r):
                    X=[]
                    for i in p:
                        for j in p:
                            if i!=j:
                                X.append(yijlt[i][j][l][t])
                    if len(X) != 0:
                        model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=X, val=[1]*len(X))],senses = ["L"], rhs = [len(X)-1])
                            
                            
                
    
    # Constraint Set 3 
    # Capacity 
    
    for i in S:
        for l in  range(L):
            for t in range(P):
                X=[]
                for j in S:
                    if i != j:
                        for c in Com:
                            X.append(vijltc[i][j][l][t][c])
                            model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=X, val=[1]*len(X))],senses = ["L"], rhs = [ql]);
    
    
    # Constraint Set 4
    # Time constraint
    
    for l in  range(L):
        for t in range(P):
            X=[]
            Y=[]
            for i in S:
                for j in S:
                    if i != j:
                        X.append(yijlt[j][i][l][t])
                        Y.append(h1[i][j] +g[i])
            model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=X, val=Y)],senses = ["L"], rhs = [rh*sh]);
    
    
    # Constraint 17
    # Relation between movement and load variables
                            
    for i in S:
        for j in A:
            if i != j :
                for l in range(L):
                    for t in range(P):
                        X=[]
                        Y=[]
                        X.append(yijlt[i][j][l][t])
                        Y.append(1)
                        for c in Com:
                            X.append(vijltc[i][j][l][t][c])
                            Y.append(-1)
                            model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=X, val=Y)],senses = ["L"], rhs = [0]);                   
    
    # Constraint 18
    # Relation between movement and load variables
    
    for i in S:
        for j in A:
            if i != j :
                for l in range(L):
                    for t in range(P):
                        X=[]
                        Y=[]
                        X.append(yijlt[i][j][l][t])
                        Y.append(100000)
                        for c in Com:
                            X.append(vijltc[i][j][l][t][c])
                            Y.append(-1)
                            model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=X, val=Y)],senses = ["G"], rhs = [0]);   
    
    j=0
    for i in S:
        for t in range(P):
            if i !=j:
                for l in range(L):
                    X=[]
                    for c in Com:
                        X.append(vijltc[i][j][l][t][c])
                    model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=X, val=[1]*len(X))],senses = ["E"], rhs = [0]);
                            
    ##Constaint 20
    # Supply Coordination
    for i in S:
        if i !=0:
            for t in range(P):
                for c in Com:
                    X=[]
                    Y=[]
                    for j in S:
                        if i !=j:
                            for l in range(L):
                                X.append(vijltc[j][i][l][t][c])
                                X.append(vijltc[i][j][l][t][c])
                                Y.append(1)
                                Y.append(-1)
                    X.append(Ditc[i][t][c])
                    Y.append(-1)
                    model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=X, val=Y)],senses = ["E"], rhs = [0]);
                    
    ### Constraint 21
    # Supply Coordination
    for i in A:
        for t in range(1,P):
            for c in Com:
                X=[]
                Y=[]
                Z=0
                X.append(Sitc[i][t][c])
                Y.append(1)
                X.append(Sitc[i][t-1][c])
                Y.append(-1)
                for j in W:
                    if i !=j:
                        for k in range(K*H):
                            Z=Z-zijktc[i][j][k][t][c]
                Z=Z-ditc[i][t][c]
                X.append(Ditc[i][t][c])    
                Y.append(-1)
                model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=X, val=Y)],senses = ["E"], rhs = [Z]);
         
    t=0
    for i in A:
        for c in Com:
            X=[]
            Y=[]
            Z=0
            d
            X.append(Sitc[i][t][c])
            Y.append(1)
            for j in W:
                if i !=j:
                    for k in range(K*H):
                        Z = Z -zijktc[i][j][k][t][c]
            X.append(Ditc[i][t][c]) 
            Z=Z-ditc[i][t][c]
            Y.append(-1)
            model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=X, val=Y)],senses = ["E"], rhs = [Z]);
            
            
    i=0
    for t in range(1,P):
        for c in Com:
            X=[]
            Y=[]
            X.append(Sitc[0][t][c])
            Y.append(1)
            X.append(Sitc[0][t-1][c])
            Y.append(-1)
            for j in S:
                if i != j:
                    for l in range(L):
                        X.append(vijltc[i][j][l][t][c])
                        Y.append(1)
            model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=X, val=Y)],senses = ["E"], rhs = [Btc[t][c]]);
            
    # Amount stored in the helicopter depot in first period less than supply in first period
    i=0
    t=0
    for c in Com:
        X=[]
        Y=[]
        X.append(Sitc[i][t][c])
        Y.append(1)
        for j in A:
            if i != j:
                for l in range(L):
                    X.append(vijltc[i][j][l][t][c])
                    Y.append(1)
        model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=X, val=Y)],senses = ["E"], rhs = [Btc[t][c]]);
            
    
    model.parameters.timelimit.set(300)
    
    try:
        model.solve()
    except CplexSolverError, e:
        print "Exception raised during solve: " + e
    else:
        solution = model.solution 
        
        
    return solution, xijktc, zijktc, ditc, diktc, rikt