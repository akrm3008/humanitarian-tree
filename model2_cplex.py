# In the second model we take the advantage of the fact that the secondary model is almost tree and can be approximated
# a tree. Hence, in this code we  first covert the secondary network to a minimum spanning tree
# Next, we build the second model in which the routing constraints have been chnaged and reduced as they take advantage 
# of the tree structure of the secondary network 
# The model is coded in Cplex and used to solve the case study 
# Reduction in computational time is expected


# importing the modules 
import cplex
from cplex.exceptions import CplexSolverError
import itertools
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
import pandas as pd



# Create a new (empty) model and populate it below.
model = cplex.Cplex()

model.objective.set_sense(model.objective.sense.maximize) 
                                                                                                   
# Initialising list to store variables

xijkt= []
yijlt= []
zijktc= []
vijltc= []
diktc= []
ditc=[]
Ditc= []
Sitc=[]
rikt=[]

# Getting Parameters

N=117
P=6
L=6
K=100
H=7
ql=1000
qk=30
rk=96
sk=2
rh=8 
sh=100
V=[x for x in range(117)]
S = [0,1,15,22,44,50,78,99]
W=[x for x in range(117) if x!=0 and x!=1 and x!=15 and x!=22 and x!=44 and x!=50 and x!=78 and x!=90]
A=[1,15,22,44,50,78,99]
Com=[0,1]
B= [x for x in range(117) if x!=0]






# Distance Matrix for primary network

df = pd.read_excel("C:/Users/S K Khare/Desktop/important stuff to save/rescue/Nepal/Remote Acess Operations/models being used now/test_data3.xlsx",sheetname="distance_matrix_HLZ",na_values=None)
np_array1 = df.as_matrix()
distance_matrix_HLZ=np_array1.tolist()
for i in range(len(distance_matrix_HLZ)):
    for j in range(len(distance_matrix_HLZ[i])):
        distance_matrix_HLZ[i][j]= int(distance_matrix_HLZ[i][j])
for i in range(len(distance_matrix_HLZ )):
    for j in range(len(distance_matrix_HLZ[i])):
        if distance_matrix_HLZ[i][j]==0:
            distance_matrix_HLZ[i][j]=None
        if i>j:
            distance_matrix_HLZ[i][j]=distance_matrix_HLZ[j][i]


# Distance matrix for secondary network
            
df = pd.read_excel("C:/Users/S K Khare/Desktop/important stuff to save/rescue/Nepal/Remote Acess Operations/models being used now/test_data3.xlsx",sheetname="distance_matrix_tree",na_values=None)
np_array2 = df.as_matrix()


# Getting the Minimum Spanning Tree of the secondary network using Krushkal's MST algorihm.

X = csr_matrix(np_array2)
Tcsr = minimum_spanning_tree(X)
distance_matrix_tree= Tcsr.toarray().astype(int)
distance_matrix_tree= distance_matrix_tree.tolist()


for i in range(len(distance_matrix_tree)):
    for j in range(len(distance_matrix_tree[i])):
        distance_matrix_tree[i][j]= int(distance_matrix_tree[i][j])
for i in range(len(distance_matrix_tree)):
    for j in range(len(distance_matrix_tree[i])):
        if distance_matrix_tree[i][j]==0:
            distance_matrix_tree[i][j]=None
        if i>j:
            distance_matrix_tree[i][j]=distance_matrix_tree[j][i]



#Getting Paramters
            
Ch= [None for x in V]
j=1
i=0
a=1
order_of_search=[]
while len(order_of_search) < N-1:
    if i!=0 and max(order_of_search)<N:
        a=max(order_of_search) + 1 
    order_of_search.append(a)
    while i < len(order_of_search) or i==0 :
        j=1
        b=order_of_search[i]
        while j < len(distance_matrix_tree[b]):
            if distance_matrix_tree[b][j] != None:
                if Ch[b]== None:
                    if Ch[j]!= None:
                        if Ch[j].count(b)==0:
                            Ch[b]=[]
                            Ch[b].append(j)
                            order_of_search.append(j)
                    else:
                        Ch[b]=[]
                        Ch[b].append(j)
                        order_of_search.append(j)
                
                elif Ch[b].count(j)==0:
                    if Ch[j]!= None:
                        if Ch[j].count(b)==0:
                            Ch[b].append(j)
                            order_of_search.append(j)
                    else:
                        Ch[b].append(j)
                        order_of_search.append(j) 
            j=j+1
        i=i+1   

for i in range(len(distance_matrix_tree)):
    for j in range(len(distance_matrix_tree[i])):
        if distance_matrix_tree[i][j]==None:
            distance_matrix_tree[i][j]=1000000        
                       

# Getting Parent nodes
Pa= [None for x in V]
for i in range(len(Ch)):
    if Ch[i] !=None:
        for j in Ch[i]:
            Pa[j]=i

# Getting subtrees         
subtrees=[None for x in V]
for i in range(len(Ch)):
    if Ch[i]!= None:
        j=0
        subtrees[i]=Ch[i]
        while j<len(subtrees[i]):
            if Ch[subtrees[i][j]] != None:
                subtrees[i]= subtrees[i]+ Ch[subtrees[i][j]]
            j=j+1
            


# Getting final distance/cost matrices

a= distance_matrix_tree
h=distance_matrix_HLZ   

        
demand_in_period=[]
for t in range(P):       
    df = pd.read_excel("C:/Users/S K Khare/Desktop/important stuff to save/rescue/Nepal/Remote Acess Operations/models being used now/test_data3.xlsx",sheetname="demand_period_"+ str(t),na_values=None)
    np_array = df.as_matrix()
    demand_in_period.append(np_array.tolist())
for t in range(P):
    for i in V: 
        for c in Com:
            demand_in_period[t][i][c]= int(demand_in_period[t][i][c])
            if demand_in_period[t][i][c]==0:
                demand_in_period[t][i][c]=None
       
m=[[[None for x in Com]for x in range(P)]for x in V]
for i in V:
    for c in Com:
        for t in range(P):
            m[i][t][c]=demand_in_period[t][i][c]

# Service time for primary vehicles at each node it delivers to 
g=[]
for i in V:
    g.append(0.5) 

            
# Service time  for secondary vehicles at each node it delivers to 

f=[]
for i in V:
    f.append(4) 

    
Btc=[]
for t in range(P):
    Btc.append([])
    for c in Com:
         Btc[t].append(500)

Ic=[]
for c in Com:
    Ic.append(500)

uprime = 1.4 
lprime = 0.6


# Initialising decision variables
# yijlt      
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
                                model.variables.add(obj = [-0.01*h[i][j]], 
                                                        lb = [0.0], 
                                                        ub = [1.0],
                                                        types= [model.variables.type.integer], 
                                                        names = [varName])
            else:
                yijlt[i].append(None)
                
        else:
             yijlt[i].append(None)


# xijkt
for i in V:
    xijkt.append([])
    for j in V:
        if i != j:
            xijkt[i].append([])
            for k in range(K*H):
                if xijkt[i][j] != None :
                    xijkt[i][j].append([])
                    for t in range(P):
                        varName = "x."+str(i)+"."+str(j)+"."+str(k)+"."+str(t)
                        if len(xijkt[i][j])!= 0: 
                            xijkt[i][j][k].append(varName)
                            model.variables.add(obj = [-0.01*a[i][j]], 
                                                lb = [0.0], 
                                                ub = [1.0], 
                                                types= [model.variables.type.integer],
                                                names = [varName])
        else:
             xijkt[i].append(None)

# zijkt
for i in V:
    zijktc.append([])
    for j in V:
        if i != j:
            zijktc[i].append([])
            for k in range(K*H):
                if zijktc[i][j] != None :
                    zijktc[i][j].append([])
                    for t in range(P):
                        zijktc[i][j][k].append([])
                        for c in Com:
                            varName = "z."+str(i)+"."+str(j)+"."+str(k)+"."+str(t)+"."+str(c)
                            if len(zijktc[i][j][k])!= 0:                      
                                zijktc[i][j][k][t].append(varName)
                                model.variables.add(obj = [0.0], 
                                                lb = [0.0], 
                                                ub = [cplex.infinity], 
                                                types= [model.variables.type.integer],
                                                names = [varName])
        else:
             zijktc[i].append(None)    

# vijltc    
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
              
# ditc
for i in V:
    if i!=0:
        ditc.append([])
        for t in range(P):
            ditc[i].append([])
            for c in Com:
                varName = "d."+str(i)+"."+str(t)+"."+str(c)
                ditc[i][t].append(varName)
                model.variables.add(obj = [1.0], 
                                    lb = [0.0], 
                                    ub = [cplex.infinity],
                                    types= [model.variables.type.integer], 
                                    names = [varName])
    else:
        ditc.append(None)                  

# diktc                            
for i in V:
    if i!=0:
        diktc.append([])
        for k in range(K*H):
            diktc[i].append([])
            for t in range(P):
                diktc[i][k].append([])
                for c in Com:
                    varName = "d."+str(i)+"."+str(k)+"."+str(t)+"."+str(c)
                    diktc[i][k][t].append(varName)
                    model.variables.add(obj = [0.0], 
                                        lb = [0.0], 
                                        ub = [cplex.infinity],
                                        types= [model.variables.type.integer], 
                                        names = [varName])                                
    else:
        diktc.append(None)
  
# Ditc                           
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

# Sitc                               
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

# rikt         
for i in V:
    rikt.append([])
    for k in range(K*H):
        rikt[i].append([])
        for t in range(P):
            varName = "r."+str(i)+"."+str(k)+"."+str(t)
            rikt[i][k].append(varName)
            model.variables.add(obj = [0.0], 
                                lb = [0.0], 
                                ub = [1.0], 
                                types= [model.variables.type.integer],
                                names = [varName]) 
            
        
        



            
# Constraints for cycling for the primary vehicles

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
            model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=X, val=Y)],senses = ["E"], rhs = [0]);

# For a primary vehicle coming into a secondary vehicle in a period

i=0
for t in range(P):
    for l in range(L):
        X=[]
        for j in S:
            if i !=j:
                X.append(yijlt[i][j][l][t])
        model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=X, val=[1]*len(X))],senses = ["E"], rhs = [1]);



# Subtours elimination for primary vehicle 

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
                        
                        
            

# Capacity constraint for primary vehicles

for i in S:
    for l in  range(L):
        for t in range(P):
            X=[]
            for j in S:
                if i != j:
                    for c in Com:
                        X.append(vijltc[i][j][l][t][c])
                    model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=X, val=[1]*len(X))],senses = ["L"], rhs = [ql]);

# Time constraint for primary vehicles

for l in  range(L):
    for t in range(P):
        X=[]
        Y=[]
        for i in S:
            for j in S:
                if i != j:
                    X.append(yijlt[j][i][l][t])
                    Y.append(h[i][j] +g[i])
        model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=X, val=Y)],senses = ["L"], rhs = [rh*sh]);





# Supply constraint for secondary vehicles

for i in V:
    for k in  range(K*H):
        for t in range(P):
            X=[]
            for j in V:
                if i != j:
                    for c in Com:
                        X.append(zijktc[i][j][k][t][c])
                    model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=X, val=[1]*len(X))],senses = ["L"], rhs = [qk]);

for k in range(K*H):
    for t in range(P):
        X=[]
        for i in range(V[1],len(V)):
            for c in Com:
                X.append(diktc[i][k][t][c])
        model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=X, val=[1]*len(X))],senses = ["L"], rhs = [qk]);






# Time constraint for secondary vehicles 

for k in  range(K*H):
    for t in range(P):
        X=[]
        Y=[]
        for i in V:
            for j in V:
                if i != j:
                    X.append(xijkt[j][i][l][t])
                    Y.append(a[i][j] + f[i])
        model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=X, val=Y)],senses = ["L"], rhs = [rk/sk]);


# No secondary vehicle going from or to primary depot
i=0 
for j in V:
    if i!=j:
        for t in range(P):
            for k in range(K*H):
                X=[]
                X.append(xijkt[i][j][k][t])
                model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=X, val=[1]*len(X))],senses = ["E"], rhs = [0]);

# Links from parents to children

for i in W:
    if Ch[i] != None:
        for j in Ch[i]:
            for k in range(K*H):
                for t in range(P):
                    X=[]
                    Y=[]
                    X.append(xijkt[Pa[i]][i][k][t])
                    Y.append(1)
                    X.append(xijkt[i][j][k][t])
                    Y.append(-1)
                    model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=X, val=Y)],senses = ["G"], rhs = [0]);
                    
for i in B:
    for j in B:
        if i!=j:
            if xijkt[i][j][k][t]!=None:
                for k in range(K*H):
                    for t in range(P):
                        X=[]
                        Y=[]
                        X.append(xijkt[i][j][k][t])
                        Y.append(1)
                        X.append(xijkt[j][i][k][t])
                        Y.append(-1)
                        model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=X, val=Y)],senses = ["E"], rhs = [0]);

for i in W:
    for k in range(K*H):
        for t in range(P):
            X=[]
            Y=[]
            X.append(xijkt[Pa[i]][i][k][t])
            Y.append(1)
            X.append(rikt[i][k][t])
            Y.append(-1)
            model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=X, val=Y)],senses = ["G"], rhs = [0]);
        

# Relation between movement and load variables
        
for i in V:
    for j in W:
        if i != j :
            for k in range(K*H):
                for t in range(P):
                    X=[]
                    Y=[]
                    X.append(xijkt[i][j][k][t])
                    Y.append(1)
                    for c in Com:
                        X.append(zijktc[i][j][k][t][c])
                        Y.append(-1)
                    model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=X, val=Y)],senses = ["L"], rhs = [0]);  


# Relation between movement and load variables

for i in V:
    for j in W:
        if i != j :
            for k in range(K*H):
                for t in range(P):
                    X=[]
                    Y=[]
                    X.append(xijkt[i][j][k][t])
                    Y.append(100000)
                    for c in Com:
                        X.append(zijktc[i][j][k][t][c])
                        Y.append(-1)
                    model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=X, val=Y)],senses = ["G"], rhs = [0]);  

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



# A secondary will deliver to a node in a period only if the node is alloted to it


for i in W:
    for k in range(K*H):
        for t in range(P):
            for c in Com:
                X=[]
                Y=[]
                X.append(diktc[i][k][t][c])
                Y.append(1)
            X.append(rikt[i][k][t])
            Y.append(-1)
            model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=X, val=Y)],senses = ["G"], rhs = [0]); 

# A secondary will deliver to a node in a period only if the node is alloted to it

for i in W:
    for k in range(K*H):
        for t in range(P):
            for c in Com:
                X=[]
                Y=[]
                X.append(diktc[i][k][t][c])
                Y.append(1)
            X.append(rikt[i][k][t])
            Y.append(-1000000)
            model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=X, val=Y)],senses = ["L"], rhs = [0]);                 


     
for i in W:
    for t in range(P):
        for c in Com:
            for k in range(K*H):
                X=[]
                Y=[]
                for j in range(V[1],len(V)):
                    if i !=j:
                        X.append(zijktc[j][i][k][t][c])
                        X.append(zijktc[i][j][k][t][c])
                        Y.append(1)
                        Y.append(-1)
                X.append(diktc[i][k][t][c])
                Y.append(-1)
                model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=X, val=Y)],senses = ["E"], rhs = [0]);
              
               
# Coordination between primaary and secondary network 
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
                

# Coordination between primary and secodary network 
for i in A:
    for t in range(1,P):
        for c in Com:
            X=[]
            Y=[]
            X.append(Sitc[i][t][c])
            Y.append(1)
            X.append(Sitc[i][t-1][c])
            Y.append(-1)
            for j in W:
                if i !=j:
                    for k in range(K*H):
                        X.append(zijktc[i][j][k][t][c])
                        Y.append(1)
            X.append(ditc[i][t][c])
            Y.append(1)
            X.append(Ditc[i][t][c])    
            Y.append(-1)
            model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=X, val=Y)],senses = ["E"], rhs = [0]);
     
t=0
for i in A:
    for c in Com:
        X=[]
        Y=[]
        X.append(Sitc[i][t][c])
        Y.append(1)
        for j in W:
            if i !=j:
                for k in range(K*H):
                    X.append(zijktc[i][j][k][t][c])
                    Y.append(1)
        X.append(Ditc[i][t][c])
        X.append(ditc[i][t][c])    
        Y.append(-1)
        Y.append(1)
        model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=X, val=Y)],senses = ["E"], rhs = [0]);
        
        


i=0
for t in range(1,P):
    for c in Com:
        X=[]
        Y=[]
        X.append(Sitc[i][t][c])
        Y.append(1)
        X.append(Sitc[i][t-1][c])
        Y.append(-1)
        for j in A:
            if i != j:
                for l in range(L):
                    X.append(vijltc[i][j][l][t][c])
                    Y.append(1)
        model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=X, val=Y)],senses = ["E"], rhs = [Btc[t][c]]);  
        
# Amount stored in the primary depot in first period less than supply in first period
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
    
# Summming up all the deliveries at a node
for i in B:
    for t in range(P):
        for c in Com:   
            X=[]
            Y=[]
            for k in range(K*H):
                X.append(diktc[i][k][t][c])
                Y.append(1)
            X.append(ditc[i][t][c])
            Y.append(-1)
            model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=X, val=Y)],senses = ["E"], rhs = [0]);

            
                

# Vehicles entering secondary depots have load variable equal to zero
for i in W:
    for t in range(P):
        for j in S:
            if i !=j:
                for k in range(K*H):
                    X=[]
                    for c in Com:
                        X.append(zijktc[i][j][k][t][c])
                        model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=X, val=[1]*len(X))],senses = ["E"], rhs = [0]);
                    

# Vehicle entering primary depots, have load variable equal to zero
j=0
for i in S:
    for t in range(P):
        if i !=j:
            for l in range(L):
                X=[]
                for c in Com:
                    X.append(vijltc[i][j][l][t][c])
                model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=X, val=[1]*len(X))],senses = ["E"], rhs = [0]);
                        



# Demand at primary depot is zero
i=0
for t in range(P):
    X=[]
    for c in Com:
        X.append(Ditc[i][t][c])
    model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=X, val=[1]*len(X))],senses = ["E"], rhs = [0]);
                   

# Total demand 
               
for i in range(V[1],len(V)):
    for c in Com:
        for t in range(P):
            X=[]
            X.append(ditc[i][t][c])
            model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=X, val=[1]*len(X))],senses = ["L"], rhs = [m[i][t][c]]);





## One secondary vehicle serving one tree in a period
for k in range(K*H):
    for t in range(P):
         X=[]
         for i in A:
            X.append(rikt[i][k][t])
         model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=X, val=[1]*len(X))],senses = ["L"], rhs = [1]);


# Constraint to transfer a secondary vehicle via primary vehicles 
for k in range(K*H):
    for t in range(P-1):
        for i in A:
            for j in A: 
                if i!=j:
                    for l in range(L):
                        X=[]
                        Y=[]
                        X.append(yijlt[i][j][l][t+1])
                        Y.append(-1)
                    print t
                    X.append(rikt[i][k][t])
                    X.append(rikt[j][k][t+1])
                    Y.append(1)
                    Y.append(1)
                    print X
                    print Y
                    model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=X, val=Y)],senses = ["L"], rhs = [1]);
            


#model.parameters.timelimit.set(300)


# Solving the model
try:
    model.solve()
except CplexSolverError, e:
    print "Exception raised during solve: " + e
else:
    solution = model.solution   



#Conflict refiner
#model.conflict.refine(model.conflict.all_constraints())
#model.conflict.write(filename="confictin_in_DFS3")



# Code for geting values of the variables at soultuon

print "objective"
solution.get_objective_value()  

for t in range(P):
    for k in range(K*H):
        for i in V:
            for j in V:
                if i != j:
                    for c in Com:
                        varName = "z."+str(i)+"."+str(j)+"."+str(k)+"."+str(t) + "."+str(c)
                        if solution.get_values(varName) !=0:
                            print varName
                            print solution.get_values(varName)                  
    

#i=1
#for j in V:
#    if i != j:
#        for k in range(K*H):
#            for t in range(P):
#                for c in Com:
#                    varName = "z."+str(i)+"."+str(j)+"."+str(k)+"."+str(t) + "."+str(c)
##                        if solution.get_values(varName) !=0:
#                    print varName
#                    print solution.get_values(varName)   

for i in V:
    for j in V:
        if i != j:
            for l in range(L):
                for t in range(P):
                    for c in Com:
                        varName = "v."+str(i)+"."+str(j)+"."+str(l)+"."+str(t)+ "."+str(c)
                        if solution.get_values(varName) !=0:
                            print varName
                            print solution.get_values(varName)   

sumx=0
for t in range(P):
    for k in range(K*H):
        for i in V:
            for j in V:
                if i != j:
                    varName = "x."+str(i)+"."+str(j)+"."+str(k)+"."+str(t)
                    if solution.get_values(varName) !=0:
                        print varName
                        print solution.get_values(varName) 
                        sumx=sumx + solution.get_values(varName)*a[i][j]
                        


sumy=0
for i in S:
    for j in S:
        if i != j:
            for l in range(L):
                for t in range(P):
                    varName = "y."+str(i)+"."+str(j)+"."+str(l)+"."+str(t)
                    if solution.get_values(varName) !=0:
                        print varName
                        print solution.get_values(varName) 
                        sumy=sumy + solution.get_values(varName)*h[i][j]
sumd=0
for i in V:
    if i !=0 :
        for t in range(P):
            for c in Com:
                 varName = "d."+str(i)+"."+str(t)+ "."+str(c)
                 if solution.get_values(varName) != 0 : 
                    print varName 
                    print solution.get_values(varName)
                    sumd=sumd + solution.get_values(varName)
print sumd
                
for i in V:
    if i!=0:
        for k in range(K*H):
            for t in range(P):
                for c in Com:
                     varName = "d."+str(i)+"."+str(k)+"."+str(t)+ "."+str(c)
                     if solution.get_values(varName) != 0 : 
                        print varName 
                        print solution.get_values(varName)
                        

            
            
            
for i in V:
    for t in range(P):
        for c in Com:
             varName = "D."+str(i)+"."+str(t)+ "."+str(c)
             if solution.get_values(varName) != 0 : 
                print varName 
                print solution.get_values(varName)
                

for i in V:
    for t in range(P):
        for c in Com:
            varName = "S."+str(i)+"."+str(t)+ "."+str(c)
            if solution.get_values(varName) != 0 : 
                print varName 
                print solution.get_values(varName)



# Printing optimal objective
print "objective"
print solution.get_objective_value() 

# optimal primary objective 
print "sumd"
print sumd 

# optimal secondary objective 
print "sumx"
print sumx 
print "sumy"
print sumy

