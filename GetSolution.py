

def printSolution(solution, param):

        
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
                                sumy= sumy + solution.get_values(varName)*h[i][j]
                                
        
        
        
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
        
        print 'Sum of Delivered Quantities"
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
        
        
        
        
        # printng optimal objective 
        print "objective"
        print solution.get_objective_value() 
        
        
        # Optimal primary objective
        print "sumd"
        print sumd 
        
        
        # Optimal secondary objectives
        print "sumx"
        print sumx 
        print "sumy"
        print sumy

        
