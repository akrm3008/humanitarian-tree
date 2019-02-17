# Description 
The last mile delivery in humanitarian relief supply often happens on a tree or an almost-tree network (which can be approximated to a tree network). I developed a model that takes advantage of the tree structure while making the vehicle routing constraints and built an approximate solution method.

# model1_cplex.py 
First, I built a multi-period multi-modal relief delivery model incorporating a tree network for last mile delivery. We developed a mixed integer programming (MIP) formulation with the goal of minimizing the unsatisfied demand of the population and used to solve the relief routing problem of Nepal 2015 earthquake. Details of the model will uploaded with research paper.

# model2_cplex.py 
In the second model, I take the advantage of the fact that the secondary model is almost tree and can be approximated a tree. Hence, iI first coverted the secondary network to a minimum spanning tree.  Next, I built a reformulation of the original model in which the routing constraints have been changed and reduced as they take advantage  of the tree structure of the secondary network. Details of these changes The model is coded in Cplex and used to solve the case study. Reduction in computational time was expected. We found that this gave an order of magnitude reduction in computational time. Details will be uploaded with the research paper.

# Approximate_method.py
To further improve computational efficiency, I developed a heuristic solution method based on a decomposition scheme applied to the tree network formulation. This led to the Capacitated Vehicle Routing Problem on trees with split deliveries (TCVRP-SD), for which I derived a closed-form solution. This decomposition scheme resulted in a further order of magnitude reduction in computation time. (Details will be shared in the research paper) 
