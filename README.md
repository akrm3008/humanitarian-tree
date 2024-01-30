# Description 

The last mile delivery in humanitarian relief supply chains is an important topic of research in humanitarian logistics. In the literature of Operations Research, there are multiple mathematical optimizations models which optimize supply chain decisions to minimize the delivery times, costs and unsatified demand of the people in need. We found that the last mile delivery network in humanitirian supply chains is often a almost-tree graph (a graph with very small number of cycles) due to damages caused to infrastructure because of the disaster. In this research, we explored the idea of exploting the approximate tree structure of these graphs to improve the compuational efficiency of solution methods for the last mile delivery problems.

The problem was modeled using two mixed integer programming (MIP) formulation: (a) Model 1 is a general MIP formulation of the problem (b) Model 2 or the Model with tree route formulation is a MIP which utilises the structure of tree graphs to formulate contraints for vehicle routing which reduces the computational complexity of the model. We solved the two models using CPLEX and found CPLEX solved the Model 2 faster as expected. We also proved some mathematical properties of vehicle routing on trees with split deliveries and used them  a decomposotion of the model and herestics to solve the decomposed model. The heurestic reduces comptational times manifolds and and also closely approximates the optimal solution. We tested the models and solution methods using Nepal Earthquake (2015) data. 
# Files 

1. Parameters.py contains class param to initialise parameters or import data and process it to obtain parameters for both the models.
2. BuildModel.py contains the functions to build a MIP for model 1 and solve it using CPLEX.
3. BuildModelTree.py contains the functions to build a MIP for model 2 and solve it using CPLEX.
4. Heurestic.py contains function to run the heurestic we created to solve the decomposition of the model.
5. GetSolution.py contains functions to print solution.
6. main.py is the script for performing experiments.




