from paramters import param
from buildModel import cplex_model
from buildModelTree import cplex_model_tree
from heurestic import heurestic 
from getSolution import printSolution


Path  = "/users/abhinavkhare/Documents/Tree/Data"
sheetname1 = "sheet1"
sheetname2 = "sheet3"
sheetname3 = "sheet4"
sheetname4 = "sheet5"
sheetname5 = "sheet6"
sheetname6 = "sheet7"
sheetname7 = "sheet8"

parameters = param()
parameters.getParamModel1(Path, sheetname1, sheetname2, sheetname3, sheetname4, sheetname5,
                         sheetname6)

solution1 = cplex_model(parameters)
printSolution(solution1, parameters)

parameters.getParamModeltree(Path, sheetname1, sheetname2, sheetname3, sheetname7, sheetname5,
                         sheetname6)

solution2 = cplex_model_tree(parameters)
printSolution(solution2, parameters)

solution3, xijktc, zijktc, ditc, diktc, rikt = heurestic(parameters)
printSolution(solution3, parameters)