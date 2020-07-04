#!/usr/bin/env python
import os
import numpy as np
import time
import placentagen as pg
from opencmiss.iron import iron

## Definition of the shape of the mesh
volume = 4.7e5  # mm^3 #volume of placenta
thickness = 22.1  # mm #thickness of placenta (z-axis dimension)
ellipticity = 1.00  # no units #ellipticity of placenta - ratio of y to x axis dimensions

# define elements in x-, y- and z-
size_el = 10  # mm #element size
circle_prop = 0.85  # maximum proportion of z surface taken up by circular mesh - helps to get good curvature in sides of mesh
squareSizeRatio = 0.50  # proportion of cross section comprised of aquare elements

el_type = 2  # Put 1 for linear element, 2 for quadratic element

# Inlet locations
stem_file = 'stem_xy.txt'
spiral_rad = 2.3#mm
# Options
# export generated mesh for visualisation
export_mesh = False
filename_mesh = 'expected-results/placenta_mesh'

export_results = True
export_directory = 'output'

# This is a constant porosity example so need to define porosity and permeability over viscosity
porosity = 0.4
perm_over_vis = 0.8  # permiability is vicosity is

initial_velocity = 0.0

dv_bc_type = 'velocity'  # or pressure
dv_bc_value = 26.5  # mm3 per sec inlet and outlet vel of spiral artery and decidual vein

sa_bc_type = 'velocity'  # or pressure
sa_bc_value = 26.5  # mm3 per sec inlet and outlet vel of spiral artery and decidual vein

wall_bc_type = 'no_slip'  # or no_penetration

debug = False
# Generate an ellipodial mesh for simulation, w
pl_mesh = pg.gen_3d_ellipsoid_structured(size_el, volume, thickness, ellipticity, squareSizeRatio, circle_prop, el_type,
                                         debug)

# Optional - export mesh prior to solution
if (export_mesh):
    pg.export_ex_coords(pl_mesh['nodes'], 'placenta', filename_mesh, 'exnode')
    if (el_type == 1):
        pg.export_exelem_3d_linear(pl_mesh['placental_el_con'], 'placenta', filename_mesh)  # Use this for linear
    elif (el_type == 2):
        pg.export_exelem_3d_quadratic(pl_mesh['elems'], 'placenta', filename_mesh)  # Use this for quadratic
    else:
        print('Cannot export elements for visualisation as element type does not exist')
print(pl_mesh['nodes'])
vnode = pg.identify_vessel_node(pl_mesh['nodes'], pl_mesh['surface_nodes'], stem_file, 1, 2, volume, thickness,ellipticity)

numberOfDimensions = 3
numberOfComponents = 1

# Set problem parameters
(coordinateSystemUserNumber,
 regionUserNumber,
 basisUserNumber,
 generatedMeshUserNumber,
 meshUserNumber,
 decompositionUserNumber,
 geometricFieldUserNumber,
 equationsSetFieldUserNumber,
 dependentFieldUserNumber,
 materialFieldUserNumber,
 equationsSetUserNumber,
 problemUserNumber) = range(1, 13)

number_of_dimensions = 3
number_of_mesh_components = 1
total_number_of_elements = len(pl_mesh["elems"])
total_number_of_nodes = len(pl_mesh["nodes"])
mesh_component_number = 1

# Create a RC coordinate system
coordinateSystem = iron.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber)
coordinateSystem.dimension = 3
coordinateSystem.CreateFinish()

# Create a region
region = iron.Region()
region.CreateStart(regionUserNumber, iron.WorldRegion)
region.label = "DarcyRegion"
region.coordinateSystem = coordinateSystem
region.CreateFinish()

# Create a tri-quadratic Lagrange basis
basis = iron.Basis()
basis.CreateStart(basisUserNumber)
basis.TypeSet(iron.BasisTypes.LAGRANGE_HERMITE_TP)
basis.numberOfXi = 3
if (el_type == 2):
    basis.interpolationXi = [iron.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE] * 3
    numberOfNodesXi = 3
    numberOfGaussXi = 3
elif (el_type == 1):
    basis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE] * 3
    numberOfNodesXi = 2
    numberOfGaussXi = 2
basis.QuadratureNumberOfGaussXiSet([numberOfGaussXi] * 3)
basis.quadratureLocalFaceGaussEvaluate = True
basis.CreateFinish()

# Start the creation of the imported mesh in the region
mesh = iron.Mesh()
mesh.CreateStart(meshUserNumber, region, number_of_dimensions)
mesh.NumberOfComponentsSet(number_of_mesh_components)
mesh.NumberOfElementsSet(total_number_of_elements)

# Define nodes for the mesh
nodes = iron.Nodes()
nodes.CreateStart(region, total_number_of_nodes)
# Refers to nodes by their user number as described in the original mesh
nodes.UserNumbersAllSet(pl_mesh['node_list'])
nodes.CreateFinish()

elements = iron.MeshElements()
elements.CreateStart(mesh, mesh_component_number, basis)

# Set the nodes pertaining to each element
for idx, elem_details in enumerate(pl_mesh['elems']):
    elem_details = elem_details.astype('int32')
    element_nodes_array = elem_details[1:28].astype('int32')
    elements.NodesSet(idx + 1, elem_details[1:28].astype('int32'))

# Refers to elements by their user number as described in the original mesh
elements.UserNumbersAllSet(pl_mesh['elems'][:, 0])
elements.CreateFinish()

mesh.CreateFinish()

# Dedine number of computational nodes and mesh properties
numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
computationalNodeNumber = iron.ComputationalNodeNumberGet()
# Create a decomposition for the mesh
decomposition = iron.Decomposition()
decomposition.CreateStart(decompositionUserNumber, mesh)
decomposition.type = iron.DecompositionTypes.CALCULATED
decomposition.numberOfDomains = numberOfComputationalNodes
decomposition.CalculateFacesSet(True)
decomposition.CreateFinish()

# Create a field for the geometry
geometricField = iron.Field()
geometricField.CreateStart(geometricFieldUserNumber, region)
geometricField.meshDecomposition = decomposition
# Set the scaling to use
geometricField.ScalingTypeSet(iron.FieldScalingTypes.NONE)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, 1)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 2, 1)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 3, 1)
geometricField.CreateFinish()

# Update the geometric field parameters
geometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

for idx, node_locs in enumerate(pl_mesh["nodes"]):
    node_num = pl_mesh['node_list'][idx]
    [x, y, z] = node_locs[1:4]

    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1,
                                            int(node_num), 1, x)
    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1,
                                            int(node_num), 2, y)
    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1,
                                            int(node_num), 3, z)

geometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

# Create standard Darcy equations set
equationsSetField = iron.Field()
equationsSet = iron.EquationsSet()
equationsSetSpecification = [iron.EquationsSetClasses.FLUID_MECHANICS,
                             iron.EquationsSetTypes.DARCY_EQUATION,
                             iron.EquationsSetSubtypes.DARCY_BRINKMAN]
equationsSet.CreateStart(equationsSetUserNumber, region, geometricField,
                         equationsSetSpecification, equationsSetFieldUserNumber, equationsSetField)
equationsSet.CreateFinish()

# Create dependent field
dependentField = iron.Field()
equationsSet.DependentCreateStart(dependentFieldUserNumber, dependentField)
dependentField.VariableLabelSet(iron.FieldVariableTypes.U, "Dependent")
dependentField.DOFOrderTypeSet(iron.FieldVariableTypes.U, iron.FieldDOFOrderTypes.SEPARATED)
dependentField.DOFOrderTypeSet(iron.FieldVariableTypes.DELUDELN, iron.FieldDOFOrderTypes.SEPARATED)
equationsSet.DependentCreateFinish()

# Create material field
materialField = iron.Field()
equationsSet.MaterialsCreateStart(materialFieldUserNumber, materialField)
materialField.VariableLabelSet(iron.FieldVariableTypes.U, "Material")
equationsSet.MaterialsCreateFinish()

materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, porosity)
materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 2,
                                          perm_over_vis)

for knode in range(0, len(pl_mesh['node_list'])):
    materialField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1,
                                           int(knode + 1), 1, porosity)  # set porosity
    materialField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1,
                                           int(knode + 1), 2, perm_over_vis)  # set perm_over_vis

# Initialise dependent field

dependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1,
                                           initial_velocity)

# Create equations
equations = iron.Equations()
equationsSet.EquationsCreateStart(equations)
equations.sparsityType = iron.EquationsSparsityTypes.SPARSE
equations.outputType = iron.EquationsOutputTypes.NONE
equationsSet.EquationsCreateFinish()

# Create Darcy equation problem
problem = iron.Problem()
problemSpecification = [iron.ProblemClasses.FLUID_MECHANICS,
                        iron.ProblemTypes.DARCY_EQUATION,
                        iron.ProblemSubtypes.DARCY_BRINKMAN]
problem.CreateStart(problemUserNumber, problemSpecification)
problem.CreateFinish()

# Create control loops
problem.ControlLoopCreateStart()
problem.ControlLoopCreateFinish()

# Create problem solver
solver = iron.Solver()
problem.SolversCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 1, solver)
# solver.outputType = iron.SolverOutputTypes.SOLVER
solver.outputType = iron.SolverOutputTypes.NONE
solver.linearType = iron.LinearSolverTypes.ITERATIVE
solver.linearIterativeAbsoluteTolerance = 1.0E-12
solver.linearIterativeRelativeTolerance = 1.0E-12
problem.SolversCreateFinish()

## Create solver equations and add equations set to solver equations
solver = iron.Solver()
solverEquations = iron.SolverEquations()
problem.SolverEquationsCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 1, solver)
solver.SolverEquationsGet(solverEquations)
solverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()

# Create boundary conditions
boundaryConditions = iron.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)

pi = np.pi
z_radius = thickness / 2.0
x_radius = np.sqrt(volume * 3.0 / (4.0 * pi * ellipticity * z_radius))
y_radius = ellipticity * x_radius

# INLET BC
if (sa_bc_type == 'velocity'):
    for i in range(0, len(vnode['spiral_array'])):
        x = pl_mesh['nodes'][vnode['spiral_array'][i]][1]
        y = pl_mesh['nodes'][vnode['spiral_array'][i]][2]
        z = pl_mesh['nodes'][vnode['spiral_array'][i]][3]
        print(x, y, z)

        vz = -1 * sa_bc_value * np.sqrt(
            x ** 2 / x_radius ** 4 + y ** 2 / y_radius ** 4 + z ** 2 / z_radius ** 4) * x_radius ** 4 * y_radius ** 2 * z_radius ** 2 * z / (
                         x_radius ** 4 * y_radius ** 2 * z ** 2 + x_radius ** 4 * z_radius ** 4 - x_radius ** 4 * z_radius ** 2 * z ** 2 - x_radius ** 2 * z_radius ** 4 * x ** 2 + y_radius ** 2 * z_radius ** 4 * x ** 2)
        if (z == 0 and x == 0):
            vx = 0
        elif (z == 0):
            vx = -1 * sa_bc_value
        else:
            vx = z_radius ** 2 * x * vz / (x_radius ** 2 * z)
        if (y == 0):
            vy = 0
        elif (z == 0):
            vy = -1 * sa_bc_value
        else:
            vy = vz * (z_radius ** 2 / z - z_radius ** 2 * x ** 2 / (x_radius ** 2 * z) - z) / y

        boundaryConditions.SetNode(dependentField, iron.FieldVariableTypes.U, 1, 1, int(vnode['spiral_array'][i]), 1,
                                   iron.BoundaryConditionsTypes.FIXED, vx)
        boundaryConditions.SetNode(dependentField, iron.FieldVariableTypes.U, 1, 1, int(vnode['spiral_array'][i]), 2,
                                   iron.BoundaryConditionsTypes.FIXED, vy)
        boundaryConditions.SetNode(dependentField, iron.FieldVariableTypes.U, 1, 1, int(vnode['spiral_array'][i]), 3,
                                   iron.BoundaryConditionsTypes.FIXED, vz)
else:
    for i in range(0, len(vnode['spiral_array'])):
        boundaryConditions.SetNode(dependentField, iron.FieldVariableTypes.U, 1, 1, int(vnode['spiral_array'][i]), 4,
                                   iron.BoundaryConditionsTypes.FIXED, sa_bc_value)

# OUTLET BC
if (dv_bc_type == 'velocity'):
    for i in range(0, len(vnode['decidual_array'])):
        x = pl_mesh['nodes'][vnode['decidual_array'][i]][1]
        y = pl_mesh['nodes'][vnode['decidual_array'][i]][2]
        z = pl_mesh['nodes'][vnode['decidual_array'][i]][3]
        print(x, y, z)

        vz = dv_bc_value * np.sqrt(
            x ** 2 / x_radius ** 4 + y ** 2 / y_radius ** 4 + z ** 2 / z_radius ** 4) * x_radius ** 4 * y_radius ** 2 * z_radius ** 2 * z / (
                         x_radius ** 4 * y_radius ** 2 * z ** 2 + x_radius ** 4 * z_radius ** 4 - x_radius ** 4 * z_radius ** 2 * z ** 2 - x_radius ** 2 * z_radius ** 4 * x ** 2 + y_radius ** 2 * z_radius ** 4 * x ** 2)
        if (z == 0 and x == 0):
            vx = 0
        elif (z == 0):
            vx = dv_bc_value
        else:
            vx = z_radius ** 2 * x * vz / (x_radius ** 2 * z)
        if (y == 0):
            vy = 0
        elif (z == 0):
            vy = dv_bc_value
        else:
            vy = vz * (z_radius ** 2 / z - z_radius ** 2 * x ** 2 / (x_radius ** 2 * z) - z) / y
        boundaryConditions.SetNode(dependentField, iron.FieldVariableTypes.U, 1, 1, int(vnode['decidual_array'][i]), 1,
                                   iron.BoundaryConditionsTypes.FIXED, vx)
        boundaryConditions.SetNode(dependentField, iron.FieldVariableTypes.U, 1, 1, int(vnode['decidual_array'][i]), 2,
                                   iron.BoundaryConditionsTypes.FIXED, vy)
        boundaryConditions.SetNode(dependentField, iron.FieldVariableTypes.U, 1, 1, int(vnode['decidual_array'][i]), 3,
                                   iron.BoundaryConditionsTypes.FIXED, vz)
else:
    for i in range(0, len(vnode['decidual_array'])):
        boundaryConditions.SetNode(dependentField, iron.FieldVariableTypes.U, 1, 1, int(vnode['decidual_array'][i]), 4,
                                   iron.BoundaryConditionsTypes.FIXED, dv_bc_value)

for i in range(0, len(vnode['surfnode_ex_vessel'])):
    boundaryConditions.SetNode(dependentField, iron.FieldVariableTypes.U, 1, 1, int(vnode['surfnode_ex_vessel'][i]), 3,
                               iron.BoundaryConditionsTypes.FIXED_WALL, 0.0)
    boundaryConditions.SetNode(dependentField, iron.FieldVariableTypes.U, 1, 1, int(vnode['surfnode_ex_vessel'][i]), 2,
                               iron.BoundaryConditionsTypes.FIXED_WALL, 0.0)
    boundaryConditions.SetNode(dependentField, iron.FieldVariableTypes.U, 1, 1, int(vnode['surfnode_ex_vessel'][i]), 1,
                               iron.BoundaryConditionsTypes.FIXED_WALL, 0.0)

solverEquations.BoundaryConditionsCreateFinish()
start_time = time.time()
# Solve the problem
problem.Solve()
end_time = time.time()
print('Total time for solve = ' + str((end_time - start_time) / 60.0) + ' mins')

if (export_results):
    if not os.path.exists(export_directory):
        os.makedirs(export_directory)
        ## Export results
    export_file = export_directory + '/StaticDarcy'
    fields = iron.Fields()
    fields.CreateRegion(region)
    fields.NodesExport(export_file, "FORTRAN")
    fields.ElementsExport(export_file, "FORTRAN")
    fields.Finalise()

iron.Finalise()
raise SystemExit