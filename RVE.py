##################################################################
########## Script to create RVE with defects #####################
##################################################################
# This script creates a RVE with defects. The defects are spherical inclusions.
# The script reads the inclusions from a csv file. The csv file should have the following columns:
# x, y, z, phi, theta, radius
# x, y, z are the coordinates of the center of the inclusion
# phi, theta are the angles of the inclusion
# radius is the radius of the inclusion
# The script creates a cube with the inclusions. The cube is meshed with tetrahedral elements.
# The script is written for Abaqus 2020
# The script is written by Hakan Celik, March 2024
##################################################################
########## Script to create RVE with defects #####################
##################################################################

from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import csv

session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)



### Read inclusions from CSV file ###

poren_size = 5
ModelName = "RVE_Foam_{}_nm_test".format(poren_size)

filename = r"D:\Plasma\Data\{}nm_pore_list_000.csv".format(poren_size)
inclusion_list = []
with open(filename, 'r') as file:
    reader = csv.reader(file, delimiter=' ')
    next(reader)  # Skip the header
    for row in reader:
        # Process each row of the CSV file
        row = [float(item) for item in row]
        #print(row)
        inclusion_list.append(row)



### Set initial cube parameters ###

cube_width = 200.0
cube_height = 200.0
cube_length = 200.0

Inclusions = [
    #[i, lambda_1, lambda_2, lambda_3, x_i, y_i,  z_i, phi_i, theta_i],
    [1, cube_width*0.25, cube_width*0.25, cube_width*0.25, cube_width*0.5, cube_width*0.5, cube_width*0.5, 0.0, 0.0],
#   [2, cube_width/2.0, cube_height/2.0, cube_length/2.0, cube_width*0.5, cube_width*0.5, cube_width*0.5, 0.0, 0.0],
]



### Create a list with name Inclusions ###

Inclusions = []
i = 1
for inclusion in inclusion_list:
    #print (inclusion)
    x_i = inclusion [0]
    y_i = inclusion [1]
    z_i = inclusion [2]
    radius = inclusion[4]
    
    Inclusions.append([i, radius, x_i, y_i, z_i])

    i = i + 1



### Calculate the max Values for x,y,z for all inclusions ###

max_x_i = max(inclusion[2] for inclusion in Inclusions)
max_y_i = max(inclusion[3] for inclusion in Inclusions)
max_z_i = max(inclusion[4] for inclusion in Inclusions)
max_radius = max(inclusion[1] for inclusion in Inclusions)
print (max_x_i, max_y_i, max_z_i)
print (max(max_x_i, max_y_i, max_z_i))
max_corner = max(max_x_i, max_y_i, max_z_i)



### Increse the cube dimensions to ensure all inclusions fit within the cube ###

cube_width  = max_corner + 4*max_radius
cube_height = max_corner + 4*max_radius
cube_length = max_corner + 4*max_radius



### Assign the model object to the variable m to interact with the model programmatically through m ###

mdb.Model(name=ModelName, modelType=STANDARD_EXPLICIT)
m = mdb.models[ModelName]



############## Create main cube and assigne material properties #######################

# Define the material properties for the cube
def define_cube_material(m, material_name):
    m.Material(name=material_name)
    m.materials[material_name].Elastic(
        table=((200000.0, 0.3),)  # Example properties: Young's modulus and Poisson's ratio
    )
    # Add other material properties as needed

# Define the section
def create_section(m, section_name, material_name):
    m.HomogeneousSolidSection(name=section_name, material=material_name, thickness=None)

# Create the cube
m.ConstrainedSketch(name='__profile__', sheetSize=400.0)
m.sketches['__profile__'].rectangle(
    point1=(0.0, 0.0),
    point2=(cube_width, cube_height)
)

m.Part(
    dimensionality=THREE_D,
    name='Cube',
    type=DEFORMABLE_BODY
)
PC = m.parts['Cube']
PC.BaseSolidExtrude(
    depth=cube_length,
    sketch=m.sketches['__profile__']
)
del m.sketches['__profile__']

# Define material and section for the cube
cube_material_name = 'CubeMaterial'
define_cube_material(m, cube_material_name)

cube_section_name = 'CubeSection'
create_section(m, cube_section_name, cube_material_name)

# Create a region object for the cube part
region = PC.Set(cells=PC.cells, name='CubeRegion')  # Create a region using a set of cells

# Assign the section to the cube part
PC.SectionAssignment(
    sectionName=cube_section_name,
    region=region,  # Provide the Region object here
    offset=0.0,
    offsetType=MIDDLE_SURFACE,
    thicknessAssignment=FROM_SECTION
)

print('Cube material and section assigned.')

########### Cube correctly positioned within the simulation environment for subsequent operations (meshing, BC, Load applications) ############################################

m.rootAssembly.DatumCsysByDefault(CARTESIAN)
m.rootAssembly.Instance(
    dependent=OFF,
    name='Cube-1', 
    part=PC
    )

m.rootAssembly.translate(
    instanceList=('Cube-1',),
    vector=(-max_radius, -max_radius, -max_radius)
    )



############################ Create spherical part #######################################################

def create_spherical_part(m, partname, i,  radius):
    PName = partname + str(i)
    
    m.ConstrainedSketch(name='__profile__', sheetSize=200.0)
    m.sketches['__profile__'].ConstructionLine(point1=(0.0, 
        -100.0), point2=(0.0, 100.0))

    m.sketches['__profile__'].Line(point1=(0.0, radius), point2=
        (0.0, -radius))
    m.sketches['__profile__'].ArcByCenterEnds(
        center=(0.0,  0.0), 
        direction=CLOCKWISE, 
        point1=(0.0, radius),
        point2=(0.0, -radius)
        )
  
    PE = m.Part(
        dimensionality=THREE_D,
        name=PName, 
        type=DEFORMABLE_BODY
        )
    PE.BaseSolidRevolve(
        angle=360.0, 
        flipRevolveDirection=OFF, 
        sketch=m.sketches['__profile__'])
    del m.sketches['__profile__']


    PE.DatumPlaneByPrincipalPlane(
        offset=0.0,
        principalPlane=XYPLANE
        )
    PE.DatumPlaneByPrincipalPlane(
        offset=0.0,
        principalPlane=YZPLANE
        )
    PE.DatumPlaneByPrincipalPlane(
        offset=0.0,
        principalPlane=XZPLANE
        )


    PE.PartitionCellByDatumPlane(
       cells=PE.cells.findAt(((0.0,0.0, 0.0), )), 
       datumPlane=PE.datums[2]
       )

    PE.PartitionCellByDatumPlane(
        cells=PE.cells.findAt(
            ((0.0, 0.0, 0.01), ),
            ((0.0, 0.0, -0.01), ),
            ), 
        datumPlane=PE.datums[3]
        )

    PE.PartitionCellByDatumPlane(
        cells=PE.cells.findAt(
            ((0.01, 0.0, 0.01), ),
            ((-0.01, 0.0, -0.01), ),
            ((0.01, 0.0, -0.01), ),
            ((-0.01, 0.0, 0.01), ),
            ), 
        datumPlane=PE.datums[4]
        )
    


    print ('Spherical part: ' + str(i) )
    return PE



### Create instances of inclusions and assign material properties ###
instances_list = []
material_sections = []  # Store material and section assignments

for Inclusion in Inclusions:
    i = Inclusion[0]
    radius = Inclusion[1]
    x_i = Inclusion[2]
    y_i = Inclusion[3]
    z_i = Inclusion[4]

    # Create spherical part
    PI_i = create_spherical_part(m, 'PI_', i, radius)

    # Create unique material for each inclusion if needed
    material_name = 'Inclusion_'+str(i)
    m.Material(name=material_name)
    # m.materials[material_name].Elastic(
    #    table=((300000.0 + i * 10000.0, 0.25 + i * 0.01),)  # Example properties: Young's modulus and Poisson's ratio
    #)
    diffusion_coefficient = 296 * radius * 1e-9  # Example: different diffusion coefficient for different ratius
    m.materials[material_name].Diffusivity(
        table=((diffusion_coefficient,),)  # Diffusion coefficient
    )
    solubility_coefficient = 1  # Example: solubility coefficient constant
    m.materials[material_name].Solubility(
        table=((solubility_coefficient,),)  # Solubility coefficient
    )


    # Create section for the spherical inclusion
    section_name = 'InclusionSection_'+str(i)
    m.HomogeneousSolidSection(name=section_name, material=material_name, thickness=None)

    # Assign the section to the spherical inclusion part
    PI_i.SectionAssignment(
        region=(PI_i.cells,),
        sectionName=section_name
    )

    material_sections.append((PI_i, section_name))

    # Create instance of the spherical inclusion
    mPI = m.rootAssembly.Instance(
        dependent=OFF,
        name='Inclusion_' + str(i),
        part=PI_i
    )

    m.rootAssembly.translate(
        instanceList=('Inclusion_' + str(i),),
        vector=(x_i, y_i, z_i)
    )

    instances_list.append(mPI)

instances_list = tuple(instances_list)

# Combine cube and inclusions using Boolean operation
m.rootAssembly.InstanceFromBooleanMerge(
    name='CombinedGeometry',
    instances=(m.rootAssembly.instances['Cube-1'],) + instances_list,
    keepIntersections=ON,
    originalInstances=DELETE,
    domain=GEOMETRY
)

print("Cube and inclusions combined into a single geometry.")

# Reference to the combined part
combined_part = m.parts['CombinedGeometry']

# Reassign material properties and sections to the combined part
# For the cube
cube_material_name = 'CubeMaterial'
cube_section_name = 'CubeSection'
define_cube_material(m, cube_material_name)
create_section(m, cube_section_name, cube_material_name)
cube_region = combined_part.Set(cells=combined_part.cells.findAt(((cube_width/2, cube_height/2, cube_length/2),)), name='CubeRegion')
combined_part.SectionAssignment(
    sectionName=cube_section_name,
    region=cube_region,
    offset=0.0,
    offsetType=MIDDLE_SURFACE,
    thicknessAssignment=FROM_SECTION
)

# For the inclusions
for PI_i, section_name in material_sections:
    inclusion_region = combined_part.Set(cells=combined_part.cells.findAt(((PI_i.cells[0].pointOn[0]),)), name='InclusionRegion_' + section_name)
    combined_part.SectionAssignment(
        sectionName=section_name,
        region=inclusion_region,
        offset=0.0,
        offsetType=MIDDLE_SURFACE,
        thicknessAssignment=FROM_SECTION
    )

print("Material properties and sections reassigned to the combined geometry.")

### Meshing the combined geometry ###
# combined_part.seedPart(size=2.5, deviationFactor=0.1, minSizeFactor=0.1)
# combined_part.generateMesh()

# print("Meshing completed for the combined geometry.")




