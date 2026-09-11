# author: Karan Taneja
# unit system: mm-N-MPa
# dimension: plane strain
# documentation: 'Astrocytes in white matter ...'

from abaqus import *
from abaqusConstants import *
import part
import material
import section
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import regionToolset
import numpy as np
# Units: mm, kPa, N,

#######################################################################################
#######################################################################################
def Create_Quarter_Ellipse(ModelName, PartName, Dimensions):
    MajorAxis_G, MajorAxis_W, MinorAxis_G, MinorAxis_W, Depth = Dimensions

    cliCommand("""session.journalOptions.setValues(replayGeometry=COORDINATE,recoverGeometry=COORDINATE)""")

    # Define Model
    mdb.Model(name=ModelName, modelType=STANDARD_EXPLICIT)

    s1 = mdb.models[ModelName].ConstrainedSketch(name='__profile__', sheetSize=100.0)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=STANDALONE)
    s1.EllipseByCenterPerimeter(center=(0.0, 0.0), axisPoint1=(MajorAxis_G, 0.0),
        axisPoint2=(0.0, MinorAxis_G))
    p = mdb.models[ModelName].Part(name='Part-1', dimensionality=THREE_D,
        type=DEFORMABLE_BODY)
    p.BaseSolidExtrude(sketch=s1, depth=Depth)
    s1.unsetPrimaryObject()

    f, e, d1 = p.faces, p.edges, p.datums
    t = p.MakeSketchTransform(sketchPlane=f[1], sketchUpEdge=e[0],
        sketchPlaneSide=SIDE1, origin=(0.0, 0.0, 2.0))
    s = mdb.models[ModelName].ConstrainedSketch(name='__profile__',
        sheetSize=187.48, gridSpacing=4.68, transform=t)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=SUPERIMPOSE)
    p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
    s.Line(point1=(MajorAxis_G*10, 0.0), point2=(0.0, 0.0))
    s.HorizontalConstraint(entity=g[3], addUndoState=False)
    s.PerpendicularConstraint(entity1=g[2], entity2=g[3], addUndoState=False)
    s.Line(point1=(0.0, 0.0), point2=(0.0, MinorAxis_G*10))
    s.VerticalConstraint(entity=g[4], addUndoState=False)
    s.PerpendicularConstraint(entity1=g[3], entity2=g[4], addUndoState=False)
    s.CoincidentConstraint(entity1=v[2], entity2=g[2], addUndoState=False)
    s.EllipseByCenterPerimeter(center=(0.0, 0.0), axisPoint1=(MajorAxis_W, 0.0),
        axisPoint2=(0.0, MinorAxis_W))
    pickedFaces = f.getSequenceFromMask(mask=('[#2 ]', ), )
    e1, d2 = p.edges, p.datums
    p.PartitionFaceBySketch(sketchUpEdge=e1[0], faces=pickedFaces, sketch=s)
    s.unsetPrimaryObject()

    c = p.cells
    pickedCells = c.getSequenceFromMask(mask=('[#1 ]', ), )
    v1, e, d1 = p.vertices, p.edges, p.datums
    p.PartitionCellByPlaneThreePoints(point1=v1[3], point2=v1[4], point3=v1[5],
        cells=pickedCells)

    c = p.cells
    pickedCells = c.getSequenceFromMask(mask=('[#2 ]', ), )
    v2, e1, d2 = p.vertices, p.edges, p.datums
    p.PartitionCellByPlaneThreePoints(point1=v2[4], point2=v2[7],
        cells=pickedCells, point3=p.InterestingPoint(edge=e1[7], rule=MIDDLE))

    c = p.cells
    pickedCells = c.getSequenceFromMask(mask=('[#4 ]', ), )
    e, d1 = p.edges, p.datums
    pickedEdges =(e[19], )
    p.PartitionCellByExtrudeEdge(line=e[3], cells=pickedCells, edges=pickedEdges,
        sense=FORWARD)

    f1 = p.faces
    p.RemoveFaces(faceList = f1[5:6]+f1[11:13]+f1[14:15]+f1[16:18],
        deleteCells=False)
    f = p.faces
    p.RemoveFaces(faceList = f[5:6]+f[7:8]+f[12:13], deleteCells=False)


#######################################################################################
#######################################################################################
def Create_Rigid_wall(ModelName):
    # Rigid plane
    s1 = mdb.models[ModelName].ConstrainedSketch(name='__profile__', sheetSize=100.0)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=STANDALONE)

    s1.rectangle(point1=(0.0, -2.00), point2=(-0.20, 4.0))
    p = mdb.models[ModelName].Part(name='wall', dimensionality=THREE_D, type=DISCRETE_RIGID_SURFACE)
    p.BaseShell(sketch=s1)
    # Reference point
    v1, e, d1, n = p.vertices, p.edges, p.datums, p.nodes
    p.ReferencePoint(point=(-0.10, 1.0, 0.0))

#######################################################################################
#######################################################################################
def Create_Step(ModelName, Step):
    TotalTime, Mass_Scaling = Step

    mdb.models[ModelName].ExplicitDynamicsStep(name='growth', previous='Initial', timePeriod=TotalTime, massScaling=((SEMI_AUTOMATIC, MODEL, AT_BEGINNING, Mass_Scaling, 0.0, None, 0, 0, 0.0, 0.0, 0, None),),)
    mdb.models[ModelName].steps['growth'].setValues(timePeriod=TotalTime, scaleFactor=1.0, linearBulkViscosity=0.0, quadBulkViscosity=0.0, improvedDtMethod=ON)
    mdb.models[ModelName].TabularAmplitude(name='Amp-1', timeSpan=STEP, smooth=SOLVER_DEFAULT, data=((0.0, 0.0), (1.0, 1.0)))

#######################################################################################
#######################################################################################
# def Create_Material(ModelName, Materials): 
#     mdb.models[ModelName].Material(name='WHIT')
#     mdb.models[ModelName].materials['WHIT'].Density(table=((density, ), ))
#     mdb.models[ModelName].materials['WHIT'].Depvar(n=10)
#     mdb.models[ModelName].materials['WHIT'].UserMaterial(mechanicalConstants=(shear_w,lambda_w,Gctx,GwByGgr,gamma_hat_pull,progen_grow_time,pull_grow_time,MajorAxis_W,MinorAxis_W,reduction,N_gyri,alpha,scaled_thresh))

#     mdb.models[ModelName].Material(name='GRAY')
#     mdb.models[ModelName].materials['GRAY'].Density(table=((density, ), ))
#     mdb.models[ModelName].materials['GRAY'].Depvar(n=7)
#     mdb.models[ModelName].materials['GRAY'].UserMaterial(mechanicalConstants=(shear_g,lambda_g,Gctx,progen_grow_time,MajorAxis_G,MinorAxis_G))

def Create_Material(ModelName, density, White_Matter_Properties, Gray_Matter_Properties):
    mdb.models[ModelName].Material(name='WHIT')
    mdb.models[ModelName].materials['WHIT'].Density(table=((density, ), ))
    mdb.models[ModelName].materials['WHIT'].Depvar(n=10)
    mdb.models[ModelName].materials['WHIT'].UserMaterial(mechanicalConstants=White_Matter_Properties)

    mdb.models[ModelName].Material(name='GRAY')
    mdb.models[ModelName].materials['GRAY'].Density(table=((density, ), ))
    mdb.models[ModelName].materials['GRAY'].Depvar(n=7)
    mdb.models[ModelName].materials['GRAY'].UserMaterial(mechanicalConstants=Gray_Matter_Properties)

#######################################################################################
#######################################################################################
def Create_Section(ModelName, PartName, Dimensions):
    MajorAxis_G, MajorAxis_W, MinorAxis_G, MinorAxis_W, Depth = Dimensions

    p = mdb.models[ModelName].parts['Part-1']
    c = p.cells


    c_subcort = c.findAt((( MajorAxis_W/2.,0.0, Depth/2.), ))
    c_cortex  = c.findAt((( MajorAxis_G - 0.1,0.0,  Depth/2.), ))

    mdb.models[ModelName].HomogeneousSolidSection(name='Subcortex', material='WHIT', thickness=None)
    region = p.Set(cells=c_subcort, name='Subcortex')
    p.SectionAssignment(region=region, sectionName='Subcortex', offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

    mdb.models[ModelName].HomogeneousSolidSection(name='Cortex', material='GRAY', thickness=None)
    region = p.Set(cells=c_cortex, name='Cortex')
    p.SectionAssignment(region=region, sectionName='Cortex', offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

#######################################################################################
#######################################################################################           
def Create_Assembly(ModelName, PartName, InstanceName, Dimensions):
    MajorAxis_G, MajorAxis_W, MinorAxis_G, MinorAxis_W, Depth = Dimensions

    p = mdb.models[ModelName].parts['Part-1']
    p1 = mdb.models[ModelName].parts['wall']
    a = mdb.models[ModelName].rootAssembly
    a.DatumCsysByDefault(CARTESIAN)
    InstanceName = 'Part-1'

    a.Instance(name=InstanceName, part=p, dependent=ON)
    a.Instance(name='wall-1', part=p1, dependent=ON)
    a.Instance(name='wall-2', part=p1, dependent=ON)

    a.rotate(instanceList=('wall-1', ), axisPoint=(0.0, 0.0, 0.01), axisDirection=(0.0, -2.0, 0.0), angle=-90.0)
    a.rotate(instanceList=('wall-2', ), axisPoint=(0.0, 0.0, 0.01), axisDirection=(0.0, -2.0, 0.0), angle=-90.0)
    a.rotate(instanceList=('wall-2', ), axisPoint=(0.0, 0.0, 0.01), axisDirection=(0.0, 0.0, -2.0), angle=-90.0)

    a.translate(instanceList=('wall-1', ), vector=(0.0,MinorAxis_W/2, 0.0))
    a.translate(instanceList=('wall-1', ), vector=(0.0,0.0, Depth/2.-0.01-0.1))
    a.translate(instanceList=('wall-2', ), vector=(MajorAxis_W, 0.0, 0.0))
    a.translate(instanceList=('wall-2', ), vector=(0.0,0.0, Depth/2.-0.01-0.1))

def y_coord(x_coord,majoraxis,minoraxis):
    output = np.sqrt((1 - (x_coord*x_coord)/(majoraxis*majoraxis)))*minoraxis
    return output

def x_coord(y_coord,majoraxis,minoraxis):
    output = np.sqrt((1 - (y_coord*y_coord)/(minoraxis*minoraxis)))*majoraxis
    return output

####################################################################################### 
####################################################################################### 
def Create_Sets(ModelName, PartName, Dimensions):
    MajorAxis_G, MajorAxis_W, MinorAxis_G, MinorAxis_W, Depth = Dimensions

    a = mdb.models[ModelName].rootAssembly
    p = mdb.models[ModelName].parts['Part-1']
    pwall = mdb.models[ModelName].parts['wall']

    f = p.faces
    fwall = pwall.faces

    e = p.edges
    c = p.cells
    v = p.vertices

    # Define coordinate of set-node
    verts = v.findAt(((MajorAxis_G, 0.0, Depth), ))

    # Define coordinate of the sets-faces
    f_back = f.findAt(((0.1, MinorAxis_G/2., 0.),),((MajorAxis_G-0.1, y_coord(MajorAxis_G-0.1,MajorAxis_G,MinorAxis_G)-0.1, 0.),))
    f_bottom = f.findAt(((MajorAxis_W/2., 0.0, Depth/2.),),((MajorAxis_G - 0.1, 0.0, Depth/2.),))#f.findAt(((MajorAxis_G/2., 0.0, Depth/2.),))
    f_side1 = f.findAt(((0, MinorAxis_G/2., Depth/2.),), ((0, MinorAxis_G - 0.1, Depth/2.),))
    f_side2 = f.findAt(((MajorAxis_W/2., 0.0, Depth/2.),),((MajorAxis_G - 0.1, 0.0, Depth/2.),))
    f_cortex_top = f.findAt(((MajorAxis_G,0.0, Depth/2),))
    f_interface = f.findAt(((MajorAxis_W,0.0, Depth/2),))
    scale_fac_G = 0.67456461111 # Apparently, this is the scaled coordinate w.r.t major axis of grey of a node formed on the ellipse contour
    scale_fac_W = 0.662037351111

    # Define coordinate of sets-edges
    e_top = e.findAt(((MajorAxis_G*scale_fac_G,y_coord(MajorAxis_G*scale_fac_G,MajorAxis_G,MinorAxis_G), Depth),))#findAt(((Length/2., 0.0, Depth), ))
    e_interface = e.findAt(((MajorAxis_W/2,y_coord(MajorAxis_W/2,MajorAxis_W,MinorAxis_W), Depth),), ((MajorAxis_W/2,y_coord(MajorAxis_W/2,MajorAxis_W,MinorAxis_W), Depth/2),), )
    e_interface_2 = e.findAt(((MajorAxis_W*scale_fac_W,y_coord(MajorAxis_W*scale_fac_W,MajorAxis_W,MinorAxis_W), Depth), ))
    s_top = f.findAt(((MajorAxis_G*scale_fac_G,y_coord(MajorAxis_G*scale_fac_G,MajorAxis_G,MinorAxis_G), Depth/2),))
    s_cortex = f.findAt(((MajorAxis_G - 0.1,y_coord(MajorAxis_G - 0.1,MajorAxis_G - 0.1,MinorAxis_G-0.1), Depth),))
    s_subcortex = f.findAt(((MajorAxis_W/2,0.1, Depth),))
    s_side1 = f.findAt(((0, MinorAxis_G/2., Depth/2.),), ((0, MinorAxis_G - 0.1, Depth/2.),))
    s_side2 = f.findAt(((MajorAxis_W/2.,0.0, Depth/2.),),((MajorAxis_G - 0.1, 0.0, Depth/2.),))
    s_interface = f.findAt(((MajorAxis_W/2,y_coord(MajorAxis_W/2,MajorAxis_W,MinorAxis_W), Depth/2),))
    # s_wall = fwall.findAt(((-1, 10, 0.0), ))#fwall.findAt(((-100E-3, 54.1, 1.1), ))
    s_wall = fwall.findAt((((0.0 + -0.20)/2., (-2.00 + 4.0)/2., 0.0), ))


    # Define coordinate of sets-cells
    cells = c.findAt(((MajorAxis_W/2,MinorAxis_W/2, Depth/2),),((0.1,MinorAxis_G-0.1, 0.0),))

    # Assign sets-faces
    region = p.Set(faces=f_back, name='back')
    region = p.Set(faces=f_bottom, name='bottom')
    region = p.Set(faces=f_side1, name='side-1')
    region = p.Set(faces=f_side2, name='side-2')
    region = p.Set(faces=s_cortex, name='cortex_front')
    region = p.Set(faces=s_subcortex, name='subcortex_front')
    region = p.Set(faces=f_cortex_top, name='cortex_top')

    region = p.Set(faces=f_interface, name='interface')

    # Assign sets-edges
    region = p.Set(edges=e_top, name='top-edge')
    region = p.Set(edges=e_interface, name='interface')
    region = p.Set(edges=e_interface_2, name='interface-edge')

    # Assign sets-surfaces
    region = p.Surface(side1Faces=s_top, name='top_s')
    region = p.Surface(side1Faces=s_cortex, name='cortex_s')
    region = p.Surface(side1Faces=s_subcortex, name='subcortex_s')
    region = p.Surface(side1Faces=s_side1, name='side1_s')
    region = p.Surface(side1Faces=s_side2, name='side2_s')
    region = p.Surface(side1Faces=s_interface, name='interface_s')
    region = pwall.Surface(side1Faces=s_wall, name='s_wall')

    # Assign sets-cells
    region = p.Set(cells=cells, name='whole-domain')

    #Assign sets-node
    p.Set(vertices=verts, name='node')

    # Rigid body surface
    p1 = mdb.models[ModelName].parts['wall']

    r = p1.referencePoints
    refPoints=(r[2], )
    p1.Set(referencePoints=refPoints, name='ref')

#######################################################################################
#######################################################################################
def Create_Boundary_Conditions(ModelName, InstanceName):
    a = mdb.models[ModelName].rootAssembly

    # Add displacement and fixed bcs

    region = a.instances[InstanceName].sets['back']
    mdb.models[ModelName].DisplacementBC(name='back', createStepName='growth', region=region, u1=UNSET, u2=UNSET, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude='Amp-1', fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)

    region = a.instances[InstanceName].sets['bottom']
    mdb.models[ModelName].DisplacementBC(name='bottom', createStepName='growth', region=region, u1=UNSET, u2=0.0, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude='Amp-1', fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)

    region = a.instances[InstanceName].sets['side-1']
    mdb.models[ModelName].DisplacementBC(name='side-1', createStepName='growth', region=region, u1=0.0, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude='Amp-1', fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)


    region = a.instances[InstanceName].sets['cortex_front']
    mdb.models[ModelName].DisplacementBC(name='cortex_front', createStepName='growth', region=region, u1=UNSET, u2=UNSET, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude='Amp-1', fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)

    region = a.instances[InstanceName].sets['subcortex_front']
    mdb.models[ModelName].DisplacementBC(name='subcortex_front', createStepName='growth', region=region, u1=UNSET, u2=UNSET, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude='Amp-1', fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)

    region = a.instances['wall-1'].sets['ref']
    mdb.models[ModelName].EncastreBC(name='ref-left', createStepName='growth', region=region, localCsys=None)

    region = a.instances['wall-2'].sets['ref']
    mdb.models[ModelName].EncastreBC(name='ref-right', createStepName='growth', region=region, localCsys=None)
    region = a.instances['Part-1'].sets['whole-domain']
    mdb.models[ModelName].DisplacementBC(name='fix_z', 
        createStepName='Initial', region=region, u1=UNSET, u2=UNSET, u3=SET, 
        ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude='Amp-1', distributionType=UNIFORM, 
        fieldName='', localCsys=None)
    mdb.models[ModelName].boundaryConditions['back'].move('growth',
        'Initial')
    mdb.models[ModelName].boundaryConditions['bottom'].move('growth',
        'Initial')
    mdb.models[ModelName].boundaryConditions['cortex_front'].move('growth',
        'Initial')
    mdb.models[ModelName].boundaryConditions['ref-left'].move('growth',
        'Initial')
    mdb.models[ModelName].boundaryConditions['ref-right'].move('growth',
        'Initial')
    mdb.models[ModelName].boundaryConditions['side-1'].move('growth',
        'Initial')
    mdb.models[ModelName].boundaryConditions['subcortex_front'].move('growth'
        , 'Initial')
    # mdb.models[ModelName].boundaryConditions['fix_z'].move('growth'
    #     , 'Initial')

#######################################################################################
###################################################################################### 
def Create_Viscous_Pressure(ModelName, InstanceName, Viscous_Pressure):
    a = mdb.models[ModelName].rootAssembly
    region = a.instances[InstanceName].surfaces['top_s']
    mdb.models[ModelName].Pressure(name='Viscous Pressure', createStepName='growth', region=region, distributionType=VISCOUS, field='', magnitude=Viscous_Pressure, amplitude='Amp-1')

#######################################################################################
#######################################################################################    
def Create_Contact(ModelName, InstanceName):
    a = mdb.models[ModelName].rootAssembly

    mdb.models[ModelName].ContactProperty('IntProp-1')
    mdb.models[ModelName].interactionProperties['IntProp-1'].TangentialBehavior(formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF, pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, table=((10.0, ), ), shearStressLimit=None, maximumElasticSlip=FRACTION, fraction=0.005, elasticSlipStiffness=None)
    mdb.models[ModelName].interactionProperties['IntProp-1'].NormalBehavior(pressureOverclosure=HARD, allowSeparation=ON, constraintEnforcementMethod=DEFAULT)
    mdb.models[ModelName].interactionProperties['IntProp-1'].GeometricProperties(contactArea=1.0, padThickness=0.01)

    mdb.models[ModelName].ContactProperty('IntProp-2')
    mdb.models[ModelName].interactionProperties['IntProp-2'].TangentialBehavior(formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF, pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, table=((10.0, ), ), shearStressLimit=None, maximumElasticSlip=FRACTION, fraction=0.005, elasticSlipStiffness=None)
    mdb.models[ModelName].interactionProperties['IntProp-2'].NormalBehavior(pressureOverclosure=HARD, allowSeparation=ON, constraintEnforcementMethod=DEFAULT)

    region = a.instances[InstanceName].surfaces['top_s']
    mdb.models[ModelName].SelfContactExp(name='Int-1', createStepName='growth', surface=region, mechanicalConstraint=KINEMATIC, interactionProperty='IntProp-1')

    region1= a.instances['wall-1'].surfaces['s_wall']#a.Surface(side1Faces=side1Faces1, name='left-wall') 
    region2=a.instances[InstanceName].surfaces['top_s']
    mdb.models[ModelName].SurfaceToSurfaceContactExp(name ='Int-2', createStepName='Initial', main = region1, secondary = region2, mechanicalConstraint=KINEMATIC, sliding=FINITE, interactionProperty='IntProp-2', initialClearance=OMIT, datumAxis=None, clearanceRegion=None)

    region1=a.instances['wall-2'].surfaces['s_wall']#a.Surface(side2Faces=side1Faces1, name='bottom_wall')
    region2=a.instances[InstanceName].surfaces['top_s']
    mdb.models[ModelName].SurfaceToSurfaceContactExp(name ='Int-3', createStepName='Initial', main = region1, secondary = region2, mechanicalConstraint=KINEMATIC, sliding=FINITE, interactionProperty='IntProp-2', initialClearance=OMIT, datumAxis=None, clearanceRegion=None)


#######################################################################################    
#######################################################################################    
def Create_Mesh(ModelName, PartName, InstanceName, Dimensions, Mesh, minsize, maxsize):
    MajorAxis_G, MajorAxis_W, MinorAxis_G, MinorAxis_W, Depth = Dimensions
    bias, cortex_size = Mesh

    # Rigif plane mesh
    p = mdb.models[ModelName].parts['wall']
    p.seedPart(size=0.15, deviationFactor=0.1, minSizeFactor=0.1)
    p = mdb.models[ModelName].parts['wall']
    p.generateMesh()

    # Bilayer Mesh for cortex

    a = mdb.models[ModelName].rootAssembly

    p = mdb.models[ModelName].parts['Part-1']
    e = p.edges
    pickedEdges = e.findAt(((1.2097431, 2.6376844, 0.0), ), ((3.1092919, 1.1744441,
    0.1), ), ((3.2761622, 1.2435281, 0.0), ), ((1.2786489, 2.8043935, 0.1), ))
    p.seedEdgeBySize(edges=pickedEdges, size=0.025, deviationFactor=0.1,
    constraint=FINER) # size =0.025

    pickedEdges = e.findAt(((3.42, 0.0, 0.025), ), ((0.0, 2.82, 0.075), ), ((3.465,
        0.0, 0.0), ), ((3.60, 0.0, 0.025), ), ((3.555, 0.0, 0.10), ), ((0.0, 2.865,
        0.10), ), ((0.0, 3.00, 0.075), ), ((0.0, 2.955, 0.0), ))
    p.seedEdgeBySize(edges=pickedEdges, size=cortex_size, deviationFactor=0.1,
        constraint=FINER) #0.5 gave bad aspect ratios

    # For subcortex, used bias to keep it finer near the interface
    pickedEdges1 = e.findAt(((0.0, 2.325, 0.0), ), ((0.0, 0.0, 0.025), ), ((2.625,
        0.0, 0.1), ))
    pickedEdges2 = e.findAt(((0.0, 0.775, 0.1), ), ((0.875, 0.0, 0.0), ))
    p.seedEdgeByBias(biasMethod=SINGLE, end1Edges=pickedEdges1,
        end2Edges=pickedEdges2, minSize=minsize, maxSize=maxsize, constraint=FINER)
    # Unstructured mesh settings
    c = p.cells
    pickedRegions = c.findAt(((2.2580953, 0.1833849, 0.1), ))
    p.deleteMesh(regions=pickedRegions)
    c = p.cells
    pickedRegions = c.findAt(((2.2580953, 0.1833849, 0.1), ))
    p.setMeshControls(regions=pickedRegions, technique=SWEEP,
        algorithm=ADVANCING_FRONT)
    partInstances =(a.instances[InstanceName], )
    # Define element type & Create mesh
    elemType1 = mesh.ElemType(elemCode=C3D8R, elemLibrary=EXPLICIT, kinematicSplit=AVERAGE_STRAIN, secondOrderAccuracy=OFF, hourglassControl=ENHANCED, distortionControl=ON)
    elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=EXPLICIT)
    elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=EXPLICIT)
    c = p.cells
    c_all = c.findAt(((MajorAxis_W/2,MinorAxis_W/2, Depth/2),),((0.1,MinorAxis_G-0.1, 0.0),))
    pickedRegions =(c_all, )
    p.setElementType(regions=pickedRegions, elemTypes=(elemType1,elemType2,elemType3))
    p.generateMesh(regions=partInstances)

   

#######################################################################################
#######################################################################################
def Create_Output(ModelName, PartName, InstanceName):
    p = mdb.models[ModelName].parts[PartName]

    # Assign desired output to sets
    regionDef=mdb.models[ModelName].rootAssembly.allInstances[InstanceName].sets['whole-domain']
    mdb.models[ModelName].FieldOutputRequest(name='F-Output-2', createStepName='growth', variables=('S', 'PE', 'LE', 'U', 'ELSE', 'COORD','SDV','EVOL'),  numIntervals=100, region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)

    regionDef=mdb.models[ModelName].rootAssembly.allInstances[InstanceName].sets['Cortex']
    mdb.models[ModelName].FieldOutputRequest(name='F-Output-3', createStepName='growth', variables=('ELSE', 'COORD', 'U','EVOL'), numIntervals=100, region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)

    regionDef=mdb.models[ModelName].rootAssembly.allInstances[InstanceName].sets['Subcortex']
    mdb.models[ModelName].FieldOutputRequest(name='F-Output-4', createStepName='growth', variables=('ELSE', 'COORD', 'U','EVOL'), numIntervals=100, region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)

    regionDef=mdb.models[ModelName].rootAssembly.allInstances[InstanceName].sets['whole-domain']
    mdb.models[ModelName].HistoryOutputRequest(name='H-Output-1', createStepName='growth', variables=('ALLIE', 'ALLKE'), numIntervals=100, region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)

#######################################################################################
#######################################################################################
def Create_Job(ModelName, JobName):
    mdb.Job(name=JobName, model=ModelName, description='', type=ANALYSIS,
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,
        memoryUnits=PERCENTAGE, explicitPrecision=DOUBLE,
        nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF,
        contactPrint=OFF, historyPrint=OFF,
        userSubroutine='vumat_pull.f',
        scratch='', resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN,
        numDomains=12, activateLoadBalancing=False, numThreadsPerMpiProcess=1,
        multiprocessingMode=DEFAULT, numCpus=12)


    mdb.jobs[JobName].writeInput(consistencyChecking=OFF)
    try: 
        mdb.jobs[JobName].submit(consistencyChecking=OFF)
        # mdb.jobs[JobName].waitForCompletion()
    except: 
        print("Job %s crashed" % JobName)

#######################################################################################
#######################################################################################


if __name__ == '__main__':
    gamma_hat_pull_list = [0.05,0.1,0.2,0.3,0.5,1.0,2.0,3.0]
    Mass_Scaling_list = [120,120,120,120,120,120,150,180]
    Viscous_Pressure_list = [1e-4,1e-5,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4]
    Job_Name_list = ['5em2','10em2','20em2','30em2','50em2','1','2','3'] #len(gamma_hat_pull_list)
    for j in range(1):

        # ======================================================
        # Dimensions
        # ======================================================
        CT = 0.15 # Units: 1/10 mm
        MajorAxis_G = 3.6 # Units: 1/10 mm
        MinorAxis_G = 3.0 # Units: 1/10 mm
        MajorAxis_W = MajorAxis_G - CT
        MinorAxis_W = MinorAxis_G - CT

        Depth = 0.1 # Units: 1/10 mm
        Dimensions = [MajorAxis_G, MajorAxis_W, MinorAxis_G, MinorAxis_W, Depth]
        # ======================================================
        # Mesh Parameters
        # ======================================================
        bias = 3                                                        # bias through width of subcortex
        ecortex = 7                                                     # number of elements in cortex                                                
        cortex_size = CT/ecortex # size of element in cortex
        Mesh = [bias,cortex_size]
        # For biased subcortex mesh
        minsize = 0.025
        maxsize = 0.2 
        # ======================================================
        # Step Parameters
        # ======================================================
        TotalTime = 1
        Mass_Scaling = Mass_Scaling_list[j] # Mass scaling necessary for stable simulations
        Step = [TotalTime, Mass_Scaling]

        # ======================================================
        # Material Parameters
        # ======================================================
        # CorticalMaterial = 'GRAY'
        # SubCortMaterial = 'WHIT'
        gamma_1 = 1.0 # gamma parameter for pushing effect due to progenitors
        T_1 = 0.1 # intial growth time less than which progenitors expand in WM only
        T_2 = 0.45 # time after which pull effect comes into play in WM @ P6
        scaled_thresh = 0.415  # parameter for phi growth rate function -> This is used in vumat_pull.f to calculate scaled_thresh which is the same as \bar{\delta} and \delta in Table 1.
        mu_g = 0.002 # N/mm2 shear modulus for cortex (Gray matter)
        mu_w = 0.001 # N/mm2  shear modulus for subcortex (White matter)
        lambda_w = 0.1 # N/mm2 
        lambda_g = 0.2 # N/mm2 
        G_GM = 1.5 # growth rate constant for cortex in simulation time 0 to 1, which translates to 0.084/day in Table 1
        N_gyri = 4.0 # Number of proliferation zones
        gamma_hat_pull = gamma_hat_pull_list[j] # gamma parameter for pulling effect due to astrocytes
        b_tilde = 0.4 # parameter \tilde{b} in 1/10 mm units
        alpha = 0.4 # standard deviation of gauss growth rate function
        delta = 0.415 # scaled threshold for Gauss function

        density = 1e-9 # N/mm3, For explicit simulation, a very small density is given to make it a quasi-static problem.

        White_Matter_Properties = (mu_w, lambda_w, G_GM, gamma_1, gamma_hat_pull, T_1, T_2,
                                       MajorAxis_W, MinorAxis_W, b_tilde, N_gyri,
                                       alpha, delta)
        Gray_Matter_Properties = (mu_g, lambda_g, G_GM, T_1, MajorAxis_G, MinorAxis_G)

        # ======================================================
        # Boundary Condition to keep KE < 10% of Total Energy
        # ======================================================
        Viscous_Pressure = Viscous_Pressure_list[j]
        # ======================================================
            # Model
        # ======================================================
        ModelName = "QuarterEllipse_pull_gammahat{}".format(Job_Name_list[j])
        PartName = 'Part-1'
        InstanceName = 'Part-1'
        JobName = ModelName
        # ======================================================
        # Call Functions
        # ======================================================
        Create_Quarter_Ellipse(ModelName, PartName, Dimensions)
        Create_Rigid_wall(ModelName)
        Create_Step(ModelName, Step)
        Create_Material(ModelName, density, White_Matter_Properties, Gray_Matter_Properties)
        Create_Section(ModelName, PartName, Dimensions)
        Create_Assembly(ModelName, PartName, InstanceName, Dimensions)
        Create_Sets(ModelName, PartName, Dimensions)
        Create_Boundary_Conditions(ModelName, InstanceName)
        Create_Viscous_Pressure(ModelName, InstanceName, Viscous_Pressure)
        Create_Contact(ModelName, InstanceName)   
        Create_Mesh(ModelName, PartName, InstanceName, Dimensions, Mesh, minsize, maxsize)
        Create_Output(ModelName, PartName, InstanceName)
        Create_Job(ModelName, JobName)
