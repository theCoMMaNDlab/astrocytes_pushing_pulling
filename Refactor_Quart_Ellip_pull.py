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
    """Quarter-ellipse bilayer: cortex (gray matter) around subcortex (white
    matter), extruded through the depth.  Same calls as the original recorded
    macro; every dimension comes from Dimensions.
      1. full outer ellipse, extruded by Depth
      2. partition the face at z = Depth with the +x axis, +y axis and the
         inner (gray/white) ellipse
      3. split the cells with the two axis planes and the extruded interface
      4. remove the faces of the other three quadrants
    The face/edge/vertex indices and masks in steps 2-4 depend on the
    construction order, not on the sizes, as long as the major axis is along x.
    """
    MajorAxis_G, MajorAxis_W, MinorAxis_G, MinorAxis_W, Depth = Dimensions

    # ---- input checks (the recorded index picks assume this layout) ----
    if not (0.0 < MajorAxis_W < MajorAxis_G and 0.0 < MinorAxis_W < MinorAxis_G):
        raise ValueError('Create_Quarter_Ellipse: need 0 < white-matter axes < gray-matter axes, '
                         'got %s' % (Dimensions, ))
    if MinorAxis_G > MajorAxis_G or MinorAxis_W > MajorAxis_W:
        raise ValueError('Create_Quarter_Ellipse: the major axis must lie along x '
                         '(MajorAxis >= MinorAxis); the vertex/edge indices below rely on it')
    if Depth <= 0.0:
        raise ValueError('Create_Quarter_Ellipse: Depth must be > 0, got %g' % Depth)

    # ---- sketch settings: the values of the original recording ----
    # None of them is a dimension of the part, but they are where the sketch
    # constraint solver starts (the +x line end is snapped onto the ellipse), so
    # changing them can move vertices by round-off amounts.  That was enough to
    # change the free subcortex mesh (+32 elements) and make the gamma = 0.05 run
    # abort at t = 0.976.  Keep them as they are.  They only need the part to fit
    # inside the sketch sheet, i.e. axes below ~40 mm.
    SHEET_1 = 100.0     # sheet of the ellipse sketch
    SHEET_2 = 187.48    # sheet of the partition sketch
    GRID    = 4.68      # grid of the partition sketch (display only)
    EXT     = 10.0      # the cutting lines run from the origin to EXT x the axes
    ORIGIN_Z = 2.0      # sketch origin of the recording (projected onto the face)

    cliCommand("""session.journalOptions.setValues(replayGeometry=COORDINATE,recoverGeometry=COORDINATE)""")

    mdb.Model(name=ModelName, modelType=STANDARD_EXPLICIT)
    model = mdb.models[ModelName]

    # ---- 1. full outer ellipse, extruded through the depth ----
    s1 = model.ConstrainedSketch(name='__profile__', sheetSize=SHEET_1)
    s1.setPrimaryObject(option=STANDALONE)
    s1.EllipseByCenterPerimeter(center=(0.0, 0.0),
                                axisPoint1=(MajorAxis_G, 0.0),
                                axisPoint2=(0.0, MinorAxis_G))
    p = model.Part(name=PartName, dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p.BaseSolidExtrude(sketch=s1, depth=Depth)
    s1.unsetPrimaryObject()

    # ---- 2. sketch on the face at z = Depth: +x axis, +y axis, inner ellipse ----
    f, e = p.faces, p.edges
    t = p.MakeSketchTransform(sketchPlane=f[1], sketchUpEdge=e[0],
                              sketchPlaneSide=SIDE1, origin=(0.0, 0.0, ORIGIN_Z))
    s = model.ConstrainedSketch(name='__profile__', sheetSize=SHEET_2,
                                gridSpacing=GRID, transform=t)
    g, v = s.geometry, s.vertices
    s.setPrimaryObject(option=SUPERIMPOSE)
    p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
    s.Line(point1=(EXT * MajorAxis_G, 0.0), point2=(0.0, 0.0))    # along +x
    s.HorizontalConstraint(entity=g[3], addUndoState=False)
    s.PerpendicularConstraint(entity1=g[2], entity2=g[3], addUndoState=False)
    s.Line(point1=(0.0, 0.0), point2=(0.0, EXT * MinorAxis_G))    # along +y
    s.VerticalConstraint(entity=g[4], addUndoState=False)
    s.PerpendicularConstraint(entity1=g[3], entity2=g[4], addUndoState=False)
    s.CoincidentConstraint(entity1=v[2], entity2=g[2], addUndoState=False)
    s.EllipseByCenterPerimeter(center=(0.0, 0.0),                  # gray/white interface
                               axisPoint1=(MajorAxis_W, 0.0),
                               axisPoint2=(0.0, MinorAxis_W))
    pickedFaces = p.faces.getSequenceFromMask(mask=('[#2 ]', ), )   # = f[1], the z = Depth face
    p.PartitionFaceBySketch(sketchUpEdge=p.edges[0], faces=pickedFaces, sketch=s)
    s.unsetPrimaryObject()

    # ---- 3. split the cells: the two axis planes, then the interface ----
    pickedCells = p.cells.getSequenceFromMask(mask=('[#1 ]', ), )
    v1 = p.vertices
    p.PartitionCellByPlaneThreePoints(point1=v1[3], point2=v1[4], point3=v1[5],
                                      cells=pickedCells)

    pickedCells = p.cells.getSequenceFromMask(mask=('[#2 ]', ), )
    v2, e1 = p.vertices, p.edges
    p.PartitionCellByPlaneThreePoints(point1=v2[4], point2=v2[7], cells=pickedCells,
                                      point3=p.InterestingPoint(edge=e1[7], rule=MIDDLE))

    pickedCells = p.cells.getSequenceFromMask(mask=('[#4 ]', ), )
    e = p.edges
    p.PartitionCellByExtrudeEdge(line=e[3], cells=pickedCells, edges=(e[19], ),
                                 sense=FORWARD)

    # ---- 4. remove the faces of the other three quadrants ----
    f1 = p.faces
    p.RemoveFaces(faceList=f1[5:6] + f1[11:13] + f1[14:15] + f1[16:18],
                  deleteCells=False)
    f = p.faces
    p.RemoveFaces(faceList=f[5:6] + f[7:8] + f[12:13], deleteCells=False)


#######################################################################################
def Quarter_Ellipse_Points(Dimensions, phi=np.pi / 4.0):
    """In-plane (x, y) pick points for the quarter-ellipse part.
 
    All come from Dimensions, so they follow any change of axes or cortex
    thickness (CT = MajorAxis_G - MajorAxis_W).  Points on or inside the
    curved boundaries sit at the parametric angle phi, i.e. (a*cos(phi),
    b*sin(phi)), away from the symmetry planes, so none of them lands on an
    edge shared by two faces.  'gm' uses the ellipse halfway through the
    cortex, which is inside the cortex for any thickness.
    Combine with a z level:  pt = xyz(P['gm'], Depth).
    """
    aG, aW, bG, bW, Depth = [float(d) for d in Dimensions]
    cs, sn = np.cos(phi), np.sin(phi)
    aM, bM = 0.5 * (aG + aW), 0.5 * (bG + bW)     # mid-cortex ellipse
    return {
        'wm':        (0.5 * aW * cs, 0.5 * bW * sn),  # inside the subcortex
        'gm':        (aM * cs, bM * sn),              # inside the cortex
        'pial':      (aG * cs, bG * sn),              # on the outer (pial) surface
        'interface': (aW * cs, bW * sn),              # on the gray/white interface
        'wm_y0':     (0.5 * aW, 0.0),                 # y = 0 plane, subcortex part
        'gm_y0':     (aM, 0.0),                       # y = 0 plane, cortex part
        'wm_x0':     (0.0, 0.5 * bW),                 # x = 0 plane, subcortex part
        'gm_x0':     (0.0, bM),                       # x = 0 plane, cortex part
        'tip_x':     (aG, 0.0),                       # outer vertex on the x axis
    }

def xyz(xy, z):
    """(x, y) + z -> a plain-float point for findAt."""
    return (float(xy[0]), float(xy[1]), float(z))

def _find(entities, points, label):
    """findAt for several points.  findAt itself only prints a warning and
    skips a point that hits nothing; this stops the script instead."""
    found = entities.findAt(*[(pt, ) for pt in points])
    if len(found) != len(points):
        raise RuntimeError('%s: %d of %d pick points hit nothing: %s'
                           % (label, len(points) - len(found), len(points), points))
    return found



#######################################################################################
#######################################################################################
def Create_Step(ModelName, Step):
    TotalTime, Mass_Scaling = Step

    mdb.models[ModelName].ExplicitDynamicsStep(name='growth', previous='Initial', timePeriod=TotalTime, massScaling=((SEMI_AUTOMATIC, MODEL, AT_BEGINNING, Mass_Scaling, 0.0, None, 0, 0, 0.0, 0.0, 0, None),),)
    mdb.models[ModelName].steps['growth'].setValues(timePeriod=TotalTime, scaleFactor=1.0, linearBulkViscosity=0.0, quadBulkViscosity=0.0, improvedDtMethod=ON)
    mdb.models[ModelName].TabularAmplitude(name='Amp-1', timeSpan=STEP, smooth=SOLVER_DEFAULT, data=((0.0, 0.0), (TotalTime, 1.0)))


def Create_Material(ModelName, density, White_Matter_Properties, Gray_Matter_Properties):
    mdb.models[ModelName].Material(name='WHIT')
    mdb.models[ModelName].materials['WHIT'].Density(table=((density, ), ))
    mdb.models[ModelName].materials['WHIT'].Depvar(n=10)
    mdb.models[ModelName].materials['WHIT'].UserMaterial(mechanicalConstants=White_Matter_Properties)

    mdb.models[ModelName].Material(name='GRAY')
    mdb.models[ModelName].materials['GRAY'].Density(table=((density, ), ))
    mdb.models[ModelName].materials['GRAY'].Depvar(n=10)
    mdb.models[ModelName].materials['GRAY'].UserMaterial(mechanicalConstants=Gray_Matter_Properties)

#######################################################################################
#######################################################################################
def Create_Section(ModelName, PartName, Dimensions):
    Depth = Dimensions[4]
    P = Quarter_Ellipse_Points(Dimensions)
 
    model = mdb.models[ModelName]
    p = model.parts[PartName]
    c = p.cells
 
    c_subcort = _find(c, [xyz(P['wm'], Depth / 2.)], 'Subcortex cell')
    c_cortex  = _find(c, [xyz(P['gm'], Depth / 2.)], 'Cortex cell')
 
    model.HomogeneousSolidSection(name='Subcortex', material='WHIT', thickness=None)
    region = p.Set(cells=c_subcort, name='Subcortex')
    p.SectionAssignment(region=region, sectionName='Subcortex', offset=0.0,
                        offsetType=MIDDLE_SURFACE, offsetField='',
                        thicknessAssignment=FROM_SECTION)
 
    model.HomogeneousSolidSection(name='Cortex', material='GRAY', thickness=None)
    region = p.Set(cells=c_cortex, name='Cortex')
    p.SectionAssignment(region=region, sectionName='Cortex', offset=0.0,
                        offsetType=MIDDLE_SURFACE, offsetField='',
                        thicknessAssignment=FROM_SECTION)

#######################################################################################
#######################################################################################
def Create_Rigid_wall(ModelName, Dimensions, reach=2.5, margin=0.5, width=0.2):
    """One flat rigid wall, used twice (wall-1 on x = 0, wall-2 on y = 0).
    In its own sketch it runs along y from -margin to reach * (largest outer
    axis) and is `width` wide, centred on x = 0, so both instances get exactly
    the same length and Create_Assembly only rotates and shifts them.
      reach  : how far the wall runs along its symmetry plane, in units of the
               largest outer axis.  The pial ends move out about 2x by t = 1,
               so 2.5 leaves room (9.0 mm for the 3.6 mm axis).
      margin : how far it runs past the origin, mm
      width  : out-of-plane width, mm; must be larger than Depth
    """
    MajorAxis_G, MajorAxis_W, MinorAxis_G, MinorAxis_W, Depth = Dimensions
    L_end = reach * max(MajorAxis_G, MinorAxis_G)
    if width <= Depth:
        raise ValueError('Create_Rigid_wall: width (%g) must exceed Depth (%g)' % (width, Depth))

    model = mdb.models[ModelName]
    s1 = model.ConstrainedSketch(name='__profile__', sheetSize=max(100.0, 4.0 * L_end))
    s1.setPrimaryObject(option=STANDALONE)
    s1.rectangle(point1=(-0.5 * width, -margin), point2=(0.5 * width, L_end))
    p = model.Part(name='wall', dimensionality=THREE_D, type=DISCRETE_RIGID_SURFACE)
    p.BaseShell(sketch=s1)
    s1.unsetPrimaryObject()
    p.ReferencePoint(point=(0.0, 0.5 * (L_end - margin), 0.0))      # centre of the wall


#######################################################################################
#######################################################################################
def Create_Assembly(ModelName, PartName, InstanceName, Dimensions, gap=0.01):
    """Brain instance plus the two walls.  Each wall sits `gap` outside its
    symmetry plane, centred on mid-depth, with its contact side (the sketch
    normal, side 1) facing the brain:
      wall-1 : plane x = -gap, runs along +y, normal +x
      wall-2 : plane y = -gap, runs along +x, normal +y
    Both cover the same range, -margin ... reach * largest axis.
    """
    Depth = Dimensions[4]
    model = mdb.models[ModelName]
    a = model.rootAssembly
    a.DatumCsysByDefault(CARTESIAN)

    a.Instance(name=InstanceName, part=model.parts[PartName], dependent=ON)
    a.Instance(name='wall-1', part=model.parts['wall'], dependent=ON)
    a.Instance(name='wall-2', part=model.parts['wall'], dependent=ON)

    O = (0.0, 0.0, 0.0)
    # wall-1: its length (sketch y) stays along +y, its normal (+z) turns to +x
    a.rotate(instanceList=('wall-1', ), axisPoint=O, axisDirection=(0.0, 1.0, 0.0), angle=90.0)
    a.translate(instanceList=('wall-1', ), vector=(-gap, 0.0, 0.5 * Depth))
    # wall-2: its length turns to +x, then its normal (+z) turns to +y
    a.rotate(instanceList=('wall-2', ), axisPoint=O, axisDirection=(0.0, 0.0, 1.0), angle=-90.0)
    a.rotate(instanceList=('wall-2', ), axisPoint=O, axisDirection=(1.0, 0.0, 0.0), angle=-90.0)
    a.translate(instanceList=('wall-2', ), vector=(0.0, -gap, 0.5 * Depth))

def y_coord(x_coord,majoraxis,minoraxis):
    output = np.sqrt((1 - (x_coord*x_coord)/(majoraxis*majoraxis)))*minoraxis
    return output

def x_coord(y_coord,majoraxis,minoraxis):
    output = np.sqrt((1 - (y_coord*y_coord)/(minoraxis*minoraxis)))*majoraxis
    return output


####################################################################################### 
####################################################################################### 
def Create_Sets(ModelName, PartName, Dimensions):
    Depth = Dimensions[4]
    P = Quarter_Ellipse_Points(Dimensions)
    z0, zm, z1 = 0.0, Depth / 2., Depth       # back face, mid-depth, front face
 
    model = mdb.models[ModelName]
    p = model.parts[PartName]
    pwall = model.parts['wall']
    f, e, c, v = p.faces, p.edges, p.cells, p.vertices
 
    # ---- vertex ----
    verts = _find(v, [xyz(P['tip_x'], z1)], 'node')
 
    # ---- faces ----
    f_back       = _find(f, [xyz(P['wm'], z0), xyz(P['gm'], z0)], 'back')         # z = 0
    f_bottom     = _find(f, [xyz(P['wm_y0'], zm), xyz(P['gm_y0'], zm)], 'bottom') # y = 0
    f_side1      = _find(f, [xyz(P['wm_x0'], zm), xyz(P['gm_x0'], zm)], 'side-1') # x = 0
    f_side2      = f_bottom                                                       # y = 0 (same faces)
    f_cortex_top = _find(f, [xyz(P['pial'], zm)], 'cortex_top')                   # pial surface
    f_interface  = _find(f, [xyz(P['interface'], zm)], 'interface face')          # gray/white interface
    s_cortex     = _find(f, [xyz(P['gm'], z1)], 'cortex_front')                   # z = Depth, cortex
    s_subcortex  = _find(f, [xyz(P['wm'], z1)], 'subcortex_front')                # z = Depth, subcortex
 
    # ---- edges (curves at z = Depth) ----
    e_top         = _find(e, [xyz(P['pial'], z1)], 'top-edge')
    e_interface   = _find(e, [xyz(P['interface'], z1)], 'interface edge')
    e_interface_2 = e_interface
 
    # ---- surfaces ----
    s_top       = f_cortex_top
    s_side1     = f_side1
    s_side2     = f_side2
    s_interface = f_interface
    s_wall      = pwall.faces            # the wall part is a single face
    if len(s_wall) != 1:
        raise RuntimeError('wall part has %d faces, expected 1' % len(s_wall))
 
    # ---- cells ----
    cells = _find(c, [xyz(P['wm'], zm), xyz(P['gm'], zm)], 'whole-domain')
 
    # ---- face sets ----
    p.Set(faces=f_back, name='back')
    p.Set(faces=f_bottom, name='bottom')
    p.Set(faces=f_side1, name='side-1')
    p.Set(faces=f_side2, name='side-2')
    p.Set(faces=s_cortex, name='cortex_front')
    p.Set(faces=s_subcortex, name='subcortex_front')
    p.Set(faces=f_cortex_top, name='cortex_top')
    # The old script also named this face set 'interface' and then overwrote it
    # with the edge set below.  It now has its own name; 'interface' stays the
    # edge set, as before.
    p.Set(faces=f_interface, name='interface-face')
 
    # ---- edge sets ----
    p.Set(edges=e_top, name='top-edge')
    p.Set(edges=e_interface, name='interface')
    p.Set(edges=e_interface_2, name='interface-edge')
 
    # ---- surfaces ----
    p.Surface(side1Faces=s_top, name='top_s')
    p.Surface(side1Faces=s_cortex, name='cortex_s')
    p.Surface(side1Faces=s_subcortex, name='subcortex_s')
    p.Surface(side1Faces=s_side1, name='side1_s')
    p.Surface(side1Faces=s_side2, name='side2_s')
    p.Surface(side1Faces=s_interface, name='interface_s')
    pwall.Surface(side1Faces=s_wall, name='s_wall')
 
    # ---- cell set ----
    p.Set(cells=cells, name='whole-domain')
 
    # ---- vertex set ----
    p.Set(vertices=verts, name='node')
 
    # ---- rigid-wall reference point (the wall part has exactly one) ----
    refPoints = tuple(pwall.referencePoints.values())
    if len(refPoints) != 1:
        raise RuntimeError('wall part has %d reference points, expected 1' % len(refPoints))
    pwall.Set(referencePoints=refPoints, name='ref')
 


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
    # region = a.instances['Part-1'].sets['whole-domain']
    # mdb.models[ModelName].DisplacementBC(name='fix_z', 
    #     createStepName='Initial', region=region, u1=UNSET, u2=UNSET, u3=SET, 
    #     ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude='Amp-1', distributionType=UNIFORM, 
    #     fieldName='', localCsys=None)
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
    mdb.models[ModelName].interactionProperties['IntProp-1'].GeometricProperties(contactArea=1.0, padThickness=0.02)

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
def Create_Z_Mesh(ModelName, PartName, InstanceName, Dimensions, Mesh, minsize, maxsize,
                  wall_seed=0.15, length_ratio=0.30):
    """Mesh the quarter ellipse with one element through the depth; same seeds
    as the original script, with every pick point computed from Dimensions.
      pial curve (z = 0)                  -> size minsize; the structured cortex
                                             copies this count to the interface
      cortex edges on x = 0 and y = 0     -> size cortex_size (= CT/ecortex)
      subcortex edges on x = 0 and y = 0  -> biased, minsize at the interface,
                                             maxsize at the centre
      edges through the depth             -> exactly 1 element
    After meshing, the node spacing on each subcortex edge is checked; an edge
    whose fine end is not at the interface is flipped and the part remeshed.
    (Abaqus replaces Python's sum() with its own, so sum() is not used.)"""
    aG, aW, bG, bW, D = Dimensions
    cortex_size = Mesh[1]
    k = 0.5 ** 0.5                                 # points on the curves at 45 deg
    aM, bM = 0.5 * (aG + aW), 0.5 * (bG + bW)      # ellipse halfway through the cortex
    zs = (0.0, D)
    model = mdb.models[ModelName]
    inst = model.rootAssembly.instances[InstanceName]

    def pick(arr, pts, what):                      # findAt that stops if a point hits nothing
        found = arr.findAt(*[((float(x), float(y), float(z)), ) for x, y, z in pts])
        if len(found) != len(pts):
            raise RuntimeError('Create_Z_Mesh: %s: %d of %d points hit nothing'
                               % (what, len(pts) - len(found), len(pts)))
        return found

    # ---- rigid wall (flat, so any element size works) ----
    pw = model.parts['wall']
    pw.seedPart(size=wall_seed, deviationFactor=0.1, minSizeFactor=0.1)
    pw.generateMesh()

    p = model.parts[PartName]
    e, c = p.edges, p.cells

    # ---- pial curve, back face only (z = 0), exactly as the original script ----
    # The sweep copies it to z = Depth.  Seeding the z = Depth curve as well gave
    # a different free mesh in the subcortex (+32 elements).
    pial = pick(e, [(aG * k, bG * k, 0.0)], 'pial curve z=0')
    p.seedEdgeBySize(edges=pial, size=minsize, deviationFactor=0.1, constraint=FINER)

    # ---- across the cortex ----
    gm = pick(e, [(aM, 0.0, z) for z in zs] + [(0.0, bM, z) for z in zs], 'cortex radial edges')
    p.seedEdgeBySize(edges=gm, size=cortex_size, deviationFactor=0.1, constraint=FINER)

    # ---- through the subcortex, biased ----
    # (point on edge, interface end, end1/end2 group used by the original script)
    wm = [((0.5 * aW, 0.0, 0.0), (aW, 0.0, 0.0), 2), ((0.5 * aW, 0.0, D), (aW, 0.0, D), 1),
          ((0.0, 0.5 * bW, 0.0), (0.0, bW, 0.0), 1), ((0.0, 0.5 * bW, D), (0.0, bW, D), 2)]
    wm_edges = [pick(e, [pt], 'subcortex radial edge')[0] for pt, _, _ in wm]
    group = [g for _, _, g in wm]

    def seed_bias():
        ends = {}
        for g, key in ((1, 'end1Edges'), (2, 'end2Edges')):
            sel = tuple([ed for ed, gi in zip(wm_edges, group) if gi == g])
            if sel:
                ends[key] = sel
        p.seedEdgeByBias(biasMethod=SINGLE, minSize=minsize, maxSize=maxsize,
                         constraint=FINER, **ends)

    def wrong_bias():                              # edges whose fine end is not at the interface
        bad = []
        for i, ed in enumerate(wm_edges):
            fine = wm[i][1]
            xs = sorted([n.coordinates for n in ed.getNodes()],
                        key=lambda q: (q[0] - fine[0]) ** 2 + (q[1] - fine[1]) ** 2)
            h = lambda a_, b_: ((xs[a_][0] - xs[b_][0]) ** 2 + (xs[a_][1] - xs[b_][1]) ** 2) ** 0.5
            ok = len(xs) < 3 or h(0, 1) < h(-1, -2)
            print('  subcortex edge to (%.3f, %.3f) z=%g: %d elements, %.4f at interface, '
                  '%.4f at centre  %s' % (fine[0], fine[1], fine[2], len(xs) - 1,
                  h(0, 1), h(-1, -2), 'ok' if ok else '<-- wrong direction'))
            if not ok:
                bad.append(i)
        return bad

    seed_bias()

    # ---- exactly one element through the depth ----
    tol, zEdges = 1.0e-3 * D, []
    for edg in p.edges:
        vid = edg.getVertices()
        if len(vid) == 2:
            q1, q2 = p.vertices[vid[0]].pointOn[0], p.vertices[vid[1]].pointOn[0]
            if abs(q1[0] - q2[0]) < tol and abs(q1[1] - q2[1]) < tol and abs(abs(q1[2] - q2[2]) - D) < tol:
                zEdges.append(edg)
    print('z-direction edges found: %d (expected 5)' % len(zEdges))
    p.seedEdgeByNumber(edges=tuple(zEdges), number=1, constraint=FIXED)

    # ---- mesh controls, element types ----
    c_wm  = pick(c, [(0.5 * aW * k, 0.5 * bW * k, 0.5 * D)], 'subcortex cell')
    c_all = pick(c, [(0.5 * aW * k, 0.5 * bW * k, 0.5 * D), (aM * k, bM * k, 0.5 * D)], 'cells')
    p.deleteMesh(regions=c_wm)
    p.setMeshControls(regions=c_wm, technique=SWEEP, algorithm=ADVANCING_FRONT)
    p.setElementType(regions=(c_all, ), elemTypes=(
        mesh.ElemType(elemCode=C3D8R, elemLibrary=EXPLICIT, kinematicSplit=AVERAGE_STRAIN,
                      secondOrderAccuracy=OFF, hourglassControl=ENHANCED,
                      distortionControl=ON, lengthRatio=length_ratio),
        mesh.ElemType(elemCode=C3D6, elemLibrary=EXPLICIT),
        mesh.ElemType(elemCode=C3D4, elemLibrary=EXPLICIT)))

    # ---- mesh, check the bias, fix it if needed ----
    p.generateMesh(regions=(inst, ))
    try:
        bad = wrong_bias()
        if bad:
            print('flipping the bias on %d subcortex edge(s) and remeshing' % len(bad))
            for i in bad:
                group[i] = 3 - group[i]
            p.deleteMesh(regions=c_all)
            seed_bias()
            p.generateMesh(regions=(inst, ))
            if wrong_bias():
                print('WARNING: subcortex bias still wrong after flipping')
    except Exception as exc:
        print('bias check skipped (%s)' % exc)

    try:
        n_pial = len(pick(e, [(aG * k, bG * k, 0.0)], 'pial')[0].getNodes()) - 1
        n_intf = len(pick(e, [(aW * k, bW * k, 0.0)], 'interface')[0].getNodes()) - 1
        print('divisions: pial %d, interface %d' % (n_pial, n_intf))
    except Exception:
        pass
    print('elements generated: %d, nodes: %d' % (len(p.elements), len(p.nodes)))


#######################################################################################
#######################################################################################
def Create_Output(ModelName, PartName, InstanceName):
    p = mdb.models[ModelName].parts[PartName]

    # ---- field output (unchanged) ----
    regionDef = mdb.models[ModelName].rootAssembly.allInstances[InstanceName].sets['whole-domain']
    mdb.models[ModelName].fieldOutputRequests['F-Output-1'].setValues(
        variables=('S','U','V','LE','SDV','CSTRESS','CDISP','ELASE','IVOL'),
        numIntervals=100)
    mdb.models[ModelName].FieldOutputRequest(name='F-Output-2', createStepName='growth',
        variables=('S','PE','LE','U','ELSE','COORD','SDV','EVOL'),
        numIntervals=100, region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)
    regionDef = mdb.models[ModelName].rootAssembly.allInstances[InstanceName].sets['Cortex']
    mdb.models[ModelName].FieldOutputRequest(name='F-Output-3', createStepName='growth',
        variables=('ELSE','COORD','U','EVOL'),
        numIntervals=100, region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)
    regionDef = mdb.models[ModelName].rootAssembly.allInstances[InstanceName].sets['Subcortex']
    mdb.models[ModelName].FieldOutputRequest(name='F-Output-4', createStepName='growth',
        variables=('ELSE','COORD','U','EVOL'),
        numIntervals=100, region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)

    # ---- whole-model energy balance ----
    mdb.models[ModelName].HistoryOutputRequest(name='H-Energy',
        createStepName='growth',
        variables=('ALLIE','ALLKE','ALLSE','ALLAE','ALLDC','ALLVD','ALLWK','ALLCD','ETOTAL'),
        numIntervals=1000, region=MODEL)

    # ---- per-region energies (optional, for cortex vs subcortex split) ----
    regionDef = mdb.models[ModelName].rootAssembly.allInstances[InstanceName].sets['whole-domain']
    mdb.models[ModelName].HistoryOutputRequest(name='H-Output-1',
        createStepName='growth', variables=('ALLIE','ALLKE','ALLAE'),
        numIntervals=1000, region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)

#######################################################################################
#######################################################################################
def Create_Job(ModelName, JobName):
    mdb.Job(name=JobName, model=ModelName, description='', type=ANALYSIS,
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,
        memoryUnits=PERCENTAGE, explicitPrecision=DOUBLE,
        nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF,
        contactPrint=OFF, historyPrint=OFF,
        userSubroutine='vumat_pull_sigm.f',
        scratch='', resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN,
        numDomains=16, activateLoadBalancing=False, numThreadsPerMpiProcess=1,
        multiprocessingMode=DEFAULT, numCpus=16)


    mdb.jobs[JobName].writeInput(consistencyChecking=OFF)
    try: 
        mdb.jobs[JobName].submit(consistencyChecking=OFF)
        mdb.jobs[JobName].waitForCompletion()
    except: 
        print("Job %s crashed" % JobName)

#######################################################################################
#######################################################################################


if __name__ == '__main__':
    gamma_hat_pull_list = [0.05,0.5,1.0,2.0,3.0]
    Mass_Scaling_list = [120,120,120,120,120]
    Viscous_Pressure_list = [1e-4,1e-4,1e-4,1e-4,1e-4]
    Job_Name_list = ['5em2','50em2','1','2','3']  #len(gamma_hat_pull_list)
    for j in range(len(gamma_hat_pull_list)):

        # ======================================================
        # Dimensions
        # ======================================================
        CT = 0.3 # Units: 1/10 mm
        MajorAxis_G = 3.6 # Units: 1/10 mm
        MinorAxis_G = 3.0 # Units: 1/10 mm
        MajorAxis_W = MajorAxis_G - CT
        MinorAxis_W = MinorAxis_G - CT

        Depth = 0.02 # Units: 1/10 mm
        Dimensions = [MajorAxis_G, MajorAxis_W, MinorAxis_G, MinorAxis_W, Depth]
        # ======================================================
        # Mesh Parameters
        # ======================================================
        bias = 3                                                        # bias through width of subcortex
        ecortex = 7                                                     # number of elements in cortex                                                
        cortex_size = CT/ecortex # size of element in cortex
        Mesh = [bias,cortex_size]
        # For biased subcortex mesh
        minsize = 0.02
        maxsize = 0.2 
        # ======================================================
        # Step Parameters
        # ======================================================
        DAYS = 18.0
        TotalTime = DAYS
        Mass_Scaling = Mass_Scaling_list[j] # Mass scaling necessary for stable simulations
        Step = [TotalTime, Mass_Scaling]

        # ======================================================
        # Material Parameters
        # ======================================================
        # CorticalMaterial = 'GRAY'
        # SubCortMaterial = 'WHIT'
        gamma_1 = 1.25 # gamma parameter for pushing effect due to progenitors
        T_1 = 2.0 # intial growth time less than which progenitors expand in WM only
        T_2 = 8.0 # time after which pull effect comes into play in WM @ P6
        scaled_thresh = 0.415  # parameter for phi growth rate function -> This is used in vumat_pull.f to calculate scaled_thresh which is the same as \bar{\delta} and \delta in Table 1.
        mu_g = 0.002 # N/mm2 shear modulus for cortex (Gray matter)
        mu_w = 0.001 # N/mm2  shear modulus for subcortex (White matter)
        lambda_w = 0.1 # N/mm2 
        lambda_g = 0.2 # N/mm2 
        G_GM = 2.0/DAYS # growth rate constant for cortex in simulation time 0 to 1, which translates to 0.084/day in Table 1
        N_gyri = 4.0 # Number of proliferation zones
        gamma_hat_pull = gamma_hat_pull_list[j] # gamma parameter for pulling effect due to astrocytes
        b_tilde = 0.4 # parameter \tilde{b} in mm units
        alpha = 0.4 # standard deviation of gauss growth rate function
        delta = 0.415 # scaled threshold for Gauss function

        density = 1e-9*DAYS**2 # N/mm3, For explicit simulation, a very small density is given to make it a quasi-static problem.

        White_Matter_Properties = (mu_w, lambda_w, G_GM, gamma_1, gamma_hat_pull, T_1, T_2,
                                       MajorAxis_W, MinorAxis_W, b_tilde, N_gyri,
                                       alpha, delta)
        Gray_Matter_Properties = (mu_g, lambda_g, G_GM, T_1, MajorAxis_G, MinorAxis_G)

        # ======================================================
        # Boundary Condition to keep KE < 10% of Total Energy
        # ======================================================
        Viscous_Pressure = Viscous_Pressure_list[j]*DAYS
        # ======================================================
            # Model
        # ======================================================
        ModelName = "Logis_pull_Z1_{}".format(Job_Name_list[j])
        PartName = 'Part-1'
        InstanceName = 'Part-1'
        JobName = ModelName
        # ======================================================
        # Call Functions
        # ======================================================
        Create_Quarter_Ellipse(ModelName, PartName, Dimensions)
        Create_Rigid_wall(ModelName,Dimensions)
        Create_Step(ModelName, Step)
        Create_Material(ModelName, density, White_Matter_Properties, Gray_Matter_Properties)
        Create_Section(ModelName, PartName, Dimensions)
        Create_Assembly(ModelName, PartName, InstanceName, Dimensions)
        Create_Sets(ModelName, PartName, Dimensions)
        Create_Boundary_Conditions(ModelName, InstanceName)
        Create_Viscous_Pressure(ModelName, InstanceName, Viscous_Pressure)
        Create_Contact(ModelName, InstanceName)   
        Create_Z_Mesh(ModelName, PartName, InstanceName, Dimensions, Mesh, minsize, maxsize)
        Create_Output(ModelName, PartName, InstanceName)
        Create_Job(ModelName, JobName)
