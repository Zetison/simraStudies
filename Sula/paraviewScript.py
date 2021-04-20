from os.path import expanduser
import sys
import numpy as np
from scipy.spatial.transform import Rotation as Rot
import copy
#from pyproj import Transformer
from matplotlib import cm
#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()
 
home = expanduser("~")
sys.path.insert(1, home+'/kode/paraUtils')
from sources import coffeeFilter
from utils import *
ImportPresets(filename=home+'/kode/colormaps/SINTEF1.json')
ImportPresets(filename=home+'/kode/colormaps/IPPC.json')
ImportPresets(filename=home+'/kode/colormaps/google.json')
#import splipy
#print(sys.executable)
#print(sys.path)
#pluginsPath = home+'/programs/paraview_build/lib/paraview-5.8/plugins/'
#pluginsPath =  '/usr/lib/ParaView-5.8.1-osmesa-MPI-Linux-Python3.7-64bit/lib/paraview-5.8/plugins/'
#LoadPlugin(pluginsPath+'/SurfaceLIC/SurfaceLIC.so', remote=False, ns=globals())
fileName = 'SED_fileName'
caseName = 'SED_caseName'
outputPath = home+'/results/simra/'+caseName+'/'
topoRes = 'SED_topoRes' 
#topoRes = '10m'
topologyFileName = outputPath+caseName+topoRes+'.vts'
textureFileName_topo = outputPath+caseName+topoRes+'_topo.png'
textureFileName_NIB = outputPath+caseName+topoRes+'.png'
windCase = 1
noSteps = 201
finalTime = 11.7006
sqrtTKE_max = 5.0
h_max = 4000
height_max = 1600 #max height of height colormap
angleOfAttack_West = 3.0
angleOfAttack_East = 3.0
runwayWidth = 150.0
runwayHeight = 30.0
runwayHeight = 10.0
bridgeBaseHeight=50.0
plotRunway               = 0 
plotTakeOffLines         = 0 
plotBridgeAndSensors     = 1 
plotUniverseBackground   = 0 

plotLIC                  = 1 
plotStreamLines          = 1 
plotVolumeRendering      = 1 
plotIPPCmapsHorizontal   = 1 
plotIPPCmapsVertical     = 1 
plotOverTime             = 1

makeVideo                = 0
saveScreenShots          = 1
useTransparentBackground = 0

plotError                = plotLIC 
bridge = SED_BRIDGE
bridgeHeight = 20
# Use i.e. norgeskart.no and "Vis koordinater" to get UTM coordinates (NORD,ØST,height)
#runwayEndWest = latLonUTM("623323.41","0060532.89",35.3)
#runwayEndEast = latLonUTM("623351.26","0060740.31",49.2)
mastLoc = np.genfromtxt('mastLoc.csv', delimiter=',')
if bridge == 1:
    runwayEndWest = mastLoc[0][:3]
    runwayEndEast = mastLoc[1][:3]
elif bridge == 2:
    runwayEndWest = mastLoc[2][:3]
    runwayEndEast = mastLoc[3][:3]

#runwayEndWest = np.array([346520,6920740,4.4])
#runwayEndEast = np.array([351140,6922074,3.6])
#runwayEndWest = np.array([549213-0.4, 6943810+0.905,0.0])
#runwayEndEast = np.array([549213,6943810,0.0])
SAVE_HIST = 20

viewSize = [1920, 1080]
scalarBarLength = 0.26
streamLines_z = 10
horizontalSlice_z = 10
glyphScale1 = 30
if windCase == 1:
    glyphScale2 = 0.03

timeStepAnnotation       = makeVideo
# Extract UTM-33 coordinates from rwyCoords.png file (available at https://ais.avinor.no/no/AIP/View/95/2020-11-05-AIRAC/html/index-no-NO.html). The coordinates in rwyCoords.png are given in Lat-lon Degrees Minutes seconds format. Here, "THR COORD" is the start of the RWY (runway), and "RWY end COORD" is its end
runwayColor      = [0.0,0.0,0.0]
takeOffLineColor = [0.0,0.0,0.0]

runwayColor      = [1.0,0.0,0.0]

P_sw = np.array([25000,6927800,streamLines_z])
P_se = np.array([64999.58,6927800,streamLines_z])
P_nw = np.array([25000,6967800,streamLines_z])
P_ne = np.array([64999.5156,6967800.5,streamLines_z])
if windCase == 1:
    streamTracerPoint1 = P_nw
    streamTracerPoint2 = P_se

runwayDir = runwayEndEast-runwayEndWest
runwayDir = runwayDir/np.linalg.norm(runwayDir)
runwayNormal = np.cross(runwayDir,[0,0,1])
runwayNormal = runwayNormal/np.linalg.norm(runwayNormal)
runwayCenter = (runwayEndEast+runwayEndWest)/2
runwayLength = np.linalg.norm(runwayEndWest-runwayEndEast)
runwayRotVector = np.array([0,0,np.degrees(np.arctan2(runwayDir[1],runwayDir[0]))])

# get animation scene
animationScene1 = GetAnimationScene()
LoadPalette(paletteName='WhiteBackground')
renderView1 = GetActiveViewOrCreate('RenderView')
renderView1.OrientationAxesVisibility = 0

if plotUniverseBackground:
    #convert -crop 50%x100% lightBlueMarble2.jpg converted_lBM2.jpg
    import universeBackground 

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# create a new 'XML Unstructured Grid Reader'
simraPVDresults = PVDReader(registrationName=fileName, FileName=outputPath+fileName+'.pvd')
#simraPVDresults = PVDReader(registrationName=fileName, FileName=outputPath+'cont_met.pvd')
animationScene1.UpdateAnimationUsingDataTimeSteps()
animationScene1.GoToLast()
# update animation scene based on data timesteps

# get active view
layout1 = GetLayoutByName("Layout #1")

# create a new 'Calculator'
calculator1 = Calculator(registrationName='Calculator1', Input=simraPVDresults)
calculator1.Function = 'sqrt(tk)'
calculator1.ResultArrayName = 'sqrtTKE'
calculatorIPPC = Calculator(registrationName='CalculatorIPPC', Input=calculator1)
calculatorIPPC.Function = 'sqrt(tk)'
calculatorIPPC.ResultArrayName = 'sqrtTKE_IPPC'

sqrtTKELUT = GetColorTransferFunction('sqrtTKE')
sqrtTKEPWF = GetOpacityTransferFunction('sqrtTKE')
sqrtTKE_IPPCLUT = GetColorTransferFunction('sqrtTKE_IPPC')

slice1 = Slice(registrationName='Slice1', Input=calculatorIPPC)
slice1.SliceType = 'Plane'
slice1.SliceType.Normal = runwayNormal 
slice1.SliceType.Origin = runwayCenter
Hide3DWidgets(proxy=slice1.SliceType)

proj_u = Calculator(registrationName='Projected u', Input=calculatorIPPC)
proj_u.Function = 'u_X*iHat+u_Y*jHat+u_Z*kHat-(u_X*'+str(runwayNormal[0])+'+u_Y*'+str(runwayNormal[1])+'+u_Z*'+str(runwayNormal[2])+')*('+str(runwayNormal[0])+'*iHat+'+str(runwayNormal[1])+'*jHat+'+str(runwayNormal[2])+'*kHat)'
proj_u.ResultArrayName = 'u_proj'

planeSlice = Plane(registrationName='Plane slice')
planeSliceLength = 43e3
planeSliceHeight = h_max
planeSlice.Origin = np.append(runwayCenter[:2]-runwayDir[:2]*15e3,np.zeros(1))
planeSlice.Point1 = planeSlice.Origin + planeSliceLength*np.append(runwayDir[:2],np.zeros(1))
planeSlice.Point2 = planeSlice.Origin + planeSliceHeight*np.array([0,0,1])
planeSlice.XResolution = np.round(planeSliceLength/1500).astype('int')
planeSlice.YResolution = np.round(10*planeSliceHeight/1500).astype('int')
resampledSlice = ResampleWithDataset(registrationName='proj_u Resampled', SourceDataArrays=proj_u, DestinationMesh=planeSlice)
resampledSlice.CellLocator = 'Static Cell Locator'

topology = XMLStructuredGridReader(registrationName='Topology',FileName=topologyFileName)
topology.PointArrayStatus = ['TCoords_']
layout1 = GetLayout()

calculator2 = Calculator(registrationName='Calculator2', Input=topology)
calculator2.Function = 'coordsZ'
calculator2.ResultArrayName = 'height'
heightLUT = GetColorTransferFunction('height')
heightPWF = GetOpacityTransferFunction('height')

sliceTopology = Slice(registrationName='Slice topology', Input=calculator2)
sliceTopology.SliceType = 'Plane'
sliceTopology.SliceType.Normal = runwayNormal 
sliceTopology.SliceType.Origin = runwayCenter
Hide3DWidgets(proxy=sliceTopology.SliceType)

def rotTran(obj,rot,tran,regName='Intermediate transform',rotateFirst=True):
    if isinstance(obj,(list,np.ndarray)):
        if rotateFirst:
            R = Rot.from_euler('z', rot[2], degrees=True)
            return tran+R.apply(obj)
        else:
            R = Rot.from_euler('z', rot[2], degrees=True)
            return R.apply(obj+tran)
    else:
        if rotateFirst:
            transform = Transform(registrationName=regName,Input=obj)
            transform.TransformAllInputVectors = True
            transform.Transform = 'Transform'
            transform.Transform.Translate = tran
            transform.Transform.Rotate = rot 
        else:
            objTemp = Transform(registrationName='Intermediate translation',Input=obj)
            objTemp.TransformAllInputVectors = True
            objTemp.Transform = 'Transform'
            objTemp.Transform.Translate = tran
            transform = Transform(registrationName=regName,Input=objTemp)
            transform.TransformAllInputVectors = True
            transform.Transform.Rotate = rot 
    return transform

box = Box()
box.Center = [0.0,0.0,0.0]
box.XLength = runwayLength 
box.YLength = runwayWidth
box.ZLength = runwayHeight
runway = rotTran(box,runwayRotVector,runwayCenter,'Runway')

takeOffLineEast = PolyLineSource()
takeOffLineEastLength = 22e3
takeOffEast = copy.deepcopy(runwayDir)
takeOffEast[2] = np.tan(np.radians(angleOfAttack_East))
takeOffLineEast.Points = np.append(runwayEndEast,runwayEndEast+takeOffLineEastLength*takeOffEast,axis=0)

takeOffLineWest = PolyLineSource()
takeOffLineWestLength = 13.9e3
takeOffWest = -copy.deepcopy(runwayDir)
takeOffWest[2] = np.tan(np.radians(angleOfAttack_East))
takeOffLineWest.Points = np.append(runwayEndWest,runwayEndWest+takeOffLineWestLength*takeOffWest,axis=0)
def showRunway(renderView,runway=runway):
    if plotRunway:
        # Display runway and takeoff visualization
        runwayDisplay = Show(runway, renderView, 'GeometryRepresentation')
        runwayDisplay.Representation = 'Surface'
        runwayDisplay.DiffuseColor = runwayColor 

def annotateTimeStep(obj,renderview,location='UpperLeftCorner'):
    pythonAnnotation = PythonAnnotation(registrationName='Time annotation', Input=obj)
    pythonAnnotation.ArrayAssociation = 'Point Data'
    pythonAnnotation.Expression = '"Step: %d\\nTime: %0.2fs" % ((time_index+1)*'+str(SAVE_HIST)+', time_value)'
    pythonAnnotationDisplay = Show(pythonAnnotation, renderview, 'TextSourceRepresentation')
    pythonAnnotationDisplay.WindowLocation = location 
    pythonAnnotationDisplay.FontSize = 5

labels = ['SulaNW (Kvitneset) - 1','SulaNW (Kvitneset) - 2','SulaNW (Kvitneset) - 3',
          'SulaNE (Trælbodneset) - 1','SulaNE (Trælbodneset) - 2','SulaNE (Trælbodneset) - 3',
          'SulaSW (Langeneset) - 1','SulaSW (Langeneset) - 2','SulaSW (Langeneset) - 3','SulaSW (Langeneset) - 4',
          'SulaSE (Kårsteinen) - 1','SulaSE (Kårsteinen) - 2','SulaSE (Kårsteinen) - 3',
          'Midspan bridge 1', 'Midspan bridge 2']
mastLabel = ['SulaNW (Kvitneset)', 'SulaNE (Trælboneset)', 'SulaSW (Langeneset)', 'SulaSE (Kårsteinen)']
points = np.genfromtxt('POINT_FILE', delimiter=' ', skip_header=True)
noMasts = mastLoc.shape[0]
noPoints = points.shape[0]
noBridges = 2
idxMap = [0,0,0,1,1,1,2,2,2,2,3,3,3]
pointSources = [''] * noPoints
lineSource2 = [''] * (noPoints-2)
for i in range(noPoints):
    pointSources[i] = PointSource(registrationName=labels[i])
    pointSources[i].Center = points[i]
    Hide3DWidgets(proxy=pointSources[i])
    if i < noPoints-2:
        lineSource2[i] = Line(registrationName=labels[i])
        lineSource2[i].Point1 = [mastLoc[idxMap[i]][0],mastLoc[idxMap[i]][1],points[i][2]]
        lineSource2[i].Point2 = points[i]
        lineSource2[i].Resolution = 1
        Hide3DWidgets(proxy=lineSource2[i])
splineSource = [''] * noBridges
bridgeNames = ['Bridge north','Bridge south']
for i in range(noBridges):
    splineSource[i] = SplineSource(registrationName=bridgeNames[i])
    splineSource[i].ParametricFunction.Points = [mastLoc[2*i,0],mastLoc[2*i,1],bridgeBaseHeight,points[i-2,0],points[i-2,1],points[i-2,2],mastLoc[2*i+1,0],mastLoc[2*i+1,1],bridgeBaseHeight]
    Hide3DWidgets(proxy=splineSource[i])

lineSource = [''] * noMasts
for i in range(noMasts):
    lineSource[i] = Line(registrationName=mastLabel[i])
    lineSource[i].Point1 = [mastLoc[i][0],mastLoc[i][1],mastLoc[i][2]]
    lineSource[i].Point2 = [mastLoc[i][0],mastLoc[i][1],mastLoc[i][2]+mastLoc[i][3]]
    lineSource[i].Resolution = 1
    Hide3DWidgets(proxy=lineSource[i])
    
def showMasts(renderView):
    if plotBridgeAndSensors:
        pointSourceDisplay = [''] * noPoints
        lineSource2Display = [''] * (noPoints-2)
        for i in range(0,noPoints):
            pointSourceDisplay[i] = Show(pointSources[i], renderView, 'GeometryRepresentation')
            pointSourceDisplay[i].RenderPointsAsSpheres = 1
            pointSourceDisplay[i].PointSize = 10.0
            pointSourceDisplay[i].AmbientColor = [1.0, 0.0, 0.0]
            pointSourceDisplay[i].DiffuseColor = [1.0, 0.0, 0.0]
            if i < noPoints-2:
                lineSource2Display[i] = Show(lineSource2[i], renderView, 'GeometryRepresentation')
                lineSource2Display[i].RenderLinesAsTubes = 1
                lineSource2Display[i].LineWidth = 2.0 
                lineSource2Display[i].AmbientColor = [1.0, 0.0, 0.0]
                lineSource2Display[i].DiffuseColor = [1.0, 0.0, 0.0]

        splineSourceDisplay = [''] * noBridges
        for i in range(noBridges):
            splineSourceDisplay[i] = Show(splineSource[i], renderView, 'GeometryRepresentation')
            splineSourceDisplay[i].RenderLinesAsTubes = 1
            splineSourceDisplay[i].LineWidth = 7.0 
            splineSourceDisplay[i].AmbientColor = [1.0, 0.0, 0.0]
            splineSourceDisplay[i].DiffuseColor = [1.0, 0.0, 0.0]

        lineSourceDisplay = [''] * noMasts
        for i in range(noMasts):
            lineSourceDisplay[i] = Show(lineSource[i], renderView, 'GeometryRepresentation')
            lineSourceDisplay[i].RenderLinesAsTubes = 1
            lineSourceDisplay[i].LineWidth = 5.0 
            lineSourceDisplay[i].AmbientColor = [1.0, 0.0, 0.0]
            lineSourceDisplay[i].DiffuseColor = [1.0, 0.0, 0.0]



####################################################################################
## Layout 1 - Surface LIC plots
# create a new 'Clip'
if plotLIC:
    topologyDisplay = Show(topology, renderView1, 'StructuredGridRepresentation')
    topologyDisplay.Representation = 'Surface'
    topologyDisplay.Texture = CreateTexture(textureFileName_NIB)
    showRunway(renderView1)
    showMasts(renderView1)
    # current camera placement for renderView1
    #renderView1.InteractionMode = '2D'
    CreateLayout('Layout #2')
    layout2 = GetLayoutByName("Layout #2")
    
    # Create a new 'Render View'
    renderView2 = CreateView('RenderView')
    renderView2.AxesGrid = 'GridAxes3DActor'
    renderView2.StereoType = 'Crystal Eyes'
    renderView2.CameraFocalDisk = 1.0
    AssignViewToLayout(view=renderView2, layout=layout2, hint=0)
    
    # create a new 'Slice'
    slice2 = Slice(registrationName='Slice2', Input=calculator1)
    slice2.SliceType = 'Plane'
    slice2.HyperTreeGridSlicer = 'Plane'
    slice2.SliceOffsetValues = [0.0]
    slice2.SliceType.Normal = [0.0, 0.0, 1.0]
    slice2.SliceType.Origin = np.append(runwayCenter[:2],horizontalSlice_z)
    Hide3DWidgets(proxy=slice2.SliceType)
    
    # show data in view
    slice1Display = Show(slice1, renderView1, 'UnstructuredGridRepresentation')
    ColorBy(slice1Display, ('POINTS', 'sqrtTKE'))
    slice1Display.SetRepresentationType('Surface LIC')
    slice1Display.ColorMode = 'Multiply'
    slice1Display.EnhanceContrast = 'Color Only'
    slice1Display.HighColorContrastEnhancementFactor = 0.3
    slice1Display.LICIntensity = 0.8
    slice1Display.MapModeBias = 0.2
    slice1Display.RescaleTransferFunctionToDataRange(True, False)
    slice1Display.SetScalarBarVisibility(renderView2, True)# set active source

    sqrtTKELUTColorBar = GetScalarBar(sqrtTKELUT, renderView1)
    sqrtTKELUTColorBar.WindowLocation = 'UpperRightCorner'
    sqrtTKELUTColorBar.Title = 'Turbulence, $\sqrt{k}$'
    sqrtTKELUTColorBar.ComponentTitle = ''
    sqrtTKELUTColorBar.ScalarBarLength = scalarBarLength

    # show data in view
    slice2Display = Show(slice2, renderView2, 'UnstructuredGridRepresentation')
    ColorBy(slice2Display, ('POINTS', 'sqrtTKE'))
    slice2Display.SetRepresentationType('Surface LIC')
    slice2Display.ColorMode = 'Multiply'
    slice2Display.EnhanceContrast = 'Color Only'
    slice2Display.HighColorContrastEnhancementFactor = 0.3
    slice2Display.LICIntensity = 0.8
    slice2Display.MapModeBias = 0.2
    slice2Display.RescaleTransferFunctionToDataRange(True, False)
    slice2Display.SetScalarBarVisibility(renderView2, True)# set active source
    
    sqrtTKELUTColorBar_2 = GetScalarBar(sqrtTKELUT, renderView2)
    sqrtTKELUTColorBar_2.WindowLocation = 'UpperRightCorner'
    sqrtTKELUTColorBar_2.Title = 'Turbulence, $\sqrt{k}$'
    sqrtTKELUTColorBar_2.ComponentTitle = ''
    sqrtTKELUTColorBar_2.ScalarBarLength = scalarBarLength

    topologyDisplay = Show(topology, renderView2, 'UnstructuredGridRepresentation')
    topologyDisplay.Representation = 'Surface'
    topologyDisplay.Texture = CreateTexture(textureFileName_NIB)

    # current camera placement for renderView2
    renderView2.OrientationAxesVisibility = 0
    showRunway(renderView2)
    showMasts(renderView2)
 
####################################################################################
## Layout 3 - Stream lines
if plotStreamLines:
    CreateLayout('Layout #3')
    layout3 = GetLayoutByName("Layout #3")
    
    # Create a new 'Render View'
    renderView3 = CreateView('RenderView')
    renderView3.AxesGrid = 'GridAxes3DActor'
    renderView3.StereoType = 'Crystal Eyes'
    renderView3.CameraFocalDisk = 1.0
    renderView3.OrientationAxesVisibility = 0
    AssignViewToLayout(view=renderView3, layout=layout3, hint=0)
    
    topologyDisplay = Show(topology, renderView3, 'UnstructuredGridRepresentation')
    topologyDisplay.Representation = 'Surface'
    topologyDisplay.Texture = CreateTexture(textureFileName_NIB)
    
    # create a new 'Stream Tracer'
    streamTracer1 = StreamTracer(registrationName='StreamTracer1', Input=calculator1, SeedType='Line')
    streamTracer1.MaximumStreamlineLength = 100000.0
    streamTracer1.SeedType.Point1 = streamTracerPoint1
    streamTracer1.SeedType.Point2 = streamTracerPoint2
    streamTracer1.SeedType.Resolution = 200
    streamTracer1.MaximumStepLength = 0.5
    streamTracer1.Vectors = ['POINTS', 'u']
    renderView3.Update()

    # show data in view
    streamTracer1Display = Show(streamTracer1, renderView3, 'GeometryRepresentation')
    streamTracer1Display.Representation = 'Surface'
    streamTracer1Display.ColorArrayName = ['POINTS', 'sqrtTKE']
    streamTracer1Display.RescaleTransferFunctionToDataRange(False, True)
    streamTracer1Display.SetScalarBarVisibility(renderView3, True)
    streamTracer1Display.SetRepresentationType('Wireframe')
    streamTracer1Display.RenderLinesAsTubes = 1
    streamTracer1Display.SetRepresentationType('Surface')
    streamTracer1Display.LineWidth = 3.0
    Hide3DWidgets(proxy=streamTracer1.SeedType)

    sqrtTKELUTColorBar_3 = GetScalarBar(sqrtTKELUT, renderView3)
    sqrtTKELUTColorBar_3.Title = 'Turbulence, $\sqrt{k}$'
    sqrtTKELUTColorBar_3.Orientation = 'Vertical'
    sqrtTKELUTColorBar_3.WindowLocation = 'UpperRightCorner'
    sqrtTKELUTColorBar_3.ComponentTitle = ''
    sqrtTKELUTColorBar_3.ScalarBarLength = scalarBarLength
    
    showRunway(renderView3)
    showMasts(renderView3)
    


####################################################################################
## Layout 4 - Volume rendering
if plotVolumeRendering:
    CreateLayout('Layout #4')
    layout4 = GetLayoutByName("Layout #4")
    
    # Create a new 'Render View'
    renderView4 = CreateView('RenderView')
    renderView4.AxesGrid = 'GridAxes3DActor'
    renderView4.StereoType = 'Crystal Eyes'
    renderView4.CameraFocalDisk = 1.0
    renderView4.OrientationAxesVisibility = 0
    AssignViewToLayout(view=renderView4, layout=layout4, hint=0)
    
    topologyDisplay = Show(topology, renderView4, 'UnstructuredGridRepresentation')
    topologyDisplay.Representation = 'Surface'
    topologyDisplay.Texture = CreateTexture(textureFileName_NIB)
    
    # create a new 'Stream Tracer'
    calculator1Display = Show(calculator1, renderView4, 'UnstructuredGridRepresentation')
    calculator1Display.Representation = 'Volume'
    calculator1Display.ColorArrayName = ['POINTS', 'sqrtTKE']

    showRunway(renderView4)
    showMasts(renderView4)
    
    sqrtTKELUTColorBar_4 = GetScalarBar(sqrtTKELUT, renderView4)
    sqrtTKELUTColorBar_4.Title = 'Turbulence, $\sqrt{k}$'
    sqrtTKELUTColorBar_4.Orientation = 'Vertical'
    sqrtTKELUTColorBar_4.WindowLocation = 'UpperRightCorner'
    sqrtTKELUTColorBar_4.ComponentTitle = ''
    sqrtTKELUTColorBar_4.ScalarBarLength = scalarBarLength
    if timeStepAnnotation:
        annotateTimeStep(calculator1,renderView4,location='UpperLeftCorner')


if plotIPPCmapsHorizontal:
    CreateLayout('Layout #5')
    layout5 = GetLayoutByName("Layout #5")
    
    # Create a new 'Render View'
    renderView5 = CreateView('RenderView')
    renderView5.AxesGrid = 'GridAxes3DActor'
    renderView5.StereoType = 'Crystal Eyes'
    renderView5.CameraFocalDisk = 1.0
    renderView5.OrientationAxesVisibility = 0
    AssignViewToLayout(view=renderView5, layout=layout5, hint=0)

    showRunway(renderView5)
    
    if plotTakeOffLines:
        takeOffLineEastDisplay = Show(takeOffLineEast, renderView5, 'GeometryRepresentation')
        takeOffLineEastDisplay.DiffuseColor = takeOffLineColor 
        takeOffLineEastDisplay.LineWidth = 2.0
        takeOffLineWestDisplay = Show(takeOffLineWest, renderView5, 'GeometryRepresentation')
        takeOffLineWestDisplay.DiffuseColor = takeOffLineColor 
        takeOffLineWestDisplay.LineWidth = 2.0

    calculator2Display = Show(calculator2, renderView5, 'StructuredGridRepresentation')
    calculator2Display.Representation = 'Surface'
    calculator2Display.Texture = CreateTexture(textureFileName_topo)
    HideScalarBarIfNotNeeded(heightLUT, renderView5)
    
    contour1 = Contour(Input=calculator2)
    contour1.ContourBy = ['POINTS', 'height']
    contour1.Isosurfaces = [0.0001]
    contour1.PointMergeMethod = 'Uniform Binning'
    contour1Display = Show(contour1, renderView5, 'GeometryRepresentation')
    ColorBy(contour1Display, None)
    contour1Display.DiffuseColor = [0.552941176470588, 0.149019607843137, 0.129411764705882]

    # Create an extended cone at which to plot results
    vtkName = outputPath+'temp'
    coffeeFilter(runwayEndWest,runwayEndEast,h_max,phi1=angleOfAttack_West,phi2=angleOfAttack_East,n=400,name=vtkName)
    extendedCone = LegacyVTKReader(registrationName='ExtendedCone', FileNames=[vtkName+'.vtk'])
    extendedConeDisplay = Show(extendedCone, renderView5, 'StructuredGridRepresentation')
    contour2 = Contour(Input=extendedCone)
    Hide(extendedCone,renderView5)
    contour2.ContourBy = ['POINTS', 'z']
    contour2.Isosurfaces = np.array([500,1000,1500,2000,2500,3000,3500,4000,4500])*0.3048 #convert to meter
    contour2.PointMergeMethod = 'Uniform Binning'
    contour2Display = Show(contour2, renderView5, 'GeometryRepresentation')
    ColorBy(contour2Display, None)
    contour2Display.DiffuseColor = [0.286274509803922, 0.317647058823529, 0.968627450980392]
    
    resampleWithDataset1 = ResampleWithDataset(SourceDataArrays=calculatorIPPC, DestinationMesh=extendedCone)
    resampleWithDataset1.CellLocator = 'Static Cell Locator'
    resampleWithDataset1Display = Show(resampleWithDataset1, renderView5, 'StructuredGridRepresentation')
    resampleWithDataset1Display.ColorArrayName = ['POINTS', 'sqrtTKE_IPPC']
    resampleWithDataset1Display.Ambient = 1.0
    resampleWithDataset1Display.Diffuse = 0.0

    sqrtTKELUTColorBar_5 = GetScalarBar(sqrtTKE_IPPCLUT, renderView5)
    sqrtTKELUTColorBar_5.Title = 'Turbulence, $\sqrt{k}$'
    sqrtTKELUTColorBar_5.Orientation = 'Vertical'
    sqrtTKELUTColorBar_5.WindowLocation = 'UpperRightCorner'
    sqrtTKELUTColorBar_5.ComponentTitle = ''
    sqrtTKELUTColorBar_5.ScalarBarLength = 0.13 
    sqrtTKELUTColorBar_5.AutomaticLabelFormat = 1
    sqrtTKELUTColorBar_5.UseCustomLabels = 1
    sqrtTKELUTColorBar_5.CustomLabels = [2, 3, 4, 5]
    sqrtTKELUTColorBar_5.AddRangeLabels = 0
    
    contour3 = Contour(Input=resampleWithDataset1)
    contour3.ContourBy = ['POINTS', 'sqrtTKE_IPPC']
    contour3.Isosurfaces = [2.0,3.0,4.0]
    contour3.PointMergeMethod = 'Uniform Binning'
    contour3Display = Show(contour3, renderView5, 'GeometryRepresentation')
    try:
        ColorBy(contour3Display, None)
    except:
        print('No turbulence above 2 to be shown')
    
    contour3Display.DiffuseColor = [0.0, 0.0, 0.0]

    calculatorConeProj = Calculator(Input=resampleWithDataset1)
    calculatorConeProj.ResultArrayName = 'projCone'
    z_proj = 5000
    calculatorConeProj.Function = 'coordsX*iHat+coordsY*jHat+'+str(z_proj)+'*kHat'
    calculatorConeProj.CoordinateResults = 1
    plane1 = Plane()
    plane1.Origin = [P_sw[0],P_sw[1],z_proj]
    plane1.Point1 = [P_se[0],P_se[1],z_proj]
    plane1.Point2 = [P_nw[0],P_nw[1],z_proj]
    plane1.XResolution = np.round(np.linalg.norm(P_sw[:2]-P_se[:2])/1000).astype('int')
    plane1.YResolution = np.round(np.linalg.norm(P_sw[:2]-P_nw[:2])/1000).astype('int')
    resampleWithDataset2 = ResampleWithDataset(SourceDataArrays=calculatorConeProj, DestinationMesh=plane1)
    resampleWithDataset2.CellLocator = 'Static Cell Locator'
    glyph1 = Glyph(Input=resampleWithDataset2, GlyphType='Arrow')
    glyph1.OrientationArray = ['POINTS', 'u']
    glyph1.ScaleArray = ['POINTS', 'u']
    glyph1.ScaleFactor = glyphScale1
    glyph1.GlyphTransform = 'Transform2'
    glyph1.GlyphMode = 'Every Nth Point'
    glyph1.Stride = 1
    glyph1Display = Show(glyph1, renderView5, 'GeometryRepresentation')
    glyph1Display.DiffuseColor = [0.0, 0.0, 0.0]
    ColorBy(glyph1Display, None)

    renderView5.InteractionMode = '2D'
    renderView5.CameraPosition = [54279.27032496591, 6950488.075249355, 25242.965782087373]
    renderView5.CameraFocalPoint = [51365.64987931158, 6949878.21067271, -393.12364711190645]
    renderView5.CameraViewUp = [-0.0026685028441370807, 0.999720760110162, -0.023479371740545256]
    renderView5.CameraParallelScale = 11840.135737817194
    
    
if plotIPPCmapsVertical:
    CreateLayout('Layout #6')
    layout6 = GetLayoutByName("Layout #6")
    # Create a new 'Render View'
    renderView6 = CreateView('RenderView')
    renderView6.AxesGrid = 'GridAxes3DActor'
    renderView6.StereoType = 'Crystal Eyes'
    renderView6.CameraFocalDisk = 1.0
    renderView6.OrientationAxesVisibility = 0
    AssignViewToLayout(view=renderView6, layout=layout6, hint=0)

    def SItoFTvsNM(obj,regName='Intermediate calculator'):
        rotObj = rotTran(obj,-runwayRotVector,-runwayCenter,rotateFirst=False)
        if isinstance(rotObj,(list,np.ndarray)): 
            return [rotObj[0]/1852, rotObj[1]/1852, 0.3048*rotObj[2]/100]
        else:
            result = Calculator(registrationName=regName, Input=rotObj)
            result.Function = '(coordsX*iHat+coordsY*jHat)/1852 + 0.3048*coordsZ*kHat/100'
            result.ResultArrayName = 'height'
            result.CoordinateResults = 1
        return result

    sliceTopologyTransform = SItoFTvsNM(sliceTopology,regName='Slice topology transformed')
    sliceTopologyTransformDisplay = Show(sliceTopologyTransform, renderView6, 'StructuredGridRepresentation')
    sliceTopologyTransformDisplay.LineWidth = 3.0
    if False:
        ColorBy(sliceTopologyTransformDisplay, 'height')
    else:
        sliceTopologyTransformDisplay.DiffuseColor = [0.0, 0.0, 0.0]
        ColorBy(sliceTopologyTransformDisplay, None)

    resampledSliceTransform = SItoFTvsNM(resampledSlice,regName='Slice1 transformed')
    slice1Transform = SItoFTvsNM(slice1,regName='Slice1 transformed')
    slice1TransformDisplay = Show(slice1Transform, renderView6, 'StructuredGridRepresentation')
    slice1TransformDisplay.ColorArrayName = ['POINTS', 'sqrtTKE_IPPC']
    slice1TransformDisplay.Ambient = 1.0
    slice1TransformDisplay.Diffuse = 0.0

    glyph2 = Glyph(registrationName='Glyphs',Input=resampledSliceTransform, GlyphType='Arrow')
    glyph2.OrientationArray = ['POINTS', 'u_proj']
    glyph2.ScaleArray = ['POINTS', 'u_proj']
    glyph2.ScaleFactor = glyphScale2
    glyph2.GlyphTransform = 'Transform2'
    glyph2.GlyphMode = 'Every Nth Point'
    glyph2.Stride = 1
    glyph2Display = Show(glyph2, renderView6, 'GeometryRepresentation')
    glyph2Display.DiffuseColor = [0.0, 0.6, 0.0]
    ColorBy(glyph2Display, None)
    
    sqrtTKELUTColorBar_6 = GetScalarBar(sqrtTKE_IPPCLUT, renderView6)
    sqrtTKELUTColorBar_6.Title = 'Turbulence, $\sqrt{k}$'
    sqrtTKELUTColorBar_6.Orientation = 'Vertical'
    sqrtTKELUTColorBar_6.WindowLocation = 'UpperRightCorner'
    sqrtTKELUTColorBar_6.ComponentTitle = ''
    sqrtTKELUTColorBar_6.ScalarBarLength = 0.13 
    sqrtTKELUTColorBar_6.AutomaticLabelFormat = 1
    sqrtTKELUTColorBar_6.UseCustomLabels = 1
    sqrtTKELUTColorBar_6.CustomLabels = [2, 3, 4, 5]
    sqrtTKELUTColorBar_6.AddRangeLabels = 0
    
    contour4 = Contour(registrationName='Contour IPPC vertical', Input=slice1Transform)
    contour4.ContourBy = ['POINTS', 'sqrtTKE_IPPC']
    contour4.Isosurfaces = [2.0,3.0,4.0]
    contour4.PointMergeMethod = 'Uniform Binning'
    contour4Display = Show(contour4, renderView6, 'GeometryRepresentation')
    try:
        ColorBy(contour4Display, None)
    except:
        print('No turbulence above 2 to be shown')
    
    contour4Display.DiffuseColor = [0.0, 0.0, 0.0]
    
    runwayTransform = SItoFTvsNM(runway,'Runway')
    showRunway(renderView6,runwayTransform)

    if plotTakeOffLines:
        takeOffLineEastTransform = SItoFTvsNM(takeOffLineEast,'East take off')
        takeOffLineEastDisplay = Show(takeOffLineEastTransform, renderView6, 'GeometryRepresentation')
        takeOffLineEastDisplay.DiffuseColor = takeOffLineColor 
        takeOffLineEastDisplay.LineWidth = 2.0
        takeOffLineWestTransform = SItoFTvsNM(takeOffLineWest, 'West take off')
        takeOffLineWestDisplay = Show(takeOffLineWestTransform, renderView6, 'GeometryRepresentation')
        takeOffLineWestDisplay.DiffuseColor = takeOffLineColor 
        takeOffLineWestDisplay.LineWidth = 2.0
    
    # Add 2D axes
    planeSliceTransform = SItoFTvsNM(planeSlice,'dummy')
    axes1 = Axes()
    Origin = SItoFTvsNM(np.array(planeSlice.Origin))
    Point1 = SItoFTvsNM(np.array(planeSlice.Point1))
    Point2 = SItoFTvsNM(np.array(planeSlice.Point2))
    axes1.Origin = Origin
    axes1.ScaleFactor = 10.0
    axes1Display = Show(axes1, renderView6, 'GeometryRepresentation')
    axes1Display.DataAxesGrid.GridAxesVisibility = 1
    axes1Display.DataAxesGrid.XTitle = 'Distance in NM'
    axes1Display.DataAxesGrid.YTitle = 'Distance in NM'
    axes1Display.DataAxesGrid.ZTitle = 'Height in 1000ft'
    axes1Display.DataAxesGrid.ShowGrid = 0
    axes1Display.DataAxesGrid.UseCustomBounds = 1
    axes1Display.DataAxesGrid.CustomBounds = [Origin[0], Point1[0], Origin[1], Origin[1], Origin[2], Point2[2]]
    axes1Display.DataAxesGrid.XLabelFontSize = 14
    axes1Display.DataAxesGrid.YLabelFontSize = 14
    axes1Display.DataAxesGrid.ZLabelFontSize = 14
    axes1Display.DataAxesGrid.XTitleFontSize = 16
    axes1Display.DataAxesGrid.YTitleFontSize = 16
    axes1Display.DataAxesGrid.ZTitleFontSize = 16
    axes1Display.DataAxesGrid.FacesToRender = 63
    axes1Display.DataAxesGrid.AxesToLabel = 5 # min-X and min-Z

    renderView6.InteractionMode = '2D'
    renderView6.CameraPosition = [3.081842464299491, -42.21395628114839, 5.605991802855324]
    renderView6.CameraFocalPoint = [3.081842464299491, -2.4408295139118075e-06, 5.605991802855324]
    renderView6.CameraViewUp = [0.0, 0.0, 1.0]
    renderView6.CameraViewAngle = 29.0
    renderView6.CameraParallelScale = 7.450788788748424 
    
if plotOverTime:
    CreateLayout('Layout #7')
    layout7 = GetLayoutByName("Layout #7")
    quartileChartView1 = CreateView('QuartileChartView')
    quartileChartView1.BottomAxisTitle = 'Iteration number'
    quartileChartView1.LeftAxisTitle = 'TKE $[m^2/s^2]$'
    AssignViewToLayout(view=quartileChartView1, layout=layout7, hint=0)

    cpcsv = [''] * noPoints
    cpcsvDisplay = [''] * noPoints
    colors = cm.jet(range(256))[0::256//noPoints]
    for i in range(0,noPoints):

        cpcsv[i] = CSVReader(registrationName=labels[i], FileName=[outputPath+fileName+'_Point'+str(i+1)+'.csv'])
        cpcsvDisplay[i] = Show(cpcsv[i], quartileChartView1, 'XYChartRepresentation')
        cpcsvDisplay[i].UseIndexForXAxis = 1
        #cpcsvDisplay.XArrayName = ''
        #cpcsvDisplay.SeriesLineStyle = ['TKE', '0']
        #cpcsvDisplay[i].SeriesMarkerSize = ['TKE', '8']
        #cpcsvDisplay[i].SeriesMarkerStyle = ['TKE', '4']
        cpcsvDisplay[i].SeriesColor = ['TKE', str(colors[i][0]), str(colors[i][1]), str(colors[i][2])]
        cpcsvDisplay[i].SeriesLabel = ['TKE', labels[i]]
        cpcsvDisplay[i].SeriesVisibility = ['TKE']
    
    saveScreenShot(quartileChartView1,outputPath+fileName+'TKE',saveScreenShots)
    

if plotIPPCmapsHorizontal:
    heightLUT.AutomaticRescaleRangeMode = "Never"
    heightLUT.RescaleOnVisibilityChange = 0
    heightLUT.EnableOpacityMapping = 0
    heightLUT.RescaleTransferFunction(0.0, height_max)
    heightLUT.ApplyPreset('google', True)

# get color transfer function/color map for 'sqrtTKE'
sqrtTKELUT.AutomaticRescaleRangeMode = "Never"
sqrtTKELUT.RescaleOnVisibilityChange = 0
sqrtTKELUT.EnableOpacityMapping = 0
sqrtTKELUT.RescaleTransferFunction(0.0, sqrtTKE_max)
sqrtTKELUT.ApplyPreset('SINTEF1', True)

if plotIPPCmapsHorizontal or plotIPPCmapsVertical:
    # get color transfer function/color map for 'sqrtTKE'
    sqrtTKE_IPPCLUT.AutomaticRescaleRangeMode = "Never"
    sqrtTKE_IPPCLUT.RescaleOnVisibilityChange = 0
    sqrtTKE_IPPCLUT.EnableOpacityMapping = 1
    sqrtTKE_IPPCLUT.ApplyPreset('IPPC', True)
    sqrtTKE_IPPCLUT.RescaleTransferFunction(2.0, sqrtTKE_max)
    sqrtTKE_IPPCPWF = GetOpacityTransferFunction('sqrtTKE_IPPC')
    sqrtTKE_IPPCPWF.ApplyPreset('IPPC', True)
    sqrtTKE_IPPCPWF.RescaleTransferFunction(2.0, sqrtTKE_max)

# get opacity transfer function/opacity map for 'sqrtTKE'
sqrtTKEPWF.AllowDuplicateScalars = 1
sqrtTKEPWF.ScalarRangeInitialized = 1
sqrtTKEPWF.RescaleTransferFunction(0.0, sqrtTKE_max)
sqrtTKEPWF.Points = [0.0, 0.0, 0.5, 0.0, 1.6, 0.0, 0.5, 0.0, 4.0, 1.0, 0.5, 0.0]

renderView1.InteractionMode = '3D'
renderView1.CameraPosition = [42315.178435332025, 6913455.09641479, 39925.03392197488]
renderView1.CameraFocalPoint = [47321.88134561954, 6944851.971259325, -590.056714551883]
renderView1.CameraViewUp = [0.10754399354149799, 0.7793741712263841, 0.6172602293023044]
renderView1.CameraParallelScale = 34572.810251475086
renderView1.CameraViewAngle = 30.0

color = 'white'
RenderAllViews()

if plotLIC:
    insertSINTEFlogo(renderView1,color)
    saveScreenShot(renderView1,outputPath+fileName+'surfaceLICside_bridge'+str(bridge),saveScreenShots)
    
    copyCamera(renderView1,renderView2)
    insertSINTEFlogo(renderView2,color)
    saveScreenShot(renderView2,outputPath+fileName+'surfaceLICtop_bridge'+str(bridge),saveScreenShots)

if plotStreamLines:
    insertSINTEFlogo(renderView3,color)
    copyCamera(renderView1,renderView3)
    saveScreenShot(renderView3,outputPath+fileName+'streamTracer_bridge'+str(bridge),saveScreenShots)
        
if plotVolumeRendering:
    insertSINTEFlogo(renderView4,color)
    copyCamera(renderView1,renderView4)
    saveScreenShot(renderView4,outputPath+fileName+'volumeRendering_bridge'+str(bridge),saveScreenShots)
    saveAnimation(renderView4,outputPath+fileName+'volumeRendering_bridge'+str(bridge),noSteps,makeVideo)

if plotIPPCmapsHorizontal:
    insertSINTEFlogo(renderView5,'white')
    saveScreenShot(renderView5,outputPath+fileName+'IPPC_horizontal_bridge'+str(bridge),saveScreenShots)
    ColorBy(calculator2Display, None)
    Hide(contour1, renderView5)
    saveScreenShot(renderView5,outputPath+fileName+'IPPC_horizontal_topo4_bridge'+str(bridge),saveScreenShots)

if plotIPPCmapsVertical:
    insertSINTEFlogo(renderView6,'blue')
    saveScreenShot(renderView6,outputPath+fileName+'IPPC_vertical_bridge'+str(bridge),saveScreenShots)

if plotError:
    slice1Display.Representation = 'Surface'
    ColorBy(slice1Display, ('CELLS', 'Continuous global L2-projection |u^*-u^h|_H1'))
    HideScalarBarIfNotNeeded(sqrtTKELUT, renderView1)
    continuousglobalL2projectionuuh_H1PWF = GetOpacityTransferFunction('ContinuousglobalL2projectionuuh_H1')
    continuousglobalL2projectionuuh_H1LUT = GetColorTransferFunction('ContinuousglobalL2projectionuuh_H1')
    continuousglobalL2projectionuuh_H1LUT.AutomaticRescaleRangeMode = "Never"
    continuousglobalL2projectionuuh_H1LUT.RescaleOnVisibilityChange = 0
    continuousglobalL2projectionuuh_H1LUT.EnableOpacityMapping = 0
    slice1Display.RescaleTransferFunctionToDataRange(False, True)
    continuousglobalL2projectionuuh_H1LUTColorBar = GetScalarBar(continuousglobalL2projectionuuh_H1LUT, renderView1)
    continuousglobalL2projectionuuh_H1LUTColorBar.WindowLocation = 'AnyLocation'
    continuousglobalL2projectionuuh_H1LUTColorBar.ScalarBarLength = 0.33000000000000007
    continuousglobalL2projectionuuh_H1LUTColorBar.Position = [0.9340269406943105, 0.09444444444444428]
    continuousglobalL2projectionuuh_H1LUTColorBar.ScalarBarLength = 0.3300000000000001
    continuousglobalL2projectionuuh_H1LUTColorBar.TitleColor = [1.0, 1.0, 1.0]
    continuousglobalL2projectionuuh_H1LUTColorBar.LabelColor = [1.0, 1.0, 1.0]
    continuousglobalL2projectionuuh_H1LUTColorBar.Title = 'Continuous global $L^2$-projection $|u^*-u^h|_{H^1}$'
    continuousglobalL2projectionuuh_H1LUT.RescaleTransferFunction(0.0, 100.0)
    continuousglobalL2projectionuuh_H1PWF.RescaleTransferFunction(0.0, 100.0)
    continuousglobalL2projectionuuh_H1LUT.ApplyPreset('SINTEF1', True)
    renderView1.CameraPosition = [39953.1497599849, 6936202.267042985, 2614.4419848900225]
    renderView1.CameraFocalPoint = [39143.54568842529, 6947093.529147031, 516.8749792426419]
    renderView1.CameraViewUp = [-0.012278300117821047, 0.18822195283001805, 0.9820497644310452]
    renderView1.CameraParallelScale = 28350.540722069076
    saveScreenShot(renderView1,outputPath+fileName+'surfaceLICside_bridge'+str(bridge)+'_Error',saveScreenShots)
