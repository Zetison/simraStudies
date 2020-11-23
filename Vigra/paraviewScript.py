from os.path import expanduser
import sys
import numpy as np
from scipy.spatial.transform import Rotation as Rot
import copy
#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()
 
home = expanduser("~")
sys.path.insert(1, home+'/kode/paraUtils')
from SINTEFlogo import insertSINTEFlogo
from sources import coffeeFilter
from utils import latLonUTM
outputPath = home+'/results/simra/Vigra_coarse/'
#topoRes = '50m'
topoRes = '10m'
topologyFileName = outputPath+'VigraFree'+topoRes+'.vts'
textureFileName_topo = outputPath+'VigraFree10m_topo.png'
textureFileName_NIB = outputPath+'VigraFree10m.png'
windCase = 4
noSteps = 201
finalTime = 11.7006
caseName = 'VigraFree10m_'+str(windCase)
sqrtTKE_max = 5.0
h_max = 4000
height_max = 1000 #max height of height colormap
angleOfAttack_West = 3.0
angleOfAttack_East = 3.4
runwayWidth = 150.0
runwayEndWest = latLonUTM("623323.41","0060532.89",35.3)
runwayEndEast = latLonUTM("623351.26","0060740.31",49.2)
SAVE_HIST = 20

viewSize = [1920, 1080]
scalarBarLength = 0.26
streamLines_z = 100
horizontalSlice_z = 300
glyphScale1 = 30
if windCase == 1:
	glyphScale2 = 0.05
elif windCase == 2:
	glyphScale2 = 0.1
elif windCase == 3:
	glyphScale2 = 0.05
elif windCase == 4:
	glyphScale2 = 0.1
elif windCase == 5:
	glyphScale2 = 0.1


plotRunwayStuff          = 1 

plotLIC                  = 0
plotStreamLines          = 0
plotVolumeRendering      = 1 
plotIPPCmapsHorizontal   = 0 
plotIPPCmapsVertical     = 0 

makeVideo                = 1
saveScreenShots          = 1
useTransparentBackground = 0

timeStepAnnotation       = makeVideo
# Extract UTM-33 coordinates from rwyCoords.png file (available at https://ais.avinor.no/no/AIP/View/95/2020-11-05-AIRAC/html/index-no-NO.html). The coordinates in rwyCoords.png are given in Lat-lon Degrees Minutes seconds format. Here, "THR COORD" is the start of the RWY (runway), and "RWY end COORD" is its end

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

# get the time-keeper
timeKeeper1 = GetTimeKeeper()
fileName = outputPath+caseName+'.pvd'

# create a new 'XML Unstructured Grid Reader'
simraPVDresults = PVDReader(registrationName=caseName, FileName=fileName)
simraPVDresults.PointArrays = ['ps', 'pts', 'tk', 'u']

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
topology.TimeArray = 'None'
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

def annotateTimeStep(obj,renderview,location='UpperLeftCorner'):
	pythonAnnotation = PythonAnnotation(registrationName='Time annotation', Input=obj)
	pythonAnnotation.ArrayAssociation = 'Point Data'
	pythonAnnotation.Expression = '"Step: %d\\nTime: %0.2fs" % ((time_index+1)*'+str(SAVE_HIST)+', time_value)'
	pythonAnnotationDisplay = Show(pythonAnnotation, renderview, 'TextSourceRepresentation')
	pythonAnnotationDisplay.WindowLocation = location 
	pythonAnnotationDisplay.FontSize = 5

####################################################################################
## Layout 1 - Surface LIC plots
# create a new 'Clip'
if plotLIC:
	topologyDisplay = Show(topology, renderView1, 'StructuredGridRepresentation')
	topologyDisplay.Representation = 'Surface'
	topologyDisplay.Texture = CreateTexture(textureFileName_NIB)
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
	slice1Display.NumberOfSteps = 100
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
	slice2Display.NumberOfSteps = 100
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
	streamTracer1.SeedType.Resolution = 200
	streamTracer1.Vectors = ['POINTS', 'u']
	streamTracer1.MaximumStepLength = 0.5
	streamTracer1.MaximumStreamlineLength = 100000.0
	P_sw = [30613.073672987,6954997.078414851,streamLines_z]
	P_se = [70493.555053226,6956794.646711984,streamLines_z]
	P_nw = [29466.912132806,6982238.693045641,streamLines_z]
	P_ne = [69366.258652280,6983800.313949561,streamLines_z]
	if windCase == 1:
		streamTracer1.SeedType.Point1 = P_sw
		streamTracer1.SeedType.Point2 = P_nw
	elif windCase == 2 or windCase == 4 or windCase == 5:
		streamTracer1.SeedType.Point1 = P_sw
		streamTracer1.SeedType.Point2 = P_se
	else:
		streamTracer1.SeedType.Point1 = P_nw
		streamTracer1.SeedType.Point2 = P_ne
	
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

	sqrtTKELUTColorBar_4 = GetScalarBar(sqrtTKELUT, renderView4)
	sqrtTKELUTColorBar_4.Title = 'Turbulence, $\sqrt{k}$'
	sqrtTKELUTColorBar_4.Orientation = 'Vertical'
	sqrtTKELUTColorBar_4.WindowLocation = 'UpperRightCorner'
	sqrtTKELUTColorBar_4.ComponentTitle = ''
	sqrtTKELUTColorBar_4.ScalarBarLength = scalarBarLength
	if timeStepAnnotation:
		annotateTimeStep(calculator1,renderView4,location='UpperLeftCorner')


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
box.ZLength = 30.0
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

	if plotRunwayStuff:
		# Display runway and takeoff visualization
		runwayDisplay = Show(runway, renderView5, 'GeometryRepresentation')
		runwayDisplay.Representation = 'Surface'
		runwayDisplay.DiffuseColor = [0.0,0.0,0.0]
		takeOffLineEastDisplay = Show(takeOffLineEast, renderView5, 'GeometryRepresentation')
		takeOffLineEastDisplay.DiffuseColor = [0.0,0.0,0.0]
		takeOffLineEastDisplay.LineWidth = 2.0
		takeOffLineWestDisplay = Show(takeOffLineWest, renderView5, 'GeometryRepresentation')
		takeOffLineWestDisplay.DiffuseColor = [0.0,0.0,0.0]
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
	coffeeFilter(runwayEndWest,runwayEndEast,h_max/3,phi1=angleOfAttack_West,phi2=angleOfAttack_East,n=400,name=vtkName)
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
	planeLength = 34e3
	planeHeight = 25e3
	plane1.Origin = [31843,6955380,z_proj]
	plane1.Point1 = [plane1.Origin[0]+planeLength,plane1.Origin[1],z_proj]
	plane1.Point2 = [plane1.Origin[0],plane1.Origin[1]+planeHeight,z_proj]
	plane1.XResolution = np.round(planeLength/1000).astype('int')
	plane1.YResolution = np.round(planeHeight/1000).astype('int')
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
	renderView5.CameraPosition = [51708.919054119346,6968532.526518099,25105.82928050565]
	renderView5.CameraFocalPoint = [48795.29860846502,6967922.661941454,-530.2601486936218]
	renderView5.CameraViewUp = [0.0, 1.0, 0.0]
	renderView5.CameraParallelScale = 9785.23614695636
	
	
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
	
	if plotRunwayStuff:
		# Display runway and takeoff visualization
		runwayTransform = SItoFTvsNM(runway,'Runway')
		runwayDisplay = Show(runwayTransform, renderView6, 'GeometryRepresentation')
		runwayDisplay.Representation = 'Surface'
		runwayDisplay.DiffuseColor = [0.0,0.0,0.0]
		takeOffLineEastTransform = SItoFTvsNM(takeOffLineEast,'East take off')
		takeOffLineEastDisplay = Show(takeOffLineEastTransform, renderView6, 'GeometryRepresentation')
		takeOffLineEastDisplay.DiffuseColor = [0.0,0.0,0.0]
		takeOffLineEastDisplay.LineWidth = 2.0
		takeOffLineWestTransform = SItoFTvsNM(takeOffLineWest, 'West take off')
		takeOffLineWestDisplay = Show(takeOffLineWestTransform, renderView6, 'GeometryRepresentation')
		takeOffLineWestDisplay.DiffuseColor = [0.0,0.0,0.0]
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
	sqrtTKE_IPPCLUT.EnableOpacityMapping = 0
	sqrtTKE_IPPCLUT.RescaleTransferFunction(2.0, sqrtTKE_max)
	sqrtTKE_IPPCLUT.ApplyPreset('IPPC', True)
	sqrtTKE_IPPCLUT.EnableOpacityMapping = 1

# get opacity transfer function/opacity map for 'sqrtTKE'
sqrtTKEPWF.AllowDuplicateScalars = 1
sqrtTKEPWF.ScalarRangeInitialized = 1
sqrtTKEPWF.RescaleTransferFunction(0.0, sqrtTKE_max)
sqrtTKEPWF.Points = [0.0, 0.0, 0.5, 0.0, 1.6, 0.0, 0.5, 0.0, 4.0, 1.0, 0.5, 0.0]

renderView1.CameraPosition = [45915.71585360529, 6938099.25941025, 22536.227838368617]
renderView1.CameraFocalPoint = [51563.52982665786, 6966996.095785035, -614.0106624878072]
renderView1.CameraViewUp = [0.08910050644944713, 0.612082952061688, 0.7857579522638646]
renderView1.CameraParallelScale = 25143.74802298738
color = 'blue'
RenderAllViews()
def saveScreenShot(renderView,name):
	if saveScreenShots:
		SaveScreenshot(outputPath+caseName+name+'.png', renderView,
				FontScaling='Scale fonts proportionally',
				TransparentBackground=useTransparentBackground,
				ImageResolution=viewSize,
				ImageQuality=100)

def saveAnimation(renderView,name):
	if makeVideo:
		renderView.ViewSize = viewSize
		animationScene1 = GetAnimationScene()
		SaveAnimation(outputPath+caseName+name+'.ogv', renderView,
				FontScaling='Scale fonts proportionally',
				OverrideColorPalette='',
				StereoMode='No change',
				TransparentBackground=0,
				ImageQuality=100,
				FrameRate=15,
				ImageResolution=viewSize,
				FrameWindow=[0, noSteps-1])

def copyCamera(renderview1,renderview2):
	renderview2.CameraPosition = renderview1.CameraPosition
	renderview2.CameraFocalPoint = renderview1.CameraFocalPoint
	renderview2.CameraViewUp = renderview1.CameraViewUp
	renderview2.CameraParallelScale = renderview1.CameraParallelScale

if plotLIC:
	insertSINTEFlogo(renderView1,color)
	saveScreenShot(renderView1,'surfaceLICside')
	insertSINTEFlogo(renderView2,color)
	copyCamera(renderView1,renderView2)
	saveScreenShot(renderView2,'surfaceLICtop')

if plotStreamLines:
	insertSINTEFlogo(renderView3,color)
	copyCamera(renderView1,renderView3)
	saveScreenShot(renderView3,'streamTracer')
		
if plotVolumeRendering:
	insertSINTEFlogo(renderView4,color)
	copyCamera(renderView1,renderView4)
	saveScreenShot(renderView4,'volumeRendering')
	saveAnimation(renderView4,'volumeRendering')

if plotIPPCmapsHorizontal:
	insertSINTEFlogo(renderView5,'white')
	saveScreenShot(renderView5,'IPPC_horizontal')
	ColorBy(calculator2Display, None)
	Hide(contour1, renderView5)
	saveScreenShot(renderView5,'IPPC_horizontal_topo4')

if plotIPPCmapsVertical:
	insertSINTEFlogo(renderView6,'blue')
	saveScreenShot(renderView6,'IPPC_vertical')
