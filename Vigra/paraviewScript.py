from os.path import expanduser
import sys
import numpy as np
#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()
 
home = expanduser("~")
sys.path.insert(1, home+'/kode/paraUtils')
from SINTEFlogo import insertSINTEFlogo
from sources import extCone

outputPath = home+'/results/simra/Vigra/'
windCase = 2
caseName = 'VigraFree10m_'+str(windCase)
sqrtTKE_max = 5.0
h_max = 4000
angleOfAttack = 4.6
viewSize = [1920, 1080]
useTransparentBackground = True
scalarBarLength = 0.26
plotLIC = False
plotStreamLines = False
plotVolumeRendering = False
plotIPPCmaps = True
makeVideo = False
saveScreenShots = False
NumberOfFrames = 201
runwayEndWest = np.array([42890.95,6967738.00,10])
runwayEndEast = np.array([45123.04,6968439.01,10])
runwayCenter = (runwayEndEast+runwayEndWest)/2
# get animation scene
animationScene1 = GetAnimationScene()
LoadPalette(paletteName='WhiteBackground')
renderView1 = GetActiveViewOrCreate('RenderView')
renderView1.OrientationAxesVisibility = 0

# get the time-keeper
timeKeeper1 = GetTimeKeeper()
fileName = home+'/results/simra/Vigra/'+caseName+'.pvd'

# create a new 'XML Unstructured Grid Reader'
huntHill = PVDReader(registrationName=caseName, FileName=fileName)
huntHill.PointArrays = ['ps', 'pts', 'tk', 'u']

animationScene1.UpdateAnimationUsingDataTimeSteps()
animationScene1.GoToLast()
# update animation scene based on data timesteps

# get active view
layout1 = GetLayoutByName("Layout #1")

# create a new 'Calculator'
calculator1 = Calculator(registrationName='Calculator1', Input=huntHill)

# Properties modified on calculator1
calculator1.Function = 'sqrt(tk)'
calculator1.ResultArrayName = 'sqrtTKE'

sqrtTKELUT = GetColorTransferFunction('sqrtTKE')
sqrtTKERWF = GetOpacityTransferFunction('sqrtTKE')

# create a new 'Slice'
slice1 = Slice(registrationName='Slice1', Input=calculator1)
slice1.SliceType = 'Plane'
slice1.HyperTreeGridSlicer = 'Plane'
slice1.SliceOffsetValues = [0.0]
slice1.HyperTreeGridSlicer.Origin = [0.11449999999999994, 0.0, 0.8015]
slice1.SliceType.Normal = [0.30300368726820465, -0.9529862428683145, 0.0]
slice1.SliceType.Origin = [48671.06678998141, 6969600.636836766, 1223.1745535236325]
# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=slice1.SliceType)

# create a new 'XML Structured Grid Reader'
topology = XMLStructuredGridReader(registrationName='Topology',FileName=[outputPath+'VigraFree10m.vts'])
topology.TimeArray = 'None'
topology.PointArrayStatus = ['TCoords_']
layout1 = GetLayout()

####################################################################################
## Layout 1 - Surface LIC plots
# create a new 'Clip'
if plotLIC:
	topologyDisplay = Show(topology, renderView1, 'StructuredGridRepresentation')
	topologyDisplay.Representation = 'Surface'
	topologyDisplay.Texture = CreateTexture(outputPath+'VigraFree10m.png')
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
	slice2.HyperTreeGridSlicer.Origin = [0.11449999999999994, 0.0, 0.8015]
	slice2.SliceType.Normal = [0.0, 0.0, 1.0]
	slice2.SliceType.Origin = [49979.89939835307, 6969396.944194966, 300]
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
	topologyDisplay.Texture = CreateTexture(outputPath+'VigraFree10m.png')

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
	topologyDisplay.Texture = CreateTexture(outputPath+'VigraFree10m.png')
	
	# create a new 'Stream Tracer'
	#streamTracer1 = StreamTracer(registrationName='StreamTracer1', Input=calculator1, SeedType='Point Cloud')
	#streamTracer1.SeedType.NumberOfPoints = 2
	streamTracer1 = StreamTracer(registrationName='StreamTracer1', Input=calculator1, SeedType='Line')
	streamTracer1.SeedType.Resolution = 20
	streamTracer1.Vectors = ['POINTS', 'u']
	streamTracer1.MaximumStepLength = 0.5
	streamTracer1.MaximumStreamlineLength = 100000.0
	#streamTracer1.SeedType.Center = [54575.67373728964, 6967974.3741878215, 374.2404115051399]
	#streamTracer1.SeedType.Radius = 2000
	z = 100.0
	P_sw = [30613.073672987,6954997.078414851,z]
	P_se = [70493.555053226,6956794.646711984,z]
	P_nw = [29466.912132806,6982238.693045641,z]
	P_ne = [69366.258652280,6983800.313949561,z]
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
	

plotVolumeRendering
# get color transfer function/color map for 'sqrtTKE'
sqrtTKELUT.AutomaticRescaleRangeMode = "Never"
sqrtTKELUT.RescaleOnVisibilityChange = 0
sqrtTKELUT.EnableOpacityMapping = 0
sqrtTKELUT.ApplyPreset('SINTEF1', True)
sqrtTKELUT.RescaleTransferFunction(0.0, sqrtTKE_max)

# get opacity transfer function/opacity map for 'sqrtTKE'
sqrtTKERWF.AllowDuplicateScalars = 1
sqrtTKERWF.ScalarRangeInitialized = 1
sqrtTKERWF.RescaleTransferFunction(0.0, sqrtTKE_max)
sqrtTKERWF.Points = [0.0, 0.0, 0.5, 0.0, 1.6, 0.0, 0.5, 0.0, 4.0, 1.0, 0.5, 0.0]

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
	topologyDisplay.Texture = CreateTexture(outputPath+'VigraFree10m.png')
	
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

if plotIPPCmaps:
	CreateLayout('Layout #5')
	layout5 = GetLayoutByName("Layout #5")
	
	# Create a new 'Render View'
	renderView5 = CreateView('RenderView')
	renderView5.AxesGrid = 'GridAxes3DActor'
	renderView5.StereoType = 'Crystal Eyes'
	renderView5.CameraFocalDisk = 1.0
	renderView5.OrientationAxesVisibility = 0
	AssignViewToLayout(view=renderView5, layout=layout5, hint=0)
	
	#topologyDisplay = Show(topology, renderView5, 'UnstructuredGridRepresentation')
	#topologyDisplay.Representation = 'Surface'
	#topologyDisplay.Texture = CreateTexture(outputPath+'VigraFree10m.png')
	
	vtkName = outputPath+'temp'
	extCone(runwayEndWest,runwayEndEast,h_max,phi=angleOfAttack,n=200,name=vtkName)
	extendedCone = LegacyVTKReader(registrationName='ExtendedCone', FileNames=[vtkName+'.vtk'])
	resampleWithDataset1 = ResampleWithDataset(SourceDataArrays=calculator1, DestinationMesh=extendedCone)
	resampleWithDataset1Display = Show(resampleWithDataset1, renderView5, 'StructuredGridRepresentation')
	resampleWithDataset1Display.Opacity = 0.51
	renderView5.InteractionMode = '2D'
	renderView5.CameraPosition = [46920.61544565432, 6968698.369576645, 25646.089429199277]
	renderView5.CameraFocalPoint = runwayCenter
	renderView5.CameraViewUp = [0.0, 1.0, 0.0]
	renderView5.CameraParallelScale = 8086.972022277982

	sqrtTKELUTColorBar_5 = GetScalarBar(sqrtTKELUT, renderView5)
	sqrtTKELUTColorBar_5.Title = 'Turbulence, $\sqrt{k}$'
	sqrtTKELUTColorBar_5.Orientation = 'Vertical'
	sqrtTKELUTColorBar_5.WindowLocation = 'UpperRightCorner'
	sqrtTKELUTColorBar_5.ComponentTitle = ''
	sqrtTKELUTColorBar_5.ScalarBarLength = scalarBarLength
	


# get color transfer function/color map for 'sqrtTKE'
sqrtTKELUT.AutomaticRescaleRangeMode = "Never"
sqrtTKELUT.RescaleOnVisibilityChange = 0
sqrtTKELUT.EnableOpacityMapping = 0
sqrtTKELUT.ApplyPreset('SINTEF1', True)
sqrtTKELUT.RescaleTransferFunction(0.0, sqrtTKE_max)

# get opacity transfer function/opacity map for 'sqrtTKE'
sqrtTKERWF.AllowDuplicateScalars = 1
sqrtTKERWF.ScalarRangeInitialized = 1
sqrtTKERWF.RescaleTransferFunction(0.0, sqrtTKE_max)
sqrtTKERWF.Points = [0.0, 0.0, 0.5, 0.0, 1.6, 0.0, 0.5, 0.0, 4.0, 1.0, 0.5, 0.0]

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
		renderView2.ViewSize = viewSize
		animationScene1 = GetAnimationScene()
		slice2Display.SetRepresentationType('Surface')
		SaveAnimation(outputPath+caseName+name+'.ogv', renderView,
				FontScaling='Scale fonts proportionally',
				OverrideColorPalette='',
				StereoMode='No change',
				TransparentBackground=0,
				ImageQuality=100,
				FrameRate=15,
				ImageResolution=viewSize,
				FrameWindow=[0, NumberOfFrames-1])

def copyCamera(renderview1,renderview2):
	renderview2.CameraPosition = renderview1.CameraPosition
	renderview2.CameraFocalPoint = renderview1.CameraFocalPoint
	renderview2.CameraViewUp = renderview1.CameraViewUp
	renderview2.CameraParallelScale = renderview1.CameraParallelScale

if plotLIC:
	insertSINTEFlogo(renderView1,color)
	insertSINTEFlogo(renderView2,color)
	saveScreenShot(renderView1,'surfaceLICside')
	saveScreenShot(renderView2,'surfaceLICtop')
	copyCamera(renderView1,renderView2)
	saveAnimation(renderView2,'surface')

if plotStreamLines:
	insertSINTEFlogo(renderView3,color)
	copyCamera(renderView1,renderView3)
	saveScreenShot(renderView3,'streamTracer')
		
if plotVolumeRendering:
	insertSINTEFlogo(renderView4,color)
	copyCamera(renderView1,renderView4)
	saveScreenShot(renderView4,'volumeRendering')
	saveAnimation(renderView4,'volumeRendering')

