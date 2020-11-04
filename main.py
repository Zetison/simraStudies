from os.path import expanduser
import sys
import numpy as np
#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()
 
home = expanduser("~")
sys.path.insert(1, home+'/kode/paraUtils')
outputPath = home+'/results/simra/HuntHill/'
saveScreenShots = True
importData = False
u_inf = 0.4366 # is the freestream velocity of the fluid
p_inf = 100000-7.0636 # is the static pressure in the freestream (i.e. remote from any disturbance)
rho_inf = 1.225 # is the freestream fluid density (Air at sea level and 15 °C is 1.225kg/m^3)
z_lnPltLvl = 0.6 # z coordinate of 1D plots
plotLIC = True
plot1Dcurves = True
plotStreamLines = True
# get animation scene
animationScene1 = GetAnimationScene()

# get the time-keeper
timeKeeper1 = GetTimeKeeper()
noVTUfiles = 20
fileNames = [''] * noVTUfiles
for i in range(0,noVTUfiles):
	fileNames[i] = home+'/kode/simraStudies/HuntHill/HuntHill-'+str(i+1)+'.vtu'
fileName = home+'/results/simra/HuntHill/cont.pvd'

# create a new 'XML Unstructured Grid Reader'
huntHill = PVDReader(registrationName='HuntHill', FileName=fileName)
#huntHill = XMLUnstructuredGridReader(registrationName='HuntHill', FileName=fileNames)
#huntHill.CellArrayStatus = []
huntHill.PointArrays = ['ps', 'pts', 'u']
#huntHill.PointArrayStatus = ['ps', 'pT', 'pts', 'u', 'u']

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
layout1 = GetLayout()

# create a new 'Calculator'
calculator1 = Calculator(registrationName='Calculator1', Input=huntHill)

# Properties modified on calculator1
calculator1.ResultArrayName = 'Cp'
denominator = 0.5*rho_inf*u_inf**2
calculator1.Function = '(ps-'+str(p_inf)+')/'+str(denominator)

# get color transfer function/color map for 'Cp'
cpLUT = GetColorTransferFunction('Cp')
cpPWF = GetOpacityTransferFunction('Cp')

# create a new 'Calculator'
calculator2 = Calculator(registrationName='Calculator2', Input=calculator1)
calculator2.ResultArrayName = 'SpeedUp'
calculator2.Function = 'u_X/'+str(u_inf)

# create a new 'Slice'
slice1 = Slice(registrationName='Slice1', Input=calculator2)
slice1.SliceType = 'Plane'
slice1.HyperTreeGridSlicer = 'Plane'
slice1.SliceOffsetValues = [0.0]
slice1.HyperTreeGridSlicer.Origin = [0.11449999999999994, 0.0, 0.8015]
slice1.SliceType.Normal = [0.0, 1.0, 0.0]
slice1.SliceType.Origin = [0.0, 0.0, 0.0]

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=slice1.SliceType)
animationScene1.GoToLast()
LoadPalette(paletteName='WhiteBackground')
renderView1.OrientationAxesVisibility = 0

# create a new 'Clip'
clip1 = Clip(registrationName='Clip1', Input=slice1)
clip1.ClipType = 'Plane'
clip1.HyperTreeGridClipper = 'Plane'
clip1.Value = 0.5223809449277922
clip1.HyperTreeGridClipper.Origin = [0.11449999999999994, 0.0, 0.8015]
clip1.ClipType.Normal = [0.0, 0.0, 1.0]
clip1.ClipType.Origin = [0.0, 0.0, 0.8015]
	
# get color legend/bar for uLUT in view renderView1
uLUT = GetColorTransferFunction('u')
uLUT.ApplyPreset('Parula', True)
uLUT.InvertTransferFunction()
uPWF = GetOpacityTransferFunction('u')
uPWF.ApplyPreset('Parula', True)

# create a new 'XML Unstructured Grid Reader'
mapvtu = XMLUnstructuredGridReader(registrationName='map.vtu', FileName=[home+'/kode/simraStudies/map.vtu'])
mapvtu.TimeArray = 'None'

####################################################################################
## Layout 1 - Surface LIC plots
# create a new 'Clip'
if plotLIC:
	clip2 = Clip(registrationName='Clip2', Input=clip1)
	clip2.ClipType = 'Plane'
	clip2.HyperTreeGridClipper = 'Plane'
	clip2.Value = 0.5223809449277922
	
	# init the 'Plane' selected for 'ClipType'
	clip2.ClipType.Origin = [0.11449999999999994, 0.0, 0.40075000000000005]
	clip2.HyperTreeGridClipper.Origin = [0.11449999999999994, 0.0, 0.40075000000000005]
	clip2.ClipType.Origin = [-0.95, 0.0, 0.0]
	clip2.ClipType.Normal = [-1.0, 0.0, 0.0]
	
	clip2Display = Show(clip2, renderView1, 'UnstructuredGridRepresentation')
	ColorBy(clip2Display, ('POINTS', 'u', 'X'))
	clip2Display.Representation = 'Surface'
	clip2Display.RescaleTransferFunctionToDataRange(True, False)
	clip2Display.SetScalarBarVisibility(renderView1, True)
	clip2Display.SetRepresentationType('Surface LIC')
	clip2Display.NumberOfSteps = 100
	clip2Display.ColorMode = 'Multiply'
	clip2Display.EnhanceContrast = 'Color Only'
	clip2Display.MapModeBias = 0.25
	clip2Display.HighColorContrastEnhancementFactor = 0.3
	Hide3DWidgets(proxy=clip2.ClipType)

	# current camera placement for renderView1
	renderView1.InteractionMode = '2D'
	renderView1.CameraPosition = [0.4886131554367844, -1.5010035486303253, 0.39706868660766104]
	renderView1.CameraFocalPoint = [0.4886131554367844, 0.0, 0.39706868660766104]
	renderView1.CameraViewUp = [0.0, 0.0, 1.0]
	renderView1.CameraParallelScale = 0.38848830515199595
	
	# split cell
	layout1.SplitVertical(0, 0.5)
	
	# Create a new 'Render View'
	renderView2 = CreateView('RenderView')
	renderView2.AxesGrid = 'GridAxes3DActor'
	renderView2.StereoType = 'Crystal Eyes'
	renderView2.CameraFocalDisk = 1.0
	AssignViewToLayout(view=renderView2, layout=layout1, hint=2)
	
	# create a new 'Slice'
	slice2 = Slice(registrationName='Slice2', Input=calculator2)
	slice2.SliceType = 'Plane'
	slice2.HyperTreeGridSlicer = 'Plane'
	slice2.SliceOffsetValues = [0.0]
	slice2.HyperTreeGridSlicer.Origin = [0.11449999999999994, 0.0, 0.8015]
	slice2.SliceType.Normal = [0.0, 0.0, 1.0]
	slice2.SliceType.Origin = [0.0, 0.0, 0.2]
	Hide3DWidgets(proxy=slice2.SliceType)
	
	# create a new 'Clip'
	clip3 = Clip(registrationName='Clip3', Input=slice2)
	clip3.ClipType = 'Plane'
	clip3.HyperTreeGridClipper = 'Plane'
	clip3.Value = 0.6420483725261659
	clip3.HyperTreeGridClipper.Origin = [0.11449999999999994, 0.0, 0.2]
	clip3.ClipType.Origin = [0.0, 0.5, 0.0]
	clip3.ClipType.Normal = [0.0, 1.0, 0.0]
	
	clip4 = Clip(registrationName='Clip4', Input=clip3)
	clip4.ClipType = 'Plane'
	clip4.HyperTreeGridClipper = 'Plane'
	clip4.Value = 0.6420483725261659
	clip4.HyperTreeGridClipper.Origin = [0.11449999999999994, -0.3225, 0.2]
	clip4.ClipType.Origin = [0.0, -0.5, 0.0]
	clip4.ClipType.Normal = [0.0, -1.0, 0.0]
	
	# show data in view
	clip4Display = Show(clip4, renderView2, 'UnstructuredGridRepresentation')
	ColorBy(clip4Display, ('POINTS', 'u', 'X'))
	clip4Display.SetRepresentationType('Surface LIC')
	clip4Display.NumberOfSteps = 100
	clip4Display.ColorMode = 'Multiply'
	clip4Display.EnhanceContrast = 'Color Only'
	clip4Display.HighColorContrastEnhancementFactor = 0.3
	clip4Display.LICIntensity = 0.8
	clip4Display.MapModeBias = 0.2
	clip4Display.RescaleTransferFunctionToDataRange(True, False)
	clip4Display.SetScalarBarVisibility(renderView2, True)# set active source
	Hide3DWidgets(proxy=clip4.ClipType)
	
	uLUTColorBar = GetScalarBar(uLUT, renderView1)
	uLUTColorBar.ScalarBarLength = 0.6
	uLUTColorBar.WindowLocation = 'UpperRightCorner'
	uLUTColorBar.Title = 'x-component of u'
	uLUTColorBar.ComponentTitle = ''
	uLUTColorBar_1 = GetScalarBar(uLUT, renderView2)
	uLUTColorBar_1.ScalarBarLength = 0.6
	uLUTColorBar_1.WindowLocation = 'UpperRightCorner'
	uLUTColorBar_1.Title = 'x-component of u'
	uLUTColorBar_1.ComponentTitle = ''

	# current camera placement for renderView2
	renderView2.OrientationAxesVisibility = 0
	renderView2.InteractionMode = '2D'
	renderView2.CameraPosition = [0.32366589551106706, -0.021501017925559134, 1.9741519036963915]
	renderView2.CameraFocalPoint = [0.32366589551106706, -0.021501017925559134, 0.8015]
	renderView2.CameraParallelScale = 0.44436115213884014
	
	if saveScreenShots:
		SaveScreenshot(outputPath+'surfaceLICside.png', renderView1,
				FontScaling='Scale fonts proportionally',
				ImageQuality=100)
		SaveScreenshot(outputPath+'surfaceLICtop.png', renderView2,
				FontScaling='Scale fonts proportionally',
				ImageQuality=100)
	 
 
####################################################################################
## Layout 2 - 1D plots
if plot1Dcurves:
	CreateLayout('Layout #2')
	layout2 = GetLayoutByName("Layout #2")
	layout2.SplitVertical(0, 0.5)
	layout2.SplitHorizontal(1, 0.5)
	layout2.SplitHorizontal(2, 0.5)
	layout2.SplitHorizontal(5, 0.5)
	layout2.SplitHorizontal(6, 0.5)
	layout2.SplitHorizontal(11, 0.5)
	layout2.SplitHorizontal(12, 0.5)
	layout2.SplitHorizontal(13, 0.5)
	layout2.SplitHorizontal(14, 0.5)
	
	# Create a new 'Line Chart View'
	lineChartView1 = CreateView('XYChartView')
	lineChartView1.BottomAxisTitle = 'x'
	lineChartView1.LeftAxisTitle = 'Cp'
	lineChartView1.ShowRightAxisGrid = 0
	lineChartView1.SortByXAxis = 0
	AssignViewToLayout(view=lineChartView1, layout=layout2, hint=0)
	
	# create a new 'Plot Over Line'
	plotOverLine1 = PlotOverLine(registrationName='PlotOverLine1', Input=calculator2, Source='Line')
	plotOverLine1.Source.Point1 = [-1.374, 0, z_lnPltLvl]
	plotOverLine1.Source.Point2 = [1.603, 0, z_lnPltLvl]
	
	# show data in view
	plotOverLine1Display = Show(plotOverLine1, lineChartView1, 'XYChartRepresentation')
	plotOverLine1Display.CompositeDataSetIndex = [0]
	plotOverLine1Display.UseIndexForXAxis = 0
	plotOverLine1Display.XArrayName = 'Points_X'
	plotOverLine1Display.SeriesLabelPrefix = ''
	plotOverLine1Display.SeriesLabel = ['Cp','z = '+str(z_lnPltLvl)]
	plotOverLine1Display.SeriesVisibility = ['Cp']
	
	resampleWithDataset1 = ResampleWithDataset(registrationName='ResampleWithDataset1', SourceDataArrays=calculator2, DestinationMesh=mapvtu)
	plotOnIntersectionCurves1 = PlotOnIntersectionCurves(registrationName='PlotOnIntersectionCurves1', Input=resampleWithDataset1)
	plotOnIntersectionCurves1.SliceType = 'Plane'
	plotOnIntersectionCurves1.SliceType.Normal = [0.0, 1.0, 0.0]
	plotOnIntersectionCurves1Display = Show(plotOnIntersectionCurves1, lineChartView1, 'XYChartRepresentation')
	plotOnIntersectionCurves1Display.CompositeDataSetIndex = [0]
	plotOnIntersectionCurves1Display.UseIndexForXAxis = 0
	plotOnIntersectionCurves1Display.XArrayName = 'Points_X'
	plotOnIntersectionCurves1Display.SeriesLabelPrefix = ''
	plotOnIntersectionCurves1Display.SeriesLabel = ['Cp', 'At surface']
	plotOnIntersectionCurves1Display.SeriesColor = ['Cp', '0.2', '1.0', '0.0']
	plotOnIntersectionCurves1Display.SeriesVisibility = ['Cp']
	
	if importData:
		cpcsv = CSVReader(registrationName='cp_experiment.csv', FileName=[home+'/kode/simraStudies/cp.csv'])
		cpcsvDisplay = Show(cpcsv, lineChartView1, 'XYChartRepresentation')
		cpcsvDisplay.UseIndexForXAxis = 0
		cpcsvDisplay.XArrayName = 'x_exp'
		cpcsvDisplay.SeriesVisibility = ['Cp_exp']
		cpcsvDisplay.SeriesLineStyle = ['Cp_exp', '0']
		cpcsvDisplay.SeriesMarkerSize = ['Cp_exp', '8']
		cpcsvDisplay.SeriesMarkerStyle = ['Cp_exp', '4']
		cpcsvDisplay.SeriesColor = ['Cp_exp', '0.569', '0', '0.85']
		cpcsvDisplay.SeriesLabel = ['Cp_exp', 'Experiment']
	
	# create a new 'Plot Over Line'
	plotOverLine2 = PlotOverLine(registrationName='PlotOverLine2', Input=calculator2, Source='Line')
	plotOverLine2.Source.Point1 = plotOverLine1.Source.Point1
	plotOverLine2.Source.Point2 = plotOverLine1.Source.Point2
	Hide3DWidgets(proxy=plotOverLine2.Source)
	
	# Create a new 'Line Chart View'
	lineChartView2 = CreateView('XYChartView')
	lineChartView2.BottomAxisTitle = 'x'
	lineChartView2.LeftAxisTitle = 'Speed-up'
	AssignViewToLayout(view=lineChartView2, layout=layout2, hint=4)
	
	# show data in view
	plotOverLine2Display = Show(plotOverLine2, lineChartView2, 'XYChartRepresentation')
	plotOverLine2Display.CompositeDataSetIndex = [0]
	plotOverLine2Display.UseIndexForXAxis = 0
	plotOverLine2Display.SeriesVisibility = ['SpeedUp']
	plotOverLine2Display.SeriesLabel = ['SpeedUp','Simra']
	plotOverLine2Display.XArrayName = 'Points_X'
	
	if importData:
		cpcsv2 = CSVReader(registrationName='speed_up_experiment.csv', FileName=[home+'/kode/simraStudies/speed_up.csv'])
		cpcsvDisplay2 = Show(cpcsv2, lineChartView2, 'XYChartRepresentation')
		cpcsvDisplay2.UseIndexForXAxis = 0
		cpcsvDisplay2.XArrayName = 'x_exp'
		cpcsvDisplay2.SeriesVisibility = ['SpeedUp_exp']
		cpcsvDisplay2.SeriesLineStyle = ['SpeedUp_exp', '0']
		cpcsvDisplay2.SeriesMarkerSize = ['SpeedUp_exp', '8']
		cpcsvDisplay2.SeriesMarkerStyle = ['SpeedUp_exp', '4']
		cpcsvDisplay2.SeriesColor = ['SpeedUp_exp', '0.569', '0', '0.85']
		cpcsvDisplay2.SeriesLabel = ['SpeedUp_exp', 'Wind tunnel results']
	
	noProfiles = 8
	hints = [23, 24, 25, 26, 27, 28, 29, 30]
	plotOverLine3 = [''] * noProfiles
	lineChartView3 = [''] * noProfiles
	plotOverLine3Display = [''] * noProfiles
	x = np.linspace(-1,1.5,noProfiles)
	for i in range(0,noProfiles):
		plotOverLine3[i] = PlotOverLine(registrationName='PlotOverLine3_'+str(i), Input=calculator2,Source='Line')
		plotOverLine3[i].Source.Point1 = [x[i], 0.0, 0.0]
		plotOverLine3[i].Source.Point2 = [x[i], 0.0, 1.603]
		
		# Create a new 'Line Chart View'
		lineChartView3[i] = CreateView('XYChartView')
		lineChartView3[i].BottomAxisTitle = 'U/U_inf'
		lineChartView3[i].LeftAxisTitle = 'z'
		AssignViewToLayout(view=lineChartView3[i], layout=layout2, hint=hints[i])
		
		# show data in view
		plotOverLine3Display[i] = Show(plotOverLine3[i], lineChartView3[i], 'XYChartRepresentation')
		plotOverLine3Display[i].CompositeDataSetIndex = [0]
		plotOverLine3Display[i].UseIndexForXAxis = 0
		Hide3DWidgets(proxy=plotOverLine3[i].Source)
		
		plotOverLine3Display[i].XArrayName = 'SpeedUp'
		plotOverLine3Display[i].SeriesLabel = ['Points_Z', 'Simra']
		plotOverLine3Display[i].SeriesVisibility = ['Points_Z']
		plotOverLine3Display[i].SeriesColor = ['Point_Z', '0.0', '0.0', '0.0']
		
		if saveScreenShots:
			SaveScreenshot(outputPath+'Cp_'+str(i)+'.png', lineChartView3[i],
					FontScaling='Scale fonts proportionally',
					ImageQuality=100)
		
	if saveScreenShots:
		SaveScreenshot(outputPath+'Cp.png', lineChartView1,
				FontScaling='Scale fonts proportionally',
				ImageQuality=100)
		SaveScreenshot(outputPath+'SpeedUp.png', lineChartView2,
				FontScaling='Scale fonts proportionally',
				ImageQuality=100)
	
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
	
	mapvtuDisplay = Show(mapvtu, renderView3, 'UnstructuredGridRepresentation')
	mapvtuDisplay.Representation = 'Surface'
	mapvtuDisplay.AmbientColor = [0.6196078431372549, 0.6549019607843137, 0.7137254901960784]
	mapvtuDisplay.DiffuseColor = [0.6196078431372549, 0.6549019607843137, 0.7137254901960784]
	
	# create a new 'Stream Tracer'
	streamTracer1 = StreamTracer(registrationName='StreamTracer1', Input=calculator2, SeedType='Point Cloud')
	streamTracer1.Vectors = ['POINTS', 'u']
	streamTracer1.MaximumStepLength = 0.1
	streamTracer1.MaximumStreamlineLength = 4.0
	streamTracer1.SeedType.Center = [-1.473, 0.0, 0.0]
	streamTracer1.SeedType.Radius = 0.3
	streamTracer1.SeedType.NumberOfPoints = 200
	
	# show data in view
	streamTracer1Display = Show(streamTracer1, renderView3, 'GeometryRepresentation')
	streamTracer1Display.Representation = 'Surface'
	streamTracer1Display.ColorArrayName = ['POINTS', 'u']
	streamTracer1Display.RescaleTransferFunctionToDataRange(False, True)
	streamTracer1Display.SetScalarBarVisibility(renderView3, True)
	streamTracer1Display.SetRepresentationType('Wireframe')
	streamTracer1Display.RenderLinesAsTubes = 1
	streamTracer1Display.SetRepresentationType('Surface')
	streamTracer1Display.LineWidth = 3.0
	Hide3DWidgets(proxy=streamTracer1.SeedType)
	
	# current camera placement for renderView3
	renderView3.InteractionMode = 'Selection'
	renderView3.CameraPosition = [2.807134777489552, 1.6654291165153725, 1.2776926672144338]
	renderView3.CameraFocalPoint = [0.555600531530124, 0.14044477759713817, -0.08468961891830855]
	renderView3.CameraViewUp = [-0.35962741945592824, -0.26761818123475595, 0.8938951998126354]
	renderView3.CameraParallelScale = 2.0418274902645424
	
	if saveScreenShots:
		SaveScreenshot(outputPath+'streamTracer.png', renderView3,
				FontScaling='Scale fonts proportionally',
				TransparentBackground=1,
				ImageQuality=100)

RenderAllViews()
