from os.path import expanduser,exists
import sys
import numpy as np

import vtk
from splipy.io import G2
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy
import vtk.util.numpy_support as vtknp
import splipy.surface_factory as sfac
from splipy import SplineObject,BSplineBasis

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()
 
home = expanduser("~")
sys.path.insert(1, home+'/kode/paraUtils')
from utils import *
ref = 0
outputPath = 'SED_resultsFolder/'
saveScreenShots = True
importData = False
u_inf = 10.0 # is the freestream velocity of the fluid
p_inf = 101325 # is the static pressure in the freestream (i.e. remote from any disturbance)
rho_inf = 1.3 # is the freestream fluid density (Air at sea level and 15 °C is 1.225kg/m^3)
z_lnPltLvl = 0.6 # z coordinate of 1D plots
sqrtTKE_max = 3.0
SAVE_HIST = 10
noSteps=120
noSteps=30
noSteps=60
color = 'blue' # color of SINTEF logo
frameRate = 5
printLogo           = 0
plotLIC             = 1
plot1Dcurves        = 0 
plotStreamLines     = 0 
plotVolumeRendering = 1 
plotyplus           = 0 
plotOverTime        = 0 
plotMesh            = 0 
makeVideo           = 0
viewSizeSlice=[1920,520]
viewSizeSlice2=[2*1920,2*520]
scalarBarLength = 0.2
scalarBarLength2 = 0.92
# get animation scene
animationScene1 = GetAnimationScene()
#LoadPlugin(home+'/programs/paraview_build/lib/paraview-5.8/plugins/SurfaceLIC/SurfaceLIC.so', remote=False, ns=globals())
# get the time-keeper
timeKeeper1 = GetTimeKeeper()
fileName = outputPath+'SED_FILENAME.pvd'
CGL2fieldsMax = [1e-2, 1e-2, 1e-1, 1e-13, 1e-4, 1e-3]
CGL2fields_idx = [0,1,2,3,4,5]
#CGL2fields_idx = [1]
CGL2fields = ['Continuous global L2-projection |(u^*,p^*,pT^*)-(u^h,p^h,pT^h)|',\
              'Continuous global L2-projection |u^*-u^h|_H1',\
              'Continuous global L2-projection |div u^*|_L2',\
              'Continuous global L2-projection |pT^*-pT^h|_H1',\
              'Continuous global L2-projection |p^*-p^h|_L2',\
              'Continuous global L2-projection |s^*-s^h|_L2']
CGL2fields_LaTeX = ['Continuous global $L^2$-projection $|(u^*,p^*,pT^*)-(u^h,p^h,pT^h)|$',\
                    'Continuous global $L^2$-projection $|u^*-u^h|_{H^1}$',\
                    'Continuous global $L^2$-projection $|\\nabla\cdot u^*|_{L^2}$',\
                    'Continuous global $L^2$-projection $|pT^*-pT^h|_{H^1}$',\
                    'Continuous global $L^2$-projection $|p^*-p^h|_{L^2}$',\
                    'Continuous global $L^2$-projection $|s^*-s^h|_{L^2}$']
CGL2fieldsFileNames = ['CGL2error_uppT_L2',\
              'CGL2error_u_H1',\
              'CGL2error_div_u_L2',\
              'CGL2error_pT_H1',\
              'CGL2error_p_L2',\
              'CGL2error_s_L2']

# create a new 'XML Unstructured Grid Reader'
huntHill = PVDReader(registrationName='SED_FILENAME', FileName=fileName)
#huntHill = XMLUnstructuredGridReader(registrationName='fileName', FileName=fileNames)
#huntHill.CellArrayStatus = []
huntHill.PointArrays = ['ps', 'pts', 'u', 'tk', 'y+','vtef', 'terrain distance']
#huntHill.PointArrayStatus = ['ps', 'pT', 'pts', 'u', 'u']

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
layout1 = GetLayout()

# create a new 'Calculator'
calculator0 = Calculator(registrationName='CalculatorTKE', Input=huntHill)
calculator0.Function = 'sqrt(tk)'
calculator0.ResultArrayName = 'sqrtTKE'

sqrtTKELUT = GetColorTransferFunction('sqrtTKE')
sqrtTKEPWF = GetOpacityTransferFunction('sqrtTKE')
CGL2LUT = GetColorTransferFunction(CGL2fields[CGL2fields_idx[0]])
CGL2PWF = GetOpacityTransferFunction(CGL2fields[CGL2fields_idx[0]])
# create a new 'Calculator'
calculator1 = Calculator(registrationName='Calculator1', Input=calculator0)
calculator1.ResultArrayName = 'Cp'
denominator = 0.5*rho_inf*u_inf**2
calculator1.Function = '(ps-'+str(p_inf)+')/'+str(denominator)

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
slice1.Crinkleslice = 1

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=slice1.SliceType)
animationScene1.GoToLast()
LoadPalette(paletteName='WhiteBackground')
renderView1.OrientationAxesVisibility = 0

# create a new 'XML Unstructured Grid Reader'
topography = XMLStructuredGridReader(registrationName='Topography', FileName=[outputPath+'SED_FILENAME.vts'])
topography.TimeArray = 'None'

####################################################################################
## Layout 1 - Surface LIC plots
# create a new 'Clip'
if plotLIC:
    # split cell
    layout1.SplitVertical(0, 0.5)
    
    slice1Display = Show(slice1, renderView1, 'UnstructuredGridRepresentation')
    ColorBy(slice1Display, ('POINTS', 'u', 'X'))
    slice1Display.Representation = 'Surface'
    slice1Display.RescaleTransferFunctionToDataRange(True, False)
    slice1Display.SetScalarBarVisibility(renderView1, True)
    slice1Display.SetRepresentationType('Surface LIC')
    slice1Display.NumberOfSteps = 100
    slice1Display.ColorMode = 'Multiply'
    slice1Display.EnhanceContrast = 'Color Only'
    slice1Display.MapModeBias = 0.25
    slice1Display.HighColorContrastEnhancementFactor = 0.3
    slice1Display.RenderLinesAsTubes = 1
    slice1Display.EdgeColor = [1.0, 1.0, 1.0]
    slice1Display.LineWidth = 0.5
    slice1Display.PointSize = 1.0

    # current camera placement for renderView1
    renderView1.InteractionMode = '2D'
    renderView1.CameraPosition = [0.2880072459488367, -1.5010035486303253, 0.44080598448416514]
    renderView1.CameraFocalPoint = [0.2880072459488367, 0.0, 0.44080598448416514]
    renderView1.CameraViewUp = [0.0, 0.0, 1.0]
    renderView1.CameraParallelScale = 0.4517288466281378
    
    # get color legend/bar for uLUT in view renderView1
    uLUT = GetColorTransferFunction('u')
    uLUT.ApplyPreset('Parula', True)
    uLUT.InvertTransferFunction()
    uLUT.RescaleTransferFunction(-3.0,14)
    uLUT.AutomaticRescaleRangeMode = "Never"
    uLUT.RescaleOnVisibilityChange = 0
    uLUT.EnableOpacityMapping = 0
    
    uLUTColorBar = GetScalarBar(uLUT, renderView1)
    uLUTColorBar.ScalarBarLength = scalarBarLength2
    uLUTColorBar.WindowLocation = 'UpperRightCorner'
    uLUTColorBar.Title = 'x-component of u'
    uLUTColorBar.ComponentTitle = ''
   
    if printLogo:
        insertSINTEFlogo(renderView1,color)
    annotateTimeStep(calculator0,renderView1,'UpperLeftCorner', SAVE_HIST,color=[1.0, 1.0, 1.0])
    saveScreenShot(renderView1,outputPath+'surfaceLICside',saveScreenShots,viewSize=viewSizeSlice)
    ColorBy(slice1Display, ('CELLS', None))
    slice1Display.SetScalarBarVisibility(renderView1, False)
    HideScalarBarIfNotNeeded(uLUT, renderView1)
    slice1Display.SetRepresentationType('Surface With Edges')
    for i in CGL2fields_idx: 
        field = CGL2fields[i]
        ColorBy(slice1Display, ('CELLS', field))
        HideScalarBarIfNotNeeded(CGL2LUT, renderView1)
        slice1Display.RescaleTransferFunctionToDataRange(True, False)
        slice1Display.SetScalarBarVisibility(renderView1, True)
        CGL2LUT = GetColorTransferFunction(field)
        CGL2LUTColorBar = GetScalarBar(CGL2LUT, renderView1)
        CGL2LUTColorBar.Title = CGL2fields_LaTeX[i]
        CGL2LUTColorBar.Orientation = 'Vertical'
        CGL2LUTColorBar.WindowLocation = 'UpperRightCorner'
        CGL2LUTColorBar.ComponentTitle = ''
        CGL2LUTColorBar.ScalarBarLength = scalarBarLength2
        #CGL2LUTColorBar.TitleFontSize = 20
        #CGL2LUTColorBar.LabelFontSize = 20
        CGL2LUT.RescaleTransferFunction(0.0, CGL2fieldsMax[i])
        CGL2LUT.ApplyPreset('SINTEF1', True)
        saveScreenShot(renderView1,outputPath+'sliceSide_'+CGL2fieldsFileNames[i],saveScreenShots,viewSize=viewSizeSlice2)
        saveAnimation(renderView1,outputPath+'sliceSide_'+CGL2fieldsFileNames[i],noSteps,makeVideo,viewSize=viewSizeSlice2,frameRate=frameRate)
    
    
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
    slice2.SliceType.Normal = [0.0, 0.0, 1.0]
    slice2.SliceType.Origin = [0.0, 0.0, 0.2]
    Hide3DWidgets(proxy=slice2.SliceType)
    
    # show data in view
    slice2Display = Show(slice2, renderView2, 'UnstructuredGridRepresentation')
    ColorBy(slice2Display, ('POINTS', 'u', 'X'))
    slice2Display.Representation = 'Surface'
    slice2Display.SetRepresentationType('Surface LIC')
    slice2Display.NumberOfSteps = 100
    slice2Display.ColorMode = 'Multiply'
    slice2Display.EnhanceContrast = 'Color Only'
    slice2Display.HighColorContrastEnhancementFactor = 0.3
    slice2Display.LICIntensity = 0.8
    slice2Display.MapModeBias = 0.2
    slice2Display.RescaleTransferFunctionToDataRange(True, False)
    slice2Display.SetScalarBarVisibility(renderView2, True)# set active source
    slice2Display.RenderLinesAsTubes = 1
    slice2Display.EdgeColor = [1.0, 1.0, 1.0]
    slice2Display.LineWidth = 0.5
    slice2Display.PointSize = 1.0
    
    uLUTColorBar_1 = GetScalarBar(uLUT, renderView2)
    uLUTColorBar_1.ScalarBarLength = scalarBarLength2
    uLUTColorBar_1.WindowLocation = 'UpperRightCorner'
    uLUTColorBar_1.Title = 'x-component of u'
    uLUTColorBar_1.ComponentTitle = ''

    # current camera placement for renderView2
    renderView2.InteractionMode = '2D'
    renderView2.CameraPosition = [0.2938256768798605, -0.06985340760893753, 1.9741519036963915]
    renderView2.CameraFocalPoint = [0.2938256768798605, 0.0, 0.8015]
    renderView2.CameraParallelScale = 0.4517288466281378
    
    if printLogo:
        insertSINTEFlogo(renderView2,color)

    annotateTimeStep(calculator0,renderView2,'UpperLeftCorner', SAVE_HIST,color=[1.0, 1.0, 1.0])
    saveScreenShot(renderView2,outputPath+'surfaceLICtop',saveScreenShots,viewSize=viewSizeSlice)
    ColorBy(slice2Display, ('CELLS', None))
    slice2Display.SetScalarBarVisibility(renderView2, False)
    HideScalarBarIfNotNeeded(uLUT, renderView2)
    slice2Display.SetRepresentationType('Surface With Edges')
    slice2.Crinkleslice = 1
    for i in CGL2fields_idx: 
        field = CGL2fields[i]
        ColorBy(slice2Display, ('CELLS', field))
        HideScalarBarIfNotNeeded(CGL2LUT, renderView2)
        slice2Display.RescaleTransferFunctionToDataRange(True, False)
        slice2Display.SetScalarBarVisibility(renderView2, True)
        CGL2LUT = GetColorTransferFunction(field)
        CGL2LUTColorBar = GetScalarBar(CGL2LUT, renderView2)
        CGL2LUTColorBar.Title = CGL2fields_LaTeX[i]
        CGL2LUTColorBar.Orientation = 'Vertical'
        CGL2LUTColorBar.WindowLocation = 'UpperRightCorner'
        CGL2LUTColorBar.ComponentTitle = ''
        CGL2LUTColorBar.ScalarBarLength = scalarBarLength2
        #CGL2LUTColorBar.TitleFontSize = 20
        #CGL2LUTColorBar.LabelFontSize = 20    
        CGL2LUT.RescaleTransferFunction(0.0, CGL2fieldsMax[i])
        CGL2LUT.ApplyPreset('SINTEF1', True)
        saveScreenShot(renderView2,outputPath+'sliceTop_'+CGL2fieldsFileNames[i],saveScreenShots,viewSize=viewSizeSlice2)
        saveAnimation(renderView2,outputPath+'sliceTop_'+CGL2fieldsFileNames[i],noSteps,makeVideo,viewSize=viewSizeSlice2,frameRate=frameRate)

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
    
    # create a new 'Plot Over Line'
    plotOverLine1 = PlotOverLine(registrationName='PlotOverLine1', Input=calculator2, Source='Line')
    plotOverLine1.Source.Point1 = [-1.374, 0, z_lnPltLvl]
    plotOverLine1.Source.Point2 = [1.603, 0, z_lnPltLvl]
    #Hide3DWidgets(proxy=plotOverLine1.Source)
    
    # Create a new 'Line Chart View'
    lineChartView1 = CreateView('XYChartView')
    lineChartView1.BottomAxisTitle = 'x'
    lineChartView1.LeftAxisTitle = 'Cp'
    AssignViewToLayout(view=lineChartView1, layout=layout2, hint=0)
    
    # show data in view
    plotOverLine1Display = Show(plotOverLine1, lineChartView1, 'XYChartRepresentation')
    plotOverLine1Display.CompositeDataSetIndex = [0]
    plotOverLine1Display.UseIndexForXAxis = 0
    plotOverLine1Display.SeriesVisibility = ['Cp']
    plotOverLine1Display.SeriesLabel = ['Cp','z = '+str(z_lnPltLvl)]
    plotOverLine1Display.XArrayName = 'Points_X'
    plotOverLine1Display.SeriesPlotCorner = ['Cp', '0']
    plotOverLine1Display.SeriesLineStyle = ['Cp', '1']
    plotOverLine1Display.SeriesLineThickness = ['Cp', '2']
    plotOverLine1Display.SeriesMarkerStyle = ['Cp', '0']
    plotOverLine1Display.SeriesMarkerSize = ['Cp', '4']

    calculatorTop = Calculator(registrationName='Elevated topography', Input=topography)
    calculatorTop.Function = 'coordsX*iHat+coordsY*jHat+(coordsZ+0.002)*kHat'
    calculatorTop.CoordinateResults = 1

    resampleWithDataset1 = ResampleWithDataset(registrationName='ResampleWithDataset1', SourceDataArrays=calculator2, DestinationMesh=calculatorTop)
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
    
    saveScreenShot(lineChartView1,outputPath+'Cp',saveScreenShots,viewSize=[955,522])
    
    # Create a new 'Line Chart View'
    lineChartView2 = CreateView('XYChartView')
    lineChartView2.BottomAxisTitle = 'x'
    lineChartView2.LeftAxisTitle = 'Speed-up'
    AssignViewToLayout(view=lineChartView2, layout=layout2, hint=4)
    
    # show data in view
    plotOverLine2Display = Show(plotOverLine1, lineChartView2, 'XYChartRepresentation')
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
    
    saveScreenShot(lineChartView2,outputPath+'SpeedUp',saveScreenShots,viewSize=[955,522])
    
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
        
        saveScreenShot(lineChartView3[i],outputPath+'Cp_'+str(i),saveScreenShots,viewSize=[238,522])
     
    
    
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
    
    topographyDisplay = Show(topography, renderView3, 'UnstructuredGridRepresentation')
    topographyDisplay.Representation = 'Surface'
    topographyDisplay.AmbientColor = [0.6196078431372549, 0.6549019607843137, 0.7137254901960784]
    topographyDisplay.DiffuseColor = [0.6196078431372549, 0.6549019607843137, 0.7137254901960784]
    
    # create a new 'Stream Tracer'
    streamTracer1 = StreamTracer(registrationName='StreamTracer1', Input=calculator0, SeedType='Point Cloud')
    streamTracer1.Vectors = ['POINTS', 'u']
    streamTracer1.MaximumStepLength = 0.1
    streamTracer1.MaximumStreamlineLength = 4.0
    streamTracer1.SeedType.Center = [-1.473, 0.0, 0.0]
    streamTracer1.SeedType.Radius = 0.3
    streamTracer1.SeedType.NumberOfPoints = 200
    
    # show data in view
    streamTracer1Display = Show(streamTracer1, renderView3, 'GeometryRepresentation')
    streamTracer1Display.Representation = 'Surface'
    streamTracer1Display.ColorArrayName = ['POINTS', 'sqrtTKE']
    streamTracer1Display.SetRepresentationType('Surface')
    streamTracer1Display.RenderLinesAsTubes = 1
    streamTracer1Display.LineWidth = 3.0
    
    Hide3DWidgets(proxy=streamTracer1.SeedType)
    
    sqrtTKELUTColorBar_3 = GetScalarBar(sqrtTKELUT, renderView3)
    sqrtTKELUTColorBar_3.Title = 'Turbulence, $\sqrt{k}$'
    sqrtTKELUTColorBar_3.Orientation = 'Vertical'
    sqrtTKELUTColorBar_3.WindowLocation = 'UpperRightCorner'
    sqrtTKELUTColorBar_3.ComponentTitle = ''
    sqrtTKELUTColorBar_3.ScalarBarLength = scalarBarLength
    annotateTimeStep(calculator0,renderView3,'UpperLeftCorner', SAVE_HIST)

    renderView3.CameraPosition = [1.877390847356428, -2.5021796320272767, 1.5784369498753288]
    renderView3.CameraFocalPoint = [0.45113064046067347, -0.3493730968624169, -0.028508921216712002]
    renderView3.CameraViewUp = [0.0, 0.0, 1.0]
    renderView3.CameraParallelScale = 2.0418274902645424    

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
    
    topographyDisplay = Show(topography, renderView4, 'UnstructuredGridRepresentation')
    topographyDisplay.Representation = 'Surface'
    topographyDisplay.AmbientColor = [0.6196078431372549, 0.6549019607843137, 0.7137254901960784]
    topographyDisplay.DiffuseColor = [0.6196078431372549, 0.6549019607843137, 0.7137254901960784]
    
    # create a new 'Stream Tracer'
    calculator0Display = Show(calculator0, renderView4, 'UnstructuredGridRepresentation')
    calculator0Display.Representation = 'Volume'
    calculator0Display.ColorArrayName = ['POINTS', 'sqrtTKE']

    sqrtTKELUTColorBar_4 = GetScalarBar(sqrtTKELUT, renderView4)
    sqrtTKELUTColorBar_4.Title = 'Turbulence, $\sqrt{k}$'
    sqrtTKELUTColorBar_4.Orientation = 'Vertical'
    sqrtTKELUTColorBar_4.WindowLocation = 'UpperRightCorner'
    sqrtTKELUTColorBar_4.ComponentTitle = ''
    sqrtTKELUTColorBar_4.ScalarBarLength = scalarBarLength
    annotateTimeStep(calculator0,renderView4,'UpperLeftCorner', SAVE_HIST)
    
    renderView4.CameraPosition = [1.877390847356428, -2.5021796320272767, 1.5784369498753288]
    renderView4.CameraFocalPoint = [0.45113064046067347, -0.3493730968624169, -0.028508921216712002]
    renderView4.CameraViewUp = [0.0, 0.0, 1.0]
    renderView4.CameraParallelScale = 2.0418274902645424    
    
## Layout 5 - Plot y+
if plotyplus:
    CreateLayout('Layout #5')
    layout5 = GetLayoutByName("Layout #5")
    renderView5 = CreateView('RenderView')
    renderView5.AxesGrid = 'GridAxes3DActor'
    renderView5.StereoType = 'Crystal Eyes'
    renderView5.CameraFocalDisk = 1.0
    renderView5.OrientationAxesVisibility = 0
    renderView5.InteractionMode = '2D'
    renderView5.CameraPosition = [0.20510506812461812, -0.025711924662068383, -3.422733425372001]
    renderView5.CameraFocalPoint = [0.20510506812461812, -0.025711924662068383, 1.53621]
    renderView5.CameraParallelScale = 2.2737436287800374
    AssignViewToLayout(view=renderView5, layout=layout5, hint=0)


    yplusFileNames = []
    i = 0
    while exists(outputPath+'yplus'+str(i)+'.vtk'):
        yplusFileNames += [outputPath+'yplus'+str(i)+'.vtk']
        i += 1

    yplus = LegacyVTKReader(registrationName='yplus', FileNames=yplusFileNames)   
    
    yplusDisplay = Show(yplus, renderView5, 'StructuredGridRepresentation')
    yplusDisplay.Representation = 'Surface'
    ColorBy(yplusDisplay, ('POINTS', 'y+'))
    yLUT = GetColorTransferFunction('y')
    yLUT.ApplyPreset('SINTEF', True)
    yLUT.RescaleTransferFunction(0.0,300)
    yLUT.AutomaticRescaleRangeMode = "Never"
    yLUT.RescaleOnVisibilityChange = 0
    yLUT.EnableOpacityMapping = 0

    yLUTColorBar = GetScalarBar(yLUT, renderView5)
    yLUTColorBar.Title = '$y_p^+$'
    yLUTColorBar.Orientation = 'Vertical'
    yLUTColorBar.WindowLocation = 'UpperRightCorner'
    yLUTColorBar.ComponentTitle = ''
    yLUTColorBar.ScalarBarLength = scalarBarLength
    annotateTimeStep(calculator0,renderView5,'UpperLeftCorner', SAVE_HIST)
    
    renderView5.InteractionMode = '3D'
    renderView5.CameraPosition = [0.11449999999999994, 0.0, 4.222279913606204]
    renderView5.CameraFocalPoint = [0.11449999999999994, 0.0, 0.11843098355000001]
    renderView5.CameraParallelScale = 1.881671065342127
    saveScreenShot(renderView5,outputPath+'yplus',saveScreenShots)

if plotOverTime:
    CreateLayout('Layout #6')
    layout6 = GetLayoutByName("Layout #6")
    quartileChartView1 = CreateView('QuartileChartView')
    quartileChartView1.BottomAxisTitle = 'Iteration number'
    quartileChartView1.LeftAxisTitle = 'TKE $[m^2/s^2]$'
    AssignViewToLayout(view=quartileChartView1, layout=layout6, hint=0)
    noPoints = 5
    cpcsv = [''] * noPoints
    cpcsvDisplay = [''] * noPoints
    colors = [[0.8,0.54,0.25],[0.6,0.2,0.1],[0.4,0.5,0.7],[0.2,0.6,0.1],[0.0,0.6,0.8]]
    for i in range(0,noPoints):
        cpcsv[i] = CSVReader(registrationName='SED_FILENAME_Point'+str(i+1), FileName=[outputPath+'SED_FILENAME_Point'+str(i+1)+'.csv'])
        cpcsvDisplay[i] = Show(cpcsv[i], quartileChartView1, 'XYChartRepresentation')
        cpcsvDisplay[i].UseIndexForXAxis = 1
        #cpcsvDisplay.XArrayName = ''
        #cpcsvDisplay.SeriesLineStyle = ['TKE', '0']
        #cpcsvDisplay[i].SeriesMarkerSize = ['TKE', '8']
        #cpcsvDisplay[i].SeriesMarkerStyle = ['TKE', '4']
        cpcsvDisplay[i].SeriesColor = ['TKE', str(colors[i][0]), str(colors[i][1]), str(colors[i][2])]
        cpcsvDisplay[i].SeriesLabel = ['TKE', 'Point at ('+str((i+1)/5)+',0,0.1)']
        cpcsvDisplay[i].SeriesVisibility = ['TKE']
    
    saveScreenShot(quartileChartView1,outputPath+'TKE',saveScreenShots)
    #pointSource1 = PointSource(registrationName='PointSource1')
    #pointSource1.Center = [0.5, 0.0, 0.1]
    #resampleWithDataset1 = ResampleWithDataset(registrationName='ResampleWithDataset1', SourceDataArrays=calculator2, DestinationMesh=pointSource1)
    #plotDataOverTime1 = PlotDataOverTime(registrationName='PlotDataOverTime1', Input=resampleWithDataset1)
    #plotDataOverTime1Display = Show(plotDataOverTime1, quartileChartView1, 'QuartileChartRepresentation')
    #plotDataOverTime1Display.SeriesVisibility = ['sqrtTKE (stats)']
    #quartileChartView1.LeftAxisTitle = '$\\sqrt{k}$'

    #quartileChartView1.BottomAxisTitle = 'Time [s]'
    #quartileChartView1.LeftAxisTitle = '$\\sqrt{k}$ [m/s]'
    #quartileChartView1.HideTimeMarker = 1
    #saveScreenShot(quartileChartView1,outputPath+'sqrtTKE',saveScreenShots)

if plotMesh:
    CreateLayout('Layout #7')
    layout7 = GetLayoutByName("Layout #7")
    renderView7 = CreateView('RenderView')
    renderView7.AxesGrid = 'GridAxes3DActor'
    renderView7.StereoType = 'Crystal Eyes'
    renderView7.CameraFocalDisk = 1.0
    renderView7.OrientationAxesVisibility = 0
    AssignViewToLayout(view=renderView7, layout=layout7, hint=0)

    huntHillDisplay = Show(huntHill, renderView7, 'UnstructuredGridRepresentation')
    huntHillDisplay.SetRepresentationType('Surface With Edges')
    huntHillDisplay.EdgeColor = [0.0, 0.0, 0.0]
    huntHillDisplay.RenderLinesAsTubes = 1
    renderView7.CameraPosition = [0.11449999999999994, 0.0, 5.920225294179859]
    renderView7.CameraFocalPoint = [0.11449999999999994, 0.0, 0.79806919]
    renderView7.InteractionMode = '2D'
    renderView7.CameraParallelScale = 1.150455376463919
    saveScreenShot(renderView7,outputPath+'sliceTop_mesh',saveScreenShots,viewSize=[1700,1080])

    Hide(huntHill)
    sliceMeshDisplay = Show(slice1, renderView7, 'UnstructuredGridRepresentation')
    sliceMeshDisplay.SetRepresentationType('Surface With Edges')
    sliceMeshDisplay.ColorMode = 'Multiply'
    sliceMeshDisplay.EnhanceContrast = 'Color Only'
    sliceMeshDisplay.MapModeBias = 0.25
    sliceMeshDisplay.HighColorContrastEnhancementFactor = 0.3
    sliceMeshDisplay.RenderLinesAsTubes = 1
    sliceMeshDisplay.EdgeColor = [0.0, 0.0, 0.0]
    sliceMeshDisplay.LineWidth = 0.5
    sliceMeshDisplay.PointSize = 1.0
    ColorBy(sliceMeshDisplay, None)
    renderView7.CameraPosition = [0.11450004577636719, -3.1793793691718686, 0.7980866522993892]
    renderView7.CameraFocalPoint = [0.11450004577636719, 0.0, 0.7980866522993892]
    renderView7.CameraViewUp = [0.0, 0.0, 1.0]
    renderView7.InteractionMode = '2D'
    renderView7.CameraParallelScale = 0.839423899387907
    saveScreenShot(renderView7,outputPath+'sliceSide_mesh',saveScreenShots,viewSize=[1920,1080])

if plotStreamLines or plotVolumeRendering:
    # get color transfer function/color map for 'sqrtTKE'
    sqrtTKELUT.AutomaticRescaleRangeMode = "Never"
    sqrtTKELUT.RescaleOnVisibilityChange = 0
    sqrtTKELUT.EnableOpacityMapping = 0
    sqrtTKELUT.RescaleTransferFunction(0.0, sqrtTKE_max)
    sqrtTKELUT.ApplyPreset('SINTEF1', True)
    
    # get opacity transfer function/opacity map for 'sqrtTKE'
    sqrtTKEPWF.AllowDuplicateScalars = 1
    sqrtTKEPWF.ScalarRangeInitialized = 1
    sqrtTKEPWF.RescaleTransferFunction(0.0, sqrtTKE_max)
    sqrtTKEPWF.Points = [0.0, 0.0, 0.5, 0.0, 1.6, 0.0, 0.5, 0.0, 4.0, 1.0, 0.5, 0.0]

if plotStreamLines:
    if printLogo:
        insertSINTEFlogo(renderView3,color)
    saveScreenShot(renderView3,outputPath+'streamTracer',saveScreenShots)
    saveAnimation(renderView3,outputPath+'streamTracer',noSteps,makeVideo)

if plotVolumeRendering:
    if printLogo:
        insertSINTEFlogo(renderView4,color)

    saveScreenShot(renderView4,outputPath+'volumeRendering',saveScreenShots)
    saveAnimation(renderView4,outputPath+'volumeRendering',noSteps,makeVideo)

RenderAllViews()
