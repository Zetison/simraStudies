from os.path import expanduser
import sys
from glob import glob
import numpy as np
from scipy.spatial.transform import Rotation as Rot
import copy
from pathlib import Path
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
simraResultsFolder = home+'/results/simra/Sula/'
outputPath = simraResultsFolder+caseName+'/'
topoRes = 'SED_topoRes' 

# Set the time
openFoamResultsFolder = home+'/results/openfoam/Sula/'
runAll = True 
if runAll:
    fileNamesOrg = [
    openFoamResultsFolder+'2020111906_OF_unStdy.pvd',
    openFoamResultsFolder+'2020111906_OF_stdy.pvd',
    openFoamResultsFolder+'2020111906_OF_stdyIso.pvd',
    simraResultsFolder+'met/cont_met.pvd',
    ]
    metFiles = glob(simraResultsFolder+"met_new/2*.pvd")

    fileNamesOrg.extend(metFiles)

    caseNamesOrg = [
    '2020111906_OF_unStdy',
    '2020111906_OF_stdy',
    '2020111906_OF_stdyIso',
    '2020111906_met',
    ]
    caseNamesOrg.extend([f[-17:-3] for f in metFiles])
    #indices = [2,5,6,7,8]
    indices = range(9) 
    indices = [0]
    fileNames = [fileNamesOrg[i] for i in indices]
    caseNames = [caseNamesOrg[i] for i in indices]
else:
    fileNames = [outputPath+fileName+'.pvd']
    caseNames = [caseName]

colorLogo = 'white' # Default color of SINTEF logo
originx=-200.0
originy=6899800.0
topographyFileName   = simraResultsFolder+'Sula'+topoRes+'.vts'
textureFileName_topo = simraResultsFolder+'Sula'+topoRes+'_topo.png'
textureFileName_NIB  = simraResultsFolder+'Sula'+topoRes+'.png'
windCase = 2 
velProfileMax_Z = 120
z_proj = 5000 # Project cone to this height
finalTime = 11.7006
sqrtTKE_max = 5.0
h_max = 4000
height_max = 1600 #max height of height colormap
angleOfAttack_West = 3.0
angleOfAttack_East = 3.0
runwayWidth = 150.0
runwayHeight = 30.0
runwayHeight = 10.0
bridgeBaseHeight = 50.0
plotRunway               = 0 
plotTakeOffLines         = 0 
plotBridgeAndSensors     = 1 
plotUniverseBackground   = 0 

plotLIC                  = 0  
plotStreamLines          = 0
plotVolumeRendering      = 0 
plotIPPCmapsHorizontal   = 0
plotIPPCmapsVertical     = 0
plotOverTime             = not runAll
plotVelocityProfiles     = runAll 

makeVideo                = 0
saveScreenShots          = 1 
useTransparentBackground = 0
plotError                = 0 

bridge = SED_BRIDGE
bridgeHeight = 20
# Use i.e. norgeskart.no and "Vis koordinater" to get UTM coordinates (NORD,ØST,height)
#runwayEndWest = latLonUTM("623323.41","0060532.89",35.3)
#runwayEndEast = latLonUTM("623351.26","0060740.31",49.2)
mastLoc = np.genfromtxt('../mastLoc.csv', delimiter=',')
if bridge == 1:
    runwayEndWest = mastLoc[0][:3]
    runwayEndEast = mastLoc[1][:3]
elif bridge == 2:
    runwayEndWest = mastLoc[2][:3]
    runwayEndEast = mastLoc[3][:3]

#runwayEndWest = np.array([346520-originx,6920740-originy,4.4])
#runwayEndEast = np.array([351140-originx,6922074-originy,3.6])
#runwayEndWest = np.array([549213-0.4-originx, 6943810+0.90-originy,0.0])
#runwayEndEast = np.array([549213-originx,6943810-originy,0.0])
SAVE_HIST = 20

viewSize = [1920, 1080]
scalarBarLength = 0.26
streamLines_z = 10
horizontalSlice_z = 10
glyphScale1 = 30
if windCase == 1:
    glyphScale2 = 0.03
else:
    glyphScale2 = 0.03

timeStepAnnotation       = makeVideo
# Extract UTM-33 coordinates from rwyCoords.png file (available at https://ais.avinor.no/no/AIP/View/95/2020-11-05-AIRAC/html/index-no-NO.html). The coordinates in rwyCoords.png are given in Lat-lon Degrees Minutes seconds format. Here, "THR COORD" is the start of the RWY (runway), and "RWY end COORD" is its end
runwayColor      = [0.0,0.0,0.0]
takeOffLineColor = [0.0,0.0,0.0]

runwayColor      = [1.0,0.0,0.0]

P_sw = np.array([25000-originx,6927800-originy,streamLines_z])
P_se = np.array([64999.58-originx,6927800-originy,streamLines_z])
P_nw = np.array([25000-originx,6967800-originy,streamLines_z])
P_ne = np.array([64999.5156-originx,6967800.5-originy,streamLines_z])
if windCase == 1:
    streamTracerPoint1 = P_nw
    streamTracerPoint2 = P_se
else:
    streamTracerPoint1 = [25749-originx,6951933-originy,streamLines_z]
    streamTracerPoint2 = [39478-originx,6966147-originy,streamLines_z]

runwayDir = runwayEndEast-runwayEndWest
runwayDir = runwayDir/np.linalg.norm(runwayDir)
runwayNormal = np.cross(runwayDir,[0,0,1])
runwayNormal = runwayNormal/np.linalg.norm(runwayNormal)
runwayCenter = (runwayEndEast+runwayEndWest)/2
runwayLength = np.linalg.norm(runwayEndWest-runwayEndEast)
runwayRotVector = np.array([0,0,np.degrees(np.arctan2(runwayDir[1],runwayDir[0]))])

LoadPalette(paletteName='WhiteBackground')
renderView1 = GetActiveViewOrCreate('RenderView')
renderView1.OrientationAxesVisibility = 0

if plotUniverseBackground:
    #convert -crop 50%x100% lightBlueMarble2.jpg converted_lBM2.jpg
    import universeBackground 


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

def annotateDateStep(obj,renderview,RegistrationName,location='LowerRightCorner'):
    annotateTimeFilter1 = AnnotateTimeFilter(registrationName=RegistrationName, Input=obj)
    annotateTimeFilter1.Format = '%10d'
    annotateTimeFilter1Display = Show(annotateTimeFilter1, renderview, 'ChartTextRepresentation')
    annotateTimeFilter1Display.FontSize = 8 
    return annotateTimeFilter1Display

def annotateTimeStep(obj,renderview,RegistrationName='Time annotation', location='UpperLeftCorner'):
    pythonAnnotation = PythonAnnotation(registrationName=RegistrationName, Input=obj)
    pythonAnnotation.ArrayAssociation = 'Point Data'
    pythonAnnotation.Expression = '"%10d" % time_value'
 #   pythonAnnotationDisplay = Show(pythonAnnotation, renderview, 'TextSourceRepresentation')
    pythonAnnotationDisplay = Show(pythonAnnotation, renderview, 'ChartTextRepresentation')
    #pythonAnnotationDisplay.WindowLocation = location 
    pythonAnnotationDisplay.FontSize = 5

labels = ['SulaNW (Kvitneset) - 1','SulaNW (Kvitneset) - 2','SulaNW (Kvitneset) - 3',
          'SulaNE (Trælbodneset) - 1','SulaNE (Trælbodneset) - 2','SulaNE (Trælbodneset) - 3',
          'SulaSW (Langeneset) - 1','SulaSW (Langeneset) - 2','SulaSW (Langeneset) - 3','SulaSW (Langeneset) - 4',
          'SulaSE (Kårsteinen) - 1','SulaSE (Kårsteinen) - 2','SulaSE (Kårsteinen) - 3',
          'Midspan bridge 1', 'Midspan bridge 2']
mastLabel = ['SulaNW (Kvitneset)', 'SulaNE (Trælboneset)', 'SulaSW (Langeneset)', 'SulaSE (Kårsteinen)']
points = np.genfromtxt('../POINT_FILE', delimiter=' ', skip_header=True)
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
        lineSource2[i] = Line(registrationName=labels[i]+' - boom')
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

mastLine = [''] * noMasts
for i in range(noMasts):
    mastLine[i] = Line(registrationName=mastLabel[i]+' line')
    mastLine[i].Point1 = [mastLoc[i][0],mastLoc[i][1],-100]
    mastLine[i].Point2 = [mastLoc[i][0],mastLoc[i][1],velProfileMax_Z]
    mastLine[i].Resolution = 1000
    Hide3DWidgets(proxy=mastLine[i])
    
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

def setVisibility(timestepvalues,visibility,totSteps):
    #visibilityTrack = GetAnimationTrack('Visibility',proxy=reader1)
    visibilityTrack = GetAnimationTrack('Visibility')
    #endStep = startStep+noStep+1
    #timeStart = startStep/totSteps
    #timeEnd = endStep/totSteps
    keyframeArr = []
    for i in range(totSteps):
        keyframe0 = CompositeKeyFrame()
        keyframe0.Interpolation = 'Boolean'
        keyframe0.KeyTime = timestepValues[i] 
        keyframe0.KeyValues = [visibility[i]]
        keyframeArr.append(keyframe0)
    #keyframe1 = CompositeKeyFrame()
    #keyframe1.Interpolation = 'Boolean'
    #keyframe1.KeyTime = timeStart
    #keyframe1.KeyValues = [1]
    #keyframe2 = CompositeKeyFrame()
    #keyframe2.KeyTime = timeEnd 
    #keyframe2.KeyValues = [0]
    visibilityTrack.KeyFrames = keyframeArr
    #if startStep == 0:
    #else:
    #    visibilityTrack.KeyFrames = [keyframe0, keyframe1, keyframe2]

topography = XMLStructuredGridReader(registrationName='Topography',FileName=topographyFileName)
topography.PointArrayStatus = ['TCoords_']
layout1 = GetLayout()

calculator2 = Calculator(registrationName='Calculator2', Input=topography)
calculator2.Function = 'coordsZ'
calculator2.ResultArrayName = 'height'
heightLUT = GetColorTransferFunction('height')
heightPWF = GetOpacityTransferFunction('height')

sliceTopography = Slice(registrationName='Slice topography', Input=calculator2)
sliceTopography.SliceType = 'Plane'
sliceTopography.SliceType.Normal = runwayNormal 
sliceTopography.SliceType.Origin = runwayCenter
Hide3DWidgets(proxy=sliceTopography.SliceType)
    
noFiles = len(fileNames)
pvdResults = [''] * noFiles
timestepValues = []
timeStamps = [''] * noFiles
for i_f in range(0,noFiles):
    # get the time-keeper
    pvdResults[i_f] = PVDReader(registrationName=fileNames[i_f], FileName=fileNames[i_f])
    timeStamps[i_f] = pvdResults[i_f].TimestepValues
    timestepValues.extend(timeStamps[i_f])

timestepValues = np.sort(np.unique(np.array(timestepValues)).astype(int))
totSteps = len(timestepValues)
visibility = [''] * noFiles
for i_f in range(noFiles):
    visibility[i_f] = np.isin(timestepValues,pvdResults[i_f].TimestepValues)


if plotVelocityProfiles:
    hints = [3,4,5,6]
    dataTypes = ['filtered','raw','rawMid']
    colorsData = [[0,0,0],[0.4,0.4,0.4],[0.8,0.8,0.8]] 
    mastNames = ['Kvitneset', 'Traelboneset','Langeneset','Kaarsteinen']
    layoutNames = ['VelocityProfiles', 'WindDirProfiles', 'AoAprofiles']
    xArrayNames = ['u_mag', 'u_mag', 'meandir']
    bottomAxisRangeMinimums = [0.0, 0.0, -90.0]
    bottomAxisRangeMaximums = [30.0, 360.0, 90.0]
    bottomAxisTitles = ['$u$ [m/s] (magnitude)', 'Wind dir', 'Angle of Attack']
    noPlots = len(layoutNames)
    noDataTypes = len(dataTypes)
    vpcsv = [''] * noMasts
    plotData1 = [''] * noMasts
    for i in range(0,noMasts):
        vpcsv[i] = [''] * noDataTypes 
        plotData1[i] = [''] * noDataTypes 
        for j in range(3):
            fileName = simraResultsFolder+'measurements/'+dataTypes[j]+'/'
            fileNames = [fileName+'VelocityProfile_'+mastNames[i]+'_'+str(time)+'.vtu' for time in timestepValues]
            writePVD(fileName+mastNames[i],fileNames,timestepValues)
            vpcsv[i][j] = PVDReader(registrationName=mastNames[i], FileName=fileName+mastNames[i]+'.pvd')
            plotData1[i][j] = PlotData(registrationName=mastNames[i]+'_'+dataTypes[j], Input=vpcsv[i][j])

    quartileChartView = [''] * noPlots 
    layouts1D = [''] * noPlots 
    vpcsvDisplay = [''] * noPlots 
    for i_l in range(noPlots):
        layouts1D[i_l] = CreateLayout(layoutNames[i_l])
        layouts1D[i_l].SplitHorizontal(0, 0.5)
        layouts1D[i_l].SplitHorizontal(1, 0.5)
        layouts1D[i_l].SplitHorizontal(2, 0.5)
        quartileChartView[i_l] = [''] * noMasts
        vpcsvDisplay[i_l] = [''] * noMasts
        for i in range(0,noMasts):
            quartileChartView[i_l][i] = CreateView('QuartileChartView')
            quartileChartView[i_l][i].BottomAxisTitle = bottomAxisTitle
            quartileChartView[i_l][i].LeftAxisTitle = '$z$ [m]'
            quartileChartView[i_l][i].ChartTitle = mastLabel[i]
            quartileChartView[i_l][i].LeftAxisUseCustomRange = 1
            quartileChartView[i_l][i].LeftAxisRangeMinimum = 0.0
            quartileChartView[i_l][i].LeftAxisRangeMaximum = velProfileMax_Z
            quartileChartView[i_l][i].BottomAxisUseCustomRange = 1
            quartileChartView[i_l][i].BottomAxisRangeMinimum = bottomAxisRangeMinimums[i_l]
            quartileChartView[i_l][i].BottomAxisRangeMaximum = bottomAxisRangeMaximums[i_l]
            quartileChartView[i_l][i].LegendLocation = 'TopRight'
            quartileChartView[i_l][i].LegendFontSize = 8 
            quartileChartView[i_l][i].ViewSize = [viewSize[0]//4, viewSize[1]]
            AssignViewToLayout(view=quartileChartView[i_l][i], layout=layouts1D[i_l], hint=hints[i])
            vpcsvDisplay[i_l][i] = [''] * noDataTypes
            for j in range(3):
                vpcsvDisplay[i_l][i][j] = Show(plotData1[i][j], quartileChartView[i_l][i], 'XYChartRepresentation')
                vpcsvDisplay[i_l][i][j].UseIndexForXAxis = 0
                vpcsvDisplay[i_l][i][j].XArrayName = xArrayNames[i_l] 
                vpcsvDisplay[i_l][i][j].SeriesLineStyle = ['coordsZ', '1']
                vpcsvDisplay[i_l][i][j].SeriesMarkerSize = ['coordsZ', '8']
                vpcsvDisplay[i_l][i][j].SeriesMarkerStyle = ['coordsZ', '4']
                vpcsvDisplay[i_l][i][j].SeriesColor = ['coordsZ', str(colorsData[j][0]), str(colorsData[j][1]), str(colorsData[j][2])]
                vpcsvDisplay[i_l][i][j].SeriesLabel = ['coordsZ', 'Exp. '+dataTypes[j]]
                vpcsvDisplay[i_l][i][j].SeriesVisibility = ['coordsZ']

calculator0 = [''] * noFiles
calculator1 = [''] * noFiles
calculator1Display = [''] * noFiles
calculatorIPPC = [''] * noFiles
slice1 = [''] * noFiles
slice2 = [''] * noFiles
slice1Display = [''] * noFiles
slice2Display = [''] * noFiles
proj_u = [''] * noFiles
planeSlice = [''] * noFiles
resampledSlice = [''] * noFiles
streamTracer1 = [''] * noFiles
streamTracer1Display = [''] * noFiles
resampledAtMast = [''] * noFiles
resampledAtMastDisplay = [''] * noFiles
contour1 = [''] * noFiles
contour1Display = [''] * noFiles
contour2 = [''] * noFiles
contour2Display = [''] * noFiles
contour3 = [''] * noFiles
contour3Display = [''] * noFiles
contour4 = [''] * noFiles
contour4Display = [''] * noFiles
glyph1 = [''] * noFiles
glyph1Display = [''] * noFiles
glyph2 = [''] * noFiles
glyph2Display = [''] * noFiles
resampleWithDataset1 = [''] * noFiles
resampleWithDataset1Display = [''] * noFiles
annotateTimeFilter1Display = [''] * noFiles
resampleWithDataset2 = [''] * noFiles
calculatorConeProj = [''] * noFiles
planeSliceTransform = [''] * noFiles
slice1Transform = [''] * noFiles
slice1TransformDisplay = [''] * noFiles
resampledSliceTransform = [''] * noFiles
annotateTimeFilter1 = [''] * noFiles
cpcsv = [''] * noFiles
cpcsvDisplay = [''] * noFiles
for i_f in range(noFiles):
    colorsCases = cm.jet(range(256))[0::256//noFiles]

    # get the time-keeper
    timeKeeper1 = GetTimeKeeper()
    
    # get active view
    layout1 = GetLayoutByName("Layout #1")
    
    # create a new 'Calculator'
    calculator0[i_f] = Calculator(registrationName='Calculator0', Input=pvdResults[i_f])
    calculator0[i_f].Function = 'u_X*'+str(runwayNormal[0])+'+u_Y*'+str(runwayNormal[1])+'+u_Z*'+str(runwayNormal[2])
    calculator0[i_f].ResultArrayName = 'dot(n,u)'
    
    # create a new 'Calculator'
    calculator1[i_f] = Calculator(registrationName='Calculator1', Input=calculator0[i_f])
    calculator1[i_f].Function = 'sqrt(tk)'
    calculator1[i_f].ResultArrayName = 'sqrtTKE'
    calculatorIPPC[i_f] = Calculator(registrationName='CalculatorIPPC', Input=calculator1[i_f])
    calculatorIPPC[i_f].Function = 'sqrt(tk)'
    calculatorIPPC[i_f].ResultArrayName = 'sqrtTKE_IPPC'
    
    sqrtTKELUT = GetColorTransferFunction('sqrtTKE')
    sqrtTKEPWF = GetOpacityTransferFunction('sqrtTKE')
    sqrtTKE_IPPCLUT = GetColorTransferFunction('sqrtTKE_IPPC')
    
    slice1[i_f] = Slice(registrationName='Slice1', Input=calculatorIPPC[i_f])
    slice1[i_f].SliceType = 'Plane'
    slice1[i_f].SliceType.Normal = runwayNormal 
    slice1[i_f].SliceType.Origin = runwayCenter
    Hide3DWidgets(proxy=slice1[i_f].SliceType)
    
    proj_u[i_f] = Calculator(registrationName='Projected u', Input=calculatorIPPC[i_f])
    proj_u[i_f].Function = 'u_X*iHat+u_Y*jHat+u_Z*kHat-(u_X*'+str(runwayNormal[0])+'+u_Y*'+str(runwayNormal[1])+'+u_Z*'+str(runwayNormal[2])+')*('+str(runwayNormal[0])+'*iHat+'+str(runwayNormal[1])+'*jHat+'+str(runwayNormal[2])+'*kHat)'
    proj_u[i_f].ResultArrayName = 'u_proj'
   
    planeSlice[i_f] = Plane(registrationName='Plane slice')
    planeSliceLength = 43e3
    planeSliceHeight = h_max
    planeSlice[i_f].Origin = np.append(runwayCenter[:2]-runwayDir[:2]*15e3,np.zeros(1))
    planeSlice[i_f].Point1 = planeSlice[i_f].Origin + planeSliceLength*np.append(runwayDir[:2],np.zeros(1))
    planeSlice[i_f].Point2 = planeSlice[i_f].Origin + planeSliceHeight*np.array([0,0,1])
    planeSlice[i_f].XResolution = np.round(planeSliceLength/1500).astype('int')
    planeSlice[i_f].YResolution = np.round(10*planeSliceHeight/1500).astype('int')
    resampledSlice[i_f] = ResampleWithDataset(registrationName='proj_u[i_f] Resampled', SourceDataArrays=proj_u[i_f], DestinationMesh=planeSlice[i_f])
    resampledSlice[i_f].CellLocator = 'Static Cell Locator'
   
   ####################################################################################
   ## Layout 1 - Surface LIC plots
   # create a new 'Clip'
    if plotLIC:
        if i_f == 0:
            topographyDisplay = Show(topography, renderView1, 'StructuredGridRepresentation')
            topographyDisplay.Representation = 'Surface'
            topographyDisplay.Texture = CreateTexture(textureFileName_NIB)
            showRunway(renderView1)
            showMasts(renderView1)
            
            sqrtTKELUTColorBar = GetScalarBar(sqrtTKELUT, renderView1)
            sqrtTKELUTColorBar.WindowLocation = 'UpperRightCorner'
            sqrtTKELUTColorBar.Title = 'Turbulence, $\sqrt{k}$'
            sqrtTKELUTColorBar.ComponentTitle = ''
            sqrtTKELUTColorBar.ScalarBarLength = scalarBarLength
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
        
            sqrtTKELUTColorBar_2 = GetScalarBar(sqrtTKELUT, renderView2)
            sqrtTKELUTColorBar_2.WindowLocation = 'UpperRightCorner'
            sqrtTKELUTColorBar_2.Title = 'Turbulence, $\sqrt{k}$'
            sqrtTKELUTColorBar_2.ComponentTitle = ''
            sqrtTKELUTColorBar_2.ScalarBarLength = scalarBarLength
        
            topographyDisplay = Show(topography, renderView2, 'StructuredGridRepresentation')
            topographyDisplay.Representation = 'Surface'
            topographyDisplay.Texture = CreateTexture(textureFileName_NIB)
        
            # current camera placement for renderView2
            renderView2.OrientationAxesVisibility = 0
            showRunway(renderView2)
            showMasts(renderView2)
     
        # create a new 'Slice'
        slice2[i_f] = Slice(registrationName='Slice2', Input=calculator1[i_f])
        slice2[i_f].SliceType = 'Plane'
        slice2[i_f].HyperTreeGridSlicer = 'Plane'
        slice2[i_f].SliceOffsetValues = [0.0]
        slice2[i_f].SliceType.Normal = [0.0, 0.0, 1.0]
        slice2[i_f].SliceType.Origin = np.append(runwayCenter[:2],horizontalSlice_z)
        Hide3DWidgets(proxy=slice2[i_f].SliceType)
        
        # show data in view
        slice1Display[i_f] = Show(slice1[i_f], renderView1, 'GeometryRepresentation')
        ColorBy(slice1Display[i_f], ('POINTS', 'sqrtTKE'))
        slice1Display[i_f].SetRepresentationType('Surface LIC')
        slice1Display[i_f].ColorMode = 'Multiply'
        slice1Display[i_f].EnhanceContrast = 'Color Only'
        slice1Display[i_f].HighColorContrastEnhancementFactor = 0.3
        slice1Display[i_f].LICIntensity = 0.8
        slice1Display[i_f].MapModeBias = 0.2
        slice1Display[i_f].RescaleTransferFunctionToDataRange(True, False)
        slice1Display[i_f].SetScalarBarVisibility(renderView2, True)# set active source
    
    
        # show data in view
        slice2Display[i_f] = Show(slice2[i_f], renderView2, 'GeometryRepresentation')
        ColorBy(slice2Display[i_f], ('POINTS', 'sqrtTKE'))
        slice2Display[i_f].SetRepresentationType('Surface LIC')
        slice2Display[i_f].ColorMode = 'Multiply'
        slice2Display[i_f].EnhanceContrast = 'Color Only'
        slice2Display[i_f].HighColorContrastEnhancementFactor = 0.3
        slice2Display[i_f].LICIntensity = 0.8
        slice2Display[i_f].MapModeBias = 0.2
        slice2Display[i_f].RescaleTransferFunctionToDataRange(True, False)
        slice2Display[i_f].SetScalarBarVisibility(renderView2, True)# set active source
        
        if i_f == 0:
            insertSINTEFlogo(renderView1,colorLogo)
            insertSINTEFlogo(renderView2,colorLogo)
    ####################################################################################
    ## Layout 3 - Stream lines
    if plotStreamLines:
        if i_f == 0:
            CreateLayout('Layout #3')
            layout3 = GetLayoutByName("Layout #3")
            
            # Create a new 'Render View'
            renderView3 = CreateView('RenderView')
            renderView3.AxesGrid = 'GridAxes3DActor'
            renderView3.StereoType = 'Crystal Eyes'
            renderView3.CameraFocalDisk = 1.0
            renderView3.OrientationAxesVisibility = 0
            AssignViewToLayout(view=renderView3, layout=layout3, hint=0)
            
            topographyDisplay = Show(topography, renderView3, 'StructuredGridRepresentation')
            topographyDisplay.Representation = 'Surface'
            topographyDisplay.Texture = CreateTexture(textureFileName_NIB)
        
            sqrtTKELUTColorBar_3 = GetScalarBar(sqrtTKELUT, renderView3)
            sqrtTKELUTColorBar_3.Title = 'Turbulence, $\sqrt{k}$'
            sqrtTKELUTColorBar_3.Orientation = 'Vertical'
            sqrtTKELUTColorBar_3.WindowLocation = 'UpperRightCorner'
            sqrtTKELUTColorBar_3.ComponentTitle = ''
            sqrtTKELUTColorBar_3.ScalarBarLength = scalarBarLength
            
            showRunway(renderView3)
            showMasts(renderView3)
        
    
        # create a new 'Stream Tracer'
        streamTracer1[i_f] = StreamTracer(registrationName='StreamTracer1', Input=calculator1[i_f], SeedType='Line')
        streamTracer1[i_f].MaximumStreamlineLength = 100000.0
        streamTracer1[i_f].SeedType.Point1 = streamTracerPoint1
        streamTracer1[i_f].SeedType.Point2 = streamTracerPoint2
        streamTracer1[i_f].SeedType.Resolution = 200
        streamTracer1[i_f].MaximumStepLength = 0.5
        streamTracer1[i_f].Vectors = ['POINTS', 'u']
        renderView3.Update()
    
        # show data in view
        streamTracer1Display[i_f] = Show(streamTracer1[i_f], renderView3, 'GeometryRepresentation')
        streamTracer1Display[i_f].Representation = 'Surface'
        streamTracer1Display[i_f].ColorArrayName = ['POINTS', 'sqrtTKE']
        streamTracer1Display[i_f].RescaleTransferFunctionToDataRange(False, True)
        streamTracer1Display[i_f].SetScalarBarVisibility(renderView3, True)
        streamTracer1Display[i_f].SetRepresentationType('Wireframe')
        streamTracer1Display[i_f].RenderLinesAsTubes = 1
        streamTracer1Display[i_f].SetRepresentationType('Surface')
        streamTracer1Display[i_f].LineWidth = 3.0
        Hide3DWidgets(proxy=streamTracer1[i_f].SeedType)
    
        if i_f == 0:
            insertSINTEFlogo(renderView3,colorLogo)
    
   ####################################################################################
   ## Layout 4 - Volume rendering
    if plotVolumeRendering:
        if i_f == 0:
            CreateLayout('Layout #4')
            layout4 = GetLayoutByName("Layout #4")
            
            # Create a new 'Render View'
            renderView4 = CreateView('RenderView')
            renderView4.AxesGrid = 'GridAxes3DActor'
            renderView4.StereoType = 'Crystal Eyes'
            renderView4.CameraFocalDisk = 1.0
            renderView4.OrientationAxesVisibility = 0
            AssignViewToLayout(view=renderView4, layout=layout4, hint=0)
            
            topographyDisplay = Show(topography, renderView4, 'StructuredGridRepresentation')
            topographyDisplay.Representation = 'Surface'
            topographyDisplay.Texture = CreateTexture(textureFileName_NIB)
            
            showRunway(renderView4)
            showMasts(renderView4)
            
            sqrtTKELUTColorBar_4 = GetScalarBar(sqrtTKELUT, renderView4)
            sqrtTKELUTColorBar_4.Title = 'Turbulence, $\sqrt{k}$'
            sqrtTKELUTColorBar_4.Orientation = 'Vertical'
            sqrtTKELUTColorBar_4.WindowLocation = 'UpperRightCorner'
            sqrtTKELUTColorBar_4.ComponentTitle = ''
            sqrtTKELUTColorBar_4.ScalarBarLength = scalarBarLength
        
        # create a new 'Stream Tracer'
        calculator1Display[i_f] = Show(calculator1[i_f], renderView4, 'StructuredGridRepresentation')
        #calculator1Display[i_f].Representation = 'Volume'
        #calculator1Display.ColorArrayName = ['POINTS', 'sqrtTKE']
        ColorBy(calculator1Display[i_f], ('POINTS', 'sqrtTKE'))
        calculator1Display[i_f].SetRepresentationType('Volume')
    
        if timeStepAnnotation and i_f == 0:
            annotateTimeStep(calculator1[i_f],renderView4,location='UpperLeftCorner')
    
        if i_f == 0:
            insertSINTEFlogo(renderView4,colorLogo)
   
    if plotIPPCmapsHorizontal:
        if i_f == 0:
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

            # Create an extended cone at which to plot results
            vtkName = outputPath+'temp'
            coffeeFilter(runwayEndWest,runwayEndEast,h_max,phi1=angleOfAttack_West,phi2=angleOfAttack_East,n=400,name=vtkName)
            extendedCone = LegacyVTKReader(registrationName='ExtendedCone', FileNames=[vtkName+'.vtk'])
        
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

            plane1 = Plane()
            plane1.Origin = [P_sw[0],P_sw[1],z_proj]
            plane1.Point1 = [P_se[0],P_se[1],z_proj]
            plane1.Point2 = [P_nw[0],P_nw[1],z_proj]
            plane1.XResolution = np.round(np.linalg.norm(P_sw[:2]-P_se[:2])/1000).astype('int')
            plane1.YResolution = np.round(np.linalg.norm(P_sw[:2]-P_nw[:2])/1000).astype('int')
            renderView5.InteractionMode = '2D'
            renderView5.CameraPosition = [49425.69655934961, 50403.11453917931, 25824.099005668097]
            renderView5.CameraFocalPoint = [46512.076113695286, 49793.249962534115, 188.00957646880653]
            renderView5.CameraViewUp = [-0.002668502844137075, 0.999720760110162, -0.023479371740545246]
            renderView5.CameraParallelScale = 11840.135737817194
            
        contour1[i_f] = Contour(Input=calculator2)
        contour1[i_f].ContourBy = ['POINTS', 'height']
        contour1[i_f].Isosurfaces = [0.0001]
        contour1[i_f].PointMergeMethod = 'Uniform Binning'
        contour1Display[i_f] = Show(contour1[i_f], renderView5, 'GeometryRepresentation')
        ColorBy(contour1Display[i_f], None)
        contour1Display[i_f].DiffuseColor = [0.552941176470588, 0.149019607843137, 0.129411764705882]
    
        contour2[i_f] = Contour(Input=extendedCone)
        contour2[i_f].ContourBy = ['POINTS', 'z']
        contour2[i_f].Isosurfaces = np.array([500,1000,1500,2000,2500,3000,3500,4000,4500])*0.3048 #convert to meter
        contour2[i_f].PointMergeMethod = 'Uniform Binning'
        contour2Display[i_f] = Show(contour2[i_f], renderView5, 'GeometryRepresentation')
        ColorBy(contour2Display[i_f], None)
        contour2Display[i_f].DiffuseColor = [0.286274509803922, 0.317647058823529, 0.968627450980392]
        
        resampleWithDataset1[i_f] = ResampleWithDataset(SourceDataArrays=calculatorIPPC[i_f], DestinationMesh=extendedCone)
        resampleWithDataset1[i_f].CellLocator = 'Static Cell Locator'
        resampleWithDataset1Display[i_f] = Show(resampleWithDataset1[i_f], renderView5, 'StructuredGridRepresentation')
        resampleWithDataset1Display[i_f].ColorArrayName = ['POINTS', 'sqrtTKE_IPPC']
        resampleWithDataset1Display[i_f].Ambient = 1.0
        resampleWithDataset1Display[i_f].Diffuse = 0.0
    
        contour3[i_f] = Contour(Input=resampleWithDataset1[i_f])
        contour3[i_f].ContourBy = ['POINTS', 'sqrtTKE_IPPC']
        contour3[i_f].Isosurfaces = [2.0,3.0,4.0]
        contour3[i_f].PointMergeMethod = 'Uniform Binning'
        contour3Display[i_f] = Show(contour3[i_f], renderView5, 'GeometryRepresentation')
        try:
            ColorBy(contour3Display[i_f], None)
        except:
            print('No turbulence above 2 to be shown')
        
        contour3Display[i_f].DiffuseColor = [0.0, 0.0, 0.0]
    
        calculatorConeProj[i_f] = Calculator(Input=resampleWithDataset1[i_f])
        calculatorConeProj[i_f].ResultArrayName = 'projCone'
        calculatorConeProj[i_f].Function = 'coordsX*iHat+coordsY*jHat+'+str(z_proj)+'*kHat'
        calculatorConeProj[i_f].CoordinateResults = 1
        resampleWithDataset2[i_f] = ResampleWithDataset(SourceDataArrays=calculatorConeProj[i_f], DestinationMesh=plane1)
        resampleWithDataset2[i_f].CellLocator = 'Static Cell Locator'
        glyph1[i_f] = Glyph(Input=resampleWithDataset2[i_f], GlyphType='Arrow')
        glyph1[i_f].OrientationArray = ['POINTS', 'u']
        glyph1[i_f].ScaleArray = ['POINTS', 'u']
        glyph1[i_f].ScaleFactor = glyphScale1
        glyph1[i_f].GlyphTransform = 'Transform2'
        glyph1[i_f].GlyphMode = 'Every Nth Point'
        glyph1[i_f].Stride = 1
        glyph1Display[i_f] = Show(glyph1[i_f], renderView5, 'GeometryRepresentation')
        glyph1Display[i_f].DiffuseColor = [0.0, 0.0, 0.0]
        ColorBy(glyph1Display[i_f], None)
    
        if i_f == 0:
            insertSINTEFlogo(renderView5,'blue')
        
    if plotIPPCmapsVertical:
        if i_f == 0:
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
        
            sliceTopographyTransform = SItoFTvsNM(sliceTopography,regName='Slice topography transformed')
            sliceTopographyTransformDisplay = Show(sliceTopographyTransform, renderView6, 'StructuredGridRepresentation')
            sliceTopographyTransformDisplay.LineWidth = 3.0
            if False:
                ColorBy(sliceTopographyTransformDisplay, 'height')
            else:
                sliceTopographyTransformDisplay.DiffuseColor = [0.0, 0.0, 0.0]
                ColorBy(sliceTopographyTransformDisplay, None)

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
            planeSliceTransform[i_f] = SItoFTvsNM(planeSlice[i_f],'dummy')
            axes1 = Axes()
            Origin = SItoFTvsNM(np.array(planeSlice[i_f].Origin))
            Point1 = SItoFTvsNM(np.array(planeSlice[i_f].Point1))
            Point2 = SItoFTvsNM(np.array(planeSlice[i_f].Point2))
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

        resampledSliceTransform[i_f] = SItoFTvsNM(resampledSlice[i_f],regName='Slice1 transformed')
        slice1Transform[i_f] = SItoFTvsNM(slice1[i_f],regName='Slice1 transformed')
        slice1TransformDisplay[i_f] = Show(slice1Transform[i_f], renderView6, 'StructuredGridRepresentation')
        slice1TransformDisplay[i_f].ColorArrayName = ['POINTS', 'sqrtTKE_IPPC']
        slice1TransformDisplay[i_f].Ambient = 1.0
        slice1TransformDisplay[i_f].Diffuse = 0.0
    
        glyph2[i_f] = Glyph(registrationName='Glyphs',Input=resampledSliceTransform[i_f], GlyphType='Arrow')
        glyph2[i_f].OrientationArray = ['POINTS', 'u_proj']
        glyph2[i_f].ScaleArray = ['POINTS', 'u_proj']
        glyph2[i_f].ScaleFactor = glyphScale2
        glyph2[i_f].GlyphTransform = 'Transform2'
        glyph2[i_f].GlyphMode = 'Every Nth Point'
        glyph2[i_f].Stride = 1
        glyph2Display[i_f] = Show(glyph2[i_f], renderView6, 'GeometryRepresentation')
        glyph2Display[i_f].DiffuseColor = [0.0, 0.6, 0.0]
        ColorBy(glyph2Display[i_f], None)
        
        contour4[i_f] = Contour(registrationName='Contour IPPC vertical', Input=slice1Transform[i_f])
        contour4[i_f].ContourBy = ['POINTS', 'sqrtTKE_IPPC']
        contour4[i_f].Isosurfaces = [2.0,3.0,4.0]
        contour4[i_f].PointMergeMethod = 'Uniform Binning'
        contour4Display[i_f] = Show(contour4[i_f], renderView6, 'GeometryRepresentation')
        try:
            ColorBy(contour4Display[i_f], None)
        except:
            print('No turbulence above 2 to be shown')
        
        contour4Display[i_f].DiffuseColor = [0.0, 0.0, 0.0]
        
        if i_f == 0:
            insertSINTEFlogo(renderView6,'blue')
        
    if plotOverTime:
        if i_f == 0:
            CreateLayout('Layout #7')
            layout7 = GetLayoutByName("Layout #7")
            quartileChartView1 = CreateView('QuartileChartView')
            quartileChartView1.BottomAxisTitle = 'Iteration number'
            quartileChartView1.LeftAxisTitle = 'TKE $[m^2/s^2]$'
            AssignViewToLayout(view=quartileChartView1, layout=layout7, hint=0)
            colors = cm.jet(range(256))[0::256//noPoints]
        
        cpcsv[i_f] = [''] * noPoints
        cpcsvDisplay[i_f] = [''] * noPoints
        
        for i in range(0,noPoints):
            cpcsv[i_f][i] = CSVReader(registrationName=labels[i], FileName=[outputPath+fileName+'_Point'+str(i+1)+'.csv'])
            cpcsvDisplay[i_f][i] = Show(cpcsv[i_f][i], quartileChartView1, 'XYChartRepresentation')
            cpcsvDisplay[i_f][i].UseIndexForXAxis = 1
            cpcsvDisplay[i_f][i].SeriesColor = ['TKE', str(colors[i][0]), str(colors[i][1]), str(colors[i][2])]
            cpcsvDisplay[i_f][i].SeriesLabel = ['TKE', labels[i]]
            cpcsvDisplay[i_f][i].SeriesVisibility = ['TKE']
        
        saveScreenShot(quartileChartView1,outputPath+fileName+'TKE',saveScreenShots)
        
    if plotVelocityProfiles:
        resampledAtMast[i_f] = [''] * noMasts
        resampledAtMastDisplay[i_f] = [''] * noMasts
        if i_f == 0:
            annotateTimeFilter1Display = [''] * noMasts

        for i in range(0,noMasts):
            resampledAtMast[i_f][i] = ResampleWithDataset(registrationName='calculator1 resampled '+mastLabel[i], SourceDataArrays=calculator1[i_f], DestinationMesh=mastLine[i])
            resampledAtMastDisplay[i_f][i] = Show(resampledAtMast[i_f][i], quartileChartView[i_l][i], 'XYChartRepresentation')
            resampledAtMastDisplay[i_f][i].UseIndexForXAxis = 0
            resampledAtMastDisplay[i_f][i].XArrayName = 'u_Magnitude'
            resampledAtMastDisplay[i_f][i].SeriesLineStyle = ['Points_Z', '1']
            resampledAtMastDisplay[i_f][i].SeriesColor = ['Points_Z', str(colorsCases[i_f][0]), str(colorsCases[i_f][1]), str(colorsCases[i_f][2])]
            resampledAtMastDisplay[i_f][i].SeriesLabel = ['Points_Z', caseNames[i_f]]
            resampledAtMastDisplay[i_f][i].SeriesVisibility = ['Points_Z']
            setVisibility(timestepValues,visibility[i_f],totSteps)
            if i_f == 0:
                annotateTimeFilter1Display[i] = annotateTimeStep(resampledAtMast[i_f][i],quartileChartView[i_l][i],'timeFilter '+mastLabel[i],location='UpperLeftCorner')
        
        
    
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
    renderView1.CameraPosition = [42315.178435332025-originx, 6913455.09641479-originy, 39925.03392197488]
    renderView1.CameraFocalPoint = [47321.88134561954-originx, 6944851.971259325-originy, -590.056714551883]
    renderView1.CameraViewUp = [0.10754399354149799, 0.7793741712263841, 0.6172602293023044]
    renderView1.CameraParallelScale = 34572.810251475086
    renderView1.CameraViewAngle = 30.0
    
    RenderAllViews()
    # get animation scene
    animationScene1 = GetAnimationScene()
    animationScene1.UpdateAnimationUsingDataTimeSteps()
    animationScene1.GoToFirst()
    # update animation scene based on data timesteps
    for i in range(totSteps):
        if not visibility[i_f][i]:
            continue

        fileNameT = caseNames[i_f]+'_'+str(i)+'_'

        if plotLIC:
            saveScreenShot(renderView1,outputPath+fileNameT+'surfaceLICside_bridge'+str(bridge),saveScreenShots)
            
            copyCamera(renderView1,renderView2)
            saveScreenShot(renderView2,outputPath+fileNameT+'surfaceLICtop_bridge'+str(bridge),saveScreenShots)
        
        if plotStreamLines:
            copyCamera(renderView1,renderView3)
            saveScreenShot(renderView3,outputPath+fileNameT+'streamTracer_bridge'+str(bridge),saveScreenShots)
                
        if plotVolumeRendering:
            copyCamera(renderView1,renderView4)
            saveScreenShot(renderView4,outputPath+fileNameT+'volumeRendering_bridge'+str(bridge),saveScreenShots)
            #saveAnimation(renderView4,outputPath+fileNameT+'volumeRendering_bridge'+str(bridge),noSteps,makeVideo)
        
        if plotIPPCmapsHorizontal:
            contour1Display[i_f] = Show(contour1[i_f], renderView5)
            ColorBy(calculator2Display, ('POINTS', 'height'))
            heightLUT.RescaleTransferFunction(0.0, height_max)
            saveScreenShot(renderView5,outputPath+fileNameT+'IPPC_horizontal_bridge'+str(bridge),saveScreenShots)
            
            Hide(contour1[i_f], renderView5)
            ColorBy(calculator2Display, None)
            saveScreenShot(renderView5,outputPath+fileNameT+'IPPC_horizontal_topo4_bridge'+str(bridge),saveScreenShots)
        
        if plotIPPCmapsVertical:
            saveScreenShot(renderView6,outputPath+fileNameT+'IPPC_vertical_bridge'+str(bridge),saveScreenShots)
        
        if plotVelocityProfiles and i_f == noFiles-1:
            for i_l in range(noPlots):
                saveScreenShot(layouts1D[i_l],outputPath+profileNames[i_l]+'_'+str(timestepValues[i]),saveScreenShots,saveAllViews=True)

        if plotError:
            slice1Display[i_f].Representation = 'Surface'
            ColorBy(slice1Display[i_f], ('CELLS', 'Continuous global L2-projection |u^*-u^h|_H1'))
            HideScalarBarIfNotNeeded(sqrtTKELUT, renderView1)
            continuousglobalL2projectionuuh_H1PWF = GetOpacityTransferFunction('ContinuousglobalL2projectionuuh_H1')
            continuousglobalL2projectionuuh_H1LUT = GetColorTransferFunction('ContinuousglobalL2projectionuuh_H1')
            continuousglobalL2projectionuuh_H1LUT.AutomaticRescaleRangeMode = "Never"
            continuousglobalL2projectionuuh_H1LUT.RescaleOnVisibilityChange = 0
            continuousglobalL2projectionuuh_H1LUT.EnableOpacityMapping = 0
            slice1Display[i_f].RescaleTransferFunctionToDataRange(False, True)
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
            renderView1.CameraPosition = [39953.1497599849-originx, 6936202.267042985-originy, 2614.4419848900225]
            renderView1.CameraFocalPoint = [39143.54568842529-originx, 6947093.529147031-originy, 516.8749792426419]
            renderView1.CameraViewUp = [-0.012278300117821047, 0.18822195283001805, 0.9820497644310452]
            renderView1.CameraParallelScale = 28350.540722069076
            saveScreenShot(renderView1,outputPath+fileNameT+'surfaceLICside_bridge'+str(bridge)+'_Error',saveScreenShots)

        animationScene1.GoToNext()

    if plotLIC:
        if not i_f == noFiles-1: 
            Hide(slice1[i_f], renderView1)
            Hide(slice2[i_f], renderView2)
    if plotStreamLines:
        if not i_f == noFiles-1: 
            Hide(streamTracer1[i_f], renderView3)
    if plotVolumeRendering:
        if not i_f == noFiles-1: 
            Hide(calculator1[i_f], renderView4)
    if plotIPPCmapsHorizontal:
        if not i_f == noFiles-1: 
            Hide(resampleWithDataset1[i_f], renderView5)
            Hide(contour1[i_f], renderView5)
            Hide(contour2[i_f], renderView5)
            Hide(contour3[i_f], renderView5)
            Hide(glyph1[i_f], renderView5)
    if plotIPPCmapsVertical:
        if not i_f == noFiles-1: 
            Hide(slice1Transform[i_f], renderView6)
            Hide(contour4[i_f], renderView6)
            Hide(glyph2[i_f], renderView6)
    if plotError:
        if not i_f == noFiles-1: 
            Hide(slice1[i_f], renderView1)
