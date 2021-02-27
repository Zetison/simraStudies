from pyproj import Transformer
import numpy as np
proj32 = "+proj=utm +zone=32K, +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
proj33 = "+proj=utm +zone=33K, +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
transformer = Transformer.from_crs(proj32,proj33)

height = 70 # Height of bridge at middle point
Boom1 = np.array([6.1,6.1,4.4,3.6])
Sensorh = np.array([[92.5,71.5,44.5],[76.8, 48.3, 27.3],[94.8, 75.0, 50.0, 27.0],[62.8, 40.0, 13.4]])
BoomOrient = np.array([[72,74,74],[289, 290, 290],[81, 81, 81, 81],[223, 223, 223]])
CoordUTM32 = np.array([[6924741, 345142],[6925267, 348347],[6920740, 346520],[6922074, 351140]])
f = open("POINT_FILE", "w")
noPoints = sum([len(listElem) for listElem in Sensorh])+len(Sensorh)//2
f.write("%d\n" % noPoints)
for i in range(0,Sensorh.shape[0]):
    for j in range(0,len(Sensorh[i])):
        alpha = np.radians(BoomOrient[i][j])
        x, y = transformer.transform(CoordUTM32[i,1]+Boom1[i]*np.sin(alpha),CoordUTM32[i,0]+Boom1[i]*np.cos(alpha)) # This is probably not the exact formula due to the UTM transformation?
        f.write("%f %f %f\n" % (x,y,Sensorh[i][j]))

for i in range(0,CoordUTM32.shape[0]//2):
    x1, y1 = transformer.transform(CoordUTM32[2*i,1],CoordUTM32[2*i,0])
    x2, y2 = transformer.transform(CoordUTM32[2*i+1,1],CoordUTM32[2*i+1,0])
    f.write("%f %f %f\n" % ((x1+x2)/2,(y1+y2)/2,height))

f.close()
