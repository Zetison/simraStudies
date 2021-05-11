from pyproj import Transformer
import pandas as pd
import click
import numpy as np
@click.command()
@click.option('--originx', default=0.0, type=float, help='x-coordinate of origin')
@click.option('--originy', default=0.0, type=float, help='y-coordinate of origin')
@click.option('--date', default='2020-11-19 06:00', type=str, help='Date for data extraction')
def main(originx,originy,date):
    proj32 = "+proj=utm +zone=32K, +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
    proj33 = "+proj=utm +zone=33K, +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
    transformer = Transformer.from_crs(proj32,proj33)

    height = 70 # Height of bridge at middle point
    # Data for 'SulaNW (Kvitneset)', 'SulaNE (Trælboneset)', 'SulaSW (Langeneset)', 'SulaSE (Kårsteinen)'
    Boom1 = np.array([6.1,6.1,4.4,3.6])
    Sensorh = np.array([[92.5,71.5,44.5],[76.8, 48.3, 27.3],[94.8, 75.0, 50.0, 27.0],[62.8, 40.0, 13.4]],dtype=object)
    BoomOrient = np.array([[72,74,74],[289, 290, 290],[81, 81, 81, 81],[223, 223, 223]],dtype=object)
    CoordUTM32 = np.array([[6924741, 345142],[6925267, 348347],[6920740, 346520],[6922074, 351140]])
    CoordUTM32 = np.array([[6924741, 345142],[6925263.08,348346.72],[6920736.2,346525.95],[6922076.51,351139.85]])
    masth = np.array([100.5, 78.0,97.0,63.0])
    #mastb = np.array([6.0,6.0,14.0,12.0])
    mastb = np.array([6.0,14.0,6.0,12.0])
    mastb = np.array([8.0,12.2,1.8,13.9])
    f = open("POINT_FILE", "w")
    noPoints = sum([len(listElem) for listElem in Sensorh])+len(Sensorh)//2
    f.write("%d\n" % noPoints)
    for i in range(Sensorh.shape[0]):
        for j in range(len(Sensorh[i])):
            alpha = np.radians(BoomOrient[i][j])
            x, y = transformer.transform(CoordUTM32[i,1]+Boom1[i]*np.sin(alpha),CoordUTM32[i,0]+Boom1[i]*np.cos(alpha)) # This is probably not the exact formula due to the UTM transformation?
            f.write("%f %f %f\n" % (x-originx,y-originy,mastb[i]+Sensorh[i][j]))

    for i in range(0,CoordUTM32.shape[0]//2):
        x1, y1 = transformer.transform(CoordUTM32[2*i,1],CoordUTM32[2*i,0])
        x2, y2 = transformer.transform(CoordUTM32[2*i+1,1],CoordUTM32[2*i+1,0])
        f.write("%f %f %f\n" % ((x1+x2)/2-originx,(y1+y2)/2-originy,height))

    f.close()

    f = open("mastLoc.csv", "w")
    for i in range(len(masth)):
        x, y = transformer.transform(CoordUTM32[i,1],CoordUTM32[i,0])
        f.write("%f,%f,%f,%f\n" % (x-originx,y-originy,mastb[i],masth[i]))

    f.close()
    
    mastNames = ['Kvitneset', 'Traelboneset','Langeneset','Kaarsteinen']
    for i in range(len(masth)):
        df = pd.DataFrame()
        for j in range(len(Sensorh[i])):
            z = np.floor(Sensorh[i][j]).astype(int)
            filename = 'measurements/'+'10hz_'+mastNames[i]+'_60mnturbulence_statistics_'+str(z)+'_202011.csv'
            df_all = pd.read_csv(filename)
            df = df.append(df_all[df_all.date==date])
        df = df.reset_index()
        df['coordsZ'] = mastb[i]+Sensorh[i]
        csvName = 'VelocityProfile_'+mastNames[i]+'_'+str(round(float(date[12])-6))+'.csv'
        print(csvName)
        print(date)
        df.to_csv(csvName)




if __name__ == '__main__':
    main()
