###################################################
### Output functions
###################################################
def lp_objects_output_ascii(fn, OBJ):

    print('Writing ascii output to: '+fn)

    #fmt = '%7.2f%8.2f%7.1f%7.1f%20.1f   %04d%02d%02d%02d      %7.2f%7.2f%7.1f  %7.2f\n' ;
    fmt = '%7.2f%8.2f%7.1f%7.1f%20.1f   %04d%02d%02d%02d\n'

    year=2011
    month=11
    day=1
    hour=0

    file = open(fn, 'w')

    for ii in range(len(OBJ['lon'])):

        print(fmt % (OBJ['lat'][ii], OBJ['lon'][ii],
                OBJ['y'][ii], OBJ['x'][ii],
                OBJ['area'][ii], year, month, day, hour ) )

        file.write(fmt % (OBJ['lat'][ii], OBJ['lon'][ii],
                OBJ['y'][ii], OBJ['x'][ii],
                OBJ['area'][ii], year, month, day, hour ) ) #, ...
                #stats1(iii).Eccentricity, ...
                #stats1(iii).MajorAxisLength/stats1(iii).MinorAxisLength,...
                #stats1(iii).Orientation,thisVolRain ) ) ;

    file.close()


def lp_objects_output_netcdf():
    pass
