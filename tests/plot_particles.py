import pandas as pd
import numpy as np


def read_pout (filetxt):
    """ Read output text file from ercore fortran version """

    buf=open(filetxt)
    nrel = int(buf.readline().strip())
    dat = None
    x0s = []; y0s = []
    for irel in range(nrel):
        nbuff,x0,y0,rstart,rend = buf.readline().strip().split()
        print 'nbuff,x0,y0,rstart,rend  = ', nbuff,x0,y0,rstart,rend
        x0s += [float(x0)]
        y0s += [float(y0)]
            
    while True:
        line=buf.readline()
        if not line: break
        aux = line.split()
        if len(aux) != 3:
            print 'Error in line:', aux
        else:
            # start reading a release and timestamp
            rel, time,npart = map(float,aux)
            rel = int(rel)
            npart = int(npart)
            print 'Reading release %d at time %.3f with %d particles' % (rel,time,npart)
            #import pdb;pdb.set_trace()
            for i in range(npart):
                aux = buf.readline()
                if not aux: 
                    print 'should not break here. file incomplete' 
                    break
                #print aux
                x,y = map(float, aux.strip().split()[:2])
                rec = np.array([time,rel,x,y])
                dat = rec if dat is None else np.vstack((dat,rec))

    for rel in range(nrel,0,-1):
        rec = np.array([0.,float(rel),x0s[rel-1],y0s[rel-1]])
        dat = np.vstack((rec,dat))

    df = pd.DataFrame(dat, columns=['time', 'release', 'x', 'y']) 
    df.set_index(['time', 'release'], inplace=True)
    return df


def read_release_txt (filein, rel=1):
    """ Read release output file from new ercore """
    from ercore import ncep2dt
    df = pd.read_table(filein)
    time = map(ncep2dt, df['Time'])
    time = map(lambda x: x.strftime('%Y-%m-%d %H:%M:%S'),time)
    release = np.ones(len(time))*rel
    df['time'] = time
    df['release'] = release
    df.set_index(['time', 'release'], inplace=True)
    return df


def plot_frames(df=None,
                filein=None, 
                colors=['b', 'r', 'g', 'k'],
                usemap=True,
                lims=None,
                dmer=2, dpar=2,
                itout=1,
                fileplotpref='plt'):

    import os
    from mpl_toolkits.basemap import Basemap
    import matplotlib.pyplot as plt
    plt.ion()

    if not df:
        if not filein:
            sys.exit('needs pandas df or pickle filein')
        else:
            df = pd.read_pickle(filein)
    
    times = df.index.get_level_values('time').unique()
    releases = df.index.get_level_values('release').unique()

    if not lims:
        lims = [np.floor(df['x'].min()/5)*5, np.ceil(df['x'].max()/5)*5, 
                np.floor(df['y'].min()/5)*5, np.ceil(df['y'].max()/5)*5]

    print 'times    = ', times
    print 'releases = ', releases
    print 'lims     = ', lims

    m = Basemap(projection='merc', llcrnrlon=lims[0],urcrnrlon=lims[1],llcrnrlat=lims[2],urcrnrlat=lims[3], resolution='l')

    fig = plt.figure()

    relx = df.loc[df.index.get_level_values('time') == times[0]]['x'].values
    rely = df.loc[df.index.get_level_values('time') == times[0]]['y'].values
    x0,y0 = m(relx,rely)

    nt = len(times)
    its = np.arange(0,nt,itout)

    for it in its:
        time = times[it]
        dft = df.loc[(df.index.get_level_values('time') == time)]
        fig.clear()
        plt.axis(lims)
        plt.title('time = %s (%d/%d)' % (time, it, nt))
        if usemap:
            m.drawcoastlines()
            m.fillcontinents()
            m.drawmeridians(np.arange(lims[0],lims[1],dmer),labels=[1,1,0,1])
            m.drawparallels(np.arange(lims[2],lims[3],dpar),labels=[1,1,0,1])
            m.plot(x0,y0, 'k+', mew=2)
        else:
            plt.plot(relx,rely, 'k+', mew=2)
        
        for i,rel in enumerate(releases):
            dftr= dft.loc[dft.index.get_level_values('release') == rel]
            if usemap:
                x,y = m(dftr['x'].values, dftr['y'].values)
                m.scatter(x,y, c=colors[i], edgecolor='')
            else:
                plt.scatter(dftr['x'], dftr['y'], c=colors[i], edgecolor='')
        plt.waitforbuttonpress()
        fileplot = fileplotpref +'_it%03i.png' % it
        plt.savefig(fileplot)


def plot_json(filein,lims, old=False,
                colors=['b', 'r', 'g', 'k'],
                usemap=True,
                dmer=2, dpar=2,
                itout=1,
                fileplotpref='plt'):

    import yaml,re, datetime
    from mpl_toolkits.basemap import Basemap
    import matplotlib.pyplot as plt
    plt.ion()

    print 'Reading ', filein
    msg = open(filein).read()
    msg2 = re.sub(r'([a-zA-Z_]+)', r'"\1"', msg) if old else msg
    pout = yaml.load(msg2)

    print 'Processing times'
    times = map(lambda x : datetime.datetime.strptime(x, '%Y/%m/%d %H:%M:%S'), pout.keys())
    times.sort()
    times = map(lambda x: x.strftime('%Y/%m/%d %H:%M:%S'), times)

    # releases = df.index.get_level_values('release').unique()

    print 'times    = ', times
    # print 'releases = ', releases
    print 'lims     = ', lims

    m = Basemap(projection='merc', llcrnrlon=lims[0],urcrnrlon=lims[1],llcrnrlat=lims[2],urcrnrlat=lims[3], resolution='l')

    fig = plt.figure()

    # relx = df.loc[df.index.get_level_values('time') == times[0]]['x'].values
    # rely = df.loc[df.index.get_level_values('time') == times[0]]['y'].values
    # x0,y0 = m(relx,rely)

    nt = len(times)
    its = np.arange(0,nt,itout)


    for it in its:
        time = times[it]
        fig.clear()
        plt.axis(lims)
        plt.title('time = %s' % time)
        if usemap:
            m.drawcoastlines()
            m.fillcontinents()
            m.drawmeridians(np.arange(lims[0],lims[1],dmer),labels=[1,1,0,1])
            m.drawparallels(np.arange(lims[2],lims[3],dpar),labels=[1,1,0,1])

        for poli in pout[time]: 
            print 'poli = ', poli
            lon,lat = poli['x'],poli['y']
            if usemap:
                x,y = m(lon,lat)
                plt.fill(x,y, 'k', alpha=0.3, label=poli['val'])
            else:
                plt.fill(lon,lat, 'k', alpha=0.3, label=poli['val'])
        plt.legend()

        # origin
        for poli in pout[times[0]]: 
            lon,lat = poli['x'],poli['y']
            if usemap:
                x,y = m(lon,lat)
                plt.fill(x,y, 'k', alpha=0.3)
            else:
                plt.fill(lon,lat, 'k', alpha=0.3)
        plt.waitforbuttonpress()
        fileplot = fileplotpref +'_it%03i.png' % it
        plt.savefig(fileplot)


if __name__ == "__main__":

    ## old fortran ercore ##
    # df = read_pout ('../eri_euro_20161002_00z/out/nsea.pout.txt')
    # print df.dtypes
    # df.to_pickle('../eri_euro_20161002_00z/out/nsea.pout.pickle')

    # plot_frames(filein='../eri_euro_20161002_00z/out/nsea.pout.pickle',
    #         usemap = True,
    #         lims = [-5,10,50,60],
    #         dmer = 2,
    #         dpar = 2,
    #         itout = 4,
    #         fileplotpref = 'plt_old_20161002_00z.pout')

    # xxx

    # plot_json(filein='../eri_euro_20161002_00z/out/f3.json',
    #             old=True,
    #             lims = [-5,10,50,60], 
    #             colors=['b', 'r', 'g', 'k'],
    #             usemap=True,
    #             dmer=2, dpar=2,
    #             itout=4,
    #             fileplotpref='plt_f3_old.json')

    # xxx

    ## new ercore ##

    # ts = read_release_txt('rel1.curr_only.txt')
    # print ts.dtypes
    # ts.to_pickle('rel1.curr_only.pickle')

    # plot_frames(filein='rel1.curr_only.pickle',
    #         usemap = True,
    #         lims = [-5,10,50,60],
    #         dmer = 2,
    #         dpar = 2,
    #         itout = 4,
    #         fileplotpref = 'plt_rel1.curr_only')

    # ts = read_release_txt('rel1.curr_wind.txt')
    # print ts.dtypes
    # ts.to_pickle('rel1.curr_wind.pickle')

    # plot_frames(filein='rel1.curr_wind.pickle',
    #         usemap = True,
    #         lims = [-5,10,50,60],
    #         dmer = 2,
    #         dpar = 2,
    #         itout = 4,
    #         fileplotpref = 'plt_rel1.curr_wind')

    # ts = read_release_txt('rel1.curr_wind_diff.txt')
    # print ts.dtypes
    # ts.to_pickle('rel1.curr_wind_diff.pickle')

    # plot_frames(filein='rel1.curr_wind_diff.pickle',
    #         usemap = True,
    #         lims = [-5,10,50,60],
    #         dmer = 2,
    #         dpar = 2,
    #         itout = 4,
    #         fileplotpref = 'plt_rel1.curr_wind_diff')


    # import sys,os,glob
    # filelist = glob.glob(sys.argv[1])
    # fileoutpref = 'rels_curr_wind_diff10_windage'
    # filelist = glob.glob('rel*_curr_wind_diff10_windage.txt')
    # print filelist

    # fileoutpref = 'new_20161002_0000z'
    # filelist = glob.glob('nsea/er20161002_0000z/out/*.txt')
    # print filelist

    # # import pdb;pdb.set_trace()
    # df = None
    # for irel,filein in enumerate(filelist):
    #     ts = read_release_txt(filein,rel=irel+1)
    #     df = ts if df is None else df.append(ts)

    # filepickle = fileoutpref+'.pickle'
    # print 'Saving to', filepickle
    # df.to_pickle(filepickle)

    # print 'Plotting to ', 'plt_'+fileoutpref
    # plot_frames(filein=filepickle,
    #         usemap = True,
    #         lims = [-5,10,50,60],
    #         dmer = 2,
    #         dpar = 2,
    #         itout = 4,
    #         fileplotpref = 'plt_'+fileoutpref)


    # plot_json(filein='nsea/er20161002_0000z/out/f3.json',
    #             old=False,
    #             lims = [-5,10,50,60], 
    #             colors=['b', 'r', 'g', 'k'],
    #             usemap=True,
    #             dmer=2, dpar=2,
    #             itout=4,
    #             fileplotpref='plt_f3.json')


    plot_json(filein='nsea3/er20161002_0000z/out/f3.json',
                old=False,
                lims = [-5,10,50,60], 
                colors=['b', 'r', 'g', 'k'],
                usemap=True,
                dmer=2, dpar=2,
                itout=4,
                fileplotpref='plt_f3_nsea3.json')


