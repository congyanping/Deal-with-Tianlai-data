import os
import numpy as np
import h5py
import healpy as hp
import aipy as a
# import iuwts
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import signal
import math
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import FK5
fl = './map_full.hdf5'
with h5py.File(fl, 'r') as f:
    c = f['map'][0, 0, :]
    nside = hp.npix2nside(c.shape[0])
    #print 'nside',nside
    #print 'c.shape[0]',c.shape[0]
    #print 'c',c
    #flux = 10 # Jy
    flux = 5 # Jy
    #flux = 2.5 # Jy
    frequency = 750 # MHz
    catalog = 'nvss'
    # catalog = 'wenss'
    dat = 'ra_dec_jy_%.1f_%.1f_%s.npz' % (flux, frequency, catalog)
    # load data from file if exist
    if os.path.exists(dat):
        rdj = np.load(dat)
        ras = rdj['ras']
        decs = rdj['decs']
        jys = rdj['jys']
    else:
        src = '%f/%f' % (flux, frequency / 1.0e3)
        srclist, cutoff, catalogs = a.scripting.parse_srcs(src, catalog)
        cat = a.src.get_catalog(srclist, cutoff, catalogs)
        nsrc = len(cat) # number of sources in cat
        ras = [ np.degrees(cat.values()[i]._ra) for i in range(nsrc) ]
        decs = [ np.degrees(cat.values()[i]._dec) for i in range(nsrc) ]
        jys = [ cat.values()[i].get_jys() for i in range(nsrc) ]
        

        print 'ras',ras
        print 'decs',decs
        print 'jys',jys
        # select sources
        inds = np.where(np.array(decs)>-15.0)[0]
        ras = np.array(ras)[inds]
        decs = np.array(decs)[inds]
        jys = np.array(jys)[inds]

        # save to file
        np.savez(dat, ras=ras, decs=decs, jys=jys)

#above select source upper 15 degree containing in the file.
    #print 'zip(ras,decs,jys)',zip(ras,decs,jys.tolist()),len(zip(ras,decs,jys.tolist()))
    pix = []
    pixx = []
    x = []
    y = []
    z = []
    ras_1 =[]
    decs_1 = []
    coordinate_ra = []
    ra_dec = []
    galactic_la = []
    for ra, dec, jy in zip(ras, decs, jys):
        theta, phi = np.radians(90.0-dec), np.radians(ra)
        decla = np.radians(dec)
        angle = np.radians(44) - decla
        angle = np.abs(angle)
        vec = hp.pixelfunc.ang2vec(theta, phi)
        #print 'vec',vec
        
        pix_number = hp.vec2pix(nside,vec[0],vec[1],vec[2])
        
        #radius = np.radians(5.0)
        
        radius = min(np.radians(5.0), np.radians(3.0 + (jy - 3.0)/100.0))
        pix = hp.query_disc(nside,vec,radius,inclusive=False)










        cp = np.ones_like(c,dtype = np.float)
        #cp[pix] = np.abs(c[pix])#cong
        sigma = 7933.36414712
        cpp = c.tolist()
        for num,val in enumerate(cpp):
            if val > sigma:
                cp[num]=c[num]
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        """
        #cong
        x0,y0,z0 = hp.pix2vec(nside,pix_number)
        x1,y1,z1 = hp.pix2vec(nside,pix)
        #print 'pix_numberi_vector',type(pix_number_vector),pix_number_vector,len(pix_number_vector)
        distance = (x1-x0)**2 + (y1-y0)**2 + (z1-z0)**2
        print 'distance',distance,distance.shape
        print x1.shape,y1.shape,z1.shape
        print x0,y0,z0
        print 'nside',nside
        n = [ni for ni in range(x1.shape[0])]
        index = []
        for value in sorted(zip(distance.tolist(),n))[0:66]:
            index.append(value[1])
        index = np.array(index)
        print 'index',index,len(index)
        x1 = x1[index]
        y1 = y1[index]
        z1 = z1[index]
        better_pix = hp.vec2pix(nside,x1,y1,z1)
        print 'better_pix',better_pix
        cp[better_pix] = c[better_pix]
        tem = max(cp[better_pix])

        """
        #N=78 #for radius =5.0
        #N=46 #for radius = 3.0
        #N = 24 #for radius = 1.5
        #N = 16 #for radius =1.0
        cp = cp[pix].tolist()
        N = len(cp)
        dimention = round(np.sqrt(N)) + 1
        dimention = int(dimention)
        #print 'len(cp)',len(cp)
        cp =np.array(cp)
        

         
        while not cp.shape[0] == dimention**2:
            if cp.shape[0] == dimention**2:
                #print 'fullfill the condition' 
                #cp = cp.reshape(N,N)
                pass 
            elif cp.shape[0] > dimention**2:
                #print 'Not fullfill the condition'
                cp = cp.tolist()
                cp.pop()
                cp = np.array(cp)
            else:
                cp = cp.tolist()
                cp.append(0)
                cp = np.array(cp)
                #print 'Not fullfill the condition'
        #print 'cp.shape',cp.shape
        

        """
        #Gaussian filtering

        import scipy.optimize as opt 
        import math
        def FittingFunction((x, y), b, a, x0, y0 ):
            g=np.exp(-b*( pow(x,2) + pow(y,2) ))*np.cos(a*pow(x-x0, 2)-a*pow(y-y0, 2) )
            return g.ravel()
        def twoD_Gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):

            xo = float(xo)
            yo = float(yo)
            #sigma_y = sigma_x
            a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
            b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
            c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
            g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo)+ c*((y-yo)**2)))
            return g.ravel()
        def twoD_Gaussian_2((x, y), amplitude, xo, yo, sigma_x, sigma_y,theta, offset):

            xo = float(xo)
            yo = float(yo)    
            #sigma_y = sigma_x
            a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
            b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
            c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
            g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) + c*((y-yo)**2)))
            return g
        def twoD_2Gaussian((x,y),a0,x0,y0,sx0,sy0,theta0,offset0,a1,x1,y1,sx1,sy1,theta1,offset1):
            
            xo = float(x0)
            yo = float(y0)    
            #sy0 = sx0
            a0 = (np.cos(theta0)**2)/(2*sx0**2) + (np.sin(theta0)**2)/(2*sy0**2)
            b0 = -(np.sin(2*theta0))/(4*sx0**2) + (np.sin(2*theta0))/(4*sy0**2)
            c0 = (np.sin(theta0)**2)/(2*sx0**2) + (np.cos(theta0)**2)/(2*sy0**2)
            g0 = offset0 + a0*np.exp( - (a0*((x-xo)**2) + 2*b0*(x-xo)*(y-yo)+ c0*((y-yo)**2)))

            x1 = float(x1)
            y1 = float(y1)    
            #sy1 = sx1
            a1 = (np.cos(theta1)**2)/(2*sx1**2) + (np.sin(theta1)**2)/(2*sy1**2)
            b1 = -(np.sin(2*theta1))/(4*sx1**2) + (np.sin(2*theta1))/(4*sy1**2)
            c1 = (np.sin(theta1)**2)/(2*sx1**2) + (np.cos(theta1)**2)/(2*sy1**2)
            g1 = offset1 + a1*np.exp( - (a1*((x-x1)**2) + 2*b1*(x-x1)*(y-y1)+ c1*((y-y1)**2)))
            return (g0+g1).ravel()
        def twoD_2Gaussian_2((x,y),a0,x0,y0,sx0,sy0,theta0,offset0,a1,x1,y1,sx1,sy1,theta1,offset1):
            
            xo = float(x0)
            yo = float(y0)
            #sy0=sx0
            a0 = (np.cos(theta0)**2)/(2*sx0**2) + (np.sin(theta0)**2)/(2*sy0**2)
            b0 = -(np.sin(2*theta0))/(4*sx0**2) + (np.sin(2*theta0))/(4*sy0**2)
            c0 = (np.sin(theta0)**2)/(2*sx0**2) + (np.cos(theta0)**2)/(2*sy0**2)
            g0 = offset0 + a0*np.exp( - (a0*((x-xo)**2) + 2*b0*(x-xo)*(y-yo)+ c0*((y-yo)**2)))

            x1 = float(x1)
            y1 = float(y1)    
            #sy1 = sx1
            a1 = (np.cos(theta1)**2)/(2*sx1**2) + (np.sin(theta1)**2)/(2*sy1**2)
            b1 = -(np.sin(2*theta1))/(4*sx1**2) + (np.sin(2*theta1))/(4*sy1**2)
            c1 = (np.sin(theta1)**2)/(2*sx1**2) + (np.cos(theta1)**2)/(2*sy1**2)
            g1 = offset1 + a1*np.exp( - (a1*((x-x1)**2) + 2*b1*(x-x1)*(y-y1)+ c1*((y-y1)**2)))
            return g0+g1
        def round_gaussian((x,y),a,x0,y0,sigma):
            g = a*np.exp((-(x-x0)**2-(y-y0)**2)/2*sigma**2)
            return g.ravel()
        
        def round_gaussian_2((x,y),a,x0,y0,sigma):
            g = a*np.exp((-(x-x0)**2-(y-y0)**2)/2*sigma**2)
            return g
        def round_2gaussian((x,y),a0,x0,y0,sigma_x0,a1,x1,y1,sigma_x1):
            g0 = a0*np.exp((-(x-x0)**2)/sigma_x0**2+(-(y-y0)**2)/sigma_x0**2)
            g1 = a1*np.exp((-(x-x1)**2)/sigma_x1**2+(-(y-y1)**2)/sigma_x1**2)
            g = g0+g1
            return g.ravel()

        def round_2gaussian_2((x,y),a0,x0,y0,sigma_x0,a1,x1,y1,sigma_x1):
            g0 = a0*np.exp((-(x-x0)**2)/sigma_x0**2+(-(y-y0)**2)/sigma_x0**2)
            g1 = a1*np.exp((-(x-x1)**2)/sigma_x1**2+(-(y-y1)**2)/sigma_x1**2)
            g = g0+g1
            return g
        x_1=np.linspace(0,dimention-1, dimention)
        y_1=np.linspace(0, dimention-1, dimention)
        x_1,y_1 = np.meshgrid(x_1,y_1)
        #initial_gaussian =(3.65688487e+03,5.01102366e+02,-1.93017373e+02,3.87846837e+04,3.53361452e+00,1.36910262e+03)
        ini_gau = (1.00288930e+05,5.01020455e+03,2.49410677e+03,1.64272669e+01,2.81258581e+04,2.14745072e+02,-4.02323482e+02)

        try:

            popt, pcov = opt.curve_fit(twoD_Gaussian, (x_1, y_1), cp, bounds=(0, dimention),method = 'trf')
            #popt, pcov = opt.curve_fit(round_gaussian, (x_1, y_1), cp)
            #print 'popt',popt
            #print 'pcov',pcov
            tem = twoD_Gaussian_2((popt[1],popt[2]),*popt)
            if np.log10(tem/jy)<0.1:
                tem = np.nan
                Jy = np.nan
                print '1_Gaussian_log(tem/jy):',np.log10(tem/jy),'(Tianlai,NVSS):',(tem,jy),'(ra,dec):',(ra,dec)
                ras_1.append(ra)
                decs_1.append(dec)
            #print 'tem',tem
            result = twoD_Gaussian((x_1,y_1),*popt)
            #print 'dimention',dimention
            fig,ax =plt.subplots(1)
            ax.imshow(result.reshape(dimention,dimention),extent=(x_1.min(),x_1.max(),y_1.min(),y_1.max()))
            titile = '1_%s_%s_%s.png'%(ra,dec,dimention)
            plt.savefig('./map_gaussian/'+titile)
            plt.close()
        except RuntimeError:
            #ini_gau = (3.90000000e+01,4.34266667e+01,5.33302406e+01,1.53363409e-02,1.17762486e+01,6.89884036e+01,3.57737789e+01,3.83361060e+01,4.38255376e+01,5.67417730e+01,6.81352527e-03,1.24354117e+01,7.57593673e+01,3.57725875e+01) 
            popt ,pcov = opt.curve_fit(twoD_2Gaussian,(x_1,y_1),cp,bounds = (0,dimention),method = 'dogbox')
            #print '2popt',popt
            tem = max(twoD_2Gaussian_2((popt[1],popt[2]),*popt),twoD_2Gaussian_2((popt[5],popt[6]),*popt))
            if RuntimeWarning:
                print '2Gaussian_tem',tem
                print '2Gaussian_Jy',jy
                print '(ra,dec)',(ra,dec)
            result = twoD_2Gaussian((x_1,y_1),*popt)
            fig,ax =plt.subplots(1)
            ax.imshow(result.reshape(dimention,dimention),extent=(x_1.min(),x_1.max(),y_1.min(),y_1.max()))
            titile = '2_%s_%s.png'%(ra,dec)
            plt.savefig('./map_gaussian/'+titile)
            plt.close()
        """     


        
        """ 
        #gaussian convolution
        mean,sigma = 1.0, 3.5
        mean,sigma = 0.00001,0.000001
        s = np.random.normal(mean,sigma,(N,N))
        tem = signal.convolve2d(cp,s,'valid')
        print 'tem.shape',tem.shape
        """
        
        pixx += hp.query_disc(nside, vec, radius, inclusive=False).tolist()
        #tem = np.median(cp[pix])
        #tem = np.sum(cp)

        """ 
        tem = np.abs(tem)
        """

        #use mr xuelei chen opinion
        tem = np.sum(cp)
        #because there are some pixel's value is zero
        if np.log10(tem/jy) < 0.1:
            tem = np.nan
            jy = np.nan




        
        #tem = np.median(np.array(cp)) 
        #tem = np.sum(tem)/(tem.shape[0]*tem.shape[1])
        #print (jy,tem)
        #x.append(jy.tolist()[0]*(np.abs(np.cos(phi))*np.sin(theta))**2.5*0.1+1)
        x.append(jy.tolist()[0])
        #x.append(jy.tolist()[0]*np.abs(np.sin(theta+phi)))
        #y.append(np.abs(tem)*1/np.cos(angle))
        y.append(tem)
        z.append(angle)
        coordinate_ra.append(np.radians(ra))
        ra_dec.append([ra,dec])
        sc = SkyCoord(ra = ra,dec = dec,unit = 'deg',frame = FK5,equinox='J1980.0')
        gl = sc.transform_to(frame = 'galactic').b
        galactic_la.append(np.radians(gl.value))
    
    """
    i = [i for i in xrange(len(y))]
    a = zip(y,x,i,ra_dec)
    print 'len(a)',len(a)
    a = sorted(a)
    #a = np.array(a)
    print 'len(y)',len(y)
    
    for i,j in enumerate(a):
        x[i]=j[1]
        y[i]=j[0]
    
    
    ras_1 =[]
    decs_1 = []
    for value in a:
        for threshold in [170,29,162,165,2,138,187,180,27,144,133,88,157,102,202,128,80,131,18,129,174,116,113,10,34]:
            if value[2]==threshold:
                ras_1.append(value[3][0])
                decs_1.append(value[3][1])
                print (value[0],value[1],value[2],value[3])
    ras = np.array(ras_1)
    decs = np.array(decs_1)
    #x=x[1:-1]
    xx=[]
    yy=[]
    zz=[]
    #for i in [11,72,49,50,39,2,138,162,187,180,142,144,133,143,157,102,202,128,80,131,18,129,174]:
    for i in [170,29,162,165,2,138,187,180,27,144,133,88,157,102,202,128,80,131,18,129,174,116,113,10,34]:
        xx.append(x[i])
        yy.append(y[i])
        zz.append(z[i])
    """
    xx = x
    yy = y
    zz = z

    
    import matplotlib.pyplot as plt 
    plt.figure(1)
    #result = np.log(x)/np.log(y)
    #for i in range(len(y)):
        #if y[i]>13:
            #y[i]=y[i]-4*2
            #pass
        #else:
            #y[i]=y[i]+4*2
            #pass
    #a = zip(yy,xx)
    #a = sorted(a)
    #for i,j in enumerate(a):
    #    xx[i]=j[1]
    #    yy[i]=j[0]
    y = np.array(yy)
    x = np.array(xx)
    z = np.array(zz)
    c_d = np.array(coordinate_ra)
    g_l = np.array(galactic_la)
    #x = x[y<500]
    #y = y[y<500]
    """
    from mpl_toolkits.mplot3d import axes3d
    import matplotlib.pyplot as plt
    from matplotlib import cm
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #X, Y, Z = axes3d.get_test_data(0.05)
    ax.scatter(x, y, z)
    #cset = ax.contour(x, y, z, cmap=cm.coolwarm)
    #ax.clabel(cset, fontsize=9, inline=1)
    """
    idx = np.isfinite(x) & np.isfinite(y) 
    plt.plot(z,np.log10(y/x),color = 'r',linestyle='',marker = 'p',MarkerSize=1.0,label = 'point source')
    #plt.plot(g_l,np.log10(y/x),color = 'r',linestyle='',marker = 'p',MarkerSize= 2.5,label = 'point source')
    #plt.plot(c_d,np.log10(y/x),color = 'g',linestyle='',marker = '*',label = 'point source')
    #plt.plot(c_d,np.log10(y/x),color = 'g',linestyle='',marker = '*',label = 'point source')
    #plt.plot(result,'ro')
    #optimize stage
    #z1=np.polyfit(z[idx],np.abs(np.log10(y[idx]/x[idx])),1)
    z1=np.polyfit(y[idx],np.abs(np.log10(y[idx]/x[idx])),3)
    
    yvals = np.polyval(z1,x)
    print 'yvals',yvals
    #yvals = np.polyval(z1,c_d)
    
    #plt.plot(x,yvals,color ='b',linestyle = '-',marker='',label='polyfit curve')
    #plt.plot(x,yvals+0.4,color = 'g',linestyle = '-.',marker = '',label = 'upper limit')
    #plt.plot(x,yvals-0.6,color = 'g',linestyle = '-.',marker = '',label = 'lower limit')
    
    #plt.plot(c_d,yvals,color ='b',linestyle = '-',marker='',label='polyfit curve')
    #plt.plot(c_d,yvals+0.4,color = 'g',linestyle = '-.',marker = '',label = 'upper limit')
    #plt.plot(c_d,yvals-0.6,color = 'g',linestyle = '-.',marker = '',label = 'lower limit')
    import time
    mean =''
    sigma = ''
    plt.legend(loc='upper left')
    #title = '(flux,fre)=%s'%((flux,frequency))
    title = 'point source'
    #plt.ylim(0,4)
    #plt.xlim(0,80)

    plt.title(title)
    plt.ylabel('log10(Tianlai/NVSS)')
    #plt.xlabel('Zenith angle')
    #plt.xlabel('NVSS')
    #plt.xlabel('Tianlai')
    plt.xlabel('Galactic latitude')
    #plt.xlabel('right ascension')
    title = '%s.png'%(time.ctime())
    plt.savefig('./map/'+title)
    plt.close()
    

    pixx = list(set(pixx)) # unique pixels
    from pylab import cm
    fig = plt.figure(1)
    cmap = cm.Blues_r
    #cmap = cm.hot
    cmap.set_under('w')
    
    """
    fig = plt.figure(1)
    hp.mollview(np.abs(c), fig=1, title='c_map', cmap=cmap, min=0, max=2000)
    hp.graticule(verbose=False)
    fig.savefig('c_map.png')
    plt.close()
    
    """
    
    cp = np.zeros_like(c)
    cp[pixx] = c[pixx]

    ras = np.array(ras_1)
    decs = np.array(decs_1)
    hp.mollview(np.abs(cp), fig=1, title='c_ps',min = 0,max = 2000)
    hp.projscatter(ras, decs, lonlat=True, s=100, facecolors='none', edgecolors='r', alpha=1.0, linewidth=1.0)
    hp.graticule(verbose=False)
    fig.savefig('c_ps.png')
    plt.close()
    
    """  
    c[pix] = 0
    fig = plt.figure(1)
    hp.mollview(np.abs(c), fig=1, title='c_hole', cmap=cmap, min=0, max=2000)
    hp.graticule(verbose=False)
    fig.savefig('c_hole.png')
    plt.close()
    """

    """ 
    # # median filtering
    c = np.abs(c)
    # c_median = np.zeros_like(c)
    c_median = c
    radius = np.radians(30)
    for i in range(c.shape[0]):
        if i in c.shape[0] / 10 * np.arange(10):
            print '%d of %d...' % (i, c.shape[0])
        # if c[i] < 400:
        #     c_median[i] = c[i]
        #     continue
        colat, lon = hp.pix2ang(nside, i)
        # if np.pi/2 < np.abs(lon-np.pi) and np.abs(lon-np.pi) < np.pi and colat > np.pi/6 and c[i] > 400:
        if np.pi/2 < np.abs(lon-np.pi) and np.abs(lon-np.pi) < np.pi and colat > np.pi/6 and c[i] > (1200.0 - 1800.0*colat/np.pi):
            vec = hp.pix2vec(nside, i)
            pix = hp.query_disc(nside, vec, radius, inclusive=False)
            vals = (c[pix])[c[pix] > 0]
            #the value of pix index in c and above zero
            # vals = c[pix]
            med = np.median(vals)
            c_median[i] = med
            # abs_diff = np.abs(vals - med)
            # mad = np.median(abs_diff) / 0.6745
            # if np.abs(c[i] - med) < 2.5 * mad:
            #     c_median[i] = c[i]
            # else:
            #     c_median[i] = med


    fig = plt.figure(1)
    hp.mollview(np.abs(c_median), fig=1, title='c_median', cmap=cmap, min=0, max=2000)
    hp.graticule(verbose=False)
    fig.savefig('c_median.png')
    plt.close()

    fig = plt.figure(1)
    hp.mollview(np.abs(cp) + np.abs(c_median), fig=1, title='c_new', cmap=cmap, min=0, max=2000)
    hp.graticule(verbose=False)
    fig.savefig('c_new.png')
    plt.close()

    fig = plt.figure(1)
    hp.mollview(np.abs(cp) + np.abs(c_median), fig=1, title='c_new_nvss', cmap=cmap, min=0, max=2000)
    hp.projscatter(ras, decs, lonlat=True, s=150, facecolors='none', edgecolors='w', alpha=1.0, linewidth=1.0)
    hp.graticule(verbose=False)
    fig.savefig('c_new_nvss.png')
    plt.close()
    """
