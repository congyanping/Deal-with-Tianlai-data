"""RFI statistics.

Inheritance diagram
-------------------

.. inheritance-diagram:: Stats
   :parts: 2

"""

from datetime import datetime
import numpy as np
import timestream_task
from tlpipe.utils.path_util import output_path
import tlpipe.plot
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import MaxNLocator, AutoMinorLocator
from caput import mpiutil
from tlpipe.container.raw_timestream import RawTimestream
from tlpipe.container.timestream import Timestream
import h5py
class Stats(timestream_task.TimestreamTask):
    """RFI statistics.

    Analysis of RFI distributions along time and frequency.

    """

    params_init = {
                    'excl_auto': False, # exclude auto-correclation
                    'plot_stats': True, # plot RFI statistics
                    'fig_name': 'stats/stats',
                    'rotate_xdate': False, # True to rotate xaxis date ticks, else half the number of date ticks
                  }

    prefix = 'rs_'

    def process(self, ts):

        excl_auto = self.params['excl_auto']
        plot_stats = self.params['plot_stats']
        fig_prefix = self.params['fig_name']
        rotate_xdate = self.params['rotate_xdate']
        tag_output_iter = self.params['tag_output_iter']

        ts.redistribute('baseline')
        
        if isinstance(ts, RawTimestream):
            func = ts.bl_data_operate
        elif isinstance(ts,Timestream):
            func = ts.pol_and_bl_data_operate
        show_progress = self.params['show_progress']
        progress_step = self.params['progress_step']
        
        func(self.realvsimag,full_data=True,show_progress = show_progress,progress_step=progress_step,keep_dist_axis=False)

        if ts.local_vis_mask.ndim == 3: # RawTimestream
            if excl_auto:
                bl = ts.local_bl
                vis_mask = ts.local_vis_mask[:, :, bl[:, 0] != bl[:, 1]].copy()
                vis = ts.local_vis[:, :, bl[:, 0] != bl[:, 1]].copy()
                rvi_vis_mask = ts.local_vis_mask[:,:,np.abs(bl[:,0]-bl[:,1])==4]
                rvi_vis = ts.local_vis[:,:,np.abs(bl[:,0]-bl[:,1]) == 4]
            else:
                vis_mask = ts.local_vis_mask.copy()
                vis = ts.local_vis.copy()
                #vis_copy = vis[vis_mask].copy()
            nt, nf, lnb = vis_mask.shape

        #princinple operate for timestream data by congyanping
        elif ts.local_vis_mask.ndim == 4: # Timestream
            # suppose masks are the same for all 4 pols
            if excl_auto:
                bl = ts.local_bl.copy()
		feed_pos = ts['feedpos'][...].copy()
                if True:
                    #bl[:,0][bl[:,1]>31] = 100
                    #bl[:,1][bl[:,1]>31] = 100
                    vis_mask = ts.local_vis_mask[:, :, 0, np.abs(bl[:, 0] - bl[:, 1])==6*3].copy()
                    #bl_label = np.vstack((bl[:,0][np.abs(bl[:,0]-bl[:,1])==6],bl[:,1][np.abs(bl[:,0]-bl[:,1])==6]))
                    bl_label = bl[np.abs(bl[:,0] - bl[:,1]) == 6*3 , :]
                    print 'bl_label.shape',bl_label.shape
                    vis = ts.local_vis[:,:,0,np.abs(bl[:,0] - bl[:,1])==6*3].copy()

                    rvi_vis_mask = ts.local_vis_mask[:,:,0,np.abs(bl[:,0]-bl[:,1])==6*3].copy()
                    #rvi_vis_mask = ts.local_vis_mask[:,:,0,bl[:,0] == bl[:,1]]
                    print "@"*20
                    print "if np.abs(bl[:,0]-bl[:,1] = 2 ",rvi_vis_mask.shape
                    print "@"*20
                    rvi_vis = ts.local_vis[:,:,0,np.abs(bl[:,0] - bl[:,1])==6*3].copy()
                    print 'rvi_vis.shape',rvi_vis.shape
            else:
                vis_mask = ts.local_vis_mask[:, :, 0].copy()
                vis = ts.local_vis[:,:,0].copy()
                #vis_copy = vis[vis_mask].copy()
            nt, nf, lnb = vis_mask.shape
        else:
            raise RuntimeError('Incorrect vis_mask shape %s' % ts.local_vis_mask.shape)

        # total number of bl
        nb = mpiutil.allreduce(lnb, comm=ts.comm)

        # un-mask ns-on positions
        if 'ns_on' in ts.iterkeys():
            vis_mask[ts['ns_on'][:]] = False
            rvi_vis_mask[ts['ns_on'][:]] = False




        # statistics along time axis
        time_mask = np.sum(vis_mask, axis=(1, 2)).reshape(-1, 1)
        # gather local array to rank0
        time_mask = mpiutil.gather_array(time_mask, axis=1, root=0, comm=ts.comm)
        if mpiutil.rank0:
            time_mask = np.sum(time_mask, axis=1)

        # statistics along time axis
        freq_mask = np.sum(vis_mask, axis=(0, 2)).reshape(-1, 1)
        # gather local array to rank0
        freq_mask = mpiutil.gather_array(freq_mask, axis=1, root=0, comm=ts.comm)
        if mpiutil.rank0:
            freq_mask = np.sum(freq_mask, axis=1)

        if plot_stats and mpiutil.rank0:
            time_fig_name = '%s_%s.png' % (fig_prefix, 'time')
            if tag_output_iter:
                time_fig_name = output_path(time_fig_name, iteration=self.iteration)
            else:
                time_fig_name = output_path(time_fig_name)

            # plot time_mask
            plt.figure()
            fig, ax = plt.subplots()
            x_vals = np.array([ datetime.fromtimestamp(s) for s in ts['sec1970'][:] ])
            xlabel = '%s' % x_vals[0].date()
            x_vals = mdates.date2num(x_vals)
            #ax.plot(x_vals, 100*time_mask/np.float(nf*nb))
            ax.bar(x_vals,100*time_mask/np.float(nf*nb),width=0.0001,align='center')
            ax.xaxis_date()
            date_format = mdates.DateFormatter('%H:%M')
            ax.xaxis.set_major_formatter(date_format)
            if rotate_xdate:
                # set the x-axis tick labels to diagonal so it fits better
                fig.autofmt_xdate()
            else:
                # reduce the number of tick locators
                locator = MaxNLocator(nbins=6)
                ax.xaxis.set_major_locator(locator)
                ax.xaxis.set_minor_locator(AutoMinorLocator(2))

            ax.set_xlabel(xlabel)
            ax.set_ylabel(r'RFI (%)')
            plt.savefig(time_fig_name)
            plt.close()

            freq_fig_name = '%s_%s.png' % (fig_prefix, 'freq')
            deal_hist_fig_name = '%s_%s.png' %(fig_prefix,'dhist')
            hist_fig_name = '%s_%s.png' % (fig_prefix, 'hist')
            rvi = '%s_%s.png' % (fig_prefix, 'rvi0')
            rvi_norm = '%s_%s.png' % (fig_prefix, 'rvi0_norm')
            
            if tag_output_iter:
                freq_fig_name = output_path(freq_fig_name, iteration=self.iteration)
                hist_fig_name = output_path(hist_fig_name, iteration=self.iteration)
                deal_hist_fig_name = output_path(deal_hist_fig_name, iteration=self.iteration)
                rvi = output_path(rvi,iteration=self.iteration)
                rvi_norm = output_path(rvi_norm,iteration=self.iteration)
            else:
                freq_fig_name = output_path(freq_fig_name)
                hist_fig_name = output_path(hist_fig_name)
                deal_hist_fig_name = output_path(deal_hist_fig_name)
                rvi = output_path(rvi)
                rvi_norm = output_path(rvi_norm)

            # plot freq_mask
            plt.figure()
            #plt.plot(ts.freq[:], 100*freq_mask/np.float(nt*nb))
            plt.bar(ts.freq[:], 100*freq_mask/np.float(nf*nb),width=0.1,align='center')
            plt.xlabel(r'$\nu$ / MHz')
            plt.ylabel(r'RFI (%)')
            plt.savefig(freq_fig_name)
            plt.close()
            
            #plot hist
            
            vis_copy = vis[vis_mask].copy()
            print '#'*20
            print 'vis_copy',vis_copy.shape
            print '#'*20
            plt.figure()

            #n, bins, patches = plt.hist(np.abs(vis_copy), 100, log=True,facecolor = 'yellow', alpha=0.75,label = "Total")
            #plt.legend(loc = 'best')
            #plt.xlabel('value of data')
            #plt.ylabel('statistic of number') 
            #plt.grid(True)
            #plt.savefig(hist_fig_name)
            #plt.close()

            #plot removel uncertain data around ns source
            data = ts['ns_on'][:].copy()
            data = list(data)
            data_before_ind = ts['ns_on'][:].copy()
            for ind,val in enumerate(data):
                if val == True:
                    data_before_ind[ind] = False
                    data_before_ind[ind-5] = True
                    



            index =[]
            for i,j in enumerate(data):
                if j == True:
                    index.append(i)
            for ii in range(len(index)/5):
                vis_mask[index[5*ii]-5:index[5*ii]+9+1] = False
            vis_copy_1 = vis[vis_mask].copy()
                         
            print '#'*20
            print 'vis_copy_1',vis_copy_1.shape
            print '#'*20
	    
	    """
            #plt.figure()
            #n, bins, patches = plt.hist(np.abs(vis_copy_1), 100, log=True,facecolor='g', alpha=0.75,label = "removol near noise source")
            plt.legend(loc = 'best')
            plt.xlabel("Vaule of RFI's data")
            plt.ylabel("Numerical statistic")
            plt.grid(True)
            plt.savefig(deal_hist_fig_name)
            plt.close()
            
            
            
            #ff = h5py.File('value_12-18_index1_with_noise_plot.hdf5','w') 
            
            #plot
            n = 0
            #n = 1
            plt.figure(figsize=(8,8))
            plt.grid(linestyle = "--")      
            ax = plt.gca()
            ax.spines['top'].set_visible(False)  
            ax.spines['right'].set_visible(False) 
            #f,ax = plt.subplots()
            #rvi_vis = rvi_vis[rvi_vis_mask]
            j = 0
            markers = ['o', '*', '+', 'x','s', 'p']*100
            col = ['b','r','g','c','y','m','k','b','r','g','c','y','m','k']*100
            print 'before rvis_vis.shape',rvi_vis.shape
            #add_value = np.zeros_like(rvi_vis[:,:,0][rvi_vis_mask[:,:,0]])
            #add_value = np.zeros_like(rvi_vis[:,8,0][rvi_vis_mask[:,8,0]])
            #add_value = np.zeros_like(rvi_vis[data,2,0])
            #print 'add_value.shape',add_value.shape
            for i in [0,2,3,4,5,7,8,9,10]:      #rvi_vis.shape[2]):
                j+=1
                #value = rvi_vis[:,8,i][rvi_vis_mask[:,8,0]]
                #value = rvi_vis[:,[0,1,2,3,4,5,6,10,11,12],i][rvi_vis_mask[:,[0,1,2,3,4,5,6,10,11,12],0]]
                value = rvi_vis[data,2,i]
                #signal_value = rvi_vis[data_before_ind,2,i]
                #value = value - signal_value
                # gather local array to rank0
                # gather local array to rank0
                #value = value.reshape(-1,1)
                #value = mpiutil.gather_array(value, axis=1, root=0, comm=ts.comm)
                #if mpiutil.rank0:
                    #value = np.sum(value, axis=1)

                print 'value.shape_after',value.shape
                print 'value',value
                #add_value = np.vstack((add_value,value))
                #add_value = np.hstack((add_value,value))
                plt.scatter(value.real, value.imag, s = 5, marker=markers[i],c=col[i],alpha=0.8,label=str(i)+'-baseline')
            #ff.create_dataset("value", data = add_value)
            #ff.create_dataset("rvi_vis_mask",data = rvi_vis_mask[:,:,0])
            #ff.close()
            """


#PLOT AND SAVE DATA

            plt.legend(loc = 'best')
            print '@'*30
            print 'j', j
            plt.title("noise source",fontsize=12,fontweight='bold')
            plt.xlabel("real axis")
            plt.ylabel("imag axis")
            plt.grid(True)
            plt.savefig(rvi)
            plt.close()

            #1 plot normalize
            n = 0
            plt.figure(figsize=(8,8))
            plt.grid(linestyle = "--")      
            ax = plt.gca()
            ax.spines['top'].set_visible(False)  
            ax.spines['right'].set_visible(False) 
            markers = ['o', '*', '+', 'x','s', 'p']*100
            col = ['b','r','g','c','y','m','k']*100
            #add_value_new = np.zeros_like(rvi_vis[data,2,0])
            add_value_raw = np.zeros_like(rvi_vis[data,2,0])
	    #value_raw = rvi_vis[data,1,0]
	    add_value_raw = np.vstack((add_value_raw,rvi_vis[data,0,0]))
	    add_value_raw = np.vstack((add_value_raw,rvi_vis[data,1,0]))
            add_value_raw = np.vstack((add_value_raw,rvi_vis[data,2,0]))
            add_value_raw = np.vstack((add_value_raw,rvi_vis[data,3,0]))
            
            add_signal_value = np.zeros_like(rvi_vis[data_before_ind,2,0])
	    add_signal_value = np.vstack((add_signal_value,rvi_vis[data_before_ind,0,0]))
	    add_signal_value = np.vstack((add_signal_value,rvi_vis[data_before_ind,1,0]))
	    add_signal_value = np.vstack((add_signal_value,rvi_vis[data_before_ind,2,0]))
	    add_signal_value = np.vstack((add_signal_value,rvi_vis[data_before_ind,3,0]))
	    """
            for i in range(rvi_vis.shape[2]):  #rvi_vis.shape[2]-10)
                #value = rvi_vis[:,:,i][rvi_vis_mask[:,:,0]]
                #value = rvi_vis[:,[0,1,2,3,4,12,13,14,15,16,17],i][rvi_vis_mask[:,[0,1,2,3,4,12,13,14,15,16,17],0]]
                value_raw = rvi_vis[data,1,i]
                signal_value = rvi_vis[data_before_ind,1,i]
                #value_new = value_raw - signal_value
                # gather local array to rank0
                #add_value_new = np.vstack((add_value_new,value_new))
                add_value_raw = np.vstack((add_value_raw,value_raw))
                add_signal_value = np.vstack((add_signal_value,signal_value))
                #plt.scatter(value.real/np.abs(np.sum(value)), value.imag/np.abs(np.sum(value)), s = 3, marker=markers[i],c=col[i],alpha=0.8,label='feed'+str(bl_label[i][0])+'and'+str(bl_label[i][1]))
                print 'np.angle(value_raw[500:550])',np.angle(value_raw[500:550])
		plt.scatter(np.arange(value_raw.shape[0]), np.angle(value_raw), s = 2, marker=markers[i],c=col[i],alpha=0.8,label='feed'+str(bl_label[i][0])+'and'+str(bl_label[i][1]))
	    """
            import time
	    a = time.ctime().split(' ')
            T = a[-1]+'_'+a[1]+'_'+a[2]+'_'+a[3]
	    f = h5py.File( T+'.hdf5','w')
            f.create_dataset('value_raw',data = add_value_raw)
            f.create_dataset('signal_value',data = add_signal_value)
            #f.create_dataset('value_new',data = add_value_new)
	    f.create_dataset('bl_label',data = bl_label)
	    f.create_dataset('feed_pos',data = feed_pos)
            f.close()


            plt.legend(loc = 'best')
            plt.xlabel("real axis")
            plt.ylabel("imag axis")
            #plt.xlim((-0.001,0.001))
            #plt.ylim((-0.0005,0.0005))
            plt.title("noise source",fontsize=12,fontweight='bold')
            plt.grid(True)
            plt.savefig(rvi_norm)
            plt.savefig(output_path('stats_rvi0_norm.svg'),format = 'svg')
            plt.close()
            
            #2 different frequency point 
            plt.figure(figsize=(8,8))
            plt.grid(linestyle = "--")      
            ax = plt.gca()
            ax.spines['top'].set_visible(False)  
            ax.spines['right'].set_visible(False) 
            markers = ['o', '*', '+', 'x','s', 'p']*100
            col = ['b','r','g','c','y','m','k']*100
            for i in [0,1,2,3,4,5,6]:  #rvi_vis.shape[2]-10)
                value = rvi_vis[data,1,i]
                plt.scatter(value.real/np.abs(np.sum(value)), value.imag/np.abs(np.sum(value)), s = 3, marker=markers[i],c=col[i],alpha=0.8,label='feed'+str(bl_label[i][0])+'and'+str(bl_label[i][1]))

            plt.legend(loc = 'best')
            plt.xlabel("real axis")
            plt.ylabel("imag axis")
            #plt.xlim((-0.001,0.001))
            plt.title("noise source",fontsize=12,fontweight='bold')
            plt.grid(True)
            plt.savefig(output_path('value=rvi_vis[data,1,i].pdf'))
            plt.savefig(output_path('value=rvi_vis[data,1,i].svg'), format = 'svg')
            plt.close()

            #3 different baseline seems have diffuse point
            plt.figure(figsize=(8,8))
            plt.grid(linestyle = "--")      
            ax = plt.gca()
            ax.spines['top'].set_visible(False)  
            ax.spines['right'].set_visible(False) 
            markers = ['o', '*', '+', 'x','s', 'p']*100
            col = ['b','b','b','r','g','g','g','c','y','m','k']*100
            for i in [0,3,4,7,10]:  #rvi_vis.shape[2]-10)
                value = rvi_vis[data,2,i]
                plt.scatter(value.real/np.abs(np.sum(value)), value.imag/np.abs(np.sum(value)), s = 3, marker=markers[i],c=col[i],alpha=0.8,label='feed'+str(bl_label[i][0])+'and'+str(bl_label[i][1]))

            plt.legend(loc = 'best')
            plt.xlabel("real axis")
            plt.ylabel("imag axis")
            #plt.xlim((-0.001,0.001))
            plt.title("noise source",fontsize=12,fontweight='bold')
            plt.grid(True)
            plt.savefig(output_path('for_i_in_[0,3,4,7,10].pdf'))
            plt.savefig(output_path('for_i_in_[0,3,4,7,10].svg'), format = 'svg')
            plt.close()
            """
	    #have RFI's frequency noise plot
            plt.figure(figsize=(8,8))
            plt.grid(linestyle = "--")      
            ax = plt.gca()
            ax.spines['top'].set_visible(False)  
            ax.spines['right'].set_visible(False) 
            markers = ['o', '*', '+', 'x','s', 'p']*100
            col = ['b','r','g','c','y','m','k']*100
            for i in [0,1,2,3,4,5,6]:  #rvi_vis.shape[2]-10)
                value = rvi_vis[data,8,i]
                plt.scatter(value.real/np.abs(np.sum(value)), value.imag/np.abs(np.sum(value)), s = 3, marker=markers[i],c=col[i],alpha=0.8,label='feed'+str(bl_label[i][0])+'and'+str(bl_label[i][1]))

            plt.legend(loc = 'best')
            plt.xlabel("real axis")
            plt.ylabel("imag axis")
            #plt.xlim((-0.001,0.001))
            plt.title("noise source",fontsize=12,fontweight='bold')
            plt.grid(True)
            plt.savefig(output_path('value=rvi_vis[data,8,i]_high_rfi_contaminate.pdf'))
            plt.savefig(output_path('value=rvi_vis[data,8,i]_high_rfi_contaminate.svg'), format = 'svg')
            plt.close()
	    """
        return super(Stats, self).process(ts)
    def realvsimag(self,vis,vis_mask,li,gi,tf,ts,**kwargs):
        pass
        """
        print 'rvi'*20
        print 'vis.shape',vis.shape
        print 'vis_mask.shape',vis_mask.shape
        print 'li',type(li),li
        print 'gi',type(gi),gi
        print 'tf',type(tf)
        print 'ts',type(ts)
        """
