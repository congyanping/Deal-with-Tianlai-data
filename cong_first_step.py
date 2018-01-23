import numpy as np
import warnings
import timestream_task
from tlpipe.container.raw_timestream import RawTimestream
from tlpipe.container.timestream import Timestream
import scipy.signal as signal
import h5py
import sys
import math
import time
class Test(timestream_task.TimestreamTask):
    params_init={}
    prefix = 'te_'
    def process(self,ts):
        ts.bl_data_operate(self.main)
        return super(Test,self).process(ts)
    def func1(self,x=0,y=0):
        chaM_0 = 2000
        listt=[]
        for i in xrange(min(vis.shape[0],vis.shape[1])/2*max(vis.shape[0],vis.shape[1]/2)):
            i=chaM_0/rho**np.log2(i)
            listt.append(i)
        return listt
    def func2(self,x,y):
        pass
    def main(self,vis,vis_mask,li,gi,bl,ts,**kwargs):
        self.vis = vis
        self.vis_mask = vis_mask
        v = 700.0
        deltaV = 0.3333
        print "vis.shape",vis.shape
        for chaM in self.func1:
            for i in xrange(vis.shape[0]):
                for j in xrange(vis.shape[1]):
                    if v+(i-j)*deltaV < 780 and i < vis.shape[0]:
                        if abs(self.vis[v+(i-j)*dealtV,i]) > chaM:
                            self.vis_mask[i,v+(i-j)*dealtv] = True

