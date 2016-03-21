"""
Created on 17.08.2013

:author: Sebastian Illing
:contact: sebastian.illing@met.fu-berlin.de

Copyright (C) 2014  Sebastian Illing
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>
"""
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as N
from mpl_toolkits.basemap import Basemap, cm

from filehandler import FileHandler


class Plotter(object):

    colorDict =  {'red':   ((0., 0, 0),  (0.35, 0, 0),  (0.45, 1, 1), (0.5, 1, 1), (0.55, 1, 1), (0.60, 0.1, 0.1),
                            (0.65, 1.0, 1.0), (0.89, 1, 1), (1, 0.5, 0.5)),
                  'green': ((0., 0, 0), (0.125, 0, 0), (0.375, 1, 1), (0.45, 1, 1), (0.5, 1, 1), (0.55, 1, 1),
                            (0.60, 1, 1), (0.65, 1, 1), (0.91, 0, 0), (1, 0, 0)),
                  'blue':  ((0., 0.5, 0.5), (0.11, 1, 1), (0.34, 1, 1), (0.45, 1, 1), (0.5, 1, 1), (0.55, 1, 1),
                            (0.56, 0, 0), (1, 0, 0))}
    
    colorDict2 = {'blue': ((0.0, 0.96862745285034180, 0.96862745285034180),
                           (0.2, 0.78039216995239258, 0.78039216995239258),
                           (0.4, 0.50980395078659058, 0.50980395078659058), 
                           (0.6, 0.30196079611778259, 0.30196079611778259), 
                           (0.8, 0.16862745583057404, 0.16862745583057404),
                           (1.0, 0.12156862765550613, 0.12156862765550613),),
                            
                  'green': ((0, 0.9686274528503418, 0.9686274528503418),
                            (0.2, 0.85882353782653809, 0.85882353782653809),
                            (0.4, 0.64705884456634521, 0.64705884456634521),
                            (0.6, 0.37647059559822083, 0.37647059559822083),
                            (0.8, 0.094117648899555206, 0.094117648899555206),
                            (1.0, 0.0, 0.0)),
                            
                  'red': ((0, 0.9686274528503418, 0.9686274528503418),
                          (0.2, 0.99215686321258545, 0.99215686321258545),
                          (0.4, 0.95686274766921997, 0.95686274766921997),
                          (0.6, 0.83921569585800171, 0.83921569585800171),
                          (0.8, 0.69803923368453979, 0.69803923368453979),
                          (1.0, 0.40392157435417175, 0.40392157435417175))
                  }
    
    colorDict2 = {'blue': [(0.0, 0.3803921639919281, 0.3803921639919281),
                           (0.2, 0.67450982332229614, 0.67450982332229614),
                           (0.4, 0.76470589637756348, 0.76470589637756348),
                  
                  (0.6, 0.87058824300765991, 0.87058824300765991),
                  (0.8, 0.94117647409439087, 0.94117647409439087),
                  (1, 0.9686274528503418, 0.9686274528503418)],
                            
                  'green': [(0.0, 0.18823529779911041, 0.18823529779911041),
                            (0.2, 0.40000000596046448, 0.40000000596046448),
                            (0.4, 0.57647061347961426, 0.57647061347961426),
                            (0.6, 0.77254903316497803, 0.77254903316497803),
                            (0.8, 0.89803922176361084, 0.89803922176361084),
                            (1.0, 0.9686274528503418, 0.9686274528503418), ],
                            
                  'red': [
                          (0.0, 0.019607843831181526, 0.019607843831181526),
                          (0.2, 0.12941177189350128, 0.12941177189350128),
                          (0.4, 0.26274511218070984, 0.26274511218070984),
                          (0.6, 0.57254904508590698, 0.57254904508590698),
                          (0.8, 0.81960785388946533, 0.81960785388946533),
                          (1.0, 0.9686274528503418, 0.9686274528503418)]}

    @staticmethod
    def plotField(fileName, vmin, vmax, colormap='goddard', output_folder='/', lonlatbox=None, region='global'):
        """
        Plot any field variable
        
        :param fileName: filepath
        :param vmin: min value for colorbar
        :param vmax: max value for colorbar
        """
        fig, ax = plt.subplots(figsize=(12, 7), dpi=500)
        if colormap == 'goddard':
            my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap', Plotter.colorDict2, 256)
        else:
            my_cmap = plt.cm.RdBu_r
        
        file_values = FileHandler.openNetCDFFile(fileName)
        mVar = file_values['variable']
        lon = file_values['lon']
        lat = file_values['lat']

        step = (lon[1]-lon[0])/2.
        
        lon = N.insert(lon, lon.shape[0], lon[-1]+2*step)
        lon = N.insert(lon, lon.shape[0], lon[-1]+2*step)
        
        lat = N.insert(lat, lat.shape[0], lat[-1]+2*step)
        lat = N.insert(lat, lat.shape[0], lat[-1]+2*step)

        mVarF = mVar[:, 0:2]
        mVar = N.concatenate((mVar, mVarF), axis=1)
        
        mVarF = mVar[0:2, :]
        mVar = N.concatenate((mVar, mVarF), axis=0)

        def plus(lon):
            return lon-step
        
        if lonlatbox is None:
            if lon[0] < 0 and lon[1] < 0:
                lonlatbox = '-180,180,-90,90'
            else:
                lonlatbox = '0,360,-90,90'
        lon = N.array(map(plus, lon))
        lat = N.array(map(plus, lat))
        
        lonlatbox = map(int, lonlatbox.split(','))
        # check which lat is greater
        if lonlatbox[2] > lonlatbox[3]:
            lonlatbox[2], lonlatbox[3] = lonlatbox[3], lonlatbox[2]  # switch position on lat entries
        
        if lonlatbox[0] > lonlatbox[1]:
            lonlatbox[1] += 360
        
        step = 0
        
        if region == 'Arctic':
            m = Basemap(projection='npstere', boundinglat=lonlatbox[2], lon_0=0)
            parallels = N.arange(0., 90., 10.)
            meridians = N.arange(0., 360., 10.)
        elif region == 'Antarctica':
            m = Basemap(projection='spstere', boundinglat=lonlatbox[3], lon_0=0)
            parallels = N.arange(-90., 0., 10.)
            meridians = N.arange(0., 360., 10.)
        else:
            m = Basemap(llcrnrlon=lonlatbox[0], llcrnrlat=lonlatbox[2], urcrnrlon=lonlatbox[1], urcrnrlat=lonlatbox[3])
            parallels = N.arange(-90., 120., 30.)
            meridians = N.arange(0., 420., 60.)

        def divi(x):
            return float(x)/10
        
        colorSteps = map(divi, range(int(vmin*10), int((vmax*10)+1), 1))
        
        if vmax == 0.5:
            colorSteps = [-0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0, 
                          0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]
            
        if vmax == 0. and vmin == -0.5:
            colorSteps = [-0.5, -0.475, -0.45, -0.425, -0.4, -0.375, -0.35, -0.325, -0.3, -0.275, -0.25, -0.225, 
                          -0.2, -0.175, -0.15, -0.125, -0.1, -0.075, -0.05, -0.025, 0]
            
        if vmax == 2. and vmin == 0:
            colorSteps = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
                          1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]
        
        if vmin == -0.6:
            colorSteps = [-0.6, -0.55, -0.5, -0.45, -0.4,-0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0,
                          0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6]
        
        if vmin == -0.1 and vmax == 0:
            colorSteps = [-0.1, -0.09, -0.08, -0.07, -0.06, -0.05, -0.04, -0.03, -0.02, -0.01, 0.0]
            
        if vmin == -0.2 and vmax == 0:
            colorSteps = [-0.2, -0.18, -0.16, -0.14, -0.12, -0.1, -0.08, -0.06, -0.04, -0.02, 0.0]
        
        colorTicks = colorSteps[0::2]
        my_cmap.set_bad("grey")  # set missing value color
        maskedArray = N.ma.masked_outside(mVar, -0.8e20, 0.8e20)
        # discrete colormap
        norm = mpl.colors.BoundaryNorm(colorSteps, my_cmap.N)
        x, y = lon,lat
        x, y = N.meshgrid(lon, lat)
        cs = m.pcolormesh(x, y, maskedArray, cmap=my_cmap, norm=norm, linewidth=0, rasterized=True)
        cs.set_edgecolor('face')
        cb = m.colorbar(cs, "right", size="5%", pad='5%', ticks=colorTicks)
        m.drawcoastlines(ax=ax)  
        m.drawparallels(parallels, labels=[1, 0, 0, 0])  # draw parallels
        m.drawmeridians(meridians, labels=[0, 0, 0, 1])  # draw meridians
        
        plt.title(Plotter.__getTitle(fileName))
        plt.text(lonlatbox[0]+(lon[1]-lon[0])/2, lonlatbox[2]+(lat[1]-lat[0])/2, 'MurCSS')
        return m

    @staticmethod
    def plotVerticalProfile(fileName, vmin, vmax, colormap='goddard', lonlatbox=None):
        """
        Plots vertical profile of zonal mean file
        """
        file_values = FileHandler.openNetCDFFile(fileName, 'plev')
        data = N.flipud(file_values['variable'])
        lat = file_values['lat']
        plev = N.flipud(file_values['plev'])
        
        fig, ax1 = plt.subplots(figsize=(12, 7), dpi=500)
        colorSteps = N.linspace(-1, 1, 21)
        if vmax == 0.5:
            colorSteps = [-0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0, 
                          0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]   
        if vmax == 0. and vmin == -0.5:
            colorSteps = [-0.5, -0.475, -0.45, -0.425, -0.4, -0.375, -0.35, -0.325, -0.3, -0.275, -0.25, -0.225, 
                          -0.2, -0.175, -0.15, -0.125, -0.1, -0.075, -0.05, -0.025, 0]     
        if vmax == 2. and vmin == 0:
            colorSteps = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
                          1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]
        if vmin == -0.6:
            colorSteps = [-0.6, -0.55, -0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0,
                          0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6]
        if vmin == -0.1 and vmax == 0:
            colorSteps = [-0.1, -0.09, -0.08, -0.07, -0.06, -0.05, -0.04, -0.03, -0.02, -0.01, 0.0]
        if vmin == -0.2 and vmax == 0:
            colorSteps = [-0.2, -0.18, -0.16, -0.14, -0.12, -0.1, -0.08, -0.06, -0.04, -0.02, 0.0]
        my_cmap = plt.cm.RdBu_r
        my_cmap.set_bad("grey")  # set missing value color
        maskedArray = N.ma.masked_outside(data, -0.8e20, 0.8e20)
        # discrete colormap
        norm = mpl.colors.BoundaryNorm(colorSteps, my_cmap.N)
        cs = ax1.imshow(maskedArray, interpolation="nearest", cmap=my_cmap, norm=norm)
        cbar = fig.colorbar(cs, ax=ax1, orientation='vertical', ticks=colorSteps[0::2])
        ax1.set_xlabel('Latitude [degrees]')
        ax1.set_ylabel("Pressure [Pa]")
        plt.xticks(range(1, len(lat), 4), lat[1::4])
        plt.yticks(range(0, len(plev), 2), plev[0::2])
        plt.title(Plotter.__getTitle(fileName))
        return ax1
        
    @staticmethod
    def plotLeadtimeseries(resultList, flag_list, plot_list, input_flags=None):
        colors = ['green', 'blue', 'red']

        def getFn(s):
            return s.split('/')[-1]

        if input_flags is None:
            input_flags = ['']*len(flag_list)
        fig = plt.figure(figsize=(17, 13))
        for i, needle in enumerate(plot_list):
            plt.subplot(2, 1, i+1)
            for j, flag in enumerate(flag_list):
                files_to_plot = list()
                for part in resultList:
                    search_needle = flag+'_'+needle[0]
                    tmpfiles = [s for s in part if input_flags[j] in s]
                    files_to_plot.append([s for s in tmpfiles if search_needle in getFn(s)][0])
                labels = list()
                plot_values = list()
                for fn in sorted(files_to_plot, key=lambda x: int(getFn(x).split('_')[0]+getFn(x).split('_')[1])):
                    tmp_label = getFn(fn).split('_')
                    if tmp_label[0] == tmp_label[1]:
                        tmp_label = tmp_label[0]
                    else:
                        tmp_label = tmp_label[0]+'-'+tmp_label[1]
                    labels.append(tmp_label)
                    tmp_values = FileHandler.openNetCDFFile(fn, 'var')
                    plot_values.append(tmp_values)
                x_val = range(1, len(files_to_plot)+1)
                plt.scatter(x_val, plot_values, color=colors[j])
                plt.plot(x_val, plot_values, color=colors[j])
                
                # get min and max plot values
                try:
                    min_val = min(min_val, min(plot_values))
                    max_val = max(max_val, max(plot_values))
                except:
                    min_val = min(plot_values)
                    max_val = max(plot_values)
                
            plt.axis([0, len(files_to_plot)+1, min(min_val, needle[2][0]), max(max_val, needle[2][1])])
            plt.ylabel(needle[1])
            plt.xlabel('Leadtimes')
            plt.xticks(x_val, labels)
            if i==0:
                plt.legend(flag_list,bbox_to_anchor=(0., 1.05, 1., .102), mode='expand',loc=3,ncol=1, borderaxespad=0.)

    @staticmethod
    def plotLeadtimeseriesSign(resultList, flag_list, plot_list):
        colors = ['green', 'blue', 'red']

        def getFn(s):
            return s.split('/')[-1]
        fig = plt.figure(figsize=(17, 13))
        for i, needle in enumerate(plot_list):
            plt.subplot(2, 1, i+1)
            for j, flag in enumerate(flag_list):
                files_to_plot = list()
#                for part in resultList:
                search_needle = flag+'_'+needle[0]
                files_to_plot += [s for s in resultList if search_needle in getFn(s)]
                labels = list()
                plot_values = list()
                min_plot_values = list()
                max_plot_values = list()
                for fn in sorted(files_to_plot, key=lambda x: int(getFn(x).split('_')[0]+getFn(x).split('_')[1])):
                    tmp_label = getFn(fn).split('_')
                    if tmp_label[0] == tmp_label[1]:
                        tmp_label = tmp_label[0]
                    else:
                        tmp_label = tmp_label[0]+'-'+tmp_label[1]
                    labels.append(tmp_label)
                    tmp_values = FileHandler.openNetCDFFile(fn, 'var')
                    plot_values.append(tmp_values)
                    tmp_min_value = FileHandler.openNetCDFFile(fn+'_bootstrap_min_val', 'var')
                    min_plot_values.append(N.abs(tmp_values-tmp_min_value))
                    tmp_max_value = FileHandler.openNetCDFFile(fn+'_bootstrap_max_val', 'var')
                    max_plot_values.append(N.abs(tmp_values-tmp_max_value))
                x_val = range(1, len(files_to_plot)+1)
                plt.errorbar(x_val, plot_values, yerr=[min_plot_values, max_plot_values], color=colors[j])
                plt.scatter(x_val, plot_values, color=colors[j])

                # get min and max plot values
                try:
                    min_val = min(min_val, min(plot_values))
                    max_val = max(max_val, max(plot_values))
                except:
                    min_val = min(plot_values)
                    max_val = max(plot_values)
                
            plt.axis([0, len(files_to_plot)+1, min(min_val, needle[2][0]), max(max_val, needle[2][1])])
            plt.ylabel(needle[1])
            plt.xlabel('Leadtimes')
            plt.xticks(x_val, labels)

            if i==0:
                plt.legend(flag_list, bbox_to_anchor=(0., 1.05, 1., .102), mode='expand', loc=3, ncol=1, borderaxespad=0.)

    @staticmethod    
    def saveFig(output_folder, fn):
        
        plt.savefig(output_folder+fn+'.eps', format='eps')
        
    @staticmethod
    def addCrosses(map, sig_mask_x, sig_mask_y, marker='x', color='k', size=9):
        # check for half crosses and append
        for i, x in enumerate(sig_mask_x):
            if x == -180:
                sig_mask_x.append(180)
                sig_mask_y.append(sig_mask_y[i])
        map.scatter(sig_mask_x, sig_mask_y, size, marker=marker, color=color)
        
    @staticmethod
    def addCrossesXY(map, fn, marker='x', color='k', size=9):
        sign_mask = FileHandler.openNetCDFFile(fn, 'var')
        sign_mask = N.flipud(sign_mask)
        sig_x = list()
        sig_y = list()
        for (x, y), value in N.ndenumerate(sign_mask):
            if value == 1:
                sig_x.append(y)
                sig_y.append(x)

        map.scatter(sig_x, sig_y, size, marker=marker, color=color)
    
    @staticmethod
    def addCrossesFile(map, fn, marker='x', color='k', size=9):
        sign_mask = FileHandler.openNetCDFFile(fn)
        lon = sign_mask['lon']
        lon = lon + (lon[1]-lon[0])/2
        lat = sign_mask['lat']
        lat = lat + (lat[1]-lat[0])/2
        sign_mask = sign_mask['variable']
        sig_x = list()
        sig_y = list()
        for (x, y), value in N.ndenumerate(sign_mask):
            if value == 1:
                sig_x.append(lon[y])
                sig_y.append(lat[x])

        map.scatter(sig_x, sig_y, size, marker=marker, color=color)

    @staticmethod    
    def __getTitle(fn):    
        
        fn = fn.split('/')[-1]
        fn = fn.split('_')        
        title = ''
        length = 0
        for part in fn:
            title += part+'_'
            length += len(part)
            if length >= 60 or part == 'vs':
                title += '\n'
                length = 0
        return title
        
        
if __name__ == '__main__':
    print 'test'
    fn = '/var/autofs/net/home/illing/workspace/murcss/integration/metrics/1_1_tos_baseline1_output_mpi-esm-lr_decs4e_1990-1992_correlation.nc_masked'
    m=Plotter.plotField(fn, -1, 1, 'RdBu')
    Plotter.addCrosses(m, [-180, 150], [10, 10])
    Plotter.saveFig('', 'test')
    plt.show()      
