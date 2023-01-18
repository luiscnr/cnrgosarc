import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import os
import json


class ArcAnalysis():

    def __init__(self, arc_opt, verbose, avg_file, resample_folder):
        self.verbose = verbose
        self.arc_opt = arc_opt
        self.avg_file = avg_file
        self.resample_folder = resample_folder
        self.list_images = {}

        # for graphics
        self.fontsizemultiple = 12
        self.line_size_multiple = 1
        self.marker_size_multiple = 5

    def check_n_overlappping(self):
        dataset = Dataset(self.avg_file, 'r')
        varnum = np.array(dataset.variables['sum_weights'][:])
        narray = np.unique(varnum)
        narray = narray[narray >= 2]
        narray = narray.astype('int32')
        ln = list(narray[:])
        return ln

    def check_overlapping_index(self, noverlap):
        dir_out = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/INTEGRATED/2019/175'
        dataset = Dataset(self.avg_file, 'r')
        varnum = np.array(dataset.variables['sum_weights'][:])
        dataset.close()
        combinations_done = np.zeros(varnum.shape)
        self.get_info_images_from_resampleinfo()

        all_images = list(self.list_images.keys())
        nall = len(all_images)

        for idx in range(nall):
            image = all_images[idx]
            indices = self.get_indices_overlap(image, varnum, noverlap)
            ntotimage = indices.sum()
            if self.list_images[image]['n_over'] == ntotimage:
                continue
            if ntotimage == 0:
                continue
            print('----------------------------------', idx, '/', nall,
                  '------------------------------------------------')
            # ymin, ymax, xmin, xmax = self.get_limits_overlaping_images(image)
            ymin, ymax, xmin, xmax = self.get_limits_from_list(image)
            index_combination = 0
            # valid_mask2_list = self.get_valid_masks(idx,varnum,noverlap)
            mask1, mask2_list = self.get_valid_masks(idx, varnum, noverlap)

            while self.list_images[image]['n_over'] < ntotimage:
                print('->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ', image, idx, '/', nall, '->',
                      self.list_images[image]['n_over'], '/', ntotimage)

                ngranulesc = 1
                index_combination = index_combination + 1
                imagesc = [image]
                oindices = None

                for ido in range(idx + 1, nall):
                    oimage = all_images[ido]
                    if self.list_images[oimage]['n_over'] == ntotimage:
                        continue
                    if not self.check_overlapping_images(image, oimage):
                        continue
                    idmask2 = ido - (idx + 1)
                    # valid_mask2 = valid_mask2_list[idmask2]
                    # if not valid_mask2:
                    #     continue
                    # oindices_here, newoindices = self.get_indices_overlap_images(image, oimage, varnum, noverlap,
                    #                                                              oindices, combinations_done)

                    mask2 = mask2_list[idmask2]
                    if mask2 is not None:
                        oindices_here, newoindices = self.get_indices_overlap_imagesv2(image, mask1, mask2, varnum,
                                                                                       noverlap, oindices,
                                                                                       combinations_done)
                    if oindices_here is None:
                        continue
                    if newoindices:
                        oindices = oindices_here
                    noi = len(oindices[0])  # noi is always greater than zero
                    if noi == 0:
                        print('AQUI NO DEBERIA LLEGAR')
                    ngranulesc = ngranulesc + 1
                    imagesc.append(oimage)
                    print('NOI ES: ', noi, ' ngranules', ngranulesc, '/', noverlap, 'ido', ido)
                    if ido == 72 and ngranulesc == 4 and idx == 44:
                        ngranulesc = 5
                    if ngranulesc == noverlap:  # new combination
                        combination_idx = f'C_{noverlap}_{idx}_{index_combination}'
                        dir_combination = os.path.join(dir_out, combination_idx)
                        if not os.path.exists(dir_combination):
                            os.mkdir(dir_combination)
                        print('Combination: ', combination_idx)
                        filec = os.path.join(dir_combination, f'{combination_idx}.json')
                        if not os.path.exists(filec):
                            combination = {
                                'granules': imagesc,
                                'ngranules': ngranulesc,
                                'npixels': noi,
                                'ymin': ymin,
                                'xmin': xmin,
                                'ymax': ymax,
                                'xmax': xmax
                            }
                            with open(filec, "w") as outfile:
                                json.dump(combination, outfile, indent=1)
                        filei = os.path.join(dir_combination, f'{combination_idx}_indices.csv')
                        if not os.path.exists(filei):
                            self.save_incices_as_csv(filei, oindices)
                        for imagec in imagesc:
                            self.list_images[imagec]['n_over'] = self.list_images[imagec]['n_over'] + noi
                            print(imagec, '->', noi, '->>', self.list_images[imagec]['n_over'])

                        combinations_done[ymin:ymax, xmin:xmax][oindices] = 1

                        # restart param for nex combinations
                        break

            # print(image,'->',len(indices[0]),len(indices[1]))

        # indices = np.where(varnum==2)

        # varsensorindices = varsensor[varnum == noverlap]
        # varsensorindices_unique = np.unique(varsensorindices)
        # print(varsensorindices_unique)
        #
        # rrs = np.array(dataset.variables['RRS510'][:])
        # rrs = rrs[np.logical_and(varnum == 2, varsensor == 0)]
        # # rrs = rrs[varsensor==0]
        # print(len(rrs))
        # print(np.mean(rrs))

    def get_indices_overlap(self, image, varnum, noverlap):
        ymin, ymax, xmin, xmax = self.get_limits_from_list(image)
        varnum_here = varnum[ymin:ymax, xmin:xmax]
        filer = os.path.join(self.resample_folder, image)
        dataset = Dataset(filer, 'r')
        mask = np.array(dataset.variables['mask'][ymin:ymax, xmin:xmax])
        indices = np.logical_and(mask == 1, varnum_here == noverlap)
        dataset.close()
        return indices

    def get_valid_masks(self, idx, varnum, noverlap):
        all_images = list(self.list_images.keys())
        nall = len(all_images)
        image = all_images[idx]
        ymin, ymax, xmin, xmax = self.get_limits_from_list(image)
        varnum_here = varnum[ymin:ymax, xmin:xmax]
        file1 = os.path.join(self.resample_folder, image)
        dataset1 = Dataset(file1, 'r')
        mask1 = np.array(dataset1.variables['mask'][ymin:ymax, xmin:xmax])
        dataset1.close()
        valid_mask2list = []
        mask2list = []
        for ido in range(idx + 1, nall):
            print('Getting mask: ', ido, '/', nall)
            image2 = all_images[ido]
            file2 = os.path.join(self.resample_folder, image2)
            dataset2 = Dataset(file2, 'r')
            mask2 = np.array(dataset2.variables['mask'][ymin:ymax, xmin:xmax])
            indices = np.logical_and(mask1 == 1, np.logical_and(varnum_here == noverlap, mask2 == 1))

            if len(indices[0]) > 0:
                valid_mask2list.append(True)
                mask2[mask2 == -999] = 0
                mask2 = mask2.astype(np.bool)
                mask2list.append(mask2)
            else:
                valid_mask2list.append(False)
                mask2list.append(None)
            dataset2.close()
        return mask1, mask2list

    def get_indices_overlap_imagesv2(self, image1, mask1, mask2, varnum, noverlap, oindices, combinations_done):
        ymin, ymax, xmin, xmax = self.get_limits_from_list(image1)
        varnum_here = varnum[ymin:ymax, xmin:xmax]
        combi_here = combinations_done[ymin:ymax, xmin:xmax]

        if oindices is None:
            indices = np.where(np.logical_and(np.logical_and(mask1 == 1, mask2 == 1),
                                              np.logical_and(combi_here == 0, varnum_here == noverlap)))
            if len(indices[0]) == 0:
                return None, False
            else:
                return indices, True
        else:
            nindices = len(oindices[0])
            mask2indices = mask2[oindices]
            nmask = (mask2indices == 1).sum()
            # print('NMASK es: ',nmask)

            if nmask == nindices:
                return oindices, False
            elif 0 < nmask < nindices:
                # indicest = np.where(np.logical_and(mask1 == 1, np.logical_and(mask2 == 1, varnum_here == noverlap)))
                indicesr = []
                indicesc = []
                for idx in range(len(oindices[0])):
                    rt = oindices[0][idx]
                    ct = oindices[1][idx]
                    if mask2[rt, ct] == 1:
                        indicesr.append(rt)
                        indicesc.append(ct)
                indices = (indicesr, indicesc)
                return indices, True
            else:
                return None, False

    def get_indices_overlap_images(self, image1, image2, varnum, noverlap, oindices, combinations_done):
        ymin, ymax, xmin, xmax = self.get_limits_from_list(image1)
        varnum_here = varnum[ymin:ymax, xmin:xmax]
        combi_here = combinations_done[ymin:ymax, xmin:xmax]
        file1 = os.path.join(self.resample_folder, image1)
        dataset1 = Dataset(file1, 'r')
        mask1 = np.array(dataset1.variables['mask'][ymin:ymax, xmin:xmax])
        dataset1.close()
        file2 = os.path.join(self.resample_folder, image2)
        dataset2 = Dataset(file2, 'r')
        mask2 = np.array(dataset2.variables['mask'][ymin:ymax, xmin:xmax])
        dataset2.close()

        if oindices is None:
            indices = np.where(np.logical_and(np.logical_and(mask1 == 1, mask2 == 1),
                                              np.logical_and(combi_here == 0, varnum_here == noverlap)))
            if len(indices[0]) == 0:
                return None, False
            else:
                return indices, True
        else:
            nindices = len(oindices[0])
            mask2indices = mask2[oindices]
            nmask = (mask2indices == 1).sum()
            # print('NMASK es: ',nmask)

            if nmask == nindices:
                return oindices, False
            elif 0 < nmask < nindices:
                # indicest = np.where(np.logical_and(mask1 == 1, np.logical_and(mask2 == 1, varnum_here == noverlap)))
                indicesr = []
                indicesc = []
                for idx in range(len(oindices[0])):
                    rt = oindices[0][idx]
                    ct = oindices[1][idx]
                    if mask2[rt, ct] == 1:
                        indicesr.append(rt)
                        indicesc.append(ct)
                indices = (indicesr, indicesc)
                return indices, True
            else:
                return None, False

    def get_limits_from_list(self, image):
        ymin = self.list_images[image]['ymin']
        ymax = self.list_images[image]['ymax']
        xmin = self.list_images[image]['xmin']
        xmax = self.list_images[image]['xmax']
        return ymin, ymax, xmin, xmax

    def get_limits_overlaping_images(self, image):
        ymin, ymax, xmin, xmax = self.get_limits_from_list(image)
        for imagel in self.list_images:
            if image == imagel:
                continue
            if self.check_overlapping_images(image, imagel):
                yminl, ymaxl, xminl, xmaxl = self.get_limits_from_list(imagel)
                if yminl < ymin:
                    ymin = yminl
                if ymaxl > ymax:
                    ymax = ymaxl
                if xminl < xmin:
                    xmin = xminl
                if xmaxl > xmax:
                    xmax = xmaxl
        return ymin, ymax, xmin, xmax

    def check_overlapping_images(self, image, imagel):
        ymin, ymax, xmin, xmax = self.get_limits_from_list(image)
        yminl, ymaxl, xminl, xmaxl = self.get_limits_from_list(imagel)
        overlap = True
        if yminl > ymax:
            overlap = False
        if ymaxl < ymin:
            overlap = False
        if xminl > xmax:
            overlap = False
        if xmaxl < xmin:
            overlap = False
        return overlap

    def save_incices_as_csv(self, file, indices):
        indicesr = indices[0]
        indicesc = indices[1]
        noi = len(indicesr)
        f = open(file, 'w')
        for idx in range(noi):
            rindex = indicesr[idx]
            cindex = indicesc[idx]
            line = f'{rindex};{cindex}'
            f.write(line)
            f.write('\n')

        f.close()

    def get_info_images_from_resampleinfo(self):

        finfo = os.path.join(self.resample_folder, 'ResampleInfo.csv')
        if not os.path.exists(finfo):
            return
        file1 = open(finfo, 'r')
        lines = file1.readlines()
        file1.close()
        for idx in range(1, len(lines)):
            line = lines[idx]
            vals = line.split(';')
            fname = vals[0].strip()
            fname = fname[:-5] + '_resampled.nc'
            ymin = int(vals[11].strip())
            ymax = int(vals[12].strip())
            xmin = int(vals[13].strip())
            xmax = int(vals[14].strip())
            self.list_images[fname] = {
                'idx': idx,
                'ymin': ymin,
                'ymax': ymax,
                'xmin': xmin,
                'xmax': xmax,
                'n_over': 0
            }

    def compute_average_spectra(self, folderc):
        # integrate_folder = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/INTEGRATED/2019/175'
        olci_l2_bands = [400, 412.5, 442.5, 490, 510, 560, 620, 665, 673.75, 681.25, 708.75, 753.75, 778.75]
        nbands = len(olci_l2_bands)
        average_variables_all = []
        for wl in olci_l2_bands:
            wls = str(wl)
            wls = wls.replace('.', '_')
            bname = f'RRS{wls}'
            average_variables_all.append(bname)
        prenames = ['Sensor', 'Time', 'OZA', 'OAA', 'SAA', 'SZA']
        name = os.path.basename(folderc)
        findices = os.path.join(folderc, f'{name}_indices.csv')
        indicesr = []
        indicesc = []
        nindices = 0
        with open(findices) as f:
            for line in f:
                sline = line.split(';')
                indicesr.append(int(sline[0].strip()))
                indicesc.append(int(sline[1].strip()))
                nindices = nindices + 1
        indices = (indicesr, indicesc)
        finfo = os.path.join(folderc, f'{name}.json')
        with open(finfo) as j:
            info = json.load(j)

        granules = info['granules']
        ymin = info['ymin']
        ymax = info['ymax']
        xmin = info['xmin']
        xmax = info['xmax']
        npixels = info['npixels']
        ngranules = len(granules)
        output_spectra = np.zeros((nindices, nbands))
        avg_spectra_pixel = np.zeros((nindices, nbands))
        avg_spectra = np.zeros((ngranules + 1, nbands))
        max_spectra = np.zeros(nbands)
        prevalues = []
        columns = []
        from datetime import datetime as dt
        time_min = dt.now()
        time_max = dt(2016, 1, 1)
        angle_min = np.zeros((4, 1))
        angle_max = np.zeros((4, 1))

        from matplotlib import pyplot as ptl

        if ngranules == 12:
            nfil = 4
            ncol = 3
            self.fontsizemultiple = 8
            self.line_size_multiple = 0.5

        maxYvalue = 0
        fig, axs = ptl.subplots(nfil, ncol)

        for idx in range(ngranules):
            granule = granules[idx]
            granule_slist = granule.split('_')
            time = dt.strptime(granule_slist[7], '%Y%m%dT%H%M%S')
            if time < time_min:
                time_min = time
            if time > time_max:
                time_max = time
            prevalues_here = [granule_slist[0], granule_slist[7]]
            print(granule, '-----------------------------------------------------')
            fgranule = os.path.join(self.resample_folder, granule)
            dataset = Dataset(fgranule)
            for iband in range(nbands):
                band = average_variables_all[iband]
                varhere = np.array(dataset.variables[band][ymin:ymax, xmin:xmax])
                varherei = varhere[indices]
                avg_spectra[idx, iband] = np.mean(varherei)
                output_spectra[:, iband] = varherei[:]
                max_spectra[iband] = np.max(varherei)
                avg_spectra_pixel[:, iband] = avg_spectra_pixel[:, iband] + varherei[:]

            max_spectra = (max_spectra * 100) / 2
            max_value = (np.ceil(np.max(max_spectra)) * 2) / 100
            if max_value > maxYvalue:
                maxYvalue = max_value

            name_granule = granule[:-3]
            fspectra = os.path.join(folderc, f'Spectra_{name_granule}.csv')
            self.save_spectra_csv(output_spectra, average_variables_all, fspectra, None, None)

            angles = ['OZA', 'OAA', 'SAA', 'SZA']
            for iangle in range(len(angles)):
                angle = angles[iangle]
                varhere = np.array(dataset.variables[angle][ymin:ymax, xmin:xmax])
                varherei = varhere[indices]
                varmean = np.mean(varherei)
                if iangle == 0:
                    oza = varmean
                if idx == 0:
                    angle_min[iangle] = varmean
                    angle_max[iangle] = varmean
                else:
                    if varmean < angle_min[iangle]:
                        angle_min[iangle] = varmean
                    if varmean > angle_max[iangle]:
                        angle_max[iangle] = varmean
                prevalues_here.append(str(varmean))

            title = f'{granule_slist[0]}_{granule_slist[7]}_{oza:.1f}'
            xdata = np.array(olci_l2_bands)
            row = np.floor(idx / ncol)
            col = idx - (row * ncol)
            row = int(row)
            col = int(col)
            print(idx, row, col)
            self.plot_spectra_granule(axs, output_spectra, xdata, row, col, title)

            prevalues.append(prevalues_here)
            columns.append(title)
            dataset.close()

        # save figure with multiple granules
        nticks = int((maxYvalue / 0.02) + 1)
        yticks = np.linspace(0, maxYvalue, nticks)
        xticks = [400, 500, 600, 700, 800]
        for fil in range(nfil):
            for col in range(ncol):
                axs[fil, col].set_ylim(0, maxYvalue)
                axs[fil, col].set_yticks(yticks)
                axs[fil, col].set_xticks(xticks)
                axs[fil, col].grid(b=True, which='major', axis='y', color='gray', linestyle='--')
                axs[fil, col].tick_params(axis='both', which='major', labelsize=self.fontsizemultiple)
        plt.gcf().tight_layout()
        fspectraimg = os.path.join(folderc, f'SpectraGranule_{name}.jpg')
        plt.savefig(fspectraimg, dpi=300)
        plt.close()

        # save average spectra as csv
        fspectra = os.path.join(folderc, f'AvgSpectra_{name}.csv')
        diftime = (time_max - time_min).total_seconds() / 3600
        angle_sum = angle_max - angle_min
        prevalues_here = ['S3', str(diftime), str(np.float(angle_sum[0])), str(np.float(angle_sum[1])),
                          str(np.float(angle_sum[2])), str(np.float(angle_sum[3]))]
        prevalues.append(prevalues_here)
        for iband in range(nbands):
            avg_spectra[ngranules, iband] = np.mean(avg_spectra[0:ngranules, iband])
        self.save_spectra_csv(avg_spectra, average_variables_all, fspectra, prevalues, prenames)

        # save figure with average spectra
        self.fontsizemultiple = 12
        self.line_size_multiple = 0.80
        fig2, axs2 = ptl.subplots(2, 1)
        xdata = np.array(olci_l2_bands)
        avg_spectra_pixel = avg_spectra_pixel / ngranules
        title = f'Avg. Spectra by pixel: {ngranules} N. Pixels: {npixels} Dif. Time: {diftime:.2f}'
        self.plot_spectra_granule(axs2, avg_spectra_pixel, xdata, 1, -1, title)

        columns.append(f'Average_{diftime:.2f}')
        title = f'Avg. Spectra by granule N. Granules: {ngranules} N. Pixels: {npixels} Dif. Time: {diftime:.2f}'
        self.plot_spetra(axs2, avg_spectra, xdata, columns, prevalues, prenames, 0, -1, title)
        # dirimage = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/SpectraExamples

        for fil in range(2):
            axs2[fil].set_ylim(0, maxYvalue)
            axs2[fil].set_yticks(yticks)
            #axs2[fil].set_xticks(xticks)
            axs2[fil].grid(b=True, which='major', axis='y', color='gray', linestyle='--')
            axs2[fil].tick_params(axis='both', which='major', labelsize=self.fontsizemultiple)
        #plt.gcf().legend(columns,loc=8,ncols=ncol)
        plt.gcf().tight_layout()
        fout = os.path.join(folderc, f'AvgSpectra_{name}.jpg')
        plt.savefig(fout, dpi=300)
        plt.close()

    def plot_spectra_granule(self, axs, output_spectra, xdata, row, col, title):

        # selection axes
        if row >= 0 and col >= 0:
            axes_here = axs[row, col]
        if row >= 0 and col == -1:
            axes_here = axs[row]

        # plotting single spectra
        nspectra = output_spectra.shape[0]
        for idx in range(nspectra):
            axes_here.plot(xdata, output_spectra[idx, :], color='gray', lw=self.line_size_multiple)

        # plotting average and std
        avg = np.mean(output_spectra, axis=0)
        std = np.std(output_spectra, axis=0)
        avgminus = avg - std
        avgplus = avg + std
        lavgsize = self.line_size_multiple * 2

        axes_here.plot(xdata, avg, color='black', lw=lavgsize)
        axes_here.plot(xdata, avgminus, color='black', ls='--', lw=self.line_size_multiple)
        axes_here.plot(xdata, avgplus, color='black', ls='--', lw=self.line_size_multiple)
        if title is not None:
            axes_here.set_title(title)
            axes_here.title.set_size(self.fontsizemultiple)

    def plot_spetra(self, axs, avg_spectra, xdata, columnsh, prevalues, prenames, row, col, title):
        from matplotlib import pyplot as plt
        import pandas as pd

        if row >= 0 and col >= 0:
            axes_here = axs[row, col]
        if row >= 0 and col == -1:
            axes_here = axs[row]

        # sorting dataframe
        ngranules = len(columnsh)
        all_avg_spectra = avg_spectra[0:ngranules - 1, :]
        all_columns = columnsh[0:ngranules - 1]
        dfprevvalues = pd.DataFrame(prevalues[0:ngranules - 1], columns=prenames, index=all_columns)
        dfavg = pd.DataFrame(all_avg_spectra, columns=xdata, index=all_columns)
        df_prev_avg = dfprevvalues.join(dfavg)
        df_prev_avg = df_prev_avg.sort_values(['Sensor', 'OZA'])
        colrrs = df_prev_avg.columns[6:19]
        df_fin = df_prev_avg[colrrs]
        #print(df_fin.index)
        df = df_fin.transpose()

        # plotin single avg averages
        df.plot(ax=axes_here, lw=self.line_size_multiple, marker='.', markersize=self.line_size_multiple)

        # ploting average
        avg_avg_spectra = avg_spectra[ngranules - 1, :]
        avg_size = self.line_size_multiple * 2
        axes_here.plot(xdata, avg_avg_spectra, lw=avg_size, marker='.', markersize=self.marker_size_multiple,
                       color='black')
        #axes_here.legend = columnsh
        legend_sorted = df_fin.index.tolist()
        legend_sorted.append(columnsh[ngranules-1])
        axes_here.legend(legend_sorted,ncol=1,bbox_to_anchor=(0.0,0.0,1,1),fontsize=5,fancybox=False)

        if title is not None:
            axes_here.set_title(title)
            axes_here.title.set_size(self.fontsizemultiple)
        # pl axs[row,col].title.set_size(self.fontsizemultiple)t.gcf().tight_layout()
        # plt.savefig(fout, dpi=300)
        # plt.close()

    def save_spectra_csv(self, spectra, name_spectra, filecsv, prevalues, prenames):
        lines = []
        if name_spectra is not None:
            line = ';'.join(name_spectra)
            if prenames is not None:
                linep = ';'.join(prenames)
                line = f'{linep};{line}'
            lines.append(line)

        nspectra = spectra.shape[0]
        for idx in range(nspectra):
            spectrum = spectra[idx][:].tolist()
            spectrum = [str(x) for x in spectrum]
            line = ';'.join(spectrum)
            if prevalues is not None:
                # print(prevalues[idx])
                linep = ';'.join(prevalues[idx])
                line = f'{linep};{line}'
            lines.append(line)
        f1 = open(filecsv, 'w')
        for line in lines:
            f1.write(line)
            f1.write('\n')
        f1.close()
