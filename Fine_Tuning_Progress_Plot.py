# -*- coding: utf-8 -*-
"""
@author: Pordon, W.M.J.
Created on Thu Jun  2 22:27:04 2022

Make plots of properties of pre-training data-set,
                            sampled data-set w/o fine-tuning (~pre-training),
                            sampled data-set w/ 1 epoch of fine-tuning,
                            ...
                            sampled data-set w/ 100 epochs of fine-tuning.
                            fine-tuning data-set
                            
Peusdocode:
----------
import dataset_comparison_plot functions
state files and names vars:
    
run dataset_comparison_plot functions                            


"""
from modlamp.analysis import GlobalAnalysis
from modlamp.core import count_aas
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
# from colour import Color
import numpy as np
import os
# import pickle
# from dataset_comparison_plot import read_store_in_lib, plot_properties


class GlobalAnalysisFineTuningProcess(GlobalAnalysis):
    """Same as GlobalAnalysis, but with adapted plot_summary method to allow
    more plots and more data-sets in the plots. Also, the method edit_distances 
    added.
    
    Methods
    -------
    edit_distances
        Calculate the edit mean edit distances between the fine-tuning data-
        set and all the other data-sets.   
    plot_summary
        Method to generate a visual summary of different characteristics of 
        the given library. The class methods are used with their standard
        options.
    """
    def __init__(self, library, names):
    #     GlobalAnalysis.__init__(self, *args)
        super().__init__(library, names)
        self.edit_dist_minima_per_sublib = None
    
    
    def edit_distances(self, filename_edit_dists=r'C:\Users\s139188\OneDrive - TU Eindhoven\Documents\01 TUe\.BMT 8\Q4 BEP\LSTM_peptides\edit_distances_minimum'):
        """Calculate the edit mean edit distances between the fine-tuning data-
        set and all the other data-sets.        

        Returns
        -------
        self.mean_eds : ndarray
            The mean edit distances between the fine-tuning data-
            set and all the other data-sets, i.e. len = n_datasets - 1.        

        """
        # Append the seq length range to the file name
        # filename_edit_dists = filename_edit_dists + \
        #     f'seqlengthrange_{seq_length_min}_{seq_length_max}' + '.ob'
            
        # # If file exists, don't compute everything anew
        # if os.path.isfile(filename_edit_dists):
        #     print(f'\n\t{filename_edit_dists} found, reading...', end=' ')
        #     with open(f'{filename_edit_dists}', 'rb') as inf:
        #         self.edit_dist_minima_per_sublib = pickle.load(inf)
        #     # self.edit_dist_minima_per_sublib = np.loadtxt(f'./{filename_edit_dists}')
        #     print('Reading min edit dists file done!')
            
        # # Else: compute edit distances
        # else:  
        #     print(f"\n\tFile '{filename_edit_dists}' not found! -->")
        print()
        print('----------------------------------')
        print('Calculating edit distances...')
        from nltk import edit_distance
        # library = self.library
        # max_sublib_length = len(max(self.library, key=len))
        # print(f'{max_sublib_length = }')
        
        # Initialize:
        FT_data = self.library[-1]
        sublibs_no_FT = self.library[:-1]
        # eventual shape (len(self.library), (len(sublibrary_i))):
        edit_dist_minima_per_sublib = [
            [] for _ in range(len(sublibs_no_FT))]  
        # edit_dist_min_mean = []  # shape (len(self.library)-1 )
        for isublib, sublib in enumerate(self.library):
            # exclude FT data i.e. the last sub-library
            if not isublib == len(self.library)-1: 
                # Track progress:
                print(f'\tSublib {isublib+1}/{len(sublibs_no_FT)} ')
                #       f'(len {len(sublibs_no_FT[isublib])})')
                # edit_dist_minima_per_sublib.append(list())
                for seq in sublib:
                    seq_edit_dist_lst = []
                    for seq_FT in FT_data:
                        seq_edit_dist_lst.append(edit_distance(seq, seq_FT))
                    # calc minimum edit distance of sequence with FT sequence    
                    edit_dist_minima_per_sublib[isublib].append(
                        min(seq_edit_dist_lst))  # len = #seqs in sublib
        print('Done calculatings edit dists!')
        self.edit_dist_minima_per_sublib = edit_dist_minima_per_sublib
            
            # # save edit distances
            # with open(filename_edit_dists, 'wb') as fp:
            #     pickle.dump(self.edit_dist_minima_per_sublib, fp)
            
    
    
    def plot_summary(self, filename=None, colors=None, plot=True, 
                     gradient_color='Blues'):
        """Method to generate a visual summary of different characteristics of 
        the given library. The class methods are used with their standard options.
        
        :param filename: {str} path to save the generated plot to.
        :param colors: {str / list} color or list of colors to use for plotting. e.g. '#4E395D', 'red', 'k'
        :param plot: {boolean} whether the plot should be created or just the features are calculated
        :param gradient_color: {str} Matplotlib Colormap, e.g. 'Blues', 'Greens' (other options at https://matplotlib.org/stable/tutorials/colors/colormaps.html#sequential)
        :return: visual summary (plot) of the library characteristics (if ``plot=True``).
        :Example:
        
        >>> g = GlobalAnalysis([seqs1, seqs2, seqs3])  # seqs being lists / arrays of sequences
        >>> g.plot_summary()
        
        .. image:: ../docs/static/summary.png
            :height: 600px
        """
        # calculate all global properties
        self.calc_len()
        self.calc_aa_freq(plot=False)
        self.calc_charge(ph=7.4, amide=True)
        self.calc_H()
        self.calc_uH()
        
        if plot:
           
            # plot settings
            plt.figure(figsize=(210/10, 297/10), dpi=30)  # A4=(210,297)
            fig, axes = plt.subplots(nrows=6, ncols=1, figsize=(25, 35), 
                gridspec_kw={'height_ratios': [2,1,1,1,1,2]})
            (ax1, ax2, ax3, ax4, ax5, ax6) = axes
            fontsize_titles, fontsize_xy_labels, ticklabel_size = 36., 28., 22.
            # Put seq-length range in title if it is not 0 to infinity:
            if not (seq_length_min == 0 and seq_length_max == 0):
                plt.suptitle(
                    'Fine-Tuning Progress\n'
                    f'seq-length range: {seq_length_min} - {seq_length_max}',
                    fontweight='bold', fontsize=fontsize_titles
                    )
            else: 
                plt.suptitle('Fine-Tuning Progress', fontweight='bold', 
                             fontsize=fontsize_titles)
                
            labels = self.libnames
            num = len(labels)
                
            # Create color gradient; with black for FT data set:
            colors = plt.get_cmap(gradient_color)(np.linspace(0.375, 1, num))
            black = np.array([[0., 0., 0., 1.]])
            colors[-1] = black
               
            for a in [ax1, ax2, ax3, ax4, ax5, ax6]:
                # only left and bottom axes, no box
                # a.spines['right'].set_visible(False)
                a.spines['top'].set_visible(False)
                a.xaxis.set_ticks_position('bottom')
                a.yaxis.set_ticks_position('left')
                a.tick_params(labelright=True, right=True, 
                              labelsize=ticklabel_size)

               
            # =================================================================
            # 1 length box plot ###############################################
            # =================================================================
            print("Creating plot 1 (length box plot)...", end=' ')
            box = ax1.boxplot(self.len, notch=0, vert=1, patch_artist=True)
            plt.setp(box['whiskers'], color='black')
            plt.setp(box['medians'], linestyle='-', linewidth=1.5, 
                     color='black')
            for p, patch in enumerate(box['boxes']):
                patch.set(facecolor=colors[p], edgecolor='black', alpha=0.8)
            ax1.set_ylabel('Sequence Length', fontweight='bold', fontsize=fontsize_xy_labels)
            ax1.set_xticks([x + 1 for x in range(len(labels))])
            if num > 30:
                ax1.set_xticklabels(labels, fontweight='bold', rotation=90)
            else:
                ax1.set_xticklabels(labels, fontweight='bold')  # rotation=90)
            # ax1.set_xlabel('Fine-Tuning Epoch', fontweight='bold', 
            #                fontsize=fontsize_xy_labels)    
            print("Done!")
            
               
            # =================================================================
            # 2 AA bar plot ###################################################
            # =================================================================
            print("Creating plot 2 (AA bar plot)...", end=' ')
            d_aa = count_aas('')
            hands = [
                mpatches.Patch(label=labels[i], facecolor=colors[i], alpha=0.8) 
                for i in range(num)
                ]
            w = .9 / num  # bar width
            offsets = np.arange(start=-w, step=w, stop=num * w)  # bar offsets if many libs
            if self.shapes:  # if the library consists of different sized sub libraries
                reversed_colors = [c for c in reversed(colors)]
                for i, l in enumerate(reversed(self.aafreq)):
                    for a in range(20):
                        ax2.bar(a - offsets[i], l[a], w, 
                                color=reversed_colors[i], alpha=0.8)
            else:
                for a in range(20):
                    ax2.bar(a, self.aafreq[0][a], w, color=colors[0], alpha=0.8)
            ax2.set_xlim([-1., 20.])
            ax2.set_ylim([0, 1.05 * np.max(self.aafreq)])
            ax2.set_xticks(range(20))
            ax2.set_xticklabels(d_aa.keys(), fontweight='bold')
            ax2.set_ylabel('Fraction', fontweight='bold', fontsize=fontsize_xy_labels)
            ax2.set_xlabel('Amino Acids', fontweight='bold', fontsize=fontsize_xy_labels)
            # Legend in 2nd subplot ####
            # ax2.legend(handles=hands, labels=labels, 
            #             ncol=(len(lib)//4 + len(lib)%4),  # 4 per column
            #             # bbox_to_anchor=(0.6, 1.8),
            #             loc = 'upper right'
            #             )
            # Legend above 1st subplot ####
            # ax1.legend(handles=hands, labels=labels, 
            #             ncol=(len(lib)//4 + len(lib)%4),  # 4 per column
            #             bbox_to_anchor=(0.6, 2.5), prop={'size': 3}
            #             )
            ### Pseudocode ###
            # legend: horizontally:
            # [lightestblue] PT = Pre-Training Data Set
            # [lightestblue+1] PT-S = Post-Training Sampled
            # [lightblue+2] - [black-1] Sampled After Fine-Tuning Epoch # 
            # [black] FT = Fine-Tuning Data Set
            ###
            # lightestblue = colors[0]        &  'PT = Pre-Training Data Set\t'
            # lightestblue_plus1 = colors[1]  &  'PT-S = Post-Training Sampled\t'
            # lightestblue_plus2 = colors[2]  &  u'\u22C5\u22C5\u22C5'  # or if that doesnt work: '-' 
            # black_min1 = colors[-2]         &  'Sampled After Fine-Tuning Epoch #\t'
            # black = colors[-1]              &  'FT = Fine-Tuning Data Set'
            # ax1.legend(handles=)
            
            # fig.legend(handles, ['oral', 'physa'], bbox_to_anchor=(2, 0),loc = 'lower right')
            # ax2.legend(bbox_to_anchor=(0, 1, 1, 0), loc="lower left", mode="expand", ncol=2)
            print("Done!")
            
            # =================================================================
            # 3 hydophobicity violin plot #####################################
            # =================================================================
            print("Creating plot 3 (hydrophobicity violin plot)...", end=' ')
            for i, l in enumerate(self.H):
                vplot = ax3.violinplot(l, positions=[i + 1], widths=0.5, 
                                       showmeans=True, showmedians=False)
                # crappy adaptions of violin dictionary elements
                vplot['cbars'].set_edgecolor('black')
                vplot['cmins'].set_edgecolor('black')
                vplot['cmeans'].set_edgecolor('black')
                vplot['cmaxes'].set_edgecolor('black')
                vplot['cmeans'].set_linestyle('--')
                for pc in vplot['bodies']:
                    pc.set_facecolor(colors[i])
                    pc.set_alpha(0.8)
                    pc.set_edgecolor('black')
                    pc.set_linewidth(1.5)
                    pc.set_alpha(0.7)
                    pc.set_label(labels[i])
            ax3.set_xticks([x + 1 for x in range(len(labels))])
            if num > 30:
                ax3.set_xticklabels(labels, fontweight='bold', rotation=90)
            else:
                ax3.set_xticklabels(labels, fontweight='bold')
            ax3.set_ylabel('Global \nHydrophobicity', fontweight='bold', 
                           fontsize=fontsize_xy_labels)
            # ax3.set_xlabel('Fine-Tuning Epoch', fontweight='bold', 
            #                fontsize=fontsize_xy_labels)    
            print("Done!")
            
            # =================================================================   
            # 4 hydrophobic moment violin plot ################################
            # =================================================================
            print("Creating plot 4 (hydrophobic moment violin plot)...",
                  end=' ')
            for i, l in enumerate(self.uH):
                vplot = ax4.violinplot(
                    l, positions=[i + 1], widths=0.5, showmeans=True, showmedians=False)
                # crappy adaptions of violin dictionary elements
                vplot['cbars'].set_edgecolor('black')
                vplot['cmins'].set_edgecolor('black')
                vplot['cmeans'].set_edgecolor('black')
                vplot['cmaxes'].set_edgecolor('black')
                vplot['cmeans'].set_linestyle('--')
                for pc in vplot['bodies']:
                    pc.set_facecolor(colors[i])
                    pc.set_alpha(0.8)
                    pc.set_edgecolor('black')
                    pc.set_linewidth(1.5)
                    pc.set_alpha(0.7)
                    pc.set_label(labels[i])
            ax4.set_xticks([x + 1 for x in range(len(labels))])
            if num > 30:
                ax4.set_xticklabels(labels, fontweight='bold', rotation=90)
            else:
                ax4.set_xticklabels(labels, fontweight='bold')
            ax4.set_ylabel('Global Hydro-\nphobic Moment', fontweight='bold', fontsize=fontsize_xy_labels)
            # ax4.set_xlabel('Fine-Tuning Epoch', fontweight='bold', 
                           # fontsize=fontsize_xy_labels)    
            print("Done!")
            
            # =================================================================
            # 5 charge histogram ##############################################
            # =================================================================
            print("Creating plot 5 (charge histogram)...", end=' ')
            if self.shapes:  # if the library consists of different sized sub libraries
                bwidth = 1. / len(self.shapes)
                for i, c in enumerate(self.charge):
                    counts, bins = np.histogram(c, range=[-5, 20], bins=25)
                    ax5.bar(bins[1:] + i * bwidth, counts / np.max(counts), 
                            bwidth, color=colors[i], label=labels[i], alpha=0.8)
            else:
                ax5.hist(
                    self.charge / np.max(self.charge), 25, alpha=0.8, 
                    align='left', rwidth=0.95, histtype='bar', label=labels,
                    color=colors[:num])
            ax5.set_xlabel('Global Charge $\it{(pH: 7.4 ,  amide: true)}$', 
                           fontweight='bold', fontsize=fontsize_xy_labels)
            # ax5.set_xlabel('pH: 7.4 ,  \namide: true', 
            #                fontsize=fontsize_xy_labels)
            ax5.set_ylabel('Fraction', fontweight='bold', fontsize=fontsize_xy_labels)
            # ax5.title.set_text('pH: 7.4 ,  \namide: true')
            # ax5.legend(loc='best')
            print("Done!")
            
            # =================================================================
            # 6 edit distance #################################################
            # =================================================================
            print("Creating plot 6 (edit distance)...", end=' ')
            # Calculate min edit dists
            self.edit_distances()
            box = ax6.boxplot(self.edit_dist_minima_per_sublib, notch=False, vert=1, 
                               patch_artist=True)
            plt.setp(box['whiskers'], color='black')
            plt.setp(box['medians'], linestyle='-', linewidth=1.5, 
                     color='black')
            for p, patch in enumerate(box['boxes']):
                patch.set(facecolor=colors[p], edgecolor='black', alpha=0.8)
            ax6.set_ylabel('Minimum Edit Distance', fontweight='bold', 
                            fontsize=fontsize_xy_labels)
            ax6.set_xticks([x + 1 for x in range(len(labels)-1)])
            if num > 30:
                ax6.set_xticklabels(labels[:-1], fontweight='bold', rotation=90)
            else:
                ax6.set_xticklabels(labels[:-1], fontweight='bold')
            # ax6.set_xlabel('Fine-Tuning Epoch', fontweight='bold', 
            #                 fontsize=fontsize_xy_labels)            
            
            print("Done!")
            
            
            # More figure settings & saving
            fig.tight_layout()  # to prevent overlap of subplots
            if filename:
                print(f"Saving plot as {filename}...", end=" ")
                plt.savefig(filename, dpi=200)
                print("Done!")
                print('File location:')
                print(f'/{filename}')
            else:
                plt.show()
            
            
            
            # 7 3D plot ###############################################
            # print("\nCreating plot 6 (3D plot)...", end=" ")
            # # plot settings
            # fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(10,10))
            # (ax6) = axes
            # plt.suptitle('3D-plot', fontweight='bold', fontsize=fontsize_titles)
            # labels = self.libnames
            # if not colors:
            #     colors = ['#FA6900', '#69D2E7', '#542437', '#53777A', '#CCFC8E', '#9CC4E4']
            # num = len(labels)
            # # plot           
            # ax7.spines['left'].set_visible(False)
            # ax7.spines['bottom'].set_visible(False)
            # ax7.set_xticks([])
            # ax7.set_yticks([])
            # ax7 = fig.add_subplot(1, 1, 1, projection='3d')
            # for i, l in enumerate(range(num)):
            #     xt = self.H[l]  # find all values in x for the given target
            #     yt = self.charge[l]  # find all values in y for the given target
            #     zt = self.uH[l]  # find all values in y for the given target
            #     ax7.scatter(xt, yt, zt, c=colors[l], alpha=.8, s=25, label=labels[i])
            
            # ax7.set_xlabel('H', fontweight='bold', fontsize=fontsize_xy_labels)
            # ax7.set_ylabel('Charge', fontweight='bold', fontsize=fontsize_xy_labels)
            # ax7.set_zlabel('uH', fontweight='bold', fontsize=fontsize_xy_labels)
            # data_c = [item for sublist in self.charge for item in sublist]  # flatten charge data into one list
            # data_H = [item for sublist in self.H for item in sublist]  # flatten H data into one list
            # data_uH = [item for sublist in self.uH for item in sublist]  # flatten uH data into one list
            # ax7.set_xlim([np.min(data_H), np.max(data_H)])
            # ax7.set_ylim([np.min(data_c), np.max(data_c)])
            # ax7.set_zlim([np.min(data_uH), np.max(data_uH)])
            # # ax7.legend(loc='best')
            
            # print("Done!")
            # if filename:
            #     plt.savefig('3D_plot_properties_comparison', dpi=200)
            # else:
            #     plt.show()
               

def read_store_in_lib(files=None, seq_length_min=0, seq_length_max=0):
    """Read sequence files, preprocess and store in ndarray.
    
    Parameters
    ----------
    files : list[str]
        Filenames/paths in a list.
    seq_length_min, seq_length_max : int, default=0
        Minimum and maximum sequence length to be processed; remove all
        sequences that are shorter than the min, and all that are longer than
        the max. If seq_length_max = 0, no maximum length is used.
    
    Returns
    -------
    lib : np.ndarray
        NumPy ndarray where each row is a sub-library
    """
    from LSTM_peptides import SequenceHandler
    data = SequenceHandler()
    
    ### Read files, store in lib
    lib = []
    for file in files:
        data.load_sequences(file)
        lib.append(data.sequences)  
    
    print(f'n of sublibs: {len(lib)}')
    print('Lengths of sublibs:', end=' ')
    for sublib in lib:
        print(len(sublib), end=' ')
    
    ### Preprocessing
    # Remove the sequence(s) with the 'B' - deleting by index not possible
    # because sequences are scrambled by SequenceHandler.
    lib[:] = [[seq for seq in sublib if 'B' not in seq] for sublib in lib]
    print("\nSequences with 'B's removed!")
    print("Lengths of sublibs *AFTER REMOVING 'B's*:")
    for sublib in lib:
        print(len(sublib), end=', ')
    
    
    # Filter based on length
    if not seq_length_min == 0 and not seq_length_max == 0:
        print('\nFiltering based on min and max seq length', seq_length_min, 
              'and', seq_length_max, '...', end=' ')
        lib_filtered = []
        lib_filtered[:] = [[seq for seq in sublib if len(seq) > seq_length_min 
                            and len(seq) <= seq_length_max] 
                           for sublib in lib
                           # exclude training, non-FT sampled and FT data-set 
                           # from filtering:
                           if not sublib in (lib[0], lib[1], lib[-1])]  
        # add the training, non-FT sampled and FT data-set again
        lib_filtered.insert(0, lib[1])
        lib_filtered.insert(0, lib[0])
        lib_filtered.append(lib[-1])
        print('Done filtering!')
        print("\nLengths of sublibs *AFTER FILTERING BASED OF LENGTH* (and removing 'B's):")
        for sublib in lib_filtered:
            print(len(sublib), end=', ')
    
        return lib_filtered

    else:
        return lib


def plot_properties2(lib, names, output_filename, colors):
    """Make plots using an adapted version of
    modlamp.analysis.GlobalAnalysis.plot_summary(). This
    version plots 6 subplots (original 6 - 3D plot + min edit distance) as 6x1 
    instead of 6 subplots as 3x2.
    
    Parameters
    ----------
    lib : list[str]
        The sequences. CANNOT be a NumPy ndarray.
    names : list[str]
        The names for the plots.
    output_filename : str
        Output file with be saved with this name.
    colors : str / list
        color or list of colors to use for plotting. e.g. '#4E395D', 'red', 'k'
    
    """
    # from modlamp.analysis import GlobalAnalysis
    analysis = GlobalAnalysisFineTuningProcess(library=lib, names=names)
    analysis.plot_summary(filename=output_filename, colors=colors)





if __name__ == "__main__":
    # ==========================================================================
    # Parameters
    # ==========================================================================
    seq_length_min, seq_length_max = 0, 0   # set to 0, 0 for no min or max
    print('Sequence length range:', seq_length_min, '-', seq_length_max)
    
    output_filename = 'Fine-Tuning_Process'
    
    if not (seq_length_min == 0 and seq_length_max == 0):
        output_filename = output_filename + \
            '__seq_length_range_{seq_length_min}-{seq_length_max}'
        
    # Add version number based on existence of previous version
    i = 1
    while os.path.exists(f'./{output_filename}_v{i}.png'):
        i += 1
    output_filename = output_filename + '_v' + str(i) + '.png'
    print('File will be saved as:', output_filename)
            
    n_epochs_FT = 80
    n_intermediate_plots = 40  # Best to leave at <=10 because of width
    
    # Total figure size, and label font sizes in plot_summary
    
    step = int(n_epochs_FT / n_intermediate_plots)
    epochs_plotted = range(step, n_epochs_FT+1, step)

    pre_trained_model_path = r'C:\Users\s139188\OneDrive - TU Eindhoven\Documents\01 TUe\.BMT 8\Q4 BEP\LSTM_peptides\mueller167v2\checkpoint\model_epoch_166.hdf5'
    # D1 = Pre-Training dataset from Mueller et al. (n=1554)
    pre_training_data_set = r'C:\Users\s139188\OneDrive - TU Eindhoven\Documents\01 TUe\.BMT 8\Q4 BEP\LSTM_peptides\training_sequences_noC.csv'
    # P1 = sampled peptides from pre-trained model (i.e. on the pre-training dataset) (n=~1000)
    sampled_peptides_pretrained = r'C:\Users\s139188\OneDrive - TU Eindhoven\Documents\01 TUe\.BMT 8\Q4 BEP\LSTM_peptides\mueller167v2\sampled_sequences_temp1.25.csv'
    # D2 = FT data from Portugal (n=8)
    finetuning_data_set = 'C:/Users/s139188/OneDrive - TU Eindhoven/Documents/01 TUe/.BMT 8/Q4 BEP/LSTM_peptides/Finetuning_dataset_8_peptides.csv'
    # P2 = pre-trained & fine-tuned samples (n=1000)
    # P2_path = [[f'./finetune_and_sample_100_epoch_nr_{epoch_FT}'] for epoch_FT in n_epochs_FT]

    ### List of lists of Filepaths and corresponding names
    filepaths_and_names = [
        [fr'C:\Users\s139188\OneDrive - TU Eindhoven\Documents\01 TUe\.BMT 8\Q4 BEP\LSTM_peptides\0 finetune_and_sample_TRUE lr 0.005 valsplit 0\finetune_and_sample_1000_epoch_nr_{epoch_plotted}\sampled_sequences_temp1.25.csv',
         f'{epoch_plotted:02d}']
        for epoch_plotted in epochs_plotted
        ]  # FT'd
    filepaths_and_names.insert(0, [pre_training_data_set, 'PT-S'])  # Non-FT
    filepaths_and_names.insert(0, [sampled_peptides_pretrained, 'PT'])  # Pre-training
    filepaths_and_names.append([finetuning_data_set, 'FT'])  # FT
    # filepaths = [path_and_name_list[0] for path_and_name_list in filepaths_and_names]
    # print(filepaths_and_names)

    # =========================================================================
    # Running
    # =========================================================================
    # Read files and store in lib
    print()
    print("Reading files...")
    print('----------------')
    lib = read_store_in_lib(
        files=[path_and_name_list[0] for path_and_name_list in filepaths_and_names]  # filepaths
        , seq_length_min=seq_length_min, seq_length_max=seq_length_max)
    print("Reading files done!\n")
    # Plot properties
    print("Plotting...\n===============")
    plot_properties2(
        lib, 
        names=[path_and_name_list[1] for path_and_name_list in filepaths_and_names],  # names
        output_filename=output_filename, colors=None
        )