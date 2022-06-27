# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 12:23:47 2022

@author: s139188

Assessment of the quality of sampled peptides, based on (mean) edit distance
and fraction of basic AAs.

PSEUDOCODE
read files -> :lib: list of lists of str
calculate AA fractions of the FT set per peptide, then compare with all peptides
    in epoch 28-3
FT_peptides = [[seq, score], ...]
FT_peptides_sorted = 

"""
from modlamp.analysis import GlobalAnalysis
from modlamp.core import count_aas
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import os
import pickle
# from dataset_comparison_plot import read_store_in_lib, plot_properties
import collections



def read_store_in_lib(files=None):
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
    
    ### Taking out FT_data and put the rest in one big lib
    FT_seqs = lib[-1]
    import itertools
    seqs = list(itertools.chain.from_iterable(lib[:-1]))
    
    return seqs, FT_seqs



def count_aas(seq, scale='relative'):
    """Function to count the amino acids occuring in a given sequence.

    :param seq: {str} amino acid sequence
    :param scale: {'absolute' or 'relative'} defines whether counts or frequencies are given for each AA
    :return: {dict} dictionary with amino acids as keys and their counts in the sequence as values.
    """
    if seq == '':  # error if len(seq) == 0
        seq = ' '
    aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    scl = 1.
    if scale == 'relative':
        scl = len(seq)
    aa = {a: (float(seq.count(a)) / scl) for a in aas}
    aa = collections.OrderedDict(sorted(list(aa.items())))
    return aa

# %%
def rank(seqs, FT_seqs):
    """
    params
    ------
    :seqs: list of lists of str
    :FT_seqs: list of lists of str
    
    returns
    -------
    :seqs_w_scores: list[tups(seq, score)]
    
    
    internal vars
    -------------
    :seqs_w_freqs_all_AA:
    :FT_seqs_w_freqs_all_AA:
    
    :seqs_w_K_and_R:
    :FT_seqs_w_K_and_R:
    
        
    
    Steps:
    ------
        calculate AA fractions of the FT set per peptide, then compare with all 
            peptides in epoch 28-32
        FT_peptides = [[seq, score], ...]
        FT_peptides_sorted = 
    """
    print('\n----------------------------------')
    print('Calculating AA fractions...')
    # aa = count_aas(seq)
    # K, R = aa['K'], aa['R']
    # %%
    # Calculate sum of K + R fractions of the FT set per peptide
    # FT_seqs = ['abcab', 'adfgijafsdgik']
    FT_seqs_freqs_KR = []
    for FT_seq in FT_seqs:
        FT_seqs_w_freqs_all_AA = count_aas(FT_seq)
        FT_seqs_freqs_KR.append((
            FT_seq_freqs_all_AA['K'], FT_seq_freqs_all_AA['R'])
            )
    
    
    FT_seqs_w_KR = list(zip(FT_seqs, *FT_seqs_freqs_KR))  # [(FTseq, fracK, fracR), ...]
    #%%
    # Calculate sum of K + R fractions of the rest per peptide
    gen_seqs_freqs_KR = []
    for gen_seq in gen_seqs:
        gen_seq_freqs_all_AA = count_aas(gen_seq)
        gen_seqs_freqs_KR.append((
            gen_seq_freqs_all_AA['K'], gen_seq_freqs_all_AA['R'])
            )
    
    gen_seqs_w_KR = list(zip(gen_seqs, *gen_seqs_freqs_KR))  # [(genseq, fracK, fracR), ...]
    
    
    
    #%%
    

    def formula(gen_seq_i, FT_seq_i):
        sum(abs(seqs_w_KR[gen_seq_i][1] - FT_seqs_w_KR[FT_seq_i][1]),
            abs(seqs_w_KR[gen_seq_i][2] - FT_seqs_w_KR[FT_seq_i][2])
            )

    
    ft_seqs, gen_seqs = list(), list()
    
    seqs_and_scores = list()
    for gen_seq_i in range(len(gen_seqs)):
    	max_score = -1
    	for FT_seq_i in range(len(FT_seqs)):
    		score = formula(gen_seq_i, FT_seq_i)
    		max_score = max(score, max_score)
            
    	gen_seq = seqs_w_KR[gen_seq_i][0]
    	seqs_and_scores.append((max_score, gen_seq))
    
    seqs_and_scores.sort()
    
    
        
    #     for FT_seqs_freq_KR in FT_seqs_freqs_KR:
    #         for aa in aas:
    #            # Calc abs diff in aafreq between seq&seqFT
    #             aa = count_aas(seq)
    #             K, R = aa['K'], aa['R']
    #         # find smallest diff between seq&seqFT
    # for sublib in lib:
    #     for seq in sublib:
            
    
    print('\n----------------------------------')
    print('Calculating edit distances...')
    from nltk import edit_distance
    
    # Initialize:
    FT_data = lib[-1]
    sublibs_no_FT = lib[:-1]
    edit_dist_per_sublib = [[] for _ in range(len(sublibs_no_FT))]  
    
    # Calculate:
        
    for isublib, sublib in enumerate(self.library):
        # exclude FT data i.e. the last sub-library
        if not isublib == len(self.library)-1: 
            # Track progress:
            print(f'\tSublib {isublib+1}/{len(sublibs_no_FT)} ')
            for seq in sublib:
                for seq_FT in FT_data:
                    # calc edit distances of sequence with FT sequence: 
                    edit_dist_per_sublib[isublib].append(
                        edit_distance(seq, seq_FT))
    eds = edit_dist_per_sublib
    print('Done calculatings edit dists!')
    return FT_peptides_sorted
a = rank()
# def edit_distances(lib):
#     """Calculate the edit mean edit distances between the fine-tuning data-
#     set and all the other data-sets.        
    
#     Params
#     -----
#     lib : list[lists['str']]
    
#     Returns
#     -------
#     eds : ndarray
#         The mean edit distances between the fine-tuning data-
#         set and all the other data-sets, i.e. len = n_datasets - 1.        

#     """
#     print('\n----------------------------------')
#     print('Calculating edit distances...')
#     from nltk import edit_distance
    
#     # Initialize:
#     FT_data = lib[-1]
#     sublibs_no_FT = lib[:-1]
#     edit_dist_per_sublib = [[] for _ in range(len(sublibs_no_FT))]  
    
#     # Calculate:
#     for isublib, sublib in enumerate(self.library):
#         # exclude FT data i.e. the last sub-library
#         if not isublib == len(self.library)-1: 
#             # Track progress:
#             print(f'\tSublib {isublib+1}/{len(sublibs_no_FT)} ')
#             for seq in sublib:
#                 for seq_FT in FT_data:
#                     # calc edit distances of sequence with FT sequence: 
#                     edit_dist_per_sublib[isublib].append(
#                         edit_distance(seq, seq_FT))
#     eds = edit_dist_per_sublib
#     print('Done calculatings edit dists!')
#     return eds



def rank(lib, )


class GlobalAnalysisScoring(GlobalAnalysis):
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
        self.edit_dist_per_sublib = None
    # def get_aa_freqs():
    #     return super().calc_aa_freqs()
    
    
            
            # # save edit distances
            # with open(filename_edit_dists, 'wb') as fp:
            #     pickle.dump(self.edit_dist_minima_per_sublib, fp)
        
    # def plot_summary(self, filename=None, colors=None, plot=True, 
    #          gradient_color='Blues'):
    #     """Method to generate a visual summary of different characteristics of 
    #     the given library. The class methods are used with their standard options.
        
    #     :param filename: {str} path to save the generated plot to.
    #     :param colors: {str / list} color or list of colors to use for plotting. e.g. '#4E395D', 'red', 'k'
    #     :param plot: {boolean} whether the plot should be created or just the features are calculated
    #     :param gradient_color: {str} Matplotlib Colormap, e.g. 'Blues', 'Greens' (other options at https://matplotlib.org/stable/tutorials/colors/colormaps.html#sequential)
    #     :return: visual summary (plot) of the library characteristics (if ``plot=True``).
    #     :Example:
        
    #     >>> g = GlobalAnalysis([seqs1, seqs2, seqs3])  # seqs being lists / arrays of sequences
    #     >>> g.plot_summary()
        
    #     .. image:: ../docs/static/summary.png
    #     :height: 600px
    #     """
    #     # calculate all global properties
    #     self.calc_len()
    #     self.calc_aa_freq(plot=False)
    #     self.calc_charge(ph=7.4, amide=True)
    #     self.calc_H()
    #     self.calc_uH()
        
        if plot:
           
            # =================================================================
            # Plot settings ###################################################
            # =================================================================
            
            plt.figure(figsize=(210/10, 297/10), dpi=30)  # A4=(210,297)
            fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(25, 35), 
                gridspec_kw={'height_ratios': [1,1]})
            (ax1, ax2) = axes
            fontsize_titles, fontsize_xy_labels, ticklabel_size = 36., 28., 22.
            plt.suptitle('Sampled Peptides Scores', fontweight='bold', 
                             fontsize=fontsize_titles)
                
            labels = self.libnames
            num = len(labels)
                
            # Create color gradient; with black for FT data set:
            colors = plt.get_cmap(gradient_color)(np.linspace(0.375, 1, num))
            black = np.array([[0., 0., 0., 1.]])
            colors[-1] = black
               
            for a in [ax1, ax2]:
                # only left and bottom axes, no box
                # a.spines['right'].set_visible(False)
                a.spines['top'].set_visible(False)
                a.xaxis.set_ticks_position('bottom')
                a.yaxis.set_ticks_position('left')
                a.tick_params(labelright=True, right=True, 
                              labelsize=ticklabel_size)
            
            # =================================================================
            # 1 edit distance #################################################
            # =================================================================
            print("Creating plot 1 (edit distance)...", end=' ')
            # Calculate min edit dists
            self.edit_distances()
            box = ax1.boxplot(self.edit_dist_per_sublib, notch=False, vert=1,
                              patch_artist=True)
            plt.setp(box['whiskers'], color='black')
            plt.setp(box['medians'], linestyle='-', linewidth=1.5, 
                     color='black')
            for p, patch in enumerate(box['boxes']):
                patch.set(facecolor=colors[p], edgecolor='black', alpha=0.8)
            ax1.set_ylabel('Minimum Edit Distance', fontweight='bold', 
                            fontsize=fontsize_xy_labels)
            ax1.set_xticks([x + 1 for x in range(len(labels)-1)])
            if num > 30:
                ax1.set_xticklabels(labels[:-1], fontweight='bold', rotation=90)
            else:
                ax1.set_xticklabels(labels[:-1], fontweight='bold')
            # ax6.set_xlabel('Fine-Tuning Epoch', fontweight='bold', 
            #                 fontsize=fontsize_xy_labels)            
            
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
            
            AAis = [7, 9, 15]
            AAs_keys = ['H', 'K', 'R']
            self.aafreqpos = [self.aafreq[AAi] for AAi in AAis]
            
            if self.shapes:  # if the library consists of different sized sub libraries
                reversed_colors = [c for c in reversed(colors)]
                for i, l in enumerate(reversed(self.aafreqpos)):
                    for a in range(20):
                        ax2.bar(a - offsets[i], l[a], w, 
                                color=reversed_colors[i], alpha=0.8)
            else:
                for a in AAis:  # H, K, R
                    ax2.bar(a, self.aafreqpos[0][a], w, color=colors[0], alpha=0.8)
            ax2.set_xlim([-1., 20.])
            ax2.set_ylim([0, 1.05 * np.max(self.aafreqpos)])
            ax2.set_xticks(range(20))
            ax2.set_xticklabels(AAs_keys, fontweight='bold')
            ax2.set_ylabel('Fraction', fontweight='bold', 
                           fontsize=fontsize_xy_labels)
            ax2.set_xlabel('Basic Amino Acids', fontweight='bold', 
                           fontsize=fontsize_xy_labels)
            print("Done!")
            
            # =================================================================
            # More figure settings & saving ###################################
            # =================================================================
            fig.tight_layout()  # to prevent overlap of subplots
            if filename:
                print(f"Saving plot as {filename}...", end=" ")
                plt.savefig(filename, dpi=200)
                print("Done!")
                print('File location:')
                print(f'/{filename}')
            else:
                plt.show()
               








if __name__ == "__main__":
    # ==========================================================================
    # Parameters
    # ==========================================================================
    output_filename = 'Scoring_Sampled_Peptides'
    # Add version number based on existence of previous version
    i = 1
    while os.path.exists(f'./{output_filename}_v{i}.png'):
        i += 1
    output_filename = output_filename + '_v' + str(i) + '.png'
    print('File will be saved as:', output_filename)
            
    # n_epochs_FT = 80
    epochs_ranked = range(28, 33)
    # n_intermediate_plots = 40  # Best to leave at <=10 because of width
    
    # Total figure size, and label font sizes in plot_summary
    
    # step = int(n_epochs_FT / n_intermediate_plots)
    # epochs_plotted = range(step, n_epochs_FT+1, step)
    

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
        [fr'C:\Users\s139188\OneDrive - TU Eindhoven\Documents\01 TUe\.BMT 8\Q4 BEP\LSTM_peptides\0 finetune_and_sample_TRUE lr 0.005 valsplit 0\finetune_and_sample_1000_epoch_nr_{epoch}\sampled_sequences_temp1.25.csv',
         f'{epoch:02d}']
        for epoch in epochs_ranked
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
        )
    print("\nReading files done!\n")
    # Plot properties
    print("Plotting...\n===============")
    rank(lib, 
        names=[path_and_name_list[1] for path_and_name_list in filepaths_and_names],  # names
        output_filename=output_filename
        )