# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 19:29:54 2022

@author: s139188

Find which sequences from the sampled sequences are already present in the 
thpdb: the database containing the FDA-approved protein and peptide 
therapeutics.
"""

import pandas as pd
import os

def read_db(db_file=r'C:\Users\s139188\OneDrive - TU Eindhoven\Documents\01 TUe\.BMT 8\Q4 BEP\thpdb.csv'):
    """Read database"""
    df_db = pd.read_csv(db_file, usecols=['Peptide Sequence'], encoding = "ISO-8859-1")
    print(f'Database file read, containing {len(df_db)} seqs.')
    return df_db


def read_seqs(sampled_seqs_file):
    """Read sampled sequences"""
    df_seqs = pd.read_csv(sampled_seqs_file)
    print(f'Sequences file read, containing {len(df_sampled)} seqs.')
    return df_seqs


def _save_common_seqs(path, filename, counter, common_seqs):
    # Save common sequences in the folder of the sampled sequences
    with open(path, 'w') as f:
        if counter == 0:
            f.write('No common sequences found!')
        else:
            f.write('\n')
            f.write(f'{counter} common sequences found!')
            f.write(f'i.e. {counter/len(df_sampled)*100}% of the sampled'
                    'sequences are in the THPDB.')
        f.write('\n')
        f.write("Common sequences:\n=======================\n")
        for common_seq in common_seqs:
            f.write(common_seq)
            f.write('\n')
    

def find_common_seqs_db(df_db, df_seqs, sampled_seqs_file):
    """
    Find common sequences between thpdb and another csv containing peptides.
    Difference with find_common_seqs_auto is that this function does not look
    for *equal* seqs but for seqs that are partly equal
    """
    print()
    print('=======================')
    print("Common sequences:")
    print('=======================')
    counter = 0
    common_seqs = []
    for sampled_seq in df_seqs:
        for db_seq in df_db:
            if sampled_seq in db_seq:
                print('sampled seq:', sampled_seq)
                print('db seq:     ', db_seq)
                same_pep = input('Is this the same peptide? [y/n]')
                if same_pep in ['y', 'yes', 'Yes', 'Y']:
                    print('-> understood input "yes". Saving...', end=' ')
                    counter += 1
                    common_seqs.append(sampled_seq)
                    print('Done!')
                elif same_pep in ['n', 'no', 'No', 'NO']:
                    print('-> understood input "no".')
    if counter == 0:
        print('No common sequences found!')
    else:
        print()
        print(f'{counter} common sequences found!')
        print(f'i.e. {counter/len(df_seqs)*100}% of the sampled sequences'
              'are in the THPDB.')
    
    # Save common seqs
    _save_common_seqs(path=os.path.dirname(sampled_seqs_file),
                      filename='seqs_in_thpdb.txt',
                      counter=counter, common_seqs=common_seqs)


def find_common_seqs_auto(df_seqs_1, df_seqs_2, sampled_seqs_file, 
                          dataset_name):
    """
    Find common sequences between 2 sequences csvs.
    :dataset_name: name of the dataset which the generated peptides are compared to. This will be appended to the name of the saved file.
    """
    common_seqs = []
    counter = 0
    for s in df_seqs_1:
        for s2 in df_seqs_2:
            if s == s2:
                counter += 1
                common_seqs.append(s)
    _save_common_seqs(path=os.path.dirname(sampled_seqs_file), 
                      filename='common_seqs_' + dataset_name,
                      counter=counter, common_seqs=common_seqs)
    

def main(db_file, sampled_seqs_file):
    df_db = read_db(db_file=r'C:\Users\s139188\OneDrive - TU Eindhoven\Documents\01 TUe\.BMT 8\Q4 BEP\thpdb.csv')
    df_sampled = read_sampled_seqs(sampled_seqs_file=r'C:\Users\s139188\OneDrive - TU Eindhoven\Documents\01 TUe\.BMT 8\Q4 BEP\LSTM_peptides\0 finetune_and_sample_TRUE lr 0.005 valsplit 0\finetune_and_sample_1000_epoch_nr_80\sampled_sequences_temp1.25.csv')
    find_common_seqs(df_db, df_sampled, sampled_seqs_file)
    
if __name__ == "__main__":
    main(db_file=r'C:\Users\s139188\OneDrive - TU Eindhoven\Documents\01 TUe\.BMT 8\Q4 BEP\thpdb.csv',
         sampled_seqs_file = r'C:\Users\s139188\OneDrive - TU Eindhoven\Documents\01 TUe\.BMT 8\Q4 BEP\LSTM_peptides\0 finetune_and_sample_TRUE lr 0.005 valsplit 0\finetune_and_sample_1000_epoch_nr_80\sampled_sequences_temp1.25.csv')
            
