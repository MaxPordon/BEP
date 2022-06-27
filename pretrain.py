# -*- coding: utf-8 -*-

import os

os.system("python LSTM_peptides.py"
          " --mode pretrain"
          " --name mueller167"
          " --dataset training_sequences_noC.csv"
          " --layers 2"
          " --neurons 256"
          " --epochs 167"
          " --refs ''"   # '' instead of False because of argparse's str interpretation
          " --sample 1000"
          " --cv 5"
          )