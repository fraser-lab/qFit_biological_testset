#!/usr/bin/env python

import pandas as pd
import argparse

def parse_log(qfit_log, pre_log, pdb):
    rval = pd.DataFrame()
    log = open(qfit_log, 'r')
    rval.loc[1,'PDB'] = pdb
    for line in log:
        if line.startswith('Final R-work'):
            rval.loc[1,'Rwork_qFit'] = line.split('=')[1][1:6]
            rval.loc[1,'Rfree_qFit'] = line.split('=')[2][1:6]
    log = open(pre_log, 'r')
    for line in log:
        if line.startswith('Final R-work'):
            rval.loc[1,'Rwork_pre'] = line.split('=')[1][1:6]
            rval.loc[1,'Rfree_pre'] = line.split('=')[2][1:6]
    rval.to_csv(pdb + '_rvalues.csv', index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('qFit_Log_File')
    parser.add_argument('pre_Log_File')
    parser.add_argument('PDB')
    args = parser.parse_args()
    parse_log(args.qFit_Log_File, args.pre_Log_File, args.PDB)
