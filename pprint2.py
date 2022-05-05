##############################################################################
# PPRInt2.0 is developed for predicting the RNA-interacting residues in protein  #
# using primary structure. It is developed by Prof G. P. S. Raghava's group.  #
# Please cite: PPRInt2.0						      #
# ############################################################################
import os
import pickle
import numpy as np
import pandas as pd
import sys
import csv
import re
import glob
import time
from time import sleep
#from tqdm import tqdm
import argparse
import warnings
from argparse import RawTextHelpFormatter
import tensorflow as tf
from keras.models import load_model
import uuid
warnings.filterwarnings('ignore')


parser = argparse.ArgumentParser(description='Please provide following arguments',formatter_class=RawTextHelpFormatter)
## Read Arguments from command
parser.add_argument("-i", "--input", type=str, required=True, help="Input: Protein sequences in FASTA format")
parser.add_argument("-o", "--output",type=str, help="Output: File for saving results by default outfile.csv")
parser.add_argument("-j", "--job",type=int, choices = [1,2,3,4], help="Feature Type:\n1: Amino Acid Binary (AAB) profile\n2: Physico-chemical based binary (PCB) profile\n3: Position-Specific Scoring Matrix (PSSM) profile\n4: Hybrid (AAB+PCB+PSSM), by default 3 (PSSM)")
parser.add_argument("-t","--threshold", type=float, help="Threshold: Value between 0 to 1 by default 0.5")
parser.add_argument("-p", "--path",type=str, help="Path: Please provide the path of python which has all libraries installed")
args = parser.parse_args()

def readseq(file):
    with open(file) as f:
        records = f.read()
    records = records.split('>')[1:]
    seqid = []
    seq = []
    for fasta in records:
        array = fasta.split('\n')
        name, sequence = array[0].split()[0], re.sub('[^ARNDCQEGHILKMFPSTWYV-]', '', ''.join(array[1:]).upper())
        seqid.append(name)
        seq.append(sequence)
    if len(seqid) == 0:
        f=open(file,"r")
        data1 = f.readlines()
        for each in data1:
            seq.append(each.replace('\n',''))
        for i in range (1,len(seq)+1):
            seqid.append("Seq_"+str(i))
    for i in seq:
        if 'B' in i:
            print('\nError: The input sequences contain non-natural amino acids. Kindly check the sequence.\n')
            sys.exit()
        if 'J' in i:
            print('\nError: The input sequences contain non-natural amino acids. Kindly check the sequence.\n')
            sys.exit()
        if 'O' in i:
            print('\nError: The input sequences contain non-natural amino acids. Kindly check the sequence.\n')
            sys.exit()
        if 'U' in i:
            print('\nError: The input sequences contain non-natural amino acids. Kindly check the sequence.\n')
            sys.exit()
        if 'Z' in i:
            print('\nError: The input sequences contain non-natural amino acids. Kindly check the sequence.\n')
            sys.exit()
        if 'X' in i:
            print('\nError: The input sequences contain non-natural amino acids. Kindly check the sequence.\n')
            sys.exit()
    df4 = pd.DataFrame(seq)
    df5 = pd.DataFrame(seqid)
    return df4,df5
def aab(file,out):
    std = list('ACDEFGHIKLMNPQRSTVWYX')
    df = pd.read_csv(file, header = None)
    uu = []
    for ss in df[0]:
         uu.append(len(ss))
    zz = df.iloc[:,0]
    f = open(out, mode='w')
    sys.stdout = f
    A=('1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0')
    C=('0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0')
    D=('0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0')
    E=('0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0')
    F=('0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0')
    G=('0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0')
    H=('0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0')
    I=('0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0')
    K=('0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0')
    L=('0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0')
    M=('0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0')
    N=('0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0')
    P=('0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0')
    Q=('0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0')
    R=('0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0')
    S=('0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0')
    T=('0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0')
    V=('0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0')
    W=('0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0')
    Y=('0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0')
    X=('0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1')
    for mm in range (1,358):
        print("F_"+str(mm),end=',')
    print("")
    for i in range(0,len(zz)):
        for j in zz[i]:
            if j == "A":
                print(''.join(A), end = ',')
            if j == "C":
                print(''.join(C), end = ',')
            if j == "D":
                print(''.join(D), end = ',')
            if j == "E":
                print(''.join(E), end = ',')
            if j == "F":
                print(''.join(F), end = ',')
            if j == "G":
                print(''.join(G), end = ',')
            if j == "H":
                print(''.join(H), end = ',')
            if j == "I":
                print(''.join(I), end = ',')
            if j == "K":
                print(''.join(K), end = ',')
            if j == "L":
                print(''.join(L), end = ',')
            if j == "M":
                print(''.join(M), end = ',')
            if j == "N":
                print(''.join(N), end = ',')
            if j == "P":
                print(''.join(P), end = ',')
            if j == "Q":
                print(''.join(Q), end = ',')
            if j == "R":
                print(''.join(R), end = ',')
            if j == "S":
                print(''.join(S), end = ',')
            if j == "T":
                print(''.join(T), end = ',')
            if j == "V":
                print(''.join(V), end = ',')
            if j == "W":
                print(''.join(W), end = ',')
            if j == "Y":
                print(''.join(Y), end = ',')
            if j == "X":
                print(''.join(X), end = ',')
        print("")
    f.truncate()
def pcp_bin(file,out):
    std = list('ACDEFGHIKLMNPQRSTVWYX')
    df = pd.read_csv(file, header = None)
    uu = []
    for ss in df[0]:
         uu.append(len(ss))
    zz = df.iloc[:,0]
    f = open(out, mode='w')
    sys.stdout = f
    A=('0,0,1,0,1,1,0,0,0,0,1,1,0,0,0,0,1,0,0,1,0,0,1,1,0')
    C=('0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,1,0,1,0,1,0,0,1,1,0')
    D=('0,1,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,0')
    E=('0,1,0,0,0,0,0,0,1,0,0,0,0,1,0,0,1,0,0,0,1,0,0,0,1')
    F=('0,0,1,0,1,0,0,1,0,0,1,1,0,0,0,0,0,1,0,1,0,0,0,0,1')
    G=('0,0,1,0,1,1,0,0,0,0,1,0,0,1,0,0,0,0,1,1,0,0,1,1,0')
    H=('1,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,1,0,0,0,0,1,0,0,1')
    I=('0,0,1,0,1,1,0,0,0,0,1,1,0,0,0,0,0,1,0,1,0,0,0,0,1')
    K=('1,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1')
    L=('0,0,1,0,1,1,0,0,0,0,1,1,0,0,0,0,1,0,0,1,0,0,0,0,1')
    M=('0,0,1,0,1,0,0,0,0,0,1,1,0,0,0,1,1,0,0,0,0,1,0,0,1')
    N=('0,0,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,1,0,1,0,0,1,0')
    P=('0,0,1,0,1,1,1,0,0,0,1,1,1,0,0,0,0,0,1,0,1,1,0,1,0')
    Q=('0,0,1,1,0,0,0,0,0,0,1,0,0,1,0,0,1,0,0,0,1,0,0,0,1')
    R=('1,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,1')
    S=('0,0,1,1,0,0,0,0,0,0,1,0,0,1,1,0,0,0,1,0,0,1,1,1,0')
    T=('0,0,1,1,0,0,0,0,0,0,1,1,0,1,1,0,0,1,0,0,0,1,0,1,0')
    V=('0,0,1,0,1,1,0,0,0,0,1,1,0,0,0,0,0,1,0,1,0,0,0,1,0')
    W=('0,0,1,0,1,0,0,1,0,0,1,1,0,0,0,0,0,1,0,1,0,0,0,0,1')
    Y=('0,0,1,1,0,0,0,1,0,0,1,0,0,0,0,0,0,1,0,0,0,1,0,0,1')
    X=('0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0')
    for mm in range (358,783):
        print("F_"+str(mm),end=',')
    print("")
    for i in range(0,len(zz)):
        for j in zz[i]:
            if j == "A":
                print(''.join(A), end = ',')
            if j == "C":
                print(''.join(C), end = ',')
            if j == "D":
                print(''.join(D), end = ',')
            if j == "E":
                print(''.join(E), end = ',')
            if j == "F":
                print(''.join(F), end = ',')
            if j == "G":
                print(''.join(G), end = ',')
            if j == "H":
                print(''.join(H), end = ',')
            if j == "I":
                print(''.join(I), end = ',')
            if j == "K":
                print(''.join(K), end = ',')
            if j == "L":
                print(''.join(L), end = ',')
            if j == "M":
                print(''.join(M), end = ',')
            if j == "N":
                print(''.join(N), end = ',')
            if j == "P":
                print(''.join(P), end = ',')
            if j == "Q":
                print(''.join(Q), end = ',')
            if j == "R":
                print(''.join(R), end = ',')
            if j == "S":
                print(''.join(S), end = ',')
            if j == "T":
                print(''.join(T), end = ',')
            if j == "V":
                print(''.join(V), end = ',')
            if j == "W":
                print(''.join(W), end = ',')
            if j == "Y":
                print(''.join(Y), end = ',')
            if j == "X":
                print(''.join(X), end = ',')
        print("")
    f.truncate()
def pssm_1(x):
        if type(x) is str:
            return x
        elif x <= -32768:
             return 0
        elif x:
            return (1/(1+(2.7182)**(-x)))
        else:
            return
def pssm_2(file):
    pd.options.display.float_format = '{:.2e}'.format
    df=file
    df1 = df
    df2 = df1.applymap(pssm_1)
    return df2
def pssm(inputfile,outputfile):
    if os.path.exists('envfile'):
        with open('envfile', 'r') as file:
            data = file.readlines()
        output = []
        for line in data:
            if not "#" in line:
                output.append(line)
        if len(output)==3:
            paths = []
            for i in range (0,len(output)):
                paths.append(output[i].split(':')[1].replace('\n',''))
            blastpgp = paths[0]
            blastdb = paths[1]
            makemat = paths[2]
        if os.path.isfile(blastpgp) and os.access(blastpgp, os.R_OK):
            print('The provided directory for blastpgp is correct and readable.')
        else:
            print("########################################################################################################################")
            print("Error: Either 'blastbgp' file is missing from the provided directory in the 'envfile', or not readable. Kindly check.", file=sys.stderr)
            print("########################################################################################################################")
            sys.exit()
        if os.path.isfile(makemat) and os.access(makemat, os.R_OK):
            print('The provided directory for makemat is correct and readable.')
        else:
            print("########################################################################################################################")
            print("Error: Either 'makemat' file is missing from the provided directory in the 'envfile', or not readable. Kindly check.", file=sys.stderr)
            print("########################################################################################################################")
            sys.exit()
        if (glob.glob(blastdb+".pin")) and (glob.glob(blastdb+".psq")) and (glob.glob(blastdb+".phr")):
            print('The provided directory for blast database is correct and readable.')
        else:
            dbfiles = blastdb.split('/')[-1]
            print("##############################################################################################################################################################################################")
            print("Error: Either the files for BLAST database are missing from the provided directory in the 'envfile', or not readable. Please provide the files with extension of", dbfiles+".pin,", dbfiles+".phr,", dbfiles+".psq in the provided directory.", "Kindly check.", file=sys.stderr)
            print("##############################################################################################################################################################################################")
            sys.exit()

    else:
        print("####################################################################################")
        print("Error: Please provide the '{}', which comprises paths for BLASTPGP and MAKEMAT".format('envfile'), file=sys.stderr)
        print("####################################################################################")
        sys.exit()
    ss = []
    file_list = []
    df1,dff12 = readseq(inputfile)
    aa = []
    for i in df1[0]:
        aa.append(len(i))
    for i in range(0,len(df1)):
        ss.append('>seq_'+str(i+1))
    df1['seq'] = ss
    df1 = df1[['seq',0]]
    for ii in range(0,len(df1)):
        name_file = df1['seq'][ii].replace('>','')+'.fasta'
        file_out = df1['seq'][ii].replace('>','')+'.pssmout'
        file_list.append(df1['seq'][ii].replace('>','')+'.mtx')
        df1.iloc[ii].to_csv(name_file, index=None, header=False, sep='\n')
        filename, file_extension = os.path.splitext(name_file)
        filename_o, file_extension_o = os.path.splitext(file_out)
        S =  blastpgp + ' -d ' +  blastdb + ' -i ' + name_file + ' -j 3 -C ' + file_out
        os.system(S)
        os.rename(file_out, filename_o + '.chk')
        outputfile1 = filename_o + ".chk"
        temp1 =filename_o + ".sn"
        C2 = 'echo {} > {}'.format(name_file,temp1)
        os.system(C2)
        temp2 = filename_o + ".pn"
        C1 = 'echo {} > {}'.format(outputfile1,temp2)
        os.system(C1)
        P = makemat + ' -P ' + filename_o
        os.system(P)
    dir = '.'
    tt = []
    for i in file_list:
        ss = []
        fp = open(i)
        all_line = fp.readlines()
        ss.append(all_line[14:])
        uu = []
        for j in all_line[14:]:
            uu.append(j.replace('\n','').replace('  ',','))
        np.savetxt(i+'_temp',uu, fmt="%s")
        df11 = pd.read_csv(i+'_temp', header=None, sep=",")
        col = [1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,23]
        df13 = df11.iloc[:,col].reset_index(drop=True)
        df12 = pssm_2(df13)
        for i in range(0,len(df12)):
            tt.extend(df12.loc[i])
    bb = []
    cc = 0
    for i in range(0,len(aa)):
        bb.append(tt[cc:cc+21*(aa[i])])
        cc = cc+21*(aa[i])
    df123 = pd.DataFrame(bb)
    df456 = df123.fillna('NA')
    df456.to_csv(outputfile,index=None,header=False,float_format='%.2e')
    allfiles = os.listdir(dir)
    for item in allfiles :
        if item.endswith(".aux") or item.endswith(".sn") or item.endswith(".pn") or item.endswith(".mn") or item.endswith(".mtx_temp") or item.endswith(".mtx") or item.endswith(".chk") or item.endswith(".fasta") or item.endswith(".readseq"):
            os.remove(os.path.join(dir, item))
def pssm2pat(seqfile,infile,nom):
    std = list('ACDEFGHIKLMNPQRSTVWYX')
    df = seqfile
    a = df[0].values
    pp = []
    tt = []
    yy = []
    df1 = infile
    for i in range(0,len(df1)):
        ss = df1.values.tolist()[i][:(len(a[i])*21)]
        uu = [0]*21*(int((nom-1)/2))
        mm = [0]*21*(int((nom-1)/2))
        uu.extend(ss)
        uu.extend(mm)
        pp.append(uu)
    for i in range(0,len(df1)):
        for j in range(0,(len(pp[i])-(int(nom)*21))+1,21):
            tt.append(pp[i][j:j+(21*int(nom))])
    for i in range(1,nom+1):
        for j in std:
            yy.append(j+str(i))
    df2 = pd.DataFrame(tt)
    df2.columns = yy
    return df2
def seq2motif(data_1):
    df = data_1
    cc = []
    for i in range(0,len(df)):
        ss = 'XXXXXXXX'+df[0][i]+'XXXXXXXX'
        for j in range(0,len(ss)):
            if len(ss[j:j+17])==17:
                cc.append(ss[j:j+17])
    df6 = pd.DataFrame(cc)
    return df6
def matrix_gen(data_train):
    matrix_input=np.array(data_train)
    matrix_input.shape
    matrix_3d_form_input= np.zeros([len(data_train),1,len(data_train.columns)], dtype=float)
    matrix_3d_form_input.shape
    for i in range(matrix_input.shape[0]):
        row=matrix_input[i]
        matrix_3d_form_input[i]=np.array(row)
    return matrix_3d_form_input
def pred(data_1,model_1):
    X1 = matrix_gen(data_1)
    model = load_model(model_1)
    train_predictions_weighted = model.predict(X1,batch_size=64)
    df5 = pd.DataFrame(train_predictions_weighted)
    return df5
def pred_hy(data_1,data_2,data_3,model_1):
    X1 = matrix_gen(data_1)
    X2 = matrix_gen(data_2)
    X3 = matrix_gen(data_3)
    model = load_model(model_1)
    train_predictions_weighted = model.predict([X1,X2,X3],batch_size=64)
    df5 = pd.DataFrame(train_predictions_weighted)
    return df5
def result_maker(data_1,data_2,data_3,data_4,num):
    df = data_1
    df1 = data_2
    df2 = data_3
    cc = []
    dd = []
    uu = []
    mm = []
    count = 0
    for i in range(0,len(df1)):
        if df1[0][i] >= num:
            cc.append('+')
        else:
            cc.append('-')
    df1['pre'] = cc
    for i in range(0,len(df)):
        dd.append(len(df[0][i]))
    for i in dd:
        ss = df1['pre'][count:count+i].values
        uu.append(''.join(ss))
    for i in range(0,len(df2)):
        mm.append('>'+df2[0][i])
    df3 = pd.DataFrame()
    df3['ID'] = mm
    df3['Seq'] = df[0]
    df3['Pred'] = uu
    df3.to_csv(data_4,index=None, header=False, sep="\n")
('##############################################################################')
print('# This program PPRInt2.0 is developed for predicting RNA-interacting resdiues  #')
print('# in protein sequence, developed by Prof G. P. S. Raghava group. #')
print('# ############################################################################')


#variables to hold input and output files
inFile= args.input

# Output file

if args.output == None:
     outFile= "outfile.csv"
else:
     outFile = args.output

if args.threshold == None:
        threshold = 0.5
else:
        threshold= float(args.threshold)

# Job Type
if args.job == None:
        Job = int(3)
else:
        Job = int(args.job)

# Python Path
if args.path == None:
        ppath = 'python3'
else:
        ppath = args.path


print('Summary of Parameters:')
print('Input File: ',inFile,'; Threshold: ', threshold,'; Job Type: ',Job)
print('Output File: ',outFile,)

path_test=inFile
path_output=outFile

"""Calling composition py files"""
#======================= Prediction Module start from here =====================
file_gen = str(uuid.uuid4())
s70,s71 = readseq(path_test)
s71.to_csv('seqids.csv',index=None,header=None)
if Job == 1:
    print('==== Initiating the process of Prediction using amino acid binaray profiles as input features: please wait ...')
    dff = seq2motif(s70)
    dff.to_csv(file_gen, index=None, header=False)
    aab(file_gen,'temp.bin')
    df_aab = pd.read_csv('temp.bin')
    df_aab_1 = df_aab.iloc[:,:-1]
    res = pred(df_aab_1,'prog/AAB_model.h5')
    result_maker(s70,res,s71,path_output,threshold)
    os.remove('temp.bin')
    os.remove('seqids.csv')
    os.remove(file_gen)
if Job ==2:
    print('==== Initiating the process of Prediction using physico-chemical properties based  binaray profiles as input features: please wait ...')
    dff = seq2motif(s70)
    dff.to_csv(file_gen, index=None, header=False)
    pcp_bin(file_gen,'temp.pcb')
    df_pcb = pd.read_csv('temp.pcb')
    df_pcb_1 = df_pcb.iloc[:,:-1]
    res = pred(df_pcb_1,'prog/PCB_model.h5')
    result_maker(s70,res,s71,path_output,threshold)
    os.remove('temp.pcb')
    os.remove('seqids.csv')
    os.remove(file_gen)
if Job ==3:
    print('==== Initiating the process of Prediction using position specific scoring matrix profiles as input features: please wait ...')
    pssm(path_test,'temp.pssm')
    df_pssm = pd.read_csv('temp.pssm', header=None)
    df_pssm_1 = pssm2pat(s70,df_pssm,17)
    res = pred(df_pssm_1,'prog/PSSM_model.h5')
    result_maker(s70,res,s71,path_output,threshold)
    os.remove('temp.pssm')
    os.remove('seqids.csv')
    os.remove(file_gen)
if Job==4:
    print('==== Initiating the process of Prediction using hybrid of all features (AAB+PCB+PSSM)s as input features: please wait ...')
    dff = seq2motif(s70)
    dff.to_csv(file_gen, index=None, header=False)
    aab(file_gen,'temp.bin')
    pcp_bin(file_gen,'temp.pcb')
    pssm(path_test,'temp.pssm')
    df_aab = pd.read_csv('temp.bin')
    df_pcb = pd.read_csv('temp.pcb')
    df_pssm = pd.read_csv('temp.pssm', header=None)
    df_aab_1 = df_aab.iloc[:,:-1]
    df_pcb_1 = df_pcb.iloc[:,:-1]
    df_pssm_1 = pssm2pat(s70,df_pssm,17)
    res = pred_hy(df_aab_1,df_pcb_1,df_pssm_1,'prog/Hybrid_model.h5')
    result_maker(s70,res,s71,path_output,threshold)
    os.remove('temp.bin')
    os.remove('temp.pcb')
    os.remove('temp.pssm')
    os.remove('seqids.csv')
    os.remove(file_gen)

print('\n======= Thanks for using PPRInt2.0. Your results are stored in file :',outFile,' =====\n\n')
print('Please cite: PPRInt2.0\n')
print("Et voila. The process is done!!!!!")
