import numpy as np
import fmqa
import dimod
import os
import subprocess
from keras.preprocessing import sequence
from keras.models import model_from_json
from keras.utils.np_utils import to_categorical
from typing import Tuple, Dict
import collections
import torch
from botorch.utils.multi_objective.hypervolume import Hypervolume
from itertools import groupby
from functools import reduce
import random
from parrot import py_predictor as ppp
import pandas as pd


#Loading classifier models
my_predictor_1 = ppp.Predictor('.../nanobody_MOQA/classif_models/network_final_classif1.pt', dtype='sequence')
my_predictor_2 = ppp.Predictor('.../nanobody_MOQA/classif_models/network_final_classif2.pt', dtype='sequence')
my_predictor_3 = ppp.Predictor('.../nanobody_MOQA/classif_models/network_final_classif3.pt', dtype='sequence')
my_predictor_4 = ppp.Predictor('.../nanobody_MOQA/classif_models/network_final_classif4.pt', dtype='sequence')
my_predictor_5 = ppp.Predictor('.../nanobody_MOQA/classif_models/network_final_classif5.pt', dtype='sequence')

#Pareto front identification
def is_pareto_efficient(costs):
    """
    Find the pareto-efficient points
    :param costs: An (n_points, n_costs) array
    :return: A (n_points, ) boolean array, indicating whether each point is Pareto efficient
    """
    is_efficient = np.ones(costs.shape[0], dtype = bool)
    for i, c in enumerate(costs):
        if is_efficient[i]:
            is_efficient[is_efficient] = np.any(costs[is_efficient]>c, axis=1)
            is_efficient[i] = True
    return is_efficient

#Elimination of sequences with non-AA characters
def seq_fix():
    lines_d = open('.../nanobody_MOQA/model_output/binary/decoded.txt','r')
    Lines_d = lines_d.readlines()
    seq_d=[]
    for i in Lines_d:
        seq_d.append(i.replace(" ", "").rstrip())

    matching_d1=[s for s in seq_d if "<pad>" in s]
    matching_d2=[s for s in seq_d if len(s)==0]
    matching_d3=[s for s in seq_d if "<unk>" in s]

    for ii in range(len(matching_d1)):
        for i in range(len(seq_d)):
            if seq_d[i]==matching_d1[ii]:
                seq_d[i]="SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS"

    for ii in range(len(matching_d2)):
        for i in range(len(seq_d)):
            if seq_d[i]==matching_d2[ii]:
                seq_d[i]="SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS"

    for ii in range(len(matching_d3)):
        for i in range(len(seq_d)):
            if seq_d[i]==matching_d3[ii]:
                seq_d[i]="SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS"

    return seq_d

#Objective function evaluation
def objectives(seq_d):

    #Whole VHH sequence construction
    seq_d_VHH=[]
    for i in seq_d:
        seq_d_VHH.append('QVQLVESGGGSVQAGGSLRLSCTAS'+i[0:13]+'WFRQAPGQEREAVA'+i[13:-16]+'ADSVKGRFTISRDNAKNTVTLQMNNLKPEDTAIYYCAA'+i[-16:]+'WGQGTQVTVSS')

    #VHH sequence containing repetitive amino acics identification
    rep_AA_cutoff=5

    counts=[]
    for i in range(len(seq_d_VHH)):
        grouped_L=[]
        grouped_L = [(k, sum(1 for i in g)) for k,g in groupby(seq_d_VHH[i])]
        counts.append(grouped_L)

    score_AA=[]
    for i in counts:
        rep_count=[]
        for j in range(len(i)):
            rep_count.append(i[j][1])
    
        rep_count_binary=[]
        for i in rep_count:
            if i >= rep_AA_cutoff:
                score=0
            else:
                score=1
            rep_count_binary.append(score)
    
        score_AA.append(reduce(lambda x, y: x*y, rep_count_binary))

    #Sequence with CDR1+CDR2+CDR3 length other than 39AA identification
    score_length=[]
    for i in range(len(seq_d)):
     
        if len(seq_d[i])!= 39:
            score=0
        else:
            score=1
    
        score_length.append(score)

    #Sequences with above constraints violated will be penalized by assigning them lowest possible values
    classifiers_ave_ref=0
    classifiers_std_ref=1
    solubility_ref=0
    
    #Sequence evaluation via classifier model ensemble
    classifier_1=[]
    classifier_2=[]
    classifier_3=[]
    classifier_4=[]
    classifier_5=[]
    for i in range(len(seq_d)):
        classifier_1.append(my_predictor_1.predict(seq_d[i])[0])
        classifier_2.append(my_predictor_2.predict(seq_d[i])[0])
        classifier_3.append(my_predictor_3.predict(seq_d[i])[0])
        classifier_4.append(my_predictor_4.predict(seq_d[i])[0])
        classifier_5.append(my_predictor_5.predict(seq_d[i])[0])

    classifiers_ave=[]
    classifiers_std=[]
    for i in range(len(seq_d)):
        if score_AA[i] == 0 or score_length[i] == 0:
            classifiers_ave.append(classifiers_ave_ref)
            classifiers_std.append(classifiers_std_ref)
        else:
            classifiers_ave.append(np.mean([classifier_1[i],classifier_2[i],classifier_3[i],classifier_4[i],classifier_5[i]]))
            classifiers_std.append(np.std([classifier_1[i],classifier_2[i],classifier_3[i],classifier_4[i],classifier_5[i]]))
        
    #Solubility via NetSolP

    oo = open('.../nanobody_MOQA/model_output/binary/seq_VHH_sol.fasta','a')
    for i in range(len(seq_d_VHH)):
        oo.write('>'+str(i) + '\n')
        oo.write(str(seq_d_VHH[i]) + '\n')
    oo.close()

    os.system('python .../netsolp-1.0.ALL/predict.py --MODELS_PATH .../netsolp-1.0.ALL/models --FASTA_PATH .../nanobody_MOQA/model_output/binary/seq_VHH_sol.fasta --OUTPUT_PATH .../nanobody_MOQA/model_output/binary/seq_VHH_sol_preds.csv --MODEL_TYPE Distilled --PREDICTION_TYPE S')

    solubility = pd.read_csv('.../nanobody_MOQA/model_output/binary/seq_VHH_sol_preds.csv', usecols = ['predicted_solubility'], low_memory = True)

    solubility_temp=[]
    for i in solubility.to_numpy():
        solubility_temp.append(i[0])
    
    solubility=[]
    for i in range(len(seq_d)):
        if score_AA[i] == 0 or score_length[i] == 0:
            solubility_c=solubility_ref
        else:
            solubility_c=solubility_temp[i]
        solubility.append(solubility_c)

    os.system('rm .../nanobody_MOQA/model_output/binary/seq_VHH_sol.fasta')
    os.system('rm .../nanobody_MOQA/model_output/binary/seq_VHH_sol_preds.csv')

    return np.array(classifiers_ave), -1.*np.array(classifiers_std), np.array(solubility)

#Hypervolume calculation
def ParetoV(classifiers_solubility):
    
    front_classifiers_solubility_0=is_pareto_efficient(classifiers_solubility)
    
    front_classifiers_solubility=[]
    front_index_classifiers_solubility_original=[]
    for i in range(len(front_classifiers_solubility_0)):
        if front_classifiers_solubility_0[i]== True:
            front_classifiers_solubility.append(classifiers_solubility[i])
            front_index_classifiers_solubility_original.append(i)

    temp0 = []
    for i in front_classifiers_solubility:
        temp0.append(i)
    
    x0 = torch.from_numpy(np.array(temp0))
    
    classifiers_ave_ref=0
    classifiers_std_ref=1
    solubility_ref=0
    
    xref = torch.from_numpy(np.array([classifiers_ave_ref,-1.*classifiers_std_ref,solubility_ref]))
    hv = Hypervolume(xref)
    pv0 = hv.compute(x0)
    
    oo = open('.../nanobody_MOQA/model_output/binary/hv_score.txt','a')
    oo.write(str(pv0) + '\n')
    oo.close()

#Non-dominated sorting procedure with pre-defined number of layers
def NSPareto(classifiers_solubility,classifiers_solubility_original,population_size,NLayers):

    SLayers=np.ones(population_size)*(-1)*10

    for ii in range(NLayers):
    
        front_classifiers_solubility_0=is_pareto_efficient(classifiers_solubility)
        
        front_classifiers_solubility=[]
        front_index_classifiers_solubility=[]
        for i in range(len(front_classifiers_solubility_0)):
            if front_classifiers_solubility_0[i]== True:
                front_classifiers_solubility.append(classifiers_solubility[i])
                front_index_classifiers_solubility.append(i)
    
        index_out=[]
        for i in front_classifiers_solubility:
            index = np.where(i == classifiers_solubility_original)
            index_out.append(index[0][0])
    
        for i in index_out:
            SLayers[i]=1./(ii+1.)
    
        classifiers_solubility = np.delete(classifiers_solubility, front_index_classifiers_solubility, axis=0)
    
    return SLayers


#Decode initial binary vector set
os.system('python3 decode_file.py experiment_configs/binary.json')

c0 = [x.split(' ')[:] for x in open('.../nanobody_MOQA/model_output/binary/vectors.txt').readlines()]

vectors_all=[]
for k in range(len(c0)):
    vectors=[]
    for i in c0[k]:
        vectors.append(int(float(i)))
    vectors_all.append(vectors)

vectors_all=np.array(vectors_all)

#Get rid of wrong sequences
seq_d_all=seq_fix()

#Calculate objective functions
scores_all=objectives(seq_d_all)
scores_all_classifiers_ave=scores_all[0]
scores_all_classifiers_std=scores_all[1]
scores_all_solubility=scores_all[2]

#Output best metrics in the initial set
oo = open('.../nanobody_MOQA/model_output/binary/all_points_best.txt','a')
oo.write(str(np.max(scores_all_classifiers_ave)) + ' ' + str(np.min(-1.*scores_all_classifiers_std)) + ' ' + str(np.max(scores_all_solubility)) + '\n')
oo.close()

classifiers_solubility=[]
for i in range(len(scores_all_classifiers_ave)):
    classifiers_solubility.append([scores_all_classifiers_ave[i],scores_all_classifiers_std[i],scores_all_solubility[i]])

classifiers_solubility_original=np.array(classifiers_solubility)
classifiers_solubility=classifiers_solubility_original

NLayers=20 #number of non-dominated layers to consider

population_size=len(scores_all_classifiers_ave)

#Calculate initial Pareto hypervolume
ParetoV(classifiers_solubility)

#Non-dominated sorting for the initial population
SLayers=NSPareto(classifiers_solubility,classifiers_solubility_original,population_size,NLayers)
scores=(-1.)*SLayers

#FM training
model = fmqa.FMBQM.from_data(vectors_all, scores)
#Simulated annealing for sampling
sampler = dimod.samplers.SimulatedAnnealingSampler()

#Sampling/Evaluation/Training
for _ in range(200):
    
    #Sample
    res = sampler.sample(model, num_reads=10)
    
    vectors_all = np.r_[vectors_all, res.record['sample']]
    vectors_sample = np.r_[res.record['sample']]
        
    vectors_sample_out=[]
    for i in vectors_sample:
        vectors_sample_out.append(''.join(str(i)[1:-1].splitlines()))
    
    os.system('rm .../nanobody_MOQA/model_output/binary/vectors.txt')
    oo = open('.../nanobody_MOQA/model_output/binary/vectors.txt','a')
    for i in vectors_sample_out:
        oo.write(str(i) + '\n')
    oo.close()

    os.system('rm .../nanobody_MOQA/model_output/binary/decoded.txt')
    os.system('python .../nanobody_MOQA/decode_file.py experiment_configs/binary.json')
    seq_d=seq_fix()
    
    seq_d_all = np.r_[seq_d_all, seq_d]
    
    #Calculate objective functions for sampled sequences
    scores_sample=objectives(seq_d)
    
    scores_sample_classifiers_ave=scores_sample[0]
    scores_all_classifiers_ave = np.r_[scores_all_classifiers_ave, scores_sample_classifiers_ave]
    
    scores_sample_classifiers_std=scores_sample[1]
    scores_all_classifiers_std = np.r_[scores_all_classifiers_std, scores_sample_classifiers_std]
    
    scores_sample_solubility=scores_sample[2]
    scores_all_solubility = np.r_[scores_all_solubility, scores_sample_solubility]

    classifiers_solubility=[]
    for i in range(len(scores_all_classifiers_ave)):
        classifiers_solubility.append([scores_all_classifiers_ave[i],scores_all_classifiers_std[i],scores_all_solubility[i]])

    classifiers_solubility_original=np.array(classifiers_solubility)
    classifiers_solubility=classifiers_solubility_original

    #Calculate Pareto hypervolume for the updated population
    ParetoV(classifiers_solubility)
    
    #Non-dominated sorting procedure with gradual layer number reduction
    STEP_GEOM=0.982
    ii= _
    if ii==0:
        NLayerss=20.
    else:
        NLayerss=STEP_GEOM*NLayerss
    
    NLayers=round(NLayerss)
    population_size=len(scores_all_classifiers_ave)
    
    SLayers=NSPareto(classifiers_solubility,classifiers_solubility_original,population_size,NLayers)
    
    scores_all=(-1.)*SLayers
    
    #Output best metrics for the updated population
    oo = open('.../nanobody_MOQA/model_output/binary/all_points_best.txt','a')
    oo.write(str(np.max(scores_all_classifiers_ave)) + ' ' + str(np.min(-1.*scores_all_classifiers_std)) + ' ' + str(np.max(scores_all_solubility)) + '\n')
    oo.close()

    #FM re-training
    model.train(vectors_all, scores_all)
