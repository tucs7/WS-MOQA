# nanobody_MOQA


# Install

* Create a new environment

```
  conda env create -n moqa -f moqa_env.yml
```

* Install fmqa package that provides a trainable binary quadratic mode within the folder from

```
  https://github.com/tsudalab/fmqa
```

# Usage

* Train binary autoencoder

```
python train.py experiment_configs/binary.json
```

* Prepare ensemble of classifier models (5 datasets --> 5 models). We utilize PARROT package from

```
https://github.com/idptools/parrot
```

* For each dataset model training is performed as follows

```
parrot-train data_set_*.tsv seq_class_model.pt --datatype sequence --classes 2 -nl 1 -lr 0.001 -hs 10 -b 16 --epochs 25 --include-figs --probabilistic-classification --set-fractions 0.8 0.1 0.1
```

* Alternatively you may skip the step above and use already trained classifier model ensemble located within /classif_models

* For solubility evaluation we use NetSolP

```
https://github.com/tvinet/NetSolP-1.0
```  

* To sample run

```
python sampler.py
```

* Note, that the default sampler is simulated anealing based in order to use quantum annealer you need a passcode from DWave

```
https://cloud.dwavesys.com/leap/signup/
```

* When you get the passcode in the sample.py replace

```
sampler = dimod.samplers.SimulatedAnnealingSampler()
```

with 

```
bqm = dimod.BinaryQuadraticModel(model)
sampler = EmbeddingComposite(DWaveSampler(endpoint='https://cloud.dwavesys.com/sapi', token='YOUR PASSCODE', solver='Advantage_system4.1'))
```

and in the loop part after model training add

```
bqm = dimod.BinaryQuadraticModel(model)
```

See also

```
https://github.com/tucs7/MOQA
```
