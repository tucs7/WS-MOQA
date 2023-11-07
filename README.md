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

* Before proceeding to sampling ensemble of classifier models has to be prepared (5 datasets --> 5 models). Here we utilize PARROT package from

```
https://github.com/idptools/parrot
```

* For each dataset model training is performed as follows

```
parrot-train data_set_*.tsv seq_class_model.pt --datatype sequence --classes 2 -nl 1 -lr 0.001 -hs 10 -b 16 --epochs 25 --include-figs --probabilistic-classification --set-fractions 0.8 0.1 0.1
```

* Alternatively you can skip the step above and use already trained classifier model ensemble located within /classif_models

* For solubility evaluation we use NetSolP

```
https://github.com/tvinet/NetSolP-1.0
```  

* To sample run

```
python sampler.py
```

For additional details see

```
https://github.com/tucs7/MOQA
```
