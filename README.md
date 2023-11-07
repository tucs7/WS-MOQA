# nanobody_MOQA


# Usage

* Binary autoencoder training

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

* For sampling run

```
python sampler.py
```

For additional details see

```
https://github.com/tucs7/MOQA
```
