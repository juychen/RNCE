# Historical results
Integrative pharmacogenomics to infer large-scale drug taxonomy 


## Apr 01 Latest

Data           |AUROC (ATC)|AUPRC (ATC)|AUROC (DTA)|AUPRC (DTA)| Dir
---------------|:---------:|----------:|----------:|----------:|-----------------
CTRPV          | 0.866     | 0.620     | 0.909     | 0.477     | 20190401
NCI60          | 0.897     | 0.554     | 0.881     | 0.573     | 20190401

```
Previous version are lce after mod jaccard. This version is modjaccard then lce.
We may explore more and discuss
```

## Mar 13 Latest

Data           |AUROC (ATC)|AUPRC (ATC)|AUROC (DTA)|AUPRC (DTA)| Dir
---------------|:---------:|----------:|----------:|----------:|-----------------
CTRPV          | 0.912     | 0.665     | 0.894     | 0.411     | WedMar131101572019

### Notes:

```
Using krnn_lce in the KRNNRR part, set k2 at a bit larger than 4 (8-16)
```


## Mar 04 Latest

Data           |AUROC (ATC)|AUPRC (ATC)|AUROC (DTA)|AUPRC (DTA)| Dir
---------------|:---------:|----------:|----------:|----------:|-----------------
CTRPV          | 0.905     | 0.684     | 0.891     | 0.442     | MonMar041718452019
NCI60          | 0.902     | 0.557     | 0.883     | 0.563     | MonMar042127162019

### Notes:

```
constPerturbationLayerDistSca replace the compeleted verion of senS layer with one has na.
Apply gaussian kernel in the pertlayer
Note that the NCI60 is in a different parameter and feature modeling
```

## Usorted History Notes 

```
# D:\rws\eDNF\Output\ MonMar041707112019
# constPerturbationLayerDistSca
# Apply gaussian kernel in the pertlayer
# ATC (0.916,0.684) D-TARGET (0.885,0.412)



# integrateStrctSensPert_post_highsetkrnn 192 4
# D:\rws\eDNF\Output\ MonMar041504192019
# ATC (0.854,0.611) D-TARGET (0.895,0.421)

# D:\rws\eDNF\Output\ D:\rws\eDNF\Output\ MonMar041444322019
# Use post krnn re inter
# ATC (0.907,0.676) D-TARGET (0.880,0.398)

## Carefull, below this line, possible bug in krnn code

# D:\rws\eDNF\Output\ ThuFeb211718412019
# Use post krnn lowset 192,32
# ATC (0.868,0.621) D-TARGET (0.911,0.428)


# D:\rws\eDNF\Output\ WedFeb201124432019
# integrateStrctSensPert_post_krnn
# k1,k2 = 192, 8 use the completion, ejaccard of perafmat, 
# ctrpv2 (ROC,PRC): 
# ATC (0.893,0.649) D-TARGET (0.901,0.412)

# Record 20190128 D:\rws\eDNF\Output\ MonJan281642382019
# Add the range01 in all Affmat and the post intergration part
# integrateStrctSensPert_post_krnn
# k1,k2 = 192, 128 use the completion, ejaccard of perafmat, 
# ctrpv2 (ROC,PRC): 
# ATC (0.893,0.649) D-TARGET (0.901,0.412)

# Record 20190124 ThuJan242039522019
#integrtStrctSensPert <- integrateStrctSensPert_post_krnn(sensAffMat,strcAffMat, pertAffMat,k1=192,k2=128)
# Use krnn lambda  discard sca, direct using jacard then lce 
# The contextual are snf with the original weight
# k1,k2 = 192, 128 use the completion, ejaccard of perafmat, 
# ctrpv2 (ROC,PRC): 
# ATC (0.896,0.644) D-TARGET (0.894,0.411)
# IF we add the o-1 range to the original layers, (0.899,0.65) (0.894, 0.407) D:\rws\eDNF\Output\ MonJan281600292019


# Record 20190123
#integrtStrctSensPert <- integrateStrctSensPert_post_krnn(sensAffMat,strcAffMat, pertAffMat,k1=192,k2=128)
# Use krnn lambda  discard sca, direct using jacard then lce 
# The contextual are snf with the original weight
# k1,k2 = 192, 128 use the completion, manhatan of perafmat, 
# ctrpv2 (ROC,PRC): 
# ATC (0.896,0.644) D-TARGET (0.894,0.411)

# Record 20190123
#integrtStrctSensPert <- integrateStrctSensPert_post_krnn(sensAffMat,strcAffMat, pertAffMat,k1=192,k2=128)
# Use knn lambda  discard sca, direct using jacard then lce 
# The contextual are snf with the original weight
# k1,k2 = 192, 128 use the completion, manhatan of perafmat, 
# ctrpv2 (ROC,PRC): 
# ATC (0.893,0.638) D-TARGET (0.894,0.410)


# Record 20190122
# Use knn lambda ten discard sca, direct using jacard then lce 
# k1,k2 = 192, 128 use the completion, manhatan of perafmat, norm mean, lambda 0
# ctrpv2 (ROC,PRC): 
# ATC (0.922,0.83) D-TARGET (0.822,0.125)

# Record 20190122
# Use krnn lambda ten discard sca, direct using jacard then lce 
# k1,k2 = 128, 64 use the completion, manhatan of perafmat, norm mean, lambda 0
# ctrpv2 (ROC,PRC): 
# ATC (0.914,0.8) D-TARGET (0.813,0.113)


# Record 20190104
# Use post combing dnf(dnf,sca)
# k1,k2 = 128, 64 use the completion, manhatan of perafmat, norm mean, lambda 0
# ctrpv2 (ROC,PRC): 
# ATC (0.829,0.598) D-TARGET (0.905,0.432)

# Record 20190104
# Use post combing dnf(dnf,sca)
# k1,k2 = 128, 64 use the completion, manhatan of perafmat, norm softmax, lambda 0
# ctrpv2 (ROC,PRC): 
# ATC (0.841,0.588) D-TARGET (0.902,0.435)

# Record 20190104
# Use post combing dnf(dnf,lce)
# k1,k2 = 128, 64 use the completion, manhatan of perafmat, norm softmax, lambda 0
# ctrpv2 (ROC,PRC): 
# ATC (0.841,0.586) D-TARGET (0.895,0.421)

# Record 20190104 
# k1,k2 = 128, 64 use the completion, manhatan of perafmat, lce/cnf softmax, lambda 0
# ctrpv2 (ROC,PRC): 
# ATC (0.892,0.752) D-TARGET (0.871,0.327)
# 
# k1,k2 = 128, 64 use the completion, manhatan of perafmat, lce/cnf norm mean, lambda 0
# ctrpv2 (ROC,PRC): 
# ATC (0.888,0.733) D-TARGET (0.880,0.347)
#
# k1,k2 = 128, 64 use the completion, manhatan of perafmat, lce/cnf norm mean, lambda 0
# ctrpv2 (ROC,PRC): 
# ATC (0.891,0.748) D-TARGET (0.878,0.332)

# Record before 20190104 
# k1,k2 = 50/100 can enhance ATC to AROC = 0.85, APRC=0.68
# Record k1,k2 = 100/50 can also have good performance
```