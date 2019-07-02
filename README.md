# icde19-iim
Code release of "Learning Individual Models for Imputation"

Code release of ["Learning Individual Models for Imputation" (ICDE 19)](https://ieeexplore.ieee.org/document/8731351).
The description of code files are listed below:

- `IIM.java`: Algorithm 1,2,3 in the paper. Use IIM algorithm to impute missing values.
- `Database.java`: the class for Relation(database instance) indicating a relation.
- `LocalCluster.java`: the class for Individual model for a certain tuple with a given $\ell$.

Datasets
----------
We use following nine datasets in the paper (details are in Experiment Section):

- ASF, CCS, CCPP and SN in [UCI](https://archive.ics.uci.edu/ml/datasets.php)
- CA, DA, MAM and HEP in [KEEL](https://sci2s.ugr.es/keel/datasets.php)
- Siemens

The schema of the data file is shown in `data/asf.csv`

Attention

- The example dataset is `data/asf.csv` which is clean and `data/asf1_0.1miss.csv` which exists 150 injected missing on Attribute A2.
- Due to the resolution of computing, using incremental algorithm may have a little difference.

Parameters
----------
The input and output of **IIM** algorithm is:

Method

```
setParams(lparams, K);
mainIIM(isInc, equalSigma, isRKNN);
```

Input:

```
int[] lparams = {1, 101, 5}; // LBEGIN, LMAX, INTERVAL(STEP), the range of adaptive learning on the number of learning neighbors(Algorithm 3)
int K = 10;  // the number of learning and imputing neighbors

isInc = True; // use incremental computing
equalSigma = True; // ignore the difference sigma when do combination (impute)
isRKNN = True; // use weighted distance when do combination (impute)
```

Output

```
Map<Position, Cell> repairedCells
```

Library
----------
[jama.jar](http://math.nist.gov/javanumerics/jama/) is used to implement the code. You can get it by maven.

Citation
----------
If you use this code for your research, please consider citing:

```
@inproceedings{DBLP:conf/icde/ZhangSSW19,
  author    = {Aoqian Zhang and
               Shaoxu Song and
               Yu Sun and
               Jianmin Wang},
  title     = {Learning Individual Models for Imputation},
  booktitle = {{ICDE}},
  pages     = {160--171},
  publisher = {{IEEE}},
  year      = {2019}
}
```