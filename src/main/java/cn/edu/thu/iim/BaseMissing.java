package cn.edu.thu.iim;

import Jama.Matrix;
import cn.edu.thu.iim.entity.Cell;
import cn.edu.thu.iim.entity.Database;
import cn.edu.thu.iim.entity.Position;
import cn.edu.thu.iim.entity.RegModel;
import cn.edu.thu.iim.entity.Tuple;
import cn.edu.thu.iim.entity.KnnPair;
import cn.edu.thu.iim.entity.LocalCluster;
import cn.edu.thu.iim.util.ComparatorKnnPair;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by Aoqian Zhang on 2019/6/30.
 * E-mail address is zaqthss2009@gmail.com
 * All Rights Reserved.
 *
 * Base Class
 *
 * @author Aoqian Zhang
 */
public class BaseMissing {
  protected Database db;
  protected List<Tuple> tpList;
  protected List<Cell> cells; // missing cells
  protected Map<Position, Cell> cellMap;

  protected double algTime = 0;

  protected double[][] dbVals;
  protected int[][] flags;
  protected double EPSILON = 0.01;
  protected double BOUNDMAX = 2E4;
  protected double BOUNDMIN = -2E4;
  public static double ALPHA = 0.001;

  protected List<Integer> rowIndexList;                 // all index list
  protected List<Integer> misRowIndexList;              // tuple as unit to impute
  protected List<Integer> comRowIndexList;              // tuples with no missing
  protected Map<Integer, List<Integer>> misRowAttrIndexMap;
  protected Map<Integer, List<Integer>> misAttrRowIndexMap;
  protected List<Integer> stableList;
  protected List<Integer> misAttrList;                  // each attrIndex here can be

  public BaseMissing(Database db) {
    this.db = db;
    this.cells = db.getCells();
  }

  public double getAlgTime() {
    return algTime;
  }

  protected void setCellMap() {
    cellMap = new HashMap<>();
    int missNum = cells.size();

    for (int i = 0; i < missNum; ++i) {
      Cell cell = cells.get(i);
      cellMap.put(cell.getPosition(), cell);
    }
  }

  /**
   * fill all the indexList
   */
  protected void initVals() {
    int size = db.getLength();
    int attrNum = db.getAttrNum();
    int[] missNums = new int[attrNum];

    dbVals = new double[size][attrNum];
    flags = db.getFlags();
    Tuple tp;

    rowIndexList = new ArrayList<>();
    misRowIndexList = new ArrayList<>();
    comRowIndexList = new ArrayList<>();

    // initialize the lists
    misRowAttrIndexMap = new HashMap<>();
    misAttrRowIndexMap = new HashMap<>();
    for (int i = 0; i < attrNum; ++i) {
      misAttrRowIndexMap.put(i, new ArrayList<>());
    }

    stableList = new ArrayList<>();
    misAttrList = new ArrayList<>();

    // find the missing attributes and collect numbers
    boolean hasMis;
    for (int i = 0; i < size; ++i) {
      tp = tpList.get(i);
      double[] datas = tp.getAllData();
      hasMis = false;

      rowIndexList.add(i);
      for (int j = 0; j < attrNum; ++j) {
        dbVals[i][j] = datas[j];

        if (flags[i][j] == 0) {
          hasMis = true;
          if (misRowAttrIndexMap.containsKey(i)) {
            misRowAttrIndexMap.get(i).add(j);
          } else {
            misRowIndexList.add(i);
            List<Integer> attrList = new ArrayList<>();
            attrList.add(j);
            misRowAttrIndexMap.put(i, attrList);
          }
          missNums[j]++;

          misAttrRowIndexMap.get(j).add(i);
        }
      }

      if (!hasMis) {
        comRowIndexList.add(i);   // complete tuples
      }
    }

    assert misRowIndexList.retainAll(db.getMisRowIndexList());
    db.setComRowIndexList(comRowIndexList);
    db.setDbVals(dbVals);

    for (int attrIndex = 0; attrIndex < attrNum; ++attrIndex) {
      if (missNums[attrIndex] == 0) {
        stableList.add(attrIndex);
      } else {
        misAttrList.add(attrIndex);
      }
    }
  }

  /**********************
   * Regression methods
   */

  protected Matrix getLassoMatrix(int attrXNum) {
    // lasso
    double[][] alphaE = new double[attrXNum][attrXNum];
    for (int i = 0; i < attrXNum; ++i) {
      alphaE[i][i] = ALPHA;
    }
    return new Matrix(alphaE);
  }

  /**
   * learn the parameter can use incremental computation for x and y
   * In case the matrix is singular
   */
  protected Matrix learnParamsOLS(Matrix xMatrix, Matrix yMatrix) {
    int attrXNum = xMatrix.getColumnDimension();

    Matrix xMatrixT = xMatrix.transpose();
    Matrix aMatrix = xMatrixT.times(xMatrix);
    Matrix lassoMatrix = getLassoMatrix(attrXNum);
    aMatrix = aMatrix.plus(lassoMatrix);
    Matrix bMatrix = xMatrixT.times(yMatrix);

    Matrix middleMatrix = aMatrix.inverse();
    Matrix phi = middleMatrix.times(bMatrix);

    return phi;
  }

  protected Matrix learnParamsOLSInc(Matrix xMatrix, Matrix yMatrix,
      Matrix[] middleMatrices, boolean isFirst) {
    int attrXNum = xMatrix.getColumnDimension();

    Matrix xMatrixT = xMatrix.transpose();

    Matrix aMatrix = xMatrixT.times(xMatrix);
    Matrix bMatrix = xMatrixT.times(yMatrix);
    if (isFirst) {
      Matrix lassoMatrix = getLassoMatrix(attrXNum);
      aMatrix = aMatrix.plus(lassoMatrix);
      middleMatrices[0] = aMatrix;
      middleMatrices[1] = bMatrix;
    } else {
      middleMatrices[0].plusEquals(aMatrix);
      middleMatrices[1].plusEquals(bMatrix);
    }

    Matrix middleMatrix = middleMatrices[0].inverse();
    Matrix phi = middleMatrix.times(middleMatrices[1]);

    return phi;
  }

  /**
   * Can handle the parameter with constants with the real size, for
   * LocalRegression
   *
   * @param xMatrix there may contains all zeros lines, which denote ignoring them
   * @param size the real size which should be considered
   */
  protected double calcModelResidual(double[] phis, int attrY, Matrix xMatrix,
      Matrix yMatrix, int size) {
    double sigma = 0;

    double[][] x = xMatrix.getArray();
    double[][] y = yMatrix.getArray();

    double estimate, residual;
    for (int i = 0; i < size; ++i) {
      estimate = 0;
      for (int j = 0; j < phis.length; ++j) {
        estimate += phis[j] * x[i][j];
      }
      residual = estimate - y[i][0];
      sigma += residual * residual;
    }

    sigma = sigma / size;
    return sigma;
  }

  /**
   * for each l in lIndexList
   *
   * @param rowIndex rowIndex
   * @param regModel regModel(attrXs -> attrY)
   * @param lIndexList selected l to learn
   * @return List<LocalCluster> for each l in lIndexList
   */
  protected List<LocalCluster> learnLocalRegression(int rowIndex, RegModel regModel,
      List<Integer> lIndexList, boolean equalSigma) {
    List<LocalCluster> clusterList = new ArrayList<>();
    LocalCluster cluster;

    int maxL = lIndexList.get(lIndexList.size() - 1);
    maxL = Math.min(maxL, comRowIndexList.size());

    double[] wholeDistances = new double[db.getLength()];
    int[] knnIndexes = new int[maxL];
    double[] knnDistances = new double[maxL];

    int[] tmpAttrXs = regModel.getAttrXs();
    int attrY = regModel.getAttrY();
    int attrXNum = tmpAttrXs.length;
    int[] usedWholeAttrs = new int[attrXNum + 1];
    for (int i = 0; i < attrXNum; ++i) {
      usedWholeAttrs[i] = tmpAttrXs[i];
    }
    usedWholeAttrs[attrXNum] = attrY;
    int attrX;
    int columnSize = attrXNum + 1;

    calcDistance(rowIndex, wholeDistances, tmpAttrXs, false);
    findCompleteKnn(wholeDistances, knnIndexes, knnDistances);

    for (int ell : lIndexList) {
      double[] phis = new double[columnSize];
      double[][] x = new double[ell][columnSize];
      double[][] y = new double[ell][1];

      for (int ki = 0; ki < ell; ++ki) {
        int rIndex = knnIndexes[ki];

        for (int j = 0; j < attrXNum; ++j) {
          attrX = tmpAttrXs[j];
          x[ki][j + 1] = dbVals[rIndex][attrX];
        }
        x[ki][0] = 1;
        y[ki][0] = dbVals[rIndex][attrY];
      }

      boolean isSingular = false;

      Matrix lxMatrix = new Matrix(x);
      Matrix lyMatrix = new Matrix(y);
      Matrix phi = null;
      // will not be singular using lasso
      try {
        phi = learnParamsOLS(lxMatrix, lyMatrix);
      } catch (Exception e) {
        isSingular = true;
      }
      if (ell == 1) {
        isSingular = true;
      }
      if (isSingular) {
        phis = setSingularPhis(lyMatrix, columnSize);
      } else {
        for (int i = 0; i < columnSize; ++i) {
          phis[i] = phi.get(i, 0);
        }
      }
//      phis = tunePhis(phis, rowIndex, attrY);
      phis = tunePhis(phis, lyMatrix);

      cluster = new LocalCluster(rowIndex, ell);
//      cluster.setKnnIndexes(knnIndexes);
//      cluster.setKnnDistances(knnDistances);
      cluster.setPhis(phis);
      if (!equalSigma) {
        double sigma2 = calcModelResidual(phis, attrY, lxMatrix, lyMatrix, ell);
        cluster.setSigma2(sigma2);
      }
      clusterList.add(cluster);
    }
    return clusterList;
  }

  /**
   * for each l in lIndexList
   *
   * @param rowIndex rowIndex
   * @param regModel regModel(attrXs -> attrY)
   * @param lIndexList selected l to learn
   * @return List<LocalCluster> for each l in lIndexList
   */
  protected List<LocalCluster> learnLocalRegressionIncremental(int rowIndex, RegModel regModel,
      List<Integer> lIndexList, boolean equalSigma) {
    List<LocalCluster> clusterList = new ArrayList<>();
    LocalCluster cluster;

    int maxL = lIndexList.get(lIndexList.size() - 1);
    maxL = Math.min(maxL, comRowIndexList.size());

    double[] wholeDistances = new double[db.getLength()];
    int[] knnIndexes = new int[maxL];
    double[] knnDistances = new double[maxL];

    int[] tmpAttrXs = regModel.getAttrXs();
    int attrY = regModel.getAttrY();
    int attrXNum = tmpAttrXs.length;
    int[] usedWholeAttrs = new int[attrXNum + 1];
    for (int i = 0; i < attrXNum; ++i) {
      usedWholeAttrs[i] = tmpAttrXs[i];
    }
    usedWholeAttrs[attrXNum] = attrY;
    int attrX;
    int columnSize = attrXNum + 1;

    calcDistance(rowIndex, wholeDistances, tmpAttrXs, false);
    findCompleteKnn(wholeDistances, knnIndexes, knnDistances);

    double[][] x = new double[maxL][columnSize];
    double[][] y = new double[maxL][1];

    for (int ki = 0; ki < maxL; ++ki) {
      int rIndex = knnIndexes[ki];

      for (int j = 0; j < attrXNum; ++j) {
        attrX = tmpAttrXs[j];
        x[ki][j + 1] = dbVals[rIndex][attrX];
      }
      x[ki][0] = 1;
      y[ki][0] = dbVals[rIndex][attrY];
    }

    Matrix lxMatrix = new Matrix(x);
    Matrix lyMatrix = new Matrix(y);
    Matrix[] middleMatrices = new Matrix[2];
    boolean isFirst = true;
    for (int li = 0; li < lIndexList.size(); ++li) {
      int ell = lIndexList.get(li);
      int beginRowIndex = isFirst ? 0 : lIndexList.get(li - 1);
      int endRowIndex = ell - 1;
      double[] phis = new double[columnSize];

      boolean isSingular = false;
      Matrix phi = null;
      // will not be singular using lasso
      try {
        phi = learnParamsOLSInc(
            lxMatrix.getMatrix(beginRowIndex, endRowIndex, 0, columnSize - 1),
            lyMatrix.getMatrix(beginRowIndex, endRowIndex, 0, 0),
            middleMatrices, isFirst);
      } catch (Exception e) {
        isSingular = true;
      }
      if (ell == 1) {
        isSingular = true;
      }
      if (isSingular) {
        phis = setSingularPhis(lyMatrix, columnSize);
//        System.out.println("is Singular");
      } else {
        for (int i = 0; i < columnSize; ++i) {
          phis[i] = phi.get(i, 0);
        }
      }
//      phis = tunePhis(phis, rowIndex, attrY);
      phis = tunePhis(phis, lyMatrix);

      cluster = new LocalCluster(rowIndex, ell);
//      cluster.setKnnIndexes(knnIndexes);
//      cluster.setKnnDistances(knnDistances);
      cluster.setPhis(phis);
      if (!equalSigma) {
        double sigma = Math.sqrt(calcModelResidual(phis, attrY, lxMatrix, lyMatrix, ell));
        cluster.setSigma2(sigma * sigma);
      }
      clusterList.add(cluster);

      if (isFirst) {
        isFirst = false;
      }
    }
    return clusterList;
  }

  /*********************
   * kNN methods
   */

  /**
   * compute the distances against dbVals[rowIndex]
   * Can only compute the complete ones
   */
  protected void calcDistance(int rowIndex, double[] distances, int[] usedAttrs,
      boolean itselfSetMax) {
    calcDistanceInSet(rowIndex, distances, usedAttrs, itselfSetMax, comRowIndexList);
  }

  protected void calcDistanceInSet(int rowIndex, double[] distances, int[] usedAttrs,
      boolean itselfSetMax, List<Integer> indexList) {
    calcDistanceInSetInDB(rowIndex, distances, usedAttrs, itselfSetMax, indexList, dbVals);
  }

  protected void calcDistanceInSetInDB(int rowIndex, double[] distances, int[] usedAttrs,
      boolean itselfSetMax, List<Integer> indexList, double[][] dbVals) {
    double[] vals = dbVals[rowIndex];

    double dis;
    double sumUp, sumDown = usedAttrs.length;
    for (int cIndex : indexList) {
      sumUp = 0;

      for (int attrIndex : usedAttrs) {
        dis = vals[attrIndex] - dbVals[cIndex][attrIndex];
        sumUp += dis * dis;
      }
      if (sumDown == 0) {
        // no both observed value
        dis = Double.MAX_VALUE;
      } else if (sumUp == 0) {
        dis = EPSILON;
      } else {
        dis = Math.sqrt(sumUp / sumDown);
      }
      distances[cIndex] = dis;
    }
    // itself
    if (itselfSetMax) {
      distances[rowIndex] = Double.MAX_VALUE;
    }
  }

  /**
   * all the attributes in the kNN results must be complete can use
   * comRowIndexList
   *
   * @param distances have real distances
   * @param knnIndexes the result index
   * @param knnDistances the result distance
   */
  protected void findCompleteKnn(double[] distances, int[] knnIndexes,
      double[] knnDistances) {
    findCompleteKnnInSet(distances, knnIndexes, knnDistances, comRowIndexList);
  }

  /**
   *
   * @param distances real distances computed before
   * @param knnIndexes
   * @param knnDistances
   * @param indexList the complete kNN is found in indexList
   */
  protected void findCompleteKnnInSet(double[] distances, int[] knnIndexes,
      double[] knnDistances, List<Integer> indexList) {
    if (knnDistances.length == 0) {
      return;
    }

    // first K
    int length = knnIndexes.length;
    if (length > indexList.size()) {
      for (int i = 0; i < indexList.size(); i++) {
        int rowIndex = indexList.get(i);
        knnIndexes[i] = rowIndex;
        knnDistances[i] = distances[rowIndex];
      }
    } else {
      for (int i = 0; i < length; ++i) {
        int rowIndex = indexList.get(i);
        knnIndexes[i] = rowIndex;
        knnDistances[i] = distances[rowIndex];
      }
      int maxIndex = getMaxIndexfromK(knnDistances);
      double maxVal = knnDistances[maxIndex];

      // the rest
      double dis;
      for (int i = length; i < indexList.size(); ++i) {
        int rowIndex = indexList.get(i);
        dis = distances[rowIndex];
        if (dis < maxVal) {
          // knnIndexes[maxIndex] = i;
          knnIndexes[maxIndex] = rowIndex;
          knnDistances[maxIndex] = dis;

          maxIndex = getMaxIndexfromK(knnDistances);
          maxVal = knnDistances[maxIndex];
        }
      }
    }

    // sort
    List<KnnPair> kpList = new ArrayList<>();
    KnnPair kp;
    for (int i = 0; i < length; ++i) {
      kp = new KnnPair(knnDistances[i], knnIndexes[i]);
      kpList.add(kp);
    }
    Collections.sort(kpList, new ComparatorKnnPair());

    for (int i = 0; i < length; ++i) {
      kp = kpList.get(i);
      knnIndexes[i] = kp.getIndex();
      knnDistances[i] = kp.getDistance();
    }
  }

  protected int getMaxIndexfromK(double[] vals) {
    int index = -1;

    double max = -1;
    for (int i = 0; i < vals.length; ++i) {
      if (vals[i] > max) {
        max = vals[i];
        index = i;
      }
    }

    return index;
  }

  /*********************
   * Proposed Local Profile
   */
  protected double[] tunePhis(double[] phis, Matrix lyMatrix) {
    double[] tPhis = new double[phis.length];
    double intercept = phis[0];

    if (intercept > BOUNDMAX || intercept < BOUNDMIN) {
      System.out.println("Tune phis hit");
      int rowNum = lyMatrix.getRowDimension();
      double sum = 0;
      for (int i = 0; i < rowNum; ++i) {
        sum += lyMatrix.get(i, 0);
      }
      tPhis[0] = sum / rowNum;
    } else {
      tPhis = phis;
    }

    return tPhis;
  }

  protected double[] setSingularPhis(Matrix lyMatrix, int columnSize) {
    double[] phis = new double[columnSize];
    int rowNum = lyMatrix.getRowDimension();
    double sum = 0;
    for (int i = 0; i < rowNum; ++i) {
      sum += lyMatrix.get(i, 0);
    }
    phis[0] = sum / rowNum;
    return phis;
  }

  /**
   * The formula sigma2=EPSILON*EPSILON may not be proper
   */
  protected double getWeightedBySigma(double[] estimates, double[] sigma2s) {
    double modify = 0;

    double estimate;
    double sigma2, weight, sumSigma2 = 0;

    int kLen = estimates.length;

    for (int ki = 0; ki < kLen; ++ki) {
      estimate = estimates[ki];
      sigma2 = sigma2s[ki];

      if (sigma2 == 0) {
        sigma2 = EPSILON * EPSILON;
      }
      weight = 1 / sigma2;
      sumSigma2 += weight;
      modify += estimate * weight;
    }
    sumSigma2 = 1 / sumSigma2;
    modify = modify * sumSigma2;

    return modify;
  }

  protected double getWeightedWithoutSigma(double[] estimates) {
    double modify = 0;

    int kLen = estimates.length;
    for (double estimate : estimates) {
      modify += estimate;
    }
    modify = modify * (1 / (double) kLen);

    return modify;
  }

  protected double getWeightedByDistance(double[] estimates) {
    double modify = 0;

    int kLen = estimates.length;
    double[][] mutualDistances = new double[kLen][kLen + 1];
    double valI, valJ, sumWeight = 0;
    for (int i = 0; i < kLen; ++i) {
      valI = estimates[i];
      for (int j = 0; j < kLen; ++j) {
        if (i == j) {
          continue;
        } else if (j > i) {
          valJ = estimates[j];
          mutualDistances[i][j] = Math.abs(valI - valJ);
        } else {
          mutualDistances[i][j] = mutualDistances[j][i];
        }
        mutualDistances[i][kLen] += mutualDistances[i][j];
      }

      double weight = mutualDistances[i][kLen] == 0 ? EPSILON : mutualDistances[i][kLen];
      modify += valI / weight;
      sumWeight += 1 / weight;
    }
    modify = modify / sumWeight;

    return modify;
  }

  /**
   * @param attrY tmpAttrXs->attrY
   */
  protected double getRegEstimation(int rowIndex, int[] tmpAttrXs, int attrY,
      double[] phis) {
    double estimate = 0;

    int attrXNum = tmpAttrXs.length;
    double intercept = phis[0];

    int attrX;
    double val;
    for (int j = 0; j < attrXNum; ++j) {
      attrX = tmpAttrXs[j];
      val = dbVals[rowIndex][attrX];
      estimate += val * phis[j + 1];
    }
    estimate += intercept;

    return estimate;
  }
}
