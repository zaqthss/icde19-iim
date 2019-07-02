package cn.edu.thu.iim;

import Jama.Matrix;
import cn.edu.thu.iim.entity.Cell;
import cn.edu.thu.iim.entity.Database;
import cn.edu.thu.iim.entity.LocalCluster;
import cn.edu.thu.iim.entity.LocalKey;
import cn.edu.thu.iim.entity.Position;
import cn.edu.thu.iim.entity.RegModel;

import cn.edu.thu.iim.util.DataUtil;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Created by Aoqian Zhang on 2019/6/28.
 * E-mail address is zaqthss2009@gmail.com
 * All Rights Reserved.
 *
 * main algorithm
 *
 * @author Aoqian Zhang
 */
public class IIM extends BaseMissing {
  // IIM
  private Map<RegModel, Map<Integer, LocalCluster>> regIndividualMap; // regModel -> rowIndex -> phi
  private int ell = 50;
  private int[] lparams = {1, 200, 20};
  private List<Integer> lIndexList;
  private Map<RegModel, Set<Integer>> regNotLearnedRowIndexSetMap;

  // for impute
  private int K = 10;
  private Map<LocalKey, LocalCluster> clusterMap;

  private double determineTime;
  private double imputeTime;

  public IIM(Database db) {
    super(db);
    tpList = db.getTpList();
    setCellMap();
    regIndividualMap = new HashMap<>();
    regNotLearnedRowIndexSetMap = new HashMap<>();
  }

  public double getDetermineTime() {
    return determineTime;
  }

  public double getImputeTime() {
    return imputeTime;
  }

  public void setParams(int[] lparams, int K) {
    for (int i = 0; i < lparams.length; ++i) {
      this.lparams[i] = lparams[i];
    }
    int LBEGIN = lparams[0];
    int LMAX = lparams[1];
    int INTERVAL = lparams[2];

    lIndexList = new ArrayList<>();
    for (int l = LBEGIN; l < LMAX; l += INTERVAL) {
      lIndexList.add(l);
    }
    Collections.sort(lIndexList);

    this.K = K;
  }

  /**
   * default for isRKNN is true, i.e., use the weighted sum distance
   * @param isInc: if use incremental computing
   * @param equalSigma: use same sigma
   * @param isRKNN: use weighted sum distance
   */
  public Map<Position, Cell> mainIIM(boolean isInc, boolean equalSigma, boolean isRKNN) {
    initVals();

    if (lparams[1] > comRowIndexList.size()) {
      lparams[1] = comRowIndexList.size();
      // reset the size
      setParams(lparams, K);
    }

    regIndividualMap = new HashMap<>();

    double startTime = System.currentTimeMillis();
    rknnLearn(isInc, equalSigma);
    determineTime = System.currentTimeMillis() - startTime;

    startTime = System.currentTimeMillis();
    Map<Position, Cell> repairedCells = impute(equalSigma, isRKNN);
    imputeTime = System.currentTimeMillis() - startTime;
    algTime = imputeTime;

    return repairedCells;
  }

  /**
   * Learned phase
   * @param isInc whether incremental learning
   * @param equalSigma equal sigma or different sigma when combinaton
   */
  private void rknnLearn(boolean isInc, boolean equalSigma) {
    List<RegModel> models = DataUtil
        .findMissingRegModels(misRowIndexList, misRowAttrIndexMap, db.getAttrNum());
    System.out.println("There are " + models.size() + " models");

    int LSIZE = lIndexList.size();

    // for each model
    for (int i = 0; i < models.size(); ++i) {
      RegModel regModel = models.get(i);
      System.out.println("model + " + i + " " + regModel.toString());
      int[] tmpAttrXs = regModel.getAttrXs();
      int attrY = regModel.getAttrY();

      double[] distances = new double[db.getLength()];
      int[] knnIndexes = new int[K];
      double[] knnDistances = new double[K];

      double[] phis;

      // Learning
      // rowIndex -> List<LocalCluster> for each ell
      Map<Integer, List<LocalCluster>> rowLocalMap = new HashMap<>();
      Map<Integer, LocalCluster> rowIndividualMap = new HashMap<>();
      for (int rowIndex : comRowIndexList) {
        List<LocalCluster> clusterList;
        if (isInc) {
          clusterList = learnLocalRegressionIncremental(rowIndex, regModel, lIndexList, equalSigma);
        } else {
          clusterList = learnLocalRegression(rowIndex, regModel, lIndexList, equalSigma);
        }
        rowLocalMap.put(rowIndex, clusterList);
      }

      // voting based on reverse kNN
      double realVal, estimate, residual;
      LocalCluster iCluster;

      int[][] hitNum = new int[db.getLength()][LSIZE];
      double[][] sumRes = new double[db.getLength()][LSIZE];

      for (int rowIndex : comRowIndexList) {
        realVal = dbVals[rowIndex][attrY];

        calcDistance(rowIndex, distances, tmpAttrXs, true);
        findCompleteKnn(distances, knnIndexes, knnDistances);

        // for each neighbor
        for (int ki = 0; ki < K; ++ki) {
          int kRowIndex = knnIndexes[ki];
          double minResidual = Double.MAX_VALUE;
          int targetLi = -1;
          for (int li = 0; li < LSIZE; ++li) {
            iCluster = rowLocalMap.get(kRowIndex).get(li);

            phis = iCluster.getPhis();
            estimate = getRegEstimation(rowIndex, tmpAttrXs, attrY, phis);
            residual = Math.abs(estimate - realVal);
            sumRes[kRowIndex][li] += residual * residual;
            if (residual < minResidual) {
              minResidual = residual;
              targetLi = li;
            }
          } // end of li
          hitNum[kRowIndex][targetLi]++;
        } // end of ki
      } // end of rowIndex

      // what if a tuple is never used to be another's neighbor
      Set<Integer> notLearnedRowIndexSet = new HashSet<>();
      regNotLearnedRowIndexSetMap.put(regModel, notLearnedRowIndexSet);
      for (int rowIndex : comRowIndexList) {
        double minRes = Double.MAX_VALUE;
        int targetLi = -1;
        for (int li = 0; li < LSIZE; ++li) {
          if (sumRes[rowIndex][li] < minRes) {
            minRes = sumRes[rowIndex][li];
            targetLi = li;
          }
        }
        if (targetLi == -1) {
          targetLi = 1;
          notLearnedRowIndexSet.add(rowIndex);
        }

        rowIndividualMap.put(rowIndex, rowLocalMap.get(rowIndex).get(targetLi));
      }
      regIndividualMap.put(regModel, rowIndividualMap);
    }
    System.out.println("rknnLearn over");
  }

  /**
   * Impute phase, for each tuple,
   * @param equalSigma equal or different
   * @param isRKNN
   * @return cell map with repaired value
   */
  private Map<Position, Cell> impute(boolean equalSigma, boolean isRKNN) {
    Map<Position, Cell> repairedCells = new HashMap<>();
    int misRowNum = misRowIndexList.size();

    // whether to learn the local regression offline
    if (clusterMap == null) {
      clusterMap = new HashMap<>();
    }

    Position position;
    Cell cell;
    int rowIndex;
    double modify, estimate;
    double[] phis;

    double[][] subDistances = new double[misRowNum][db.getLength()];
    int[][] sIndexes = new int[misRowNum][K];
    double[][] sDistances = new double[misRowNum][K];

    int[] knnIndexes;

    LocalKey lKey;
    LocalCluster iCluster;

    for (int ri = 0; ri < misRowNum; ++ri) {
      rowIndex = misRowIndexList.get(ri);

      List<Integer> misList = misRowAttrIndexMap.get(rowIndex);
      int[] tmpAttrXs = DataUtil.getAttrXsFromMisList(db.getAttrNum(), misList); // identify LHS

      // may be have more than one missing
      for (int misAttrIndex : misList) {
        int sumNotLearnedMissNum = 0;
        position = new Position(rowIndex, misAttrIndex);
        cell = new Cell(position);
        RegModel regModel = new RegModel(tmpAttrXs, misAttrIndex);

        int L = ell;
        int lIndex = -1;
        if (regIndividualMap.containsKey(regModel)) {
          lIndex = 0;
        }

        // true or false does not matter, since it is incomplete
        calcDistance(rowIndex, subDistances[ri], tmpAttrXs, true);
        findCompleteKnn(subDistances[ri], sIndexes[ri], sDistances[ri]);

        knnIndexes = sIndexes[ri];

        double[] estimates = new double[K];
        double[] sigma2s = new double[K];

        for (int ki = 0; ki < K; ++ki) {
          int kRowIndex = knnIndexes[ki];

          if (regNotLearnedRowIndexSetMap.get(regModel).contains(kRowIndex)) {
            sumNotLearnedMissNum++;
          }

          if (lIndex > -1) {
            // It must have been computed
            iCluster = regIndividualMap.get(regModel).get(kRowIndex);
          } else {
            System.out.println("Only hit if no determine process");
            lKey = new LocalKey(regModel, kRowIndex);
            if (!clusterMap.containsKey(lKey)) {
              List<Integer> lList = new ArrayList<>();
              lList.add(L);
              LocalCluster tmpCluster = learnIndividualRegression(kRowIndex, regModel, lList);
              clusterMap.put(lKey, tmpCluster);
            }
            iCluster = clusterMap.get(lKey);
          }

          phis = iCluster.getPhis();
          estimate = getRegEstimation(rowIndex, tmpAttrXs, misAttrIndex, phis);

          estimates[ki] = estimate;
          if (!equalSigma) {
            sigma2s[ki] = iCluster.getSigma2();
          }
        } // end of k
        if (isRKNN) {
          modify = getWeightedByDistance(estimates);
        } else {
          if (equalSigma) {
            modify = getWeightedWithoutSigma(estimates);
          } else {
            modify = getWeightedBySigma(estimates, sigma2s);
          }
        }

        if (sumNotLearnedMissNum > 0) {
          System.out.println("AttrY = " + misAttrIndex + " has miss " + sumNotLearnedMissNum);
        }

        cell.setValue(modify);
        repairedCells.put(position, cell);
      } // end of misAttrIndex
    } // end of misRowNum

    return repairedCells;
  }

  /**
   * for each l in lIndexList, only one LocalCluster
   *
   * @param rowIndex rowIndex
   * @param regModel regModel(attrXs -> attrY)
   * @param lIndexList selected l to learn
   * @return LocalCluster with minimum residual in lIndexList
   */
  private LocalCluster learnIndividualRegression(int rowIndex, RegModel regModel,
      List<Integer> lIndexList) {
    LocalCluster cluster;

    int maxL = lIndexList.get(lIndexList.size() - 1);

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

    double minResidual = Double.MAX_VALUE;
    int targetL = -1;
    double[] targetPhis = new double[columnSize];

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
      phis = tunePhis(phis, lyMatrix);
      double sigma2 = calcModelResidual(phis, attrY, lxMatrix, lyMatrix, ell);

      if (sigma2 < minResidual && ell > 1) {
        minResidual = sigma2;
        targetL = ell;
        targetPhis = phis;
      }
    }
    cluster = new LocalCluster(rowIndex, targetL);
    cluster.setPhis(targetPhis);
    cluster.setSigma2(minResidual);

    return cluster;
  }
}