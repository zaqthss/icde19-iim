package cn.edu.thu.iim.entity;

/**
 * LocalCluster is the individual model learned for tuple:rowIndex
 *
 * @author Aoqian Zhang
 */
public class LocalCluster {
  private int rowIndex;
  private int ell;
  private double[] phis;
  private double sigma2;
  private int[] knnIndexes;
  private double[] knnDistances;

  public LocalCluster(int rowIndex, int ell) {
    this.rowIndex = rowIndex;
    this.ell = ell;
  }

  public int getRowIndex() {
    return rowIndex;
  }

  public void setRowIndex(int rowIndex) {
    this.rowIndex = rowIndex;
  }

  public int getEll() {
    return ell;
  }

  public void setEll(int ell) {
    this.ell = ell;
  }

  public double[] getPhis() {
    return phis;
  }

  public void setPhis(double[] phis) {
    this.phis = phis;
  }

  public double getSigma2() {
    return sigma2;
  }

  public void setSigma2(double sigma2) {
    this.sigma2 = sigma2;
  }

  public int[] getKnnIndexes() {
    return knnIndexes;
  }

  public void setKnnIndexes(int[] knnIndexes) {
    this.knnIndexes = knnIndexes;
  }

  public double[] getKnnDistances() {
    return knnDistances;
  }

  public void setKnnDistances(double[] knnDistances) {
    this.knnDistances = knnDistances;
  }
}
