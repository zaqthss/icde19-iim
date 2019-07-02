package cn.edu.thu.iim.entity;

/**
 * Created by Aoqian Zhang on 2019/6/30.
 * E-mail address is zaqthss2009@gmail.com
 *  All Rights Reserverd.
 *
 * Tuple class
 *
 * @author Aoqian Zhang
 */
public class Tuple {
  protected int tIndex;
  protected int attrNum;
  protected double[] values;
  protected int[] status; // status = 1, observed; status = 0, missing
  
  public Tuple(int attrNum) {
    this.attrNum = attrNum;
    values = new double[attrNum];
    status = new int[attrNum];
  }
  
  public void buildTuple(int tIndex, double[] vals) {
    settIndex(tIndex);
    
    if (vals.length != attrNum) {
      System.out.println("Inconsistent attrNum !");
    }
    
    for (int i = 0; i < attrNum; ++i) {
      this.values[i] = vals[i];
    }
  }
  
  public void setAttrNum(int attrNum) {
    this.attrNum = attrNum;
  }
  
  public void settIndex(int tIndex) {
    this.tIndex = tIndex;
  }
  
  public int gettIndex() {
    return tIndex;
  }

  public void setStatusbyIndex(int attrIndex, int state) {
    this.status[attrIndex] = state;
  }
  
  public int[] getStatus() {
    return status;
  }
  
  public double[] getAllData() {
    return values;
  }
  
  public double getDataByIndex(int attrIndex) {
    return values[attrIndex];
  }
  
  public void clear() {
    for (int i = 0; i < attrNum; ++i) {
      status[i] = 1;
    }
  }
}
