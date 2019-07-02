package cn.edu.thu.iim.entity;

import java.util.Arrays;

/**
 * Regression model
 * attrXs -> attrY, e.g., A1, A2, A3, A4, A5 -> A6
 *
 * @author Aoqian Zhang
 */
public class RegModel {
  private int[] attrXs;
  private int attrY;
  private double[] phi;

  public RegModel(int[] attrXs, int attrY) {
    setAttrXs(attrXs);
    setAttrY(attrY);
  }

  public int[] getAttrXs() {
    return attrXs;
  }

  public void setAttrXs(int[] attrXs) {
    // notice the attrXs can change this.attrXs
    this.attrXs = attrXs;
    Arrays.sort(this.attrXs);
  }

  public int getAttrY() {
    return attrY;
  }

  public void setAttrY(int attrY) {
    this.attrY = attrY;
  }

  public double[] getPhi() {
    return phi;
  }

  public void setPhi(double[] phi) {
    this.phi = phi;
  }

  @Override
  public int hashCode() {
    final int prime = 31;
    int result = 1;
    
    for (int attrX : attrXs) {
      result = prime * result + attrX;
    }
    result = prime * result + attrY;
    return result;
  }

  @Override
  public boolean equals(Object obj) {
    if (this == obj) {
      return true;
    }
    if (obj == null) {
      return false;
    }
    if (getClass() != obj.getClass()) {
      return false;
    }
    RegModel other = (RegModel) obj;
    if (attrY != other.getAttrY()) {
      return false;
    }
    if (attrXs.length != other.getAttrXs().length) {
      return false;
    }
    
    for (int i = 0; i < attrXs.length; ++i) {
      if (attrXs[i] != other.getAttrXs()[i]) {
        return false;
      }
    }
    return true;
  }
  
  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    
    sb.append("[");
    for (int attrX : attrXs) {
      sb.append(attrX).append(",");
    }
    sb.deleteCharAt(sb.length() - 1);
    sb.append("]").append("->").append(attrY);
    
    return sb.toString();
  }
}
