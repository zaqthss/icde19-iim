package cn.edu.thu.iim.entity;

/**
 * to identify the position of cell Generally, tIndex = rowIndex+1
 *
 * @author Aoqian Zhang
 */
public class Position {
  private int tIndex;
  private int attrIndex;

  public Position(int tIndex, int attrIndex) {
    this.tIndex = tIndex;
    this.attrIndex = attrIndex;
  }

  public int gettIndex() {
    return tIndex;
  }

  public void settIndex(int tIndex) {
    this.tIndex = tIndex;
  }

  public int getAttrIndex() {
    return attrIndex;
  }

  public void setAttrIndex(int attrIndex) {
    this.attrIndex = attrIndex;
  }

  @Override
  public int hashCode() {
    final int prime = 31;
    int result = 1;
    result = prime * result + attrIndex;
    result = prime * result + tIndex;
    // for (int tempAttrIndex : misList) {
    // result = prime * result + (tempAttrIndex + 1);
    // }
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
    Position other = (Position) obj;
    if (attrIndex != other.attrIndex) {
      return false;
    }
    if (tIndex != other.tIndex) {
      return false;
    }
    return true;
  }
}
