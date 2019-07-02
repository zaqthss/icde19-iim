package cn.edu.thu.iim.entity;

/**
 * Identify a local cluster
 * 
 * @author Aoqian Zhang
 *
 */
public class LocalKey {
  private RegModel regModel;
  private int rowIndex;
  
  public LocalKey(RegModel regModel, int rowIndex) {
    setRegModel(regModel);
    setRowIndex(rowIndex);
  }

  public RegModel getRegModel() {
    return regModel;
  }

  public void setRegModel(RegModel regModel) {
    this.regModel = regModel;
  }

  public int getRowIndex() {
    return rowIndex;
  }

  public void setRowIndex(int rowIndex) {
    this.rowIndex = rowIndex;
  }

  public int[] getAttrXs() {
    return regModel.getAttrXs();
  }

  public int getAttrY() {
    return regModel.getAttrY();
  }

  @Override
  public int hashCode() {
    final int prime = 31;
    int result = 1;

    result = prime * result + regModel.hashCode();
    result = prime * result + rowIndex;

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
    LocalKey other = (LocalKey) obj;
    if (regModel.getAttrY() != other.getAttrY()) {
      return false;
    }
    if (regModel.getAttrXs().length != other.getAttrXs().length) {
      return false;
    }
    int[] attrXs = regModel.getAttrXs();
    for (int i = 0; i < attrXs.length; ++i) {
      if (attrXs[i] != other.getAttrXs()[i]) {
        return false;
      }
    }
    return true;
  }
}
