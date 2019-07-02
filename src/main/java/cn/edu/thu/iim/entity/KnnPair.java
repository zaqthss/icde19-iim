package cn.edu.thu.iim.entity;

/**
 * Created by Aoqian Zhang on 2019/6/28.
 * E-mail address is zaqthss2009@gmail.com
 * All Rights Reserved.
 *
 * @author Aoqian Zhang
 */
public class KnnPair {
  private double distance;
  private int index;
  
  public KnnPair(double distance, int index) {
    setDistance(distance);
    setIndex(index);
  }

  public double getDistance() {
    return distance;
  }

  public void setDistance(double distance) {
    this.distance = distance;
  }

  public int getIndex() {
    return index;
  }

  public void setIndex(int index) {
    this.index = index;
  }
}
