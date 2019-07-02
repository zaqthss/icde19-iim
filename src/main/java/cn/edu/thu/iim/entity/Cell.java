package cn.edu.thu.iim.entity;

/**
 * Created by Aoqian Zhang on 2019/6/30.
 * E-mail address is zaqthss2009@gmail.com
 * All Rights Reserverd.
 *
 * @author Aoqian Zhang
 */
public class Cell {
  private Position position;
  private double value;

  public Cell(Position position) {
    setPosition(position);
  }

  public Position getPosition() {
    return position;
  }

  public void setPosition(Position position) {
    this.position = position;
  }

  public double getValue() {
    return value;
  }

  public void setValue(double value) {
    this.value = value;
  }
}
