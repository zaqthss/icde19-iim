package cn.edu.thu.iim.entity;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by Aoqian Zhang on 2019/6/30.
 * E-mail address is zaqthss2009@gmail.com
 * All Rights Reserverd.
 *
 * @author Aoqian Zhang
 */
public class Database {
  private String head;
  private int attrNum;
  private List<Tuple> tpList;
  private int[][] flags; // 1 for observed, 0 for missing, n*attrNum
  private double[][] dbVals;

  private List<Integer> comRowIndexList = new ArrayList<>();
  private List<Integer> misRowIndexList = new ArrayList<>();
  private List<Cell> cells = new ArrayList<>();

  public Database() {
    // TODO Auto-generated constructor stub
    tpList = new ArrayList<>();
  }

  public Database(int attrNum,
      List<Tuple> tpList) {
    this.attrNum = attrNum;

    this.tpList = tpList;
    this.flags = new int[tpList.size()][attrNum];
  }

  public String getHead() {
    return head;
  }

  public void setHead(String head) {
    this.head = head;
  }

  public int getAttrNum() {
    return attrNum;
  }

  public void setAttrNum(int attrNum) {
    this.attrNum = attrNum;
  }

  public List<Tuple> getTpList() {
    return tpList;
  }

  public void setTpList(List<Tuple> tpList) {
    this.tpList = tpList;
  }

  public void addTuple(Tuple tp) {
    this.tpList.add(tp);
  }

  public List<Integer> getComRowIndexList() {
    return comRowIndexList;
  }

  public Tuple getTupleByIndex(int tIndex) {
    return tpList.get(tIndex);
  }

  public int[][] getFlags() {
    return flags;
  }

  public void setFlags(int[][] flags) {
    this.flags = flags;
  }

  // set flags according to cells
  public void setFlags() {
    int size = tpList.size();

    flags = new int[size][attrNum];

    for (int i = 0; i < size; ++i) {
      for (int j = 0; j < attrNum; ++j) {
        flags[i][j] = 1;
      }
    }

    for (Cell cell : cells) {
      Position pos = cell.getPosition();
      flags[pos.gettIndex()][pos.getAttrIndex()] = 0;
    }
  }

  public double[][] getDbVals() {
    return dbVals;
  }

  public void setDbVals(double[][] dbVals) {
    this.dbVals = dbVals;
  }

  public int getLength() {
    return tpList.size();
  }

  public void setComRowIndexList(List<Integer> comRowIndexList) {
    this.comRowIndexList = comRowIndexList;
  }

  public List<Integer> getMisRowIndexList() {
    return misRowIndexList;
  }

  public void setMisRowIndexList(List<Integer> misRowIndexList) {
    this.misRowIndexList = misRowIndexList;
  }

  public List<Cell> getCells() {
    return cells;
  }

  public void setCells(List<Cell> cells) {
    this.cells = cells;
  }
}
