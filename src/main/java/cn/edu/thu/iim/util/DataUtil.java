package cn.edu.thu.iim.util;

import cn.edu.thu.iim.entity.Cell;
import cn.edu.thu.iim.entity.Database;
import cn.edu.thu.iim.entity.Position;
import cn.edu.thu.iim.entity.RegModel;
import cn.edu.thu.iim.entity.Tuple;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

/**
 * Created by Aoqian Zhang on 2019/6/30.
 * E-mail address is zaqthss2009@gmail.com
 * All Rights Reserved.
 *
 * Handle data
 *
 * @author Aoqian Zhang
 */
public class DataUtil {
  public static String PATH = "data/";

  /**
   * Read Data from file to get Database
   * @param input
   * @return
   */
  public static Database readData(String input) {
    Database db = new Database();
    List<Integer> misRowIndexList = new ArrayList<>();
    List<Cell> cells = new ArrayList<>();

    try {
      FileReader fr = new FileReader(DataUtil.PATH + input);
      BufferedReader br = new BufferedReader(fr);

      String line;
      String[] vals;

      String head = br.readLine();
      db.setHead(head);
      int attrNum = head.split(",").length;
      db.setAttrNum(attrNum);

      int tid = 0;
      double value;
      while ((line = br.readLine()) != null) {
        vals = line.split(",");

        double[] data = new double[attrNum];
        for (int j = 0; j < attrNum; ++j) {
          // if the last one is empty, then vals will shrink
          if (j >= vals.length || vals[j].equals("")) {
            value = Double.NaN;

            Position pos = new Position(tid, j);
            Cell cell = new Cell(pos);
            cells.add(cell);

            if (!misRowIndexList.contains(tid)) {
              misRowIndexList.add(tid);
            }
          } else {
            value = Double.parseDouble(vals[j]);
          }
          data[j] = value;
        }

        Tuple tp = new Tuple(attrNum);
        tp.buildTuple(tid, data);
        db.addTuple(tp);
        tid++;
      }

      db.setMisRowIndexList(misRowIndexList);
      db.setCells(cells);
      db.setFlags();

      br.close();
      fr.close();

    } catch (IOException e) {
      e.printStackTrace();
    }

    return db;
  }

  public static List<RegModel> findMissingRegModels(List<Integer> misRowIndexList,
      Map<Integer, List<Integer>> misRowAttrIndexMap, int attrNum) {
    List<RegModel> models = new ArrayList<>();

    RegModel regModel;
    for (int misRowIndex : misRowIndexList) {
      List<Integer> misList = misRowAttrIndexMap.get(misRowIndex);

      int[] attrXs = getAttrXsFromMisList(attrNum, misList);
      for (int mi = 0; mi < misList.size(); ++mi) {
        int attrY = misList.get(mi);
        regModel = new RegModel(attrXs, attrY);

        if (!models.contains(regModel)) {
          models.add(regModel);
        }
      }
    } // end of i

    return models;
  }

  public static int[] getAttrXsFromMisList(int attrNum, List<Integer> misList) {
    int attrXNum = attrNum - misList.size();

    int[] attrXs = new int[attrXNum];
    int curIndex = 0;

    for (int attri = 0; attri < attrNum; ++attri) {
      if (misList.contains(attri)) {
        continue;
      }
      attrXs[curIndex++] = attri;
    }
    return attrXs;
  }

  public static Map<Position, Cell> getTruth(Map<Position, Cell> repairedCells, Database truthDb) {
    Map<Position, Cell> truthCells = new HashMap<>();
    for (Entry<Position, Cell> entry : repairedCells.entrySet()) {
      Position pos = entry.getKey();

      int rowIndex = pos.gettIndex();
      int attrIndex = pos.getAttrIndex();

      double truthVal = truthDb.getTupleByIndex(rowIndex).getDataByIndex(attrIndex);
      Cell cell = new Cell(pos);
      cell.setValue(truthVal);
      truthCells.put(pos, cell);
    }

    return truthCells;
  }

  /**
   * RMS sqrt(|modify - truth|^2 / len)
   */
  public static double calcRMS(Map<Position, Cell> truthCells, Map<Position, Cell> repairedCells) {
    double cost = 0;
    double delta;

    int len = truthCells.size();
    for (Position pos: truthCells.keySet()) {
      Cell truthCell = truthCells.get(pos);
      Cell repairedCell = repairedCells.get(pos);

      delta = truthCell.getValue() - repairedCell.getValue();

      cost += delta * delta;
    }
    cost /= len;

    return Math.sqrt(cost);
  }
}
