package cn.edu.thu.iim.test;

import cn.edu.thu.iim.IIM;
import cn.edu.thu.iim.entity.Cell;
import cn.edu.thu.iim.entity.Database;
import cn.edu.thu.iim.entity.Position;
import cn.edu.thu.iim.util.DataUtil;
import java.util.List;
import java.util.Map;

/**
 * Created by Aoqian Zhang on 2019/6/28.
 * E-mail address is zaqthss2009@gmail.com
 * All Rights Reserved.
 *
 * @author Aoqian Zhang
 */

public class IIMTest {
  public static void main(String[] args) {

    String rawInput = "asf1_0.1miss.csv";
    int[] lparams = {1, 101, 5}; // LBEGIN, LMAX, INTERVAL
    int K = 10;

    // run IIM
    Database rawDb = DataUtil.readData(rawInput);
    IIM iim = new IIM(rawDb);
    iim.setParams(lparams, K);
    Map<Position, Cell> repairedCells = iim.mainIIM(true, false, true);

    // evaluate
    String truthInput = "asf.csv";
    Database truthDb = DataUtil.readData(truthInput);
    Map<Position, Cell> truthCells = DataUtil.getTruth(repairedCells, truthDb);

    double rmse = DataUtil.calcRMS(truthCells, repairedCells);

    System.out.println("rms error is " + rmse);
    System.out.println("test over");
  }
}
