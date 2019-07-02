package cn.edu.thu.iim.util;

import cn.edu.thu.iim.entity.KnnPair;
import java.util.Comparator;

/**
 * Created by Aoqian Zhang on 2019/6/30.
 * E-mail address is zaqthss2009@gmail.com
 * All Rights Reserved.
 *
 * Sort kNN
 *
 * @author Aoqian Zhang
 */
public class ComparatorKnnPair implements Comparator<KnnPair> {

  @Override
  public int compare(KnnPair kp1, KnnPair kp2) {
    // TODO Auto-generated method stub
    double distance1 = kp1.getDistance();
    double distance2 = kp2.getDistance();
    
    if (distance1 > distance2) {
      return 1;
    } else if (distance1 < distance2) {
      return -1;
    }
    
    return 0;
  }
  
}
