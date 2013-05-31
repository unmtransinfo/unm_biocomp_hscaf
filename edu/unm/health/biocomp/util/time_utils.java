package edu.unm.health.biocomp.util;

import java.io.*;
import java.util.*;
import java.text.*;

/**	Static methods for time and date processing.
	@author Jeremy J Yang
*/
public class time_utils
{
  private time_utils() {} //disable default constructor
  /////////////////////////////////////////////////////////////////////////////
  /**	Human readable time interval.
  */
  public static String timeDeltaStr(java.util.Date t_i,java.util.Date t_f)
  {
    long t_d=t_f.getTime()-t_i.getTime();
    int t_d_s = (int)((t_d/1000L)%60L);
    int t_d_m = (int)((t_d/60000L)%60);
    int t_d_h = (int)(t_d/3600000L);
    if (t_d_h>0)
      return (String.format("%02d:%02d:%02d",t_d_h,t_d_m,t_d_s));
    else if (t_d_m>0)
      return (String.format("%02d:%02d",t_d_m,t_d_s));
    else
      return (String.format("%2ds",t_d_s));
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Time of day (00:00 - 23:59)..
  */
  public static String currentTime()
  {
    Calendar cal = Calendar.getInstance(TimeZone.getDefault());
    return String.format("%02d:%02d",cal.get(Calendar.HOUR_OF_DAY),cal.get(Calendar.MINUTE));
  }
}
