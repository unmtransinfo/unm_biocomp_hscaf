package edu.unm.health.biocomp.hscaf;

import java.io.*;
import java.util.*;
import java.text.*;

/**	Static utility methods for time and date processing, etc&#46;
	@author Jeremy J Yang
*/
public class util
{
  private util() {} //disable default constructor
  /////////////////////////////////////////////////////////////////////////////
  /**	Human readable file size&#46;
  */
  public static String niceBytes(long bytes)
  {
    if (bytes<1.0e3)
      return String.format("%d bytes",bytes);
    else if (bytes<1.0e6)
      return String.format("%.1f KB",((float)bytes/1.0e3));
    else if (bytes<1.0e9)
      return String.format("%.1f MB",((float)bytes/1.0e6));
    else if (bytes<1.0e12)
      return String.format("%.1f GB",((double)bytes/1.0e9));
    else
      return String.format("%.1f TB",((double)bytes/1.0e12));
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Human readable time interval&#46;
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
}
