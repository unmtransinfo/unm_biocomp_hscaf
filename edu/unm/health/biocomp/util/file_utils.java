package edu.unm.health.biocomp.util;

import java.io.*;
import java.util.*;
import java.text.*;

/**	Static methods for file processing&#46;
	@author Jeremy J Yang
*/
public class file_utils
{
  private file_utils() {} //disable default constructor
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
}
