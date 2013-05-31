package edu.unm.health.biocomp.http;
import java.util.*;

/**	HTTP parameter container.
	@author Jeremy J Yang
*/
public class HttpParams
{
  private HashMap<String,String> params;
  public HttpParams()
  {
    params = new HashMap<String,String>();
  }
  public String getVal(String key)
  {
    if (!params.containsKey(key)) return "";
    return params.get(key);
  }
  public void setVal(String key,String val)
  {
    params.put(key,val);
  }
  public boolean isChecked(String key)
  {
    if (!params.containsKey(key)) return false;
    if (params.get(key).equals("")) return false;
    return true;
  }
  public boolean containsKey(String key)
  {
    return params.containsKey(key);
  }
  public void clear()
  {
    params.clear();
  }
}
