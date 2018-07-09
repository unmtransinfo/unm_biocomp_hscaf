package edu.unm.health.biocomp.hscaf;

import java.util.*;

import edu.unm.health.biocomp.db.*;

/**	Represents scaffold record used for storage in ScaffoldDB.
	Normally comprised of data from tables "scaffold" and "scaf2scaf".
        <br>
	@see ScaffoldDB
	@author Jeremy J Yang
*/
public class ScaffoldDBRecord
{
  private long id;
  private String cansmi;
  private HashSet<Long> childIDs = new HashSet<Long>();
  private long parentID;

  private ScaffoldDBRecord() { } //disallow default constructor

  public ScaffoldDBRecord(Long id,String cansmi) {
    this.id = id;
    this.cansmi = cansmi;
  }

  public long getParentID() { return this.parentID; }
  public long getID() { return id; }
  public String getCansmi() { return cansmi; }
  public HashSet<Long> getChildIDs() { return this.childIDs; }
  public void setParentID(Long id) { this.parentID=id; }
  public void setID(Long i) { this.id=i; }
  public void setCansmi(String c) { this.cansmi=c; }

  public void setChildIDs(ArrayList<Long> cids)
  {
    this.childIDs.clear();
    for (long id: cids) { this.childIDs.add(id); }
  }
  public void addChildID(Long id)
  {
    this.childIDs.add(id);
  }
  public String toString() { return cansmi + "\t" + id ; }
}
