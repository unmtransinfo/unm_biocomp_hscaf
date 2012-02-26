package edu.unm.health.biocomp.hscaf;

import java.io.*;
import java.util.*;

import chemaxon.formats.*;
import chemaxon.struc.*;
import chemaxon.sss.search.*;
import chemaxon.marvin.io.MolExportException;
 
/**	Represents a set of unique HierS scaffolds &#46;  Purpose of this class is
	to 
	<OL>
	<LI> Store set of scaffolds, compactly, generated for input molecular dataset, for output &#46;
	<LI> Optimize generation of scaffolds by avoiding repitition of prior analyses&#46;
	</OL>
	IDs are consecutive and start at 1&#46;
	Child scaffolds are not separate instances but rather among the set and may be shared&#46;
	Developer note: scafhash and idhash contain references to the same Scaffold objects, 
	thereby conserving memory&#46;
	<br />
	Note that the ScaffoldStore class is an alternative utilizing BerkeleyDB for
	improved performance and memory management for large datasets&#46;
	<br />
	@see edu.unm.health.biocomp.hscaf.Scaffold
	@see edu.unm.health.biocomp.hscaf.ScaffoldTree
	@see edu.unm.health.biocomp.hscaf.hier_scaffolds_utils
	@see edu.unm.health.biocomp.hscaf.hier_scaffolds
	@see edu.unm.health.biocomp.hscaf.hscaf_utils
	@author Jeremy J Yang
*/
public class ScaffoldSet
{
  private String name;
  private HashMap<String,Scaffold> scafhash;
  private HashMap<Long,Scaffold> idhash;
  /////////////////////////////////////////////////////////////////////////////
  /**	Constructor with name&#46;
  */
  public ScaffoldSet(String name)
  {
    this.name=name;
    this.scafhash = new HashMap<String,Scaffold>();
    this.idhash = new HashMap<Long,Scaffold>();
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Return number of Scaffolds in ScaffoldSet&#46;
  */
  public long getCount() { return this.scafhash.size(); }
  /////////////////////////////////////////////////////////////////////////////
  /**	True if Scaffold present with same cansmi&#46;
  */
  public boolean containsScaffold(Scaffold scaf)
  {
    return (this.scafhash.containsKey(scaf.getCansmi()));
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Returns ID of Scaffold with same cansmi if present, null otherwise&#46;
  */
  public Long getScaffoldID(Scaffold scaf)
  {
    return this.scafhash.get(scaf.getCansmi()).getID(); //null if missing
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Merges a complete [sub]ScaffoldTree defined by the Scaffold argument with the
	ScaffoldSet&#46;  New scaffolds are assigned new IDs, which is the main purpose of
	this method&#46; Existing scaffolds are recognized as such, should have IDs,
	since found during hscaf_utils.findChildScaffolds&#46;
	Initially called with new root scaffold, then
	recursively&#46; For each call,
	process all immediate child scaffolds, then recurse&#46;
	Child scaffolds present in set are same instances, not only equal&#46;
	Returns number of new scaffolds merged&#46;
  */
  public int mergeScaffoldTree(Scaffold scaf)
  {
//    //System.err.println("DEBUG: (mergeScaffoldTree) scaf: "+scaf.getCansmi());
//    int n_new=0;
//    long id=0L;
//    if (this.containsScaffold(scaf)) 
//    {
//      id=this.getScaffoldID(scaf);
//      System.err.println("DEBUG: (mergeScaffoldTree) scaf ID="+id+" contained in set.");
//      Scaffold foundscaf=this.getScaffoldByID(id);
//      scaf.setID(id); // should already be set?  But no...
//      for (Scaffold cscaf: foundscaf.getChildScaffolds())
//      {
//        System.err.println("DEBUG: (mergeScaffoldTree) addChild(cscaf) cID="+cscaf.getID());
//        scaf.addChild(cscaf);
//      }
//    }
//    else
//    {
//      id=this.getCount()+1;  //new ID [1...)
//      System.err.println("DEBUG: (mergeScaffoldTree) new scaf ID="+id+" cansmi: "+scaf.getCansmi());
//      scaf.setID(id);
//      this.addScaffold(scaf);
//      ++n_new;
//      for (Scaffold cscaf: scaf.getChildScaffolds())
//      {
//        n_new+=mergeScaffoldTree(cscaf); //recurse
//      }
//    }
//    return n_new;


// EXPERIMENT...


    //System.err.println("DEBUG: (mergeScaffoldTree) scaf: "+scaf.getCansmi());
    int n_new=0;
    long id=0L;
    if (this.containsScaffold(scaf)) 
    {
      id=this.getScaffoldID(scaf);
      //System.err.println("DEBUG: (mergeScaffoldTree) scaf ID="+id+" contained in set.");
      Scaffold foundscaf=this.getScaffoldByID(id);
      //if (foundscaf==null) System.err.println("DEBUG: (mergeScaffoldTree) foundscaf==null.");
      for (Scaffold cscaf: foundscaf.getChildScaffolds())
      {
        //System.err.println("DEBUG: (mergeScaffoldTree) addChild(cscaf) cID="+cscaf.getID());
        scaf.addChild(cscaf);
      }
      scaf.setID(id);
    }
    else
    {
      id=this.getCount()+1;  //new ID [1...)
      //System.err.println("DEBUG: (mergeScaffoldTree) new scaf ID="+id+" cansmi: "+scaf.getCansmi());
      scaf.setID(id);
      this.addScaffold(scaf);
      ++n_new;
    }

    for (Scaffold cscaf: scaf.getChildScaffolds())
    {
      long cid=0L;
      if (this.containsScaffold(cscaf)) 
      {
        cid=this.getScaffoldID(cscaf);
        cscaf.setID(cid);
      }
      else
      {
        cid=this.getCount()+1;  //new ID [1...)
        //System.err.println("DEBUG: (mergeScaffoldTree) new scaf cID="+cid+" cansmi: "+cscaf.getCansmi());
        cscaf.setID(cid);
        this.addScaffold(cscaf);
        ++n_new;
      }
      n_new+=mergeScaffoldTree(cscaf); //recurse
    }
    return n_new;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Returns Set of cansmi strings&#46;
  */
  public Set<String> keySet()
  {
    return this.scafhash.keySet();
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Add scaffold to ScaffoldSet, and compress&#46;
  */
  public void addScaffold(Scaffold scaf)
  {
    String cansmi=scaf.getCansmi();
    Long id=scaf.getID();
    //if (cansmi==null) System.err.println("DEBUG: (addScaffold) ERROR: cansmi==null.");
    //if (id==null) System.err.println("DEBUG: (addScaffold) ERROR: id==null.");
    //if (id==0) System.err.println("DEBUG: (addScaffold) ERROR: id==0.");
    scaf.compress();
    this.scafhash.put(cansmi,scaf);
    this.idhash.put(id,scaf);
  }
  /////////////////////////////////////////////////////////////////////////////
  public Scaffold getScaffoldByID(long id)
  {
    Scaffold scaf=null;
    if (this.idhash.containsKey(id))
      scaf=this.idhash.get(id);
    return scaf;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Return Scaffold with same cansmi, if present, otherwise null&#46;
  */
  public Scaffold findScaffold(Scaffold scaf)
  {
    return this.scafhash.get(scaf.getCansmi()); //null if missing
  }
  /////////////////////////////////////////////////////////////////////////////
  public void DEBUG_dump()
  {
    for (long id=1L;true;++id)
    {
      Scaffold scaf=null;
      try {
        scaf=this.getScaffoldByID(id);
      }
      catch (Exception e) { break; }
      if (scaf==null) break;
      System.err.println("Scaffold ID="+id+":");
      System.err.println("\tCansmi: "+scaf.getCansmi());
      System.err.println("\tChildIDs: "+scaf.getChildIDs());
    }
    return;
  }
}
