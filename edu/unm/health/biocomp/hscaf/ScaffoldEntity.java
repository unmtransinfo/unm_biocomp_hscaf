package edu.unm.health.biocomp.hscaf;

import java.util.*;

import com.sleepycat.persist.model.*;
import static com.sleepycat.persist.model.Relationship.*;

/**	Scaffold entity used for storage in BerkeleyDB-based ScaffoldStore.
        <br />
	Developer note: Auto-sequencing did not work correctly for cases where a database
	was re-opened, which resulted in skipping ahead (to nearest largest multiple of 100).
	So scaffold IDs are manually assigned 1 + count().
        <br />
	@see ScaffoldStore
	@author Jeremy J Yang
*/
@Entity
public class ScaffoldEntity
{
  //@PrimaryKey(sequence="ID")  //auto, 1+	... not working correctly on re-opened database, skips ahead.

  @PrimaryKey	//must assign manually as next available count+1, after [1,2,...count]
  private long id;

  @SecondaryKey(relate=ONE_TO_ONE)
  private String canSmi;

  @SecondaryKey(relate=MANY_TO_MANY)
  private HashSet<Long> childIds = new HashSet<Long>();

  @SecondaryKey(relate=MANY_TO_ONE)
  private long parentId;

  public ScaffoldEntity(Long id,String canSmi) {
    this.id = id;
    this.canSmi = canSmi;
  }

  public void setChildIds(ArrayList<Long> cIds)
  {
    this.childIds.clear();
    for (long id: cIds)
    {
      this.childIds.add(id);
    }
  }

  public void addChildId(Long Id)
  {
    //System.err.println("DEBUG: (ScaffoldEntity.addChildId) id="+Id);
    this.childIds.add(Id);
  }

  public void setParentId(Long Id) {
    this.parentId=Id;
  }

  public long getParentId() {
    return this.parentId;
  }

  @SuppressWarnings("unused")
  private ScaffoldEntity() { }

  public long getId() {
    return id;
  }

  public String getCanSmi() {
    return canSmi;
  }

  @Override
  public String toString()
  {
    return canSmi + "\t" + id ;
  }

  public HashSet<Long> getChildIds()
  {
    return this.childIds;
  }
}
