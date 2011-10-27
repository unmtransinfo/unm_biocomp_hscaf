package edu.unm.health.biocomp.hscaf;

import java.util.*;

import com.sleepycat.persist.model.*;
import static com.sleepycat.persist.model.Relationship.*;

/**	Scaffold entity used for storage in DB
	@see ScaffoldStore
	@author Oleg Ursu, Jeremy J Yang
*/
@Entity
public class ScaffoldEntity
{
  @PrimaryKey(sequence="ID")  //auto, 1+
  private int id;

  @SecondaryKey(relate=ONE_TO_ONE)
  private String canSmi;

  @SecondaryKey(relate=MANY_TO_MANY)
  private HashSet<Integer> childIds = new HashSet<Integer>();

  @SecondaryKey(relate=MANY_TO_ONE)
  private int parentId;

  public ScaffoldEntity(String canSmi,ArrayList<Integer> cIds) {
    this.canSmi = canSmi;
    for (int id: cIds) this.childIds.add(id);
  }

  @SuppressWarnings("unused")
  private ScaffoldEntity() { }

  public int getId() {
    return id;
  }

  public String getCanSmi() {
    return canSmi;
  }

  @Override
  public String toString() {
    return canSmi + "\t" + id ;
  }
}
