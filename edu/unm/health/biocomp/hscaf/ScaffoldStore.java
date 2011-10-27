package edu.unm.health.biocomp.hscaf;

import java.io.File;

import com.sleepycat.je.CursorConfig;
import com.sleepycat.je.Environment;
import com.sleepycat.je.EnvironmentConfig;
import com.sleepycat.je.Transaction;
import com.sleepycat.persist.EntityCursor;
import com.sleepycat.persist.EntityStore;
import com.sleepycat.persist.PrimaryIndex;
import com.sleepycat.persist.SecondaryIndex;
import com.sleepycat.persist.StoreConfig;
import com.sleepycat.je.DatabaseException;

/**	Initiate and provide local storage &amp; retrieval via BerkeleyDB&#46;
	For large jobs this avoids memory limits and improves performance by
	storing solved scaffolds for rapid lookup and avoiding re-calculation&#46;
	Algorithm logic is unaffected&#46;
	@author Oleg Ursu, Jeremy J Yang
 */
public class ScaffoldStore
{
  private EnvironmentConfig envConf = new EnvironmentConfig();
  private Environment env;
  private StoreConfig storeConf = new StoreConfig();
  private EntityStore scafEntityStore;
  
  protected PrimaryIndex<Integer, ScaffoldEntity> scaffoldById;
  protected SecondaryIndex<String, Integer, ScaffoldEntity> scaffoldByCanSmi;
  
  public ScaffoldStore(File location)
	throws DatabaseException
  {
    assert (!(location.exists() && location.isFile()));
    this.envConf.setAllowCreate(true);
    this.envConf.setTransactional(true);
    this.env = new Environment(location, envConf);
    this.storeConf.setAllowCreate(true);
    this.storeConf.setTransactional(true);
    this.scafEntityStore = new EntityStore(env, "scaffolds", storeConf);
    this.scaffoldById = this.scafEntityStore.getPrimaryIndex(Integer.class, ScaffoldEntity.class);
    this.scaffoldByCanSmi = this.scafEntityStore.getSecondaryIndex(scaffoldById, String.class, "canSmi");
  }
  
  public void closeAll()
	throws DatabaseException
  {
    this.scafEntityStore.close();
    this.env.cleanLog();
    this.env.close();
  }
  
  public void removeAll()
	throws DatabaseException
  {
    removeAllFromScaffStore();
  }
  
  public void removeAllFromScaffStore()
	throws DatabaseException
  {
    Transaction txn = env.beginTransaction(null, null);
    EntityCursor<ScaffoldEntity> curr = scaffoldById.entities(txn, CursorConfig.DEFAULT);
    for(ScaffoldEntity scaf = curr.first(); scaf != null; scaf = curr.next()) {
      curr.delete();
    }
    curr.close();
    txn.commit();
  }

  public int mergeScaffoldTree(ScaffoldTree scaftree)
    throws DatabaseException
  {
    int n_new=0;
    for (Scaffold scaf : scaftree.getAllScaffolds())
    {
      String canSmi = scaf.getCansmi();
      int id=0;
      if (this.scaffoldByCanSmi.contains(canSmi))
      {
        id=this.scaffoldByCanSmi.get(canSmi).getId();
        scaf.setID(id);
      }
      else
      {
        ScaffoldEntity entity = new ScaffoldEntity(canSmi,scaf.getChildIDs());
        this.scaffoldById.put(entity);
        id=entity.getId();
        scaf.setID(id);
        ++n_new;
      }
    }
    return n_new;
  }
  
}
