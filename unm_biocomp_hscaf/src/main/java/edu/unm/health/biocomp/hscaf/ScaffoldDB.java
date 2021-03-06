package edu.unm.health.biocomp.hscaf;

import java.io.*;
import java.util.*;
import java.math.*;
import java.sql.*;

import edu.unm.health.biocomp.util.db.*;

import chemaxon.formats.MolFormatException;

/**	Provides persistent scaffold storage via standard RDB.
	For large jobs this avoids memory limits and improves performance by
	storing solved scaffolds for rapid lookup and avoiding re-calculation.
	<br>
	An additional benefit of DB usage is paralellization, since
	separate processes can concurrently read and write to the same
	scaffold database by locking and unlocking records.
	<br>
	Minimal, standard RDB features required; no chemical cartridge used,
	smiles canonicalization by client.
	<br>
	PostgreSQL supported only.  Maybe should use another layer to be database
	agnostic (support MySql, etc).
	<br>
	@author Jeremy J Yang
	@see	ScaffoldDBRecord
 */
public class ScaffoldDB
{
  private Boolean stereo;
  private Boolean keep_nitro_attachments;
  private String dbhost;
  private int dbport;
  private String dbname;
  private String dbschema;
  private String dbusr;
  private String dbpw;
  private String dbtableprefix;
  private Connection dbcon;
  
  /**	@return host for database
  */
  public String getDBHost() { return this.dbhost; }
  /**	@return name for database
  */
  public String getDBName() { return this.dbname; }
  /**	@return schema for database
  */
  public String getDBSchema() { return this.dbschema; }
  /**	@return user for database connection
  */
  public String getDBUser() { return this.dbusr; }
  /**	@return port for database connection
  */
  public int getDBPort() { return this.dbport; }
  /**	@return table prefix for database
  */
  public String getDBTableprefix() { return this.dbtableprefix; }
  /**	@return Connection for database
  */
  public Connection getDBConnection() { return this.dbcon; }

  /////////////////////////////////////////////////////////////////////////////
  /**	Constructor creates or opens existing database at location specified.<br/>
	<p>
	<code>CREATE SCHEMA IF NOT EXISTS <i>DBSCHEMA</i><br/>
	Schema:<br>
	<code>CREATE TABLE <i>DBSCHEMA.DBTABLEPREFIX</i>scaffold (<br/>
	&nbsp;	id INTEGER PRIMARY KEY,<br/>
	&nbsp;	scafsmi VARCHAR(512) NOT NULL,<br/>
	&nbsp;	scaftree VARCHAR(2048))<br/>
	<br>
	CREATE TABLE <i>DBSCHEMA.DBTABLEPREFIX</i>scaf2scaf (<br/>
	&nbsp;	parent_id INTEGER,<br/>
	&nbsp;	child_id INTEGER )<br/>
	</code>

	@param	host	database host
	@param	port	database port
	@param	name	database name
	@param	schema	database schema
	@param	usr	database connection user
	@param	pw	database connection password
	@param	tableprefix	database table name prefix
	@param	stereo	flag, whether scaffolds are stereo
	@param	keep_nitro_attachments	flag, whether to keep nitrogens on rings
	@param	create_new	flag indicating whether to create new database
  */
  public ScaffoldDB(
	String host,
	int port,
	String name,
	String schema,
	String usr,
	String pw,
	String tableprefix,
	Boolean stereo,
	Boolean keep_nitro_attachments,
	Boolean create_new)
	throws SQLException
  {
    this.dbhost=host;
    this.dbport=port;
    this.dbname=name;
    this.dbschema=schema;
    this.dbusr=usr;
    this.dbpw=pw;
    this.dbtableprefix=tableprefix;
    this.stereo=stereo;
    this.keep_nitro_attachments=keep_nitro_attachments;

    this.dbcon=pg_utils.DBConnect(this.dbhost,this.dbport,this.dbname,this.dbusr,this.dbpw);
    //System.err.println("DEBUG: "+pg_utils.ServerStatusTxt(this.dbcon));

    boolean ok;
    if (create_new)
    {
       createNewDB();
    }
    ok=checkValid();
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Creates new database at location specified.
	Scafsmi column intentionally not specified UNIQUE, to allow re-canonicalization
	by database cartridge, which could result in collisions (unfortunate
	but better than a fatal duplicate-key error).
  */
  private boolean createNewDB()
	throws SQLException
  {
    boolean ok=true;
    String sql="";

    this.dbcon.setAutoCommit(true);  //Should be true by default but let's be sure.

    //sql="CREATE SCHEMA IF NOT EXISTS "+this.dbschema; //should work in PG 9.3
    try {ok&=pg_utils.Execute(this.dbcon,"CREATE SCHEMA "+this.dbschema);} catch (Exception e) {} //Ignore; assume schema exists

    sql="CREATE TABLE "+this.dbschema+"."+this.dbtableprefix+"scaffold"
      +"(id INTEGER PRIMARY KEY,scafsmi VARCHAR(512) NOT NULL,scaftree VARCHAR(2048))";
    ok&=pg_utils.Execute(this.dbcon,sql);

    sql="CREATE TABLE "+this.dbschema+"."+this.dbtableprefix+"scaf2scaf"
      +"(parent_id INTEGER NOT NULL,child_id INTEGER NOT NULL)";
    ok&=pg_utils.Execute(this.dbcon,sql);

    sql="CREATE INDEX "+this.dbtableprefix+"_scaffold_id_idx ON "+this.dbschema+"."+this.dbtableprefix+"scaffold (id)";
    ok&=pg_utils.Execute(this.dbcon,sql);

    sql="CREATE INDEX "+this.dbtableprefix+"_scaf2scaf_parent_id_idx ON "+this.dbschema+"."+this.dbtableprefix+"scaf2scaf (parent_id)";
    ok&=pg_utils.Execute(this.dbcon,sql);

    return ok;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Re-index modified tables periodically for faster response.
  */
  public boolean reindex()
	throws SQLException
  {
    boolean ok=true;
    if (this.dbcon==null) return false;
    String sql="REINDEX TABLE "+this.dbschema+"."+this.dbtableprefix+"scaffold";
    ok&=pg_utils.Execute(this.dbcon,sql);

    sql="REINDEX TABLE "+this.dbschema+"."+this.dbtableprefix+"scaf2scaf";
    ok&=pg_utils.Execute(this.dbcon,sql);
    return ok;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Minimal sanity check whether DB exists, is accessible and valid schema.
  */
  public boolean checkValid()
	throws SQLException
  {
    boolean ok=true;
    if (this.dbcon==null) return false;
    String sql="SELECT COUNT(scafsmi) FROM "+this.dbschema+"."+this.dbtableprefix+"scaffold";
    ResultSet rset=pg_utils.ExecuteSql(this.dbcon,sql); // Check ok?  How?
    ok&=(rset.next());
    rset.getStatement().close();

    sql="SELECT COUNT(parent_id) FROM "+this.dbschema+"."+this.dbtableprefix+"scaf2scaf";
    rset=pg_utils.ExecuteSql(this.dbcon,sql); // Check ok?  How?
    ok&=(rset.next());
    rset.getStatement().close();
    return ok;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	@return number of scaffolds in ScaffoldDB
  */
  public long count()
	throws SQLException
  {
    long n=0;
    String sql="SELECT COUNT(scafsmi) FROM "+this.dbschema+"."+this.dbtableprefix+"scaffold";
    ResultSet rset=pg_utils.ExecuteSql(this.dbcon,sql); // Check ok?  How?
    if (rset.next())
      n=rset.getInt(1);
    rset.getStatement().close();
    return n;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	@return information about contents of ScaffoldDB
  */
  public String info()
	throws SQLException
  {
    String txt="";
    txt+=(this.getDBUser()+"@");
    txt+=(this.getDBHost()+":"+this.getDBPort()+":"+this.getDBName()+":"+this.getDBSchema());
    txt+=(":("+this.getDBTableprefix()+"scaffold,"+this.getDBTableprefix()+"scaf2scaf)");
    txt+=("[n_scaf="+count()+"]\n");
    return txt;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	@param	id	scaffold ID
	@return	{@link ScaffoldDBRecord} for specified scaffold ID 
  */
  public ScaffoldDBRecord getScaffoldByID(long id)
	throws SQLException
  {
    ScaffoldDBRecord scafrec = null;
    String sql="SELECT scafsmi FROM "+this.dbschema+"."+this.dbtableprefix+"scaffold WHERE id="+id;
    ResultSet rset=pg_utils.ExecuteSql(this.dbcon,sql); // Check ok?  How?
    if (!rset.next()) return null; //ERROR

    scafrec = new ScaffoldDBRecord(id,rset.getString("scafsmi")); 
    rset.getStatement().close();

    sql="SELECT child_id FROM "+this.dbschema+"."+this.dbtableprefix+"scaf2scaf WHERE parent_id="+id;
    rset=pg_utils.ExecuteSql(this.dbcon,sql); // Check ok?  How?
    while (rset.next())
    {
      scafrec.addChildID(rset.getLong("child_id"));
    }
    return scafrec;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	@param	smi	scaffold SMILES
	@return	{@link ScaffoldDBRecord} for specified smiles
  */
  public ScaffoldDBRecord getScaffoldBySmiles(String smi)
	throws SQLException,MolFormatException,IOException
  {
    String cansmi=hscaf_utils.cansmi(smi,this.isStereo());
    return getScaffoldByCansmi(cansmi);
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	@param	cansmi	scaffold SMILES, already canonical
	@return	{@link ScaffoldDBRecord} for specified smiles
  */
  public ScaffoldDBRecord getScaffoldByCansmi(String cansmi)
	throws SQLException,MolFormatException,IOException
  {
    String sql="SELECT id FROM "+this.dbschema+"."+this.dbtableprefix+
	"scaffold WHERE scafsmi='"+pg_utils.QuoteString(cansmi)+"'";
    ResultSet rset=pg_utils.ExecuteSql(this.dbcon,sql); // Check ok?  How?
    if (!rset.next()) return null; //NOT FOUND
    return getScaffoldByID(rset.getLong("id"));
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Check whether specified ID is in database.
  */
  public boolean containsScaffoldByID(Long id)
	throws SQLException
  {
    String sql="SELECT scafsmi FROM "+this.dbschema+"."+this.dbtableprefix+"scaffold"+" WHERE id="+id;
    ResultSet rset=pg_utils.ExecuteSql(this.dbcon,sql);
    boolean ok=rset.next();
    rset.getStatement().close();
    return ok;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Check whether scaffold specified by smiles is in database.
	@param	smi	scaffold SMILES
  */
  public boolean containsScaffoldBySmiles(String smi)
	throws SQLException,MolFormatException,IOException
  {
    String cansmi=hscaf_utils.cansmi(smi,this.isStereo());
    return containsScaffoldByCansmi(cansmi);
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Check whether scaffold specified by canonical smiles is in database.
	@param	cansmi	scaffold canonical SMILES
  */
  public boolean containsScaffoldByCansmi(String cansmi)
	throws SQLException
  {
    String sql="SELECT id FROM "+this.dbschema+"."+this.dbtableprefix+
	"scaffold WHERE scafsmi='"+pg_utils.QuoteString(cansmi)+"'";
    ResultSet rset=pg_utils.ExecuteSql(this.dbcon,sql);
    boolean ok=rset.next();
    rset.getStatement().close();
    return ok;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	@param	cansmi	scaffold canonical SMILES
	@return scaffold ID specified by canonical smiles
  */
  public long scaffoldCansmi2ID(String cansmi)
	throws SQLException
  {
    long id=0;
    String sql="SELECT id FROM "+this.dbschema+"."+this.dbtableprefix+
	"scaffold WHERE scafsmi='"+pg_utils.QuoteString(cansmi)+"'";
    ResultSet rset=pg_utils.ExecuteSql(this.dbcon,sql);
    boolean ok=rset.next();
    if (!ok) return 0; //NOT FOUND
    id=rset.getLong("id");
    rset.getStatement().close();
    return id;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Write {@link ScaffoldDBRecord} to database.
	Should assure scaffold not already known.
  */
  public boolean addScaffold(ScaffoldDBRecord scafrec)
	throws SQLException
  {
    boolean ok=true;
    this.dbcon.setAutoCommit(true);  //Should be true by default but let's be sure.
    String sql="INSERT INTO "+this.dbschema+"."+this.dbtableprefix+"scaffold (id,scafsmi,scaftree)"
	+" VALUES ("+scafrec.getID()+",'"+scafrec.getCansmi()+"','"+getScaffoldString(scafrec)+"')";
    pg_utils.Execute(this.dbcon,sql);
    return ok;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Update database with revised {@link ScaffoldDBRecord},
	matched by ID.
  	@param scafrec	 revised {@link ScaffoldDBRecord} with pre-existing ID
  */
  public boolean updateScaffoldByID(ScaffoldDBRecord scafrec)
	throws SQLException
  {
    boolean ok=true;
    this.dbcon.setAutoCommit(true);  //Should be true by default but let's be sure.
    String sql="UPDATE "+this.dbschema+"."+this.dbtableprefix+"scaffold"
	+" SET scafsmi='"+scafrec.getCansmi()+"',scaftree='"+getScaffoldString(scafrec)+"'"
	+" WHERE id="+scafrec.getID();
    pg_utils.Execute(this.dbcon,sql); // Check ok?  How?

    sql="DELETE FROM "+this.dbschema+"."+this.dbtableprefix+"scaf2scaf"
	+" WHERE parent_id="+scafrec.getID() ;
    pg_utils.Execute(this.dbcon,sql);

    HashSet<Long> chids = scafrec.getChildIDs();
    for (Iterator<Long> itr=chids.iterator(); itr.hasNext(); )
    {
      long chid=itr.next();
      sql="INSERT INTO "+this.dbschema+"."+this.dbtableprefix+"scaf2scaf"
	+" (parent_id,child_id)  VALUES ("+scafrec.getID()+","+chid+")" ;
      pg_utils.Execute(this.dbcon,sql); // Check ok?  How?
    }

    return ok;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Close database.
  */
  public void close()
	throws SQLException
  {
    this.dbcon.close();
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Delete all data from scaffold databases.  Tables and schema remain.
  */
  public void removeAll()
	throws SQLException
  {
    this.dbcon.setAutoCommit(true);  //Should be true by default but let's be sure.
    String sql="DELETE FROM "+this.dbschema+"."+this.dbtableprefix+"scaffold";
    pg_utils.Execute(this.dbcon,sql); // Check ok?  How?

    sql="DELETE FROM "+this.dbschema+"."+this.dbtableprefix+"scaf2scaf";
    pg_utils.Execute(this.dbcon,sql); // Check ok?  How?
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Destroy databases.  Tables are dropped (not schema).
  */
  public void destroyDB()
	throws SQLException
  {
    this.dbcon.setAutoCommit(true);  //Should be true by default but let's be sure.
    String sql="DROP TABLE "+this.dbschema+"."+this.dbtableprefix+"scaffold CASCADE";
    pg_utils.Execute(this.dbcon,sql); // Check ok?  How?
    sql="DROP TABLE "+this.dbschema+"."+this.dbtableprefix+"scaf2scaf CASCADE";
    pg_utils.Execute(this.dbcon,sql); // Check ok?  How?
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Merges a fully populated ScaffoldTree defined by the Scaffold
	argument with the ScaffoldDB.  

	New scaffolds are assigned new IDs.  Existing scaffolds are
	recognized as such.  Initially called with new root scaffold,
	then recursively.  For each call, process all immediate child
	scaffolds, then assign ChildIDs and ParentID to ScaffoldDBRecord.

	Note that a new unknown parent may have known and/or new child scaffolds.	
	Hence this recursive method must handle the known-scaffold case.
  */
  public int mergeScaffoldTree(Scaffold scaf)
	throws SQLException
  {
    int n_new=0;
    String cansmi=scaf.getCansmi();
    //System.err.println("DEBUG: (mergeScaffoldTree) scafsmi="+cansmi);
    long id=0L;
    // First either find existing or create new ScaffoldDBRecord.
    if (this.containsScaffoldByCansmi(cansmi))
    {
      id=this.scaffoldCansmi2ID(cansmi);
      //System.err.println("DEBUG: (mergeScaffoldTree) scaf found; id="+id);
    }
    else
    {
      id=this.count()+1;	//new ID (consecutive)
      this.addScaffold(new ScaffoldDBRecord(id,cansmi));
      ++n_new;
    }
    scaf.setID(id);
    ScaffoldDBRecord scafrec=this.getScaffoldByID(id);

    // For each child, find existing or create new ScaffoldDBRecord, then recurse.
    for (Scaffold cscaf : scaf.getChildScaffolds())
    {
      cansmi=cscaf.getCansmi();
      long cid=0L;
      if (this.containsScaffoldByCansmi(cansmi))
      {
        cid=this.scaffoldCansmi2ID(cansmi);
      }
      else
      {
        cid=this.count()+1;	//new ID (consecutive)
        this.addScaffold(new ScaffoldDBRecord(cid,cansmi));
        ++n_new;
      }
      ScaffoldDBRecord cscafrec=getScaffoldByID(cid); //existing or just added
      cscaf.setID(cid);
      scafrec.addChildID(cid);
      this.updateScaffoldByID(scafrec); //update
      //System.err.println("DEBUG: (mergeScaffoldTree) ID="+id+" scafrec.getChildIDs() "+scafrec.getChildIDs());
      //System.err.println("DEBUG: (mergeScaffoldTree) CID="+cid+" cscafrec.setParentID("+id+") ...");
      cscafrec.setParentID(id);
      this.updateScaffoldByID(cscafrec); //update
      //System.err.println("DEBUG: (mergeScaffoldTree) CID="+cid+" ID="+id+" cscafrec.getParentID()"+cscafrec.getParentID());
      n_new+=mergeScaffoldTree(cscaf);
    }
    return n_new;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Populates an incomplete ScaffoldTree from known root Scaffold.
	Note that since root scaffold is known, all child scaffolds must also
	be known.

	For known scaffold, obtain ID from store, and find child scaffolds
	in store to populate ScaffoldTree.
  */
  public int populateScaffoldTree(Scaffold scaf)
	throws SQLException
  {
    int n_scaf=0;
    //System.err.println("DEBUG: (populateScaffoldTree); scaf.getID() = "+scaf.getID());
    if (!this.containsScaffoldByID(scaf.getID()))
    {
      //System.err.println("DEBUG: (populateScaffoldTree) aaack! ID not found...");
      return 0; //Should not happen.
    }
    ScaffoldDBRecord scafrec = this.getScaffoldByID(scaf.getID());
    HashSet<Long> chids = scafrec.getChildIDs();
    //System.err.println("DEBUG: (populateScaffoldTree) chids.size() = "+chids.size());
    for (Iterator<Long> itr=chids.iterator(); itr.hasNext(); )
    {
      long chid=itr.next();
      //if (chid==0) System.err.println("DEBUG: (populateScaffoldTree) ERROR; chid = "+chid);
      ScaffoldDBRecord cscafrec = this.getScaffoldByID(chid);
      String cansmi = cscafrec.getCansmi();
      //System.err.println("DEBUG: (populateScaffoldTree) chid = "+chid+" ; smi = "+cansmi);
      Scaffold cscaf = null;
      try { cscaf = new Scaffold(cansmi,this.isStereo(),this.isKeep_nitro_attachments()); }
      catch (MolFormatException e) { System.err.println("DEBUG: aaack! (populateScaffoldTree)"); } //Should not happen.
      cscaf.setID(chid);
      scaf.addChild(cscaf);
      ++n_scaf;
      n_scaf+=populateScaffoldTree(cscaf);      //recurse
    }
    //System.err.println("DEBUG: (populateScaffoldTree) leaving... ; n_scaf = "+n_scaf);
    return n_scaf;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Generates string representing the hierarchical scaffold sub tree
        rooted by this scaffold. Same format as Scaffold.subTreeAsString().
	e.g. "1:(2,3)" or "1:(2:(3,4,5),6:(4,7))"
  */
  public String getScaffoldString(ScaffoldDBRecord scafrec)
	throws SQLException
  {
    if (scafrec==null) return "";
    String str=""+scafrec.getID();
    if (scafrec.getChildIDs().size()>0)
    {
      //System.err.println("DEBUG: scafrec.getChildIDs().size()="+scafrec.getChildIDs().size());
      str+=":(";
      long i=0L;
      HashSet<Long> chids = scafrec.getChildIDs();
      for (Iterator<Long> itr=chids.iterator(); itr.hasNext(); )
      {
        long id=itr.next();
        if (i>0) str+=",";
        ScaffoldDBRecord cscafrec=this.getScaffoldByID(id);
        //if (cscafrec==null) System.err.println("DEBUG: id="+id+", cscafrec==null.");
        str+=this.getScaffoldString(cscafrec); //recurse
        ++i;
      }
      str+=")";
    }
    return str;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Is stereo flag on for scaffolds in this ScaffoldDB.
  */
  public Boolean isStereo()
  {
    return this.stereo;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Is keep_nitro_attachments flag on for scaffolds in this ScaffoldDB.
  */
  public Boolean isKeep_nitro_attachments()
  {
    return this.keep_nitro_attachments;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Dump all scaffolds in this ScaffoldDB to file; format is 
	legal smiles file and normal hier_scaffolds output:
	(SMILES ID SCAFTREE).  
	@return number of scaffolds written
  */
  public long dumpToFile(File fout,int verbose)
        throws IOException,SQLException
  {
    long nscaf=0;
    if (verbose>0)
    {
      System.err.println("Dumping scaffold DB to file.");
      System.err.println(this.info());
    }
    PrintWriter fout_writer=new PrintWriter(new BufferedWriter(new FileWriter(fout,false))); //overwrite

    String sql="SELECT id,scafsmi,scaftree FROM "+this.dbschema+"."+this.dbtableprefix+"scaffold";
    ResultSet rset=pg_utils.ExecuteSql(this.dbcon,sql);
    while (rset.next())
    {
      fout_writer.write(rset.getString("scafsmi")+" "+rset.getLong("id")+" "+rset.getString("scaftree")+"\n");
    }
    rset.getStatement().close();
    fout_writer.close();
    return nscaf;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	For testing only.
	<br/>
	CREATE ROLE www WITH LOGIN PASSWORD 'foobar';
	GRANT CONNECT ON DATABASE mydb TO www ;
	GRANT CREATE ON DATABASE mydb TO www ;
	(Then use schema "public".)
  */
  public static void main(String[] args)
  {
    if (args.length<8)
    {
      System.err.println("ERROR: syntax: ScaffoldDB DBHOST DBPORT DBNAME DBSCHEMA DBUSR DBPW TABLEPREFIX OFILE");
      System.exit(1);
    }

    Boolean stereo=false;
    Boolean keep_nitro_attachments=false;
    Boolean create_new=true;
    ScaffoldDB scafdb = null;
    try {
      scafdb = new ScaffoldDB(
	args[0], Integer.parseInt(args[1]), args[2], args[3], args[4], args[5], args[6],
	stereo, keep_nitro_attachments,
	create_new);

      boolean ok = scafdb.checkValid();
      System.err.println("DEBUG: checkValid(): "+ok);
    }
    catch (SQLException e) { System.err.println("SQLException: "+e.getMessage()); }

    if (scafdb!=null)
    {
      try {
        ScaffoldDBRecord scafrec = new ScaffoldDBRecord(1L,"N1CCC1");
        scafdb.addScaffold(scafrec);
        System.err.println("DEBUG: info(): "+scafdb.info());
      }
      catch (SQLException e) { System.err.println("SQLException: "+e.getMessage()); }
    }

    if (scafdb!=null)
    {
      try {
        File fout = new File(args[7]);
        int verbose=1;
        long n_scaf = scafdb.dumpToFile(fout,verbose);
        System.err.println("DEBUG: "+n_scaf+" scaffolds dumped to: "+args[7]);
        scafdb.destroyDB();
      }
      catch (SQLException e) { System.err.println("SQLException: "+e.getMessage()); }
      catch (IOException e) { System.err.println("IOException: "+e.getMessage()); }
    }
    return;
  }
}
