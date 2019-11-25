package edu.unm.health.biocomp.hscaf;

import java.io.*;
import java.util.*;
import java.util.regex.*;
import java.text.DateFormat;
import java.sql.*;

import com.sleepycat.je.DatabaseException; // Some versions ".db." instead of ".je."?

import chemaxon.formats.*;
import chemaxon.sss.search.*;
import chemaxon.struc.*;
import chemaxon.marvin.io.MolExportException;

import edu.unm.health.biocomp.util.db.*;
import edu.unm.health.biocomp.util.*;

/**	HierS hierarchical scaffolds application.  Run with
	no args for command-line help.
	<br>
	For each input molecule, scaffolds are perceived, including the hierarchical
	scaffold tree, and a set of unique scaffolds for the dataset is generated, with IDs assigned
	sequentially.
	Output molecule file annotated with scafIDs.
	Output scaffold file annotated with scafIDs and hierarchy.
	<br>
	@author Jeremy J Yang
*/
public class hier_scaffolds
{
  private static int verbose=0;
  private static String ifile=null;
  private static String ofile=null;
  private static Boolean keep_nitro_attachments=false;
  private static Boolean stereo=false;
  private static String ofile_scaf=null;
  private static Boolean append2ofile=false;
  private static Boolean inc_mol=false;
  private static Boolean inc_link=false;
  private static Boolean inc_chain=false;
  private static Boolean inc_scaf=false;
  private static Boolean scaflist_title=false;
  private static Boolean scaflist_append2title=false;
  private static String scaflist_sdtag="HSCAF_SCAFLIST";
  private static String scaftree_sdtag="HSCAF_SCAFTREE";
  private static Boolean show_js=false;

  private static Boolean bdb=false;
  private static Boolean bdb_keep=false;
  private static Boolean bdb_resume=false;
  private static Boolean bdb_predelete=false;
  private static String bdb_dump=null;

  private static Boolean rdb=false;
  private static Boolean rdb_keep=false;
  private static Boolean rdb_resume=false;
  private static Boolean rdb_predelete=false;
  private static String rdb_dump=null;
  private static String rdb_host=null;
  private static Integer rdb_port=5432; //PGSql default
  private static String rdb_name=null;
  private static String rdb_schema="public";
  private static String rdb_user=null;
  private static String rdb_pw=null;
  private static String rdb_tableprefix="";
  private static Integer rdb_reindex_per=5000;

  private static Boolean ddb=false;
  private static Boolean ddb_keep=false;
  private static Boolean ddb_resume=false;
  private static Boolean ddb_predelete=false;
  private static String ddb_dump=null;

  private static int maxatoms=100;
  private static int maxrings=10;
  private static int nmax=0;
  private static int nskip=0;
  //private static String db_preload=null;
  private static String SMIFMT_HSCAF="cxsmiles:u-L-l-e-d-D-p-R-f-w";
  private static String SMIFMT_OUT="smiles:+n-a-H"; //Kekule, H-implicit
  private static String SMIFMT_TSV="smiles:+n-a-HT*"; //Kekule, H-implicit, noHeader
  private static String bdb_dir="/tmp/hscaf";
  private static String ddb_dir="/tmp/hscaf";

  private static void Help(String msg)
  {
    System.err.println(msg+"\n"
      +"hier_scaffolds - HierS scaffolds analyzer\n"
      +"\n"
      +"From the UNM Translational Informatics Division\n"
      +"\n"
      +"usage: hier_scaffolds [options]\n"
      +"\n"
      +"required:\n"
      +"    -i IFILE .................. input file of molecules\n"
      +"  or\n"
      +"    -bdb_dump OFILE_SCAF ...... dump existing BerkeleyDB at BDBDIR\n"
      +"  or\n"
      +"    -rdb_dump OFILE_SCAF ...... dump existing RDB\n"
      +"\n"
      +"options:\n"
      +"  Algorithm:\n"
      +"    -keep_nitro_attachments ... atoms single bonded to ring N remain in scaffold\n"
      +"    -stereo ................... stereo scaffolds (default non-stereo)\n"
      +"  I/O:\n"
      +"    -o OFILE .................. one (possibly) disconnected mol of parts per input mol (.smi|.sdf)\n"
      +"    -out_scaf OFILE_SCAF ...... unique scafs numbered sequentially (new + db) (.smi|.sdf)\n"
      +"    -inc_mol .................. include mols in output (-o)\n"
      +"    -inc_scaf ................. include scaffolds in output (-o)\n"
      +"    -inc_link ................. include linkers in output (-o)\n"
      +"    -inc_chain ................ include sidechains in output (-o)\n"
      +"    -scaflist_title  .......... scaf/link/chain list written to title\n"
      +"    -scaflist_append2title .... scaf/link/chain list appended to title\n"
      +"    -scaflist_sdtag SDTAG ..... scaf list written to SD dataitem\n"
      +"    -append2ofile ............. output mols appended to -o OFILE (resume mode)\n"
      +"  RDB:\n"
      +"    -rdb ...................... use RDB for storage and performance\n"
      +"    -rdb_keep ................. keep RDB tables after job completion\n"
      +"    -rdb_predelete ............ initially delete RDB (if existing)\n"
      +"    -rdb_resume ............... resume job using existing RDB\n"
      +"    -rdb_host DBHOST .......... RDB host ["+rdb_host+"]\n"
      +"    -rdb_port DBPORT .......... RDB port ["+rdb_port+"]\n"
      +"    -rdb_name DBNAME .......... RDB db name\n"
      +"    -rdb_schema DBSCHEMA ...... RDB db schema ["+rdb_schema+"]\n"
      +"    -rdb_user DBUSER .......... RDB db user ["+rdb_user+"]\n"
      +"    -rdb_pw DBPW .............. RDB db password\n"
      +"    -rdb_tableprefix PREFIX ... RDB db table prefix\n"
      +"    -rdb_reindex_per N ........ RDB reindex once per N mols ["+rdb_reindex_per+"]\n"
      +"  BerkeleyDB:\n"
      +"    -bdb ...................... use BerkeleyDB for storage and performance\n"
      +"    -bdb_dir BDBDIR ........... scratch dir for BerkeleyDB files ["+bdb_dir+"]\n"
      +"    -bdb_keep ................. keep BerkeleyDB files after job completion\n"
      +"    -bdb_predelete ............ initially delete BerkeleyDB at BDBDIR (if existing)\n"
      +"    -bdb_resume ............... resume job using existing BerkeleyDB\n"
      +"  DerbyDB:\n"
      +"    -ddb ...................... use DerbyDB for storage and performance\n"
      +"    -ddb_dir DDBDIR ........... scratch dir for DerbyDB files ["+ddb_dir+"]\n"
      +"    -ddb_keep ................. keep DerbyDB files after job completion\n"
      +"    -ddb_predelete ............ initially delete DerbyDB at DDBDIR (if existing)\n"
      +"    -ddb_resume ............... resume job using existing DerbyDB\n"
      +"  Misc:\n"
      +"    -maxatoms MAX ............. max atom count of input mol ["+maxatoms+"]\n"
      +"    -maxrings MAX ............. max ring count of input mol ["+maxrings+"]\n"
      +"    -show_js .................. show junctions (as pseudoatoms) -- for debugging, visualizing\n"
      +"    -nmax NMAX ................ quit after NMAX molecules\n"
      +"    -nskip NSKIP .............. skip NSKIP molecules\n"
      //+"    -db_preload SCAFS ......... preload db from scaffold file for incremental analysis\n"
      +"    -v ........................ verbose\n"
      +"    -vv ....................... very verbose\n"
      +"    -vvv or -debug ............ very very verbose (slows process)\n"
      +"    -h ........................ this help\n");
    if (verbose>1)
      System.err.println("JChem version: "+com.chemaxon.version.VersionInfo.getVersion());
    System.exit(1);
  }
  /////////////////////////////////////////////////////////////////////////////
  private static void ParseCommand(String args[])
  {
    for (int i=0;i<args.length;++i)
    {
      if (args[i].equals("-i")) ifile=args[++i];
      else if (args[i].equals("-o")) ofile=args[++i];
      else if (args[i].equals("-out_scaf")) ofile_scaf=args[++i];
      else if (args[i].equals("-inc_mol")) inc_mol=true;
      else if (args[i].equals("-inc_scaf")) inc_scaf=true;
      else if (args[i].equals("-inc_link")) inc_link=true;
      else if (args[i].equals("-inc_chain")) inc_chain=true;
      else if (args[i].equals("-scaflist_title")) scaflist_title=true;
      else if (args[i].equals("-scaflist_append2title")) scaflist_append2title=true;
      else if (args[i].equals("-append2ofile")) append2ofile=true;
      else if (args[i].equals("-show_js")) show_js=true;
      else if (args[i].equals("-keep_nitro_attachments")) keep_nitro_attachments=true;
      else if (args[i].equals("-stereo")) stereo=true;
      else if (args[i].equals("-scaflist_sdtag")) scaflist_sdtag=args[++i];
      else if (args[i].equals("-maxatoms") || args[i].equals("-maxmol")) maxatoms=Integer.parseInt(args[++i]);
      else if (args[i].equals("-maxrings")) maxrings=Integer.parseInt(args[++i]);

      else if (args[i].equals("-rdb")) rdb=true;
      else if (args[i].equals("-rdb_keep")) rdb_keep=true;
      else if (args[i].equals("-rdb_resume")) rdb_resume=true;
      else if (args[i].equals("-rdb_predelete")) rdb_predelete=true;
      else if (args[i].equals("-rdb_dump")) rdb_dump=args[++i];
      else if (args[i].equals("-rdb_host")) rdb_host=args[++i];
      else if (args[i].equals("-rdb_port")) rdb_port=Integer.parseInt(args[++i]);
      else if (args[i].equals("-rdb_name")) rdb_name=args[++i];
      else if (args[i].equals("-rdb_schema")) rdb_schema=args[++i];
      else if (args[i].equals("-rdb_user")) rdb_user=args[++i];
      else if (args[i].equals("-rdb_pw")) rdb_pw=args[++i];
      else if (args[i].equals("-rdb_tableprefix")) rdb_tableprefix=args[++i];
      else if (args[i].equals("-rdb_reindex_per")) rdb_reindex_per=Integer.parseInt(args[++i]);

      else if (args[i].equals("-bdb")) bdb=true;
      else if (args[i].equals("-bdb_dir")) bdb_dir=args[++i];
      else if (args[i].equals("-bdb_keep")) bdb_keep=true;
      else if (args[i].equals("-bdb_resume")) bdb_resume=true;
      else if (args[i].equals("-bdb_predelete")) bdb_predelete=true;
      else if (args[i].equals("-bdb_dump")) bdb_dump=args[++i];

      else if (args[i].equals("-ddb")) ddb=true;
      else if (args[i].equals("-ddb_dir")) ddb_dir=args[++i];
      else if (args[i].equals("-ddb_keep")) ddb_keep=true;
      else if (args[i].equals("-ddb_resume")) ddb_resume=true;
      else if (args[i].equals("-ddb_predelete")) ddb_predelete=true;
      else if (args[i].equals("-ddb_dump")) ddb_dump=args[++i];

      else if (args[i].equals("-nmax")) nmax=Integer.parseInt(args[++i]);
      else if (args[i].equals("-nskip")) nskip=Integer.parseInt(args[++i]);
      //else if (args[i].equals("-db_preload")) db_preload=args[++i];
      else if (args[i].equals("-v")) verbose=1;
      else if (args[i].equals("-vv")) verbose=2;
      else if (args[i].equals("-vvv") || args[i].equals("-debug")) verbose=3;
      else if (args[i].equals("-h")) Help("");
      else Help("Unknown option: "+args[i]);
    }
  }
  /////////////////////////////////////////////////////////////////////////////
  public static void main(String[] args)
    throws IOException, DatabaseException
  {
    ParseCommand(args);

    /// Mode of operation:
    if (bdb_dump!=null)
    {
      if (ifile!=null || ofile!=null || ofile_scaf!=null)
        Help("ERROR: -i, -o, -o_scaf incompatible with -bdb_dump.");
      if (!bdb)
        Help("ERROR: -bdb_dump requires -bdb.");
      if (bdb_resume || bdb_predelete)
        Help("ERROR: -bdb_resume, -bdb_predelete incompatible with -bdb_dump.");
      if (verbose>0) System.err.println("Dumping scaffold store at: "+bdb_dir);

      File bdbDir = new File(bdb_dir);
      ScaffoldStore scafstore = new ScaffoldStore(bdbDir, stereo, keep_nitro_attachments);
      File fout = new File(bdb_dump);
      long nscaf = scafstore.dumpToFile(fout, verbose);
      System.err.println("Scaffold store at: "+bdb_dir+" dumped to: "+bdb_dump);
      System.exit(0);
    }
    else if (rdb_dump!=null)
    {
      if (ifile!=null || ofile!=null || ofile_scaf!=null)
        Help("ERROR: -i, -o, -o_scaf incompatible with -rdb_dump.");
      if (!rdb)
        Help("ERROR: -rdb_dump requires -rdb.");
      if (rdb_resume || rdb_predelete)
        Help("ERROR: -rdb_resume, -rdb_predelete incompatible with -rdb_dump.");

      if (verbose>0) System.err.println("Dumping scaffold DB: "+rdb_host);
      ScaffoldDB scafdb = null;
      try {
        scafdb = new ScaffoldDB(rdb_host, rdb_port, rdb_name, rdb_schema, rdb_user, rdb_pw, rdb_tableprefix,
                stereo, keep_nitro_attachments, false);
        File fout = new File(rdb_dump);
        long nscaf = scafdb.dumpToFile(fout, verbose);
        System.err.println("Scaffold DB: "+rdb_host+":"+rdb_port+"/"+rdb_name+" dumped: "+rdb_dump);
      }
      catch (SQLException e) {
        System.err.println("SQLException: "+e.getMessage());
        System.exit(1);
      }
      System.exit(0);
    }

    if (ifile==null) Help("Input file required.");

    if (verbose>1)
      System.err.println("JChem version: "+com.chemaxon.version.VersionInfo.getVersion());
      //System.err.println("JChem version: "+chemaxon.jchem.version.VersionInfo.JCHEM_VERSION); //pre-v6.3

    if (bdb_resume && bdb_predelete)
      Help("ERROR: -bdb_resume and -bdb_predelete incompatible.");
    if (bdb_resume && !bdb)
      Help("ERROR: -bdb_resume requires -bdb.");
    if ((bdb_resume||bdb_keep||bdb_predelete) && !bdb)
      Help("ERROR: -bdb_resume|-bdb_keep|-bdb_predelete require -bdb.");
    if (rdb_resume && rdb_predelete)
      Help("ERROR: -rdb_resume and -rdb_predelete incompatible.");
    if (rdb_resume && !rdb)
      Help("ERROR: -rdb_resume requires -rdb.");
    if ((rdb_host!=null ||rdb_name!=null) && !rdb)
      Help("ERROR: -rdb_host|-rdb_name require -rdb.");
    if ((rdb_user!=null ||rdb_pw!=null) && !rdb)
      Help("ERROR: -rdb_user|-rdb_pw require -rdb.");
    if ((rdb_resume||rdb_keep||rdb_predelete) && !rdb)
      Help("ERROR: -rdb_resume|-rdb_keep|-rdb_predelete require -rdb.");
    if (append2ofile && !bdb_resume && !rdb_resume)
      Help("ERROR: -append2ofile requires -bdb_resume or -rdb_resume.");
    //if (db_preload!=null && !rdb && !bdb)
    //  Help("ERROR: -db_preload requires -rdb or -bdb.");
    //if (db_preload!=null && (rdb_resume || bdb_resume))
    //  Help("ERROR: -db_preload incompatible w/ -rdb_resume, -bdb_resume.");

    if (!(new File(ifile).exists())) Help("Non-existent input file: "+ifile);
    MolImporter molReader = new MolImporter(ifile);

    MolExporter molWriter=null;
    if (ofile!=null)
    {
      if (!inc_mol && !inc_scaf && !inc_link && !inc_chain)
        Help("-o requires one or more of: -inc_mol, inc_scaf, inc_link, inc_chain");
      String ofmt=MFileFormatUtil.getMostLikelyMolFormat(ofile);
      if (ofmt.equals("smiles")) ofmt=SMIFMT_TSV;
      molWriter=new MolExporter(new FileOutputStream(ofile, append2ofile), ofmt);
    }
    else
    {
      molWriter=new MolExporter(System.out, SMIFMT_TSV);
    }
    MolExporter molWriter_scaf=null;
    if (ofile_scaf!=null)
    {
      String ofmt=MFileFormatUtil.getMostLikelyMolFormat(ofile_scaf);
      if (ofmt.equals("smiles")) ofmt=SMIFMT_OUT;
      molWriter_scaf=new MolExporter(new FileOutputStream(ofile_scaf), ofmt);
    }
    else
    {
      if (ofile==null)
        molWriter_scaf=molWriter;
      else
        molWriter_scaf=new MolExporter(System.out, SMIFMT_OUT);
    }

    ScaffoldStore scafstore=null; // ScaffoldStore stores global unique scafs in BerkeleyDB.
    File bdbDir=null;

    ScaffoldDB scafdb=null; // ScaffoldDB stores global unique scafs in RDB.
    Connection rdb_con=null;

    long n_db_scaf_init=0L; // Used by bdb or rdb.

    if (bdb)
    {
      bdbDir = new File(bdb_dir);
      if (!(bdbDir.exists() && bdbDir.isDirectory()))
      {
        try {
          bdbDir.mkdir();
          bdbDir.setWritable(true, true);
        }
        catch (Exception e)
        { System.err.println("ERROR: DBDIR creation failed ("+bdb_dir+"): "+e.getMessage()); }
      }
      scafstore = new ScaffoldStore(bdbDir, stereo, keep_nitro_attachments);

      if (bdb_predelete)
      {
        if (verbose>0)
          System.err.println("Deleting and reinitializing BerkeleyDB at \""+bdb_dir+"\".");
        scafstore.destroyDB();
        scafstore = new ScaffoldStore(bdbDir, stereo, keep_nitro_attachments);
      }

      // Check for existing contents of scafstore:
      n_db_scaf_init=scafstore.count();
      if (n_db_scaf_init>0)
      {
        if (!bdb_resume)
          Help("ERROR: Non-empty scratch DB at \""+bdb_dir+"\" (n_scafs="+n_db_scaf_init+") requires -bdb_resume or -bdb_predelete.");
        else
          System.err.println("Non-empty scratch DB at \""+bdb_dir+"\" (n_scafs="+n_db_scaf_init+").");
        if (verbose>0)
          System.err.println("scafstore.info:\n"+scafstore.info(verbose));
      }
    }
    else if (rdb)
    {
      try {
        rdb_con=pg_utils.DBConnect(rdb_host, rdb_port, rdb_name, rdb_user, rdb_pw);
        if (rdb_con==null)
          Help("RDB connection failed.");
        if (verbose>1)
        {
          System.err.println("PostgreSQL connection ok ("+rdb_host+":"+rdb_port+"/"+rdb_name+").");
          System.err.println(pg_utils.ServerStatusTxt(rdb_con));
        }
        if (rdb_resume)
        {
          scafdb = new ScaffoldDB(rdb_host, rdb_port, rdb_name, rdb_schema, rdb_user, rdb_pw, rdb_tableprefix,
            stereo, keep_nitro_attachments, false);
        }
        else
        {
          if (rdb_predelete)
          {
            if (hier_scaffolds_utils.CheckExistsDB(rdb_host, rdb_port, rdb_name, rdb_schema, rdb_user, rdb_pw, rdb_tableprefix))
            {
              if (verbose>0)
                System.err.println("Deleting and re-creating DB at "+rdb_host+".");
              scafdb = new ScaffoldDB(rdb_host, rdb_port, rdb_name, rdb_schema, rdb_user, rdb_pw, rdb_tableprefix,
                stereo, keep_nitro_attachments, false);
              scafdb.destroyDB();
            }
            else
            {
              if (verbose>0)
                System.err.println("No DB to delete at "+rdb_host+".");
            }
          }
          scafdb = new ScaffoldDB(rdb_host, rdb_port, rdb_name, rdb_schema, rdb_user, rdb_pw, rdb_tableprefix,
            stereo, keep_nitro_attachments, true);
          if (verbose>0)
            System.err.println("Scaffold DB created: "+rdb_host+":"+rdb_port+"/"+rdb_name+" (tableprefix: \""+rdb_tableprefix+"\")");
        }
        if (scafdb==null || !scafdb.checkValid()) Help("ERROR: scafdb initialization failed.");

        // Check for existing contents of scafdb:
        n_db_scaf_init=scafdb.count();
        if (n_db_scaf_init>0)
        {
          if (!rdb_resume)
            Help("ERROR: Non-empty RDB (n_scafs="+n_db_scaf_init+") requires -rdb_resume or -rdb_predelete.");
          else
            System.err.println("Non-empty RDB (n_scafs="+n_db_scaf_init+").");
          if (verbose>0)
            System.err.println("scafdb.info:\n"+scafdb.info());
        }
      }
      catch (SQLException e) { Help("SQLException: "+e.getMessage()); }
      if (verbose>0)
        System.err.println("database: "+rdb_host+":"+rdb_port+":"+rdb_name+":"+rdb_schema);
    }

    // If not -bdb, not -rdb, ScaffoldSet used to store global unique scafs.
    ScaffoldSet scafset = null;
    if (!bdb && !rdb)
    {
      if (verbose>0)
        System.err.println("NOTE: No persistent scaffold database used.  Slower for large datasets.");
      scafset = new ScaffoldSet("scaffold set from: "+ifile);
    }

    Molecule mol;
    int n_mol=0;
    int n_err=0;
    int n_mol_toobig=0;
    int n_mol_toomanyrings=0;
    int n_mol_frag=0;
    int n_total_scaf=0;
    java.util.Date t_0 = new java.util.Date();
    if (verbose>0)
      System.err.println(DateFormat.getDateTimeInstance().format(t_0));
    java.util.Date t_i = t_0;
    int n_chunk=100;
    for (n_mol=0;true;)
    {
      boolean ok=true;
      try { mol=molReader.read(); }
      catch (MolFormatException e)
      {
        System.err.println(e.getMessage());
        ++n_err;
        continue;
      }
      if (mol==null) break; //EOF
      ++n_mol;
      if (nskip>0 && n_mol<=nskip) continue;

      Molecule outmol = new Molecule();
      if (inc_mol)
      {
        outmol.setName(mol.getName());
        outmol.fuse(mol.cloneMolecule(), false);
      }
      String molname=mol.getName();
      if (verbose>1)
      {
        System.err.println(""+n_mol+". "+molname);
        try { System.err.println("\t"+MolExporter.exportToFormat(mol, SMIFMT_HSCAF)); }
        catch (MolExportException e) { System.err.println(e.getMessage()); }
        catch (IOException e) { System.err.println(e.getMessage()); }
      }
      if (mol.getAtomCount()>maxatoms)
      {
        if (verbose>1)
          System.err.println("Warning: skipping mol; natoms="+mol.getAtomCount() +">"+ +maxatoms+" ["+n_mol+"] "+molname);
        ++n_mol_toobig;
        ok=false;
      }
      if (mol.getFragCount(MoleculeGraph.FRAG_BASIC)>1)
      {
        if (verbose>1)
          System.err.println("Warning: multi-frag mol; analyzing largest frag only: ["+n_mol+"] "+molname);
        ++n_mol_frag;
        mol=hier_scaffolds_utils.LargestPart(mol);
      }
      int ring_count=hier_scaffolds_utils.RawRingsystemCount(mol);
      if (verbose>1)
      {
        System.err.println("\traw_ringsys_count: "+ring_count);
      }
      if (maxrings>0 && ring_count>maxrings)
      {
        if (verbose>1)
          System.err.println("Warning: skipping mol; nrings="+ring_count +">"+ +maxrings+" ["+n_mol+"] "+molname);
        ++n_mol_toomanyrings;
        ok=false;
      }
      if (!ok)
      {
        if (ofile!=null) { if (!WriteMol(molWriter, outmol)) ++n_err; }
        continue;
      }
      java.util.Date t_this_0 = new java.util.Date();
      ScaffoldTree scaftree=null;
      //System.err.println("DEBUG: ScaffoldTree() next...");
      try {
        if (rdb)
          scaftree = new ScaffoldTree(mol, stereo, keep_nitro_attachments, scafdb); //scafIDs assigned.
        else if (bdb)
          scaftree = new ScaffoldTree(mol, stereo, keep_nitro_attachments, scafstore); //scafIDs assigned.
        else 
          scaftree = new ScaffoldTree(mol, stereo, keep_nitro_attachments, scafset); //scafIDs assigned.
      }
      catch (Exception e) {
        ++n_err;
        //System.err.println(e.getMessage());
        System.err.println("DEBUG: ScaffoldTree() Exception: "+e.toString());
        ok=false;
      }
      if (!ok)
      {
        if (ofile!=null) { if (!WriteMol(molWriter, outmol)) ++n_err; }
        continue;
      }
      if (verbose>1)
      {
        System.err.print("\tn_scaf="+scaftree.getScaffoldCount());
        System.err.print("\tn_link="+scaftree.getLinkerCount());
        System.err.println("\tn_chain="+scaftree.getSidechainCount());
      }
      n_total_scaf+=scaftree.getScaffoldCount();

      ArrayList<Long> scaflist = new ArrayList<Long>(); //scafIDs for this mol
      int n_scaf=0;
      for (Scaffold scaf: scaftree.getScaffolds())
      {
        ++n_scaf;
        scaflist.add(scaf.getID());
        if (verbose>2)
        {
          System.err.println("\tscaf: "+n_scaf+". "+(show_js?scaf.getJsmi():scaf.getCansmi()));
          //System.err.println("\tDEBUG: "+hscaf_utils.ScaffoldSmarts(scaf));
          System.err.print("\t\tID="+(scaf.getID()));
          System.err.print(" cIDs="+scaf.getChildIDs());
          System.err.print(scaf.isRoot()?" (root)":"");
          System.err.print(scaf.isLeaf()?" (leaf)":"");
          System.err.println("");
        }
        if (inc_scaf)
          outmol.fuse(scaf.cloneMolecule(), false);
      }
      if (verbose>1)
      {
        System.err.println("\t"+scaftree.toString());
      }

      Collections.sort(scaflist);
      String scaflist_str="S:";
      for (long id: scaflist) { scaflist_str+=""+id+","; }
      scaflist_str=scaflist_str.replaceFirst(",$", "");

      if (scaflist_title) outmol.setName(scaflist_str);
      else if (scaflist_append2title)
        outmol.setName(mol.getName()+" "+scaflist_str);
      outmol.setProperty(scaflist_sdtag, scaflist_str);
      outmol.setProperty(scaftree_sdtag, scaftree.toString());

      int n_link=0;
      for (Linker link: scaftree.getLinkers())
      {
        ++n_link;
        if (verbose>2)
        {
          System.err.println("\tlinker: "+n_link+". "+(show_js?link.getJsmi():link.getSmi()));
        }
        if (inc_link)
          outmol.fuse(link.cloneMolecule(), false);
      }
      int n_chain=0;
      for (Sidechain chain: scaftree.getSidechains())
      {
        ++n_chain;
        if (verbose>2)
          System.err.println("\tsidechain: "+n_chain+". "+(show_js?chain.getJsmi():chain.getSmi()));
        if (inc_chain)
          outmol.fuse(chain.cloneMolecule(), false);
      }
      if (verbose>1)
        System.err.println("\tt_this: "+time_utils.TimeDeltaStr(t_this_0, new java.util.Date()));
      if (ofile!=null) { if (!WriteMol(molWriter, outmol)) ++n_err; }
      if (verbose>0 && n_mol%n_chunk==0)
      {
        System.err.print(" mols: "+n_mol+"; errors: "+n_err+": tod: "+time_utils.CurrentTime());
        java.util.Date t_j = new java.util.Date();
        System.err.print("; elapsed: "+time_utils.TimeDeltaStr(t_0, t_j)+"; dt: "+time_utils.TimeDeltaStr(t_i, t_j));
        System.err.println("; @: \""+molname+"\"");
        t_i=t_j;
      }
      if (nmax>0 && n_mol==(nmax+nskip)) break;
      if (rdb && ((n_mol-nskip)%rdb_reindex_per)==0)
      {
        if (verbose>0)
          System.err.println("REINDEX-ing database: "+rdb_host+":"+rdb_port+":"+rdb_name+":"+rdb_schema);
        try { scafdb.reindex(); }
        catch (SQLException e) { System.err.println("SQLException: "+e.getMessage()); }
      }
    }
    int n_scaf_unique=0;
    if (bdb)
    {
      for (long scaf_id=1L;true;++scaf_id)
      {
        String scafsmi=null;
        ScaffoldStoreEntity scent=null;
        try {
          scent=scafstore.scaffoldById.get(scaf_id);
          scafsmi=scent.getCanSmi();
        }
        catch (Exception e) { break; }
        if (scafsmi==null) break;
        ++n_scaf_unique;
        if (ofile_scaf!=null)
        {
          Molecule scafmol=null;
          try { scafmol=MolImporter.importMol(scafsmi, "smiles:"); }
          catch (Exception e) {
          //System.err.println(e.getMessage());
          System.err.println(e.toString());
            ++n_err;
          }
          if (scafmol==null)
          {
            System.err.println("error: importMol smiles: \""+scafsmi+"\"");
            ++n_err;
            scafmol = new Molecule();
          }
          scafmol.setName(""+scaf_id+" "+scafstore.getScaffoldString(scent));
          molWriter_scaf.write(scafmol);
        }
      }
    }
    else if (rdb)
    {
      String sql="SELECT id,scafsmi,scaftree FROM "+scafdb.getDBSchema()+"."+scafdb.getDBTableprefix()+"scaffold ORDER BY id ASC";
      try {
        ResultSet rset=pg_utils.ExecuteSql(rdb_con, sql);
        while (rset.next())
        {
          String scafsmi=rset.getString("scafsmi");
          if (scafsmi==null) continue;
          long scaf_id=rset.getLong("id");
          String scaftree=rset.getString("scaftree");
          ++n_scaf_unique;
          if (ofile_scaf!=null)
          {
            Molecule scafmol=null;
            try { scafmol=MolImporter.importMol(scafsmi, "smiles:"); }
            catch (Exception e) {
              //System.err.println(e.getMessage());
              System.err.println(e.toString());
              ++n_err;
            }
            if (scafmol==null)
            {
              System.err.println("error: importMol smiles: \""+scafsmi+"\"");
              ++n_err;
              scafmol = new Molecule();
            }
            scafmol.setName(""+scaf_id+" "+scaftree);
            molWriter_scaf.write(scafmol);
          }
        }
        rset.getStatement().close();
      }
      catch (SQLException e) {
        System.err.println("SQLException: "+e.getMessage());
        ++n_err;
      }
    }
    else
    {
      for (int scaf_id=1;scaf_id<=scafset.getCount();++scaf_id)
      {
        Scaffold scaf=scafset.getScaffoldByID(scaf_id);
        scaf.decompress();
        scaf.setName(""+scaf_id+" "+scaf.subTreeAsString());
        molWriter_scaf.write(scaf);
        ++n_scaf_unique;
      }
    }

    long n_db_scaf=0L;
    if (bdb)
    {
      if (verbose>0)
        System.err.println("scafstore.info:\n"+scafstore.info(verbose));
      n_db_scaf=scafstore.count();
      if (bdb_keep)
      {
        if (verbose>0)
          System.err.println("NOTE: scratch db kept at: "+bdb_dir+" N = "+n_db_scaf);
        scafstore.closeAll();
      }
      else
      {
        scafstore.destroyDB();
      }
    }
    else if (rdb)
    {
      try {
        if (verbose>0)
          System.err.println("scafdb.info:\n"+scafdb.info());
        n_db_scaf=scafdb.count();
        if (rdb_keep)
        {
          if (verbose>0)
            System.err.println("NOTE: RDB kept at: "+rdb_host+" N = "+n_db_scaf);
          scafdb.close();
        }
        else
        {
          scafdb.destroyDB();
        }
      }
      catch (SQLException e) {
        System.err.println("SQLException: "+e.getMessage());
        ++n_err;
      }
    }
    molWriter_scaf.close();
    molWriter.close();
    System.err.println("Total elapsed time: "+time_utils.TimeDeltaStr(t_0, new java.util.Date()));
    if (verbose>0)
      System.err.println(DateFormat.getDateTimeInstance().format(new java.util.Date()));
    System.err.println("Total input mols: "+n_mol);
    if (nskip>0)
    {
      System.err.println("Skipped mols: "+nskip);
      System.err.println("Processed mols: "+(n_mol-nskip));
    }
    System.err.println("Scaffolds found: "+n_total_scaf);
    if (bdb || rdb)
      if (n_db_scaf>n_db_scaf_init)
        System.err.println("Unique scaffolds found this run: "+(n_db_scaf-n_db_scaf_init));
    System.err.println("Total unique scaffolds: "+n_scaf_unique);
    System.err.println("Oversized mols rejected: "+n_mol_toobig);
    System.err.println("Too-many-rings mols rejected: "+n_mol_toomanyrings);
    System.err.println("Multi-fragment mols (largest part analyzed only): "+n_mol_frag);
    System.err.println("Errors: "+n_err);
    System.exit(0);
  }
  /////////////////////////////////////////////////////////////////////////////
  private static boolean WriteMol(MolExporter molWriter, Molecule mol)
  {
    try { molWriter.write(mol); }
    catch (Exception e) {
      //System.err.println(e.getMessage());
      System.err.println(e.toString());
      return false;
    }
    return true;
  }
}
