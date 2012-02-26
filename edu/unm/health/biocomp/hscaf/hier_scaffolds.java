package edu.unm.health.biocomp.hscaf;

import java.io.*;
import java.util.*;
import java.util.regex.*;
import java.text.DateFormat;

import com.sleepycat.je.DatabaseException; // Some versions ".db." instead of ".je."?

import chemaxon.formats.*;
import chemaxon.sss.search.*;
import chemaxon.struc.*;
import chemaxon.marvin.io.MolExportException;

import edu.unm.health.biocomp.hscaf.ScaffoldTree;
 
/// to do: [ ] When using ScaffoldStore, use also ScaffoldSet, for selected (big, ugly) scaffolds.

/**	HierS hierarchical scaffolds application&#46;  Run program with
	no args for command-line help&#46;
	<br />
	@author Jeremy J Yang
	@see edu.unm.health.biocomp.hscaf.ScaffoldTree
	@see edu.unm.health.biocomp.hscaf.Scaffold
	@see edu.unm.health.biocomp.hscaf.Linker
	@see edu.unm.health.biocomp.hscaf.Sidechain
	@see edu.unm.health.biocomp.hscaf.hier_scaffolds_utils
	@see edu.unm.health.biocomp.hscaf.ScaffoldSet
	@see edu.unm.health.biocomp.hscaf.ScaffoldStore
	@see edu.unm.health.biocomp.hscaf.ScaffoldEntity
*/
public class hier_scaffolds
{
  private static void help(String msg)
  {
    System.err.println(msg+"\n"
      +"hier_scaffolds - HierS scaffolds analyzer\n"
      +"\n"
      +"From the University of New Mexico Division of Biocomputing\n"
      +"(http://biocomp.health.unm.edu/).\n"
      +"\n"
      +"usage: hier_scaffolds [options]\n"
      +"  required:\n"
      +"    -i IFILE .................. input file of molecules\n"
      +"    -o OFILE .................. one (possibly) disconnected mol of parts per input mol\n"
      +"  options:\n"
      +"    -out_scaf OUTSCAFS ........ unique scafs numbered sequentially (new + db)\n"
      +"    -inc_mol .................. include mols in output (-o)\n"
      +"    -inc_scaf ................. include scaffolds in output (-o)\n"
      +"    -inc_link ................. include linkers in output (-o)\n"
      +"    -inc_chain ................ include sidechains in output (-o)\n"
      +"    -scaflist_title  .......... scaf/link/chain list written to title\n"
      +"    -scaflist_append2title .... scaf/link/chain list appended to title\n"
      +"    -scaflist_sdtag SDTAG ..... scaf list written to SD dataitem\n"
      +"    -maxatoms MAX ............. max atom count of input mol [default=100]\n"
      +"    -maxrings MAX ............. max ring count of input mol [no default]\n"
      +"    -show_js .................. show junctions (as pseudoatoms) -- for debugging, visualizing\n"
      +"    -keep_nitro_attachments ... atoms single bonded to ring N remain in scaffold\n"
      +"    -stereo ................... stereo scaffolds (default non-stereo)\n"
      +"    -dbdir DBDIR .............. scratch dir for db files [default=/tmp/hscaf]\n"
      +"    -keepdb ................... keep scratch db files\n"
      +"    -destroyexistingdb ........ initially destroy existing db at DBDIR\n"
      +"    -nodb ..................... do not use BerkeleyDB for storage and performance\n"
      +"    -resume ................... resume job using kept scratch db files\n"
      +"    -nmax NMAX ................ quit after NMAX molecules\n"
      +"    -nskip NSKIP .............. skip NSKIP molecules\n"
      +"    -v ........................ verbose\n"
      +"    -vv ....................... very verbose\n"
      +"    -vvv or -debug ............ very very verbose (slows process)\n"
      +"    -h ........................ this help\n");
    System.exit(1);
  }
  private static int verbose=0;
  private static String ifile=null;
  private static String ofile=null;
  private static String ofile_scaf=null;
  private static int maxatoms=100;
  private static int maxrings=100;
  private static Boolean inc_mol=false;
  private static Boolean inc_link=false;
  private static Boolean inc_chain=false;
  private static Boolean inc_scaf=false;
  private static String scaflist_sdtag="HSCAF_SCAFLIST";
  private static Boolean scaflist_title=false;
  private static Boolean scaflist_append2title=false;
  private static Boolean show_js=false;
  private static Boolean keep_nitro_attachments=false;
  private static Boolean stereo=false;
  private static Boolean keepdb=false;
  private static Boolean resume=false;
  private static Boolean destroyexistingdb=false;
  private static Boolean nodb=false;
  private static int nmax=0;
  private static int nskip=0;
  private static String smifmt="cxsmiles:u-L-l-e-d-D-p-R-f-w";
  private static String dbdir="/tmp/hscaf";

  /////////////////////////////////////////////////////////////////////////////
  private static void parseCommand(String args[])
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
      else if (args[i].equals("-show_js")) show_js=true;
      else if (args[i].equals("-keep_nitro_attachments")) keep_nitro_attachments=true;
      else if (args[i].equals("-stereo")) stereo=true;
      else if (args[i].equals("-scaflist_sdtag")) scaflist_sdtag=args[++i];
      else if (args[i].equals("-maxatoms") || args[i].equals("-maxmol")) maxatoms=Integer.parseInt(args[++i]);
      else if (args[i].equals("-maxrings")) maxrings=Integer.parseInt(args[++i]);
      else if (args[i].equals("-dbdir")) dbdir=args[++i];
      else if (args[i].equals("-keepdb")) keepdb=true;
      else if (args[i].equals("-resume")) resume=true;
      else if (args[i].equals("-destroyexistingdb")) destroyexistingdb=true;
      else if (args[i].equals("-nodb")) nodb=true;
      else if (args[i].equals("-nmax")) nmax=Integer.parseInt(args[++i]);
      else if (args[i].equals("-nskip")) nskip=Integer.parseInt(args[++i]);
      else if (args[i].equals("-v")) verbose=1;
      else if (args[i].equals("-vv")) verbose=2;
      else if (args[i].equals("-vvv") || args[i].equals("-debug")) verbose=3;
      else if (args[i].equals("-h")) help("");
      else help("Unknown option: "+args[i]);
    }
  }
  /////////////////////////////////////////////////////////////////////////////
  public static void main(String[] args)
    throws IOException,DatabaseException
  {
    parseCommand(args);
    if (ifile==null) help("Input file required.");
    if (!(new File(ifile).exists())) help("Non-existent input file: "+ifile);
    MolImporter molReader = new MolImporter(ifile);

    MolExporter molWriter=null;
    if (ofile!=null)
    {
      if (!inc_mol && !inc_scaf && !inc_link && !inc_chain)
        help("-o requires one or more of: -inc_mol, inc_scaf, inc_link, inc_chain");
      String ofmt=MFileFormatUtil.getMostLikelyMolFormat(ofile);
      if (ofmt.equals("smiles")) ofmt="smiles:+n-a"; //Kekule for compatibility
      molWriter=new MolExporter(new FileOutputStream(ofile),ofmt);
    }
    if (resume && destroyexistingdb)
      help("ERROR: -resume and -destroyexisting incompatible.");
    if (resume && nodb)
      help("ERROR: -resume and -nodb incompatible.");

    MolExporter molWriter_scaf=null;
    if (ofile_scaf!=null)
    {
      String ofmt=MFileFormatUtil.getMostLikelyMolFormat(ofile_scaf);
      if (ofmt.equals("smiles")) ofmt="smiles:+n-a"; //Kekule for compatibility
      molWriter_scaf=new MolExporter(new FileOutputStream(ofile_scaf),ofmt);
    }

    if (verbose>1)
      System.err.println("JChem version: "+chemaxon.jchem.version.VersionInfo.JCHEM_VERSION);

    // ScaffoldStore used to store global unique scafs.
    File dbDir = new File(dbdir);
    assert(!(dbDir.exists() && dbDir.isDirectory()));
    ScaffoldStore scafstore=null;
    if (!nodb)
      scafstore = new ScaffoldStore(dbDir,stereo,keep_nitro_attachments);

    if (destroyexistingdb && scafstore!=null)
    {
      if (verbose>0)
        System.err.println("Destroying scratch DB at \""+dbdir+"\".");
      scafstore.destroyDB();
      scafstore = new ScaffoldStore(dbDir,stereo,keep_nitro_attachments);
    }

    // Check for existing contents of scafstore:
    long n_db_scaf_init=0L;
    if (scafstore!=null)
    {
      n_db_scaf_init=scafstore.count();
      if (n_db_scaf_init>0)
      {
        if (!resume)
          help("ERROR: Non-empty scratch DB at \""+dbdir+"\" (n_scafs="+n_db_scaf_init+") requires -resume or -destroyexistingdb.");
        else
          System.err.println("Non-empty scratch DB at \""+dbdir+"\" (n_scafs="+n_db_scaf_init+").");
        if (verbose>0)
          System.err.println("scafstore.info:\n"+scafstore.info(verbose));
      }
    }

    // If -nodb, ScaffoldSet used to store global unique scafs.
    ScaffoldSet scafset = null;
    if (nodb)
    {
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

      String molname=mol.getName();
      if (verbose>1)
      {
        System.err.println(""+n_mol+". "+molname);
        try { System.err.println("\t"+mol.exportToFormat(smifmt)); }
        catch (MolExportException e) { System.err.println(e.getMessage()); }
      }
      if (mol.getAtomCount()>maxatoms)
      {
        if (verbose>1)
          System.err.println("Warning: skipping mol; natoms="+mol.getAtomCount() +">"+ +maxatoms+" ["+n_mol+"] "+molname);
        ++n_mol_toobig;
        ok=false;
      }
      if (mol.getFragCount()>1)
      {
        System.err.println("Warning: multi-frag mol; analyzing largest frag only: ["+n_mol+"] "+molname);
        ++n_mol_frag;
        mol=hier_scaffolds_utils.largestPart(mol);
      }
      int ring_count=hier_scaffolds_utils.rawRingsystemCount(mol);
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
      Molecule outmol = new Molecule();
      if (inc_mol)
        outmol.fuse(mol.cloneMolecule(),true);
      if (!ok)
      {
        if (ofile!=null) { if (!writeMol(molWriter,outmol)) ++n_err; }
        continue;
      }
      java.util.Date t_this_0 = new java.util.Date();
      ScaffoldTree scaftree=null;
      try {
        scaftree = new ScaffoldTree(mol,stereo,keep_nitro_attachments,scafset,scafstore); //scafIDs assigned.
      }
      catch (Exception e) {
        ++n_err;
        System.err.println(e.getMessage());
        ok=false;
      }
      if (!ok)
      {
        if (ofile!=null) { if (!writeMol(molWriter,outmol)) ++n_err; }
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
          System.err.print("\t\tID="+(scaf.getID()));
          System.err.print(" cIDs="+scaf.getChildIDs());
          System.err.print(scaf.isRoot()?" (root)":"");
          System.err.print(scaf.isLeaf()?" (leaf)":"");
          System.err.println("");
        }
        if (inc_scaf)
          outmol.fuse(scaf.cloneMolecule(),true);
      }
      if (verbose>1)
      {
        System.err.println("\t"+scaftree.toString());
      }

      //if (scafset!=null) scafset.DEBUG_dump();

      Collections.sort(scaflist);
      String scaflist_str="S:";
      for (long id: scaflist) { scaflist_str+=""+id+","; }
      scaflist_str=scaflist_str.replaceFirst(",$","");

      if (scaflist_title) outmol.setName(scaflist_str);
      else if (scaflist_append2title)
        outmol.setName(mol.getName()+" "+scaflist_str);
      outmol.setProperty(scaflist_sdtag,scaflist_str);

      int n_link=0;
      for (Linker link: scaftree.getLinkers())
      {
        ++n_link;
        if (verbose>2)
        {
          System.err.println("\tlinker: "+n_link+". "+(show_js?link.getJsmi():link.getSmi()));
        }
        if (inc_link)
          outmol.fuse(link.cloneMolecule(),true);
      }
      int n_chain=0;
      for (Sidechain chain: scaftree.getSidechains())
      {
        ++n_chain;
        if (verbose>2)
          System.err.println("\tsidechain: "+n_chain+". "+(show_js?chain.getJsmi():chain.getSmi()));
        if (inc_chain)
          outmol.fuse(chain.cloneMolecule(),true);
      }
      if (verbose>1)
        System.err.println("\tt_this: "+util.timeDeltaStr(t_this_0,new java.util.Date()));
      if (ofile!=null) { if (!writeMol(molWriter,outmol)) ++n_err; }
      if (verbose>0 && n_mol%n_chunk==0)
      {
        System.err.print(" mols: "+n_mol+"; errors: "+n_err+": ");
        java.util.Date t_j = new java.util.Date();
        System.err.print("elapsed time: "+util.timeDeltaStr(t_0,t_j)+"; ");
        System.err.println("last "+n_chunk+" mols: "+util.timeDeltaStr(t_i,t_j));
        t_i=t_j;
      }
      if (nmax>0 && n_mol==(nmax+nskip)) break;
    }
    int n_scaf_unique=0;
    if (scafstore!=null)
    {
      for (long scaf_id=1L;true;++scaf_id)
      {
        String scafsmi=null;
        ScaffoldEntity scent=null;
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
          try { scafmol=MolImporter.importMol(scafsmi,"smiles:"); }
          catch (Exception e) {
          System.err.println(e.getMessage());
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
    else if (nodb)
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
    if (scafstore!=null)
    {
      if (verbose>0)
        System.err.println("scafstore.info:\n"+scafstore.info(verbose));
      n_db_scaf=scafstore.count();
      if (keepdb)
      {
        if (verbose>0)
          System.err.println("NOTE: scratch db kept at: "+dbdir+" N = "+n_db_scaf);
        scafstore.closeAll();
      }
      else
      {
        scafstore.destroyDB();
      }
    }
    molWriter_scaf.close();
    molWriter.close();
    System.err.println("Total elapsed time: "+util.timeDeltaStr(t_0,new java.util.Date()));
    if (verbose>0)
      System.err.println(DateFormat.getDateTimeInstance().format(new java.util.Date()));
    System.err.println("Total input mols: "+n_mol);
    if (nskip>0)
    {
      System.err.println("Skipped mols: "+nskip);
      System.err.println("Processed mols: "+(n_mol-nskip));
    }
    System.err.println("Scaffolds found: "+n_total_scaf);
    if (scafstore!=null)
      if (n_db_scaf>n_db_scaf_init)
        System.err.println("Unique scaffolds found this run: "+(n_db_scaf-n_db_scaf_init));
    System.err.println("Total unique scaffolds: "+n_scaf_unique);
    System.err.println("Oversized mols rejected: "+n_mol_toobig);
    System.err.println("Too-many-rings mols rejected: "+n_mol_toomanyrings);
    System.err.println("Multi-fragment mols (largest part analyzed only): "+n_mol_frag);
    System.err.println("Errors: "+n_err);
  }
  /////////////////////////////////////////////////////////////////////////////
  private static boolean writeMol(MolExporter molWriter,Molecule mol)
  {
    try {
      molWriter.write(mol);
    }
    catch (Exception e) {
      System.err.println(e.getMessage());
      return false;
    }
    return true;
  }
}
