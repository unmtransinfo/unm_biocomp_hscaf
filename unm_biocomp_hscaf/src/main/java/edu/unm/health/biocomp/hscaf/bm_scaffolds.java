package edu.unm.health.biocomp.hscaf;

import java.io.*;
import java.util.*;
import java.util.regex.*;
import java.text.DateFormat;
import java.sql.*;

import chemaxon.formats.*;
import chemaxon.sss.search.*;
import chemaxon.struc.*;
import chemaxon.marvin.io.MolExportException;

import edu.unm.health.biocomp.util.db.*;
import edu.unm.health.biocomp.util.*;


/**	Bemis-Murko scaffolds application.  Run program with
	no args for command-line help.
	<br>
	@author Jeremy J Yang
	@see edu.unm.health.biocomp.hscaf.BMScaffold
*/
public class bm_scaffolds
{
  private static int verbose=0;
  private static String ifile=null;
  private static String ofile=null;
  private static Boolean keep_nitro_attachments=false;
  private static Boolean stereo=false;
  private static String ofile_scaf=null;
  private static Boolean show_js=false;

  private static int maxatoms=100;
  private static int maxrings=10;
  private static int nmax=0;
  private static int nskip=0;
  private static String smifmt="cxsmiles:u-L-l-e-d-D-p-R-f-w-a"; //Kekule

  private static void Help(String msg)
  {
    System.err.println(msg+"\n"
      +"bm_scaffolds - Bemis-Murko scaffolds analyzer\n"
      +"\n"
      +"From the University of New Mexico Translational Informatics Division\n"
      +"(http://medicine.unm.edu/informatics/).\n"
      +"\n"
      +"usage: bm_scaffolds [options]\n"
      +"\n"
      +"required:\n"
      +"    -i IFILE .................. input file of molecules\n"
      +"\n"
      +"options:\n"
      +"  Algorithm:\n"
      +"    -keep_nitro_attachments ... atoms single bonded to ring N remain in scaffold\n"
      +"    -stereo ................... stereo scaffolds (default non-stereo)\n"
      +"  I/O:\n"
      +"    -o OFILE .................. input mol and BM scaffold SMILES\n"
      +"  Misc:\n"
      +"    -maxatoms MAX ............. max atom count of input mol ["+maxatoms+"]\n"
      +"    -maxrings MAX ............. max ring count of input mol ["+maxrings+"]\n"
      +"    -show_js .................. show junctions (as pseudoatoms) -- for debugging, visualizing\n"
      +"    -nmax NMAX ................ quit after NMAX molecules\n"
      +"    -nskip NSKIP .............. skip NSKIP molecules\n"
      +"    -v ........................ verbose\n"
      +"    -vv ....................... very verbose\n"
      +"    -h ........................ this help\n");
    System.exit(1);
  }
  /////////////////////////////////////////////////////////////////////////////
  private static void ParseCommand(String args[])
  {
    for (int i=0;i<args.length;++i)
    {
      if (args[i].equals("-i")) ifile=args[++i];
      else if (args[i].equals("-o")) ofile=args[++i];
      else if (args[i].equals("-show_js")) show_js=true;
      else if (args[i].equals("-keep_nitro_attachments")) keep_nitro_attachments=true;
      else if (args[i].equals("-stereo")) stereo=true;
      else if (args[i].equals("-maxatoms") || args[i].equals("-maxmol")) maxatoms=Integer.parseInt(args[++i]);
      else if (args[i].equals("-maxrings")) maxrings=Integer.parseInt(args[++i]);
      else if (args[i].equals("-nmax")) nmax=Integer.parseInt(args[++i]);
      else if (args[i].equals("-nskip")) nskip=Integer.parseInt(args[++i]);
      else if (args[i].equals("-v")) verbose=1;
      else if (args[i].equals("-vv")) verbose=2;
      else if (args[i].equals("-vvv") || args[i].equals("-debug")) verbose=3;
      else if (args[i].equals("-h")) Help("");
      else Help("Unknown option: "+args[i]);
    }
  }
  /////////////////////////////////////////////////////////////////////////////
  public static void main(String[] args)
    throws IOException
  {
    ParseCommand(args);

    if (ifile==null) Help("Input file required.");

    if (!(new File(ifile).exists())) Help("Non-existent input file: "+ifile);
    MolImporter molReader = new MolImporter(ifile);

    MolExporter molWriter=null;
    if (ofile!=null)
    {
      String ofmt=MFileFormatUtil.getMostLikelyMolFormat(ofile);
      if (ofmt.equals("smiles")) ofmt="smiles:+n-a-H"; //Kekule for compatibility, H-implicitized
      molWriter=new MolExporter(new FileOutputStream(ofile,false),ofmt);
    }
    else
    {
      molWriter=new MolExporter(System.out,"smiles:+n-a-H");
    }
    //if (verbose>1)
      //System.err.println("JChem version: "+com.chemaxon.version.VersionInfo.getVersion()); //v6.3+
      //System.err.println("JChem version: "+chemaxon.jchem.version.VersionInfo.JCHEM_VERSION);

    Molecule mol;
    Molecule outmol = new Molecule();
    int n_mol=0;
    int n_err=0;
    int n_mol_toobig=0;
    int n_mol_frag=0;
    int n_total_scaf=0;
    java.util.Date t_0 = new java.util.Date();
    if (verbose>0)
      System.err.println(DateFormat.getDateTimeInstance().format(t_0));
    java.util.Date t_i = t_0;
    int n_chunk=1000;
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
        try { System.err.println("\t"+MolExporter.exportToFormat(mol,smifmt)); }
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
      outmol.clear();
      outmol.setName(mol.getName());
      outmol.fuse(mol.cloneMolecule(),true);
      if (!ok)
      {
        if (ofile!=null) { if (!WriteMol(molWriter,outmol)) ++n_err; }
        continue;
      }

      BMScaffold scaf = new BMScaffold(mol,stereo,keep_nitro_attachments);
      
      String scafsmi=scaf.getSmi(); //Kekule, for compatibility
      if (scafsmi!=null && !scafsmi.isEmpty())
      {
        outmol.setName(mol.getName()+" "+scafsmi);
        outmol.setProperty("SCAFSMI",scafsmi);
        ++n_total_scaf;
      }

      if (verbose>1)
      {
        System.err.println("\tscaf: "+(show_js?scaf.getJsmi():scaf.getSmi()));
        //System.err.println("\tDEBUG: "+hscaf_utils.ScaffoldSmarts(scaf));
        System.err.println("");
      }

      if (ofile!=null) { if (!WriteMol(molWriter,outmol)) ++n_err; }
      if (verbose>0 && n_mol%n_chunk==0)
      {
        System.err.print(" mols: "+n_mol+"; errors: "+n_err+": tod: "+time_utils.CurrentTime());
        java.util.Date t_j = new java.util.Date();
        System.err.print("; elapsed: "+time_utils.TimeDeltaStr(t_0,t_j)+"; dt: "+time_utils.TimeDeltaStr(t_i,t_j));
        System.err.println("; @: \""+molname+"\"");
        t_i=t_j;
      }
      if (nmax>0 && n_mol==(nmax+nskip)) break;
    }

    molWriter.close();
    System.err.println("Total elapsed time: "+time_utils.TimeDeltaStr(t_0,new java.util.Date()));
    if (verbose>0)
      System.err.println(DateFormat.getDateTimeInstance().format(new java.util.Date()));
    System.err.println("Total input mols: "+n_mol);
    if (nskip>0)
    {
      System.err.println("Skipped mols: "+nskip);
      System.err.println("Processed mols: "+(n_mol-nskip));
    }
    System.err.println("Scaffolds found: "+n_total_scaf);
    System.err.println("Oversized mols rejected: "+n_mol_toobig);
    System.err.println("Multi-fragment mols (largest part analyzed only): "+n_mol_frag);
    System.err.println("Errors: "+n_err);
    System.exit(0);
  }
  /////////////////////////////////////////////////////////////////////////////
  private static boolean WriteMol(MolExporter molWriter,Molecule mol)
  {
    try { molWriter.write(mol); }
    catch (Exception e) {
      System.err.println(e.getMessage());
      return false;
    }
    return true;
  }
}
