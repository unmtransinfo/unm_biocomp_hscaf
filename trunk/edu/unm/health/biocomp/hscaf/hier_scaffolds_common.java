package edu.unm.health.biocomp.hscaf;

import java.io.*;
import java.util.*;
import java.util.regex.*;

import chemaxon.formats.*;
import chemaxon.sss.search.*;
import chemaxon.struc.*;
import chemaxon.marvin.io.MolExportException;

import edu.unm.health.biocomp.hscaf.ScaffoldTree;
 
/**	HierS hierarchical scaffolds application for max common scaffold search&#46; 
	Run program with no args for command-line help&#46;
	<br />
	@author Jeremy J Yang
	@see edu.unm.health.biocomp.hscaf.ScaffoldTree
	@see edu.unm.health.biocomp.hscaf.Scaffold
	@see edu.unm.health.biocomp.hscaf.Linker
	@see edu.unm.health.biocomp.hscaf.Sidechain
	@see edu.unm.health.biocomp.hscaf.ScaffoldSet
	@see edu.unm.health.biocomp.hscaf.hier_scaffolds_utils
*/
public class hier_scaffolds_common
{
  private static void help(String msg)
  {
    System.err.println(msg+"\n"
      +"hier_scaffolds_common - search for common HierS scaffolds\n"
      +"\n"
      +"From the University of New Mexico Division of Biocomputing\n"
      +"(http://biocomp.health.unm.edu/).\n"
      +"\n"
      +"usage: hier_scaffolds_common [options]\n"
      +"  required:\n"
      +"    -q IFILE .................. query molecule\n"
      +"    -db IFILE ................. db molecule[s]\n"
      +"    -o OFILE .................. common scaffold[s] ordered by size\n"
      +"  options:\n"
      +"    -out_scaf OUTSCAFS ........ unique scafs numbered sequentially\n"
      +"    -maxmol MAX ............... max size/atoms of input mol [default=100]\n"
      +"    -show_js .................. show junction points (as pseudoatoms) -- for debugging, visualizing\n"
      +"    -keep_nitro_attachments ... atoms single bonded to ring N remain in scaffold\n"
      +"    -stereo ................... stereo scaffolds (default is non-stereo)\n"
      +"    -v ........................ verbose\n"
      +"    -vv ....................... very verbose\n"
      +"    -h ........................ this help\n");
    System.exit(1);
  }
  private static int verbose=0;
  private static String ifileQ=null;
  private static String ifileN=null;
  private static String ofile=null;
  private static String ofile_scaf=null;
  private static int maxmol=100;
  private static Boolean show_js=false;
  private static Boolean keep_nitro_attachments=false;
  private static Boolean stereo=false;
  private static String smifmt="cxsmiles:u-L-l-e-d-D-p-R-f-w";

  /////////////////////////////////////////////////////////////////////////////
  private static void parseCommand(String args[])
  {
    for (int i=0;i<args.length;++i)
    {
      if (args[i].equals("-q")) ifileQ=args[++i];
      else if (args[i].equals("-db")) ifileN=args[++i];
      else if (args[i].equals("-o")) ofile=args[++i];
      else if (args[i].equals("-out_scaf")) ofile_scaf=args[++i];
      else if (args[i].equals("-show_js")) show_js=true;
      else if (args[i].equals("-keep_nitro_attachments")) keep_nitro_attachments=true;
      else if (args[i].equals("-stereo")) stereo=true;
      else if (args[i].equals("-maxmol")) maxmol=Integer.parseInt(args[++i]);
      else if (args[i].equals("-v")) verbose=1;
      else if (args[i].equals("-vv")) verbose=2;
      else if (args[i].equals("-d")) verbose=2;
      else if (args[i].equals("-h")) help("");
      else help("Unknown option: "+args[i]);
    }
  }
  /////////////////////////////////////////////////////////////////////////////
  public static void main(String[] args)
    throws IOException
  {
    parseCommand(args);
    if (ifileQ==null) help("Input query file required.");
    if (!(new File(ifileQ).exists())) help("Non-existent input file: "+ifileQ);
    MolImporter molReaderQ = new MolImporter(ifileQ);
    if (ifileN==null) help("Input db file required.");
    if (!(new File(ifileN).exists())) help("Non-existent input file: "+ifileN);
    MolImporter molReaderN = new MolImporter(ifileN);

    MolExporter molWriter=null;
    if (ofile!=null)
    {
      String ofmt=MFileFormatUtil.getMostLikelyMolFormat(ofile);
      if (ofmt.equals("smiles")) ofmt="smiles:+n-a"; //Kekule for compatibility
      molWriter=new MolExporter(new FileOutputStream(ofile),ofmt);
    }
    MolExporter molWriter_scaf=null;
    if (ofile_scaf!=null)
    {
      String ofmt=MFileFormatUtil.getMostLikelyMolFormat(ofile_scaf);
      if (ofmt.equals("smiles")) ofmt="smiles:+n-a"; //Kekule for compatibility
      molWriter_scaf=new MolExporter(new FileOutputStream(ofile_scaf),ofmt);
    }

    if (verbose>0)
      System.err.println("JChem version: "+chemaxon.jchem.version.VersionInfo.JCHEM_VERSION);

    //First analyze query mol.
    Molecule molQ=null;
    try { molQ=molReaderQ.read(); }
    catch (MolFormatException e) {
      System.err.println(e.getMessage());
      System.err.println("No query mol; quitting.");
      System.exit(1);
    }
    if (molQ.getAtomCount()>maxmol)
    {
      System.err.println("Oversized query mol; natoms>"+maxmol+"; quitting.");
      System.exit(1);
    }
    String molnameQ=molQ.getName();
    if (molQ.getFragCount()>1)
    {
      System.err.println("Warning: multi-frag query mol; analyzing largest frag only: "+molnameQ);
      molQ=hier_scaffolds_utils.largestPart(molQ);
    }
    ScaffoldTree scaftreeQ=null;
    try {
      scaftreeQ = new ScaffoldTree(molQ,stereo,keep_nitro_attachments);
    }
    catch (Exception e) {
      System.err.println(e.getMessage());
      System.err.println("HierS failed on query mol; quitting.");
    }

    // ScaffoldSet used to store global unique scafs.
    ScaffoldSet scafset = new ScaffoldSet("scaffold set from: "+ifileQ+" and "+ifileN);

    Molecule mol;
    int n_mol=0;
    int n_mcscaf=0;
    int n_err=0;
    int n_mol_toobig=0;
    int n_mol_frag=0;
    int n_total_scaf=0;
    for (n_mol=0;true;)
    {
      boolean ok=true;
      try { mol=molReaderN.read(); }
      catch (MolFormatException e) {
        System.err.println(e.getMessage());
        ++n_err;
        continue;
      }
      if (mol==null) break; //EOF
      ++n_mol;

      String molname=mol.getName();
      if (verbose>0)
      {
        System.err.println(""+n_mol+". "+molname);
        try { System.err.println("\t"+mol.exportToFormat(smifmt)); }
        catch (MolExportException e) { System.err.println(e.getMessage()); }
      }

      if (mol.getAtomCount()>maxmol)
      {
        System.err.println("Warning: skipping mol; natoms="+mol.getAtomCount() +">"+ +maxmol+" ["+n_mol+"] "+molname);
        ++n_mol_toobig;
        ok=false;
      }
      if (mol.getFragCount()>1)
      {
        System.err.println("Warning: multi-frag mol; analyzing largest frag only: ["+n_mol+"] "+molname);
        ++n_mol_frag;
        mol=hier_scaffolds_utils.largestPart(mol);
      }

      Molecule outmol = new Molecule();	//largest common scaf[s]

      if (!ok)
      {
        if (ofile!=null) { if (!writeMol(molWriter,outmol)) ++n_err; }	//empty placholder mol
        continue;
      }

      ScaffoldTree scaftree=null;
      try {
        scaftree = new ScaffoldTree(mol,stereo,keep_nitro_attachments,scafset);
      }
      catch (Exception e) {
        ++n_err;
        System.err.println(e.getMessage());
        ok=false;
      }
      if (verbose>0)
      {
        System.err.print("\tn_scaf="+scaftree.getScaffoldCount());
        System.err.print("\tn_link="+scaftree.getLinkerCount());
        System.err.println("\tn_chain="+scaftree.getSidechainCount());
      }
      n_total_scaf+=scaftree.getScaffoldCount();


      Scaffold mcscaf = hier_scaffolds_utils.maxCommonScaffold(scaftreeQ,scaftree);
      float sim = hier_scaffolds_utils.commonScaffoldTanimoto(scaftreeQ,scaftree,mcscaf);

      if (mcscaf!=null)
      {
        outmol.fuse(mcscaf.cloneMolecule(),true);
        if (show_js)
          System.err.println("max common scaf: "+mcscaf.getJsmi());
        else
          System.err.println("max common scaf: "+mcscaf.getSmi());
        ++n_mcscaf;
      }
      else
        outmol.fuse(new Molecule(),true);	//empty placeholder mol
      if (verbose>0)
        System.err.println(String.format("\tscaffold Tanimoto similarity = %.2f",sim));

      if (ofile!=null) { if (!writeMol(molWriter,outmol)) ++n_err; }
    }
    System.err.println("Total mols: "+n_mol);
    System.err.println("Max common scafs found: "+n_mcscaf);
    System.err.println("Errors: "+n_err);
    System.err.println("Total scaffolds found: "+n_total_scaf);
    System.err.println("Total unique scaffolds found: "+scafset.getCount());
    System.err.println("Oversized mols rejected: "+n_mol_toobig);
    System.err.println("Multi-fragment mols (largest part analyzed only): "+n_mol_frag);
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
