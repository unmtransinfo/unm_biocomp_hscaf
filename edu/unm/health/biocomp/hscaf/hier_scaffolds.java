package edu.unm.health.biocomp.hscaf;

import java.io.*;
import java.util.*;
import java.util.regex.*;

import chemaxon.formats.*;
import chemaxon.sss.search.*;
import chemaxon.struc.*;
import chemaxon.marvin.io.MolExportException;

import edu.unm.health.biocomp.hscaf.ScaffoldTree;
 
/**	HierS hierarchical scaffolds application&#46;  Run program with
	no args for command-line help&#46;
	<br />
	@author Jeremy J Yang
	@see edu.unm.health.biocomp.hscaf.ScaffoldTree
	@see edu.unm.health.biocomp.hscaf.Scaffold
	@see edu.unm.health.biocomp.hscaf.Linker
	@see edu.unm.health.biocomp.hscaf.Sidechain
	@see edu.unm.health.biocomp.hscaf.hier_scaffolds_utils
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
      +"    -i IFILE\n"
      +"    -o OFILE .................. one (possibly) disconnected mol of parts per input mol\n"
      +"  options:\n"
      +"    -out_scaf <outscafs> ...... unique scafs numbered sequentially\n"
      +"    -inc_mol      ............. include mols in output (-o)\n"
      +"    -inc_scaf ................. include scaffolds in output (-o)\n"
      +"    -inc_link   ............... include linkers in output (-o)\n"
      +"    -inc_chain ................ include sidechains in output (-o)\n"
      +"    -scaflist_title  .......... scaf/link/chain list written to title\n"
      +"    -scaflist_append2title  ... scaf/link/chain list appended to title\n"
      +"    -scaflist_sdtag <sdtag> ... scaf list written to SD dataitem\n"
      +"    -maxmol <max>   ........... max size/atoms of input mol [default=100]\n"
      +"    -show_js         .......... show junctions (as pseudoatoms) -- for debugging, visualizing\n"
      +"    -keep_nitro_attachments ... atoms single bonded to ring N remain in scaffold\n"
      +"    -stereo ................... stereo scaffolds (default non-stereo)\n"
      +"    -v ........................ verbose\n"
      +"    -vv ....................... very verbose\n"
      +"    -h ........................ this help\n");
    System.exit(1);
  }
  private static int verbose=0;
  private static String ifile=null;
  private static String ofile=null;
  private static String ofile_scaf=null;
  private static int maxmol=100;
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
  private static String smifmt="cxsmiles:u-L-l-e-d-D-p-R-f-w";

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
    if (ifile==null) help("Input file required.");
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

    MolExporter molWriter_scaf=null;
    if (ofile_scaf!=null)
    {
      String ofmt=MFileFormatUtil.getMostLikelyMolFormat(ofile_scaf);
      if (ofmt.equals("smiles")) ofmt="smiles:+n-a"; //Kekule for compatibility
      molWriter_scaf=new MolExporter(new FileOutputStream(ofile_scaf),ofmt);
    }

    if (verbose>0)
      System.err.println("JChem version: "+chemaxon.jchem.version.VersionInfo.JCHEM_VERSION);

    HashMap<String,Integer> scaf_usmis_global = new HashMap<String,Integer>();
    Molecule mol;
    int n_mol=0;
    int n_err=0;
    int n_mol_toobig=0;
    int n_mol_frag=0;
    int n_total_scaf=0;
    for (n_mol=0;true;)
    {
      boolean ok=true;
      try { mol=molReader.read(); }
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

      Molecule outmol = new Molecule();
      if (inc_mol)
        outmol.fuse(mol.cloneMolecule(),true);

      if (!ok)
      {
        if (ofile!=null) { if (!writeMol(molWriter,outmol)) ++n_err; }
        continue;
      }

      ScaffoldTree scaftree=null;
      try {
        scaftree = new ScaffoldTree(mol,keep_nitro_attachments,stereo);
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
      hier_scaffolds_utils.assignScaffoldIDs(scaf_usmis_global,scaftree);

      ArrayList<Integer> scaflist = new ArrayList<Integer>(); //for this mol
      int n_scaf=0;
      for (Scaffold scaf: scaftree.getScaffolds())
      {
        ++n_scaf;
        scaflist.add(scaf_usmis_global.get(scaf.getCansmi()));
        if (verbose>0)
        {
          System.err.println("\tscaf: "+n_scaf+". "+(show_js?scaf.getJsmi():scaf.getCansmi()));
          if (verbose>1)
          {
            System.err.print("\t\tID="+(scaf.getID()));
            System.err.print(" r["+(scaf.isRoot()?"x":" ")+"]");
            System.err.print(" l["+(scaf.isLeaf()?"x":" ")+"]");
            System.err.print(" pID="+(scaf.getParentScaffold()==null?" ":scaf.getParentScaffold().getID()));
            System.err.println(" nc="+scaf.getChildCount());
          }
        }
        if (inc_scaf)
          outmol.fuse(scaf.cloneMolecule(),true);
      }
      n_total_scaf+=n_scaf;
      if (verbose>0)
      {
        System.err.println("\t"+hier_scaffolds_utils.scafTreeAsString(scaftree.getRootScaffold()));
      }

      Collections.sort(scaflist);
      String scaflist_str="S:";
      for (int id: scaflist) { scaflist_str+=""+id+","; }
      scaflist_str=scaflist_str.replaceFirst(",$","");

      if (scaflist_title) outmol.setName(scaflist_str);
      else if (scaflist_append2title)
        outmol.setName(mol.getName()+" "+scaflist_str);
      outmol.setProperty(scaflist_sdtag,scaflist_str);

      int n_link=0;
      for (Linker link: scaftree.getLinkers())
      {
        ++n_link;
        if (verbose>0)
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
        if (verbose>0)
          System.err.println("\tsidechain: "+n_chain+". "+(show_js?chain.getJsmi():chain.getSmi()));
        if (inc_chain)
          outmol.fuse(chain.cloneMolecule(),true);
      }
      if (ofile!=null) { if (!writeMol(molWriter,outmol)) ++n_err; }
    }
    if (ofile_scaf!=null)
    {
      //Must output scafids as assigned and in order.
      HashMap<Integer,String> scaf_usmis_global_rev = new HashMap<Integer,String>();
      for (String scafsmi: scaf_usmis_global.keySet())
      {
        int scaf_id=scaf_usmis_global.get(scafsmi);
        scaf_usmis_global_rev.put(scaf_id,scafsmi);
      }

      for (int scaf_id=1;scaf_id<=scaf_usmis_global.size();++scaf_id)
      {
        Molecule scafmol=null;
        String scafsmi=scaf_usmis_global_rev.get(scaf_id);
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
        scafmol.setName(""+scaf_id);
        molWriter_scaf.write(scafmol);
      }
    }
    System.err.println("Total mols: "+n_mol);
    System.err.println("Errors: "+n_err);
    System.err.println("Total scaffolds found: "+n_total_scaf);
    System.err.println("Total unique scaffolds found: "+scaf_usmis_global.size());
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
