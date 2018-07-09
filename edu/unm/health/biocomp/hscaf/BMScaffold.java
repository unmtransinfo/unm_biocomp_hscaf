package edu.unm.health.biocomp.hscaf;

import java.io.*;
import java.util.*;

import chemaxon.formats.*;
import chemaxon.struc.*;
import chemaxon.sss.search.*;
import chemaxon.marvin.io.MolExportException;
 
/**	Represents the Bemis-Murko (BM) framework, for use in analyses not
	requiring a hierarchical scaffold (HierS) model.
	All molecules have zero or one BM framework.
	<br>
	@see edu.unm.health.biocomp.hscaf.bm_scaffolds
	@author Jeremy J Yang
*/
public class BMScaffold extends Molecule implements Comparable<Object>
{
  /** JChem format used for getCansmi() , for canonical scaffold smiles.
      E.g. cansmi=MolExporter.exportToFormat(mol,BMScaffold.CANSMIFMT)
  */
  public static final String CANSMIFMT="cxsmiles:u-L-l-e-d-D-p-R-f-w-H+a_gen+0"; //fmt for getCansmi()
  public static final String CANSMIFMT_STEREO="cxsmiles:u-L-l-e-d-D-p-R-f-w-H+a_gen"; //fmt for getCansmi()
  private static String smifmt="cxsmiles:u-L-l-e-d-D-p-R-f-w-a"; //fmt for getSmi()
  private static String cxsmifmt="cxsmiles:+L+l+e+d+D+p+R+f+w"; //fmt for getJsmi()
  private String cansmi; //cansmiles for equality comparison
  private String smi; //Kekule smiles (good for export)
  private Long id; //unique ID for dataset scope
  private boolean keep_nitro_attachments;
  private boolean stereo;
  /////////////////////////////////////////////////////////////////////////////
  private BMScaffold() {} //disallow default constructor
  /////////////////////////////////////////////////////////////////////////////
  /**	When initialized from a molecule, sidechains are removed
	and explicit hydrogens are added at each junction bond.
	If molecule has no rings, return empty scaffold.
  */
  public BMScaffold(Molecule mol,boolean stereo,boolean keep_nitro_attachments)
    throws MolFormatException
  {
    if (!stereo)
    {
      this.smifmt+="+0"; //non-stereo
      this.cxsmifmt+="+0"; //non-stereo
    }
    this.id=0L;
    this.cansmi=null; // Lazily evaluated
    this.smi=null; // Lazily evaluated
    this.keep_nitro_attachments=keep_nitro_attachments;
    this.stereo=stereo;
    int [][] sssr = mol.getSSSR();
    boolean ok=(sssr.length>0);

    if (ok)
    {
      mol.clonecopy(this);
      ok=this.dearomatize(); //avoid auto-H-adjustments by JChem
    }
    if (ok)
    {
      hscaf_utils.rmSideChains(this,keep_nitro_attachments);
      this.aromatize(MoleculeGraph.AROM_GENERAL);
    }
    ok=this.isLegal(); //benzene
    if (!ok) this.clear();
  }
  /////////////////////////////////////////////////////////////////////////////
  /**   
  */
  public void clear()
  {
    super.clear();
    this.cansmi=null;
    this.smi=null;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**   Initialized molecule from smiles and call usual constructor.
  */
  public BMScaffold(String smiles,Boolean stereo,Boolean keep_nitro_attachments)
    throws MolFormatException
  {
    this(new Molecule(),stereo,keep_nitro_attachments);
    try { MolImporter.importMol(smiles.getBytes(),"smiles:",null,this); }
    catch (MolFormatException e) { } // should not happen!
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Copy constructor.
  */
  public BMScaffold(BMScaffold scaf)
  {
    scaf.clonecopy(this);
    this.cansmi=scaf.cansmi;
    this.smi=scaf.smi;
    this.id=scaf.id;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	For sorting by size.
  */
  public int compareTo(Object o) throws ClassCastException
  {
    return ((BMScaffold)o).getAtomCount()-this.getAtomCount();
  }
  /////////////////////////////////////////////////////////////////////////////
  public boolean equals(BMScaffold scaf2)
  {
    return (this.getCansmi().equals(scaf2.getCansmi()));
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Reject (1) benzene, and (2) empty scaffolds .
	Note that benzene is identified by canonical smiles, unfortunately 
	hard-coded.
  */
  public boolean isLegal()
  {
    if (this.isEmpty()) this.decompress(); // isEmpty() from Molecule
    if (this.getAtomCount()==0) return false;
    if (this.getCansmi().equals("c1ccccc1")) return false; //benzene disallowed
    if (this.getCansmi().equals("C1=CC=CC=C1")) return false; //benzene disallowed
    return true;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Get unique ID for dataset scope.
  */
  public long getID() { return this.id; }
  /**	Set unique ID for dataset scope.
  */
  public void setID(long id2) { this.id=id2; }
  /////////////////////////////////////////////////////////////////////////////
  /**	Get canonical SMILES for this scaffold (as used for equality comparison).
	Canonical SMILES are generated using JChem format Scaffold.CANSMIFMT.
	Lazy evaluation.
  */
  public String getCansmi()
  {
    if (this.cansmi==null)
    {
      try { this.cansmi=MolExporter.exportToFormat(this,(this.stereo?this.CANSMIFMT_STEREO:this.CANSMIFMT)); }
      catch (MolExportException e) { System.err.println(e.getMessage()); this.cansmi=""; }
      catch (IOException e) { System.err.println(e.getMessage()); this.cansmi=""; }
    }
    return this.cansmi;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Get Kekule SMILES for this scaffold.  Recommended for export
	as resulting SMILES will be more universally compatible with
	other software. Lazy evaluation.
  */
  public String getSmi()
  {
    if (this.smi==null)
    {
      this.decompress();
      try { this.smi=MolExporter.exportToFormat(this,this.smifmt); }
      catch (MolExportException e) { System.err.println(e.getMessage()); this.smi=""; }
      catch (IOException e) { System.err.println(e.getMessage()); this.smi=""; }
    }
    return this.smi;
  }
  ///////////////////////////////////////////////////////////////////////////
  /**	Get SMILES for this scaffold with junction atoms as pseudo-atoms.
	Note this is in ChemAxon cxsmiles format e.g.
	"*C1CC1* |$p_J;;;;p_J$|".
  */
  public String getJsmi()
  {
    this.decompress();
    Molecule scafmol=this.cloneMolecule();
    hscaf_utils.replaceJHydrogensWithJPseudoatoms(scafmol);
    String jsmi=null;
    try { jsmi=MolExporter.exportToFormat(scafmol,cxsmifmt); }
    catch (MolExportException e) { jsmi=""; }
    catch (IOException e) { jsmi=""; }
    return jsmi;
  }
  ///////////////////////////////////////////////////////////////////////////
  /**	For storage; delete molecule object and retain canonical SMILES.
  	Canonical SMILES used since it is an identifier.
  */
  public void compress()
  {
    this.smi=null;
    if (this.cansmi==null) this.cansmi=this.getCansmi();
    if (!this.cansmi.isEmpty() && !this.isEmpty()) this.clear();
  }
  ///////////////////////////////////////////////////////////////////////////
  /**	Reconstruct molecule object from stored canonical SMILES.
  */
  public void decompress()
  {
    if (!this.isEmpty()) return; //already decompressed
    if (this.cansmi==null) return; //empty scaffold
    try { MolImporter.importMol(this.cansmi.getBytes(),"smiles:",null,this); }
    catch (MolFormatException e) { } // should not happen!
  }
  ///////////////////////////////////////////////////////////////////////////
  public Boolean isStereo()
  {
    return this.stereo;
  }
  ///////////////////////////////////////////////////////////////////////////
  public Boolean isKeep_nitro_attachments()
  {
    return this.keep_nitro_attachments;
  }
  ///////////////////////////////////////////////////////////////////////////
}
