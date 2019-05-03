package edu.unm.health.biocomp.hscaf;

import java.io.*;
import java.util.*;

import chemaxon.formats.*;      //MolImporter
import chemaxon.struc.*; //Molecule

/**     Test program.
*/
public class hscaf_scaf2smarts
{
  /////////////////////////////////////////////////////////////////////////////
  public static void main(String args[])
        throws IOException
  {
    if (args.length != 1)
    {
      System.err.println("syntax: hscaf_scaf2smarts IFILE");
      System.err.println("\tIFILE must contain scaffolds only (smiles or other formats).");
      System.exit(1);
    }
    String ifile=args[0];


    int n_all=0;
    int n_out=0;
    int n_err=0;

    MolImporter molReader = new MolImporter(ifile);

    Molecule mol;
    while ((mol=molReader.read())!=null)
    {
      ++n_all;
      String sma = hscaf_utils.ScaffoldMol2Smarts(mol);
      if (!sma.isEmpty())
      {
        System.out.println(sma);
        ++n_out;
      }
    }
    System.err.println("hscaf_scaf2smarts: n_all: "+n_all+"; n_out: "+n_out+"; n_err: "+n_err);
  }
}

