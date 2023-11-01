package edu.unm.health.biocomp.hscaf;

import java.util.*;
import java.io.IOException;
import java.sql.SQLException;

import junit.framework.*; //Test, TestCase, TestSuite

import chemaxon.formats.*; //MolFormatException
import chemaxon.marvin.io.MolExportException;
import chemaxon.sss.search.SearchException;
import chemaxon.struc.Molecule;

import edu.unm.health.biocomp.hscaf.*; //ScaffoldSet, ScaffoldTree, ScaffoldException

/**
 * Unit test for simple App.
 */
public class AppTest 
    extends TestCase
{
    /**
     * Create the test case
     *
     * @param testName name of the test case
     */
    public AppTest( String testName )
    {
        super( testName );
    }

    /**
     * @return the suite of tests being tested
     */
    public static Test suite()
    {
        return new TestSuite( AppTest.class );
    }

    /**
     * Rigourous Test :-)
     */
    public void testApp()
    {
        assertTrue( true );
    }

    /**
     * Better Test :-)
     */
    public void testScaf()
	throws Exception
	//throws IOException, SQLException, MolFormatException, MolExportException, SearchException, ScaffoldException
    {
      String smi = "COc1cc2c(ccnc2cc1)C(O)C4CC(CC3)C(C=C)CN34"; //quinine
      Molecule mol = MolImporter.importMol(smi, "smiles:");
      ScaffoldSet scafset = new ScaffoldSet("test");

      boolean stereo = false;
      boolean keep_nitro_attachments = false;
      ScaffoldTree scaftree = new ScaffoldTree(mol, stereo, keep_nitro_attachments, scafset);
      assertTrue( scaftree.getScaffoldCount() > 0 );
    }
}
