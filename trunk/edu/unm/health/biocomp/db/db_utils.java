package edu.unm.health.biocomp.db;

import java.io.*;
import java.util.*;
import java.sql.*; //Connection

//import org.postgresql.Driver; //Driver
//import com.mysql.jdbc.Driver; //Driver
//import oracle.jdbc.driver.OracleDriver; //OracleDriver

/**	
*/
public class db_utils
{
  public static void test_connection(String dbtype,String dbhost,Integer dbport,String dbname,String dbusr,String dbpw)
    throws Exception
  {
    DBCon dbcon = new DBCon(dbtype,dbhost,dbport,dbname,dbusr,dbpw);
    System.err.print("DEBUG: Attempting "+dbtype+" connection: "+dbhost+":"+dbport+":"+dbname+" ("+dbusr+")...");
    if (dbcon==null)
    {
      System.err.println(" FAILED.");
      return;
    }
    else
    {
      System.err.println(" OK.");
      System.err.println(dbcon.serverStatusTxt());
    }
  }

  public static void main(String[] args)
    throws IOException
    ,Exception
  {
    if (args.length < 6)
    {
      System.err.println("ERROR: args required: <DBTYPE> <DBHOST> <DBPORT> <DBNAME> <DBUSR> <DBPW>");
      System.exit(1);
    }
    String dbtype=args[0];
    String dbhost=args[1];
    Integer dbport=Integer.parseInt(args[2]);
    String dbname=args[3];
    String dbusr=args[4];
    String dbpw=args[5];

    test_connection(dbtype,dbhost,dbport,dbname,dbusr,dbpw);

  }
}
