package edu.unm.health.biocomp.db;

import java.io.*;
import java.util.*;
import java.sql.*;

import oracle.jdbc.driver.OracleDriver;

/**	Static utility methods for Oracle databases.
	<br />
	Uses Oracle JDBC driver (oracle.jdbc.driver.OracleDriver).
	<br />
	@author Jeremy J Yang
	@see java.sql.DriverManager
	@see oracle.jdbc.driver.OracleDriver
*/
public class ora_utils
{
  /////////////////////////////////////////////////////////////////////////////
  /**	Return text with server status information.
  */
  public static String serverStatusTxt(Connection dbcon)
  {
    String txt="";
    try {
      ResultSet rset=executeSql(dbcon,"SELECT jchem_core_pkg.getenvironment() FROM dual");
      if (rset.next())
        txt+=(rset.getString(1));
      rset.getStatement().close();
    }
    catch (Exception e)
    {
      txt+=("error: "+e.getMessage());
    }
    return txt;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**   Executes SQL statement. Normally use this function.
  */
  public static ResultSet executeSql(Connection dbcon,String sql)
      throws SQLException
  {
    Statement stmt=dbcon.createStatement();
    ResultSet rset=stmt.executeQuery(sql);
    return rset;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**   Executes SQL statement. Use this function for non-queries,
        UPDATE, INSERT, DELETE, CREATE, REINDEX, etc.
        Note that autocommit is normally true so there is no
        need to call commit() directly.
  */
  public static boolean execute(Connection dbcon,String sql)
      throws SQLException
  {
    Statement stmt=dbcon.createStatement();
    boolean ok=stmt.execute(sql);
    stmt.close();
    return ok;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Return Oracle connection.
  */
  public static Connection dbConnect(String dbhost,Integer dbport,String dbid,String dbusr,String dbpw)
    throws SQLException
  {
    DriverManager.registerDriver(new oracle.jdbc.driver.OracleDriver());
    Connection dbcon=DriverManager.getConnection("jdbc:oracle:thin:@"+dbhost+":"+dbport+":"+dbid,dbusr,dbpw);
    return dbcon;
  }
  /////////////////////////////////////////////////////////////////////////////
}
