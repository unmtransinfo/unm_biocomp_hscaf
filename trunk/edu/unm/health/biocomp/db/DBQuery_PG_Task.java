package edu.unm.health.biocomp.db;

import java.io.*;
import java.util.*;
import java.sql.*;
import java.util.concurrent.*;

import edu.unm.health.biocomp.threads.*;
 
/**	Callable task for processing PostgreSql search&#46;
	@author Jeremy J Yang
*/
public class DBQuery_PG_Task
	implements Callable<Boolean>
{
  private String dbhost;
  private Integer dbport;
  private String dbid;
  private String dbusr;
  private String dbpw;
  private Connection dbcon;
  private String sql;
  private String errtxt;
  private Integer n_max;
  private ResultSet rset;
  private int n_total;
  private int n_done;
  public TaskStatus taskstatus;
  private java.util.Date t0;
  public DBQuery_PG_Task(String dbhost,Integer dbport,String dbid,String dbusr,String dbpw,String sql,Integer n_max)
  {
    this.dbhost=dbhost;
    this.dbport=dbport;
    this.dbid=dbid;
    this.dbusr=dbusr;
    this.dbpw=dbpw;
    this.sql=sql;
    this.n_max=n_max;
    this.taskstatus=new Status(this);
    this.n_total=0;
    this.n_done=0;
    this.t0 = new java.util.Date();
    this.rset=null;
    this.errtxt="";
  }
  /////////////////////////////////////////////////////////////////////////
  public synchronized ResultSet getRSet() { return this.rset; }
  public synchronized String getErrtxt() { return this.errtxt; }
  public synchronized Boolean call()
  {
    try {
      Connection dbcon=pg_utils.dbConnect(this.dbhost,this.dbport,this.dbid,this.dbusr,this.dbpw);
      this.rset=pg_utils.executeSql(dbcon,sql);
      dbcon.close();
    }
    catch (SQLException e)
    {
      this.errtxt=e.getMessage();
      return false;
    }
    return true;
  }
  class Status implements TaskStatus
  {
    private DBQuery_PG_Task task;
    public Status(DBQuery_PG_Task task) { this.task=task; }
    public String status()
    {
      long t=(new java.util.Date()).getTime()-t0.getTime();
      int m=(int)(t/60000L);
      int s=(int)((t/1000L)%60L);
      String statstr=("["+String.format("%02d:%02d",m,s)+"]");
      if (task.n_done>0)
        statstr+=(String.format(" %5d;",task.n_done));
      if (task.n_total>0)
        statstr+=(String.format(" %.0f%%",100.0f*task.n_done/task.n_total));
      return statstr;
    }
  }
}
