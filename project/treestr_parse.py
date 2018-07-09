#!/usr/bin/env python
#
### Input: 23:(45,56:(3),17:(5,99,4))
### Output:
#	{
#	  23:
#	    {
#	      45:{},
#	      56:{ 3:{} },
#	      17:{ 5:{},99:{},4:{} }
#	    }
#	}
#
import sys,os,re

#############################################################################
def ErrorExit(msg):
  print >>sys.stderr, msg
  sys.exit(1)

#############################################################################
def str2tree(str):
  print >>sys.stderr, ("DEBUG: before str: \"%s\""%str)

  str=re.sub(r'\s','',str) ## remove spaces

  # 1st look for all terminal sets of leaves, e.g. "(5,99,4)"
  # "17:(5,99,4)" -> "{17:{5,99,4}}"

  rob = re.compile(r'^(.*)(\d+):\(([\d,]+)\)(.*)$')
  while True:
    m=rob.search(str)
    if not m: break
    s1=m.group(1)
    s2=m.group(2)
    s3=m.group(3)
    s4=m.group(4)
    print >>sys.stderr, ("DEBUG: s2: \"%s\"  s3: \"%s\""%(s2,s3))
    s3=re.sub(r',',r':{},',s3)
    str = '%s%s:{%s:{}}%s'%(s1,s2,s3,s4)

  print >>sys.stderr, ("DEBUG: after str: \"%s\""%str)

  str = re.sub(r'(\d)([,\)])',r'\1:{}\2',str)

  str = re.sub(r'\(','{',str)
  str = re.sub(r'\)','}',str)

  str = '{'+str+'}'

  print >>sys.stderr, ("DEBUG: final str: \"%s\""%str)

  tree = eval(str)
  return tree

#############################################################################
def tree2str(tree,level=0):
  if not tree:
    return ''
  else:
    str=''
    for key,val in tree.items():
      indent=('\t'*level)
      str+=("%s%s:\n%s{%s\n%s}"%(indent,key,indent,tree2str(val,level+1),indent))
    return str

#############################################################################
if __name__=='__main__':

  treestr = "23:(45,56:(3),17:(5,99,4))"
  tree = str2tree(treestr)
  print "treestr: %s"%treestr
  print "   tree: %s"%tree

  print tree2str(tree)
