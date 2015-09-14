#!/usr/bin/env python

import sys
from itertools import permutations


def getDimNames(ndims):  
  dim_names = ['i', 'j', 'k', 'l', 'm', 'n']
  return dim_names[0:ndims]
  
def getDimPerms(dim_names):
  return permutations(dim_names)

def getDimPermNames(dim_names):
  PERM = permutations(dim_names)
  l = []
  for p in PERM:
    l.append("".join(p).upper())
  return l
  
def getEnumNames(ndims):
  dim_names = getDimNames(ndims)
  perm_names = getDimPermNames(dim_names)
  enum_names = map( lambda a: "LAYOUT_"+a, perm_names)
  return enum_names
  
def getEnumName(PERM):
  return "LAYOUT_" + "".join(PERM).upper();
 

def writeEnumDecl(ndims_list):
  for ndims in ndims_list:
    # Get names of each permutation    
    enum_names = getEnumNames(ndims)
    
    # Write an enum for each permutation
    for enum in enum_names:
      print "    struct %s {};" % enum    
    continue
    
    #print "    enum LAYOUT%dD {"%ndims
    #print "      " + ",\n      ".join(enum_names)
    #print "    };"
    #print ""

  print ""

def writeLayoutDecl(ndims_list):

  for ndims in ndims_list:
    dim_names = getDimNames(ndims)
  
    print "    template<typename L>"
    print "    struct Layout%dd {" % ndims
    
    # Define constructor
    args = map(lambda a: "int n"+a, dim_names)
    argstr = ", ".join(args)
    expstr = ""
    if ndims == 1:
      expstr = " explicit"
    print "        inline%s Layout%dd(%s);" % (expstr, ndims, argstr)
    
    # Define () Operator
    args = map(lambda a: "int "+a, dim_names)
    argstr = ", ".join(args)
    print "        inline int operator()(%s) const;" % argstr
    
    # Define toIndices function
    args = map(lambda a: "int &"+a, dim_names)
    argstr = ", ".join(args)
    print "        inline void toIndices(int linear, %s) const;" % argstr
        
    # Add local variables
    print ""
    args = map(lambda a: "int const size_"+a, dim_names)
    for arg in args:
      print "        %s;" % arg
      
    # Add stride variables
    print ""
    args = map(lambda a: "int const stride_"+a, dim_names)
    for arg in args:
      print "        %s;" % arg
      
    # Close out struct decl
    print "    };"
    print ""
  
 
  
def writeLayoutImpl(ndims_list):

  for ndims in ndims_list:
    dim_names = getDimNames(ndims)
    
    print ""
    print "/******************************************************************"
    print " *  Implementation for Layout%dD" % ndims
    print " ******************************************************************/"
    print ""
                
    # Loop over each permutation specialization
    perms = getDimPerms(dim_names)
    for perm in perms:
      # get enumeration name
      enum = getEnumName(perm)
    
      # Define constructor
      args = map(lambda a: "int n"+a, dim_names)
      argstr = ", ".join(args)    
      print "      template<>"
      print "      inline Layout%dd<%s>::Layout%dd(%s):" % (ndims, enum, ndims, argstr)    
      
      # initialize size of each dim
      args = map(lambda a: "size_%s(n%s)"%(a,a), dim_names)
      
      # initialize stride of each dim      
      for i in range(0,ndims):
        remain = perm[i+1:]
        if len(remain) > 0:
          remain = map(lambda a: "n"+a, remain)
          stride = "stride_%s(%s)" % ( perm[i],  "*".join(remain) )          
        else:
          stride = "stride_%s(1)" % ( perm[i] )
        args.append(stride)
          
      # output all initializers
      argstr = ", ".join(args)
      print "        %s" % argstr                    
      print "      {"
      print "      }"
      print ""

      # Define () Operator   
      args = map(lambda a: "int "+a, dim_names)
      argstr = ", ".join(args)   
      idxparts = []
      for i in range(0,ndims):
        remain = perm[i+1:]        
        if len(remain) > 0:
          idxparts.append("%s*stride_%s" % (perm[i], perm[i]))
        else:
          idxparts.append(perm[i])
      idx = " + ".join(idxparts)  

      print "      template<>"
      print "      inline int Layout%dd<%s>::operator()(%s) const {" % (ndims, enum, argstr)
      print "        return(" + idx + ");"
      print "      }"
      print ""
                 
      # Print out all of the linear->indices functions      
      args = map(lambda a: "int &"+a, dim_names)
      argstr = ", ".join(args)

      print "      template<>"
      print "      inline void Layout%dd<%s>::toIndices(int linear, %s) const {" % (ndims, enum, argstr)
      for i in range(0, ndims):
        idx = perm[i]
        prod = "*".join(map(lambda a: "size_%s"%a, perm[i+1 : ndims]))
        if prod != '':
          print "        %s = linear / (%s);" % (idx, prod)
          print "        linear -= %s*(%s);" % (idx, prod)
      print "        %s = linear;" % (perm[ndims-1]) 
      print "      }"
      print ""          


def writeViewDecl(ndims_list):

  for ndims in ndims_list:
    dim_names = getDimNames(ndims)
  
    print "    template<typename T, typename L>"
    print "    struct View%dd {" % ndims
    
    # Define constructor
    args = map(lambda a: "int n"+a, dim_names)
    argstr = ", ".join(args)
    print "        inline View%dd(T *data_ptr, %s);" % (ndims, argstr)
    
    # Define () Operator (const)
    args = map(lambda a: "int "+a, dim_names)
    argstr = ", ".join(args)
    print "        inline T const &operator()(%s) const;" % argstr
    
    # Define () Operator
    args = map(lambda a: "int "+a, dim_names)
    argstr = ", ".join(args)
    print "        inline T &operator()(%s);" % argstr
            
    # Add local variables
    print ""
    print "        Layout%dd<L> const layout;" % ndims
    print "        T *data;"
    print "    };"
    print ""


def writeViewImpl(ndims_list):

  for ndims in ndims_list:
    dim_names = getDimNames(ndims)
    
    print ""
    print "/******************************************************************"
    print " *  Implementation for View%dD" % ndims
    print " ******************************************************************/"
    print ""
                
    # Define constructor
    args = map(lambda a: "int n"+a, dim_names)
    argstr = ", ".join(args)    
    print "      template<typename T, typename L>"
    print "      inline View%dd<T,L>::View%dd(T *data_ptr, %s):" % (ndims, ndims, argstr)    
    args = map(lambda a: "n"+a, dim_names)
    argstr = ", ".join(args)
    print "        layout(%s)," % argstr
    print "        data(data_ptr)"
    print "      {"
    print "      }"
    print ""

    # Define () Operator (const)
    if True:
      args = map(lambda a: "int "+a, dim_names)
      argstr = ", ".join(args)      
      print "      template<typename T, typename L>"
      print "      inline T const &View%dd<T,L>::operator()(%s) const {" % (ndims, argstr)
      argstr = ", ".join(dim_names)
      print "        return(data[layout(%s)]);" % argstr
      print "      }"
      print ""

    # Define () Operator   
    args = map(lambda a: "int "+a, dim_names)
    argstr = ", ".join(args)      
    print "      template<typename T, typename L>"
    print "      inline T &View%dd<T,L>::operator()(%s){" % (ndims, argstr)
    argstr = ", ".join(dim_names)
    print "        return(data[layout(%s)]);" % argstr
    print "      }"
    print ""
               
  
 


# ACTUAL SCRIPT ENTRY:
print """//AUTOGENERATED BY genLayout.py
  
#ifndef __DOMAIN_LAYOUT_H__
#define __DOMAIN_LAYOUT_H__


"""

ndims_list = range(1,4+1)

# Dump all declarations (with documentation, etc)
writeEnumDecl(ndims_list)
writeLayoutDecl(ndims_list)
writeViewDecl(ndims_list)
#writeTLayoutDecl(ndims_list)

# Dump all implementations and specializations
writeLayoutImpl(ndims_list)
writeViewImpl(ndims_list)
#writeTLayoutImpl(ndims_list)

print """
  
#endif
"""

