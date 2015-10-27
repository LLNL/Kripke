#!/usr/bin/env python

import sys
from itertools import permutations
from lperm import *

def writeForallBase(ndims):

  dim_names = getDimNames(ndims)

  print ""
  print "/******************************************************************"
  print " *  OpenMP Parallel Region forall%d()" % ndims
  print " ******************************************************************/"
  print ""
  print "#ifdef _OPENMP"
  print ""
  
  args = map(lambda a: "typename T"+a.upper(), dim_names)
  argstr = ", ".join(args)    
  print "    template<typename POLICY, %s, typename BODY>" % (argstr)
  
  args = map(lambda a: "T%s const &is_%s"%(a.upper(), a), dim_names)
  argstr = ", ".join(args)
  print "    RAJA_INLINE void forall%d(Forall%d_OMP_Parallel_Tag, %s, BODY const &body){" % (ndims, ndims, argstr)
    
  print "      typedef typename POLICY::NextPolicy NextPolicy;"
  print "      typedef typename POLICY::NextPolicy::PolicyTag NextPolicyTag;"
  print "      // create OpenMP Parallel Region"
  print "#pragma omp parallel"
  print "      {"
  print "        // execute the next policy"
  
  args = map(lambda a: "T"+a.upper(), dim_names)
  argstr = ", ".join(args)    
  
  args = map(lambda a: "is_"+a, dim_names)
  argstr2 = ", ".join(args)  
  print "        forall%d<NextPolicy, %s, BODY>(NextPolicyTag(), %s, body);" % (ndims, argstr, argstr2)
  print "      }"
  print "    }"
  print ""
  print "#endif"
  print ""



  print ""
  print "/******************************************************************"
  print " *  Tiling Policy for forall%d()" % ndims
  print " ******************************************************************/"
  print ""
  
  args = map(lambda a: "typename T"+a.upper(), dim_names)
  argstr = ", ".join(args)    
  print "    template<typename POLICY, %s, typename BODY>" % (argstr)
  
  args = map(lambda a: "T%s const &is_%s"%(a.upper(), a), dim_names)
  argstr = ", ".join(args)
  print "    RAJA_INLINE void forall%d(Forall%d_Tile_Tag, %s, BODY const &body){" % (ndims, ndims, argstr)
    
  print "      typedef typename POLICY::NextPolicy NextPolicy;"
  print "      typedef typename POLICY::NextPolicy::PolicyTag NextPolicyTag;"
  for d in dim_names:
    print "      typedef typename POLICY::Tile%s Tile%s;" % (d.upper(), d.upper())
  print ""
  print "      // execute the next policy"
  
  indent = ""
  close_paren = []
  for d in dim_names:
    print "%s      forall_tile(Tile%s(), is_%s, [=](auto is_%s%s){" % (indent, d.upper(), d, d, d)    
    close_paren.append(indent + "      });")
    indent += "  "
    
  # call body with tiled index sets
  args = map(lambda a: "is_%s%s"%(a,a), dim_names)
  argstr2 = ", ".join(args)  
  print "%s      forall%d<NextPolicy>(NextPolicyTag(), %s, body);" % (indent, ndims, argstr2)

  # close forall parenthesis
  close_paren.reverse()
  for c in close_paren:
    print c

  print "    }"
  print ""
  print ""



    
  print ""
  print "/******************************************************************"
  print " *  Execute policy, forall%d()" % ndims
  print " ******************************************************************/"
  print ""
  
  args = map(lambda a: "typename T"+a.upper(), dim_names)
  argstr = ", ".join(args)    
  print "    template<typename POLICY, %s, typename BODY>" % (argstr)
  
  args = map(lambda a: "T%s const &is_%s"%(a.upper(), a), dim_names)
  argstr = ", ".join(args)
  print "    RAJA_INLINE void forall%d(Forall%d_Execute_Tag, %s, BODY const &body){" % (ndims, ndims, argstr)
  
  args = map(lambda a: "T"+a.upper(), dim_names)
  argstr = ", ".join(args)    
  
  args = map(lambda a: "is_"+a, dim_names)
  argstr2 = ", ".join(args)
  print "      typedef typename POLICY::LoopOrder L;"
  print "      forall%d_permute<POLICY, %s, BODY>(L(), %s, body);" % (ndims, argstr, argstr2)  
  print "    }"  
  print ""

  print ""
  print "/******************************************************************"
  print " *  User interface, forall%d()" % ndims
  print " ******************************************************************/"
  print ""

  args = map(lambda a: "typename Idx%s=int"%a.upper(), dim_names)
  idxstr = ", ".join(args)    
  args = map(lambda a: "typename T"+a.upper(), dim_names)
  argstr = ", ".join(args)    
  print "    template<typename POLICY, %s, %s, typename BODY>" % (idxstr, argstr)
  
  args = map(lambda a: "T%s const &is_%s"%(a.upper(), a), dim_names)
  argstr = ", ".join(args)
  print "    RAJA_INLINE void forall%d(%s, BODY const &body){" % (ndims, argstr)
  
  args = map(lambda a: "T"+a.upper(), dim_names)
  argstr = ", ".join(args)    
  
  args = map(lambda a: "is_"+a, dim_names)
  argstr2 = ", ".join(args)
  print "      typedef typename POLICY::PolicyTag PolicyTag;"
  print "      forall%d<POLICY, %s>(PolicyTag(), %s, " % (ndims, argstr, argstr2)
  
  args = map(lambda a: "int %s"%a, dim_names)
  argstr = ", ".join(args)
  print "        [=](%s){" % argstr
  
  args = map(lambda a: "Idx%s(%s)"%(a.upper(), a), dim_names)
  argstr = ", ".join(args)
  print "          body(%s);" % argstr
  print "        }"
  print "      );"
  print "    }"  
  print ""


def writeForallPolicy(ndims):

  print ""
  print "/******************************************************************"
  print " *  Policy base class, forall%d()" % ndims
  print " ******************************************************************/"
  print ""

  dim_names = getDimNames(ndims)

  print "    // Interchange-loops and Execute (Base-case for all policies)" 
  print "    struct Forall%d_Execute_Tag {};" % ndims
  args = map(lambda a: "typename POL_"+(a.upper()), dim_names)
  argstr = ", ".join(args)
  print "    template<typename LOOP_ORDER, %s>" % argstr
  print "    struct Forall%d_Execute {" % (ndims)
  print "      typedef Forall%d_Execute_Tag PolicyTag;" % ndims
    
  args = map(lambda a: "end_"+a, dim_names)
  argstr = ", ".join(args)
  print "      typedef LOOP_ORDER LoopOrder;"
  for dim in dim_names:
    print "      typedef POL_%s Policy%s;" % (dim.upper(), dim.upper())
  print "    };"
  print ""
  
  print "    // Begin OpenMP Parallel Block"
  print "    struct Forall%d_OMP_Parallel_Tag {};" % ndims
  print "    template<typename NEXT>"
  print "    struct Forall%d_OMP_Parallel {" % ndims
  print "      typedef Forall%d_OMP_Parallel_Tag PolicyTag;" % ndims
  print "      typedef NEXT NextPolicy;"
  print "    };"
  print ""
  
  print "    // Tiling Policy"
  print "    struct Forall%d_Tile_Tag {};" % ndims
  args = map(lambda a: "typename TILE_"+(a.upper()), dim_names)
  argstr = ", ".join(args)
  print "    template<%s, typename NEXT>" % argstr
  print "    struct Forall%d_Tile {" % ndims
  print "      typedef Forall%d_Tile_Tag PolicyTag;" % ndims
  print "      typedef NEXT NextPolicy;"
  for dim in dim_names:
    print "      typedef TILE_%s Tile%s;" % (dim.upper(), dim.upper())
  print "    };"
  print ""
  
def writeForallPermutations(ndims):

  dim_names = getDimNames(ndims)
  
  print ""
  print "/******************************************************************"
  print " *  Permutations layer for forall%d()" % ndims
  print " ******************************************************************/"
  print ""
  
    
  # Loop over each permutation specialization
  perms = getDimPerms(dim_names)
  for perm in perms:
    # get enumeration name
    enum = getEnumName(perm)
  
    # print function declaration
    args = map(lambda a: "typename T"+a.upper(), dim_names)
    argstr = ", ".join(args)    
    print "      template<typename POLICY, %s, typename BODY>" % argstr
    
    args = map(lambda a: "T%s const &is_%s"%(a.upper(), a), dim_names)
    argstr = ", ".join(args)    
    print "      RAJA_INLINE void forall%d_permute(%s, %s, BODY const &body){" % (ndims, enum, argstr)    
    
    # Create executor
    args = map(lambda a: "typename POLICY::Policy"+(a.upper()), perm)
    polstr = ", ".join(args)
    args = map(lambda a: "T"+(a.upper()), perm)
    setstr = ", ".join(args)
    print "        Forall%dExecutor<%s, %s> exec;" % (ndims, polstr, setstr)
    
    # Call executor
    args = map(lambda a: "is_"+a, perm)
    setstr = ", ".join(args)
    args = map(lambda a: "int "+a, perm)
    idxstr = ", ".join(args)
    print "        exec(%s, RAJA_LAMBDA(%s){" % (setstr, idxstr)
    argstr = ", ".join(dim_names)  # NOT PERMUTED!
    print "          body(%s);" % argstr
    print "        });"
    print "      }"
         
    print ""
    

def writeForallExecutor(ndims):

  print ""
  print "/******************************************************************"
  print " *  Default Executor for forall%d()" % ndims
  print " ******************************************************************/"
  print ""

  dim_names = getDimNames(ndims)

  args = map(lambda a: "typename POLICY_"+(a.upper()), dim_names)
  polstr = ", ".join(args)
  args = map(lambda a: "typename T"+(a.upper()), dim_names)
  setstr = ", ".join(args)
  print "    template<%s, %s>" % (polstr, setstr)
  print "    class Forall%dExecutor {" % (ndims)
  print "      public:  "
  
  # Create default executor
  args = map(lambda a: "T%s const &is_%s"%(a.upper(), a), dim_names)
  argstr = ", ".join(args)  
  print "        template<typename BODY>"
  print "        inline void operator()(%s, BODY const &body) const {" % argstr
  print "          RAJA::forall<POLICY_I>(is_i, RAJA_LAMBDA(int i){"
  if ndims == 2:  # 2 dimension termination case:
    print "            RAJA::forall<POLICY_J>(is_j, RAJA_LAMBDA(int j){"
  else: # more than 2 dimensions, we just peel off the outer loop, and call an N-1 executor
    args = map(lambda a: "is_"+a, dim_names[1:])
    setstr = ", ".join(args)
    args = map(lambda a: "int "+a, dim_names[1:])
    idxstr = ", ".join(args)  
    print "            exec(%s, RAJA_LAMBDA(%s){" % (setstr, idxstr)
  
  argstr = ", ".join(dim_names)  
  print "              body(%s);" % argstr
  print "            });"
  print "          });"
  print "        }"
    
  # More than 2 dims: create nested ForallNExecutor  
  if ndims > 2:
    print ""
    args = map(lambda a: "POLICY_"+(a.upper()), dim_names[1:])
    polstr = ", ".join(args)
    args = map(lambda a: "T"+(a.upper()), dim_names[1:])
    argstr = ", ".join(args)
    print "      private:"
    print "        Forall%dExecutor<%s, %s> exec;" % (ndims-1, polstr, argstr)
  print "    };"
  print ""    


def writeForallOpenMP(ndims):

  dim_names = getDimNames(ndims)
  
  print ""
  print "/******************************************************************"
  print " *  OpenMP Auto-Collapsing Executors for forall%d()" % ndims
  print " ******************************************************************/"
  print ""
  print "#ifdef _OPENMP"
  print ""

  for omp_policy in ['omp_parallel_for_exec', 'omp_for_nowait_exec']:
    for depth in range(2,ndims+1):
    
      remainder_ndims = ndims - depth

      polargs = []
      setargs = []
      args =  map(lambda a: "typename POLICY_"+a.upper(), dim_names[depth:])
      args.extend(map(lambda a: "typename T"+a.upper(), dim_names[depth:]))
      argstr = ", ".join(args)   
      print "    // OpenMP Executor with collapse(%d) for %s" % (depth, omp_policy)
      print "    template<%s>" % argstr
      
      args =  map(lambda a: "RAJA::"+omp_policy, range(0,depth))
      args.extend(map(lambda a: "POLICY_"+a.upper(), dim_names[depth:]))
      args.extend(map(lambda a: "RAJA::RangeSegment", range(0,depth)))
      args.extend(map(lambda a: "T"+a.upper(), dim_names[depth:]))
      argstr = ", ".join(args)   
      print "    class Forall%dExecutor<%s> {" % (ndims, argstr)
      print "      public:  "
      
      # Create collapse(depth) executor function
      print "        template<typename BODY>"
      
      args = map(lambda a: "RAJA::RangeSegment const &is_"+ a, dim_names[0:depth])
      args.extend(map(lambda a: "T%s const &is_%s"%(a.upper(),a), dim_names[depth:ndims]))
      argstr = ", ".join(args)  
      print "        inline void operator()(%s, BODY const &body) const {" % argstr
  #    print "          printf(\"collapse(%d)\\n\");" % depth
      
      # get begin and end indices each collapsed RangeSegment
      for a in dim_names[0:depth]:
        print "          int const %s_start = is_%s.getBegin();" % (a,a)
        print "          int const %s_end   = is_%s.getEnd();" % (a,a)
        print ""
      
      # Generate nested collapsed for loops
      if omp_policy == 'omp_parallel_for_exec':
        print "#pragma omp parallel for schedule(static) collapse(%d)" % depth
      elif omp_policy == 'omp_for_nowait_exec':
        print "#pragma omp for schedule(static) collapse(%d) nowait" % depth
      indent = ""
      for d in dim_names[0:depth]:
        print "          %sfor(int %s = %s_start;%s < %s_end;++ %s){" % (indent, d, d, d, d, d)
        indent += "  "
      
      # No more inner loops, so call the loop body directly
      if remainder_ndims == 0:
        argstr = argstr = ", ".join(dim_names)
        print "          %sbody(%s);" % (indent, argstr)
      
      # Just one inner loop, so issue a RAJA::forall
      elif remainder_ndims == 1:      
        d = dim_names[depth]
        print "          %sRAJA::forall<POLICY_%s>(is_%s, RAJA_LAMBDA(int %s){" % (indent, d.upper(), d, d)
        argstr = argstr = ", ".join(dim_names)
        print "          %s  body(%s);" % (indent, argstr)
        print "          %s});" % (indent)
      
      # More than one inner loop, so call an inner executor
      else:      
        #    exec(is_j, is_k, is_l, RAJA_LAMBDA(int j, int k, int l){
        #      body(i, j, k, l);
        #    });
        
        args = map(lambda a: "is_"+a, dim_names[depth:])
        setstr = ", ".join(args)
        args = map(lambda a: "int "+a, dim_names[depth:])
        argstr = ", ".join(args)      
        print "          %sexec(%s, RAJA_LAMBDA(%s){" % (indent, setstr, argstr)
        argstr = argstr = ", ".join(dim_names)
        print "          %s  body(%s);" % (indent, argstr)
        print "          %s});" % (indent)
      
      # Close out collapsed loops
      argstr = "";
      for d in range(0,depth):
        argstr += "} "
      print "          %s" % argstr
      print "        }"
      
        
      # More than 2 dims: create nested ForallNExecutor
      if remainder_ndims >= 2:
        print ""
        args = map(lambda a: "POLICY_"+(a.upper()), dim_names[depth:])
        polstr = ", ".join(args)
        args = map(lambda a: "T"+(a.upper()), dim_names[depth:])
        argstr = ", ".join(args)
        print "      private:"
        print "        Forall%dExecutor<%s, %s> exec;" % (ndims-depth, polstr, argstr)
      print "    };"
      print ""    
      
  print ""
  print "#endif // _OPENMP"
  print ""

ndims = int(sys.argv[1])


# ACTUAL SCRIPT ENTRY:
print """//AUTOGENERATED BY genForallN.py
  
#ifndef __DOMAIN_FORALL%d_H__
#define __DOMAIN_FORALL%d_H__

#include<RAJA/RAJA.hxx>
#include<Domain/Tile.h>

""" % (ndims, ndims)


# Create the policy struct so the user can define loop policies
writeForallPolicy(ndims)

# Create the default executor
writeForallExecutor(ndims)

# Create the OpenMP collapse() executors
writeForallOpenMP(ndims)

# Create all permutation MUX's 
writeForallPermutations(ndims)

# Dump out the base function that the user calls directly
writeForallBase(ndims)

print """
  
#endif
"""

