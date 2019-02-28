/**
 * @file HorvathKernel.h
 * @author Benoit Gauzere <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Mon Mar 19 2012
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 */

#ifndef __HORVATHKERNEL_H__
#define __HORVATHKERNEL_H__
#include "pandore.h"
#include "GraphKernel.h"
class HorvathKernel  : public GraphKernel
{
public:
  double operator()(pandore::Collection* c1, pandore::Collection* c2);
};

#endif // __HORVATHKERNEL_H__
