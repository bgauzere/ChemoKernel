/*
 * @file StereoName.cpp
 *
 * @author Pierre-Anthony Grenier <pierre-anthony.grenier@ensicaen.fr>
 *   
 * @version     0.0.1 - Ven Sep 13 2013
 * 
 */

#include <vector>
#include <iostream>
#include "StereoName.hpp"
using namespace std;

StereoName::StereoName(const StereoName & s)
{
  this->clear();
  for(unsigned int i=0;i<s._fromList.size();i++)
    this->_fromList.push_back(s._fromList[i]);
  for(unsigned int i=0;i<s._ringClosureList.size();i++)
    this->_ringClosureList.push_back(s._ringClosureList[i]);
  for(unsigned int i=0;i<s._atomTypeList.size();i++)
    this->_atomTypeList.push_back(s._atomTypeList[i]);
  for(unsigned int i=0;i<s._bondTypeList.size();i++)
    this->_bondTypeList.push_back(s._bondTypeList[i]);
 for(unsigned int i=0;i<s._ringClosureBondTypeList.size();i++)
    this->_ringClosureBondTypeList.push_back(s._ringClosureBondTypeList[i]);
 for(unsigned int i=0;i<s._stereoDoubleCarbonList.size();i++)
    this->_stereoDoubleCarbonList.push_back(s._stereoDoubleCarbonList[i]);
  for(unsigned int i=0;i<s._ringClosureStereoDoubleCarbonList.size();i++)
    this->_ringClosureStereoDoubleCarbonList.push_back(s._ringClosureStereoDoubleCarbonList[i]);
   for(unsigned int i=0;i<s._saturatedCarbonList.size();i++)
    this->_saturatedCarbonList.push_back(s._saturatedCarbonList[i]);
}

void StereoName::clear()
{
   _fromList.clear();
   _ringClosureList.clear();
   _atomTypeList.clear();
   _bondTypeList.clear();
   _ringClosureBondTypeList.clear();
   _stereoDoubleCarbonList.clear();
   _ringClosureStereoDoubleCarbonList.clear();
   _saturatedCarbonList.clear();
}

void StereoName::getName(std::vector<int> * name)
{
  name->clear();
  name->insert( name->end(), _fromList.begin(), _fromList.end());
  name->insert( name->end(), _ringClosureList.begin(), _ringClosureList.end());
  name->insert( name->end(), _atomTypeList.begin(), _atomTypeList.end());  
  name->insert( name->end(), _bondTypeList.begin(), _bondTypeList.end());  
  name->insert( name->end(), _ringClosureBondTypeList.begin(), _ringClosureBondTypeList.end());
  name->insert( name->end(), _stereoDoubleCarbonList.begin(), _stereoDoubleCarbonList.end());
  name->insert( name->end(), _ringClosureStereoDoubleCarbonList.begin(), _ringClosureStereoDoubleCarbonList.end());
  name->insert( name->end(), _saturatedCarbonList.begin(), _saturatedCarbonList.end());
}

bool StereoName::operator<(StereoName const &s)
{
  for(unsigned int i=0;i<_fromList.size();i++)
    {
      if(this->_fromList[i]<s._fromList[i])
	return true;
      else if(this->_fromList[i]>s._fromList[i])
	return false;
    }
  for(unsigned int i=0;i<_ringClosureList.size();i++)
    {
      if(this->_ringClosureList[i]<s._ringClosureList[i])
	return true;
      else if(this->_ringClosureList[i]>s._ringClosureList[i])
	return false;
    }
  for(unsigned int i=0;i<_atomTypeList.size();i++)
    {
      if(this->_atomTypeList[i]<s._atomTypeList[i])
	return true;
      else if(this->_atomTypeList[i]>s._atomTypeList[i])
	return false;
    }
  for(unsigned int i=0;i<_bondTypeList.size();i++)
    {
      if(this->_bondTypeList[i]<s._bondTypeList[i])
	return true;
      else if(this->_bondTypeList[i]>s._bondTypeList[i])
	return false;
    }
  for(unsigned int i=0;i<_ringClosureBondTypeList.size();i++)    
    {
      if(this->_ringClosureBondTypeList[i]<s._ringClosureBondTypeList[i])
	return true;
      else if(this->_ringClosureBondTypeList[i]>s._ringClosureBondTypeList[i])
	return false;
    }

  for(unsigned int i=0;i<_stereoDoubleCarbonList.size();i++)
    {
      if(this->_stereoDoubleCarbonList[i]<s._stereoDoubleCarbonList[i])
	return true;
      else if(this->_stereoDoubleCarbonList[i]>s._stereoDoubleCarbonList[i])
	return false;
    }

  for(unsigned int i=0;i<_ringClosureStereoDoubleCarbonList.size();i++)
    {
      if(this->_ringClosureStereoDoubleCarbonList[i]<s._ringClosureStereoDoubleCarbonList[i])
	return true;
      else if(this->_ringClosureStereoDoubleCarbonList[i]>s._ringClosureStereoDoubleCarbonList[i])
	return false;
    }

  for(unsigned int i=0;i<_saturatedCarbonList.size();i++)
    {
      if(this->_saturatedCarbonList[i]<s._saturatedCarbonList[i])
	return true;
      else if(this->_saturatedCarbonList[i]>s._saturatedCarbonList[i])
	return false;
    }

     return false;
}

int StereoName::nbDifference(StereoName const &s)
{
  int nbDiff=0;
  for(unsigned int i=0;i<s._fromList.size();i++)
    if(this->_fromList[i]!=s._fromList[i])
      nbDiff++;
  for(unsigned int i=0;i<s._ringClosureList.size();i++)
    if(this->_ringClosureList[i]!=s._ringClosureList[i])
      nbDiff++;
  for(unsigned int i=0;i<s._atomTypeList.size();i++)
    if(this->_atomTypeList[i]!=s._atomTypeList[i])
      nbDiff++;
  for(unsigned int i=0;i<s._bondTypeList.size();i++)
    if(this->_bondTypeList[i]!=s._bondTypeList[i])
      nbDiff++;
  for(unsigned int i=0;i<s._ringClosureBondTypeList.size();i++)
    if(this->_ringClosureBondTypeList[i]!=s._ringClosureBondTypeList[i])
      nbDiff++;
  for(unsigned int i=0;i<s._stereoDoubleCarbonList.size();i++)
    if(this->_stereoDoubleCarbonList[i]!=s._stereoDoubleCarbonList[i])
      nbDiff++;
  for(unsigned int i=0;i<s._ringClosureStereoDoubleCarbonList.size();i++)
    if(this->_ringClosureStereoDoubleCarbonList[i]!=s._ringClosureStereoDoubleCarbonList[i])
      nbDiff++;
  for(unsigned int i=0;i<s._saturatedCarbonList.size();i++)
    if(this->_saturatedCarbonList[i]!=s._saturatedCarbonList[i])
      nbDiff++;
  return nbDiff;
}

bool StereoName::difInStereoCarbon(StereoName const &s)
{
  for(unsigned int i=0;i<s._saturatedCarbonList.size();i++)
    if(this->_saturatedCarbonList[i]!=s._saturatedCarbonList[i])
      return true;
  return false;
}

bool StereoName::difInStereoDoubleBond(StereoName const &s)
{
  for(unsigned int i=0;i<s._ringClosureStereoDoubleCarbonList.size();i++)
    if(this->_ringClosureStereoDoubleCarbonList[i]!=s._ringClosureStereoDoubleCarbonList[i])
      return true;
  for(unsigned int i=0;i<s._stereoDoubleCarbonList.size();i++)
    if(this->_stereoDoubleCarbonList[i]!=s._stereoDoubleCarbonList[i])
      return true;
  return false;
}

int StereoName::numDiffCarbon(StereoName const &s)
{
  int pos=-1;
  for(unsigned int i=0;i<s._saturatedCarbonList.size() && pos==-1;i++)
    if(this->_saturatedCarbonList[i]!=s._saturatedCarbonList[i])
      pos=i;
  return pos;
}

int StereoName::numDiffDoubleBond(StereoName const &s)
{
  int pos=-1;
  for(unsigned int i=0;i<s._ringClosureStereoDoubleCarbonList.size() && pos==-1;i++)
    if(this->_ringClosureStereoDoubleCarbonList[i]!=s._ringClosureStereoDoubleCarbonList[i])
      pos=_ringClosureList[i*2];
  for(unsigned int i=0;i<s._stereoDoubleCarbonList.size() && pos==-1;i++)
    if(this->_stereoDoubleCarbonList[i]!=s._stereoDoubleCarbonList[i])
      pos=i;
  return pos;
}

void StereoName::showName()
{
  cout<<"From List:";
  for(unsigned int i=0;i<_fromList.size();i++)
    cout<<" "<<_fromList[i];
  cout<<endl<<"Ring Closure List:";
  for(unsigned int i=0;i<_ringClosureList.size();i++)
    cout<<" "<<_ringClosureList[i];
  cout<<endl<<"Atom Type List:";
  for(unsigned int i=0;i<_atomTypeList.size();i++)
    cout<<" "<<_atomTypeList[i];
  cout<<endl<<"Bond Type List:";
  for(unsigned int i=0;i<_bondTypeList.size();i++)
    cout<<" "<<_bondTypeList[i];
  cout<<endl<<"Ring Closure Bond Type List:";
  for(unsigned int i=0;i<_ringClosureBondTypeList.size();i++)    
        cout<<" "<<_ringClosureBondTypeList[i];
  cout<<endl<<"Stereo Double Carbon List:";
  for(unsigned int i=0;i<_stereoDoubleCarbonList.size();i++)
            cout<<" "<<_stereoDoubleCarbonList[i];
  cout<<endl<<"Ring Closure Stereo Double Carbon List:";
  for(unsigned int i=0;i<_ringClosureStereoDoubleCarbonList.size();i++)
    cout<<" "<<_ringClosureStereoDoubleCarbonList[i];
  cout<<endl<<"Saturated Carbon List:";
  for(unsigned int i=0;i<_saturatedCarbonList.size();i++)
        cout<<" "<<_saturatedCarbonList[i];
  cout<<endl;
}
