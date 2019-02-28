/*
 * @file StereoName.hpp
 *
 * @author Pierre-Anthony Grenier <pierre-anthony.grenier@ensicaen.fr>
 *
 * @version     0.0.1 - Ven Sep 13 2013
 *
 */

#ifndef __STEREONAME_H__
#define __STEREONAME_H__

class StereoName
{
  std::vector<int> _fromList;
  std::vector<int> _ringClosureList;
  std::vector<int> _atomTypeList;
  std::vector<int> _bondTypeList;
  std::vector<int> _ringClosureBondTypeList;
  std::vector<int> _stereoDoubleCarbonList;
  std::vector<int> _ringClosureStereoDoubleCarbonList;
  std::vector<int> _saturatedCarbonList;

public:
  StereoName(){};
  StereoName(const StereoName & c);

  void clear();

  void getName(std::vector<int> * name);
  
  void pushFL(int val){_fromList.push_back(val);};
  void pushRCL(int val){_ringClosureList.push_back(val);};
  void pushATL(int val){_atomTypeList.push_back(val);};
  void pushBTL(int val){_bondTypeList.push_back(val);};
  void pushRCBTL(int val){_ringClosureBondTypeList.push_back(val);};
  void pushSDCL(int val){_stereoDoubleCarbonList.push_back(val);};
  void pushRCSDCL(int val){_ringClosureStereoDoubleCarbonList.push_back(val);};
  void pushSCL(int val){_saturatedCarbonList.push_back(val);};
  
  int fromList(int i){return _fromList[i];};

  bool operator<(StereoName const &s);

  int nbDifference(StereoName const &s);

  bool difInStereoCarbon(StereoName const &s);

  bool difInStereoDoubleBond(StereoName const &s);

  int numDiffCarbon(StereoName const &s);

  int numDiffDoubleBond(StereoName const &s);

  void showName();

};

#endif // __STEREONAME_H__
