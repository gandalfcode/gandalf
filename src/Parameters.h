//=============================================================================
//  Parameters.h
//=============================================================================


#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_


#include <map>
#include <string>


//=============================================================================
//  Class Parameters
/// \brief   Class for continaing all simuilation parameter information.
/// \details Class for continaing all simuilation parameter information.
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=============================================================================
class Parameters
{
 public:

  Parameters();
  ~Parameters();
  Parameters(const Parameters&);

  void ReadParamsFile(std::string);
  void ParseLine(std::string);
  void SetDefaultValues(void);
  std::string GetParameter (std::string);
  void SetParameter(std::string , std::string);
  void PrintParameters(void);
  void RecordParametersToFile(void);
  std::string TrimWhiteSpace(std::string);
  

  std::map <std::string,int> intparams;
  std::map <std::string,float> floatparams;
  std::map <std::string,std::string> stringparams;

};

#endif
