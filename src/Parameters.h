// ============================================================================
// PARAMETERES.H
// ============================================================================


#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_


#include <map>
#include <string>


// ============================================================================
// CLASS PARAMETERS
// ============================================================================
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
  void trim2(std::string&);
  

  std::map <std::string,int> intparams;
  std::map <std::string,float> floatparams;
  std::map <std::string,std::string> stringparams;

};



#endif
