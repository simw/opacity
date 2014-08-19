#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <map>
#include <string>
#include <vector>
#include <list>

/**
  Parameters: reads a file, retrieves and stores the settings in a map structure
  When we want to retrieve that setting, a call to GetParameters with the relevant
  ParamID returns a list of strings containing the relevant parameters.

  The file format: the file can contain anything. This object only looks at lines
  between "# Begin settings for " + StringID, and "# End settings for " + StringID
  where StringID was supplied to the constructor of the object.
  Between these lines, blank lines and lines beginning with # (comment lines) are ignored.
  The rest are taken as "ParamID param1 param2 ...", space separated values.

  ParseInputFile returns an integer. If this is > 0, then this specifies the number of
  settings lines that have been read. If < 0, this gives an error code.
  = 0 is also an error - no settings read in.

  TODO: extend so that multiple files can be read in, where subsequent files are
  specified by a line "+ ./furtherparams.txt" line. Must check for looping references.
  This would be useful for eg multiple jet energies and types, but where everything
  else is the same. One common settings file, one specialized.
  
  @author Simon Wicks <simon_wicks@yahoo.com>
*/
class Parameters
{
public:
  Parameters( std::string StringID_ );

  /**
  Read in a file which contains all the settings for the run
  @param FileName The filename to process
  @return Integer error code
  */
  int ParseInputFile( std::string FileName );
  /// Helper function for ParseInputFile - does one line supplied as LineString
  int ParseSingleLine( std::string LineString );
  ///
  bool IsParamInList( std::string ParamID );
  /// Get a specific parameter in the map, identified by ParamID
  std::list<std::string> GetParametersString( std::string ParamID );
  /// Get the parameter set, but then convert to a vector of doubles
  std::vector<double> GetParametersDouble( std::string ParamID );
  /// Get the parameter set, but then convert to a vector of longs
  std::vector<long> GetParametersLong( std::string ParamID );
  /** See whether all the parameters read in from file have been accessed by the program
     @return Gives 0 if all parameters used, otherwise gives the number of unaccessed
     parameters (ie >0)
  **/
  int CheckForUnaccessedParameters();
  /** Write our parameters to file, for record and possible reopening later
    @param FileName The filename of the file to write to
    @param append A true / false flag on whether to erase the current contents of the file
    of to append to it.
    @return Integer error code, 0 if ok
  */
  int WriteToFile( std::string FileName, bool append );

private:
  /// The unique ID to identify the settings section in the file
  std::string StringID;
  /// The map in which to store ParamID, list of parameters
  std::map<std::string, std::list<std::string> > TheParameters;
  /// The map in which to store whether ParamID has been accessed
  std::map<std::string, bool> AccessedList;
};

#endif
