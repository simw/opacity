#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include "parameters.h"

using namespace std;

Parameters::Parameters( string StringID_ )
{
  // Get ride of any whitespace around the given string
  boost::trim( StringID_ );
  StringID = StringID_;
}

int Parameters::ParseInputFile( string FileName )
{
  // First, we open up the file, check its existence
  ifstream FileIn;
  FileIn.open( FileName.c_str(), ios::in );  
  if ( FileIn.fail() )
  {
    cerr << "Parameters, unable to open file " << FileName << " for reading." << endl;
    return -20;
  }

  // Each line, we'll read into LineReadIn
  // Each line of settings will be counted in NumberOfLines
  string LineReadIn; long NumberOfLines = 0;

  // The string that identifies the beginning of the settings section
  string BeginString = "# Begin settings for " + StringID;
  // The string that identifies the end of the settings section
  string EndString = "# End setting for " + StringID;

  // Indicator whether we have found the settings section yet
  bool FoundStart = false;
  // Run through the file until we find the settings section
  do
  {
    if ( FileIn.eof() )
    {
      cerr << "Parameters, reached end of file at line before "; 
      cerr << "finding settings section: " << StringID << endl;
      return -21;
    }
    getline( FileIn, LineReadIn );
    boost::trim( LineReadIn );
    if ( LineReadIn == BeginString )
      FoundStart = true;
  }
  while ( !FoundStart );

  // Now we read each line, check whether it begins with a #, if not then we
  // pass it on to ParseSingleLine
  // When we parse each line, we will get a return value from the function ParseSingleLine
  int SingleLineReturnValue;

  while (! FileIn.eof() )
  {
    getline ( FileIn, LineReadIn );
    // First, trim the line of whitespace at each end
    boost::trim( LineReadIn );

    // If we've found the ending phrase, exit the while loop
    if ( LineReadIn == EndString )
      break;

    // If it's a blank line, or starts with a comment character #, skip this line
    if ( !(LineReadIn.size() == 0 || LineReadIn[0] == '#') )
    {
      SingleLineReturnValue = ParseSingleLine( LineReadIn );
      if ( SingleLineReturnValue != 0 )
      {
        // We've hit an error in the ParseSingleLine function
        FileIn.close();
        return SingleLineReturnValue;
      }
      ++NumberOfLines;
    }
  }

  // Now we're done, we close the file
  FileIn.close();
  // Everything has worked fine. Return a success code - the number of settings lines parsed
  return NumberOfLines;
}

int Parameters::ParseSingleLine( string LineString )
{
  string ParamID;
  list<string> ParametersIn;

  // Use the boost library to split the string 'LineString' into sections, with
  // a space delimiter
  boost::split( ParametersIn, LineString, boost::is_any_of(" "));

  if ( ParametersIn.size() <= 1 )
  {
    // The length of a line must be at least a defining token and a parameter setting
    // value. If the length is <=1, then the format is wrong
    cerr << "Parameters, line length is too short: " << LineString << endl;
    return -3;
  }

  // The first element is the ParamID. Get it and then remove it from the Parameters list
  ParamID = ParametersIn.front();
  ParametersIn.pop_front();

  // First, we check whether we've already used this parameter
  map<string, list<string> >::iterator i = TheParameters.find( ParamID );
  if ( i != TheParameters.end() )
  {
    cerr << "Parameters, duplicate parameter read in: " << ParamID << endl;
    return -4;
  }

  // Now we have all the information to set up a 'parameter'
  // Put this into our map
  TheParameters.insert( pair<string,list<string> >( ParamID, ParametersIn ) );
  AccessedList.insert( pair<string,bool>( ParamID, false ) );

  // All is ok, return ok
  return 0;
}

bool Parameters::IsParamInList( string ParamID )
{
  // We use the 'map' structure to turn the ParamID into a list of parameters
  map<string, list<string> >::iterator i = TheParameters.find( ParamID );

  // If we're at the end of the parameters map, then we haven't found it
  if ( i == TheParameters.end() )
  {
    return false;
  }

  return true;
}

list<string> Parameters::GetParametersString( string ParamID )
{
  // We use the 'map' structure to turn the ParamID into a list of parameters
  map<string, list<string> >::iterator i = TheParameters.find( ParamID );

  // If we can't find the ParamID in the map, then we return an error code -1
  if ( i == TheParameters.end() )
  {
    list<string> EmptyList;
    std::cerr << "ParamParser, cannot find ParamID in the map: " << ParamID << std::endl;
    return EmptyList;
  }

  // We've found it. Set accessed to true (we don't mind multiple accesses)
  map<string, bool>::iterator j = AccessedList.find( ParamID );
  j->second = true;

  // Ok, we've found it, now we get the list<string> from the iterator i, and return it
  return (i->second);
}

vector<double> Parameters::GetParametersDouble( string ParamID )
{
  list<string> StringList = GetParametersString( ParamID );

  vector<double> ReturnVector;
  double ConvertedString;

  for ( list<string>::iterator i = StringList.begin();
    i != StringList.end(); ++i )
  {
    ConvertedString = boost::lexical_cast<double>( *i );
    ReturnVector.push_back( ConvertedString );
  }

  // We've found it. Set accessed to true (we don't mind multiple accesses)
  map<string, bool>::iterator j = AccessedList.find( ParamID );
  j->second = true;

  return ReturnVector;
}

vector<long> Parameters::GetParametersLong( string ParamID )
{
  list<string> StringList = GetParametersString( ParamID );

  vector<long> ReturnVector;
  long ConvertedString;

  for ( list<string>::iterator i = StringList.begin();
    i != StringList.end(); ++i )
  {
    ConvertedString = boost::lexical_cast<long>( *i );
    ReturnVector.push_back( ConvertedString );
  }

  // We've found it. Set accessed to true (we don't mind multiple accesses)
  map<string, bool>::iterator j = AccessedList.find( ParamID );
  j->second = true;

  return ReturnVector;
}

int Parameters::CheckForUnaccessedParameters()
{
  // How many unaccessed parameters have we found so far
  int ReturnValue = 0;

  // Now iterate over all parameters in the map
  // (the ParamIDs in AccessedList should be the same as in TheParameters)
  for ( map<string, bool>::iterator i = AccessedList.begin();
        i != AccessedList.end(); ++i )
  {
    if ( i->second != true )
    {
      ++ReturnValue;
      cerr << "Unused parameter: " << i->first << "\n";
    }
  }

  return ReturnValue;
}

int Parameters::WriteToFile( string FileName, bool Append )
{
  // First, open the file
  std::ofstream FileOut;
  if ( Append )
  {
    FileOut.open( FileName.c_str(), ios::app );
  }
  else
  {
    FileOut.open( FileName.c_str(), ios::trunc );
  }
  
  // Check for errors
  if ( FileOut.fail() )
  {
    cerr << "Parameters, unable to open file " << FileName << " for reading." << endl;
    return -20;
  }

  // First, write out our ID string
  FileOut << "# Begin settings for " << StringID << "\n";

  // Now, go through the map and write out all the settings!
  map<string, list<string> >::const_iterator i;
  list<string>::const_iterator j;
  for ( i = TheParameters.begin(); i != TheParameters.end(); ++i )
  {
    // First, we have the ParamID string
    FileOut << i->first << " ";
    // Now we go through each setting in the list
    for ( j = (i->second).begin(); j != (i->second).end(); ++j )
    {
      FileOut << *j << " ";
    }
    FileOut << "\n";
  }

  // Now, write out our ending ID string
  FileOut << "# End settings for " << StringID << "\n";

  // We're done
  FileOut.close();

  return 0;
}
