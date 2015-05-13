/** parse input file for 2D macro solver, inlcudes
 * boost json file parsing
 * **/

#include <iostream>
#include "types.h"

#include <string>
#include <cstring>
#include <iostream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>

void parse_input(string input_file, Input *in)
{
  std::ifstream file(input_file.c_str());
  boost::property_tree::ptree pt;

  if (file)
    {
      std::stringstream buffer;
      buffer << file.rdbuf();

      file.close();
      boost::property_tree::read_json(buffer, pt);
    }
  try
    {
      BOOST_FOREACH(const boost::property_tree::ptree::value_type & v, pt.get_child("parameter.CoarseScaleModel"))
      {
        if (v.second.get<std::string>("id") == "element dim x")
          {
            in->elemDimX = v.second.get<int>("value");
            CkPrintf("set element dim x:                %d\n", in->elemDimX);
          }
        if (v.second.get<std::string>("id") == "element dim y")
          {
            in->elemDimY = v.second.get<int>("value");
            CkPrintf("set element dim y:                %d\n", in->elemDimY);
          }
        if (v.second.get<std::string>("id") == "element dim z")
          {
            in->elemDimZ = v.second.get<int>("value");
            CkPrintf("set element dim z:                %d\n", in->elemDimZ);
          }
        if (v.second.get<std::string>("id") == "block dim x")
          {
            in->blockDimX = v.second.get<int>("value");
            CkPrintf("set block dim x                   %d\n", in->blockDimX);
          }
        if (v.second.get<std::string>("id") == "block dim y")
          {
            in->blockDimY = v.second.get<int>("value");
            CkPrintf("set block dim y                   %d\n", in->blockDimY);
          }
        if (v.second.get<std::string>("id") == "block dim z")
          {
            in->blockDimZ = v.second.get<int>("value");
            CkPrintf("set block dim z                   %d\n", in->blockDimZ);
          }
        if (v.second.get<std::string>("id") == "max timesteps")
          {
            in->maxTimesteps = v.second.get<int>("value");
            CkPrintf("max timesteps:                    %d\n", in->maxTimesteps);
          }
      }
    }
  catch (std::exception const& e)
    {
      std::cerr << e.what() << std::endl;
      CkExit();
    }
}
