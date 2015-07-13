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
        if (v.second.get<std::string>("id") == "count")
          {
            in->coarseCount = v.second.get<int>("value");
            CkPrintf("count                             %d\n", in->coarseCount);
          }
        if (v.second.get<std::string>("id") == "use adaptive sampling")
          {
            in->useAdaptiveSampling = v.second.get<int>("value");
            CkPrintf("use adaptive sampling             %d\n", in->useAdaptiveSampling);
          }
      }
    }
  catch (std::exception const& e)
    {
      std::cerr << e.what() << std::endl;
      CkExit();
    }
}
