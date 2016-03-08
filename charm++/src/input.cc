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
  } else {
      std::cerr << "Could not open json file " << input_file << std::endl;
      CkExit();
 }
  try
    {
      BOOST_FOREACH(const boost::property_tree::ptree::value_type & v, pt.get_child("parameter.CoarseScaleModel"))
      {
        if (v.second.get<std::string>("id") == "type")
          {
            in->coarseType = v.second.get<int>("value");
            CkPrintf("coarse model type                 %d\n", in->coarseType);
          }
        if (v.second.get<std::string>("id") == "count")
          {
            in->coarseCount = v.second.get<int>("value");
            CkPrintf("coarse model count                %d\n", in->coarseCount);
          }
        if (v.second.get<std::string>("id") == "use adaptive sampling")
          {
            in->useAdaptiveSampling = v.second.get<int>("value");
            CkPrintf("use adaptive sampling             %d\n", in->useAdaptiveSampling);
          }
        if (v.second.get<std::string>("id") == "stop time")
          {
            in->stopTime = v.second.get<Real_t>("value");
            CkPrintf("stop time                         %e\n", in->stopTime);
          }
        if (v.second.get<std::string>("id") == "file parts")
          {
            in->file_parts = v.second.get<int>("value");
            CkPrintf("file parts                        %e\n", in->file_parts);
          }
        if (v.second.get<std::string>("id") == "visit data interval")
          {
            in->visit_data_interval = v.second.get<int>("value");
            CkPrintf("visit data interval               %e\n", in->visit_data_interval);
          }
        if (v.second.get<std::string>("id") == "edge elems")
          {
            in->edgeElems = v.second.get<int>("value");
            CkPrintf("Number of Edge ELems               %e\n", in->edgeElems);
          }
        if (v.second.get<std::string>("id") == "height elems")
          {
            in->heightElems = v.second.get<int>("value");
            CkPrintf("Number of Height ELems               %e\n", in->heightElems);
          }
      }

      BOOST_FOREACH(const boost::property_tree::ptree::value_type & v, pt.get_child("parameter.FineScaleModel"))
      {
        if (v.second.get<std::string>("id") == "type")
          {
            in->fineType = v.second.get<int>("value");
            CkPrintf("fine model type                 %d\n", in->fineType);
          }
      }

      BOOST_FOREACH(const boost::property_tree::ptree::value_type & v, pt.get_child("parameter.NearestNeighborSearch"))
      {
        if (v.second.get<std::string>("id") == "type")
          {
            in->nnsType = v.second.get<int>("value");
            CkPrintf("nns type                        %d\n", in->nnsType);
          }
        if (v.second.get<std::string>("id") == "count")
          {
            in->nnsCount = v.second.get<int>("value");
            CkPrintf("nns count                       %d\n", in->nnsCount);
          }
        if (v.second.get<std::string>("id") == "point dimension")
          {
            in->pointDim = v.second.get<int>("value");
            CkPrintf("nns point dimension             %d\n", in->pointDim);
          }
        if (v.second.get<std::string>("id") == "number of trees")
          {
            in->numTrees = v.second.get<int>("value");
            CkPrintf("nns number of trees             %d\n", in->numTrees);
          }
       }

       BOOST_FOREACH(const boost::property_tree::ptree::value_type & v, pt.get_child("parameter.Interpolate"))
      {
        if (v.second.get<std::string>("id") == "type")
          {
            in->interpType = v.second.get<int>("value");
            CkPrintf("interpolate type                 %d\n", in->interpType);
          }
        if (v.second.get<std::string>("id") == "count")
          {
            in->interpCount = v.second.get<int>("value");
            CkPrintf("interpolate count                %d\n", in->interpCount);
          }
       }

       BOOST_FOREACH(const boost::property_tree::ptree::value_type & v, pt.get_child("parameter.Evaluate"))
      {     
        if (v.second.get<std::string>("id") == "type")
          {
            in->evalType = v.second.get<int>("value");
            CkPrintf("evaluate type                 %d\n", in->evalType);
          } 
        if (v.second.get<std::string>("id") == "count")
          {
            in->evalCount = v.second.get<int>("value");
            CkPrintf("evaluate count                %d\n", in->evalCount);
          } 
       }  

       BOOST_FOREACH(const boost::property_tree::ptree::value_type & v, pt.get_child("parameter.DBInterface"))
      { 
        if (v.second.get<std::string>("id") == "type")
          {
            in->dbType = v.second.get<int>("value");
            CkPrintf("DB type                 %d\n", in->dbType);
          }
        if (v.second.get<std::string>("id") == "count")
          {
            in->dbCount = v.second.get<int>("value");
            CkPrintf("DB count                %d\n", in->dbCount);
          }
        if (v.second.get<std::string>("id") == "remote")
          {
            in->dbRemote = v.second.get<int>("value");
            CkPrintf("Remote DB                %d\n", in->dbRemote);
          }
       }

    }
  catch (std::exception const& e)
    {
      std::cerr << e.what() << std::endl;
      CkExit();
    }
}
