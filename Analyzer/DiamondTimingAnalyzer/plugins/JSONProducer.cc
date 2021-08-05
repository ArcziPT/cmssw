#include "JSONProducer.h"

namespace pt = boost::property_tree;

void JSON::save(const CTPPSGeometry& geom,
	const DiamondTimingCalibration& calib,
	std::map<ChannelKey, double>& res,
	std::string output_file_name){

	std::ofstream output_file(output_file_name.c_str());
	if (!output_file.is_open()){
		std::cout << "ERROR: Can't open output file " << output_file_name << std::endl;
		return;
	}

	std::map<PlaneKey, pt::ptree> plane_node_map;
	for(auto it = geom.beginSensor(); it != geom.endSensor(); ++it){
		if (!CTPPSDiamondDetId::check(it->first))
			continue;
		
		const CTPPSDiamondDetId detid(it->first);
		ChannelKey key(detid);

		//create channel node
		pt::ptree ch_node, param_node;
		for(auto &param : calib.parameters(key)){
			pt::ptree base_node;
			base_node.put("", param);
			param_node.push_back(std::make_pair("", base_node));
		}

		ch_node.put("channel", key.channel);
		ch_node.put("time_offset", calib.timeOffset(key));
		ch_node.put("time_precision", res[key]);
		ch_node.add_child("param", param_node);

		plane_node_map[key.planeKey].push_back(std::make_pair("", ch_node));
	}

	std::map<std::pair<int, int>, pt::ptree> stations_node_map;
	for(auto& it : plane_node_map){
		auto& key = it.first;
		auto& plane = stations_node_map[{key.sector, key.station}];

		pt::ptree base_node;
		base_node.put("plane", key.plane);
		base_node.add_child("Channels", it.second);
		plane.push_back(std::make_pair("", base_node));
	}

	std::map<int, pt::ptree> sector_node_map;
	for(auto& it : stations_node_map){
		auto sec = it.first.first;
		auto st = it.first.second;
		auto& station = sector_node_map[sec];

		pt::ptree base_node;
		base_node.put("station", st);
		base_node.add_child("Planes", it.second);
		station.push_back(std::make_pair("", base_node));
	}

	pt::ptree sectors;
	for(auto& it : sector_node_map){
		auto sec = it.first;

		pt::ptree base_node;
		base_node.put("sector", sec);
		base_node.add_child("Stations", it.second);
		sectors.push_back(std::make_pair("", base_node));
	}
	
	pt::ptree parameters;
	parameters.add_child("Sectors", sectors);
	
	pt::ptree root;
	root.put("formula", "[0]/(exp((x-[1])/[2])+1)+[3]");
	root.add_child("Parameters", parameters);

	pt::write_json(output_file, root);
}
