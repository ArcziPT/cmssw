#include "DiamondTimingCalibration.h"

std::ostream& operator<<(std::ostream& os, const DiamondTimingCalibration& data){
    os<<data.calib;
    return os;
}

std::ostream& operator<<(std::ostream& os, const ChannelKey& key){
    os<<"{"<<key.planeKey<<", "<<key.channel<<"}";
    return os;
}

std::ostream& operator<<(std::ostream& os, const PlaneKey& planeKey){
    os<<"{"<<planeKey.sector<<", "<<planeKey.station<<", "<<planeKey.plane<<"}";
    return os;
}