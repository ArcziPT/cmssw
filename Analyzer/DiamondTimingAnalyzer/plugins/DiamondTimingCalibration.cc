#include "DiamondTimingCalibration.h"

std::ostream& operator<<(std::ostream& os, const DiamondTimingCalibration& data){
    os<<data.calib;
    return os;
}

std::ostream& operator<<(std::ostream& os, const ChannelKey& key){
    os<<"{planeKey="<<key.planeKey<<", channel="<<key.channel<<"}";
    return os;
}

std::ostream& operator<<(std::ostream& os, const PlaneKey& planeKey){
    os<<"{sector="<<planeKey.sector<<", station="<<planeKey.station<<", plane="<<planeKey.plane<<"}";
    return os;
}