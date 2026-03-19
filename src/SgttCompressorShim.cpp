#include "compressor/Compressor.hpp"

namespace sgtt::compressor {

std::string Compressor::getFileName()
{
   return getName();
}

std::string Compressor::getInfo()
{
   return "{}";
}

void Compressor::printInfo()
{
}

} // namespace sgtt::compressor
