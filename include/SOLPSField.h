#ifndef SOLPS_FIELD_H
#define SOLPS_FIELD_H

#include <string>
#include <vector>

struct SOLPSIndexedValue
{
    int i;
    int j;
    int component;
    double value;
};

struct SOLPSIndexedField
{
    std::string header;
    std::vector<SOLPSIndexedValue> values;
};

SOLPSIndexedField ReadSOLPSIndexedField(
    const std::string &path, const std::string &required_header_text);

#endif
