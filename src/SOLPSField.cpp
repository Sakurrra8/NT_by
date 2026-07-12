#include "SOLPSField.h"

#include <cmath>
#include <fstream>
#include <stdexcept>

SOLPSIndexedField ReadSOLPSIndexedField(
    const std::string &path, const std::string &required_header_text)
{
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Cannot open SOLPS indexed field: " + path);

    SOLPSIndexedField field;
    if (!std::getline(input, field.header))
        throw std::runtime_error("SOLPS indexed field is empty: " + path);
    if (field.header.find(required_header_text) == std::string::npos)
        throw std::runtime_error(
            "Unexpected SOLPS field header in " + path + ": " + field.header);

    SOLPSIndexedValue row{};
    while (input >> row.i >> row.j >> row.component >> row.value)
    {
        if (!std::isfinite(row.value))
            throw std::runtime_error("Non-finite value in SOLPS indexed field: " + path);
        field.values.push_back(row);
    }
    if (!input.eof())
        throw std::runtime_error("Malformed row in SOLPS indexed field: " + path);
    if (field.values.empty())
        throw std::runtime_error("No data rows in SOLPS indexed field: " + path);
    return field;
}
