#include <sstream>
#include <fstream>
#include "NumberWithUnits.hpp"
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <stdexcept>
using namespace std;
constexpr double tolerance = 0.0001;
constexpr int increment = 5;

namespace ariel
{

    void NumberWithUnits::initTbl(std::vector<std::string> &units, Tbl &tbl)
    {
        vector<pair<string, int>> tbl_units;
        

        for (unsigned int i = 0; i < units.size(); i += increment)
        {
            std::string unit1 = units.at(i + 1);
            std::string unit2 = units.at(i + 4);

            if (tbl.umap.count(units.at(i + 1)) == 0 && tbl.umap.count(units.at(i + 4)) == 0)
            {
                insert2units(units, tbl, i);
            }
            else if (tbl.umap.count(units.at(i + 1)) == 0 || tbl.umap.count(units.at(i + 4)) == 0)
            {
                insert1unit(units, tbl, i);
            }
        }
    }

    NumberWithUnits::NumberWithUnits(double num, string const &unit)
    {
        std::unordered_map<std::string, uint>::const_iterator got = tbl.umap.find(unit);
        if (got == tbl.umap.end())
        {
            throw std::runtime_error("unit is not exist");
        }
        this->num = num;
        this->unit = unit;
    }
   
    void NumberWithUnits::read_units(std::ifstream &file)
    {
        if (file.is_open())
        {
            string s;
            while (file >> s)
            {
                units.push_back(s);
            }
        }
        initTbl(units, tbl);
    }
    
 
    uint NumberWithUnits::findIndex(Tbl &tbl, string const &unit)
    {
        uint line = tbl.umap.at(unit);
        for (size_t i = 0; i < tbl.matrix.at(line).size(); i++)
        {
            if (tbl.matrix.at(line).at(i).first == unit)
            {
                return i;
            }
        }
        return 0;
    }

    double NumberWithUnits::connect(Tbl &tbl, string const &left_unit, string const &right_unit)
    {
        if (left_unit == right_unit)
        {
            return 1;
        }
        double conversion = 0.0;
        unordered_map<string, uint> existing_units;
        existing_units[left_unit] = 0;

        if (tbl.umap.at(left_unit) == tbl.umap.at(right_unit))
        {
            uint line = tbl.umap.at(left_unit);
            vector<pair<string, vector<pair<string, double>>>> &ver_vec = tbl.matrix.at(line);

            uint ind = findIndex(tbl, left_unit);
            vector<pair<string, double>> &vec = ver_vec.at(ind).second; //KM_VEC

            int vec_length = vec.size();
            for (unsigned int i = 0; i < ver_vec.size(); i++)
            {
                string temp_unit = vec.at(i).first; // looking for the same unit
                if (temp_unit == right_unit)
                {
                    return 1 / vec.at(i).second;
                }
                conversion = vec.at(i).second;

                uint temp_ind = findIndex(tbl, temp_unit);
                vector<pair<string, double>> &vec_temp = ver_vec.at(temp_ind).second;
                unsigned int temp_length = vec_temp.size();

                for (unsigned int j = 0; j < temp_length; j++)
                {
                    string temp_name = vec_temp.at(j).first;
                    if (existing_units.find(temp_name) == existing_units.end())
                    {
                        conversion *= vec_temp.at(j).second;
                        vec.push_back(make_pair(temp_name, conversion));
                        existing_units[temp_name] = i;
                        vec_length++;
                        conversion /= vec_temp.at(j).second;
                    }
                }
            }
            for (unsigned int i = 0; i < vec.size(); i++)
            {
                if (vec.at(i).first == right_unit)
                {
                    return 1 / vec.at(i).second;
                }
            }
        } 
        else
        {
            const std::exception ex;
            throw(ex);
        }

        return 1 / conversion;
    }

   
    void NumberWithUnits::insert1unit(vector<string> &units, Tbl &tbl, uint i)
    {
        uint pointer = 0;
        double conv1 = stod(units.at(i + 3));
        string unit1 = units.at(i + 1);
        string unit2 = units.at(i + 4);
        pair<string, uint> left_unit;
        pair<string, uint> right_unit;
        pair<string, vector<pair<string, double>>> add_to_existing_vec; 
        vector<pair<string, double>> add_to_pair;                       

        if (tbl.umap.count(units.at(i + 4)) == 0)
        {
            pointer = tbl.umap.at(unit1); 
        }

        else
        {
            unit1 = units.at(i + 4);
            unit2 = units.at(i + 1);
            pointer = tbl.umap.at(unit1);
            conv1 = 1 / stod(units.at(i + 3));
        }

        pair<string, double> curr_conversion(unit1, 1 / conv1); 
        add_to_pair.push_back(curr_conversion);                 

       
        add_to_existing_vec.first = unit2;
        add_to_existing_vec.second = add_to_pair;

        tbl.matrix.at(pointer).push_back(add_to_existing_vec); 

        left_unit.first = unit2;
        left_unit.second = pointer;

        uint ind = 0;
        for (unsigned int i = 0; i < tbl.matrix.at(pointer).size(); i++)
        {
            if (tbl.matrix.at(pointer).at(i).first == unit1)
            {
                ind = i;
            }
        }

        pair<string, double> left_conversion(unit2, conv1); 
        tbl.matrix.at(pointer).at(ind).second.push_back(left_conversion);
        tbl.umap.insert(left_unit);
    }
   
    void NumberWithUnits::insert2units(vector<string> &units, Tbl &tbl, uint i)
    {
        double conv1 = stod(units.at(i + 3));
        double conv2 = 1 / stod(units.at(i + 3));
        string unit1 = units.at(i + 1);
        string unit2 = units.at(i + 4);

        vector<pair<string, vector<pair<string, double>>>> final_vec;

        vector<pair<string, double>> right_inside_units;
        vector<pair<string, double>> left_inside_units;

        pair<string, double> conversion;

        pair<string, double> right_conversion(units.at(i + 1), conv2); 

        pair<string, double> left_conversion(units.at(i + 4), conv1); 
        left_inside_units.push_back(left_conversion);
        right_inside_units.push_back(right_conversion);

        pair<string, vector<pair<string, double>>> left_pair(units.at(i + 1), left_inside_units);
        pair<string, vector<pair<string, double>>> right_pair(units.at(i + 4), right_inside_units);

        final_vec.push_back(left_pair);
        final_vec.push_back(right_pair);

        tbl.matrix.emplace_back(final_vec);

        pair<string, int> left_unit(units.at(i + 1), tbl.matrix.size() - 1);
        pair<string, int> right_unit(units.at(i + 4), tbl.matrix.size() - 1);
        tbl.umap.insert(left_unit);
        tbl.umap.insert(right_unit);
    }

    NumberWithUnits operator+(const NumberWithUnits &f1, const NumberWithUnits &f2)
    {
        double d = NumberWithUnits::connect(tbl, f1.unit, f2.unit);
        NumberWithUnits a(f1.num + f2.num * d, f1.unit);
        return a;
    }

    NumberWithUnits& NumberWithUnits::operator+=(const NumberWithUnits &f2)
    {
        double d = NumberWithUnits::connect(tbl, unit, f2.unit);
        num = num + f2.num * d;
        return *this;
    }

    NumberWithUnits operator+(const NumberWithUnits &f1)
    {
        return f1;
    }
    NumberWithUnits operator-(const NumberWithUnits &f1, const NumberWithUnits &f2)
    {
        double d = NumberWithUnits::connect(tbl, f1.unit, f2.unit);
        NumberWithUnits a(f1.num - f2.num * d, f1.unit);
        return a;
    }
    NumberWithUnits& NumberWithUnits::operator-=(const NumberWithUnits& f2)
    {
        double d = NumberWithUnits::connect(tbl, unit, f2.unit);
        num = num - f2.num * d;
        return *this;
    }

    NumberWithUnits operator-(const NumberWithUnits &f1)
    {
        NumberWithUnits a{-1 * f1.num, f1.unit};
        return a;
    }
    bool NumberWithUnits::operator<(const NumberWithUnits &f) const
    {
        double d = NumberWithUnits::connect(tbl, unit, f.unit);
        return !(*this == f) && num < f.num * d;
    }
    bool NumberWithUnits::operator<=(const NumberWithUnits &f) const
    {
        double d = NumberWithUnits::connect(tbl, unit, f.unit);
        return num < f.num * d || num == f.num * d;
    }
    bool NumberWithUnits::operator>(const NumberWithUnits &f) const
    {
        double d = NumberWithUnits::connect(tbl, unit, f.unit);
        return num > f.num * d;
    }
    bool NumberWithUnits::operator>=(const NumberWithUnits &f) const
    {
        double d = NumberWithUnits::connect(tbl, unit, f.unit);
        return num > f.num * d || num == f.num * d;
    }
    bool NumberWithUnits::operator==(const NumberWithUnits &f) const
    {
        double d = NumberWithUnits::connect(tbl, unit, f.unit);
        return ((-1) * tolerance < num - f.num * d && num - f.num * d < tolerance);
    }
    bool NumberWithUnits::operator!=(const NumberWithUnits &f) const
    {
        return !(*this == f);
    }
    NumberWithUnits &NumberWithUnits::operator++()
    {
        this->num += 1;
        return *this;
    } 
    NumberWithUnits NumberWithUnits::operator++(int)
    {
        NumberWithUnits temp = *this;
        ++(this->num);
        return temp;
    } 

    NumberWithUnits &NumberWithUnits::operator--()
    {
        this->num -= 1;
        return *this;
    } 
    NumberWithUnits NumberWithUnits::operator--(int)
    {
        NumberWithUnits temp = *this;
        --(this->num);
        return temp;
    } 
   
    NumberWithUnits operator*(const NumberWithUnits &f1, const double &f2)
    {
        NumberWithUnits a{f1.num * f2, f1.unit};
        return a;
    }
    NumberWithUnits operator*(const double &f1, const NumberWithUnits &f2)
    {
        NumberWithUnits a{f1 * f2.num, f2.unit};
        return a;
    }

    std::istream &operator>>(std::istream &is, NumberWithUnits &f)
    {
        double num = 0;
        string unit;
        char leftS = '.', rsign = '.';

        is >> skipws >> num >> leftS >> unit;
        if (unit.at(unit.length() - 1) == ']')
        {
            unit = unit.substr(0, unit.length() - 1);
        }
        else
        {
            is >> skipws >> rsign;
        }

        std::unordered_map<std::string, uint>::const_iterator got = tbl.umap.find(unit);
        if (got == tbl.umap.end())
        {
            throw std::runtime_error("this unit is not exist");
        }

        f.num = num;
        f.unit = unit;
        return is;
    }

    std::ostream &operator<<(std::ostream &os, const NumberWithUnits &f)
    {
        return os << f.num << "[" << f.unit << "]";
    }
    
}
