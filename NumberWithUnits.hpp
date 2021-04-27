#ifndef NUMBERWITHUNITS_H
#define NUMBERWITHUNITS_H
#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <sstream>
#include <stdexcept>
using namespace std;

namespace ariel
{

    struct Tbl
    {
        vector<vector<pair<string, vector<pair<string, double>>>>> matrix;
        unordered_map<string, uint> umap;
        unordered_map<string, uint> lomap;
    };

    static vector<string> units;
    static Tbl tbl;

    class NumberWithUnits
    {

    public:
        string unit;
        double num;
        NumberWithUnits(double num, string const &unit);

        NumberWithUnits();


        static void read_units(std::ifstream &file);

        NumberWithUnits& operator+=(const NumberWithUnits &f2);
        NumberWithUnits& operator-=(const NumberWithUnits &f2);

        friend NumberWithUnits operator-(const NumberWithUnits &f1);
        friend NumberWithUnits operator+(const NumberWithUnits &f1);

        friend NumberWithUnits operator+(const NumberWithUnits &f1, const NumberWithUnits &f2);
        friend NumberWithUnits operator-(const NumberWithUnits &f1, const NumberWithUnits &f2);

        bool operator<(const NumberWithUnits &f) const;
        bool operator>(const NumberWithUnits &f) const;

        bool operator<=(const NumberWithUnits &f) const;
        bool operator>=(const NumberWithUnits &f) const;
        
        bool operator==(const NumberWithUnits &f) const;
        bool operator!=(const NumberWithUnits &f) const;

        NumberWithUnits &operator++();   // prefix: ++a
        NumberWithUnits operator++(int); // postfix: a++

        NumberWithUnits &operator--();   // prefix: --a
        NumberWithUnits operator--(int); // postfix: a--

        friend NumberWithUnits operator*(const NumberWithUnits &f1, const double &f2);
        friend NumberWithUnits operator*(const double &f1, const NumberWithUnits &f2);

        friend std::ostream &operator<<(std::ostream &os, const NumberWithUnits &f);
        friend std::istream &operator>>(std::istream &is, NumberWithUnits &f);

        static void initTbl(vector<string> &units, Tbl &tbl);
        static uint findIndex(Tbl &tbl, string const &unit);
        static double connect(Tbl &tbl, string const &left_unit, std::string const &right_unit);
        static void insert2units(vector<string> &units, Tbl &tbl, uint i);
        static void insert1unit(vector<string> &units, Tbl &tbl, uint i);
    };
}

#endif
