/*
 * SCC_GeometricEntity.h
 *
 *  Created on: Dec 29, 2021
 *      Author: anderson
 */
/*
#############################################################################
#
# Copyright 2021- Chris Anderson
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# For a copy of the GNU General Public License see
# <http://www.gnu.org/licenses/>.
#
#############################################################################
*/

//
// Virtual base class for 2D geometric entities
//

#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>

#include "SCC_GeoType.h"

#ifndef SCC_GEOENTITY_
#define SCC_GEOENTITY_
namespace SCC
{
class  GeometricEntity
{

public :

//
//  Constructors
//
    GeometricEntity(){};
    GeometricEntity(const GeometricEntity&){};
//
//  Destructor
//
    virtual  ~GeometricEntity(){};
//
//  output functions
//
    friend std::ostream&  operator <<(std::ostream& out_stream, const GeometricEntity& A)
    {A.output(out_stream);return out_stream;};
    virtual void  output(std::ostream&) const {};

    friend std::istream&  operator >>(std::istream& in_stream, GeometricEntity& A)
    {A.input(in_stream);return in_stream;};
    
    virtual void input(std::istream&) {};

//
//  Geometric Entity Member Functions
//
    virtual SCC::GeoType getEntityType() const = 0;
    virtual int  interiorExteriorTest(double, double) const {return 0;};
    virtual int  boundaryTest(double, double) const {return 0;};
    virtual int  boundaryTest(double, double, double, double) const {return 0;};
    virtual void  scaleBy(double){};
    virtual void  rotateBy(double){};
    virtual void  translateBy(double, double){};
    virtual void  translateTo(double, double){};
    virtual double  getMinX() const {return 0.0;};
    virtual double  getMaxX() const {return 0.0;};
    virtual double  getMinY() const {return 0.0;};
    virtual double  getMaxY() const {return 0.0;};
    virtual double  getTotalArcLength() const {return 0.0;};
    virtual int  getParametricCoordinate(double& /* coordinate */, double, double) const {return 0;};
    virtual int  getParametricCoordinate(double& /* coordinate */, double, double, double, double) const {return 0;};
    virtual double  getXcoordinate(double) const {return 0.0;};
    virtual double  getYcoordinate(double) const {return 0.0;};
    virtual int  getInteriorPoint(double& x, double& y) const {x = 0.0; y = 0.0; return 0;};
    virtual int  getSegmentIntersection(double& intersectPoint, double,
				  double, double, double) const {intersectPoint = 0.0; return 0;};
    virtual int  getSegmentIntersection(double& intersectPoint, double,
				  double, double, double, double, double) const {intersectPoint = 0.0; return 0;};
    virtual void  getUnitNormal(double, double& n_x, double& n_y) const {n_x =0.0; n_y = 0.0;};
    virtual void  getUnitTangent(double, double& t_x, double& t_y) const {t_x =0.0; t_y = 0.0;};

	virtual double  getDistanceToBoundary(double /* x */, double /* y */) const {return 0.0;};

//
//  Creation
//
    virtual GeometricEntity*  newDuplicateEntity() const {GeometricEntity* A = nullptr; return A;};
//
//  Equality/Inequality
//
    virtual bool  operator == (const GeometricEntity&)const {return 0;};
    virtual bool  operator != (const GeometricEntity&)const {return 0;};
    virtual void getConstructorData(std::vector<double>&, std::vector<long>&, std::string&) const {};
    virtual bool compareConstructorData(std::vector<double>&, std::vector<long>&, std::string&)const {return 0;};
//
//  Error Handling Functions
//
    static void  parametricCoordinateError()
    {
    throw std::runtime_error("\n GeometricEntity error : getParametricCoordinate arguments not on entity  \n");
    };

    static std::string getGeoTypeString(const SCC::GeoType& geoType)
    {
    std::string geoTypeString;

    switch (geoType)
    {
    case SCC::GeoType::CIRCLE       : geoTypeString = "SCC_GeoType_CIRCLE"; break;
    case SCC::GeoType::POLYGON      : geoTypeString = "SCC_GeoType_POLYGON"; break;
    case SCC::GeoType::XY_RECTANGLE : geoTypeString = "SCC_GeoType_XY_RECTANGLE"; break;
    case SCC::GeoType::COMBINED     : geoTypeString = "SCC_GeoType_COMBINED"; break;
    default:
    geoTypeString = "SCC_GeoType_NONE"; break;
    }

    return geoTypeString;
    }

};
}// namespace SCC
#endif


 
