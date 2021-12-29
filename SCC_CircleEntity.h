/*
 * SCC_CircleEntity.h
 *
 *  Created on: Dec. 28, 2021
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

#include <cmath>
#include <iostream>
#include "SCC_GeometricEntity.h"

#ifndef SCC_CIRCLE_ENTITY_
#define SCC_CIRCLE_ENTITY_

#ifndef SCC_PI
#define SCC_PI 3.14159265358979323846
#endif

#define LINE_TOLERANCE 1.0e-08

namespace SCC
{
class CircleEntity : public GeometricEntity
{

private :

//
//  internal representation data
//
    double        center_x;
    double        center_y;
    double          radius;
    double     start_theta;
    int        orientation;
    double  line_tolerance;

public:

CircleEntity(): GeometricEntity()
{
	this->center_x       = 0.0;
	this->center_y       = 0.0;
	this->radius         = 0.0;
	this->start_theta    = 0.0;
	this->orientation    = 1;
	this->line_tolerance = LINE_TOLERANCE;
}

CircleEntity(const CircleEntity& A): GeometricEntity(A)
{
	this->center_x       = A.center_x;
	this->center_y       = A.center_y;
	this->radius         = A.radius;
	this->start_theta    = A.start_theta;
	this->orientation    = A.orientation;
	this->line_tolerance = A.line_tolerance;
}

CircleEntity(double x_center, double y_center, double r): GeometricEntity()
{
	this->center_x       = x_center;
	this->center_y       = y_center;
	this->radius         = r;
    this->start_theta    = 0.0;
	this->orientation    = 1;
	this->line_tolerance = LINE_TOLERANCE;
}

//
//********************************************************************************
//                    DESTRUCTOR
//********************************************************************************
//
~CircleEntity(){}

//
//********************************************************************************
//                    ASSIGNMENT
//********************************************************************************
//
CircleEntity& operator =( const CircleEntity& A)
{
    this->center_x       = A.center_x;
    this->center_y       = A.center_y;
    this->radius         = A.radius;
    this->start_theta    = A.start_theta;
	this->orientation    = A.orientation;
	this->line_tolerance = A.line_tolerance;

	return *this;
}

//********************************************************************************
//                    INITIALIZATION
//********************************************************************************

void  initialize()
{
	center_x = 0.0;
	center_y = 0.0;
	radius = 0.0;
	start_theta = 0.0;
	orientation = 1;
	line_tolerance = 1.0e-8;

}
void initialize(const CircleEntity& A)
{
	center_x        = A.center_x;
	center_y        = A.center_y;
	radius          = A.radius;
	start_theta     = A.start_theta;
	orientation     = A.orientation;
	line_tolerance  = A.line_tolerance;
}

void initialize(double x_center, double y_center, double r)
{
	center_x = x_center;
	center_y = y_center;
	radius = r;
	start_theta = 0.0;
	orientation = 1;
	line_tolerance = 1.0e-8;
}

//
//********************************************************************************
//                    Creation
//********************************************************************************
//
GeometricEntity*  newDuplicateEntity() const
{
	CircleEntity* R = new CircleEntity(*this);
	return  R;
}
//
//********************************************************************************
//                    Equality/Inequality
//********************************************************************************
//
bool operator ==(const GeometricEntity &A) const
{
	if(A.getEntityType() != getEntityType()) return false;
	std::vector<double> D;
	std::vector<long>   L;
	std::string         C;
	A.getConstructorData(D,L,C);
	bool flag = compareConstructorData(D,L,C);
	return flag;
}
bool operator !=(const GeometricEntity &A) const
{
    return !(*this == A);
}

bool  compareConstructorData(std::vector<double>& D, std::vector<long>& L, std::string&) const
{
   bool flag = true;;
   if(D[0] != center_x)    flag = false;
   if(D[1] != center_y)    flag = false;
   if(D[2] != radius)      flag = false;
   if(D[3] != start_theta) flag = false;
   if(L[0] != orientation) flag = false;
   return flag;
}
void getConstructorData(std::vector<double>& D, std::vector<long>& L, std::string&)  const
{
    D.resize(4);
    L.resize(1);

	D[0] = center_x;
	D[1] = center_y;
	D[2] = radius;
	D[3] = start_theta;
	L[0] = orientation;
}

//
//******************************************************************************
//                    MEMBER_FUNCTIONS
//******************************************************************************
//
//
//
//******************************************************************************
//                       ACCESS
//******************************************************************************
//
//
    SCC::GeoType getEntityType() const {return SCC::GeoType::CIRCLE;};

	double  getRadius()        const {return radius;}
    double  getXcenter()       const {return center_x;}
    double  getYcenter()       const {return center_y;}
    double  getStartingTheta() const {return start_theta;}
    int     getOrientation()   const {return orientation;}


    void setRadius(double val)         {radius = val;}
    void setXcenter(double val)        {center_x = val;}
    void setYcenter(double val)        {center_y = val;}

    void setCenter(double xPos, double yPos)
    {center_x = xPos; center_y = yPos;}

    void setStartingTheta(double val)  {start_theta = val;}
    void setOrientation(int val)       {orientation = val;}

//
//******************************************************************************
//                    INPUT/OUTPUT
//******************************************************************************
//
void output(std::ostream& out_stream) const
{
   out_stream << "[BEGIN_ENTITY]\n";
   out_stream <<  GeometricEntity::getGeoTypeString(SCC::GeoType::CIRCLE) << "\n";
   out_stream << "[ENTITY_DATA]\n";
   out_stream << center_x << "  " << center_y <<  '\n';
   out_stream << radius   <<  '\n';
   out_stream << orientation << '\n';
   out_stream << "[END_ENTITY]";
}
friend std::ostream&  operator <<(std::ostream& out_stream, const SCC::CircleEntity& A)
{
	A.output(out_stream);
    return(out_stream);
}

friend std::istream&  operator >>(std::istream& in_stream, SCC::CircleEntity& A)
{
	A.input(in_stream);
    return(in_stream);
}

void input(std::istream& in_stream)
{
	std::string str;
	char delim = '\n';
    std::getline(in_stream,str,delim);
    std::getline(in_stream,str,delim);
    inputData(in_stream);
}

void inputData(std::istream& in_stream)
{
	std::string str;
	char delim = '\n';
    std::getline(in_stream,str,delim);
    in_stream >> center_x;
    in_stream >> center_y;
    in_stream >> radius;
    in_stream >> orientation;
    std::getline(in_stream,str,delim);
    line_tolerance = 1.0e-8;
}


//
//********************************************************************************
//                    MEMBER_FUNCTIONS
//********************************************************************************
//

int interiorExteriorTest(double xTest, double yTest) const
/*
@begin_doc@
@title    @ interiorExteriorTest()
@Purpose  @ Returns +1 if test point is interior to the circle and
-1 if point is exterior to the circle, otherwise returns 0;
@end_doc  @
*/
{
	double rTest;
	double relTolerance;

	relTolerance = line_tolerance*radius;

	rTest   = std::sqrt((xTest - center_x)*(xTest - center_x) +
			  (yTest - center_y)*(yTest - center_y));

	if(rTest < radius - relTolerance) return   1;
	if(rTest > radius + relTolerance) return  -1;
	return 0;
}
/*
@begin_doc@
@title    @ getDistanceToBoundary()
@Purpose  @ Returns signed distance to the boundary.(+) for interior distances
(-1) for exterior distances;
@end_doc  @
*/
double getDistanceToBoundary(double x, double y) const
{
	double rTest;

	rTest   = std::sqrt((x - center_x)*(x - center_x) +
			  (y - center_y)*(y - center_y));
	return radius - rTest;
}

int boundaryTest(double xTest, double yTest) const
/*
@begin_doc@
@title    @ boundaryTest(double xTest, double yTest)
@Purpose  @ Returns 1 if test point is on the boundary of the
circle. Returns 0 otherwise.
@end_doc  @ 
*/
{
	double rTest;
	double relTolerance;

	relTolerance = line_tolerance*radius;

	rTest   = std::sqrt((xTest - center_x)*(xTest - center_x) +
			  (yTest - center_y)*(yTest - center_y));

	if((rTest > radius - relTolerance)&&
	   (rTest < radius + relTolerance)) return 1;

	return 0;
}
int  boundaryTest(double s1, double s2, double xTest,
double yTest) const
/*
@begin_doc@
@title    @ boundaryTest(double xTest, double yTest)
@Purpose  @ Returns 1 if the test point is on the boundary of the
circle within the parametric coordinate range [s1, s2].
Returns 0 otherwise.
@end_doc  @ 
*/
{
	double rTest;
	double relTolerance;

	relTolerance = line_tolerance*radius;

	rTest   = sqrt((xTest - center_x)*(xTest - center_x) +
			  (yTest - center_y)*(yTest - center_y));

	if((rTest > radius - relTolerance)&&
	   (rTest < radius + relTolerance))
	{
	double s;
	double theta = atan2(yTest - center_y,xTest - center_x);
	if(theta < 0.0) theta += 2.0*SCC_PI;

	if(theta < start_theta)
	{s = radius*(theta + start_theta);}
	else
	{s = radius*(theta - start_theta);}

	if((s1 <= s)&&(s <= s2)) return 1;
	}
	return 0;
}

void scaleBy(double alpha)
/*
@begin_doc@
@title    @ scaleBy(double alpha)
@Purpose  @ 
@end_doc  @ 
*/
{
	radius = alpha*radius;
}

void  rotateBy(double theta)
/*
@begin_doc@
@title    @ rotateBy(double theta)
@Purpose  @ 
@end_doc  @ 
*/
{
	start_theta = start_theta + theta;
}

void  translateBy(double x, double y)
/*
@begin_doc@
@title    @ translateBy(double x, double y)
@Purpose  @ 
@end_doc  @ 
*/
{
	center_x = center_x + x;
	center_y = center_y + y;
}

void translateTo(double x, double y)
/*
@begin_doc@
@title    @ translateTo(double #include "camgeoimpexp.h"x, double y)
@Purpose  @ 
@end_doc  @ 
*/
{
	center_x = x;
	center_y = y;
}

double 	getMinX() const
/*
@begin_doc@
@title    @ getMinX()
@Purpose  @ 
@end_doc  @ 
*/
{
	return center_x - radius;
}

double getMaxX() const
/*
@begin_doc@
@title    @ getMaxX()
@Purpose  @ 
@end_doc  @ 
*/
{
	return center_x + radius;
}

double getMinY() const
/*
@begin_doc@
@title    @ getMinY()
@Purpose  @ 
@end_doc  @ 
*/
{
	return center_y - radius;
}

double 	getMaxY() const
/*
@begin_doc@
@title    @ getMaxY()
@Purpose  @ 
@end_doc  @ 
*/
{
	return center_y + radius;
}

double	getTotalArcLength() const
/*
@begin_doc@
@title    @ getTotalArcLength()
@Purpose  @ 
@end_doc  @ 
*/
{
	return 2.0*SCC_PI*radius;
}

int	getParametricCoordinate(double& s, double x, double y) const
/*
@begin_doc@
@title    @ getParametricCoordinate(double& s, double x, double y)
@Purpose  @ This routine gets the parametric coordinate of the
point (x,y). If the point is not on the curve the routine returns
0, otherwise the routine returns a value of 1, and the variable s
is assigned to the parametric coordinate value.
@end_doc  @
*/
{
	if(boundaryTest(x,y) == 0)
	{return 0;}

	double theta = std::atan2(y - center_y,x - center_x);
	if(theta < 0.0) theta += 2.0*SCC_PI;

	if(theta < start_theta)
	{s = radius*(theta + start_theta);}
	else
	{s = radius*(theta - start_theta);}

	return 1;
}

int	getParametricCoordinate(double& s, double s1, double s2, double x, double y) const
/*
@begin_doc@
@title    @ getParametricCoordinate(double& s, double s1, double s2, double x, double y)
@Purpose  @ This routine gets the parametric coordinate of the
point (x,y) when the cooresponding parametric coordinate of
the point is between s1 and s2.
If the point is not on the curve within the parametric coordinate
range [s1,s2] the routine returns
0, otherwise the routine returns a value of 1, and the variable s
is assigned to the parametric coordinate value.
@end_doc  @
*/
{
	if(boundaryTest(s1,s2,x,y) == 0)
	{return 0;}

	double theta = std::atan2(y - center_y,x - center_x);
	if(theta < 0.0) theta += 2.0*SCC_PI;

	if(theta < start_theta)
	{s = radius*(theta + start_theta);}
	else
	{s = radius*(theta - start_theta);}

	return 1;
}

double getXcoordinate(double s) const
/*
@begin_doc@
@title    @ getXcoordinate(double s)
@Purpose  @ 
@end_doc  @ 
*/
{
	if((s > radius*2.0*SCC_PI)||(s < 0))
	{
	throw std::runtime_error("\n CircleEntity error : getXCoordinate() argument out of range \n");
    }
	return center_x + radius*(cos((s/radius) + start_theta));
}

double getYcoordinate(double s) const
/*
@begin_doc@
@title    @ getYcoordinate(double s)
@Purpose  @ 
@end_doc  @ 
*/
{
	if((s > radius*2.0*SCC_PI)||(s < 0))
	throw std::runtime_error("\n CircleEntity error : getYCoordinate() argument out of range \n");
	return center_y + radius*(sin((s/radius) + start_theta));
}

int getInteriorPoint(double& x, double& y) const
/*
@begin_doc@
@title    @ getInteriorPoint(double& x, double& y)
@Purpose  @ 
@end_doc  @ 
*/
{
	x = center_x;
    y = center_y;
	return 1;
}

int getSegmentIntersection(double& intersectPoint, double x_a,
			  double y_a, double x_b, double y_b) const
/*
@begin_doc@
@title    @ getSegmentIntersection(double& intersectPoint, double x_1,
			  double y_1, double x_2, double y_2)
@Purpose  @ 
@end_doc  @ 
*/
{
	double x_ba = x_b  - x_a;
	double y_ba = y_b  - y_a;
	double x_ac = x_a  - center_x;
	double y_ac = y_a  - center_y;

	double a = x_ba*x_ba + y_ba*y_ba;
	double b = 2.0*(x_ac*x_ba + y_ac*y_ba);
	double c = x_ac*x_ac + y_ac*y_ac - radius*radius;

	double s_1; double s_2;
	double x_p; double y_p;
	int returnFlag = 0;

	if(getQuadraticRoots(a, b, c, s_1, s_2) == 0) return 0;

	if((s_1 >=0)&&(s_1 <= 1.0))
	{
	x_p = x_a + s_1*(x_b - x_a);
	y_p = y_a + s_1*(y_b - y_a);
	if(getParametricCoordinate(intersectPoint, x_p, y_p) == 0)
	GeometricEntity::parametricCoordinateError();
	returnFlag = 1;
	}

	if((s_2 >=0)&&(s_2 <= 1.0))
	{
	if(returnFlag == 1)
	{
	    if(s_2 < s_1)
		{
		x_p = x_a + s_2*(x_b - x_a);
		y_p = y_a + s_2*(y_b - y_a);
		if(getParametricCoordinate(intersectPoint, x_p, y_p) == 0)
		GeometricEntity::parametricCoordinateError();
		}
	}
	else
	{
	  x_p = x_a + s_2*(x_b - x_a);
	  y_p = y_a + s_2*(y_b - y_a);
	  if(getParametricCoordinate(intersectPoint, x_p, y_p) == 0)
	  GeometricEntity::parametricCoordinateError();
	  returnFlag = 1;
	}
	}

	return returnFlag;
}

int getSegmentIntersection(double& intersectPoint, double s_1,
			  double s_2, double x_1, double y_1, double x_2, double y_2) const
/*
@begin_doc@
@title    @ getSegmentIntersection(double& intersectPoint, double s_1,
		      double s_2, double x_1, double y_1, double x_2, double y_2)
@Purpose  @ 
@end_doc  @ 
*/
{
	if(getSegmentIntersection(intersectPoint, x_1,y_1, x_2, y_2) == 1)
	{
    if((intersectPoint >= s_1)&&(intersectPoint <= s_2)) return 1;
	}
	return 0;
}

void getUnitNormal(double s, double& n_x, double& n_y) const
/*
@begin_doc@
@title    @ getNormal(double s, double& n_x, double& n_y)
@Purpose  @ 
@end_doc  @ 
*/
{
	if((s > radius*2.0*SCC_PI)||(s < 0))
	{std::cerr << "Error : getUnitNormal parametric argument out of range "
		  <<std:: endl; exit(1);}

	n_x = cos((s/radius) + start_theta);
	n_y = sin((s/radius) + start_theta);
}

void getUnitTangent(double s, double& t_x, double& t_y) const
/*
@begin_doc@
@title    @ getTangent(double s, double& t_x, double& t_y)
@Purpose  @ 
@end_doc  @ 
*/
{
	if((s > radius*2.0*SCC_PI)||(s < 0))
	{std::cerr << "Error : getUnitTangent parametric argument out of range "
		  << std::endl; exit(1);}

	t_x = -sin((s/radius) + start_theta);
	t_y =  cos((s/radius) + start_theta);
}
//
//********************************************************************************
//                    INTERNAL MEMBER_FUNCTIONS
//********************************************************************************
//
int  getQuadraticRoots(double a, double b, double c, double& r_1, double& r_2) const
/*
@begin_doc@
@title    @ getQuadraticRoots(double a, double b, double c, double& r_1, double r_2)
@Purpose  @ 
@end_doc  @ 
*/
//
// Added b== 0 and c==0 code 1/15/98 CRA
//
{
    double discriminant = (b*b - 4.0*a*c);
    if( (discriminant < 0)||
        ((a == 0.0)&&(b == 0))
      ) return 0;

	 if((b == 0)&&(c==0))
	 {
	  r_1 = 0.0;
	  r_2 = 0.0;
	  return 1;
	 }
	 
    if(a == 0.0)
    {
	  r_1 = -c/b;
	  r_2 = -c/b;
	  return 1;
    }

	 double sqrtD;
	 if(discriminant == 0.0)
	 {sqrtD = 0.0;}
	 else
	 {sqrtD  = sqrt(discriminant);}
	 
    if(b >= 0)
    {
        r_1 = (-2.0*c)/(b + sqrtD);
        r_2 = ((-b) - sqrtD)/(2.0*a);
    }
    else
    {
        r_1 = ((-b) + sqrtD)/(2.0*a);
        r_2 = (-2.0*c)/(b - sqrtD);
    }
    return 1;
}

};
} //namespace SCC

#undef SCC_PI
#undef LINE_TOLERANCE
#endif

//
//********************************************************************************
//                     CPP File End
//********************************************************************************
//

 
