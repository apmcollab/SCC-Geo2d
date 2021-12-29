/*
 * SCC_XYrectangleEntity.h
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


#include <cmath>
#include <iostream>
#include "SCC_GeometricEntity.h"

#ifndef SCC_XYRECTANGLE_ENTITY__
#define SCC_XYRECTANGLE_ENTITY__

#ifndef SCC_PI
#define SCC_PI 3.14159265358979323846
#endif

#define LINE_TOLERANCE 1.0e-08

namespace SCC
{
class XYrectangleEntity : public GeometricEntity
{

private :

//
//  internal representation data
//
    double  x_a;
    double  y_a;
    double  x_b;
    double  y_b;
    int     orientation;
    double  line_tolerance;

public:


//
//********************************************************************************
//                    CONSTRUCTORS
//********************************************************************************
//
//
//********************************************************************************
//
//********************************************************************************
//

XYrectangleEntity(): GeometricEntity()
{
	x_a = 0.0;
	y_a = 0.0;
	x_b = 0.0;
	y_b = 0.0;
	orientation    = 1;
	line_tolerance = LINE_TOLERANCE;
}

XYrectangleEntity( const XYrectangleEntity& A): GeometricEntity(A)
{

	x_a = A.x_a;
	y_a = A.y_a;
	x_b = A.x_b;
	y_b = A.y_b;
	orientation    = A.orientation;
	line_tolerance = A.line_tolerance;

}

XYrectangleEntity(double x_1, double y_1, double x_2, double y_2)
{
	if((x_1 >= x_2)||(y_1 >= y_2))
	{
	throw std::runtime_error("\n XYrectrangleEntity error : Invalid constructor specification  \n");
	}
	x_a = x_1;
	y_a = y_1;
	x_b = x_2;
	y_b = y_2;
    orientation    = 1;
	line_tolerance = LINE_TOLERANCE;
}
//
//********************************************************************************
//                    DESTRUCTOR
//********************************************************************************
//
~XYrectangleEntity()
{}

//
//********************************************************************************
//                    ASSIGNMENT
//********************************************************************************
//
XYrectangleEntity&  operator =( const XYrectangleEntity& A)
{
//
//  Default Assignment  : Memberwise Assignment
//
	x_a = A.x_a;
	y_a = A.y_a;
	x_b = A.x_b;
	y_b = A.y_b;
	orientation    = A.orientation;
	line_tolerance = A.line_tolerance;

    return *this;

}
//
//********************************************************************************
//                    INITIALIZATION
//********************************************************************************
//
void  initialize()
{
	x_a = 0.0;
	y_a = 0.0;
	x_b = 0.0;
	y_b = 0.0;
	orientation    = 1;
	line_tolerance = 1.0e-08;
}
void  initialize(const XYrectangleEntity& A)
{
	x_a = A.x_a;
	y_a = A.y_a;
	x_b = A.x_b;
	y_b = A.y_b;
	orientation    = A.orientation;
	line_tolerance = A.line_tolerance;
}
void  initialize(double x_1, double y_1, double x_2, double y_2)
{
	if((x_1 >= x_2)||(y_1 >= y_2))
	{
	throw std::runtime_error("\n XYrectrangleEntity error : Invalid constructor specification  \n");
    }

	x_a = x_1;
	y_a = y_1;
	x_b = x_2;
	y_b = y_2;
    orientation    = 1;
    line_tolerance = 1.0e-08;
}
//
//******************************************************************************
//                    Creation
//******************************************************************************
//
GeometricEntity*  newDuplicateEntity() const
{
	XYrectangleEntity* R = new XYrectangleEntity(*this);
	return  R;
}
//
//******************************************************************************
//                    Access
//******************************************************************************
//
    SCC::GeoType getEntityType() const {return SCC::GeoType::XY_RECTANGLE;};

    double   getLowerLeftXpoint()  const {return x_a;};
    double   getLowerLeftYpoint()  const {return y_a;};
    double   getUpperRightXpoint() const {return x_b;};
    double   getUpperRightYpoint() const {return y_b;};
    int getOrientation() const {return orientation;};

    int  getVertexCount() const {return 4;};
//
//
//******************************************************************************
//                    OUTPUT
//******************************************************************************
//
void  output(std::ostream& out_stream) const
{
   out_stream << "[BEGIN_ENTITY]\n";
   out_stream <<  GeometricEntity::getGeoTypeString(SCC::GeoType::XY_RECTANGLE) + "\n";
   out_stream << "[ENTITY_DATA]\n";
   out_stream << x_a << "   " <<  y_a  <<  '\n';
   out_stream << x_b << "   " <<  y_b  <<  '\n';
   out_stream << orientation << '\n';
   out_stream << "[END_ENTITY]\n";
}
friend std::ostream&  operator <<(std::ostream& out_stream, const SCC::XYrectangleEntity& A)
{
	A.output(out_stream);
	return(out_stream);
}

friend std::istream&  operator >>(std::istream& in_stream, SCC::XYrectangleEntity& A)
{
	 A.input(in_stream);
    return(in_stream);
}

void  input(std::istream& in_stream)
{
	std::string str;
	char delim = '\n';
    std::getline(in_stream,str,delim);
    std::getline(in_stream,str,delim);
    inputData(in_stream);
}

void  inputData(std::istream& in_stream)
{
	std::string str;
	char delim = '\n';
    std::getline(in_stream,str,delim);
    in_stream >> x_a;
    in_stream >> y_a;
    in_stream >> x_b;
    in_stream >> y_b;
    in_stream >> orientation;
    std::getline(in_stream,str,delim);

    line_tolerance = 1.0e-08;
}

//
//********************************************************************************
//                    Equality/Inequality
//********************************************************************************
//
bool  operator ==(const GeometricEntity &A) const
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
bool compareConstructorData(std::vector<double>& D, std::vector<long>& L, std::string&) const
{
   bool flag = true;
   if(D[0] != x_a)         flag = false;
   if(D[1] != y_a)         flag = false;
   if(D[2] != x_b)         flag = false;;
   if(D[3] != y_b)         flag = false;
   if(L[0] != orientation) flag = false;
   return flag;
}
void getConstructorData(std::vector<double>& D, std::vector<long>& L, std::string&)  const
{
	D.resize(4);
	L.resize(1);
	D[0] = x_a;
	D[1] = y_a;
	D[2] = x_b;
	D[3] = y_b;
	L[0] = orientation;
}
//
//********************************************************************************
//                    MEMBER_FUNCTIONS
//********************************************************************************
//

int	  interiorExteriorTest(double xTest, double yTest) const
/*
@begin_doc@
@title    @ interiorExteriorTest(double xTest, double yTest)
@Purpose  @ Returns +1 if test point is interior to the rectangle and
-1 if point is exterior to the rectangle, otherwise returns 0;
@end_doc  @ 
*/
{
	if(
		((xTest >= x_a)&&(xTest <= x_b))&&
		((yTest >= y_a)&&(yTest <= y_b))
	  )
	{
	  if(
		((xTest >= x_a)&&(xTest <= x_b))&&
		((yTest == y_a)||(yTest == y_b))
		) return 0;
	  if(
		((xTest == x_a)||(xTest == x_b))&&
		((yTest >= y_a)&&(yTest <= y_b))
		) return 0;

	  return 1;
	}


  return -1;
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
//
//  Compute the minimum distance to each of the sides
//
    double d[4];

    d[0] = computeSegmentDistance(x, y, x_a, y_a, x_b, y_a);
    d[1] = computeSegmentDistance(x, y, x_b, y_a, x_b, y_b);
    d[2] = computeSegmentDistance(x, y, x_b, y_b, x_a, y_b);
    d[3] = computeSegmentDistance(x, y, x_a, y_b, x_a, y_a);

    double dMin = d[0];
    
    int i;
    for(i = 1; i <= 3; i++)
    {
    dMin = (dMin < d[i]) ? dMin : d[i];
    }

    int intExtTest = interiorExteriorTest(x,y);

    if(intExtTest > 0)      return  dMin; // inside
    else if(intExtTest < 0) return -dMin; // outside
   
	return 0.0;
}


int	  boundaryTest(double xTest, double yTest) const
/*
@begin_doc@
@title    @ boundaryTest(double xTest, double yTest)
@Purpose  @ Returns 1 if test point is on the boundary of the
rectangle. Returns 0 otherwise.
@end_doc  @
*/
{
	if(
		((xTest >= x_a)&&(xTest <= x_b))&&
		((yTest == y_a)||(yTest == y_b))
	  ) return 1;

	if(
		((xTest == x_a)||(xTest == x_b))&&
		((yTest >= y_a)&&(yTest <= y_b))
	  ) return 1;

	return 0;
}

int	  boundaryTest(double s1, double s2, double xTest, double yTest) const
/*
@begin_doc@
@title    @ boundaryTest(double s1, double s2, double xTest, double yTest)
@Purpose  @ Returns 1 if test point is on the boundary of the
rectangle and whose coordinate is in the range [s1, s2].
Returns 0 otherwise.
@end_doc  @
*/
{
	double sx = (x_b - x_a);
	double sy = (y_b - y_a);

	double distance;

	if(
		((xTest >= x_a)&&(xTest <= x_b))&&
		((yTest == y_a)||(yTest == y_b))
	  )
	{
	if(yTest == y_a) distance = (xTest - x_a);
	if(xTest == x_b) distance = (yTest - y_a) + sx;
	if(yTest == y_b) distance = (x_b - xTest) + sx + sy;
	if(xTest == x_a) distance = (y_b - yTest) + sx + sy + sx;

	if(distance == 2.0*(sx + sy)) distance = 0;
	if((s1 <= distance)&&(distance<= s2) )return 1;
	}

	if(
		((xTest == x_a)||(xTest == x_b))&&
		((yTest >= y_a)&&(yTest <= y_b))
	  )
	{
	if(yTest == y_a) distance = (xTest - x_a);
	if(xTest == x_b) distance = (yTest - y_a) + sx;
	if(yTest == y_b) distance = (x_b - xTest) + sx + sy;
	if(xTest == x_a) distance = (y_b - yTest) + sx + sy + sx;

	if(distance == 2.0*(sx + sy)) distance = 0;
	if((s1 <= distance)&&(distance<= s2) )return 1;
	}

	return 0;
}

void	  scaleBy(double alpha)
/*
@begin_doc@
@title    @ scaleBy(double alpha)
@Purpose  @
@end_doc  @
*/
{
	double x_cent = (x_a + x_b)/2.0;
	double y_cent = (y_a + y_b)/2.0;

	x_a = x_cent + alpha*(x_a - x_cent);
	y_a = y_cent + alpha*(y_a - y_cent);
	x_b = x_cent + alpha*(x_b - x_cent);
	y_a = y_cent + alpha*(y_b - y_cent);
}

void rotateBy(double theta)
/*
@begin_doc@
@title    @ rotateBy(double theta)
@Purpose  @ 
@end_doc  @ 
*/
{
	double x_cent = (x_a + x_b)/2.0;
	double y_cent = (y_a + y_b)/2.0;
	double x_d    = x_b - x_cent;
	double y_d    = y_b - y_cent;


	if((           0 < theta)&&( theta <=  .5* SCC_PI))
	{
	x_a = x_cent - y_d;
	y_a = y_cent - x_d;
	x_b = x_cent + y_d;
	y_b = y_cent + x_d;
	}
	if((     SCC_PI  < theta)&&( theta <= 1.5* SCC_PI))
	{
	x_a = x_cent - y_d;
	y_a = y_cent - x_d;
	x_b = x_cent + y_d;
	y_b = y_cent + x_d;
	}

}

void translateBy(double x, double y)
/*
@begin_doc@
@title    @ translateBy(double x, double y)
@Purpose  @ 
@end_doc  @ 
*/
{
	x_a = x_a + x;
	y_a = y_a + y;
	x_b = x_b + x;
	y_b = y_b + y;
}

void  	  translateTo(double x, double y)
/*
@begin_doc@
@title    @ translateTo(double x, double y)
@Purpose  @ Translates the center of the rectangle to the input
point.
@end_doc  @ 
*/
{
	double x_trans = x - (x_a + x_b)/2.0;
	double y_trans = y - (y_a + y_b)/2.0;
	translateBy(x_trans,y_trans);
}

double getMinX() const
/*
@begin_doc@
@title    @ getMinX()
@Purpose  @ 
@end_doc  @ 
*/
{
	return x_a;
}

double getMaxX()const
/*
@begin_doc@
@title    @ getMaxX()
@Purpose  @ 
@end_doc  @ 
*/
{
	return x_b;
}

double getMinY() const
/*
@begin_doc@
@title    @ getMinY()
@Purpose  @ 
@end_doc  @ 
*/
{
	return y_a;
}

double getMaxY() const
/*
@begin_doc@
@title    @ getMaxY()
@Purpose  @ 
@end_doc  @ 
*/
{
	return y_b;
}

double getTotalArcLength() const
/*
@begin_doc@
@title    @ getTotalArcLength()
@Purpose  @ 
@end_doc  @ 
*/
{
	return (2.0*(x_b - x_a) + 2.0*(y_b - y_a));
}

int getParametricCoordinate(double& s, double x, double y) const
/*
@begin_doc@
@title    @ getParametricCoordinate(double& s, double x, double y)
@Purpose  @ If the point (x,y) is on the curve then the routine
returns the value 1 and the argument s is set to the parametric
coordinate. The routine returns a 0 if the point is not on the
curve.
@end_doc  @
*/
{
	if(boundaryTest(x,y) == 0) return 0;

	double s1 = (x_b - x_a);
	double s2 = (y_b - y_a);

	double distance;

	if(y == y_a) distance = (x - x_a);
	if(x == x_b) distance = (y - y_a) + s1;
	if(y == y_b) distance = (x_b - x) + s1 + s2;
	if(x == x_a) distance = (y_b - y) + s1 + s2 + s1;

	if(distance == 2.0*(s1 + s2)) distance = 0;
	s = distance;
	return 1;
}

int getParametricCoordinate(double& s, double s1, double
s2, double x, double y) const
/*
@begin_doc@
@title    @ getParametricCoordinate(double& s, double x, double y)
@Purpose  @ If the point (x,y) is on the curve with parametric
coordinate in the range [s1, s2] then the routine
returns the value 1 and the argument s is set to the parametric
coordinate. The routine returns a 0 if the point is not on the
curve.
@end_doc  @
*/
{
	if(boundaryTest(s1,s2,x,y) == 0) return 0;

	double sx = (x_b - x_a);
	double sy = (y_b - y_a);

	double distance;

	if(y == y_a) distance = (x - x_a);
	if(x == x_b) distance = (y - y_a) + sx;
	if(y == y_b) distance = (x_b - x) + sx + sy;
	if(x == x_a) distance = (y_b - y) + sx + sy + sx;

	if(distance == 2.0*(sx + sy)) distance = 0;
	s = distance;
	return 1;
}
double	  getXcoordinate(double s)const
/*
@begin_doc@
@title    @ getXcoordinate(double s)
@Purpose  @ 
@end_doc  @ 
*/
{
	double d1 = (x_b - x_a);
	double d2 = d1 + (y_b - y_a);
	double d3 = d2 + (x_b - x_a);
	double d4 = d3 + (y_b - y_a);

	double returnValue;

    if((s > d4)||(s < 0))
	{
	throw std::runtime_error("\n XYrectrangleEntity error : getXCoordinate() argument out of range  \n");
	}

	if(s <= d1)              returnValue =  x_a + s;
	if((d1 < s)&&(s <= d2))  returnValue = x_b;
	if((d2 < s)&&(s <= d3))  returnValue = x_b - (s - d2);
	if((d3 < s)&&(s <= d4))  returnValue = x_a;

	return returnValue;
}

double	  getYcoordinate(double s) const
/*
@begin_doc@
@title    @ getYcoordinate(double s)
@Purpose  @ 
@end_doc  @ 
*/
{
	double d1 = (x_b - x_a);
	double d2 = d1 + (y_b - y_a);
	double d3 = d2 + (x_b - x_a);
	double d4 = d3 + (y_b - y_a);

	double returnValue;

	if((s > d4)||(s < 0))
	{
	throw std::runtime_error("\n XYrectrangleEntity error : getXCoordinate() argument out of range  \n");
	}

	if(s <= d1)              returnValue = y_a;
	if((d1 < s)&&(s <= d2))  returnValue = y_a + (s - d1);
	if((d2 < s)&&(s <= d3))  returnValue = y_b;
	if((d3 < s)&&(s <= d4))  returnValue = y_b - (s - d3);

    return returnValue;
}

double  getVertexCoordinate(long vertexIndex) const
{
	double d[4];
	long vertexCount = getVertexCount();

	if((vertexIndex <= 0)||(vertexIndex > vertexCount))
	{throw std::runtime_error("\n XYrectrangleEntity error : vertex index out of range \n");}

	d[0] = 0.0;
	d[1] = (x_b - x_a);
	d[2] = d[1] + (y_b - y_a);
	d[3] = d[2] + (x_b - x_a);

    return d[vertexIndex - 1];
}

int  getInteriorPoint(double& x, double& y)  const
/*
@begin_doc@
@title    @ getInteriorPoint(double& x, double& y)
@Purpose  @ 
@end_doc  @ 
*/
{
	x = (x_a + x_b)/2.0;
    y = (y_a + y_b)/2.0;
	return 1;
}

int getSegmentIntersection(double& intersectPoint,
double x_1, double y_1, double x_2, double y_2) const
/*
@begin_doc@
@title    @ getSegmentIntersection(double& intersectPoint, double x_1,
			  double y_1, double x_2, double y_2)
@Purpose  @ 
@end_doc  @
*/
{
	double U_a, V_a, U_b, V_b;

	double s           = 2.0;
	double t           = 2.0;
	double s_intersect = 2.0;
	double t_intersect;
	int intersectFlag  = 0;

	U_a = x_a;
	U_b = x_b;
	V_a = y_a;
	V_b = y_a;
	if(intersectSegments(s, x_1, y_1, x_2, y_2, t, U_a, V_a, U_b, V_b) == 1)
	{s_intersect = s; t_intersect = t*(x_b - x_a); intersectFlag = 1;}

	U_a = x_b; U_b = x_b;
	V_a = y_a; V_b = y_b;
	if(intersectSegments(s, x_1, y_1, x_2, y_2, t, U_a, V_a, U_b, V_b) == 1)
	{ if(s <= s_intersect)
	  {
	  s_intersect = s;
	  t_intersect = t*(y_b - y_a) + (x_b - x_a);
	  }
	  intersectFlag = 1;
	}

	U_a = x_b; U_b = x_a;
	V_a = y_b; V_b = y_b;
	if(intersectSegments(s, x_1, y_1, x_2, y_2, t, U_a, V_a, U_b, V_b) == 1)
	{ if(s <= s_intersect)
	  {
	  s_intersect = s;
	  t_intersect = t*(x_b - x_a) + (y_b - y_a) + (x_b - x_a);
	  }
	  intersectFlag = 1;
	}

	U_a = x_a; U_b = x_a;
	V_a = y_b; V_b = y_a;
	if(intersectSegments(s, x_1, y_1, x_2, y_2, t, U_a, V_a, U_b, V_b) == 1)
	{ if(s <= s_intersect)
	  {
	  s_intersect = s;
	  t_intersect = t*(y_b - y_a) + (y_b - y_a) + 2.0*(x_b - x_a);
	  }
	  intersectFlag = 1;
	}

	if(intersectFlag == 1)
	{
	intersectPoint = t_intersect;
	}
	return intersectFlag;
}

int  getSegmentIntersection(double& intersectPoint,
double s_1, double s_2, double x_1, double y_1, double x_2, double y_2)  const
/*
@begin_doc@
@title    @ getSegmentIntersection(double& intersectPoint, double s_1,
			  double s_2, double x_1, double y_1, double x_2, double y_2)
@Purpose  @ 
@end_doc  @ 
*/
{
	if(getSegmentIntersection(intersectPoint,x_1, y_1, x_2, y_2) == 1)
	{
	if((s_1 <= intersectPoint)&&(intersectPoint <= s_2)) return 1;
	}
	return 0;

}

void  getUnitNormal(double s, double& n_x, double& n_y) const
/*
@begin_doc@
@title    @ getUnitNormal(double s, double& n_x, double& n_y)
@Purpose  @ 
@end_doc  @ 
*/
{
	double d1 = (x_b - x_a);
	double d2 = d1 + (y_b - y_a);
	double d3 = d2 + (x_b - x_a);
	double d4 = d3 + (y_b - y_a);

	double oneOverSqrt2 = 1.0/sqrt(2.0);

	if((s > d4)||(s < 0))
	{
	throw std::runtime_error("\n XYrectrangleEntity error : getUnitNormal argument out of range \n");
	}



	if(s == 0.0)             {n_x =-oneOverSqrt2; n_y =-oneOverSqrt2;}
	else if((0.0 < s)&&(s <  d1)) {n_x = 0.0;          n_y =-1.0;}
	else if(s == d1)              {n_x = oneOverSqrt2; n_y =-oneOverSqrt2;}
	else if((d1  < s)&&(s <  d2)) {n_x = 1.0;          n_y = 0.0;}
	else if(s == d2)              {n_x = oneOverSqrt2; n_y =  oneOverSqrt2;}
	else if((d2  < s)&&(s <  d3)) {n_x = 0.0;          n_y = 1.0;}
	else if(s == d3)              {n_x =-oneOverSqrt2; n_y =  oneOverSqrt2;}
	else if((d3  < s)&&(s <  d4)) {n_x =-1.0;          n_y = 0.0;}
	else if(s == d4)              {n_x =-oneOverSqrt2; n_y =-oneOverSqrt2;}
}

void  getUnitTangent(double s, double& t_x, double& t_y) const
/*
@begin_doc@
@title    @ getUnitTangent(double s, double& t_x, double& t_y)
@Purpose  @ 
@end_doc  @ 
*/
{
	double d1 = (x_b - x_a);
	double d2 = d1 + (y_b - y_a);
	double d3 = d2 + (x_b - x_a);
	double d4 = d3 + (y_b - y_a);

	double oneOverSqrt2 = 1.0/sqrt(2.0);

	if((s > d4)||(s < 0))
	{throw std::runtime_error("\n XYrectrangleEntity error : getUnitTangent argument out of range \n");}


	if(s == 0.0)                  {t_x = oneOverSqrt2; t_y =-oneOverSqrt2;}
	else if((0.0 < s)&&(s <  d1)) {t_x = 1.0;          t_y = 0.0;}
	else if(s == d1)              {t_x = oneOverSqrt2; t_y = oneOverSqrt2;}
	else if((d1  < s)&&(s <  d2)) {t_x = 0.0;          t_y = 1.0;}
	else if(s == d2)              {t_x =-oneOverSqrt2; t_y = oneOverSqrt2;}
	else if((d2  < s)&&(s <  d3)) {t_x = -1.0;         t_y = 0.0;}
	else if(s == d3)              {t_x =-oneOverSqrt2; t_y =-oneOverSqrt2;}
	else if((d3  < s)&&(s <  d4)) {t_x = 0.0;          t_y =-1.0;}
	else if(s == d4)              {t_x = oneOverSqrt2; t_y =-oneOverSqrt2;}
}

int  intersectSegments(double& s, double X_a, double Y_a,
double X_b, double Y_b, double& t, double U_a, double V_a, double U_b, double V_b) const
/*
@begin_doc@
@title    @ intersectSegments(double& s, double x_a, double y_a,
double x_b, double y_b, double& t, double u_a, double v_a, double u_b, double v_b)
@Purpose  @
This routine returns the intersection point of an
"interrogator" segment (X_a, Y_a)<->(X_b,Y_b) with a "target"
segment (U_a, V_a)<->(U_b,V_b). The interrogator segment is
parameterized by s, and the target segment is parameterized by t.

The notion of interrogator and target segments is introduced
to handle the case when the segments overlap one another. In
such a case, the result is always one endpoint of the target
segment; the intersection point being defined as the
first point of the target which is hit by the interrogator.
In the case that the interrogator contains the target segment,
the endpoint closest to (X_a, Y_a), the origin of the
interrogator, is used.

This routine sets the values of s and t to be the parametric
coordinate so that the intersection point on each line is given
by
	(X_a,Y_a) + s*(X_b - X_a, Y_b - Y_a)  (interrogator)

	(U_a,V_a) + t*(U_b - U_a, V_b - V_a)  (target)

Things to do : Add line tolerance correctness

@end_doc  @ 
*/
{
	double s_0; double s_1;
	double x_ba = (X_b - X_a);
	double y_ba = (Y_b - Y_a);
	double u_ba = (U_b - U_a);
	double v_ba = (V_b - V_a);
	double determinant = (u_ba*y_ba)-(x_ba*v_ba);
	double area;

	if(std::abs(determinant) > 1.0e-09)
	{
	s = (-v_ba*(U_a - X_a) + u_ba*(V_a - Y_a))/determinant;
	t = (-y_ba*(U_a - X_a) + x_ba*(V_a - Y_a))/determinant;
	if(
	   ((0.0 <= s)&&(s <= 1.0)) && ((0.0 <= t)&&(t <= 1.0))
	   ) return 1;
	else
	{s = 0.0; t = 0.0; return 0;}

	}
//
//  Parallel segments - do work if co-linear
//
	area = (U_b*Y_a - X_a*V_b) - (U_a*Y_a - X_a*V_a) +
		   (U_a*V_b - U_b*V_a);

	if(std::abs(area) <= 1.0e-09)
	{

//
//  Compute interrogator values for target endpoints
//
	if(std::abs(x_ba) > std::abs(y_ba))
	{
	 s_0 =  (U_a - X_a)/x_ba;
	 s_1 = ((U_a - X_a) + u_ba)/x_ba;
	}
	else
	{
	 s_0 = ( V_a - Y_a)/y_ba;
	 s_1 = ((V_a - Y_a) + v_ba)/y_ba;
	}

	if(((0.0 <= s_0)&&(s_0 <= 1.0))&&((s_1 < 0)||(s_1 > 1.0)))
	  {s = s_0; t = 0.0; return 1;}
	else if(((0.0 <= s_1)&&(s_1 <= 1.0))&&((s_0 < 0)||(s_0 > 1.0)))
	  {s = s_1; t = 1.0; return 1;}
	else if(((0.0 <= s_0)&&(s_0 <= 1.0))&&((0.0 <= s_1)&&(s_1 <= 1.0)))
	{
	  if(s_0 <= s_1)
	  {s = s_0; t = 0.0; return 1;}
	  else
	  {s = s_1; t = 1.0; return 1;}
	}
	}


	return 0;
}

/**
This routine computes the distance to a segment from (xA,yA) to (xB,yB). 
It does this by minimizing the square of the distance from the 
test point to the line through the two points as a function of the 
parametric coordinate. 

*/
double computeSegmentDistance(double x, double y,
       double xA, double yA, double xB, double yB) const
{
//
//	Find s* = value at which the distance between x and 
//  line though segment endpoins is minimum.
//
    double dxAB     = xB - xA;
    double dyAB     = yB - yA;

    double dxAX     = x  - xA;
    double dyAY     = y  - yA;

    double xStar; double yStar;
    double sStar;

    double dA; double dB;

    double dABnorm2 = dxAB*dxAB + dyAB*dyAB;
//
//  If endpoints of segment coincide to line tolerance, then return 
//  distance from (x,y) to (xA,yA).
//
    if(dABnorm2 < line_tolerance*line_tolerance) return sqrt(dxAX*dxAX + dyAY*dyAY);

    sStar  = (dxAX*dxAB + dyAY*dyAB)/dABnorm2;

    if((sStar > 0)&&(sStar < 1.0))
    {
        xStar = xA + sStar*dxAB;
        yStar = yA + sStar*dyAB;
        return sqrt((x-xStar)*(x-xStar) + (y-yStar)*(y-yStar));
    }
    else
    {
        dA = sqrt(dxAX*dxAX + dyAY*dyAY);
        dB = sqrt((x-xB)*(x-xB) + (y-yB)*(y-yB));
        if(dA < dB) return dA;
        return dB;
    }

    return 0.0;
}

};
} //namespace SCC

#undef SCC_PI
#undef LINE_TOLERANCE
#endif

 
