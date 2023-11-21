/*
 * SCC_EntityGraphics.h
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
#include <vector>

#include "CppGraphics/CAMgraphics.h"


#include "SCC_XYrectangleEntity.h"
#include "SCC_PolygonEntity.h"
#include "SCC_CircleEntity.h"
#include "SCC_CombinedEntity.h"
/*
#############################################################################
#
# Copyright 1997-2019 Chris Anderson
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
#ifndef ENTITY_GRAPHICS_
#define ENTITY_GRAPHICS_

#ifndef SCC_PI
#define SCC_PI 3.14159265358979323846
#endif


namespace SCC
{
class EntityGraphics : public CAMgraphics
{
public :

EntityGraphics(): CAMgraphics() {};

virtual ~EntityGraphics(){};

void draw(const XYrectangleEntity &R)
{
	double  x[5];
	double  y[5];
	double x_a    =  R.getLowerLeftXpoint();
	double y_a    =  R.getLowerLeftYpoint();
	double x_b    =  R.getUpperRightXpoint();
	double y_b    =  R.getUpperRightYpoint();
	x[0] = x_a; y[0] = y_a;
	x[1] = x_b; y[1] = y_a;
	x[2] = x_b; y[2] = y_b;
	x[3] = x_a; y[3] = y_b;
	x[4] = x_a; y[4] = y_a;
	plot(x,y,5);
}

void fill(const XYrectangleEntity &R)
{
	double  x[5];
	double  y[5];
	double x_a    =  R.getLowerLeftXpoint();
	double y_a    =  R.getLowerLeftYpoint();
	double x_b    =  R.getUpperRightXpoint();
	double y_b    =  R.getUpperRightYpoint();
	x[0] = x_a; y[0] = y_a;
	x[1] = x_b; y[1] = y_a;
	x[2] = x_b; y[2] = y_b;
	x[3] = x_a; y[3] = y_b;
	x[4] = x_a; y[4] = y_a;
	region(x,y,5);
}


void draw(const PolygonEntity &P)
{
    std::vector<double> x = P.getXvertices();
    std::vector<double> y = P.getYvertices();
    long n = P.getSideCount();
    if(n != 0)
    {
    	plot(&x[0],&y[0],n+1);
    }
}

void fill(const PolygonEntity &P)
{
    std::vector<double> x = P.getXvertices();
    std::vector<double> y = P.getYvertices();
    long n = P.getSideCount();
    if(n != 0)
    {
    	region(&x[0],&y[0],n+1);
    }
}

void draw(const CircleEntity &C)
{
	double center_x = C.getXcenter();
	double center_y = C.getYcenter();
	double radius   = C.getRadius();
	long i;
	long nTheta = 50;
	double dTheta = (2.0*SCC_PI)/double(nTheta);
	double* x = new double[nTheta +1];
	double* y = new double[nTheta +1];
	for(i = 1; i <= 51; i++)
	{
	x[i-1] =  center_x + radius*cos(dTheta*double((i-1)));
	y[i-1] =  center_y + radius*sin(dTheta*double((i-1)));
	}

	plot(x,y,nTheta + 1);

	delete [] x;
	delete [] y;
}

void fill(const CircleEntity &C)
{
	double center_x = C.getXcenter();
	double center_y = C.getYcenter();
	double radius   = C.getRadius();
	long i;
	long nTheta = 50;
	double dTheta = (2.0*SCC_PI)/double(nTheta);
	double* x = new double[nTheta +1];
	double* y = new double[nTheta +1];
	for(i = 1; i <= 51; i++)
	{
	x[i-1] =  center_x + radius*cos(dTheta*double((i-1)));
	y[i-1] =  center_y + radius*sin(dTheta*double((i-1)));
	}

	region(x,y,nTheta + 1);

	delete [] x;
	delete [] y;
}

void draw(const CombinedEntity& E)
{
	CircleEntity       C;
	XYrectangleEntity  R;
	PolygonEntity      P;
	for(long i = 0; i < E.getEntityCount(); i++)
	{
     if(E[i].getEntityType() == SCC::GeoType::CIRCLE)
     {
        draw(static_cast<CircleEntity&>(E[i]));
     }
     else if(E[i].getEntityType() == SCC::GeoType::XY_RECTANGLE)
     {
        draw(static_cast<XYrectangleEntity&>(E[i]));
     }
     else if(E[i].getEntityType() == SCC::GeoType::POLYGON)
     {
        draw(static_cast<PolygonEntity&>(E[i]));
     }
     else if(E[i].getEntityType() == SCC::GeoType::COMBINED)
     {
        draw(static_cast<CombinedEntity&>(E[i]));
     }
	}
}

void fill(const CombinedEntity& E)
{
	CircleEntity       C;
	XYrectangleEntity  R;
	PolygonEntity      P;
	for(long i = 0; i < E.getEntityCount(); i++)
	{
     if(E[i].getEntityType() == SCC::GeoType::CIRCLE)
     {
        fill(static_cast<CircleEntity&>(E[i]));
     }
     else if(E[i].getEntityType() == SCC::GeoType::XY_RECTANGLE)
     {
        fill(static_cast<XYrectangleEntity&>(E[i]));
     }
     else if(E[i].getEntityType() == SCC::GeoType::POLYGON)
     {
        fill(static_cast<PolygonEntity&>(E[i]));
     }
     else if(E[i].getEntityType() == SCC::GeoType::COMBINED)
     {
        fill(static_cast<CombinedEntity&>(E[i]));
     }
	}
}


};
}

#undef SCC_PI
#endif 
