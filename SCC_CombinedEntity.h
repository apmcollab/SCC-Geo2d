/*
 * SCC_CombinedEntity.h
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
#include "SCC_GeometricEntity.h"
#include "SCC_CircleEntity.h"
#include "SCC_XYrectangleEntity.h"
#include "SCC_PolygonEntity.h"

#ifndef SCC_COMBINEDENTITY_
#define SCC_COMBINEDENTITY_

namespace SCC
{
class CombinedEntity : public GeometricEntity
{
private :

//
//  Internal representation data using base class pointers
//
//
    std::vector<GeometricEntity*> entityList;
    std::vector<std::string>  entityNameList;
    long  entityCount;

public :

CombinedEntity()
{
	entityList.clear();
	entityNameList.clear();
	entityCount = 0;
}

CombinedEntity(const CombinedEntity& E)
{
	initialize(E);
}

//
//********************************************************************************
//                    DESTRUCTOR
//********************************************************************************
//

virtual ~CombinedEntity()
{
	deleteData();
}

//********************************************************************************
//                    INITIALIZATION
//********************************************************************************


void initialize()
{
	deleteData();
}

void initialize(const CombinedEntity& E)
{
    initialize();
    if(E.entityList.size() == 0) {return;}

    entityCount = (long)E.entityList.size();

    entityList.resize(entityCount, nullptr);
    entityNameList.resize(entityCount,"");

    for(long i = 0; i < entityCount; i++)
    {
    entityList[i]     = E.entityList[i]->newDuplicateEntity();
    entityNameList[i] = E.entityNameList[i];
    }
}

//
//********************************************************************************
//                    ASSIGNMENT
//********************************************************************************
//


CombinedEntity& operator =(const CombinedEntity& A)
{
	this->initialize(A);
	return *this;
}

//
//********************************************************************************
//                    Creation
//********************************************************************************
//
GeometricEntity*  newDuplicateEntity() const
{
	CombinedEntity* R = new CombinedEntity(*this);
	return  R;
}

//
//********************************************************************************
//                    OUTPUT
//********************************************************************************
//
friend std::ostream&  operator  <<(std::ostream& out_stream, const SCC::CombinedEntity& A)
{
   out_stream << "[BEGIN_ENTITY]" << "\n";
   out_stream <<  GeometricEntity::getGeoTypeString(SCC::GeoType::COMBINED) << "\n";
   out_stream << "[ENTITY_DATA]" << "\n";

   for(long i = 0; i < A.getEntityCount(); i++)
   {
     out_stream << A[i];
   }
   out_stream << "[END_ENTITY]\n";
   return(out_stream);
}


friend std::istream& operator >>(std::istream& in_stream, SCC::CombinedEntity& A)
{
    CircleEntity         C;
    XYrectangleEntity    R;
    PolygonEntity        P;

	std::string str;
	char delim = '\n';
    std::getline(in_stream,str,delim); // [BEGIN_ENTITY]
    std::getline(in_stream,str,delim); // SCC_GeoType_COMBINED
    std::getline(in_stream,str,delim); // [ENTITY_DATA]

    while(!(in_stream.eof()))
    {
    std::getline(in_stream,str,delim); // [BEGIN_ENTITY]
    if(str.compare("[END_ENTITY]") == 0) break;

    std::getline(in_stream,str,delim); // SCC_GeoType_XXX
    if(str.compare("SCC_GeoType_CIRCLE") == 0)
    {
        C.initialize();
    	C.inputData(in_stream);
    	A.addEntity("SCC_GeoType_CIRCLE",C);
        std::getline(in_stream,str,delim);
    }
    else if(str.compare("SCC_GeoType_XY_RECTANGLE") == 0)
    {
        R.initialize();
    	R.inputData(in_stream);
    	A.addEntity("SCC_GeoType_XY_RECTANGLE",R);
        std::getline(in_stream,str,delim);
    }
    else if(str.compare("SCC_GeoType_POLYGON") == 0)
    {
        P.initialize();
    	P.inputData(in_stream);
    	A.addEntity("SCC_GeoType_POLYGON",P);
    	std::getline(in_stream,str,delim);
    }


    }


    /*
    std::getline(in_stream,str,delim); // [ENTITY_DATA]

    in_stream.getline(nameInput,256);   // clear eol
    in_stream.getline(nameInput,256);   // Name of entry
    if(strcmp(nameInput,"[END_ENTITY]")==0) return(in_stream);

    in_stream >> lineInput;
    while(strcmp(lineInput,"[BEGIN_ENTITY]")==0)
    {
    in_stream >> lineInput;  // class name

    if(strcmp(lineInput,"CAMcircleEntity") == 0)
    {
      C.inputData(in_stream);
      A.addEntity(nameInput,C);
    }
    if(strcmp(lineInput,"CAMxyRectangleEntity") == 0)
    {
      R.inputData(in_stream);
      A.addEntity(nameInput,R);
    }
    if(strcmp(lineInput,"CAMpolygonEntity") == 0)
    {
      P.inputData(in_stream);
      A.addEntity(nameInput,P);
    }

    in_stream.getline(nameInput,256);    // clear eol
    in_stream.getline(nameInput,256);
    if(strcmp(nameInput,"[END_ENTITY]")==0) return(in_stream);
	 in_stream >> lineInput;
    }
    */
    return(in_stream);
}

//
//********************************************************************************
//                    Equality/Inequality
//********************************************************************************
//



bool  equals(const CombinedEntity &E) const
{
	if(entityCount != E.entityCount) return false;

	bool flag = true;
	for(int i = 0; i < entityCount; i++)
	{
	if(*entityList[i] !=  *E.entityList[i]) flag = false;
	}
	return flag;
}
bool notEquals(const CombinedEntity &E) const
{
	return !(*this == E);
}

//
//******************************************************************************
//                       ACCESS
//******************************************************************************
//
//

SCC::GeoType getEntityType() const {return SCC::GeoType::COMBINED;};

long getEntityCount() const
{
	return entityCount;
}
void  addEntity(const GeometricEntity& E)
{
	entityCount++;
    entityList.push_back(E.newDuplicateEntity());
    entityNameList.push_back(GeometricEntity::getGeoTypeString(E.getEntityType()));
}


void  addEntity(const std::string& entityName, const GeometricEntity& E)
{
	entityCount++;
    entityList.push_back(E.newDuplicateEntity());
    entityNameList.push_back(entityName);
}


GeometricEntity& operator[](std::string& entityName) const
{
    long i = 0;
    long i_hit = -1;

	while((i < entityCount)&&(i_hit < 0))
	if(!(entityName.compare(entityNameList[i]) == 0))
    { i++;}
    else
    {i_hit = i;}

	if(i > entityCount-1)
	{
	throw std::runtime_error("\n CombinedEntity error : Invalid CombinedEntity string index in operator[] \n");
	}

	return *(entityList[i]);
}

GeometricEntity&  operator[](long  entityIndex) const
{
	if((entityIndex < 0)||(entityIndex > entityCount - 1))
	{
	throw std::runtime_error("\n CombinedEntity error : Invalid CombinedEntity integer index in operator[] \n");
	}

	return *entityList[entityIndex];
}

std::string getIndexName(long index) const
{
	if((index < 0)||(index > entityCount - 1))
	{
	throw std::runtime_error("\n CombinedEntity error : Invalid CombinedEntity integer index in operator[] \n");
	}

	return entityNameList[index];
}

void  setIndexName(long index, const std::string& newName)
{
	if((index < 0)||(index > entityCount - 1))
	{
	throw std::runtime_error("\n CombinedEntity error : Invalid CombinedEntity integer index in operator[] \n");
	}

	entityNameList[index] = newName;
}


void scaleBy(double alpha)
/*
@begin_doc@
@title    @ scaleBy(double alpha)
@Purpose  @
@end_doc  @
*/
{
	for(size_t i = 0; i < entityList.size(); i++)
	{
	entityList[i]->scaleBy(alpha);
	}
}

void  rotateBy(double theta)
/*
@begin_doc@
@title    @ rotateBy(double theta)
@Purpose  @
@end_doc  @
*/
{
	for(size_t i = 0; i < entityList.size(); i++)
	{
	entityList[i]->rotateBy(theta);
	}
}

void  translateBy(double x, double y)
/*
@begin_doc@
@title    @ translateBy(double x, double y)
@Purpose  @
@end_doc  @
*/
{
	for(size_t i = 0; i < entityList.size(); i++)
	{
	entityList[i]->translateBy(x,y);
	}
}

void translateTo(double x, double y)
/*
@begin_doc@
@title    @ translateTo(double #include "camgeoimpexp.h"x, double y)
@Purpose  @
@end_doc  @
*/
{
	for(size_t i = 0; i < entityList.size(); i++)
	{
	entityList[i]->translateTo(x,y);
	}
}

int interiorExteriorTest(double xTest, double yTest) const
//
//  Returns +1 if test point is interior to at least one of the elements
//  in the combined entity, -1 if exterior to all elements, 0 if on one
//  of the elements.
//
{
    long k;
    long interiorExteriorFlag = -1;

    if(entityCount == 0) return interiorExteriorFlag;

    else
    {
    interiorExteriorFlag = operator[](0).interiorExteriorTest(xTest,yTest);
    k = 1;
    while((interiorExteriorFlag == -1)&&(k < entityCount))
    {
    interiorExteriorFlag = operator[](k).interiorExteriorTest(xTest,yTest);
    k++;
    }
    }

    return interiorExteriorFlag;
}

/*
@begin_doc@
@title    @ getDistanceToBoundary()
@Purpose  @ Returns signed minimal distance to the boundary of the collection of
entities. This routine is only defined for non-overlapping collections of entities.
If a point is inside one of the entites, then it returns the (+) distance to it's boundary,
otherwise it returns (-) the minimum distance to all of the entities. If the test
point is interior to more than one entity, an error message is generated and exit(1) is
called.
@end_doc  @
*/
double getDistanceToBoundary(double x, double y) const
{
    long k;
    if(entityCount == 0) return 0.0;

    long interiorCount        =  0;
    long interiorIndex        =  0;
    long interiorExteriorFlag = -1;

    for(k = 0; k < entityCount; k++)
    {
        interiorExteriorFlag = operator[](k).interiorExteriorTest(x,y);
        if(interiorExteriorFlag > 0)
        {
        interiorCount++;
        interiorIndex = k;
        }
    }


    if(interiorCount > 1)
    {
    std::string errMsg;
    errMsg = "\n CombinedEntity error : combinedEntity contains overlapping entities \n";
    errMsg = errMsg + "      : getDistanceToBoundary(...) is not implemented for  \n";
    errMsg = errMsg + "      : points inside overlapping entities.";
    throw std::runtime_error(errMsg);
    }


    if(interiorCount == 1)
    {
    return operator[](interiorIndex).getDistanceToBoundary(x,y);
    }

    //
    // point exterior to all entities, so take the -minimum |distance|
    // to the entities.
    //
    double dMin;

    dMin = std::abs(operator[](0).getDistanceToBoundary(x,y));
    double dTest;

    for(k = 1; k < entityCount; k++)
    {
    dTest = std::abs(operator[](k).getDistanceToBoundary(x,y));
    dMin = (dMin < dTest) ? dMin : dTest;
    }

    return -dMin;
}

/*
@begin_doc@
@title    @ getExteriorDistanceToBoundary()
@Purpose  @ For a point exterior to a collection of entities, this routine returns
(-) the minimum distance to all of the entities. For a point interior to one or
more entities this routine returns the value 1.0. This routine is well defined for
overlapping entities.
@end_doc  @
*/
double getExteriorDistanceToBoundary(double x, double y) const
{
    long k;
    if(entityCount == 0) return 0.0;

    long interiorCount        =  0;
    long interiorIndex        =  0;
    long interiorExteriorFlag = -1;

    for(k = 0; k < entityCount; k++)
    {
        interiorExteriorFlag = operator[](k).interiorExteriorTest(x,y);
        if(interiorExteriorFlag > 0)
        {
        interiorCount++;
        interiorIndex = k;
        }
    }

    if(interiorCount > 1)
    {
    return 1.0;
    }

    double dMin;
    if(interiorCount == 1)
    {
    dMin = operator[](interiorIndex).getDistanceToBoundary(x,y);
    if(dMin > 0.0){return 1.0; }
    else          {return dMin;}
    }

    //
    // point exterior to all entities, so take the -minimum |distance|
    // to the entities.
    //

    dMin = std::abs(operator[](0).getDistanceToBoundary(x,y));
    double dTest;

    for(k = 1; k < entityCount; k++)
    {
    dTest = std::abs(operator[](k).getDistanceToBoundary(x,y));
    dMin = (dMin < dTest) ? dMin : dTest;
    }

    return -dMin;
}


void  getBoundingBox(double& bMinX, double& bMaxX,
double& bMinY, double& bMaxY)
{
    long entityCount = getEntityCount();

    if(entityCount == 0)
    {
    bMinX = 0.0; bMinY = 0.0;
    bMaxX = 0.0; bMaxY = 0.0;
    return;
    }

    long k;

    double minX; double maxX;
    double minY; double maxY;

    bMinX = operator[](0).getMinX();
    bMaxX = operator[](0).getMaxX();
    bMinY = operator[](0).getMinY();
    bMaxY = operator[](0).getMaxY();

    for(k = 1; k < entityCount; k++)
    {
    minX = operator[](k).getMinX();
    maxX = operator[](k).getMaxX();
    minY = operator[](k).getMinY();
    maxY = operator[](k).getMaxY();
    bMinX = (bMinX > minX) ? minX : bMinX;
    bMaxX = (bMaxX < maxX) ? maxX : bMaxX;
    bMinY = (bMinY > minY) ? minY : bMinY;
    bMaxY = (bMaxY < maxY) ? maxY : bMaxY;
    }

    return;
}


private:

void deleteData()
{
 	if(entityList.size() != 0)
	{
		for(size_t i = 0; i < entityList.size(); i++)
		{
			delete entityList[i];
		}
	}
	entityList.clear();
	entityNameList.clear();
	entityCount = 0;
}

};

} // Namespace SCC




#endif /* SCC_COMBINEDENTITY_ */
