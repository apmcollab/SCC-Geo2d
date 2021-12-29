/*
 * SCC_LinkUtiilty.h
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

#ifndef SCC_LINK_
#define SCC_LINK_

namespace SCC
{
class LinkUtility
{
public :
    
int onLink(double xa, double ya, double xb, double yb,
double xp, double yp, double tol, double& s, double& d) const
{
//
// Returns 1 if the point (xp, yp) is "on" the link
// from (xa,ya) to (xb,yb). ("on" is defined as the point
// being within the tolerance region for the link).
//
// s is the normalized parametric coordinate of the point with the
// the segments parameterization being (xa,ya)*(1-s) + (xb,yb)*s
//
// d is the distance of (xp,yp) to the point on the line with
// parametric coordinate s (which will not necessarily be zero
// because of non-zero tolerance).
//
// Returns 0 if the point is not on the segment.
//
    double xab  = (xa - xb);
    double yab  = (ya - yb);

	s = ((xa - xp)*xab + (ya - yp)*yab)/(xab*xab + yab*yab);

    double xpl  =     xp - (xa*(1.0 - s) + s*xb);
    double ypl  =     yp - (ya*(1.0 - s) + s*yb);
    double  d2  =     xpl*xpl + ypl*ypl;
    double tol2 =     tol*tol;

    if(((0.0 <= s)&&(s <= 1.0))&&(d2 < tol2))
    {d = std::sqrt(d2); return 1;}

    double xpb = (xp - xb);
    double ypb = (yp - yb);

    if((1.0 < s)&&((xpb*xpb + ypb*ypb) < tol2))
    {d = std::sqrt(xpb*xpb + ypb*ypb); return 1;}

    double xpa = (xp - xa);
    double ypa = (yp - ya);

    if((s <  0.0)&&((xpa*xpa + ypa*ypa) < tol2))
    {d = std::sqrt(xpa*xpa + ypa*ypa); return 1;}

    return 0;
}

int onConnectedLinks(const double* X, const double* Y, long nLinks,
double xp, double yp, double tol, long& linkIndex, double& s) const
{
//
// Returns 1 if the point (xp, yp) is on the collection of
// connected links whose vertex coordinates are given by the
// arrays (X[0], Y[0]) i = 0 .. nLinks.
//
// The arrays X and Y are of size nLinks
//
// The link index linkIndex is set to the link upon which
// the point (xp, yp) lies, and the value of s is set to the
// normalized parametric coordinate of (xp,yp) on that link.
//
// Returns 0 if the point is not on the link collection
//
    double d;
    long indexA, indexB;
    double sA,sB;
    double dA,dB;


    int  hits  = 0;
    long i     = 0;

    double sa, sb;
 //
 // Check for Hits :
 //
    while((i < nLinks)&&(hits <= 2))
    {
    if(onLink(X[i],Y[i], X[i+1], Y[i+1],xp, yp,tol, s,d))
    {
    hits++;
    if(hits == 1)
    {indexA = i+1; sA     = s; dA = d;}
    if(hits == 2)
    {indexB = i+1; sB     = s; dB = d;}
    }
    i++;
    }
    if(hits == 1)
    {
       linkIndex = indexA;
       s            = sA;
       return 1;
    }
    if(hits == 2)
    {
       if
       (((0.0 <= sA)&&(sA <= 1.0))
       &&
       ((0.0 <= sB)&&(sB <= 1.0)))
       {
       if(dA == dB)
       {
         if(sA < sB){ linkIndex = indexA; s = sA;}
         else       { linkIndex = indexB; s = sB;}
       }
       else if(dA < dB) { linkIndex = indexA; s = sA;}
       else { linkIndex = indexB; s = sB;}
       return 1;
       }
       else if
       (((0.0 <= sA)&&(sA <= 1.0))
       &&
       ((sB < 0.0)||(sB > 1.0)))
       { linkIndex = indexA; s = sA; return 1;}
       else if
       (((0.0 <= sB)&&(sB <= 1.0))
       &&
       ((sA < 0.0)||(sA > 1.0)))
       { linkIndex = indexB; s = sB; return 1;}
       else
       {
       if(sA < 0.0) sa = -sA;
       if(sB < 0.0) sb = -sB;
       if(sA > 1.0) sa =  sA - 1.0;
       if(sB > 1.0) sb =  sB - 1.0;
       if(sa <  sb){ linkIndex = indexA; s = sA; return 1;}
       else { linkIndex = indexB; s = sB; return 1;}
       }
    }
    return 0;
}

int linkIntersection(const double* Qx, const double* Qy,const double* Rx , const double* Ry, double tol, double& s) const
{
//
// Returns 1 if the link defined by Q == link from
// (Qx[0], Qy[0]) to (Qx[1], Qy[1]) intersects the link defined
// by R == link from (Rx[0], Ry[0]) to (Rx[1], Ry[1])
//
// The notion of intersection is defined by the existence of a
// value of s_Q in [0,1] such that the point on Q with that
// parametric coordinate lies "on" the link R.
//
//  This intersection function is a directed intersection; if
//  the segments overlap, then the intersection point is that point
//  on Q which is closest to the base point of Q, i.e.
//  (Qx[0], Qy[0]).
//
//  (A) Check for overlap conditions
//  (B) Check for non-overlap intersection
//
    double d;

    double sQ0 = 0.0;
    double sQ1 = 0.0;
    double sR0 = 0.0;
    double sR1 = 0.0;
    int Q1on   = 0;
    int R0on   = 0;
    int R1on   = 0;

	if(onLink(Rx[0],Ry[0],Rx[1],Ry[1],Qx[0],Qy[0],tol,sQ0,d))
    {
    s = 0.0;
    return 1;
    }
    else
    {
     if(onLink(Rx[0],Ry[0],Rx[1],Ry[1],Qx[1],Qy[1],tol,sQ1,d))
     {Q1on = 1; sQ1 = 1.0;}

     if(onLink(Qx[0],Qy[0],Qx[1],Qy[1],Rx[0],Ry[0],tol,sR0,d))
     {R0on = 1;}

     if(onLink(Qx[0],Qy[0],Qx[1],Qy[1],Rx[1],Ry[1],tol,sR1,d))
     {R1on = 1;}

     if(R0on && (!(R1on))){ s = normalize(sR0); return 1;}
     else
     if(R1on && (!(R0on))){ s = normalize(sR1); return 1;}
     else
     if(R0on && R1on)
     {s = sR1; if (sR0 <= sR1) s = normalize(sR0); return 1;}
     else if(Q1on)
     {s = 1.0; return 1;}
    }
//
//  No link overlap --- so see if there is a link intersection
//
	double A[2][2];
    double B[2];

    A[0][0] =  (Qx[1] - Qx[0]);
    A[0][1] = -(Rx[1] - Rx[0]);
    A[1][0] =  (Qy[1] - Qy[0]);
    A[1][1] = -(Ry[1] - Ry[0]);

    B[0] = Rx[0] - Qx[0];
    B[1] = Ry[0] - Qy[0];

    double rScale0 = std::sqrt(A[0][0]*A[0][0] + A[0][1]*A[0][1]);
    double rScale1 = std::sqrt(A[1][0]*A[1][0] + A[1][1]*A[1][1]);

    if((rScale0 < tol) ||( rScale1 < tol))
    {
    /* then both lines parallel to on or other axis so exit */
    return 0;
    }

    rScale0 =  std::abs(A[0][0])/rScale0;
    rScale1 =  std::abs(A[1][0])/rScale1;

    if(rScale0 < rScale1) // pivot
    {
    A[1][0] =  (Qx[1] - Qx[0]);
    A[1][1] = -(Rx[1] - Rx[0]);
    A[0][0] =  (Qy[1] - Qy[0]);
    A[0][1] = -(Ry[1] - Ry[0]);

    B[1] = Rx[0] - Qx[0];
    B[0] = Ry[0] - Qy[0];
    }

    A[1][1] = A[1][1] - (A[1][0]/A[0][0])*A[0][1];
    B[1]    = B[1]    - (A[1][0]/A[0][0])*B[0];

    if(std::abs(A[1][1])== 0.0)
    {
    /* then parallel lines  so exit */
    return 0;
    }

    double s1 =  B[1]/A[1][1];
    double s0 = (B[0] - A[0][1]*s1)/A[0][0];

    if(((s1 > 0.0)&&(s1 < 1.0))&&((s0 > 0.0)&&(s0 < 1.0)))
    {s = s0; return 1;}

    return 0;
    }

    private :

    double normalize(double sVal) const
    {
	double s = sVal;
	if(s <= 0.0) s = 0.0;
	else if(s >= 1.0) s = 1.0;
	return s;
	}
};
} //Namespace SCC
#endif

