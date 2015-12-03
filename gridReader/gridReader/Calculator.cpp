//
//  Calculator.cpp
//  gridReader
//
//  Created by Franz Neubert on 02/11/15.
//  Copyright © 2015 Franz Neubert. All rights reserved.
//

#include <stdio.h>
#include <math.h>
#include <vtkPoints.h>
#include "Calculator.h"



void Calculator::calcCrossProduct(double *x, double *y, double *crossProduct) {
    double c1[3];
    double c2[3];
    
    c1[0] = x[1] * y[2];
    c1[1] = x[2] * y[0];
    c1[2] = x[0] * y[1];
    
    c2[0] = x[2] * y[1];
    c2[1] = x[0] * y[2];
    c2[2] = x[1] * y[0];
    
    Calculator::subtractVectors(c2, c1, crossProduct);
}

double Calculator::calcDotProduct(double *x, double *y) {
    return x[0]*y[0] + x[1]*y[1] + x[2]*y[2] ;
}

void Calculator::subtractVectors(double *x, double *y, double *result) {
    result[0] = y[0] - x[0];
    result[1] = y[1] - x[1];
    result[2] = y[2] - x[2];
}

float Calculator::calcVectorLength(double *vector) {
    // std::cout << "vector length: " << sqrt(pow(vector[0], 2) + pow(vector[1], 2) + pow(vector[2], 2)) << std::endl;
    return sqrt(pow(vector[0], 2) + pow(vector[1], 2) + pow(vector[2], 2));
}

void Calculator::calcMidPoint(double *A, double *B, double *midpoint) {
    midpoint[0] = (A[0]+B[0])/2;
    midpoint[1] = (A[1]+B[1])/2;
    midpoint[2] = (A[2]+B[2])/2;
}

double r8_acos ( double c )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ACOS computes the arc cosine function, with argument truncation.
//
//  Discussion:
//
//    If you call your system ACOS routine with an input argument that is
//    outside the range [-1.0, 1.0 ], you may get an unpleasant surprise.
//    This routine truncates arguments outside the range.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double C, the argument, the cosine of an angle.
//
//    Output, double R8_ACOS, an angle whose cosine is C.
//
{
    const double r8_pi = 3.141592653589793;
    double value;
    
    if ( c <= -1.0 )
    {
        value = r8_pi;
    }
    else if ( 1.0 <= c )
    {
        value = 0.0;
    }
    else
    {
        value = acos ( c );
    }
    return value;
}


double r8vec_dot ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DOT computes the dot product of a pair of R8VEC's in ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], A2[N], the two vectors to be considered.
//
//    Output, double R8VEC_DOT, the dot product of the vectors.
//
{
    int i;
    double value;
    
    value = 0.0;
    for ( i = 0; i < n; i++ )
    {
        value = value + a1[i] * a2[i];
    }
    
    return value;
}






void tetrahedron_edges ( double tetra[3*4], double ab[], double ac[],
                        double ad[], double bc[], double bd[], double cd[] )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_EDGES returns the edges of a tetrahedron.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 May 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double TETRA[3*4], the tetrahedron vertices.
//
//    Output, double AB[3], AC[3], AD[3], BC[3], BD[3], CD[3], the edges.
//
{
    int i;
    //
    //  Compute the vectors that represent the sides.
    //
    for ( i = 0; i < 3; i++ )
    {
        ab[i] = tetra[i+1*3] - tetra[i+0*3];
        ac[i] = tetra[i+2*3] - tetra[i+0*3];
        ad[i] = tetra[i+3*3] - tetra[i+0*3];
        bc[i] = tetra[i+2*3] - tetra[i+1*3];
        bd[i] = tetra[i+3*3] - tetra[i+1*3];
        cd[i] = tetra[i+3*3] - tetra[i+2*3];
    }
    
    return;
}

double r8vec_angle_3d ( double u[], double v[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_ANGLE_3D computes the angle between two vectors in 3D.
//
//  Modified:
//
//    07 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double U[3], V[3], the vectors.
//
//    Output, double ANGLE, the angle between the two vectors.
//
{
    double angle;
    double angle_cos;
    double u_norm;
    double uv_dot;
    double v_norm;
    
    uv_dot = r8vec_dot ( 3, u, v );
    
    u_norm = sqrt ( r8vec_dot ( 3, u, u ) );
    
    v_norm = sqrt ( r8vec_dot ( 3, v, v ) );
    
    angle_cos = uv_dot / u_norm / v_norm;
    
    angle = r8_acos ( angle_cos );
    
    return angle;
}

double *r8vec_cross_3d ( double v1[3], double v2[3] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_CROSS_3D computes the cross product of two R8VEC's in 3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double V1[3], V2[3], the coordinates of the vectors.
//
//    Output, double R8VEC_CROSS_3D[3], the cross product vector.
//
{
    double *v3;
    
    v3 = new double[3];
    
    v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
    v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
    v3[2] = v1[0] * v2[1] - v1[1] * v2[0];
    
    return v3;
}

double *tetrahedron_dihedral_angles ( double tetra[] )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_DIHEDRAL_ANGLES computes dihedral angles of a tetrahedron.
//
//  Modified:
//
//    08 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, real ( kind = 8 ) TETRA(3,4), the vertices of the tetrahedron,
//    which can be labeled as A, B, C and D.
//
//    Output, double TETRAHEDRON_DIHEDRAL_ANGLES[6], the dihedral angles
//    along the axes AB, AC, AD, BC, BD and CD, respectively.
//
{
    double ab[3];
    double *abc_normal;
    double *abd_normal;
    double ac[3];
    double *acd_normal;
    double ad[3];
    double *angle;
    double bc[3];
    double *bcd_normal;
    double bd[3];
    double cd[3];
    int i;
    const double r8_pi = 3.141592653589793;
    
    tetrahedron_edges ( tetra, ab, ac, ad, bc, bd, cd );
    
    abc_normal = r8vec_cross_3d ( ac, ab );
    abd_normal = r8vec_cross_3d ( ab, ad );
    acd_normal = r8vec_cross_3d ( ad, ac );
    bcd_normal = r8vec_cross_3d ( bc, bd );
    
    angle = new double[6];
    
    angle[0] = r8vec_angle_3d ( abc_normal, abd_normal );
    angle[1] = r8vec_angle_3d ( abc_normal, acd_normal );
    angle[2] = r8vec_angle_3d ( abd_normal, acd_normal );
    angle[3] = r8vec_angle_3d ( abc_normal, bcd_normal );
    angle[4] = r8vec_angle_3d ( abd_normal, bcd_normal );
    angle[5] = r8vec_angle_3d ( acd_normal, bcd_normal );
    
    for ( i = 0; i < 6; i++ )
    {
        angle[i] = r8_pi - angle[i];
    }
    
    delete [] abc_normal;
    delete [] abd_normal;
    delete [] acd_normal;
    delete [] bcd_normal;
    
    return angle;
}

double *tetrahedron_solid_angles ( double tetra[] )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_SOLID_ANGLES computes solid angles of a tetrahedron.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double TETRA[3*4], the vertices of the tetrahedron.
//
//    Output, double TETRAHEDRON_SOLID_ANGLES[4], the solid angles.
//
{
    double *angle;
    double *dihedral_angles;
    const double r8_pi = 3.141592653589793;
    
    dihedral_angles = tetrahedron_dihedral_angles ( tetra );
    
    angle = new double[4];
    
    angle[0] = dihedral_angles[0]
    + dihedral_angles[1]
    + dihedral_angles[2] - r8_pi;
    
    angle[1] = dihedral_angles[0]
    + dihedral_angles[3]
    + dihedral_angles[4] - r8_pi;
    
    angle[2] = dihedral_angles[1]
    + dihedral_angles[3]
    + dihedral_angles[5] - r8_pi;
    
    angle[3] = dihedral_angles[2]
    + dihedral_angles[4]
    + dihedral_angles[5] - r8_pi;
    
    delete [] dihedral_angles;
    
    return angle;
}

double Calculator::calcSolidAngle(vtkCell *tetra, vtkIdType origin) {
    vtkPoints *points = tetra->GetPoints();
    double tetraCoords[3*4];
    points->GetPoint(tetra->GetPointId(0), &tetraCoords[0]);
    points->GetPoint(tetra->GetPointId(1), &tetraCoords[3]);
    points->GetPoint(tetra->GetPointId(2), &tetraCoords[6]);
    points->GetPoint(tetra->GetPointId(3), &tetraCoords[9]);
    
    double *angles;
    angles = tetrahedron_solid_angles(tetraCoords);
    
    return tetrahedron_solid_angles(tetraCoords)[sizeof(double)*origin];
}




