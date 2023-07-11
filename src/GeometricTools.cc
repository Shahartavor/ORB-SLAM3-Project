/**
* This file is part of ORB-SLAM3
*
* Copyright (C) 2017-2021 Carlos Campos, Richard Elvira, Juan J. Gómez Rodríguez, José M.M. Montiel and Juan D. Tardós, University of Zaragoza.
* Copyright (C) 2014-2016 Raúl Mur-Artal, José M.M. Montiel and Juan D. Tardós, University of Zaragoza.
*
* ORB-SLAM3 is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
* License as published by the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* ORB-SLAM3 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
* the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License along with ORB-SLAM3.
* If not, see <http://www.gnu.org/licenses/>.
*/

#include "GeometricTools.h"
#include "Thirdparty/PolyToEigen/include/Poly.h"
#include "KeyFrame.h"

namespace ORB_SLAM3
{

Eigen::Matrix3f GeometricTools::ComputeF12(KeyFrame* &pKF1, KeyFrame* &pKF2)
{
    Sophus::SE3<float> Tc1w = pKF1->GetPose();
    Sophus::Matrix3<float> Rc1w = Tc1w.rotationMatrix();
    Sophus::SE3<float>::TranslationMember tc1w = Tc1w.translation();

    Sophus::SE3<float> Tc2w = pKF2->GetPose();
    Sophus::Matrix3<float> Rc2w = Tc2w.rotationMatrix();
    Sophus::SE3<float>::TranslationMember tc2w = Tc2w.translation();

    Sophus::Matrix3<float> Rc1c2 = Rc1w * Rc2w.transpose();
    Eigen::Vector3f tc1c2 = -Rc1c2 * tc2w + tc1w;

    Eigen::Matrix3f tc1c2x = Sophus::SO3f::hat(tc1c2);

    const Eigen::Matrix3f K1 = pKF1->mpCamera->toK_();
    const Eigen::Matrix3f K2 = pKF2->mpCamera->toK_();

    return K1.transpose().inverse() * tc1c2x * Rc1c2 * K2.inverse();
}

bool GeometricTools::Triangulate(Eigen::Vector3f &x_c1, Eigen::Vector3f &x_c2,
Eigen::Matrix<float,3,4> &Tc1w ,Eigen::Matrix<float,3,4> &Tc2w ,   Eigen::Matrix3f &mK_1, Eigen::Matrix3f &mK_2,
Eigen::Matrix<float,3,3> &Rcw1 ,Eigen::Matrix<float,3,3> &Rcw2,Eigen::Vector3f &tcw1,
Eigen::Vector3f &tcw2 ,Eigen::Vector3f &x3D)
{
  Eigen::VectorXf point1(2);
  point1 << x_c1[0], x_c1[1];
  Eigen::VectorXf point2(2);
 point2 << x_c2[0], x_c2[1];
  //std::cout << "new 8 parameters function" << std::endl;
Triangulation::Poly p(Tc1w,Tc2w,mK_1,mK_2,Rcw1,Rcw2,tcw1,tcw2);

    Eigen::Vector3f result = p.triangulate(point1,point2);
     if(result.isZero())
    {
        std::cout << "error in triangulate" << std::endl;
            return false;
    }
    x3D = result;
    return  true;
    
}
    bool GeometricTools::Triangulate(Eigen::Vector3f &x_c1, Eigen::Vector3f &x_c2,
                                     Eigen::Matrix<float,3,4> &Tc1w ,Eigen::Matrix<float,3,4> &Tc2w ,   Eigen::Matrix3f &mK_1, Eigen::Matrix3f &mK_2,
                                     Eigen::Matrix<float,3,3> &Rcw1 ,Eigen::Matrix<float,3,3> &Rcw2,Eigen::Vector3f &tcw1,
                                     Eigen::Vector3f &tcw2 ,Eigen::Vector3f &x3D,Eigen::Matrix3f &F)
    {
        Eigen::VectorXf point1(2);
        point1 << x_c1[0], x_c1[1];
        Eigen::VectorXf point2(2);
        point2 << x_c2[0], x_c2[1];
        //std::cout << "new 8 parameters function" << std::endl;
        Triangulation::Poly p(Tc1w,Tc2w,mK_1,mK_2,Rcw1,Rcw2,tcw1,tcw2,F);

        Eigen::Vector3f result = p.triangulate(point1,point2);
        if(result.isZero())
        {
            std::cout << "error in triangulate" << std::endl;
            return false;
        }
        x3D = result;
        return  true;

    }
    bool GeometricTools::Triangulate(Eigen::Vector3f &x_c1, Eigen::Vector3f &x_c2,Eigen::Matrix<float,3,4> &Tc1w ,Eigen::Matrix<float,3,4> &Tc2w , Eigen::Vector3f &x3D)
    {
        Eigen::Matrix4f A;
        A.block<1,4>(0,0) = x_c1(0) * Tc1w.block<1,4>(2,0) - Tc1w.block<1,4>(0,0);
        A.block<1,4>(1,0) = x_c1(1) * Tc1w.block<1,4>(2,0) - Tc1w.block<1,4>(1,0);
        A.block<1,4>(2,0) = x_c2(0) * Tc2w.block<1,4>(2,0) - Tc2w.block<1,4>(0,0);
        A.block<1,4>(3,0) = x_c2(1) * Tc2w.block<1,4>(2,0) - Tc2w.block<1,4>(1,0);

        Eigen::JacobiSVD<Eigen::Matrix4f> svd(A, Eigen::ComputeFullV);

        Eigen::Vector4f x3Dh = svd.matrixV().col(3);

        if(x3Dh(3)==0)
            return false;

        // Euclidean coordinates
        x3D = x3Dh.head(3)/x3Dh(3);

        return true;
    }
/*
bool GeometricTools::Triangulate(Eigen::Vector3f &x_c1, Eigen::Vector3f &x_c2,Eigen::Matrix<float,3,4> &Tc1w ,Eigen::Matrix<float,3,4> &Tc2w , Eigen::Vector3f &x3D)
{
    Eigen::Matrix4f A;
    std::cout << "old geometric tools " << std::endl;  
    //We got a match of 2 points from triangulation - x_c1 and x_c2.
    //x_c1 is point1 =(y,x) , x_c2 is point2 =(y',x')
    //Tc1w is camera matrix1= p, Tc2w is camera matrix2= p'

    //first row of A: x*p3 (row number 2 in p) - p1 (row number 0 in p)
    A.block<1,4>(0,0) = x_c1(0) * Tc1w.block<1,4>(2,0) - Tc1w.block<1,4>(0,0);
    //second row of A: y*p3 (row number 2 in p) - p2 (row number 1 in p)
    A.block<1,4>(1,0) = x_c1(1) * Tc1w.block<1,4>(2,0) - Tc1w.block<1,4>(1,0);

    //third row of A: x'*p3' (row number 2 in p') - p1' (row number 0 in p')
    A.block<1,4>(2,0) = x_c2(0) * Tc2w.block<1,4>(2,0) - Tc2w.block<1,4>(0,0);
    //fourth row of A: y'*p3' (row number 2 in p') - p2' (row number 1 in p')
    A.block<1,4>(3,0) = x_c2(1) * Tc2w.block<1,4>(2,0) - Tc2w.block<1,4>(1,0);
    
    //A matrix is:
    //  [  y*p3 - p2    ]
    //  |  p1 - x*p3    |
    //  |  y'*p3'-p2'   |
    //  [  p1' - x'p3'  ]
    Eigen::JacobiSVD<Eigen::Matrix4f> svd(A, Eigen::ComputeFullV);
    //JacobiSVD= U * sigma * V^t where U and V 
    //Soultion x is the column of V corresponding to the smallest signular value  =col3.
    
    Eigen::Vector4f x3Dh = svd.matrixV().col(3);

    if(x3Dh(3)==0)
        return false;

    // Euclidean coordinates

    //divide by Z to return Euclidean coordinates
    x3D = x3Dh.head(3)/x3Dh(3);

    return true;
}
*/

} //namespace ORB_SLAM
