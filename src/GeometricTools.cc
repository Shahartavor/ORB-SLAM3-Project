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

bool GeometricTools::Triangulate(Eigen::Vector3f &x_c1, Eigen::Vector3f &x_c2,Eigen::Matrix<float,3,4> &Tc1w ,Eigen::Matrix<float,3,4> &Tc2w , Eigen::Vector3f &x3D)
{
    Eigen::Matrix4f A;  
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

    bool GeometricTools::TriangulatePoly(Eigen::Vector3f &x_c1, Eigen::Vector3f &x_c2,Eigen::Matrix<float,3,4> &Tc1w ,Eigen::Matrix<float,3,4> &Tc2w , Eigen::Vector3f &x3D) {

        //We got a match of 2 points from triangulation - x_c1 and x_c2.
        //u = x_c1 is point1 =(y,x) , u' = x_c2 is point2 =(y',x')
        //Tc1w is camera matrix1= p, Tc2w is camera matrix2= p'

        // params is polynomial 6th degrees
        //1.calculate 2 epoipolar lines: e =essentialMatrix*x
        //2.function calcluate distance from point to line : x1,y1 - point , a,b,c-parameters of line
        //    void shortest_distance(float x1, float y1, float a, float b,
        //                           float c)
        //    {
        //        float d = fabs((a * x1 + b * y1 + c))
        //                  / (sqrt(a * a + b * b));
        //        cout << "Perpendicular distance is, " << d << endl;
        //        return;
        //    }
        //3.params = distance of u from u epipole + distance of u' from u' epipole
        //4. poly:
        std::vector<double> Poly::PreparePolyCoeffs(const Poly::PolyParams &params) const {
            double a, b, c, d, e, f;
            std::tie(a, b, c, d, e, f) = params;
            std::vector<double> result = {
                    // a*b*d^2 - b^2*c*d
                    a * b * d * d - b * b * c * d,
                    // - 2*b^2*d^2*f^2*z + a^2*d^2*z - d^4*f^4*z - b^2*c^2*z - b^4*z
                    -2 * b * b * d * d * f * f + a * a * d * d - d * d * d * d * f * f * f * f - b * b * c * c -
                    b * b * b * b,
                    // - 4*b^2*c*d*f^2*z^2 - 2*b^2*c*d*e^2*z^2 - 4*a*b*d^2*f^2*z^2 + 2*a*b*d^2*e^2*z^2 + a^2*c*d*z^2 - 4*c*d^3*f^4*z^2 - a*b*c^2*z^2 - 4*a*b^3*z^2
                    -4 * b * b * c * d * f * f - 2 * b * b * c * d * e * e - 4 * a * b * d * d * f * f +
                    2 * a * b * d * d * e * e + a * a * c * d - 4 * c * d * d * d * f * f * f * f - a * b * c * c -
                    4 * a * b * b * b,
                    // - 8*a*b*c*d*f^2*z^3 - 6*c^2*d^2*f^4*z^3 - 2*b^2*c^2*f^2*z^3 - 2*a^2*d^2*f^2*z^3 - 2*b^2*c^2*e^2*z^3 + 2*a^2*d^2*e^2*z^3 - 6*a^2*b^2*z^3
                    -8 * a * b * c * d * f * f - 6 * c * c * d * d * f * f * f * f - 2 * b * b * c * c * f * f -
                    2 * a * a * d * d * f * f - 2 * b * b * c * c * e * e + 2 * a * a * d * d * e * e -
                    6 * a * a * b * b,
                    // - 4*a^2*c*d*f^2*z^4 + 2*a^2*c*d*e^2*z^4 - 4*a*b*c^2*f^2*z^4 + a*b*d^2*e^4*z^4 - 2*a*b*c^2*e^2*z^4 - b^2*c*d*e^4*z^4 - 4*c^3*d*f^4*z^4 - 4*a^3*b*z^4
                    -4 * a * a * c * d * f * f + 2 * a * a * c * d * e * e - 4 * a * b * c * c * f * f +
                    a * b * d * d * e * e * e * e - 2 * a * b * c * c * e * e - b * b * c * d * e * e * e * e -
                    4 * c * c * c * d * f * f * f * f - 4 * a * a * a * b,
                    // a^2*d^2*e^4*z^5 - 2*a^2*c^2*f^2*z^5 - b^2*c^2*e^4*z^5 - c^4*f^4*z^5 - a^4*z^5
                    a * a * d * d * e * e * e * e - 2 * a * a * c * c * f * f - b * b * c * c * e * e * e * e -
                    c * c * c * c * f * f * f * f - a * a * a * a,
                    // a^2*c*d*e^4*z^6 - a*b*c^2*e^4*z^6
                    a * a * c * d * e * e * e * e - a * b * c * c * e * e * e * e
            };
            result.resize(FindPolynomialOrder(result) + 1);
            double max_coeff = *std::max_element(result.begin(), result.end());
            if (max_coeff > 0)
            {
                std::transform(result.begin(), result.end(), result.begin(), [&max_coeff](double& a){ return a / max_coeff; });
            }
            else if (max_coeff < 0)
            {
                double min_coeff = *std::min_element(result.begin(), result.end());
                std::transform(result.begin(), result.end(), result.begin(), [&min_coeff](double& a){ return a / min_coeff; });
            }
            return result;
        }
    }
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

} //namespace ORB_SLAM
