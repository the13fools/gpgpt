  
  #include "UtilsMisc.h"


  #include <Eigen/Core>
  #include <Eigen/Dense>



//   Eigen::Matrix4d rstar_from_r(Eigen::Matrix2d r)
//     {
//       Eigen::Matrix4d ret;
//       double a = r(0,0);
//       double b = r(0,1);
//       double c = r(1,0);
//       double d = r(1,1);
//       ret << a*a, a*b, a*b, b*b,
//              a*c, b*c, a*d, b*d,
//              a*c, b*c, a*d, b*d,
//              c*c, c*d, c*d, d*d;

//       // ret << r(0,0), r(0, 1), r(1,0), r(1,1); // could maybe also use v1.resize(1, 4); possibly faster
//       return ret;
//     }

    
    // Eigen::Vector4d rstar_xcomp_from_r(Eigen::Matrix2d r)
    // {
    //   Eigen::Vector4d ret;
    //   double a = r(0,0);
    //   double b = r(0,1);
    //   double c = r(1,0);
    //   double d = r(1,1);
    //   ret << a*a, a*b, a*b, b*b;

    //   // ret << r(0,0), r(0, 1), r(1,0), r(1,1); // could maybe also use v1.resize(1, 4); possibly faster
    //   return ret;
    // }







void GN_proj_to_rank_1(Eigen::Matrix2d p, Eigen::Vector2d& v)
{        
    // std::cout << "Newtons method to do rank-1 projection test" << std::endl;

    // Eigen::Matrix2d p = targ*targ.transpose();
        // std::cout << " p " << p << std::endl;
        double a = p(0,0);
        double b = p(1,0);
        double c = p(0,1);
        double d = p(1,1);
        Eigen::Vector2d grad;
        Eigen::Matrix2d hess;
        // Eigen::Vector2d v = Eigen::Vector2d::Random()*10;
        // v = .1*v + targ;

        Eigen::Vector2d v_prev;
        int iter;
        
        for(int i = 0; i < 500; i++)
        {
            // std::cout << v.transpose() << std::endl;
            v_prev = v;

            grad(0) = -4 * v(0) * (a - v(0)*v(0)) - 2 * v(1) * (b - v(0) * v(1)) - 2 * v(1) * (c - v(0) * v(1));
            grad(1) = -2 * v(0) * (b - v(0) * v(1)) - 2 * v(0) * (c - v(0) * v(1)) - 4 * v(1) * (d - v(1)*v(1));
            hess(0,0) = 8 * v(0)*v(0) - 4 * (a - v(0)*v(0)) + 4 * v(1)*v(1); 
            hess(0,1) = 4 * v(0) * v(1) - 2 * (b - v(0) * v(1)) - 2 * (c - v(0) * v(1));
            hess(1,0) = 4 * v(0) * v(1) - 2 * (b - v(0) * v(1)) - 2 * (c - v(0) * v(1));
            hess(1,1) = 4 * v(0)*v(0) + 8 * v(1)*v(1) - 4 * (d - v(1)*v(1));
            
            hess.colPivHouseholderQr().solve(-grad);
            Eigen::Matrix2d inv = hess.inverse();
            v += hess.colPivHouseholderQr().solve(-grad); // .5 * (inv * grad);
            // v += inv*(-grad-hess*v);

            if ( ( v - v_prev ).squaredNorm() < 1e-8 )
            {
                // std::cout << "break at " << i << std::endl;
                break;
            }
      
        }

        std::cout << "p \n " << p << " v \n" << v << v.transpose() << std::endl;
        if(v.norm() < 1e-4)
        {
            v = Eigen::Vector2d::Random();
        }
}