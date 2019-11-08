/*
**    TP CPE Lyon
**    Copyright (C) 2015 Damien Rohmer
**
**    This program is free software: you can redistribute it and/or modify
**    it under the terms of the GNU General Public License as published by
**    the Free Software Foundation, either version 3 of the License, or
**    (at your option) any later version.
**
**   This program is distributed in the hope that it will be useful,
**    but WITHOUT ANY WARRANTY; without even the implied warranty of
**    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**    GNU General Public License for more details.
**
**    You should have received a copy of the GNU General Public License
**    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "mesh_parametric_cloth.hpp"

#include "../lib/common/error_handling.hpp"
#include <cmath>

namespace cpe
{


void mesh_parametric_cloth::update_force()
{

    int const Nu = size_u();
    int const Nv = size_v();
    int const N_total = Nu*Nv;
    ASSERT_CPE(static_cast<int>(force_data.size()) == Nu*Nv , "Error of size");


    //Gravity
    static vec3 const g (0.0f,0.0f,-9.81f);
    vec3 const g_normalized = g/N_total;
    for(int ku=0 ; ku<Nu ; ++ku)
    {
        for(int kv=0 ; kv<Nv ; ++kv)
        {
            force(ku,kv) =  g_normalized;

        }
    }

    //*************************************************************//
    // TO DO, Calculer les forces s'appliquant sur chaque sommet
    //*************************************************************//
    float Lo = 1./(Nu - 1.);
    float L1 = (float) 1.0f * sqrt(2)/(Nu-1.);
    float k = 10;
    float k2 = 8;
    vec3 voisins[12] = {vec3(1,0,Lo),vec3(-1,0,Lo),vec3(0,1,Lo),vec3(0,-1,Lo),
                      vec3(1,1,L1),vec3(-1,1,L1),vec3(1,-1, L1),vec3(-1,-1, L1),
                       vec3(2,0,2*Lo),vec3(-2,0,2*Lo),vec3(0,2,2*Lo),vec3(0,-2,2*Lo)};

    vec3 v;
    vec3 u;
    float K;
    for(int ku=0 ; ku<Nu; ++ku)
    {
        for(int kv=0 ; kv<Nv ; ++kv)
        {
           for(int ik = 0 ; ik < 12; ik++){
               v = voisins[ik];
               if(ku + v.x() < 0 || ku + v.x()>Nu - 1 || kv + v.y()<0 || kv + v.y()>Nv - 1)
                   continue;
               K = k;
               if(ik > 7)
                   K = k2;
               u = vertex(ku, kv) - vertex(ku + v.x(), kv + v.y());
               force(ku, kv) += K * (v.z() - norm(u))*u/norm(u);
           }

         }

    }
    //*************************************************************//

}


void mesh_parametric_cloth::integration_step(float const dt)
{
    ASSERT_CPE(speed_data.size() == force_data.size(),"Incorrect size");
    ASSERT_CPE(static_cast<int>(speed_data.size()) == size_vertex(),"Incorrect size");



    int const Nu = size_u();
    int const Nv = size_v();
    //*************************************************************//
    // TO DO: Calculer l'integration numerique des positions au cours de l'intervalle de temps dt.
    //*************************************************************//
    //
    // ...
    //
    //
    //
    //*************************************************************//



    force (0,0) = vec3(0,0,0);
    force(Nu - 1, 0) = vec3(0,0,0);
    float mu = 0.01f;
    vec3 c = vec3(0.5f,0.05f,-1.1f);
    float r = 0.198;
    vec3 normal = vec3(0.0f,0.0f,1.0f);
    vec3 p0 = vec3(-0.5f,-1.0f,-1.1f);
    float m1 = 0;
    float m2 = 1;
    vec3 test = vec3(0,0,0);
    for(int ku=0 ; ku<Nu  ; ++ku)
    {
        for(int kv=0 ; kv<Nv   ; ++kv)
        {

           speed(ku,kv) = (1 - mu*dt) *speed(ku,kv) + dt * force(ku, kv);
           vertex(ku, kv) += dt * speed(ku, kv);
           if(vertex(ku, kv).y() >= 1.1f){
               speed(ku,kv) = (1 - mu*dt) *speed(ku,kv);
               vertex(ku, kv) += dt * speed(ku, kv);
            }

           if(dot(vertex(ku, kv)- p0, normal) <= 0){
               speed(ku,kv) = -mu*dot(speed(ku,kv),normal) * normal +
                       (speed(ku,kv)- dot(speed(ku,kv),normal)* normal);
               vertex(ku, kv) = vertex(ku,kv) - dot(vertex(ku, kv)- p0, normal) *normal;
           }

           if(norm(vertex(ku, kv) - c) < r + 1 || norm(vertex(ku, kv) - c) > r - 1 )
               test = vertex(ku,kv);
           if(norm(vertex(ku, kv) - c) < r){
               std::cout << "test" << std::endl;
               vec3 u = vertex(ku, kv) - c;
               speed(ku,kv) += 1/(m1+m2)*(m2*dot(vec3(0,0,0),u)-1/2*(m1+3*m2)*dot(speed(ku,kv),u))*u;
               vertex(ku, kv) = vertex(ku,kv) - dot(vertex(ku, kv) - test , u) *u;
           }


        }
    }

    //security check (throw exception if divergence is detected)
    static float const LIMIT=30.0f;
    for(int ku=0 ; ku<Nu ; ++ku)
    {
        for(int kv=0 ; kv<Nv ; ++kv)
        {
            vec3 const& p = vertex(ku,kv);

            if( norm(p) > LIMIT )
            {
                throw exception_divergence("Divergence of the system",EXCEPTION_PARAMETERS_CPE);
            }
        }
    }

}


void mesh_parametric_cloth::set_plane_xy_unit(int const size_u_param,int const size_v_param)
{
    mesh_parametric::set_plane_xy_unit(size_u_param,size_v_param);

    int const N = size_u()*size_v();
    speed_data.resize(N);
    force_data.resize(N);
}

vec3 const& mesh_parametric_cloth::speed(int const ku,int const kv) const
{
    ASSERT_CPE(ku >= 0 , "Value ku ("+std::to_string(ku)+") should be >=0 ");
    ASSERT_CPE(ku < size_u() , "Value ku ("+std::to_string(ku)+") should be < size_u ("+std::to_string(size_u())+")");
    ASSERT_CPE(kv >= 0 , "Value kv ("+std::to_string(kv)+") should be >=0 ");
    ASSERT_CPE(kv < size_v() , "Value kv ("+std::to_string(kv)+") should be < size_v ("+std::to_string(size_v())+")");

    int const offset = ku + size_u()*kv;

    ASSERT_CPE(offset < static_cast<int>(speed_data.size()),"Internal error");

    return speed_data[offset];
}

vec3& mesh_parametric_cloth::speed(int const ku,int const kv)
{
    ASSERT_CPE(ku >= 0 , "Value ku ("+std::to_string(ku)+") should be >=0 ");
    ASSERT_CPE(ku < size_u() , "Value ku ("+std::to_string(ku)+") should be < size_u ("+std::to_string(size_u())+")");
    ASSERT_CPE(kv >= 0 , "Value kv ("+std::to_string(kv)+") should be >=0 ");
    ASSERT_CPE(kv < size_v() , "Value kv ("+std::to_string(kv)+") should be < size_v ("+std::to_string(size_v())+")");

    int const offset = ku + size_u()*kv;

    ASSERT_CPE(offset < static_cast<int>(speed_data.size()),"Internal error");

    return speed_data[offset];
}

vec3 const& mesh_parametric_cloth::force(int const ku,int const kv) const
{
    ASSERT_CPE(ku >= 0 , "Value ku ("+std::to_string(ku)+") should be >=0 ");
    ASSERT_CPE(ku < size_u() , "Value ku ("+std::to_string(ku)+") should be < size_u ("+std::to_string(size_u())+")");
    ASSERT_CPE(kv >= 0 , "Value kv ("+std::to_string(kv)+") should be >=0 ");
    ASSERT_CPE(kv < size_v() , "Value kv ("+std::to_string(kv)+") should be < size_v ("+std::to_string(size_v())+")");

    int const offset = ku + size_u()*kv;

    ASSERT_CPE(offset < static_cast<int>(force_data.size()),"Internal error");

    return force_data[offset];
}

vec3& mesh_parametric_cloth::force(int const ku,int const kv)
{
    ASSERT_CPE(ku >= 0 , "Value ku ("+std::to_string(ku)+") should be >=0 ");
    ASSERT_CPE(ku < size_u() , "Value ku ("+std::to_string(ku)+") should be < size_u ("+std::to_string(size_u())+")");
    ASSERT_CPE(kv >= 0 , "Value kv ("+std::to_string(kv)+") should be >=0 ");
    ASSERT_CPE(kv < size_v() , "Value kv ("+std::to_string(kv)+") should be < size_v ("+std::to_string(size_v())+")");

    int const offset = ku + size_u()*kv;

    ASSERT_CPE(offset < static_cast<int>(force_data.size()),"Internal error");

    return force_data[offset];
}


}
