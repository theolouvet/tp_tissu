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

#pragma once

#ifndef MESH_PARAMETRIC_HPP
#define MESH_PARAMETRIC_HPP

#include "mesh_basic.hpp"
#include <math.h>
#include <vector>


namespace cpe
{

/** Special mesh class dealing with parametric (u,v) coordinates */
class mesh_parametric : public mesh_basic
{
public:

    mesh_parametric();

    /** Init the mesh to a basic plane of size (Nu,Nv) */
    void set_plane_xy_unit(int size_u_param,int size_v_param);

    /** The number of samples in the u direction */
    int size_u() const;
    /** The number of samples in the v direction */
    int size_v() const;

    vec3 vertex(int ku,int kv) const;
    vec3& vertex(int ku,int kv);
    vec3 normal(int ku,int kv) const;
    vec3& normal(int ku,int kv);
    vec3 color(int ku,int kv) const;
    vec3& color(int ku,int kv);
    vec2 texture_coord(int ku,int kv) const;
    vec2& texture_coord(int ku,int kv);

    /** Check if the mesh is valid */
    bool valid_mesh() const;


private:

    /** Internal storage of the u size */
    int size_u_data;
    /** Internal storage of the v size */
    int size_v_data;

};

}


#endif
