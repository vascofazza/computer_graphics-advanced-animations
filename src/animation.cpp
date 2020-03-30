#include "animation.h"
#include "tesselation.h"
#include "tetgen.h"
#include <tuple>
#include <unordered_set>
#include <unordered_map>
#include "ThreadPool.h"

// compute the frame from an animation
frame3f animate_compute_frame(FrameAnimation* animation, int time) {
    // grab keyframe interval
    auto interval = 0;
    for(auto t : animation->keytimes) if(time < t) break; else interval++;
    interval--;
    // get translation and rotation matrices
    auto t = float(time-animation->keytimes[interval])/float(animation->keytimes[interval+1]-animation->keytimes[interval]);
    auto m_t = translation_matrix(animation->translation[interval]*(1-t)+animation->translation[interval+1]*t);
    auto m_rz = rotation_matrix(animation->rotation[interval].z*(1-t)+animation->rotation[interval+1].z*t,z3f);
    auto m_ry = rotation_matrix(animation->rotation[interval].y*(1-t)+animation->rotation[interval+1].y*t,y3f);
    auto m_rx = rotation_matrix(animation->rotation[interval].x*(1-t)+animation->rotation[interval+1].x*t,x3f);
    // compute combined xform matrix
    auto m = m_t * m_rz * m_ry * m_rx;
    // return the transformed frame
    return transform_frame(m, animation->rest_frame);
}

// update mesh frames for animation
void animate_frame(Scene* scene) {
    // YOUR CODE GOES HERE ---------------------
    // foreach mesh
    for(auto mesh : scene->meshes)
    {
        // if not animation, continue
        if(not mesh->animation)
            continue;
        // update frame
        mesh->frame = animate_compute_frame(mesh->animation, scene->animation->time);
    }
    // foreach surface
    for(auto surf : scene->surfaces){
        // if not animation, continue
        if(not surf->animation)
            continue;
        // update frame
        surf->frame = animate_compute_frame(surf->animation, scene->animation->time);
        // update the _display_mesh
        surf->_display_mesh->frame = surf->frame;
    }
}

// skinning scene
void animate_skin(Scene* scene) {
    // YOUR CODE GOES HERE ---------------------
    // foreach mesh
    for(auto mesh : scene->meshes)
    {
        // if no skinning, continue
        if(not mesh->skinning) continue;
        // foreach vertex index
        for(int i = 0; i < mesh->pos.size(); i++) //ogni vertice e' assegnato ad un bone o piu a seconda dell idx
        {
            // set pos/norm to zero
            mesh->pos[i] = zero3f;
            mesh->norm[i] = zero3f;
            // for each bone slot (0..3)
            for(int j = 0; j < 4; j++)
            {
                // get bone weight and index
                auto idx = mesh->skinning->bone_ids[i][j];
                auto w = mesh->skinning->bone_weights[i][j];
                // if index < 0, continue
                if(idx < 0) continue; //se il bone e' attivo
                // grab bone xform
                auto xform = mesh->skinning->bone_xforms[scene->animation->time][idx];
                // update position and normal
                mesh->pos[i] += transform_point(xform, mesh->skinning->rest_pos[i])*w;
                mesh->norm[i] += transform_normal(xform, mesh->skinning->rest_norm[i])*w;
            }
            // normalize normal
            mesh->norm[i] = normalize(mesh->norm[i]);
        }
    }
}


const static float thick = 0.05f;

static std::unordered_multimap<vec3i, int> multimap;

#define DEBUG_ false

vec3i s_hash(vec3f pos)
{
    int x = pos.x/thick*2;
    int y = pos.y/thick*2;
    int z = pos.z/thick*2;
    return vec3i(x,y,z);
}

void to_triangle_mesh(Mesh* mesh)
{
    multimap.clear();
    if(mesh->triangle.empty())
    {
        mesh->triangle.clear();
        for(auto q : mesh->quad)
        {
            vec3i t1 = {q.x,q.y,q.z};
            vec3i t2 = {q.z,q.w,q.x};
            mesh->triangle.push_back(t1);
            mesh->triangle.push_back(t2);
        }
    }
    for(int t = 0; t < mesh->triangle.size(); t++)
    {
        for(int i = 0; i < 3; i++)
        {
            vec3i key = s_hash(mesh->pos[mesh->triangle[t][i]]);
//            auto key = std::make_tuple(t1.x,t1.y,t1.z);
            auto ret = multimap.equal_range(key);
            for (auto it=ret.first; it!=ret.second; ++it)
            {
                if(it->second == t)
                    continue;
            }
            multimap.insert(std::pair<vec3i, int>(key, t));
        }
    }
}
static std::unordered_set<int> collisions;
//ThreadPool tp(10);
std::pair<vec3i, vec3f> self_collision(Mesh* mesh, int point, vec3f pos)
{
    float T = thick;
    static const float sT = sqrt(thick);
    collisions.clear();
    vector<vec3i> visited;
//    for(int i = -1; i < 2; i++)
//        for(int j = -1; j < 2; j++)
//            for(int k = -1; k < 2; k++)
            {
//                if(abs(i)+abs(j)+abs(k) > 1)
//                    T = sT;
//                else
//                    T = thick;
                //Se mi muovo in diagonale
                auto key = s_hash(pos);//+ vec3f(T*i, T*j, T*k));
                for(auto last_key : visited)
                    if(last_key == key)
                        continue;
                visited.push_back(key);
                auto ret = multimap.equal_range(key);
                for (auto it=ret.first; it!=ret.second; ++it)
                {
                    collisions.insert(it->second);
                }
            }
    std::unordered_set<int>::iterator it;
    for (it=collisions.begin(); it!=collisions.end(); ++it)
    {
        auto tri = mesh->triangle[(*it)];
        if(point == tri.x || point == tri.y || point == tri.z)
            continue;
        //controlliamo la thickness
        auto x4 = pos;//mesh->pos[point];
        auto x43 = x4-mesh->pos[tri.z];
        // compute face normal
        auto u = mesh->pos[tri.z]-mesh->pos[tri.x];
        auto v = mesh->pos[tri.y]-mesh->pos[tri.x];
        auto n = cross(u,v);
        auto dd = abs(dot(x43,n));
        //            message("%f\n",dd);
        if(dd < thick)
        {
            double oneOver4ASquared= 1.0 / dot( n, n );
            vec3f w = x4 - mesh->pos[tri.x];
            //                message("OK\n");
            //project point to triangle
            auto w1 = dot(cross(u,w),n)*oneOver4ASquared;//dot(cross(mesh->pos[tri.y]-mesh->pos[tri.x], mesh->pos[point]-mesh->pos[tri.x]), n)/(dot(n,n));//dot((mesh->pos[tri.x]- mesh->pos[tri.z]), x43);
            auto w2 = dot(cross(w,v),n)*oneOver4ASquared;//dot(cross(mesh->pos[point]-mesh->pos[tri.x],mesh->pos[tri.z]-mesh->pos[tri.x]), n)/(dot(n,n));//dot((mesh->pos[tri.y]- mesh->pos[tri.z]), x43);
            auto w3 = 1-w2-w1;
            if(w1+w2+w3 <= 1 && w1 > 0 && w2 > 0 && w3 > 0 && w1 < 1 && w2 < 1 && w3 < 1)
                return std::pair<vec3i, vec3f>(tri, vec3f(w1,w2,w3));
        }
    }
    return std::pair<vec3i, vec3f>(zero3i, zero3f);
}

// particle simulation
void simulate(Scene* scene) {
    // YOUR CODE GOES HERE ---------------------
    // for each mesh
    auto dt = ((float)scene->animation->dt)/scene->animation->simsteps;
    for(auto mesh : scene->meshes)
    {
        // skip if no simulation
        if(not mesh->simulation) continue;
        // compute time per step
        //
        // foreach simulation steps
        for(auto step : range(scene->animation->simsteps))
        {
                // compute extenal forces (gravity)
            for(int i : range(mesh->simulation->force.size()))
            {
                auto gravity = mesh->simulation->mass[i]*scene->animation->gravity;
                mesh->simulation->force[i] = gravity; //non accumulare, ripartiamo ad ogni frame
            }
                
                // for each spring, compute spring force on points
            for(auto spring : mesh->simulation->springs)
            {
                    // compute spring distance and length
                int p1 = spring.ids.x;
                int p2 = spring.ids.y;
                float ls = length(mesh->pos[p2]-mesh->pos[p1]);
                vec3f ds = normalize(mesh->pos[p2]-mesh->pos[p1]);
                vec3f vs = mesh->simulation->vel[p2]-mesh->simulation->vel[p1];
                    // compute static force
                vec3f f_1_static = spring.ks*(ls-spring.restlength)*ds;
                    // accumulate static force on points
                mesh->simulation->force[p1] += f_1_static;
                mesh->simulation->force[p2] += -f_1_static;
                    // compute dynamic force
                vec3f f_1_dynamic = spring.kd*(dot(vs,ds))*ds;
                    // accumulate dynamic force on points
                mesh->simulation->force[p1] += f_1_dynamic;
                mesh->simulation->force[p2] += -f_1_dynamic;
            }
                // newton laws
            for(int i : range(mesh->simulation->init_pos.size()))
            {
                    // if pinned, skip
                if(mesh->simulation->pinned[i]) continue;
                    // acceleration
                auto a = mesh->simulation->force[i]/mesh->simulation->mass[i];
                    // update velocity and positions using Euler's method
                auto new_vel = mesh->simulation->vel[i] + a*dt;//viene inizializzata nel reset
                auto new_pos = mesh->pos[i] + mesh->simulation->vel[i]*dt + ((a*powf(dt,2))/2.f); //la inizializzazione viene fatta nel reset
                std::pair<vec3i, vec3f> tri = mesh->self_collision? self_collision(mesh, i, new_pos) :pair<vec3i, vec3f>(zero3i, zero3f);
                if(tri.first == zero3i)
                {
                    mesh->simulation->vel[i] = new_vel;
                    mesh->pos[i] = new_pos;
#if DEBUG_
                    mesh->texcoord[i] = {0,0};
#endif
                }
                else
                {
                    
                    //e se fossero piccole sfere?
                    
                    //centro sfera
                    auto c = mesh->pos[tri.first[0]]*tri.second[0] + mesh->pos[tri.first[1]]*tri.second[1] + mesh->pos[tri.first[2]]*tri.second[2];
                    //normal
                    auto n = normalize(new_pos-c);
                    
                    //punto di intersezione
                    auto p = 2*thick*n+c;
                    
                    //posizione invariata;
                    mesh->norm[i] = n;
                    auto v = new_vel;
                    auto vp = v-(dot(n,v))*n;
                    auto vo = -(dot(n,v))*n;
                    auto d_p = scene->animation->bounce_dump.x;
                    auto d_o = scene->animation->bounce_dump.y;
                    mesh->simulation->vel[i] = vp*(1-d_p)+vo*(1-d_o);
                    mesh->pos[i] = new_pos;
#if DEBUG_
                    mesh->texcoord[i] = {1,0};
#endif
                    
                    /*
                    auto dir = normalize(mesh->pos[i]-new_pos);//one3f;//-normalize(mesh->simulation->vel[i]);
                    auto Ic = mesh->simulation->mass[i]*dir/2.f;
                    auto I = (2*Ic)/(1+tri.second[0]*tri.second[0] + tri.second[1]*tri.second[1] + tri.second[2] * tri.second[2]);
                    for(int t = 0; t < 3; t++)
                    {
                        mesh->simulation->vel[tri.first[t]] += tri.second[t]*(I/mesh->simulation->mass[tri.first[t]])*dir;
//                        mesh->pos[tri.first[t]] += mesh->pos[tri.first[t]] + mesh->simulation->vel[tri.first[t]]*dt + ((a*powf(dt,2))/2.f);
                    }
//                     message("LOL\n");
                    mesh->simulation->vel[i] -= (I/mesh->simulation->mass[i])*dir; //dir?
//                    mesh->pos[i] = mesh->pos[i] + mesh->simulation->vel[i]*dt + ((a*powf(dt,2))/2.f);
                    mesh->pos[i] = new_pos;
                    */
                }
                for(auto o_mesh : scene->surfaces)
                {
                    // compute inside tests
                    bool inside = false;
                    auto p = mesh->pos[i];
                    auto n = zero3f;
                    // if quad
                    if(o_mesh->isquad)
                    {
                        // compute local poisition
                        auto p_local = transform_point_inverse(o_mesh->frame, p); //il punto dal punto di vista del piano
                        // perform inside test
                        if(p_local.z < 0 && -o_mesh->radius < p_local.x && p_local.x < o_mesh->radius &&
                           -o_mesh->radius < p_local.y && p_local.y < o_mesh->radius)
                        {
                            inside = true;
                            // if inside, set position and normal
                            n = o_mesh->frame.z;
                            p = transform_point(o_mesh->frame, vec3f(p_local.x, p_local.y, 0.f));
                        }
                    }
                    // else sphere
                    else
                    {
                        // inside test
                        if(length(p-o_mesh->frame.o) < o_mesh->radius)
                        {
                            // if inside, set position and normal
                            inside = true;
                            p = o_mesh->radius*normalize(p-o_mesh->frame.o)+o_mesh->frame.o;
                            n = normalize(p-o_mesh->frame.o);
                        }
                    }
                    // if inside
                    if(inside)
                    {
                        // set particle position
                        mesh->pos[i] = p;
                        mesh->norm[i] = n;
                        // update velocity
                        auto v = mesh->simulation->vel[i];
                        auto vp = v-(dot(n,v))*n;
                        auto vo = -(dot(n,v))*n;
                        auto d_p = scene->animation->bounce_dump.x;
                        auto d_o = scene->animation->bounce_dump.y;
                        mesh->simulation->vel[i] = mesh->friction ? zero3f : vp*(1-d_p)+vo*(1-d_o);
                    }
                }
            }
        }
        // smooth normals if it has triangles or quads
        if(mesh->triangle.size() || mesh->quad.size())
        {
            smooth_normals(mesh);
        }
    }
}

// scene reset
void animate_reset(Scene* scene) {
    scene->animation->time = 0;
    for(auto mesh : scene->meshes) {
        if(mesh->animation) {
            mesh->frame = mesh->animation->rest_frame;
        }
        if(mesh->skinning) {
            mesh->pos = mesh->skinning->rest_pos;
            mesh->norm = mesh->skinning->rest_norm;
        }
        if(mesh->simulation) {
            mesh->pos = mesh->simulation->init_pos;
            mesh->simulation->vel = mesh->simulation->init_vel;
            mesh->simulation->force.resize(mesh->simulation->init_pos.size());
        }
    }
}

void Tetrahedralize(Mesh* mesh)
{
    tetgenio in, out;
    tetgenio::facet *f;
    tetgenio::polygon *p;
    
    // All indices start from 1.
    in.firstnumber = 0;
    
    in.numberofpoints = mesh->pos.size();
    in.pointlist = new REAL[in.numberofpoints * 3];
    for (int i = 0; i < in.numberofpoints; i++) {
        in.pointlist[i * 3]     = mesh->pos[i].x;
        in.pointlist[i * 3 + 1] = mesh->pos[i].y;
        in.pointlist[i * 3 + 2] = mesh->pos[i].z;
    }
    
    for(auto q : mesh->quad)
    {
        vec3i t1 = {q.x,q.y,q.z};
        vec3i t2 = {q.z,q.w,q.x};
        mesh->triangle.push_back(t1);
        mesh->triangle.push_back(t2);
    }
    
    in.numberoffacets = mesh->triangle.size();
    in.facetlist = new tetgenio::facet[in.numberoffacets];
    in.facetmarkerlist = new int[in.numberoffacets];
    
    // Facet 1. The leftmost facet.
    for(int i = 0 ; i < in.numberoffacets; i++)
    {
        f = &in.facetlist[i];
        f->numberofpolygons = 1;
        f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
        f->numberofholes = 0;
        f->holelist = NULL;
        p = &f->polygonlist[0];
        p->numberofvertices = 3;
        p->vertexlist = new int[p->numberofvertices];
        p->vertexlist[0] = mesh->triangle[i].x;
        p->vertexlist[1] = mesh->triangle[i].y;
        p->vertexlist[2] = mesh->triangle[i].z;
//        p->vertexlist[3] = mesh->quad[i].w;
    }
    
    // Output the PLC to files 'barin.node' and 'barin.poly'.
//    in.save_nodes("/tmp/barin");
//    in.save_poly("/tmp/barin");
    
    // Tetrahedralize the PLC. Switches are chosen to read a PLC (p),
    //   do quality mesh generation (q) with a specified quality bound
    //   (1.414), and apply a maximum volume constraint (a0.1).
    
    tetrahedralize("pq1.414a0.001", &in, &out);
    
    // Output mesh to files 'barout.node', 'barout.ele' and 'barout.face'.
//    out.save_nodes("/tmp/barout");
//    out.save_elements("/tmp/barout");
//    out.save_faces("/tmp/barout");
    
    mesh->pos.clear();
    vector<vec4i> old_quad = mesh->quad;//.clear();
    mesh->quad.clear();
    mesh->triangle.clear();
    //mesh->norm.clear();
    mesh->line.clear();
    
    mesh->simulation = new MeshSimulation();
    
    for(int i = 0 ; i < out.numberofpoints*3; i+=3)
    {
        auto p = vec3f(out.pointlist[i+0],out.pointlist[i+1],out.pointlist[i+2]);
        mesh->pos.push_back(p);
        mesh->simulation->init_pos.push_back(p);
        mesh->simulation->init_vel.push_back(zero3f);
        mesh->simulation->mass.push_back(0.005f);
        mesh->simulation->pinned.push_back(0);
    }
    for(int i = 0; i < out.numberoftrifaces*3; i +=3)
    {
//        mesh->line.push_back(vec2i(out.trifacelist[i+0], out.trifacelist[i+1]));
//        mesh->line.push_back(vec2i(out.trifacelist[i+1], out.trifacelist[i+2]));
//        mesh->line.push_back(vec2i(out.trifacelist[i+2], out.trifacelist[i+0]));
       mesh->triangle.push_back(vec3i(out.trifacelist[i+0],out.trifacelist[i+1],out.trifacelist[i+2]));
    }
    for(int i = 0; i < out.numberoftetrahedra*4; i +=4)
    {
//        mesh->line.push_back(vec2i(out.tetrahedronlist[i+0], out.tetrahedronlist[i+1]));
//        mesh->line.push_back(vec2i(out.tetrahedronlist[i+1], out.tetrahedronlist[i+2]));
//        mesh->line.push_back(vec2i(out.tetrahedronlist[i+2], out.tetrahedronlist[i+3]));
//        mesh->line.push_back(vec2i(out.tetrahedronlist[i+3], out.tetrahedronlist[i+0]));
        mesh->quad.push_back(vec4i(out.tetrahedronlist[i+0],out.tetrahedronlist[i+1],out.tetrahedronlist[i+2],out.tetrahedronlist[i+3]));
    }
    EdgeMap map = EdgeMap(mesh->triangle, mesh->quad);
    for(auto edge : map.edges())
    {
        mesh->simulation->springs.push_back({edge, length(mesh->pos[edge.x]-mesh->pos[edge.y]),400,1});
    }
    smooth_normals(mesh);
    //mesh->triangle.clear();
    mesh->quad.clear();
    //mesh->quad = old_quad;
}

Mesh* generate_cloth()
{
    Mesh* mesh = new Mesh();
    mesh->self_collision = true;
    mesh->mat->double_sided = true;
    mesh->simulation = new MeshSimulation();
    int step = 4;
    int size = 100;
    int round = size*2/step;
    for(int i = -size; i < size; i+= step)
        for(int j = -size; j < size; j+= step)
        {
            vec3f p = vec3f(i/(float)size,0.03f,j/(float)size);
            mesh->pos.push_back(p);
            int current = mesh->pos.size()-1;
            mesh->texcoord.push_back(zero2f);
            mesh->simulation->force.push_back(zero3f);
            mesh->simulation->init_pos.push_back(p);
            mesh->simulation->init_vel.push_back(zero3f);
            mesh->simulation->vel.push_back(zero3f);
            mesh->simulation->mass.push_back(0.0004f);
            mesh->simulation->pinned.push_back(false);
            if(j > -size && i > -size)
            {
                mesh->triangle.push_back(vec3i(current, current-round, current-1));
                mesh->triangle.push_back(vec3i(current-round, current-round-1, current-1));
                vec2i edge = vec2i(current, current-round-1);
                mesh->simulation->springs.push_back({edge, length(mesh->pos[edge.x]-mesh->pos[edge.y]),100,10});
            }
//            else if(j == 99 && i < 99)
//                mesh->triangle.push_back(vec3i(current, current+39, current+40));
//            else if(i == 99 && j < 99)
//                mesh->triangle.push_back(vec3i(current, current+1, current-39));
        }
    EdgeMap map = EdgeMap(mesh->triangle, mesh->quad);
    for(auto edge : map.edges())
    {
        mesh->simulation->springs.push_back({edge, length(mesh->pos[edge.x]-mesh->pos[edge.y]),100,10});
    }
    smooth_normals(mesh);
    return mesh;
}

bool first = true;
// scene update
void animate_update(Scene* scene) {
    if(first && scene->meshes[0]->self_collision){scene->meshes.clear(); scene->meshes.push_back(generate_cloth()); first = !first;}
    if(scene->meshes[0]->self_collision)to_triangle_mesh(scene->meshes[0]);
    if(scene->animation->time >= scene->animation->length-1) {
        if(scene->animation->loop) animate_reset(scene);
        else return;
    } else scene->animation->time ++;
    animate_frame(scene);
    if(not scene->animation->gpu_skinning) animate_skin(scene);
    simulate(scene);
}
