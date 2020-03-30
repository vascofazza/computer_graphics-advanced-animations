#version 120

attribute vec3 vertex_pos;          // vertex position (in mesh coordinate frame)
attribute vec3 vertex_norm;         // vertex normal   (in mesh coordinate frame)
attribute vec2 vertex_texcoord;     // vertex texture coordinate

uniform mat4 mesh_frame;            // mesh frame (as a matrix)
uniform mat4 camera_frame_inverse;  // inverse of the camera frame (as a matrix)
uniform mat4 camera_projection;     // camera projection

varying vec3 pos;                   // [to fragment shader] vertex position (in world coordinate)
varying vec3 norm;                  // [to fragment shader] vertex normal (in world coordinate)
varying vec2 texcoord;              // [to fragment shader] vertex texture coordinate

uniform bool skin_enabled;                  // skinning
uniform mat4 skin_bone_xforms[48];          // bone xform
attribute vec4 vertex_skin_bone_ids;        // skin bone indices
attribute vec4 vertex_skin_bone_weights;    // skin weights
attribute vec3 rest_pos;

// main function
void main() {
    // apply skinning if necessary
    if(skin_enabled) {
        // YOUR CODE GOES HERE ---------------------
        // (only for extra credit)
        // this is a placeholder
        pos = vec3(0,0,0);
        norm = vec3(0,0,0);
        for(int j = 0; j < 4; j++)
        {
            // get bone weight and index
            int idx = int(vertex_skin_bone_ids[j]);
            float w = vertex_skin_bone_weights[j];
            // if index < 0, continue
            if(idx < 0) continue; //se il bone e' attivo
            // grab bone xform
            mat4 xform = skin_bone_xforms[idx];
            // update position and normal
            vec4 npos = vec4(rest_pos,1); //in realta' non serve rest_pos, non aggiorno veramente le pos
            vec4 tmp =(xform*npos);
            pos+= (tmp.xyz/tmp.w)*w;
            norm += (xform*vec4(vertex_norm,0)).xyz*w;
            //pos += (xform*npos).xyz*w; //transform_point(xform, mesh->skinning->rest_pos[i])*w;
            //norm += transform_normal(xform, mesh->skinning->rest_norm[i])*w;
        }
        norm = normalize(norm);
    } else {
        pos = vertex_pos;
        norm = vertex_norm;
    }
    // compute pos and normal in world space and set up variables for fragment shader (use mesh_frame)
    pos = (mesh_frame * vec4(pos,1)).xyz / (mesh_frame * vec4(pos,1)).w;
    norm = (mesh_frame * vec4(norm,0)).xyz;
    // copy texture coordinates down
    texcoord = vertex_texcoord;
    // project vertex position to gl_Position using mesh_frame, camera_frame_inverse and camera_projection
    gl_Position = camera_projection * camera_frame_inverse * vec4(pos,1);
}
