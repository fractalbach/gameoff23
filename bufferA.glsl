#define KEY_W 87
#define KEY_A 65
#define KEY_S 83
#define KEY_D 68

#define KEY_LEFT  37
#define KEY_UP    38
#define KEY_RIGHT 39
#define KEY_DOWN  40

const float velocity = 15.0;

const vec2 min_offset_wasd = vec2(-1.0, 0.0);
const vec2 max_offset_wasd = vec2( 1.0, 1.0);

const vec2 min_offset_arrows = vec2(-1.0, -0.2);
const vec2 max_offset_arrows = vec2( 1.0,  1.0);


void handleMemory( out vec4 fragColor, in vec2 fragCoord, float player_scale, bool player_enemy_collision )
{
    ivec2 coord = ivec2(floor(fragCoord));
    
    if ( iTime < 0.5 ) { 
        fragColor = vec4(0.); 
        return; 
    }
    
    // Movement and Scaling using W,A,S,D
    if ( coord == ivec2(0,0) )
    {
        vec2 offset = texelFetch( iChannel0, ivec2(0,0), 0 ).xy;
        vec2 left   = texelFetch( iChannel1, ivec2(KEY_A,0), 0 ).x * vec2(-0.1,  0.0);
        vec2 up     = texelFetch( iChannel1, ivec2(KEY_W,0), 0 ).x * vec2( 0.0,  0.1);
        vec2 right  = texelFetch( iChannel1, ivec2(KEY_D,0), 0 ).x * vec2( 0.1,  0.0);
        vec2 down   = texelFetch( iChannel1, ivec2(KEY_S,0), 0 ).x * vec2( 0.0, -0.1);
        
        offset += (left + up + right + down) * velocity * iTimeDelta;
        offset = clamp( offset, min_offset_wasd, max_offset_wasd );
        
        fragColor.xy = offset;
    }
    
    // Moving the camera around using arrow keys
    else if ( coord == ivec2(0,1) )
    {
        vec2 offset = texelFetch( iChannel0, ivec2(0,1), 0 ).xy;
        vec2 left   = texelFetch( iChannel1, ivec2(KEY_LEFT ,0), 0 ).x * vec2(-0.1,  0.0);
        vec2 up     = texelFetch( iChannel1, ivec2(KEY_UP   ,0), 0 ).x * vec2( 0.0,  0.1);
        vec2 right  = texelFetch( iChannel1, ivec2(KEY_RIGHT,0), 0 ).x * vec2( 0.1,  0.0);
        vec2 down   = texelFetch( iChannel1, ivec2(KEY_DOWN ,0), 0 ).x * vec2( 0.0, -0.1);
        
        offset += (left + up + right + down) * velocity * iTimeDelta;
        offset = clamp( offset, min_offset_arrows, max_offset_arrows );
        
        fragColor.xy = offset;
    }
    
    else if ( coord == MEMORY_START_TIME )
    {
        float start_time = texelFetch( iChannel0, MEMORY_START_TIME, 0 ).x;
        
        if ( player_enemy_collision ) {
            fragColor.x = iTime;
        } else {
            fragColor.x = start_time;
        }
    }
    
    else if ( coord == MEMORY_SCORE )
    {
        float score = texelFetch( iChannel0, MEMORY_SCORE, 0 ).x;
        if ( player_enemy_collision ) {
            score = 0.0;
        }
        
        float start_time = texelFetch( iChannel0, MEMORY_START_TIME, 0 ).x;
        
        float new_score = score += (iTimeDelta / 2.0) * player_scale * 2.0;
        fragColor.x = new_score;
        
    }
    
    else if ( coord == MEMORY_HIGH_SCORE )
    {
        float high_score = texelFetch( iChannel0, MEMORY_HIGH_SCORE, 0 ).x;
        float score = texelFetch( iChannel0, MEMORY_SCORE, 0 ).x;
        fragColor.x = max( high_score, score );
    }
}



//================================================
// HEY, YOU! TRY ADJUSTING THESE CONSTANTS :D
//------------------------------------------------

// Higher numbers --> Higher quality
// Lower numbers  --> Better Performance

const int MAX_SAMPLES = 2; // initial number of rays born into the world, for each pixel.
const int MAX_DEPTH   = 6; // maximum number of reflections or refractions, for each ray.


//================================================
// Structures
//------------------------------------------------

struct ray {
    vec3 origin;    // point in world-space where the ray begins
    vec3 direction; // unit-vector for the direction of the ray
};

struct hit_result {
    // bool exists        // true if the ray actually hit something
    vec3  p;            // point in world-space at which the hit occurred
    vec3  n;            // normal vector (could point either outside or inside)
    float t;            // time at which the ray hit
    bool  front_facing; // normal is pointing outwards
};

struct sphere {
    float radius;
    vec3  center;
    uint  material;
};

struct camera {
    float aspect_ratio;
    int image_width;
    int image_height;
    vec3 center;
    vec3 pixel00_loc;
    vec3 pixel_delta_u;
    vec3 pixel_delta_v;
};

struct material {
    uint type; // enumerated value corresponding to unique material types
    float albedo;
    float attenuation;
    float fuzz;
};

//================================================
// Constants
//------------------------------------------------

const float FLOAT_BIG  = 3.4e38; // not the maximum possible float, but it's up there.
const vec3 ALBEDO_METAL_GOLD = vec3(0.8, 0.6, 0.2);

// Material Enumerations
const uint MAT_LAMBERTIAN = 1u;
const uint MAT_METAL = 2u;
const uint MAT_DIELECTRIC = 3u;
const uint MAT_CHECKER = 4u;

//================================================
// Variables
//------------------------------------------------

// Need to make these variable for the moment so that they can be 
// accessed by the random number generator, but there's gotta be a better way
int current_sample = 0;

// CAMERA
camera cam;


//================================================
// The World
//------------------------------------------------

const int MAX_SPHERES = 17;
const int ENEMY_INDEX_BEGIN = 5;
const int NUM_ENEMIES = 8;

sphere all_spheres[MAX_SPHERES] = sphere[MAX_SPHERES]
(
    sphere( -0.2 , vec3( 0.0, -0.3, 2.0), MAT_DIELECTRIC ),    // the sphere you control
    sphere( 0.2 ,  vec3(1.5, -0.3 , 1.0), MAT_DIELECTRIC ),     // another sphere you may control later (off the right atm)
    sphere( 200.0 , vec3(0., -200.5 , 0. ), MAT_CHECKER ), // the ground
    sphere( 1.0,  vec3(-0.8, 3.0, 20.0), MAT_METAL ), // the sun
    sphere( -1.0, vec3(-0.8, 3.0, 20.0), MAT_DIELECTRIC ), // the moon
    
    // some spheres that will approach us through the tunnel
    sphere( 0.2, vec3(0.3), MAT_LAMBERTIAN ),
    sphere( 0.2, vec3(0.3), MAT_METAL ),
    sphere( 0.2, vec3(0.3), MAT_LAMBERTIAN ),
    sphere( 0.2, vec3(0.3), MAT_METAL ),
    sphere( 0.2, vec3(0.3), MAT_LAMBERTIAN ),
    sphere( 0.2, vec3(0.3), MAT_METAL ),
    sphere( 0.2, vec3(0.3), MAT_LAMBERTIAN ),
    sphere( 0.2, vec3(0.3), MAT_METAL ),
    
    // extras
    sphere( -5.5 , vec3(- 8.0, 2.3 , 15.0 ), MAT_DIELECTRIC ), // the left globe
    sphere( -5.5 , vec3(  8.0, 2.3 , 15.0 ), MAT_DIELECTRIC ), // the right globe
    
    sphere( 2.0 , vec3(-6.0, 1.0 ,  6.0 ), MAT_METAL ), // the left globe
    sphere( 2.0 , vec3( 6.0, 1.0 ,  6.0 ), MAT_METAL ) // the right globe
    
    // sphere( 0.15 , vec3(0.0, 0.0 , 0.0 ), MAT_LAMBERTIAN ), // head
    // sphere( 0.05 , vec3(0.0, 0.3 , 4.0 ), MAT_DIELECTRIC ), // right arm
    // sphere( 0.05 , vec3(0.0, 0.3 , 4.0 ), MAT_DIELECTRIC ) // left arm
);

#define SPHERE_SUN 3
#define SPHERE_MOON 4
#define SPHERE_HEAD 15
#define SPHERE_ARM_LEFT 16
#define SPHERE_ARM_RIGHT 17

//================================================
// Update the World
//------------------------------------------------

void update_world()
{
    float C = cos(iTime/PIE);
    float S = sin(iTime/PIE);
    
    // Control of the Primary Sphere
    // all_spheres[0].center.x = (iMouse.x / iResolution.x - 0.5) * 5.0;
    // all_spheres[0].center.y = (iMouse.y / iResolution.y - 0.5) / 2.0 - 0.15;
    all_spheres[0].center.z = 2.0;
    // all_spheres[0].radius   = ((iMouse.y / iResolution.y) / 2.0 + 0.2);
    
    vec2 input_wasd = texelFetch( iChannel0, MEMORY_WASD, 0 ).xy;
    all_spheres[0].center.x = 2.0 * input_wasd.x;
    all_spheres[0].center.y = ( input_wasd.y ) / 1.0 - 0.21;
    all_spheres[0].radius   = ( input_wasd.y ) / 1.0 + 0.4;
    
    
    // all_spheres[1].center.x = (iMouse.x / iResolution.x - 0.5) * 5.0 + 0.5;
    all_spheres[1].center.x = all_spheres[0].center.x;
    all_spheres[1].center.y = all_spheres[0].center.y + 1.5 * abs(all_spheres[0].radius);
    all_spheres[1].center.z = 2.0;
    all_spheres[1].radius   = all_spheres[0].radius;
    
    // all_spheres[SPHERE_HEAD].center = all_spheres[1].center + vec3( 0.0, all_spheres[0].radius + all_spheres[SPHERE_HEAD].radius - 0.05, 0.0 );

    /*
    all_spheres[3].center.x = S;
    all_spheres[3].center.z = -2.-cos(iTime/PIE);
    all_spheres[1].center.y = 1. + C;
    all_spheres[2].cent er.y = 1. + S;
    */

    // update the sun and moon
    all_spheres[SPHERE_SUN].center.x = 4.0 * cos(iTime);
    all_spheres[SPHERE_SUN].center.y = 4.0 * sin(iTime);
   
    all_spheres[SPHERE_MOON].center.x = 4.0 * cos(iTime + PI);
    all_spheres[SPHERE_MOON].center.y = 4.0 * sin(iTime + PI);
    
    for (int i = ENEMY_INDEX_BEGIN; i < (ENEMY_INDEX_BEGIN + NUM_ENEMIES); i++) {
        float zi = 50.0;
        float zf = -1.0;
        float tf = 2.5;
        float rate = (zf - zi) / tf;
        float newTime = iTime + ( float(i)/float(NUM_ENEMIES) * tf );
        float t = mod(newTime, tf);
        all_spheres[i].center.z = zi + rate * t;
        // all_spheres[i].center.x = ( float(i%5)/5. - 0.5 ) * 6.0;
        all_spheres[i].center.x = ( hash12( uvec2( newTime/tf, i) ) * 2.0 - 1.0 ) * 3.0;
        all_spheres[i].center.y = -0.25;
    }
}


//================================================
// Utils
//------------------------------------------------

// point in world-space where ray r would be at time t
vec3 at( ray r, float t ) {
    return r.origin + t * r.direction;
}

bool near_zero(in vec3 v) {
    return all(lessThan(v,vec3(1e-8)));
}

bool sphere_sphere_collides( sphere s1, sphere s2 ) {
    float rsquared = ( abs(s1.radius) + abs(s2.radius) );
    rsquared *= rsquared;
    vec3 delta = s2.center - s1.center;
    return dot(delta,delta) < rsquared;
}

// RANDOM

// uniformly random on the surface of a sphere
// produces normal vectors as well
// requires 2 randomly generated numbers between (0.0, 1.0)
/*
vec3 uniform_sphere_area(vec2 rnds) {
    vec2 u = rnds;
    float phi = 6.28318530718*u.x;
    float rho_c = 2.0 * u.y - 1.0;
    float rho_s = sqrt(1.0 - (rho_c * rho_c));
    return vec3(rho_s * cos(phi), rho_s * sin(phi), rho_c);
}

    float r1 = rnds[0];
    float r2 = rnds[1];
    float theta = 2. * PI * r1;
    float phi = acos(1. - 2. * r2);
    float x = sin(phi) * cos(theta);
    float y = sin(phi) * sin(theta);
    float z = cos(phi);
    return vec3(x,y,z);
*/
vec3 uniform_sphere_area(int sample_number, int depth, vec2 fragCoord) {
    int i = 0;
    while (true) {
        vec3 p = hash33(uvec3(iFrame+i, sample_number*MAX_SAMPLES + depth, fragCoord.x + iResolution.x * fragCoord.y));
        p = 2. * p - 1.;
        if (dot(p,p) < 1.)
            return p;
    }
}

// returns a randomly generated normal vector
// on the outward surface of a sphere
vec3 random_on_hemisphere(vec3 normal, int sample_number, int depth, vec2 fragCoord) {
    vec3 on_unit_sphere = 2. * normalize(uniform_sphere_area(sample_number, depth, fragCoord)) + 1.;
    if ( dot( on_unit_sphere, normal ) > 0. ) 
        return on_unit_sphere;
    else
        return -on_unit_sphere;
}

//================================================
// Textures
//------------------------------------------------

// color of a checkered texture, given some uv's and a point in world space.
// However.. the uv's are actually ignored in this current implementation.
vec3 color_checkers( vec2 uv, vec3 p )
{
    const float CHECKER_SCALE = 1.0 / 1.32;
    const vec3 CHECKER_COLOR_1 = vec3( 0.1, 0.3, 0.2 );
    const vec3 CHECKER_COLOR_2 = vec3( 0.9, 0.9, 0.9 );
    ivec3 ip = ivec3( floor( CHECKER_SCALE * p ) );
    bool isEven = ( ip.x + ip.y + ip.z ) % 2 == 0;
    return isEven ? CHECKER_COLOR_1 : CHECKER_COLOR_2;
}

vec3 sky_color(ray r) {
    vec3 unit = normalize( r.direction );
    float a = 0.5 * ( unit.y + 1. );
    vec3 color1 = vec3( 1.0, 1.0, 1.0 );
    vec3 color2 = vec3( 0.5, 0.7, 1.0 );
    vec3 result = mix( color1, color2, a );
    
    // special sky experiment! I want to add a starfield eventually.
    // if ( r.direction.y > cos(iTime/PIE) && r.direction.x < (cos(iTime/PIE)) && r.direction.z < sin(iTime/PIE) ) { result = vec3(1.0, 0.5, 0.5); }
    // iChannelResolution.x * r.direction;
    // convert cartesian to spherical coordinates
    
    float theta = acos( r.direction.z );
    float phi = sign( r.direction.y ) * acos( r.direction.x / length(r.direction.xy) );
    
    float theta1 = theta / PI;
    float phi1 = phi / PIE;
    color2 = vec3(cos(theta1 + iTime*2.), 1., 1.);
    float trippy = sin(10.*abs(r.direction.z / r.direction.x) + 5.*iTime) * 0.5 + 0.5;
    
    color1 = vec3(0.9-trippy, 0.5 , 0.7);
    
    result = mix(color1, color2, a);
    
    return color1;
}


//================================================
// Path Tracing
//------------------------------------------------


vec3 tunnel1(in vec2 pt) {
    pt = 3.0*(pt.xy / iResolution.xy - 0.5)*vec2(iResolution.x/iResolution.y,1);
    float rInv = 1./length(pt);
    pt = pt * rInv - vec2(rInv + iTime,0.5);
    return (texture(iChannel0,pt)*rInv/2.).xyz;
}

vec3 tunnel_reverse(in vec2 pt) {
    pt = 3.0*(pt.xy / iResolution.xy - 0.5)*vec2(iResolution.x/iResolution.y,1);
    float rInv = 1./length(pt);
    pt = pt * rInv - vec2(rInv - iTime,0.5);
    return (texture(iChannel0,pt)*rInv/2.).xyz;
}

// Intersection test between a ray and a sphere. If there's a hit, 
// this function will populate the hit_result and then return true.
bool hit_sphere(in ray r, in sphere s, float tmin, float tmax, out hit_result res) {
    vec3 oc = r.origin - s.center;
    float a = dot(r.direction, r.direction);
    float half_b = dot(oc, r.direction);
    float c = dot(oc, oc) - s.radius * s.radius;
    float d = half_b * half_b - a * c;
    if (d <= 0.) { return false; }
    
    float sqrtd = sqrt(d);
    float root = (-half_b - sqrtd) / a;
    if (root <= tmin || tmax <= root) {
        root = (-half_b + sqrtd) / a;
        if (root <= tmin || tmax <= root)
            return false;
    }
    res.t = root;
    res.p = at(r, root);
    res.n = normalize( (res.p - s.center) / s.radius );
    res.front_facing = dot(r.direction, res.n) < 0.;
    res.n = res.front_facing ? res.n : -res.n;
    return true;
}


bool scatter_metal(in ray r_in, in hit_result rec, in vec3 albedo, inout vec3 attenuation, inout ray scattered, vec2 fragCoord) {
    const float fuzz = 0.1; // 0.3; // TODO: make this a material parameter
    vec3 reflected = reflect(normalize(r_in.direction), rec.n);
    
    // vec3 uniform_sphere_area(int sample_number, int depth, vec2 fragCoord) {
    vec3 vec_random = uniform_sphere_area(current_sample, int(fragCoord.x), fragCoord);
    
    scattered = ray(rec.p, reflected + fuzz*vec_random);
    attenuation = albedo;
    return true;
}


bool scatter_dielectric(in ray r_in, in hit_result rec, in vec3 albedo, inout vec3 attenuation, inout ray scattered) {
    attenuation = vec3(1.0, 1.0, 1.0);
    const float ir = 1.5; // index of refraction, TODO: make this a material parameter
    float refraction_ratio = rec.front_facing ? (1.0/ir) : ir;
    vec3 unit_direction = normalize(r_in.direction);

    float cos_theta = min(dot(-unit_direction, rec.n), 1.0);
    float sin_theta = sqrt(1.0 - cos_theta*cos_theta);

    bool cannot_refract = refraction_ratio * sin_theta > 1.0;
    vec3 direction;

    if (cannot_refract)
        direction = reflect(unit_direction, rec.n);
    else
        direction = refract(unit_direction, rec.n, refraction_ratio);

    scattered = ray(rec.p, direction);

    return true;
}


vec3 ray_color( ray r, in vec2 fragCoord, int sample_number, bool player_turns_red ) {

    vec3 final_color = vec3(1., 1., 1.);
    int depth=0;
    float reflectance = 0.8;
    
    for (depth=0; depth<MAX_DEPTH; depth++)
    {
        hit_result best_result;
        sphere sphere_hit;
        int sphere_hit_index;
        float tmin = 0.001;
        float tmax = FLOAT_BIG;
        bool hit_anything = false;
        vec3 col = vec3(0.0, 0.0, 0.0);

        // Check for intersections with spheres
        for (int i=0; i<MAX_SPHERES; i++)
        {
            // if (i==0 && tan(hash12(uvec2(iFrame/12, fragCoord.x + iResolution.x*fragCoord.y))/PIE) < 0.4) { continue; } // HACKY TRANSPARENCY
            // if ( i==0 && ( int(iResolution.x*fragCoord.y)%2==0 ) ) { continue; } // HACKY TRANSPARENCY
            // if ( i==0 && ( int(fragCoord.x)%2 == 0 && int(fragCoord.y)%2 == 0 ) ) { continue; } // HACKY TRANSPARENCY
            sphere s = all_spheres[i];
            hit_result result;
            bool hit = hit_sphere( r, s, tmin, tmax, result );
            if (hit) {
                hit_anything = true;
                tmax = result.t;
                best_result = result;
                sphere_hit = s;
                sphere_hit_index = i;
            }
        }

        if (hit_anything)
        {
            // Calculate where the next ray should go
            float some_random = hash13(uvec3(iFrame, depth*sample_number, fragCoord.y*fragCoord.x));
            float another_random = hash13(uvec3(iFrame, depth*sample_number, fragCoord.y*fragCoord.x));
            
            ray scattered;
            vec3 attenuation;
            
            if (sphere_hit.material == MAT_METAL && scatter_metal(r, best_result, ALBEDO_METAL_GOLD, attenuation, scattered, fragCoord)) {
                col = attenuation;
                r = scattered;
            } 
            else if (sphere_hit.material == MAT_DIELECTRIC && scatter_dielectric(r, best_result, ALBEDO_METAL_GOLD, attenuation, scattered)) {
                col = attenuation;
                r = scattered;
            }
            else { 
                // lambertian is the default right now
                // Determine the color (TODO: do this better)
                if  (sphere_hit.material == MAT_CHECKER ) {
                    col = color_checkers( vec2(0.0), best_result.p );
                } else {
                    col = 0.5 * (normalize(best_result.n) + 1.); // Normal Vector Debug Color
                }
                
                r.direction = best_result.n + random_on_hemisphere(best_result.n, sample_number, depth, fragCoord); // Lambertian
                if ( near_zero(r.direction) ) {
                    r.direction = best_result.n;
                }
                r.origin = best_result.p + 0.01 * r.direction;
            }

            // apply red color if the player sphere has been hit recently
            if ( player_turns_red && (sphere_hit_index == 0 || sphere_hit_index == 1) ) {
                col = vec3( 1.0, 0.0, 0.0 );
            }
            
            // Fully Apply Color?
            final_color *= reflectance * col;
            
            // final_color = final_color + col * pow(reflectance, float(depth+1));
        }
        else
        {
            col = sky_color(r);
            // if ( tan(hash12(uvec2(iFrame/12, fragCoord.x + iResolution.x*fragCoord.y))/PIE) < 0.4 ) { col *= 0.9 ; } //HACKY COLOR CHANGE
            // starfield(col, vec2((r.direction.x+1.)/2.*iResolution.x, (r.direction.y+1.)/2.*iResolution.y));
            final_color *= reflectance * col; // pow(reflectance, float(depth+1));
            // final_color = col;
            break;
        }
    }
    
    if (depth==0) {
        // starfield(final_color, fragCoord);
        // final_color = tunnel1(fragCoord) * ( tan(hash12(uvec2(iFrame/12, fragCoord.x+iResolution.x*fragCoord.y))/PIE)<0.5 ? 0.4 : 1. );
        final_color = sky_color(r);
    } else {
        final_color = final_color * pow(reflectance, float(depth));
    }
    
    return final_color;
}


//================================================
// TEXT
//------------------------------------------------
makeStr1f(printHighScore) _T _o _p _SPA _S _c _o _r _e _COL _SPA _dec(i,1) _end
makeStr1f(printScore) _S _c _o _r _e _COL _SPA _dec(i,1) _end
makeStr1f(printScale) _S _c _a _l _e _COL _SPA _dec(i,2) _end


//================================================
// Main Image
//------------------------------------------------

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{

    // UPDATE WORLD
    update_world();

    // Collision Check
    bool player_enemy_collision = false;
    for ( int i = ENEMY_INDEX_BEGIN; i < (ENEMY_INDEX_BEGIN + NUM_ENEMIES); i++ )
    {
        if ( sphere_sphere_collides( all_spheres[0], all_spheres[i] ) ) {
            player_enemy_collision = true;
            break;
        }
    }
    
    // Some variables used to turn the player red if there has been a collision recently
    float start_time = texelFetch( iChannel0, MEMORY_START_TIME, 0 ).x;
    bool player_turns_red = ( iTime - start_time ) < 0.5;
    
    // Camera
    cam.center  = vec3(  0.0,  0.5, -1.0 ); // Point camera is looking from
    vec3 lookat = vec3(  0.0,  0.5,  2.0 ); // Point camera is looking at

    // Camera Hack for fun
    vec2 input_arrow_keys = texelFetch( iChannel0, ivec2(0,1), 0 ).xy;
    cam.center.xy += 5.0 * input_arrow_keys.xy;
    // lookat.xy = lookat.xy + input_arrow_keys.xy;
    // cam.center.y = max(cam.center.y, 0.0);
        

    vec3 vup = vec3(0.0, 1.0,  0.0);     // Camera-relative "up" direction
    float vfov = 45.;
    float focal_length = length(cam.center - lookat);
    float theta = radians(vfov);
    float h = tan(theta/2.); 

    // Viewport Calculations
    cam.aspect_ratio = iResolution.x / iResolution.y;
    float viewport_height = 2. * h * focal_length;
    float viewport_width = viewport_height * cam.aspect_ratio;
    
    // Calculate the u,v,w unit basis vectors for the camera coordinate frame.
    vec3 w = normalize(cam.center - lookat);
    vec3 u = normalize(cross(vup, w));
    vec3 v = cross(w, u);
    
    // vectors across the viewport
    vec3 viewport_u = viewport_width * -u;
    vec3 viewport_v = viewport_height * v;
    
    // horizontal and vertical delta vectors from pixel to pixel
    cam.pixel_delta_u = viewport_u / iResolution.x;
    cam.pixel_delta_v = viewport_v / iResolution.y;
    
    // calculate location of the upper left pixel
    vec3 viewport_upper_left = cam.center - (focal_length * w) - viewport_u/2. - viewport_v/2.;
    cam.pixel00_loc = viewport_upper_left + 0.5 * (cam.pixel_delta_u + cam.pixel_delta_v);
    
    // AIM AND LOAD THE RAY-CANNON
    vec3 pixel_center = cam.pixel00_loc + (fragCoord.x * cam.pixel_delta_u) + (fragCoord.y * cam.pixel_delta_v);
    vec3 ray_direction = pixel_center - cam.center;
    ray r = ray(cam.center, ray_direction);
    
    // FIRE THE RAY CANNONS!!!!
    vec3 col = vec3(0.,0.,0.);
    for (int i=0; i<MAX_SAMPLES; i++) {
        r.origin = cam.center;
        vec2 offset = hash23( uvec3(iFrame, i, fragCoord.x*fragCoord.y) )/2.;
        vec3 pixel_sample_square = offset.x * cam.pixel_delta_u + offset.y * cam.pixel_delta_v;
        vec3 pixel_sample = pixel_center + pixel_sample_square;
        r.direction = pixel_sample - r.origin;
        col += ray_color(r, fragCoord, i, player_turns_red);
    }
    col = col / float(MAX_SAMPLES);
    
    // col = ray_color(r, fragCoord);
    

    // Time varying pixel color
    // vec3 col = 0.5 + 0.5*cos(iTime+uv.xyx+vec3(0,2,4));
    
    
    // Gamma Correction
    col = sqrt(col);


    // Font on the top
    float player_scale = 1.0 + texelFetch( iChannel0, MEMORY_WASD, 0 ).y;
    float score = texelFetch( iChannel0, MEMORY_SCORE, 0 ).x;
    float high_score = texelFetch( iChannel0, MEMORY_HIGH_SCORE, 0 ).x;

    
    vec2 uv = fragCoord / iResolution.y;
    const float fontsize = 14.0;
    uv *= fontsize;        // Scale font with font_size
    uv.y -= fontsize - 1.; // Start drawing from the top
    // TEXT HIGH SCORE
    col -= vec3( 0.2, 0.2, 0.2 ) * printHighScore( uv + vec2( -0.05, 0.05 ), high_score ); // shadow text
    col += vec3( 0.8, 1.0, 0.8 ) * printHighScore( uv, high_score ); // foreground text
    uv.y++;
    // TEXT SCORE 
    vec3 score_col1 = player_turns_red ? vec3( 1.0, 0.0, 0.0 ) : vec3( 0.8, 1.0, 0.8 );
    col -= vec3( 0.2, 0.2, 0.2 ) * printScore( uv + vec2( -0.05, 0.05 ), score ); // shadow text
    col += score_col1 * printScore( uv, score ); // foreground text
    uv.y++;
    // TEXT PLAYER SCALE
    col -= vec3( 0.2, 0.2, 0.2 ) * printScale( uv + vec2( -0.05, 0.05 ), player_scale ); // shadow text
    col += vec3( 0.8, 1.0, 0.8 ) * printScale( uv, player_scale ); // foreground text

    // Prevent out-of-bounds colors
    col = clamp(col, 0., 1.);

    // Output to screen
    fragColor = vec4( col, 1.0 );
    
    handleMemory( fragColor, fragCoord, player_scale, player_enemy_collision );
    // fragColor = texelFetch( iChannel0, ivec2(fragCoord), 0 );
}



