//============================================================================
// Constants.
//============================================================================
const float PI = 3.1415926536;

const vec3 BACKGROUND_COLOR = vec3( 0.1, 0.2, 0.6 );

// Vertical field-of-view angle of camera. In radians.
const float FOVY = 50.0 * PI / 180.0;

// Use this for avoiding the "epsilon problem" or the shadow acne problem.
const float DEFAULT_TMIN = 10.0e-4;

// Use this for tmax for non-shadow ray intersection test.
const float DEFAULT_TMAX = 10.0e6;

// Equivalent to number of recursion levels (0 means ray-casting only).
// We are using iterations to replace recursions.
const int NUM_ITERATIONS = 2;

// Constants for the scene objects.
const int NUM_LIGHTS = 2;
const int NUM_MATERIALS = 9;
const int NUM_SPHERES = 22;
const int NUM_CUBES = 1;


//============================================================================
// Define new struct types.
//============================================================================
struct Ray_t {
    vec3 o;  // Ray Origin.
    vec3 d;  // Ray Direction. A unit vector.
};

struct Plane_t {
    // The plane equation is Ax + By + Cz + D = 0.
    float A, B, C, D;
    int materialID;
};

struct Sphere_t {
    vec3 center;
    float radius;
    int materialID;
};


struct Cube_t {
    vec3 center;
    int materialID;
    vec3 size;
};


struct Light_t {
    vec3 position;  // Point light 3D position.
    vec3 I_a;       // For Ambient.
    vec3 I_source;  // For Diffuse and Specular.
};

struct Material_t {
    vec3 k_a;   // Ambient coefficient.
    vec3 k_d;   // Diffuse coefficient.
    vec3 k_r;   // Reflected specular coefficient.
    vec3 k_rg;  // Global reflection coefficient.
    float n;    // The specular reflection exponent. Ranges from 0.0 to 128.0.
};

//============================================================================
// Global scene data.
//============================================================================

Plane_t Table;
Sphere_t Sphere[NUM_SPHERES];
Cube_t Cube[NUM_CUBES];
Light_t Light[NUM_LIGHTS];
Material_t Material[NUM_MATERIALS];

/////////////////////////////////////////////////////////////////////////////
// Initializes the scene.
/////////////////////////////////////////////////////////////////////////////

void InitScene() {
    // Snooker Table Plane
    Table.A = 0.0;
    Table.B = 1.0;
    Table.C = 0.0;
    Table.D = 0.0;
    Table.materialID = 8;
    
    // Black Ball
    Sphere[0].center = vec3(-15.0, 0.5, 0.0);
    Sphere[0].radius = 0.5;
    Sphere[0].materialID = 0;
    
    // Red Ball * 15
    for (int i = 1; i <= 15; i++) {
        Sphere[i].radius = 0.5;
        Sphere[i].materialID = 1;
    }
    Sphere[1].center = vec3(-12.0, 0.5, 2.0);
    
    Sphere[2].center = vec3(-12.0, 0.5, 1.0);
    
    Sphere[3].center = vec3(-12.0, 0.5, 0.0);
    
    Sphere[4].center = vec3(-12.0, 0.5, -1.0);
    
    Sphere[5].center = vec3(-12.0, 0.5, -2.0);
    
    Sphere[6].center = vec3(-11.13, 0.5, 1.5);
    
    Sphere[7].center = vec3(-11.14, 0.5, 0.5);
    
    Sphere[8].center = vec3(-11.14, 0.5, -0.5);
    
    Sphere[9].center = vec3(-11.14, 0.5, -1.5);
    
    Sphere[10].center = vec3(-10.26, 0.5, 1.0);
    
    Sphere[11].center = vec3(-10.28, 0.5, 0.0);
    
    Sphere[12].center = vec3(-10.28, 0.5, -1.0);
    
    Sphere[13].center = vec3(-9.39, 0.5, 0.5);
    
    Sphere[14].center = vec3(-9.39, 0.5, -0.5);
    
    Sphere[15].center = vec3(-8.52, 0.5, 0.0);
    
    // Purple Ball
    Sphere[16].center = vec3(-7.52, 0.5, 0.0);
    Sphere[16].radius = 0.5;
    Sphere[16].materialID = 2;
    
    // Blue Ball
    Sphere[17].center = vec3(0.0, 0.5, 0.0);
    Sphere[17].radius = 0.5;
    Sphere[17].materialID = 3;
    
    // Green Ball
    Sphere[18].center = vec3(11.0, 0.5, 2.5);
    Sphere[18].radius = 0.5;
    Sphere[18].materialID = 4;
    
    // Brown Ball
    Sphere[19].center = vec3(11.0, 0.5, 0.0);
    Sphere[19].radius = 0.5;
    Sphere[19].materialID = 5;
    
    // Yellow Ball
    Sphere[20].center = vec3(11.0, 0.5, -2.5);
    Sphere[20].radius = 0.5;
    Sphere[20].materialID = 6;
    
    // White Ball
    Sphere[21].center = vec3(3.0, 0.5, 4.5);
    Sphere[21].radius = 0.5;
    Sphere[21].materialID = 7;
    
    //cube1
    Cube[0].center = vec3(0.0, 0.5, -10.0);
    Cube[0].size = vec3(50.0, 1.0, 1.0);
    Cube[0].materialID = 5;
    
    // Black Plastic Material
    Material[0].k_d = vec3(0.0, 0.0, 0.0);
    Material[0].k_a = 0.2 * Material[0].k_d;
    Material[0].k_r = vec3(1.0, 1.0, 1.0);
    Material[0].k_rg = 0.5 * Material[0].k_r;
    Material[0].n = 128.0;
    
    // Red Plastic Material
    Material[1].k_d = vec3(0.89, 0.09, 0.05);
    Material[1].k_a = 0.2 * Material[1].k_d;
    Material[1].k_r = vec3(1.0, 1.0, 1.0);
    Material[1].k_rg = 0.5 * Material[1].k_r;
    Material[1].n = 128.0;
    
    // Purple Plastic Material
    Material[2].k_d = vec3(0.86, 0.44, 0.84);
    Material[2].k_a = 0.2 * Material[2].k_d;
    Material[2].k_r = vec3(1.0, 1.0, 1.0);
    Material[2].k_rg = 0.5 * Material[2].k_r;
    Material[2].n = 128.0;
    
    // Blue Plastic Material
    Material[3].k_d = vec3(0.0, 0.0, 1.0);
    Material[3].k_a = 0.2 * Material[3].k_d;
    Material[3].k_r = vec3(1.0, 1.0, 1.0);
    Material[3].k_rg = 0.5 * Material[3].k_r;
    Material[3].n = 128.0;
    
    // Green Plastic Material
    Material[4].k_d = vec3(0.13, 0.55, 0.13);
    Material[4].k_a = 0.2 * Material[4].k_d;
    Material[4].k_r = vec3(1.0, 1.0, 1.0);
    Material[4].k_rg = 0.5 * Material[4].k_r;
    Material[4].n = 128.0;
    
    // Brown Plastic Material
    Material[5].k_d = vec3(0.78, 0.38, 0.08);
    Material[5].k_a = 0.2 * Material[5].k_d;
    Material[5].k_r = vec3(1.0, 1.0, 1.0);
    Material[5].k_rg = 0.5 * Material[5].k_r;
    Material[5].n = 128.0;
    
    // Yellow Plastic Material
    Material[6].k_d = vec3(1.0, 1.0, 0.0);
    Material[6].k_a = 0.2 * Material[6].k_d;
    Material[6].k_r = vec3(1.0, 1.0, 1.0);
    Material[6].k_rg = 0.5 * Material[6].k_r;
    Material[6].n = 128.0;
    
    // White Plastic Material
    Material[7].k_d = vec3(1.0, 1.0, 1.0);
    Material[7].k_a = 0.2 * Material[7].k_d;
    Material[7].k_r = vec3(1.0, 1.0, 1.0);
    Material[7].k_rg = 0.5 * Material[7].k_r;
    Material[7].n = 128.0;
    
    // Table Material
    Material[8].k_d = vec3(0.42, 0.56, 0.14);
    Material[8].k_a = 0.2 * Material[7].k_d;
    Material[8].k_r = vec3(1.0, 1.0, 1.0);
    Material[8].k_rg = 0.5 * Material[7].k_r;
    Material[8].n = 128.0;
    
    // Light 0
    Light[0].position = vec3(9.0 * cos(iTime), 15.0, 9.0 * sin(iTime));
    Light[0].I_a = vec3(0.1, 0.1, 0.1);
    Light[0].I_source = vec3(1.0, 1.0, 1.0);
    
    // Light 1
    Light[1].position = vec3(-18.0, 15.0, -9.0);
    Light[1].I_a = vec3(0.1, 0.1, 0.1);
    Light[1].I_source = vec3(1.0, 1.0, 1.0);
}

/////////////////////////////////////////////////////////////////////////////
// Computes intersection between a plane and a ray.
// Returns true if there is an intersection where the ray parameter t is
// between tmin and tmax, otherwise returns false.
// If there is such an intersection, outputs the value of t, the position
// of the intersection (hitPos) and the normal vector at the intersection
// (hitNormal).
/////////////////////////////////////////////////////////////////////////////
bool IntersectPlane( in Plane_t pln, in Ray_t ray, in float tmin, in float tmax,
                     out float t, out vec3 hitPos, out vec3 hitNormal )
{
    vec3 N = vec3( pln.A, pln.B, pln.C );
    float NRd = dot( N, ray.d );
    float NRo = dot( N, ray.o );
    float t0 = (-pln.D - NRo) / NRd;
    if ( t0 < tmin || t0 > tmax ) return false;

    // We have a hit -- output results.
    t = t0;
    hitPos = ray.o + t0 * ray.d;
    hitNormal = normalize( N );
    return true;
}

/////////////////////////////////////////////////////////////////////////////
// Computes intersection between a plane and a ray.
// Returns true if there is an intersection where the ray parameter t is
// between tmin and tmax, otherwise returns false.
/////////////////////////////////////////////////////////////////////////////
bool IntersectPlane( in Plane_t pln, in Ray_t ray, in float tmin, in float tmax )
{
    vec3 N = vec3( pln.A, pln.B, pln.C );
    float NRd = dot( N, ray.d );
    float NRo = dot( N, ray.o );
    float t0 = (-pln.D - NRo) / NRd;
    if ( t0 < tmin || t0 > tmax ) return false;
    return true;
}

/////////////////////////////////////////////////////////////////////////////
// Computes intersection between a sphere and a ray.
// Returns true if there is an intersection where the ray parameter t is
// between tmin and tmax, otherwise returns false.
// If there is one or two such intersections, outputs the value of the
// smaller t, the position of the intersection (hitPos) and the normal
// vector at the intersection (hitNormal).
/////////////////////////////////////////////////////////////////////////////
bool IntersectSphere( in Sphere_t sph, in Ray_t ray, in float tmin, in float tmax,
                      out float t, out vec3 hitPos, out vec3 hitNormal )
{
    vec3 translateRay = ray.o;
    translateRay -= sph.center;
    float a = 1.0;
    float b = 2.0 * dot(ray.d, translateRay);
    float c = dot(translateRay, translateRay) - sph.radius * sph.radius;
    float d = b * b - 4.0 * a * c;
    if (d == 0.0) {
        float t0 = - b / (2.0 * a);
        if ( t0 < tmin || t0 > tmax ) return false;
        t = t0;
        hitPos = ray.o + t * ray.d;
        hitNormal = normalize(translateRay + t * ray.d);
        return true;
    }
    if (d < 0.0) {
        return false;
    }
    if (d > 0.0) {
        float t1 = (-b + sqrt(d)) / (2.0 * a);
        float t2 = (-b - sqrt(d)) / (2.0 * a);
        bool flag1 = true;
        bool flag2 = true;
        if ( t1 < tmin || t1 > tmax ){
            flag1 = false;
        }
        if ( t2 < tmin || t2 > tmax ){
            flag2 = false;
        }
        if (flag1 == false && flag2 == false) {
            return false;
        }
        if (flag1 == false && flag2 == true) {
            t = t2;
            hitPos = ray.o + t * ray.d;
            hitNormal = normalize(translateRay + t * ray.d);
            return true;
        }
        if (flag1 == true && flag2 == false) {
            t = t1;
            hitPos = ray.o + t * ray.d;
            hitNormal = normalize(translateRay + t * ray.d);
            return true;
        }
        if (flag1 == true && flag2 == true) {
            t = min(t1, t2);
            hitPos = ray.o + t * ray.d;
            hitNormal = normalize(translateRay + t * ray.d);
            return true;
        }
    }
}

/////////////////////////////////////////////////////////////////////////////
// Computes intersection between a sphere and a ray.
// Returns true if there is an intersection where the ray parameter t is
// between tmin and tmax, otherwise returns false.
/////////////////////////////////////////////////////////////////////////////
bool IntersectSphere( in Sphere_t sph, in Ray_t ray, in float tmin, in float tmax )
{
    vec3 translateRay = ray.o;
    translateRay -= sph.center;
    float a = 1.0;
    float b = 2.0 * dot(ray.d, translateRay);
    float c = dot(translateRay, translateRay) - sph.radius * sph.radius;
    float d = b * b - 4.0 * a * c;
    if (d == 0.0) {
        float t0 = - b / (2.0 * a);
        if ( t0 < tmin || t0 > tmax ) return false;
        return true;
    }
    if (d < 0.0) {
        return false;
    }
    if (d > 0.0) {
        float t1 = (-b + sqrt(d)) / (2.0 * a);
        float t2 = (-b - sqrt(d)) / (2.0 * a);
        bool flag1 = true;
        bool flag2 = true;
        if ( t1 < tmin || t1 > tmax ){
            flag1 = false;
        }
        if ( t2 < tmin || t2 > tmax ){
            flag2 = false;
        }
        if (flag1 == false && flag2 == false) {
            return false;
        }
        return true;
    }
}

bool IntersectCube( in Cube_t cube, in Ray_t ray, in float tmin, in float tmax,
                      out float t, out vec3 hitPos, out vec3 hitNormal )
{
    vec3 minBound = cube.center - (cube.size / 2.0);
    vec3 maxBound = cube.center + (cube.size / 2.0);

    vec3 t1 = (minBound - ray.o) / ray.d;
    vec3 t2 = (maxBound - ray.o) / ray.d;
    
    vec3 tsmall = min(t1, t2);
    vec3 tbig = max(t1, t2);
    
    float tstart = max(max(tsmall.x, tsmall.y), max(tsmall.y, tsmall.z));
    float tend = min(min(tbig.x, tbig.y), min(tbig.y, tbig.z));
    
    if(tend < tstart || tstart > tmax || tend < tmin)
        return false;
        
    t = tstart < tmin ? tend : tstart;
    hitPos = ray.o + ray.d * t;
    vec3 d = abs(hitPos - cube.center) / (cube.size / 2.0);
    float maxD = max(max(d.x, d.y), d.z);

    if(maxD == d.x)
        hitNormal = vec3(sign(hitPos.x - cube.center.x), 0.0, 0.0);
    else if(maxD == d.y)
        hitNormal = vec3(0.0, sign(hitPos.y - cube.center.y), 0.0);
    else
        hitNormal = vec3(0.0, 0.0, sign(hitPos.z - cube.center.z));
    hitNormal = normalize(hitNormal);

    return true;
}


bool IntersectCube( in Cube_t cube, in Ray_t ray, in float tmin, in float tmax)
{
    vec3 minBound = cube.center - (cube.size / 2.0);
    vec3 maxBound = cube.center + (cube.size / 2.0);

    vec3 t1 = (minBound - ray.o) / ray.d;
    vec3 t2 = (maxBound - ray.o) / ray.d;
    
    vec3 tsmall = min(t1, t2);
    vec3 tbig = max(t1, t2);
    
    float tstart = max(max(tsmall.x, tsmall.y), max(tsmall.y, tsmall.z));
    float tend = min(min(tbig.x, tbig.y), min(tbig.y, tbig.z));
    
    if(tend < tstart || tstart > tmax || tend < tmin)
        return false;

    return true;

}

/////////////////////////////////////////////////////////////////////////////
// Computes (I_a * k_a) + k_shadow * I_source * [ k_d * (N.L) + k_r * (R.V)^n ].
// Input vectors L, N and V are pointing AWAY from surface point.
// Assume all vectors L, N and V are unit vectors.
/////////////////////////////////////////////////////////////////////////////
vec3 PhongLighting( in vec3 L, in vec3 N, in vec3 V, in bool inShadow,
                    in Material_t mat, in Light_t light )
{
    if ( inShadow ) {
        return light.I_a * mat.k_a;
    }
    else {
        vec3 R = reflect( -L, N );
        float N_dot_L = max( 0.0, dot( N, L ) );
        float R_dot_V = max( 0.0, dot( R, V ) );
        float R_dot_V_pow_n = ( R_dot_V == 0.0 )? 0.0 : pow( R_dot_V, mat.n );

        return light.I_a * mat.k_a +
               light.I_source * (mat.k_d * N_dot_L + mat.k_r * R_dot_V_pow_n);
    }
}


/////////////////////////////////////////////////////////////////////////////
// Casts a ray into the scene and returns color computed at the nearest
// intersection point. The color is the sum of light from all light sources,
// each computed using Phong Lighting Model, with consideration of
// whether the interesection point is being shadowed from the light.
// If there is no interesection, returns the background color, and outputs
// hasHit as false.
// If there is intersection, returns the computed color, and outputs
// hasHit as true, the 3D position of the intersection (hitPos), the
// normal vector at the intersection (hitNormal), and the k_rg value
// of the material of the intersected object.
/////////////////////////////////////////////////////////////////////////////
vec3 CastRay( in Ray_t ray,
              out bool hasHit, out vec3 hitPos, out vec3 hitNormal, out vec3 k_rg )
{
    // Find whether and where the ray hits some object.
    // Take the nearest hit point.

    bool hasHitSomething = false;
    float nearest_t = DEFAULT_TMAX;   // The ray parameter t at the nearest hit point.
    vec3 nearest_hitPos;              // 3D position of the nearest hit point.
    vec3 nearest_hitNormal;           // Normal vector at the nearest hit point.
    int nearest_hitMatID;             // MaterialID of the object at the nearest hit point.

    float temp_t;
    vec3 temp_hitPos;
    vec3 temp_hitNormal;
    bool temp_hasHit;
    

    if (IntersectPlane(Table, ray, DEFAULT_TMIN, nearest_t, temp_t, temp_hitPos, temp_hitNormal) == true) {
        temp_hasHit = true;
        if (temp_t < nearest_t) {
            nearest_t = temp_t;
            nearest_hitPos = temp_hitPos;
            nearest_hitNormal = temp_hitNormal;
            nearest_hitMatID = Table.materialID;
            hasHitSomething = temp_hasHit;
        }
    }
    for (int i = 0; i < NUM_SPHERES; i++) {
        if (IntersectSphere(Sphere[i], ray, DEFAULT_TMIN, nearest_t, temp_t, temp_hitPos, temp_hitNormal) == true) {
            temp_hasHit = true;
            if (temp_t < nearest_t) {
                nearest_t = temp_t;
                nearest_hitPos = temp_hitPos;
                nearest_hitNormal = temp_hitNormal;
                nearest_hitMatID = Sphere[i].materialID;
                hasHitSomething = temp_hasHit;
            }
        }
    }
     for (int i = 0; i < NUM_CUBES; i++) {
        if(IntersectCube( Cube[i], ray, DEFAULT_TMIN, DEFAULT_TMAX,
                      temp_t, temp_hitPos, temp_hitNormal ) && temp_t < nearest_t) {
                         hasHitSomething = true;
                         nearest_t = temp_t;
                         nearest_hitPos = temp_hitPos;
                         nearest_hitNormal = temp_hitNormal;
                         nearest_hitMatID = Cube[i].materialID;
                      }
    }


    // One of the output results.
    hasHit = hasHitSomething;
    if ( !hasHitSomething ) return BACKGROUND_COLOR;

    vec3 I_local = vec3( 0.0 );  // Result color will be accumulated here.

    for (int i = 0; i < NUM_LIGHTS; i++) {
        bool inShadow = false;
        Ray_t shadowRay;
        shadowRay.o = nearest_hitPos;
        shadowRay.d = normalize(Light[i].position - shadowRay.o);
        if (IntersectPlane(Table, shadowRay, DEFAULT_TMIN, length(Light[i].position - shadowRay.o)) == true) {
            inShadow = true;
        }
        for (int j = 0; j < NUM_SPHERES; j++) {
            if (IntersectSphere(Sphere[j], shadowRay, DEFAULT_TMIN, length(Light[i].position - shadowRay.o)) == true) {
                inShadow = true;
            }
        }
        for (int j = 0; j < NUM_CUBES; j++) {
            bool inShadow = IntersectCube( Cube[j], shadowRay, DEFAULT_TMIN, length(Light[i].position - shadowRay.o));
            if(inShadow) {
                inShadow = true;
                break;
            }
        }
        I_local += PhongLighting(shadowRay.d, nearest_hitNormal, -ray.d, inShadow, Material[nearest_hitMatID], Light[i]);
    }

    // Populate output results.
    hitPos = nearest_hitPos;
    hitNormal = nearest_hitNormal;
    k_rg = Material[nearest_hitMatID].k_rg;

    return I_local;
}

void CalcMove(in int sph) {
    vec3 holePos = vec3(18.0, 0.5, -9.0);
    vec3 centerPos = Sphere[sph].center;
    vec3 TargetMove = holePos - centerPos;
    float TargetMovex = abs(TargetMove.x);
    float TargetMovez = abs(TargetMove.z);
    float TargetMoveLength = sqrt(TargetMovex * TargetMovex + TargetMovez * TargetMovez);
    float S = TargetMovez / TargetMoveLength;
    float C = TargetMovex / TargetMoveLength;
    vec3 arrivePos = vec3(centerPos.x - 1.0 * C, 0.5, centerPos.z + 1.0 * S);
    vec3 HitMove = arrivePos - Sphere[21].center;
    float HitMovex = abs(HitMove.x);
    float HitMovez = abs(HitMove.z);
    float HitMoveLength = sqrt(HitMovex * HitMovex + HitMovez * HitMovez);
    float s = HitMovez / HitMoveLength;
    float c = HitMovex / HitMoveLength;
    float T = 3.0;
    float t = mod(iTime, T);
    float Speed = (HitMoveLength + 2.0 * TargetMoveLength) / T;
    float t1 = HitMoveLength / Speed;
    float xDir = 1.0;
    if (Sphere[sph].center.x < Sphere[21].center.x) {
        xDir = -1.0;
    }
    float zDir = 1.0;
    if (Sphere[sph].center.z < Sphere[21].center.z) {
        zDir = -1.0;
    }
    if (t <= t1) {
        
        Sphere[21].center.x += xDir * Speed * t * c;
        Sphere[21].center.z += zDir * Speed * t * s;
    }
    else {
        Sphere[21].center = arrivePos;
        float delta = t - t1;
        Sphere[sph].center.x += Speed * 0.5 * delta * C;
        Sphere[sph].center.z -= Speed * 0.5 * delta * S;
        vec3 dir = HitMove - 0.5 * TargetMove;
        float dirLength = sqrt(dir.x * dir.x + dir.z * dir.z);
        float ss = dir.z / dirLength;
        float cc = dir.x / dirLength;
        Sphere[21].center.x += Speed * 0.15 * delta * cc;
        Sphere[21].center.z += Speed * 0.15 * delta * ss;
    }
}

/////////////////////////////////////////////////////////////////////////////
// Execution of fragment shader starts here.
// 1. Initializes the scene.
// 2. Compute a primary ray for the current pixel (fragment).
// 3. Trace ray into the scene with NUM_ITERATIONS recursion levels.
/////////////////////////////////////////////////////////////////////////////
void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    InitScene();
    
    float T = 3.0;
    int turn = int(mod(iTime / T, 4.0));
    
    CalcMove(17 + turn);

    // Scale pixel 2D position such that its y coordinate is in [-1.0, 1.0].
    vec2 pixel_pos = (2.0 * fragCoord.xy - iResolution.xy) / iResolution.y;

    // Position the camera.
    vec3 cam_pos = vec3( Sphere[21].center.x - 20.0, 10.0, Sphere[21].center.z + 15.0 );
    vec3 cam_lookat = Sphere[21].center;
    vec3 cam_up_vec = vec3( 0.0, 1.0, 0.0 );

    // Set up camera coordinate frame in world space.
    vec3 cam_z_axis = normalize( cam_pos - cam_lookat );
    vec3 cam_x_axis = normalize( cross(cam_up_vec, cam_z_axis) );
    vec3 cam_y_axis = normalize( cross(cam_z_axis, cam_x_axis));

    // Create primary ray.
    float pixel_pos_z = -1.0 / tan(FOVY / 2.0);
    Ray_t pRay;
    pRay.o = cam_pos;
    pRay.d = normalize( pixel_pos.x * cam_x_axis  +  pixel_pos.y * cam_y_axis  +  pixel_pos_z * cam_z_axis );


    // Start Ray Tracing.
    // Use iterations to emulate the recursion.

    vec3 I_result = vec3( 0.0 );
    vec3 compounded_k_rg = vec3( 1.0 );
    Ray_t nextRay = pRay;

    for ( int level = 0; level <= NUM_ITERATIONS; level++ )
    {
        bool hasHit;
        vec3 hitPos, hitNormal, k_rg;

        vec3 I_local = CastRay( nextRay, hasHit, hitPos, hitNormal, k_rg );

        I_result += compounded_k_rg * I_local;

        if ( !hasHit ) break;

        compounded_k_rg *= k_rg;

        nextRay = Ray_t( hitPos, normalize( reflect(nextRay.d, hitNormal) ) );
    }

    fragColor = vec4( I_result, 1.0 );
}