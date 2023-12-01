//================================================
// Constants
//------------------------------------------------
const float PI = 3.18309886183790671537767526745028724;
const float PIE = 6.36619772367581343075535053490057448;

#define MEMORY_WASD       ivec2(0,0)
#define MEMORY_ARROWS     ivec2(0,1)
#define MEMORY_START_TIME ivec2(1,0)
#define MEMORY_SCORE      ivec2(1,1)
#define MEMORY_HIGH_SCORE ivec2(2,0)

//================================================
// Some Intersections
//------------------------------------------------

// plane degined by p (p.xyz must be normalized)
float plaIntersect( in vec3 ro, in vec3 rd, in vec4 p )
{
    return -(dot(ro,p.xyz)+p.w)/dot(rd,p.xyz);
}


//================================================
// Starfield
//------------------------------------------------
// "SmallStars (138 chars)"
// by P_Malin
// https://www.shadertoy.com/view/Ml2XDt
// copied on 11/18/2023
//
void starfield( out vec4 f, vec2 p )
{
    p=p/2e3-.2;
    float b = ceil(atan(p.x, p.y) * 6e2), h = cos(b), z = h / dot(p,p);
    f = exp(fract(z + h * b + iDate.wwww) * -1e2) / z;
}
void starfield( out vec3 f, vec2 p )
{
    p=p/2e3-.2;
    float b = ceil(atan(p.x, p.y) * 6e2), h = cos(b), z = h / dot(p,p);
    f = exp(fract(z + h * b + iDate.www) * -1e2) / z;
}


//================================================
// Slit scan tunnel
//------------------------------------------------
// by roywig
// copied on 11/18/2023
// Inspired by physical slit-scan photography, where you dolly the
// camera backwards while you drag the image across a slit. This "slit" is circular.

// Sound borrowed from https://www.shadertoy.com/view/MdfXWX



//================================================
// Random Number Generator
//------------------------------------------------
//Quality hashes collection
//by nimitz 2018 (twitter: @stormoid)

//The MIT License
//Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#if 1
//Modified from: iq's "Integer Hash - III" (https://www.shadertoy.com/view/4tXyWN)
uint baseHash(uvec3 p)
{
    p = 1103515245U*((p.xyz >> 1U)^(p.yzx));
    uint h32 = 1103515245U*((p.x^p.z)^(p.y>>3U));
    return h32^(h32 >> 16);
}

//Modified from: iq's "Integer Hash - III" (https://www.shadertoy.com/view/4tXyWN)
//Faster than "full" xxHash and good quality
uint baseHash(uvec2 p)
{
    p = 1103515245U*((p >> 1U)^(p.yx));
    uint h32 = 1103515245U*((p.x)^(p.y>>3U));
    return h32^(h32 >> 16);
}

uint baseHash(uint p)
{
    p = 1103515245U*((p >> 1U)^(p));
    uint h32 = 1103515245U*((p)^(p>>3U));
    return h32^(h32 >> 16);
}
#else
//XXHash32 based (https://github.com/Cyan4973/xxHash)
uint baseHash(uvec3 p)
{
	const uint PRIME32_2 = 2246822519U, PRIME32_3 = 3266489917U;
	const uint PRIME32_4 = 668265263U, PRIME32_5 = 374761393U;
	uint h32 =  p.z + PRIME32_5 + p.x*PRIME32_3;
	h32 = PRIME32_4*((h32 << 17) | (h32 >> (32 - 17)));
	h32 += p.y * PRIME32_3;
	h32 = PRIME32_4*((h32 << 17) | (h32 >> (32 - 17))); //Initial testing suggests this line could be omitted for extra perf
    h32 = PRIME32_2*(h32^(h32 >> 15));
    h32 = PRIME32_3*(h32^(h32 >> 13));
    return h32^(h32 >> 16);
}

//XXHash32 based (https://github.com/Cyan4973/xxHash)
//Slower, higher quality
uint baseHash(uvec2 p)
{
    const uint PRIME32_2 = 2246822519U, PRIME32_3 = 3266489917U;
	const uint PRIME32_4 = 668265263U, PRIME32_5 = 374761393U;
    uint h32 = p.y + PRIME32_5 + p.x*PRIME32_3;
    h32 = PRIME32_4*((h32 << 17) | (h32 >> (32 - 17))); //Initial testing suggests this line could be omitted for extra perf
    h32 = PRIME32_2*(h32^(h32 >> 15));
    h32 = PRIME32_3*(h32^(h32 >> 13));
    return h32^(h32 >> 16);
}

uint baseHash(uint p)
{
	const uint PRIME32_2 = 2246822519U, PRIME32_3 = 3266489917U;
	const uint PRIME32_4 = 668265263U, PRIME32_5 = 374761393U;
	uint h32 = p + PRIME32_5;
	h32 = PRIME32_4*((h32 << 17) | (h32 >> (32 - 17))); //Initial testing suggests this line could be omitted for extra perf
    h32 = PRIME32_2*(h32^(h32 >> 15));
    h32 = PRIME32_3*(h32^(h32 >> 13));
    return h32^(h32 >> 16);
}
#endif

//---------------------3D input---------------------
float hash13(uvec3 x)
{
    uint n = baseHash(x);
    return float(n)*(1.0/float(0xffffffffU));
}

vec2 hash23(uvec3 x)
{
    uint n = baseHash(x);
    uvec2 rz = uvec2(n, n*48271U); //see: http://random.mat.sbg.ac.at/results/karl/server/node4.html
    return vec2((rz.xy >> 1) & uvec2(0x7fffffffU))/float(0x7fffffff);
}

vec3 hash33(uvec3 x)
{
    uint n = baseHash(x);
    uvec3 rz = uvec3(n, n*16807U, n*48271U); //see: http://random.mat.sbg.ac.at/results/karl/server/node4.html
    return vec3((rz >> 1) & uvec3(0x7fffffffU))/float(0x7fffffff);
}

vec4 hash43(uvec3 x)
{
    uint n = baseHash(x);
    uvec4 rz = uvec4(n, n*16807U, n*48271U, n*69621U); //see: http://random.mat.sbg.ac.at/results/karl/server/node4.html
    return vec4((rz >> 1) & uvec4(0x7fffffffU))/float(0x7fffffff);
}

//---------------------2D input---------------------

float hash12(uvec2 x)
{
    uint n = baseHash(x);
    return float(n)*(1.0/float(0xffffffffU));
}

vec2 hash22(uvec2 x)
{
    uint n = baseHash(x);
    uvec2 rz = uvec2(n, n*48271U);
    return vec2((rz.xy >> 1) & uvec2(0x7fffffffU))/float(0x7fffffff);
}

vec3 hash32(uvec2 x)
{
    uint n = baseHash(x);
    uvec3 rz = uvec3(n, n*16807U, n*48271U);
    return vec3((rz >> 1) & uvec3(0x7fffffffU))/float(0x7fffffff);
}

vec4 hash42(uvec2 x)
{
    uint n = baseHash(x);
    uvec4 rz = uvec4(n, n*16807U, n*48271U, n*69621U); //see: http://random.mat.sbg.ac.at/results/karl/server/node4.html
    return vec4((rz >> 1) & uvec4(0x7fffffffU))/float(0x7fffffff);
}

//---------------------1D input---------------------
float hash11(uint x)
{
    uint n = baseHash(x);
    return float(n)*(1.0/float(0xffffffffU));
}

vec2 hash21(uint x)
{
    uint n = baseHash(x);
    uvec2 rz = uvec2(n, n*48271U); //see: http://random.mat.sbg.ac.at/results/karl/server/node4.html
    return vec2((rz.xy >> 1) & uvec2(0x7fffffffU))/float(0x7fffffff);
}

vec3 hash31(uint x)
{
    uint n = baseHash(x);
    uvec3 rz = uvec3(n, n*16807U, n*48271U); //see: http://random.mat.sbg.ac.at/results/karl/server/node4.html
    return vec3((rz >> 1) & uvec3(0x7fffffffU))/float(0x7fffffff);
}

vec4 hash41(uint x)
{
    uint n = baseHash(x);
    uvec4 rz = uvec4(n, n*16807U, n*48271U, n*69621U); //see: http://random.mat.sbg.ac.at/results/karl/server/node4.html
    return vec4((rz >> 1) & uvec4(0x7fffffffU))/float(0x7fffffff);
}



//================================================
// Fonts
//------------------------------------------------
// by kishimisu
// https://www.shadertoy.com/view/dsGXDt

#define FONT_TEXTURE iChannel2
#define CHAR_SPACING 0.44
#define makeStr(func_name) float func_name(vec2 u) { _print 
#define makeStr1i(func_name) float func_name(vec2 u, int i) { _print
#define makeStr1f(func_name) float func_name(vec2 u, float i) { _print
#define makeStr2f(func_name) float func_name(vec2 u, float i, float j) { _print
#define makeStrXX(func_name) float func_name(vec2 u, ...) { _print
#define _end    ); return d; }
#define _ch(i)  _ 65+int(i)
#define _dig(i) _ 48+int(i)
#define _dec(x, dec) ); d += _decimal(FONT_TEXTURE, u, x, dec); (0
#define _SPA    ); u.x -= CHAR_SPACING; (0
#define _EXC  _ 33 // " ! "
#define _DBQ  _ 34 // " " "
#define _NUM  _ 35 // " # "
#define _DOL  _ 36 // " $ "
#define _PER  _ 37 // " % "
#define _AMP  _ 38 // " & "
#define _QUOT _ 39 // " ' "
#define _LPR  _ 40 // " ( "
#define _RPR  _ 41 // " ) "
#define _MUL  _ 42 // " * "
#define _ADD  _ 43 // " + "
#define _COM  _ 44 // " , "
#define _SUB  _ 45 // " - "
#define _DOT  _ 46 // " . "
#define _DIV  _ 47 // " / "
#define _COL  _ 58 // " : "
#define _SEM  _ 59 // " ; "
#define _LES  _ 60 // " < "
#define _EQU  _ 61 // " = "
#define _GRE  _ 62 // " > "
#define _QUE  _ 63 // " ? "
#define _AT   _ 64 // " @ "
#define _LBR  _ 91 // " [ "
#define _ANTI _ 92 // " \ "
#define _RBR  _ 93 // " ] "
#define _UND  _ 95 // " _ "
#define _A _ 65
#define _B _ 66
#define _C _ 67
#define _D _ 68
#define _E _ 69
#define _F _ 70
#define _G _ 71
#define _H _ 72
#define _I _ 73
#define _J _ 74
#define _K _ 75
#define _L _ 76
#define _M _ 77
#define _N _ 78
#define _O _ 79
#define _P _ 80
#define _Q _ 81
#define _R _ 82
#define _S _ 83
#define _T _ 84
#define _U _ 85
#define _V _ 86
#define _W _ 87
#define _X _ 88
#define _Y _ 89
#define _Z _ 90
#define _a _ 97
#define _b _ 98
#define _c _ 99
#define _d _ 100
#define _e _ 101
#define _f _ 102
#define _g _ 103
#define _h _ 104
#define _i _ 105
#define _j _ 106
#define _k _ 107
#define _l _ 108
#define _m _ 109
#define _n _ 110
#define _o _ 111
#define _p _ 112
#define _q _ 113
#define _r _ 114
#define _s _ 115
#define _t _ 116
#define _u _ 117
#define _v _ 118
#define _w _ 119
#define _x _ 120
#define _y _ 121
#define _z _ 122
#define _0 _ 48
#define _1 _ 49
#define _2 _ 50
#define _3 _ 51
#define _4 _ 52
#define _5 _ 53
#define _6 _ 54
#define _7 _ 55
#define _8 _ 56
#define _9 _ 57
#define _print  float d = 0.; (u.x += CHAR_SPACING
#define _       ); u.x -= CHAR_SPACING; d += _char(FONT_TEXTURE, u,

// Print character
float _char(sampler2D s, vec2 u, int id) {
    vec2 p = vec2(id%16, 15. - floor(float(id)/16.));
         p = (u + p) / 16.;
         u = step(abs(u-.5), vec2(.5));
    return texture(s, p).r * u.x * u.y;
}

// Floating point debug
float _decimal(sampler2D FONT_TEXTURE, inout vec2 u, float n, int decimals) {
    float d = 0., N = 1.; // d is the final color, N the number of digits before the decimal

    if (n < 0.) {  // If the number is negative
        n *= -1.;  // Make it positive
        (0 _SUB ); // Print a minus sign
    }
    
    // Calculate the number of digits before the decimal point
    for (float x = n; x >= 10.; x /= 10.) N++;

    // Print the digits before the decimal point
    for (float i = 0.; i < N; i++) {        
        float magnitude = pow(10., N-i-1.);
        float leftDigit = floor(n / magnitude);
        n -= leftDigit * magnitude;
        
        (0 _dig(leftDigit) );
    }

    if (decimals > 0) {
        (0 _DOT ); // Print a dot
    }
    
    // Print the digits after the decimal point
    for (int i = 0; i < decimals; i++) {
        float firstDecimal = floor((n - floor(n)) * 10.);
        n *= 10.;
        
        (0 _dig(firstDecimal) );
    }
    
    return d;
}

//================================================
// Other
//------------------------------------------------
