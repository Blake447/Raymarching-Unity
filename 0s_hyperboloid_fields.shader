Shader "Nave/Raymarching/Hyperbolic Fields"
{
	Properties
	{
		_AO("Ambient Occlusion", Range(0, 5)) = 1.0
		_Color("Color", Color) = (1,1,1,1)
		_Reflect("Reflectivity", Range(0,1)) = 0
		_Matte("Matte", Range(0,1)) = 0
		_Cubemap("Cubemap", Cube) = "black" {}
		_Epsilon("Epsilon", float) = .01
		_Delta("Delta", range(0,1)) = .02
		_LX("Light X position", range(-1,1)) = 1
		_LY("Light Y position", range(0,1)) = 1
		_LZ("Light Z position", range(-1,1)) = 1
	}
	SubShader
	{
		Tags { "RenderType"="Transparent" "Queue"="Transparent-2" }
		LOD 100
		Cull Front
		ZWrite Off

		Pass
		{
			CGPROGRAM
			#pragma vertex vert
			#pragma fragment frag
			
			#include "UnityCG.cginc"
			#define PI 3.14159265
			#define ROOT_HALF .707106781
			struct appdata
			{
				float4 vertex : POSITION;
			};

			struct v2f
			{
				float4 vertex : SV_POSITION;
				float3 wpos : TEXCOORD0;
				float3 wcenter : TEXCOORD1;
			};

			struct fragOut
			{
				float4 color : SV_Target;
				float depth : SV_Depth;
			};
			
			sampler2D_float _CameraDepthTexture;
			samplerCUBE _Cubemap;
			float _AO;
			float _Reflect;
			fixed4 _Color;
			float _Epsilon;
			float _Delta;
			float _Matte;
			float _LX;
			float _LY;
			float _LZ;

			// Beginning of a transcribed library from GLS:
			float safeAbs(float x)
			{
				if (x == 0)
				{
					return 0;
				}
				return abs(x);

			}

			float sgn(float x)
			{
				return step(0, x) * 2 - 1;
			}

			float2 sgn(float2 v)
			{
				return float2(step(0, v.x) * 2 - 1, step(0, v.y) * 2 - 1);
			}

			float lengthSqr(float3 x)
			{
				return dot(x, x);
			}

			float vmax(float2 v)
			{
				return max(v.x, v.y);
			}

			float vmax(float3 v)
			{
				return max(max(v.x, v.y), v.z);
			}

			float vmax(float4 v)
			{
				return max(max(max(v.x, v.y), v.z), v.w);
			}

			float vmin(float2 v)
			{
				return min(v.x, v.y);
			}

			float vmin(float3 v)
			{
				return min(min(v.x, v.y), v.z);
			}

			float vmin(float4 v)
			{
				return min(min(min(v.x, v.y), v.z), v.w);
			}

			void pR(inout float2 p, float a)
			{
				p = cos(a)*p + sin(a)*float2(p.y, -p.x);
			}

			void pR45(inout float2 p)
			{
				p = (p + float2(p.y, -p.x))*ROOT_HALF;
			}

			void pR45(inout float2 p, float2 offset)
			{
				p = (p-offset + float2(p.y-offset.y, -p.x+offset.x))*ROOT_HALF + offset;
			}

			float pMod1(inout float p, float size)
			{
				float halfsize = size * 0.5;
				float c = floor((p + halfsize) / size);
				p = fmod(abs(p) + halfsize, size) - halfsize;
				return c;

			}

			float pModSingle1(inout float p, float size)
			{
				float halfsize = size * .5;
				float c = floor((p + halfsize) / size);
				if (p >= 0)
					p = fmod(p + halfsize, size) - halfsize;
				return c;
			}

			float pModPolar(inout float2 p, float repetitions) {
				float angle = 2 * PI / repetitions;
				float r = length(p);
				float a = (atan2(p.x, p.y) + PI) + angle * 0.5;
				float c = floor(a / angle);
				a = fmod(a, angle) - angle * 0.5;
				p = float2(cos(a), sin(a))*r;
				if (abs(c) >= (repetitions/2)) c = abs(c);
				return c;
			}

			float pModPolar(inout float2 p, float repetitions, float phase) {
				float angle = 2 * PI / repetitions;
				float r = length(p);
				float a = (atan2(p.x, p.y) + PI + (2 * PI * phase / 360)) + angle * 0.5;
				float c = floor(a / angle);
				a = fmod(a, angle) - angle * 0.5;
				p = float2(cos(a), sin(a))*r;
				if (abs(c) >= (repetitions / 2)) c = abs(c);
				return c;
			}

			void pModMirror(inout float p)
			{
				p = abs(p);
			}

			float fOpUnionChamfer(float a, float b, float r)
			{
				return min(min(a, b), (a - r + b)*ROOT_HALF);
			}

			float fOpIntersectionChamfer(float a, float b, float r)
			{
				return max(max(a, b), (a - r + b)*ROOT_HALF);
			}

			float fOpDifferenceChamfer(float a, float b, float r)
			{
				return fOpIntersectionChamfer(a, -b, r);
			}

			float fOpUnionRound(float a, float b, float r)
			{
				float2 u = max(float2(r + a, r + b), float2(0.0, 0.0));
				return max(r, min(a, b)) - length(u);
			}

			float fOpIntersectionRound(float a, float b, float r)
			{
				float2 u = max(float2(r + a, r + b), float2(0.0, 0.0));
				return min(r, min(a, b)) + length(u);
			}

			float fOpDifferenceRound(float a, float b, float r)
			{
				return fOpIntersectionRound(a, -b, r);
			}

			float fOpUnionSoft(float a, float b, float r)
			{
				float e = max(r - abs(a - b), 0);
				return min(a, b) - e * e*0.25 / r;
			}
			
			float fOpUnionStairs(float a, float b, float r, float n)
			{
				float s = r / n;
				float u = b - r;
				return min(min(a, b), 0.5 * (u + a + abs((fmod(u - a + s, 2 * s)) - s)));
			}

			float fOpDifferenceStairs(float a, float b, float r, float n)
			{
				return -fOpUnionStairs(-a, b, r, n);
			}

			// This one doesn't quite work all the time and i'm not sure why, so be careful with it
			float fOpPipe(float a, float b, float r)
			{
				return length(float2(a, b)) - r;
			}

			float fOpGroove(float a, float b, float ra, float rb)
			{
				return max(a, min(a + ra, rb - abs(b)));
			}

			float fOpTongue(float a, float b, float ra, float rb)
			{
				return min(a, max(a - ra, abs(b) - rb));
			}

			float fPlane(float3 p, float3 n, float3 c)
			{
				float dist = dot(p, n) + c;
				return dist;
			}

			float fSphere(float3 p, float r, float3 c) 
			{
				return length(p - c) - r;
			}
			
			float fBox(float3 p, float3 b)
			{
				float3 d = abs(p) - b;
				return max(max(d.x, d.y), d.z);
			}

			float fBox(float3 p, float3 b, float3 c)
			{
				float3 d = abs(p-c) - b;
				return max(max(d.x, d.y), d.z);
			}

			float fPillar(float3 p, float3 c, float r, float height)
			{
				float d = length(p.xz - c.xz) - r;
				d = max(d, abs(p.y - c.y) - height);
				return d;
			}

			// End of functions transcribed from GLSL for unity. These functions from here on out are personally derived
			float2 intersectionOfTwoLines(float2 A, float sa, float2 B, float sb)
			{
				float buffer = (sa*A.x - A.y - sb * B.x + B.y) / (sa - sb);
				float2 intersection = float2(buffer, sa*(buffer - A.x) + A.y);
				return intersection;
			}

			// SDF for a paraboloid, including scale support (though not the best option for optimization purposes)
			float fAltParaboloid(float3 Point, float3 Center, float scale)
			{

				// Transform into polar coordinates
				float h = sqrt(dot(Point.xz - Center.xz, Point.xz - Center.xz))/scale;
				float k = (Point.y-Center.y)/scale;

				// Calculate point A
				float2 A = float2(h, k);

				// Calculate the point on the curve at h
				float2 r_of_h = float2(h, h*h);

				// Calculate the slope of the curve at h (differentiating x^2 gives us 2c)
				float t_of_h = 2 * h;

				// Calculate the slope of the normal of tangent line through (h,k)
				float n_of_h = -1 / t_of_h;

				// Calculate the intersection of the line normal to the tangent through (h,k) and the tangent line itself through (h, f(h))
				float2 t_int_n = intersectionOfTwoLines(r_of_h, t_of_h, A, n_of_h);

				// Calcuate two points for a chord line for underestimating the distance from inside the paraboloid
				float2 P = float2(h, h*h);
				float2 Q = float2(sqrt(k), k);

				// calculate the slope of the chord line and its normal
				float s_of_P2Q = (P.y - Q.y) / (P.x - Q.x);
				float n_of_P2Q = -1 / s_of_P2Q;

				// intersect the chord line with the normal line that goes through our point in question
				float2 R = intersectionOfTwoLines(P, s_of_P2Q, A, n_of_P2Q);

				// calculate the distance from the point in question to the tangent line
				float dist = distance(A, t_int_n);

				// If we are above the paraboloid
				if (k > r_of_h.y)
				{
					// calculate the minimum distance to the chord line instead
					dist = -distance(A, R);
				}

				// Return the distance
				return dist * scale;
			}

			// SDF for the hyperboloid of Two Sheets.
			float fAltHyperboloidTwoSheet(float3 Point, float3 Center, float scale)
			{
				// Transform into polar coordinates
				float h = sqrt(dot(Point.xz - Center.xz, Point.xz - Center.xz))/scale;
				float k = safeAbs((Point.y-Center.y)/scale);

				// Calculate point A
				float2 A = float2(h, k);

				// Calculate the point on the curve at h
				float2 r_of_h = float2(h, sqrt(h*h+1));

				// Calculate the slope of the curve at h. differentiating y = sqrt(x^2 + 1) gives x / sqrt((x^2 + 1))
				float t_of_h = h / r_of_h.y;

				// Calculate the slope of the normal of tangent line through (h,k). 
				float n_of_h = -1 / t_of_h;

				// Calculate the intersection of the line normal to the tangent through (h,k) and the tangent line itself through (h, f(h))
				float2 t_int_n = intersectionOfTwoLines(r_of_h, t_of_h, A, n_of_h);

				// Calculate the points for the chord line
				float2 P = r_of_h;
				float2 Q = float2(sqrt(k*k-1), k);

				// Calculate the slope oof the chord line and its normal
				float s_of_P2Q = (P.y - Q.y) / (P.x - Q.x);
				float n_of_P2Q = -1 / s_of_P2Q;

				
				// intersect the chord line with the normal line that goes through our point in question
				float2 R = intersectionOfTwoLines(P, s_of_P2Q, A, n_of_P2Q);

				// calculate the distance from the point in question to the tangent line
				float dist = distance(A, t_int_n);

				// If we are above the paraboloid
				if (k > r_of_h.y)
				{
					// calculate the minimum distance to the chord line instead
					dist = -distance(A, R);
				}

				// Return the distance
				return dist * scale;
			}

			// SDF for a hyperboloid of one sheet. Note that the only difference really is that the h and k are switched. 
			float fAltHyperboloid(float3 Point, float3 Center, float scale)
			{
				// Transform into polar coordinates
				float k = sqrt(dot(Point.xz - Center.xz, Point.xz - Center.xz))/scale;
				float h = safeAbs((Point.y-Center.y)/scale);

				// Calculate point A
				float2 A = float2(h, k);

				// Calculate the point on the curve at h
				float2 r_of_h = float2(h, sqrt(h*h+1));

				// Calculate the slope of the curve at h. differentiating y = sqrt(x^2 + 1) gives x / sqrt((x^2 + 1))
				float t_of_h = h / r_of_h.y;

				// Calculate the slope of the normal of tangent line through (h,k)
				float n_of_h = -1 / t_of_h;

				// Calculate the intersection of the line normal to the tangent through (h,k) and the tangent line itself through (h, f(h))
				float2 t_int_n = intersectionOfTwoLines(r_of_h, t_of_h, A, n_of_h);

				// Calculate the points for the chord line
				float2 P = r_of_h;
				float2 Q = float2(sqrt(k*k-1), k);

				// Calculate the slope oof the chord line and its normal
				float s_of_P2Q = (P.y - Q.y) / (P.x - Q.x);
				float n_of_P2Q = -1 / s_of_P2Q;

				
				// intersect the chord line with the normal line that goes through our point in question
				float2 R = intersectionOfTwoLines(P, s_of_P2Q, A, n_of_P2Q);

				// calculate the distance from the point in question to the tangent line
				float dist = distance(A, t_int_n);

				// If we are above the paraboloid
				if (k > r_of_h.y)
				{
					// calculate the minimum distance to the chord line instead
					dist = -distance(A, R);
				}

				// Return the distance
				return -dist * scale;
			}

			// Versions for when scale is not provided
			float fAltParaboloid(float3 Point, float3 Center)
			{

				// Transform into polar coordinates
				float h = sqrt(dot(Point.xz - Center.xz, Point.xz - Center.xz));
				float k = (Point.y-Center.y);

				// Calculate point A
				float2 A = float2(h, k);

				// Calculate the point on the curve at h
				float2 r_of_h = float2(h, h*h);

				// Calculate the slope of the curve at h
				float t_of_h = 2 * h;

				// Calculate the slope of the normal of tangent line through (h,k)
				float n_of_h = -1 / t_of_h;

				// Calculate the intersection of the line normal to the tangent through (h,k) and the tangent line itself through (h, f(h))
				float2 t_int_n = intersectionOfTwoLines(r_of_h, t_of_h, A, n_of_h);

				float2 P = float2(h, h*h);
				float2 Q = float2(sqrt(k), k);

				float s_of_P2Q = (P.y - Q.y) / (P.x - Q.x);
				float n_of_P2Q = -1 / s_of_P2Q;

				float2 R = intersectionOfTwoLines(P, s_of_P2Q, A, n_of_P2Q);




				float dist = distance(A, t_int_n);

				if (k > r_of_h.y)
				{
					dist = -distance(A, R);
				}

				return dist;
			}

			float fAltHyperboloidTwoSheet(float3 Point, float3 Center)
			{
				// Transform into polar coordinates
				float h = sqrt(dot(Point.xz - Center.xz, Point.xz - Center.xz));
				float k = safeAbs((Point.y-Center.y));

				// Calculate point A
				float2 A = float2(h, k);

				// Calculate the point on the curve at h
				float2 r_of_h = float2(h, sqrt(h*h+1));

				// Calculate the slope of the curve at h
				float t_of_h = h / r_of_h.y;

				// Calculate the slope of the normal of tangent line through (h,k)
				float n_of_h = -1 / t_of_h;

				// Calculate the intersection of the line normal to the tangent through (h,k) and the tangent line itself through (h, f(h))
				float2 t_int_n = intersectionOfTwoLines(r_of_h, t_of_h, A, n_of_h);

				float2 P = r_of_h;
				float2 Q = float2(sqrt(k*k-1), k);

				float s_of_P2Q = (P.y - Q.y) / (P.x - Q.x);
				float n_of_P2Q = -1 / s_of_P2Q;

				float2 R = intersectionOfTwoLines(P, s_of_P2Q, A, n_of_P2Q);

				float dist = distance(A, t_int_n);

				if (k > r_of_h.y)
				{
					dist = -distance(A, R);
				}

				return dist;
			}

			float fAltHyperboloid(float3 Point, float3 Center)
			{
				// Transform into polar coordinates
				float k = sqrt(dot(Point.xz - Center.xz, Point.xz - Center.xz));
				float h = safeAbs((Point.y-Center.y));

				// Calculate point A
				float2 A = float2(h, k);

				// Calculate the point on the curve at h
				float2 r_of_h = float2(h, sqrt(h*h+1));

				// Calculate the slope of the curve at h
				float t_of_h = h / r_of_h.y;

				// Calculate the slope of the normal of tangent line through (h,k)
				float n_of_h = -1 / t_of_h;

				// Calculate the intersection of the line normal to the tangent through (h,k) and the tangent line itself through (h, f(h))
				float2 t_int_n = intersectionOfTwoLines(r_of_h, t_of_h, A, n_of_h);

				float2 P = r_of_h;
				float2 Q = float2(sqrt(k*k-1), k);

				float s_of_P2Q = (P.y - Q.y) / (P.x - Q.x);
				float n_of_P2Q = -1 / s_of_P2Q;

				float2 R = intersectionOfTwoLines(P, s_of_P2Q, A, n_of_P2Q);

				float dist = distance(A, t_int_n);

				if (k > r_of_h.y)
				{
					dist = -distance(A, R);
				}

				return -dist;
			}

			// Now the rest is Neen and Naves work, with maybe a small modification to allow reflectivity with cube maps,
			// and then clipping the raymarched shader if the distance gets too far awa

			float DE(float3 p)
			{

				if (dot(p, p) < 1951)
				{

					float3 q = p;
					pModPolar(q.xz, 2);
					
					float3 cPar = float3(0, 0, 0);
					float d = fAltHyperboloidTwoSheet(p, cPar + float3(2, 0, 2), .5);
					d = min(d,fAltHyperboloid(p, cPar + float3(-2, 0, -2), .5));
					d = min(d, fAltParaboloid(p, cPar + float3(2, 0, -2)));
					d = min(d, fAltParaboloid(p*float3(1, -1, 1), cPar + float3(-2, 0, 2)));

					d = max(d, fBox(p, float3(20, 2, 20)));

					return d;
				}
				return -1.0;
			}

			float3 calcNormal(float3 p)
			{
				float2 e = float2(1.0, -1.0) * .0005;
				return normalize(
					e.xyy * DE(p + e.xyy) +
					e.yyx * DE(p + e.yyx) +
					e.yxy * DE(p + e.yxy) +
					e.xxx * DE(p + e.xxx));
			}
		
			
			v2f vert (appdata v)
			{
				v2f o;
				o.vertex = UnityObjectToClipPos(v.vertex);
				o.wpos = mul(unity_ObjectToWorld, v.vertex).xyz;
				o.wcenter = mul(unity_ObjectToWorld, float4(0, 0, 0, 1)).xyz;
				return o;
			}
			
			fragOut frag (v2f i)
			{
				//Heyy Neen, I'm taking the time to comment this for you :)
				fragOut f;

				//Alright so this just defines the maximum amount of steps we'll use, I probably don't
				//have to explain something as simple as this but whatever :D
				float maxStep = 64;
				float4 clipPos = float4(0,0,0,0);
				float clipDepth = -1;
				float3 normal = float3(0, 1, 0);

				fixed4 colorHit = float4(0,0,0,1);

				// Calculate world space ray direction and origin
				float3 raydir = normalize(i.wpos - _WorldSpaceCameraPos.xyz);
				float3 raypos = _WorldSpaceCameraPos.xyz - i.wcenter;

				//alright, time to raytrace, simple for-loop with maxStep iterations
				for (int stp = 0; stp < maxStep; stp++) {
					//alright, here we simply get the distance to the nearest object calculated by our amazing DE
					//or "Distance Estimation" function, genius, you can look into the comments on that one further above
					float d = DE(raypos);
					//now if the distance is super small, that means we hit something, I tried checking against <= 0.0 but
					//that made everything noisy and stuff so we'll just use something super tiny here
					if (d < 0)
					{
						break;
					}
					//if (d <= 0.008 + .04*stp*.015625) {
					if (d <= (0.008 + .04*stp*.015625 + .001*dot(raypos, raypos)*0.01625)*_Epsilon) {



						//Now if we did hit something, we just return white times my magic ambient occlusion formular... Ohhhh :O
						//it's super simple tho, the main core is (stp / maxStep) which basically gives us the ratio of how many
						//steps it took to get here, if we hit something on the first step, then stp/maxStep is gonna be super small
						//but if it was the last step of the for loop then stp/maxStep is basically almost 1...
						//so if we hit something early it's gonna be close to 0 and if it's super far away or we needed a lot of steps
						//it's gonna be close to 1... then we just invert that with a "1 - bla" so, far is 1 and near is 0 and then
						//we take all of that and take that to a power of something which has a slider so you can kinda play around with
						//the "intensity" a little bit... Oh yes, also, if a ray takes multiple steps to get somewhere, not only does that
						//mean it may be far away, but it could also mean the surface it hit was more complex to get to, that's why you
						//see the spheres having some small gradient torwards the edges, which looks cool!
						clipPos = mul(UNITY_MATRIX_VP, float4(raypos + i.wcenter, 1.0));
						clipDepth = clipPos.z / clipPos.w;
		
						normal = calcNormal(raypos);
						break;
					}
					
					//oh yes, also, if we didn't hit something we just add the direction times our minimum distance to the nearest
					//object to the current ray pos and take the next step
					raypos += raydir * d;
				}
				

				if (dot(raypos, raypos) > 1950|| stp > maxStep-1)
				{
					clip(-1);
				}

				float3 lightVector = normalize(float3(_LX, _LY, _LZ));
				//also look, if it went through all steps and didn't hit anything then just return pure black... it must have either
				//went on for infinity or found itself in a super complex crack of some fractal surface or whatever
				colorHit = texCUBE(_Cubemap, reflect(raydir, normal));
				//colorHit = _Color*(dot(normal, lightVector) + 1)*.2;
				float4 matte = _Color * (dot(normal, lightVector) + 1)*.4;
				
				float4 reflective = colorHit * float4(1, 1, 1, 1) * saturate(pow(1 - stp / maxStep, _AO))*_Color * 4;
				//float4 matte = float4(1, 1, 1, 1) * saturate(pow(1 - stp / maxStep, _AO))*_Color * 4;
					//f.color = float4(0, 0, 0, 1) + colorHit* float4(1, 1, 1, 1) * saturate(pow(1 - stp / maxStep, _AO))*_Color*4;
				f.color = (reflective * _Reflect + matte * (1 - _Reflect));
				f.depth = clipDepth;
				return f;
			}
			ENDCG
		}
	}
}
