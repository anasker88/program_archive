// 単純なレイトレーシングの雛形

struct Ray
{
    vec3 org;
    vec3 dir;
};

struct Hit
{
    float distanceToHitpoint;
    // 交差点での法線ベクトル
    vec3 normal;
};

// 各種パラメータの例
float FilmWidth(){return iResolution.x/100.;}
float FilmHeight(){return iResolution.y/100.;}
float FilmDistance(){return 8.;}

vec3 CameraFrom(){return vec3(5.,2.,3.);}
vec3 CameraTo(){return vec3(.2,.7,.2);}
vec3 CameraUp(){return vec3(0.,1.,0.);}

float LargeFloat(){return 1e+6;}

// 正規直交基底を計算する関数の例
void createOrthoNormalBasis(
    vec3 from,vec3 to,vec3 up,
    out vec3 u,out vec3 v,out vec3 w,out vec3 e
)
{
    // TODO: ベクトル正規化normalize()や外積cross()を用いて実装する。
    vec3 fromTo=to-from;
    w=normalize(fromTo);
    vec3 upCrossW=cross(up,w);
    u=normalize(upCrossW);
    v=cross(w,u);
    e=from;
}

vec3 convertToCameraCoordinateSystem(vec2 pixelCoordinate)
{
    float filmX=(pixelCoordinate.x/iResolution.x-.5)*FilmWidth();
    float filmY=(pixelCoordinate.y/iResolution.y-.5)*FilmHeight();
    return vec3(filmX,filmY,-FilmDistance());
}

Ray generateCameraRay(
    vec2 pixelCoordinate
)
{
    // TODO: 以下を実装する。
    // 1. ピクセル座標をカメラ座標系に変換
    vec3 cameraCoordinate=convertToCameraCoordinateSystem(pixelCoordinate);
    // 2. カメラパラメータからカメラ座標系の正規直交基底を計算。
    vec3 u,v,w,e;
    createOrthoNormalBasis(CameraFrom(),CameraTo(),CameraUp(),u,v,w,e);
    // 3. ピクセル座標を基底を用いてワールド座標系に変換
    vec3 worldCoordinate=e+cameraCoordinate.x*u+cameraCoordinate.y*v+cameraCoordinate.z*w;
    // 4. カメラレイを計算。
    Ray cameraRay;
    cameraRay.org=e;
    cameraRay.dir=normalize(worldCoordinate-e);
    return cameraRay;
}

bool intersectToSphere(
    vec3 center,float radius,Ray ray,
    out Hit hit
)
{
    // TODO: レイと球の交差判定を実装する。
    // 二次方程式の解の計算に帰着する。
    // 解が両方とも正の場合に交差していると判定する。
    vec3 oc=ray.org-center;
    float a=dot(ray.dir,ray.dir);
    float b=2.*dot(oc,ray.dir);
    float c=dot(oc,oc)-radius*radius;
    float discriminant=b*b-4.*a*c;
    
    if(discriminant>0.){
        float t1=(-b-sqrt(discriminant))/(2.*a);
        float t2=(-b+sqrt(discriminant))/(2.*a);
        if(t1>0.&&t1<t2){
            hit.distanceToHitpoint=t1;
        }else if(t2>0.){
            hit.distanceToHitpoint=t2;
        }else{
            hit.distanceToHitpoint=LargeFloat();
            return false;
        }
        hit.normal=normalize(ray.org+hit.distanceToHitpoint*ray.dir-center);
        return true;
    }
    hit.distanceToHitpoint=LargeFloat();
    return false;
    
}

bool intersect(Ray ray,out Hit hit)
{
    hit.distanceToHitpoint=LargeFloat();
    
    // TODO: intersectToSphere を用いて具体的な球との交差判定を行う。
    intersectToSphere(vec3(0.,0.,0.),1.,ray,hit);
    
    return hit.distanceToHitpoint<LargeFloat();
}

vec3 shade(Ray ray,Hit hit)
{
    // TODO: なんらかのシェーディングを行う。
    vec3 lightDir=normalize(vec3(1.,1.,1.));
    float brightness=max(dot(hit.normal,lightDir),0.);
    return vec3(brightness);
}

void mainImage(out vec4 fragColor,in vec2 fragCoord)
{
    Ray ray=generateCameraRay(fragCoord);
    
    Hit hit;
    if(intersect(ray,hit))
    {
        fragColor=vec4(shade(ray,hit),0.);
    }
    else
    {
        fragColor=vec4(0.);
    }
}
