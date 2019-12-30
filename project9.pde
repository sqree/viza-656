//
// VIZA 656: Image Synthesis
// Name: Anne Fu
// Project 9
//

float specHighlight (PVector znorm, PVector light, PVector viewer, float specStrength) {
  PVector reflect = (PVector.mult(znorm, 2.0*(znorm.dot(PVector.mult(light, -1))))).add(light);
  reflect.normalize();
  float specLight = ((max(viewer.dot(reflect), 0.0)-0.7)/(1-0.7))*specStrength;
  if (specLight>0.9) {
    specLight = 1;
  } else if (specLight<0) {
    specLight = 0;
  }
  return specLight;
}

class light {
  PVector Pl;
  PVector lUp, lN0, lN1, lN2; //for areaLight
  PVector normL; //for spotlight
  float lAngle; //for spotlight
  float lS0, lS1, lM, lN; // size & sampling size for arealight
  color colorL;
  String lightType;
  light (String typeL, PVector lightCenter, PVector lightU, PVector lightN, float angleL, float lightSX, float lightSY, float lightX, float lightY, color lightC) {
    lightType = typeL; //point, spot, direct
    Pl = lightCenter; //light position
    lUp = lightU; //up vector for arealight
    normL = lightN; //light direction for spot light & view direction for area light
    lN2 = normL.normalize(); //local coordinates calculations for area light
    lN0 = lUp.cross(normL);
    lN1 = lN2.cross(lN0); 
    lAngle = angleL; //angle for spot light
    lS0 = lightSX; //size of arealight
    lS1 = lightSY;
    lM = lightX; //sublight sections
    lN = lightY;
    colorL = lightC;
  }
  PVector lightRay(PVector hitPoint) {
    PVector lRay;
    switch(lightType) {
    case "direct":
      lRay = normL;
      break;
    default:
      lRay = PVector.sub(hitPoint, Pl);
      break;
    }
    lRay.normalize();
    return lRay;
  }
  float castShadow (PVector hitPoint, PVector lightRay, String sMode, geometry[] sceneGeo, PVector hitNorm) {
    PVector rayL = lightRay;
    float cosTheta = 0;
    if (lightType.equals("direct")) {
      Pl = new PVector(hitPoint.x-lightRay.x*400, hitPoint.y-lightRay.y*400, hitPoint.z-lightRay.z*400);
    }
    //check for shadow        
    float th = sceneGeo[0].intersect(rayL, Pl, sceneGeo[0].rad);
    for (int j=1; j<scene.length; j=j+1) {
      if (th>sceneGeo[j].intersect(rayL, Pl, sceneGeo[j].rad)) {
        th = sceneGeo[j].intersect(rayL, Pl, sceneGeo[j].rad);
      }
    }
    if ((hitPoint.x-Pl.x)/rayL.x>th+0.5 || (lightType.equals("spot") && normL.dot(rayL)<lAngle)) {
      switch(sMode) {
      case "modedR":
        //d/R calculation
        PVector Pray = hitPoint;
        float tempT = 1;
        float dR = 1;
        //trace along light ray
        while (PVector.sub(Pray, Pl).mag() > 1) {
          Pray = PVector.sub(hitPoint, PVector.mult(rayL, tempT));
          tempT = tempT + 1;
          for (int j=1; j<3; j = j+1) {
            //increment if inside sphere -- I should probably use the implicit function for these but since there are only planes and spheres in this scene...
            if (PVector.sub(Pray, scene[j].centerPoint).mag() < scene[j].rad) {
              dR = dR + 1;
            }
          }
        }
        //normalize shadow
        cosTheta = 3*(dR/PVector.sub(Ph, lightPos).mag());
        break;
      case "modeArea":
        //arealight calculation
        PVector tempLightPos = Pl;
        PVector tempLight = lightRay;
        float tempCount = 0;
        float randomOffset = random(0, lM-1); //randomizing sublight sequence; since M=N can use same variable
        for (int j=0; j<lM; j = j+1) {
          for (int k=0; k<lN; k = k+1) {
            //X, Y coordinates essentially simplify to randomized sublight section + random offset for random point within sublight section
            tempLightPos = PVector.add(PVector.add(tempLightPos, PVector.mult(lN0, lS0*(((j+randomOffset)% lM)/lM) + random(0, lS0/lM))), PVector.mult(lN1, lS1*(((k+randomOffset)%lN)/lN)+ random(0, lS1/lN)));
            tempLight = PVector.sub(hitPoint, tempLightPos);
            tempLight.normalize();
            //check for intersection between light ray and objects in scene
            th = scene[0].intersect(tempLight, tempLightPos, scene[0].rad);
            for (int l=1; l<scene.length; l=l+1) {
              if (th>scene[l].intersect(tempLight, tempLightPos, scene[l].rad)) {
                th = scene[l].intersect(tempLight, tempLightPos, scene[l].rad);
              }
            }
            //increment shadow counter if intersection found
            if ((hitPoint.x-tempLightPos.x)/tempLight.x>th+0.5) {
              tempCount = tempCount+1;
            }
          }
        }
        //normalize shadow
        cosTheta = tempCount/(lM*lN);
        break;
      default:
        cosTheta = hitNorm.dot(PVector.mult(rayL, -1));
        //cosTheta = 0.5;
        break;
      }
    }
    return cosTheta;
  }
}

class geometry {
  PVector centerPoint; 
  PVector pNormal;
  String geoType;
  color geoColor;
  color geoDarkColor;
  float rad;
  PVector p0, p1, p2;
  boolean refr, refl;
  geometry (PVector position, PVector planeNormal, String type, color objectColor, color darkColor, float radius, PVector pt0, PVector pt1, PVector pt2, boolean trs, boolean ref) {
    centerPoint = position; 
    geoType = type;
    geoColor = objectColor;
    geoDarkColor = darkColor;
    switch(geoType) {
    case "sphere":
      rad = radius;
      pNormal = new PVector(0, 0, 0);
      p0 = p1 = p2 = new PVector(0,0,0);
      refr = trs;
      refl = ref;
      break;
    case "plane":
      rad = Float.NaN;
      pNormal = planeNormal;
      p0 = p1 = p2 = new PVector(0,0,0);
      refr = trs;
      refl = ref;
      break;
    case "triangle":
      rad = Float.NaN;
      pNormal = planeNormal;
      p0 = pt0;
      p1 = pt1;
      p2 = pt2;
      refr = trs;
      refl = ref;
    default:
      break;
    }
  } 
  PVector normal (PVector hitPoint) {
    PVector geoNorm = new PVector(0, 0, 0);
    switch(geoType) {
    case "sphere":
      geoNorm = PVector.sub(hitPoint, centerPoint);
      break;
    case "plane":
      geoNorm = pNormal;
      break;
    case "triangle":
      geoNorm = pNormal;
      break;
    default:
      break;
    }
    geoNorm.normalize();
    return geoNorm;
  }
  float intersect (PVector ray, PVector startPoint, float radius) {
    float intersectPoint=0;
    switch(geoType) {
    case "sphere":
      float b = (ray.dot(PVector.sub(centerPoint, startPoint)));
      float c = (PVector.sub(centerPoint, startPoint)).dot(PVector.sub(centerPoint, startPoint))-radius*radius;
      intersectPoint = b - sqrt(b*b-c);
      if(intersectPoint<0) {
        intersectPoint = b + sqrt(b*b-c);
      }
      break;
    case "plane":
      intersectPoint = (pNormal.dot(PVector.sub(centerPoint, startPoint)))/(pNormal.dot(ray));
      break;
    case "triangle":
      //calculate hit point with plane coinciding with triangle
      intersectPoint = (pNormal.dot(PVector.sub(centerPoint, startPoint)))/(pNormal.dot(ray));
      PVector tPh = PVector.add(startPoint, PVector.mult(ray, intersectPoint));
      //check if hit point inside triangle by calculating area of sub-triangles
      PVector A = PVector.sub(p1,p0).cross(PVector.sub(p2,p0));
      PVector A0 = PVector.sub(tPh,p1).cross(PVector.sub(tPh,p2));
      PVector A1 = PVector.sub(tPh,p2).cross(PVector.sub(tPh,p0));
      PVector A2 = PVector.sub(tPh,p0).cross(PVector.sub(tPh,p1));
      float u,v,w;
      u = v = w = 0;
      if(max(abs(A.x),abs(A.y),abs(A.z)) == abs(A.x)) {
        u = A1.x/A.x;
        v = A2.x/A.x;
        w = A0.x/A.x;
      }
      else if(max(abs(A.x),abs(A.y),abs(A.z)) == abs(A.y)) {
        u = A1.y/A.y;
        v = A2.y/A.y;
        w = A0.y/A.y;
      }
      else if(max(abs(A.x),abs(A.y),abs(A.z)) == abs(A.z)) {
        u = A1.z/A.z;
        v = A2.z/A.z;
        w = A0.z/A.z;
      }
      //check if areas valid
      if((u<0 || u>1)||(v<0 || v>1)||(1-u-v<0 || 1-u-v>1)|| Float.isNaN(u) || Float.isNaN(v)) {
        intersectPoint = Float.MAX_VALUE;
      }
      break;
    default:
      break;
    }
    if (Float.isNaN(intersectPoint) || intersectPoint<0) { //check for invalid intersection
      intersectPoint = Float.MAX_VALUE;
    }
    return intersectPoint;
  }
}

PImage img, nimg, envimg, wimg, bimg, dimg;
PImage[] timg;
float x, y, z;
float sX, sY;
float th;
float sphereStart;
color lightColor, darkColor, objectColor, specLightColor;
PVector znorm;
PVector light, lightPos;
PVector nL; //for spot, area light
float spotAngle; //for spot light
PVector Pe, Ph, P0;
PVector vView, vUp;
PVector n0, n1, n2;
PVector[][] vertex; 
int totalImageSize;
float ambient, diffuse, specStrength, spec;
color[][] rawColor;
String mode, modeM, mapM; //shadow calculation mode hard, d/R, area, mode/map type
String lightT; //point, spot, direct

String imgName;
String imgPath, nimgPath;

PShape cub,cub1;
PShape tetra;
PShape dodeca;
PImage texture;
color iri;

geometry[] scene;
light[] lightSet;

void setup() {
  size(500, 500, P3D);  
  //camera set up
  Pe = new PVector(50, 0, 300);
  vView = new PVector(0, 0, 1);
  vView.normalize();
  vUp = new PVector (0, 1, 0);
  n2 = vView.normalize();
  n0 = (vView.cross(vUp)).normalize();
  n1 = n0.cross(n2);
  P0 = PVector.add(Pe, PVector.mult(n2, 50));
  lightPos = new PVector(-50, 100, 100);

  //light setup
  nL = new PVector(-0.5, -1, 0.5); //for direct, spot, area
  nL.normalize();
  spotAngle = 0.85;
  ambient = 0.5;
  specStrength = 0.9;
  specLightColor = color(255, 253, 239);

  mode = "reflect";
  lightT = "";
  modeM = "default";
  mapM = "";
  
  cub=loadShape("cube.obj");
  cub1=loadShape("cube1.obj");
  tetra = loadShape("tetrahedron.obj");
  dodeca = loadShape("dodecahedron.obj");
}

void mouseClicked() {
  image(nimg, 0, 0);
  lightPos.x = mouseX-img.width/2;
  lightPos.y = -(mouseY-img.width/2);
  lightPos.z = 150;
  nL = new PVector(-(mouseX-img.width/2), (mouseY-img.width/2), -200);
  nL.normalize();
  redraw();
  print(lightPos,"\n");
}

void keyPressed() {
  if (key == CODED) {
    if (keyCode == UP) {
      Pe.y = Pe.y - 5;
    } else if (keyCode == DOWN) {
      Pe.y = Pe.y + 5;
    } else if (keyCode == RIGHT) {
      Pe.x = Pe.x - 5;
    } else if (keyCode == LEFT) {
      Pe.x = Pe.x + 5;
    } else if (keyCode == ALT) {
      Pe.z = Pe.z - 5;
    } else if (keyCode == CONTROL) {
      Pe.z = Pe.z + 5;
    }
  } else if (key == 'd') {
    mode = "modedR";
    if (lightT.equals("spot")) {
      lightT = "point";
    }
  } else if (key == 'a') {
    mode = "modeArea";
    lightT = "point";
  } else if (key == 'h') {
    mode = "";
  } else if (key == 'l') {
    lightT = "direct";
    mode = "";
  } else if (key == 's') {
    lightT = "spot";
    mode = "";
  } else if (key == 'n') {
    mapM = "normal";
    if(!modeM.equals("default")){
      modeM = "default";
    }
  } else if (key == 'p') {
    modeM = "refract";
    mapM = "";
    lightPos = new PVector(-9,19,150);
  } else if (key == 'z') {
    modeM = "reflect";
    mapM = "";
  } else if (key == 'v') {
    modeM = "envMap";
    mapM = "";
  } else if (key == 'i') {
    modeM = "iridescent";
    mapM = "";
  } else if (key == 'y') {
    modeM = "ior";
    mapM = "";
  }
  redraw();
}

void draw() {

  PVector nPlane0 = new PVector(0, 1, 0);
  nPlane0.normalize();

  PVector nPlane1 = new PVector(0, 0, 1);
  nPlane1.normalize();

  lightSet = new light[1];
  lightSet[0] = new light(lightT, lightPos, new PVector(0, 1, 0), nL, 0.85, 20.0, 20.0, 4.0, 4.0, color(255, 247, 170));
  
  vertex = new PVector[2*cub.getChildCount()+dodeca.getChildCount()+tetra.getChildCount()][3];
  PVector[][] vt = new PVector[2*cub.getChildCount()+dodeca.getChildCount()+tetra.getChildCount()][3];
    
  //get all face vertices from obj -- not currently pulling normals
  //cube
  for(int i = 0;i<cub.getChildCount();i=i+1) {
    for(int j=0; j<cub.getChild(i).getVertexCount(); j=j+1) {
      PVector translate = new PVector(-50,5,40);
      //get face vertices
      vertex[i][j] = PVector.add(translate,PVector.mult(cub.getChild(i).getVertex(j),70));
      //get texture coordinates
      vt[i][j] = new PVector(cub.getChild(i).getTextureU(j), cub.getChild(i).getTextureV(j));
    }
  }
  
  nimgPath = "./images/output.png";
  nimg = createImage(2000, 2000, RGB);
  totalImageSize = nimg.width * nimg.height;
  nimg.loadPixels();
  
  //PImage textureN = loadImage("./images/rubix_nrm.png");
  
  envimg = loadImage("./images/mars.png");
  bimg = loadImage("./images/crystal.jpg");
  wimg = loadImage("./images/iri.png");
  rawColor = new color[nimg.width][nimg.height];
  
  for (int i = 0; i < totalImageSize-1; i = i + 1) {

    x = (i+1)%nimg.width;
    y = nimg.height-floor(((i+1)/nimg.width))-1;

    //shoot ray
    PVector nPe = new PVector((x-nimg.width/2-Pe.x)/5-P0.x, (y-nimg.height/2-Pe.y)/5-P0.y, Pe.z/5-P0.z);
    nPe.normalize();

    scene = new geometry[cub.getChildCount()+3];
    
    if(modeM.equals("envMap")) {
      scene = new geometry[cub.getChildCount()+2];
    }

    PVector point1;
    PVector point2;
    PVector point3;
    PVector pNormal;
    PVector pCenter;
    
    //create triangle geometry for obj faces
    //cube
    for (int j = 0;j<cub.getChildCount();j=j+1) {
      point1 = vertex[j][0];
      point2 = vertex[j][1];
      point3 = vertex[j][2];
      pNormal = PVector.sub(point2, point1).cross(PVector.sub(point3, point1));
      pNormal.normalize();
      pCenter = point1;
      scene[j+1] = new geometry(pCenter, pNormal, "triangle", color(26, 154, 139), color(134, 167, 174), 0, point1, point2, point3, false, false);
    }
    
    
    scene[cub.getChildCount()+1] = new geometry(new PVector(-50, 10, 120), new PVector(0, 0, 0), "sphere", color(152, 165, 73), color(78, 38, 43), 20,new PVector(0,0,0), new PVector(0,0,0), new PVector(0,0,0),false, false);
    scene[0] = new geometry(new PVector(30, 50, 0), new PVector(0, 0, 0), "sphere", color(152, 165, 73), color(115, 58, 47), 40,new PVector(0,0,0), new PVector(0,0,0), new PVector(0,0,0),false, false);
    
    if(!modeM.equals("envMap")) {
      scene[cub.getChildCount()+2] = new geometry(new PVector(0, -30, 0), nPlane0, "plane", color(147, 188, 190), color(134, 167, 174), 0 ,new PVector(0,0,0), new PVector(0,0,0), new PVector(0,0,0),false, false);
    }
    if(modeM.equals("reflect")){
      scene[0].refl = true;
      for (int j = 0;j<cub.getChildCount();j=j+1) {
        scene[j+1].refl = true;
      }
    }
    if(modeM.equals("refract")){
      scene[cub.getChildCount()+1].refr = true;
      for (int j = 0;j<cub.getChildCount();j=j+1) {
        scene[j+1].refr = true;
      }
    }
    if(modeM.equals("envMap") || mapM.equals("normal")) {
      scene[0].refl = true;
    }
    if(modeM.equals("ior")) {
      scene[0].refr = true;
    }

    float tTemp = scene[0].intersect(nPe, Pe, scene[0].rad)+1;
    //calculate plane hit point
    Ph = PVector.add(Pe, PVector.mult(nPe, tTemp));
    znorm = scene[0].normal(Ph);
    objectColor = scene[0].geoColor;
    darkColor = scene[0].geoDarkColor;
    
    for (int j = 0; j<scene.length; j=j+1) {
      if (tTemp>scene[j].intersect(nPe, Pe, scene[j].rad)) {
        tTemp = scene[j].intersect(nPe, Pe, scene[j].rad);
        Ph = PVector.add(Pe, PVector.mult(nPe, scene[j].intersect(nPe, Pe, scene[j].rad)));
        znorm = scene[j].normal(Ph);
        
        iri = color(0,0,0);
        
        //for reflection defined by normal map
        if(j==0 && mapM.equals("normal")) {
          PVector sphereN0 = new PVector(1, 0, 0);
          PVector sphereN1 = new PVector(0, 1, 0);
          PVector sphereN2 = new PVector(0, 0, 1);

          float sphereX = sphereN0.dot(PVector.sub(Ph, scene[j].centerPoint))/scene[j].rad;
          float sphereY = sphereN1.dot(PVector.sub(Ph, scene[j].centerPoint))/scene[j].rad;
          float sphereZ = sphereN2.dot(PVector.sub(Ph, scene[j].centerPoint))/scene[j].rad;
          
          float v = acos(sphereZ)/PI;
          float u = acos(sphereY/sqrt(1-sphereZ*sphereZ))/(2*PI);
          //flip u if x negative
          if(sphereX<0) {
            u = 1-u;
          }
          //normal map
          color cM = bimg.get((int)(u*(bimg.width)),(int)(v*(bimg.height)));
          //convert color to normal vector
          PVector nM = new PVector((red(cM)/256)*2-1,(green(cM)/256)*2-1,(blue(cM)/256)*2-1);
          nM.normalize();
          //average normals
          znorm = PVector.add(nM, znorm);
          znorm.normalize();
        }
        objectColor = scene[j].geoColor;
        darkColor = scene[j].geoDarkColor;
        
        light = PVector.sub(Ph,lightPos);
        light.normalize();
        float none = 1;
        float ntwo = 1.3;
        
        //for iridescence and ior defined by texture map
        if(j==0 && (modeM.equals("iridescent")||modeM.equals("ior"))) {
          float irix, iriy;
          irix = iriy = -light.dot(znorm);
          iri = wimg.get((int)(irix*(wimg.width)),(int)(iriy*(wimg.height)));
          if(modeM.equals("iridescent")) {
            objectColor = color(100,100,100);
            objectColor = color(irix*red(iri)+red(objectColor),irix*green(iri)+green(objectColor),irix*blue(iri)+blue(objectColor));
          }
          if(modeM.equals("ior")){
          //set IOR as a function of red in iridescence texture
          ntwo = 1+red(iri)/255;
          }
        }
      
        
        //reflect/refract calculation          
        if(scene[j].refr || scene[j].refl) {
          //calculate reflect vector
          PVector Re = (PVector.mult(znorm,2*znorm.dot(PVector.mult(nPe,-1)))).add(nPe);
          Re.normalize();
          
          if (scene[j].refr && PVector.mult(light,-1).dot(znorm)<0) {
            znorm = znorm.mult(-1);
          }
          
          //calculate transmit vector  
          float n = none/ntwo;
          float ci = max(-znorm.dot(light),0);
          float c2 = sqrt(1 - (n*n*(1-ci*ci)));
          
          PVector Tr = PVector.mult(light, n).add(PVector.mult(znorm,((n)*ci-c2)));
          Tr.normalize();
          
          //initializing variables for reflect
          boolean reflect = false;
          PVector Pr = PVector.add(Ph,PVector.mult(Re,0.1)); //offset to prevent self intersects
          PVector rnorm = znorm;
          float rTemp = 0;
          
          //initializing variables for transmit
          boolean transmit = false;
          PVector Ptr = PVector.add(Ph,PVector.mult(Tr,0.01)); //offset to prevent self intersects
          PVector tnorm = znorm;
          float trTemp = 0;
                      
          float countr = 0;
          int hitIndex = j;
          float Reff = 0;
          
          if(scene[j].refl) {
            Reff = 1;
            reflect = true;
          }
          if(scene[j].refr) {
            transmit = true;
            reflect = true;
            float Schlick = ((none-ntwo)/(none+ntwo))*((none-ntwo)/(none+ntwo));          
            Reff = Schlick + (1-Schlick)*((1-ci)*(1-ci)*(1-ci)*(1-ci)*(1-ci));
          }

          //while hit point is reflective/transmittive, continue iterating
          while((reflect || transmit) && countr<20) {
            //tracking reflect
            rTemp = Float.MAX_VALUE;
            //initialize reflectColor
            color reflectColor = scene[0].geoColor;
            reflect = false;
            
            //tracking refract
            trTemp = Float.MAX_VALUE;
            //initialize refractColor
            color refractColor = scene[0].geoColor;
            transmit = false;

            for(int k = 0; k<scene.length; k = k+1) {
              if(k!=hitIndex) { //ignoring self intersects for now
                //reflect intersect calculation
                if (rTemp>scene[k].intersect(Re,Pr,scene[k].rad) && Reff>0) {
                  rTemp = scene[k].intersect(Re,Pr,scene[k].rad);
                  rnorm = PVector.mult(scene[k].normal(Pr),1);
                  reflectColor = scene[k].geoColor;
                  hitIndex = k;
                  reflect = false;
                  transmit = false;
                  
                  //if hit reflective object
                  if(scene[k].refl){
                    reflect = true;
                    hitIndex = k;
                  }
                  //if hit transmittive object
                  if(scene[k].refr){                      
                    transmit = true;
                    hitIndex = k;
                  }
                }
                //transmit intersect calculation
                if (trTemp>scene[k].intersect(Tr,Ptr,scene[k].rad) && (1-Reff)>0) {
                  trTemp = scene[k].intersect(Tr,Ptr,scene[k].rad);
                  tnorm = scene[k].normal(Ptr);
                  refractColor = scene[k].geoColor;
                  transmit = false;
                  reflect = false;
                  
                  if(Float.isNaN(Tr.x)) {
                    print(Tr,"\n");
                  }
                  
                  //if hit reflective object
                  if(scene[k].refl){
                    reflect = true;
                    hitIndex = k;
                  }
                  //if hit refractive object
                  if(scene[k].refr){
                    transmit = true;
                    hitIndex = k;
                  }
                }
              }
            }
            if(rTemp==Float.MAX_VALUE) {
              reflectColor = color(0,0,0);
              if(modeM.equals("envMap")) {
                float sphereX = n2.dot(PVector.mult(Re,-1));
                float sphereY = n0.dot(PVector.mult(Re,-1));
                float sphereZ = n1.dot(PVector.mult(Re,-1));
                //u,v conversion
                float v = acos(sphereZ)/PI;
                float u = acos(sphereY/sqrt(1-sphereZ*sphereZ))/(2*PI);
                //flip u if x negative
                if(sphereX<0) {
                  u = 1-u;
                }
                float tS0 = envimg.width;
                float tS1 = envimg.height;
                reflectColor = envimg.get((int)(u*tS0),(int)(v*tS1));
              }
            }
            if(trTemp==Float.MAX_VALUE) {
              refractColor = color(0,0,0);
              if(modeM.equals("envMap")) {
                float sphereX = n2.dot(PVector.mult(Tr,-1));
                float sphereY = n0.dot(PVector.mult(Tr,-1));
                float sphereZ = n1.dot(PVector.mult(Tr,-1));
                //u,v conversion
                float v = acos(sphereZ)/PI;
                float u = acos(sphereY/sqrt(1-sphereZ*sphereZ))/(2*PI);
                //flip u if x negative
                if(sphereX<0) {
                  u = 1-u;
                }
                float tS0 = envimg.width;
                float tS1 = envimg.height;
                refractColor = envimg.get((int)(u*tS0),(int)(v*tS1));
              }
            }
            
            //recalculate hitpoint for recursion
            Pr = PVector.add(Pr, PVector.mult(Re,rTemp));
            Ptr = PVector.add(Ptr,PVector.mult(Tr,trTemp));
            
            if(reflect) {
              Re = (PVector.mult(scene[hitIndex].normal(Pr),2*scene[hitIndex].normal(Pr).dot(PVector.mult(Re,-1)))).add(Re);
              Re.normalize();
              Pr = Pr.mult(1.01);
              Reff = 1;
            }
            if(transmit) {
              light = PVector.sub(Ptr,lightPos);
              light.normalize();
              ci = max(tnorm.dot(light),0);
              c2 = sqrt(1 - n*n*(1-ci*ci));
              //if hitting back face, flip
              if (light.dot(tnorm)>0) {
                //n = 1/n;
                ci = max(tnorm.dot(light),0);
                c2 = sqrt(1 - n*n*(1-ci*ci));
                tnorm = tnorm.mult(-1);
              }
              Tr = PVector.mult(light, n).add(PVector.mult(tnorm,((n)*ci-c2)));
              Tr.normalize();
              Reff = 0;
            }
            objectColor = color(Reff*red(reflectColor)+(1-Reff)*red(refractColor),Reff*green(reflectColor)+(1-Reff)*green(refractColor),Reff*blue(reflectColor)+(1-Reff)*blue(refractColor));
            countr = countr+1;
          }
          //shading reflected point            
          specStrength = 0.9;
          //change hitpoint/norm for refraction
          if(Reff < 1) {
            Pr = Ptr;
            rnorm = tnorm;
            specStrength=0;
          }
          light = PVector.sub(Pr,lightPos);
          light.normalize();
          
          float cosTheta=0;
          float lightCounter=0;
          spec = 0;
          diffuse = 0;
          //shading & shadow casting for each light in scene (only 1 right now)
          for (int l=0; l<lightSet.length; l=l+1) {
            light = lightSet[l].lightRay(Pr);
            light.normalize();
      
            //shading calculation
            diffuse = diffuse + max(rnorm.dot((PVector.mult(light, -1))), 0);
            spec = spec + specHighlight(rnorm, light, PVector.mult(Re, -1), specStrength);
      
            //shadow calculation
            cosTheta = cosTheta+lightSet[l].castShadow(Pr, lightSet[l].lightRay(Pr), mode, scene, rnorm);
            lightCounter = lightCounter+1;
          }
      
          //average values for multiple lights
          diffuse = diffuse/lightCounter;
          spec = spec/lightCounter;
          cosTheta = cosTheta/lightCounter;
          
          if((rTemp == Float.MAX_VALUE && reflect) || (trTemp == Float.MAX_VALUE && transmit)) {
            cosTheta = 0;
            diffuse = 1-ambient;
          }
          
          if(modeM.equals("refract")) {
            cosTheta=0;
          }
                                           
          if (cosTheta>0) {
            spec = 0;
          }
      
          //shadow calculated color, blocked for env map
          if(!modeM.equals("envMap")){
            objectColor = color((1-cosTheta)*(red(objectColor))+cosTheta*47, (1-cosTheta)*(green(objectColor))+cosTheta*27, (1-cosTheta)*(blue(objectColor))+cosTheta*24);
            objectColor = color(spec*red(specLightColor)+(1-spec)*(diffuse)*red(objectColor), spec*green(specLightColor)+(1-spec)*(diffuse)*green(objectColor), spec*blue(specLightColor)+(1-spec)*(diffuse)*blue(objectColor));
          }
        }
      }
    }
    
    light = PVector.sub(Ph,lightPos);
    light.normalize();

    
    float cosTheta=0;
    float lightCounter=0;
    specStrength = 0.9;
    spec = 0;
    diffuse = 0;
    //shading & shadow casting for each light in scene (only 1 right now)
    for (int j=0; j<lightSet.length; j=j+1) {
      light = lightSet[j].lightRay(Ph);
      light.normalize();

      //shading calculation
      diffuse = diffuse + max(znorm.dot((PVector.mult(light, -1))), 0);
      spec = spec + specHighlight(znorm, light, PVector.mult(nPe, -1), specStrength);

      //shadow calculation
      cosTheta = cosTheta+lightSet[j].castShadow(Ph, lightSet[j].lightRay(Ph), mode, scene, znorm);
      lightCounter = lightCounter+1;
    }

    //average values for multiple lights
    diffuse = diffuse/lightCounter;
    spec = spec/lightCounter;
    cosTheta = cosTheta/lightCounter;

    //shadow calculated color
    objectColor = color((1-cosTheta)*(red(objectColor))+cosTheta*47, (1-cosTheta)*(green(objectColor))+cosTheta*27, (1-cosTheta)*(blue(objectColor))+cosTheta*24);

    //block spec if in shadow
    if (cosTheta>0) {
      spec = 0;
    }
  
    //environment map
    if (tTemp==Float.MAX_VALUE) {
      objectColor = darkColor = color(0,0,0);
      if(modeM.equals("envMap")) {
        float sphereX = n2.dot(PVector.mult(nPe,-1));
        float sphereY = n0.dot(PVector.mult(nPe,-1));
        float sphereZ = n1.dot(PVector.mult(nPe,-1));
        //u,v conversion
        float v = acos(sphereZ)/PI;
        float u = acos(sphereY/sqrt(1-sphereZ*sphereZ))/(2*PI);
        //flip u if x negative
        if(sphereX<0) {
          u = 1-u;
        }
        float tS0 = envimg.width;
        float tS1 = envimg.height;
        objectColor = darkColor = envimg.get((int)(u*tS0),(int)(v*tS1));
      }
      spec = 0;
      cosTheta = 0;
      diffuse = 1-ambient;
    }
    
    nimg.pixels[i] = color(spec*red(specLightColor)+(1-spec)*(diffuse+ambient)*red(objectColor), spec*green(specLightColor)+(1-spec)*(diffuse+ambient)*green(objectColor), spec*blue(specLightColor)+(1-spec)*(diffuse+ambient)*blue(objectColor));
    rawColor[(int) x][(int) y] = nimg.pixels[i];
  }
  
  imgPath = "./images/rasterized.png";
  img = createImage(500, 500, RGB);
  totalImageSize = img.width * img.height;
  img.loadPixels();

  for (int i = 0; i < totalImageSize-img.width; i = i + 1) {
    x = (i+1)%img.width;
    y = img.height-floor(((i+1)/img.width))-1;

    int subX = 4*(int) x;
    int subY = 4*(int) y;
    color avg = rawColor[subX][subY];
    float avgr = red(avg);
    float avgg = green(avg);
    float avgb = blue(avg);

    for (int k=1; k<4; k = k+1) {
      for (int j=1; j<4; j = j+1) {
        color a = rawColor[subX+j][subY+k];
        float red = red(a);
        float green = green(a);
        float blue = blue(a);
        avgr = (avgr + red);
        avgg = (avgg + green);
        avgb = (avgb + blue);
      }
    }
    img.pixels[(i)] = color(avgr/9, avgg/9, avgb/9);
  }

  image(img, 0, 0);
  save(imgPath);
  noLoop();
}
