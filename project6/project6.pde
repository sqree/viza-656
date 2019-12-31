//
// VIZA 656: Image Synthesis
// Name: Anne Fu
// Project 6
//

float specHighlight (PVector znorm, PVector light, PVector viewer, float specStrength) {
  PVector reflect = (PVector.mult(znorm,2.0*(znorm.dot(PVector.mult(light,-1))))).add(light);
  reflect.normalize();
  float specLight = ((max(viewer.dot(reflect), 0.0)-0.9)/(1-0.9))*specStrength;
  if (specLight>1) {
    specLight = 1;
  }
  else if (specLight<0) {
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
  float castShadow (PVector hitPoint, PVector lightRay, String sMode, geometry[] sceneGeo) {
    PVector rayL = lightRay;
    float cosTheta = 0;
    if(lightType.equals("direct")) {
      Pl = new PVector(hitPoint.x-lightRay.x*400, hitPoint.y-lightRay.y*400, hitPoint.z-lightRay.z*400);      
    }
    //check for shadow        
    float th = sceneGeo[0].intersect(rayL, Pl, sceneGeo[0].rad);
    for(int j=1; j<scene.length; j=j+1) {
      if (th>sceneGeo[j].intersect(rayL,Pl,sceneGeo[j].rad)) {
        th = sceneGeo[j].intersect(rayL,Pl,sceneGeo[j].rad);
      }
    }
    if((hitPoint.x-Pl.x)/rayL.x>th+0.5 || (lightType.equals("spot") && normL.dot(rayL)<lAngle)) {
      switch(sMode) {
        case "modedR":
        //d/R calculation
        PVector Pray = hitPoint;
        float tempT = 1;
        float dR = 1;
        //trace along light ray
        while (PVector.sub(Pray, Pl).mag() > 1) {
          Pray = PVector.sub(hitPoint, PVector.mult(rayL,tempT));
          tempT = tempT + 1;
          for (int j=2; j<5; j = j+1) {
            //increment if inside sphere -- I should probably use the implicit function for these but since there are only planes and spheres in this scene...
            if(PVector.sub(Pray, scene[j].centerPoint).mag() < scene[j].rad) {
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
        float randomOffsetM = random(0,lM-1); //randomizing sublight sequence
        float randomOffsetN = random(0,lN-1);
        for (int j=0; j<lM; j = j+1) {
          for (int k=0; k<lN; k = k+1) {
            //X, Y coordinates essentially simplify to randomized sublight section + random offset for random point within sublight section
            tempLightPos = PVector.add(PVector.add(tempLightPos, PVector.mult(lN0,lS0*(((j+randomOffsetM)% lM)/lM) + random(0,lS0/lM))),PVector.mult(lN1,lS1*(((k+randomOffsetN)%lN)/lN)+ random(0,lS1/lN)));
            tempLight = PVector.sub(hitPoint, tempLightPos);
            tempLight.normalize();
            //check for intersection between light ray and objects in scene
            th = scene[0].intersect(tempLight, tempLightPos, scene[0].rad);
            for(int l=1; l<scene.length; l=l+1) {
              if (th>scene[l].intersect(tempLight,tempLightPos,scene[l].rad)) {
                th = scene[l].intersect(tempLight,tempLightPos,scene[l].rad);
              }
            }
            //increment shadow counter if intersection found
            if((hitPoint.x-tempLightPos.x)/tempLight.x>th+0.5) {
              tempCount = tempCount+1;
            }
          }
        }
        //normalize shadow
        cosTheta = znorm.dot(PVector.mult(rayL, -1))*tempCount/(lM*lN);
        break;
        default:
        cosTheta = znorm.dot(PVector.mult(rayL, -1));
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
  geometry (PVector position, PVector planeNormal, String type, color objectColor, color darkColor, float radius) {
    centerPoint = position; 
    geoType = type;
    geoColor = objectColor;
    geoDarkColor = darkColor;
    switch(geoType){
      case "sphere":
      rad = radius;
      pNormal = new PVector(0,0,0);
      break;
      case "plane":
      rad = Float.NaN;
      pNormal = planeNormal;
      break;
      default:
      break;
    }
  } 
  PVector normal (PVector hitPoint) {
    PVector geoNorm = new PVector(0,0,0);
    switch(geoType){
      case "sphere":
        geoNorm = PVector.sub(hitPoint, centerPoint);
      break;
      case "plane":
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
    switch(geoType){
      case "sphere":
        float b = (ray.dot(PVector.sub(centerPoint,startPoint)));
        float c = (PVector.sub(centerPoint,startPoint)).dot(PVector.sub(centerPoint,startPoint))-radius*radius;
        intersectPoint = b - sqrt(b*b-c);
      break;
      case "plane":
        intersectPoint = (pNormal.dot(PVector.sub(centerPoint,startPoint)))/(pNormal.dot(ray));
      break;
      default:
      break;
    }
    if(Float.isNaN(intersectPoint) || intersectPoint<0) { //check for invalid intersection
      intersectPoint = Float.MAX_VALUE;
    }
    return intersectPoint;
  }
}

PImage img, nimg;
float x, y, z;
float sX, sY;
float th;
color lightColor, darkColor, objectColor, specLightColor;
PVector znorm;
PVector light, lightPos;
PVector nL; //for spot, area light
float spotAngle; //for spot light
PVector Pe, Ph, P0;
PVector vView, vUp;
PVector n0, n1, n2;
int totalImageSize;
float ambient, diffuse, specStrength, spec;
color[][] rawColor;
String mode; //shadow calculation mode hard, d/R, area
String lightT; //point, spot, direct

String imgName;
String imgPath;
String nimgPath;

geometry[] scene;
light[] lightSet;

void setup() {
  size(500,500);  
  //camera set up
  Pe = new PVector(0, 0, 200);
  vView = new PVector(0,0,1);
  vView.normalize();
  vUp = new PVector (0,1,0);
  n2 = vView.normalize();
  n0 = (vView.cross(vUp)).normalize();
  n1 = n0.cross(n2);
  P0 = PVector.add(Pe,PVector.mult(n2, 50));
  
    
  lightPos = new PVector(150,100,30);
  //print(lightPos,"\n");
  
  //light setup
  nL = new PVector(-0.5,-1,0.5); //for direct, spot, area
  nL.normalize();
  spotAngle = 0.85;
  ambient = 0.2;
  specStrength = 0.7;
  specLightColor = color(255,247,170);
  
  mode = "";
  lightT = "";
}

void mouseClicked() {
  image(nimg,0,0);
  lightPos.x = mouseX-img.width/2;
  lightPos.y = -(mouseY-img.width/2);
  lightPos.z = 150;
  nL = new PVector(-(mouseX-img.width/2), (mouseY-img.width/2), -200);
  nL.normalize();
  redraw();
}

void keyPressed() {
  if (key == CODED) {
    if (keyCode == UP) {
      Pe.y = Pe.y - 5;
    } 
    else if (keyCode == DOWN) {
      Pe.y = Pe.y + 5;
    }
        else if (keyCode == RIGHT) {
      Pe.x = Pe.x - 5;
    }
    else if (keyCode == LEFT) {
      Pe.x = Pe.x + 5;
    }
    else if (keyCode == ALT) {
      Pe.z = Pe.z - 5;
    }
    else if (keyCode == CONTROL) {
      Pe.z = Pe.z + 5;
    }
  }
  else if (key == 'd') {
    mode = "modedR";
    if(lightT.equals("spot")) {
      lightT = "point";
    }
  }
  else if (key == 'a') {
    mode = "modeArea";
    lightT = "point";
  }
  else if (key == 'h') {
    mode = "";
  }
  else if (key == 'l') {
    lightT = "direct";
    mode = "";
  }
  else if (key == 's') {
    lightT = "spot";
    mode = "";
  }
  redraw();
}

void draw() {
  
  PVector nPlane0 = new PVector(0,1,0);
  nPlane0.normalize();
    
  PVector nPlane1 = new PVector(0,0,1);
  nPlane1.normalize();
      
  scene = new geometry[5];
  scene[0] = new geometry(new PVector(0,-30,0), nPlane0, "plane", color(147,188,190), color(134,167,174), 0);
  scene[1] = new geometry(new PVector(0,0,-50), nPlane1, "plane", color(210,215,218), color(158,163,157), 0);
  scene[2] = new geometry(new PVector(90,10,0), new PVector(0,0,0), "sphere", color(152, 165, 73), color(78,38,43), 40);
  scene[3] = new geometry(new PVector(20,10,30), new PVector(0,0,0), "sphere", color(152, 165, 73), color(78,38,43), 40);
  scene[4] = new geometry(new PVector(-50,10,60), new PVector(0,0,0), "sphere", color(152, 165, 73), color(78,38,43), 40);
  
  lightSet = new light[1];
  lightSet[0] = new light(lightT, lightPos, new PVector(0,1,0), nL, 0.85, 40.0, 40.0, 4.0, 4.0, color(255,247,170));
    
  nimgPath = "./images/output.png";
  nimg = createImage(2000,2000,RGB);
  totalImageSize = nimg.width * nimg.height;
  nimg.loadPixels();
  
  rawColor = new color[nimg.width][nimg.height];

  for (int i = 0; i < totalImageSize-1; i = i + 1) {

    x = (i+1)%nimg.width;
    y = nimg.height-floor(((i+1)/nimg.width))-1;
    
    //shoot ray
    PVector nPe = new PVector((x-nimg.width/2-Pe.x)/5-P0.x, (y-nimg.height/2-Pe.y)/5-P0.y, Pe.z/5-P0.z);
    nPe.normalize();
    
    //initialize for hitpoint calculation
    float tTemp = scene[0].intersect(nPe, Pe, scene[0].rad);
    Ph = PVector.add(Pe, PVector.mult(nPe, tTemp));
    znorm = scene[0].normal(Ph);
    objectColor = scene[0].geoColor;
    darkColor = scene[0].geoDarkColor;
    for(int j=1; j<scene.length; j=j+1) {
      if (tTemp>scene[j].intersect(nPe,Pe,scene[j].rad)) {
        tTemp = scene[j].intersect(nPe,Pe,scene[j].rad);
        Ph = PVector.add(Pe, PVector.mult(nPe,tTemp));
        znorm = scene[j].normal(Ph);
        objectColor = scene[j].geoColor;
        darkColor = scene[j].geoDarkColor;
      }
    }
    znorm.normalize();
    float cosTheta=0;
    float lightCounter=0;
    spec = 0;
    diffuse = 0;
    //shading & shadow casting for each light in scene (only 1 right now)
    for(int j=0; j<lightSet.length; j=j+1) {
      specLightColor = lightSet[j].colorL;
      light = lightSet[j].lightRay(Ph);
      light.normalize();
      
      //shading calculation
      diffuse = diffuse + max(znorm.dot((PVector.mult(light,-1))),0);
      spec = spec + specHighlight(znorm,light,PVector.mult(nPe,-1),specStrength);
      
      //shadow calculation
      cosTheta = cosTheta+lightSet[j].castShadow(Ph, lightSet[j].lightRay(Ph), mode, scene);
      lightCounter = lightCounter+1;
    }
    
    //average values for multiple lights
    diffuse = diffuse/lightCounter;
    spec = spec/lightCounter;
    cosTheta = cosTheta/lightCounter;
    
    //shadow calculated color
    objectColor = color((1-cosTheta)*(red(objectColor))+cosTheta*68,(1-cosTheta)*(green(objectColor))+cosTheta*61,(1-cosTheta)*(blue(objectColor))+cosTheta*105);
    
    //block spec if in shadow
    if(cosTheta>0) {
      spec = 0;
    }
    
    //final color mixing
    nimg.pixels[i] = color((spec*red(specLightColor)+(1-spec)*(diffuse+ambient)*red(objectColor)+(1-diffuse-ambient)*red(darkColor)),(spec*green(specLightColor)+(1-spec)*(diffuse+ambient)*green(objectColor)+(1-diffuse-ambient)*green(darkColor)),(spec*blue(specLightColor)+(1-spec)*(diffuse+ambient)*blue(objectColor)+(1-diffuse-ambient)*blue(darkColor)));
    rawColor[(int) x][(int) y] = nimg.pixels[i];
  }
  
  imgPath = "./images/rasterized.png";
  img = createImage(500,500,RGB);
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

    for (int k=1; k<4; k = k+1){
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
    img.pixels[(i)] = color(avgr/9,avgg/9,avgb/9);
  }

  image(img,0,0);
  save(imgPath);
  noLoop();
}
