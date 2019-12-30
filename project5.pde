//
// VIZA 656: Digital Image
// Name: Anne Fu
// Project 5
//

float specHighlight (PVector znorm, PVector light, PVector viewer,PVector Pl, PVector Phit, float specStrength) {
  PVector reflect = (PVector.mult(znorm, 2.0*(znorm.dot(PVector.mult(light, -1))))).add(light);
  reflect.normalize();
  float specLight = ((max(viewer.dot(reflect), 0.0)-0.9)/(1-0.9))*specStrength;
  if (specLight>0.3) {
    specLight = 1;
  } else if (specLight<0) {
    specLight = 0;
  }
  
  //for rectangular highlight
  if (rectSpec) {
    reflect = PVector.sub(viewer, PVector.mult(znorm,2*viewer.dot(znorm)));
    
    float localX = PVector.sub(Phit,Pl).dot(PVector.sub(n0Light, PVector.mult(n2Light, reflect.dot(n0Light)/reflect.dot(n2Light))));
    float localY = PVector.sub(Phit,Pl).dot(PVector.sub(n1Light, PVector.mult(n2Light, reflect.dot(n1Light)/reflect.dot(n2Light))));
    
    if(max(abs(localX)/s0,abs(localY)/s1)<1 && reflect.z<0) {
      float intensity = 1-max(abs(localX)/s0,abs(localY)/s1);
      intensity = (intensity-0.3)/(1-0.3);
      if(intensity>0.7) {
        intensity = 1;
      }
      else if(intensity<0) {
        intensity = 0;
      }
      specLight = 1*intensity;
    }
    else {
      specLight = 0;
    }
  }
  return specLight;
}

PImage img;
PImage nimg;
float x;
float y;
float sX, sY;
int pixelx;
int pixely;
float z;
color objectColor;
color specLightColor;
PVector znorm;
PVector light;
PVector lightPos;
PVector Pe;
PVector Ph;
PVector vView;
PVector vUp;
PVector n0, n1, n2;
PVector sphereCenter;
PVector n0Light;
PVector n1Light;
PVector n2Light;
PVector nL;
int totalImageSize;
float ambient;
float diffuse;
float specStrength;
float spec;
float spotAngle;
float s0, s1;

color[][] rawColor;
String imgName;
String imgPath;
String nimgPath;
boolean rectSpec;

void setup() {
  size(500, 500);

  lightPos = new PVector(-20, 50, 40);
  
  //spotlight
  nL = new PVector(1,-1,-1);
  nL.normalize();
  spotAngle = 0.2;
  
  //arealight
  s0 = 40;
  s1 = 30;

  n0Light = new PVector(1, 0, 0);
  n1Light = new PVector(0, 1, 0);
  n2Light = new PVector(0,0,-1);

  ambient = 0.3;
  specStrength = 0.7;
  
  //set up camera
  Pe = new PVector(0, 10, 200);
  vView = new PVector(0, 0, 1);
  vView.normalize();
  vUp = new PVector (0, 1, 0);
  n2 = vView;
  n0 = (vView.cross(vUp)).normalize();
  n1 = n0.cross(n2);

  specLightColor = color(255, 247, 144);
  rectSpec = false; 
}

void mouseClicked() {
  image(nimg, 0, 0);
  lightPos.x = mouseX-img.width/2;
  lightPos.y = -(mouseY-img.width/2);
  redraw();
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
    } else if (keyCode == SHIFT) {
      rectSpec = !rectSpec;
    }
  }
  redraw();
}

void draw() {

  //create planes  
  PVector plane0 = new PVector(0, -30, 0);
  PVector nPlane0 = new PVector(0, 1, 0);
  nPlane0.normalize();

  PVector plane1 = new PVector(40, 0, 0);
  PVector nPlane1 = new PVector(-1, 0, 0);
  nPlane1.normalize();

  PVector plane2 = new PVector(0, 0, -80);
  PVector nPlane2 = new PVector(0, 0, 1);
  nPlane2.normalize();

  //create sphere
  sphereCenter = new PVector(0, 10, 0);
  float r = 35;
  //spherical border
  float rBorder = 37;

  //assign colors
  color sphere0Color = color (118, 161, 177);
  color plane2Color = color(160, 179, 124);
  color plane1Color = color(242, 215, 72);
  color plane0Color = color(191, 47, 39);

  nimgPath = "./images/output.png";
  nimg = createImage(2000, 2000, RGB);
  totalImageSize = nimg.width * nimg.height;
  nimg.loadPixels();

  rawColor = new color[nimg.width][nimg.height];

  for (int i = 0; i < totalImageSize-1; i = i + 1) {

    x = (i+1)%nimg.width;
    y = nimg.height-floor(((i+1)/nimg.width))-1;
    pixelx = (int) x;
    pixely = floor(((i+1)/nimg.width));

    //create ray
    PVector nPe = new PVector((x-nimg.width/2)/13-Pe.x, (y-nimg.height/2)/13-Pe.y, 100-Pe.z);
    nPe.normalize();

    //solve for plane
    float th0 = (nPlane0.dot(PVector.sub(plane0, Pe)))/(nPlane0.dot(nPe));
    float th1 = (nPlane1.dot(PVector.sub(plane1, Pe)))/(nPlane1.dot(nPe));
    float th3 = (nPlane2.dot(PVector.sub(plane2, Pe)))/(nPlane2.dot(nPe));
    //solve for sphere
    float b = (nPe.dot(PVector.sub(sphereCenter, Pe)));
    float c = (PVector.sub(sphereCenter, Pe)).dot(PVector.sub(sphereCenter, Pe))-r*r;
    float cBorder = (PVector.sub(sphereCenter, Pe)).dot(PVector.sub(sphereCenter, Pe))-rBorder*rBorder;
    float th2 = b - sqrt(b*b-c);
    //for sphere border
    float thBorder = b - sqrt(b*b-cBorder);

    if (!Float.isNaN(th2) && (th2<th0 || th0<0) && (th2<th1 || th1<0) && (th2<th3 || th3<0)) {
      //calculate sphere normal
      Ph = PVector.add(Pe, PVector.mult(nPe, th2));
      znorm = PVector.sub(Ph, sphereCenter);
      objectColor = sphere0Color;
      //print(th0,th2,"\n");
    } else if ((th2>th0 || Float.isNaN(th2)) && (th0<th1 || th1<0) && (th0<th3 || th3<0) && th0>0) {
      Ph = PVector.add(Pe, PVector.mult(nPe, th0));
      znorm = nPlane0;
      objectColor = plane0Color;
    } else if ((th0>th1 || th0<0) && (th2>th1 || Float.isNaN(th2)) && (th1<th3 || th3<0) && th1>0) {
      Ph = PVector.add(Pe, PVector.mult(nPe, th1));
      znorm = nPlane1;
      objectColor = plane1Color;
    } else if ((th0>th3 || th0<0) && (th1>th3 || th1<0) && (th2>th3 || Float.isNaN(th2)) && th3>0) {
      Ph = PVector.add(Pe, PVector.mult(nPe, th3));
      znorm = nPlane2;
      objectColor = plane2Color;
    } else {
      znorm = new PVector(0, 0, 0);
      Ph = lightPos;
      objectColor = color(0, 0, 0);
    }
    
    light = PVector.sub(Ph, lightPos);
    
    light.normalize();
    znorm.normalize();

    diffuse = max(znorm.dot((PVector.mult(light, -1))), 0);
    spec = specHighlight(znorm, light, PVector.mult(nPe, -1), lightPos, Ph, specStrength);

    //check if outside spotlight
    if (nL.dot(light)<spotAngle) {
      diffuse=0;
      spec=0;
    }
    else {
      diffuse = diffuse*((nL.dot(light)-spotAngle));
      spec = spec*((nL.dot(light)-spotAngle));
    }
    //border color
    if (!Float.isNaN(thBorder) && Float.isNaN(th2)) {
      spec = 0;
      diffuse = 1-ambient;
      objectColor = specLightColor;
    }
    nimg.pixels[i] = color(spec*red(specLightColor)+(1-spec)*(diffuse+ambient)*red(objectColor), spec*green(specLightColor)+(1-spec)*(diffuse+ambient)*green(objectColor), spec*blue(specLightColor)+(1-spec)*(diffuse+ambient)*blue(objectColor));

    rawColor[(int) x][(int) y] = nimg.pixels[i];
  }

  //resize for rasterization
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

  //Draw the image to the screen at coordinate (0,0)
  image(img, 0, 0);
  save(imgPath);
  noLoop();
}
