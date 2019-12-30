//
// VIZA 656: Digital Image
// Name: Anne Fu
// Project 4
//

PImage antialias(PImage vimg) {
  //antialiasing
  int totalImageSize = vimg.width * vimg.height;
  for (int t=1; t<3 ; t = t+1){
    for (int i = 0; i < totalImageSize; i = i + 1) {
      x = (i+1)%vimg.width;
      y = floor(((i+1)/vimg.width));
      color avg = vimg.get(x, y);
      float avgr = red(avg);
      float avgg = green(avg);
      float avgb = blue(avg);
      //get color values of the pixels surrounding
      for(int j = -2; j < 2; j = j+1) {
        for(int k = -2; k < 2; k = k+1) {
          color a = vimg.get(x+j, y+k);
          float r = red(a);
          float g = green(a);
          float b = blue(a);
          if (abs(avgr-r) > 10 || abs(avgg - g) > 10 || abs(avgb - b) > 10){
            avgr = sqrt((avgr*avgr + r*r)/2);
            avgg = sqrt((avgg*avgg + g*g)/2);
            avgb = sqrt((avgb*avgb + b*b)/2);
          }
        }
      }
      vimg.pixels[i] = color(avgr,avgg,avgb);
    }
  }
  vimg.updatePixels();
  return vimg;
}

color blur(int x, int y, int glossy, PImage vimg) {
  color avg = vimg.get(x, y);
  float avgr = red(avg);
  float avgg = green(avg);
  float avgb = blue(avg);

  //get color values of the pixels surrounding
  for(int j = -1*glossy; j < glossy; j = j+1) {
    for(int k = -1*glossy; k < glossy; k = k+1) {
        int blurx = (x+floor(random(0,5*glossy)));
        int blury = (y+floor(random(0,5*glossy)));
        color a = vimg.get(blurx, blury);
        float r = red(a);
        float g = green(a);
        float b = blue(a);
        avgr = sqrt((avgr*avgr + r*r)/2);
        avgg = sqrt((avgg*avgg + g*g)/2);
        avgb = sqrt((avgb*avgb + b*b)/2);
    }
  }
  return color(floor(avgr),floor(avgg), floor(avgb));
}

PVector inshadow(int x, int y, float[][] heightValues, PVector light, int maxx, int maxy, float[][] hAngle) {
  float pHeight = heightValues[x][y];
  int currentx = x;
  int currenty = y;
  int t = 1;
  int hCount = 1;
  int shade = 0;
  float R = pHeight;
  float shadowRatio = hAngle[x][y];
  while(currentx>0 && currentx<maxx && currenty>0 && currenty<maxy){
    currentx = (int) (currentx - light.x*t+0.5);
    currenty = (int) (currenty - light.y*t+0.5);
    pHeight = (pHeight - light.z*t);
    if(currentx>0 && currentx<maxx && currenty>0 && currenty<maxy) {
      if (pHeight < heightValues[currentx][currenty]) {
        shade=1;
        R = R-light.z*t;        
      }
    }
    t = t+1;
  }
  shadowRatio = (10/R);
  if (shadowRatio>1) {
    shadowRatio = 1;
  }
  else if (shadowRatio<0.9) {
    shadowRatio = 0.8;
  }
  PVector shadowInfo = new PVector(shade,shadowRatio);
  return shadowInfo;
}

float specHighlight (PVector znorm, PVector light, PVector viewer, float specStrength) {
  reflect = (PVector.mult(znorm,2.0*(znorm.dot(PVector.mult(light,-1))))).add(light);
  reflect.normalize();
  float specLight = pow(max(viewer.dot(reflect), 0.0), 32)*specStrength;
  //float specLight = (max(viewer.dot(reflect), 0.0)-0.95)*20;
  if (specLight>1) {
    specLight = 1;
  }
  else if (specLight<0) {
    specLight = 0;
  }
  return specLight;
}

PImage img;
PImage nimg;
PImage diffuseDark;
PImage diffuseLight;
PImage backgroundImage;
int x;
int y;
int pixelx;
int pixely;
int glossy;
float z;
color lightColor;
color objectColor;
color imgColor;
color reflectColor;
color refractColor;
color specLightColor;
PVector norm;
PVector znorm;
PVector light;
PVector viewer;
PVector reflect;
PVector gradient;
PVector blah;
PVector bg;
PVector Re;
PVector transmit;
float n;
int h;
int totalImageSize;
float d;
float ambient;
float diffuse;
float specStrength;
float spec;
float zratio;
float gradientx, gradienty;
float gradientxz,gradientyz,gradientxo,gradientyo;
float depthRatio;
float n1,n2;
float ci, c2;
float[][] zValues;
float[][] cosTheta;
PFont f;
String input = "";
String userinput = "";
boolean mapType;
boolean objType;
boolean done;
String mapTemp;
String imgName;
String imgPath;
String nimgPath;
String [] imgNameParse;

void setup() {
  size(800,800);
  f = createFont("Arial",16);
  mapType = false;
  objType = false;
  done = false;
  
  light = new PVector(1,-1,-1);
  light.normalize();
  viewer = new PVector(0,0,1);
  viewer.normalize();
  specStrength = 0.7;
  ambient = 0.7; 
  bg = new PVector(0,0,0);
  lightColor = color(239, 100, 74);
  specLightColor = color(248,193,70);
  //distance
  d = 256;
  //glossiness
  glossy = 0;
  n1 = 1;
  n2 = 1.5;
  n = n1/n2;
}

void mousePressed() {
  if (mouseButton == LEFT) {
  } 
  else if (mouseButton == RIGHT) {
    mapType = false;
    objType = false;
    done = false;
    redraw();
  } 
}

void mouseClicked() {
  if(done) {
  image(nimg,0,0);
  light.x = nimg.width/2-mouseX;
  light.y = -(nimg.height/2-mouseY);
  light.z = -mouseX%(nimg.width/2);
  light.normalize();
  input = userinput;
  keyPressed();
  }
}

void keyPressed() {
  if (key == CODED) {
    if (keyCode == UP) {
      glossy = glossy + 5;
    } 
    else if (keyCode == DOWN) {
      glossy = glossy - 5;
    }
    input = userinput;
  } 
  if (key != '\n' && key != CODED ) {
    // As long as 'enter' isn't hit
    // Each character typed by the user is added to the end of the String variable.
    input = input + key; 
    loop();
  } 
  else {
    userinput = input;
    // A String can be cleared by setting it equal to ""
    if (input.equals("normal") || input.equals("depth")) {
      mapTemp = input;
      mapType = true;
      loop();
    }
    else if (input.equals("sphere")) {
      userinput = mapTemp + " " + input;
      //objectColor = color(random(0,255),random(0,255),random(0,255));
      objType = true;
    }
    else if ((input.substring(0,5)).equals("image")) {
      userinput = mapTemp + " " + input.substring(0,5);
      imgName = input.substring(6,input.length());
      objType = true;
    }

    //clear input  
    input = ""; 
 
    switch(userinput) {
      
      case "depth image":
        imgPath = "./images/" + imgName;
        imgNameParse = split(imgName, '.');
        nimgPath = "./images/" + imgNameParse[0] + "_shaded." + imgNameParse[1];
        
        img = loadImage(imgPath);
        nimg = createImage(img.width,img.height,RGB);
        
        diffuseDark = loadImage("./images/" + imgNameParse[0] + "1." + imgNameParse[1]);
        diffuseLight = loadImage("./images/" + imgNameParse[0] + "2." + imgNameParse[1]);

    
        // The pixels are stored  in an array that is totalImageSize long.
        totalImageSize = img.width * img.height;
        // Make sure the pixels are stored properly in the array
        img.loadPixels();
        nimg.loadPixels();
        
        zValues = new float[nimg.width][nimg.height];
        cosTheta = new float[nimg.width][nimg.height];
        // depth map
        for (int i = 0; i < totalImageSize-1; i = i + 1) {
          x = (i+1)%img.width;
          y = img.height-floor(((i+1)/img.width))-1;
    
          pixelx = (i+1)%img.width;
          pixely = floor(((i+1)/img.width));
          objectColor = diffuseDark.get(pixelx,pixely);
          lightColor = diffuseLight.get(pixelx,pixely);
          
          z = red(img.get(pixelx,pixely));
          zValues[x][y] = z;
          depthRatio = z/255;
      
          h = 8;
          
          float gradientxavg = 0;
          float gradientyavg = 0;
          float count=0;
          //calculate gradient
          for (int j = 1; j < h; j = j + 1) {
            gradientxz = red(img.get(pixelx-j,pixely));
            gradientyz = red(img.get(pixelx,pixely+j));
            gradientxo = red(img.get(pixelx+j,pixely));
            gradientyo = red(img.get(pixelx,pixely-j));
            gradientx = (gradientxo - gradientxz)/2;
            gradienty = (gradientyo - gradientyz)/2;
            gradientxavg = gradientxavg+gradientx;
            gradientyavg = gradientyavg+gradienty;
            count = count+1;          
          }
          gradientx = gradientxavg/count;
          gradienty = gradientyavg/count;
          PVector gradientz = new PVector(1, 0, gradientx);
          PVector gradiento = new PVector(0, 1, gradienty);
          gradient = gradientz.cross(gradiento);
       
          gradient.normalize();
          znorm = gradient;
          znorm.normalize();
                   
          diffuse = max(znorm.dot((PVector.mult(light,-1))),0);
          
          cosTheta[x][y] = diffuse;
      
          //adding spec
          spec = specHighlight(znorm,light,viewer,specStrength);
          
          nimg.pixels[i] = color(spec*red(specLightColor)+(1-spec)*(ambient+diffuse)*red(lightColor)*red(objectColor)/255, spec*green(specLightColor)+(1-spec)*(ambient+diffuse)*green(lightColor)*green(objectColor)/255, spec*blue(specLightColor)+(1-spec)*(ambient+diffuse)*blue(lightColor)*blue(objectColor)/255);

          //mask out background
          if (z == 0) {
            nimg.pixels[i] = diffuseDark.get(pixelx,pixely);
          }
        }
        for (int i = 0; i < totalImageSize-1; i = i + 1) {
          x = (i+1)%nimg.width;
          y = nimg.height-floor(((i+1)/nimg.width))-1;
          z = zValues[x][y];
          float hitShadow;
          hitShadow = inshadow(x, y, zValues, light, nimg.width, nimg.height, cosTheta).x;
          if (hitShadow==1) {
            nimg.pixels[i] = color(red(nimg.pixels[i])/1.5,green(nimg.pixels[i])/2,blue(nimg.pixels[i])/2);
          }
        }
      
        done = true;
        break;
      
      case "depth sphere":      
        //simple shaded sphere

        nimgPath = "./images/depth_sphere.png";
        nimg = createImage(500,500,RGB);
        
        totalImageSize = nimg.width * nimg.height;
        // Make sure the pixels are stored properly in the array
        nimg.loadPixels();
        
        zValues = new float[nimg.width][nimg.height];
        cosTheta = new float[nimg.width][nimg.height];
        
        for (int i = 0; i < totalImageSize-1; i = i + 1) {
          x = (i+1)%nimg.width;
          y = nimg.height-floor(((i+1)/nimg.width))-1;
          pixelx = x;
          pixely = floor(((i+1)/nimg.width));

          if ((x-floor(nimg.width/2))*(x-floor(nimg.width/2)) + (y-floor(nimg.height/2))*(y-floor(nimg.height/2)) < 10000) {
            z = (sqrt(10000-(x-floor(0.5*nimg.width))*(x-floor(0.5*nimg.width))-(y-floor(0.5*nimg.height))*(y-floor(0.5*nimg.height))));

            zValues[x][y] = z;
            
            znorm = new PVector(x-250,y-250,z);
            znorm.normalize();
      
            objectColor = color(184,43,67);
            
            //adding diffuse
            diffuse = max(znorm.dot((PVector.mult(light,-1))),0);
            
            cosTheta[x][y] = diffuse;
      
            //adding spec
            spec = specHighlight(znorm,light,viewer,specStrength);

            nimg.pixels[i] = color(spec*red(specLightColor)+(1-spec)*(ambient+diffuse)*red(lightColor)*red(objectColor)/255, spec*green(specLightColor)+(1-spec)*(ambient+diffuse)*green(lightColor)*green(objectColor)/255, spec*blue(specLightColor)+(1-spec)*(ambient+diffuse)*blue(lightColor)*blue(objectColor)/255);

          }
          else {
            z = 0;
            cosTheta[x][y] = max(new PVector(0,0,1).dot(PVector.mult(light,-1)),0);
            nimg.pixels[i] = color(134,36,66);
          }
        }
        //color shadow regions 
        for (int i = 0; i < totalImageSize-1; i = i + 1) {
          x = (i+1)%nimg.width;
          y = nimg.height-floor(((i+1)/nimg.width))-1;
          z = zValues[x][y];
          float hitShadow;
          hitShadow = inshadow(x, y, zValues, light, nimg.width, nimg.height, cosTheta).x;
          if (hitShadow==1) {
            nimg.pixels[i] = color(red(nimg.pixels[i])*inshadow(x, y, zValues, light, nimg.width, nimg.height, cosTheta).y,green(nimg.pixels[i])*inshadow(x, y, zValues, light, nimg.width, nimg.height, cosTheta).y,blue(nimg.pixels[i])*inshadow(x, y, zValues, light, nimg.width, nimg.height, cosTheta).y);
          }
        }
        done = true;
        //print("done\n");
        break;
            
    default:
      //redraw();      
      break;
    }
    
    if (!done) {
      loop();
    }
    else {
    //Update when you're done messing with the pixels.
    nimg.updatePixels();
    //run through antialias if using gen sphere
    if (userinput.equals("normal sphere") || userinput.equals("depth sphere")) {
      nimg = antialias(nimg);
    }
    redraw();
    }
  }
}

//Draw is run last by the program. 
void draw() {
  background(255);
  int indent = 25;
  
  // Set the font and fill for text
  textFont(f);
  fill(0);
  //enter type of map to be used, depth/normal
  if(!mapType) {
    text("Click in this window and type. \nHit enter to save. ", indent, 40);
    text("'normal' or 'depth'? : " + input,indent,190);
    noLoop();
  }
  //if generated sphere chosen
  else if(!objType || !done) {
    text("Click in this window and type. \nHit enter to save. ", indent, 40);
    text("'sphere' or 'image'? : " + input,indent,190);
    noLoop();
  }
  else {
    // Draw the image to the screen at coordinate (0,0)
    image(nimg,0,0);
    save(nimgPath);
  }
}
