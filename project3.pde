//
// VIZA 654: Digital Image
// Name: NAME GOES HERE!
// Project X
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
  /*if (avgr == 0) {
        print(x,y,x%vimg.width,y%vimg.height,avg,"\n");
  }*/

  //get color values of the pixels surrounding
  for(int j = -1*glossy; j < glossy; j = j+1) {
    for(int k = -1*glossy; k < glossy; k = k+1) {
        int blurx = (x+floor(random(0,5*glossy)));
        int blury = (y+floor(random(0,5*glossy)));
        color a = vimg.get(blurx, blury);
        float r = red(a);
        float g = green(a);
        float b = blue(a);
      //if (blurx >0 && blurx<vimg.width && blury>0 && blury<vimg.height){
        avgr = sqrt((avgr*avgr + r*r)/2);
        avgg = sqrt((avgg*avgg + g*g)/2);
        avgb = sqrt((avgb*avgb + b*b)/2);
        //println(avgr,avgg,avgb);
      //}
    }
  }
  //print(avgr,avgg,avgb);
  return color(floor(avgr),floor(avgg), floor(avgb));
}

color cubemap_cast(PVector position, PVector pixelPos, PVector ray, float d, PImage cubemap) {
  color castColor = color(0,0,0);
  float x = position.x;
  float y = position.y;
  float z = position.z;
  float testx, testy, testz, t;
  int pixelx = (int)pixelPos.x;
  int pixely = (int)pixelPos.y;
  float dx, dy, dz;
  dx=dy=dz=0;
  testx=testy=testz=0;
  t=0;
  
  //cubemap image  
  //print("in");
  int cubesize = cubemap.width/4;

  while(testx<d&&testy<d&&testz<d) {
    testx = abs(x+ray.x*t);
    testy = abs(y+ray.y*t);
    testz = abs(z+ray.z*t);
    t=t+1;
  }
    
  if(max(testx,testy,testz) == testz){
    if (ray.z<0) {
      castColor = blur(3*cubesize/2-((pixelx)+floor((d-abs(z))*ray.x/ray.z)-floor(bg.x))%(cubesize/2), ((pixely)+floor((d-abs(z))*ray.y/ray.z)+floor(bg.y))%(cubesize/2)+3*cubesize/2, glossy, cubemap);
    }
    else {
      castColor = blur(7*cubesize/2-(((pixelx)+floor((d-z)*ray.x/ray.z)+floor(bg.x))%(cubesize/2)), ((pixely)-floor((d-z)*ray.y/ray.z)-floor(bg.y))%(cubesize/2)+3*cubesize/2, glossy, cubemap);
    }              
  }
  else if(max(testx,testy,testz) == testx){
    if (ray.x<0) {
      castColor = blur(((pixelx)+floor((d-abs(x))*ray.z/ray.x)+floor(bg.x))%(cubesize/2)+cubesize/2, ((pixely)+floor((d-abs(x))*ray.y/ray.x)+floor(bg.y))%(cubesize/2)+3*cubesize/2, glossy, cubemap);
    }
    else {
      castColor = blur(((pixelx)+floor((d-abs(x))*ray.z/ray.x)+floor(bg.x))%(cubesize/2)+5*cubesize/2, ((pixely)-floor((d-abs(x))*ray.y/ray.x)-floor(bg.y))%(cubesize/2)+3*cubesize/2, glossy, cubemap);
    }
  }
  else if(max(testx,testy,testz) == testy){
    if (ray.y<0) {
      castColor = blur(((pixelx)-floor((d+(y))*ray.x/ray.y)-floor(bg.x))%(cubesize/2)+3*cubesize/2, ((pixely)-floor((d+(y))*ray.z/ray.y)-floor(bg.y))%(cubesize/2)+cubesize*5/2, glossy, cubemap);
    }
    else {
      castColor = blur(((pixelx)+floor((d-(y))*ray.x/ray.y)+floor(bg.x))%(cubesize/2)+3*cubesize/2, ((pixely)-floor((d-(y))*ray.z/ray.y)-floor(bg.y))%(cubesize/2)+cubesize/2, glossy, cubemap);
    }
  }
  return castColor;
}


float specHighlight (PVector znorm, PVector light, PVector viewer, float specStrength) {
  reflect = (PVector.mult(znorm,2.0*(znorm.dot(PVector.mult(light,-1))))).add(light);
  reflect.normalize();
  float specLight = pow(max(viewer.dot(reflect), 0.0), 32)*specStrength;
  //print(spec,"\n");
  return specLight;
}

// Declaring a variable of type PImage
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
  // The Size Command can only take in actual numbers, no variables.
  size(800,800);
  // Make a new instance of a PImage by loading an image file
  // Useful methods for the PImage type can be found at https://processing.org/reference/PImage.html
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
  //environment map
  diffuseDark = loadImage("./images/skybox_texture.jpg");
  lightColor = color(255, 255, 255);
  specLightColor = color(255,255,255);
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
  bg.x = nimg.width/2-mouseX;
  bg.y = -(nimg.height/2-mouseY);
  //light.z = mouseX%(nimg.width/2);
  //light.normalize();
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
        //println("generating shaded image based off depth map " + imgName);
        imgPath = "./images/" + imgName;
        imgNameParse = split(imgName, '.');
        //println(imgNameParse[0], imgNameParse[1]);
        nimgPath = "./images/" + imgNameParse[0] + "_shaded." + imgNameParse[1];
        
        img = loadImage(imgPath);
        nimg = createImage(img.width,img.height,RGB);
        
        //diffuseDark = loadImage("./images/" + imgNameParse[0] + "1." + imgNameParse[1]);
        //diffuseLight = loadImage("./images/" + imgNameParse[0] + "2." + imgNameParse[1]);

    
        // The pixels are stored  in an array that is totalImageSize long.
        totalImageSize = img.width * img.height;
        // Make sure the pixels are stored properly in the array
        img.loadPixels();
        nimg.loadPixels();
        // Go through the array pixel by pixel from left to right then top to bottom
        // A 3x3 PImage would be ordered:
        // 1 2 3
        // 4 5 6
        // 7 8 9
     
        // depth map
        for (int i = 0; i < totalImageSize; i = i + 1) {
          x = (i+1)%img.width;
          y = img.height-floor(((i+1)/img.width))-1;
    
          pixelx = (i+1)%img.width;
          pixely = floor(((i+1)/img.width));
          objectColor = diffuseDark.get(pixelx,pixely);
          //lightColor = diffuseLight.get(pixelx,pixely);
          
          z = red(img.get(pixelx,pixely));
          depthRatio = z/255;
      
          h = 3;
          
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
          
          //normal vector
          znorm = gradient;
          znorm.normalize();
          
          //calculate reflect vector
          Re = (PVector.mult(znorm,2*znorm.dot(viewer))).sub(viewer);
          Re.normalize();
          
          //calculate transmit vector            
          ci = max(znorm.dot(light),0);
          c2 = sqrt(1 - n*n*(1-ci*ci));
          
          transmit = PVector.mult(light, n).add(PVector.mult(znorm,((n)*ci-c2)));
          transmit.normalize();
          
          //calculate Fresnel     
          //angle of transmittance
          float ct = abs(znorm.dot(transmit));
          
          //check for total reflection
          float Frpar = ((n2*ci-n1*ct)/(n2*ci+n1*ct))*((n2*ci-n1*ct)/(n2*ci+n1*ct));
          float Frper = ((n1*ct-n2*ci)/(n1*ct+n2*ci))*((n1*ct-n2*ci)/(n1*ct+n2*ci));
          
          float Reff = 1-0.5*(Frpar+Frper);

          PVector temppos = new PVector(x-img.width/2,y-img.height/2,0);
          PVector temppospixel = new PVector(pixelx-img.width/2,img.height-y+1);

          //cast color
          reflectColor = cubemap_cast(temppos, temppospixel, Re, d, diffuseDark);                      
          refractColor = blur(((pixelx)+floor((90-abs(z/10))*transmit.x/transmit.z)+floor(bg.x))+700, ((pixely)-floor((90-abs(z/10))*transmit.y/transmit.z)-floor(bg.y))+512, glossy, diffuseDark);
          objectColor = color(floor((red(reflectColor))*(Reff)+(red(refractColor))*(1-Reff)),floor((green(reflectColor))*(Reff)+(green(refractColor))*(1-Reff)),floor((blue(reflectColor))*(Reff)+(blue(refractColor))*(1-Reff)));

          //adding diffuse
          diffuse = max(znorm.dot((PVector.mult(light,-1))),0);
          //println(i,gradient,gradientz.cross(gradiento));
      
          //adding spec
          spec = specHighlight(znorm,light,viewer,specStrength);
          
          nimg.pixels[i] = color(spec*red(specLightColor)+(1-spec)*(ambient+diffuse)*red(lightColor)*red(objectColor)/255, spec*green(specLightColor)+(1-spec)*(ambient+diffuse)*green(lightColor)*green(objectColor)/255, spec*blue(specLightColor)+(1-spec)*(ambient+diffuse)*blue(lightColor)*blue(objectColor)/255);

          //mask out background
          if (z == 0) {
            nimg.pixels[i] = diffuseDark.get(pixelx-250+diffuseDark.width/2,pixely-250+diffuseDark.height/2);
          }
        }
        done = true;
        break;
      
      case "normal image":
        imgPath = "./images/" + imgName;
        imgNameParse = split(imgName, '.');
        nimgPath = "./images/" + imgNameParse[0] + "_shaded." + imgNameParse[1];
        img = loadImage(imgPath);
        nimg = createImage(img.width,img.height,RGB);
    
        // The pixels are stored  in an array that is totalImageSize long.
        totalImageSize = img.width * img.height;
        // Make sure the pixels are stored properly in the array
        img.loadPixels();
        nimg.loadPixels();
        
        //initialize znorm
        znorm = new PVector(0,0,0);
        
        // normal map
        for (int i = 0; i < totalImageSize; i = i + 1) {
          pixelx = (i+1)%img.width;
          pixely = floor(((i+1)/img.width));
          objectColor = diffuseDark.get(pixelx,pixely);
          //lightColor = diffuseLight.get(pixelx,pixely);

          znorm.x = (red(img.get(pixelx,pixely))/256)*2-1;
          znorm.y = (green(img.get(pixelx,pixely))/256)*2-1;
          znorm.z = (blue(img.get(pixelx,pixely))/256)*2-1;
          
          //normal vector
          znorm.normalize();
          
          //calculate reflect vector
          Re = (PVector.mult(znorm,2*znorm.dot(viewer))).sub(viewer);
          Re.normalize();
          
          //calculate transmit vector            
          ci = max(znorm.dot(light),0);
          c2 = sqrt(1 - n*n*(1-ci*ci));
          
          transmit = PVector.mult(light, n).add(PVector.mult(znorm,((n)*ci-c2)));
          transmit.normalize();
          
          //calculate Fresnel          
          //angle of transmittance
          float ct = abs(znorm.dot(transmit));
          
          float Frpar = ((n2*ci-n1*ct)/(n2*ci+n1*ct))*((n2*ci-n1*ct)/(n2*ci+n1*ct));
          float Frper = ((n1*ct-n2*ci)/(n1*ct+n2*ci))*((n1*ct-n2*ci)/(n1*ct+n2*ci));
          
          float Reff = 1-0.5*(Frpar+Frper);

          PVector temppos = new PVector(x-img.width/2,y-img.height/2,0);
          PVector temppospixel = new PVector(pixelx-img.width/2,img.height-y+1);

          //cast color
          reflectColor = cubemap_cast(temppos, temppospixel, Re, d, diffuseDark);                      
          refractColor = blur(((pixelx)+floor((100-abs(z))*transmit.x/transmit.z)+floor(bg.x))+512, ((pixely)-floor((100-abs(z))*transmit.y/transmit.z)-floor(bg.y))+512, glossy, diffuseDark);
          objectColor = color(floor((red(reflectColor))*(Reff)+(red(refractColor))*(1-Reff)),floor((green(reflectColor))*(Reff)+(green(refractColor))*(1-Reff)),floor((blue(reflectColor))*(Reff)+(blue(refractColor))*(1-Reff)));

          //adding diffuse
          diffuse = max(znorm.dot((PVector.mult(light,-1))),0);
      
          //adding spec
          spec = specHighlight(znorm,light,viewer,specStrength);
          
          nimg.pixels[i] = color(spec*red(specLightColor)+(1-spec)*(ambient+diffuse)*red(lightColor)*red(objectColor)/255, spec*green(specLightColor)+(1-spec)*(ambient+diffuse)*green(lightColor)*green(objectColor)/255, spec*blue(specLightColor)+(1-spec)*(ambient+diffuse)*blue(lightColor)*blue(objectColor)/255);
  
          //masking out backgound
          if (red(img.get(pixelx,pixely)) == 0) {
            nimg.pixels[i] = diffuseDark.get(pixelx-250+diffuseDark.width/2,pixely-250+diffuseDark.height/2);
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
        
        for (int i = 0; i < totalImageSize; i = i + 1) {
          x = (i+1)%nimg.width;
          y = nimg.height-floor(((i+1)/nimg.width))-1;
          pixelx = x;
          pixely = floor(((i+1)/nimg.width));

          if ((x-floor(nimg.width/2))*(x-floor(nimg.width/2)) + (y-floor(nimg.height/2))*(y-floor(nimg.height/2)) < 10000) {
            z = (floor(sqrt(10000-(x-floor(0.5*nimg.width))*(x-floor(0.5*nimg.width))-(y-floor(0.5*nimg.height))*(y-floor(0.5*nimg.height)))));
            
            //normal vector
            znorm = new PVector(x-250,y-250,z);
            znorm.normalize();
      
            //calculate reflect vector
            Re = (PVector.mult(znorm,2*znorm.dot(viewer))).sub(viewer);
            Re.normalize();
            
            //calculate transmit vector            
            ci = max(znorm.dot(light),0);
            c2 = sqrt(1 - n*n*(1-ci*ci));
            
            transmit = PVector.mult(light, n).add(PVector.mult(znorm,((n)*ci-c2)));
            transmit.normalize();
            
            //calculate Fresnel
            //angle of transmittance
            float ct = abs(znorm.dot(transmit));
            
            //check for total reflection
            float Frpar = ((n2*ci-n1*ct)/(n2*ci+n1*ct))*((n2*ci-n1*ct)/(n2*ci+n1*ct));
            float Frper = ((n1*ct-n2*ci)/(n1*ct+n2*ci))*((n1*ct-n2*ci)/(n1*ct+n2*ci));
            
            float Reff = 1-0.5*(Frpar+Frper);
            //transmission only
            Reff = 1;

            PVector temppos = new PVector(x-250,y-250,z);
            PVector temppospixel = new PVector(pixelx-250,pixely-250);
            //objectColor = cubemap_cast(temppos, temppospixel, Re, d);
            
            //cast color
            reflectColor = cubemap_cast(temppos, temppospixel, Re, d, diffuseDark);
            
            temppos = new PVector(x-250,y-250,-z);            
            
            //transmit vector cast differently bc i'm too lazy to figure out non-centered cube mapping 
            refractColor = blur(((pixelx)+floor((100-abs(z))*transmit.x/transmit.z)+floor(bg.x))+512, ((pixely)-floor((100-abs(z))*transmit.y/transmit.z)-floor(bg.y))+512, glossy, diffuseDark);

            objectColor = color(floor((red(reflectColor))*(1-Reff)+(red(refractColor))*(Reff)),floor((green(reflectColor))*(1-Reff)+(green(refractColor))*(Reff)),floor((blue(reflectColor))*(1-Reff)+(blue(refractColor))*(Reff)));

            //adding diffuse
            diffuse = max(znorm.dot((PVector.mult(light,-1))),0);
      
            //adding spec
            spec = specHighlight(znorm,light,viewer,specStrength);

            nimg.pixels[i] = color(spec*red(specLightColor)+(1-spec)*(ambient+diffuse)*red(lightColor)*red(objectColor)/255, spec*green(specLightColor)+(1-spec)*(ambient+diffuse)*green(lightColor)*green(objectColor)/255, spec*blue(specLightColor)+(1-spec)*(ambient+diffuse)*blue(lightColor)*blue(objectColor)/255);

          }
          else {
            nimg.pixels[i] = diffuseDark.get(pixelx+diffuseDark.width/4,pixely-250+diffuseDark.height/2);
          }
        }
        done = true;
        break;
        
      case "normal sphere":
        //normal map
        nimgPath = "./images/normal_sphere.png";
        nimg = createImage(500,500,RGB);
        
        totalImageSize = nimg.width * nimg.height;
        // Make sure the pixels are stored properly in the array

        nimg.loadPixels();
  
        //sphere code
        for (int i = 0; i < totalImageSize; i = i + 1) {
          x = (i+1)%nimg.width;
          y = nimg.height-floor(((i+1)/nimg.width))-1;
          pixelx = x;
          pixely = floor(((i+1)/nimg.width));

          if ((x-floor(nimg.width/2))*(x-floor(nimg.width/2)) + (y-floor(nimg.height/2))*(y-floor(nimg.height/2)) < 10000) {
            z = (floor(sqrt(10000-(x-floor(0.5*nimg.width))*(x-floor(0.5*nimg.width))-(y-floor(0.5*nimg.height))*(y-floor(0.5*nimg.height)))));
            
            //normal vector
            znorm = new PVector((x-250),(y-250),z);
            znorm.normalize();
            
            //calculate reflect vector
            Re = (PVector.mult(znorm,2*znorm.dot(viewer))).sub(viewer);
            Re.normalize();         

            //calculate transmit vector            
            ci = max(znorm.dot(light),0);
            c2 = sqrt(1 - n*n*(1-ci*ci));
            
            transmit = PVector.mult(light, n).add(PVector.mult(znorm,((n)*ci-c2)));
            transmit.normalize();
            
            //calculate Fresnel            
            //angle of transmittance
            float ct = abs(znorm.dot(transmit));
            
            float Frpar = ((n2*ci-n1*ct)/(n2*ci+n1*ct))*((n2*ci-n1*ct)/(n2*ci+n1*ct));
            float Frper = ((n1*ct-n2*ci)/(n1*ct+n2*ci))*((n1*ct-n2*ci)/(n1*ct+n2*ci));
            
            float Reff = 1-0.5*(Frpar+Frper);
            //transmission only
            Reff = 1;

            PVector temppos = new PVector(x-250,y-250,0);
            PVector temppospixel = new PVector(pixelx-250,pixely-250);
            //objectColor = cubemap_cast(temppos, temppospixel, Re, d);
            reflectColor = cubemap_cast(temppos, temppospixel, Re, d, diffuseDark);                        
            refractColor = blur(((pixelx)+floor((100-abs(z))*transmit.x/transmit.z)+floor(bg.x))+512, ((pixely)-floor((100-abs(z))*transmit.y/transmit.z)-floor(bg.y))+512, glossy, diffuseDark);
            objectColor = color(floor((red(reflectColor))*(1-Reff)+(red(refractColor))*(Reff)),floor((green(reflectColor))*(1-Reff)+(green(refractColor))*(Reff)),floor((blue(reflectColor))*(1-Reff)+(blue(refractColor))*(Reff)));

            //adding diffuse
            diffuse = max(znorm.dot((PVector.mult(light,-1))),0);
      
            //adding spec
            spec = specHighlight(znorm,light,viewer,specStrength);
            
            nimg.pixels[i] = color((ambient+spec+diffuse)*red(lightColor)*red(objectColor)/255, (ambient+spec+diffuse)*green(lightColor)*green(objectColor)/255, (ambient+spec+diffuse)*blue(lightColor)*blue(objectColor)/255);
          }
          else {
            nimg.pixels[i] = diffuseDark.get(pixelx-250+diffuseDark.width/2,pixely-250+diffuseDark.height/2);
          }
        }
      
      done = true;
      //println("u got in");
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
