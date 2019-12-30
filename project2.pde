//
// VIZA 656: Digital Image
// Name: Anne Fu
// Project 2
//

//antialias function for generated spheres
PImage antialias(PImage vimg) {
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

//blur function for glossy reflection
color blur(int x, int y, int glossy, PImage vimg) {
  color avg = vimg.get(x, y);
  float avgr = red(avg);
  float avgg = green(avg);
  float avgb = blue(avg);

  //get color values of the pixels surrounding
  for(int j = -1*glossy; j < glossy; j = j+1) {
    for(int k = -1*glossy; k < glossy; k = k+1) {
      color a = vimg.get(x+floor(random(0,glossy)), y+floor(random(0,glossy)));
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

//function for calculating specular highlight
float specHighlight (PVector znorm, PVector light, PVector viewer, float specStrength) {
  reflect = (PVector.mult(znorm,2.0*znorm.dot(PVector.mult(light,-1)))).sub(light);
  reflect.normalize();
  spec = pow(max(viewer.dot(reflect), 0.0), 64)*specStrength;
  return spec;
}

PImage img;
PImage nimg;
PImage diffuseDark;
PImage diffuseLight;
int x, y;
int pixelx, pixely;
int glossy;
float z;
color lightColor;
color objectColor;
PVector norm;
PVector znorm;
PVector light;
PVector viewer;
PVector reflect;
PVector gradient;
PVector bg;
PVector Re;
float n;
int h;
int totalImageSize;
float d;
float ambient, diffuse, spec;
float specStrength;
float zratio;
float gradientx, gradienty, gradientxz,gradientyz,gradientxo,gradientyo;
float depthRatio;
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
  ambient = 0.5; 
  bg = new PVector(0,0,0);
  //environment map
  diffuseDark = loadImage("./images/google-alphabet.jpg");
  lightColor = color(255, 255, 255);
  //distance
  d = 150;
  //glossiness factor; increase to blur more
  glossy = 0;
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
    input = input + key; 
    loop();
  } 
  else {
    userinput = input;
    if (input.equals("normal") || input.equals("depth")) {
      mapTemp = input;
      mapType = true;
      loop();
    }
    else if (input.equals("sphere")) {
      userinput = mapTemp + " " + input;
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
        //load map
        imgPath = "./images/" + imgName;
        imgNameParse = split(imgName, '.');
        //load path for rendered image
        nimgPath = "./images/" + imgNameParse[0] + "_shaded." + imgNameParse[1];
        
        img = loadImage(imgPath);
        nimg = createImage(img.width,img.height,RGB);
    
        totalImageSize = img.width * img.height;

        img.loadPixels();
        nimg.loadPixels();
     
        for (int i = 0; i < totalImageSize; i = i + 1) {
          //calculating equivalent x, y pos for Processing image structure
          pixelx = (i+1)%img.width;
          pixely = floor(((i+1)/img.width));
                    
          z = red(img.get(pixelx,pixely));
          depthRatio = z/255;
          
          //set h for how far in each direction for gradient calculation, must be greater than 1
          h = 2;
          
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
          
          //calculate reflection vector
          Re = (PVector.mult(znorm,2*znorm.dot(viewer))).sub(viewer);
          Re.normalize();
          
          //cast color
          //objectColor = diffuseDark.get(pixelx+floor((d-z)*Re.x)+floor(bg.x)+1000,pixely-floor((d-z)*Re.y)-floor(bg.y));
          
          objectColor = blur((pixelx-img.width/2+floor((d-z/2)*Re.x/Re.z)+floor(bg.x))%(diffuseDark.width/2)+(diffuseDark.width/2), (pixely-img.height/2-floor((d-z/2)*Re.y/Re.z)-floor(bg.y))%(diffuseDark.height/2)+(diffuseDark.height/2), glossy, diffuseDark);
          
          //mask out black areas of map
          if (z == 0) {
            objectColor = color(0,0,0);
          }
      
          //adding diffuse
          diffuse = max(znorm.dot((PVector.mult(light,-1))),0);
      
          //adding spec
          spec = specHighlight(znorm,light,viewer,specStrength);
          
          nimg.pixels[i] = color((ambient+spec+diffuse)*red(lightColor)*red(objectColor)/255, (ambient+spec+diffuse)*green(lightColor)*green(objectColor)/255, (ambient+spec+diffuse)*blue(lightColor)*blue(objectColor)/255);
        }
        done = true;
        break;
      
      case "normal image":
        //load map
        imgPath = "./images/" + imgName;
        imgNameParse = split(imgName, '.');
        //load path for rendered image
        nimgPath = "./images/" + imgNameParse[0] + "_shaded." + imgNameParse[1];
        
        img = loadImage(imgPath);
        nimg = createImage(img.width,img.height,RGB);
    
        totalImageSize = img.width * img.height;

        img.loadPixels();
        nimg.loadPixels();
        
        //initialize znorm
        znorm = new PVector(0,0,0);
        
        // normal map
        for (int i = 0; i < totalImageSize; i = i + 1) {
          pixelx = (i+1)%img.width;
          pixely = floor(((i+1)/img.width));
          objectColor = diffuseDark.get(pixelx,pixely);

          znorm.x = (red(img.get(pixelx,pixely))/256)*2-1;
          znorm.y = (green(img.get(pixelx,pixely))/256)*2-1;
          znorm.z = (blue(img.get(pixelx,pixely))/256)*2-1;         
         
          znorm.normalize();
         
          //calculate reflection vector
          Re = (PVector.mult(znorm,2*znorm.dot(viewer))).sub(viewer);
          Re.normalize();
          
          //cast color          
          objectColor = blur((pixelx-img.width/2+floor((d)*Re.x/Re.z)+floor(bg.x))%(diffuseDark.width/2)+(diffuseDark.width/2), (pixely-img.height/2-floor((d)*Re.y/Re.z)-floor(bg.y))%(diffuseDark.height/2)+(diffuseDark.height/2), glossy, diffuseDark);

          //mask out black areas of map
          if (red(img.get(pixelx,pixely)) == 0) {
            objectColor = color(0,0,0);
          }
      
          //adding diffuse
          diffuse = max(znorm.dot((PVector.mult(light,-1))),0);
      
          //adding spec
          spec = specHighlight(znorm,light,viewer,specStrength);
          
          nimg.pixels[i] = color((ambient+spec+diffuse)*red(lightColor)*red(objectColor)/255, (ambient+spec+diffuse)*green(lightColor)*green(objectColor)/255, (ambient+spec+diffuse)*blue(lightColor)*blue(objectColor)/255);
        }
        done = true;
        break;
        
      case "depth sphere":      
        nimgPath = "./images/depth_sphere.png";
        nimg = createImage(500,500,RGB);
        
        totalImageSize = nimg.width * nimg.height;

        nimg.loadPixels();

        for (int i = 0; i < totalImageSize; i = i + 1) {
          x = (i+1)%nimg.width;
          y = nimg.height-floor(((i+1)/nimg.width))-1;
          pixelx = x;
          pixely = floor(((i+1)/nimg.width));

          //create centered sphere w/radius 100
          if ((x-floor(nimg.width/2))*(x-floor(nimg.width/2)) + (y-floor(nimg.height/2))*(y-floor(nimg.height/2)) < 10000) {
            //calculate z
            z = (floor(sqrt(10000-(x-floor(0.5*nimg.width))*(x-floor(0.5*nimg.width))-(y-floor(0.5*nimg.height))*(y-floor(0.5*nimg.height)))));
                  
            znorm = new PVector(x-250,y-250,z);
            znorm.normalize();      
            
            //calculate reflection vector
            Re = (PVector.mult(znorm,2*znorm.dot(viewer))).sub(viewer);
            Re.normalize();
            

            //cast color
            objectColor = blur((pixelx-250+floor((d-z)*Re.x/Re.z)+floor(bg.x))%(diffuseDark.width/2)+diffuseDark.width/2, (pixely-250-floor((d-z)*Re.y/Re.z)-floor(bg.y))%(diffuseDark.height/2)+diffuseDark.height/2, glossy, diffuseDark);
                        
            //for cube map
            /*if (abs(Re.z)>abs(Re.x) && abs(Re.z)>abs(Re.y)) {
              //z is max
              if (Re.z<0) {
                objectColor = blur(((pixelx-250)-floor((d-z)*Re.x/Re.z)-floor(bg.x))%400+1200, ((pixely-250)+floor((d-z)*Re.y/Re.z)+floor(bg.y))%400+1200, glossy, diffuseDark);
                if(abs((pixely-250)+floor((d-z)*Re.y/Re.z)+floor(bg.y))>400) {
                  objectColor = blur(((pixelx-250)-floor((d-z)*Re.x/Re.z)-floor(bg.x))/400+1200, ((pixely-250)+floor((d-z)*Re.y/Re.z)+floor(bg.y))/400+1200, glossy, diffuseDark);
                  if (Re.y<0) {
                    objectColor = blur(((pixelx-250)-floor((d+(y-250))*Re.x/Re.y)-floor(bg.x))%400+1200, ((pixely-250)-floor((d+(y-250))*Re.z/Re.y)-floor(bg.y))%400+2000, glossy, diffuseDark);
                  }
                  else {
                    objectColor = blur(((pixelx-250)+floor((d-(y-250))*Re.x/Re.y)+floor(bg.x))%400+1200, ((pixely-250)-floor((d-(y-250))*Re.z/Re.y)-floor(bg.y))%400+400, glossy, diffuseDark);
                  }
                }
                if(abs((pixelx-250)-floor((d-z)*Re.x/Re.z)+floor(bg.x))>400) {
                  objectColor = blur(((pixelx-250)-floor((d-z)*Re.x/Re.z)-floor(bg.x))/400+1200, ((pixely-250)+floor((d-z)*Re.y/Re.z)+floor(bg.y))/400+1200, glossy, diffuseDark);
                  if (Re.x<0) {
                    objectColor = blur(((pixelx-250)+floor((d+(x-250))*Re.z/Re.x)+floor(bg.x))%400+400, ((pixely-250)+floor((d+(x-250))*Re.y/Re.x)+floor(bg.y))%400+1200, glossy, diffuseDark);
                  }
                  else {
                    objectColor = blur(((pixelx-250)+floor((d-(x-250))*Re.z/Re.x)+floor(bg.x))%400+2000, ((pixely-250)-floor((d-(x-250))*Re.y/Re.x)-floor(bg.y))%400+1200, glossy, diffuseDark);
                  }
                }
              }
              else {
                objectColor = blur(2800-(((pixelx-250)+floor((d-z)*Re.x/Re.z)+floor(bg.x))%400), ((pixely-250)-floor((d-z)*Re.y/Re.z)-floor(bg.y))%400+1200, glossy, diffuseDark);
              }
            }
            else if (abs(Re.x)>abs(Re.y) && abs(Re.x)>abs(Re.z)) {
              //x is max
              if (Re.x<0) {
                objectColor = blur(((pixelx-250)+floor((d+(x-250))*Re.z/Re.x)+floor(bg.x))%400+400, ((pixely-250)+floor((d+(x-250))*Re.y/Re.x)+floor(bg.y))%400+1200, glossy, diffuseDark);
              }
              else {
                objectColor = blur(((pixelx-250)+floor((d-(x-250))*Re.z/Re.x)+floor(bg.x))%400+2000, ((pixely-250)-floor((d-(x-250))*Re.y/Re.x)-floor(bg.y))%400+1200, glossy, diffuseDark);
              }
            }
            else {
              //y is max
              if (Re.y<0) {
                objectColor = blur(((pixelx-250)-floor((d+(y-250))*Re.x/Re.y)-floor(bg.x))%400+1200, ((pixely-250)-floor((d+(y-250))*Re.z/Re.y)-floor(bg.y))%400+2000, glossy, diffuseDark);
              }
              else {
                objectColor = blur(((pixelx-250)+floor((d-(y-250))*Re.x/Re.y)+floor(bg.x))%400+1200, ((pixely-250)-floor((d-(y-250))*Re.z/Re.y)-floor(bg.y))%400+400, glossy, diffuseDark);
              }
            }*/

            /*if(Re.z<0) {
              objectColor = blur((pixelx+floor((d+z)*Re.x/Re.z)+floor(bg.x))%(diffuseDark.width/2)+1200, (pixely+floor((d+z)*Re.y/Re.z)-floor(bg.y))%(diffuseDark.height/2)+700, glossy, diffuseDark);
            }*/

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
        break;
        
      case "normal sphere":
        nimgPath = "./images/normal_sphere.png";
        nimg = createImage(500,500,RGB);
        
        totalImageSize = nimg.width * nimg.height;

        nimg.loadPixels();
  
        //create centered sphere w/radius 100
        for (int i = 0; i < totalImageSize; i = i + 1) {
          x = (i+1)%nimg.width;
          y = nimg.height-floor(((i+1)/nimg.width))-1;
          pixelx = x;
          pixely = floor(((i+1)/nimg.width));

          if ((x-floor(nimg.width/2))*(x-floor(nimg.width/2)) + (y-floor(nimg.height/2))*(y-floor(nimg.height/2)) < 10000) {
            z = (floor(sqrt(10000-(x-floor(0.5*nimg.width))*(x-floor(0.5*nimg.width))-(y-floor(0.5*nimg.height))*(y-floor(0.5*nimg.height)))));
            znorm = new PVector((x-250)/z,(y-250)/z,1);
            znorm.normalize();
            
            //calculate reflection vector
            Re = (PVector.mult(znorm,2*znorm.dot(viewer))).sub(viewer);
            Re.normalize();
            
            //cast color            
            objectColor = blur((pixelx-250+floor((d)*Re.x/Re.z)+floor(bg.x))%(diffuseDark.width/2)+diffuseDark.width/2, (pixely-250-floor((d)*Re.y/Re.z)-floor(bg.y))%(diffuseDark.height/2)+diffuseDark.height/2, glossy, diffuseDark);

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
      break;
    
    default:
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

void draw() {
  background(255);
  int indent = 25;
  
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
    image(nimg,0,0);
    save(nimgPath);
  }
}
