1) Run .exe.
Note: Majority of shadow implementation is in the castshadow() function in the light class.

2) Scene should automatically load with 2 planes & 3 spheres. Default is point light with hard shadow. The shaded image will also be saved to the ./images folder as ‘rasterized.png’.

3) The light position can be updated by clicking around the screen. UP, DOWN, RIGHT, & LEFT arrows can be used to pivot camera up, down, right & left, respectively, ALT & CTRL to zoom in & out.

Lighting modes:
'L' = direction
'P' = point
'S' = spot*
(*) Only working with hard shadows; will automatically default to 'P' if shadow mode changed to soft area or d/R.

3) Shadow modes:
'A' = soft area*
'D' = soft d/R
'H' = hard cosTheta
(*) May take longer to load. Only working with point light mode; will automatically default to 'H' if light mode changed to direction or spot.

Only the latest image will be saved. 
