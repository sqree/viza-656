Triangle meshes

1) Run .exe.
Note: Majority of added .obj support is in the geometry class (triangle) and first half of the draw function (prior to shadow casting).

2) Scene should automatically load with 2 planes, 2 cubes (1 textured), a tetrahedron, and a dodecahedron. Default is point light with hard shadow. The shaded image will also be saved to the ./images folder as ‘rasterized.png’.

Disclaimer: All the below functionality still exists (theoretically, I haven’t tested them out), but will increase render time significantly.

----------

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
