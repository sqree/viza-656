1) Run .exe.

2) Scene should automatically load with 1 planes, 3 spheres, & 1 infinite sphere. Default is uv-wrapped texture. The shaded image will be saved to the ./images folder as ‘rasterized.png’.

3) The light position can be updated by clicking around the screen. UP, DOWN, RIGHT, & LEFT arrows can be used to pivot camera up, down, right & left, respectively, ALT & CTRL to zoom in & out.

Lighting modes:
'L' = direction
'P' = point
'S' = spot*
(*) Only working with hard shadows; will automatically default to 'P' if shadow mode changed to soft area or d/R.

Shadow modes:
'A' = soft area*
'D' = soft d/R
'H' = hard cosTheta
(*) May take longer to load. Only working with point light mode; will automatically default to 'H' if light mode changed to direction or spot.

Texturing modes:
’N’ = normal
‘W’ = projection (parallel default)
‘E’ = perspective (*)
‘P’ = procedural texture
‘Z’ = default (uv wrapped)
(*) Must be pressed in projection mode (‘W’) or scene will not update correctly.

Only the latest image will be saved. 
