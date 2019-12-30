1) Run .exe
2) Initial prompt will ask whether "normal" or "depth" information will be provided. Select "normal" if the input image will be a normal map, "depth" if the input image will be a depth map.

To see the shaded generated sphere, select "depth" for this prompt.

Hit return.

3) Second prompt will ask whether a sphere will be generated or an input image will be provided. If an input image is being provided, type "image" followed by a space, then the full file name for the map image (e.g. image normal map.png).

The reference images for dark diffuse and light diffuse should have the same name as the map image followed by 1 or 2, respectively (e.g. normalmap1.png, normalmap2.png). All input images should be stored in the ./images folder.

To see the shaded generated sphere, select "sphere" for this prompt (*).

Hit return.

(*) if "normal" + "sphere" is selected, the output image will be the normal map for the generated sphere.

4) The window will be refreshed and display the shaded image once calculations are complete. The shaded image will also be saved to the ./images folder as the map image name + "_shaded" (e.g. normalmap_shaded.png).

5) The background position can be updated by clicking around the screen. Glossiness can be set using the up & down arrow keys (*). Only the latest image will be saved.

(*) the higher the glossiness, the longer the image will take to update.

6) Right-click to return to the first prompt.

Included depth map:
depthmap.png

Included normal map:
wiki.png
