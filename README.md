# Installation 
To use this bundle, first download the whole repository and save it in a desired directory. Next, open ChimeraX and run the following command:
```
devel build Path/To/Source/Code
```
and
```
devel install Path/To/Source/Code
```
and
```
devel clean Path/To/Source/Code
```

# Usage
## Illustrate point clouds 
Both methods covered in thid bundle, rely on converting density maps to points clouds. For this pupose we used two methods, Topology Representing Networks (TRN) and Centroidal Voronoi Tesselation (CVT). For illustrations purposes you can use ```ot_alignment show_points``` command to show the point_clouds on chimera. The basic usage is:
```
ot_alignment show_points [model_id]
```
Also, these are some optional keywords:
* `n`: Number of sampled points. the default value is `500`.
* `size`: Radius of each sphere representing the points. The default value is `2`.
* `random_seed`: Rnadom seed used in methods.
* `thresh`: Thresholding value used on the maps. Highly recommend to use the default level value that chimera sets for illustration of the map. The default valu is `0.0`.
* `color`: Color used for the spheres that represent points. The default value is `#ff5a5a`.
* `sampling_method`: The moethod used for sampling point clouds. Only valid options are `cvt` and `trn` with default to `trn`.

Here are some examples of usage of this command:
``` ot_alignment show_points #1 sampling_method cvt n 50 random_seed 1 color #55aaff thresh 0.1```
``` ot_alignment show_points #2```

## Perform AlignOT
To use the AlignOT method for alignment of two density maps, first open two maps and put them in the desired initial position. Next, use the `ot_alignment AlignOT` command. The basic usage is:
```
ot_alignment AlignOT [model_id1] [model_id2]
```
In result this command moves **model_id2** to align with **model_id1**(fixed). Also, these are some optional keywords:
* `n`: Number of sampled points. the default value is `500`.
* `thresh`: Thresholding value used on the maps. Highly recommend to use the default level value that chimera sets for illustration of the map. The default valu is `0.0`.
* `max_iter`: Maximum number of iterations, before stoping the algorithm. the default value is `500`.
* `random_seed`: Random seed used in methods.
* `lr`: The learning rate used in the SGD-like algorithm. The default value is `0.0001`.
* `reg`: The regularization parameter used for computing entropy-regularized wasserstein distance. Teh default value is `100`.

Here are some examples of usage of this command:
```ot_alignment AlignOT #1 #2 random_seed 1 max_iter 100 thresh 0.1```
```ot_alignment AlignOT #3 #4```

## Perform EMPOT
To use the EMPOT method for alignment of two density maps, first open two maps and put them in the desired initial position. Next, use the `ot_alignment EMPOT` command. The basic usage is:
```
ot_alignment EMPOT [model_id1] [model_id2]
```
In result this command moves **model_id2** to align with **model_id1**(fixed). Also, these are some optional keywords:
* `n`: Number of sampled points. the default value is `500`.
* `thresh`: Thresholding value used on the maps. Highly recommend to use the default level value that chimera sets for illustration of the map. The default valu is `0.0`.
* `num`: Number of alignments done, before choosing the best one. The default value is `2`.
* `random_seed`: Random seed used in methods.

Here are some examples of usage of this command:
```ot_alignment EMPOT #1 #2 random_seed 1 num 5 thresh 0.1```
```ot_alignment EMPOT #3 #4```