# Description
ot_alignment is a plugin for ChimeraX, which allows users to align density maps with methods based on the Optimal Transport (OT) theory. This bundle aggregates the methods used in [AlignOT]{https://arxiv.org/abs/2210.09361} and [EMPOT]{https://arxiv.org/abs/2311.00850} paper. 

This is the repository for distributing our software. For installing and running the software, please refer to the User Manual pdf file or see below.


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
## Generating point clouds 
Alignment methods covered in this bundle rely on converting density maps to point clouds. For this purpose, two methods can be used: (1) Topology Representing Networks (TRN) and (2) Centroidal Voronoi Tesselation (CVT). To visualize point clouds, use ```ot_alignment show_points``` command to show the point_clouds on ChimeraX. The basic usage is:
```
ot_alignment show_points [model_id]
```
with the optional keywords:
* `n`: Number of sampled points. the default value is `500`.
* `size`: Radius of each sphere representing the points. The default value is `2`.
* `random_seed`: Random seed used in methods.
* `thresh`: Thresholding value used on the maps. Highly recommend to use the default level value that ChimeraX sets for the visualization of the map. The default value is `0.0`.
* `color`: Color used for the spheres that represent points. The default value is `#ff5a5a`.
* `sampling_method`: The method used for sampling point clouds. The only valid options are `cvt` and `trn` with default to `trn`.

Here are some examples of usage of this command:
```
ot_alignment show_points #1 sampling_method cvt n 50 random_seed 1 color #55aaff thresh 0.1
```
```
ot_alignment show_points #2
```

## Performing AlignOT
To use the AlignOT method for the alignment of two density maps, first open two maps and put them in the desired initial position. Next, use the `ot_alignment AlignOT` command. The basic usage is:
```
ot_alignment AlignOT [model_id1] [model_id2]
```
As a result this command moves **model_id2** to align with **model_id1**(fixed). Also, these are some optional keywords:
* `n`: Number of sampled points. the default value is `500`.
* `thresh`: Thresholding value used on the maps. Highly recommend to use the default level value that ChimeraX sets for the visualization of the map. The default value is `0.0`.
* `max_iter`: Maximum number of iterations, before stopping the algorithm. the default value is `500`.
* `random_seed`: Random seed used in methods.
* `lr`: The learning rate used in the SGD-like algorithm. The default value is `0.0001`.
* `reg`: The regularization parameter used for computing entropy-regularized Wasserstein distance. The default value is `100`.
* `sampling_method`: The method used for sampling point clouds. The only valid options are `cvt` and `trn` with default to `trn`.
* `local_refinement`: Whether we use a final step of `fitmap` to refine results or not. The default value is `True`.

Here are some examples of usage of this command:
```
ot_alignment AlignOT #1 #2 random_seed 1 max_iter 100 thresh 0.1
```
```
ot_alignment AlignOT #3 #4
```

## Performing EMPOT
To use the EMPOT method for the alignment of two density maps, first open two maps and put them in the desired initial position. Next, use the `ot_alignment EMPOT` command. The basic usage is:
```
ot_alignment EMPOT [model_id1] [model_id2]
```
As a result this command moves **model_id2** to align with **model_id1**(fixed). Also, these are some optional keywords:
* `n`: Number of sampled points. the default value is `500`.
* `thresh`: Thresholding value used on the maps. Highly recommend to use the default level value that ChimeraX sets for the visualization of the map. The default value is `0.0`.
* `num`: Number of alignments done, before choosing the best one. The default value is `2`.
* `random_seed`: Random seed used in methods.
* `sampling_method`: The method used for sampling point clouds. The only valid options are `cvt` and `trn` with default to `trn`.
* `local_refinement`: Whether we use a final step of `fitmap` to refine results or not. The default value is `True`.

Here are some examples of usage of this command:
```
ot_alignment EMPOT #1 #2 random_seed 1 num 5 thresh 0.1
```
```
ot_alignment EMPOT #3 #4
```
