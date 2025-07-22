Visionaray Ambient Occlusion Example
------------------------------------

### Command line

```
Usage:
   hpcnrw [OPTIONS] filename

Positional options:
   filename               Input file in wavefront obj format

Options:
   -bvh=<ARG>             BVH build strategy:
      =default            - Binned SAH
      =split              - Binned SAH with spatial splits
   -camera=<ARG>          Text file with camera parameters
   -width=<ARG>           Image width
   -height=<ARG>          Image height
   -threads=<ARG>         Number of threads
   -spp=<ARG>             Number of frames to be accumulated
```

