// get river width
// dhemerson.costa@ipam.org.br

// define the years to be computed 
var years = [
  2023
  ];

// read water surface
var collection = ee.ImageCollection("projects/mapbiomas-workspace/COLECAO9/agua")
  .filter("version == '4'");
// for each year, convert into image
var recipe = ee.Image([]);  // empty recipe 
years.forEach(function(year_i) {
  // convert into image and store
   recipe = recipe.addBands(
     collection.filter(ee.Filter.eq('year', year_i))
      .mosaic()
      .remap([33], [1])
      .rename('classification_' + year_i));
});

// plot
//Map.addLayer(test, {palette:['blue'], min:1 , max:1}, 'water');

/* high-level image functions*/
var AssignDefault = function(x, dv) {
  return(typeof x !== 'undefined' ? x : dv);
};

var rwGenSR = function(MAXDISTANCE, FILL_SIZE, MAXDISTANCE_BRANCH_REMOVAL, AOI) {
   
  var grwl = ee.FeatureCollection("users/eeProject/grwl");

  // generate function based on user choice
  var tempFUN = function(imageID) {
    
    // get mapbiomas classification
    var image = recipe.select(imageID);
    
    // clip by AOI
    AOI = ee.Algorithms.If(AOI, AOI, image.geometry());
    image = image.clip(AOI);
    
    /////////////////////////////////////////////
    // derive water mask 
    var imgOut = image.rename('waterMask');


    /////////////////////////////////////////////
    // calculate river mask
    /* functions to extract river mask */
    var GetCenterline = function(clDataset, bound) {
      // filter the GRWL centerline based on area of interest
      var cl = clDataset.filterBounds(bound); 
      return(cl);
    };
    
    var ExtractChannel = function(image, centerline, maxDistance) {
      // extract the channel water bodies from the water mask, based on connectivity to the reference centerline.
      var connectedToCl = image.not().cumulativeCost({
        source: ee.Image().toByte().paint(centerline, 1).and(image), // only use the centerline that overlaps with the water mask
        maxDistance: maxDistance,
        geodeticDistance: false
      }).eq(0);
      var channel = image.updateMask(connectedToCl).unmask(0).updateMask(image.gte(0)).rename(['channelMask']);
      return channel;
    };
    
    var RemoveIsland = function(channel, FILL_SIZE) {
      /* fill in island as water if the size (number of pixels) of the island is smaller than FILL_SIZE */
      var fill = channel.not().selfMask().connectedPixelCount(FILL_SIZE).lt(FILL_SIZE);
      var river = channel.where(fill, ee.Image(1)).rename(['riverMask']);
      return(river);
    };
    
    var ExtractRiver = function(imgIn, clData, maxDist, minIslandRemoval) {
      var waterMask = imgIn.select('waterMask');
      var bound = waterMask.geometry();
      var cl = GetCenterline(clData, bound);
      var channelMask = ExtractChannel(waterMask, cl, maxDist);
      var riverMask = RemoveIsland(channelMask, minIslandRemoval);
      return(imgIn.addBands(channelMask).addBands(riverMask));
    };
    
    imgOut = ExtractRiver(imgOut, grwl, MAXDISTANCE, FILL_SIZE);
    //Map.addLayer(imgOut.select('channelMask'), {palette:['green'], min:1 , max:1}, 'channel');
    //Map.addLayer(imgOut.select('riverMask'), {palette:['red'], min:1 , max:1}, 'river');
    
    //////////////////////////////////////////////////////////////////////////////////////////
    var CalcDistanceMap = function(img, neighborhoodSize, scale) {
    // assign each river pixel with the distance (in meter) between itself and the closest non-river pixel
      var imgD2 = img.focal_max(1.5, 'circle', 'pixels', 2);
      var imgD1 = img.focal_max(1.5, 'circle', 'pixels', 1);
      var outline = imgD2.subtract(imgD1);
      var dpixel = outline.fastDistanceTransform(neighborhoodSize).sqrt();
      var dmeters = dpixel.multiply(scale); // for a given scale
      var DM = dmeters.mask(dpixel.lte(neighborhoodSize).and(imgD2));
      return(DM.rename(['distanceMap']));
      };
    
    var CalcGradientMap = function(image, gradMethod, scale) {
        // calculate gradient of the given image
        var dx, dy, g, k_dx, k_dy;
        
        // calculate the gradient
        if (gradMethod == 1) {
          var grad = image.gradient();
          dx = grad.select(['x']);
          dy = grad.select(['y']);
          g = dx.multiply(dx).add(dy.multiply(dy)).sqrt();
        }
        if (gradMethod == 2) {
          k_dx = ee.Kernel.fixed(3, 3,
                                [[ 1/8,  0,  -1/8],
                                  [ 2/8,  0,  -2/8],
                                  [ 1/8,  0,  -1/8]]);
          k_dy = ee.Kernel.fixed(3, 3,
                                [[ -1/8, -2/8,  -1/8],
                                  [ 0,    0,    0],
                                  [ 1/8, 2/8,   1/8]]);
          dx = image.convolve(k_dx);
          dy = image.convolve(k_dy);
          g = dx.multiply(dx).add(dy.multiply(dy)).divide(scale.multiply(scale)).sqrt();
        }
        if (gradMethod == 3) {
          k_dx = ee.Kernel.fixed(3, 1,
                                [[-0.5,  0,  0.5]]);
          k_dy = ee.Kernel.fixed(1, 3,
                                [[0.5],
                                  [0],
                                  [-0.5]]);
          dx = image.convolve(k_dx);
          dy = image.convolve(k_dy);
          g = dx.multiply(dx).add(dy.multiply(dy)).divide(scale.multiply(scale));
        }
        
        return(g);
      };
      
    var SplitKernel = function(kernel, value) {
      // recalculate the kernel according to the given foreground value
      var result = [];
      for(var r = 0; r < kernel.length; r++) {
        var row = [];
        for(var c = 0; c < kernel.length; c++) {
          row.push(kernel[r][c] == value ? 1 : 0);
        }
        result.push(row);
      }
      return(result);
    };
    
    var HitOrMiss = function(image, se1, se2) {
      // perform hitOrMiss transform
      var e1 = image.reduceNeighborhood(ee.Reducer.min(), se1);
      var e2 = image.not().reduceNeighborhood(ee.Reducer.min(), se2);
      return(e1.and(e2));
    };
              
    var Skeletonize = function(image, iterations, method) {
      // perform skeletonization
      // skeletonization implementation from Donchyts et al., 2016
        var se1w = [[2, 2, 2],
                    [0, 1, 0],
                    [1, 1, 1]];
        if(method == 2) {
          se1w = [[2, 2, 2],
                  [0, 1, 0],
                  [0, 1, 0]];
        }
        var se11 = ee.Kernel.fixed(3, 3, SplitKernel(se1w, 1));
        var se12 = ee.Kernel.fixed(3, 3, SplitKernel(se1w, 2));
        var se2w = [[2, 2, 0],
                    [2, 1, 1],
                    [0, 1, 0]];
        if(method == 2) {
          se2w = [[2, 2, 0],
                  [2, 1, 1],
                  [0, 1, 1]];
        }
        var se21 = ee.Kernel.fixed(3, 3, SplitKernel(se2w, 1));
        var se22 = ee.Kernel.fixed(3, 3, SplitKernel(se2w, 2));
        var result = image;
        for(var i = 0; i < iterations; i++) {
          for(var j=0; j<4; j++) { // rotate kernels
            result = result.subtract(HitOrMiss(result, se11, se12));
            se11 = se11.rotate(1);
            se12 = se12.rotate(1);
            result = result.subtract(HitOrMiss(result, se21, se22));
            se21 = se21.rotate(1);
            se22 = se22.rotate(1);
          }
        }
        return(result.rename(['clRaw']));
      };
      
    var CalcOnePixelWidthCenterline = function(img, GM, hGrad) {
      /***
      calculate the 1px centerline from:
      1. fast distance transform of the river banks
      2. gradient of the distance transform, mask areas where gradient greater than a threshold hGrad
      3. apply skeletonization twice to get a 1px centerline
      thresholding gradient map inspired by Pavelsky and Smith., 2008
      ***/
      // var DM = getDistanceMap(img, neighborhoodSize, scale);
      // var g = CalcGradientMap(DM, gradMethod);
      
      var imgD2 = img.focal_max(1.5, 'circle', 'pixels', 2);
      var cl = ee.Image(GM).mask(imgD2).lte(hGrad).and(img);
      // apply skeletonization twice
      var cl1px = Skeletonize(cl, 2, 1);
  
    return(cl1px);
  };
  
  var ExtractEndpoints = function(CL1px) {
    // calculate end points in the one pixel centerline
    var se1w = [[0, 0, 0],
              [2, 1, 2],
              [2, 2, 2]];
    var se11 = ee.Kernel.fixed(3, 3, SplitKernel(se1w, 1));
    var se12 = ee.Kernel.fixed(3, 3, SplitKernel(se1w, 2));
    var result = CL1px;
    // the for loop removes the identified endpoints from the imput image
    for(var i=0; i<4; i++) { // rotate kernels
      result = result.subtract(HitOrMiss(result, se11, se12));
      se11 = se11.rotate(1);
      se12 = se12.rotate(1);
    }
    var endpoints = CL1px.subtract(result);
    return(endpoints);
  };
  
  var ExtractCorners = function(CL1px) {
    /* calculate corner points in the one pixel centerline */
    var se1w = [[2, 2, 0],
              [2, 1, 1],
              [0, 1, 0]];
    var se11 = ee.Kernel.fixed(3, 3, SplitKernel(se1w, 1));
    var se12 = ee.Kernel.fixed(3, 3, SplitKernel(se1w, 2));
    var result = CL1px;
    // the for loop removes the identified corners from the imput image
    for(var i=0; i<4; i++) { // rotate kernels
      result = result.subtract(HitOrMiss(result, se11, se12));
      se11 = se11.rotate(1);
      se12 = se12.rotate(1);
    }
    var cornerPoints = CL1px.subtract(result);
    return(cornerPoints);
  };
  
  var CleanCenterline = function(cl1px, maxBranchLengthToRemove, rmCorners) {
    /*** clean the 1px centerline:
    1. remove branches
    2. remove corners to insure 1px width (optional)
    ***/
    var nearbyPoints = cl1px.mask(cl1px).reduceNeighborhood({
      reducer: ee.Reducer.count(),
      kernel: ee.Kernel.circle(1.5),
      skipMasked: true});
    var endsByNeighbors = nearbyPoints.lte(2);
    var joints = nearbyPoints.gte(4);
    var costMap = cl1px.mask(cl1px).updateMask(joints.not()).cumulativeCost({
      source: endsByNeighbors.mask(endsByNeighbors),
      maxDistance: maxBranchLengthToRemove,
      geodeticDistance: false});
    var branchMask = costMap.gte(0).unmask(0);
    var cl1Cleaned = cl1px.updateMask(branchMask.not()); // mask short branches;
    var ends = ExtractEndpoints(cl1Cleaned);
    cl1Cleaned = cl1Cleaned.updateMask(ends.not());
    if (rmCorners) {
      var corners = ExtractCorners(cl1Cleaned);
      cl1Cleaned = cl1Cleaned.updateMask(corners.not());
    }
    return(cl1Cleaned);
  };
    
    var CalculateCenterline = function(imgIn, MAXDISTANCE_BRANCH_REMOVAL) {
    
    MAXDISTANCE_BRANCH_REMOVAL = AssignDefault(MAXDISTANCE_BRANCH_REMOVAL, 500);
    
    var crs = imgIn.get('crs');
    var scale = ee.Number(30);
    var riverMask = imgIn.select('riverMask');
    var distM = CalcDistanceMap(riverMask, 256, scale);
    var gradM = CalcGradientMap(distM, 2, scale);
    var cl1 = CalcOnePixelWidthCenterline(riverMask, gradM, 0.9);
    var cl1Cleaned1 = CleanCenterline(cl1, MAXDISTANCE_BRANCH_REMOVAL, true);
    var cl1px = CleanCenterline(cl1Cleaned1, MAXDISTANCE_BRANCH_REMOVAL, false);
    
    var imgOut = imgIn
    .addBands(cl1px.rename('cleanedCL'))
    .addBands(cl1.rename('rawCL'))
    .addBands(gradM.rename('gradientMap'))
    .addBands(distM.rename('distanceMap'));
    
    return(imgOut);
  };
  
  imgOut = CalculateCenterline(imgOut, MAXDISTANCE_BRANCH_REMOVAL);
  
  /* function to calculate orthogonal directions for 1px-width centerline */
  var CalculateAngle = function(clCleaned) {
    // calculate the orthogonal direction of each pixel of the centerline
    var w3 = (ee.Kernel.fixed(9, 9, [
    [135.0, 126.9, 116.6, 104.0, 90.0, 76.0, 63.4, 53.1, 45.0],
    [143.1, 0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0, 36.9],
    [153.4, 0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0, 26.6],
    [166.0, 0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0, 14.0],
    [180.0, 0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0, 1e-5],
    [194.0, 0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0, 346.0],
    [206.6, 0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0, 333.4],
    [216.9, 0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0, 323.1],
    [225.0, 233.1,  243.4,  256.0,  270.0,  284.0,  296.6,  306.9, 315.0]]));
    var combinedReducer = ee.Reducer.sum().combine(ee.Reducer.count(), null, true);
    var clAngle = (clCleaned.mask(clCleaned)
        .rename(['clCleaned'])
        .reduceNeighborhood({
        reducer: combinedReducer,
        kernel: w3,
        inputWeight: 'kernel',
        skipMasked: true}));
    // ## mask calculating when there are more than two inputs into the angle calculation
    var clAngleNorm = (clAngle
        .select('clCleaned_sum')
        .divide(clAngle.select('clCleaned_count'))
        .mask(clAngle.select('clCleaned_count').gt(2).not()));
    // ## if only one input into the angle calculation, rotate it by 90 degrees to get the orthogonal
    clAngleNorm = (clAngleNorm
        .where(clAngle.select('clCleaned_count').eq(1), clAngleNorm.add(ee.Image(90))));
    return(clAngleNorm.rename(['orthDegree']));
  }; 
  
  // calculate orthogonal direction of the centerline
  var CalculateOrthAngle = function(imgIn) {
    var cl1px = imgIn.select('cleanedCL');
    var angle = CalculateAngle(cl1px);
    var imgOut = imgIn.addBands(angle);
    return(imgOut);
  };
  
  imgOut = CalculateOrthAngle(imgOut);
  print(imgOut);
  Map.addLayer(imgOut);
  
  /* function to calculate river width based on river mask and centerline orthogonal direction */
  var GetWidth = function(clAngleNorm, segmentInfo, endInfo, DM, crs, bound, scale, sceneID) {
    // """calculate the width of the river at each centerline pixel, measured according to the orthgonal direction of the river
    // """
    var GetXsectionEnds = function(f) {
      // calculate the ends of the cross-sectional lines
      var xc = ee.Number(f.get('x'));
      var yc = ee.Number(f.get('y'));
      var orthRad = ee.Number(f.get('angle')).divide(180).multiply(Math.PI); //convert angle in degree to radiance
    
      var halfWidth = ee.Number(f.get('toBankDistance')).multiply(1.5);
      var cosRad = halfWidth.multiply(orthRad.cos());
      var sinRad = halfWidth.multiply(orthRad.sin());
      var p1 = ee.Geometry.Point([xc.add(cosRad), yc.add(sinRad)], crs);
      var p2 = ee.Geometry.Point([xc.subtract(cosRad), yc.subtract(sinRad)], crs);
    
      var xlEnds = (ee.Feature(ee.Geometry.MultiPoint([p1, p2]), {
          'xc': xc,
          'yc': yc,
          'longitude': f.get('lon'),
          'latitude': f.get('lat'),
          'orthogonalDirection': orthRad,
          'MLength': halfWidth.multiply(2),
          'p1': p1,
          'p2': p2,
          'crs': crs,
          'image_id': sceneID
          }));
    
      return(xlEnds);
    };
    var SwitchGeometry = function(f) {
      // switch the geometry to cross-sectional buffered lines
      return(f
      .setGeometry(ee.Geometry.LineString({coords: [f.get('p1'), f.get('p2')], geodesic: false})) 
      .set('p1', null).set('p2', null)); //# remove p1 and p2
    };
    // ## convert centerline image to a featurecollection. prepare for map function
    var clPoints = (clAngleNorm.rename(['angle'])
    .addBands(ee.Image.pixelCoordinates(crs))
    .addBands(ee.Image.pixelLonLat().rename(['lon', 'lat']))
    .addBands(DM.rename(['toBankDistance']))
    .sample({
        region: bound,
        scale: scale,
        projection: null,
        factor: 1,
        dropNulls: true
    }));
    // ## calculate the cross-section lines, returning a featureCollection
    var xsectionsEnds = clPoints.map(GetXsectionEnds);
    // ## calculate the flags at the xsection line end points
    var endStat = (endInfo.reduceRegions({
        collection: xsectionsEnds,
        reducer: ee.Reducer.anyNonZero().combine(ee.Reducer.count(), null, true), //# test endpoints type 1. if in water or 2. if extends over the image bound
        scale: scale,
        crs: crs}));
    // ## calculate the width of the river and other flags along the xsection lines
    var xsections1 = endStat.map(SwitchGeometry);
    var combinedReducer = ee.Reducer.mean();
    var xsections = (segmentInfo.reduceRegions({
        collection: xsections1,
        reducer: combinedReducer,
        scale: scale,
        crs: crs}));
    return(xsections);
  };
  
  var prepExport = function(f) {
    f = f.set({
      'width': ee.Number(f.get('MLength')).multiply(f.get('channelMask')),
      'endsInWater': ee.Number(f.get('any')).eq(1),
      'endsOverEdge': ee.Number(f.get('count')).lt(2)
    });
    
    var fOut = ee.Feature(ee.Geometry.Point([f.get('longitude'), f.get('latitude')]), {})
    .copyProperties(f, null, ['any', 'count', 'MLength', 'xc', 'yc', 'channelMask']);
    
    return(fOut);
  };
  
  // export widths
  var CalculateWidth = function(imgIn) {
    var crs = ee.String('EPSG:4326');
    var scale = ee.Number(30);
    var imgId = "classification_2023";
    var bound = imgIn.select('riverMask').geometry();
    var angle = imgIn.select('orthDegree');
    var dem = ee.Image("MERIT/DEM/v1_0_3");
    var infoEnds = imgIn.select('riverMask');
    var infoExport = imgIn.select('channelMask')
    //.addBands(imgIn.select('^flag.*'))
    .addBands(dem.rename('flag_elevation'));
    var dm = imgIn.select('distanceMap');
    
    var widths = GetWidth(angle, infoExport, infoEnds, dm, crs, bound, scale, imgId).map(prepExport);
    return(widths);
  };
  
  var widthOut = CalculateWidth(imgOut);
    
    return(widthOut);
  };
  
  return(tempFUN);
};

/////////////// apply
// set AOI
var aoi = geometry;

// set parameters
var rwc = rwGenSR(
    // max distance:
    4000,
    // fill size:
    333,
    // distance brach removal:
    500,
    // aoi:
    aoi);

var imageId = "classification_2023";
var widths = rwc(imageId);

// remove the geometry before exporting the width as CSV file
widths = widths.map(function(f) {return(f.setGeometry(null))});
print(widths.first());

// export the result as a CSV file into Google drive
Export.table.toDrive({
  collection: widths,
  description: imageId,
  folder: "widths",
  fileNamePrefix: imageId,
  fileFormat: "CSV"});


