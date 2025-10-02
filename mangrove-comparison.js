var L8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2"),
    l7 = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2"),
    SRTM = ee.Image("CGIAR/SRTM90_V4");

// ===== GLOBAL VARIABLES FOR CHANGE DETECTION =====
var userGeometry = null;

// ===== IMPROVED MANGROVE DETECTION FUNCTION =====
function updateMap(year) {
  // print('updating map for year ', year);

  var currentYear = year;

  // Don't run if no area selected
  if (userGeometry === null) {
    return;
  }

  // var layersToRemove = [];
  // Map.layers().forEach(function(layer) {
  //   var name = layer.getName();
  //   if (name && name.indexOf('Mangroves') !== -1) {
  //     layersToRemove.push(layer);
  //   }
  // });
  
  // // Now remove them
  // layersToRemove.forEach(function(layer) {
  //   Map.remove(layer);
  // });

  //Center the map to the region of interest using the region shapefile
  // Map.centerObject(userGeometry,11)
  Map.setOptions('satellite')
  if (year >= 2013) {
    mapL8(year);  // Use only Landsat 8 for 2013+
  } else {
    mapL7(year);  // Use Landsat 7 for older years
  }
}

//2.1) Cloud Masking
////////////////////

//Landsat data includes a 'pixel_qa' band which can be used to create 
//     a function to mask clouds

function maskClouds(image) {
  
  // Bits 3 and 5 are cloud shadow and cloud, respectively.
    var cloudShadowBitMask = ee.Number(2).pow(3).int();
    var cloudsBitMask = ee.Number(2).pow(5).int();  
    
    // Get the pixel QA band.
    var qa = image.select('QA_PIXEL');
    
     // Both flags should be set to zero, indicating clear conditions.
    var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0).and(qa.bitwiseAnd(cloudsBitMask).eq(0)); 
  
  // Return the masked image, scaled to [0, 1].
  return image.updateMask(mask).divide(10000).copyProperties(image, ["system:time_start"]);
}

//2.2) Adding Spectral Indices
///////////////////////////////

// This function maps spectral indices for Mangrove Mapping using Landsat 8 Imagery
var addIndicesL8 = function(img) {
  // NDVI
  var ndvi = img.normalizedDifference(['SR_B5','SR_B4']).rename('NDVI');
  // NDMI (Normalized Difference Mangrove Index - Shi et al 2016 - New spectral metrics for mangrove forest identification)
  var ndmi = img.normalizedDifference(['SR_B7','SR_B3']).rename('NDMI');
  // MNDWI (Modified Normalized Difference Water Index - Hanqiu Xu, 2006)
  var mndwi = img.normalizedDifference(['SR_B3','SR_B6']).rename('MNDWI');
  // SR (Simple Ratio)
  var sr = img.select('SR_B5').divide(img.select('SR_B4')).rename('SR');
  // Band Ratio 54
  var ratio54 = img.select('SR_B6').divide(img.select('SR_B5')).rename('R54');
  // Band Ratio 35
  var ratio35 = img.select('SR_B4').divide(img.select('SR_B6')).rename('R35');
  // GCVI
  var gcvi = img.expression('(NIR/GREEN)-1',{
    'NIR':img.select('SR_B5'),
    'GREEN':img.select('SR_B3')
  }).rename('GCVI');
  return img
    .addBands(ndvi)
    .addBands(ndmi)
    .addBands(mndwi)
    .addBands(sr)
    .addBands(ratio54)
    .addBands(ratio35)
    .addBands(gcvi);
};

// This function maps spectral indices for Mangrove Mapping using Landsat 7 Imagery
var addIndicesL7 = function(img) {
  // NDVI
  var ndvi = img.normalizedDifference(['SR_B4','SR_B3']).rename('NDVI');
  // NDMI (Normalized Difference Mangrove Index - Shi et al 2016 - New spectral metrics for mangrove forest identification)
  var ndmi = img.normalizedDifference(['SR_B7','SR_B2']).rename('NDMI');
  // MNDWI (Modified Normalized Difference Water Index - Hanqiu Xu, 2006)
  var mndwi = img.normalizedDifference(['SR_B2','SR_B5']).rename('MNDWI');
  // SR (Simple Ratio)
  var sr = img.select('SR_B4').divide(img.select('SR_B3')).rename('SR');
  // Band Ratio 54
  var ratio54 = img.select('SR_B5').divide(img.select('SR_B4')).rename('R54');
  // Band Ratio 35
  var ratio35 = img.select('SR_B3').divide(img.select('SR_B5')).rename('R35');
  // GCVI
  var gcvi = img.expression('(NIR/GREEN)-1',{
    'NIR':img.select('SR_B4'),
    'GREEN':img.select('SR_B2')
  }).rename('GCVI');
  return img
    .addBands(ndvi)
    .addBands(ndmi)
    .addBands(mndwi)
    .addBands(sr)
    .addBands(ratio54)
    .addBands(ratio35)
    .addBands(gcvi);
};

//2.3) Filter Landsat data by Date and Region for Landsat 7
/////////////////////////////////////////////

function mapL7(year) {
  // print('mapping l7 for year ', year);
    // Select the desired central year here
  // var year = 2009; 
  
  // Start date will be set one year before the central year
  var startDate = (year-1)+'-01-01'; 
  
  // End date will be set to one year later than the central year.
  var endDate = (year+1)+'-12-31'; 
  
  //4.3) Apply filters and masks to Landsat 7 imagery
  ///////////////////////////////////////////////////
  
  var l7 = L7.filterDate(startDate,endDate)
  // Mask for clouds and cloud shadows
      .map(maskClouds)  //We use the same function we used for Landsat 8
  //Add the indices
      .map(addIndicesL7)
      
  //4.4) Composite the Landsat image collection
  /////////////////////////////////////////////
  
  //You can composite on a per pixel, per-band basis using .median()
  // OR with quality bands like .qualityMosaic('NDVI')
  
  // Clip SRTM data to region
  var srtmClip = SRTM.clip(userGeometry);
  
  //Mask to elevations less than 65 meters
  var elevationMask = srtmClip.lt(65);

  var L7composite = l7
                // Uses the median reducer
                .median() 
                // Clips the composite to our area of interest
                .clip(userGeometry); 
  
  //4.5) Mask to areas of low elevation and high NDVI and MNDWI
  /////////////////////////////////////////////////////////////
  
  //Used the NDVI and MNDWI bands to create masks
  var L7NDVIMask = L7composite.select('NDVI').gt(0.25);
  var L7MNDWIMask = L7composite.select('MNDWI').gt(-0.50);
  
  //Apply the masks
  var L7compositeNew = L7composite
                          .updateMask(L7NDVIMask)
                          .updateMask(L7MNDWIMask)
                          .updateMask(elevationMask) //We can use the same mask as before
                          
  //4.6) Display results
  ///////////////////////
  //Select bands and parameters for visualization
  //We use bands 4, 5, and 3 instead
  var L7visPar = {bands:['SR_B4','SR_B5','SR_B3'], min: 0.7, max: 2.4}; 
  
  //Add layer to map
  Map.addLayer(L7compositeNew.clip(userGeometry), L7visPar, 'Mangroves ' + year)
}

function mapL8(year) {
  print('mapping l8 for year ', year);
  // Select the desired central year here
  // var year = 2019; 
  
  // Start date will be set one year before the central year
  var startDate = (year-1)+'-01-01'; 
  
  // End date will be set to one year later than the central year.
  var endDate = (year+1)+'-12-31'; 
  
  //2.4) Apply filters and masks to Landsat 8 imagery
  ///////////////////////////////////////////////////
  
  var l8 = L8.filterDate(startDate,endDate)
  // Mask for clouds and cloud shadows
      .map(maskClouds)
  //Add the indices
      .map(addIndicesL8)
      
  //2.5) Composite the Landsat image collection
  /////////////////////////////////////////////
  //You can composite on a per pixel, per-band basis using .median()
  // OR with quality bands like .qualityMosaic('NDVI')
  
  var composite = l8
                // Uses the median reducer
                .median() 
                // Clips the composite to our area of interest
                .clip(userGeometry); 
  
  //2.6) Mask to areas of low elevation and high NDVI and MNDWI
  /////////////////////////////////////////////////////////////
  // Clip SRTM data to region
  var srtmClip = SRTM.clip(userGeometry);
  //Mask to elevations less than 65 meters
  var elevationMask = srtmClip.lt(65);
  
  //Used the NDVI and MNDWI bands to create masks
  var NDVIMask = composite.select('NDVI').gt(0.25);
  var MNDWIMask = composite.select('MNDWI').gt(-0.50);
  
  //Apply the masks
  var compositeNew = composite
                          .updateMask(NDVIMask)
                          .updateMask(MNDWIMask)
                          .updateMask(elevationMask)
                          
  //2.7) Display results
  ///////////////////////
  //Select bands and parameters for visualization
  var visPar = {bands:['SR_B5','SR_B6','SR_B4'], min: 0.7, max: 2.4}; 
  
  // Increase maxPixels and add bestEffort
  print('Composite pixel values:', composite.select(['SR_B5','SR_B6','SR_B4']).reduceRegion({
    reducer: ee.Reducer.minMax(),
    geometry: userGeometry,
    scale: 100,  // Larger scale to reduce pixel count
    maxPixels: 1e8,  // Much higher limit
    bestEffort: true  // Let GEE use whatever scale needed
  }));
  
  
  //Add layer to map
  print('adding landsat composite L8 layer');
  Map.addLayer(compositeNew.clip(userGeometry), visPar, 'Mangroves ' + year)
}


// Create a side panel
var sidePanel = ui.Panel({
  layout: ui.Panel.Layout.flow('vertical'),
  style: {
    width: '300px',
    position: 'bottom-left',  // options: 'top-left', 'top-right', 'bottom-left', 'bottom-right'
    padding: '8px'
  }
});

// // Create a slider for years
// var yearSlider = ui.Slider({
//   min: 2014,
//   max: 2024,
//   step: 1,
//   value: 2014,
//   style: {stretch: 'horizontal'},
// });


// Draw region button
var drawButton = ui.Button({
  label: 'Draw a region',
  onClick: function() {
    Map.drawingTools().setShape('polygon'); // can be 'rectangle', 'point', etc.
    Map.drawingTools().draw();
  }
});

var instructLabel = ui.Label('Draw a region by clicking button, select 2 years to compare using sliders, then click Compare Years to compare mangrove coverage between those years.')


var year1Slider = ui.Slider({
  min: 2000,
  max: 2024,
  step: 1,
  value: 2010,
  style: {stretch: 'horizontal'}
});

var year2Slider = ui.Slider({
  min: 2000,
  max: 2024,
  step: 1,
  value: 2020,
  style: {stretch: 'horizontal'}
});


    
var year1Label = ui.Label('Year 1: 2010');
var year2Label = ui.Label('Year 2: 2020');

year1Slider.onChange(function(year) {
  year1Label.setValue('Year 1: ' + year);
});

year2Slider.onChange(function(year) {
  year2Label.setValue('Year 2: ' + year);
});

// Update button to generate comparison
var compareButton = ui.Button({
  label: 'Compare Years',
  onClick: function() {
    var layersToRemove = [];
    Map.layers().forEach(function(layer) {
      var name = layer.getName();
      if (name && name.indexOf('Mangroves') !== -1) {
        layersToRemove.push(layer);
      }
    });
    
    layersToRemove.forEach(function(layer) {
      Map.remove(layer);
    });
    var y1 = year1Slider.getValue();
    var y2 = year2Slider.getValue();
    updateMap(y1);
    updateMap(y2);
  }
});

// Add widgets to side panel
sidePanel.add(instructLabel);
sidePanel.add(drawButton);
sidePanel.add(year1Label);
sidePanel.add(year1Slider);
sidePanel.add(year2Label);
sidePanel.add(year2Slider);
sidePanel.add(compareButton);

// Add the side panel to the map
Map.add(sidePanel);

// Enable drawing tools
Map.drawingTools().setShown(false);
Map.drawingTools().layers().reset();


Map.setCenter(0, 20, 4);

Map.drawingTools().onDraw(function() {
  var layers = Map.drawingTools().layers();
  if (layers.length() > 0) {
    userGeometry = layers.get(0).toGeometry();
    print('User drew geometry:', userGeometry);

    // Clear the drawing overlay
    Map.drawingTools().layers().reset();
    
    // Add the geometry to the map as a real layer
    var geomLayer = ui.Map.Layer(userGeometry, {color: 'red'}, 'User Geometry');
    Map.layers().set(1, geomLayer);  // slot 1 so your composite stays on top (slot 0)

  }
});




