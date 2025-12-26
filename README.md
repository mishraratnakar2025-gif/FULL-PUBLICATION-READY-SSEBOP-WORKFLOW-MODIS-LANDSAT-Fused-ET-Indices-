# FULL-PUBLICATION-READY-SSEBOP-WORKFLOW-MODIS-LANDSAT-Fused-ET-Indices-
FULL PUBLICATION-READY SSEBOP WORKFLOW   MODIS + LANDSAT + Fused ET   Indices, Uncertainty, Confidence, Quadrants, Chart



/********************************************************
  FULL PUBLICATION-READY SSEBOP WORKFLOW
  MODIS + LANDSAT + Fused ET
  Indices, Uncertainty, Confidence, Quadrants, Charts
********************************************************/

/* ================================
   1️⃣ ROI
================================ */
var roi = ee.FeatureCollection('projects/eng-charge-478806-c6/assets/GOPAD_BASIN');
Map.centerObject(roi, 7);

/* ================================
   2️⃣ SAFE HELPER FUNCTIONS
================================ */
function safeMean(ic, band){
  var img = ic.select(band).reduce(ee.Reducer.mean());
  return ee.Image(ee.Algorithms.If(img.bandNames().size().gt(0), img.rename(band), ee.Image.constant(0).rename(band)));
}
function safeStd(ic, band){
  var img = ic.select(band).reduce(ee.Reducer.stdDev());
  return ee.Image(ee.Algorithms.If(img.bandNames().size().gt(0), img.rename(band), ee.Image.constant(0).rename(band)));
}
function safeImage(img, band){
  return ee.Image(ee.Algorithms.If(img.bandNames().contains(band), img, ee.Image.constant(0).rename(band)));
}

/* ================================
   3️⃣ MODIS ET 2001–2023
================================ */
var modisET = ee.ImageCollection('MODIS/061/MOD16A2')
  .filterBounds(roi)
  .filterDate('2001-01-01','2023-12-31')
  .select('ET')
  .map(function(img){ return img.multiply(0.1).rename('ET'); });

var years = ee.List.sequence(2001,2023);
var annualET_MODIS = ee.ImageCollection(
  years.map(function(y){
    var etYear = modisET.filter(ee.Filter.calendarRange(y,y,'year'));
    return safeMean(etYear,'ET').set('year',y);
  })
);

/* ================================
   4️⃣ CHIRPS SPI-12
================================ */
var rain = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
  .filterBounds(roi)
  .filterDate('2001-01-01','2023-12-31')
  .select('precipitation')
  .map(function(img){ return img.rename('P').set('year', ee.Date(img.get('system:time_start')).get('year')); });

var annualRain = ee.ImageCollection(
  years.map(function(y){
    var ic = rain.filter(ee.Filter.eq('year', y));
    return safeMean(ic,'P').set('year',y);
  })
);

var rainMean = safeMean(annualRain,'P');
var rainStd  = safeStd(annualRain,'P').max(0.001);

var SPI = annualRain.map(function(img){
  return img.subtract(rainMean).divide(rainStd).rename('SPI').set('year',img.get('year'));
});

/* ================================
   5️⃣ Drought & Normal Years
================================ */
var droughtYears = SPI.filter(ee.Filter.lt('SPI',-1)).aggregate_array('year');
var normalYears  = SPI.filter(ee.Filter.gt('SPI',-0.5)).aggregate_array('year');

var ET_drought = annualET_MODIS.filter(ee.Filter.inList('year',droughtYears));
var ET_normal  = annualET_MODIS.filter(ee.Filter.inList('year',normalYears));

var meanET_drought = safeMean(ET_drought,'ET');
var meanET_normal  = safeMean(ET_normal,'ET').max(0.001);
var stdET_drought  = safeStd(ET_drought,'ET').max(0.001);

/* ================================
   6️⃣ ET INDICES
================================ */
meanET_drought = ee.Image(meanET_drought).rename('ET_mean_drought').clip(roi);
meanET_normal  = ee.Image(meanET_normal).rename('ET_mean_normal').clip(roi);
stdET_drought  = ee.Image(stdET_drought).rename('ET_std_drought').clip(roi);

var ET_resistance = meanET_drought.divide(meanET_normal.max(0.001)).rename('Resistance').clip(roi);
var ET_resilience = meanET_normal.divide(stdET_drought.max(0.001)).rename('Resilience').clip(roi);
var ET_RTI = stdET_drought.divide(meanET_normal.max(0.001)).rename('RTI').clip(roi);
var ET_vulnerability = ET_RTI.multiply(stdET_drought).rename('Vulnerability').clip(roi);
var ET_resilience_cond = meanET_normal.subtract(meanET_drought).divide(meanET_drought.max(0.001)).rename('Resilience_cond').clip(roi);
var ET_RTI_months = ET_RTI.multiply(12).rename('RTI_months').clip(roi);

/* ================================
   7️⃣ ET UNCERTAINTY & CONFIDENCE
================================ */
var annualET_safe = annualET_MODIS.map(function(img){ return img.toFloat().updateMask(img.mask().or(1)); });
var ET_mean_all = annualET_safe.select('ET').reduce(ee.Reducer.mean()).rename('ET_mean').clip(roi);
var ET_std_all  = annualET_safe.select('ET').reduce(ee.Reducer.stdDev()).rename('ET_std').clip(roi);
var ET_CV = ET_std_all.divide(ET_mean_all.max(0.001)).rename('ET_CV').clip(roi);

var stats = ET_CV.reduceRegion({reducer:ee.Reducer.minMax(), geometry:roi, scale:500, maxPixels:1e13});
var ET_uncertainty = ET_CV.subtract(ee.Number(stats.get('ET_CV_min')))
  .divide(ee.Number(stats.get('ET_CV_max')).subtract(ee.Number(stats.get('ET_CV_min'))).max(0.001))
  .rename('Uncertainty').toFloat().clip(roi);
var ET_confidence = ee.Image(1).subtract(ET_uncertainty).rename('Confidence').toFloat().clip(roi);

/* ================================
   8️⃣ Resistance–Resilience Quadrants
================================ */
var Rm = ET_resistance.reduceRegion({reducer:ee.Reducer.mean(), geometry:roi, scale:500, maxPixels:1e13}).getNumber('Resistance');
var Sm = ET_resilience.reduceRegion({reducer:ee.Reducer.mean(), geometry:roi, scale:500, maxPixels:1e13}).getNumber('Resilience');
var quadrant = ET_resistance.gt(Rm).multiply(2).add(ET_resilience.gt(Sm)).rename('Quadrant').clip(roi);

/* ================================
   9️⃣ MODIS–Landsat FUSION 2019–2023
================================ */
var landsat = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filterBounds(roi)
  .filterDate('2019-01-01','2023-12-31')
  .map(function(img){
    var qa = img.select('QA_PIXEL');
    var mask = qa.bitwiseAnd(1<<3).eq(0).and(qa.bitwiseAnd(1<<4).eq(0));
    var LST = img.select('ST_B10').multiply(0.00341802).add(149.0).rename('LST');
    return LST.updateMask(mask).copyProperties(img,['system:time_start']);
  });

var LSTmin = landsat.reduce(ee.Reducer.percentile([5])).rename('LSTmin');
var LSTmax = landsat.reduce(ee.Reducer.percentile([95])).rename('LSTmax');

var EF = landsat.map(function(img){
  var ef = LSTmax.subtract(img).divide(LSTmax.subtract(LSTmin)).rename('EF');
  return ef.clamp(0,1).copyProperties(img,['system:time_start']);
});

var Rn = 15;   // MJ m-2 day-1
var lambda = 2.45; // MJ kg-1
var ET_Landsat = EF.map(function(img){ return img.multiply(Rn).divide(lambda).rename('ET').copyProperties(img,['system:time_start']); });

var annualET_LS = ee.ImageCollection(
  ee.List.sequence(2019,2023).map(function(y){
    var ic = ET_Landsat.filter(ee.Filter.calendarRange(y,y,'year'));
    return safeMean(ic,'ET').multiply(365).rename('ET').set('year',y);
  })
);

var MODIS_recent = annualET_MODIS.filter(ee.Filter.inList('year',[2019,2020,2021,2022,2023]));
var fusedET = ee.ImageCollection(
  ee.List.sequence(2019,2023).map(function(y){
    var mod = MODIS_recent.filter(ee.Filter.eq('year',y)).first();
    var ls  = annualET_LS.filter(ee.Filter.eq('year',y)).first();
    return ee.Image(mod).multiply(0.6).add(ee.Image(ls).multiply(0.4)).rename('ET').set('year',y);
  })
);
var fusedET_mean = safeMean(fusedET,'ET').clip(roi);
Map.addLayer(fusedET_mean,{min:500,max:1200,palette:['white','cyan','blue']},'Fused MODIS–Landsat ET');

/* ================================
   🔹 EXPORTS
================================ */
var exportList = [
  {img:ET_resistance, name:'ET_Resistance_500m'},
  {img:ET_resilience, name:'ET_Resilience_500m'},
  {img:ET_vulnerability, name:'ET_Vulnerability_500m'},
  {img:ET_RTI_months, name:'ET_RTI_500m'},
  {img:ET_resilience_cond, name:'ET_SPI_SPEI_Resilience_500m'},
  {img:ET_confidence, name:'ET_Confidence_500m'},
  {img:quadrant, name:'ET_Quadrants_500m'},
  {img:fusedET_mean, name:'Fused_ET_MODIS_Landsat_250m'}
];
exportList.forEach(function(o){
  Export.image.toDrive({
    image: o.img,
    description: o.name,
    region: roi,
    scale:500,
    crs:'EPSG:4326',
    maxPixels:1e13
  });
});

/* ================================
   🔹 ANNUAL ET TIME SERIES CHART
================================ */
function setTimeStart(ic){ return ic.map(function(img){ return img.set('system:time_start', ee.Date.fromYMD(ee.Number(img.get('year')),1,1).millis()); }); }
var MODIS_ts = setTimeStart(annualET_MODIS);
var LS_ts    = setTimeStart(annualET_LS);
var Fused_ts = setTimeStart(fusedET);

var combinedChart = ui.Chart.image.seriesByRegion({
  imageCollection: MODIS_ts.merge(LS_ts).merge(Fused_ts),
  band: 'ET',
  regions: roi,
  reducer: ee.Reducer.mean(),
  scale: 500,
  seriesProperty: 'year'
}).setOptions({
  title:'Annual ET Comparison (MODIS, Landsat, Fused)',
  hAxis:{title:'Year'},
  vAxis:{title:'ET (mm/year)'},
  lineWidth:2,
  pointSize:4,
  colors:['blue','green','red']
});
print(combinedChart);

/* ================================
   🔹 LEGENDS
================================ */
function addLegend(title, palette, min, max) {
  var legend = ui.Panel({style:{position:'bottom-left', padding:'8px 15px', backgroundColor:'white'}});
  var legendTitle = ui.Label(title, {fontWeight:'bold', fontSize:'14px', margin:'0 0 4px 0'});
  legend.add(legendTitle);
  var lon = ee.Image.pixelLonLat().select('latitude');
  var gradient = lon.multiply((max-min)/100.0).add(min);
  var colorBar = ui.Thumbnail({image: gradient.visualize({min:min, max:max, palette:palette}), style:{stretch:'horizontal', height:'12px', margin:'0 0 4px 0'}});
  legend.add(colorBar);
  var panel = ui.Panel({layout: ui.Panel.Layout.flow('horizontal')});
  panel.add(ui.Label(min));
  panel.add(ui.Label((min+max)/2, {margin: '0 45px'}));
  panel.add(ui.Label(max));
  legend.add(panel);
  Map.add(legend);
}

addLegend('ET Resistance', ['red','yellow','green'], 0, 1);
addLegend('ET Resilience', ['red','yellow','green'], 0, 3);
addLegend('RTI (months)', ['green','yellow','orange','red'], 0, 12);
addLegend('ET Vulnerability', ['yellow','red'], 0, 10);
addLegend('SPI/SPEI Resilience', ['red','yellow','green'], 0, 3);
addLegend('ET Uncertainty', ['darkgreen','yellow','red'], 0, 1);
addLegend('ET Confidence', ['red','yellow','darkgreen'], 0, 1);
addLegend('Resistance–Resilience Quadrant', ['red','yellow','blue','darkgreen'], 0, 3);
addLegend('Fused MODIS–Landsat ET', ['white','cyan','blue'], 500, 1200);
