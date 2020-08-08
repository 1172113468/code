import numpy as np
from scipy.spatial import Voronoi,voronoi_plot_2d
from osgeo import osr
import scipy
import os
import sys
import   pandas  as pd
import math
from math import sin,cos
from osgeo import ogr,osr,gdal
import datetime
start = datetime.datetime.now()

def getlayer(filename):
    gdal.SetConfigOption("SHAPE_ENCODING", "")
    driver = ogr.GetDriverByName('ESRI Shapefile')
    dataSource = driver.Open(filename,0)
    geosrs = osr.SpatialReference()
    print(geosrs)
    if dataSource is None:
        print('File cannot be opened！')
    return dataSource

def get_features(layer):
    layer.ResetReading()
    pointss=[]
    feat = layer.GetNextFeature()
    while feat:
        geom = feat.GetGeometryRef()
        x=geom.GetX()
        y=geom.GetY()
        pointss.append((x,y))
        feat = layer.GetNextFeature()
    layer.ResetReading()
    return pointss

def setXY(XY):
    n = len(XY)
    CC = (XY[0:n-1]+XY[1:n])/2
    return CC

def setgrid(extent):
    dn = 1000
    width = extent[1] - extent[0]
    height = extent[3] - extent[2]
    Y = np.linspace(extent[2],extent[3],dn+1)

    rsize = scipy.diff(Y)[0]
    Y = setXY(Y)
    Y=np.expand_dims(Y,axis=1)
    
    
    N_x = int(np.ceil(width/rsize))
    Y=np.repeat(Y,N_x,axis=1)
    

    X = np.zeros((N_x+1,))
    X[0] = extent[0]
    for i in range(N_x):
        X[i+1] = X[i]+rsize
    X = setXY(X)
    X=np.expand_dims(X,axis=0)
    X=np.repeat(X,dn,axis=0)
    return X,Y,dn,N_x,rsize

def getpoints(filename):
    print(filename)
    pointss=[]
    dataSource=getlayer(filename)
    layer = dataSource.GetLayer(0)
    extent2 = layer.GetExtent()
    layer.ResetReading()
    feat = layer.GetNextFeature()
    while feat:
        geom = feat.GetGeometryRef()
        x=geom.GetX()
        y=geom.GetY()
        pointss.append((x,y))
        feat = layer.GetNextFeature()
    layer.ResetReading()
    dataSource.Destroy()
    return pointss,extent2
    
def cal_d(points,X,Y):
    n = len(points)
    M = np.ones(X.shape)*float("inf")
    ot = np.zeros(X.shape)
    for i in range(n):
        p = points[i]
        cg = np.sqrt((X-p[0])**2+(Y-p[1])**2)
        mask = np.less(cg,M)
        if i==0:
            M = cg*mask
        else:
            M = M*(~mask)+cg*mask
        ot = ot*(~mask) + (i+1)*mask
    return ot

def writeASCII(path,ot,extent,nrows,ncols,rsize,nodata):
    head = {}
    head[0] = 'ncols         '+str(ncols)
    head[1] = 'nrows         '+str(nrows)
    head[2] = 'xllcorner     '+str(extent[0])
    head[3] = 'yllcorner     '+str(extent[2])
    head[4] = 'cellsize      '+str(rsize)
    head[5] = 'NODATA_value  '+str(nodata)
    f = open(path,'w+')
    n = ot.shape[0]
    m = ot.shape[1]
    for i in range(len(head)):
        f.write(head[i])
        f.write('\n')

    i = n-1
    while i>=0:
        string=''
        for j in range(m):
            cg = ot[i,j]
            if j<m-1:
                string = string + str(cg) +' '
            else:
                string = string + str(cg)
        f.write(string)
        f.write('\n')
        i = i-1
    f.close()

def  get_voronoi(mask,filename,out_put,path):
    dataSource=getlayer(mask)
    layer = dataSource.GetLayer(0)
    gs = layer.GetSpatialRef()
    extent1 = layer.GetExtent()
    points,extent2 = getpoints(filename)
    x_min=min([extent1[0], extent2[0]])
    x_max=max([extent1[1], extent2[1]])
    y_min=min([extent1[2], extent2[2]])
    y_max=max([extent1[3], extent2[3]])
    extent = (x_min,x_max,y_min,y_max)
    X,Y,nrows,ncols,rsize=setgrid(extent)


    ot = cal_d(points,X,Y)

    oDriver =ogr.GetDriverByName("ESRI Shapefile")
    oDS =oDriver.CreateDataSource(out_put)
    oLayer =oDS.CreateLayer("vion_polygon", gs, ogr.wkbPolygon, [])

    oFieldID =ogr.FieldDefn("Input_FID", ogr.OFTInteger)
    oLayer.CreateField(oFieldID, 1)
    oFieldName =ogr.FieldDefn("ID", ogr.OFTInteger)
    oLayer.CreateField(oFieldName, 2)
    nodata=-9999
    ot = ot.astype('int32')
    writeASCII(path,ot,extent,nrows,ncols,rsize,nodata)
    dataset=gdal.Open(path)
    im_proj = dataset.GetProjection()
    band=dataset.GetRasterBand(1)
    gdal.Polygonize(band,None,oLayer,0)
    oDS.Destroy()

class SchoolDistrictDivision():
      def rad(self,d):
          return d * math.pi / 180.0

      def getDistance(self,lat1, lng1, lat2, lng2):
          EARTH_REDIUS = 6378.137
          radLat1 = self.rad(lat1)
          radLat2 = self.rad(lat2)
          a = radLat1 - radLat2
          b = self.rad(lng1) - self.rad(lng2)
          s = 2 * math.asin(math.sqrt(math.pow(sin(a / 2), 2) + cos(radLat1) * cos(radLat2) * math.pow(sin(b / 2), 2)))
          s = s * EARTH_REDIUS
          # print("distance=", s)
          return s

      def get_fields(self,lyr):
          field_name = [field.name for field in lyr.schema]
          print(field_name)
          C_CLASS = []
          for feature0 in lyr:
              C_FULL = []
              for field in field_name:
                  C_FULL.append(feature0.GetFieldAsString(field))
              geom0 = feature0.GetGeometryRef()
              C_FULL.append(geom0.GetEnvelope()[0])
              C_FULL.append(geom0.GetEnvelope()[2])
              C_FULL.append(geom0.ExportToWkt())
              # print(geom0.GetEnvelope()[0])
              C_CLASS.append(C_FULL)
          # print(C_CLASS)
          field_name.append("X")
          field_name.append("Y")
          field_name.append('WKT')
          df = pd.DataFrame(C_CLASS, columns=field_name)
          lyr.ResetReading()
          return df

      def get_fields1(self,lyr):
          field_name = [field.name for field in lyr.schema]
          print(field_name)
          C_CLASS = []
          for feature0 in lyr:
              C_FULL = []
              for field in field_name:
                  C_FULL.append(feature0.GetFieldAsString(field))
              geom0 = feature0.GetGeometryRef()
              C_FULL.append(geom0.Centroid().GetX())
              C_FULL.append(geom0.Centroid().GetY())
              C_FULL.append(geom0.ExportToWkt())
              C_CLASS.append(C_FULL)
          # print(C_CLASS)
          field_name.append("X")
          field_name.append("Y")
          field_name.append('WKT')
          df = pd.DataFrame(C_CLASS, columns=field_name)
          lyr.ResetReading()
          return df

      def get_shp(self,a, b, type, spatial_ref):
          # spatial_ref = osr.SpatialReference()
          # spatial_ref.ImportFromEPSG(4326)
          data = a.values
          # print(data)
          colname = list(a.columns)
          # print(colname)
          # print(data[0][0])
          tem1 = driver.CreateDataSource(b)
          # print(tem1)
          TLayer1 = tem1.CreateLayer("TestPoint", spatial_ref, type_dict[type], [])
          # print(TLayer1)
          for field in colname:
              if field != 'WKT':
                  oID = ogr.FieldDefn(field, ogr.OFTString)
                  oID.SetWidth(100)
                  TLayer1.CreateField(oID, 1)
          oDefn = TLayer1.GetLayerDefn()
          # print(oDefn)
          # print(len(data))
          for i in range(len(data)):
              # print(i)
              oFeatureTriangle = ogr.Feature(oDefn)
              for j in range(len(colname) - 1):
                  # print(j)
                  oFeatureTriangle.SetField(j, data[i, j])
              geomTriangle = ogr.CreateGeometryFromWkt(data[i][len(colname) - 1])
              oFeatureTriangle.SetGeometry(geomTriangle)
              TLayer1.CreateFeature(oFeatureTriangle)
          tem1.Destroy()


      def calculate(self,original_grid,original_school,tyson_polygon,temp_txt,new_grid,new_school,temp_result,end_result,n,field_school_name,field_school_zrs,field_grid_rks):
          global type_dict, driver, ref3

          type_dict = {1: ogr.wkbPoint, 2: ogr.wkbLineString, 3: ogr.wkbPolygon, 4: ogr.wkbMultiPoint,
                       5: ogr.wkbMultiLineString, 6: ogr.wkbMultiPolygon, 7: ogr.wkbGeometryCollection}
          driver = ogr.GetDriverByName('ESRI Shapefile')
          gdal.SetConfigOption("SHAPE_ENCODING", "")

          result = []
          col = []
          for i in range(int(n)-2):
              print(i)
              if i == 0:
                  mask =original_grid
                  filename = original_school
                  out_put = tyson_polygon
                  path = temp_txt
                  get_voronoi(mask, filename, out_put, path)

                  filename_v1 = tyson_polygon
                  dataSource1 = driver.Open(filename_v1, 0)
                  if dataSource1 is None:
                      print('File cannot be opened！')
                  layer1 = dataSource1.GetLayer(0)
                  type1 = layer1.GetGeomType()
                  ref1 = layer1.GetSpatialRef()
                  ts_df = self.get_fields(layer1)

                  filename_v2 = original_school
                  dataSource2 = driver.Open(filename_v2, 0)
                  if dataSource1 is None:
                      print('File cannot be opened！')
                  layer2 = dataSource2.GetLayer(0)
                  type2 = layer2.GetGeomType()
                  ref2 = layer2.GetSpatialRef()
                  xx_df = self.get_fields1(layer2)
                  xx_df[field_school_zrs] = xx_df[field_school_zrs].astype(float)
                  xx_df[field_school_zrs] = xx_df[field_school_zrs] * 1.29
                  ts_df.sort_values(by=['X', 'Y'], ascending=(True, True), inplace=True)
                  ts_df['sort'] = [i for i in range(len(ts_df))]
                  min_wkt = ts_df['WKT'].values[0]
                  min_value = 0
                  min_name = ''
                  col_xxdf=list(xx_df.columns)
                  name_index=col_xxdf.index(field_school_name)
                  value_index=col_xxdf.index(field_school_zrs)
                  xx_values = xx_df.values.tolist()
                  #print(xx_values)
                  for row in range(len(xx_values)):
                      if ogr.CreateGeometryFromWkt(min_wkt).Contains(
                              ogr.CreateGeometryFromWkt(xx_values[row][-1])) or ogr.CreateGeometryFromWkt(
                              min_wkt).Intersects(ogr.CreateGeometryFromWkt(xx_values[row][-1])):
                          min_value = float(xx_values[row][value_index])
                          min_name = xx_values[row][name_index]
                          break

                 
                  filename_v3 = original_grid
                  dataSource3 = driver.Open(filename_v3, 0)
                  if dataSource3 is None:
                      print('File cannot be opened！')
                  layer3 = dataSource3.GetLayer(0)
                  type3 = layer3.GetGeomType()
                  ref3 = layer3.GetSpatialRef()
                  wg_df = self.get_fields1(layer3)
                  wg_dfcopy = wg_df.copy()
                  wg_df[field_grid_rks] = wg_df[field_grid_rks].fillna(0)
                  wg_df[field_grid_rks] = wg_df[field_grid_rks].replace('', 0)
                  wg_df[field_grid_rks] = wg_df[field_grid_rks].astype(float)

                  wg_values = wg_df.values.tolist()
                  flag = []
                  for row in range(len(wg_values)):
                      if ogr.CreateGeometryFromWkt(min_wkt).Intersects(ogr.CreateGeometryFromWkt(wg_values[row][-1])):
                          flag.append(1)
                      else:
                          flag.append(0)
                  wg_df['flag'] = flag
                  select_1 = wg_df[wg_df['flag'] > 0]
                  print('min_value  RKS')
                  print(min_value)
                  print(select_1[field_grid_rks].sum())
                  print(select_1.columns)
                  print(min_value)
                  print(type(select_1[field_grid_rks].sum()))
                  select_1.sort_values(by=['Y', 'X'], ascending=(True, True), inplace=True)
                  select_1_values = select_1.values.tolist()
                  select_1_col = list(select_1.columns)
                  select_1_vindex = select_1_col.index(field_grid_rks)
 
                  s = 0.0
                  print(select_1.columns)
                  flag_xx = []
                  for row in range(len(select_1_values)):
                      if s <= int(min_value-min_value*0.1):
                          flag_xx.append(min_name)
                          s = s + select_1_values[row][select_1_vindex]
                      else:
                          flag_xx.append(0)
                  select_1['flag_xx'] = flag_xx
                  select_2 = select_1[select_1['flag_xx'] == min_name]

                  tem_wkt = select_2['WKT']
                  del select_2['WKT']
                  select_2['WKT'] = tem_wkt
                  col = list(select_2.columns)
                  list1 = select_2.values.tolist()
                  for a in list1:
                      result.append(a)
                
                  rest_xx = xx_df[~xx_df[field_school_name].isin([min_name])]
                  tem_wkt = rest_xx['WKT']
                  del rest_xx['WKT']
                  rest_xx['WKT'] = tem_wkt
                  del rest_xx['X']
                  del rest_xx['Y']
                  self.get_shp(rest_xx, new_school, type2, ref3)

                  rest_wg = wg_df[~wg_df['WKT'].isin(select_2['WKT'].values.tolist())]
                  tem_wkt = rest_wg['WKT']
                  del rest_wg['WKT']
                  rest_wg['WKT'] = tem_wkt
                  del rest_wg['X']
                  del rest_wg['Y']
                  self.get_shp(rest_wg, new_grid, type3, ref3)

                  print(min_name)
                  dataSource1.Destroy()
                  dataSource2.Destroy()
                  dataSource3.Destroy()

              else:
                  
                  mask = new_grid
                  filename = new_school
                  out_put = tyson_polygon  #
                  path = temp_txt  
                  get_voronoi(mask, filename, out_put, path)

                  filename_v1 = tyson_polygon
                  dataSource1 = driver.Open(filename_v1, 0)
                  if dataSource1 is None:
                      print('File cannot be opened！')
                  layer1 = dataSource1.GetLayer(0)
                  type1 = layer1.GetGeomType()
                  ref1 = layer1.GetSpatialRef()
                  ts_df = self.get_fields(layer1)

                  filename_v2 = new_school
                  dataSource2 = driver.Open(filename_v2, 0)
                  if dataSource1 is None:
                      print('File cannot be opened！')
                  layer2 = dataSource2.GetLayer(0)
                  type2 = layer2.GetGeomType()
                  ref2 = layer2.GetSpatialRef()
                  xx_df = self.get_fields1(layer2)
                  xx_df[field_school_zrs] = xx_df[field_school_zrs].astype(float)
                  # xx_df['ZRS'] = xx_df['ZRS'] * 1.2
                  ts_df.sort_values(by=['X', 'Y'], ascending=(True, True), inplace=True)
                  ts_df['sort'] = [i for i in range(len(ts_df))]
                  min_wkt = ts_df['WKT'].values[0]
                  min_value = 0
                  min_name = ''
                  col_xxdf = list(xx_df.columns)
                  name_index = col_xxdf.index(field_school_name)
                  value_index = col_xxdf.index(field_school_zrs)
                  xx_values = xx_df.values.tolist()
                  print(xx_values)
                  for row in range(len(xx_values)):
                      if ogr.CreateGeometryFromWkt(min_wkt).Contains(
                              ogr.CreateGeometryFromWkt(xx_values[row][-1])) or ogr.CreateGeometryFromWkt(
                              min_wkt).Intersects(ogr.CreateGeometryFromWkt(xx_values[row][-1])):
                          min_value = xx_values[row][value_index]
                          min_name = xx_values[row][name_index]
                          break

                  
                  filename_v3 = new_grid
                  dataSource3 = driver.Open(filename_v3, 0)
                  if dataSource3 is None:
                      print('File cannot be opened！')
                  layer3 = dataSource3.GetLayer(0)
                  type3 = layer3.GetGeomType()
                  ref3 = layer3.GetSpatialRef()
                  wg_df = self.get_fields1(layer3)
                  wg_dfcopy = wg_df.copy()
                  wg_df[field_grid_rks] = wg_df[field_grid_rks].fillna(0)
                  wg_df[field_grid_rks] = wg_df[field_grid_rks].replace('', 0)
                  wg_df[field_grid_rks] = wg_df[field_grid_rks].astype(float)

                  wg_values = wg_df.values.tolist()
                  flag = []
                  sum_RKS = 0.0
                  for row in range(len(wg_values)):
                      if ogr.CreateGeometryFromWkt(min_wkt).Intersects(ogr.CreateGeometryFromWkt(wg_values[row][-1])):
                          flag.append(1)
                      else:
                          flag.append(0)
                  wg_df['flag'] = flag
                  select_1 = wg_df[wg_df['flag'] > 0]
                  sum_RKS = select_1[field_grid_rks].sum()
                  boder = 0.001
                  while sum_RKS < min_value-min_value*0.1:
                      print('while')
                      print(min_value, sum_RKS)
                      flag1 = []
                      wg_df['flag'] = None
                      min_wkt = ogr.CreateGeometryFromWkt(min_wkt).Buffer(boder).ExportToWkt()
                      for row in range(len(wg_values)):
                          if ogr.CreateGeometryFromWkt(min_wkt).Intersects(
                                  ogr.CreateGeometryFromWkt(wg_values[row][-1])):
                              flag1.append(1)
                          else:
                              flag1.append(0)
                      wg_df['flag'] = flag1
                      select_1 = wg_df[wg_df['flag'] > 0]
                      sum_RKS = select_1[field_grid_rks].sum()
                      boder = boder + 0.0001
                  print('min_value  RKS')
                  print(min_value)
                  print(select_1[field_grid_rks].sum())

                  print(select_1.columns)
                  print(min_value)
                  print(type(select_1[field_grid_rks].sum()))
                  select_1.sort_values(by=['Y', 'X'], ascending=(True, True), inplace=True)
                  select_1_values = select_1.values.tolist()
                  select_1_col = list(select_1.columns)
                  select_1_vindex = select_1_col.index(field_grid_rks)
                  
                  s = 0.0
                  print(select_1.columns)
                  flag_xx = []
                  print('min_value')
                  print(min_value)
                  print(type(min_value))
                  for row in range(len(select_1_values)):
                      if s <= min_value:
                          flag_xx.append(min_name)
                          s = s + select_1_values[row][select_1_vindex]
                      else:
                          flag_xx.append(0)
                  select_1['flag_xx'] = flag_xx
                  select_2 = select_1[select_1['flag_xx'] == min_name]

                  tem_wkt = select_2['WKT']
                  del select_2['WKT']
                  select_2['WKT'] = tem_wkt
                  list1 = select_2.values.tolist()
                  for a in list1:
                      result.append(a)
                  

                  rest_xx = xx_df[~xx_df[field_school_name].isin([min_name])]
                  tem_wkt = rest_xx['WKT']
                  del rest_xx['WKT']
                  rest_xx['WKT'] = tem_wkt
                  del rest_xx['X']
                  del rest_xx['Y']
                  dataSource2.Destroy()
                  self.get_shp(rest_xx, new_school, type2, ref3)

                  rest_wg = wg_df[~wg_df['WKT'].isin(select_2['WKT'].values.tolist())]
                  tem_wkt = rest_wg['WKT']
                  del rest_wg['WKT']
                  rest_wg['WKT'] = tem_wkt
                  del rest_wg['X']
                  del rest_wg['Y']
                  dataSource3.Destroy()
                  self.get_shp(rest_wg, new_grid, type3, ref3)
                  dataSource1.Destroy()
                  print(min_name)
                  result_end1 = pd.DataFrame(None, columns=col)
                  result_end1 = pd.DataFrame(result, columns=col)
                  self.get_shp(result_end1, temp_result, type_dict[1], ref3)

          result_end = pd.DataFrame(result, columns=col)
          del result_end['flag']

          # REST
          filename_v4 = new_school
          dataSource4 = driver.Open(filename_v4, 0)
          if dataSource4 is None:
              print('File cannot be opened！')
          layer4 = dataSource4.GetLayer(0)
          type4 = layer4.GetGeomType()
          ref4 = layer4.GetSpatialRef()
          xx_df = self.get_fields1(layer4)
          xx_df[field_school_zrs] = xx_df[field_school_zrs].astype(float)
          xx_df.sort_values(by=['X', 'Y'], ascending=(True, True), inplace=True)

         
          filename_v5 = new_grid
          dataSource5 = driver.Open(filename_v5, 0)
          if dataSource5 is None:
              print('File cannot be opened！')
          layer5 = dataSource5.GetLayer(0)
          type5 = layer5.GetGeomType()
          ref5 = layer5.GetSpatialRef()
          wg_df = self.get_fields1(layer5)
          wg_df[field_grid_rks] = wg_df[field_grid_rks].fillna(0)
          wg_df[field_grid_rks] = wg_df[field_grid_rks].replace('', 0)
          wg_df[field_grid_rks] = wg_df[field_grid_rks].astype(float)

          xx_xy = xx_df[['X', 'Y',field_school_name, field_school_zrs]].values.tolist()

          all_list = []
          # ONE
          dis = []
          WG_XY = wg_df[['X', 'Y', field_grid_rks]].values.tolist()
          for loc in WG_XY:
              dis1 = self.getDistance(loc[1], loc[0], xx_xy[0][1], xx_xy[0][0])
              dis.append(dis1)
          wg_df['dis'] = dis
          wg_df.sort_values(by=['dis'], ascending=True, inplace=True)
          wg_dis = wg_df['dis'].values.tolist()
          WG_XY = wg_df[['X', 'Y', field_grid_rks]].values.tolist()

          flag1 = []
          sum_RKS = 0
          for row in range(len(wg_dis)):
              if sum_RKS <= xx_xy[0][3]:
                  sum_RKS = sum_RKS + WG_XY[row][2]
                  flag1.append(xx_xy[0][2])
              else:
                  flag1.append(0)
          wg_df['flag_xx'] = flag1
          select_11 = wg_df[wg_df['flag_xx'] == xx_xy[0][2]]
          wg_df1 = wg_df[~wg_df['flag_xx'].isin(select_11['flag_xx'].values.tolist())]

          # TWO
          dis2 = []
          WG_XY1 = wg_df1[['X', 'Y', field_grid_rks]].values.tolist()
          for loc in WG_XY1:
              dis1 = self.getDistance(loc[1], loc[0], xx_xy[1][1], xx_xy[1][0])
              dis2.append(dis1)
          wg_df1['dis'] = dis2
          wg_df1.sort_values(by=['dis'], ascending=True, inplace=True)
          wg_dis1 = wg_df1['dis'].values.tolist()
          flag2 = []
          sum_RKS2 = 0
          WG_XY2 = wg_df1[['X', 'Y', field_grid_rks]].values.tolist()
          print(WG_XY2)
          for row in range(len(wg_dis1)):
              if sum_RKS2 <= xx_xy[1][3]:
                  sum_RKS2 = sum_RKS2 + WG_XY2[row][2]
                  flag2.append(xx_xy[1][2])
              else:
                  flag2.append('rest_net')
          wg_df1['flag_xx'] = flag2
          select_22 = wg_df1

          for a in select_11.values.tolist():
              all_list.append(a)
          for b in select_22.values.tolist():
              all_list.append(b)
          result1 = pd.DataFrame(all_list, columns=list(select_11.columns))
          del result1['dis']
          del result1['flag']
          tem_wkt = result1['WKT']
          del result1['WKT']
          result1['WKT'] = tem_wkt
          all_result = pd.concat([result_end, result1], axis=0)
          X = [ogr.CreateGeometryFromWkt(g).Centroid().GetX() for g in all_result['WKT'].values.tolist()]
          Y = [ogr.CreateGeometryFromWkt(g).Centroid().GetY() for g in all_result['WKT'].values.tolist()]
          all_result['X'] = X
          all_result['Y'] = Y
          self.get_shp(all_result, end_result, type_dict[1], ref3)
          end = datetime.datetime.now()
          print("start_time:" + start.strftime("%Y-%m-%d %H:%M:%S"))
          print("end_time:" + end.strftime("%Y-%m-%d %H:%M:%S"))

if __name__ == '__main__' :
    original_grid = 'net.shp'#Input data,the data format is shpfile,this data represents the number of students in each unit,which contains the field of field_grid_rks.
    original_school = 'school.shp'#Input data,the data format is shpfile,this data represents the student capacity in each school,which contains the field of field_school_zrs and ield_school_name.
    tyson_polygon = 'voronoi_polygon.shp'#Output data,the data format is shpfile,this data is automatically generated during the calculation process.
    temp_txt = 'temp.txt'#Output data,the data format is txt,this data is automatically generated during the calculation process.
    new_grid = 'wg.shp'#Output data,the data format is shpfile,this data is automatically generated during the calculation process.
    new_school = 'xx.shp'#Output data,the data format is shpfile,this data is automatically generated during the calculation process.
    temp_result = 'result_temp.shp'#Output data,the data format is shpfile,this data is automatically generated during the calculation process.
    end_result = 'result.shp'#Output data,the data format is shpfile,this data is the target data.
    field_school_name = 'NAME'#The field of school name.
    field_school_zrs = 'ZRS'#The field of student capacity in each school.
    field_grid_rks = 'RKS'#The field of students number in each unit.
    sdd=SchoolDistrictDivision()
    sdd.calculate(original_grid,original_school,tyson_polygon,temp_txt,new_grid,new_school,temp_result,end_result,n,field_school_name,field_school_zrs,field_grid_rks)






