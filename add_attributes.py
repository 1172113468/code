import   pandas  as pd
import math
from osgeo import ogr,osr,gdal
#显示所有列
pd.set_option('display.max_columns', None)
#显示所有行
pd.set_option('display.max_rows', None)
#设置value的显示长度为100，默认为50
pd.set_option('max_colwidth',100)

class  add_attributes():
    def __init__(self, input_shpfile, group_id, mutifield_sum, output_csv,output_shpfile,threshold=90):
        self.input_shpfile = input_shpfile
        self.group_id = group_id
        self.mutifield_sum = mutifield_sum
        self.output_csv=output_csv
        self.output_shpfile = output_shpfile
        self.threshold = float(threshold)
    def get_fields(self,lyr):
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

    def get_convex(self,a):
        a["K"] = a["Y"] / a["X"]
        a.sort_values(by=[self.group_id])
        grouped = a.groupby(self.group_id)
        poly = []
        poly1 = []
        new_list = []
        for name, group in grouped:
            # print(name)
            group_sort = group.sort_values(by=['K', 'X'], ascending=[True, True])
            # group_list=group_sort.values.tolist()
            group_sort["cen_x"] = group_sort["X"].mean()
            group_sort["cen_y"] = group_sort["Y"].mean()
            group_sort["x_cen"] = group_sort["X"] - group_sort["cen_x"]
            group_sort["y_cen"] = group_sort["Y"] - group_sort["cen_y"]
            group_list = group_sort[["x_cen", "y_cen"]].values.tolist()
            # print(group_list)
            atan = []
            temp = []
            for gp_list in group_list:
                atan.append(math.atan2(gp_list[1], gp_list[0]))
            b = pd.DataFrame(atan, columns=["atan_tan"])
            group_sort_re = group_sort.reset_index()
            # print(group_sort_re)
            b1 = b
            # print(b1)
            c = pd.concat([group_sort_re, b1], axis=1)
            # print(c)
            c1 = c.sort_values(by=["atan_tan", 'X'], ascending=[True, True])
            c2 = c1[["X", "Y"]].values.tolist()
            boundary_value = self.threshold
            min_value = boundary_value - boundary_value*0.05
            max_value = boundary_value + boundary_value*0.05
            c1[[self.mutifield_sum]] = c1[[self.mutifield_sum]].astype(float)
            z = c1[self.mutifield_sum].sum()
            print(z)
            n = math.ceil(z / boundary_value)
            n = int(n * 5)
            v = c1[self.mutifield_sum].values.tolist()
            label = []
            label_1 = []
            i = 0
            if i <= len(group):
                for m in range(n):
                    sum_len = float(0)
                    while (sum_len <= max_value) and i < len(group):
                        # print(i)
                        if sum_len >= min_value:
                            break
                        else:
                            sum_len = sum_len + v[i]
                            label.append(name + '_' + str(m))
                            label_1.append(m)
                            i = i + 1
                c1['class1'] = label
                c1['class2'] = label_1
                max_label = int(c1['class2'].max())
                max_label_df = c1[c1['class2'] == max_label]
                sum_max = float(max_label_df[self.mutifield_sum].sum())
                if sum_max <= boundary_value * 0.45:
                    c1['class2'] = c1['class2'].replace(max_label, max_label - 1)
                    c1['class1'] = c1['class1'].replace(str(name) + '_' + str(max_label),
                                                        str(name) + '_' + str(max_label - 1))
                print('sum')
                print(type(float(max_label_df[self.mutifield_sum].sum())))
                print(name)
                print(type(int(max_label)))
                print(max_label_df)
            c1_list = c1.values.tolist()
            for list1 in c1_list:
                new_list.append(list1)
        new_df = pd.DataFrame(new_list, columns=c1.columns)
        return new_df

    def get_aopt(self,a):
        a["K"] = a["Y"] / a["X"]
        a.sort_values(by=[self.groupid])
        grouped = a.groupby(self.groupid)
        poly = []
        poly1 = []
        for name, group in grouped:
            group_sort = group.sort_values(by=['K', 'X'], ascending=[True, True])
            # group_list=group_sort.values.tolist()
            group_sort["cen_x"] = group_sort["X"].mean()
            group_sort["cen_y"] = group_sort["Y"].mean()
            group_sort["x_cen"] = group_sort["X"] - group_sort["cen_x"]
            group_sort["y_cen"] = group_sort["Y"] - group_sort["cen_y"]
            group_list = group_sort[["x_cen", "y_cen"]].values.tolist()
            # print(group_list)
            atan = []
            temp = []
            for gp_list in group_list:
                atan.append(math.atan2(gp_list[1], gp_list[0]))
            b = pd.DataFrame(atan, columns=["atan_tan"])
            group_sort_re = group_sort.reset_index()
            # print(group_sort_re)
            b1 = b
            # print(b1)
            c = pd.concat([group_sort_re, b1], axis=1)
            # print(c)
            c1 = c.sort_values(by=["atan_tan", 'X'], ascending=[True, True])
            c2 = c1[["X", "Y"]].values.tolist()
            # print(c2)
            str1 = "POLYGON(("
            for row1 in range(len(c2)):
                if row1 < len(c2) - 1:
                    str1 = str1 + str(c2[row1][0]) + " " + str(c2[row1][1]) + ","
                else:
                    str1 = str1 + str(c2[row1][0]) + " " + str(c2[row1][1]) + "," + str(c2[0][0]) + " " + str(
                        c2[0][1]) + "))"
            temp.append(str1)
            temp.append(name)
            poly.append(str1)
            poly1.append(temp)

        return poly, poly1

    def get_shp(self,a,b,type, spatial_ref):
        # spatial_ref = osr.SpatialReference()
        # spatial_ref.ImportFromEPSG(4326)
        data = a.values
        # print(data)
        print('a.columns')
        print(a.columns)
        colname = list(a.columns)
        print()
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

    def write_result(self,Data, out_root):
        file_write_obj = open(out_root, 'w')
        for key in Data:
            print(key)
            file_write_obj.writelines(key)
            file_write_obj.write('\n')
    def main(self):
        global type_dict, driver

        type_dict = {1: ogr.wkbPoint, 2: ogr.wkbLineString, 3: ogr.wkbPolygon, 4: ogr.wkbMultiPoint,
                     5: ogr.wkbMultiLineString, 6: ogr.wkbMultiPolygon, 7: ogr.wkbGeometryCollection}
        driver = ogr.GetDriverByName('ESRI Shapefile')
        gdal.SetConfigOption("SHAPE_ENCODING", "")

        filename_v = self.input_shpfile
        dataSource1 = driver.Open(filename_v, 0)
        if dataSource1 is None:
            print('文件打不开！')
        layer1 = dataSource1.GetLayer(0)
        a = self.get_fields(layer1)
        # print(a)
        type1 = layer1.GetGeomType()
        ref1 = layer1.GetSpatialRef()
        d = self.get_convex(a)

       


i