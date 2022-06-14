def gson = GsonTools.getInstance(true)
def json = new File("/path/to/file/cell_segmentation_geo.json").text

def type = new com.google.gson.reflect.TypeToken<List<qupath.lib.objects.PathObject>>() {}.getType()
def deserializedAnnotations = gson.fromJson(json, type)
addObjects(deserializedAnnotations)   
