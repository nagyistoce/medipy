def configure_file(source, destination, **kwargs) :
    data = open(source).read()
    
    for key, value in kwargs.items() :
        data = data.replace("@{0}@".format(key), value)
    
    open(destination, "w").write(data)

def get_mangled_type_name(type):
    dictionary = {
        "unsigned short" : "US",
        "unsigned char" : "UC"
    }
    
    return dictionary[type]