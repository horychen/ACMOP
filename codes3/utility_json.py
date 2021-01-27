import json
import jsonpickle


def is_jsonable(x):
    try:
        json.dumps(x)
        return True
    except:
        return False

def build_json_fname(fname, suffix):
    if fname[:2] != '__':
        return '__'+fname+suffix+'.json'
    else:
        return fname

def to_json(obj, fname, suffix='', save_here=''):

    # remove not json serilizable content
    keys_to_delete = []
    for key, value in obj.__dict__.items():
        if not is_jsonable(value):
            keys_to_delete.append(key)
            obj.__dict__[key] = None
    print('keys_to_delete:', keys_to_delete)

    # write to file
    with open(save_here+build_json_fname(fname, suffix), 'w') as f:
        json.dump(obj.__dict__, f, indent=4)        

def recursive_dictionarify(obj, keys_not_jsonable):
    if type(obj) == type({}):
        the_big_dict = obj
    else:
        the_big_dict = obj.__dict__

    for key1, value1 in the_big_dict.items():

        if not is_jsonable(the_big_dict[key1]):

            print('DEBUG', value1)

            if type(value1) == type({}):
                the_big_dict[key1], keys_not_jsonable = recursive_dictionarify(value1, keys_not_jsonable)
            else:
                the_big_dict[key1] = value1.__dict__
            keys_not_jsonable.append(obj.__class__.__name__ + '-' + key1)

            if type(obj) == type({}):
                obj = the_big_dict
            else:
                obj.__dict__.update(the_big_dict)
            return obj, keys_not_jsonable

def to_json_recursively(obj, fname, suffix=''):
    s = json.dumps( json.loads( jsonpickle.encode(obj)), indent=4 )
    print('[to_json_recursively]', type( s ))
    # print( s )

    print('[to_json_recursively]', obj.spec_geometry_dict['DriveW_CurrentAmp'])

    with open(build_json_fname(fname, suffix), 'w') as f:
        f.write( s )

    # Obsolete
    # keys_not_jsonable = []
    # obj, keys_not_jsonable = recursive_dictionarify(obj, keys_not_jsonable)
    # with open(build_json_fname(fname, suffix), 'w') as f:
    #     json.dump(obj.__dict__, f, indent=4)

def from_json_recursively(fname, suffix=''):
    with open(build_json_fname(fname, suffix), 'r') as f:
        read_obj = jsonpickle.decode(f.read())
    return read_obj

def from_json(obj, fname, suffix=''):
    with open(build_json_fname(fname, suffix), 'r') as f:
        d = json.load(f)
        obj.__dict__.update(d)

class tempClass(object):
    """docstring for tempClass"""
    def __init__(self, arg):
        super(tempClass, self).__init__()
        self.arg = arg
        self.digit = 100

if __name__ == '__main__':
    import varname

    C = tempClass('C-string')
    B = tempClass(C)
    A = tempClass(B)
    # from varname import nameof
    # print( nameof(A) )
    # quit()

    # Option 1
    if False:
        count = 0
        print(A.__dict__, '\n\n\n')
        keys_not_jsonable = []
        while not is_jsonable(A.__dict__):
            if count>10:
                break
            count += 1
            A, keys_not_jsonable = recursive_dictionarify(A, keys_not_jsonable)
            print(count, A.__dict__)
            print(is_jsonable(A.__dict__))
        print(keys_not_jsonable)

    # Option 2: jsonpickle, see https://stackoverflow.com/questions/3768895/how-to-make-a-class-json-serializable
    import json
    import jsonpickle
    print(varname.nameof(A), ':')
    print(json.dumps(json.loads(jsonpickle.encode(A)), indent=4))

    # obj = jsonpickle.decode(file.read())
    obj = jsonpickle.decode(jsonpickle.encode(A))
    print(dir(obj))
    print(obj.arg)
    print(obj.arg.arg)
    print(obj.arg.arg.arg)

    # # Option 3: jsonS, see
    # import jsons
    # a_dict = jsons.dump(A, indent=4)
    # print(a_dict)

    # a_str = jsons.dumps(A, indent=4)
    # print(a_str)
    

    # Option 4: No third-party codes
        # ust add to_json method to your class like this:

        # def to_json(self):
        #   return self.message # or how you want it to be serialized
        # And add this code (from this answer), to somewhere at the top of everything:

        # from json import JSONEncoder

        # def _default(self, obj):
        #     return getattr(obj.__class__, "to_json", _default.default)(obj)

        # _default.default = JSONEncoder().default
        # JSONEncoder.default = _default
        # This will monkey-patch json module when it's imported so JSONEncoder.default() automatically checks for a special "to_json()" method and uses it to encode the object if found.

        # Just like Onur said, but this time you don't have to update every json.dumps() in your project.

        # share  edit  follow    
        # edited Apr 27 '18 at 2:13

        # Samuel Liew♦
        # 61.8k4040 gold badges127127 silver badges205205 bronze badges
        # answered Aug 4 '16 at 10:27

        # Fancy John
        # 33.4k22 gold badges2020 silver badges2424 bronze badges
        # 6
        # Big thanks! This is the only answer that allows me to do what I want: be able to serialize an object without changing the existing code. The other methods mostly do not work for me. The object is defined in a third-party library, and the serialization code is third-party too. Changing them will be awkward. With your method, I only need to do TheObject.to_json = my_serializer. – Yongwei Wu Oct 11 '17 at 13:12 


