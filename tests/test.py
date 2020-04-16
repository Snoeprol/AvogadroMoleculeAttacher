import os
import glob
path = os.path.curdir
#path = os.path.join('src','input','molecule')
path = os.path.join('input','molecule')
for filename in glob.glob(os.path.join(path, '*.xyz')):
    print(filename)
print(glob.glob(os.path.join(path, '*.xyz'))[0])
#print(os.listdir(os.getcwd()))
