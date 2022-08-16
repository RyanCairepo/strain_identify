import sys,re
from PIL import Image

img = Image.open(sys.argv[1])
eps_name = re.sub('\.png',".eps",sys.argv[1])
print(eps_name)
#img.save(eps_name)