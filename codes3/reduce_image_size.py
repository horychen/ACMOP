import os
from PIL import Image
import PIL
PIL.Image.MAX_IMAGE_PIXELS = 933120000

for file in os.listdir('./'):
    if 'Figure_selected_optimal_design' in file and 'png' in file:
        print(file)

        foo = Image.open(file)
        print(foo.size)

        foo = foo.resize( (int(foo.size[0]*0.01),int(foo.size[1]*0.01) ), Image.ANTIALIAS)
        foo.save('resized' + file, quality=300)

        # # The saved downsized image size is 24.8kb
        # foo.save("path\\to\\save\\image_scaled_opt.jpg",optimize=True,quality=95)

