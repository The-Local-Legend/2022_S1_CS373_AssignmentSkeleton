from cProfile import label
import math
import sys
from pathlib import Path
import statistics
from matplotlib import pyplot
from matplotlib.patches import Rectangle

# import our basic, light-weight png reader library
import imageIO.png
class Queue:
    def __init__(self):
        self.items = []

    def isEmpty(self):
        return self.items == []

    def enqueue(self, item):
        self.items.insert(0,item)

    def dequeue(self):
        return self.items.pop()

    def size(self):
        return len(self.items)
# this function reads an RGB color png file and returns width, height, as well as pixel arrays for r,g,b
def readRGBImageToSeparatePixelArrays(input_filename):

    image_reader = imageIO.png.Reader(filename=input_filename)
    # png reader gives us width and height, as well as RGB data in image_rows (a list of rows of RGB triplets)
    (image_width, image_height, rgb_image_rows, rgb_image_info) = image_reader.read()

    print("read image width={}, height={}".format(image_width, image_height))

    # our pixel arrays are lists of lists, where each inner list stores one row of greyscale pixels
    pixel_array_r = []
    pixel_array_g = []
    pixel_array_b = []

    for row in rgb_image_rows:
        pixel_row_r = []
        pixel_row_g = []
        pixel_row_b = []
        r = 0
        g = 0
        b = 0
        for elem in range(len(row)):
            # RGB triplets are stored consecutively in image_rows
            if elem % 3 == 0:
                r = row[elem]
            elif elem % 3 == 1:
                g = row[elem]
            else:
                b = row[elem]
                pixel_row_r.append(r)
                pixel_row_g.append(g)
                pixel_row_b.append(b)

        pixel_array_r.append(pixel_row_r)
        pixel_array_g.append(pixel_row_g)
        pixel_array_b.append(pixel_row_b)

    return (image_width, image_height, pixel_array_r, pixel_array_g, pixel_array_b)


# a useful shortcut method to create a list of lists based array representation for an image, initialized with a value
def createInitializedGreyscalePixelArray(image_width, image_height, initValue = 0):

    new_array = [[initValue for x in range(image_width)] for y in range(image_height)]
    return new_array

def computeRGBToGreyscale(pixel_array_r, pixel_array_g, pixel_array_b, image_width, image_height):
    
    greyscale_pixel_array = createInitializedGreyscalePixelArray(image_width, image_height)
    
    for i in range(len(greyscale_pixel_array)):
        for j in range(len(greyscale_pixel_array[i])):
            greyscale_pixel_array[i][j] =  round((0.299 * pixel_array_r[i][j]) + (0.587 * pixel_array_g[i][j]) + (0.114 * pixel_array_b[i][j]))
    
    return greyscale_pixel_array
def contrastStretch(pix_array, image_width, image_height):
    copyarray = pix_array[:]
    mymin = min([min(r) for r in copyarray]) #lowest pixel in image 
    mymax = max([max(r) for r in copyarray]) #highest pixel in image
    print(mymin)
    print(mymax)
    #contrast = 255 / (mymax - mymin)
    #print(contrast)
    for i in range(len(copyarray)):
        for j in range(len(copyarray[i])):
            copyarray[i][j] = round((copyarray[i][j] - mymin) * (255 / (mymax - mymin))) # the math thingy
            if copyarray[i][j] > 255:
                copyarray[i][j] = 255
    return copyarray 

def standardDev(pix_array, image_width, image_height):
    copyarray = createInitializedGreyscalePixelArray(image_width, image_height)
    for i in range(2, image_height - 2): #loop through rows without border
        for j in range(2, image_width - 2): #loop through columns without border
            neighbourhood = [] # array for the 5x5
            for row in range(i - 2, i + 3): #range of -2 + 2 around the current center pixel
                for col in range(j - 2, j + 3): #range of -2 + 2 around the current center pixel
                    neighbourhood.append(pix_array[row][col]) #add to neighbourhood
            standev = statistics.pstdev(neighbourhood) #standard deviation
            #avrg = statistics.mean(neighbourhood) 
            copyarray[i][j] = standev
    print(neighbourhood)
    return copyarray

def binaryimage(pix_array, image_width, image_height):
    copyarray = createInitializedGreyscalePixelArray(image_width, image_height)
    for i in range(image_height):
        for j in range(image_width):
            if pix_array[i][j] > 150:
                copyarray[i][j] = 255
            else:
                copyarray[i][j] = 0
    return copyarray

def erosion(pix_array, image_width, image_height):
    copyarray = createInitializedGreyscalePixelArray(image_width, image_height)
    for i in range(image_height):
        for j in range(image_width):
            try:    
                count = 0
                if pix_array[i][j] > 1:
                    for k in range(i - 1, i + 2):
                        for l in range(j - 1, j + 2):
                            if pix_array[k][l] > 1:
                                count += 1
                    if count == 9:
                        copyarray[i][j] = 255
                    else:
                        copyarray[i][j] = 0
                else:
                    copyarray[i][j] = 0
            except(IndexError):
                continue
    return copyarray

def dilation(pix_array, image_width, image_height):
    copyarray = createInitializedGreyscalePixelArray(image_width, image_height)
    for i in range(image_height):
        for j in range(image_width):
            try:
                if pix_array[i][j] > 1:
                    for k in range(i - 1, i + 2):
                        for l in range(j - 1, j + 2):
                            copyarray[k][l] = 255
            except(IndexError):
                continue
    return copyarray

def connections(pix_array, image_width, image_height):
    #debug = open("debug.txt", "w")
    current_label = 1
    copyarray = createInitializedGreyscalePixelArray(image_width, image_height) #Empty array of 0s
    visitedarray = createInitializedGreyscalePixelArray(image_width, image_height)
    for i in range(image_height):
        for j in range(image_width):
            if pix_array[i][j] > 1 and visitedarray[i][j] != True:
                thequeue = Queue()
                thequeue.enqueue((i, j))
                #debug.write(str(thequeue.size()))
                while not thequeue.isEmpty():
                    primeindexes = thequeue.dequeue()
                    #debug.writelines(str(primeindexes))
                    copyarray[primeindexes[0]][primeindexes[1]] = current_label
                    try: #Left pixel
                        if visitedarray[primeindexes[0]][primeindexes[1] - 1] != True and pix_array[primeindexes[0]][primeindexes[1] - 1] > 1:
                            visitedarray[primeindexes[0]][primeindexes[1] - 1] = True
                            thequeue.enqueue((primeindexes[0], primeindexes[1] - 1))
                            
                        else:
                            visitedarray[primeindexes[0]][primeindexes[1] - 1] = True
                    except(IndexError):
                        pass
                    try: #Right pixel
                        if visitedarray[primeindexes[0]][primeindexes[1] + 1] != True and pix_array[primeindexes[0]][primeindexes[1] + 1] > 1:
                            visitedarray[primeindexes[0]][primeindexes[1] + 1] = True
                            thequeue.enqueue((primeindexes[0], primeindexes[1] + 1))
                        else:
                            visitedarray[primeindexes[0]][primeindexes[1] + 1] = True
                    except(IndexError):
                        pass
                    try: #Upper pixel
                        if visitedarray[primeindexes[0] + 1][primeindexes[1]] != True and pix_array[primeindexes[0] + 1][primeindexes[1]] > 1:
                            visitedarray[primeindexes[0] + 1][primeindexes[1]] = True
                            thequeue.enqueue((primeindexes[0] + 1, primeindexes[1]))
                        else:
                            visitedarray[primeindexes[0] + 1][primeindexes[1]] = True
                    except(IndexError):
                        pass
                    try:
                        if visitedarray[primeindexes[0] - 1][primeindexes[1]] != True and pix_array[primeindexes[0] - 1][primeindexes[1]] > 1:
                            visitedarray[primeindexes[0] - 1][primeindexes[1]]
                            thequeue.enqueue((primeindexes[0] - 1, primeindexes[1]))
                        else:
                            visitedarray[primeindexes[0] - 1][primeindexes[1]] = True
                    except(IndexError):
                        pass
                    #debug.writelines(str(thequeue.items))
                current_label += 1
    dictionary = dict()
    for i in range(image_height):
        for j in range(image_width):
            if copyarray[i][j] > 0 and copyarray[i][j] not in dictionary.keys():
                dictionary[copyarray[i][j]] = 1
            elif copyarray[i][j] > 0:
                dictionary[copyarray[i][j]] += 1
    region = 0
    highest = 0
    for i in dictionary:
        if dictionary[i] > highest:
            region = i
            highest = dictionary[i]
    #file = open("connect.txt", "w")
    #file.write(str(copyarray))
    #file.close()
    #debug.close()
    print(dictionary)
    #print(max(dictionary.values()))
    print(region)
    tl = (image_width, image_height)
    br = (0, 0)
    for i in range(image_height):
        for j in range(image_width):
            if copyarray[i][j] == region:
                if j < tl[0]:
                    tl = (j, tl[1])
                    if i < tl[1]:
                        tl = (tl[0], i)
                if j > br[0]:
                    br = (j, br[1])
                    if i > br[1]:
                        br = (br[0], i)
    print(tl, br)
    return (tl, br)

# This is our code skeleton that performs the license plate detection.
# Feel free to try it on your own images of cars, but keep in mind that with our algorithm developed in this lecture,
# we won't detect arbitrary or difficult to detect license plates!
def main():

    command_line_arguments = sys.argv[1:]

    SHOW_DEBUG_FIGURES = True

    # this is the default input image filename
    input_filename = "numberplate1.png"

    if command_line_arguments != []:
        input_filename = command_line_arguments[0]
        SHOW_DEBUG_FIGURES = False

    output_path = Path("output_images")
    if not output_path.exists():
        # create output directory
        output_path.mkdir(parents=True, exist_ok=True)

    output_filename = output_path / Path(input_filename.replace(".png", "_output.png"))
    if len(command_line_arguments) == 2:
        output_filename = Path(command_line_arguments[1])


    # we read in the png file, and receive three pixel arrays for red, green and blue components, respectively
    # each pixel array contains 8 bit integer values between 0 and 255 encoding the color values
    (image_width, image_height, px_array_r, px_array_g, px_array_b) = readRGBImageToSeparatePixelArrays(input_filename)

    # setup the plots for intermediate results in a figure
    fig1, axs1 = pyplot.subplots(2, 2)
    #axs1[0, 0].set_title('Input red channel of image')
    #axs1[0, 0].imshow(px_array_r, cmap='gray')
    #axs1[0, 1].set_title('Input green channel of image')
    #axs1[0, 1].imshow(px_array_g, cmap='gray')
    #axs1[1, 0].set_title('Input blue channel of image')
    #axs1[1, 0].imshow(px_array_b, cmap='gray')


    # STUDENT IMPLEMENTATION here

    px_array = computeRGBToGreyscale(px_array_r, px_array_b, px_array_g, image_width, image_height)
    px_arraycopy = px_array[:]
    contrastarray = contrastStretch(px_arraycopy, image_width, image_height)
    
    standarddevarray = standardDev(contrastarray, image_width, image_height)
    standarddevarray = contrastStretch(standarddevarray, image_width, image_height)
    thresholdarray = binaryimage(standarddevarray, image_width, image_height)
    dilationarray = dilation(thresholdarray, image_width, image_height)
    for i in range(3):
        dilationarray = dilation(dilationarray, image_width, image_height)
    erosionarray = erosion(dilationarray, image_width, image_height)
    for i in range(2):
        erosionarray = erosion(erosionarray, image_width, image_height)
    tl, br = connections(erosionarray, image_width, image_height)
    axs1[0, 1].set_title('Standard Deviation')
    axs1[0, 1].imshow(standarddevarray, cmap='gray')
    axs1[1, 0].set_title('Erosion')
    axs1[1, 0].imshow(erosionarray, cmap='gray')
    axs1[0, 0].set_title('Thresholding')
    axs1[0, 0].imshow(thresholdarray, cmap='gray')
    # compute a dummy bounding box centered in the middle of the input image, and with as size of half of width and height
    center_x = image_width / 2.0
    center_y = image_height / 2.0
    bbox_min_x = center_x - image_width / 4.0
    bbox_max_x = center_x + image_width / 4.0
    bbox_min_y = center_y - image_height / 4.0
    bbox_max_y = center_y + image_height / 4.0





    # Draw a bounding box as a rectangle into the input image
    axs1[1, 1].set_title('Final image of detection')
    axs1[1, 1].imshow(px_array, cmap='gray')
    rect = Rectangle((bbox_min_x, bbox_min_y), bbox_max_x - bbox_min_x, bbox_max_y - bbox_min_y, linewidth=1,
                     edgecolor='g', facecolor='none')
    axs1[1, 1].add_patch(rect)



    # write the output image into output_filename, using the matplotlib savefig method
    extent = axs1[1, 1].get_window_extent().transformed(fig1.dpi_scale_trans.inverted())
    pyplot.savefig(output_filename, bbox_inches=extent, dpi=600)

    if SHOW_DEBUG_FIGURES:
        # plot the current figure
        pyplot.show()


if __name__ == "__main__":
    main()