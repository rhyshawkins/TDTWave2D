
import sys
import argparse
import numpy

def model_constant(amplitude, scale, nx, ny):

    return amplitude


MODELS = {'Constant' : model_constant}

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-W', '--width', type = int, default = 32, help = 'Image width')
    parser.add_argument('-H', '--height', type = int, default = 32, help = 'Image height')

    parser.add_argument('-o', '--output', type = str, required = True, help = 'Output file')

    parser.add_argument('-A', '--amplitude', type = float, default = 0.0, help = 'Amplitude')
    parser.add_argument('-S', '--scale', type = float, default = 1.0, help = 'Scale')
    
    parser.add_argument('-m', '--model', type = str, default = 'Constant', help = 'Model name')
                        
    args = parser.parse_args()

    if (args.width <= 0) or ((args.width & (args.width - 1)) != 0):
        print 'Width must be positive power of 2'
        sys.exit(-1)
        
    if (args.height <= 0) or ((args.height & (args.height - 1)) != 0):
        print 'Height must be positive power of 2'
        sys.exit(-1)

    image = numpy.zeros((args.height, args.width))

    
    if not MODELS.has_key(args.model):
        print 'Unknown model %s, select from:' % args.model
        print MODELS.keys
        sys.exit(-1)

    func = MODELS[args.model]

    f = open(args.output, 'w')

    f.write('%d %d\n' % (args.height, args.width))
    
    for j in range(args.height):

        ny = -1.0 + 2.0 * (float(j) + 0.5)/float(args.height)

        for i in range(args.width):

            
            nx = -1.0 + 2.0 * (float(i) + 0.5)/float(args.width)

            f.write('%15.9f ' % func(args.amplitude, args.scale, nx, ny))

        f.write('\n')
            
    f.close()
    

    
    
