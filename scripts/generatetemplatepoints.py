
import argparse
import numpy

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-N', '--observations', type = int, default = 100, help = 'No. observations')
    parser.add_argument('--minx', type = float, default = -1.0, help = 'Min x')
    parser.add_argument('--miny', type = float, default = -1.0, help = 'Min y')
    parser.add_argument('--maxx', type = float, default = 1.0, help = 'Max x')
    parser.add_argument('--maxy', type = float, default = 1.0, help = 'Max y')

    parser.add_argument('-o', '--output', type = str, required = True, help = 'Output file')

    args = parser.parse_args()

    f = open(args.output, 'w')

    f.write('%15.9f %15.9f\n%15.9f %15.9f\n' % (args.minx, args.maxy, args.miny, args.maxy))

    f.write('%d\n' % args.observations)

    for i in range(args.observations):

        x = numpy.random.uniform() * (args.maxx - args.minx) + args.minx
        y = numpy.random.uniform() * (args.maxy - args.miny) + args.miny

        f.write('%15.9f %15.9f 0.0 1.0\n' % (x, y))

    f.close()
    

    
    
