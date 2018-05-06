
import argparse
import numpy

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-N', '--observations', type = int, default = 100, help = 'No. observations')
    parser.add_argument('--minx', type = float, default = -1.0, help = 'Min x')
    parser.add_argument('--miny', type = float, default = -1.0, help = 'Min y')
    parser.add_argument('--maxx', type = float, default = 1.0, help = 'Max x')
    parser.add_argument('--maxy', type = float, default = 1.0, help = 'Max y')
    parser.add_argument('--samples', type = int, default = 10, help = 'Samples per path')

    parser.add_argument('-o', '--output', type = str, required = True, help = 'Output file')

    args = parser.parse_args()

    f = open(args.output, 'w')

    f.write('%15.9f %15.9f\n%15.9f %15.9f\n' % (args.minx, args.maxy, args.miny, args.maxy))

    f.write('%d\n' % args.observations)

    for i in range(args.observations):

        sx = numpy.random.uniform() * (args.maxx - args.minx) + args.minx
        sy = numpy.random.uniform() * (args.maxy - args.miny) + args.miny

        dx = numpy.random.uniform() * (args.maxx - args.minx) + args.minx
        dy = numpy.random.uniform() * (args.maxy - args.miny) + args.miny

        px = numpy.linspace(sx, dx, args.samples)
        py = numpy.linspace(sy, dy, args.samples)

        f.write('0.0 1.0 %d\n' % args.samples)
        for j in range(args.samples):
            f.write('%15.9f %15.9f\n' % (px[j], py[j]))

    f.close()
    

    
    
