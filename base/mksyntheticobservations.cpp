
#include <stdio.h>
#include <stdlib.h>

#include <getopt.h>

#include "genericinterface.hpp"
#include "tdtwave2dimage.hpp"
#include "tdtwave2dutil.hpp"

#include "rng.hpp"

static char short_options[] = "i:I:o:O:N:s:h";

static struct option long_options[] = {
  {"input-image", required_argument, 0, 'i'},
  {"input-points", required_argument, 0, 'I'},

  {"output", required_argument, 0, 'o'},
  {"output-true", required_argument, 0, 'O'},

  {"noise", required_argument, 0, 'N'},
  {"seed", required_argument, 0, 's'},
  
  {"help", no_argument, 0, 'h'},
  {0, 0, 0, 0}
};



static void usage(const char *pname);

  
int main(int argc, char *argv[])
{
  int c;
  int option_index;

  char *input_image;
  char *input_points;

  char *output_file;
  char *output_true;

  double noise;

  int seed;

  //
  // Defaults
  //
  input_image = nullptr;
  input_points = nullptr;

  output_file = nullptr;
  output_true = nullptr;

  seed = 983;

  noise = 1.0;
  
  
  while (true) {

    c = getopt_long(argc, argv, short_options, long_options, &option_index);
    if (c == -1) {
      break;
    }

    switch (c) {

    case 'i':
      input_image = optarg;
      break;

    case 'I':
      input_points = optarg;
      break;

    case 'o':
      output_file = optarg;
      break;

    case 'O':
      output_true = optarg;
      break;

    case 'N':
      noise = atof(optarg);
      if (noise <= 0.0) {
	fprintf(stderr, "error: noise must be positive\n");
	return -1;
      }
      break;

    case 's':
      seed = atoi(optarg);
      break;

    case 'h':
    default:
      usage(argv[0]);
      return -1;
    }
  }

  if (input_image == nullptr) {
    fprintf(stderr, "error: required input image parameter missing\n");
    return -1;
  }

  if (input_points == nullptr) {
    fprintf(stderr, "error: required input points parameter missing\n");
    return -1;
  }

  if (output_file == nullptr) {
    fprintf(stderr, "error: required output file parameter missing\n");
    return -1;
  }
      

  tdtwave2dimage image;
  if (!image.load(input_image)) {
    fprintf(stderr, "error: failed to load image file\n");
    return -1;
  }


  return 0;
}


static void usage(const char *pname)
{
  fprintf(stderr,
          "usage: %s [options]\n"
          "where options is one or more of:\n"
          "\n"
	  " -i|--input-image <filename>        Input true image\n"
	  " -I|--input-points <filename>         Input flight path\n"
	  "\n"
	  " -o|--output <filename>             Output file to write (required)\n"
	  " -O|--output-true <filename>        Output true observations to file\n"
	  "\n"
	  " -n|--noise <float>                 Std dev of Gaussian noise\n"
	  " -s|--seed <int>                    Random seed\n"
	  "\n"
	  " -h|--help                          Usage information\n"
	  "\n",
	  pname);
}

  
