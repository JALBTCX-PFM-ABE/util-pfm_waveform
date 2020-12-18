
/*********************************************************************************************

    This is public domain software that was developed by or for the U.S. Naval Oceanographic
    Office and/or the U.S. Army Corps of Engineers.

    This is a work of the U.S. Government. In accordance with 17 USC 105, copyright protection
    is not available for any work of the U.S. Government.

    Neither the United States Government, nor any employees of the United States Government,
    nor the author, makes any warranty, express or implied, without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE, or assumes any liability or
    responsibility for the accuracy, completeness, or usefulness of any information,
    apparatus, product, or process disclosed, or represents that its use would not infringe
    privately-owned rights. Reference herein to any specific commercial products, process,
    or service by trade name, trademark, manufacturer, or otherwise, does not necessarily
    constitute or imply its endorsement, recommendation, or favoring by the United States
    Government. The views and opinions of authors expressed herein do not necessarily state
    or reflect those of the United States Government, and shall not be used for advertising
    or product endorsement purposes.

*********************************************************************************************/

#include "pfm_waveform.h"

#include "version.h"


/*****************************************************************************

    Program:    pfm_waveform

    Purpose:    Extracts and analyzes waveforms for HOF points in an area of the
                input PFM.

    Programmer: Jan C. Depner

    Date:       11/03/10

*****************************************************************************/


OPTIONS                options;


void usage ()
{
  fprintf (stderr, "\nUsage: pfm_waveform PFM_FILE AREA_FILE\n");
  fprintf (stderr, "\nWhere:\n\n");
  fprintf (stderr, "\tPFM_FILE = PFM file (required)\n\n");
  fprintf (stderr, "\tAREA_FILE = Area file (required)\n\n");
  fprintf (stderr, "\t\tThe area file names must have a .ARE extension\n");
  fprintf (stderr, "\t\tfor ISS60 type area files, a .are extension for generic area files, or\n");
  fprintf (stderr, "\t\ta .afs extension for Army Corps area files.\n\n");
  fflush (stderr);
}



int32_t get_waveforms (char *path, double *polygon_x, double *polygon_y, int32_t polygon_count,
                        int32_t start, int32_t end, int32_t total, FILE *txt_fp, char *areafile __attribute__ ((unused)))
{
  FILE                   *data_fp, *wave_fp;
  POS_OUTPUT_T           pos;
  WAVE_HEADER_T          wave_header;
  WAVE_DATA_T            wave_data;
  int64_t                new_stamp, data_timestamp;
  int32_t                i;
  uint8_t                good_rec;
  static HYDRO_OUTPUT_T  hof;
  static int32_t         count = 0, percent = 0, old_percent = -1, good_count = 0, bad_count = 0;
  char                   wave_file[512], pos_file[512];
  static FILE            *pos_fp = NULL;
  /*
  static int32_t         sum_count = 0;
  static double          sum = 0.0;
  */

  uint8_t process_waveforms (HYDRO_OUTPUT_T *, WAVE_HEADER_T *, WAVE_DATA_T *, FILE *, int32_t);

  strcpy (wave_file, path);
  strcpy (&wave_file[strlen (wave_file) - 4], ".inh");
  fprintf(stderr,"%s %s %d %s\n",__FILE__,__FUNCTION__,__LINE__,path);

  wave_fp = open_wave_file (wave_file);

  if (wave_fp == NULL) return (0);

  wave_read_header (wave_fp, &wave_header);


  if ((data_fp = open_hof_file (path)) == NULL)
    {
      perror (path);
      exit (-1);
    }

  for (i = start ; i <= end ; i++)
    {
      /*  Find the record based on the timestamp from the hof file.  */

      hof_read_record (data_fp, i, &hof);
      data_timestamp = hof.timestamp;


      if (hof.correct_depth != -998.0)
        {
          if (wave_read_record (wave_fp, i, &wave_data))
            {
              /*  Find, open, and read the sbet (reprocessed) or pos file information.  */

              if (pos_fp != NULL)
                {
                  fclose (pos_fp);
                  pos_fp = NULL;
                }

              if (!get_pos_file (path, pos_file))
                {
                  fprintf (stderr, "Unable to find pos/sbet file for hof file %s\n", path);
                  fflush (stderr);
                }
              else
                {
                  if ((pos_fp = open_pos_file (pos_file)) == NULL)
                    {
                      fprintf (stderr, "Unable to open pos/sbet file for hof file %s\n", path);
                      fflush (stderr);
                    }
                  else
                    {
                      /*  Get the attitude data for this shot.  */

                      new_stamp = pos_find_record (pos_fp, &pos, data_timestamp);

                      if (!new_stamp) 
                        {
                          fprintf (stderr, "\n\nUnable to get timestamp ");
                          fprintf (stderr, "%"PRId64, data_timestamp);
                          fprintf (stderr, " for pos/sbet file %s\n", pos_file);
                          fprintf (stderr, "This usually indicates that the above pos/sbet file is FUBAR or the name is incorrect!\n");
                          fprintf (stderr, "Make sure the file name conforms to the naming convention (_YYMMDD_NNNN.out or .pos) and\n");
                          fprintf (stderr, "check the start and end times of this file (dump_pos) against the data in the HOF/TOF/IMG files.\n\n\n");

                          bad_count++;

                          if (bad_count > 100) exit (-1);
                        }
                      else
                        {
                          good_rec = NVFalse;
                          /*
                          if (hof.abdc > 70 && hof.correct_sec_depth != -998.0 && hof.kgps_sec_elev > 0.0)
                            {
                              fprintf (stderr, "%s %s %d %d %d %f %f %f\n", __FILE__, __FUNCTION__, __LINE__, hof.bot_bin_first, hof.bot_bin_second,
                                       hof.kgps_res_elev, hof.kgps_sec_elev, (hof.kgps_res_elev - hof.kgps_sec_elev) /
                                       (double) (hof.bot_bin_first - hof.bot_bin_second));
                              sum += (hof.kgps_res_elev - hof.kgps_sec_elev) / (double) (hof.bot_bin_first - hof.bot_bin_second);
                              sum_count++;
                            }
                          */

                          if (hof.abdc > 70 || (hof.correct_sec_depth != -998.0 && hof.sec_abdc > 70))
                            {
                              /*  Assume GCS was right if it picked two returns.  */

                              if (hof.correct_depth == -998.0 || hof.correct_sec_depth == -998.0)
                                {
                                  /*  Make sure we're inside the area we specified.  */

                                  if (inside_polygon2 (polygon_x, polygon_y, polygon_count, hof.longitude, hof.latitude))
                                    {
                                      /*  Get the altitude and other stuff.  */

                                      /*fprintf (txt_fp, "%.11f, %.11f, %f, %f, %d, %f\n", hof.latitude, hof.longitude, hof.kgps_res_elev, pos.altitude,
                                        hof.bot_bin_first, (pos.altitude - hof.kgps_res_elev) / (double) hof.bot_bin_first);*/

                                      good_rec = NVTrue;
                                      good_count++;
                                    }
                                }
                            }

                          /*55213 & 26990 & 26706*/
                          if (good_rec && (i == 54150)/* || i == 26990)*/) process_waveforms (&hof, &wave_header, &wave_data, txt_fp, i);
                        }
                    }
                }
            }
        }

      count++;
      percent = ((float) count / (float) total) * 100.0;
      if (percent != old_percent)
        {
          old_percent = percent;
          fprintf (stderr, "%03d%% processed            \r", percent);
          fflush (stderr);
        }
    }


  fclose (wave_fp);
  fclose (data_fp);

  /* 0.106795
  fprintf(stderr,"%s %s %d %d %f\n",__FILE__,__FUNCTION__,__LINE__,sum_count, sum / (double) sum_count);*/

  return (good_count);
}



int32_t main (int32_t argc, char **argv)
{
  FILE                   *txt_fp;
  int32_t                i, j, k, polygon_count = 0, x_start, y_start, width, height, total, percent = 0, old_percent = -1, recnum, icount = 0;
  int16_t                type;
  NV_I32_COORD2          coord;
  char                   pfm_file[512], areafile[512], txt_file[512], path[512];
  char                   c;
  double                 polygon_x[200], polygon_y[200];
  LIST_NUM               list[10000];
  BIN_RECORD             bin;
  DEPTH_RECORD           *depth;
  extern char            *optarg;
  extern int             optind;


  printf ("\n\n %s \n\n\n", VERSION);


  while ((c = getopt (argc, argv, "n")) != EOF)
    {
      switch (c)
        {
        case 'n':
          /*  Placeholder  */
          break;

        default:
          usage ();
          exit (-1);
          break;
        }
    }


  /* Make sure we got the mandatory file name arguments.  */

  if (optind >= argc)
    {
      usage ();
      exit (-1);
    }


  strcpy (pfm_file, argv[optind]);
  strcpy (areafile, argv[optind + 1]);
 

  get_area_mbr (areafile, &polygon_count, polygon_x, polygon_y, &options.mbr);

  strcpy (options.open_args.list_path, argv[optind]);

  options.open_args.checkpoint = 0;
  options.pfm_handle = open_existing_pfm_file (&options.open_args);

  if (options.pfm_handle < 0) pfm_error_exit (pfm_error);


  x_start = 0;
  y_start = 0;
  width = options.open_args.head.bin_width;
  height = options.open_args.head.bin_height;

  if (options.mbr.min_y > options.open_args.head.mbr.max_y || options.mbr.max_y < options.open_args.head.mbr.min_y ||
      options.mbr.min_x > options.open_args.head.mbr.max_x || options.mbr.max_x < options.open_args.head.mbr.min_x)
    {
      fprintf (stderr, "\n\nSpecified area is completely outside of the PFM bounds!\n\n");
      exit (-1);
    }


  /*  Match to nearest cell.  */

  x_start = NINT ((options.mbr.min_x - options.open_args.head.mbr.min_x) / options.open_args.head.x_bin_size_degrees);
  y_start = NINT ((options.mbr.min_y - options.open_args.head.mbr.min_y) / options.open_args.head.y_bin_size_degrees);
  width = NINT ((options.mbr.max_x - options.mbr.min_x) / options.open_args.head.x_bin_size_degrees);
  height = NINT ((options.mbr.max_y - options.mbr.min_y) / options.open_args.head.y_bin_size_degrees);


  /*  Adjust to PFM bounds if necessary.  */

  if (x_start < 0) x_start = 0;
  if (y_start < 0) y_start = 0;
  if (x_start + width > options.open_args.head.bin_width) width = options.open_args.head.bin_width - x_start;
  if (y_start + height > options.open_args.head.bin_height) height = options.open_args.head.bin_height - y_start;


  /*  Open the output file.  */

  strcpy (txt_file, pfm_basename (areafile));
  strcpy (&txt_file[strlen (txt_file) - 4], ".pts");

  if ((txt_fp = fopen (txt_file, "w")) == NULL) 
    {
      perror (txt_file);
      exit (-1);
    }


  /*  Set all hits to false.  */

  for (i = 0 ; i < 10000 ; i++)
    {
      list[i].hit = NVFalse;
      list[i].start = UINT32_MAX;
      list[i].end = 0;
    }


  /*  Double loop over height & width of area  */

  for (i = y_start ; i < y_start + height ; i++)
    {
      coord.y = i;
      for (j = x_start ; j < x_start + width ; j++)
        {
          coord.x = j;

          read_bin_record_index (options.pfm_handle, coord, &bin);

          if (bin.num_soundings)
            {
              /*  Get file numbers and min and max record numbers in file so we can figure out which   */
              /*  waveforms to retrieve.  */

              if (!read_depth_array_index (options.pfm_handle, coord, &depth, &recnum))
                {
                  for (k = 0 ; k < recnum ; k++)
                    {
                      if (!(depth[k].validity & PFM_DELETED))
                        {
                          int32_t m = depth[k].file_number;

                          if (m > 10000)
                            {
                              fprintf (stderr, "\n\nFile number out of bounds - %d\n\n", depth[k].file_number);
                              exit (-1);
                            }
                          list[m].hit = NVTrue;
                          if (depth[k].ping_number < list[m].start) list[m].start = depth[k].ping_number;
                          if (depth[k].ping_number > list[m].end) list[m].end = depth[k].ping_number;
                        }
                    }
                  free (depth);
                }
            }
        }

      percent = NINT (((float) (i - y_start) / (float) height) * 100.0);
      if (percent != old_percent)
        {
          fprintf (stderr, "%03d%% read\r", percent);
          fflush (stderr);
          old_percent = percent;
        }
    }


  total = 0;
  for (i = 0 ; i < 10000 ; i++)
    {
      if (list[i].hit) total += ((list[i].end - list[i].start) + 1);
    }


  for (i = 0 ; i < 10000 ; i++)
    {
      if (list[i].hit)
        {
          read_list_file (options.pfm_handle, (int16_t) i, path, &type);


          /*  Check for HOF data type.  */

          if (type == PFM_CHARTS_HOF_DATA)
            icount += get_waveforms (path, polygon_x, polygon_y, polygon_count, list[i].start, list[i].end, total, txt_fp, areafile);
        }
    }

  fprintf (stderr, "Extracted %d waveforms\n\n", icount);
  fflush (stderr);

  close_pfm_file (options.pfm_handle);
  fclose (txt_fp);


  return (0);
}
