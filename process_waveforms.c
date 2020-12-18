
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

#define DEBUG

/***************************************************************************\
*                                                                           *
*   Module Name:        process_waveforms                                   *
*                                                                           *
*   Programmer(s):      Jan C. Depner                                       *
*                                                                           *
*   Date Written:       November 2010                                       *
*                                                                           *
*   Purpose:            Process the APD and PMT waveforms to try to detect  *
*                       small objects near the selected return(s).          *
*                                                                           *
*   Arguments:          hof_record     - the HOF record                     *
*                       wave_header    - wave file header                   *
*                       wave_data      - the wave file header               *
*                                                                           *
\***************************************************************************/

uint8_t process_waveforms (HYDRO_OUTPUT_T *hof, WAVE_HEADER_T *wave_header, WAVE_DATA_T *wave_data, FILE *txt_fp __attribute__ ((unused)),
                           int32_t rec __attribute__ ((unused)))
{
  int32_t        i, j, k, rise, start_run, drop = 0, p_start_data[2], p_end_data[2], p_data_rise[2], p_data_run[2], p_count, first_drop,
                 threshold_count, a_start_data[2], a_end_data[2], a_data_rise[2], a_data_run[2], a_count, run_req, start_loc, end_loc;


  /*  Don't mess with shoreline depth swapped or shallow water algorithm data.  */

  if (hof->abdc == 72 || hof->sec_abdc == 72 || hof->abdc == 74 || hof->sec_abdc == 74) return (NVFalse);


  /*  Hard wired to 6 for the moment.  */

  run_req = 6;


  /*  Initialize the run and slope variables for the PMT data.  */

  rise = 0;
  drop = 0;
  start_run = 0;
  threshold_count = 0;
  start_loc = 0;
  end_loc = 0;
  p_start_data[0] = p_start_data[1] = 0;
  p_end_data[0] = p_end_data[1] = 0;
  p_data_rise[0] = p_data_rise[1] = 0;
  p_data_run[0] = p_data_run[1] = 0;
  p_count = 0;


  /*  Initialize the first_drop index.  We don't want to start searching for runs until we've cleared the surface  */
  /*  return.  This is not how Optech does it but I'm not really interested in very shallow water for this.  */

  first_drop = 0;


  /*  Loop through the PMT data looking for runs of increasing value that exceed "run_req".  Skip the first 20  */
  /*  bins so that we don't start looking in the noisy section prior to the surface return.  */

  for (i = 20 ; i < wave_header->pmt_size ; i++)
    {
      /*  If we get 10 points within 15 of the ac zero offset we're done.  */

      if (wave_data->pmt[i] - wave_header->ac_zero_offset[PMT] < 15)
        {
          threshold_count++;
          if (threshold_count > 10) break;
        }


      /*  If the value is increasing...   */

      if (first_drop && (wave_data->pmt[i] - wave_data->pmt[i - 1] > 0))
        {
          if (!start_loc) start_loc = i;


          /*  Increment the rise count.  */

          rise++;


          /*  If we have not already started a run and the rise count is greater than "run_req", start a new run.  */

          if (!start_run && rise > run_req) start_run = start_loc;


          /*  Set the drop count to 0.  */

          drop = 0;
        }


      /*  If the value is decreasing  */

      else if (wave_data->pmt[i] - wave_data->pmt[i - 1] <= 0)
        {
          if (!drop) end_loc = i;


          /*  Increment the drop counter.  */

          drop++;


          /*  If we have five consecutive drops...  */

          if (drop >= 5)
            {
              if (!first_drop)
                {
                  first_drop = i;
                }
              else
                {
                  /*  If we have a run going of more than "run_req" points (qualifying run) ...  */

                  if (start_run)
                    {
#ifdef DEBUG
                      fprintf(stderr,"%s %s %d %d %d\n",__FILE__,__FUNCTION__,__LINE__,start_run,i);
#endif
                      /*  Save the start, end, run, and rise values for this run.  */

                      p_start_data[p_count] = start_run;
                      p_end_data[p_count] = end_loc;
                      p_data_run[p_count] = end_loc - start_loc + 1;
                      p_data_rise[p_count] = wave_data->pmt[end_loc] - wave_data->pmt[start_run];


                      /*  Increment the qualifying run counter.  */

                      p_count++;


                      /*  If we have encountered 2 qualifying runs we can stop looking at the data.  */

                      if (p_count == 2) break;
                    }


                  /*  Zero out the rise count and start_run to get ready for the next qualifying run.  */

                  rise = 0;
                  start_run = 0;
                }
            }
          start_loc = 0;
        }
    }


  /*  If we incremented to p_count == 2 but never got another start_run we need to decrement p_count.  */

  if (p_count == 2 && !start_run) p_count--;



  /*  Hard wired to 6 for the moment.  */

  run_req = 6;


  /*  Initialize the run and slope variables for the APD data.  */

  rise = 0;
  drop = 0;
  start_run = 0;
  threshold_count = 0;
  start_loc = 0;
  end_loc = 0;
  a_start_data[0] = a_start_data[1] = 0;
  a_end_data[0] = a_end_data[1] = 0;
  a_data_rise[0] = a_data_rise[1] = 0;
  a_data_run[0] = a_data_run[1] = 0;
  a_count = 0;


  /*  Initialize the first_drop index.  We don't want to start searching for runs until we've cleared the surface  */
  /*  return.  This is not how Optech does it but I'm not really interested in very shallow water for this.  */

  first_drop = 0;


  /*  Loop through the APD data looking for runs of increasing value that exceed "run_req".  Skip the first 20  */
  /*  bins so that we don't start looking in the noisy section prior to the surface return.  */

  for (i = 20 ; i < wave_header->apd_size ; i++)
    {
      /*  If we get 20 points under the ac zero offset we're done.  */

      if (wave_data->apd[i] - wave_header->ac_zero_offset[APD] < 0)
        {
          threshold_count++;
          if (threshold_count > 20) break;
        }


      /*  If the value is increasing...  */

      if (first_drop && (wave_data->apd[i] - wave_data->apd[i - 1] > 0))
        {
          if (!start_loc) start_loc = i;


          /*  Increment the rise count.  */

          rise++;


          /*  If we have not already started a run and the rise count is greater than "run_req", start a new run.  */

          if (!start_run && rise > run_req) start_run = start_loc;


          /*  Set the drop count to 0.  */

          drop = 0;
        }


      /*  If the value is decreasing  */

      else if (wave_data->apd[i] - wave_data->apd[i - 1] <= 0)
        {
          if (!drop) end_loc = i;


          /*  Increment the drop counter.  */

          drop++;


          /*  If we have five consecutive drops...  */

          if (drop >= 5)
            {
              if (!first_drop)
                {
                  first_drop = i;
                }
              else
                {
                  /*  If we have a run going of more than "run_req" points (qualifying run) ...  */

                  if (start_run)
                    {
                      /*  Save the start, end, run, and rise values for this run.  */

                      a_start_data[a_count] = start_run;
                      a_end_data[a_count] = end_loc;
                      a_data_run[a_count] = end_loc - start_loc + 1;
                      a_data_rise[a_count] = wave_data->apd[end_loc] - wave_data->apd[start_run];


                      /*  Increment the qualifying rise counter.  */

                      a_count++;


                      /*  If we have encountered 2 qualifying runs we can stop looking at the data.  */

                      if (a_count == 2) break;
                    }


                  /*  Zero out the rise count and start_run to get ready for the next qualifying run.  */

                  rise = 0;
                  start_run = 0;
                }
            }
          start_loc = 0;
        }
    }


  /*  If we incremented to a_count == 2 but never got another start_run we need to decrement a_count.  */

  if (a_count == 2 && !start_run) a_count--;


#ifdef DEBUG
  fprintf(stderr,"%s %s %d %d %d %d\n",__FILE__,__FUNCTION__,__LINE__,p_count,a_count,hof->bot_bin_used_pmt);
#endif
  /*  Loop through the qualifying PMT runs.  */
  /*
  uint8_t        tagged;
  tagged = NVFalse;
  for (j = 0 ; j < p_count ; j++)
    {
      prev_slope = 0.0;
#ifdef DEBUG
      fprintf(stderr,"%s %s %d %d %d\n",__FILE__,__FUNCTION__,__LINE__,p_data_rise[j],p_data_run[j]);
#endif
      for (i = p_start_data[j] + 1 ; i <= p_end_data[j] ; i++)
        {
          if ((i - (p_start_data[j] + 1)) >= 3)
            {
              slope = (wave_data->pmt[i] - wave_data->pmt[i - 3]) / 3.0;

#ifdef DEBUG
              fprintf (stderr,"%s %s %d %d %d %d %d %f\n",__FILE__,__FUNCTION__,__LINE__,j,i,wave_data->pmt[i - 1], wave_data->pmt[i], slope);
#endif

              if (!tagged && (p_end_data[j] - i) > 10 && (i - p_start_data[j]) > 5 && (prev_slope - slope) > 2.1)
                {
                  fprintf (txt_fp, "%.11f,%.11f,%f\n", hof->latitude, hof->longitude, 10000000.0 + (float) rec);
                  tagged = NVTrue;
                }

              prev_slope = slope;  
            }
        }
    }


  for (j = 0 ; j < a_count ; j++)
    {
      prev_slope = 0.0;
#ifdef DEBUG
      fprintf(stderr,"%s %s %d %d %d\n",__FILE__,__FUNCTION__,__LINE__,a_data_rise[j],a_data_run[j]);
#endif
      for (i = a_start_data[j] + 1 ; i <= a_end_data[j] ; i++)
        {
          if ((i - (a_start_data[j] + 1)) >= 3)
            {
              slope = (wave_data->apd[i] - wave_data->apd[i - 3]) / 3.0;

#ifdef DEBUG
              fprintf (stderr,"%s %s %d %d %d %d %d %f\n",__FILE__,__FUNCTION__,__LINE__,j,i,wave_data->apd[i - 1], wave_data->apd[i], slope);
#endif

              if (!tagged && (a_end_data[j] - i) > 10 && (i - a_start_data[j]) > 5 && (prev_slope - slope) > 2.1)
                {
                  fprintf (txt_fp, "%.11f,%.11f,%f\n", hof->latitude, hof->longitude, 20000000.0 + (float) rec);
                  tagged = NVTrue;
                }

              prev_slope = slope;
            }
        }
    }

  tagged = NVFalse;

  */

  for (j = 0 ; j < p_count ; j++)
    {
      float *first_diff, *second_diff;
      int32_t size_first = p_end_data[j] - p_start_data[j];
      int32_t size_second = p_end_data[j] - p_start_data[j] - 1;
      first_diff = (float *) calloc (size_first, sizeof (float));
      second_diff = (float *) calloc (size_second, sizeof (float));

      for (i = p_start_data[j] + 1, k = 0 ; i <= p_end_data[j] ; i++, k++)
        {
          first_diff[k] = wave_data->pmt[i] - wave_data->pmt[i - 1];
          fprintf(stderr, "%s %s %d %d %d %d %f\n",__FILE__,__FUNCTION__,__LINE__,k,wave_data->pmt[i - 1],wave_data->pmt[i],first_diff[k]);
        }

      for (i = 1, k = 0 ; i < size_first ; i++, k++)
        {
          second_diff[k] = first_diff[i] - first_diff[i - 1];
          fprintf(stderr, "%s %s %d %d %f\n",__FILE__,__FUNCTION__,__LINE__,k,second_diff[k]);
        }

      free (first_diff);
      free (second_diff);
    }

  for (j = 0 ; j < a_count ; j++)
    {
      float *first_diff, *second_diff;
      int32_t size_first = a_end_data[j] - a_start_data[j];
      int32_t size_second = a_end_data[j] - a_start_data[j] - 1;
      first_diff = (float *) calloc (size_first, sizeof (float));
      second_diff = (float *) calloc (size_second, sizeof (float));

      for (i = a_start_data[j] + 1, k = 0 ; i <= a_end_data[j] ; i++, k++)
        {
          first_diff[k] = wave_data->apd[i] - wave_data->apd[i - 1];
          fprintf(stderr, "%s %s %d %d %d %d %f\n",__FILE__,__FUNCTION__,__LINE__,k,wave_data->apd[i - 1],wave_data->apd[i],first_diff[k]);
        }

      for (i = 1, k = 0 ; i < size_first ; i++, k++)
        {
          second_diff[k] = first_diff[i] - first_diff[i - 1];
          fprintf(stderr, "%s %s %d %d %f\n",__FILE__,__FUNCTION__,__LINE__,k,second_diff[k]);
        }

      free (first_diff);
      free (second_diff);
    }

  return (NVTrue);
}
