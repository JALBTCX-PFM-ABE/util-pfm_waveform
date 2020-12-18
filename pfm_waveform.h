
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

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <proj_api.h>

/*
#include <gdal.h>
#include <cpl_string.h>
#include <ogr_api.h>
#include <ogr_srs_api.h>
*/

#include "FileHydroOutput.h"
#include "FilePOSOutput.h"
#include "FileWave.h"

#include "nvutility.h"

#include "pfm.h"
#include "pfm_extras.h"


typedef struct
{
  uint8_t       hit;
  uint32_t      start;
  uint32_t      end;
} LIST_NUM;


typedef struct
{
  double        altitude;
  PFM_OPEN_ARGS open_args;
  int32_t       pfm_handle;
  NV_F64_XYMBR  mbr;
} OPTIONS;
