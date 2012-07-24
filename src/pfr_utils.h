/*
 *  pfr_utils.h
 *  pfr
 *
 *  Created by Eric C. Anderson on 8/2/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

extern int XNumOfaType(analtype aType);
extern int NumAofaType(analtype aType);
extern int XdimOfaType(analtype aType);
extern int **SetXKeys(analtype aType);
extern int IdxFromPopName(char *Str, pfr_forback_data *D);
extern int Xstate(int a, analtype aType);