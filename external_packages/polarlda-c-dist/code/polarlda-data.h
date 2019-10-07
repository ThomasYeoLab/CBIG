// (C) Copyright 2018, Nanbo Sun (nanbosun@u.nus.edu), Xiuming Zhang

// This file is part of POLARLDA-C.

// POLARLDA-C is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your
// option) any later version.

// POLARLDA-C is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.

// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA

#ifndef POLARLDA_DATA_H
#define POLARLDA_DATA_H

#include <stdio.h>
#include <stdlib.h>

#include "polarlda.h"

#define OFFSET 0;                  // offset for reading data

corpus* read_data(char* data_filename);
int max_corpus_length(corpus* c);

#endif
